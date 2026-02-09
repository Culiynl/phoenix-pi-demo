import React, { useState, useMemo } from 'react';
import MolstarViewer from './MolstarViewer';
import FileUploadPanel from './FileUploadPanel';
import GeminiAdvantageShowcase from './GeminiAdvantageShowcase';
import {
    ResponsiveContainer, AreaChart, Area, XAxis, YAxis, Tooltip, CartesianGrid
} from 'recharts';
import { db } from '../firebase';
import { doc, onSnapshot } from "firebase/firestore";
import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';

const LiveWorkspace = ({ missionId, activeView = 'summary', setActiveView }) => {
    const [missionData, setMissionData] = useState(null);
    const [chartMetric, setChartMetric] = useState('pic50'); // pic50, mw, logp, qed, sa
    const [leaderboardView, setLeaderboardView] = useState('table'); // table, cards

    React.useEffect(() => {
        if (!missionId) return;
        const missionRef = doc(db, "missions", missionId);
        const unsubscribe = onSnapshot(missionRef, (docSnap) => {
            if (docSnap.exists()) {
                setMissionData(docSnap.data());
            }
        });
        return () => unsubscribe();
    }, [missionId]);

    // Multi-objective sorting logic (matches backend)
    const sortedCandidates = useMemo(() => {
        const cands = missionData?.chembl_summary?.top_candidates || [];
        const weights = missionData?.active_weights || { qed: 1.0, pIC50: 1.0, logp: 0, mw: 0.01, sascore: 0.5 };

        return [...cands].map(c => {
            if (!c) return null;
            const score = (
                ((c.qed || 0) * (weights.qed || 1.0) * 10) +
                ((c.pIC50 || 0) * (weights.pIC50 || 1.0)) +
                ((c.logp || 0) * (weights.logp || 0)) +
                ((500 - (c.mw || 500)) * (weights.mw || 0.01)) +
                ((10 - (c.sascore || 5)) * (weights.sascore || 0.5))
            );
            return { ...c, weighted_score: score };
        }).filter(Boolean).sort((a, b) => b.weighted_score - a.weighted_score);
    }, [missionData?.chembl_summary?.top_candidates, missionData?.active_weights]);

    const renderBellCurve = (data, color, label) => {
        if (!data || data.length === 0) return (
            <div style={{ height: '180px', display: 'flex', alignItems: 'center', justifyContent: 'center', color: '#444', fontSize: '11px' }}>
                Waiting for data grounding...
            </div>
        );
        return (
            <ResponsiveContainer width="100%" height={180}>
                <AreaChart data={data}>
                    <defs>
                        <linearGradient id={`grad-${label}`} x1="0" y1="0" x2="0" y2="1">
                            <stop offset="5%" stopColor={color} stopOpacity={0.8} />
                            <stop offset="95%" stopColor={color} stopOpacity={0} />
                        </linearGradient>
                    </defs>
                    <CartesianGrid strokeDasharray="3 3" stroke="#222" />
                    <XAxis dataKey="bin" stroke="#666" fontSize={10} />
                    <YAxis hide />
                    <Tooltip
                        contentStyle={{ backgroundColor: '#111', border: '1px solid #333', fontSize: '10px' }}
                        itemStyle={{ color: color }}
                    />
                    <Area type="monotone" dataKey="count" stroke={color} fillOpacity={1} fill={`url(#grad-${label})`} />
                </AreaChart>
            </ResponsiveContainer>
        );
    };

    const PropertyBar = ({ label, value, max, color = "#4285F4" }) => (
        <div className="prop-bar-container">
            <div className="prop-label"><span>{label}</span> <span>{(value || 0).toFixed(2)}</span></div>
            <div className="prop-track">
                <div className="prop-fill" style={{ width: `${Math.min(100, ((value || 0) / max) * 100)}%`, background: color }} />
            </div>
        </div>
    );

    return (
        <div className="live-workspace">
            <div className="workspace-header">
                <div className="tabs">
                    <button className={activeView === 'summary' ? 'active' : ''} onClick={() => setActiveView('summary')}>SUMMARY</button>
                    <button className={activeView === 'research' ? 'active' : ''} onClick={() => setActiveView('research')}>LITERATURE</button>
                    <button className={activeView === '3d' ? 'active' : ''} onClick={() => setActiveView('3d')}>STRUCTURAL</button>
                    <button className={activeView === 'stats' ? 'active' : ''} onClick={() => setActiveView('stats')}>GROUNDING</button>
                    <button className={activeView === 'retro' ? 'active' : ''} onClick={() => setActiveView('retro')}>SYNTHESIS</button>
                    <button className={activeView === 'logs' ? 'active' : ''} onClick={() => setActiveView('logs')}>LOGS</button>
                </div>
            </div>

            <div className="workspace-content">
                {activeView === 'summary' && (
                    <div className="view-panel">
                        <div className="summary-dashboard">
                            <div className="summary-header">
                                <h3>Mission Report: {missionData?.target_name || missionId}</h3>
                                <div className="status-badge">Phase: Step {missionData?.current_step || 1}</div>
                            </div>

                            <GeminiAdvantageShowcase missionData={missionData} />

                            <FileUploadPanel missionId={missionId} />

                            <div className="summary-grid">
                                <div className="summary-card">
                                    <h4>Primary Target</h4>
                                    <div className="card-val">{missionData?.chembl_summary?.target_name || "Scanning..."}</div>
                                    <p style={{ color: '#888', fontSize: '12px' }}>ChEMBL ID: {missionData?.chembl_summary?.target_id || "N/A"}</p>
                                </div>
                                <div className="summary-card best-candidate-profile">
                                    <h4>Best Candidate Profile</h4>
                                    {missionData?.best_candidate_profile ? (
                                        <div className="profile-layout">
                                            <div className="profile-img">
                                                <img src={`http://localhost:8000/api/render?smiles=${encodeURIComponent(missionData.best_candidate_profile.smiles)}`} alt="best cand" />
                                            </div>
                                            <div className="profile-details">
                                                <div className="profile-score">Score: {missionData.best_candidate_profile.score?.toFixed(2)}</div>
                                                <div className="profile-status">{missionData.best_candidate_profile.status}</div>
                                                <div className="profile-props">
                                                    <PropertyBar label="QED" value={missionData.best_candidate_profile.properties?.qed || 0} max={1} color="#4CAF50" />
                                                    <PropertyBar label="SA" value={10 - (missionData.best_candidate_profile.properties?.sa || 5)} max={10} color="#FFC107" />
                                                    <PropertyBar label="pIC50" value={missionData.best_candidate_profile.properties?.pIC50 || 0} max={12} color="#FF8042" />
                                                </div>
                                                <p className="profile-rationale">{missionData.best_candidate_profile.rationale}</p>
                                            </div>
                                        </div>
                                    ) : (
                                        <div className="profile-placeholder">Running Pareto-front analysis on grounding set...</div>
                                    )}
                                </div>

                                <div className="summary-card wide final-synthesis-card">
                                    <h4>Final Mission Synthesis</h4>
                                    {missionData?.final_summary ? (
                                        <div className="final-summary-box">
                                            <ReactMarkdown
                                                remarkPlugins={[remarkGfm]}
                                                disallowedElements={['a']}
                                                unwrapDisallowed={true}
                                            >
                                                {missionData?.final_report || ""}
                                            </ReactMarkdown>
                                        </div>
                                    ) : (
                                        <div className="synthesis-pending">
                                            <p>The mission is currently in Phase {missionData?.current_step || 1}. The final comprehensive report will be generated and published here once the synthesis and verification phase is completed (Step 7).</p>
                                            <div className="dynamic-progress-bar">
                                                <div className="progress-fill" style={{ width: `${((missionData?.current_step || 1) / 7) * 100}%` }} />
                                            </div>
                                        </div>
                                    )}
                                </div>
                            </div>
                        </div>
                    </div>
                )}

                {activeView === 'research' && (
                    <div className="view-panel research-split">
                        <div className="research-sources">
                            <h3>Grounding Sources (ArXiv)</h3>
                            <div className="source-list">
                                {(missionData?.research_papers || []).map((p, i) => (
                                    <div key={i} className="source-card">
                                        <h5>{p.title}</h5>
                                        <p className="abstract-peek">{p.summary?.substring(0, 150)}...</p>
                                        <div className="source-footer">
                                            <span className="author-tag">{p.authors?.[0] || "Unknown"} et al.</span>
                                            <a href={p.url} target="_blank" rel="noreferrer" className="read-btn">PDF ↗</a>
                                        </div>
                                    </div>
                                ))}
                                {(!missionData?.research_papers || missionData.research_papers.length === 0) && (
                                    <div className="internal-research-status">
                                        {missionData?.internal_research_papers?.length > 0 ? (
                                            <p className="internal-msg">
                                                <span className="dot pulse yellow" /> Grounding logic cached {missionData.internal_research_papers.length} sources.
                                                Summary will be published once Step 1 research satisfies the mission strategy.
                                            </p>
                                        ) : (
                                            <p className="empty-msg">No papers found yet.</p>
                                        )}
                                    </div>
                                )}
                            </div>
                        </div>
                        <div className="research-synthesis">
                            <h3>Intelligence Summary</h3>
                            <div className="summary-box">
                                {missionData?.research_summary ? (
                                    <ReactMarkdown
                                        remarkPlugins={[remarkGfm]}
                                        disallowedElements={['a']}
                                        unwrapDisallowed={true}
                                    >
                                        {missionData?.research_summary || "Synthesis will appear after literature review..."}
                                    </ReactMarkdown>
                                ) : (
                                    <div className="synthesis-placeholder">
                                        {missionData?.internal_research_summary ? (
                                            <div className="staging-indicator">
                                                <div className="spinner-small" />
                                                <p>Analyzing and refining research intelligence internally...</p>
                                            </div>
                                        ) : (
                                            <p>Initiating literature grounding search...</p>
                                        )}
                                    </div>
                                )}
                            </div>
                        </div>
                    </div>
                )}

                {activeView === '3d' && (
                    <div className="view-panel" style={{ height: '100%', display: 'flex', flexDirection: 'column' }}>
                        <div className="viewer-header">
                            <h3>Structural Analysis</h3>
                            <div className="viewer-controls">
                                {missionData?.active_pdb && <span className="pdb-badge">Protein: {missionData.active_pdb}</span>}
                                {missionData?.last_docking_pose && <span className="complex-badge">Ligand Bound</span>}
                            </div>
                        </div>
                        <div className="viewer-main">
                            {missionData?.active_pdb ? (
                                <MolstarViewer
                                    pdbId={missionData.active_pdb}
                                    ligandUrl={missionData?.last_docking_pose?.ligand_url}
                                    receptorUrl={missionData?.last_docking_pose?.receptor_url}
                                />
                            ) : (
                                <div className="empty-viewer">Structure not selected.</div>
                            )}
                        </div>
                    </div>
                )}

                {activeView === 'stats' && (
                    <div className="view-panel stats-layout">
                        <div className="stats-main">
                            <div className="grounding-header">
                                <h3>Property Distributions</h3>
                                <div className="chart-controls">
                                    <select value={chartMetric} onChange={(e) => setChartMetric(e.target.value)}>
                                        <option value="pic50">pIC50 (Potency)</option>
                                        <option value="mw">Molecular Weight</option>
                                        <option value="logp">LogP (Lipophilicity)</option>
                                        <option value="qed">QED (Drug-likeness)</option>
                                        <option value="sascore">SA Score (Synthetic Accessibility)</option>
                                    </select>
                                </div>
                            </div>
                            <div className="main-chart-box">
                                {chartMetric === 'pic50' && renderBellCurve(missionData?.chembl_summary?.pic50_bins, "#FF8042", "pic50")}
                                {chartMetric === 'mw' && renderBellCurve(missionData?.chembl_summary?.mw_bins, "#4285F4", "mw")}
                                {chartMetric === 'logp' && renderBellCurve(missionData?.chembl_summary?.logp_bins, "#82ca9d", "logp")}
                                {chartMetric === 'qed' && renderBellCurve(missionData?.chembl_summary?.qed_bins, "#8884d8", "qed")}
                                {chartMetric === 'sascore' && renderBellCurve(missionData?.chembl_summary?.sascore_bins, "#ffc658", "sa")}
                                <div className="chart-meta-footer">
                                    Dataset: {missionData?.chembl_summary?.count || 0} molecules from ChEMBL {missionData?.chembl_summary?.target_id}
                                </div>
                            </div>
                        </div>

                        <div className="potency-leaderboard wider">
                            <div className="leaderboard-header">
                                <h3>Weighted Leaderboard</h3>
                                <div className="view-switcher">
                                    <button className={leaderboardView === 'table' ? 'active' : ''} onClick={() => setLeaderboardView('table')}>Table</button>
                                    <button className={leaderboardView === 'cards' ? 'active' : ''} onClick={() => setLeaderboardView('cards')}>Cards</button>
                                </div>
                            </div>

                            <div className="leaderboard-scroll">
                                {leaderboardView === 'table' ? (
                                    <table>
                                        <thead>
                                            <tr>
                                                <th>Molecule</th>
                                                <th>Properties</th>
                                                <th>Weighted Score</th>
                                            </tr>
                                        </thead>
                                        <tbody>
                                            {(sortedCandidates || []).map((cand, i) => (
                                                <tr key={i} className="potency-row">
                                                    <td className="mol-cell">
                                                        <div className="mol-preview-mini">
                                                            {cand.smiles ? (
                                                                <img
                                                                    src={`http://localhost:8000/api/render?smiles=${encodeURIComponent(cand.smiles)}`}
                                                                    alt="mol"
                                                                    onError={(e) => { e.target.src = "https://via.placeholder.com/44?text=Mol"; }}
                                                                />
                                                            ) : (
                                                                <div className="mol-placeholder">NO SMILES</div>
                                                            )}
                                                        </div>
                                                        <div>
                                                            <div style={{ fontWeight: 800 }}>{cand.molecule_id}</div>
                                                            <div style={{ fontSize: '9px', color: '#666' }}>{cand.smiles?.substring(0, 20)}...</div>
                                                        </div>
                                                    </td>
                                                    <td>
                                                        <div className="prop-dots">
                                                            <span title={`pIC50: ${cand.pIC50?.toFixed(2)}`} style={{ color: '#FF8042' }}>P: {cand.pIC50?.toFixed(1)}</span>
                                                            <span title={`QED: ${cand.qed?.toFixed(2)}`} style={{ color: '#4CAF50' }}>Q: {cand.qed?.toFixed(2)}</span>
                                                            <span title={`MW: ${cand.mw}`} style={{ color: '#4285F4' }}>MW: {cand.mw}</span>
                                                        </div>
                                                    </td>
                                                    <td className="score-cell">
                                                        <div className="score-val">{cand.weighted_score.toFixed(1)}</div>
                                                    </td>
                                                </tr>
                                            ))}
                                            {(sortedCandidates || []).length === 0 && (
                                                <tr><td colSpan="3" className="empty-msg">No candidates found in grounding set.</td></tr>
                                            )}
                                        </tbody>
                                    </table>
                                ) : (
                                    <div className="card-grid">
                                        {(sortedCandidates || []).map((cand, i) => (
                                            <div key={i} className="cand-card">
                                                <div className="cand-img">
                                                    {cand.smiles ? (
                                                        <img
                                                            src={`http://localhost:8000/api/render?smiles=${encodeURIComponent(cand.smiles)}`}
                                                            alt="mol"
                                                            onError={(e) => { e.target.src = "https://via.placeholder.com/140?text=Structure+Unavailable"; }}
                                                        />
                                                    ) : (
                                                        <div className="mol-placeholder">NO SMILES</div>
                                                    )}
                                                </div>
                                                <div className="cand-info">
                                                    <div className="cand-title">{cand.molecule_id}</div>
                                                    <div className="cand-score">Rank Score: {cand.weighted_score.toFixed(1)}</div>
                                                    <div className="cand-props">
                                                        <PropertyBar label="pIC50" value={cand.pIC50} max={12} color="#FF8042" />
                                                        <PropertyBar label="QED" value={cand.qed} max={1} color="#4CAF50" />
                                                    </div>
                                                </div>
                                            </div>
                                        ))}
                                        {(sortedCandidates || []).length === 0 && <p className="empty-msg">No candidates found.</p>}
                                    </div>
                                )}
                            </div>
                        </div>
                    </div>
                )}

                {activeView === 'retro' && (
                    <div className="view-panel">
                        <h3>Retrosynthesis Path (AiZynthFinder)</h3>
                        <div className="retro-flow">
                            {(missionData?.retro_data?.routes || []).map((route, i) => (
                                <div key={i} className="route-block">
                                    <div className="route-header">
                                        <h4>Route {i + 1}</h4>
                                        <div className="route-score-badge">State Score: {route.score}</div>
                                    </div>
                                    <div className="route-viz">
                                        {route.image_url ? (
                                            <img src={route.image_url} alt="pathway" className="pathway-img" />
                                        ) : (
                                            <div className="steps-container">
                                                {(route.steps || []).map((s, j) => (
                                                    <div key={j} className="retro-step">
                                                        <div className="step-img-box">
                                                            {s.reactant_smiles?.[0] ? (
                                                                <img src={`http://localhost:8000/api/render?smiles=${encodeURIComponent(s.reactant_smiles[0])}`} alt="molecule" />
                                                            ) : (
                                                                <div className="mol-placeholder" />
                                                            )}
                                                        </div>
                                                        <div className="step-arrow">→</div>
                                                    </div>
                                                ))}
                                            </div>
                                        )}
                                    </div>
                                </div>
                            ))}
                            {(!missionData?.retro_data?.routes || missionData.retro_data.routes.length === 0) && (
                                <p className="empty-msg">Synthesizability analysis will run in Step 7.</p>
                            )}
                        </div>
                    </div>
                )}

                {activeView === 'logs' && (
                    <div className="view-panel">
                        <h3>Detailed Tool Operation History</h3>
                        <div className="logs-table-wrapper">
                            <table>
                                <thead>
                                    <tr>
                                        <th>Time</th>
                                        <th>Tool Name</th>
                                        <th>Arguments</th>
                                        <th>Status</th>
                                    </tr>
                                </thead>
                                <tbody>
                                    {(missionData?.tool_history || []).map((log, i) => (
                                        <tr key={i}>
                                            <td style={{ color: '#666', fontSize: '11px' }}>{new Date(log.timestamp).toLocaleTimeString()}</td>
                                            <td style={{ color: '#4285F4', fontWeight: 'bold' }}>{log.tool}</td>
                                            <td><code style={{ fontSize: '10px', color: '#888' }}>{JSON.stringify(log.args)}</code></td>
                                            <td><code className="status-pill">COMPLETED</code></td>
                                        </tr>
                                    ))}
                                </tbody>
                            </table>
                        </div>
                    </div>
                )}
            </div>

            <style>{`
                .live-workspace { display: flex; flex-direction: column; height: 100%; background: #000; color: #fff; font-family: 'Inter', sans-serif; }
                .workspace-header { padding: 5px 20px; background: #0a0a0a; border-bottom: 2px solid #1a1a1a; }
                .tabs { display: flex; gap: 20px; }
                .tabs button { background: none; border: none; color: #444; font-size: 10px; font-weight: 800; letter-spacing: 1px; padding: 12px 0; cursor: pointer; text-transform: uppercase; border-bottom: 2px solid transparent; transition: all 0.2s; }
                .tabs button.active { color: #4285F4; border-color: #4285F4; }
                .workspace-content { flex-grow: 1; padding: 25px; overflow-y: auto; }
                
                .summary-dashboard { display: flex; flex-direction: column; gap: 20px; }
                .summary-header { display: flex; justify-content: space-between; align-items: center; }
                .summary-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }
                .summary-card { background: #0d0d0d; border: 1px solid #1a1a1a; padding: 25px; border-radius: 12px; }
                .summary-card.wide { grid-column: span 2; }
                .card-val { font-size: 28px; font-weight: 800; color: #4285F4; margin: 10px 0; line-height: 1.1; }
                
                .best-candidate-profile { border-left: 4px solid #4285F4; }
                .best-candidate-profile .profile-layout { display: flex; gap: 20px; align-items: center; }
                .profile-img { background: #fff; padding: 10px; border-radius: 8px; width: 140px; height: 140px; border: 1px solid #333; }
                .profile-img img { width: 100%; height: 100%; object-fit: contain; }
                .profile-details { flex-grow: 1; }
                .profile-score { font-size: 20px; font-weight: 800; color: #4CAF50; }
                .profile-status { font-size: 10px; background: #333; display: inline-block; padding: 2px 8px; border-radius: 4px; margin: 5px 0; text-transform: uppercase; letter-spacing: 1px; color: #bbb; }
                .profile-props { margin-top: 10px; }
                .profile-rationale { font-size: 11px; color: #888; margin-top: 10px; line-height: 1.4; border-top: 1px solid #222; padding-top: 8px; }

                .prop-bar-container { margin-bottom: 8px; }
                .prop-label { display: flex; justify-content: space-between; font-size: 9px; color: #666; margin-bottom: 2px; font-weight: bold; }
                .prop-track { height: 4px; background: #222; border-radius: 2px; }
                .prop-fill { height: 100%; border-radius: 2px; transition: width 0.5s ease-out; }

                .final-synthesis-card { background: linear-gradient(145deg, #0d0d0d 0%, #151515 100%); position: relative; overflow: hidden; }
                .final-synthesis-card::after { content: "SYNTHESIS"; position: absolute; top: -10px; right: -10px; font-size: 60px; font-weight: 900; color: rgba(255,255,255,0.02); pointer-events: none; }
                .final-summary-box { color: #ddd; font-size: 14px; line-height: 1.6; padding: 10px; }
                .synthesis-pending { padding: 40px; text-align: center; color: #666; }
                .dynamic-progress-bar { height: 6px; background: #1a1a1a; border-radius: 3px; margin-top: 20px; overflow: hidden; }
                .dynamic-progress-bar .progress-fill { height: 100%; background: #4285F4; transition: width 1s ease-in-out; box-shadow: 0 0 10px #4285F4; }

                .stats-layout { display: grid; grid-template-columns: 1fr 450px; gap: 25px; align-items: start; }
                .main-chart-box { background: #0d0d0d; border: 1px solid #1a1a1a; padding: 25px; border-radius: 12px; }
                .chart-controls select { background: #111; border: 1px solid #333; color: #fff; padding: 6px 12px; border-radius: 6px; font-size: 12px; cursor: pointer; }
                .chart-meta-footer { font-size: 10px; color: #444; margin-top: 15px; text-align: right; }

                .potency-leaderboard { background: #0d0d0d; border: 1px solid #1a1a1a; border-radius: 12px; padding: 20px; }
                .leaderboard-header { display: flex; justify-content: space-between; align-items: center; margin-bottom: 15px; }
                .view-switcher { display: flex; background: #000; padding: 3px; border-radius: 8px; border: 1px solid #222; }
                .view-switcher button { background: none; border: none; color: #444; padding: 6px 16px; font-size: 10px; cursor: pointer; border-radius: 6px; font-weight: bold; transition: all 0.2s; }
                .view-switcher button.active { background: #4285F4; color: #fff; }

                .leaderboard-scroll { max-height: 600px; overflow-y: auto; }
                .leaderboard-scroll table { width: 100%; font-size: 11px; border-collapse: collapse; }
                .leaderboard-scroll th { border-bottom: 2px solid #1a1a1a; color: #666; padding: 10px; text-align: left; }
                .potency-row td { padding: 12px 10px; border-bottom: 1px solid #111; }
                .research-split { display: grid; grid-template-columns: 350px 1fr; gap: 25px; height: 100%; overflow: hidden; }
                .research-sources { display: flex; flex-direction: column; background: #0a0a0a; border: 1px solid #1a1a1a; border-radius: 12px; padding: 20px; overflow: hidden; }
                .source-list { flex-grow: 1; overflow-y: auto; display: flex; flex-direction: column; gap: 12px; margin-top: 15px; }
                .source-card { background: #111; border: 1px solid #222; padding: 15px; border-radius: 8px; border-left: 3px solid #4285F4; }
                .source-card h5 { font-size: 11px; margin-bottom: 8px; color: #fff; line-height: 1.3; }
                .abstract-peek { font-size: 10px; color: #666; line-height: 1.4; margin-bottom: 10px; }
                .source-footer { display: flex; justify-content: space-between; align-items: center; }
                .author-tag { font-size: 9px; color: #4285F4; font-weight: bold; background: rgba(66, 133, 244, 0.1); padding: 2px 6px; border-radius: 4px; }
                .read-btn { font-size: 9px; color: #fff; text-decoration: none; border: 1px solid #333; padding: 2px 8px; border-radius: 4px; transition: all 0.2s; }
                .read-btn:hover { background: #222; border-color: #4285F4; }

                .research-synthesis { display: flex; flex-direction: column; background: #0d0d0d; border: 1px solid #1a1a1a; border-radius: 12px; padding: 25px; overflow: hidden; }
                .summary-box { flex-grow: 1; overflow-y: auto; color: #ddd; font-size: 14px; line-height: 1.8; }
                .synthesis-placeholder { height: 100%; display: flex; align-items: center; justify-content: center; color: #444; text-align: center; }

                .mol-cell { display: flex; align-items: center; gap: 10px; }
                .mol-preview-mini { background: #fff; width: 44px; height: 44px; border-radius: 6px; padding: 3px; border: 1px solid #333; display: flex; align-items: center; justify-content: center; overflow: hidden; }
                .mol-preview-mini img { max-width: 100%; max-height: 100%; object-fit: contain; }
                .prop-dots { display: flex; gap: 10px; font-size: 10px; font-weight: bold; }
                .score-cell .score-val { color: #4285F4; font-weight: 800; font-size: 18px; }

                .card-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 15px; }
                .cand-card { background: #111; border: 1px solid #222; border-radius: 10px; overflow: hidden; transition: transform 0.2s; }
                .cand-card:hover { transform: translateY(-5px); border-color: #4285F4; }
                .cand-img { background: #fff; height: 160px; display: flex; align-items: center; justify-content: center; padding: 15px; border-bottom: 1px solid #222; }
                .cand-img img { width: 100%; height: 100%; object-fit: contain; }
                .cand-info { padding: 15px; }
                .cand-title { font-size: 11px; font-weight: bold; color: #4285F4; margin-bottom: 6px; font-family: monospace; }
                .cand-score { font-size: 18px; font-weight: 900; color: #fff; margin-bottom: 12px; }

                .route-block { background: #0d0d0d; border: 1px solid #1a1a1a; padding: 25px; border-radius: 12px; margin-bottom: 25px; }
                .route-header { display: flex; justify-content: space-between; align-items: center; margin-bottom: 20px; }
                .route-score-badge { background: #222; padding: 5px 12px; border-radius: 6px; font-size: 11px; color: #4CAF50; font-weight: bold; }
                .pathway-img { width: 100%; border-radius: 10px; border: 1px solid #333; background: #fff; }

                @keyframes spin { 0% { transform: rotate(0deg); } 100% { transform: rotate(360deg); } }
                .spinner-small { width: 24px; height: 24px; border: 3px solid #1a1a1a; border-top: 3px solid #4285F4; border-radius: 50%; animation: spin 0.8s linear infinite; margin: 0 auto; }
                .empty-msg { text-align: center; color: #444; padding: 40px; font-size: 12px; }
            `}</style>
        </div>
    );
};

export default LiveWorkspace;
