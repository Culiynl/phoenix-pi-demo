import React, { useState, useMemo } from 'react';
import MolstarViewer from './MolstarViewer';
import {
    ResponsiveContainer, AreaChart, Area, XAxis, YAxis, Tooltip, CartesianGrid,
    BarChart, Bar, Cell
} from 'recharts';
import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';

const LiveWorkspace = ({ missionId, missionData, activeView = 'summary', setActiveView }) => {
    const [chartMetric, setChartMetric] = useState('pic50'); // pic50, mw, logp, qed, sa
    const [leaderboardView, setLeaderboardView] = useState('table'); // table, cards
    const [localPdbData, setLocalPdbData] = useState(null);
    const [activeRouteIndex, setActiveRouteIndex] = useState(0);

    React.useEffect(() => {
        if (!missionId) return;
        setActiveRouteIndex(0); // Reset route when switching missions
    }, [missionId]);

    React.useEffect(() => {
        if (missionData?.active_pdb) setLocalPdbData(null);
    }, [missionData?.active_pdb]);

    const handleLocalPdbUpload = (e) => {
        const file = e.target.files[0];
        if (!file) return;

        const reader = new FileReader();
        reader.onload = (event) => {
            setLocalPdbData(event.target.result); // Pass string directly to viewer
        };
        reader.readAsText(file);
    };

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

    // Case-insensitive property lookup — backend mixes casing (QED vs qed, MW vs mw, etc.)
    const getProp = (props, ...keys) => {
        if (!props) return undefined;
        // Try exact keys first, then lowercase comparison
        for (const key of keys) {
            if (props[key] !== undefined && props[key] !== null) return props[key];
        }
        const lower = keys.map(k => k.toLowerCase());
        for (const [k, v] of Object.entries(props)) {
            if (lower.includes(k.toLowerCase()) && v !== undefined && v !== null) return v;
        }
        return undefined;
    };

    const PropertyBar = ({ label, value, max, color = "#EA4335" }) => (
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
                    <button className={activeView === 'analysis' ? 'active' : ''} onClick={() => setActiveView('analysis')}>MOLECULE ANALYSIS</button>
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

                            {/* Integrated Chat / Message Phoenix removed for clarity */}

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
                                                    <PropertyBar label="QED" value={getProp(missionData.best_candidate_profile.properties, 'qed', 'QED') || 0} max={1} color="#4CAF50" />
                                                    <PropertyBar label="SA" value={10 - (getProp(missionData.best_candidate_profile.properties, 'sa', 'SA', 'sascore', 'SAScore') || 5)} max={10} color="#FFC107" />
                                                    <PropertyBar label="pIC50" value={getProp(missionData.best_candidate_profile.properties, 'pIC50', 'pic50', 'pIC50_mean') || 0} max={12} color="#FF8042" />
                                                </div>
                                                <p className="profile-rationale">{missionData.best_candidate_profile.rationale}</p>
                                            </div>
                                        </div>
                                    ) : (
                                        <div className="profile-placeholder">Running Pareto-front analysis on grounding set...</div>
                                    )}
                                </div>

                                <div className="summary-card wide plot-card">
                                    <h4>Potency vs Drug-likeness (Champion Search)</h4>
                                    <div className="discovery-plot-container" style={{ height: '300px' }}>
                                        <ResponsiveContainer width="100%" height="100%">
                                            <AreaChart data={missionData?.chembl_summary?.pic50_bins || []}>
                                                <CartesianGrid strokeDasharray="3 3" stroke="#222" />
                                                <XAxis dataKey="bin" stroke="#71717a" fontSize={10} label={{ value: 'pIC50', position: 'insideBottom', offset: -5, fill: '#71717a' }} />
                                                <YAxis stroke="#71717a" fontSize={10} />
                                                <Tooltip contentStyle={{ background: '#000', border: '1px solid #333' }} />
                                                <Area type="monotone" dataKey="count" stroke="#ef4444" fill="#ef4444" fillOpacity={0.2} />
                                            </AreaChart>
                                        </ResponsiveContainer>
                                    </div>
                                </div>

                                <div className="summary-card wide final-synthesis-card">
                                    <h4>Final Mission Synthesis</h4>
                                    {(missionData?.final_summary || missionData?.final_report) ? (
                                        <div className="final-summary-box">
                                            <ReactMarkdown
                                                remarkPlugins={[remarkGfm]}
                                                disallowedElements={['a']}
                                                unwrapDisallowed={true}
                                            >
                                            </ReactMarkdown>
                                            <div style={{ marginTop: '24px' }}>
                                                <button
                                                    className="pdf-export-btn"
                                                    onClick={() => window.open(`http://localhost:8000/api/mission/${missionId}/export/pdf`, '_blank')}
                                                >
                                                    <span className="pdf-icon">⬇</span> Export Final Report (PDF)
                                                </button>
                                            </div>
                                        </div>
                                    ) : (
                                        <div className="synthesis-pending">
                                            <p>The mission is currently in Phase {missionData?.current_step || 1}. The final comprehensive report will be generated and published here once the synthesis and verification phase is completed (Step 7).</p>
                                            <div className="dynamic-progress-bar">
                                                <div className="progress-fill" style={{ width: `${((missionData?.current_step || 1) / 7) * 100}%` }} />
                                            </div>
                                            <div style={{ marginTop: '24px' }}>
                                                <button
                                                    className="pdf-export-btn"
                                                    onClick={() => window.open(`http://localhost:8000/api/mission/${missionId}/export/pdf`, '_blank')}
                                                >
                                                    <span className="pdf-icon">⬇</span> Export Final Report (PDF)
                                                </button>
                                            </div>
                                        </div>
                                    )}
                                </div>
                            </div>
                        </div>
                    </div>
                )}


                {activeView === 'analysis' && (
                    <div className="view-panel">
                        <div className="analysis-dashboard">
                            <div className="analysis-header">
                                <h3>Best Candidate Analysis</h3>
                                {missionData?.best_candidate_profile?.smiles && (
                                    <div className="header-actions">
                                        <button
                                            className="primary-btn download-btn"
                                            onClick={() => window.open(`http://localhost:8000/api/mission/${missionId}/download/pdb`, '_blank')}
                                        >
                                            Download PDB Complex
                                        </button>
                                        <button
                                            className="secondary-btn download-btn"
                                            onClick={() => {
                                                const element = document.createElement("a");
                                                const file = new Blob([missionData.best_candidate_profile.smiles + " " + missionData.best_candidate_profile.molecule_id], { type: 'text/plain' });
                                                element.href = URL.createObjectURL(file);
                                                element.download = `${missionData.best_candidate_profile.molecule_id || "molecule"}.sdf`;
                                                document.body.appendChild(element);
                                                element.click();
                                            }}
                                        >
                                            Download SDF (Simulated)
                                        </button>
                                    </div>
                                )}
                            </div>

                            {/* Interaction Distribution Chart */}
                            {missionData?.best_candidate_profile?.interactions && (
                                <div className="interaction-chart-card card wide" style={{ marginTop: '20px' }}>
                                    <h4>Structural Interaction Distribution</h4>
                                    <div style={{ height: '240px' }}>
                                        <ResponsiveContainer width="100%" height="100%">
                                            <BarChart data={(() => {
                                                const types = {};
                                                (missionData.best_candidate_profile.interactions.plip || []).forEach(it => {
                                                    types[it.type] = (types[it.type] || 0) + 1;
                                                });
                                                (missionData.best_candidate_profile.interactions.key_residues || []).forEach(it => {
                                                    types[it.interaction || 'Contact'] = (types[it.interaction || 'Contact'] || 0) + 1;
                                                });
                                                return Object.entries(types).map(([name, count]) => ({ name, count }));
                                            })()}>
                                                <CartesianGrid strokeDasharray="3 3" stroke="#1a1a1a" vertical={false} />
                                                <XAxis dataKey="name" stroke="#71717a" fontSize={10} />
                                                <YAxis stroke="#71717a" fontSize={10} allowDecimals={false} />
                                                <Tooltip contentStyle={{ background: '#000', border: '1px solid #333', fontSize: '11px' }} />
                                                <Bar dataKey="count" radius={[4, 4, 0, 0]}>
                                                    {(() => {
                                                        const colors = {
                                                            'Hydrophobic': '#aaa',
                                                            'H-Bond': '#4CAF50',
                                                            'Salt Bridge': '#FF8042',
                                                            'Pi-Stacking': '#9C27B0',
                                                            'Contact': '#EA4335',
                                                            'van der Waals': '#71717a'
                                                        };
                                                        return Object.keys(colors).map((key, index) => (
                                                            <Cell key={`cell-${index}`} fill={colors[key] || '#555'} />
                                                        ));
                                                    })()}
                                                </Bar>
                                            </BarChart>
                                        </ResponsiveContainer>
                                    </div>
                                </div>
                            )}

                            <div className="analysis-grid">
                                <div className="mol-structure-card">
                                    <h4>Structure 2D</h4>
                                    <div className="structure-img-large">
                                        {missionData?.best_candidate_profile?.smiles ? (
                                            <img src={`http://localhost:8000/api/render?smiles=${encodeURIComponent(missionData.best_candidate_profile.smiles)}`} alt="Structure" />
                                        ) : <div className="placeholder">Analysis Pending...</div>}
                                    </div>
                                    <div className="smiles-box">
                                        <code>{missionData?.best_candidate_profile?.smiles || "No SMILES data"}</code>
                                    </div>
                                </div>

                                <div className="mol-props-card">
                                    <h4>ADMET & Properties</h4>
                                    <div className="props-list">
                                        {missionData?.best_candidate_profile?.properties ? (
                                            Object.entries(missionData.best_candidate_profile.properties).map(([k, v]) => (
                                                <div key={k} className="prop-row">
                                                    <span className="key">{k.toUpperCase()}</span>
                                                    <span className="val">{typeof v === 'number' ? v.toFixed(3) : String(v)}</span>
                                                </div>
                                            ))
                                        ) : (
                                            <p>No property data available.</p>
                                        )}

                                        {/* Dynamic ADMET Data */}
                                        {missionData?.best_candidate_profile?.admet && Object.entries(missionData.best_candidate_profile.admet).length > 0 ? (
                                            Object.entries(missionData.best_candidate_profile.admet).map(([k, v]) => (
                                                <div key={k} className="prop-row">
                                                    <span className="key">{k.replace(/_/g, " ")}</span>
                                                    <span className="val good">{typeof v === 'number' ? v.toFixed(2) : v.toString()}</span>
                                                </div>
                                            ))
                                        ) : (
                                            <div style={{ marginTop: '15px' }}>
                                                <small style={{ color: '#666' }}>ADMET Analysis pending...</small>
                                            </div>
                                        )}
                                    </div>
                                </div>

                                {/* 🔗 New Interaction Analysis Results */}
                                {missionData?.best_candidate_profile?.interactions && (
                                    <div className="interaction-card">
                                        <h4>Structural Interaction Mapping (PLIP)</h4>
                                        <div className="interaction-list">
                                            {missionData.best_candidate_profile.interactions.plip?.map((it, idx) => (
                                                <div key={idx} className="interaction-row">
                                                    <span className={`type-tag ${it.type.toLowerCase().replace(' ', '-')}`}>{it.type}</span>
                                                    <span className="residue">{it.residue}</span>
                                                    <span className="dist">{it.distance} Å</span>
                                                </div>
                                            ))}
                                            {missionData.best_candidate_profile.interactions.key_residues?.map((it, idx) => (
                                                <div key={`key-${idx}`} className="interaction-row key-res">
                                                    <span className="type-tag contact">Contact</span>
                                                    <span className="residue bold">{it.residue}</span>
                                                    <span className="atom">{it.prot_atom} - {it.lig_atom}</span>
                                                    <span className="dist">{it.distance} Å</span>
                                                </div>
                                            ))}
                                        </div>
                                        <div className="interaction-logic">
                                            <h5>Strategy Rationale</h5>
                                            <p>{missionData.best_candidate_profile.rationale}</p>
                                        </div>
                                    </div>
                                )}
                            </div>
                        </div>
                        <style>{`
                            .analysis-dashboard { display: flex; flex-direction: column; gap: 20px; }
                            .analysis-header { display: flex; justify-content: space-between; align-items: center; border-bottom: 1px solid #222; padding-bottom: 15px; }
                            .header-actions { display: flex; gap: 10px; }
                            .download-btn { font-size: 11px; padding: 6px 12px; }
                            .primary-btn { background: #4285F4; color: white; border: none; padding: 8px 16px; border-radius: 6px; cursor: pointer; font-weight: bold; }
                            .secondary-btn { background: #333; color: white; border: none; padding: 8px 16px; border-radius: 6px; cursor: pointer; font-weight: bold; }
                            .success-btn { background: #4CAF50; color: white; border: none; padding: 8px 16px; border-radius: 6px; cursor: pointer; font-weight: bold; }
                            .analysis-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 30px; }
                            .mol-structure-card { background: #111; padding: 20px; border-radius: 12px; border: 1px solid #222; }
                            .structure-img-large { background: #fff; border-radius: 8px; height: 300px; display: flex; align-items: center; justify-content: center; margin: 15px 0; }
                            .structure-img-large img { max-width: 100%; max-height: 100%; }
                            .smiles-box { background: #000; padding: 10px; border-radius: 6px; font-family: monospace; color: #EA4335; word-break: break-all; font-size: 11px; }
                            
                            .mol-props-card { background: #111; padding: 20px; border-radius: 12px; border: 1px solid #222; }
                            .props-list { display: flex; flex-direction: column; gap: 10px; margin-top: 15px; }
                            .prop-row { display: flex; justify-content: space-between; border-bottom: 1px solid #222; padding-bottom: 8px; }
                            .prop-row .key { color: #888; font-size: 12px; font-weight: bold; }
                            .prop-row .val { color: #ddd; font-weight: bold; }
                            .val.good { color: #4CAF50; }
                            .val.safe { color: #EA4335; }
                        `}</style>
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
                    <div className="view-panel-fixed">
                        <div className="viewer-header-clean">
                            <div className="header-left">
                                <h3>Structural Analysis</h3>
                                {!missionData?.active_pdb && (
                                    <div className="upload-control">
                                        <label htmlFor="pdb-upload" className="action-btn-blue">
                                            {localPdbData ? "Change Local PDB" : "Upload Local PDB"}
                                        </label>
                                        <input
                                            id="pdb-upload" type="file" accept=".pdb,.cif"
                                            onChange={handleLocalPdbUpload} style={{ display: 'none' }}
                                        />
                                    </div>
                                )}
                            </div>

                            <div className="header-right">
                                {missionData?.active_pdb ? (
                                    <span className="pdb-badge">Official: {missionData.active_pdb}</span>
                                ) : localPdbData ? (
                                    <span className="pdb-badge local-glow">Local Preview Active</span>
                                ) : (
                                    <span className="pdb-badge pending">No Structure Loaded</span>
                                )}
                            </div>
                        </div>

                        <div className="viewer-main-constrained">
                            <MolstarViewer
                                pdbId={missionData?.active_pdb}
                                customPdbData={missionData?.active_pdb ? null : localPdbData}
                                ligandUrl={missionData?.last_docking_pose?.ligand_url}
                                receptorUrl={missionData?.last_docking_pose?.receptor_url}
                                visualState={missionData?.visual_state}
                            />
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
                                {chartMetric === 'mw' && renderBellCurve(missionData?.chembl_summary?.mw_bins, "#EA4335", "mw")}
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
                                                            <span title={`MW: ${cand.mw}`} style={{ color: '#EA4335' }}>MW: {cand.mw}</span>
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
                        <div className="view-header-flex">
                            <h3>Retrosynthesis Path (AiZynthFinder)</h3>
                            {missionData?.retro_data?.routes?.length > 1 && (
                                <div className="route-switcher">
                                    <span className="switcher-label">Switch Route:</span>
                                    <div className="switcher-buttons">
                                        {missionData.retro_data.routes.map((r, idx) => (
                                            <button
                                                key={idx}
                                                className={`switcher-btn ${activeRouteIndex === idx ? 'active' : ''}`}
                                                onClick={() => setActiveRouteIndex(idx)}
                                            >
                                                Path {idx + 1} <span className="conf">({(r.score * 100).toFixed(0)}%)</span>
                                            </button>
                                        ))}
                                    </div>
                                </div>
                            )}
                        </div>

                        <div className="retro-flow-vertical">
                            {(() => {
                                const retroData = missionData?.retro_data;
                                const route = retroData?.routes?.[activeRouteIndex || 0];

                                // Case 1: AIZynth returned an error
                                if (retroData && !route && (retroData.error || retroData.routes?.length === 0)) {
                                    return (
                                        <div className="retro-empty-state">
                                            <div className="retro-placeholder-icon" style={{ color: '#f87171' }}>⚠</div>
                                            <p style={{ color: '#f87171', marginBottom: '10px' }}>
                                                AiZynthFinder could not generate graphical routes for this molecule.
                                            </p>
                                            {retroData.error && (
                                                <pre style={{ fontSize: '10px', color: '#666', background: '#0a0a0a', padding: '12px', borderRadius: '6px', maxHeight: '120px', overflowY: 'auto', textAlign: 'left' }}>
                                                    {retroData.error}
                                                </pre>
                                            )}
                                            {/* Show AI-reasoned synthesis from the final report if available */}
                                            {missionData?.final_summary && (
                                                <div style={{ marginTop: '20px', background: '#111', border: '1px solid #991b1b', borderRadius: '8px', padding: '20px', textAlign: 'left' }}>
                                                    <h4 style={{ color: '#f87171', fontSize: '12px', marginBottom: '10px', textTransform: 'uppercase', letterSpacing: '1px' }}>
                                                        ✦ AI Synthetic Feasibility Assessment
                                                    </h4>
                                                    <p style={{ color: '#aaa', fontSize: '13px', lineHeight: '1.7' }}>
                                                        {/* Extract synthesis section from final summary */}
                                                        {missionData.final_summary.split('\n').filter(l =>
                                                            l.toLowerCase().includes('synth') ||
                                                            l.toLowerCase().includes('coupling') ||
                                                            l.toLowerCase().includes('reaction') ||
                                                            l.toLowerCase().includes('precursor') ||
                                                            l.toLowerCase().includes('retro')
                                                        ).join('\n') || 'See Final Mission Synthesis tab for detailed synthetic feasibility assessment.'}
                                                    </p>
                                                </div>
                                            )}
                                        </div>
                                    );
                                }

                                // Case 2: No retro data yet
                                if (!route) {
                                    return (
                                        <div className="retro-empty-state">
                                            <div className="retro-placeholder-icon">⌬</div>
                                            <p>Retrosynthetic analysis will be generated automatically in Step 7.</p>
                                        </div>
                                    );
                                }


                                return (
                                    <div className="route-block-modern">
                                        <div className="route-header-modern">
                                            <div className="route-index">ROUTE {(activeRouteIndex || 0) + 1}</div>
                                            <div className="route-score-pill">Confidence: {(route.score * 100).toFixed(1)}%</div>
                                        </div>

                                        <div className="route-visual-container">
                                            {route.image_url ? (
                                                <div className="full-route-image-box">
                                                    <img src={route.image_url} alt={`Route ${(activeRouteIndex || 0) + 1}`} className="full-route-img" />
                                                    <div className="img-caption">Complete Synthetic Pathway (AiZynthFinder Native)</div>
                                                </div>
                                            ) : (
                                                <div className="placeholder-box">Visual Graph Unavailable</div>
                                            )}
                                        </div>

                                        <div className="route-details-modern">
                                            <h5>Synthetic Step Details (JSON Data)</h5>
                                            <table className="retro-details-table">
                                                <thead>
                                                    <tr>
                                                        <th>Step</th>
                                                        <th>Reaction Type</th>
                                                        <th>Starting Materials</th>
                                                        <th>Est. Yield</th>
                                                    </tr>
                                                </thead>
                                                <tbody>
                                                    {(route.steps || []).map((s, j) => (
                                                        <tr key={j}>
                                                            <td>{j + 1}</td>
                                                            <td className="rxn-column">{s.reaction_name || "Organic Transformation"}</td>
                                                            <td>
                                                                <div className="smmi-list">
                                                                    {(s.reactant_smiles || []).length > 0 ? (
                                                                        s.reactant_smiles.map((rsmi, k) => (
                                                                            <code key={k} title={rsmi}>{rsmi.substring(0, 15)}...</code>
                                                                        ))
                                                                    ) : (
                                                                        <span className="sm-empty">{s.reactants || "Commercial Precursor"}</span>
                                                                    )}
                                                                </div>
                                                            </td>
                                                            <td>{s.yield || "85%"}</td>
                                                        </tr>
                                                    ))}
                                                </tbody>
                                            </table>

                                            <div className="route-summary-stats">
                                                <div className="stat-pill"><b>Total Reactions:</b> {route.steps?.length || 0}</div>
                                                <div className="stat-pill"><b>Starting Materials:</b> {
                                                    // Count all unique reactants across all steps that aren't products of other steps
                                                    // For simplicity, we'll use the pre-calculated score or sum reactant_smiles
                                                    route.scores?.["number of pre-cursors"] ||
                                                    Array.from(new Set(route.steps?.flatMap(s => s.reactant_smiles || []))).length
                                                }</div>
                                                <div className="stat-pill"><b>Path Difficulty:</b> {route.score > 0.8 ? "Easy" : "Challenging"}</div>
                                            </div>
                                        </div>
                                    </div>
                                );
                            })()}
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
                                            <td style={{ color: '#EA4335', fontWeight: 'bold' }}>{log.tool}</td>
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
                .tabs button.active { color: #EA4335; border-color: #EA4335; }
                .workspace-content { flex-grow: 1; padding: 25px; overflow-y: auto; }
                
                .summary-dashboard { display: flex; flex-direction: column; gap: 20px; }
                .summary-header { display: flex; justify-content: space-between; align-items: center; }
                .summary-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }
                .summary-card { background: #0d0d0d; border: 1px solid #1a1a1a; padding: 25px; border-radius: 12px; }
                .summary-card.wide { grid-column: span 2; }
                .card-val { font-size: 28px; font-weight: 800; color: #EA4335; margin: 10px 0; line-height: 1.1; }
                
                .best-candidate-profile { border-left: 4px solid #EA4335; }
                .best-candidate-profile .profile-layout { display: flex; gap: 20px; align-items: center; }
                .profile-img { background: #fff; padding: 10px; border-radius: 8px; width: 640px; height: 140px; border: 1px solid #333; }
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
                .dynamic-progress-bar .progress-fill { height: 100%; background: #991b1b; transition: width 1s ease-in-out; box-shadow: 0 0 10px #991b1b; }

                /* PDF Export Button */
                .pdf-export-btn {
                    display: inline-flex;
                    align-items: center;
                    gap: 8px;
                    background: linear-gradient(135deg, #7f1d1d 0%, #991b1b 100%);
                    color: #fff;
                    border: 1px solid #7f1d1d;
                    padding: 11px 22px;
                    border-radius: 8px;
                    font-family: 'Inter', sans-serif;
                    font-size: 13px;
                    font-weight: 700;
                    letter-spacing: 0.5px;
                    cursor: pointer;
                    box-shadow: 0 0 16px rgba(153, 27, 27, 0.4), 0 2px 8px rgba(0,0,0,0.4);
                    transition: all 0.25s cubic-bezier(0.4, 0, 0.2, 1);
                    text-transform: uppercase;
                }
                .pdf-export-btn:hover {
                    background: linear-gradient(135deg, #991b1b 0%, #b91c1c 100%);
                    box-shadow: 0 0 28px rgba(153, 27, 27, 0.6), 0 4px 12px rgba(0,0,0,0.5);
                    transform: translateY(-2px);
                }
                .pdf-export-btn:active { transform: translateY(0px); box-shadow: 0 0 10px rgba(153,27,27,0.3); }
                .pdf-icon { font-size: 16px; }


                .stats-layout { display: grid; grid-template-columns: 1fr 450px; gap: 25px; align-items: start; }
                .main-chart-box { background: #0d0d0d; border: 1px solid #1a1a1a; padding: 25px; border-radius: 12px; }
                .chart-controls select { background: #111; border: 1px solid #333; color: #fff; padding: 6px 12px; border-radius: 6px; font-size: 12px; cursor: pointer; }
                .chart-meta-footer { font-size: 10px; color: #444; margin-top: 15px; text-align: right; }

                .potency-leaderboard { background: #0d0d0d; border: 1px solid #1a1a1a; border-radius: 12px; padding: 20px; }
                .leaderboard-header { display: flex; justify-content: space-between; align-items: center; margin-bottom: 15px; }
                .view-switcher { display: flex; background: #000; padding: 3px; border-radius: 8px; border: 1px solid #222; }
                .view-switcher button { background: none; border: none; color: #444; padding: 6px 16px; font-size: 10px; cursor: pointer; border-radius: 6px; font-weight: bold; transition: all 0.2s; }
                .view-switcher button.active { background: #EA4335; color: #fff; }

                .leaderboard-scroll { max-height: 600px; overflow-y: auto; }
                .leaderboard-scroll table { width: 100%; font-size: 11px; border-collapse: collapse; }
                .leaderboard-scroll th { border-bottom: 2px solid #1a1a1a; color: #666; padding: 10px; text-align: left; }
                .potency-row td { padding: 12px 10px; border-bottom: 1px solid #111; }
                .research-split { display: grid; grid-template-columns: 350px 1fr; gap: 25px; height: 100%; overflow: hidden; }
                .research-sources { display: flex; flex-direction: column; background: #0a0a0a; border: 1px solid #1a1a1a; border-radius: 12px; padding: 20px; overflow: hidden; }
                .source-list { flex-grow: 1; overflow-y: auto; display: flex; flex-direction: column; gap: 12px; margin-top: 15px; }
                .source-card { background: #111; border: 1px solid #222; padding: 15px; border-radius: 8px; border-left: 3px solid #EA4335; }
                .source-card h5 { font-size: 11px; margin-bottom: 8px; color: #fff; line-height: 1.3; }
                .abstract-peek { font-size: 10px; color: #666; line-height: 1.4; margin-bottom: 10px; }
                .source-footer { display: flex; justify-content: space-between; align-items: center; }
                .author-tag { font-size: 9px; color: #EA4335; font-weight: bold; background: rgba(234, 67, 53, 0.1); padding: 2px 6px; border-radius: 4px; }
                .read-btn { font-size: 9px; color: #fff; text-decoration: none; border: 1px solid #333; padding: 2px 8px; border-radius: 4px; transition: all 0.2s; }
                .read-btn:hover { background: #222; border-color: #EA4335; }

                .research-synthesis { display: flex; flex-direction: column; background: #0d0d0d; border: 1px solid #1a1a1a; border-radius: 12px; padding: 25px; overflow: hidden; }
                .summary-box { flex-grow: 1; overflow-y: auto; color: #ddd; font-size: 14px; line-height: 1.8; }
                .synthesis-placeholder { height: 100%; display: flex; align-items: center; justify-content: center; color: #444; text-align: center; }

                .mol-cell { display: flex; align-items: center; gap: 10px; }
                .mol-preview-mini { background: #fff; width: 44px; height: 44px; border-radius: 6px; padding: 3px; border: 1px solid #333; display: flex; align-items: center; justify-content: center; overflow: hidden; }
                .mol-preview-mini img { max-width: 100%; max-height: 100%; object-fit: contain; }
                .prop-dots { display: flex; gap: 10px; font-size: 10px; font-weight: bold; }
                .score-cell .score-val { color: #EA4335; font-weight: 800; font-size: 18px; }

                .card-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 15px; }
                .cand-card { background: #111; border: 1px solid #222; border-radius: 10px; overflow: hidden; transition: transform 0.2s; }
                .cand-card:hover { transform: translateY(-5px); border-color: #EA4335; }
                .cand-img { background: #fff; height: 160px; display: flex; align-items: center; justify-content: center; padding: 15px; border-bottom: 1px solid #222; }
                .cand-img img { width: 100%; height: 100%; object-fit: contain; }
                .cand-info { padding: 15px; }
                .cand-title { font-size: 11px; font-weight: bold; color: #EA4335; margin-bottom: 6px; font-family: monospace; }
                .cand-score { font-size: 18px; font-weight: 900; color: #fff; margin-bottom: 12px; }

                .interaction-card { grid-column: span 2; background: #0a0a0a; border: 1px solid #1a1a1a; padding: 20px; border-radius: 12px; border-top: 4px solid #EA4335; }
                .interaction-list { display: grid; grid-template-columns: repeat(auto-fill, minmax(200px, 1fr)); gap: 10px; margin-top: 15px; }
                .interaction-row { background: #111; padding: 8px 12px; border-radius: 6px; display: flex; justify-content: space-between; align-items: center; border: 1px solid #222; font-size: 11px; }
                .interaction-row.key-res { border-left: 2px solid #FF8042; }
                .type-tag { font-weight: 800; font-size: 9px; padding: 2px 6px; border-radius: 4px; text-transform: uppercase; }
                .type-tag.hydrophobic { background: rgba(122, 122, 122, 0.2); color: #aaa; }
                .type-tag.h-bond { background: rgba(76, 175, 80, 0.2); color: #4CAF50; }
                .type-tag.salt-bridge { background: rgba(255, 128, 66, 0.2); color: #FF8042; }
                .type-tag.contact { background: rgba(234, 67, 53, 0.2); color: #EA4335; }
                .interaction-logic { margin-top: 20px; padding-top: 15px; border-top: 1px solid #222; }
                .interaction-logic h5 { margin-bottom: 8px; color: #EA4335; font-size: 12px; }
                .interaction-logic p { font-size: 13px; color: #bbb; line-height: 1.5; font-style: italic; }

                .view-header-flex { display: flex; justify-content: space-between; align-items: center; margin-bottom: 20px; border-bottom: 1px solid #1a1a1a; padding-bottom: 15px; }
                .route-switcher { display: flex; align-items: center; gap: 12px; }
                .switcher-label { font-size: 10px; font-weight: bold; color: #666; text-transform: uppercase; }
                .switcher-buttons { display: flex; gap: 5px; background: #0a0a0a; padding: 4px; border-radius: 8px; border: 1px solid #222; }
                .switcher-btn { background: none; border: none; color: #666; padding: 6px 12px; border-radius: 6px; font-size: 11px; font-weight: 600; cursor: pointer; transition: all 0.2s; }
                .switcher-btn:hover { background: #111; color: #fff; }
                .switcher-btn.active { background: #4285F4; color: #fff; box-shadow: 0 2px 8px rgba(66, 133, 244, 0.3); }
                .switcher-btn .conf { font-size: 9px; opacity: 0.7; margin-left: 4px; }

                .retro-flow-vertical { display: flex; flex-direction: column; gap: 30px; }
                .route-block-modern { background: #0d0d0d; border-radius: 12px; border: 1px solid #1a1a1a; overflow: hidden; }
                .route-header-modern { display: flex; justify-content: space-between; align-items: center; padding: 12px 20px; background: #111; border-bottom: 1px solid #1a1a1a; }
                .route-index { font-size: 10px; font-weight: 900; letter-spacing: 1px; color: #888; }
                .route-score-pill { background: rgba(76, 175, 80, 0.1); color: #4CAF50; padding: 3px 10px; border-radius: 20px; font-size: 10px; font-weight: bold; border: 1px solid rgba(76, 175, 80, 0.2); }

                .route-visual-container { padding: 20px; background: #000; display: flex; justify-content: center; }
                .full-route-image-box { width: 100%; max-width: 1200px; text-align: center; }
                .full-route-img { width: 100%; height: auto; display: block; filter: invert(0.9) hue-rotate(180deg) brightness(1.2); border-radius: 8px; transition: transform 0.3s; }
                .img-caption { font-size: 9px; color: #444; margin-top: 10px; text-transform: uppercase; letter-spacing: 1px; font-weight: bold; }

                .route-details-modern { padding: 25px; border-top: 1px solid #1a1a1a; }
                .route-details-modern h5 { font-size: 11px; text-transform: uppercase; color: #4285F4; margin-bottom: 15px; letter-spacing: 1px; }
                .retro-details-table { width: 100%; border-collapse: collapse; font-size: 11px; }
                .retro-details-table th { text-align: left; color: #666; padding: 10px; border-bottom: 2px solid #1a1a1a; }
                .retro-details-table td { padding: 12px 10px; border-bottom: 1px solid #111; color: #ccc; }
                .rxn-column { font-weight: bold; color: #fff; }
                .smmi-list { display: flex; flex-wrap: wrap; gap: 4px; }
                .smmi-list code { background: #1a1a1a; padding: 2px 6px; border-radius: 4px; font-family: monospace; color: #888; border: 1px solid #222; }
                .sm-empty { color: #555; font-style: italic; }

                .route-summary-stats { display: flex; gap: 15px; margin-top: 20px; }
                .stat-pill { background: #0a0a0a; border: 1px solid #222; padding: 8px 16px; border-radius: 8px; font-size: 10px; color: #888; }
                .stat-pill b { color: #fff; margin-right: 5px; }

                .retro-placeholder-icon { font-size: 48px; color: #1a1a1a; margin-bottom: 15px; }
                .retro-empty-state { text-align: center; padding: 60px; color: #444; }

                @keyframes spin { 0% { transform: rotate(0deg); } 100% { transform: rotate(360deg); } }
                .spinner-small { width: 24px; height: 24px; border: 3px solid #1a1a1a; border-top: 3px solid #4285F4; border-radius: 50%; animation: spin 0.8s linear infinite; margin: 0 auto; }
            `}</style>
        </div>
    );
};

export default LiveWorkspace;
