import React, { useState, useEffect, useRef } from 'react';
import { useParams, useNavigate } from 'react-router-dom';
import { LineChart, Line, XAxis, YAxis, Tooltip, ResponsiveContainer, CartesianGrid } from 'recharts';

const FeudalOptimization = () => {
    const { projectId } = useParams();
    const navigate = useNavigate();
    const logEndRef = useRef(null); 
    
    // State
    const [smiles, setSmiles] = useState("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O");
    const [loading, setLoading] = useState(false);
    const [history, setHistory] = useState([]); 
    const [winner, setWinner] = useState(null);
    const [logs, setLogs] = useState([]); 
    // ADDED similarity to weights
    const [weights, setWeights] = useState({ affinity: 1.0, qed: 1.0, sa: 1.0, similarity: 1.0 });

    useEffect(() => {
        logEndRef.current?.scrollIntoView({ behavior: "smooth" });
    }, [logs]);

    const addLog = (msg, type = 'info') => {
        const timestamp = new Date().toLocaleTimeString();
        setLogs(prev => [...prev, { msg, type, time: timestamp }]);
    };

    const runOptimization = async () => {
        setLoading(true);
        setLogs([]);
        addLog("Initializing Feudal RL Engine...", "info");

        try {
            await fetch('http://localhost:8000/api/feudal/optimize', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ project_id: projectId, smiles: smiles })
            });
            addLog("Optimization task sent to background worker.", "success");
            startPolling();
        } catch (e) {
            addLog("FAILED to connect to backend.", "error");
            setLoading(false);
        }
    };

    const startPolling = () => {
        const interval = setInterval(async () => {
            try {
                // 1. Fetch live logs and history
                const statusRes = await fetch(`http://localhost:8000/api/feudal/status/${projectId}`);
                const statusData = await statusRes.json();
                
                if (statusData.logs) setLogs(statusData.logs);
                if (statusData.history) setHistory(statusData.history);

                // 2. Fetch Pareto Winner (for the card)
                const winnerRes = await fetch(`http://localhost:8000/api/feudal/pareto-rank`, {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ project_id: projectId, weights: weights })
                });
                
                if (winnerRes.ok) {
                    const winnerData = await winnerRes.json();
                    setWinner(winnerData);
                    
                    // Stop polling if the backend log says complete
                    const lastLog = statusData.logs[statusData.logs.length - 1];
                    if (lastLog?.msg.includes("complete")) {
                        setLoading(false);
                        clearInterval(interval);
                    }
                }
            } catch (e) { console.error("Polling error", e); }
        }, 3000);
    };

    return (
        <div className="main-content">
            <div className="docking-container-wrapper" style={{ transform: 'translateX(-60px)' }}>
                
                <div className="docking-header">
                    <button className="primary-btn" onClick={() => navigate(`/dashboard/${projectId}`)}>← Back</button>
                    <h2 style={{ margin: 0 }}>Feudal RL Optimization: {projectId}</h2>
                </div>

                <div className="docking-grid">
                    {/* LEFT PANEL */}
                    <div className="card" style={{ margin: 0, padding: '24px', display: 'flex', flexDirection: 'column', gap: '20px' }}>
                        <div>
                            <h3 style={{ fontSize: '0.7rem', color: 'var(--text-dim)', marginBottom: '10px' }}>GENERATIVE SETTINGS</h3>
                            <input type="text" className="technical-input" value={smiles} onChange={(e) => setSmiles(e.target.value)} />
                        </div>

                        <div>
                            <h3 style={{ fontSize: '0.7rem', color: 'var(--text-dim)', marginBottom: '10px' }}>PARETO WEIGHTS</h3>
                            <div className="weight-sliders" style={{ display: 'flex', flexDirection: 'column', gap: '12px' }}>
                                {['affinity', 'qed', 'sa', 'similarity'].map(key => (
                                    <div key={key}>
                                        <div style={{ display: 'flex', justifyContent: 'space-between', fontSize: '0.75rem', marginBottom: '4px' }}>
                                            <span style={{ textTransform: 'uppercase', color: '#ccc' }}>{key}</span>
                                            <span style={{ color: 'var(--primary)' }}>{weights[key]}</span>
                                        </div>
                                        <input 
                                            type="range" min="0" max="2" step="0.1" 
                                            style={{ width: '100%', accentColor: 'var(--primary)' }}
                                            value={weights[key]} 
                                            onChange={(e) => setWeights({...weights, [key]: parseFloat(e.target.value)})}
                                        />
                                    </div>
                                ))}
                            </div>
                        </div>

                        <button className="primary-btn" style={{ background: 'var(--primary)', color: 'white', border: 'none' }} 
                                onClick={runOptimization} disabled={loading}>
                            {loading ? "Optimizing..." : "Start Feudal Loop"}
                        </button>

                        <div>
                            <h3 style={{ fontSize: '0.7rem', color: 'var(--text-dim)', marginBottom: '10px' }}>SYSTEM LOGS</h3>
                            <div className="feudal-log" style={{ height: '180px', fontSize: '0.7rem', overflowY: 'auto', background: '#000', border: '1px solid #222' }}>
                                {logs.map((log, i) => (
                                    <div key={i} style={{ marginBottom: '4px' }}>
                                        <span style={{ color: '#555' }}>[{log.time}]</span>{' '}
                                        <span style={{ color: log.type === 'error' ? '#f87171' : log.type === 'success' ? '#4ade80' : '#3b82f6' }}>{log.msg}</span>
                                    </div>
                                ))}
                                <div ref={logEndRef} />
                            </div>
                        </div>
                    </div>

                    {/* RIGHT PANEL */}
                    <div style={{ display: 'flex', flexDirection: 'column', gap: '24px' }}>
                        <div className="card" style={{ margin: 0, height: '300px' }}>
                            <h3 style={{ fontSize: '0.8rem', color: 'var(--text-dim)', marginBottom: '15px' }}>BINDING AFFINITY TRAJECTORY</h3>
                            <ResponsiveContainer width="100%" height="90%">
                                <LineChart data={history?.length > 0 ? history : [{step:0, score:0}]}>
                                    <CartesianGrid strokeDasharray="3 3" stroke="#222" />
                                    <XAxis dataKey="step" stroke="#71717a" />
                                    <YAxis reversed stroke="#71717a" domain={['auto', 'auto']} />
                                    <Tooltip contentStyle={{ background: '#0f0f0f', border: '1px solid #353535' }} />
                                    <Line type="monotone" dataKey="score" stroke="#3b82f6" strokeWidth={3} />
                                </LineChart>
                            </ResponsiveContainer>
                        </div>

                        <div className="card" style={{ margin: 0, flex: 1, padding: '24px' }}>
                            <h3 style={{ fontSize: '0.8rem', color: 'var(--text-dim)', marginBottom: '15px' }}>PARETO OPTIMAL CANDIDATE</h3>
                            {winner ? (
                                <div style={{ display: 'flex', gap: '40px', alignItems: 'center' }}>
                                    <div style={{ width: '250px', height: '250px', background: '#fff', borderRadius: '8px', padding: '10px', display: 'flex', alignItems: 'center', justifyContent: 'center' }}>
                                        <img 
                                            src={`http://localhost:8000/api/render?smiles=${encodeURIComponent(winner.smiles)}`} 
                                            alt="Winner Molecule" 
                                            style={{ maxWidth: '100%', maxHeight: '100%' }} 
                                        />
                                    </div>
                                    <div style={{ flex: 1 }}>
                                        <div style={{ fontSize: '2rem', fontWeight: 'bold', color: winner.score < 0 ? '#4ade80' : '#f87171' }}>
                                            {winner.score} kcal/mol
                                            {winner.score === 0 && <span style={{fontSize: '0.8rem', display: 'block', color: '#f87171'}}>⚠️ Vina Execution Failed</span>}
                                        </div>
                                        
                                        <div className="feudal-log" style={{ height: 'auto', width: '350px', background: 'transparent', border: 'none', padding: 0, marginTop: '15px' }}>
                                            <div style={{display: 'flex', justifyContent: 'space-between', marginBottom: '8px', fontSize: '1rem'}}>
                                                <span>QED Score:</span> <strong>{winner.qed}</strong>
                                            </div>
                                            <div style={{display: 'flex', justifyContent: 'space-between', marginBottom: '8px', fontSize: '1rem'}}>
                                                <span>SA Score:</span> <strong>{winner.sa}</strong>
                                            </div>
                                            <div style={{display: 'flex', justifyContent: 'space-between', marginBottom: '8px', fontSize: '1rem'}}>
                                                <span>Similarity:</span> <strong style={{color: 'var(--primary)'}}>{(winner.similarity * 100).toFixed(1)}%</strong>
                                            </div>
                                            <div style={{borderTop: '1px solid #333', marginTop: '10px', paddingTop: '10px', fontSize: '0.8rem', color: '#666'}}>
                                                SMILES: <span style={{wordBreak: 'break-all'}}>{winner.smiles}</span>
                                            </div>
                                        </div>
                                    </div>
                                </div>
                            ) : (
                                <div className="stats-placeholder">Waiting for evolution data...</div>
                            )}
                        </div>
                    </div>
                </div>
            </div>
        </div>
    );
};

export default FeudalOptimization;