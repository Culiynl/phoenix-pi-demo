import React, { useState, useEffect, useRef } from 'react';
import { useParams, useNavigate } from 'react-router-dom';
import { LineChart, Line, XAxis, YAxis, Tooltip, ResponsiveContainer, CartesianGrid } from 'recharts';

const FeudalOptimization = () => {
    const { projectId } = useParams();
    const navigate = useNavigate();
    const logEndRef = useRef(null); // For auto-scrolling
    
    // State
    const [smiles, setSmiles] = useState("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O");
    const [loading, setLoading] = useState(false);
    const [history, setHistory] = useState([]); 
    const [winner, setWinner] = useState(null);
    const [logs, setLogs] = useState([]); // NEW: System Logs
    const [weights, setWeights] = useState({ affinity: 1.0, qed: 1.0, sa: 1.0 });

    // Auto-scroll logs to bottom
    useEffect(() => {
        logEndRef.current?.scrollIntoView({ behavior: "smooth" });
    }, [logs]);

    const addLog = (msg, type = 'info') => {
        const timestamp = new Date().toLocaleTimeString();
        setLogs(prev => [...prev, { msg, type, time: timestamp }]);
    };

    const runOptimization = async () => {
        setLoading(true);
        setLogs([]); // Clear old logs
        addLog("Initializing Feudal RL Engine...", "info");
        addLog(`Target Project: ${projectId}`, "info");

        try {
            const res = await fetch('http://localhost:8000/api/feudal/optimize', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ project_id: projectId, smiles: smiles })
            });
            const data = await res.json();
            
            addLog("Optimization task sent to background worker.", "success");
            addLog("Starting live log stream...", "warn");

            // Start Polling for results/logs
            startPolling();
        } catch (e) {
            addLog("FAILED to connect to backend engine.", "error");
        }
    };

    const startPolling = () => {
        const interval = setInterval(async () => {
            try {
                // We fetch the latest history/logs from the backend
                const res = await fetch(`http://localhost:8000/api/feudal/pareto-rank`, {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ project_id: projectId, weights: weights })
                });
                
                if (res.ok) {
                    const winnerData = await res.json();
                    setWinner(winnerData);
                    addLog(`Update: Best score found: ${winnerData.score} kcal/mol`, "success");
                    
                    // If we found a winner, we can stop polling (or keep going until steps are done)
                    // For now, let's stop loading once we get a valid response
                    setLoading(false);
                    clearInterval(interval);
                }
            } catch (e) {
                // Continue polling until backend writes the file
            }
        }, 3000);
    };

    return (
        <div className="main-content">
            <div className="docking-container-wrapper" style={{ transform: 'translateX(-40px)' }}>
                
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
                                <div className="weight-sliders" style={{ display: 'flex', flexDirection: 'column', gap: '12px', marginBottom: '20px' }}>
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

                        {/* LOGGING CONSOLE */}
                        <div>
                            <h3 style={{ fontSize: '0.7rem', color: 'var(--text-dim)', marginBottom: '10px' }}>SYSTEM LOGS</h3>
                            <div className="feudal-log" style={{ height: '200px', fontSize: '0.7rem', overflowY: 'auto', background: '#000', border: '1px solid #222' }}>
                                {logs.length === 0 && <div style={{color: '#444'}}>No active tasks...</div>}
                                {logs.map((log, i) => (
                                    <div key={i} style={{ marginBottom: '4px' }}>
                                        <span style={{ color: '#555' }}>[{log.time}]</span>{' '}
                                        <span style={{ color: log.type === 'error' ? '#f87171' : log.type === 'success' ? '#4ade80' : '#3b82f6' }}>
                                            {log.msg}
                                        </span>
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

                        <div className="winner-stats">
    <div style={{ fontSize: '1.8rem', fontWeight: 'bold', color: winner?.score < 0 ? '#4ade80' : '#f87171', marginBottom: '10px' }}>
        {winner?.score} kcal/mol
        {winner?.score === 0 && <span style={{fontSize: '0.8rem', display: 'block', color: '#f87171'}}>⚠️ Docking Failed</span>}
    </div>
    
    <div className="feudal-log" style={{ height: 'auto', width: '300px', padding: '15px' }}>
        <div style={{display: 'flex', justifyContent: 'space-between', marginBottom: '5px'}}>
            <span>QED Score:</span> <strong>{winner?.qed}</strong>
        </div>
        <div style={{display: 'flex', justifyContent: 'space-between', marginBottom: '5px'}}>
            <span>SA Score:</span> <strong>{winner?.sa}</strong>
        </div>
        <div style={{display: 'flex', justifyContent: 'space-between', marginBottom: '5px'}}>
            <span>Similarity:</span> <strong style={{color: 'var(--primary)'}}>{(winner?.similarity * 100).toFixed(1)}%</strong>
        </div>
        <div style={{borderTop: '1px solid #222', marginTop: '10px', paddingTop: '10px', fontSize: '0.7rem', color: '#888'}}>
            Step Identified: {winner?.step}
        </div>
    </div>
</div>
                    </div>
                </div>
            </div>
        </div>
    );
};

export default FeudalOptimization;