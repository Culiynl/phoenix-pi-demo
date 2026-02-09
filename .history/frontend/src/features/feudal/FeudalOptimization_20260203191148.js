import React, { useState, useEffect } from 'react';
import { useParams, useNavigate } from 'react-router-dom';
import { LineChart, Line, XAxis, YAxis, Tooltip, ResponsiveContainer, CartesianGrid } from 'recharts';

const FeudalOptimization = () => {
    const { projectId } = useParams();
    const navigate = useNavigate();
    
    // Safety: Always initialize history as an empty array []
    const [smiles, setSmiles] = useState("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O");
    const [loading, setLoading] = useState(false);
    const [history, setHistory] = useState([]); 
    const [winner, setWinner] = useState(null);
    const [weights, setWeights] = useState({ affinity: 1.0, qed: 1.0, sa: 1.0 });

    const runOptimization = async () => {
        setLoading(true);
        try {
            const res = await fetch('http://localhost:8000/api/feudal/optimize', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ project_id: projectId, smiles: smiles })
            });
            const data = await res.json();

            // Safety Fix 1: Use a fallback if data.history doesn't exist yet
            const newHistory = data.history || [];
            setHistory(newHistory);
            
            // Safety Fix 2: Only call winner calculation if we actually got history
            if (newHistory.length > 0) {
                calculateWinner(newHistory);
            } else {
                alert("Optimization started in background. Please wait a moment and refresh/re-rank.");
            }
        } catch (e) {
            console.error("Optimization error", e);
            alert("Optimization Task Failed");
        }
        setLoading(false);
    };

    const calculateWinner = async (overrideHistory) => {
        try {
            const res = await fetch('http://localhost:8000/api/feudal/pareto-rank', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ project_id: projectId, weights: weights })
            });
            const data = await res.json();
            setWinner(data);
        } catch (e) {
            console.error("Pareto ranking failed", e);
        }
    };

    return (
        <div className="main-content">
            <div className="docking-container-wrapper" style={{ transform: 'translateX(-40px)' }}>
                
                <div className="docking-header">
                    <button className="primary-btn" onClick={() => navigate(`/dashboard/${projectId}`)}>
                        ← Back
                    </button>
                    <h2 style={{ margin: 0 }}>Feudal RL Optimization: {projectId}</h2>
                    <div style={{ width: '100px' }} />
                </div>

                <div className="docking-grid">
                    <div className="card" style={{ margin: 0, padding: '24px', display: 'flex', flexDirection: 'column' }}>
                        <h3 style={{ fontSize: '0.7rem', color: 'var(--text-dim)', marginBottom: '15px' }}>GENERATIVE SETTINGS</h3>
                        
                        <label style={{ fontSize: '0.7rem', color: 'var(--text-dim)' }}>INITIAL LEAD (SMILES)</label>
                        <input 
                            type="text" 
                            className="technical-input" 
                            style={{ marginBottom: '20px', marginTop: '5px' }}
                            value={smiles}
                            onChange={(e) => setSmiles(e.target.value)}
                        />

                        <h3 style={{ fontSize: '0.7rem', color: 'var(--text-dim)', marginBottom: '10px' }}>PARETO WEIGHTS</h3>
                        <div className="weight-sliders" style={{ display: 'flex', flexDirection: 'column', gap: '15px', marginBottom: '20px' }}>
                            {['affinity', 'qed', 'sa'].map(key => (
                                <div key={key}>
                                    <div style={{ display: 'flex', justifyContent: 'space-between', fontSize: '0.8rem' }}>
                                        <span style={{ textTransform: 'uppercase' }}>{key}</span>
                                        <span>{weights[key]}</span>
                                    </div>
                                    <input 
                                        type="range" min="0" max="2" step="0.1" 
                                        style={{ width: '100%' }}
                                        value={weights[key]} 
                                        onChange={(e) => setWeights({...weights, [key]: parseFloat(e.target.value)})}
                                    />
                                </div>
                            ))}
                        </div>

                        <button className="primary-btn" style={{ width: '100%', background: 'var(--primary)', color: 'white', border: 'none' }} 
                                onClick={runOptimization} disabled={loading}>
                            {loading ? "Optimizing (Running Vina)..." : "Start Feudal Loop"}
                        </button>
                        
                        {/* Safety Fix 3: Optional Chaining ?. */}
                        {history?.length > 0 && (
                            <button className="primary-btn" style={{ width: '100%', marginTop: '10px' }} 
                                    onClick={() => calculateWinner()}>
                                Recalculate Pareto Winner
                            </button>
                        )}
                    </div>

                    <div style={{ display: 'flex', flexDirection: 'column', gap: '24px' }}>
                        
                        <div className="card" style={{ margin: 0, height: '350px' }}>
                            <h3 style={{ fontSize: '0.8rem', color: 'var(--text-dim)', marginBottom: '15px' }}>BINDING AFFINITY TRAJECTORY</h3>
                            <ResponsiveContainer width="100%" height="90%">
                                {/* Safety Fix 4: Fallback data array */}
                                <LineChart data={history || []}>
                                    <CartesianGrid strokeDasharray="3 3" stroke="#222" />
                                    <XAxis dataKey="step" stroke="#71717a" />
                                    <YAxis reversed stroke="#71717a" domain={['auto', 'auto']} />
                                    <Tooltip contentStyle={{ background: '#0f0f0f', border: '1px solid #353535' }} />
                                    <Line type="monotone" dataKey="score" stroke="#3b82f6" strokeWidth={3} dot={{ r: 6 }} />
                                </LineChart>
                            </ResponsiveContainer>
                        </div>

                        <div className="card" style={{ margin: 0, flex: 1 }}>
                            <h3 style={{ fontSize: '0.8rem', color: 'var(--text-dim)', marginBottom: '15px' }}>PARETO OPTIMAL CANDIDATE</h3>
                            {winner ? (
                                <div style={{ display: 'flex', gap: '40px', alignItems: 'center' }}>
                                    <div className="stats-placeholder" style={{ width: '250px', height: '250px', background: '#fff' }}>
                                        {/* Safety Fix 5: Optional chaining on smiles */}
                                        <img 
                                            src={`http://localhost:8000/api/render?smiles=${encodeURIComponent(winner?.smiles || "")}`} 
                                            alt="Molecule" 
                                            style={{ width: '100%', height: '100%', objectFit: 'contain' }}
                                        />
                                    </div>
                                    <div className="winner-stats">
                                        <div style={{ fontSize: '1.5rem', fontWeight: 'bold', color: '#4ade80', marginBottom: '10px' }}>
                                            Score: {winner?.score} kcal/mol
                                        </div>
                                        <div className="feudal-log" style={{ height: 'auto', width: '300px' }}>
                                            <div>QED Score: {winner?.qed}</div>
                                            <div>SA Score: {winner?.sa}</div>
                                            <div>Step Identified: {winner?.step}</div>
                                            <div style={{marginTop: '10px', color: 'var(--primary)'}}>SMILES:</div>
                                            <div style={{fontSize: '0.7rem', wordBreak: 'break-all'}}>{winner?.smiles}</div>
                                        </div>
                                    </div>
                                </div>
                            ) : (
                                <div className="stats-placeholder">Run optimization to identify Pareto winner</div>
                            )}
                        </div>
                    </div>
                </div>
            </div>
        </div>
    );
};

export default FeudalOptimization;