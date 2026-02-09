import React, { useState } from 'react';
import { LineChart, Line, XAxis, YAxis, Tooltip, ResponsiveContainer } from 'recharts';

const FeudalDashboard = ({ projectId }) => {
    const [weights, setWeights] = useState({ affinity: 1.0, qed: 1.0, sa: 1.0 });
    const [history, setHistory] = useState([]);
    const [winner, setWinner] = useState(null);

    const runOptimization = async () => {
        const res = await fetch('/api/feudal/optimize', {
            method: 'POST',
            body: JSON.stringify({ project_id: projectId, smiles: "..." })
        });
        const data = await res.json();
        setHistory(data.history);
    };

    return (
        <div className="feudal-panel">
            <div className="controls card">
                <h3>Pareto Weights</h3>
                <label>Affinity: <input type="range" step="0.1" min="0" max="2" 
                    value={weights.affinity} onChange={e => setWeights({...weights, affinity: e.target.value})} /></label>
                <label>QED: <input type="range" step="0.1" min="0" max="2" 
                    value={weights.qed} onChange={e => setWeights({...weights, qed: e.target.value})} /></label>
                <button className="primary-btn" onClick={runOptimization}>Start Feudal Evolution</button>
            </div>

            <div className="grid-2">
                <div className="card">
                    <h4>Optimization Trajectory (Vina Score)</h4>
                    <ResponsiveContainer width="100%" height={300}>
                        <LineChart data={history}>
                            <XAxis dataKey="step" />
                            <YAxis reversed />
                            <Tooltip />
                            <Line type="monotone" dataKey="score" stroke="#3b82f6" strokeWidth={3} />
                        </LineChart>
                    </ResponsiveContainer>
                </div>

                <div className="card">
                    <h4>Manager Goal (Stepwise Vector)</h4>
                    <pre className="feudal-log">
                        {history.map(h => `Step ${h.step}: [${h.goal?.map(g => g.toFixed(2))}]`).join('\n')}
                    </pre>
                </div>
            </div>

            {winner && (
                <div className="card winner-display">
                    <h3>🏆 Pareto Winner</h3>
                    <div className="mol-info">
                        <img src={`/api/render/${winner.smiles}`} alt="molecule" />
                        <div>
                            <p>Affinity: {winner.score}</p>
                            <p>QED: {winner.qed}</p>
                            <p>SA Score: {winner.sa}</p>
                        </div>
                    </div>
                </div>
            )}
        </div>
    );
};