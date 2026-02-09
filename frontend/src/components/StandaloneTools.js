import React, { useState } from 'react';

const StandaloneTools = () => {
    const [activeTab, setActiveTab] = useState('docking');
    const [dockingInput, setDockingInput] = useState({ pdb: '', smiles: '' });
    const [dockingResult, setDockingResult] = useState(null);
    const [loading, setLoading] = useState(false);

    const runDocking = async () => {
        setLoading(true);
        try {
            const response = await fetch('http://localhost:8000/api/docking', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    pdb_id: dockingInput.pdb,
                    ligand_smiles: dockingInput.smiles
                })
            });
            const result = await response.json();
            setDockingResult(result);
        } catch (error) {
            setDockingResult({ error: error.message });
        } finally {
            setLoading(false);
        }
    };

    return (
        <div className="standalone-tools">
            <h3>🔬 Standalone Tools</h3>
            <p className="tools-desc">Quick analysis tools powered by AutoDock Vina, RDKit, and more</p>

            <div className="tool-tabs">
                <button
                    className={activeTab === 'docking' ? 'active' : ''}
                    onClick={() => setActiveTab('docking')}
                >
                    Molecular Docking
                </button>
                <button
                    className={activeTab === 'properties' ? 'active' : ''}
                    onClick={() => setActiveTab('properties')}
                >
                    Molecule Properties
                </button>
                <button
                    className={activeTab === 'similarity' ? 'active' : ''}
                    onClick={() => setActiveTab('similarity')}
                >
                    Similarity Search
                </button>
            </div>

            {activeTab === 'docking' && (
                <div className="tool-panel">
                    <h4>AutoDock Vina Docking</h4>
                    <input
                        type="text"
                        placeholder="PDB ID (e.g., 1A2B)"
                        value={dockingInput.pdb}
                        onChange={(e) => setDockingInput({ ...dockingInput, pdb: e.target.value })}
                    />
                    <input
                        type="text"
                        placeholder="Ligand SMILES"
                        value={dockingInput.smiles}
                        onChange={(e) => setDockingInput({ ...dockingInput, smiles: e.target.value })}
                    />
                    <button onClick={runDocking} disabled={loading}>
                        {loading ? 'Docking...' : 'Run Docking'}
                    </button>

                    {dockingResult && (
                        <div className="result-box">
                            {dockingResult.error ? (
                                <div className="error">{dockingResult.error}</div>
                            ) : (
                                <div className="success">
                                    <strong>Binding Affinity:</strong> {dockingResult.affinity} kcal/mol
                                    <br />
                                    <strong>Best Pose:</strong> {dockingResult.pose}
                                </div>
                            )}
                        </div>
                    )}
                </div>
            )}

            {activeTab === 'properties' && (
                <div className="tool-panel">
                    <h4>Calculate Molecular Properties</h4>
                    <input type="text" placeholder="Enter SMILES" />
                    <button>Calculate</button>
                    <div className="info-text">
                        Calculates: MW, LogP, TPSA, QED, Lipinski violations
                    </div>
                </div>
            )}

            {activeTab === 'similarity' && (
                <div className="tool-panel">
                    <h4>Tanimoto Similarity Search</h4>
                    <input type="text" placeholder="Query SMILES" />
                    <input type="text" placeholder="Target SMILES (comma-separated)" />
                    <button>Search</button>
                    <div className="info-text">
                        Finds similar molecules using fingerprint similarity
                    </div>
                </div>
            )}

            <style>{`
                .standalone-tools {
                    background: #0d0d0d;
                    border: 1px solid #1a1a1a;
                    border-radius: 12px;
                    padding: 20px;
                    margin-bottom: 20px;
                }
                .standalone-tools h3 {
                    margin: 0 0 8px 0;
                    color: #F4B400;
                    font-size: 16px;
                }
                .tools-desc {
                    color: #888;
                    font-size: 12px;
                    margin: 0 0 15px 0;
                }
                .tool-tabs {
                    display: flex;
                    gap: 8px;
                    margin-bottom: 15px;
                    border-bottom: 1px solid #222;
                    padding-bottom: 10px;
                }
                .tool-tabs button {
                    background: none;
                    border: none;
                    color: #888;
                    padding: 8px 12px;
                    border-radius: 6px;
                    cursor: pointer;
                    font-size: 12px;
                    transition: all 0.2s;
                }
                .tool-tabs button:hover {
                    background: #1a1a1a;
                    color: #fff;
                }
                .tool-tabs button.active {
                    background: #F4B400;
                    color: #000;
                    font-weight: bold;
                }
                .tool-panel {
                    background: #111;
                    border-radius: 8px;
                    padding: 15px;
                }
                .tool-panel h4 {
                    margin: 0 0 12px 0;
                    color: #fff;
                    font-size: 13px;
                }
                .tool-panel input {
                    width: 100%;
                    background: #0a0a0a;
                    border: 1px solid #222;
                    border-radius: 6px;
                    padding: 10px;
                    color: #fff;
                    font-size: 12px;
                    margin-bottom: 10px;
                }
                .tool-panel input::placeholder {
                    color: #555;
                }
                .tool-panel button {
                    width: 100%;
                    background: #F4B400;
                    border: none;
                    border-radius: 6px;
                    padding: 10px;
                    color: #000;
                    font-weight: bold;
                    font-size: 12px;
                    cursor: pointer;
                    transition: background 0.2s;
                }
                .tool-panel button:hover:not(:disabled) {
                    background: #ffcc00;
                }
                .tool-panel button:disabled {
                    opacity: 0.5;
                    cursor: not-allowed;
                }
                .result-box {
                    margin-top: 12px;
                    padding: 12px;
                    border-radius: 6px;
                    font-size: 12px;
                }
                .result-box .success {
                    background: rgba(76, 175, 80, 0.1);
                    border: 1px solid #4CAF50;
                    color: #4CAF50;
                }
                .result-box .error {
                    background: rgba(244, 67, 54, 0.1);
                    border: 1px solid #f44336;
                    color: #f44336;
                }
                .info-text {
                    margin-top: 10px;
                    font-size: 11px;
                    color: #666;
                    font-style: italic;
                }
            `}</style>
        </div>
    );
};

export default StandaloneTools;
