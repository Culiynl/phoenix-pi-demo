import React, { useState } from 'react';
import MolstarViewer from './MolstarViewer';

const StandaloneTools = () => {
    const [activeTab, setActiveTab] = useState('docking');
    const [dockingInput, setDockingInput] = useState({ pdb: '', smiles: '' });
    const [dockingResult, setDockingResult] = useState(null);
    const [loading, setLoading] = useState(false);

    const [pdbSearchResults, setPdbSearchResults] = useState(null);

    const runPdbSearch = async () => {
        if (!dockingInput.pdb) return;
        setLoading(true);
        try {
            // const res = await fetch(`http://localhost:8000/api/pdb/search?q=${encodeURIComponent(dockingInput.pdb)}`);
            const res = await fetch(`http://localhost:8000/api/pdb/search?q=${encodeURIComponent(dockingInput.pdb)}`);
            const data = await res.json();
            setPdbSearchResults(data);
        } catch (e) {
            console.error(e);
            alert("Search failed");
        } finally {
            setLoading(false);
        }
    };

    const selectPdb = (id) => {
        // We need to communicate up to the parent viewer, but StandaloneTools is currently standalone.
        // The user asked to "add mol star visualization and also maybe pair it with RCSB".
        // StandaloneTools usually is just a side panel. 
        // If we want to Load it, we might need a callback prop 'onLoadPdb'.
        // For now, let's just pre-fill the docking input or alert.
        // Actually, looking at the user request: "tools section you need to add the mol star visualization"
        // If the Tools section is separate, it might need its own viewer?
        // Or if this component is inside LiveWorkspace, it might have access to setMissionData?
        // Based on file context, StandaloneTools seems unrelated to LiveWorkspace state.
        // I will assume for now we just prefill for docking or copy to clipboard, 
        // OR add a MolstarViewer right here in the tools panel?
        // User said: "tools section you need to add the mol star visualization"
        // Let's notify the user via a callback if possible, or add a viewer here.
        // Let's add an activePdb state here and show a viewer if selected.
        setDockingInput({ ...dockingInput, pdb: id });
        setActivePdb(id);
    };

    const [activePdb, setActivePdb] = useState(null);

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
            <h3>Standalone Tools</h3>
            <p className="tools-desc">Quick analysis tools powered by AutoDock Vina, RDKit, and more</p>

            <div className="tool-tabs">
                <button
                    className={activeTab === 'docking' ? 'active' : ''}
                    onClick={() => setActiveTab('docking')}
                >
                    Molecular Docking
                </button>
                <button
                    className={activeTab === 'proteins' ? 'active' : ''}
                    onClick={() => setActiveTab('proteins')}
                >
                    Protein Search
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

            {activeTab === 'proteins' && (
                <div className="tool-panel">
                    <h4>RCSB PDB Search</h4>
                    {/* ... search inputs ... */}

                    {activePdb && (
                        <div style={{ height: '400px', margin: '20px 0', border: '1px solid #333', borderRadius: '8px', overflow: 'hidden' }}>
                            <MolstarViewer pdbId={activePdb} />
                        </div>
                    )}

                    <div className="search-bar-row">
                        <input
                            type="text"
                            placeholder="Protein Name (e.g. EGFR, Kinase)"
                            value={dockingInput.pdb}
                            onChange={(e) => setDockingInput({ ...dockingInput, pdb: e.target.value })}
                            onKeyDown={(e) => e.key === 'Enter' && runPdbSearch()}
                        />
                        <button onClick={runPdbSearch} disabled={loading} style={{ width: '80px', marginBottom: '10px' }}>
                            {loading ? '...' : 'Search'}
                        </button>
                    </div>

                    {pdbSearchResults && (
                        <div className="result-list-scroll">
                            {pdbSearchResults.length === 0 && <div className="info-text">No results found.</div>}
                            {pdbSearchResults.map((res) => (
                                <div key={res.id} className="pdb-result-item" onClick={() => selectPdb(res.id)}>
                                    <div className="pdb-id-badge">{res.id}</div>
                                    <div className="pdb-info">
                                        <div className="pdb-title">{res.title}</div>
                                        <div className="pdb-meta">{res.organism} • {res.resolution || "N/A"}Å • {res.date}</div>
                                    </div>
                                    <button className="load-btn">LOAD</button>
                                </div>
                            ))}
                        </div>
                    )}
                </div>
            )}

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
                    color: #ffffffff;
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
                    background: #3b82f6;
                    color: #ffffffff;
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
                    background: #3b82f6;
                    border: none;
                    border-radius: 6px;
                    padding: 10px;
                    color: #ffffffff;
                    font-weight: bold;
                    font-size: 12px;
                    cursor: pointer;
                    transition: background 0.2s;
                }
                .tool-panel button:hover:not(:disabled) {
                    background: #3b83f677;
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
                
                .result-list-scroll { max-height: 300px; overflow-y: auto; display: flex; flex-direction: column; gap: 8px; margin-top: 10px; }
                .pdb-result-item { display: flex; align-items: center; gap: 10px; background: #0a0a0a; border: 1px solid #222; padding: 10px; border-radius: 8px; cursor: pointer; transition: all 0.2s; }
                .pdb-result-item:hover { background: #151515; border-color: #4285F4; }
                .pdb-id-badge { background: #4285F4; color: #fff; font-weight: 900; font-size: 12px; padding: 4px 8px; border-radius: 4px; width: 50px; text-align: center; }
                .pdb-info { flex-grow: 1; overflow: hidden; }
                .pdb-title { font-size: 11px; font-weight: bold; color: #ddd; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
                .pdb-meta { font-size: 10px; color: #666; margin-top: 2px; }
                .load-btn { background: none; border: 1px solid #333; color: #4285F4; font-size: 9px; padding: 4px 8px; width: auto; }
                .search-bar-row { display: flex; gap: 10px; }
            `}</style>
        </div>
    );
};

export default StandaloneTools;
