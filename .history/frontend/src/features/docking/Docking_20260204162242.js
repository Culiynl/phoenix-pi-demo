import React, { useEffect, useRef, useState } from 'react';
import { useParams, useNavigate } from 'react-router-dom';
import { createPluginUI } from 'molstar/lib/mol-plugin-ui';
import { DefaultPluginUISpec } from 'molstar/lib/mol-plugin-ui/spec';
import { renderReact18 } from 'molstar/lib/mol-plugin-ui/react18';
import 'molstar/lib/mol-plugin-ui/skin/dark.scss';
const highlightPocket = async (residueIdsString) => {
    if (!pluginRef.current) return;
    const plugin = pluginRef.current;

    // residueIdsString looks like "A_107 A_108 A_110"
    const selection = residueIdsString.split(' ').map(id => {
        const [auth_asym_id, auth_seq_id] = id.split('_');
        return { auth_asym_id, auth_seq_id: parseInt(auth_seq_id) };
    });

    // Clear previous highlights and apply new ones
    // This uses Mol* Script to create a selection
    const data = plugin.managers.structure.hierarchy.current.structures[0].cell.obj.data;
    const sel = Script.getStructureSelection(Q => Q.struct.generator.atomGroups({
        'residue-test': Q.core.set.has([
            Q.set(...selection.map(s => s.auth_seq_id)),
            Q.ammp('auth_seq_id')
        ]),
        'chain-test': Q.core.set.has([
            Q.set(...selection.map(s => s.auth_asym_id)),
            Q.ammp('auth_asym_id')
        ])
    }), data);

    plugin.managers.structure.selection.fromSelection(sel);
    plugin.managers.camera.focusSelection(sel);
};
const Docking = () => {
    const { projectId } = useParams();
    const navigate = useNavigate();
    const parentRef = useRef(null);
    const pluginRef = useRef(null);
    
    // UI State
    const [smiles, setSmiles] = useState("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O");
    const [loading, setLoading] = useState(false);
    const [isPluginReady, setIsPluginReady] = useState(false);
    const [pdbSearch, setPdbSearch] = useState(""); // FIX: Define state
    const [pockets, setPockets] = useState([]);      // FIX: Define state
    
    // Docking Results State
    const [poses, setPoses] = useState([]); // Array of { id, energy }
    const [activePose, setActivePose] = useState(0);
    const [lastResultUrls, setLastResultUrls] = useState(null); // Stores { receptor_url, ligand_url }

    // 1. Initialize Mol*
    useEffect(() => {
        if (!parentRef.current || pluginRef.current) return;

        async function init() {
            const spec = DefaultPluginUISpec();
            
            // Force a very clean, dark UI
            spec.layout = {
                initialShowControls: false,
                initialShowRemoteState: false,
                initialShowSequence: false,
                initialShowLog: false,
                initialShowLeftPanel: false,
            };
            spec.components = {
                remoteState: 'none',
                sequence: 'none',
                log: 'none',
                leftPanel: 'none'
            };

            try {
                const ctx = await createPluginUI({
                    target: parentRef.current,
                    spec: spec,
                    render: renderReact18 
                });
                pluginRef.current = ctx;
                setIsPluginReady(true);
            } catch (error) {
                console.error("Mol* Init Error:", error);
            }
        }
        init();

        return () => {
            if (pluginRef.current) {
                pluginRef.current.dispose();
                pluginRef.current = null;
            }
        };
    }, []);

    const handleFetchPdb = async () => {
        if (!pdbSearch) return alert("Enter a PDB ID");
        setLoading(true);
        try {
            const res = await fetch(`http://localhost:8000/api/pdb/fetch/${pdbSearch}`);
            const data = await res.json();
            if (data.status === "success") {
                alert(`Downloaded ${data.filename}. You can now run P2Rank.`);
                // Automatically load it into the viewer
                loadPdbIntoViewer(`/static/proteins/${data.filename}`);
            }
        } catch (e) { alert("Failed to fetch PDB"); }
        setLoading(false);
    };
    const handleRunP2Rank = async () => {
        setLoading(true);
        try {
            // Trigger the backend CLI
            await fetch(`http://localhost:8000/api/p2rank/run?pdb_filename=${pdbSearch}.pdb&project_id=${projectId}`, {
                method: 'POST'
            });
            
            // Poll for results (Simulated delay, in production use the logging system)
            setTimeout(async () => {
                const res = await fetch(`http://localhost:8000/api/p2rank/results?pdb_filename=${pdbSearch}.pdb&project_id=${projectId}`);
                const data = await res.json();
                setPockets(data);
                setLoading(false);
            }, 5000);
        } catch (e) { 
            alert("P2Rank failed to start");
            setLoading(false); 
        }
    };
    const loadPdbIntoViewer = async (url) => {
        if (!pluginRef.current) return;
        const plugin = pluginRef.current;
        await plugin.clear();
        const data = await plugin.builders.data.download({ url: `http://localhost:8000${url}` });
        const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
        await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default');
    };

    const highlightPocket = async (residueIdsString) => {
        if (!pluginRef.current || !residueIdsString) return;
        const plugin = pluginRef.current;

        // Parse "A_107 A_108" into Molstar selection
        const segments = residueIdsString.split(' ');
        const selection = segments.map(s => {
            const [chain, resIdx] = s.split('_');
            return { auth_asym_id: chain, auth_seq_id: parseInt(resIdx) };
        });

        // Use Molstar Script to highlight and zoom
        const data = plugin.managers.structure.hierarchy.current.structures[0].cell.obj.data;
        const sel = Script.getStructureSelection(Q => Q.struct.generator.atomGroups({
            'residue-test': Q.core.set.has([Q.set(...selection.map(s => s.auth_seq_id)), Q.ammp('auth_seq_id')]),
            'chain-test': Q.core.set.has([Q.set(...selection.map(s => s.auth_asym_id)), Q.ammp('auth_asym_id')])
        }), data);

        plugin.managers.structure.selection.fromSelection(sel);
        plugin.managers.camera.focusSelection(sel);
    };
    // 2. Synchronize Mol* View -> React Table
    // If user clicks the "Model Selector" arrows in Mol*, update the highlighted row in our table
    useEffect(() => {
        if (!pluginRef.current) return;
        
        const sub = pluginRef.current.state.data.events.cell.stateUpdated.subscribe(e => {
            const newIndex = e.cell.transform.params?.modelIndex;
            if (typeof newIndex === 'number' && newIndex !== activePose) {
                setActivePose(newIndex);
            }
        });

        return () => sub.unsubscribe();
    }, [isPluginReady, activePose]);

    // 3. Load Receptor and Ligand into Mol*
    const loadStructure = async (receptorUrl, ligandUrl, poseIndex) => {
        if (!pluginRef.current) return;
        const plugin = pluginRef.current;
        
        // Clear previous viewer state
        await plugin.clear();

        const baseUrl = 'http://localhost:8000';
        
        try {
            // LOAD RECEPTOR (Protein)
            const recData = await plugin.builders.data.download({ url: `${baseUrl}${receptorUrl}` });
            const recTraj = await plugin.builders.structure.parseTrajectory(recData, 'pdb');
            const recModel = await plugin.builders.structure.createModel(recTraj);
            const recStruct = await plugin.builders.structure.createStructure(recModel);
            
            await plugin.builders.structure.representation.addRepresentation(recStruct, {
                type: 'cartoon', 
                color: 'uniform', 
                colorParams: { value: 0x3b82f6 } // Blue color
            });

            // LOAD LIGAND (Specific Pose)
            const ligData = await plugin.builders.data.download({ url: `${baseUrl}${ligandUrl}` });
            const ligTraj = await plugin.builders.structure.parseTrajectory(ligData, 'pdb');
            
            // This isolates the specific MODEL index from the PDB file
            const ligModel = await plugin.builders.structure.createModel(ligTraj, { modelIndex: poseIndex });
            const ligStruct = await plugin.builders.structure.createStructure(ligModel);

            await plugin.builders.structure.representation.addRepresentation(ligStruct, {
                type: 'ball-and-stick', 
                color: 'element-symbol'
            });

            // Re-center camera on the new docking site
            plugin.managers.camera.reset();

        } catch (e) {
            console.error("Mol* Visualization Error:", e);
        }
    };

    // 4. Handle Backend Docking Run
    const handleRunDocking = async () => {
        if (!projectId) return alert("Please select a project first.");
        
        setLoading(true);
        try {
            const response = await fetch('http://localhost:8000/api/docking/run', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    project_id: String(projectId),
                    smiles: smiles.trim()
                })
            });

            const data = await response.json();
            
            if (data.success) {
                setLastResultUrls(data);
                
                // Map real scores from backend
                const realPoses = data.scores.map((s, index) => ({
                    id: index,
                    energy: s.toFixed(2)
                }));
                
                setPoses(realPoses);
                setActivePose(0);
                
                // Initial 3D load
                loadStructure(data.receptor_url, data.ligand_url, 0);
            } else {
                alert("Docking Failed: " + data.error);
            }
        } catch (e) {
            console.error(e);
            alert("Network Error");
        }
        setLoading(false);
    };

    // 5. Switch Pose via Table Click
    const selectPose = (idx) => {
        if (loading) return;
        setActivePose(idx);
        if (lastResultUrls) {
            loadStructure(lastResultUrls.receptor_url, lastResultUrls.ligand_url, idx);
        }
    };

    return (
        <div className="main-content">
            <div className="docking-container-wrapper">
                
                {/* Header Row */}
                <div className="docking-header">
                    <button className="primary-btn" onClick={() => navigate(`/dashboard/${projectId}`)}>
                        ← Back to Dashboard
                    </button>
                    <h2 style={{ margin: 0, fontSize: '1.4rem' }}>Molecular Docking: {projectId}</h2>
                    <div style={{ width: '130px' }} /> {/* Symmetry Spacer */}
                </div>

                <div className="docking-grid">
                    <div className="card" style={{ gridColumn: 'span 2' }}>
                        <h3>PROTEIN & POCKET PREPARATION</h3>
                        <div style={{ display: 'flex', gap: '10px', marginBottom: '20px' }}>
                            <input 
                                type="text" 
                                placeholder="Search RCSB PDB (e.g. 6OD6)" 
                                className="technical-input"
                                onChange={(e) => setPdbSearch(e.target.value)}
                            />
                            <button className="primary-btn" onClick={() => handleFetchPdb()}>Import from RCSB</button>
                            <button className="primary-btn" onClick={() => handleRunP2Rank()}>Run P2Rank Pocket Detection</button>
                        </div>

                        {/* Results Table */}
                        {pockets.length > 0 && (
                            <div className="technical-table-container" style={{maxHeight: '300px', overflow: 'auto'}}>
                                <table>
                                    <thead>
                                        <tr>
                                            <th>Rank</th>
                                            <th>Score</th>
                                            <th>Probability</th>
                                            <th>Center (X, Y, Z)</th>
                                            <th>Actions</th>
                                        </tr>
                                    </thead>
                                    <tbody>
                                        {pockets.map(p => (
                                            <tr key={p.rank}>
                                                <td>{p.rank}</td>
                                                <td>{p.score}</td>
                                                <td>{p.probability}</td>
                                                <td style={{fontSize: '0.7rem'}}>{p.center_x}, {p.center_y}, {p.center_z}</td>
                                                <td>
                                                    <button 
                                                        className="primary-btn" 
                                                        onClick={() => highlightPocket(p.residue_ids)}
                                                    >
                                                        View in Mol*
                                                    </button>
                                                </td>
                                            </tr>
                                        ))}
                                    </tbody>
                                </table>
                            </div>
                        )}
                    </div>
                    
                    {/* Left: Input & Poses */}
                    <div className="card" style={{ margin: 0, padding: '24px', display: 'flex', flexDirection: 'column' }}>
                        <h3 style={{ fontSize: '0.7rem', color: 'var(--text-dim)', letterSpacing: '1px', marginBottom: '15px' }}>CONTROLS</h3>
                        
                        <div style={{ background: '#111', padding: '12px', borderRadius: '8px', marginBottom: '20px', border: '1px solid #222' }}>
                            <div style={{ color: '#4ade80', fontSize: '0.8rem' }}>• 3D Engine: Ready</div>
                            <div style={{ color: 'var(--text-dim)', fontSize: '0.8rem' }}>• Project ID: {projectId}</div>
                        </div>

                        <label style={{ fontSize: '0.7rem', color: 'var(--text-dim)', marginBottom: '5px' }}>LIGAND (SMILES)</label>
                        <textarea 
                            className="technical-input" 
                            style={{ 
                                background: '#000', border: '1px solid var(--border)', color: 'white', 
                                padding: '12px', borderRadius: '6px', resize: 'none', height: '80px', 
                                marginBottom: '20px', width: '100%', boxSizing: 'border-box',
                                fontFamily: 'monospace', fontSize: '0.85rem'
                            }}
                            value={smiles}
                            onChange={(e) => setSmiles(e.target.value)}
                        />

                        <button 
                            className="primary-btn" 
                            style={{ width: '100%', padding: '12px' }} 
                            onClick={handleRunDocking} 
                            disabled={loading || !isPluginReady}
                        >
                            {loading ? "Running Vina..." : "Start Docking Calculation"}
                        </button>

                        {poses.length > 0 && (
                            <div style={{ marginTop: '30px' }}>
                                <h3 style={{ fontSize: '0.7rem', color: 'var(--text-dim)', marginBottom: '10px' }}>DOCKING POSES</h3>
                                <div className="pose-table-container">
                                    {poses.map((p) => (
                                        <div 
                                            key={p.id} 
                                            className={`pose-row ${activePose === p.id ? 'active' : ''}`}
                                            onClick={() => selectPose(p.id)}
                                        >
                                            <span className="pose-label">#{p.id + 1}</span>
                                            <span style={{color: '#eee'}}>Pose Result</span>
                                            <span className="pose-energy">{p.energy}</span>
                                        </div>
                                    ))}
                                </div>
                            </div>
                        )}
                    </div>

                    {/* Right: Mol* Viewer */}
                    <div className="molstar-container" ref={parentRef} />

                </div>
            </div>
        </div>
    );
};

export default Docking;