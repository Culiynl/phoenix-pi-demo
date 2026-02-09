import React, { useEffect, useRef, useState } from 'react';
import { useParams, useNavigate } from 'react-router-dom';
import { createPluginUI } from 'molstar/lib/mol-plugin-ui';
import { DefaultPluginUISpec } from 'molstar/lib/mol-plugin-ui/spec';
import { renderReact18 } from 'molstar/lib/mol-plugin-ui/react18';
import 'molstar/lib/mol-plugin-ui/skin/dark.scss';
import { Script } from 'molstar/lib/mol-script/script';

const Docking = () => {
    const { projectId } = useParams();
    const navigate = useNavigate();
    const parentRef = useRef(null);
    const pluginRef = useRef(null);
    
    // --- State: Search & Import ---
    const [pdbSearch, setPdbSearch] = useState("");
    const [searchResults, setSearchResults] = useState([]);
    const [isSearching, setIsSearching] = useState(false);
    const [loading, setLoading] = useState(false);
    
    // --- State: P2Rank & Pockets ---
    const [pockets, setPockets] = useState([]);
    const [dockingCenter, setDockingCenter] = useState(null); // Calculated from P2Rank
    
    // --- State: Vina Docking ---
    const [smiles, setSmiles] = useState("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O");
    const [poses, setPoses] = useState([]);
    const [activePose, setActivePose] = useState(0);
    const [lastResultUrls, setLastResultUrls] = useState(null);
    const [isPluginReady, setIsPluginReady] = useState(false);
    const [p2rankLogs, setP2rankLogs] = useState([]);
    const p2LogEndRef = useRef(null);
    
    useEffect(() => {
        p2LogEndRef.current?.scrollIntoView({ behavior: "smooth" });
    }, [p2rankLogs]);

    const handleRunP2Rank = async () => {
        if (!pdbSearch) return alert("Select or Import a protein first.");
        setLoading(true);
        setPockets([]);
        setP2rankLogs([{ msg: "Initializing P2Rank process...", type: "info" }]);
        
        const filename = pdbSearch; 

        try {
            await fetch(`http://localhost:8000/api/p2rank/run?pdb_filename=${filename}&project_id=${projectId}`, { 
                method: 'POST' 
            });

            // POLLING LOOP
            const pollInterval = setInterval(async () => {
                try {
                    // 1. Fetch Logs
                    const logRes = await fetch(`http://localhost:8000/api/p2rank/status/${projectId}`);
                    const logs = await logRes.json();
                    setP2rankLogs(logs);

                    // 2. Try fetching results
                    const res = await fetch(`http://localhost:8000/api/p2rank/results?pdb_filename=${filename}&project_id=${projectId}`);
                    if (res.status === 200) {
                        const data = await res.json();
                        setPockets(data);
                        if (data.length > 0) {
                            setDockingCenter({ x: data[0].center_x, y: data[0].center_y, z: data[0].center_z });
                            highlightPocket(data[0].residue_ids);
                        }
                        clearInterval(pollInterval);
                        setLoading(false);
                    }
                    
                    // Stop if we see "complete" or "failed" in logs
                    const lastLog = logs[logs.length - 1]?.msg.toLowerCase();
                    if (lastLog?.includes("complete") || lastLog?.includes("failed")) {
                        clearInterval(pollInterval);
                        setLoading(false);
                    }
                } catch (err) { console.log("Polling logs..."); }
            }, 2000);

        } catch (e) { 
            setLoading(false); 
            alert("Server error starting P2Rank.");
        }
    };

    const getFormat = (url) => {
        if (url.toLowerCase().endsWith('.cif') || url.toLowerCase().endsWith('.bcif')) return 'mmcif';
        return 'pdb';
    };
    // 1. Initialize Molstar
    useEffect(() => {
        if (!parentRef.current || pluginRef.current) return;
        async function init() {
            const spec = DefaultPluginUISpec();
            spec.layout = {
                initialShowControls: false,
                initialShowRemoteState: false,
                initialShowSequence: false,
                initialShowLog: false,
                initialShowLeftPanel: false,
            };
            spec.components = { remoteState: 'none', sequence: 'none', log: 'none' };

            try {
                const ctx = await createPluginUI({
                    target: parentRef.current,
                    spec: spec,
                    render: renderReact18 
                });
                pluginRef.current = ctx;
                setIsPluginReady(true);
            } catch (error) {
                console.error("Molstar Init Error:", error);
            }
        }
        init();
        return () => { if (pluginRef.current) { pluginRef.current.dispose(); pluginRef.current = null; } };
    }, []);

    // 2. RCSB Search & Import Logic
    const handlePdbSearch = async () => {
        if (!pdbSearch) return;
        setIsSearching(true);
        try {
            const res = await fetch(`http://localhost:8000/api/pdb/search?q=${encodeURIComponent(pdbSearch)}`);
            const data = await res.json();
            setSearchResults(data || []);
        } catch (e) { alert("Search failed."); }
        setIsSearching(false);
    };

    const handlePreview = (pdbId) => {
        const url = `https://files.rcsb.org/download/${pdbId.toUpperCase()}.pdb`;
        loadPdbFromUrl(url);
    };

    const handleImport = async (pdbId) => {
        setLoading(true);
        try {
            const res = await fetch(`http://localhost:8000/api/pdb/fetch/${pdbId.toUpperCase()}`);
            if (!res.ok) throw new Error("Structure not found");
            
            const data = await res.json();
            // data.filename is "9BA8.cif"
            setPdbSearch(data.filename); 
            alert(`Imported ${data.filename}`);
            loadPdbIntoViewer(`/static/proteins/${data.filename}`);
        } catch (e) {
            alert(e.message);
        }
        setLoading(false);
    };

    // 3. Visualization Helpers
    const loadPdbFromUrl = async (url) => {
        if (!pluginRef.current) return;
        const plugin = pluginRef.current;
        const format = getFormat(url); // Dynamically detect format
        
        try {
            await plugin.clear();
            const data = await plugin.builders.data.download({ url }, { state: { isGhost: true } });
            const trajectory = await plugin.builders.structure.parseTrajectory(data, format);
            await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default');
        } catch (e) {
            console.error("Preview failed", e);
            alert("RCSB Preview Error: Structure might not be released yet or format is unsupported.");
        }
    };

    const loadPdbIntoViewer = async (relativeUrl) => {
        if (!pluginRef.current) return;
        const plugin = pluginRef.current;
        const format = getFormat(relativeUrl);
        
        try {
            await plugin.clear();
            const data = await plugin.builders.data.download({ url: `http://localhost:8000${relativeUrl}` });
            const trajectory = await plugin.builders.structure.parseTrajectory(data, format);
            await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default');
        } catch (e) {
            console.error("Viewer failed", e);
        }
    };

    // 4. P2Rank Logic
    const handleRunP2Rank = async () => {
        if (!pdbSearch) return alert("Please import a protein first.");
        setLoading(true);
        
        // Ensure we use the exact filename (9BA8.cif)
        const filename = pdbSearch; 

        try {
            await fetch(`http://localhost:8000/api/p2rank/run?pdb_filename=${filename}&project_id=${projectId}`, { 
                method: 'POST' 
            });

            const poll = setInterval(async () => {
                const res = await fetch(`http://localhost:8000/api/p2rank/results?pdb_filename=${filename}&project_id=${projectId}`);
                if (res.status === 200) {
                    const data = await res.json();
                    setPockets(data);
                    if (data.length > 0) {
                        setDockingCenter({ x: data[0].center_x, y: data[0].center_y, z: data[0].center_z });
                    }
                    clearInterval(poll);
                    setLoading(false);
                }
            }, 3000);
        } catch (e) { setLoading(false); }
    };

    const highlightPocket = async (residueIdsString) => {
        if (!pluginRef.current || !residueIdsString) return;
        const plugin = pluginRef.current;
        const segments = residueIdsString.split(' ').map(s => {
            const [chain, resIdx] = s.split('_');
            return { auth_asym_id: chain, auth_seq_id: parseInt(resIdx) };
        });

        const data = plugin.managers.structure.hierarchy.current.structures[0].cell.obj.data;
        const sel = Script.getStructureSelection(Q => Q.struct.generator.atomGroups({
            'residue-test': Q.core.set.has([Q.set(...segments.map(s => s.auth_seq_id)), Q.ammp('auth_seq_id')]),
            'chain-test': Q.core.set.has([Q.set(...segments.map(s => s.auth_asym_id)), Q.ammp('auth_asym_id')])
        }), data);

        plugin.managers.structure.selection.fromSelection(sel);
        plugin.managers.camera.focusSelection(sel);
    };

    // 5. Vina Docking Logic
    const handleRunDocking = async () => {
        setLoading(true);
        try {
            const response = await fetch('http://localhost:8000/api/docking/run', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ 
                    project_id: String(projectId), 
                    smiles: smiles.trim(),
                    center_override: dockingCenter 
                })
            });
            const data = await response.json();
            if (data.success) {
                setLastResultUrls(data);
                setPoses(data.scores.map((s, i) => ({ id: i, energy: s.toFixed(2) })));
                setActivePose(0);
                loadComplex(data.receptor_url, data.ligand_url, 0);
            }
        } catch (e) { alert("Docking Error"); }
        setLoading(false);
    };

    const loadComplex = async (recUrl, ligUrl, idx) => {
        const plugin = pluginRef.current;
        await plugin.clear();
        const baseUrl = 'http://localhost:8000';
        
        const recData = await plugin.builders.data.download({ url: `${baseUrl}${recUrl}` });
        const recStruct = await plugin.builders.structure.createStructure(await plugin.builders.structure.createModel(await plugin.builders.structure.parseTrajectory(recData, 'pdb')));
        await plugin.builders.structure.representation.addRepresentation(recStruct, { type: 'cartoon', color: 'uniform', colorParams: { value: 0x3b82f6 } });

        const ligData = await plugin.builders.data.download({ url: `${baseUrl}${ligUrl}` });
        const ligStruct = await plugin.builders.structure.createStructure(await plugin.builders.structure.createModel(await plugin.builders.structure.parseTrajectory(ligData, 'pdb'), { modelIndex: idx }));
        await plugin.builders.structure.representation.addRepresentation(ligStruct, { type: 'ball-and-stick', color: 'element-symbol' });
        
        plugin.managers.camera.reset();
    };

    return (
        <div className="main-content">
            <div className="docking-container-wrapper">
                
                {/* Header */}
                <div className="docking-header" style={{ justifyContent: 'space-between' }}>
                    <div style={{ display: 'flex', alignItems: 'center', gap: '20px' }}>
                        <button className="primary-btn" onClick={() => navigate(`/dashboard/${projectId}`)}>← Back</button>
                        <h2 style={{ margin: 0 }}>Structure Preparation & Pocket Analysis</h2>
                    </div>
                    {dockingCenter && (
                        <div style={{ background: '#111', padding: '8px 15px', borderRadius: '8px', border: '1px solid #222', color: '#4ade80', fontSize: '0.8rem' }}>
                            🎯 Active Docking Center: {dockingCenter.x.toFixed(1)}, {dockingCenter.y.toFixed(1)}, {dockingCenter.z.toFixed(1)}
                        </div>
                    )}
                </div>

                {/* Protein Search & P2Rank Results (Stacked vertically at the top) */}
                <div className="card">
                    <h3 style={{ fontSize: '0.7rem', color: 'var(--text-dim)', marginBottom: '15px' }}>RCSB PDB DATABASE & POCKET DETECTION</h3>
                    <div className="action-group" style={{ marginBottom: '20px' }}>
                        <input 
                            type="text" 
                            placeholder="Search keyword or ID (e.g. Alzheimer's, 6OD6)" 
                            className="technical-input" 
                            value={pdbSearch} 
                            onChange={(e) => setPdbSearch(e.target.value)} 
                            onKeyDown={(e) => e.key === 'Enter' && handlePdbSearch()} 
                        />
                        <button className="primary-btn" style={{ background: 'var(--primary)', color: 'white' }} onClick={handlePdbSearch}>Search PDB</button>
                        <button className="primary-btn" onClick={handleRunP2Rank} disabled={loading || !pdbSearch}>
                            {loading ? "Analyzing..." : "Run P2Rank Analysis"}
                        </button>
                    </div>

                    {/* RCSB Search Results Table */}
                    {searchResults.length > 0 && (
                        <div className="technical-table-container" style={{ maxHeight: '300px', marginBottom: '20px' }}>
                            <table className="technical-table">
                                <thead><tr><th>ID</th><th>Description</th><th>Resolution</th><th>Actions</th></tr></thead>
                                <tbody>
                                    {searchResults.map(r => (
                                        <tr key={r.id}>
                                            <td style={{ fontWeight: 'bold', color: 'var(--primary)' }}>{r.id}</td>
                                            <td style={{ fontSize: '0.8rem' }}>{r.title}</td>
                                            <td>{r.resolution ? `${r.resolution} Å` : 'N/A'}</td>
                                            <td>
                                                <div style={{ display: 'flex', gap: '10px' }}>
                                                    <button className="primary-btn" onClick={() => handlePreview(r.id)}>Preview</button>
                                                    <button className="primary-btn" style={{ background: '#4ade80', color: '#000' }} onClick={() => handleImport(r.id)}>Import</button>
                                                </div>
                                            </td>
                                        </tr>
                                    ))}
                                </tbody>
                            </table>
                        </div>
                    )}

                    {/* P2Rank Pocket Results Table */}
                    {pockets.length > 0 && (
                        <div className="technical-table-container" style={{ maxHeight: '400px' }}>
                            <table className="technical-table">
                                <thead><tr><th>Rank</th><th>Score</th><th>Prob.</th><th>Center (X,Y,Z)</th><th>Actions</th></tr></thead>
                                <tbody>
                                    {pockets.map(p => (
                                        <tr key={p.rank} style={{ background: dockingCenter?.x === p.center_x ? 'rgba(59, 130, 246, 0.15)' : 'transparent' }}>
                                            <td>#{p.rank}</td>
                                            <td style={{ color: '#4ade80' }}>{p.score}</td>
                                            <td>{(p.probability * 100).toFixed(1)}%</td>
                                            <td style={{ fontSize: '0.7rem', fontFamily: 'monospace' }}>{p.center_x.toFixed(2)}, {p.center_y.toFixed(2)}, {p.center_z.toFixed(2)}</td>
                                            <td>
                                                <div style={{ display: 'flex', gap: '8px' }}>
                                                    <button className="primary-btn" onClick={() => highlightPocket(p.residue_ids)}>Focus</button>
                                                    <button className="primary-btn" style={{ background: '#fff', color: '#000' }} 
                                                        onClick={() => setDockingCenter({ x: p.center_x, y: p.center_y, z: p.center_z })}>
                                                        Target
                                                    </button>
                                                </div>
                                            </td>
                                        </tr>
                                    ))}
                                </tbody>
                            </table>
                        </div>
                    )}
                </div>
                <div className="card" style={{ gridColumn: 'span 2' }}>
                    <div style={{ display: 'flex', justifyContent: 'space-between' }}>
                        <h3 style={{ fontSize: '0.7rem', color: 'var(--text-dim)' }}>P2RANK LIVE CONSOLE</h3>
                        {loading && <span className="loader-small">Processing Structure...</span>}
                    </div>
                    <div className="feudal-log" style={{ height: '150px', marginTop: '10px', background: '#000', fontSize: '0.75rem' }}>
                        {p2rankLogs.map((log, i) => (
                            <div key={i} style={{ color: log.type === 'success' ? '#4ade80' : log.type === 'error' ? '#f87171' : '#aaa' }}>
                                > {log.msg}
                            </div>
                        ))}
                        <div ref={p2LogEndRef} />
                    </div>
                </div>

                {/* Main Docking Grid (Controls + Mol*) */}
                <div className="docking-grid">
                    
                    {/* Vina Controls Card */}
                    <div className="card" style={{ height: '100%', padding: '24px' }}>
                        <h3 style={{ fontSize: '0.7rem', color: 'var(--text-dim)', marginBottom: '15px' }}>DOCKING EXECUTION</h3>
                        
                        <label style={{ fontSize: '0.7rem', color: 'var(--text-dim)' }}>LIGAND (SMILES)</label>
                        <textarea 
                            className="technical-input" 
                            style={{ background: '#000', border: '1px solid var(--border)', color: 'white', padding: '12px', borderRadius: '6px', height: '100px', marginBottom: '20px', resize: 'none', fontFamily: 'monospace' }}
                            value={smiles} 
                            onChange={(e) => setSmiles(e.target.value)} 
                        />

                        <button 
                            className="primary-btn" 
                            style={{ width: '100%', padding: '15px', background: dockingCenter ? 'var(--primary)' : '#333', color: 'white' }} 
                            onClick={handleRunDocking} 
                            disabled={loading || !dockingCenter}
                        >
                            {loading ? "Running Simulation..." : "Start Auto-Centered Docking"}
                        </button>
                        {!dockingCenter && <p style={{ fontSize: '0.65rem', color: '#f87171', marginTop: '10px' }}>⚠️ Target a pocket above to enable docking.</p>}

                        {/* Docking Poses List */}
                        {poses.length > 0 && (
                            <div style={{ marginTop: '30px', flexGrow: 1, overflow: 'hidden', display: 'flex', flexDirection: 'column' }}>
                            <h3 style={{ fontSize: '0.7rem', color: 'var(--text-dim)', marginBottom: '10px' }}>DOCKING POSES</h3>
                            <div className="pose-table-container">
                                {poses.map((p) => (
                                    <div key={p.id} className={`pose-row ${activePose === p.id ? 'active' : ''}`} 
                                        onClick={() => { setActivePose(p.id); loadComplex(lastResultUrls.receptor_url, lastResultUrls.ligand_url, p.id); }}>
                                        <span className="pose-label">#{p.id + 1}</span>
                                        <span>Binding Energy</span>
                                        <span className="pose-energy">{p.energy}</span>
                                    </div>
                                ))}
                            </div>
                        </div>
                    )}
                </div>

                {/* Mol* Viewer Card */}
                <div className="molstar-container" ref={parentRef} />

            </div>
        </div>
    </div>
);
};

export default Docking;