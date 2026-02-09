import React, { useEffect, useRef, useState } from 'react';
import { useParams, useNavigate } from 'react-router-dom';
import { createPluginUI } from 'molstar/lib/mol-plugin-ui';
import { DefaultPluginUISpec } from 'molstar/lib/mol-plugin-ui/spec';
import { renderReact18 } from 'molstar/lib/mol-plugin-ui/react18';
import { Script } from 'molstar/lib/mol-script/script';
import 'molstar/lib/mol-plugin-ui/skin/dark.scss';

const Docking = () => {
    const { projectId } = useParams();
    const navigate = useNavigate();
    const parentRef = useRef(null);
    const pluginRef = useRef(null);
    
    // UI & Search State
    const [pdbSearch, setPdbSearch] = useState("");
    const [searchResults, setSearchResults] = useState([]);
    const [isSearching, setIsSearching] = useState(false);
    const [loading, setLoading] = useState(false);
    
    // Pocket & Docking Center State
    const [pockets, setPockets] = useState([]);
    const [dockingCenter, setDockingCenter] = useState(null); // {x, y, z} from P2Rank
    
    // Vina State
    const [smiles, setSmiles] = useState("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O");
    const [poses, setPoses] = useState([]);
    const [activePose, setActivePose] = useState(0);
    const [lastResultUrls, setLastResultUrls] = useState(null);

    // --- Mol* Initialization ---
    useEffect(() => {
        if (!parentRef.current || pluginRef.current) return;
        async function init() {
            const spec = DefaultPluginUISpec();
            spec.layout = { initialShowControls: false, initialShowSequence: false, initialShowLog: false };
            spec.components = { remoteState: 'none' };
            try {
                const ctx = await createPluginUI({ target: parentRef.current, spec, render: renderReact18 });
                pluginRef.current = ctx;
            } catch (error) { console.error("Mol* Init Error:", error); }
        }
        init();
        return () => { if (pluginRef.current) { pluginRef.current.dispose(); pluginRef.current = null; } };
    }, []);

    // --- PDB Search & Import ---
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

    const handleImport = async (pdbId) => {
        setLoading(true);
        try {
            const res = await fetch(`http://localhost:8000/api/pdb/fetch/${pdbId.toUpperCase()}`);
            if (res.status === 404) throw new Error("PDB not found on RCSB");
            const data = await res.json();
            alert(`Structure ${pdbId} added. Run P2Rank now.`);
            setPdbSearch(pdbId.toUpperCase());
            loadPdbIntoViewer(`/static/proteins/${pdbId.toUpperCase()}.pdb`);
        } catch (e) { alert(e.message); }
        setLoading(false);
    };

    const handlePreview = (pdbId) => {
        const url = `https://files.rcsb.org/download/${pdbId.toUpperCase()}.pdb`;
        loadPdbFromUrl(url);
    };

    const loadPdbFromUrl = async (url) => {
        if (!pluginRef.current) return;
        const plugin = pluginRef.current;
        await plugin.clear();
        try {
            const data = await plugin.builders.data.download({ url }, { state: { isGhost: true } });
            const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
            await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default');
        } catch (e) { alert("Visualizer error: File might not be available for preview."); }
    };

    const loadPdbIntoViewer = async (url) => {
        if (!pluginRef.current) return;
        await pluginRef.current.clear();
        const data = await pluginRef.current.builders.data.download({ url: `http://localhost:8000${url}` });
        const trajectory = await pluginRef.current.builders.structure.parseTrajectory(data, 'pdb');
        await pluginRef.current.builders.structure.hierarchy.applyPreset(trajectory, 'default');
    };

    // --- P2Rank & Center Logic ---
    const handleRunP2Rank = async () => {
        setLoading(true);
        const filename = `${pdbSearch.toUpperCase()}.pdb`;
        try {
            await fetch(`http://localhost:8000/api/p2rank/run?pdb_filename=${filename}&project_id=${projectId}`, { method: 'POST' });
            
            // Poll for results
            const poll = setInterval(async () => {
                const res = await fetch(`http://localhost:8000/api/p2rank/results?pdb_filename=${filename}&project_id=${projectId}`);
                if (res.ok) {
                    const data = await res.json();
                    setPockets(data);
                    // AUTO-SET CENTER to the #1 Rank pocket
                    if (data.length > 0) {
                        const top = data[0];
                        setDockingCenter({ x: top.center_x, y: top.center_y, z: top.center_z });
                    }
                    clearInterval(poll);
                    setLoading(false);
                }
            }, 3000);
        } catch (e) { setLoading(false); }
    };

    const highlightPocket = async (residueIdsString) => {
        if (!pluginRef.current || !residueIdsString) return;
        const segments = residueIdsString.split(' ').map(s => {
            const [chain, resIdx] = s.split('_');
            return { auth_asym_id: chain, auth_seq_id: parseInt(resIdx) };
        });
        const data = pluginRef.current.managers.structure.hierarchy.current.structures[0].cell.obj.data;
        const sel = Script.getStructureSelection(Q => Q.struct.generator.atomGroups({
            'residue-test': Q.core.set.has([Q.set(...segments.map(s => s.auth_seq_id)), Q.ammp('auth_seq_id')]),
            'chain-test': Q.core.set.has([Q.set(...segments.map(s => s.auth_asym_id)), Q.ammp('auth_asym_id')])
        }), data);
        pluginRef.current.managers.structure.selection.fromSelection(sel);
        pluginRef.current.managers.camera.focusSelection(sel);
    };

    // --- Vina Docking ---
    const handleRunDocking = async () => {
        setLoading(true);
        try {
            const response = await fetch('http://localhost:8000/api/docking/run', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ 
                    project_id: String(projectId), 
                    smiles: smiles.trim(),
                    center_override: dockingCenter // Pass the P2Rank center to the backend
                })
            });
            const data = await response.json();
            if (data.success) {
                setLastResultUrls(data);
                setPoses(data.scores.map((s, i) => ({ id: i, energy: s.toFixed(2) })));
                setActivePose(0);
                loadStructure(data.receptor_url, data.ligand_url, 0);
            }
        } catch (e) { alert("Docking Error"); }
        setLoading(false);
    };

    const loadStructure = async (receptorUrl, ligandUrl, poseIndex) => {
        const plugin = pluginRef.current;
        await plugin.clear();
        const baseUrl = 'http://localhost:8000';
        const recData = await plugin.builders.data.download({ url: `${baseUrl}${receptorUrl}` });
        const recTraj = await plugin.builders.structure.parseTrajectory(recData, 'pdb');
        const recStruct = await plugin.builders.structure.createStructure(await plugin.builders.structure.createModel(recTraj));
        await plugin.builders.structure.representation.addRepresentation(recStruct, { type: 'cartoon', color: 'uniform', colorParams: { value: 0x3b82f6 } });

        const ligData = await plugin.builders.data.download({ url: `${baseUrl}${ligandUrl}` });
        const ligTraj = await plugin.builders.structure.parseTrajectory(ligData, 'pdb');
        const ligStruct = await plugin.builders.structure.createStructure(await plugin.builders.structure.createModel(ligTraj, { modelIndex: poseIndex }));
        await plugin.builders.structure.representation.addRepresentation(ligStruct, { type: 'ball-and-stick', color: 'element-symbol' });
        plugin.managers.camera.reset();
    };

    return (
        <div className="main-content">
            <div className="docking-container-wrapper" style={{ transform: 'translateX(-60px)' }}>
                <div className="docking-header">
                    <button className="primary-btn" onClick={() => navigate(`/dashboard/${projectId}`)}>← Back</button>
                    <h2 style={{ margin: 0 }}>Project Docking Workspace</h2>
                </div>

                <div className="docking-grid">
                    {/* RCSB Search */}
                    <div className="card" style={{ gridColumn: 'span 2' }}>
                        <h3 style={{ fontSize: '0.7rem', color: 'var(--text-dim)' }}>RCSB PDB SEARCH</h3>
                        <div className="action-group">
                            <input type="text" placeholder="Search (e.g. Alzheimer's, 6OD6)" className="technical-input" value={pdbSearch} onChange={(e) => setPdbSearch(e.target.value)} />
                            <button className="primary-btn" onClick={handlePdbSearch}>Search</button>
                            <button className="primary-btn" onClick={handleRunP2Rank} disabled={!pdbSearch || loading}>Run P2Rank</button>
                        </div>
                        
                        {searchResults.length > 0 && (
                            <div className="technical-table-container" style={{marginTop: '15px', maxHeight: '250px'}}>
                                <table className="technical-table">
                                    <thead><tr><th>ID</th><th>Description</th><th>Actions</th></tr></thead>
                                    <tbody>
                                        {searchResults.map(r => (
                                            <tr key={r.id}>
                                                <td>{r.id}</td>
                                                <td style={{fontSize: '0.75rem'}}>{r.title}</td>
                                                <td>
                                                    <button className="primary-btn" onClick={() => handlePreview(r.id)}>Preview</button>
                                                    <button className="primary-btn" style={{background: '#4ade80', color: '#000'}} onClick={() => handleImport(r.id)}>Import</button>
                                                </td>
                                            </tr>
                                        ))}
                                    </tbody>
                                </table>
                            </div>
                        )}
                    </div>

                    {/* P2Rank Result Table */}
                    {pockets.length > 0 && (
                        <div className="card" style={{ gridColumn: 'span 2' }}>
                            <h3 style={{ fontSize: '0.7rem', color: 'var(--text-dim)' }}>DETECTED POCKETS (AUTO-CENTER READY)</h3>
                            <table className="technical-table">
                                <thead><tr><th>Rank</th><th>Score</th><th>Center (X,Y,Z)</th><th>Actions</th></tr></thead>
                                <tbody>
                                    {pockets.map(p => (
                                        <tr key={p.rank} style={{background: p.rank === 1 ? '#1a1a1a' : 'transparent'}}>
                                            <td>#{p.rank} {p.rank === 1 && "(Target)"}</td>
                                            <td>{p.score}</td>
                                            <td style={{fontSize: '0.7rem'}}>{p.center_x}, {p.center_y}, {p.center_z}</td>
                                            <td><button className="primary-btn" onClick={() => highlightPocket(p.residue_ids)}>Focus</button></td>
                                        </tr>
                                    ))}
                                </tbody>
                            </table>
                        </div>
                    )}

                    {/* Vina Controls */}
                    <div className="card">
                        <h3 style={{ fontSize: '0.7rem', color: 'var(--text-dim)' }}>VINA DOCKING</h3>
                        {dockingCenter && <div style={{fontSize: '0.7rem', color: '#4ade80', marginBottom: '10px'}}>✓ Using P2Rank Center</div>}
                        <textarea className="technical-input" value={smiles} onChange={(e) => setSmiles(e.target.value)} style={{height: '80px', marginBottom: '15px'}} />
                        <button className="primary-btn" style={{width: '100%'}} onClick={handleRunDocking} disabled={loading}>Run Docking</button>
                        
                        {poses.map(p => (
                            <div key={p.id} className={`pose-row ${activePose === p.id ? 'active' : ''}`} onClick={() => { setActivePose(p.id); loadStructure(lastResultUrls.receptor_url, lastResultUrls.ligand_url, p.id); }}>
                                <span>Pose #{p.id+1}</span><span className="pose-energy">{p.energy}</span>
                            </div>
                        ))}
                    </div>

                    <div className="molstar-container" ref={parentRef} />
                </div>
            </div>
        </div>
    );
};

export default Docking;