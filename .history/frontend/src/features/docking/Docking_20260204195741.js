import React, { useEffect, useRef, useState, useCallback } from 'react';
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
    const isInitialMount = useRef(true);
    
    // UI State
    const [pdbSearch, setPdbSearch] = useState("");
    const [searchResults, setSearchResults] = useState([]);
    const [isSearching, setIsSearching] = useState(false);
    const [loading, setLoading] = useState(false);
    
    // P2Rank & Docking State
    const [pockets, setPockets] = useState([]);
    const [dockingCenter, setDockingCenter] = useState(null);
    const [p2rankLogs, setP2rankLogs] = useState([]);
    const p2LogEndRef = useRef(null);

    // Vina State
    const [smiles, setSmiles] = useState("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O");
    const [poses, setPoses] = useState([]);
    const [activePose, setActivePose] = useState(0);
    const [lastResultUrls, setLastResultUrls] = useState(null);

    // --- Mol* Initialization (React 18 Safe) ---
    useEffect(() => {
        if (!parentRef.current || pluginRef.current) return;

        async function init() {
            const spec = DefaultPluginUISpec();
            spec.layout = { initialShowControls: false, initialShowSequence: false, initialShowLog: false };
            spec.components = { remoteState: 'none', sequence: 'none', log: 'none' };

            try {
                const ctx = await createPluginUI({ target: parentRef.current, spec, render: renderReact18 });
                pluginRef.current = ctx;
            } catch (error) { console.error("Mol* Init Error", error); }
        }
        init();

        return () => {
            if (pluginRef.current) {
                pluginRef.current.dispose();
                pluginRef.current = null;
            }
        };
    }, []);

    // --- Search Logic (Manual Trigger Only) ---
    const handlePdbSearch = async () => {
        if (!pdbSearch || isSearching) return;
        setIsSearching(true);
        try {
            const res = await fetch(`http://localhost:8000/api/pdb/search?q=${encodeURIComponent(pdbSearch)}`);
            const data = await res.json();
            setSearchResults(data || []);
        } catch (e) {
            console.error("Search failed", e);
        } finally {
            setIsSearching(false);
        }
    };

    // --- Import & Preview Logic ---
    const handleImport = async (pdbId) => {
        setLoading(true);
        try {
            const res = await fetch(`http://localhost:8000/api/pdb/fetch/${pdbId.toUpperCase()}`);
            if (!res.ok) throw new Error("Structure download failed");
            const data = await res.json();
            
            setPdbSearch(data.filename); 
            // Important: This loads the locally saved file
            loadPdbIntoViewer(`/static/proteins/${data.filename}`);
            alert(`Imported ${data.filename}`);
        } catch (e) { alert(e.message); }
        setLoading(false);
    };

    const handlePreview = (pdbId) => {
        const url = `https://files.rcsb.org/download/${pdbId.toUpperCase()}.cif`;
        loadPdbFromUrl(url);
    };

    // --- Visualizer Helpers ---
    const loadPdbFromUrl = async (url) => {
        if (!pluginRef.current) return;
        try {
            await pluginRef.current.clear();
            const data = await pluginRef.current.builders.data.download({ url }, { state: { isGhost: true } });
            const format = url.toLowerCase().endsWith('.cif') ? 'mmcif' : 'pdb';
            const trajectory = await pluginRef.current.builders.structure.parseTrajectory(data, format);
            await pluginRef.current.builders.structure.hierarchy.applyPreset(trajectory, 'default');
        } catch (e) { console.error("Preview failed", e); }
    };

    const loadPdbIntoViewer = async (relativeUrl) => {
        if (!pluginRef.current) return;
        const format = relativeUrl.toLowerCase().endsWith('.cif') ? 'mmcif' : 'pdb';
        try {
            await pluginRef.current.clear();
            const data = await pluginRef.current.builders.data.download({ url: `http://localhost:8000${relativeUrl}` });
            const trajectory = await pluginRef.current.builders.structure.parseTrajectory(data, format);
            await pluginRef.current.builders.structure.hierarchy.applyPreset(trajectory, 'default');
        } catch (e) { console.error("Viewer failed", e); }
    };

    // --- P2Rank Execution & Polling ---
    const handleRunP2Rank = async () => {
        if (!pdbSearch) return alert("Please import a protein first.");
        setLoading(true);
        setPockets([]);
        setP2rankLogs([{ msg: "Starting P2Rank Analysis...", type: "info" }]);

        try {
            await fetch(`http://localhost:8000/api/p2rank/run?pdb_filename=${pdbSearch}&project_id=${projectId}`, { method: 'POST' });

            const pollInterval = setInterval(async () => {
                try {
                    const logRes = await fetch(`http://localhost:8000/api/p2rank/status/${projectId}`);
                    const logs = await logRes.json();
                    setP2rankLogs(logs);

                    const res = await fetch(`http://localhost:8000/api/p2rank/results?pdb_filename=${pdbSearch}&project_id=${projectId}`);
                    if (res.status === 200) {
                        const data = await res.json();
                        setPockets(data);
                        if (data.length > 0) {
                            setDockingCenter({ x: data[0].center_x, y: data[0].center_y, z: data[0].center_z });
                            // Don't auto-focus immediately to prevent camera jitter during polling
                        }
                        clearInterval(pollInterval);
                        setLoading(false);
                    }
                } catch (err) { /* silent poll */ }
            }, 3000);
        } catch (e) { setLoading(false); }
    };

    // --- JSX Layout ---
    return (
        <div className="main-content">
            <div className="docking-container-wrapper">
                
                <div className="docking-header">
                    <button className="primary-btn" onClick={() => navigate(`/dashboard/${projectId}`)}>← Back</button>
                    <h2 style={{ margin: 0 }}>Docking & Structure Prep</h2>
                </div>

                {/* SEARCH SECTION */}
                <div className="card">
                    <h3 style={{ fontSize: '0.7rem', color: 'var(--text-dim)', marginBottom: '15px' }}>RCSB PDB SEARCH</h3>
                    <div className="action-group">
                        <input 
                            type="text" 
                            className="technical-input" 
                            placeholder="Search (e.g. Alzheimer's, 6OD6)" 
                            value={pdbSearch} 
                            onChange={(e) => setPdbSearch(e.target.value)}
                            onKeyDown={(e) => e.key === 'Enter' && handlePdbSearch()}
                        />
                        <button className="primary-btn" style={{background: 'var(--primary)', color: 'white'}} onClick={handlePdbSearch} disabled={isSearching}>
                            {isSearching ? 'Searching...' : 'Search Database'}
                        </button>
                    </div>

                    {searchResults.length > 0 && (
                        <div className="technical-table-container" style={{ marginTop: '15px', maxHeight: '250px' }}>
                            <table className="technical-table">
                                <thead><tr><th>ID</th><th>Description</th><th>Actions</th></tr></thead>
                                <tbody>
                                    {searchResults.map(r => (
                                        <tr key={r.id}>
                                            <td style={{ fontWeight: 'bold', color: 'var(--primary)' }}>{r.id}</td>
                                            <td style={{ fontSize: '0.8rem' }}>{r.title}</td>
                                            <td>
                                                <div style={{ display: 'flex', gap: '5px' }}>
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
                </div>

                <div className="docking-grid">
                    <div style={{ display: 'flex', flexDirection: 'column', gap: '20px' }}>
                        {/* P2RANK CONTROLS */}
                        <div className="card" style={{height: 'auto'}}>
                            <h3 style={{ fontSize: '0.7rem', color: 'var(--text-dim)' }}>POCKET DETECTION</h3>
                            <button className="primary-btn" style={{width: '100%', marginTop: '10px'}} onClick={handleRunP2Rank} disabled={loading}>
                                Run P2Rank Analysis
                            </button>
                            <div className="feudal-log" style={{ height: '120px', marginTop: '10px', fontSize: '0.7rem' }}>
                                {p2rankLogs.map((l, i) => <div key={i}>> {l.msg}</div>)}
                                <div ref={p2LogEndRef} />
                            </div>
                        </div>

                        {/* VINA CONTROLS */}
                        <div className="card">
                            <h3 style={{ fontSize: '0.7rem', color: 'var(--text-dim)' }}>DOCKING EXECUTION</h3>
                            <textarea className="technical-input" style={{height: '80px', marginTop: '10px'}} value={smiles} onChange={(e) => setSmiles(e.target.value)} />
                            <button className="primary-btn" style={{width: '100%', marginTop: '10px', background: dockingCenter ? 'var(--primary)' : '#333'}} disabled={!dockingCenter}>
                                Start Docking
                            </button>
                        </div>
                    </div>

                    <div className="molstar-container" ref={parentRef} />
                </div>
            </div>
        </div>
    );
};

export default Docking;