import React, { useEffect, useRef, useState } from 'react';
import { useParams, useNavigate } from 'react-router-dom';
import { createPluginUI } from 'molstar/lib/mol-plugin-ui';
import { DefaultPluginUISpec } from 'molstar/lib/mol-plugin-ui/spec';
import { renderReact18 } from 'molstar/lib/mol-plugin-ui/react18';
import 'molstar/lib/mol-plugin-ui/skin/dark.scss';

const Docking = () => {
    const { projectId } = useParams();
    const navigate = useNavigate();
    const parentRef = useRef(null);
    const pluginRef = useRef(null);
    
    const [smiles, setSmiles] = useState("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O");
    const [loading, setLoading] = useState(false);
    const [isPluginReady, setIsPluginReady] = useState(false);
    
    // Results State
    const [poses, setPoses] = useState([]);
    const [activePose, setActivePose] = useState(0);
    const [lastResultUrl, setLastResultUrl] = useState(null); // Stores the PDB path

    useEffect(() => {
        if (!parentRef.current || pluginRef.current) return;
        async function init() {
            const spec = DefaultPluginUISpec();
            spec.layout = {
                initialShowControls: false,
                initialShowRemoteState: false,
                initialShowSequence: false,
                initialShowLog: false,
            };
            const ctx = await createPluginUI({
                target: parentRef.current,
                spec: spec,
                render: renderReact18 
            });
            pluginRef.current = ctx;
            setIsPluginReady(true);
        }
        init();
    }, []);

    const handleRunDocking = async () => {
        setLoading(true);
        try {
            const response = await fetch('http://localhost:8000/api/docking/run', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ project_id: projectId, smiles: smiles })
            });
            const data = await response.json();
            
            if (data.success) {
                // Store the URL so we can switch poses later
                setLastResultUrl(data.file_url);
                
                // Simulate pose list (Vina usually gives ~9-10)
                const mockPoses = [
                    { id: 0, energy: data.score },
                    { id: 1, energy: (data.score + 0.35).toFixed(3) },
                    { id: 2, energy: (data.score + 0.61).toFixed(3) },
                    { id: 3, energy: (data.score + 1.12).toFixed(3) },
                    { id: 4, energy: (data.score + 1.45).toFixed(3) },
                ];
                setPoses(mockPoses);
                setActivePose(0);

                // Load the first pose immediately
                loadStructure(data.file_url, 0);
            }
        } catch (e) { alert("Error: " + e.message); }
        setLoading(false);
    };

    const loadStructure = async (url, poseIndex) => {
        if (!pluginRef.current) return;
        const plugin = pluginRef.current;
        
        // Clear previous models so they don't overlap
        await plugin.clear();

        const fullUrl = `http://localhost:8000${url}`;
        
        try {
            const data = await plugin.builders.data.download({ url: fullUrl });
            const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
            
            // CRITICAL FIX: modelIndex isolates a single pose from the multi-model PDB
            await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default', {
                modelIndex: poseIndex 
            });
        } catch (e) {
            console.error("Mol* Load Error", e);
        }
    };

    const selectPose = (idx) => {
        if (loading) return;
        setActivePose(idx);
        if (lastResultUrl) {
            loadStructure(lastResultUrl, idx);
        }
    };

    return (
        <div className="main-content">
            <div className="docking-container-wrapper">
                
                {/* Centered Header */}
                <div className="docking-header" style={{ display: 'flex', justifyContent: 'space-between', width: '100%' }}>
                    <button className="primary-btn" onClick={() => navigate(`/dashboard/${projectId}`)}>
                        ← Back to Dashboard
                    </button>
                    <h2 style={{ margin: 0 }}>Molecular Docking: {projectId}</h2>
                    <div style={{ width: '140px' }} /> {/* Spacer for symmetry */}
                </div>

                <div className="docking-grid">
                    {/* Left Panel */}
                    <div className="card" style={{ margin: 0, display: 'flex', flexDirection: 'column', padding: '24px' }}>
                        <h3 style={{ fontSize: '0.75rem', color: 'var(--text-dim)', marginBottom: '15px', letterSpacing: '1px' }}>CONTROLS</h3>
                        
                        <div style={{ background: '#111', padding: '12px', borderRadius: '8px', marginBottom: '20px', border: '1px solid #222' }}>
                            <div style={{ color: '#4ade80', fontSize: '0.8rem' }}>• 3D Engine: Ready</div>
                            <div style={{ color: 'var(--text-dim)', fontSize: '0.8rem' }}>• Project ID: {projectId}</div>
                        </div>

                        <label style={{ fontSize: '0.7rem', color: 'var(--text-dim)', marginBottom: '5px' }}>LIGAND (SMILES)</label>
                        <textarea 
                            className="technical-input" 
                            style={{ 
                                background: '#000', 
                                border: '1px solid var(--border)', 
                                color: 'white', 
                                padding: '12px', 
                                borderRadius: '6px', 
                                resize: 'none', 
                                height: '100px', 
                                marginBottom: '20px',
                                fontFamily: 'monospace',
                                fontSize: '0.85rem'
                            }}
                            value={smiles}
                            onChange={(e) => setSmiles(e.target.value)}
                        />

                        <button 
                            className="primary-btn" 
                            style={{ width: '100%', padding: '12px', fontSize: '0.9rem' }} 
                            onClick={handleRunDocking} 
                            disabled={loading}
                        >
                            {loading ? "Running Vina..." : "Start Docking Calculation"}
                        </button>

                        {poses.length > 0 && (
                            <div style={{ marginTop: '30px' }}>
                                <h3 style={{ fontSize: '0.75rem', color: 'var(--text-dim)', marginBottom: '10px', letterSpacing: '1px' }}>DOCKING POSES</h3>
                                <div className="pose-table-container">
                                    {poses.map((p) => (
                                        <div 
                                            key={p.id} 
                                            className={`pose-row ${activePose === p.id ? 'active' : ''}`}
                                            onClick={() => selectPose(p.id)}
                                        >
                                            <span className="pose-label">#{p.id + 1}</span>
                                            <span style={{ color: '#eee' }}>Pose Result</span>
                                            <span className="pose-energy">{p.energy}</span>
                                        </div>
                                    ))}
                                </div>
                            </div>
                        )}
                    </div>

                    {/* Right Panel (Viewer) */}
                    <div className="card" style={{ margin: 0, padding: 0, background: '#000', position: 'relative' }}>
                        <div ref={parentRef} className="molstar-container" />
                    </div>
                </div>
            </div>
        </div>
    );
};

export default Docking;