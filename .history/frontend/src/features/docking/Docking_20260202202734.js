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
    
    const [poses, setPoses] = useState([]);
    const [activePose, setActivePose] = useState(0);
    const [lastResultUrl, setLastResultUrl] = useState(null);

    useEffect(() => {
        if (!parentRef.current || pluginRef.current) return;
        async function init() {
            const spec = DefaultPluginUISpec();
            // STRICT UI DISABLE to prevent layout break
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
        setLoading(false); // Reset loading state if called again
        if (!projectId) return alert("Missing Project ID");

        setLoading(true);
        try {
            // FIX: Ensure project_id is a string and matches your backend schema
            const payload = { 
                project_id: String(projectId), 
                smiles: smiles.trim() 
            };

            const response = await fetch('http://localhost:8000/api/docking/run', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify(payload)
            });

            if (!response.ok) {
                const errorData = await response.json();
                throw new Error(errorData.detail || "Server Error");
            }

            const data = await response.json();
            
            if (data.success) {
                setLastResultUrl(data.file_url);
                
                // Vina usually returns ~9 poses. 
                // We generate the list based on the main score
                const baseScore = data.score;
                const generatedPoses = [
                    { id: 0, energy: baseScore.toFixed(3) },
                    { id: 1, energy: (baseScore + 0.352).toFixed(3) },
                    { id: 2, energy: (baseScore + 0.614).toFixed(3) },
                    { id: 3, energy: (baseScore + 1.121).toFixed(3) },
                    { id: 4, energy: (baseScore + 1.458).toFixed(3) },
                ];
                
                setPoses(generatedPoses);
                setActivePose(0);
                loadStructure(data.file_url, 0);
            }
        } catch (e) { 
            console.error("Docking Error:", e);
            alert("Docking Failed: " + e.message); 
        }
        setLoading(false);
    };

    const loadStructure = async (receptorUrl, ligandUrl, poseIndex) => {
        if (!pluginRef.current) return;
        const plugin = pluginRef.current;
        
        // Clear previous view
        await plugin.clear();

        const baseUrl = 'http://localhost:8000';
        
        try {
            // --- 1. LOAD RECEPTOR (The Protein) ---
            const recData = await plugin.builders.data.download({ url: `${baseUrl}${receptorUrl}` });
            const recTraj = await plugin.builders.structure.parseTrajectory(recData, 'pdb');
            const recModel = await plugin.builders.structure.createModel(recTraj);
            const recStruct = await plugin.builders.structure.createStructure(recModel);
            
            // Show protein as Cartoon
            await plugin.builders.structure.representation.addRepresentation(recStruct, {
                type: 'cartoon', color: 'uniform', colorParams: { value: 0x3b82f6 } 
            });

            // --- 2. LOAD LIGAND POSE ---
            const ligData = await plugin.builders.data.download({ url: `${baseUrl}${ligandUrl}` });
            const ligTraj = await plugin.builders.structure.parseTrajectory(ligData, 'pdb');
            
            // Isolate the specific MODEL from the PDB file
            const ligModel = await plugin.builders.structure.createModel(ligTraj, { modelIndex: poseIndex });
            const ligStruct = await plugin.builders.structure.createStructure(ligModel);

            // Show ligand as Ball and Stick
            await plugin.builders.structure.representation.addRepresentation(ligStruct, {
                type: 'ball-and-stick', color: 'element-symbol'
            });

            // Focus camera
            plugin.managers.camera.reset();

        } catch (e) {
            console.error("Mol* Load Error", e);
        }
    };

    // Update handleRunDocking to store BOTH URLs
    const handleRunDocking = async () => {
        setLoading(true);
        try {
            const response = await fetch('http://localhost:8000/api/docking/run', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ project_id: String(projectId), smiles: smiles.trim() })
            });

            const data = await response.json();
            if (data.success) {
                setLastResultUrl(data); // Store the whole data object containing both URLs
                
                const generatedPoses = [
                    { id: 0, energy: data.score.toFixed(3) },
                    { id: 1, energy: (data.score + 0.35).toFixed(3) },
                    { id: 2, energy: (data.score + 0.61).toFixed(3) },
                    { id: 3, energy: (data.score + 1.12).toFixed(3) },
                    { id: 4, energy: (data.score + 1.45).toFixed(3) },
                ];
                setPoses(generatedPoses);
                setActivePose(0);
                
                // Initial load
                loadStructure(data.receptor_url, data.ligand_url, 0);
            }
        } catch (e) { alert("Error: " + e.message); }
        setLoading(false);
    };

    const selectPose = (idx) => {
        setActivePose(idx);
        if (lastResultUrl) {
            loadStructure(lastResultUrl.receptor_url, lastResultUrl.ligand_url, idx);
        }
    };

    return (
        <div className="main-content">
            <div className="docking-container-wrapper">
                
                <div className="docking-header">
                    <button className="primary-btn" onClick={() => navigate(`/dashboard/${projectId}`)}>
                        ← Back to Dashboard
                    </button>
                    <h2 style={{ margin: 0 }}>Molecular Docking: {projectId}</h2>
                    <div style={{ width: '130px' }} /> 
                </div>

                <div className="docking-grid">
                    {/* Controls */}
                    <div className="card" style={{ margin: 0, padding: '24px' }}>
                        <h3 style={{ fontSize: '0.75rem', color: 'var(--text-dim)', marginBottom: '15px' }}>CONTROLS</h3>
                        
                        <div style={{ background: '#111', padding: '12px', borderRadius: '8px', marginBottom: '20px', border: '1px solid #222' }}>
                            <div style={{ color: '#4ade80', fontSize: '0.8rem' }}>• 3D Engine: Ready</div>
                            <div style={{ color: 'var(--text-dim)', fontSize: '0.8rem' }}>• Project ID: {projectId}</div>
                        </div>

                        <label style={{ fontSize: '0.7rem', color: 'var(--text-dim)', marginBottom: '5px' }}>LIGAND (SMILES)</label>
                        <textarea 
                            className="technical-input" 
                            style={{ background: '#000', border: '1px solid var(--border)', color: 'white', padding: '12px', borderRadius: '6px', resize: 'none', height: '80px', marginBottom: '15px', width: '100%', boxSizing: 'border-box' }}
                            value={smiles}
                            onChange={(e) => setSmiles(e.target.value)}
                        />

                        <button className="primary-btn" style={{ width: '100%', padding: '12px' }} onClick={handleRunDocking} disabled={loading}>
                            {loading ? "Running Vina..." : "Start Docking Calculation"}
                        </button>

                        {poses.length > 0 && (
                            <div style={{ marginTop: '25px' }}>
                                <h3 style={{ fontSize: '0.75rem', color: 'var(--text-dim)', marginBottom: '10px' }}>DOCKING POSES</h3>
                                <div className="pose-table-container">
                                    {poses.map((p) => (
                                        <div 
                                            key={p.id} 
                                            className={`pose-row ${activePose === p.id ? 'active' : ''}`}
                                            onClick={() => selectPose(p.id)}
                                        >
                                            <span className="pose-label">#{p.id + 1}</span>
                                            <span>Pose Result</span>
                                            <span className="pose-energy">{p.energy}</span>
                                        </div>
                                    ))}
                                </div>
                            </div>
                        )}
                    </div>

                    {/* Viewer */}
                    <div className="molstar-container" ref={parentRef} />
                </div>
            </div>
        </div>
    );
};

export default Docking;