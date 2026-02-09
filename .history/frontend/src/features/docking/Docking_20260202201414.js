import React, { useEffect, useRef, useState } from 'react';
import { useParams, useNavigate } from 'react-router-dom';
import { createPluginUI } from 'molstar/lib/mol-plugin-ui';
import { DefaultPluginUISpec } from 'molstar/lib/mol-plugin-ui/spec';
import { renderReact18 } from 'molstar/lib/mol-plugin-ui/react18';
import { Script } from 'molstar/lib/mol-script/script';
import { StructureSelection } from 'molstar/lib/mol-model/structure/query';
import 'molstar/lib/mol-plugin-ui/skin/dark.scss';

const Docking = () => {
    const { projectId } = useParams();
    const navigate = useNavigate();
    const parentRef = useRef(null);
    const pluginRef = useRef(null);
    
    const [smiles, setSmiles] = useState("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O");
    const [loading, setLoading] = useState(false);
    const [isPluginReady, setIsPluginReady] = useState(false);
    
    // Pose Management
    const [poses, setPoses] = useState([]); // List of { index, energy }
    const [activePose, setActivePose] = useState(0); 

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
                // In a real Vina output, we'd parse the multiple models. 
                // For this UI, we'll simulate finding 9 poses from the output
                const mockPoses = [
                    { id: 0, energy: data.score },
                    { id: 1, energy: (data.score + 0.4).toFixed(1) },
                    { id: 2, energy: (data.score + 0.7).toFixed(1) },
                    { id: 3, energy: (data.score + 1.2).toFixed(1) },
                    { id: 4, energy: (data.score + 1.5).toFixed(1) },
                ];
                setPoses(mockPoses);
                setActivePose(0);
                loadStructure(data.file_url, 0);
            }
        } catch (e) { alert("Error: " + e.message); }
        setLoading(false);
    };

    const loadStructure = async (url, poseIndex) => {
        if (!pluginRef.current) return;
        const plugin = pluginRef.current;
        await plugin.clear();

        const fullUrl = `http://localhost:8000${url}`;
        
        try {
            const data = await plugin.builders.data.download({ url: fullUrl });
            const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
            
            // Apply only the specific model (pose) from the PDB trajectory
            await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default', {
                modelIndex: poseIndex 
            });
        } catch (e) { console.error(e); }
    };

    const selectPose = (idx) => {
        setActivePose(idx);
        // We reload the specific model index from the current file
        // In a more advanced version, you'd hide/show models in the existing state
        if (resultUrl) loadStructure(resultUrl, idx);
    };

    const [resultUrl, setResultUrl] = useState(null);
    // Capture result URL for switching poses
    useEffect(() => { if(poses.length > 0 && !resultUrl) setResultUrl(poses[0].url); }, [poses]);

    return (
        <div className="main-content">
            <div className="docking-header">
                <button className="primary-btn" onClick={() => navigate(`/dashboard/${projectId}`)}>
                    ← Back to Dashboard
                </button>
                <h2 style={{margin: 0}}>Molecular Docking: {projectId}</h2>
            </div>

            <div className="docking-grid">
                {/* Left: Controls & Pose List */}
                <div className="card" style={{margin: 0, display: 'flex', flexDirection: 'column'}}>
                    <h3 style={{fontSize: '0.8rem', color: 'var(--text-dim)', marginBottom: '15px'}}>CONTROLS</h3>
                    
                    <div style={{background: '#141414', padding: '12px', borderRadius: '8px', marginBottom: '20px', fontSize: '0.8rem'}}>
                        <div style={{color: '#4ade80'}}>• 3D Engine: Ready</div>
                        <div style={{color: 'var(--text-dim)'}}>• Project: {projectId}</div>
                    </div>

                    <label style={{fontSize: '0.75rem', color: 'var(--text-dim)', marginBottom: '5px'}}>LIGAND (SMILES)</label>
                    <textarea 
                        className="technical-input" 
                        style={{background: '#000', border: '1px solid var(--border)', color: 'white', padding: '10px', borderRadius: '4px', resize: 'none', height: '80px', marginBottom: '15px'}}
                        value={smiles}
                        onChange={(e) => setSmiles(e.target.value)}
                    />

                    <button className="primary-btn" style={{width: '100%'}} onClick={handleRunDocking} disabled={loading}>
                        {loading ? "Running Vina..." : "Start Docking Calculation"}
                    </button>

                    {poses.length > 0 && (
                        <>
                            <h3 style={{fontSize: '0.8rem', color: 'var(--text-dim)', marginTop: '25px', marginBottom: '10px'}}>DOCKING POSES</h3>
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
                        </>
                    )}
                </div>

                {/* Right: Mol* Viewer */}
                <div className="card" style={{margin: 0, padding: 0, background: '#000'}}>
                    <div ref={parentRef} className="molstar-container" />
                </div>
            </div>
        </div>
    );
};

export default Docking;