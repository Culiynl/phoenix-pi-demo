import React, { useEffect, useRef, useState } from 'react';
import { useParams, useNavigate } from 'react-router-dom';
import { createPluginUI } from 'molstar/lib/mol-plugin-ui';
import { DefaultPluginUISpec } from 'molstar/lib/mol-plugin-ui/spec';
import { renderReact18 } from 'molstar/lib/mol-plugin-ui/react18';
import 'molstar/lib/mol-plugin-ui/skin/dark.scss';

const Docking = () => {
    const { projectId } = useParams();
    const navigate = useNavigate();
    const navigate = useNavigate();
    const parentRef = useRef(null);
    const pluginRef = useRef(null);
    
    // UI State
    const [smiles, setSmiles] = useState("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O");
    const [loading, setLoading] = useState(false);
    const [result, setResult] = useState(null);
    const [isPluginReady, setIsPluginReady] = useState(false);

    // Initialize Mol*
    useEffect(() => {
        if (!parentRef.current || pluginRef.current) return;

        async function init() {
            console.log("[Mol*] Initializing...");
            const spec = DefaultPluginUISpec();
            spec.layout = {
                initialShowControls: false,
                initialShowRemoteState: false,
                initialShowSequence: false,
                initialShowLog: false,
                controlsDisplay: 'reactive', 
            };
            spec.components = { remoteState: 'none' };

            try {
                const ctx = await createPluginUI({
                    target: parentRef.current,
                    spec: spec,
                    render: renderReact18 
                });
                pluginRef.current = ctx;
                setIsPluginReady(true);
                console.log("[Mol*] Ready.");
            } catch (error) {
                console.error("[Mol*] Initialization Failed:", error);
            }
        }
        init();
    }, []);

    // Function to call Backend
    const handleRunDocking = async () => {
        console.log("--> BUTTON CLICKED <--");

        // 1. Check Project ID
        if (!projectId) {
            alert("⚠️ Missing Project ID.\n\nPlease go to the 'Projects' page and click on a specific project to open this dashboard.");
            return;
        }

        // 2. Check Mol*
        if (!isPluginReady) {
            alert("⚠️ 3D Viewer is not ready yet.\nPlease wait a moment for Mol* to load.");
            return;
        }

        console.log("--> Prerequisites met. Sending request...");
        
        setLoading(true);
        try {
            const response = await fetch('http://localhost:8000/api/docking/run', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    project_id: projectId,
                    smiles: smiles
                })
            });

            console.log("--> Response Status:", response.status);
            const data = await response.json();
            console.log("--> Data:", data);
            
            if (data.success) {
                setResult(data);
                loadStructure(data.file_url);
            } else {
                alert("Docking Failed: " + (data.detail || data.error));
            }

        } catch (e) {
            console.error(e);
            alert("Network Error. Check console (F12) for details.");
        }
        setLoading(false);
    };

    const loadStructure = async (url) => {
        if (!pluginRef.current) return;
        const fullUrl = `http://localhost:8000${url}`;
        try {
            const plugin = pluginRef.current;
            plugin.clear();
            const data = await plugin.builders.data.download({ url: fullUrl }, { state: { isGhost: true } });
            const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
            await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default');
        } catch (e) {
            console.error("Mol* Load Error", e);
        }
    };

    // --- DEBUG STATUS ---
    // If you came from the sidebar, projectId is likely undefined.
    const isProjectIdMissing = !projectId;

    return (
        <div className="docking-container" style={{display: 'flex', gap: '20px', height: '100%'}}>
            
            {/* Left Panel: Controls */}
            <div className="card" style={{flex: '0 0 300px', height: 'fit-content'}}>
                <h2>Docking Controls</h2>
                
                {/* STATUS INDICATORS - VISUAL DEBUGGING */}
                <div style={{background: '#222', padding: '10px', borderRadius: '5px', marginBottom: '15px', fontSize: '0.85rem'}}>
                    <div style={{color: projectId ? '#4ade80' : '#f87171'}}>
                        • Project ID: {projectId || "MISSING (Select a project first)"}
                    </div>
                    <div style={{color: isPluginReady ? '#4ade80' : '#fbbf24'}}>
                        • 3D Engine: {isPluginReady ? "Ready" : "Loading..."}
                    </div>
                    <div style={{color: 'white'}}>
                        • Backend: http://localhost:8000
                    </div>
                </div>

                {isProjectIdMissing && (
                    <button 
                        className="primary-btn" 
                        style={{background: '#f59e0b', marginBottom: '15px', border: 'none', color: 'black'}}
                        onClick={() => navigate('/projects')}
                    >
                        ← Go Select a Project
                    </button>
                )}

                <div style={{marginBottom: '15px'}}>
                    <label style={{fontSize: '0.8rem', color: '#888'}}>Ligand (SMILES)</label>
                    <textarea 
                        className="technical-input" 
                        rows={4}
                        value={smiles}
                        onChange={(e) => setSmiles(e.target.value)}
                        style={{width: '100%', marginTop: '5px'}}
                    />
                </div>

                <button 
                    className="primary-btn" 
                    onClick={handleRunDocking}
                    // REMOVED DISABLED PROP FOR DEBUGGING
                    // disabled={loading || !isPluginReady || !projectId}
                    style={{width: '100%', opacity: loading ? 0.7 : 1, cursor: 'pointer'}}
                >
                    {loading ? "Calculating..." : "Run Vina Docking"}
                </button>

                {result && (
                    <div style={{marginTop: '20px', padding: '10px', background: '#1a1a1a', borderRadius: '4px'}}>
                        <h4 style={{margin: '0 0 5px 0', color: '#4ade80'}}>Success</h4>
                        <div>Score: <strong>{result.score} kcal/mol</strong></div>
                    </div>
                )}
            </div>

            {/* Right Panel: Viewer */}
            <div className="card" style={{flex: 1, padding: 0, overflow: 'hidden', display: 'flex', flexDirection: 'column'}}>
                 <div 
                    ref={parentRef} 
                    style={{ 
                        flex: 1,
                        background: '#000',
                        position: 'relative'
                    }} 
                />
            </div>
        </div>
    );
};

export default Docking;