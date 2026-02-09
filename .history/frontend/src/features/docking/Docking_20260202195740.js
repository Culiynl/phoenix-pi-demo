import React, { useEffect, useRef, useState } from 'react';
import { useParams } from 'react-router-dom'; // To get Project ID
import { createPluginUI } from 'molstar/lib/mol-plugin-ui';
import { DefaultPluginUISpec } from 'molstar/lib/mol-plugin-ui/spec';
import { renderReact18 } from 'molstar/lib/mol-plugin-ui/react18';
import 'molstar/lib/mol-plugin-ui/skin/dark.scss';

const Docking = () => {
    const { projectId } = useParams(); // Get ID from URL
    const parentRef = useRef(null);
    const pluginRef = useRef(null);
    
    // UI State
    const [smiles, setSmiles] = useState("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"); // Default Ibuprofen
    const [loading, setLoading] = useState(false);
    const [result, setResult] = useState(null);
    const [isPluginReady, setIsPluginReady] = useState(false);

    // Initialize Mol*
    useEffect(() => {
        if (!parentRef.current || pluginRef.current) return;

        async function init() {
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
            } catch (error) {
                console.error("Failed to initialize Mol*", error);
            }
        }
        init();
    }, []);

    // Function to call Backend
    const handleRunDocking = async () => {
        if (!projectId) {
            alert("No project selected. Please go to Projects and select one.");
            return;
        }
        
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

            const data = await response.json();
            
            if (data.success) {
                setResult(data);
                loadStructure(data.file_url);
            } else {
                alert("Docking Failed: " + data.detail);
            }

        } catch (e) {
            console.error(e);
            alert("Network Error");
        }
        setLoading(false);
    };

    // Helper to load URL into Mol*
    const loadStructure = async (url) => {
        if (!pluginRef.current) return;
        const plugin = pluginRef.current;
        
        plugin.clear();
        
        // Use full URL
        const fullUrl = `http://localhost:8000${url}`;
        
        try {
            const data = await plugin.builders.data.download({ url: fullUrl }, { state: { isGhost: true } });
            const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
            await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default');
        } catch (e) {
            console.error("Mol* Load Error", e);
        }
    };

    return (
        <div className="docking-container" style={{display: 'flex', gap: '20px', height: '100%'}}>
            
            {/* Left Panel: Controls */}
            <div className="card" style={{flex: '0 0 300px', height: 'fit-content'}}>
                <h2>Docking Controls</h2>
                
                <div style={{marginBottom: '10px'}}>
                    <label style={{fontSize: '0.8rem', color: '#888'}}>Project ID</label>
                    <div style={{fontWeight: 'bold'}}>{projectId || "None Selected"}</div>
                </div>

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
                    disabled={loading || !isPluginReady || !projectId}
                    style={{width: '100%'}}
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