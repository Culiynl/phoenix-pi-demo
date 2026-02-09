import React, { useEffect, useRef, useState } from 'react';
import { createPluginUI } from 'molstar/lib/mol-plugin-ui';
import { DefaultPluginUISpec } from 'molstar/lib/mol-plugin-ui/spec';
import { renderReact18 } from 'molstar/lib/mol-plugin-ui/react18';
import 'molstar/lib/mol-plugin-ui/skin/dark.scss';

const Docking = () => {
    const parentRef = useRef(null);
    const pluginRef = useRef(null); // Ref to store plugin instance immediately
    const [fileName, setFileName] = useState("No file selected");
    const [isPluginReady, setIsPluginReady] = useState(false); // State to toggle UI availability

    useEffect(() => {
        // Guard: If the plugin is already initialized, do not run again.
        // This fixes the "You are calling ReactDOMClient.createRoot()..." error.
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
            spec.components = {
                remoteState: 'none'
            };

            try {
                // Initialize Mol*
                const ctx = await createPluginUI({
                    target: parentRef.current,
                    spec: spec,
                    render: renderReact18 
                });

                // Assign to ref immediately to block subsequent init calls
                pluginRef.current = ctx;
                setIsPluginReady(true);
            } catch (error) {
                console.error("Failed to initialize Mol*", error);
            }
        }

        init();

        // Cleanup: Only dispose if the component is truly unmounting (not just re-rendering)
        return () => {
            // Optional: You can add cleanup logic here, but for Mol* in React 18,
            // it is often safer to let the instance persist during dev hot-reloads 
            // to avoid WebGL context loss, or strictly manage disposal.
            // For now, we leave the instance alive to prevent the crash.
        };
    }, []);

    const handleFileUpload = async (e) => {
        const file = e.target.files[0];
        // Use the ref here to ensure we access the active instance
        if (!file || !pluginRef.current) return;
        
        setFileName(file.name);

        const plugin = pluginRef.current;
        plugin.clear();

        const reader = new FileReader();
        reader.onload = async (event) => {
            const data = event.target.result;
            
            try {
                // CORRECT API WORKFLOW:
                // 1. Create a raw data node from the text content
                const dataNode = await plugin.builders.data.rawData({ 
                    data: data, 
                    label: file.name 
                });

                // 2. Parse that data node into a Trajectory
                // Note: We use parseTrajectory, not parse
                const trajectory = await plugin.builders.structure.parseTrajectory(dataNode, 'pdb');

                // 3. Create the visual hierarchy (Cartoon/Ball-and-Stick)
                await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default');
                
            } catch (err) {
                console.error("Failed to load PDB:", err);
                alert("Error loading PDB file. See console for details.");
            }
        };
        reader.readAsText(file);
    };

    return (
        <div className="card">
            <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '20px' }}>
                <h2>Molecular Docking</h2>
                <div style={{ display: 'flex', gap: '10px', alignItems: 'center' }}>
                    <span style={{ fontSize: '0.8rem', color: 'gray' }}>{fileName}</span>
                    <label className="primary-btn" style={{ 
                        cursor: isPluginReady ? 'pointer' : 'not-allowed', 
                        padding: '8px 16px', 
                        background: isPluginReady ? '#333' : '#999', 
                        color: '#fff', 
                        borderRadius: '4px' 
                    }}>
                        Upload PDB
                        <input 
                            type="file" 
                            accept=".pdb" 
                            onChange={handleFileUpload} 
                            style={{ display: 'none' }} 
                            disabled={!isPluginReady}
                        />
                    </label>
                </div>
            </div>

            <div 
                ref={parentRef} 
                style={{ 
                    width: '100%', 
                    height: '600px', 
                    border: '1px solid #ccc', 
                    borderRadius: '8px', 
                    overflow: 'hidden',
                    background: '#000',
                    position: 'relative'
                }} 
            />
        </div>
    );
};

export default Docking;