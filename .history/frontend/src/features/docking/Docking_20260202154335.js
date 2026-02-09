import React, { useEffect, useRef, useState } from 'react';
// Remove createRoot; Mol* handles the rendering internally via createPluginUI
import { PluginContext } from 'molstar/lib/mol-plugin/context';
import { DefaultPluginUISpec } from 'molstar/lib/mol-plugin-ui/spec';
import { createPluginUI } from 'molstar/lib/mol-plugin-ui'; // Import the UI creator helper
import 'molstar/lib/mol-plugin-ui/skin/dark.scss';
import { renderReact18 } from 'molstar/lib/mol-plugin-ui/react18'; // Optional: ensures React 18 support if your Mol* version supports it

const Docking = () => {
    const parentRef = useRef(null);
    const [pluginCtx, setPluginCtx] = useState(null);
    const [fileName, setFileName] = useState("No file selected");
    const pluginRef = useRef(null); // Keep track of the instance to prevent double-init

    useEffect(() => {
        if (!parentRef.current || pluginRef.current) return;

        async function init() {
            // 1. Configure the Viewer
            const spec = DefaultPluginUISpec();
            
            spec.layout = {
                initialShowControls: false,
                initialShowRemoteState: false,
                initialShowSequence: false,
                initialShowLog: false,
                // Keep the rest of the UI flexible
                controlsDisplay: 'reactive', 
            };
            
            // Be careful with disabling components; sometimes UI parts depend on them.
            // Safe to disable remoteState, but if errors persist, try removing this block.
            spec.components = {
                remoteState: 'none'
            };

            // 2. Initialize using createPluginUI
            // This replaces `new PluginContext`, `ctx.init()`, and `createRoot`
            try {
                const ctx = await createPluginUI({
                    target: parentRef.current,
                    spec: spec,
                    render: renderReact18 // Pass this if using React 18, otherwise Mol* defaults to ReactDOM.render
                });

                pluginRef.current = ctx;
                setPluginCtx(ctx);
            } catch (error) {
                console.error("Failed to initialize Mol*", error);
            }
        }

        init();

        // Cleanup
        return () => {
            // In React 18 Strict Mode, this might run immediately. 
            // In a real app, you might want to debounce disposal or check if component is truly unmounting.
            if (pluginRef.current) {
                pluginRef.current.dispose();
                pluginRef.current = null;
                setPluginCtx(null);
            }
        };
    }, []);

    const handleFileUpload = async (e) => {
        const file = e.target.files[0];
        if (!file || !pluginCtx) return;
        setFileName(file.name);

        // Clear previous structure
        pluginCtx.clear();

        const reader = new FileReader();
        reader.onload = async (event) => {
            const data = event.target.result;
            
            try {
                // Parse and Visualize the PDB Data
                const task = pluginCtx.builders.structure.parse({ 
                    data, 
                    format: 'pdb' 
                });
                
                const trajectory = await pluginCtx.runTask(task);
                
                // Use 'default' preset to apply cartoon/ball-and-stick representation automatically
                await pluginCtx.builders.structure.hierarchy.applyPreset(trajectory, 'default');
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
                    <span style={{ fontSize: '0.8rem', color: 'var(--text-dim)' }}>{fileName}</span>
                    <label className="primary-btn" style={{ cursor: 'pointer', padding: '8px 16px', background: '#333', color: '#fff', borderRadius: '4px' }}>
                        Upload PDB
                        <input type="file" accept=".pdb" onChange={handleFileUpload} style={{ display: 'none' }} />
                    </label>
                </div>
            </div>

            {/* Container for the Mol* Viewer */}
            <div 
                ref={parentRef} 
                style={{ 
                    width: '100%', 
                    height: '600px', 
                    border: '1px solid #ccc', 
                    borderRadius: '8px', 
                    overflow: 'hidden',
                    background: '#000',
                    position: 'relative', // Essential for absolute positioning of Mol* UI
                    zIndex: 0
                }} 
            />
        </div>
    );
};

export default Docking;