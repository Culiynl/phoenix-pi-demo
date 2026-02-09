import React, { useEffect, useRef, useState } from 'react';
import { createRoot } from 'react-dom/client'; // Use standard React 18 root
import { PluginContext } from 'molstar/lib/mol-plugin/context';
import { Plugin } from 'molstar/lib/mol-plugin-ui/plugin';
import { DefaultPluginUISpec } from 'molstar/lib/mol-plugin-ui/spec';
import 'molstar/lib/mol-plugin-ui/skin/dark.scss'; // Uses your new sass compiler

const Docking = () => {
    const parentRef = useRef(null);
    const [pluginCtx, setPluginCtx] = useState(null);
    const [fileName, setFileName] = useState("No file selected");

    useEffect(() => {
        if (!parentRef.current) return;

        let root = null;
        let ctx = null;

        async function init() {
            // 1. Configure the Viewer
            const spec = DefaultPluginUISpec();
            spec.layout = {
                initialShowControls: false,
                initialShowRemoteState: false,
                initialShowSequence: false,
                initialShowLog: false,
            };
            spec.components = {
                remoteState: 'none'
            };

            // 2. Initialize the "Brain" (PluginContext)
            ctx = new PluginContext(spec);
            await ctx.init();

            // 3. Mount the "Body" (UI) using standard React 18
            root = createRoot(parentRef.current);
            root.render(<Plugin plugin={ctx} />);

            // 4. Save context to state so we can use it for uploads
            setPluginCtx(ctx);
        }

        init();

        // Cleanup when leaving the page
        return () => {
            if (root) root.unmount();
            if (ctx) ctx.dispose();
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
                    <label className="primary-btn" style={{ cursor: 'pointer' }}>
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
                    border: '1px solid var(--border)', 
                    borderRadius: '8px', 
                    overflow: 'hidden',
                    background: '#000',
                    position: 'relative',
                    zIndex: 0
                }} 
            />
        </div>
    );
};

export default Docking;