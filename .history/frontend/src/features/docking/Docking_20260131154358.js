import React, { useEffect, useRef, useState } from 'react';
// Import from the UI core, NOT the heavy "apps" folder
import { createPluginUI } from 'molstar/lib/mol-plugin-ui';
import { DefaultPluginUISpec } from 'molstar/lib/mol-plugin-ui/spec';
import 'molstar/lib/mol-plugin-ui/skin/dark.scss'; // We can use SCSS now that sass is installed

const Docking = () => {
    const parentRef = useRef(null);
    const pluginRef = useRef(null);
    const [fileName, setFileName] = useState("No file selected");
    const [initialized, setInitialized] = useState(false);

    useEffect(() => {
        // Prevent double-initialization in React 18 Strict Mode
        if (initialized || !parentRef.current) return;

        async function init() {
            setInitialized(true); // Mark as running

            const spec = DefaultPluginUISpec();
            
            // Clean up the UI options
            spec.layout = {
                initialShowControls: false,
                initialShowRemoteState: false,
                initialShowSequence: false,
                initialShowLog: false,
            };

            spec.components = {
                remoteState: 'none' // Disable remote server features
            };

            try {
                // mount the UI into the parentRef div
                const plugin = await createPluginUI(parentRef.current, spec);
                pluginRef.current = plugin;
                
                // Force dark background
                const canvas = parentRef.current.querySelector('canvas');
                if (canvas) canvas.style.background = '#000';
                
            } catch (e) {
                console.error("Mol* Init Error:", e);
            }
        }

        init();

        // Cleanup
        return () => {
            // We generally don't dispose pluginRef here in StrictMode
            // because it can cause race conditions. We let it live 
            // since the component is main-page persistence.
        };
    }, [initialized]);

    const handleFileUpload = async (e) => {
        const file = e.target.files[0];
        if (!file || !pluginRef.current) return;
        setFileName(file.name);

        const plugin = pluginRef.current;
        plugin.clear();

        const reader = new FileReader();
        reader.onload = async (event) => {
            const data = event.target.result;
            const task = plugin.builders.structure.parse({ data, format: 'pdb' });
            const trajectory = await plugin.runTask(task);
            await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default');
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
                    zIndex: 0 /* Ensures it doesn't overlap navbar */
                }} 
            />
        </div>
    );
};

export default Docking;