import React, { useEffect, useRef, useState } from 'react';
import { PluginContext } from 'molstar/lib/mol-plugin/context';
import { DefaultPluginUISpec } from 'molstar/lib/mol-plugin-ui/spec';
import { renderReact18 } from 'molstar/lib/mol-plugin-ui/react18'; // Use renderReact18 as suggested by the error
import 'molstar/build/viewer/molstar.css'; // Updated CSS path

const Docking = () => {
    const parentRef = useRef(null);
    const pluginRef = useRef(null);
    const [fileName, setFileName] = useState("No file selected");

    useEffect(() => {
        async function init() {
            if (!parentRef.current) return;

            // 1. Define the specification
            const spec = DefaultPluginUISpec();
            spec.layout = {
                initialShowControls: false,
                initialShowRemoteState: false,
            };

            // 2. Create the Plugin Context manually
            const plugin = new PluginContext(spec);
            await plugin.init();

            // 3. Render the UI using the React 18 specific renderer
            await renderReact18(plugin, parentRef.current);
            
            pluginRef.current = plugin;
        }
        
        init();

        return () => {
            if (pluginRef.current) {
                pluginRef.current.dispose();
            }
        };
    }, []);

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
                    position: 'relative',
                    background: '#000'
                }} 
            />
        </div>
    );
};

export default Docking;