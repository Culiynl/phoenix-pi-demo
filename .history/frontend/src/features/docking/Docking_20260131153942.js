import React, { useEffect, useRef, useState } from 'react';
import { createPluginUI } from 'molstar/lib/mol-plugin-ui/react18';
import { DefaultPluginUISpec } from 'molstar/lib/mol-plugin-ui/spec';
import 'molstar/lib/mol-plugin-ui/skin/dark.css'; // Fits your dark theme

const Docking = () => {
    const parentRef = useRef(null);
    const pluginRef = useRef(null);
    const [fileName, setFileName] = useState("No file selected");

    // Initialize Mol* Plugin
    useEffect(() => {
        async function init() {
            const spec = DefaultPluginUISpec();
            // Customizing the UI to be more minimal
            spec.layout = {
                initialShowControls: false,
                initialShowRemoteState: false,
            };

            const plugin = await createPluginUI(parentRef.current, spec);
            pluginRef.current = plugin;
        }
        init();

        return () => {
            if (pluginRef.current) pluginRef.current.dispose();
        };
    }, []);

    // Handle File Upload & Visualization
    const handleFileUpload = async (e) => {
        const file = e.target.files[0];
        if (!file) return;
        setFileName(file.name);

        const plugin = pluginRef.current;
        if (!plugin) return;

        // Clear previous structures
        plugin.clear();

        // Read file as string
        const reader = new FileReader();
        reader.onload = async (event) => {
            const data = event.target.result;
            
            // Load the PDB data into Mol*
            const task = plugin.builders.structure.parse({ 
                data, 
                format: 'pdb' 
            });
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

            {/* Mol* Viewer Container */}
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
            
            <div style={{ marginTop: '15px', fontSize: '0.8rem', color: 'var(--text-dim)' }}>
                <strong>Controls:</strong> Left Click = Rotate | Right Click = Pan | Scroll = Zoom
            </div>
        </div>
    );
};

export default Docking;