import React, { useEffect, useState } from 'react';
import { Plugin } from 'molstar/lib/mol-plugin-ui/plugin';
import { PluginContext } from 'molstar/lib/mol-plugin/context';
import { DefaultPluginUISpec } from 'molstar/lib/mol-plugin-ui/spec';
import 'molstar/build/viewer/molstar.css';

const Docking = () => {
    const [plugin, setPlugin] = useState(null);
    const [fileName, setFileName] = useState("No file selected");

    useEffect(() => {
        async function init() {
            // 1. Create the spec
            const spec = DefaultPluginUISpec();
            spec.layout = {
                initialShowControls: false,
                initialShowRemoteState: false,
            };

            // 2. Initialize the context
            const ctx = new PluginContext(spec);
            await ctx.init();
            
            // 3. Store the context in state
            setPlugin(ctx);
        }

        init();

        // Cleanup on unmount
        return () => {
            if (plugin) plugin.dispose();
        };
    }, []);

    const handleFileUpload = async (e) => {
        const file = e.target.files[0];
        if (!file || !plugin) return;
        setFileName(file.name);

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

            <div style={{ 
                width: '100%', 
                height: '600px', 
                border: '1px solid var(--border)', 
                borderRadius: '8px', 
                overflow: 'hidden',
                background: '#000',
                position: 'relative'
            }}>
                {/* 
                  Instead of a Ref, we use the <Plugin /> component.
                  We only render it once the 'plugin' context is ready.
                */}
                {plugin && <Plugin plugin={plugin} />}
            </div>
        </div>
    );
};

export default Docking;