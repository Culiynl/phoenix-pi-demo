import React, { useEffect, useRef, useState } from 'react';
import { Viewer } from 'molstar/lib/apps/viewer';
import 'molstar/build/viewer/molstar.css';

const Docking = () => {
    const viewerWrapperRef = useRef(null);
    const viewerInstance = useRef(null);
    const [fileName, setFileName] = useState("No file selected");

    useEffect(() => {
        async function initViewer() {
            if (!viewerWrapperRef.current || viewerInstance.current) return;

            // Initialize the Viewer in the div
            // This API is much more stable for React 18
            const viewer = new Viewer(viewerWrapperRef.current, {
                layoutIsExpanded: false,
                layoutShowControls: false,
                layoutShowRemoteState: false,
                layoutShowSequence: false,
                layoutShowLog: false,
                viewportShowExpand: false,
                viewportShowSelectionMode: false,
                viewportShowAnimation: false,
                collapseLeftPanel: true,
                pdbProvider: 'rcsb',
                emdbProvider: 'rcsb',
            });

            viewerInstance.current = viewer;
        }

        initViewer();

        // Cleanup function
        return () => {
            if (viewerInstance.current) {
                viewerInstance.current.plugin.dispose();
                viewerInstance.current = null;
            }
        };
    }, []);

    const handleFileUpload = async (e) => {
        const file = e.target.files[0];
        if (!file || !viewerInstance.current) return;
        setFileName(file.name);

        const plugin = viewerInstance.current.plugin;
        
        // Clear existing models
        await plugin.clear();

        const reader = new FileReader();
        reader.onload = async (event) => {
            const data = event.target.result;
            
            // Using the Viewer's built-in loading logic
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

            {/* Viewer Container */}
            <div 
                ref={viewerWrapperRef} 
                style={{ 
                    width: '100%', 
                    height: '600px', 
                    border: '1px solid var(--border)', 
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