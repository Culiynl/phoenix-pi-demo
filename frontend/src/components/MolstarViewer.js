import React, { useEffect, useRef, useState, useLayoutEffect } from 'react';
import { createPluginUI } from 'molstar/lib/mol-plugin-ui';
import { renderReact18 } from 'molstar/lib/mol-plugin-ui/react18';
import { DefaultPluginUISpec } from 'molstar/lib/mol-plugin-ui/spec';
import 'molstar/lib/mol-plugin-ui/skin/light.scss';
import { Color } from 'molstar/lib/mol-util/color';

const MolstarViewer = ({ pdbId, ligandUrl, receptorUrl, customPdbData, visualState }) => {
    const parentRef = useRef(null);
    const pluginRef = useRef(null);
    const [initialized, setInitialized] = useState(false);

    useLayoutEffect(() => {
        let isCurrent = true;
        const currentParent = parentRef.current;
        async function init() {
            if (!currentParent) return;
            currentParent.innerHTML = '';

            const spec = {
                ...DefaultPluginUISpec(),
                layout: {
                    initial: {
                        isExpanded: false,
                        showControls: false,
                        regionState: { left: 'hidden', top: 'hidden', right: 'hidden', bottom: 'hidden' }
                    }
                },
                components: { remoteState: 'none' }
            };

            try {
                const ctx = await createPluginUI({
                    target: currentParent,
                    spec: spec,
                    render: renderReact18
                });

                if (!isCurrent) {
                    ctx.dispose();
                    return;
                }

                // Set Background to White as requested
                ctx.canvas3d?.setProps({
                    renderer: { ...ctx.canvas3d.props.renderer, backgroundColor: Color(0xFFFFFF) }
                });

                pluginRef.current = ctx;
                setInitialized(true);
            } catch (error) {
                if (!error.message?.includes('already added')) console.error(error);
            }
        }
        init();
        return () => {
            isCurrent = false;
            if (pluginRef.current) pluginRef.current.dispose();
            if (currentParent) currentParent.innerHTML = '';
        };
    }, []);

    useEffect(() => {
        if (!initialized || !pluginRef.current) return;

        const load = async () => {
            const plugin = pluginRef.current;
            await plugin.clear();

            try {
                // ... (existing loading logic) ...
                let data;
                // 1. Raw PDB
                if (customPdbData && !pdbId) {
                    data = await plugin.builders.data.rawData({ data: customPdbData });
                }
                // 2. RCSB PDB
                else if (pdbId) {
                    const url = `https://files.rcsb.org/download/${pdbId.toUpperCase()}.pdb`;
                    data = await plugin.builders.data.download({ url, isBinary: false });
                }
                // 3. Docking Results (Receptor + Ligand)
                else if (receptorUrl && ligandUrl) {
                    // const baseUrl = 'http://localhost:8000';
                    const baseUrl = 'http://localhost:8000';
                    // Load Receptor
                    const r = await plugin.builders.data.download({ url: `${baseUrl}${receptorUrl}`, label: 'Receptor' });
                    const rTraj = await plugin.builders.structure.parseTrajectory(r, 'pdb');
                    await plugin.builders.structure.hierarchy.applyPreset(rTraj, 'default');

                    // Load Ligand
                    const l = await plugin.builders.data.download({ url: `${baseUrl}${ligandUrl}`, label: 'Ligand (Best Pose)' });
                    const lTraj = await plugin.builders.structure.parseTrajectory(l, 'pdb');
                    const lComp = await plugin.builders.structure.hierarchy.applyPreset(lTraj, 'default');

                    // Focus on Ligand
                    plugin.managers.camera.focusLoci(lComp.representation.Loci);
                    return;
                }

                if (data) {
                    const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
                    await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default');
                }
            } catch (e) { console.warn("Load failed", e); }
        };
        load();
    }, [initialized, pdbId, ligandUrl, receptorUrl, customPdbData]);

    // Programmatic Control Effect
    useEffect(() => {
        if (!pluginRef.current || !visualState) return;

        const applyState = async () => {
            const plugin = pluginRef.current;
            // Example: { type: 'focus', target: 'ligand' } or { type: 'zoom', residue: 45 }
            if (visualState.type === 'focus' && visualState.target === 'ligand') {
                // Try to find ligand component and focus
                // Simplified: Reset camera if no specific ligand tracking yet
                plugin.managers.camera.reset();
            }
            if (visualState.type === 'zoom' && visualState.residue) {
                // Implement residue focus logic (requires structure traversal, skipping for MVP robustness)
                console.log("Zoom request for residue", visualState.residue);
            }
            if (visualState.type === 'highlight') {
                // Highlight logic
            }
        };
        applyState();
    }, [visualState]);

    return (
        <div className="molstar-wrapper-fixed">
            <div ref={parentRef} className="molstar-host-canvas" />
            <style>{`
                .molstar-wrapper-fixed {
                    width: 100%;
                    height: 100%;
                    background: #fff; /* Match requested white background */
                    position: relative;
                    overflow: hidden;
                }
                .molstar-host-canvas {
                    width: 100%;
                    height: 100%;
                }
                .msp-plugin-ui-main {
                    background: #fff !important;
                }
                /* Hide glitchy UI overlays */
                .msp-help, .msp-logo { display: none !important; }
            `}</style>
        </div>
    );
};

export default MolstarViewer;