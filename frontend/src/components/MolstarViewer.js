import React, { useEffect, useRef, useState } from 'react';
import { renderReact18 } from 'molstar/lib/mol-plugin-ui/react18';
import { DefaultPluginUISpec } from 'molstar/lib/mol-plugin-ui/spec';
import 'molstar/lib/mol-plugin-ui/skin/light.scss';

const MolstarViewer = ({ pdbId, ligandUrl, receptorUrl }) => {
    const parentRef = useRef(null);
    const pluginRef = useRef(null);
    const [initError, setInitError] = useState(null);

    // 1️⃣ Initialize Molstar ONCE (no dependencies on props)
    useEffect(() => {
        let mounted = true;

        const spec = {
            ...DefaultPluginUISpec(),
            layout: {
                initial: {
                    isExpanded: false,
                    showControls: false,
                    regionState: {
                        left: 'hidden',
                        top: 'hidden',
                        right: 'hidden',
                        bottom: 'hidden'
                    }
                }
            }
        };

        const init = async () => {
            // Wait for React to mount the DOM element
            await new Promise(resolve => setTimeout(resolve, 0));

            if (!mounted || !parentRef.current || pluginRef.current) return;

            try {
                const plugin = await renderReact18(parentRef.current, spec);
                if (mounted) {
                    pluginRef.current = plugin;
                }
            } catch (err) {
                console.error("Molstar init failed:", err);
                if (mounted) {
                    setInitError(err.message);
                }
            }
        };

        init();

        return () => {
            mounted = false;
            if (pluginRef.current) {
                try {
                    pluginRef.current.dispose();
                } catch (e) {
                    console.warn("Molstar disposal error:", e);
                }
                pluginRef.current = null;
            }
        };
    }, []); // ✅ Empty deps - init once only

    // 2️⃣ Load structures separately when props change
    useEffect(() => {
        if (!pluginRef.current) return;

        const loadStructures = async () => {
            const plugin = pluginRef.current;

            try {
                // Clear previous structures
                plugin.clear();

                if (pdbId) {
                    const url = `http://localhost:8000/static/pdb_files/${pdbId}_clean.pdb`;
                    await plugin.builders.structure.download({ url }, 'pdb');
                } else if (receptorUrl && ligandUrl) {
                    const receptor = await plugin.builders.structure.download({
                        url: `http://localhost:8000${receptorUrl}`
                    }, 'pdb');
                    const ligand = await plugin.builders.structure.download({
                        url: `http://localhost:8000${ligandUrl}`
                    }, 'pdb');

                    if (receptor && ligand) {
                        await plugin.builders.structure.representation.addRepresentation(
                            receptor.root,
                            { type: 'cartoon', color: 'chain-id' }
                        );
                        await plugin.builders.structure.representation.addRepresentation(
                            ligand.root,
                            { type: 'ball-and-stick', color: 'element-symbol' }
                        );
                    }
                }
            } catch (err) {
                console.error("Structure loading failed:", err);
            }
        };

        loadStructures();
    }, [pdbId, ligandUrl, receptorUrl]); // ✅ Reload data when props change

    return (
        <div className="molstar-wrapper">
            <div
                ref={parentRef}
                style={{
                    width: '100%',
                    height: '100%',
                    position: 'relative'
                }}
            />

            {initError && (
                <div className="init-fallback">
                    <p>3D Viewer Error</p>
                    <code>{initError}</code>
                </div>
            )}

            <style>{`
                .molstar-wrapper {
                    width: 100%;
                    height: 100%;
                    position: relative;
                    background: #000;
                    border-radius: 8px;
                    overflow: hidden;
                    border: 1px solid #1a1a1a;
                }
                .init-fallback {
                    position: absolute;
                    inset: 0;
                    display: flex;
                    flex-direction: column;
                    align-items: center;
                    justify-content: center;
                    background: rgba(0,0,0,0.9);
                    z-index: 100;
                    padding: 20px;
                    text-align: center;
                }
                .init-fallback p {
                    color: #ff4444;
                    font-weight: bold;
                    margin-bottom: 10px;
                }
                .init-fallback code {
                    font-size: 11px;
                    color: #888;
                    background: #111;
                    padding: 8px;
                    border-radius: 4px;
                }
                .msp-plugin-ui {
                    background: #000 !important;
                }
            `}</style>
        </div>
    );
};

// Wrapper to opt out of StrictMode for Molstar
const NonStrict = ({ children }) => <>{children}</>;

const MolstarViewerWrapper = (props) => (
    <NonStrict>
        <MolstarViewer {...props} />
    </NonStrict>
);

export default MolstarViewerWrapper;
