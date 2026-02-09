import React from 'react';
import StandaloneTools from '../../components/StandaloneTools';

const ToolsPage = () => {
    return (
        <div className="tools-page">
            <div className="page-header">
                <h1>🔬 Research Tools</h1>
                <p>Standalone computational chemistry tools for quick analysis</p>
            </div>

            <div className="tools-container">
                <StandaloneTools />
            </div>

            <style>{`
                .tools-page {
                    padding: 30px;
                    max-width: 1200px;
                    margin: 0 auto;
                }
                .page-header {
                    margin-bottom: 30px;
                }
                .page-header h1 {
                    margin: 0 0 10px 0;
                    color: #fff;
                    font-size: 32px;
                }
                .page-header p {
                    margin: 0;
                    color: #888;
                    font-size: 16px;
                }
                .tools-container {
                    display: grid;
                    gap: 20px;
                }
            `}</style>
        </div>
    );
};

export default ToolsPage;
