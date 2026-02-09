import React from 'react';

const GeminiAdvantageShowcase = ({ missionData }) => {
    // Only show real metrics from actual mission data
    const contextUsage = missionData?.context_tokens || 0;
    const maxTokens = 2000000;
    const contextPercent = Math.min((contextUsage / maxTokens) * 100, 100);

    const papersAnalyzed = missionData?.research_papers?.length || 0;
    const moleculesAnalyzed = missionData?.chembl_summary?.molecules_count || 0;
    const toolCallsExecuted = missionData?.thought_log?.filter(log => log.type === 'action')?.length || 0;

    // Don't show if no meaningful data yet
    if (contextUsage === 0 && papersAnalyzed === 0 && moleculesAnalyzed === 0) {
        return null;
    }

    return (
        <div className="gemini-showcase">
            <h3 className="showcase-title">⚡ Gemini 3.0 Capabilities in Action</h3>

            <div className="metrics-panel">
                {contextUsage > 0 && (
                    <div className="metric-card">
                        <h5>Long Context Usage</h5>
                        <div className="progress-bar">
                            <div className="progress-fill" style={{ width: `${contextPercent}%` }}>
                                {contextUsage.toLocaleString()} / 2M tokens
                            </div>
                        </div>
                        <p className="metric-note">
                            {papersAnalyzed > 0 && `${papersAnalyzed} papers analyzed`}
                            {moleculesAnalyzed > 0 && ` • ${moleculesAnalyzed} molecules processed`}
                        </p>
                    </div>
                )}

                {toolCallsExecuted > 0 && (
                    <div className="metric-card">
                        <h5>Autonomous Tool Calls</h5>
                        <div className="metric-value">{toolCallsExecuted}</div>
                        <p className="metric-note">Function calls executed by Gemini</p>
                    </div>
                )}

                {missionData?.current_step && (
                    <div className="metric-card">
                        <h5>Workflow Progress</h5>
                        <div className="metric-value">Step {missionData.current_step} / 7</div>
                        <p className="metric-note">Autonomous multi-step execution</p>
                    </div>
                )}
            </div>

            <style>{`
                .gemini-showcase {
                    background: linear-gradient(135deg, #0d0d0d 0%, #1a1a2e 100%);
                    border: 1px solid #2a2a3e;
                    border-radius: 16px;
                    padding: 25px;
                    margin-bottom: 25px;
                }
                .showcase-title {
                    margin: 0 0 20px 0;
                    font-size: 18px;
                    background: linear-gradient(90deg, #4285F4, #0F9D58, #F4B400);
                    -webkit-background-clip: text;
                    -webkit-text-fill-color: transparent;
                    background-clip: text;
                }
                .metrics-panel {
                    display: grid;
                    grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
                    gap: 15px;
                }
                .metric-card {
                    background: #0a0a0a;
                    border: 1px solid #1a1a1a;
                    border-radius: 8px;
                    padding: 15px;
                }
                .metric-card h5 {
                    margin: 0 0 10px 0;
                    font-size: 12px;
                    color: #888;
                    text-transform: uppercase;
                }
                .progress-bar {
                    background: #1a1a1a;
                    border-radius: 4px;
                    overflow: hidden;
                    height: 24px;
                }
                .progress-fill {
                    background: linear-gradient(90deg, #4285F4, #0F9D58);
                    height: 100%;
                    display: flex;
                    align-items: center;
                    justify-content: center;
                    font-size: 10px;
                    color: #fff;
                    font-weight: bold;
                    transition: width 0.3s ease;
                    min-width: 80px;
                }
                .metric-value {
                    font-size: 32px;
                    font-weight: bold;
                    color: #4285F4;
                    margin-bottom: 5px;
                }
                .metric-note {
                    font-size: 10px;
                    color: #666;
                    margin: 5px 0 0 0;
                }
            `}</style>
        </div>
    );
};

export default GeminiAdvantageShowcase;
