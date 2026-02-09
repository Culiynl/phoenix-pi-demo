import React, { useState } from 'react';

const FileUploadPanel = ({ missionId }) => {
    const [uploading, setUploading] = useState(false);
    const [uploadResult, setUploadResult] = useState(null);

    const handleFileUpload = async (event, fileType) => {
        const file = event.target.files[0];
        if (!file) return;

        setUploading(true);
        setUploadResult(null);

        const formData = new FormData();
        formData.append('file', file);
        formData.append('file_type', fileType);

        try {
            const response = await fetch(`http://localhost:8000/api/mission/${missionId}/upload`, {
                method: 'POST',
                body: formData
            });

            const result = await response.json();
            setUploadResult({ success: true, ...result });
        } catch (error) {
            setUploadResult({ success: false, error: error.message });
        } finally {
            setUploading(false);
        }
    };

    return (
        <div className="upload-panel">
            <h4>Multimodal Analysis</h4>
            <p className="upload-desc">Upload protein structure images or research PDFs for AI analysis</p>

            <div className="upload-grid">
                <div className="upload-card">
                    <div className="upload-icon">🧬</div>
                    <h5>Protein Structure Image</h5>
                    <p>Gemini will identify binding sites and druggable pockets</p>
                    <label className="upload-btn">
                        <input
                            type="file"
                            accept="image/*"
                            onChange={(e) => handleFileUpload(e, 'image')}
                            disabled={uploading}
                        />
                        {uploading ? 'Uploading...' : 'Choose Image'}
                    </label>
                </div>

                <div className="upload-card">
                    <div className="upload-icon">📄</div>
                    <h5>Research Paper PDF</h5>
                    <p>Extract molecules, targets, and experimental data</p>
                    <label className="upload-btn">
                        <input
                            type="file"
                            accept=".pdf"
                            onChange={(e) => handleFileUpload(e, 'pdf')}
                            disabled={uploading}
                        />
                        {uploading ? 'Uploading...' : 'Choose PDF'}
                    </label>
                </div>
            </div>

            {uploadResult && (
                <div className={`upload-result ${uploadResult.success ? 'success' : 'error'}`}>
                    {uploadResult.success ? (
                        <>
                            <strong>✓ Uploaded:</strong> {uploadResult.filename}
                            <p>File is ready for analysis. Ask PHOENIX to analyze it!</p>
                        </>
                    ) : (
                        <>
                            <strong>✗ Error:</strong> {uploadResult.error}
                        </>
                    )}
                </div>
            )}

            <style>{`
                .upload-panel {
                    background: #0d0d0d;
                    border: 1px solid #1a1a1a;
                    border-radius: 12px;
                    padding: 25px;
                    margin-bottom: 20px;
                }
                .upload-panel h4 {
                    margin: 0 0 8px 0;
                    color: #4285F4;
                    font-size: 16px;
                }
                .upload-desc {
                    color: #888;
                    font-size: 12px;
                    margin: 0 0 20px 0;
                }
                .upload-grid {
                    display: grid;
                    grid-template-columns: 1fr 1fr;
                    gap: 15px;
                }
                .upload-card {
                    background: #111;
                    border: 1px solid #222;
                    border-radius: 8px;
                    padding: 20px;
                    text-align: center;
                }
                .upload-icon {
                    font-size: 40px;
                    margin-bottom: 10px;
                }
                .upload-card h5 {
                    margin: 0 0 8px 0;
                    color: #fff;
                    font-size: 14px;
                }
                .upload-card p {
                    color: #666;
                    font-size: 11px;
                    margin: 0 0 15px 0;
                    line-height: 1.4;
                }
                .upload-btn {
                    display: inline-block;
                    background: #4285F4;
                    color: #fff;
                    padding: 8px 20px;
                    border-radius: 6px;
                    font-size: 12px;
                    font-weight: bold;
                    cursor: pointer;
                    transition: background 0.2s;
                }
                .upload-btn:hover {
                    background: #3367d6;
                }
                .upload-btn input {
                    display: none;
                }
                .upload-result {
                    margin-top: 15px;
                    padding: 12px;
                    border-radius: 6px;
                    font-size: 12px;
                }
                .upload-result.success {
                    background: rgba(76, 175, 80, 0.1);
                    border: 1px solid #4CAF50;
                    color: #4CAF50;
                }
                .upload-result.error {
                    background: rgba(244, 67, 54, 0.1);
                    border: 1px solid #f44336;
                    color: #f44336;
                }
                .upload-result strong {
                    display: block;
                    margin-bottom: 5px;
                }
                .upload-result p {
                    margin: 5px 0 0 0;
                    color: #888;
                }
            `}</style>
        </div>
    );
};

export default FileUploadPanel;
