import React, { useState, useEffect, useRef } from 'react';
import './App.css';

// Simple Inline CSS for Tabs to keep it in one file
const styles = {
    tabContainer: { display: 'flex', gap: '10px', marginBottom: '20px', borderBottom: '1px solid #333', paddingBottom: '10px' },
    tabButton: (isActive) => ({
        padding: '10px 20px',
        cursor: 'pointer',
        border: 'none',
        background: isActive ? '#3b82f6' : '#27272a',
        color: isActive ? '#fff' : '#9ca3af',
        borderRadius: '5px',
        fontWeight: 'bold'
    }),
    panel: { animation: 'fadeIn 0.3s ease-in-out' }
};

function App() {
  const [activeTab, setActiveTab] = useState('discovery'); // 'discovery', 'training', 'retro'
  
  // -- DISCOVERY STATE --
  const [disease, setDisease] = useState('');
  const [reviewData, setReviewData] = useState(null);
  const [loadingReview, setLoadingReview] = useState(false);
  const [downloadingTarget, setDownloadingTarget] = useState(null);

  // -- TRAINING STATE --
  const [pdbId, setPdbId] = useState('');
  const [selectedFile, setSelectedFile] = useState(null);
  const [trainingStatus, setTrainingStatus] = useState({
    is_running: false, progress: 0, logs: [], total_ligands: 0, best_r2: 0, viz_image: null, status_text: 'Idle'
  });
  const logsEndRef = useRef(null);

  // -- RETRO STATE --
  const [retroSmiles, setRetroSmiles] = useState('');
  const [retroResults, setRetroResults] = useState(null);
  const [loadingRetro, setLoadingRetro] = useState(false);

  // --- POLLING ---
  useEffect(() => {
    const interval = setInterval(async () => {
      try {
        const res = await fetch('http://localhost:8000/api/training/status');
        const data = await res.json();
        setTrainingStatus(data);
      } catch (e) { 
          // console.error(e); // Suppress log spam if backend is offline
      }
    }, 1000);
    return () => clearInterval(interval);
  }, []);

  useEffect(() => {
    if(activeTab === 'training') logsEndRef.current?.scrollIntoView({ behavior: "smooth" });
  }, [trainingStatus.logs, activeTab]);

  // --- ACTIONS ---

  const fetchReview = async () => {
    if (!disease) return;
    setLoadingReview(true);
    try {
      const res = await fetch('http://localhost:8000/api/review', {
        method: 'POST', headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({ disease_name: disease })
      });
      const data = await res.json();
      setReviewData(data);
    } catch (e) { alert("Error fetching review"); }
    setLoadingReview(false);
  };

  const fetchChemblData = async (targetName) => {
    setDownloadingTarget(targetName);
    try {
      const res = await fetch('http://localhost:8000/api/chembl/download', {
        method: 'POST', headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({ target_name: targetName })
      });
      const data = await res.json();
      downloadCSV(data.compounds, `${targetName}_chembl.csv`);
    } catch(e) { alert("Error finding data"); }
    setDownloadingTarget(null);
  };

  const downloadCSV = (data, filename) => {
    const headers = Object.keys(data[0]).join(",");
    const rows = data.map(obj => Object.values(obj).map(val => `"${val}"`).join(","));
    const csv = "data:text/csv;charset=utf-8," + [headers, ...rows].join("\n");
    const link = document.createElement("a");
    link.href = encodeURI(csv); link.download = filename;
    document.body.appendChild(link); link.click(); document.body.removeChild(link);
  };

  const handleTrainUpload = async () => {
    if (!selectedFile || !pdbId) return alert("Missing PDB ID or CSV.");
    const formData = new FormData();
    formData.append("file", selectedFile);
    formData.append("pdb_id", pdbId);
    await fetch("http://localhost:8000/api/train", { method: "POST", body: formData });
  };

  const runRetrosynthesis = async () => {
    if (!retroSmiles) return;
    setLoadingRetro(true);
    setRetroResults(null);
    try {
        const res = await fetch('http://localhost:8000/api/retrosynthesis', {
            method: 'POST', headers: {'Content-Type': 'application/json'},
            body: JSON.stringify({ smiles: retroSmiles })
        });
        const data = await res.json();
        
        if(data.error) {
            alert("Error: " + data.error);
        } else if(data.detail) {
            alert(data.detail);
        } else {
            setRetroResults(data);
        }
    } catch(e) { alert("Failed to run retrosynthesis"); }
    setLoadingRetro(false);
  };

  // Helper to safely extract score
  const getScore = (route) => {
      if (typeof route.score === 'number') return route.score;
      if (route.scores && typeof route.scores.score === 'number') return route.scores.score;
      return 0; // Default if not found
  };

  return (
    <div className="container">
      <header style={{marginBottom: '20px'}}>
        <h1>AI Drug Architect</h1>
        <p style={{color: '#9ca3af'}}>Target Discovery • Binding Optimization • Synthetic Route Planning</p>
      </header>

      {/* --- TABS --- */}
      <div style={styles.tabContainer}>
        <button onClick={() => setActiveTab('discovery')} style={styles.tabButton(activeTab === 'discovery')}>
          1. Discovery
        </button>
        <button onClick={() => setActiveTab('training')} style={styles.tabButton(activeTab === 'training')}>
          2. Training Engine
        </button>
        <button onClick={() => setActiveTab('retro')} style={styles.tabButton(activeTab === 'retro')}>
          3. Retrosynthesis
        </button>
      </div>

      {/* --- PANEL 1: DISCOVERY --- */}
      {activeTab === 'discovery' && (
        <div style={styles.panel} className="card">
            <h2>Literature Analysis & Target ID</h2>
            <div style={{ display: 'flex', gap: '10px' }}>
                <input 
                    type="text" 
                    placeholder="Enter Disease (e.g. Glioblastoma)" 
                    value={disease}
                    onChange={(e) => setDisease(e.target.value)}
                />
                <button onClick={fetchReview} disabled={loadingReview}>
                    {loadingReview ? "Analyzing..." : "Find Targets"}
                </button>
            </div>
            
            {reviewData && (
                <div className="results" style={{marginTop: '20px'}}>
                    <h3>Suggested Targets</h3>
                    <table>
                        <thead>
                            <tr><th>Target</th><th>Type</th><th>Druggability</th><th>Action</th></tr>
                        </thead>
                        <tbody>
                            {reviewData.targets.map((t, i) => (
                                <tr key={i}>
                                    <td style={{fontWeight:'bold', color: '#60a5fa'}}>{t.name}</td>
                                    <td>{t.target_type}</td>
                                    <td>{t.druggability}</td>
                                    <td>
                                        <button onClick={() => fetchChemblData(t.name)} disabled={downloadingTarget === t.name}>
                                            {downloadingTarget === t.name ? '...' : 'Get Data'}
                                        </button>
                                    </td>
                                </tr>
                            ))}
                        </tbody>
                    </table>
                </div>
            )}
        </div>
      )}

      {/* --- PANEL 2: TRAINING --- */}
      {activeTab === 'training' && (
        <div style={styles.panel} className="card">
            <h2>GNN-LLM Fine-Tuning</h2>
            <div style={{ display: 'flex', gap: '20px', alignItems: 'center', marginBottom: '20px' }}>
                <input type="text" placeholder="PDB ID (e.g. 5T35)" value={pdbId} onChange={(e) => setPdbId(e.target.value)} style={{ width: '120px' }} />
                <input type="file" accept=".csv" onChange={(e) => setSelectedFile(e.target.files[0])} style={{ width: 'auto' }} />
                <button onClick={handleTrainUpload} disabled={trainingStatus.is_running} 
                    style={{ background: trainingStatus.is_running ? '#27272a' : '#10b981' }}>
                    {trainingStatus.is_running ? "Running..." : "Start Loop"}
                </button>
            </div>

            <div className="metrics-grid">
                <div className="metric-box"><h4>Status</h4><span style={{fontSize:'0.8rem'}}>{trainingStatus.status_text}</span></div>
                <div className="metric-box"><h4>Val R²</h4>
                <span style={{color: '#34d399'}}>
                    {typeof trainingStatus.best_r2 === "number"
                        ? trainingStatus.best_r2.toFixed(3)
                        : "—"}
                    </span>
                </div>
                <div className="metric-box"><h4>Progress</h4><span>{trainingStatus.progress}%</span></div>
            </div>

            {/* Progress Bar */}
            <div className="progress-container" style={{height:'10px', background:'#333', borderRadius:'5px', margin:'10px 0', overflow:'hidden'}}>
                <div style={{width: `${trainingStatus.progress}%`, background:'#3b82f6', height:'100%', transition:'width 0.5s'}}></div>
            </div>

            {trainingStatus.viz_image && (
                <div style={{margin: '20px 0', textAlign: 'center'}}>
                    <img src={trainingStatus.viz_image} alt="Viz" style={{maxWidth: '100%', borderRadius: '8px', border: '1px solid #444'}} />
                </div>
            )}

            <div className="terminal-window" style={{height: '200px'}}>
                {trainingStatus.logs.map((log, i) => <div key={i} className="log-line">{log}</div>)}
                <div ref={logsEndRef} />
            </div>
        </div>
      )}

      {/* --- PANEL 3: RETROSYNTHESIS --- */}
      {activeTab === 'retro' && (
        <div style={styles.panel} className="card">
            <h2>AiZynthFinder Retrosynthesis</h2>
            <div style={{ display: 'flex', gap: '10px' }}>
                <input 
                    type="text" 
                    placeholder="Enter SMILES (e.g. CC(=O)Oc1ccccc1C(=O)O)" 
                    value={retroSmiles}
                    onChange={(e) => setRetroSmiles(e.target.value)}
                />
                <button onClick={runRetrosynthesis} disabled={loadingRetro}>
                    {loadingRetro ? "Calculating..." : "Plan Route"}
                </button>
            </div>

            {retroResults && (
                <div style={{marginTop: '20px'}}>
                    
                    {/* Status Badge */}
                    {retroResults.is_solved !== undefined && (
                        <div style={{marginBottom: '15px'}}>
                            <span className={`tag ${retroResults.is_solved ? 'tag-high' : 'tag-med'}`}>
                                {retroResults.is_solved ? "ROUTE FOUND" : "NO ROUTE"}
                            </span>
                        </div>
                    )}

                    {/* Routes List */}
                    {retroResults.routes?.map((route, i) => (
                        <div key={i} style={{background: '#18181b', padding: '15px', borderRadius: '8px', marginBottom: '20px', border: '1px solid #333'}}>
                            
                            <h4 style={{margin: '0 0 10px 0', color: '#60a5fa', borderBottom: '1px solid #333', paddingBottom: '10px'}}>
                                Route {i+1} 
                                <span style={{fontSize: '0.8rem', color: '#9ca3af', marginLeft: '10px'}}>
                                    (Score: {getScore(route).toFixed(2)})
                                </span>
                            </h4>

                            {/* Image Display */}
                            {route.image_url && (
                                <div style={{textAlign: 'center', margin: '15px 0', background: '#fff', borderRadius: '5px', padding: '5px'}}>
                                    <img 
                                        src={route.image_url} 
                                        alt={`Route ${i+1}`} 
                                        style={{maxWidth: '100%', maxHeight: '400px'}} 
                                    />
                                </div>
                            )}

                            {/* Steps List */}
                            <div style={{fontSize: '0.9rem'}}>
                                <h5 style={{marginTop: '10px', color: '#e5e7eb'}}>Reaction Steps:</h5>
                                {Array.isArray(route.steps) && route.steps.map((step, s) => {

                                    return (
                                        <div key={s} style={{marginBottom: '10px', padding: '8px', background: '#27272a', borderRadius: '4px'}}>
                                            <div style={{fontWeight: 'bold', color: '#34d399', marginBottom: '4px'}}>
                                                Step {s+1}
                                            </div>
                                            <div style={{color: '#d1d5db', wordBreak: 'break-all'}}>
                                                <span style={{color: '#9ca3af'}}>Reactants: </span> 
                                                {step.reactants}
                                            </div>
                                            <div style={{color: '#6b7280', fontSize: '0.8rem', marginTop: '2px'}}>
                                                Template: {step.reaction_smarts
                                                    ? (step.reaction_smarts.length > 50
                                                        ? step.reaction_smarts.substring(0, 50) + "..."
                                                        : step.reaction_smarts)
                                                    : "N/A"}
                                            </div>
                                        </div>
                                    );
                                })}
                            </div>

                        </div>
                    ))}
                </div>
            )}
        </div>
      )}

    </div>
  );
}

export default App;