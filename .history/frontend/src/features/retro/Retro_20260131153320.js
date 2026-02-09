import React, { useState } from 'react';
import { api } from '../../services/api';

const Retro = () => {
  const [smiles, setSmiles] = useState('');
  const [results, setResults] = useState(null);
  const [loading, setLoading] = useState(false);

  const handleRun = async () => {
    setLoading(true);
    try {
      const data = await api.runRetrosynthesis(smiles);
      setResults(data);
    } catch (e) { alert("Check backend connection"); }
    setLoading(false);
  };

  return (
  <div className="card">
    <h2>Retrosynthesis Planner</h2>
    
    <div className="action-group">
      <input 
        type="text" 
        className="technical-input"
        placeholder="Enter SMILES (e.g. CC(=O)Oc1ccccc1C(=O)O)" 
        value={smiles}
        onChange={(e) => setSmiles(e.target.value)}
      />
      <button className="primary-btn" onClick={handleRun} disabled={loading}>
        {loading ? "Calculating..." : "Plan Route"}
      </button>
    </div>


      {results && results.routes?.map((route, i) => (
        <div key={i} className="route-card" style={{background: '#1a1a1a', padding: '15px', borderRadius: '8px', marginTop: '10px'}}>
           <h4>Route {i+1} - Steps: {route.scores?.["number of reactions"]}</h4>
           {route.image_url && <img src={route.image_url} alt="route" style={{width: '100%', background: 'white', borderRadius: '4px'}} />}
        </div>
      ))}
    </div>
  );
};

export default Retro;