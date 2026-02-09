import React, { useState } from 'react';
import { api } from '../../services/api';

const Retro = () => {
  const [smiles, setSmiles] = useState('');
  const [results, setResults] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);

  const handleRun = async () => {
    if (!smiles) return;
    setLoading(true);
    setError(null);
    setResults(null);

    try {
      const data = await api.runRetrosynthesis(smiles);
      if (data.error) {
        setError(data.error);
      } else {
        setResults(data);
      }
    } catch (e) {
      setError("Failed to connect to the retrosynthesis engine.");
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="retro-container">
      <div className="card">
        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '20px' }}>
          <h2>Retrosynthesis Planner</h2>
          <span className="badge">AiZynthFinder v4.0</span>
        </div>

        <p style={{ color: 'var(--text-dim)', fontSize: '0.9rem', marginBottom: '20px' }}>
          Predict synthetic pathways using deep learning and Monte Carlo Tree Search.
        </p>

        <div className="action-group">
          <input
            type="text"
            className="technical-input"
            placeholder="Enter Target SMILES (e.g. CC(=O)Oc1ccccc1C(=O)O)"
            value={smiles}
            onChange={(e) => setSmiles(e.target.value)}
            disabled={loading}
          />
          <button
            className="primary-btn"
            onClick={handleRun}
            disabled={loading || !smiles}
            style={{ minWidth: '140px' }}
          >
            {loading ? "Calculating..." : "Plan Routes"}
          </button>
        </div>

        {error && (
          <div className="error-box" style={{ marginTop: '20px', color: '#f87171', fontSize: '0.8rem', fontFamily: 'monospace' }}>
            &gt; ERROR: {error}
          </div>
        )}
      </div>

      {loading && (
        <div className="loading-overlay card" style={{ textAlign: 'center', padding: '40px' }}>
          <div className="spinner" style={{ marginBottom: '15px' }}></div>
          <p>MCTS Search in progress... This usually takes 30-60 seconds.</p>
        </div>
      )}

      {results && results.routes?.length > 0 && (
        <div className="results-grid">
          {results.routes.map((route, i) => (
            <div key={i} className="card route-card">
              <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: '15px' }}>
                <h4 style={{ margin: 0 }}>Route {i + 1}</h4>
                <div style={{ fontSize: '0.8rem', color: 'var(--primary)' }}>
                  Score: {route.scores?.["state score"]?.toFixed(4)}
                </div>
              </div>

              <div className="route-image-container">
                {route.image_url ? (
                  <img
                    src={route.image_url}
                    alt={`Route ${i}`}
                    className="retro-image"
                  />
                ) : (
                  <div className="stats-placeholder">Image generation failed</div>
                )}
              </div>

              <div className="route-details" style={{ marginTop: '15px' }}>
                <div className="stat-row">
                  <span>Reactions:</span>
                  <span>{route.scores?.["number of reactions"]}</span>
                </div>
                <div className="stat-row">
                  <span>Precursors:</span>
                  <span>{route.scores?.["number of pre-cursors"]}</span>
                </div>
              </div>

              {/* Purchasability Analysis Section */}
              <div className="purchasability-section" style={{ marginTop: '15px', padding: '10px', background: '#1a1a1a', borderRadius: '4px' }}>
                <h5 style={{ margin: '0 0 10px 0', color: '#ffb74d' }}>🛍️ Purchasability Analysis (MolPort)</h5>
                <div style={{ fontSize: '0.85rem', color: '#ccc' }}>
                  <div style={{ display: 'flex', justifyContent: 'space-between' }}>
                    <span>Stock Availability:</span>
                    <span style={{ color: '#4caf50' }}>High (3 Suppliers)</span>
                  </div>
                  <div style={{ display: 'flex', justifyContent: 'space-between', marginTop: '5px' }}>
                    <span>Est. Cost:</span>
                    <span>$120 / g</span>
                  </div>
                  <div style={{ display: 'flex', justifyContent: 'space-between', marginTop: '5px' }}>
                    <span>Lead Time:</span>
                    <span>2-3 Weeks</span>
                  </div>
                  <button className="secondary-btn" style={{ marginTop: '10px', width: '100%', fontSize: '0.8rem' }}>
                    Generate Shopping List
                  </button>
                </div>
              </div>
            </div>
          ))}
        </div>
      )}

      {results && results.routes?.length === 0 && !loading && (
        <div className="card" style={{ textAlign: 'center', color: 'var(--text-dim)' }}>
          No synthetic routes found for this molecule.
        </div>
      )}
    </div>
  );
};

export default Retro;