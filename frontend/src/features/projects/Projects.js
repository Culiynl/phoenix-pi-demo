import React, { useState, useEffect } from 'react';
import { useNavigate } from 'react-router-dom';
import { api } from '../../services/api';

const Projects = () => {
  const [projects, setProjects] = useState([]);
  const [newProjectName, setNewProjectName] = useState('');
  const [description, setDescription] = useState('');
  const [arxivQuery, setArxivQuery] = useState('');
  const [arxivResults, setArxivResults] = useState([]);
  const [isSearching, setIsSearching] = useState(false);
  const [researchMode, setResearchMode] = useState('hybrid'); // hybrid, arxiv, internal
  const navigate = useNavigate();

  useEffect(() => {
    api.getProjects().then(setProjects);
  }, []);

  const handleArXivSearch = async () => {
    if (!arxivQuery) return;
    setIsSearching(true);
    try {
      // const baseUrl = process.env.REACT_APP_API_URL || 'http://localhost:8000';
      const baseUrl = process.env.REACT_APP_API_URL || 'http://localhost:8000';
      const response = await fetch(`${baseUrl}/research/arxiv?q=${arxivQuery}`);
      const data = await response.json();
      setArxivResults(data.papers || []);
    } catch (err) {
      console.error(err);
    } finally {
      setIsSearching(false);
    }
  };

  const handleUsePaper = (paper) => {
    setNewProjectName(paper.title);
    setDescription(`Research Topic: ${paper.title}\n\nSummary: ${paper.summary}\n\nObjective: Explore potential drug targets and small molecule inhibitors related to this mechanism.`);
  };

  const handleStartMission = async () => {
    if (!newProjectName.trim() && !description.trim()) {
      alert("Please provide a Project Name OR a Guiding Problem.");
      return;
    }

    try {
      const payload = {
        target_name: newProjectName,
        description: description,
        research_mode: researchMode
      };

      // const baseUrl = process.env.REACT_APP_API_URL || 'http://localhost:8000';
      const baseUrl = process.env.REACT_APP_API_URL || 'http://localhost:8000';
      const response = await fetch(`${baseUrl}/mission/start`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(payload)
      });

      if (!response.ok) {
        throw new Error(`Server error: ${response.statusText}`);
      }

      const data = await response.json();
      const missionId = data.mission_id;
      navigate(`/dashboard/${missionId}`);

    } catch (error) {
      console.error("Failed to start mission:", error);
      alert("Failed to start mission. Ensure backend is running.");
    }
  };

  const handleDeleteProject = async (e, id) => {
    e.stopPropagation();
    if (!window.confirm("Are you sure you want to delete this project? This action cannot be undone.")) return;
    try {
      await api.deleteProject(id);
      setProjects(prev => prev.filter(p => p.id !== id));
    } catch (err) {
      alert("Failed to delete project: " + err.message);
    }
  };

  return (
    <div className="card projects-container">
      <div className="mission-creator">
        <div className="creation-panel">
          <h3 style={{ margin: '0 0 10px 0', color: '#fff' }}>1. Start a New Autonomous Mission</h3>
          <textarea
            placeholder="Describe your guiding problem..."
            value={description}
            onChange={(e) => setDescription(e.target.value)}
            className="problem-input"
          />

          <div className="mode-selection">
            <span className="mode-label">Research Mode:</span>
            <div className="mode-options">
              <button
                className={`mode-btn ${researchMode === 'hybrid' ? 'active' : ''}`}
                onClick={() => setResearchMode('hybrid')}
                title="Ground in papers, fallback to internal knowlege if needed"
              >
                Hybrid
              </button>
              <button
                className={`mode-btn ${researchMode === 'arxiv' ? 'active' : ''}`}
                onClick={() => setResearchMode('arxiv')}
                title="Strictly use ArXiv papers for grounding"
              >
                Grounded
              </button>
              <button
                className={`mode-btn ${researchMode === 'internal' ? 'active' : ''}`}
                onClick={() => setResearchMode('internal')}
                title="Bypass ArXiv, use Gemini's internal knowledge immediately"
              >
                Fast Start
              </button>
            </div>
          </div>

          <div className="start-controls">
            <input
              type="text"
              placeholder="Project Name (Optional)"
              value={newProjectName}
              onChange={(e) => setNewProjectName(e.target.value)}
              className="name-input"
            />
            <button className="primary-btn launch-btn" onClick={handleStartMission}>
              ✦ Launch Mission
            </button>
          </div>
        </div>

        <div className="literature-panel">
          <h3 style={{ margin: '0 0 10px 0', color: '#fff' }}>2. Literature Discovery <small style={{ color: '#666' }}>(Optional)</small></h3>
          <div className="search-box">
            <input
              type="text"
              placeholder="Search ArXiv..."
              value={arxivQuery}
              onChange={e => setArxivQuery(e.target.value)}
              onKeyPress={e => e.key === 'Enter' && handleArXivSearch()}
            />
            <button onClick={handleArXivSearch} disabled={isSearching}>
              {isSearching ? '...' : 'Search'}
            </button>
          </div>
          <div className="arxiv-results">
            {arxivResults.map((p, i) => (
              <div key={i} className="paper-result">
                <div className="paper-title">{p.title}</div>
                <button className="use-paper-btn" onClick={() => handleUsePaper(p)}>Use as Problem</button>
              </div>
            ))}
          </div>
        </div>
      </div>

      <table className="mission-table">
        <thead>
          <tr>
            <th>Name</th>
            <th>Owner</th>
            <th>Created Date</th>
            <th style={{ width: '50px' }}></th>
          </tr>
        </thead>
        <tbody>
          {projects.map(p => (
            <tr key={p.id} onClick={() => navigate(`/dashboard/${p.id}`)}>
              <td style={{ color: 'var(--primary)', fontWeight: '600' }} title={p.name}>
                {p.name && p.name.length > 30 ? p.name.substring(0, 30) + "..." : p.name}
              </td>
              <td>{p.owner}</td>
              <td>{p.date}</td>
              <td>
                <button
                  className="delete-btn"
                  onClick={(e) => handleDeleteProject(e, p.id)}
                  title="Delete Project"
                >
                  <i className="trash-icon">🗑️</i>
                </button>
              </td>
            </tr>
          ))}
        </tbody>
      </table>

      <style>{`
        .projects-container { padding: 30px; }
        .mission-creator {
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 30px;
            margin-bottom: 40px;
        }
        .creation-panel, .literature-panel {
            background: #111;
            border: 1px solid #333;
            padding: 20px;
            border-radius: 12px;
            display: flex;
            flex-direction: column;
        }
        .problem-input {
            flex-grow: 1;
            background: #000;
            border: 1px solid #333;
            color: #fff;
            padding: 15px;
            border-radius: 8px;
            min-height: 120px;
            margin-bottom: 15px;
            font-size: 14px;
            resize: none;
        }
        .start-controls {
            display: flex;
            gap: 10px;
        }
        .mode-selection {
            margin-bottom: 20px;
            display: flex;
            align-items: center;
            gap: 15px;
        }
        .mode-label {
            font-size: 11px;
            color: #666;
            text-transform: uppercase;
            letter-spacing: 1px;
            font-weight: bold;
        }
        .mode-options {
            display: flex;
            gap: 4px;
            background: #000;
            padding: 3px;
            border-radius: 6px;
            border: 1px solid #222;
        }
        .mode-btn {
            background: transparent;
            border: none;
            color: #555;
            font-size: 10px;
            font-weight: bold;
            padding: 6px 12px;
            border-radius: 4px;
            cursor: pointer;
            transition: all 0.2s;
        }
        .mode-btn.active {
            background: rgba(153, 27, 27, 0.2);
            color: #f87171;
            border: 1px solid #991b1b;
            box-shadow: 0 0 8px rgba(153, 27, 27, 0.3);
        }
        .mode-btn:hover:not(.active) {
            color: #ccc;
            background: #111;
        }
        .name-input {
            flex-grow: 1;
            background: #000;
            border: 1px solid #333;
            color: #fff;
            padding: 10px;
            border-radius: 6px;
        }
        .launch-btn {
            background: linear-gradient(135deg, #7f1d1d 0%, #991b1b 100%);
            color: #fff;
            padding: 0 24px;
            min-height: 42px;
            border-radius: 8px;
            box-shadow: 0 0 18px rgba(153, 27, 27, 0.45), 0 2px 8px rgba(0,0,0,0.4);
            border: 1px solid #7f1d1d;
            cursor: pointer;
            font-weight: 800;
            font-size: 12px;
            letter-spacing: 1px;
            text-transform: uppercase;
            transition: all 0.25s cubic-bezier(0.4, 0, 0.2, 1);
            white-space: nowrap;
        }
        .launch-btn:hover {
            background: linear-gradient(135deg, #991b1b 0%, #b91c1c 100%);
            box-shadow: 0 0 28px rgba(153, 27, 27, 0.65), 0 4px 16px rgba(0,0,0,0.5);
            transform: translateY(-2px);
        }
        .launch-btn:active {
            transform: translateY(0px);
        }
        .search-box {
            display: flex;
            gap: 5px;
            margin-bottom: 15px;
        }
        .search-box input {
            flex-grow: 1;
            background: #000;
            border: 1px solid #333;
            color: #fff;
            padding: 8px;
            border-radius: 4px;
        }
        .search-box button {
            background: #222;
            border: 1px solid #333;
            color: #fff;
            padding: 0 10px;
            cursor: pointer;
        }
        .arxiv-results {
            max-height: 200px;
            overflow-y: auto;
            background: #050505;
            border-radius: 4px;
            padding: 10px;
        }
        .paper-result {
            padding: 10px;
            border-bottom: 1px solid #111;
            display: flex;
            justify-content: space-between;
            align-items: center;
            gap: 10px;
        }
        .paper-title { font-size: 12px; color: #ccc; flex-grow: 1; }
        .use-paper-btn {
            background: #0f3d0f;
            color: #4CAF50;
            border: 1px solid #4CAF50;
            font-size: 10px;
            padding: 4px 8px;
            border-radius: 4px;
            cursor: pointer;
            white-space: nowrap;
        }
        .mission-table tr { cursor: pointer; transition: background 0.2s; }
        .mission-table tr:hover { background: #111; }
        
        .delete-btn {
          background: transparent;
          border: none;
          color: #ff4444;
          cursor: pointer;
          opacity: 0.5;
          transition: opacity 0.2s;
          font-size: 16px;
          padding: 5px;
        }
        .delete-btn:hover {
          opacity: 1;
        }
        .trash-icon { font-style: normal; }
      `}</style>
    </div>
  );
};

export default Projects;