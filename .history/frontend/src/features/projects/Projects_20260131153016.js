import React, { useState, useEffect } from 'react';
import { useNavigate } from 'react-router-dom';
import { api } from '../../services/api';

const Projects = () => {
  const [projects, setProjects] = useState([]);
  const navigate = useNavigate();

  useEffect(() => {
    api.getProjects().then(setProjects);
  }, []);

  return (
    <div className="card">
      <div style={{display:'flex', justifyContent:'space-between', alignItems:'center', marginBottom: '20px'}}>
        <h2>PROJECT LIST</h2>
        <button className="primary-btn">+ New Project</button>
      </div>
      <table>
        <thead>
          <tr><th>Name</th><th>Owner</th><th>Created Date</th></tr>
        </thead>
        <tbody>
          {projects.map(p => (
            <tr key={p.id} onClick={() => navigate(`/dashboard/${p.id}`)} style={{cursor:'pointer'}}>
              <td style={{color: 'var(--primary)', fontWeight:'600'}}>{p.name}</td>
              <td>{p.owner}</td>
              <td>{p.date}</td>
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
};

export default Projects;