import React, { useState, useEffect } from 'react';
import { useNavigate } from 'react-router-dom';
import { api } from '../../services/api';

const Projects = () => {
  return (
    <div className="card">
      <div style={{display:'flex', justifyContent:'space-between', alignItems:'center'}}>
        <h2>PROJECT LIST</h2>
        <button className="primary-btn">+ New Project</button>
      </div>
      
      <table>
        <thead>
          <tr>
            <th>Name</th>
            <th>Owner</th>
            <th>Creation date</th>
          </tr>
        </thead>
        <tbody>
          {/* Apply the "magenta" class to the first row to match your wireframe */}
          <tr className="project-row-highlighted">
            <td>rwar</td>
            <td>Adrian Pham</td>
            <td>2025-05-20</td>
          </tr>
          <tr>
            <td>Alzheimer's BACE1</td>
            <td>Adrian Pham</td>
            <td>2025-05-20</td>
          </tr>
        </tbody>
      </table>
    </div>
  );
};

export default Projects;