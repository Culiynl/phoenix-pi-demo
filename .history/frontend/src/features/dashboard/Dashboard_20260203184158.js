import React from 'react';
import { useParams, useNavigate } from 'react-router-dom'; // Added useNavigate
import { PieChart, Pie, Cell, ResponsiveContainer } from 'recharts';

const Dashboard = () => {
  const { id } = useParams(); // This catches the ID from the URL (/dashboard/123)
  const navigate = useNavigate();

  const retroData = [
    { name: '1-3 steps', value: 400 },
    { name: '4-6 steps', value: 300 },
  ];
  const COLORS = ['#3b82f6', '#1e40af'];

  return (
    <div className="dashboard">
      <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
        <h2>Project Dashboard: <span style={{color: 'var(--primary)'}}>{id}</span></h2>
      </div>

      <div className="dashboard-grid">
        {/* Box 1: Top Molecule */}
        <div className="card" style={{ display: 'flex', flexDirection: 'column' }}>
          <h3 style={{fontSize: '0.8rem', color: 'var(--text-dim)'}}>TOP MOLECULE</h3>
          
          <div className="stats-placeholder" style={{flexGrow: 1, display: 'flex', alignItems: 'center', justifyContent: 'center', minHeight: '150px', background: '#000', borderRadius: '4px'}}>
             <span style={{color: '#444'}}>No Structure Loaded</span>
          </div>

          {/* NEW: Navigation Button to Docking */}
          <button 
            className="primary-btn" 
            style={{marginTop: '15px', width: '100%', background: 'var(--primary)', color: 'white', border: 'none'}}
            onClick={() => navigate(`/docking/${id}`)}
          >
            Open Molecular Docking Tool →
          </button>
        </div>

          <button className="primary-btn" style={{ background: 'var(--primary)', color: 'white' }} 
                  onClick={() => navigate(`/optimization/${id}`)}>
            AI Optimization Loop →
          </button>
        {/* Box 2: Feudal Network Stats */}
        <div className="card">
          <h3 style={{fontSize: '0.8rem', color: 'var(--text-dim)'}}>FEUDAL NETWORK STATS</h3>
          <div className="feudal-log">
            <div>&gt; Initializing Optimization...</div>
            <div>&gt; Epoch 42: Success Rate 84%</div>
            <div>&gt; Score: -9.2 kcal/mol</div>
          </div>
        </div>

        {/* Box 3: Retro Statistics */}
        <div className="card">
          <h3 style={{fontSize: '0.8rem', color: 'var(--text-dim)'}}>RETRO STATISTICS</h3>
          <div style={{height: '180px'}}>
            <ResponsiveContainer width="100%" height="100%">
              <PieChart>
                <Pie data={retroData} innerRadius={50} outerRadius={70} dataKey="value">
                  {retroData.map((entry, index) => (
                    <Cell key={`cell-${index}`} fill={COLORS[index % COLORS.length]} />
                  ))}
                </Pie>
              </PieChart>
            </ResponsiveContainer>
          </div>
        </div>
      </div>

      {/* Problem Summary (Bottom Section) */}
      <div className="card">
        <h2>PROBLEM SUMMARY</h2>
        <div className="summary-container">
           <div className="prop-list">
              <strong>DESIRED PROPERTIES:</strong><br/>
              • 200 &lt; molwt &lt; 400 <br/>
              • 1 &lt; logp &lt; 3 <br/>
              • Target: BACE1
           </div>
           <div className="stats-placeholder">
              [ ChEMBL Distribution Graph ]
           </div>
        </div>
      </div>
    </div>
  );
};

export default Dashboard;