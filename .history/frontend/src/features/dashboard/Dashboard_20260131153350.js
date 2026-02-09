import React from 'react';
import { PieChart, Pie, Cell, ResponsiveContainer } from 'recharts';

const Dashboard = () => {
  const retroData = [
    { name: '1-3 steps', value: 400 },
    { name: '4-6 steps', value: 300 },
  ];
  const COLORS = ['#3b82f6', '#1e40af'];

  return (
    <div className="dashboard">
      <h2>Project Dashboard</h2>

      <div className="dashboard-grid">
        {/* Box 1: Top Molecule */}
        <div className="card">
          <h3 style={{fontSize: '0.8rem', color: 'var(--text-dim)'}}>TOP MOLECULE</h3>
          <div className="stats-placeholder" style={{height: '180px'}}>
             [ Mol* Viewer ]
          </div>
        </div>

        {/* Box 2: Feudal Network Stats (NEW) */}
        <div className="card">
          <h3 style={{fontSize: '0.8rem', color: 'var(--text-dim)'}}>FEUDAL NETWORK STATS</h3>
          <div className="feudal-log">
            <div>&gd; Initializing Optimization...</div>
            <div>&ld; Epoch 42: Success Rate 84%</div>
            <div>&ld; Score: -9.2 kcal/mol</div>
            <div>&ld; Generating Fragments...</div>
          </div>
          <div style={{marginTop: '10px', fontSize: '0.8rem'}}>
            Success Rate: <span style={{color: '#4ade80'}}>Stable</span>
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