import { PieChart, Pie, Cell, ResponsiveContainer } from 'recharts';

const Dashboard = ({ projectName = "BACE 1" }) => {
  const data = [
    { name: 'Step 1', value: 400 },
    { name: 'Step 2', value: 300 },
  ];

  return (
    <div className="dashboard">
      <header>
        <h1>{projectName} Dashboard</h1>
      </header>

      <div className="grid-dashboard">
        {/* Top Molecule - Placeholder for Mol* */}
        <div className="card">
          <h3>Top Molecule</h3>
          <div style={{height: '200px', background: '#000', borderRadius: '8px', display: 'flex', alignItems: 'center', justifyContent: 'center'}}>
            [ Mol* Viewer Placeholder ]
          </div>
        </div>

        {/* Stats - React-Charts */}
        <div className="card">
          <h3>Retro Statistics</h3>
          <ResponsiveContainer width="100%" height={200}>
            <PieChart>
              <Pie data={data} innerRadius={60} outerRadius={80} fill="#3b82f6" dataKey="value" />
            </PieChart>
          </ResponsiveContainer>
        </div>
      </div>

      <div className="card" style={{marginTop: '20px'}}>
        <h3>Problem Summary: [Disease + Target]</h3>
        <p>Targeting BACE1 for Alzheimer's Disease</p>
        <div style={{display: 'flex', gap: '10px'}}>
            <span className="tag tag-high">Desired LogP: 1-3</span>
            <span className="tag tag-high">MW %lt; 450</span>
        </div>
      </div>
    </div>
  );
};