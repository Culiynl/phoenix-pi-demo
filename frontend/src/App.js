import React from 'react';
import { BrowserRouter, Routes, Route, Navigate } from 'react-router-dom';
import { onAuthStateChanged } from "firebase/auth";
import { auth } from "./firebase";
import Login from "./components/Login";
import Layout from './components/Layout';
import Projects from './features/projects/Projects';
import ToolsPage from './features/tools/ToolsPage';
import Dashboard from './features/dashboard/Dashboard';
import LiteratureReview from './features/research/LiteratureReview';
import Retro from './features/retro/Retro';
import './styles/theme.css';
import Docking from './features/docking/Docking';
import FeudalOptimization from './features/feudal/FeudalOptimization';

function App() {
  const [user, setUser] = React.useState(null);
  const [loading, setLoading] = React.useState(true);

  React.useEffect(() => {
    const unsubscribe = onAuthStateChanged(auth, (currentUser) => {
      setUser(currentUser);
      setLoading(false);
    });
    return () => unsubscribe();
  }, []);

  if (loading) return <div className="loading-screen">Loading PHOENIX-PI...</div>;

  if (!user) {
    return <Login onLogin={setUser} />;
  }

  return (
    <BrowserRouter>
      <Layout user={user}>
        <Routes>
          <Route path="/" element={<Navigate to="/projects" replace />} /> {/* Added 'replace' */}
          <Route path="/projects" element={<Projects />} />
          <Route path="/tools" element={<ToolsPage />} />
          <Route path="/literature" element={<LiteratureReview />} />

          {/* Allow Docking to handle specific project IDs */}
          <Route path="/docking" element={<Docking />} />
          <Route path="/docking/:projectId" element={<Docking />} />

          <Route path="/dashboard/:id" element={<Dashboard />} />
          <Route path="/optimization/:projectId" element={<FeudalOptimization />} />
          <Route path="/retro" element={<Retro />} />
        </Routes>
      </Layout>
    </BrowserRouter>
  );
}
export default App;