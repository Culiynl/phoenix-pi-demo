import React from 'react';
import { BrowserRouter, Routes, Route, Navigate } from 'react-router-dom';
import Layout from './components/Layout';
import Projects from './features/projects/Projects';
import Dashboard from './features/dashboard/Dashboard';
import LiteratureReview from './features/lit-review/lit';
import Retro from './features/retro/Retro';
import './styles/theme.css';
import Docking from './features/docking/Docking';

function App() {
  return (
    <BrowserRouter>
      <Layout>
        <Routes>
          <Route path="/" element={<Navigate to="/projects" />} />
          <Route path="/projects" element={<Projects />} />
          
          {/* Allow Docking to handle specific project IDs */}
          <Route path="/docking" element={<Docking />} />
          <Route path="/docking/:projectId" element={<Docking />} />
          
          <Route path="/dashboard/:id" element={<Dashboard />} />
          <Route path="/literature" element={<LiteratureReview />} />
          <Route path="/retro" element={<Retro />} />
        </Routes>
      </Layout>
    </BrowserRouter>
  );
}
export default App;