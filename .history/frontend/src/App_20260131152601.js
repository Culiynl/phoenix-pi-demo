import React from 'react';
import { BrowserRouter, Routes, Route, Navigate } from 'react-router-dom';
import Layout from './components/Layout';
import Projects from './features/projects/Projects';
import Dashboard from './features/dashboard/Dashboard';
import LiteratureReview from './features/lit-review/lit';
import Retro from './features/retro/Retro';
import './styles/theme.css';

function App() {
  return (
    <BrowserRouter>
      <Layout>
        <Routes>
          <Route path="/" element={<Navigate to="/projects" />} />
          <Route path="/projects" element={<Projects />} />
          <Route path="/dashboard/:id" element={<Dashboard />} />
          <Route path="/literature" element={<Literature />} />
          <Route path="/retro" element={<Retro />} />
          {/* Add Molecular Docking route here later */}
        </Routes>
      </Layout>
    </BrowserRouter>
  );
}

export default App;