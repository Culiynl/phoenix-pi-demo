import React from 'react';
import { NavLink } from 'react-router-dom';

const Layout = ({ children }) => {
  return (
    <>
      <nav className="sidebar">
        <h3 style={{fontSize: '0.8rem', color: 'var(--text-dim)'}}>FLOATING NAVBAR</h3>
        <NavLink to="/projects" className="nav-item">Projects</NavLink>
        <NavLink to="/literature" className="nav-item">Literature Review</NavLink>
        <NavLink to="/docking" className="nav-item">Molecular Docking</NavLink>
        <NavLink to="/retro" className="nav-item">Retrosynthesis</NavLink>
        
        <div style={{marginTop: 'auto', textAlign: 'center'}}>
           <div className="user-avatar">👤 Guest</div>
        </div>
      </nav>
      <main className="main-content">{children}</main>
    </>
  );
};

export default Layout;