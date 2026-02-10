import React from 'react';
import { NavLink } from 'react-router-dom';

const Layout = ({ children, user }) => {
  return (
    <div className="app-container">
      <nav className="sidebar">
        {/* Match the h3 styling from theme.css */}
        <h3>PHOENIX-PI</h3>

        <nav className="sidebar-nav" style={{ display: 'flex', flexDirection: 'column', gap: '4px' }}>
          {/* Use 'nav-item' to match your theme.css exactly */}
          <NavLink
            to="/projects"
            className={({ isActive }) => isActive ? 'nav-item active' : 'nav-item'}
          >
            <span>Projects</span>
          </NavLink>

          <NavLink
            to="/tools"
            className={({ isActive }) => isActive ? 'nav-item active' : 'nav-item'}
          >
            <span>Tools</span>
          </NavLink>

          <NavLink
            to="/literature"
            className={({ isActive }) => isActive ? 'nav-item active' : 'nav-item'}
          >
            <span>Literature</span>
          </NavLink>
        </nav>

        {/* User section styling */}
        <div className="user-section" style={{ marginTop: 'auto', paddingTop: '20px', borderTop: '1px solid var(--border)' }}>
          {user ? (
            <div style={{ display: 'flex', alignItems: 'center', gap: '12px' }}>
              <img
                src={user.photoURL || 'https://via.placeholder.com/32'}
                alt="Profile"
                style={{ width: '32px', height: '32px', borderRadius: '50%', border: '1px solid var(--border)' }}
              />
              <div style={{ display: 'flex', flexDirection: 'column', overflow: 'hidden' }}>
                <span style={{ fontSize: '0.8rem', color: 'white', whiteSpace: 'nowrap', overflow: 'hidden', textOverflow: 'ellipsis' }}>
                  {user.displayName || user.email}
                </span>
                <span style={{ fontSize: '0.65rem', color: '#4caf50' }}>● Online</span>
              </div>
            </div>
          ) : (
            <div className="nav-item">Guest User</div>
          )}
        </div>
      </nav>

      <main className="main-content">
        {children}
      </main>
    </div>
  );
};

export default Layout;