import React from 'react';
import { NavLink } from 'react-router-dom';

const Layout = ({ children, user }) => {
  return (
    <div className="layout-container">
      <nav className="sidebar">
        <h3 style={{ fontSize: '0.9rem', color: '#888', marginBottom: '20px', letterSpacing: '2px' }}>PHOENIX</h3>

        <nav className="sidebar-nav">
          <NavLink to="/projects" className={({ isActive }) => isActive ? 'nav-link active' : 'nav-link'}>
            <span className="nav-icon">📊</span>
            <span>Projects</span>
          </NavLink>
          <NavLink to="/tools" className={({ isActive }) => isActive ? 'nav-link active' : 'nav-link'}>
            <span className="nav-icon">🔬</span>
            <span>Tools</span>
          </NavLink>
          <NavLink to="/literature" className={({ isActive }) => isActive ? 'nav-link active' : 'nav-link'}>
            <span className="nav-icon">📚</span>
            <span>Literature</span>
          </NavLink>
        </nav>

        <div className="user-section">
          {user ? (
            <div className="user-profile">
              <img src={user.photoURL} alt="Profile" className="avatar-img" />
              <div className="user-info">
                <span className="user-name">{user.displayName}</span>
                <span className="user-status">Online</span>
              </div>
            </div>
          ) : (
            <div className="user-avatar">Guest</div>
          )}
        </div>
      </nav>
      <main className="main-content">{children}</main>

      <style>{`
        .layout-container { display: flex; height: 100vh; width: 100vw; overflow: hidden; }
        .sidebar { width: 240px; background: #0f0f13; border-right: 1px solid #222; padding: 20px; display: flex; flex-direction: column; }
        .nav-links { display: flex; flex-direction: column; gap: 10px; }
        .nav-item { color: #888; text-decoration: none; padding: 10px 15px; border-radius: 8px; font-size: 0.9rem; transition: all 0.2s; }
        .nav-item:hover { background: #1a1a20; color: white; }
        .nav-item.active { background: #2a2a35; color: white; border-left: 3px solid #4285F4; }
        .main-content { flex: 1; overflow: auto; background: #050505; }
        
        .user-section { margin-top: auto; padding-top: 20px; border-top: 1px solid #222; }
        .user-profile { display: flex; align-items: center; gap: 10px; }
        .avatar-img { width: 36px; height: 36px; border-radius: 50%; object-fit: cover; border: 2px solid #333; }
        .user-info { display: flex; flex-direction: column; }
        .user-name { font-size: 0.85rem; color: white; font-weight: 500; }
        .user-status { font-size: 0.7rem; color: #4caf50; }
        .user-avatar { padding: 10px; background: #222; border-radius: 50%; width: 40px; height: 40px; text-align: center; line-height: 20px; color: #888; }
      `}</style>
    </div>
  );
};

export default Layout;