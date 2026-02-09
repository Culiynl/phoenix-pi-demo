import React from 'react';
import ReactDOM from 'react-dom/client';
import App from './App';
// This connects your global styles to the app
import './styles/theme.css'; 

const root = ReactDOM.createRoot(document.getElementById('root'));
root.render(
  <React.StrictMode>
    <App />
  </React.StrictMode>
);