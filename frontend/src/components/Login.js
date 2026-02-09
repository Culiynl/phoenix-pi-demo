import React from 'react';
import { signInWithPopup } from "firebase/auth";
import { auth, googleProvider } from "../firebase";

const Login = ({ onLogin }) => {
  const handleLogin = async () => {
    try {
      const result = await signInWithPopup(auth, googleProvider);
      const user = result.user;
      console.log("Logged in user:", user);
      onLogin(user);
    } catch (error) {
      console.error("Login failed:", error);
      alert("Login failed: " + error.message);
    }
  };

  return (
    <div className="login-container">
      <div className="login-box">
        <h1>PHOENIX-PI</h1>
        <p>Autonomous Research Agent</p>
        <button className="login-btn" onClick={handleLogin}>
          Sign in with Google
        </button>
      </div>
      <style>{`
        .login-container {
          display: flex;
          justify-content: center;
          align-items: center;
          width: 100vw;
          height: 100vh;
          background: #0f0f13; /* Dark background */
          color: white;
          position: fixed;
          top: 0;
          left: 0;
        }
        .login-box {
          text-align: center;
          padding: 2rem;
          background: #1a1a20;
          border-radius: 12px;
          border: 1px solid #333;
        }
        .login-btn {
          background: #4285F4;
          color: white;
          border: none;
          padding: 10px 20px;
          border-radius: 4px;
          font-size: 16px;
          cursor: pointer;
          margin-top: 20px;
        }
        .login-btn:hover {
          background: #357ae8;
        }
      `}</style>
    </div>
  );
};

export default Login;
