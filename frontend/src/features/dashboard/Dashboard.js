import React, { useEffect, useState } from 'react';
import { useParams, useNavigate } from 'react-router-dom';
import ThoughtStream from '../../components/ThoughtStream';
import LiveWorkspace from '../../components/LiveWorkspace';
import { db } from '../../firebase';
import { doc, onSnapshot } from "firebase/firestore";
import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';

const Dashboard = () => {
  const { id } = useParams();
  const navigate = useNavigate();
  const [activeView, setActiveView] = useState('3d');
  const [summary, setSummary] = React.useState('');
  const [awaitingConfirm, setAwaitingConfirm] = useState(false);
  const [currentStep, setCurrentStep] = useState(1);
  const [isResuming, setIsResuming] = useState(false);
  const [isThinking, setIsThinking] = useState(false);
  const [autoCorrect, setAutoCorrect] = useState(false);
  const [researchMode, setResearchMode] = useState('hybrid');

  React.useEffect(() => {
    if (!id) return;
    const missionRef = doc(db, "missions", id);
    const unsubscribe = onSnapshot(missionRef, (docSnap) => {
      if (docSnap.exists()) {
        const data = docSnap.data();
        setSummary(data.research_summary || '');
        setAwaitingConfirm(data.awaiting_confirmation || false);
        setCurrentStep(data.current_step || 1);
        setIsThinking(data.is_thinking || false);
        setAutoCorrect(data.auto_correct || false);
        setResearchMode(data.research_mode || 'hybrid');
      }
    });
    return () => unsubscribe();
  }, [id]);

  // Auto-trigger resume if awaiting confirmation and auto-confirm is ON
  useEffect(() => {
    if (awaitingConfirm && autoCorrect && !isResuming && !isThinking) {
      console.log("Auto-Confirm active: Triggering resume...");
      handleResume();
    }
  }, [awaitingConfirm, autoCorrect, isResuming, isThinking]);

  const handleToggleAutoCorrect = async () => {
    try {
      const baseUrl = process.env.REACT_APP_API_URL || 'http://localhost:8000';
      const formData = new FormData();
      formData.append('mission_id', id);
      formData.append('auto_correct', !autoCorrect);

      await fetch(`${baseUrl}/api/mission/settings`, {
        method: 'POST',
        body: formData
      });
    } catch (err) {
      console.error(err);
    }
  };

  const handleResume = async () => {
    setIsResuming(true);
    try {
      const baseUrl = process.env.REACT_APP_API_URL || 'http://localhost:8000';
      const formData = new FormData();
      formData.append('mission_id', id);

      const response = await fetch(`${baseUrl}/api/mission/resume`, {
        method: 'POST',
        body: formData
      });
      if (!response.ok) throw new Error("Failed to resume");
    } catch (err) {
      console.error(err);
      alert("Failed to proceed: " + err.message);
    } finally {
      setIsResuming(false);
    }
  };

  const steps = [
    "Literature & Abstracts",
    "Planning & ID",
    "Protein/Ligand Extraction",
    "Molecular Stats",
    "Starter Suggestion",
    "RL Loop Optimization",
    "Retrosynthesis"
  ];

  return (
    <div className="mission-control-container">
      {/* Top Bar */}
      <div className="mission-header-stacked">
        <div className="header-top">
          <div className="title-area">
            <h2>MISSION: <span className="highlight">{id}</span> <small style={{ fontSize: '10px', color: '#666', marginLeft: '10px', textTransform: 'uppercase' }}>Mode: {researchMode}</small></h2>
            {isThinking && <div className="thinking-status">
              <span className="dot pulse"></span> Gemini is Thinking...
            </div>}
          </div>
          <div className="controls">
            {awaitingConfirm && !autoCorrect && (
              <button
                className="resume-btn"
                onClick={handleResume}
                disabled={isResuming}
              >
                {isResuming ? 'Resuming...' : '✦ Confirm & Proceed'}
              </button>
            )}
            <button
              className={`success-btn ${autoCorrect ? 'on' : ''}`}
              onClick={handleToggleAutoCorrect}
            >
              Auto-Confirm: {autoCorrect ? 'ON' : 'OFF'}
            </button>
            <button className="danger-btn" onClick={() => navigate('/projects')}>Abort</button>
          </div>
        </div>

        <div className="progress-bar-container">
          {steps.map((s, i) => (
            <div key={i} className={`step ${i + 1 === currentStep ? 'active' : ''} ${i + 1 < currentStep ? 'complete' : ''}`}>
              <div className="step-num">{i + 1}</div>
              <div className="step-label">{s}</div>
            </div>
          ))}
        </div>
      </div>

      {/* Main Split Layout */}
      <div className="mission-body">
        <div className="left-panel">
          <ThoughtStream missionId={id} />
        </div>
        <div className="right-panel">
          <LiveWorkspace missionId={id} activeView={activeView} setActiveView={setActiveView} />
        </div>
      </div>

      <style>{`
        .mission-control-container {
            display: flex;
            flex-direction: column;
            height: calc(100vh - 64px); 
            background: #000;
            overflow: hidden;
            color: #fff;
        }
        .research-summary-panel {
            background: #1a1a1a;
            border-bottom: 2px solid #444;
            padding: 12px;
            height: 200px; /* Fixed height to avoid squashing chat */
            overflow-y: auto;
            font-size: 12px;
            color: #ccc;
            line-height: 1.4;
        }
        .left-panel {
            width: 350px;
            min-width: 350px;
            height: 100%;
            display: flex;
            flex-direction: column;
            border-right: 1px solid #333;
        }
        .left-panel > *:last-child {
            flex-grow: 1;
            overflow: hidden;
        }
        .panel-label {
            font-size: 10px;
            color: #4285F4;
            font-weight: bold;
            margin-bottom: 8px;
            letter-spacing: 1px;
        }
        .summary-content {
            white-space: pre-wrap;
        }
        .mission-header-stacked {
            background: #111;
            border-bottom: 2px solid #333;
            padding: 15px 25px;
            display: flex;
            flex-direction: column;
            gap: 15px;
        }
        .header-top {
            display: flex;
            justify-content: space-between;
            align-items: center;
        }
        .progress-bar-container {
            display: flex;
            gap: 20px;
            justify-content: space-between;
            background: #000;
            padding: 10px 20px;
            border-radius: 8px;
            border: 1px solid #222;
        }
        .controls {
            display: flex;
            gap: 10px;
            align-items: center;
        }
        .controls button {
            padding: 8px 15px;
            border-radius: 4px;
            font-weight: bold;
            font-size: 11px;
            cursor: pointer;
            border: none;
            transition: all 0.2s;
            text-transform: uppercase;
        }
        .resume-btn {
            background: #4285F4;
            color: #fff;
            box-shadow: 0 0 15px rgba(66, 133, 244, 0.4);
            animation: pulse-border 2s infinite;
        }
        @keyframes pulse-border {
            0% { box-shadow: 0 0 5px rgba(66, 133, 244, 0.4); }
            50% { box-shadow: 0 0 20px rgba(66, 133, 244, 0.8); }
            100% { box-shadow: 0 0 5px rgba(66, 133, 244, 0.4); }
        }
        .success-btn {
            background: #1a1a1a;
            color: #444;
            border: 1px solid #333 !important;
            transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
        }
        .success-btn:hover {
            border-color: #666 !important;
            color: #888;
        }
        .success-btn.on {
            background: rgba(76, 175, 80, 0.1);
            color: #4CAF50;
            border: 1px solid #4CAF50 !important;
            box-shadow: 0 0 10px rgba(76, 175, 80, 0.2);
        }
        .success-btn.on:hover {
            background: rgba(76, 175, 80, 0.2);
            box-shadow: 0 0 15px rgba(76, 175, 80, 0.3);
        }
        .danger-btn {
            background: #3d0f0f;
            color: #F44336;
        }
        .thinking-status {
            font-size: 11px;
            color: #4285F4;
            display: flex;
            align-items: center;
            gap: 6px;
            margin-top: 5px;
        }
        .dot {
            width: 8px;
            height: 8px;
            border-radius: 50%;
            background: #4285F4;
        }
        .pulse {
            animation: pulse-animation 1.5s infinite;
        }
        @keyframes pulse-animation {
            0% { transform: scale(0.9); opacity: 0.6; }
            50% { transform: scale(1.2); opacity: 1; }
            100% { transform: scale(0.9); opacity: 0.6; }
        }
        .title-area {
            display: flex;
            flex-direction: column;
        }
        .step {
            display: flex;
            align-items: center;
            gap: 10px;
            opacity: 0.2;
            transition: all 0.3s;
        }
        .step.active {
            opacity: 1;
            transform: translateY(-2px);
        }
        .step.complete {
            opacity: 0.7;
            color: #4CAF50;
        }
        .step-num {
            width: 22px;
            height: 22px;
            border-radius: 50%;
            background: #222;
            display: flex;
            align-items: center;
            justify-content: center;
            font-size: 10px;
            font-weight: bold;
            border: 1px solid #333;
        }
        .step.active .step-num {
            background: #4285F4;
            border-color: #4285F4;
            color: white;
        }
        .step.complete .step-num {
            background: #4CAF50;
            border-color: #4CAF50;
            color: white;
        }
        .step-label {
            font-size: 10px;
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 0.8px;
        }
        .mission-header-stacked h2 {
            font-size: 20px;
            color: #fff;
            margin: 0;
            letter-spacing: -0.5px;
        }
        .highlight { color: #4285F4; }
        .mission-body {
            display: flex;
            flex-grow: 1;
            overflow: hidden;
            height: 0; /* Important for flex children to scroll */
        }
        .right-panel {
            flex-grow: 1;
            height: 100%;
            display: flex;
            flex-direction: column;
        }
        .danger-btn {
            background: #d32f2f;
            color: white;
            border: none;
            padding: 5px 10px;
            border-radius: 4px;
            margin-right: 10px;
            cursor: pointer;
        }
        .success-btn {
            background: #2e7d32;
            color: white;
            border: none;
            padding: 5px 10px;
            border-radius: 4px;
            cursor: default;
        }
      `}</style>
    </div>
  );
};

export default Dashboard;