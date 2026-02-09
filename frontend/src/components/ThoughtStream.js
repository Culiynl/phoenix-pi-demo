import React, { useEffect, useState } from 'react';
import { db } from '../firebase';
import { doc, onSnapshot } from "firebase/firestore";
import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';

const ThoughtStream = ({ missionId }) => {
    const [logs, setLogs] = useState([
        { type: 'system', content: 'Connecting to PHOENIX-PI Mission Control...', timestamp: new Date() }
    ]);
    const [isThinking, setIsThinking] = useState(false);
    const scrollRef = React.useRef(null);

    useEffect(() => {
        if (scrollRef.current) {
            scrollRef.current.scrollTop = scrollRef.current.scrollHeight;
        }
    }, [logs]);

    useEffect(() => {
        if (!missionId) return;

        const missionRef = doc(db, "missions", missionId);

        // Real-time listener
        const unsubscribe = onSnapshot(missionRef, (docSnap) => {
            if (docSnap.exists()) {
                const data = docSnap.data();
                setIsThinking(data.is_thinking || false);
                if (data.thought_log) {
                    const formattedLogs = data.thought_log.map(log => ({
                        ...log,
                        timestamp: log.timestamp ? new Date(log.timestamp) : new Date()
                    }));
                    setLogs(formattedLogs);
                }
            } else {
                // Document hasn't been created yet (or backend failed to create it)
                setLogs([{ type: 'system', content: 'Waiting for Mission Initialization...', timestamp: new Date() }]);
            }
        }, (error) => {
            console.error("Firestore error:", error);
            setLogs(prev => [...prev, { type: 'error', content: `Connection Error: ${error.message}`, timestamp: new Date() }]);
        });

        return () => unsubscribe();
    }, [missionId]);

    const [input, setInput] = useState("");

    const handleSend = async () => {
        if (!input.trim()) return;
        const msg = input.trim();
        setInput(""); // Clear input immediately

        // Optimistic UI update handled by Firestore listener eventually, 
        // but we could add a temp local log if needed.

        try {
            await fetch(`http://localhost:8000/api/mission/${missionId}/chat`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ message: msg })
            });
        } catch (e) {
            console.error("Chat failed", e);
        }
    };

    return (
        <div className="thought-stream-panel">
            <div className="panel-header">
                <h3>THOUGHT STREAM</h3>
                <div style={{ display: 'flex', gap: '8px', alignItems: 'center' }}>
                    {isThinking && <div className="thinking-badge">THINKING</div>}
                    <div className="live-indicator">
                        <span className="blink">●</span> LIVE
                    </div>
                </div>
            </div>
            <div className="logs-container" ref={scrollRef}>
                {logs.map((log, i) => (
                    <div key={i} className={`log-entry ${log.type}`}>
                        <span className="timestamp">[{log.timestamp.toLocaleTimeString()}]</span>
                        <div className="content">
                            <ReactMarkdown
                                remarkPlugins={[remarkGfm]}
                                disallowedElements={['a']}
                                unwrapDisallowed={true}
                            >
                                {log.content}
                            </ReactMarkdown>
                        </div>
                    </div>
                ))}
            </div>
            <div className="chat-input-area">
                <input
                    className="chat-input"
                    type="text"
                    placeholder="Message PHOENIX..."
                    value={input}
                    onChange={(e) => setInput(e.target.value)}
                    onKeyDown={(e) => e.key === 'Enter' && handleSend()}
                />
            </div>
            <style>{`
        .thought-stream-panel {
            background: #111;
            border-right: 1px solid #333;
            flex: 1;
            min-height: 0;
            display: flex;
            flex-direction: column;
            color: #ccc;
            font-family: 'Consolas', 'Monaco', monospace;
        }
        .panel-header {
            padding: 10px;
            border-bottom: 1px solid #333;
            display: flex;
            justify-content: space-between;
            align-items: center;
        }
        .panel-header h3 { margin: 0; font-size: 14px; color: #888; }
        .live-indicator { font-size: 10px; color: #0f0; }
        .thinking-badge {
            background: #4285F4;
            color: white;
            font-size: 9px;
            padding: 2px 6px;
            border-radius: 10px;
            font-weight: bold;
            animation: pulse-animation 1.5s infinite;
        }
        @keyframes pulse-animation {
            0% { transform: scale(0.95); opacity: 0.7; }
            50% { transform: scale(1.1); opacity: 1; }
            100% { transform: scale(0.95); opacity: 0.7; }
        }
        .blink { animation: blinker 1s linear infinite; }
        @keyframes blinker { 50% { opacity: 0; } }
        
        .logs-container {
            flex-grow: 1;
            overflow-y: auto;
            padding: 10px;
            font-size: 12px;
            padding-bottom: 60px; /* Space for input */
        }
        .log-entry { margin-bottom: 6px; line-height: 1.4; }
        .timestamp { color: #666; margin-right: 8px; }
        .log-entry.thought { color: #aaa; }
        .log-entry.system { color: #4285F4; }
        .log-entry.action { color: #0f0; }
        .log-entry.error { color: #f44; }
        .log-entry.user { color: #fff; background: #222; padding: 4px 8px; border-radius: 4px; }
        
        .content p { margin: 0 0 8px 0; }
        .content p:last-child { margin-bottom: 0; }
        .content ul, .content ol { padding-left: 20px; margin: 8px 0; }
        .content code { background: #333; padding: 2px 4px; border-radius: 3px; font-size: 11px; }
        .content strong { color: #4285F4; }
        
        .chat-input-area {
            padding: 10px;
            border-top: 1px solid #333;
            background: #000;
        }
        .chat-input {
            width: 100%;
            background: #111;
            border: 1px solid #333;
            color: white;
            padding: 8px;
            border-radius: 4px;
            font-family: monospace;
        }
        .chat-input:focus { outline: none; border-color: #4285F4; }
      `}</style>
        </div>
    );
};
export default ThoughtStream;
