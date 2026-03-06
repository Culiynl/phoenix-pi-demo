import React, { useState, useEffect, useRef } from 'react';
import { api } from '../services/api';
import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';

const ChatInterface = ({ missionId, missionData, embedded = false }) => {
    const [input, setInput] = useState('');
    const [messages, setMessages] = useState([]);
    const [sending, setSending] = useState(false);
    const endRef = useRef(null);

    // Sync logs from missionData to messages
    useEffect(() => {
        const rawLog = missionData?.log || missionData?.thought_log;
        if (!rawLog) return;

        const mapped = (rawLog || []).map((entry, i) => ({
            id: i,
            role: entry.type === 'user' ? 'user' : 'assistant',
            content: entry.content,
            type: entry.type,
            timestamp: entry.timestamp
        }));

        setMessages(mapped);
    }, [missionData]);

    useEffect(() => {
        endRef.current?.scrollIntoView({ behavior: 'smooth' });
    }, [messages]);

    const handleSend = async () => {
        if (!input.trim() || sending) return;
        const msg = input;
        setInput('');
        setSending(true);
        try {
            await api.chat(missionId, msg);
        } catch (e) {
            console.error(e);
            alert("Failed to send message");
        } finally {
            setSending(false);
        }
    };

    return (
        <div className={`chat-interface ${embedded ? 'embedded' : ''}`}>
            {embedded && <h3>Message: Phoenix</h3>}
            <div className="chat-history">
                {messages.length === 0 && <div className="empty-chat">Start the conversation...</div>}

                {messages.map((msg) => (
                    <div key={msg.id} className={`chat-msg ${msg.role} ${msg.type}`}>
                        <div className="msg-bubble">
                            {msg.type === 'thought' && <span className="thought-label">Thinking...</span>}
                            {msg.type === 'action' && <span className="action-label">Executing Tool...</span>}
                            {msg.type === 'system' && <span className="system-label">System</span>}

                            <ReactMarkdown remarkPlugins={[remarkGfm]}>
                                {msg.content}
                            </ReactMarkdown>
                        </div>
                        <span className="msg-time">{new Date(msg.timestamp).toLocaleTimeString()}</span>
                    </div>
                ))}
                {missionData?.is_thinking && (
                    <div className="chat-msg assistant thought">
                        <div className="msg-bubble">
                            <span className="typing-dot">.</span><span className="typing-dot">.</span><span className="typing-dot">.</span>
                        </div>
                    </div>
                )}
                <div ref={endRef} />
            </div>

            <div className="chat-input-area">
                <input
                    type="text"
                    placeholder="Ask Gemini to research, analyze, or plan..."
                    value={input}
                    onChange={(e) => setInput(e.target.value)}
                    onKeyDown={(e) => e.key === 'Enter' && handleSend()}
                    disabled={sending}
                />
                <button onClick={handleSend} disabled={sending || !input.trim()}>
                    {sending ? '...' : 'Send'}
                </button>
            </div>

            <style>{`
                .chat-interface { display: flex; flex-direction: column; height: 100%; border: 1px solid #222; border-radius: 12px; background: #0d0d0d; overflow: hidden; }
                .chat-interface.embedded { height: 400px; margin-bottom: 20px; border: 1px solid #333; }
                .chat-interface.embedded h3 { padding: 15px 20px; margin: 0; border-bottom: 1px solid #222; font-size: 14px; color: #4285F4; background: #111; }

                .chat-history { flex-grow: 1; overflow-y: auto; padding: 20px; display: flex; flex-direction: column; gap: 15px; }
                .empty-chat { color: #444; text-align: center; margin-top: 50px; }
                
                .chat-msg { display: flex; flex-direction: column; max-width: 80%; }
                .chat-msg.user { align-self: flex-end; align-items: flex-end; }
                .chat-msg.assistant { align-self: flex-start; align-items: flex-start; }
                
                .msg-bubble { padding: 10px 15px; border-radius: 12px; font-size: 13px; line-height: 1.5; color: #ddd; word-break: break-word; }
                .chat-msg.user .msg-bubble { background: #4285F4; color: #fff; border-bottom-right-radius: 2px; }
                .chat-msg.assistant .msg-bubble { background: #1a1a1a; border: 1px solid #333; border-bottom-left-radius: 2px; }
                
                .chat-msg.thought .msg-bubble { background: #111; border: 1px dashed #333; color: #888; font-style: italic; }
                .chat-msg.action .msg-bubble { background: #1a0f0f; border: 1px solid #422; color: #f88; }
                .chat-msg.system .msg-bubble { background: #0f1a0f; border: 1px solid #242; color: #8f8; font-family: monospace; }
                
                .thought-label, .action-label, .system-label { display: block; font-size: 9px; font-weight: bold; margin-bottom: 4px; text-transform: uppercase; opacity: 0.7; }
                
                .msg-time { font-size: 9px; color: #444; margin-top: 4px; }
                
                .chat-input-area { padding: 15px; background: #0a0a0a; border-top: 1px solid #222; display: flex; gap: 10px; }
                .chat-input-area input { flex-grow: 1; background: #111; border: 1px solid #333; color: #fff; padding: 10px 15px; border-radius: 20px; outline: none; }
                .chat-input-area input:focus { border-color: #4285F4; }
                .chat-input-area button { background: #4285F4; color: #fff; border: none; padding: 0 20px; border-radius: 20px; font-weight: bold; cursor: pointer; }
                .chat-input-area button:disabled { background: #222; color: #666; cursor: not-allowed; }
                
                .typing-dot { animation: bg-blink 1.4s infinite both; }
                .typing-dot:nth-child(2) { animation-delay: 0.2s; }
                .typing-dot:nth-child(3) { animation-delay: 0.4s; }
                @keyframes bg-blink { 0% { opacity: 0.2; } 20% { opacity: 1; } 100% { opacity: 0.2; } }
            `}</style>
        </div>
    );
};

export default ChatInterface;
