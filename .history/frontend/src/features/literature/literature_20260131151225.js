const LiteratureReview = () => {
  return (
    <div style={{display: 'flex', gap: '20px'}}>
      <div style={{flex: 2}}>
        <div className="card">Paper 1: BACE1 Inhibitors...</div>
        <div className="card">Paper 2: Blood-Brain Barrier...</div>
      </div>
      
      <div style={{flex: 1}}>
         <div className="card" style={{height: '100%'}}>
            <h3>Problem Summary</h3>
            <div style={{background: '#222', height: '150px'}}>Images/Graphs</div>
         </div>
      </div>

      {/* Floating Chatbot PhoenixGPT */}
      <div className="chatbot-window">
        <h4>PhoenixGPT 🐼</h4>
        <div className="chat-history">...</div>
        <input placeholder="Ask about pathways..." />
      </div>
    </div>
  );
};