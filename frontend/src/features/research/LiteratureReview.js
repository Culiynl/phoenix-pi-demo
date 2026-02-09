import React, { useState } from 'react';

const LiteratureReview = () => {
    const [query, setQuery] = useState('');
    const [papers, setPapers] = useState([]);
    const [loading, setLoading] = useState(false);
    const [summary, setSummary] = useState('');
    const [synthesizing, setSynthesizing] = useState(false);

    const handleSearch = async () => {
        if (!query) return;
        setLoading(true);
        setSummary(''); // Clear previous summary
        try {
            const baseUrl = process.env.REACT_APP_API_URL || 'http://localhost:8000';
            const res = await fetch(`${baseUrl}/api/research/arxiv?q=${encodeURIComponent(query)}`);
            if (!res.ok) throw new Error("Search failed");
            const data = await res.json();
            setPapers(data.papers || []);
        } catch (e) {
            console.error(e);
            alert("Failed to fetch papers");
        } finally {
            setLoading(false);
        }
    };

    const handleSynthesize = async () => {
        if (papers.length === 0) return;
        setSynthesizing(true);
        try {
            const baseUrl = process.env.REACT_APP_API_URL || 'http://localhost:8000';
            const res = await fetch(`${baseUrl}/api/research/synthesize`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ papers })
            });
            const data = await res.json();
            setSummary(data.summary);
        } catch (e) {
            console.error(e);
            alert("Synthesis failed");
        } finally {
            setSynthesizing(false);
        }
    };

    return (
        <div className="card" style={{ height: '100%', display: 'flex', flexDirection: 'column' }}>
            <h2>LITERATURE REVIEW</h2>

            <div style={{ display: 'flex', gap: '10px', marginBottom: '20px' }}>
                <input
                    type="text"
                    value={query}
                    onChange={(e) => setQuery(e.target.value)}
                    placeholder="Search ArXiv (e.g. 'BACE1 inhibitors', 'Glioblastoma targets')"
                    style={{ flexGrow: 1, padding: '10px', borderRadius: '4px', border: '1px solid #333', background: '#222', color: '#fff' }}
                    onKeyDown={(e) => e.key === 'Enter' && handleSearch()}
                />
                <button className="primary-btn" onClick={handleSearch} disabled={loading}>
                    {loading ? "Searching..." : "Search"}
                </button>
            </div>

            <div style={{ display: 'flex', gap: '20px', flexGrow: 1, overflow: 'hidden' }}>
                {/* Results List */}
                <div style={{ flex: 1, overflowY: 'auto', paddingRight: '10px', borderRight: '1px solid #333' }}>
                    {papers.map((p, i) => (
                        <div key={i} style={{ background: '#1a1a1a', padding: '15px', marginBottom: '10px', borderRadius: '4px' }}>
                            <h4 style={{ margin: '0 0 5px 0', color: '#4285F4' }}>
                                <a href={p.pdf_url} target="_blank" rel="noreferrer" style={{ color: '#4285F4', textDecoration: 'none' }}>{p.title}</a>
                            </h4>
                            <p style={{ fontSize: '12px', color: '#888', margin: '0 0 10px 0' }}>{p.authors.join(', ')} • {p.published}</p>
                            <p style={{ fontSize: '13px', lineHeight: '1.4', color: '#ccc' }}>{p.summary.substring(0, 300)}...</p>
                        </div>
                    ))}
                    {papers.length === 0 && !loading && <p style={{ color: '#666', textAlign: 'center' }}>No papers found. Try a search.</p>}
                </div>

                {/* Synthesis Panel */}
                <div style={{ flex: 1, display: 'flex', flexDirection: 'column' }}>
                    <div style={{ marginBottom: '10px', display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
                        <h3 style={{ margin: 0 }}>Gemini Synthesis</h3>
                        {papers.length > 0 && (
                            <button className="success-btn" onClick={handleSynthesize} disabled={synthesizing}>
                                {synthesizing ? "Synthesizing..." : "Synthesize Findings"}
                            </button>
                        )}
                    </div>
                    <div style={{
                        flexGrow: 1,
                        background: '#111',
                        padding: '15px',
                        borderRadius: '4px',
                        overflowY: 'auto',
                        border: '1px solid #333',
                        fontFamily: 'monospace',
                        whiteSpace: 'pre-wrap',
                        color: '#d4d4d4'
                    }}>
                        {summary || (synthesizing ? "Gemini is reading the abstracts..." : "Select papers and click Synthesize to generate a strategic summary.")}
                    </div>
                </div>
            </div>
        </div>
    );
};

export default LiteratureReview;
