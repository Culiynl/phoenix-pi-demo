const BASE_URL = 'http://localhost:8000/api';

export const api = {
  // Mock projects for now - match your "Guest" user logic
  getProjects: async () => {
    return [
      { id: '1', name: "Alzheimer's BACE1", owner: "Adrian Pham", date: "2025-05-20" },
      { id: '2', name: "Glioblastoma EGFR", owner: "Adrian Pham", date: "2025-06-12" }
    ];
  },

  // Your working retrosynthesis logic
  runRetrosynthesis: async (smiles) => {
    const res = await fetch(`${BASE_URL}/retrosynthesis`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ smiles })
    });
    if (!res.ok) throw new Error("Backend Error");
    return res.json();
  }
};