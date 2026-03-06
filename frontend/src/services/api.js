// const BASE_URL = process.env.REACT_APP_API_URL || 'http://localhost:8000/api';
const BASE_URL = 'http://localhost:8000/api';

export const api = {
  getProjects: async () => {
    const res = await fetch(`${BASE_URL}/projects`, {
      headers: { "x-ms-devtunnel-skip": "true" }
    });
    if (!res.ok) throw new Error("Failed to fetch projects");
    return res.json();
  },

  exportProject: (missionId) => {
    // Note: window.open opens a new tab. If Microsoft blocks this, 
    // the user will see the "Continue" screen in the new tab and can just click it.
    window.open(`${BASE_URL}/mission/export/${missionId}`, '_blank');
  },

  importProject: async (file) => {
    const formData = new FormData();
    formData.append('file', file);
    const res = await fetch(`${BASE_URL}/mission/import`, {
      method: 'POST',
      headers: { "x-ms-devtunnel-skip": "true" },
      body: formData
    });
    if (!res.ok) throw new Error("Import failed");
    return res.json();
  },

  getMission: async (missionId) => {
    const res = await fetch(`${BASE_URL}/mission/${missionId}`, {
      headers: { "x-ms-devtunnel-skip": "true" }
    });
    if (!res.ok) throw new Error("Failed to fetch mission");
    return res.json();
  },

  getMissionLog: async (missionId) => {
    const res = await fetch(`${BASE_URL}/mission/${missionId}/log`, {
      headers: { "x-ms-devtunnel-skip": "true" }
    });
    if (!res.ok) throw new Error("Failed to fetch logs");
    return res.json();
  },

  runRetrosynthesis: async (smiles) => {
    const res = await fetch(`${BASE_URL}/retrosynthesis`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        "x-ms-devtunnel-skip": "true"
      },
      body: JSON.stringify({ smiles })
    });
    if (!res.ok) throw new Error("Backend Error");
    return res.json();
  },

  resumeMission: async (missionId) => {
    const formData = new FormData();
    formData.append('mission_id', missionId);
    const res = await fetch(`${BASE_URL}/mission/resume`, {
      method: 'POST',
      headers: { "x-ms-devtunnel-skip": "true" },
      body: formData
    });
    if (!res.ok) throw new Error("Failed to resume");
    return res.json();
  },

  chat: async (missionId, message) => {
    const res = await fetch(`${BASE_URL}/mission/${missionId}/chat`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        "x-ms-devtunnel-skip": "true"
      },
      body: JSON.stringify({ message })
    });
    if (!res.ok) throw new Error("Chat failed");
    return res.json();
  },

  updateMissionSettings: async (missionId, autoCorrect) => {
    const formData = new FormData();
    formData.append('mission_id', missionId);
    formData.append('auto_correct', autoCorrect);
    const res = await fetch(`${BASE_URL}/mission/settings`, {
      method: 'POST',
      headers: { "x-ms-devtunnel-skip": "true" },
      body: formData
    });
    if (!res.ok) throw new Error("Failed to update settings");
    return res.json();
  },

  deleteProject: async (missionId) => {
    const res = await fetch(`${BASE_URL}/projects/${missionId}`, {
      method: 'DELETE',
      headers: { "x-ms-devtunnel-skip": "true" }
    });
    if (!res.ok) throw new Error("Delete failed");
    return res.json();
  }
};