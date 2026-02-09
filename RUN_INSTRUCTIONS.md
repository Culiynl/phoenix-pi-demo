# How to Run PHOENIX-PI

### 1. Prerequisites
- Python 3.9+
- Node.js & npm
- Firebase Project (Auth and Firestore enabled)
- Gemini API Key

### 2. Backend Setup
1. Open a terminal in `backend/`.
2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```
3. Ensure `serviceAccountKey.json` is in the `backend/` root.
4. Start the server:
   ```bash
   uvicorn api:app --reload --port 8000
   ```

### 3. Frontend Setup
1. Open a terminal in `frontend/`.
2. Install dependencies:
   ```bash
   npm install
   ```
3. Start the development server:
   ```bash
   npm start
   ```

### 4. Running a Mission
1. Go to `http://localhost:3000`.
2. Log in with Google.
3. In the sidebar, enter a research target (e.g., "Alzheimer’s BACE1").
4. Watch the **Mission Control** dashboard. The agent will begin its autonomous cycle, and you'll see its thoughts appearing in real-time in the left panel.
