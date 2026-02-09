# Quick Deployment Guide - Vertical Slice Demo

## 🎯 Goal: Minimal Working Demo Online

This is a **vertical slice** - just enough to demonstrate Gemini 3.0's capabilities without deploying everything.

---

## Option 1: Railway (Recommended - Easiest)

**Why**: Free tier, auto-deploys from GitHub, handles both frontend & backend

### Steps:

1. **Create GitHub Repo**
```bash
cd "f:\Coding Stuff\Full Model 2026 Fixed\webserver_stuff"
git init
git add .
git commit -m "Initial commit"
gh repo create phoenix-pi-demo --public --source=. --remote=origin --push
```

2. **Deploy to Railway**
- Go to https://railway.app
- Click "New Project" → "Deploy from GitHub"
- Select your `phoenix-pi-demo` repo
- Railway auto-detects FastAPI and React

3. **Add Environment Variables**
```
GOOGLE_API_KEY=your_key_here
FIREBASE_CREDENTIALS=your_firebase_json
```

4. **Done!** Railway gives you URLs:
- Backend: `https://phoenix-backend.railway.app`
- Frontend: `https://phoenix-pi.railway.app`

**Cost**: $0 (free tier: 500 hours/month)

---

## Option 2: Vercel (Frontend) + Render (Backend)

### Frontend (Vercel)
```bash
cd frontend
npm run build
vercel --prod
```

### Backend (Render)
1. Go to https://render.com
2. New Web Service → Connect GitHub repo
3. Build: `pip install -r requirements-light.txt`
4. Start: `uvicorn api:app --host 0.0.0.0 --port $PORT`

**Cost**: $0 (both have free tiers)

---

## Option 3: Google Cloud Run (What You Tried Before)

**Problem**: Heavy models (AiZynthFinder) break Cloud Run
**Solution**: Deploy lightweight version only

### Lightweight Backend (No Heavy Models)
```bash
cd backend

# Build lightweight container
gcloud builds submit --tag gcr.io/YOUR_PROJECT/phoenix-light -f Dockerfile.light

# Deploy
gcloud run deploy phoenix-backend \
    --image gcr.io/YOUR_PROJECT/phoenix-light \
    --platform managed \
    --region us-central1 \
    --allow-unauthenticated \
    --set-env-vars GOOGLE_API_KEY=your_key
```

### Frontend (Firebase Hosting)
```bash
cd frontend
npm run build

# Deploy to Firebase
firebase deploy --only hosting
```

**Cost**: ~$0 (Cloud Run free tier: 2M requests/month)

---

## 🎬 Vertical Slice Features to Deploy

### Must-Have (3 Gemini Requirements):
1. ✅ **Multimodal**: Upload protein image → Gemini analyzes
2. ✅ **Long Context**: Analyze 20+ papers simultaneously
3. ✅ **Agentic**: Autonomous 7-step workflow with function calling

### Can Skip for Demo:
- ❌ AiZynthFinder (too heavy)
- ❌ P2Rank docking (requires GPU)
- ❌ Ligand generation (large models)

### Lightweight Alternatives:
- ✅ ChEMBL queries (API-based, fast)
- ✅ Literature review (Gemini-powered)
- ✅ Molecular property analysis (RDKit, lightweight)

---

## 📝 Deployment Checklist

- [ ] Create GitHub repo
- [ ] Choose deployment platform (Railway recommended)
- [ ] Add environment variables
- [ ] Test deployed backend API
- [ ] Update frontend to use deployed backend URL
- [ ] Deploy frontend
- [ ] Test end-to-end workflow
- [ ] Record demo video

---

## 🚀 Recommended: Railway Deployment

**Single Command Deployment:**
```bash
# Install Railway CLI
npm i -g @railway/cli

# Login
railway login

# Deploy backend
cd backend
railway up

# Deploy frontend
cd ../frontend
railway up
```

Railway automatically:
- Detects Python/Node.js
- Installs dependencies
- Exposes ports
- Provides HTTPS URLs

**Total time**: ~5 minutes

---

## 🎥 Demo Video Script (2 minutes)

### 0:00-0:30 - Problem
"Traditional drug discovery tools are fragmented and can't adapt."

### 0:30-1:00 - Gemini Multimodal
Upload protein image → Gemini identifies binding sites

### 1:00-1:30 - Gemini Long Context
Load 20 papers → Gemini finds cross-paper patterns

### 1:30-2:00 - Gemini Agentic
Show autonomous workflow → Gemini pivots when ChEMBL fails

**Tools needed**: OBS Studio (free screen recorder)
