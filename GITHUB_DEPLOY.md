# Alternative Deployment: GitHub + Railway (RECOMMENDED)

## ❌ Problem: Railway CLI Upload Fails

The `railway up` command times out or fails with large projects. This is a known Railway CLI issue.

## ✅ Solution: Deploy via GitHub (Much More Reliable)

Railway works best when deploying from GitHub. Here's how:

---

## Step 1: Create GitHub Repository

```bash
# Navigate to project root
cd "f:\Coding Stuff\Full Model 2026 Fixed\webserver_stuff"

# Initialize git (if not already)
git init

# Add all files (gitignore will exclude large files)
git add .

# Commit
git commit -m "Initial commit for Railway deployment"

# Create GitHub repo and push
gh repo create phoenix-pi-demo --public --source=. --push
```

**Don't have GitHub CLI?** Install it:
```bash
winget install GitHub.cli
```

Or create repo manually:
1. Go to https://github.com/new
2. Create repo "phoenix-pi-demo"
3. Push your code:
   ```bash
   git remote add origin https://github.com/YOUR_USERNAME/phoenix-pi-demo.git
   git push -u origin main
   ```

---

## Step 2: Deploy Backend from GitHub

1. Go to https://railway.app/dashboard
2. Click on your "kind-success" project
3. Click "New Service" → "GitHub Repo"
4. Select `phoenix-pi-demo`
5. Set **Root Directory**: `backend`
6. Railway auto-detects Python and deploys!

---

## Step 3: Add Environment Variables

In Railway dashboard → Backend service → Variables:

```
GOOGLE_API_KEY=your_actual_key_here
GEMINI_MODEL=gemini-2.0-flash
```

---

## Step 4: Configure Start Command

Railway dashboard → Backend service → Settings → Start Command:

```
uvicorn api:app --host 0.0.0.0 --port $PORT
```

---

## Step 5: Deploy Frontend

1. In Railway dashboard, click "New Service" → "GitHub Repo"
2. Select same `phoenix-pi-demo` repo
3. Set **Root Directory**: `frontend`
4. Railway auto-detects React and builds!

---

## Step 6: Update Frontend API URL

After backend deploys, you'll get a URL like:
`https://kind-success-production-xxxx.up.railway.app`

Update `frontend/src/firebase.js` (or wherever you define API_URL):

```javascript
const API_URL = "https://YOUR_BACKEND_URL_HERE";
```

Commit and push:
```bash
git add .
git commit -m "Update API URL"
git push
```

Railway auto-redeploys on every push!

---

## Why This Works Better

✅ **No upload limits** - Railway pulls from GitHub
✅ **Auto-deploys** - Push to GitHub = automatic deployment
✅ **Better logs** - See build progress in dashboard
✅ **Rollbacks** - Easy to revert to previous commits
✅ **Team collaboration** - Others can contribute

---

## Quick Commands

```bash
# One-time setup
cd "f:\Coding Stuff\Full Model 2026 Fixed\webserver_stuff"
git init
git add .
git commit -m "Initial commit"
gh repo create phoenix-pi-demo --public --source=. --push

# Then use Railway dashboard to connect GitHub repo
# No more CLI needed!
```

---

## Alternative: Render.com (Even Easier)

If Railway continues to have issues, try Render:

1. Go to https://render.com
2. Click "New +" → "Web Service"
3. Connect GitHub repo
4. Set:
   - **Root Directory**: `backend`
   - **Build Command**: `pip install -r requirements-light.txt`
   - **Start Command**: `uvicorn api:app --host 0.0.0.0 --port $PORT`
5. Add environment variables
6. Deploy!

Render's free tier is very generous and more reliable than Railway CLI.
