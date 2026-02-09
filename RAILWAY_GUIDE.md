# Railway Deployment Guide - Step by Step

## ✅ You've Already Done: `railway link`

Good! You're linked to project "kind-success". Now let's deploy.

---

## Step 1: Create a Service (Choose "Empty Service")

When Railway asks "What would you like to create?", select:
- **Empty Service** (not GitHub repo)

This creates a service that you'll deploy to directly from your local files.

---

## Step 2: Deploy the Backend

You're already in the backend directory. Now run:

```bash
railway up
```

This will:
- Upload your backend code
- Auto-detect it's a Python app
- Install dependencies from `requirements.txt`
- Start the service

**Note**: Railway will use `requirements.txt` by default. If you want the lightweight version for deployment, you can temporarily rename:
```bash
# Optional: Use lightweight deps for faster deployment
mv requirements.txt requirements-full.txt
mv requirements-light.txt requirements.txt
railway up
# Then switch back
mv requirements.txt requirements-light.txt
mv requirements-full.txt requirements.txt
```

---

## Step 3: Set Environment Variables

After deployment, you need to add your API keys:

```bash
# Set Google API key
railway variables set GOOGLE_API_KEY=your_actual_key_here

# Set Firebase credentials (if needed)
railway variables set FIREBASE_CREDENTIALS='{"type":"service_account",...}'
```

Or use the Railway dashboard:
1. Go to https://railway.app/dashboard
2. Click on your "kind-success" project
3. Click on the backend service
4. Go to "Variables" tab
5. Add:
   - `GOOGLE_API_KEY`
   - `FIREBASE_CREDENTIALS` (paste your Firebase JSON)

---

## Step 4: Configure Start Command

Railway needs to know how to start your app. Create a `Procfile`:

```bash
echo "web: uvicorn api:app --host 0.0.0.0 --port $PORT" > Procfile
```

Or set it in Railway dashboard → Settings → Start Command:
```
uvicorn api:app --host 0.0.0.0 --port $PORT
```

---

## Step 5: Deploy Frontend

```bash
cd ../frontend
railway up
```

Railway will auto-detect React and:
- Run `npm install`
- Run `npm run build`
- Serve the static files

---

## Step 6: Update Frontend to Use Backend URL

After backend deploys, Railway gives you a URL like:
`https://kind-success-production.up.railway.app`

Update your frontend to use this URL:

**In `frontend/src/firebase.js` or wherever you define API_URL:**
```javascript
const API_URL = "https://kind-success-production.up.railway.app";
```

Then redeploy frontend:
```bash
railway up
```

---

## Quick Commands Summary

```bash
# Backend deployment
cd backend
railway up
railway variables set GOOGLE_API_KEY=your_key

# Frontend deployment
cd ../frontend
# Update API_URL in code first!
railway up
```

---

## Troubleshooting

### "Available options can not be empty"
- This means you need to create a new service first
- Go to Railway dashboard → New → Empty Service
- Then run `railway link` again

### Deployment fails with module errors
- Use `requirements-light.txt` instead of full requirements
- Remove heavy dependencies (AiZynthFinder, etc.)

### Port binding errors
- Make sure your start command uses `$PORT` variable
- Railway assigns a random port, you must use their variable

---

## Alternative: Use Railway Dashboard (Easier)

1. Go to https://railway.app/dashboard
2. Click "New Project"
3. Select "Deploy from GitHub repo"
4. Connect your GitHub account
5. Push your code to GitHub first:
   ```bash
   git init
   git add .
   git commit -m "Initial commit"
   gh repo create phoenix-pi --public --source=. --push
   ```
6. Select the repo in Railway
7. Railway auto-deploys on every push!

This is actually **easier** than CLI deployment.
