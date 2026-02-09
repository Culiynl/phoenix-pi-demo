# Render.com Deployment - Simple & Reliable

## ✅ Why Render is Better for This Project

- **No CLI needed** - Everything through web dashboard
- **Free tier** - 750 hours/month free
- **Auto-deploys** from GitHub
- **More reliable** than Railway CLI
- **Better logs** and debugging

---

## Step-by-Step Deployment

### Step 1: Push to GitHub (You've Already Done This!)

Your code is ready to push. Just need to create the repo:

1. Go to https://github.com/new
2. Name: `phoenix-pi-demo`
3. Public
4. Create repository
5. Then run: `git push -u origin main`

**OR** skip GitHub entirely - Render can deploy from local git!

---

### Step 2: Sign Up for Render

1. Go to https://render.com
2. Click "Get Started"
3. Sign up with GitHub (easiest) or email

---

### Step 3: Deploy Backend

1. Click **"New +"** → **"Web Service"**
2. Connect your GitHub account (or use "Public Git Repository")
3. If using GitHub: Select `phoenix-pi-demo` repo
4. If using public git: Enter `https://github.com/Culiynl/phoenix-pi-demo`

**Configure Backend:**
- **Name**: `phoenix-backend`
- **Root Directory**: `backend`
- **Runtime**: `Python 3`
- **Build Command**: `pip install -r requirements-light.txt`
- **Start Command**: `uvicorn api:app --host 0.0.0.0 --port $PORT`

**Environment Variables** (click "Advanced"):
- `GOOGLE_API_KEY` = your actual key
- `GEMINI_MODEL` = `gemini-2.0-flash`
- `PYTHON_VERSION` = `3.10.0`

**Instance Type**: Free

Click **"Create Web Service"**

Render will build and deploy! Takes ~5 minutes.

---

### Step 4: Get Backend URL

After deployment, Render gives you a URL like:
`https://phoenix-backend.onrender.com`

Copy this URL!

---

### Step 5: Update Frontend API URL

Edit `frontend/src/firebase.js` or wherever you define the API URL:

```javascript
const API_URL = "https://phoenix-backend.onrender.com";
```

Commit and push:
```bash
git add .
git commit -m "Update API URL for Render"
git push
```

---

### Step 6: Deploy Frontend

1. Click **"New +"** → **"Static Site"**
2. Select same `phoenix-pi-demo` repo

**Configure Frontend:**
- **Name**: `phoenix-frontend`
- **Root Directory**: `frontend`
- **Build Command**: `npm install && npm run build`
- **Publish Directory**: `build`

Click **"Create Static Site"**

---

### Step 7: Test Your Deployment!

Frontend URL: `https://phoenix-frontend.onrender.com`

Visit it and test:
- ✅ Create a new mission
- ✅ Upload a file
- ✅ Use standalone tools
- ✅ Check if Gemini responds

---

## 🎬 You're Done!

Your app is now live at:
- **Frontend**: https://phoenix-frontend.onrender.com
- **Backend**: https://phoenix-backend.onrender.com

**Free tier notes:**
- Services spin down after 15 min of inactivity
- First request after spin-down takes ~30 seconds
- Perfect for demos and competitions!

---

## Troubleshooting

### Build fails with "requirements not found"
- Make sure Root Directory is set to `backend` or `frontend`
- Check that `requirements-light.txt` exists

### Frontend can't reach backend
- Check CORS settings in `api.py`
- Verify API_URL in frontend code
- Check Render logs for errors

### "Out of memory" during build
- Use `requirements-light.txt` instead of full requirements
- Render free tier has 512MB RAM limit

---

## Alternative: Skip GitHub Entirely

Render can deploy from a public git URL:

1. Make your local git repo public somehow, OR
2. Use Render's "Manual Deploy" option
3. Upload a zip file of your code

But GitHub integration is recommended for auto-deploys!
