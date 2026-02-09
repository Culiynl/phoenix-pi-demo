# GitHub Authentication Quick Fix

## Problem
Git push is failing because you need to authenticate with GitHub.

## Solution: Personal Access Token (PAT)

### Step 1: Create a Personal Access Token

1. Go to **https://github.com/settings/tokens**
2. Click **"Generate new token"** → **"Generate new token (classic)"**
3. Give it a name: `phoenix-pi-deployment`
4. Set expiration: 90 days (or No expiration for demo)
5. Select scopes:
   - ✅ **repo** (all repo permissions)
6. Click **"Generate token"**
7. **COPY THE TOKEN** - you won't see it again!

### Step 2: Push with Token

When you run `git push`, it will ask for credentials:

```
Username: Culiynl
Password: <paste your token here, NOT your GitHub password>
```

**Or** set it up once:

```bash
git remote set-url origin https://Culiynl:YOUR_TOKEN@github.com/Culiynl/phoenix-pi-demo.git
git push -u origin main
```

Replace `YOUR_TOKEN` with the token you just copied.

---

## Alternative: GitHub Desktop (Easier)

1. Download **GitHub Desktop**: https://desktop.github.com/
2. Sign in with your GitHub account
3. File → Add Local Repository
4. Select: `f:\Coding Stuff\Full Model 2026 Fixed\webserver_stuff`
5. Click **"Push origin"** button

GitHub Desktop handles authentication automatically!

---

## After Successful Push

Once your code is on GitHub:

1. Go to **https://render.com**
2. Sign up with GitHub
3. New Web Service → Select `phoenix-pi-demo`
4. Deploy!

No more git commands needed - Render auto-deploys from GitHub.
