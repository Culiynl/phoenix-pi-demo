# PHOENIX-PI Deployment Guide

## Architecture Overview

```
┌──────────────────────────────────────────┐
│     Vercel (Frontend - Static)           │
│  - React app                             │
│  - Connects to Cloud Run backend         │
└──────────────────────────────────────────┘
                    │
                    ▼
┌──────────────────────────────────────────┐
│   Cloud Run (Lightweight Backend)        │
│  - Gemini orchestration                  │
│  - ChEMBL queries                        │
│  - Firestore operations                  │
│  - Multimodal analysis                   │
│  - Long-context synthesis                │
└──────────────────────────────────────────┘
                    │
                    │ (API calls for heavy compute)
                    ▼
┌──────────────────────────────────────────┐
│  Compute Engine VM (Heavy Compute)       │
│  - AiZynthFinder (retrosynthesis)        │
│  - P2Rank (docking)                      │
│  - Large ML models                       │
└──────────────────────────────────────────┘
```

## Step 1: Deploy Heavy Compute Service (Compute Engine)

### Create VM
```bash
gcloud compute instances create phoenix-compute \
    --zone=us-central1-a \
    --machine-type=n1-standard-4 \
    --image-family=debian-11 \
    --image-project=debian-cloud \
    --boot-disk-size=50GB
```

### SSH into VM and setup
```bash
gcloud compute ssh phoenix-compute --zone=us-central1-a

# Install Python and dependencies
sudo apt-get update
sudo apt-get install -y python3-pip python3-venv

# Clone repo or copy files
git clone YOUR_REPO_URL
cd webserver_stuff/backend

# Install heavy dependencies
pip3 install -r requirements.txt

# Run heavy compute service
python3 heavy_compute_service.py
```

### Get VM IP
```bash
gcloud compute instances describe phoenix-compute \
    --zone=us-central1-a \
    --format='get(networkInterfaces[0].accessConfigs[0].natIP)'
```

## Step 2: Deploy Backend to Cloud Run

### Update config with VM IP
Edit `config.py`:
```python
HEAVY_COMPUTE_URL = "http://YOUR_VM_IP:8001"
```

### Build and deploy
```bash
cd backend

# Build container
gcloud builds submit --tag gcr.io/YOUR_PROJECT_ID/phoenix-backend -f Dockerfile.light

# Deploy to Cloud Run
gcloud run deploy phoenix-backend \
    --image gcr.io/YOUR_PROJECT_ID/phoenix-backend \
    --platform managed \
    --region us-central1 \
    --allow-unauthenticated \
    --memory 2Gi \
    --timeout 300
```

### Get Cloud Run URL
```bash
gcloud run services describe phoenix-backend \
    --region us-central1 \
    --format='value(status.url)'
```

## Step 3: Deploy Frontend to Vercel

### Update frontend config
Edit `frontend/src/firebase.js`:
```javascript
const API_URL = "YOUR_CLOUD_RUN_URL";
```

### Deploy
```bash
cd frontend

# Install Vercel CLI
npm i -g vercel

# Deploy
vercel --prod
```

## Step 4: Test End-to-End

1. Visit your Vercel URL
2. Create a new mission
3. Verify:
   - Gemini orchestration works (Cloud Run)
   - ChEMBL queries work (Cloud Run)
   - Retrosynthesis works (Compute Engine)

## Cost Estimate

- **Cloud Run**: ~$0 (free tier covers demo usage)
- **Compute Engine**: ~$0.05/hour ($1.20/day if running 24/7)
- **Vercel**: $0 (free tier)

**Total**: ~$5-10/month for demo

## Optimization Tips

1. **Use preemptible VM**: Reduce Compute Engine cost by 80%
   ```bash
   --preemptible
   ```

2. **Auto-shutdown VM**: Stop when not in use
   ```bash
   gcloud compute instances stop phoenix-compute --zone=us-central1-a
   ```

3. **Use Cloud Run min instances**: Keep 0 to avoid idle costs
   ```bash
   --min-instances=0
   ```

## Troubleshooting

### Cloud Run timeout
- Increase timeout: `--timeout 600`
- Check logs: `gcloud run logs read phoenix-backend`

### Compute Engine connection refused
- Check firewall: Allow port 8001
  ```bash
  gcloud compute firewall-rules create allow-heavy-compute \
      --allow tcp:8001
  ```

### Frontend can't reach backend
- Check CORS in `api.py`
- Verify Cloud Run URL in frontend config
