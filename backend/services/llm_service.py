import json
import google.generativeai as genai
from Bio import Entrez
from config import GOOGLE_API_KEY, ENTREZ_EMAIL
from schemas import DiseaseQuery

# Setup
genai.configure(api_key=GOOGLE_API_KEY)
Entrez.email = ENTREZ_EMAIL

def analyze_full_literature(pdf_paths: list[str], target_disease: str):
    """
    PHOENIX-PI: Uses Gemini 1.5 Pro's 1M context window to analyze full PDF texts.
    """
    model = genai.GenerativeModel('gemini-3-pro-preview')
    
    # In a real implementation:
    # 1. Upload PDFs to Gemini File API
    # 2. Build prompt: "Analyze these papers and design a drug discovery mission for {target_disease}"
    # 3. Return a multi-step Plan for ResearchOrchestrator
    
    prompt = (
        f"You are PHOENIX-PI, an autonomous research agent. Analyze the full text of the "
        f"provided {len(pdf_paths)} papers and design a detailed mission to target {target_disease}. "
        "Output a multi-step JSON plan: { 'mission_plan': [ {'step': 1, 'goal': '...', 'tool': '...'} ] }"
    )
    
    # Placeholder for actual multimodal execution
    return {"plan": "Mission plan generated from full-text ArXiv PDFs (Placeholder)."}

def review_literature(query: DiseaseQuery):
    """Fetches PubMed abstracts and uses LLM to extract targets."""
    search = f"{query.disease_name} AND (drug target OR mechanism)"
    handle = Entrez.esearch(db="pubmed", term=search, retmax=query.max_papers)
    ids = Entrez.read(handle)["IdList"]
    
    if not ids:
        return {"targets": []}
    
    handle = Entrez.efetch(db="pubmed", id=ids, rettype="xml", retmode="text")
    papers = Entrez.read(handle)
    
    abstracts = []
    for p in papers['PubmedArticle']:
        try:
            abstracts.append(p['MedlineCitation']['Article']['Abstract']['AbstractText'][0])
        except:
            pass
    
    prompt = f"Extract drug targets for {query.disease_name} from these abstracts:\n" + "\n".join(abstracts[:5])
    model = genai.GenerativeModel('gemini-3-pro-preview')
    res = model.generate_content(
        prompt + "\nOutput JSON: { 'targets': [ { 'name': '...', 'target_type': '...', 'druggability': '...' } ] }"
    )
    
    try:
        clean = res.text.replace('```json','').replace('```','').strip()
        return json.loads(clean)
    except Exception:
        return {"targets": []}

def search_arxiv_papers(query: str, max_results=5):
    """Searches ArXiv for relevant papers with respectful rate limiting."""
    import arxiv
    import time
    print(f"DEBUG: Starting ArXiv search for: {query}")
    try:
        # Respectful client configuration
        client = arxiv.Client(
            page_size=100,
            delay_seconds=3.0,
            num_retries=3
        )
        search = arxiv.Search(
            query=query,
            max_results=max_results,
            sort_by=arxiv.SortCriterion.Relevance
        )
        
        results = []
        # Use an iterator approach to catch 429s during the loop
        results_iter = client.results(search)
        
        for _ in range(max_results):
            try:
                r = next(results_iter)
                results.append({
                    "title": getattr(r, 'title', 'Untitled'),
                    "authors": [getattr(a, 'name', 'Unknown') for a in getattr(r, 'authors', [])],
                    "summary": getattr(r, 'summary', ''),
                    "pdf_url": getattr(r, 'pdf_url', ''),
                    "published": r.published.strftime("%Y-%m-%d") if hasattr(r, 'published') and r.published else "Unknown"
                })
            except StopIteration:
                break
            except Exception as e:
                print(f"ERROR: ArXiv iteration error: {e}")
                break # Return whatever we have so far
                
        print(f"DEBUG: Found {len(results)} papers.")
        return results
    except Exception as e:
        print(f"ERROR: ArXiv Search Failed: {e}")
        return []

def fetch_paper_full_text(pdf_url: str):
    """Downloads and extracts text from an ArXiv PDF."""
    import requests
    import io
    from pypdf import PdfReader
    
    try:
        response = requests.get(pdf_url, timeout=30)
        if response.status_code == 200:
            f = io.BytesIO(response.content)
            reader = PdfReader(f)
            text = ""
            # Extract first 10 pages to avoid context overflow/bloat
            for i in range(min(len(reader.pages), 10)):
                text += reader.pages[i].extract_text() + "\n"
            return text
        return f"Failed to download PDF: HTTP {response.status_code}"
    except Exception as e:
        return f"Error extracting PDF: {str(e)}"

def summarize_literature_gemini(papers: list, full_texts: list = None, existing_summary: str = None):
    """
    Summarizes a list of papers using gemini-3-pro-preview.
    If existing_summary is provided, it incorporates new findings into it.
    """
    if not papers: return existing_summary or "No papers to summarize."
    
    if full_texts:
        context_parts = []
        for i, p in enumerate(papers):
            content = full_texts[i] if i < len(full_texts) else p['summary']
            context_parts.append(f"New Paper {i+1} Title: {p['title']}\nContent:\n{content[:5000]}...")
        context = "\n\n---\n\n".join(context_parts)
    else:
        context = "\n\n".join([f"New Paper Title: {p['title']}\nAbstract: {p['summary']}" for p in papers])
        
    prompt = (
        "You are an expert scientific synthesizer. Analyze the provided scientific context and merge it with the existing research summary if provided.\n"
        "GOAL: Create a single, cohesive strategic summary that captures the most important targets, mechanisms, and structural data discovered so far.\n"
        "IMPORTANT: \n"
        "1. Do NOT just append; integrate the information.\n"
        "2. Keep the output concise and high-level (approx. 4-6 sentences total).\n"
        "3. Focus on what informs a drug discovery mission."
    )
    
    if existing_summary:
        full_prompt = f"{prompt}\n\nEXISTING SUMMARY:\n{existing_summary}\n\nNEW CONTEXT TO INTEGRATE:\n{context}"
    else:
        full_prompt = f"{prompt}\n\nNEW CONTEXT:\n{context}"
    
    model = genai.GenerativeModel('gemini-3-pro-preview')
    try:
        res = model.generate_content(full_prompt)
        return res.text
    except Exception as e:
        return f"Error generating summary: {str(e)}"

def generate_mission_title(prompt: str) -> str:
    """Generates a short, catchy title/target name from a user description."""
    model = genai.GenerativeModel('gemini-3-pro-preview')
    try:
        res = model.generate_content(
            f"Generate a short (3-5 words) scientific project title for this request: '{prompt}'. "
            "Output ONLY the title, no quotes or extra text."
        )
        return res.text.strip().replace('"', '').replace("'", "")
    except:
        return "New Mission"

def extract_search_queries(prompt: str) -> list[str]:
    """Extracts key scientific terms for ArXiv/PDB searching."""
    model = genai.GenerativeModel('gemini-3-pro-preview')
    try:
        res = model.generate_content(
            f"Extract 2-3 specific scientific search queries (e.g., protein names, disease mechanisms) "
            f"from this text: '{prompt}'. Output as JSON list: ['query1', 'query2']"
        )
        text = res.text.strip()
        if text.startswith("```json"): text = text.replace("```json", "").replace("```", "")
        return json.loads(text)
    except:
        return [prompt]