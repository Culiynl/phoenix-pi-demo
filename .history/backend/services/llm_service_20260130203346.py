import json
import google.generativeai as genai
from Bio import Entrez
from config import GOOGLE_API_KEY, ENTREZ_EMAIL
from schemas import DiseaseQuery

# Setup
genai.configure(api_key=GOOGLE_API_KEY)
Entrez.email = ENTREZ_EMAIL

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
    model = genai.GenerativeModel('gemini-1.5-flash')
    res = model.generate_content(
        prompt + "\nOutput JSON: { 'targets': [ { 'name': '...', 'target_type': '...', 'druggability': '...' } ] }"
    )
    
    try:
        clean = res.text.replace('```json','').replace('```','').strip()
        return json.loads(clean)
    except:
        return {"targets": []}