"""
Long Context Service
Implements deep literature synthesis and large-scale data analysis using Gemini's 2M token context window.
"""

import google.generativeai as genai
from config import GOOGLE_API_KEY
from typing import List, Dict
import asyncio

genai.configure(api_key=GOOGLE_API_KEY)


class LongContextAnalyzer:
    """Handles long-context analysis tasks using Gemini's 2M token window."""
    
    @staticmethod
    async def deep_literature_synthesis(paper_ids: List[str], mission_id: str, max_papers: int = 50):
        """
        Synthesizes insights from up to 50 research papers using long context.
        
        Args:
            paper_ids: List of ArXiv paper IDs
            mission_id: Mission ID for logging
            max_papers: Maximum number of papers to analyze (default 50)
            
        Returns:
            dict: Comprehensive synthesis across all papers
        """
        model = genai.GenerativeModel('gemini-2.0-flash')
        
        try:
            # Fetch paper abstracts and metadata
            from services.llm_service import LLMService
            llm = LLMService()
            
            papers_text = []
            for i, pid in enumerate(paper_ids[:max_papers]):
                try:
                    # Search for the specific paper
                    results = await llm.search_arxiv_papers(pid, max_results=1)
                    if results:
                        paper = results[0]
                        papers_text.append(f"""
## Paper {i+1}: {paper.get('title', 'Unknown')}
**Authors**: {paper.get('authors', 'N/A')}
**Published**: {paper.get('published', 'N/A')}
**ArXiv ID**: {pid}

### Abstract
{paper.get('summary', 'No abstract available')}

---
""")
                except Exception as e:
                    print(f"Error fetching paper {pid}: {e}")
                    continue
            
            if not papers_text:
                return {"success": False, "error": "No papers could be fetched"}
            
            combined_text = "\n".join(papers_text)
            
            prompt = f"""You are analyzing {len(papers_text)} research papers in the drug discovery domain.

# Papers to Analyze
{combined_text}

# Analysis Tasks

Provide a comprehensive synthesis addressing:

1. **Common Therapeutic Targets**: What proteins, enzymes, or receptors are most frequently studied across these papers?

2. **Novel Drug Candidates**: Identify any promising molecules, compounds, or SMILES strings mentioned. List them with their reported activities.

3. **Emerging Trends**: What are the key trends or hot topics in this research area based on these papers?

4. **Methodological Insights**: What experimental techniques or computational methods are commonly used?

5. **Contradictions & Debates**: Are there any conflicting findings or ongoing debates in the literature?

6. **Research Gaps**: What questions remain unanswered? What future directions are suggested?

Format your response with clear section headers and bullet points."""

            response = await model.generate_content_async(prompt)
            
            synthesis_text = response.text
            
            # Save to Firestore
            from services.firebase_service import update_mission_field
            from datetime import datetime
            
            update_mission_field(mission_id, "deep_synthesis", {
                "papers_analyzed": len(papers_text),
                "synthesis": synthesis_text,
                "timestamp": datetime.now().isoformat()
            })
            
            return {
                "success": True,
                "papers_analyzed": len(papers_text),
                "synthesis": synthesis_text
            }
            
        except Exception as e:
            return {
                "success": False,
                "error": str(e)
            }
    
    @staticmethod
    async def analyze_full_chembl_dataset(target_name: str, mission_id: str):
        """
        Loads a large ChEMBL dataset into context for comprehensive SAR analysis.
        
        Args:
            target_name: Protein target name
            mission_id: Mission ID for logging
            
        Returns:
            dict: Comprehensive SAR insights
        """
        model = genai.GenerativeModel('gemini-2.0-flash-exp')
        
        try:
            from services.bio_service import BioDataManager
            
            # Fetch a larger dataset (up to 500 molecules)
            results = BioDataManager.search_chembl(target_name, max_results=500)
            
            if not results or 'molecules' not in results:
                return {"success": False, "error": "No ChEMBL data found"}
            
            molecules = results['molecules']
            
            # Format molecules for context
            mol_text = []
            for i, mol in enumerate(molecules[:500]):
                mol_text.append(f"""
Molecule {i+1}:
- SMILES: {mol.get('smiles', 'N/A')}
- pIC50: {mol.get('pIC50', 'N/A')}
- Molecular Weight: {mol.get('molwt', 'N/A')}
- LogP: {mol.get('logp', 'N/A')}
- QED: {mol.get('qed', 'N/A')}
""")
            
            combined_mol_text = "\n".join(mol_text)
            
            prompt = f"""You are analyzing {len(molecules)} molecules from ChEMBL for target: {target_name}

# Molecular Dataset
{combined_mol_text}

# Analysis Tasks

Provide a comprehensive Structure-Activity Relationship (SAR) analysis:

1. **Potency Patterns**: What structural features correlate with high pIC50 values?

2. **Physicochemical Trends**: How do MW, LogP, and QED relate to potency?

3. **Scaffold Analysis**: Identify common molecular scaffolds or core structures.

4. **Optimization Opportunities**: Based on the data, what modifications might improve potency while maintaining drug-likeness?

5. **Outliers**: Are there any unusual molecules with unexpected properties?

6. **Design Recommendations**: What would be your top 3 recommendations for designing new molecules?

Format your response with clear sections and specific examples (use SMILES when relevant)."""

            response = await model.generate_content_async(prompt)
            
            sar_analysis = response.text
            
            # Save to Firestore
            from services.firebase_service import update_mission_field
            from datetime import datetime
            
            update_mission_field(mission_id, "full_chembl_analysis", {
                "molecules_analyzed": len(molecules),
                "target": target_name,
                "sar_analysis": sar_analysis,
                "timestamp": datetime.now().isoformat()
            })
            
            return {
                "success": True,
                "molecules_analyzed": len(molecules),
                "sar_analysis": sar_analysis
            }
            
        except Exception as e:
            return {
                "success": False,
                "error": str(e)
            }
