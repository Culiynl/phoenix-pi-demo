"""
Multimodal Analysis Service
Provides Gemini-powered analysis to enhance existing tools (P2Rank, etc.)
"""

import os
import google.generativeai as genai
from config import GOOGLE_API_KEY, IMG_DIR, DATA_DIR

genai.configure(api_key=GOOGLE_API_KEY)


class MultimodalAnalyzer:
    """Enhances existing bioinformatics tools with Gemini's multimodal capabilities."""
    
    @staticmethod
    async def analyze_p2rank_results(p2rank_output: dict, pdb_id: str, mission_id: str):
        """
        Analyzes P2Rank pocket detection results and suggests next steps.
        
        Args:
            p2rank_output: P2Rank results with pocket predictions
            pdb_id: PDB structure ID
            mission_id: Mission ID for logging
            
        Returns:
            dict: Gemini's analysis and recommendations
        """
        model = genai.GenerativeModel('gemini-2.0-flash')
        
        try:
            # Format P2Rank results for Gemini
            pockets_summary = []
            for i, pocket in enumerate(p2rank_output.get('pockets', [])[:5]):
                pockets_summary.append(f"""
Pocket {i+1}:
- Score: {pocket.get('score', 'N/A')}
- Residues: {', '.join(pocket.get('residues', [])[:10])}
- Center: {pocket.get('center', 'N/A')}
""")
            
            pockets_text = "\n".join(pockets_summary)
            
            prompt = f"""You are analyzing protein binding pocket predictions from P2Rank for PDB structure {pdb_id}.

# P2Rank Results
{pockets_text}

# Analysis Tasks

Based on these pocket predictions, provide:

1. **Best Pocket Selection**: Which pocket(s) look most druggable and why?

2. **Ligand Design Strategy**: What types of molecular interactions would be favorable for each pocket?

3. **Next Steps**: Should we:
   - Proceed with virtual screening against the top pocket?
   - Search ChEMBL for known binders to similar pockets?
   - Design de novo ligands for this pocket?

4. **Risk Assessment**: Any concerns about these pockets (e.g., too shallow, too hydrophobic)?

Provide actionable recommendations."""

            response = await model.generate_content_async(prompt)
            
            analysis_text = response.text
            
            # Save to Firestore
            from services.firebase_service import update_mission_field
            from datetime import datetime
            
            update_mission_field(mission_id, "p2rank_analysis", {
                "pdb_id": pdb_id,
                "pockets_analyzed": len(p2rank_output.get('pockets', [])),
                "gemini_recommendations": analysis_text,
                "timestamp": datetime.now().isoformat()
            })
            
            return {
                "success": True,
                "analysis": analysis_text,
                "pockets_count": len(p2rank_output.get('pockets', []))
            }
            
        except Exception as e:
            return {
                "success": False,
                "error": str(e)
            }

        """
        Analyzes protein structure images to identify binding sites and druggable pockets.
        
        Args:
            image_path: Absolute path to the protein structure image
            mission_id: Mission ID for logging
            
        Returns:
            dict: Analysis results with binding sites, residues, and coordinates
        """
        model = genai.GenerativeModel('gemini-2.0-flash-exp')
        
        try:
            with open(image_path, 'rb') as f:
                image_data = f.read()
            
            prompt = """Analyze this protein structure image and provide:

1. **Binding Pocket Identification**: Identify potential druggable binding pockets visible in the structure
2. **Key Residues**: List important amino acid residues that appear to be in or near binding sites
3. **Structural Features**: Describe relevant structural features (helices, sheets, loops) near potential binding sites
4. **Drug Design Insights**: Suggest what types of molecular interactions would be favorable for ligand binding

Provide your analysis in a structured format with clear sections."""

            response = await model.generate_content_async([
                prompt,
                {"mime_type": "image/png", "data": image_data}
            ])
            
            analysis_text = response.text
            
            # Save to Firestore
            from services.firebase_service import update_mission_field
            update_mission_field(mission_id, "structure_analysis", {
                "image_path": image_path,
                "analysis": analysis_text,
                "timestamp": datetime.now().isoformat()
            })
            
            return {
                "success": True,
                "analysis": analysis_text,
                "image_path": image_path
            }
            
        except Exception as e:
            return {
                "success": False,
                "error": str(e)
            }
    
    @staticmethod
    async def analyze_research_pdf(pdf_path: str, mission_id: str):
        """
        Extracts insights from research paper PDFs.
        
        Args:
            pdf_path: Absolute path to the PDF file
            mission_id: Mission ID for logging
            
        Returns:
            dict: Extracted insights including molecules, targets, and findings
        """
        model = genai.GenerativeModel('gemini-2.0-flash-exp')
        
        try:
            # Upload PDF to Gemini
            uploaded_file = genai.upload_file(pdf_path)
            
            prompt = """Analyze this research paper and extract:

1. **Main Findings**: Key discoveries and conclusions
2. **Molecular Targets**: Proteins, enzymes, or receptors studied
3. **Drug Candidates**: Any molecules, compounds, or SMILES strings mentioned
4. **Experimental Data**: IC50, EC50, binding affinities, or other quantitative results
5. **Therapeutic Applications**: Diseases or conditions being targeted
6. **Novel Insights**: Unique contributions or unexpected findings

Format the output as structured sections with bullet points."""

            response = await model.generate_content_async([uploaded_file, prompt])
            
            analysis_text = response.text
            
            # Save to Firestore
            from services.firebase_service import update_mission_field
            from datetime import datetime
            
            update_mission_field(mission_id, "pdf_analysis", {
                "pdf_path": pdf_path,
                "analysis": analysis_text,
                "timestamp": datetime.now().isoformat()
            })
            
            return {
                "success": True,
                "analysis": analysis_text,
                "pdf_path": pdf_path
            }
            
        except Exception as e:
            return {
                "success": False,
                "error": str(e)
            }
    
    @staticmethod
    async def sketch_to_smiles(sketch_image_path: str):
        """
        Converts hand-drawn molecular structure sketches to SMILES notation.
        
        Args:
            sketch_image_path: Path to the sketch image
            
        Returns:
            dict: SMILES string and confidence
        """
        model = genai.GenerativeModel('gemini-2.0-flash-exp')
        
        try:
            with open(sketch_image_path, 'rb') as f:
                image_data = f.read()
            
            prompt = """Analyze this hand-drawn chemical structure and convert it to SMILES notation.

Instructions:
1. Identify all atoms (C, N, O, S, etc.)
2. Identify all bonds (single, double, triple, aromatic)
3. Determine stereochemistry if indicated
4. Generate the SMILES string

Respond ONLY with the SMILES string on the first line, followed by your confidence level (high/medium/low) on the second line."""

            response = await model.generate_content_async([
                prompt,
                {"mime_type": "image/png", "data": image_data}
            ])
            
            lines = response.text.strip().split('\n')
            smiles = lines[0].strip()
            confidence = lines[1].strip() if len(lines) > 1 else "unknown"
            
            return {
                "success": True,
                "smiles": smiles,
                "confidence": confidence
            }
            
        except Exception as e:
            return {
                "success": False,
                "error": str(e)
            }
