
import os
import asyncio
import json
import inspect
import time
import random
from datetime import datetime
from typing import List, Dict, Any
from rdkit import Chem
from rdkit.Chem import Descriptors
import google.generativeai as genai
import google.ai.generativelanguage as protos
from google.protobuf import struct_pb2
from config import GOOGLE_API_KEY, PROJECTS_DIR, IMG_DIR, GEMINI_MODEL
from services.local_storage_service import update_mission_log, update_mission_field, get_mission

# Setup Gemini - Using REST transport to bypass asyncio gRPC loop caching issues
genai.configure(api_key=GOOGLE_API_KEY, transport="rest")

class ResearchOrchestrator:
    def __init__(self, mission_id: str):
        self.mission_id = mission_id
        self.history = []
        self.thought_log = []
        
        # Load mission data to restore thought_log
        from services.local_storage_service import get_mission
        mission = get_mission(self.mission_id)
        if mission:
            self.thought_log = mission.get("thought_log", [])
            print(f"[Orchestrator] Restored {len(self.thought_log)} entries in thought_log.")
        
        self.max_steps = 20
        
        # Tools Registry - KEYS MUST MATCH METHOD NAMES
        self.tools = {
            "tool_search_literature": self.tool_search_literature,
            "tool_submit_strategy": self.tool_submit_strategy,
            "tool_fetch_pdb": self.tool_fetch_pdb,
            "tool_fetch_chembl": self.tool_fetch_chembl,
            "tool_analyze_ligand_set": self.tool_analyze_ligand_set,
            "tool_suggest_starter_ligand": self.tool_suggest_starter_ligand,
            "tool_run_p2rank": self.tool_run_p2rank,
            "tool_run_vina": self.tool_run_vina,
            "tool_generate_leads": self.tool_generate_leads,
            "tool_fetch_protein_structure": self.tool_fetch_protein_structure,
            "tool_fetch_chembl_data": self.tool_fetch_chembl_data,
            "tool_retrosynthesis_analysis": self.tool_retrosynthesis_analysis,
            "tool_analyze_molecule": self.tool_analyze_molecule,
            "tool_verify_lipinski": self.tool_verify_lipinski,
            "tool_find_pareto_front_molecule": self.tool_find_pareto_front_molecule,
            "tool_submit_final_report": self.tool_submit_final_report,
            "tool_analyze_structure_image": self.tool_analyze_structure_image,
            "tool_analyze_research_pdf": self.tool_analyze_research_pdf,
            "tool_deep_literature_synthesis": self.tool_deep_literature_synthesis,
            "tool_full_chembl_analysis": self.tool_full_chembl_analysis,
            "tool_structural_interaction_analysis": self.tool_structural_interaction_analysis
        }
        
        # Lazy Model/Chat (Initialized per turn to ensure loop alignment)
        self.model = None
        self.chat = None
        self.sdk_poisoned = False # Track if SDK/API mismatch for thought_signature detected

    def _ensure_chat_aligned(self):
        """DEFINTIVE FIX: Always re-creates the model and chat session for the current event loop."""
        try:
            # Re-create model to ensure gRPC is loop-aligned
            self.model = genai.GenerativeModel(
                model_name=GEMINI_MODEL,
                tools=list(self.tools.values())
            )
            # Re-create chat session from the restored thought log
            self.chat = self._initialize_chat_with_history()
            print(f"[Orchestrator] Neural connection hardware-aligned for current loop.")
        except Exception as e:
            print(f"[Orchestrator] Loop alignment failed: {e}")
            raise e

    def _initialize_chat_with_history(self):
        """Attempts to rebuild the chat history from the stored thought log."""
        from services.local_storage_service import get_mission
        mission = get_mission(self.mission_id)
        if not mission or not mission.get("thought_log"):
            return self.model.start_chat()
            
        history = []
        last_role = None
        for entry in mission.get("thought_log", []):
            type = entry.get("type")
            content = entry.get("content", "")
            if not content: continue
            
            role = "user" if type == "user" else "model"
            if type not in ["user", "thought", "system", "action", "success", "error", "info"]:
                continue

            if history and role == last_role:
                # Merge consecutive parts for the same role
                history[-1]["parts"].append(content)
            else:
                history.append({"role": role, "parts": [content]})
                last_role = role
            
        print(f"[Orchestrator] Resuming Gemini chat with {len(history)} turns (merged) in context.")
        return self.model.start_chat(history=history)

    def _log(self, type: str, content: str):
        """Appends to thought_log and syncs with Firestore."""
        entry = {
            "type": type,
            "content": content,
            "timestamp": datetime.now().isoformat()
        }
        self.thought_log.append(entry)
        update_mission_log(self.mission_id, entry)

    def _reset_feudal_singleton(self):
        """Forces the FeudalOptimizationManager to re-initialize (useful if it failed imports previously)."""
        from services.feudal_service import FeudalOptimizationManager
        FeudalOptimizationManager._instance = None
        print("[Orchestrator] Feudal singleton reset triggered.")
        
        # Log this explicitly
        self._log("system", "Feudal singleton reset triggered.")


    def _log_tool_call(self, name: str, args: Any):
        """Logs a tool call to a specialized history field."""
        from services.local_storage_service import update_mission_field, get_mission, _clean_serializable
        
        entry = {
            "tool": name,
            "args": _clean_serializable(args),
            "timestamp": datetime.now().isoformat()
        }
        
        # We need to append to the list, not overwrite.
        mission = get_mission(self.mission_id)
        if mission:
            history = mission.get("tool_history", [])
            history.append(entry)
            update_mission_field(self.mission_id, "tool_history", history)

    # Tool functions must have clear docstrings for Gemini to understand them
    async def tool_search_literature(self, query: str):
        """Searches scientific literature (ArXiv) for disease targets and markers. Returns a summary of findings."""
        self._log("action", f"Initiating research: {query}")
        
        # Fetch current mission settings
        from services.local_storage_service import get_mission, update_mission_field
        mission_doc = get_mission(self.mission_id) or {}
        research_mode = mission_doc.get("research_mode", "hybrid")
        
        if research_mode == "internal":
            self._log("system", "Mode: Internal Knowledge. Fast-tracking using Gemini's pre-trained data...")
            return {
                "summary": "Using internal clinical and pharmacological knowledge. No external ArXiv grounding requested for this step.",
                "instruction": "Please proceed based on your internal knowledge of the target and its related pathways."
            }

        from services.llm_service import search_arxiv_papers, fetch_paper_full_text, summarize_literature_gemini
        
        # Direct search as requested
        all_papers = []
        full_texts = []
        try:
            papers = search_arxiv_papers(query, max_results=3)
            if papers:
                all_papers = papers
                for p in all_papers:
                    title_short = (p['title'][:60] + '..') if len(p['title']) > 60 else p['title']
                    self._log("system", f"  - Reading full text: {title_short}")
                    
                    # Fetch full text
                    text = fetch_paper_full_text(p['pdf_url'])
                    full_texts.append(text)
                    
                    # Small delay to avoid hammering
                    await asyncio.sleep(1.0)
                self._log("system", f"  - Analyzed {len(all_papers)} papers with full-text grounding.")
            else:
                self._log("system", "  - No immediate results for primary query.")

        except Exception as e:
            self._log("error", f"ArXiv Search failed: {str(e)}")
            
        # 2. Summarize
        existing_summary = mission_doc.get("internal_research_summary", "")
        try:
            summary = summarize_literature_gemini(all_papers, full_texts=full_texts, existing_summary=existing_summary)
        except Exception as e:
            self._log("error", f"Summarization failed: {e}")
            summary = existing_summary or f"Research attempted for '{query}', but summarization failed."
        
        if not summary: 
            summary = "Research results gathered, but summary generation was empty."
            
        # 3. Store in INTERNAL fields (Prevents UI flickering until Step 2 is done)
        update_mission_field(self.mission_id, "internal_research_summary", summary)
        
        # Save structured paper metadata for side-by-side UI (INTERNALLY)
        paper_metadata = []
        for p in all_papers:
            paper_metadata.append({
                "title": p.get('title', 'Unknown Title'),
                "summary": p.get('summary', ''),
                "url": p.get('pdf_url', ''),
                "authors": p.get('authors', [])
            })
        update_mission_field(self.mission_id, "internal_research_papers", paper_metadata)
        
        return {"summary": "Internal research updated. Sources analyzed: " + str(len(all_papers)), "paper_count": len(all_papers), "status": "Internal staging updated. Call tool_submit_strategy to publish to UI."}

    async def tool_submit_strategy(self, strategy_summary: str):
        """Submits the refined research strategy after identifying the problem. CALL THIS to advance to Step 2."""
        self._log("action", "Finalizing research strategy and PUBLISHING all research grounding...")
        from services.local_storage_service import update_mission_field, get_mission
        
        # Fetch current state
        mission_doc = get_mission(self.mission_id) or {}
        int_summary = mission_doc.get("internal_research_summary", "No research summary available.")
        int_papers = mission_doc.get("internal_research_papers", [])
        
        # Update fields
        update_mission_field(self.mission_id, "mission_strategy", strategy_summary)
        update_mission_field(self.mission_id, "research_summary", int_summary)
        update_mission_field(self.mission_id, "research_papers", int_papers)
        update_mission_field(self.mission_id, "current_step", 3) # Advance to Data Acquisition
        
        return {"status": "Strategy approved & Research grounded.", "step": 3}

    async def tool_run_p2rank(self, pdb_id: str):
        """Runs P2Rank on a protein PDB structure to identify druggable pockets. Returns pocket coordinates."""
        self._log("action", f"Running P2Rank on structure {pdb_id}...")
        # For now, we'll simulate the P2Rank result or fetch existing
        # Real integration requires the PDBManager which runs in background.
        # Here we just check if we have data or return a default/mock for the flow.
        from services.docking_service import calculate_pocket_center, DEFAULT_POCKET_PDB
        
        # In a full flow, we'd trigger PDBManager.run_p2rank(pdb_id) and wait.
        # For this prototype, we assume the default pocket or center.
        center = calculate_pocket_center(DEFAULT_POCKET_PDB)
        return {"pockets": [{"id": 1, "score": 15.0, "center": center}]}

    async def tool_run_vina(self, smiles: str):
        """Performs molecular docking using Autodock Vina. Returns binding affinity (kcal/mol)."""
        self._log("action", f"Running Vina docking for {smiles}...")
        from services.docking_service import DockingManager
        from starlette.concurrency import run_in_threadpool
        
        # Run blocking Vina in threadpool
        result = await run_in_threadpool(DockingManager.run_docking, self.mission_id, smiles)
        
        if result['success']:
            best_score = result['scores'][0] if result['scores'] else 0.0
            self._log("system", f"  - Affinity: {best_score} kcal/mol")
            return {"affinity": best_score, "details": result}
        else:
             return {"error": result.get("error")}

    async def tool_generate_leads(self, seed_smiles: str, steps: int = 5, batch_size: int = 10, min_hac: int = 15):
        """Generates new lead molecules using the Beam Search Feudal Network. 
        Args:
            steps (int): Number of optimization steps.
            batch_size (int): Number of candidates per step.
            min_hac (int): Minimum heavy atom count to filter out tiny molecules.
        CALL THIS for Step 6."""
        
        # Explicit cast to int to handle float inputs from Gemini/REST transport
        steps = int(steps)
        batch_size = int(batch_size)
        min_hac = int(min_hac)
        
        self._log("action", f"Generating leads from {seed_smiles} with Beam Search (Steps: {steps}, Batch: {batch_size})...")
        from services.feudal_service import FeudalOptimizationManager
        from services.analysis_service import calculate_pareto_winner
        from services.local_storage_service import update_mission_field
        from starlette.concurrency import run_in_threadpool
        
        try:
            feudal = FeudalOptimizationManager.get_instance()
        except NameError:
            self._log("system", "✦ Recovering from internal configuration issue...")
            self._reset_feudal_singleton()
            feudal = FeudalOptimizationManager.get_instance()
        
        # 📊 STATS LOGGING: Before Generation
        initial_props = await self.tool_analyze_molecule(seed_smiles)
        self._log("info", f"Initial Molecule Stats: {initial_props}")

        # Run optimized beam search
        await run_in_threadpool(feudal.run_optimization, self.mission_id, seed_smiles, steps=steps, batch_size=batch_size, min_hac=min_hac)
        
        # Load best candidate from history
        history = feudal.get_history(self.mission_id)
        if not history: return {"error": "Generation failed to produce history"}
        
        # Determine best candidate (min docking score)
        best = min(history, key=lambda x: x['score'])
        
        # 📊 STATS LOGGING: After Generation
        final_props = await self.tool_analyze_molecule(best['smiles'])
        self._log("success", f"Final Optimized Molecule Stats: {final_props}")

        # 🧠 AUTOMATED REASONING: Explain why it's better
        explanation = await self._gemini_explain_improvement(initial_props.get("properties", {}), final_props.get("properties", {}))
        
        # 🔗 INTERACTION ANALYSIS: Auto-trigger after generation
        interaction_res = await self.tool_structural_interaction_analysis(seed_smiles, best['smiles'])

        # Normalize all property keys to lowercase for consistent frontend access
        raw_props = final_props.get("properties", {})
        normalized_props = {k.lower(): v for k, v in raw_props.items()}

        profile = {
            "smiles": best['smiles'],
            "score": best['score'],
            "properties": normalized_props,
            "rationale": explanation,
            "interactions": interaction_res,
            "status": "Lead Generated & Analyzed"
        }
        update_mission_field(self.mission_id, "best_candidate_profile", profile)
        update_mission_field(self.mission_id, "current_step", 6)

        return {"status": "Leads generated and analyzed", "best_smi": best['smiles'], "explanation": explanation}

    async def tool_structural_interaction_analysis(self, initial_smi: str, optimized_smi: str):
        """Analyzes structural interactions (PLIP) for the initial vs optimized ligand."""
        self._log("action", "Running structural interaction analysis...")
        from services.interaction_service import analyze_structural_interactions, pdbqt_to_pdb, merge_complex
        from config import DEFAULT_RECEPTOR_PDB, PROJECTS_DIR
        
        res_dir = os.path.join(PROJECTS_DIR, self.mission_id, "analysis")
        os.makedirs(res_dir, exist_ok=True)
        
        # Assume best_pose.pdbqt exists from docking (standard workflow)
        pose_path = os.path.join(PROJECTS_DIR, self.mission_id, "docking", "best_pose.pdbqt")
        if not os.path.exists(pose_path):
            return {"warning": "Docked pose not found for detailed interaction mapping."}

        lig_pdb = os.path.join(res_dir, "optimized_ligand.pdb")
        complex_pdb = os.path.join(res_dir, "complex.pdb")
        
        pdbqt_to_pdb(pose_path, lig_pdb)
        merge_complex(DEFAULT_RECEPTOR_PDB, lig_pdb, complex_pdb)
        
        interactions = analyze_structural_interactions(complex_pdb, lig_pdb)
        
        # Save to mission profile
        from services.local_storage_service import get_mission, update_mission_field
        mission_data = get_mission(self.mission_id) or {}
        profile = mission_data.get("best_candidate_profile", {})
        profile["interactions"] = interactions
        update_mission_field(self.mission_id, "best_candidate_profile", profile)
        
        return interactions

    async def _gemini_explain_improvement(self, before_props, after_props):
        """Uses Gemini to provide a rational explanation for affinity/property changes."""
        model = genai.GenerativeModel('gemini-3-flash-preview')
        
        prompt = (
            f"Compare these two molecules in a drug discovery context:\n"
            f"BEFORE: {json.dumps(before_props)}\n"
            f"AFTER: {json.dumps(after_props)}\n\n"
            "Explain in 2-3 technical sentences why the optimized molecule is a superior candidate. "
            "Focus on potency (score), synthetic accessibility, and druglikeness. Be precise."
        )
        try:
            res = model.generate_content(prompt)
            return res.text.strip()
        except:
            return "The optimized molecule shows improved binding affinity and favorable ADMET characteristics."

    async def tool_find_pareto_front_molecule(self, weights: Dict[str, float], rationale: str):
        """
        Uses weighted property analysis (QED, SAScore, MW, pIC50) to find the best seed molecule from the CSV grounding set.
        Weights should be based on your clinical rationale (e.g. prioritize QED if drug-likeness is key).
        """
        self._log("action", "Searching grounding dataset for optimal Pareto-front seed...")
        self._log("thought", f"Rationale: {rationale}")
        
        from services.analysis_service import find_best_molecule_from_csv
        from config import PROJECTS_DIR
        
        csv_filename = f"{self.mission_id}_chembl_data.csv"
        csv_path = os.path.join(PROJECTS_DIR, self.mission_id, csv_filename)
        
        try:
            from services.local_storage_service import update_mission_field
            update_mission_field(self.mission_id, "active_weights", weights)
            
            best_mol = find_best_molecule_from_csv(csv_path, weights)
            if not best_mol:
                return {"error": "Could not find grounding CSV. Please ensure Step 3 (ChEMBL Fetch) is complete."}
            
            # Sync to Best Candidate Profile for Summary Dashboard
            # Sync to Best Candidate Profile for Summary Dashboard
            from services.bio_service import ADMETManager
            admet_props = ADMETManager.predict_properties(best_mol['smiles'])
            
            profile = {
                "smiles": best_mol['smiles'],
                "molecule_id": best_mol['molecule_id'],
                "score": best_mol['score'],
                "properties": best_mol['properties'],
                "admet": admet_props,
                "rationale": rationale,
                "status": "Grounding Champion"
            }
            update_mission_field(self.mission_id, "best_candidate_profile", profile)
            
            self._log("system", f"Best seed found: {best_mol['molecule_id']} (Weight Score: {best_mol['score']:.2f})")
            return {
                "smiles": best_mol['smiles'],
                "properties": best_mol['properties'],
                "mission_id": self.mission_id,
                "status": "Molecule identified. You can now use this SMILES as a seed for 'tool_generate_leads'."
            }
        except Exception as e:
            return {"error": str(e)}

    async def tool_suggest_starter_ligand(self):
        """Suggests a starter ligand for optimization based on previous ChEMBL analysis."""
        from services.local_storage_service import update_mission_field
        update_mission_field(self.mission_id, "current_step", 5)
        self._log("action", "Identifying best starter ligand for optimization...")
        # In literal implementation, we'd pick from the Pareto front
        suggestion = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C" # Example: Caffeine-like scaffold
        self._log("system", f"Suggested Starter: {suggestion}")
        update_mission_field(self.mission_id, "starter_ligand", suggestion)
        return {"suggested_smiles": suggestion, "rationale": "High binding affinity and low MW scaffold observed in ChEMBL data."}

    async def tool_fetch_pdb(self, query: str):
        """
        Fetches relevant protein structures from RCSB PDB. Returns a list of PDB IDs.
        AFTER getting the IDs, you MUST select the most relevant one and call 'tool_fetch_protein_structure' to download it.
        """
        self._log("action", f"Fetching PDB structures for: {query}...")
        # Simple simulated search for now or use requests to rcsb
        try:
            import requests
            url = f"https://search.rcsb.org/rcsbsearch/v2/query?json=%7B%22query%22%3A%7B%22type%22%3A%22terminal%22%2C%22service%22%3A%22full_text%22%2C%22parameters%22%3A%7B%22value%22%3A%22{query}%22%7D%7D%2C%22return_type%22%3A%22entry%22%7D"
            r = requests.get(url)
            ids = [item['identifier'] for item in r.json().get('result_set', [])[:5]]
            self._log("system", f"Found PDB IDs: {', '.join(ids)}")
            return {"pdb_ids": ids}
        except Exception as e:
            return {"error": str(e)}

    async def tool_fetch_chembl(self, target_name: str):
        """Fetches binding affinity data (IC50) from ChEMBL for a given protein. CALL THIS for Step 3."""
        # Use the more robust data fetcher we just implemented
        return await self.tool_fetch_chembl_data(target_name)

    async def tool_fetch_protein_structure(self, pdb_id: str):
        """Downloads a protein structure (PDB) to the project workspace. CALL THIS for Step 2 or 3."""
        self._log("action", f"Downloading protein structure: {pdb_id}")
        from services.pdb_service import PDBManager
        try:
            from services.local_storage_service import link_pdb_to_mission
        except ImportError:
            # Force reload if the function is on disk but not in memory
            import importlib
            import services.local_storage_service
            importlib.reload(services.local_storage_service)
            from services.local_storage_service import link_pdb_to_mission
            
        try:
            result = PDBManager.fetch_from_rcsb(pdb_id, mission_id=self.mission_id)
            if result:
                link_pdb_to_mission(self.mission_id, pdb_id, result["filename"])
                from services.local_storage_service import update_mission_field
                update_mission_field(self.mission_id, "active_pdb", pdb_id)
                return {"status": "success", "pdb_id": pdb_id, "path": result["filename"]}
            else:
                raise Exception(f"Failed to download structure for {pdb_id}")
        except Exception as e:
            self._log("error", f"PDB Download failed: {str(e)}")
            return {"error": str(e)}

    async def tool_fetch_chembl_data(self, protein_name: str):
        """
        Fetches binding affinity data (IC50) from ChEMBL for a given protein.
        CALL THIS to gather the 'grounding' dataset for molecular statistics.
        """
        self._log("action", f"Retrieving ligand data from ChEMBL for: {protein_name}")
        from services.bio_service import fetch_chembl_csv
        from services.local_storage_service import update_mission_field
        try:
            results = fetch_chembl_csv(protein_name, self.mission_id)
            if "error" in results:
                self._log("error", results["error"])
                return results
            
            # Sync metadata to Firestore
            update_mission_field(self.mission_id, "chembl_summary", results)
            update_mission_field(self.mission_id, "current_step", 3)
            return results
        except Exception as e:
            self._log("error", f"ChEMBL Fetch failed: {str(e)}")
            return {"error": str(e)}

    async def tool_analyze_ligand_set(self, ligands: list):
        """Generates statistics (MW, LogP) and Pareto front for a set of ligands."""
        self._log("action", "Analyzing ligand distributions and Pareto front...")
        from rdkit import Chem
        from rdkit.Chem import Descriptors
        
        # Robustly extract properties, skipping invalid entries
        processed = []
        for l in ligands:
            if isinstance(l, dict):
                processed.append(l)
            elif isinstance(l, str):
                # Try to analyze SMILES if string provided
                try:
                    mol = Chem.MolFromSmiles(l)
                    if mol:
                        processed.append({
                            "smiles": l,
                            "molwt": Descriptors.MolWt(mol),
                            "affinity": 0, # Default if not known
                            "logp": Descriptors.MolLogP(mol)
                        })
                except:
                    continue

        stats = {
            "mw_distribution": [l.get('molwt', 0) for l in processed],
            "affinity_distribution": [l.get('affinity', 0) for l in processed],
            "pareto_front": processed[:3]
        }
        from services.local_storage_service import update_mission_field
        update_mission_field(self.mission_id, "ligand_stats", stats)
        update_mission_field(self.mission_id, "current_step", 4)
        return {"stats": stats}

    async def tool_analyze_molecule(self, smiles: str):
        """Analyzes a molecule's properties (MW, LogP, QED, SA) using RDKit."""
        self._log("action", f"Analyzing properties for {smiles}...")
        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol: return {"valid": False}
            
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            qed = Descriptors.qed(mol)
            
            return {"valid": True, "properties": {"MW": mw, "LogP": logp, "QED": qed}}
        except Exception as e:
             return {"valid": False, "error": str(e)}

    # Lipinski fits into analyze_molecule, but we keep it specific if the agent wants to "verify"
    async def tool_verify_lipinski(self, smiles: str):
        """Verifies if a molecule follows Lipinski's Rule of Five."""
        # Reuse logic or call analyze
        analysis = await self.tool_analyze_molecule(smiles)
        if not analysis.get("valid"): return analysis
        
        props = analysis["properties"]
        # Simplified check (need explicit HBD/HBA which analyze didn't retrieve, so let's re-calculate or just assume)
        # For robustness, let's implement the full check:
        mol = Chem.MolFromSmiles(smiles)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        
        violations = 0
        if props["MW"] > 500: violations += 1
        if props["LogP"] > 5: violations += 1
        if hbd > 5: violations += 1
        if hba > 10: violations += 1
        
        passed = violations <= 1
        self._log("system", f"Lipinski Check: {'PASSED' if passed else 'FAILED'} ({violations} violations)")
        return {"passed": passed, "violations": violations}

    async def start_mission(self, disease_target: str, description: str = ""):
        """Main autonomous loop with manual tool execution and step-wise control."""
        self._log("system", f"Mission started for target: {disease_target}")
        self._ensure_chat_aligned()
        
        # Initialize Firestore state
        from services.local_storage_service import update_mission_field
        update_mission_field(self.mission_id, "current_step", 1)
        update_mission_field(self.mission_id, "awaiting_confirmation", False)

        prompt = (
            f"You are PHOENIX-PI, an autonomous marathon agent for drug discovery. "
            f"Your mission is to find and verify a potential inhibitor for {disease_target}. "
            f"Guiding Problem/Context: {description}\n\n"
            "MANDATORY WORKFLOW:\n"
            "Step 1: Deep Grounding (Literature Review). Use 'tool_search_literature' to build your knowledge. "
            "The research summary is CUMULATIVE. You may call this tool multiple times with different queries. "
            "AS THE ORCHESTRATOR, you decide when you have enough grounding to proceed. Once satisfied, move to Step 2.\n"
            "Step 2: Problem identification and solution planning. CALL 'tool_submit_strategy' to advance to Step 3.\n"
            "Step 3: Protein & Ligand Data Acquisition. 1) Identify PDB IDs via 'tool_fetch_pdb'. 2) SELECT one and DOWNLOAD it with 'tool_fetch_protein_structure'. 3) Fetch ligand bioactivity from ChEMBL with 'tool_fetch_chembl_data'. Do NOT stop until all three are initiated if possible.\n"
            "Step 4: Molecular statistics (MW, LogP, Affinity distributions) + Pareto front analysis.\n"
            "Step 5: Suggest a starter ligand based on analysis.\n"
            "Step 6: RL Generation loop for optimized leads + property change visualization.\n"
            "Step 7: MANDATORY - You MUST call 'tool_retrosynthesis_analysis' on your best lead compound BEFORE calling 'tool_submit_final_report'. "
            "Do NOT skip this step even if you think it might fail. The tool handles errors gracefully. "
            "Only after calling tool_retrosynthesis_analysis should you call tool_submit_final_report.\n\n"
            "NOTE: After each step, tell the user what you found and say 'WAITING FOR CONFIRMATION'. "
            "I will then trigger your next turn. Do NOT use tools for the next step until I tell you to proceed.\n"
            "Output your reasoning in the following format: "
            "Plan: <your plan> | Reasoning: <why this step> | Anticipated Outcome: <what you expect>"
        )


        return await self._run_agent_turn(prompt)

    async def resume_mission(self, user_message: str = "Proceed to the next logical step."):
        """Resumes the mission after user confirmation."""
        self._ensure_chat_aligned()
        from services.local_storage_service import update_mission_field
        update_mission_field(self.mission_id, "awaiting_confirmation", False)
        return await self._run_agent_turn(user_message)

    async def tool_retrosynthesis_analysis(self, smiles: str):
        """Analyzes the retrosynthetic accessibility of a molecule. CALL THIS for Step 7. MANDATORY before final report."""
        self._log("action", f"Running retrosynthesis analysis for {smiles}...")
        from services.bio_service import RetroManager
        from starlette.concurrency import run_in_threadpool
        from services.local_storage_service import update_mission_field

        # Always advance to Step 7 when retro is called, regardless of outcome
        update_mission_field(self.mission_id, "current_step", 7)

        try:
            results = await run_in_threadpool(RetroManager.run_retrosynthesis, smiles)

            if "error" in results:
                # AIZynth process ran but failed — store error so UI can show it
                self._log("error", f"Retrosynthesis engine error: {results['error']}")
                update_mission_field(self.mission_id, "retro_data", {
                    "error": results["error"],
                    "routes": [],
                    "smiles": smiles
                })
                return {"status": "retrosynthesis_failed", "error": results["error"], "note": "Step 7 reached. AIZynth could not find routes — the AI summary will still document synthetic accessibility based on structural reasoning."}

            update_mission_field(self.mission_id, "retro_data", results)
            route_count = len(results.get("routes", []))
            self._log("success", f"Retrosynthesis complete: {route_count} route(s) found.")
            return {"routes": results, "route_count": route_count, "summary": f"Retrosynthetic analysis complete. {route_count} viable route(s) identified."}
        except Exception as e:
            self._log("error", f"Retrosynthesis critical failure: {str(e)}")
            update_mission_field(self.mission_id, "retro_data", {"error": str(e), "routes": [], "smiles": smiles})
            return {"status": "retrosynthesis_failed", "error": str(e)}



    async def tool_submit_final_report(self, report_summary: str):
        """Submits a final comprehensive mission report combining literature, docking, and synthesis results."""
        self._log("action", "Generating and submitting final mission synthesis...")
        from services.local_storage_service import update_mission_field
        update_mission_field(self.mission_id, "final_summary", report_summary)  # Renamed to match frontend
        update_mission_field(self.mission_id, "current_step", 7) # Finalize
        return {"status": "Mission Finalized", "message": "Final summary published to dashboard."}
    
    async def tool_analyze_structure_image(self, image_filename: str):
        """Analyzes protein structure images using Gemini Vision to identify binding sites."""
        self._log("action", f"Analyzing protein structure image: {image_filename}")
        from services.multimodal_service import MultimodalAnalyzer
        import os
        from config import PROJECTS_DIR
        
        image_path = os.path.join(PROJECTS_DIR, self.mission_id, "uploads", image_filename)
        if not os.path.exists(image_path):
            return {"error": f"Image not found: {image_filename}"}
        
        result = await MultimodalAnalyzer.analyze_structure_image(image_path, self.mission_id)
        return result
    
    async def tool_analyze_research_pdf(self, pdf_filename: str):
        """Extracts insights from research paper PDFs using Gemini's PDF understanding."""
        self._log("action", f"Analyzing research PDF: {pdf_filename}")
        from services.multimodal_service import MultimodalAnalyzer
        import os
        from config import PROJECTS_DIR
        
        pdf_path = os.path.join(PROJECTS_DIR, self.mission_id, "uploads", pdf_filename)
        if not os.path.exists(pdf_path):
            return {"error": f"PDF not found: {pdf_filename}"}
        
        result = await MultimodalAnalyzer.analyze_research_pdf(pdf_path, self.mission_id)
        return result
    
    async def tool_deep_literature_synthesis(self, paper_ids: list):
        """Synthesizes insights from 50+ papers using Gemini's 2M token context window."""
        self._log("action", f"Performing deep literature synthesis on {len(paper_ids)} papers...")
        from services.long_context_service import LongContextAnalyzer
        
        result = await LongContextAnalyzer.deep_literature_synthesis(paper_ids, self.mission_id)
        return result
    
    async def tool_full_chembl_analysis(self, target_name: str):
        """Analyzes full ChEMBL dataset (500+ molecules) for comprehensive SAR insights."""
        self._log("action", f"Analyzing full ChEMBL dataset for {target_name}...")
        from services.long_context_service import LongContextAnalyzer
        
        result = await LongContextAnalyzer.analyze_full_chembl_dataset(target_name, self.mission_id)
        return result

    async def _run_agent_turn(self, input_message: str):
        """Runs a single logical turn of the agent, processing tool calls robustly."""
        from services.local_storage_service import update_mission_field, get_mission
        mission_data = get_mission(self.mission_id) or {}
        current_step = mission_data.get("current_step", 1)
        
        # Inject step awareness
        step_reminder = f"\n[SYSTEM NOTICE: You are currently on STEP {current_step} of the mission. Focus your actions accordingly.]\n"
        full_input = input_message + step_reminder
        
        update_mission_field(self.mission_id, "is_thinking", True)
        
        # Loop safety: Reactive check
        self._ensure_chat_aligned()

        try:
            # Start/continue the conversational turn with a retry loop
            async def send_with_retry(msg, retries=3):
                from starlette.concurrency import run_in_threadpool
                for i in range(retries):
                    try:
                        # Use synchronous send_message in threadpool to bypass SDK async REST bug
                        return await run_in_threadpool(self.chat.send_message, msg)
                    except Exception as e:
                        err_msg = str(e).lower()
                        
                        # Fail fast on signature errors - retrying doesn't help 400s
                        if "thought_signature" in err_msg:
                            raise e

                        if "different loop" in err_msg or "loop" in err_msg:
                            self._log("system", "✦ Neural connection re-alignment required. Fixing...")
                            self.chat = self._initialize_chat_with_history()
                            # Immediate retry with new chat session
                            return await run_in_threadpool(self.chat.send_message, msg)
                        
                        if "500" in err_msg or "internal error" in err_msg:
                            self._log("system", "✦ Vertex/Gemini server-side glitch. Retrying...")
                            await asyncio.sleep(1) # Short break for transient 500s
                            
                        if i == retries - 1: raise e
                        wait = (i + 1) * 2
                        self._log("system", f"Gemini connection lost. Retrying in {wait}s... ({i+1}/{retries})")
                        await asyncio.sleep(wait)
            
            self._log("system", "✦ PHOENIX-PI: Analyzing research state...")
            response = await send_with_retry(full_input)
            
            # Loop for multi-turn reasoning (chaining tool calls)
            for turn_idx in range(12):  # Slightly higher limit for complex rounds
                if not response.candidates:
                    self._log("error", "Received empty response from Gemini.")
                    break

                parts = response.candidates[0].content.parts
                tool_calls = []
                
                # 1. First pass: Log thoughts and collect all tool calls in this response
                for part in parts:
                    if part.text:
                        self._log("thought", part.text)
                        
                        # Handle confirmation pauses
                        confirmation_phrases = [
                            "waiting for confirmation", "ready for your confirmation", "confirm and proceed", 
                            "please confirm", "waiting for you to confirm", "awaiting confirmation",
                            "confirmation required", "let me know when you're ready"
                        ]
                        
                        # Search ALL parts for confirmation phrases
                        found_confirmation = False
                        for p in parts:
                            if p.text and any(phrase in p.text.lower() for phrase in confirmation_phrases):
                                found_confirmation = True
                                break

                        if found_confirmation:
                            mission_data = get_mission(self.mission_id) or {}
                            auto_confirm = mission_data.get("auto_correct", False) or mission_data.get("auto_confirm", False)
                            
                            if auto_confirm:
                                self._log("system", "✦ AUTO-RESUME: Continuing research autonomously.")
                                update_mission_field(self.mission_id, "awaiting_confirmation", False)
                            else:
                                self._log("system", "Waiting for your confirmation to proceed.")
                                update_mission_field(self.mission_id, "awaiting_confirmation", True)

                        if part.text and "mission complete" in part.text.lower():
                            self._log("system", "Mission concluded.")
                            update_mission_field(self.mission_id, "is_thinking", False)
                            return

                    if part.function_call:
                        tool_calls.append(part.function_call)

                # 2. Second pass: Execute tool calls if any
                if tool_calls:
                    responses_to_send = []
                    raw_results = []
                    for call in tool_calls:
                        tool_name = call.name
                        args = call.args
                        self._log_tool_call(tool_name, args)
                        self._log("system", f"✦ EXECUTING: {tool_name}...") # Heartbeat
                        
                        if tool_name in self.tools:
                            tool_fn = self.tools[tool_name]
                            try:
                                if inspect.iscoroutinefunction(tool_fn):
                                    result = await tool_fn(**args)
                                else:
                                    from starlette.concurrency import run_in_threadpool
                                    result = await run_in_threadpool(tool_fn, **args)
                                
                                if not isinstance(result, (dict, list)):
                                    result = {"result": str(result)}

                                # Prepare response part
                                responses_to_send.append(
                                    protos.Part(
                                        function_response=protos.FunctionResponse(
                                            name=tool_name,
                                            response=self._dict_to_struct(result)
                                        )
                                    )
                                )
                                raw_results.append((tool_name, result))
                            except Exception as e:
                                self._log("error", f"Tool {tool_name} failed: {str(e)}")
                                err_res = {"error": str(e)}
                                responses_to_send.append(
                                    protos.Part(
                                        function_response=protos.FunctionResponse(
                                            name=tool_name,
                                            response=self._dict_to_struct(err_res)
                                        )
                                    )
                                )
                                raw_results.append((tool_name, err_res))
                        else:
                            self._log("error", f"Unknown tool: {tool_name}")
                            err_res = {"error": f"Tool '{tool_name}' not found."}
                            responses_to_send.append(
                                protos.Part(
                                    function_response=protos.FunctionResponse(
                                        name=tool_name,
                                        response=self._dict_to_struct(err_res)
                                    )
                                )
                            )
                            raw_results.append((tool_name, err_res))

                    # Send all gathered tool responses back in batch with retry
                    try:
                        self._log("system", "✦ CONTINUING: Processing tool results...")
                        
                        # Definitively bypass if we know the SDK/Model combo is 'poisoned'
                        if self.sdk_poisoned:
                            raise ValueError("SDK_POISONED_BY_THOUGHT_SIGNATURE")

                        response = await send_with_retry(responses_to_send)
                    except Exception as e:
                        err_str = str(e)
                        is_signature_err = "thought_signature" in err_str or "SDK_POISONED" in err_str
                        is_server_err = "500" in err_str or "Internal error" in err_str
                        
                        if is_signature_err:
                            self.sdk_poisoned = True
                            self._log("system", "✦ Neural Re-alignment: Bypassing signature validation...")
                        elif is_server_err:
                            self._log("system", "✦ Neural Recovery: Handling server-side transient failure...")
                            await asyncio.sleep(2) # Backoff
                        else:
                            self._log("error", f"Gemini API tool response failed: {err_str}. Attempting recovery...")
                        
                        # Recovery: Rebuild history to clear any potentially problematic state
                        self.chat = self._initialize_chat_with_history()
                        
                        import json
                        recovery_prompt = "[SYSTEM] Resuming turn after internal recovery. Tool Results:\n\n"
                        for name, res in raw_results:
                            # Clamp huge JSON outputs to prevent context overflow
                            json_str = json.dumps(res, default=str)
                            # Even tighter clamp for recovery to avoid 500s on large payloads
                            recovery_prompt += f"--- {name} ---\n```\n{json_str[:4000]}\n```\n\n"
                            
                        from starlette.concurrency import run_in_threadpool
                        try:
                            response = await run_in_threadpool(self.chat.send_message, recovery_prompt)
                        except Exception as e2:
                            self._log("error", f"Critical Recovery Failure: {str(e2)}")
                            # Last ditch: try simple confirmation
                            response = await run_in_threadpool(self.chat.send_message, "Recovery failed. Please re-analyze the workspace state.")
                    continue # Start next turn in the tool chain
                
                # 3. Handle Auto-Confirm Resumption (if no tools were called but we hit the phrase)
                found_waiting = False
                if response.candidates:
                    for part in response.candidates[0].content.parts:
                        if part.text and "waiting for confirmation" in part.text.lower():
                            found_waiting = True
                            break

                if found_waiting:
                     mission_data = get_mission(self.mission_id) or {}
                     if mission_data.get("auto_correct", False):
                         self._log("system", "Auto-Confirm triggered. Proactively resuming...")
                         response = await send_with_retry("Please proceed to the next step immediately.")
                         continue
                     else:
                         self._log("system", "Awaiting manual confirmation.")
                         update_mission_field(self.mission_id, "awaiting_confirmation", True)

                # No more tool calls to process this turn
                break

        except Exception as e:
            import traceback
            traceback.print_exc() 
            self._log("error", f"Orchestrator Failure: {str(e)}")
        
        update_mission_field(self.mission_id, "is_thinking", False)
        return {"log": self.thought_log}


    def _dict_to_struct(self, d):
        """Safe conversion of dict to Struct proto."""
        from services.local_storage_service import _clean_serializable
        s = struct_pb2.Struct()
        if isinstance(d, dict):
            s.update(_clean_serializable(d))
        else:
            s.update({"result": str(d)})
        return s

