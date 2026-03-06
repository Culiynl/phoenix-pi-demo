import os
import json
import subprocess
import tempfile
import textwrap

# Configuration for the environment
AIZYNTH_ENV_PYTHON = r"F:\PyMol\envs\aizynth-env\python.exe"
AIZYNTH_CONFIG = "config.yml" 
IBM_RXN_API_KEY = os.getenv("IBM_RXN_API_KEY", "apk-bba11db0846bc8203f3fe5965c8d2b20e1fdd96a57f1ef44a956e81a293fec9d")

class RXNManager:
    _instance = None
    
    @classmethod
    def get_instance(cls):
        if cls._instance is None:
            cls._instance = RXNManager()
        return cls._instance
        
    def __init__(self):
        # We store settings here, but we initialize the actual 
        # heavy objects (AiZynthApp/RXN) inside the worker process
        self.python_exe = AIZYNTH_ENV_PYTHON
        self.config_path = os.path.abspath(AIZYNTH_CONFIG)
        self.api_key = IBM_RXN_API_KEY

    def _execute_worker(self, code_block):
        """Helper to run code in the aizynth-env and capture output."""
        with tempfile.NamedTemporaryFile(delete=False, suffix=".py", mode="w", encoding="utf-8") as temp:
            temp.write(code_block)
            temp_path = temp.name
        
        try:
            # Set environment variables for the worker
            env = os.environ.copy()
            env["IBM_RXN_API_KEY"] = self.api_key
            
            result = subprocess.run(
                [self.python_exe, temp_path],
                capture_output=True, text=True, check=True, env=env
            )
            return result.stdout
        except subprocess.CalledProcessError as e:
            print(f"Error in aizynth-env worker: {e.stderr}")
            return None
        finally:
            if os.path.exists(temp_path):
                os.remove(temp_path)

    def predict_reaction(self, reactants_smiles):
        """Calls IBM RXN predict_reaction."""
        worker_code = textwrap.dedent(f"""
            import json
            import os
            from rxn4chemistry import RXN4ChemistryWrapper
            
            api_key = os.getenv("IBM_RXN_API_KEY")
            try:
                rxn = RXN4ChemistryWrapper(api_key=api_key)
                # Replicating your original call
                result = rxn.predict_reaction("{reactants_smiles}")
                print(json.dumps(result))
            except Exception as e:
                print(json.dumps({{"error": str(e)}}))
        """)
        output = self._execute_worker(worker_code)
        return json.loads(output) if output else {{"error": "Execution failed"}}

    def run_retrosynthesis(self, target_smiles):
        """
        Calls AiZynthApp exactly as specified:
        target_smiles -> search_tree() -> build_routes() -> routes.to_dict()
        """
        worker_code = textwrap.dedent(f"""
            import json
            from aizynthfinder.api import AiZynthApp
            
            try:
                az_finder = AiZynthApp(r"{self.config_path}")
                az_finder.target_smiles = "{target_smiles}"
                az_finder.search_tree()
                az_finder.build_routes()
                
                # Returns standard list of dicts
                print(json.dumps(az_finder.routes.to_dict()))
            except Exception as e:
                print(json.dumps({{"error": str(e)}}))
        """)
        output = self._execute_worker(worker_code)
        return json.loads(output) if output else {{"error": "Retrosynthesis failed"}}

    def get_reaction_name(self, reaction_smiles):
        """Heuristics to name a reaction (Local string processing)."""
        common_reactions = [
            ("Suzuki Coupling", ["[c:1][B:2]", "[c:1][B](O)O"]),
            ("Amide Coupling", ["C(=O)N", "NC(=O)"]),
            ("Esterification", ["C(=O)O[C]"]),
            ("Reductive Amination", ["[C][NH][C]", "[C][N]([C])[C]"]),
            ("Grignard Reaction", ["[C][Mg]"]),
            ("Friedel-Crafts", ["c[C]=O"]),
            ("Boc Protection", ["NC(=O)OC(C)(C)C"]),
            ("Deprotection", ["[NH2]"])
        ]
        for name, patterns in common_reactions:
            for p in patterns:
                if p in reaction_smiles:
                    return name
        return "Synthetic Step"

    def generate_mermaid_route(self, route_data):
        """
        Generates Mermaid code. 
        Note: This is run in the worker to allow use of ReactionTree objects.
        """
        route_json = json.dumps(route_data)
        worker_code = textwrap.dedent(f"""
            import json
            from aizynthfinder.reactiontree import ReactionTree
            
            # Helper inside worker
            def get_rxn_name(smiles):
                patterns = [("Suzuki Coupling", ["[c:1][B:2]", "[c:1][B](O)O"]), ("Amide Coupling", ["C(=O)N", "NC(=O)"])]
                for name, ps in patterns:
                    if any(p in smiles for p in ps): return name
                return "Synthetic Step"

            route_data = json.loads('''{route_json}''')
            actual_route = route_data[0] if isinstance(route_data, list) else route_data
            
            try:
                tree = ReactionTree.from_dict(actual_route)
                lines = ["graph LR"]
                
                def traverse(node, child_id=None):
                    node_id = f"mol_{{abs(hash(node.smiles)) % 10000}}"
                    label = (node.smiles[:15] + "...") if len(node.smiles) > 15 else node.smiles
                    color = "#99ff99" if node.in_stock else "#ffcc99"
                    lines.append(f"  {{node_id}}[\"{{label}}\"]")
                    lines.append(f"  style {{node_id}} fill:{{color}}")
                    if child_id:
                        lines.append(f"  {{node_id}} -->|{{get_rxn_name(node.smiles)}}| {{child_id}}")
                    for reaction in node.reactions:
                        for precursor in reaction.precursors:
                            traverse(precursor, node_id)
                
                traverse(tree.root)
                print("\\n".join(lines))
            except Exception as e:
                print(f"graph TD\\n  Error[\\"{{str(e)}}\\"]")
        """)
        output = self._execute_worker(worker_code)
        return output if output else "graph TD\n  Error[\"Worker Failure\"]"

    def save_route_image(self, route_data, filepath):
        """Calls ReactionTree.from_dict -> tree.to_image -> img.save()."""
        route_json = json.dumps(route_data)
        abs_filepath = os.path.abspath(filepath)
        
        worker_code = textwrap.dedent(f"""
            import json
            from aizynthfinder.reactiontree import ReactionTree
            
            route_data = json.loads('''{route_json}''')
            actual_route = route_data[0] if isinstance(route_data, list) else route_data
            
            tree = ReactionTree.from_dict(actual_route)
            img = tree.to_image(in_stock_colors={{True: "#99ff99", False: "#ffcc99"}})
            img.save(r"{abs_filepath}")
        """)
        self._execute_worker(worker_code)
        return filepath