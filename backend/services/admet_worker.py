import sys
import json
import argparse
import contextlib

# Try imports inside the script to avoid breaking if run in wrong env
try:
    from admet_ai import ADMETModel
except ImportError:
    print(json.dumps({"error": "ADMET-AI not installed in this environment"}))
    sys.exit(1)

def predict(smiles):
    try:
        # Redirect ALL internal library prints to stderr
        with contextlib.redirect_stdout(sys.stderr):
            model = ADMETModel()
            preds = model.predict(smiles=smiles)
            
        # The ONLY thing on stdout should be our JSON
        print(json.dumps(preds), flush=True)
    except Exception as e:
        # Error still goes to stdout as JSON for the caller to parse
        print(json.dumps({"error": str(e)}), flush=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--smiles", required=True, help="SMILES string to analyze")
    args = parser.parse_args()
    
    predict(args.smiles)
