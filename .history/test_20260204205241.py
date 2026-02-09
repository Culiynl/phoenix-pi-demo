import subprocess
import os

pdb_path = r"F:\Coding Stuff\Full Model 2026 Fixed\webserver_stuff\backend\data\proteins\6od6_protein.pdb"
output_dir = r"F:\Coding Stuff\Full Model 2026 Fixed\webserver_stuff\backend\data\proteins\p2rank_results\1"
prank_bat = r"F:\Coding Stuff\Full Model 2026 Fixed\webserver_stuff\backend\p2rank_2.5.1\prank.bat"

# Ensure output folder exists
os.makedirs(output_dir, exist_ok=True)

# Proper quoting
cmd = f'cmd /c ""{prank_bat}" predict -f "{pdb_path}" -o "{output_dir}""'

print(f"Running: {cmd}")

process = subprocess.Popen(
    cmd,
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
    text=True,
    shell=True
)

for line in iter(process.stdout.readline, ''):
    print(line.strip())

process.stdout.close()
process.wait()
print(f"Return code: {process.returncode}")
