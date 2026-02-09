import subprocess
import os

pdb_path = r"F:\Coding Stuff\Full Model 2026 Fixed\webserver_stuff\backend\data\proteins\6OD6.pdb"
output_dir = r"F:\p2rank_2.5.1\"
prank_bat = r"F:\p2rank_2.5.1\prank.bat"

# Make sure output folder exists
os.makedirs(output_dir, exist_ok=True)

# Wrap paths in extra quotes for Windows batch/Java
cmd = f'cmd /c ""{prank_bat}" predict -f "{pdb_path}" -o "{output_dir}" -v""'

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
