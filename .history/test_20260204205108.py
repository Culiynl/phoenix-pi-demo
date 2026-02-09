import subprocess

pdb_path = r"F:\Coding Stuff\Full Model 2026 Fixed\webserver_stuff\backend\data\proteins\6OD6_protein.pdb"
output_dir = r"F:\Coding Stuff\Full Model 2026 Fixed\webserver_stuff\backend\data\proteins\p2rank_results\1"
prank_bat = r"F:\Coding Stuff\Full Model 2026 Fixed\webserver_stuff\backend\p2rank_2.5.1\prank.bat"

cmd = f'cmd /c ""{prank_bat}" predict -f "{pdb_path}" -o "{output_dir}" -v""'

process = subprocess.Popen(
    cmd,
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
    text=True,
    shell=True  # required when passing a single command string
)

for line in iter(process.stdout.readline, ''):
    print(line.strip())

process.stdout.close()
process.wait()
