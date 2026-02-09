import subprocess
import os

def get_short_path(long_path):
    import ctypes
    buf = ctypes.create_unicode_buffer(260)
    ctypes.windll.kernel32.GetShortPathNameW(long_path, buf, 260)
    return buf.value

pdb_path = get_short_path(r"F:\Coding Stuff\Full Model 2026 Fixed\webserver_stuff\backend\data\proteins\6od6_protein.pdb")
output_dir = get_short_path(r"F:\p2rank_2.5.1")
prank_bat = get_short_path(r"F:\p2rank_2.5.1\prank.bat")

cmd = f'cmd /c "{prank_bat} predict -f {pdb_path} -o {output_dir}"'

print(cmd)
subprocess.run(cmd, shell=True)
