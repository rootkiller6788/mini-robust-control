"""Build script for mini-structured-singular-value ¡ª generates all missing files."""
import os

BASE = r"F:\nano-everything\mini-complex-control-theory\13. mini-robust-control\mini-structured-singular-value"

def write_file(relpath, content):
    path = os.path.join(BASE, relpath)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w', encoding='utf-8') as f:
        f.write(content)
    lines = content.count('\n')
    print(f"  wrote {relpath} ({lines} lines)")

print("Build script loaded. BASE =", BASE)
