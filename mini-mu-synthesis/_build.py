import os, sys
BASE = r"F:\nano-everything\mini-complex-control-theory\13. mini-robust-control\mini-mu-synthesis"

def wf(relpath, content):
    path = os.path.join(BASE, relpath)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write(content)
    lines = content.count("\n")
    print(f"  {relpath}: {lines} lines")

print("Building mini-mu-synthesis files...")
