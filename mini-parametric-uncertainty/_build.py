import os, sys
BASE = os.path.dirname(os.path.abspath(__file__))
def wf(relpath, content):
    full = os.path.join(BASE, relpath)
    os.makedirs(os.path.dirname(full), exist_ok=True)
    with open(full, "w", encoding="utf-8") as fh:
        fh.write(content)
    lc = content.count("
")
    print(f"  {relpath}: {lc} lines")
    return lc
total = 0
print("Building files via _build.py ...")
