import os
BASE = r"F:\nano-everything\mini-complex-control-theory\13. mini-robust-control\mini-parametric-uncertainty"
def wf(rel, content):
    p = os.path.join(BASE, rel)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p, "w", encoding="utf-8") as f:
        f.write(content)
    lc = content.count("\n")
    print(f"  OK {rel}: {lc} lines")
    return lc
total = 0
print("Writing source files...")
