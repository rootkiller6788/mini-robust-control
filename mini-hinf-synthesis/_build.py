
import os
BASE = os.path.dirname(os.path.abspath(__file__))

def wf(relpath, content):
    path = os.path.join(BASE, relpath)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w', encoding='utf-8', newline='
') as f:
        f.write(content)
    n = content.count('
')
    print(f'  wrote {relpath} ({n} lines)')

# Add the hinf_math.h we already wrote
# Now write the implementation files
