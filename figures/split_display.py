"""Split figures/latex_includes.tex into paper/display/{fig1..fig6,table1}.tex (one display item per file)."""
import re
from pathlib import Path
root = Path(__file__).resolve().parents[1] / "paper"
src = (root.parent / "figures" / "latex_includes.tex").read_text()
blocks = re.split(r"(?m)^(?=% === )", src)
out = root / "display"; out.mkdir(exist_ok=True)
for b in blocks:
    m = re.match(r"% === (Fig|Table) (\d+)", b)
    if m:
        name = f"{m.group(1).lower()}{m.group(2)}.tex"
        (out / name).write_text(b.rstrip() + "\n"); print(name)
