"""Fig 1: event-construction schematic (editable SVG following FigEdit layer/ID conventions -> PDF/PNG via cairosvg).
Run with the Python environment installed from requirements.txt. Inset numbers are read from results/pw_p1_joint_icc.csv (floor 5, joint6),
results/rv_v7_ceiling_icc.csv (1,000 deg/s ceiling, same cohort) and results/rv_v4_loo_summary.csv (one participant omitted)."""
import re
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
J = pd.read_csv(ROOT / "results" / "pw_p1_joint_icc.csv")
J = J[(J.floor == 5) & (J.comparison == "R1_S1S2") & (J.cohort == "joint6")].set_index("method")
N = int(J.loc["mean_u", "n"]); assert N == int(J.loc["adm_u", "n"]) == 376
f3 = lambda x: f"{x:.3f}".replace("0.", ".", 1)
ICC_U, ICC_ADM = f3(J.loc["mean_u", "icc_a1"]), f3(J.loc["adm_u", "icc_a1"])
CE = pd.read_csv(ROOT / "results" / "rv_v7_ceiling_icc.csv").query("tag == 'floor5_ceil1000'").set_index("method")
assert int(CE.loc["ceil1000_u", "n"]) == N and abs(CE.loc["adm_u", "icc_a1"] - J.loc["adm_u", "icc_a1"]) < 1e-9
ICC_CE = f3(CE.loc["ceil1000_u", "icc_a1"])
LO = pd.read_csv(ROOT / "results" / "rv_v4_loo_summary.csv").query("floor == 5").iloc[0]
assert int(LO.n) == N
WO_U, WO_ADM = f3(LO.icc_mean_u_wo_top), f3(LO.icc_adm_u_wo_top)

W, H = 1600, 1075  # printed at \textwidth = 160 mm, so 1 px = 0.1 mm and 28 px = 8 pt
OR, BL, GR, INK, GREY = "#D55E00", "#0072B2", "#009E73", "#222222", "#8C8C8C"
FONT = "DejaVu Sans, Arial, Helvetica, sans-serif"
L = {k: [] for k in ("panels", "connectors", "sections", "icons", "texts", "annotations")}


FD = "/usr/share/fonts/truetype/dejavu/DejaVuSans"
_fonts = {}


def _w(seg, size, bold, ital):
    """Advance width of seg in DejaVu Sans (used to place italic-symbol runs; cairosvg ignores tspan flow)."""
    from PIL import ImageFont
    f = FD + {(0, 0): "", (1, 0): "-Bold", (0, 1): "-Oblique", (1, 1): "-BoldOblique"}[(bold, ital)] + ".ttf"
    key = (f, size)
    if key not in _fonts: _fonts[key] = ImageFont.truetype(f if Path(f).exists() else FD + ".ttf", size)
    return _fonts[key].getlength(seg)


def text(id_, x, y, s, size=22, fill=INK, anchor="start", weight="normal", style="normal", layer="texts"):
    esc = lambda t: t.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    attrs = f'font-family="{FONT}" font-size="{size}" fill="{fill}" font-weight="{weight}"'
    if "*" not in s:
        L[layer].append(f'<text id="{id_}" x="{x}" y="{y}" {attrs} text-anchor="{anchor}" font-style="{style}">{esc(s)}</text>')
        return
    segs = [(t, i % 2 == 1) for i, t in enumerate(re.split(r"\*", s)) if t]  # *l* -> italic run
    bold = weight == "bold"
    ws = [_w(t, size, bold, it or style == "italic") for t, it in segs]
    x0 = x - {"start": 0, "middle": sum(ws) / 2, "end": sum(ws)}[anchor]
    parts = []
    for k, ((t, it), w) in enumerate(zip(segs, ws)):
        st = "italic" if it else style
        parts.append(f'<text id="{id_}-{k}" x="{x0:.1f}" y="{y}" {attrs} text-anchor="start" font-style="{st}" xml:space="preserve">{esc(t)}</text>')
        x0 += w
    L[layer].append(f'<g id="{id_}">' + "".join(parts) + "</g>")


def rect(id_, x, y, w, h, stroke=INK, fill="#FFFFFF", sw=2.5, rx=10, dash=None, layer="sections"):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    L[layer].append(f'<rect id="{id_}" x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"{d}/>')


def line(id_, pts, color=INK, sw=2.5, arrow=True, dash=None, layer="connectors"):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    m = f' marker-end="url(#arr-{color.strip("#")})"' if arrow else ""
    p = " ".join(f"{x},{y}" for x, y in pts)
    L[layer].append(f'<polyline id="{id_}" points="{p}" fill="none" stroke="{color}" stroke-width="{sw}" stroke-linejoin="round"{d}{m}/>')


SMALL, BODY, TITLE = 28, 29, 32
FIT = []  # (id, text width, available width) checked at the end


def wrap(s, size, width, bold=False):
    out, cur = [], ""
    for wd in s.split():
        t = (cur + " " + wd).strip()
        if _w(t.replace("*", ""), size, bold, 0) <= width: cur = t
        else: out.append(cur); cur = wd
    return out + [cur]


def box(id_, x, y, w, color, title, body, fill, dash=None, size=BODY, lh=36):
    lines = body if isinstance(body, list) else wrap(body, size, w - 36)
    for t in lines: FIT.append((id_, _w(t.replace('*', ''), size, 0, 0), w - 20))
    h = 50 + lh * len(lines) + 18
    rect(f"section-{id_}", x, y, w, h, stroke=color, fill=fill, dash=dash)
    FIT.append((id_, _w(title, TITLE, 1, 0), w - 20))
    text(f"label-{id_}-title", x + w / 2, y + 40, title, size=TITLE, fill=color, anchor="middle", weight="bold")
    for i, t in enumerate(lines):
        text(f"label-{id_}-{i + 1}", x + w / 2, y + 50 + lh * (i + 1) - 6, t, size=size, anchor="middle")
    return h


def fit(id_, s, size, width, bold=False):
    FIT.append((id_, _w(s.replace("*", ""), size, bold, 0), width))


# ---- stage headings (row 1)
for i, (x, s) in enumerate([(20, "Gaze signals"), (340, "Event construction")]):
    text(f"label-stage-{i + 1}", x, 38, s, size=TITLE, fill=GREY, weight="bold")

# ---- (1) gaze traces and per-eye speed
C1, C1W = 20, 280
rect("panel-traces", C1, 62, C1W, 290, stroke=GREY, sw=1.5, layer="panels")
rect("panel-speed", C1, 412, C1W, 290, stroke=GREY, sw=1.5, layer="panels")
text("label-trace-axis", C1 + 14, 98, "horizontal gaze", size=SMALL, fill=GREY); fit("trace-axis", "horizontal gaze", SMALL, C1W - 28)
X0, Y0, STEP = C1 + 22, 236, 29
xs = [X0]; ys = [Y0]; sac = []
for k in range(6):  # four progressive steps, a return sweep, then a step
    x0 = xs[-1] + STEP
    if k == 4:
        xs += [x0, x0 + 12]; ys += [ys[-1], Y0]; sac.append((x0 + 6, "rs"))
    else:
        xs += [x0, x0 + 5]; ys += [ys[-1], ys[-1] - 26]; sac.append((x0 + 3, "prog"))
xs.append(xs[-1] + STEP); ys.append(ys[-1])
for eye, col, dy in (("left", INK, 0), ("right", GREY, 9)):
    line(f"trace-{eye}", [(x, y + dy) for x, y in zip(xs, ys)], color=col, sw=3, arrow=False, layer="icons")
for i, (eye, col) in enumerate((("left eye", INK), ("right eye", GREY))):
    kx, ky = C1 + 14, 290 + 32 * i
    L["icons"].append(f'<line id="key-{i}" x1="{kx}" y1="{ky - 9}" x2="{kx + 26}" y2="{ky - 9}" stroke="{col}" stroke-width="3"/>')
    text(f"label-key-{i}", kx + 34, ky, eye, size=SMALL - 2, fill=col)
text("label-speed-axis", C1 + 14, 448, "per-eye speed", size=SMALL, fill=GREY); fit("speed-axis", "per-eye speed", SMALL, C1W - 28)
base = 680
for eye, col, scale, off in (("left", INK, 1.0, 0), ("right", GREY, 0.86, 4)):
    pts = [(C1 + 12, base)]
    for x, kind in sac:
        h = 185 if kind == "rs" else 66
        pts += [(x - 8 + off, base), (x + off, base - h * scale), (x + 8 + off, base)]
    pts.append((C1 + C1W - 12, base))
    line(f"speed-{eye}", pts, color=col, sw=2.5, arrow=False, layer="icons")
text("label-trace-time", C1 + C1W - 12, 708 - 12, "", size=SMALL, fill=GREY, anchor="end")
line("arrow-traces-to-speed", [(C1 + C1W / 2, 354), (C1 + C1W / 2, 408)])

# ---- (2) detectors, gate, event sets
BX, BW = 340, 520
y = 62
hm = box("mean-speed", BX, y, BW, OR, "Mean-speed detector", "saccades detected on the mean of the two per-eye speed magnitudes", "#FBEFE8")
ym = y + hm / 2
yg = y + hm + 46
hg = box("gate", BX, yg, BW, BL, "Admission gate", "keep a mean-speed event only if it overlaps a matched left/right pair; "
         "windows, peaks, classes unchanged", "#E8F1F8")
yp = yg + hg + 46
hp = box("per-eye", BX, yp, BW, GR, "Per-eye detector", "saccades detected in each eye and matched into left/right pairs", "#E6F5F0")
ya, ype = yg + hg / 2, yp + hp / 2
R2 = yp + hp + 110; H = R2 + 285
SPX = (C1 + C1W + BX) / 2
line("arrow-speed-to-mean", [(C1 + C1W, 557), (SPX, 557), (SPX, ym), (BX - 4, ym)])
line("arrow-speed-to-pereye", [(SPX, 557), (SPX, ype), (BX - 4, ype)])
line("arrow-mean-to-gate", [(BX + 110, y + hm), (BX + 110, yg - 4)], color=OR)
line("arrow-pereye-to-gate", [(BX + 370, yp), (BX + 370, yg + hg + 4)], color=GR, dash="7,5")
text("label-matched-pairs", BX + 384, yg + hg + 33, "pairs", size=SMALL, fill=GR)
TX, TW = 888, 224
for id_, yy, col, fill, s in (("mean", ym, OR, "#FBEFE8", ("mean-speed", "events")), ("adm", ya, BL, "#E8F1F8", ("admitted", "events")),
                              ("pereye", ype, GR, "#E6F5F0", ("per-eye", "events"))):
    rect(f"section-tag-{id_}", TX, yy - 42, TW, 84, stroke=col, fill=fill, rx=26, sw=2)
    for k, t in enumerate(s):
        text(f"label-tag-{id_}-{k}", TX + TW / 2, yy - 6 + 34 * k, t, size=BODY, fill=col, anchor="middle", weight="bold")
        fit(f"tag-{id_}-{k}", t, BODY, TW - 16, True)
    line(f"arrow-to-tag-{id_}", [(BX + BW, yy), (TX - 4, yy)], color=col)
    line(f"bracket-{id_}", [(TX + TW, yy), (1135, yy)], color=INK, sw=2, arrow=False)
line("bracket-spine", [(1135, ym), (1135, ype)], color=INK, sw=2, arrow=False)
# optional centring branch vs unfiltered bypass, applied to each event set
CX, CW, OUTX = 1190, 350, 1572
text("label-stage-3", CX - 30, 38, "Optional centring", size=TITLE, fill=GREY, weight="bold")
yc = ya - 30
line("path-unfiltered", [(1135, ya), (1160, ya), (1160, yc - 60), (OUTX, yc - 60), (OUTX, R2 - 50)], color=INK, arrow=False)
text("label-unfiltered", CX + CW / 2, yc - 74, "unfiltered", size=BODY, anchor="middle", style="italic")
line("path-to-centring", [(1160, ya), (1160, yc + 90), (CX - 4, yc + 90)], color=INK)
hc = box("centring", CX, yc, CW, INK, "Centring", ["keep events with", "|*l* − median *l*| < log 2"], "#FFFFFF", dash="9,6")
line("path-centred", [(CX + CW, yc + 90), (OUTX, yc + 90)], color=INK, arrow=False)
text("label-centred", CX + CW / 2, yc + hc + 38, "centred", size=BODY, anchor="middle", style="italic")
text("label-each-set", CX + CW / 2, yc + hc + 84, "applied to each event set", size=SMALL, fill=GREY, anchor="middle", style="italic")
fit("each-set", "applied to each event set", SMALL, OUTX - 1160 - 20)

# ---- (3) row 2: marker definition and result inset
line("path-to-marker", [(OUTX, R2 - 50), (OUTX, R2 - 50), (470, R2 - 50), (470, R2 - 4)], color=INK)
text("label-stage-4", 20, R2 - 14, "Marker per recording", size=TITLE, fill=GREY, weight="bold")
MW = 900
rect("section-marker", 20, R2, MW, 270, stroke=INK, fill="#FFFFFF")
text("label-marker-title", 44, R2 + 44, "Recordings that meet the event floors", size=TITLE, weight="bold")
fit("marker-title", "Recordings that meet the event floors", TITLE, MW - 48, True)
for k, t in enumerate(["*l* = log(left / right peak velocity)", "*A* = ½ [median *l* (return sweeps)", "        − median *l* (progressive saccades)]"]):
    text(f"label-marker-{k}", 44, R2 + 102 + 46 * k, t, size=BODY)
    fit(f"marker-{k}", t, BODY, 640)
text("label-marker-ex", 44, R2 + 248, "floor 5: ≥ 5 return sweeps and ≥ 20 progressive saccades", size=SMALL - 2, fill=GREY)
fit("marker-ex", "floor 5: ≥ 5 return sweeps and ≥ 20 progressive saccades", SMALL - 2, MW - 48)
# reading-line icon: word bars, progressive hops, return sweep
IX = 690
for i in range(4):
    rect(f"icon-word-{i}", IX + 52 * i, R2 + 132, 40, 12, stroke="none", fill="#C8C8C8", rx=3, sw=0, layer="icons")
for i in range(3):
    x0 = IX + 20 + 52 * i
    L["icons"].append(f'<path id="icon-prog-{i}" d="M{x0},{R2 + 124} q26,-22 52,0" fill="none" stroke="{INK}" stroke-width="2.5" '
                      f'marker-end="url(#arr-222222)"/>')
L["icons"].append(f'<path id="icon-return-sweep" d="M{IX + 180},{R2 + 154} Q{IX + 100},{R2 + 182} {IX + 18},{R2 + 160}" fill="none" stroke="{OR}" stroke-width="3" '
                  f'marker-end="url(#arr-D55E00)"/>')
text("label-icon-prog", IX + 100, R2 + 92, "progressive", size=SMALL, fill=INK, anchor="middle")
text("label-icon-rs", IX + 100, R2 + 208, "return sweep", size=SMALL, fill=OR, anchor="middle")
fit("icon-prog", "progressive", SMALL, 200); fit("icon-rs", "return sweep", SMALL, 210)
IN = 960
line("arrow-marker-to-inset", [(20 + MW, R2 + 135), (IN - 4, R2 + 135)])
rect("section-inset", IN, R2, 1580 - IN, 270, stroke=BL, fill="#E8F1F8", sw=4)
ICX = (IN + 1580) / 2
text("label-inset-title", ICX, R2 + 42, f"R1 S1–S2 ICC(*A*,1), *n* = {N}, floor 5", size=BODY, anchor="middle", weight="bold")
fit("inset-title", f"R1 S1–S2 ICC(A,1), n = {N}, floor 5", BODY, 1580 - IN - 30, True)
PU = "#CC79A7"
for k, (v, sub, col) in enumerate(((ICC_U, "unfiltered", OR), (ICC_ADM, "admitted", BL), (ICC_CE, "≤ 1,000°/s", PU))):
    x = ICX + 195 * (k - 1)
    text(f"label-inset-v{k}", x, R2 + 112, v, size=56, fill=col, anchor="middle", weight="bold")
    text(f"label-inset-s{k}", x, R2 + 150, sub, size=SMALL - 2, fill=col, anchor="middle")
    fit(f"inset-v{k}", v, 56, 185, True); fit(f"inset-s{k}", sub, SMALL - 2, 185)
text("label-inset-wo", ICX, R2 + 202, f"one participant omitted: {WO_U} vs {WO_ADM}", size=BODY, anchor="middle")
fit("inset-wo", f"one participant omitted: {WO_U} vs {WO_ADM}", BODY, 1580 - IN - 30)
text("label-inset-fixed", ICX, R2 + 246, "windows, peaks and classes unchanged", size=SMALL - 2, anchor="middle", style="italic")
fit("inset-fixed", "windows, peaks and classes unchanged", SMALL - 2, 1580 - IN - 30)
bad = [(i, round(a), b) for i, a, b in FIT if a > b]
assert not bad, bad

defs = "".join(f'<marker id="arr-{c.strip("#")}" viewBox="0 0 10 10" refX="8.5" refY="5" markerWidth="5" markerHeight="5" orient="auto">'
               f'<path d="M0,0 L10,5 L0,10 z" fill="{c}"/></marker>' for c in (INK, OR, BL, GR))
svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">', f"<defs>{defs}</defs>",
       f'<g id="background"><rect id="canvas" x="0" y="0" width="{W}" height="{H}" fill="#FFFFFF"/></g>']
svg += [f'<g id="{k}">' + "".join(v) + "</g>" for k, v in L.items()]
svg.append("</svg>")
out = Path(__file__).with_name("fig1_schematic.svg"); out.write_text("\n".join(svg), encoding="utf-8")
import cairosvg
cairosvg.svg2pdf(url=str(out), write_to=str(out.with_suffix(".pdf")))
cairosvg.svg2png(url=str(out), write_to=str(out.with_suffix(".png")), output_width=W)
print("Saved:", out, out.with_suffix(".pdf"), out.with_suffix(".png"), f"inset {ICC_U} / {ICC_ADM} / {ICC_CE}, omitted {WO_U} vs {WO_ADM}, n = {N}")
