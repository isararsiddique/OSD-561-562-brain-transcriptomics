"""Build the Paper 1 graphical abstract (figures/Graphical_Abstract.svg) with the Python standard library only.

Every text line is its own <text> element, because some SVG renderers (CoreSVG, used by macOS Preview and AppKit)
misplace multiply-positioned <tspan> children. Before writing, the script estimates the width of each line in
Helvetica metrics and REFUSES to write the file if a line would overflow its box.

House rules (user's figure preferences):
- every box is a rounded rectangle;
- box type is shown by Okabe-Ito fill colour, explained in a key;
- outline style (solid / dashed / dotted) and text labels are the greyscale backup.
Every number is taken from the manuscript text or research/F_signatures/F_REPORT.md.
"""
from pathlib import Path
from xml.sax.saxutils import escape

OUT = Path(__file__).resolve().parent / "Graphical_Abstract.svg"
W, H = 1600, 980
INK, MUTED = "#1a1a1a", "#555555"
KINDS = {  # fill, fill-opacity, stroke, stroke-dasharray, stroke-width
    "protocol": ("#D55E00", 0.12, "#D55E00", None, 2.5),   # Okabe-Ito vermillion, solid outline
    "flight":   ("#0072B2", 0.10, "#0072B2", "9 6", 2.5),  # Okabe-Ito blue, dashed outline
    "method":   ("#009E73", 0.12, "#009E73", None, 2.5),   # Okabe-Ito bluish green, solid outline
    "null":     ("#000000", 0.05, "#000000", "2 5", 2.0),  # black tint, dotted outline
}

# Approximate Helvetica advance widths in em units; unknown characters count as 0.60 em.
_NARROW = set("ijl.,:;'|!()[]/ ")
_WIDE = set("mwMW%@")


def text_width(s, size, bold=False):
    em = 0.0
    for ch in s:
        if ch in _NARROW:
            em += 0.28
        elif ch in _WIDE:
            em += 0.86
        elif ch.isupper() or ch.isdigit():
            em += 0.62
        else:
            em += 0.54
    return em * size * (1.06 if bold else 1.0)


parts, problems = [], []


def rect(x, y, w, h, kind, rx=12):
    fill, op, stroke, dash, sw = KINDS[kind]
    d = f' stroke-dasharray="{dash}"' if dash else ""
    cap = ' stroke-linecap="round"' if kind == "null" else ""
    parts.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{fill}" fill-opacity="{op}" '
                 f'stroke="{stroke}" stroke-width="{sw}"{d}{cap}/>')


def text(x, y, s, size=17, bold=False, italic=False, color=INK, anchor=None, limit=None, where=""):
    if limit is not None and text_width(s, size, bold) > limit:
        problems.append(f"{where}: '{s}' needs ~{text_width(s, size, bold):.0f}px > {limit}px")
    attrs = f'x="{x}" y="{y}" font-size="{size}" fill="{color}"'
    if bold:
        attrs += ' font-weight="700"'
    if italic:
        attrs += ' font-style="italic"'
    if anchor:
        attrs += f' text-anchor="{anchor}"'
    parts.append(f"<text {attrs}>{escape(s)}</text>")


def card(x, y, w, h, kind, lines, title_first=True, where=""):
    """A rounded box with lines of text; the first line is a bold heading if title_first."""
    rect(x, y, w, h, kind)
    inner = w - 32
    for i, ln in enumerate(lines):
        italic = ln.startswith("_") and ln.endswith("_")
        s = ln.strip("_") if italic else ln
        bold = title_first and i == 0
        text(x + 16, y + 30 + 23 * i, s, size=18 if bold else 17, bold=bold, italic=italic, limit=inner, where=where)
    need = 30 + 23 * (len(lines) - 1) + 18
    if need > h:
        problems.append(f"{where}: {len(lines)} lines need {need}px > box height {h}px")


# ------------------------------------------------------------------ title
text(40, 54, "Collection protocol outweighs spaceflight in the fastest brain immediate-early genes", 30, bold=True,
     limit=W - 80, where="title")
text(40, 90, "Seven NASA OSDR brain RNA-seq deposits \u00b7 RRRM-2, RR-10 and RR-23 \u00b7 83 mice \u00b7 "
     "whole animals as units, exact permutation tests", 19, color=MUTED, limit=W - 80, where="subtitle")

# ------------------------------------------------------------------ panels
PX, PW, PY, PH = [30, 422, 814, 1206], 364, 115, 620
heads = ["A design confound", "Protocol outweighs flight", "Only rapid genes move", "Detectable, not removable"]
for i, (x, head) in enumerate(zip(PX, heads), start=1):
    parts.append(f'<rect x="{x}" y="{PY}" width="{PW}" height="{PH}" rx="16" fill="#ffffff" stroke="#bdbdbd" stroke-width="1.5"/>')
    parts.append(f'<rect x="{x + 16}" y="{PY + 14}" width="32" height="32" rx="8" fill="{INK}"/>')
    text(x + 32, PY + 37, str(i), 19, bold=True, color="#ffffff", anchor="middle")
    text(x + 60, PY + 38, head, 21, bold=True, limit=PW - 76, where=f"panel {i} head")

BX, BW = [x + 16 for x in PX], PW - 32

# Panel 1 — the design confound
card(BX[0], 180, BW, 170, "protocol", ["ISS-terminal (ISS-T)", "Ketamine/xylazine overdose,", "exsanguination, thoracotomy",
                                       "Flight: cryochiller on orbit", "Ground: liquid nitrogen"], where="P1 ISS-T")
card(BX[0], 362, BW, 140, "flight", ["Live animal return (LAR)", "Gas anaesthesia on Earth",
                                     "after 24 days of recovery", "Liquid nitrogen"], where="P1 LAR")
card(BX[0], 514, BW, 206, "null", ["The key comparison", "Ground controls that never flew", "were collected both ways",
                                   "(12 per brain region), so the", "protocol effect is measured",
                                   "with no spaceflight at all."], where="P1 key")

# Panel 2 — protocol versus protocol-matched flight (rapid primary-response genes, log2)
x0, scale = PX[1] + 30, 260 / 1.2
parts.append('<defs><pattern id="hatch" width="8" height="8" patternUnits="userSpaceOnUse" patternTransform="rotate(45)">'
             '<rect width="8" height="8" fill="#ffffff"/><line x1="0" y1="0" x2="0" y2="8" stroke="#0072B2" stroke-width="3.2"/>'
             '</pattern></defs>')
text(x0, 196, "Rapid primary-response genes (log\u2082)", 16, color=MUTED, limit=PW - 46, where="P2 axis title")
for label, ytop, prot, flt in [("Cerebellum", 226, 0.70, 0.19), ("Hippocampus", 318, 1.12, 0.28)]:
    text(x0, ytop, label, 16, bold=True)
    wp, wf = round(prot * scale), round(flt * scale)
    parts.append(f'<rect x="{x0}" y="{ytop + 7}" width="{wp}" height="26" rx="4" fill="#D55E00"/>')
    text(x0 + wp + 8, ytop + 26, f"{prot:.2f}", 17, bold=True)
    parts.append(f'<rect x="{x0}" y="{ytop + 37}" width="{wf}" height="26" rx="4" fill="url(#hatch)" stroke="#0072B2" '
                 'stroke-width="2" stroke-dasharray="5 3"/>')
    text(x0 + wf + 8, ytop + 56, f"{flt:.2f}", 17)
    if x0 + wp + 8 + text_width(f"{prot:.2f}", 17, True) > PX[1] + PW - 10:
        problems.append(f"P2 {label}: bar label overflows the panel")
parts.append(f'<line x1="{x0}" y1="392" x2="{x0 + 260}" y2="392" stroke="{MUTED}" stroke-width="1.5"/>')
for v in (0, 0.5, 1.0):
    xt = x0 + round(v * scale)
    parts.append(f'<line x1="{xt}" y1="392" x2="{xt}" y2="399" stroke="{MUTED}" stroke-width="1.5"/>')
    text(xt, 416, "0" if v == 0 else f"{v:.1f}", 14, color=MUTED, anchor="middle")
parts.append(f'<rect x="{x0}" y="430" width="26" height="15" rx="3" fill="#D55E00"/>')
text(x0 + 34, 443, "Protocol effect, never flown", 15, limit=PW - 80, where="P2 legend")
parts.append(f'<rect x="{x0}" y="454" width="26" height="15" rx="3" fill="url(#hatch)" stroke="#0072B2" stroke-width="2" stroke-dasharray="5 3"/>')
text(x0 + 34, 467, "Flight effect, protocol matched", 15, limit=PW - 80, where="P2 legend")
card(BX[1], 482, BW, 112, "protocol", ["Direct exact test", "protocol > matched flight", "q = 0.034 and 0.042"], where="P2 test")
card(BX[1], 606, BW, 114, "flight", ["Three missions", "protocol-matched flight effects",
                                     "\u22120.40 to +0.28 log\u2082"], where="P2 missions")

# Panel 3 — specificity of the signal
card(BX[2], 180, BW, 118, "protocol", ["Rapid primary-response genes", "_Fos, Npas4, Egr1, Arc_",
                                       "up 0.70\u20131.12 log\u2082 with protocol"], where="P3 rapid")
card(BX[2], 310, BW, 80, "null", ["Delayed primary-response genes", "no protocol dominance"], title_first=False, where="P3 delayed")
card(BX[2], 402, BW, 80, "null", ["Secondary-response genes", "no protocol effect"], title_first=False, where="P3 SRG")
card(BX[2], 494, BW, 118, "null", ["Five gene sets from other labs", "(activity, IEG, dissociation;",
                                   "pre-registered): none match"], title_first=False, where="P3 published sets")
card(BX[2], 624, BW, 96, "method", ["Independent re-implementation", "reproduces the key numbers"], where="P3 reproduction")

# Panel 4 — detectability versus removability
card(BX[3], 180, BW, 102, "null", ["Collection sensor (30 genes)", "separates the two cohorts:",
                                   "held-out AUC 0.96\u20131.00"], where="P4 sensor")
card(BX[3], 294, BW, 128, "null", ["Genome-wide loading is only", "as reproducible as noise allows:",
                                   "cosine 0.23 / 0.52 / 0.17 vs", "ceiling 0.21 / 0.50 / 0.15"],
     title_first=False, where="P4 loading")
card(BX[3], 434, BW, 128, "protocol", ["Subtracting the signature", "amplifies other directions 8\u201315\u00d7",
                                       "and worsened prediction", "in 2 of 3 settings"], where="P4 removal")
card(BX[3], 574, BW, 146, "method", ["Gain-capped regression", "held-out error down 1.2\u20137.2%",
                                     "in 7 settings with gain \u2264 2,", "but it does not separate",
                                     "biology from nuisance"], where="P4 constrained")

# ------------------------------------------------------------------ key
key = [("protocol", "Collection protocol (solid outline)"), ("flight", "Spaceflight (dashed outline, hatched bars)"),
       ("method", "Validation and method (solid outline)"), ("null", "Context or null (dotted outline)")]
for (kind, label), kx in zip(key, [40, 420, 860, 1230]):
    rect(kx, 756, 30, 18, kind, rx=5)
    text(kx + 40, 771, label, 16)

# ------------------------------------------------------------------ recommendations
parts.append(f'<rect x="30" y="796" width="1540" height="160" rx="16" fill="#000000" fill-opacity="0.04" stroke="{INK}" stroke-width="1.5"/>')
text(54, 834, "What future spaceflight studies should do", 21, bold=True)
recs = [(54, ["Match the collection protocol", "between flight and ground", "animals by design"]),
        (560, ["Record euthanasia agent, time to", "freezing, freezing hardware and", "time of day for every animal"]),
        (1066, ["Report the gain of any post hoc", "correction alongside its results"])]
for n, (rx_, lines) in enumerate(recs, start=1):
    text(rx_, 872, str(n), 18, bold=True)
    for j, ln in enumerate(lines):
        text(rx_ + 24, 872 + 26 * j, ln, 18, limit=470, where=f"recommendation {n}")

if problems:
    raise SystemExit("Refusing to write; text does not fit:\n  " + "\n  ".join(problems))

svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" '
       'font-family="Helvetica, Arial, sans-serif">\n'
       '<title>Graphical abstract: collection protocol outweighs spaceflight in the fastest brain immediate-early genes</title>\n'
       f'<rect x="0" y="0" width="{W}" height="{H}" fill="#ffffff"/>\n' + "\n".join(parts) + "\n</svg>\n")
OUT.write_text(svg, encoding="utf-8")
print(f"wrote {OUT} ({len(svg):,} bytes, {sum(p.startswith('<text') for p in parts)} text elements)")
