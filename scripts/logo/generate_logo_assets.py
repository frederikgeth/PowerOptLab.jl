#!/usr/bin/env python3
"""Generate the PowerOptLab.jl logo assets (the "feeder flask").

Writes SVG, PNG and PDF versions of every logo variant to
docs/src/assets/logo/, plus the Documenter sidebar logo and favicon in
docs/src/assets/, from the geometry defined below.

The style matches the BMOPFTools.jl logo family (cream ground, dark-brown
outlines, Poppins, Julia's four colours as phases a/b/c and the neutral).
Text is converted to glyph outlines, so the SVGs render identically without
the font installed.

Requirements: Python packages ``fonttools`` and ``pillow``, the
``rsvg-convert`` binary (librsvg), and the Poppins TTF files (SIL Open Font
License; https://github.com/google/fonts):

    python3 scripts/logo/generate_logo_assets.py --fonts /path/to/fonts

The fonts directory must contain Poppins-Regular.ttf, Poppins-Medium.ttf,
Poppins-SemiBold.ttf and Poppins-Bold.ttf.
"""

from __future__ import annotations

import argparse
import io
import shutil
import subprocess
from pathlib import Path

from fontTools.pens.svgPathPen import SVGPathPen
from fontTools.pens.transformPen import TransformPen
from fontTools.ttLib import TTFont
from PIL import Image

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "docs" / "src" / "assets" / "logo"
DOC_ASSETS = REPO / "docs" / "src" / "assets"

NAME = "PowerOptLab"
TAGLINE = "Models and methods for uncertain four-wire distribution networks"

LIGHT = dict(
    ink="#26221F", cream="#FBF7EF", tan="#D3BC8C", wash="#EFE9DD",
    red="#CB3C33", green="#389826", purple="#9558B2", blue="#4063D8",
    muted="#8A7E6C", sub="#6B6257", ground="#F3EDE2",
)
DARK = dict(
    ink="#F2EBDD", cream="#26221F", tan="#7A6A4C", wash="#3A332B",
    red="#D66661", green="#6BAB5B", purple="#AA79C1", blue="#6682DF",
    muted="#B3A68F", sub="#C9BEAB", ground="#1A1714",
)
PAPER, RULE = "#FFFDF8", "#E4D7BE"


# --------------------------------------------------------------------------
# Text to outlines
# --------------------------------------------------------------------------

class Font:
    def __init__(self, path: Path):
        self.tt = TTFont(str(path))
        self.glyphs = self.tt.getGlyphSet()
        self.cmap = self.tt.getBestCmap()
        self.upm = self.tt["head"].unitsPerEm
        self.hmtx = self.tt["hmtx"]

    def width(self, s: str, size: float, tracking: float = 0.0) -> float:
        adv = sum(self.hmtx[self.cmap[ord(ch)]][0] for ch in s)
        return adv * size / self.upm + tracking * max(len(s) - 1, 0)

    def path(self, s: str, size: float, x: float, y: float, tracking: float = 0.0) -> str:
        scale = size / self.upm
        pen = SVGPathPen(self.glyphs)
        cursor = x
        for ch in s:
            name = self.cmap[ord(ch)]
            self.glyphs[name].draw(TransformPen(pen, (scale, 0, 0, -scale, cursor, y)))
            cursor += self.hmtx[name][0] * scale + tracking
        return pen.getCommands()


FONTS: dict[str, Font] = {}


def text(s: str, font: str, size: float, x: float, y: float, fill: str,
         anchor: str = "start", tracking: float = 0.0) -> str:
    f = FONTS[font]
    w = f.width(s, size, tracking)
    if anchor == "middle":
        x -= w / 2
    return f'<path fill="{fill}" d="{f.path(s, size, x, y, tracking)}"/>'


def wordmark(size: float, x: float, y: float, c: dict, anchor: str = "start") -> tuple[str, float]:
    """NAME in Poppins Bold plus '.jl' in Poppins Medium."""
    tracking = -size * 0.026
    w_main = FONTS["bold"].width(NAME, size, tracking)
    total = w_main + tracking + FONTS["medium"].width(".jl", size, tracking)
    if anchor == "middle":
        x -= total / 2
    svg = (text(NAME, "bold", size, x, y, c["ink"], tracking=tracking)
           + text(".jl", "medium", size, x + w_main + tracking, y, c["muted"], tracking=tracking))
    return svg, total


# --------------------------------------------------------------------------
# The mark: a lab flask with four liquid layers (phases a, b, c, neutral),
# rising bubbles and a feeder growing out of the neck. 120x120 units.
# --------------------------------------------------------------------------

FLASK = "M52 37 V58 L22 104 Q20 110 27 110 H93 Q100 110 98 104 L68 58 V37 Z"


def flask(c: dict, small: bool = False, mono: bool = False, uid: str = "f") -> str:
    ink, cream = c["ink"], c["cream"]

    def node(col):
        return (ink if small else cream) if mono else col

    def liquid(col):
        return ink if mono else col

    tan = cream if mono else c["tan"]
    ew, nr, nw, ow, rw = (7, 9, 0, 6.5, 0) if small else (4, 6, 3, 4, 3)
    p = [f'<clipPath id="{uid}-clip"><path d="{FLASK}"/></clipPath>',
         f'<path d="{FLASK}" fill="{cream}"/>',
         f'<g clip-path="url(#{uid}-clip)">'
         f'<rect x="0" y="74" width="120" height="9" fill="{liquid(c["red"])}"/>'
         f'<rect x="0" y="83" width="120" height="9" fill="{liquid(c["green"])}"/>'
         f'<rect x="0" y="92" width="120" height="9" fill="{liquid(c["purple"])}"/>'
         f'<rect x="0" y="101" width="120" height="12" fill="{liquid(c["blue"])}"/>'
         f'<path d="M0 74 H120" stroke="{ink}" stroke-width="2.5"/>'
         # One-colour version: keep the four layers distinguishable.
         + (f'<path d="M0 83 H120 M0 92 H120 M0 101 H120" stroke="{cream}" stroke-width="1.8"/>' if mono else "")
         + '</g>',
         f'<path d="{FLASK}" fill="none" stroke="{ink}" stroke-width="{ow}" stroke-linejoin="round"/>']
    if not small:
        p.append(f'<g fill="{cream}" stroke="{ink}" stroke-width="2">'
                 '<circle cx="50" cy="88" r="3"/><circle cx="68" cy="79" r="2.5"/>'
                 '<circle cx="60" cy="65" r="3.5"/><circle cx="57" cy="49" r="2.8"/></g>')
    p.append(f'<rect x="47" y="30" width="26" height="8" rx="3" fill="{tan}" stroke="{ink}" stroke-width="{rw}"/>'
             f'<path d="M60 30 V22 L44 12 M60 22 L78 14 L95 22" fill="none" stroke="{ink}" '
             f'stroke-width="{ew}" stroke-linecap="round" stroke-linejoin="round"/>')
    if not small:
        p.append(f'<circle cx="60" cy="22" r="4.5" fill="{cream}" stroke="{ink}" stroke-width="2.5"/>')
    p.append(f'<circle cx="44" cy="12" r="{nr}" fill="{node(c["red"])}" stroke="{ink}" stroke-width="{nw}"/>'
             f'<circle cx="78" cy="14" r="{nr}" fill="{node(c["green"])}" stroke="{ink}" stroke-width="{nw}"/>'
             f'<circle cx="95" cy="22" r="{nr}" fill="{node(c["purple"])}" stroke="{ink}" stroke-width="{nw}"/>')
    return "".join(p)


def tile(c: dict, small: bool = False, uid: str = "f") -> str:
    return (f'<rect x="4" y="4" width="112" height="112" rx="28" fill="{c["cream"]}" '
            f'stroke="{c["ink"]}" stroke-width="{6 if small else 4}"/>'
            f'<g transform="translate(10 9) scale(0.82)">{flask(c, small=small, uid=uid)}</g>')


def bmopftools_mark(c: dict) -> str:
    """The BMOPFTools.jl mark (B4), used as the 'built on' badge in the hero."""
    ink, cream = c["ink"], c["cream"]
    return (f'<ellipse cx="68" cy="70" rx="47" ry="36" fill="none" stroke="{c["tan"]}" stroke-width="2.5"/>'
            f'<ellipse cx="68" cy="70" rx="33" ry="25" fill="none" stroke="{c["tan"]}" stroke-width="2.5"/>'
            f'<ellipse cx="68" cy="70" rx="19" ry="14" fill="{c["wash"]}" stroke="{c["tan"]}" stroke-width="2.5"/>'
            f'<g fill="none" stroke="{ink}" stroke-width="4.5" stroke-linecap="round" stroke-linejoin="round">'
            '<path d="M24 24 L52 34 L44 52 L62 58 L60 66 L68 70"/><path d="M52 34 L80 26"/><path d="M44 52 L28 80"/></g>'
            f'<circle cx="15" cy="15" r="9" fill="{cream}" stroke="{ink}" stroke-width="3"/>'
            f'<circle cx="24" cy="24" r="9" fill="{cream}" stroke="{ink}" stroke-width="3"/>'
            f'<circle cx="52" cy="34" r="7" fill="{c["red"]}" stroke="{ink}" stroke-width="3"/>'
            f'<circle cx="80" cy="26" r="7" fill="{c["green"]}" stroke="{ink}" stroke-width="3"/>'
            f'<circle cx="28" cy="80" r="7" fill="{c["purple"]}" stroke="{ink}" stroke-width="3"/>'
            f'<circle cx="68" cy="70" r="10" fill="{c["blue"]}" stroke="{ink}" stroke-width="3"/>')


def svg_doc(w: float, h: float, body: str, title: str) -> str:
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{w:g}" height="{h:g}" '
            f'viewBox="0 0 {w:g} {h:g}" role="img" aria-label="{title}">'
            f'<title>{title}</title>{body}</svg>\n')


def place(inner: str, x: float, y: float, size: float) -> str:
    return f'<g transform="translate({x:g} {y:g}) scale({size / 120:g})">{inner}</g>'


# --------------------------------------------------------------------------
# Compositions
# --------------------------------------------------------------------------

def lockup(c: dict) -> str:
    m = 120
    wm, w_wm = wordmark(48, m + 14, 62, c)
    w = m + 14 + max(w_wm, FONTS["regular"].width(TAGLINE, 14)) + 4
    body = place(flask(c), 0, 0, m) + wm + text(TAGLINE, "regular", 14, m + 14, 89, c["sub"])
    return svg_doc(round(w), m, body, f"{NAME}.jl: {TAGLINE}")


def lockup_stacked(c: dict) -> str:
    m = 150
    _, w_wm = wordmark(46, 0, 0, c)
    w = max(m, w_wm) + 8
    wm, _ = wordmark(46, w / 2, m + 46, c, anchor="middle")
    return svg_doc(round(w), m + 62, place(flask(c), (w - m) / 2, 0, m) + wm, f"{NAME}.jl")


WIRES = [("red", 0), ("green", 14), ("purple", 28), ("blue", 42)]


def wire_band(c: dict, y0: float, x0: float, x1: float, under: float = 8, over: float = 4.5) -> str:
    out = ['<g stroke-linecap="round">']
    for col, dy in WIRES:
        y = y0 + dy
        out.append(f'<path d="M{x0:g} {y:g} H{x1:g}" stroke="{c["ink"]}" stroke-width="{under}"/>'
                   f'<path d="M{x0:g} {y:g} H{x1:g}" stroke="{c[col]}" stroke-width="{over}"/>')
    out.append("</g>")
    return "".join(out)


def card(x: float, c: dict, w: float = 232) -> str:
    return (f'<rect x="{x}" y="106" width="{w}" height="220" rx="18" fill="{c["ink"]}" opacity="0.1"/>'
            f'<rect x="{x}" y="100" width="{w}" height="220" rx="18" fill="{PAPER}" stroke="{c["ink"]}" stroke-width="3"/>')


def ellipse(cx, cy, rx, ry, rot, col):
    return (f'<ellipse cx="{cx}" cy="{cy}" rx="{rx}" ry="{ry}" transform="rotate({rot} {cx} {cy})" '
            f'fill="{col}" fill-opacity="0.18" stroke="{col}" stroke-width="2.5" stroke-dasharray="5 4"/>')


def hero_band(c: dict) -> str:
    """Study areas along the four-wire band (1600x440 units), built on BMOPFTools."""
    ink, red, green, purple, blue, tan = (c[k] for k in ("ink", "red", "green", "purple", "blue", "tan"))
    p = [f'<ellipse cx="800" cy="404" rx="740" ry="20" fill="{ink}" opacity="0.07"/>',
         wire_band(c, 204, 90, 1576)]

    # Built on BMOPFTools
    p.append(f'<circle cx="88" cy="225" r="54" fill="{c["cream"]}" stroke="{ink}" stroke-width="3"/>'
             f'<g transform="translate(44 182) scale(0.72)">{bmopftools_mark(c)}</g>')

    # Estimation: feeder with confidence ellipses
    p.append(card(170, c))
    p.append(f'<g fill="none" stroke="{ink}" stroke-width="3" stroke-linecap="round" stroke-linejoin="round">'
             '<path d="M206 138 L240 170 L276 200 L322 246"/><path d="M240 170 L304 146"/><path d="M276 200 L232 268"/></g>'
             + ellipse(322, 246, 30, 17, 35, blue) + ellipse(322, 246, 16, 9, 35, blue)
             + ellipse(240, 170, 19, 11, -30, red) + ellipse(304, 146, 17, 10, 15, green)
             + ellipse(232, 268, 17, 12, 40, purple)
             + f'<circle cx="198" cy="130" r="9" fill="{c["cream"]}" stroke="{ink}" stroke-width="3"/>'
             f'<circle cx="207" cy="139" r="9" fill="{c["cream"]}" stroke="{ink}" stroke-width="3"/>'
             f'<circle cx="276" cy="200" r="5" fill="{PAPER}" stroke="{ink}" stroke-width="2.5"/>'
             f'<g stroke="{ink}" stroke-width="2.5"><circle cx="240" cy="170" r="6" fill="{red}"/>'
             f'<circle cx="304" cy="146" r="6" fill="{green}"/><circle cx="232" cy="268" r="6" fill="{purple}"/>'
             f'<circle cx="322" cy="246" r="8" fill="{blue}"/></g>')

    # Devices: battery, EV plug, generator
    p.append(card(448, c))
    p.append(f'<rect x="481" y="146" width="20" height="9" rx="2" fill="{ink}"/>'
             f'<rect x="466" y="154" width="50" height="124" rx="9" fill="{PAPER}" stroke="{ink}" stroke-width="3.5"/>'
             f'<rect x="474" y="170" width="34" height="20" rx="3" fill="{red}"/>'
             f'<rect x="474" y="196" width="34" height="20" rx="3" fill="{green}"/>'
             f'<rect x="474" y="222" width="34" height="20" rx="3" fill="{purple}"/>'
             f'<rect x="474" y="248" width="34" height="20" rx="3" fill="{blue}"/>'
             f'<rect x="561" y="168" width="8" height="18" rx="2" fill="{ink}"/>'
             f'<rect x="579" y="168" width="8" height="18" rx="2" fill="{ink}"/>'
             f'<rect x="548" y="184" width="52" height="58" rx="14" fill="{tan}" stroke="{ink}" stroke-width="3.5"/>'
             f'<path d="M574 242 V256 C574 286 620 276 640 300" fill="none" stroke="{ink}" stroke-width="5" stroke-linecap="round"/>'
             f'<circle cx="636" cy="148" r="24" fill="{c["cream"]}" stroke="{ink}" stroke-width="3.5"/>'
             f'<path d="M620 148 Q628 134 636 148 T652 148" fill="none" stroke="{ink}" stroke-width="3.5" stroke-linecap="round"/>')

    # Inverters: inverter symbol feeding four wires
    p.append(card(726, c))
    p.append(f'<path d="M738 198 H760 M738 222 H760" stroke="{ink}" stroke-width="4" stroke-linecap="round"/>'
             f'<rect x="760" y="150" width="112" height="112" rx="12" fill="{c["cream"]}" stroke="{ink}" stroke-width="3.5"/>'
             f'<path d="M764 258 L868 154" stroke="{ink}" stroke-width="3"/>'
             f'<path d="M776 172 H808 M776 184 H808" stroke="{ink}" stroke-width="4" stroke-linecap="round"/>'
             f'<path d="M824 236 Q832 222 840 236 T856 236" fill="none" stroke="{ink}" stroke-width="4" stroke-linecap="round"/>')
    p.append(wire_band(c, 184, 872, 940, under=7, over=4))

    # Envelopes: a fan between dashed limits
    p.append(card(1004, c))
    p.append(f'<path d="M1030 210 L1212 128 L1212 292 Z" fill="{c["wash"]}"/>'
             f'<path d="M1030 210 L1212 166 L1212 254 Z" fill="{tan}" fill-opacity="0.5"/>'
             f'<g fill="none" stroke="{ink}" stroke-width="3" stroke-linecap="round" stroke-dasharray="6 5">'
             '<path d="M1030 210 L1212 128"/><path d="M1030 210 L1212 292"/></g>'
             '<g fill="none" stroke-width="4" stroke-linecap="round">'
             f'<path d="M1030 210 L1204 178" stroke="{red}"/><path d="M1030 210 L1204 200" stroke="{green}"/>'
             f'<path d="M1030 210 L1204 222" stroke="{purple}"/><path d="M1030 210 L1204 246" stroke="{blue}"/></g>'
             f'<circle cx="1030" cy="210" r="7" fill="{PAPER}" stroke="{ink}" stroke-width="3"/>')

    # Operability: voltage nose curve with margin to collapse
    p.append(card(1282, c))
    p.append(f'<path d="M1304 126 V292 H1494" fill="none" stroke="{ink}" stroke-width="3" '
             'stroke-linecap="round" stroke-linejoin="round"/>'
             f'<path d="M1314 142 C1400 150 1462 196 1468 226" fill="none" stroke="{blue}" stroke-width="4" stroke-linecap="round"/>'
             f'<path d="M1468 226 C1470 252 1420 272 1318 280" fill="none" stroke="{blue}" stroke-width="4" '
             'stroke-linecap="round" stroke-dasharray="6 6"/>'
             f'<path d="M1400 176 H1460" stroke="{green}" stroke-width="3.5" stroke-linecap="round"/>'
             f'<path d="M1452 169 L1461 176 L1452 183" fill="none" stroke="{green}" stroke-width="3.5" '
             'stroke-linecap="round" stroke-linejoin="round"/>'
             f'<circle cx="1394" cy="160" r="7" fill="{red}" stroke="{ink}" stroke-width="2.5"/>'
             f'<circle cx="1468" cy="226" r="7" fill="{purple}" stroke="{ink}" stroke-width="2.5"/>')

    for x in (425, 703, 981, 1259):
        p.append(f'<circle cx="{x}" cy="225" r="15" fill="{c["cream"]}" stroke="{ink}" stroke-width="3"/>'
                 f'<path d="M{x - 4} 218 L{x + 3} 225 L{x - 4} 232" fill="none" stroke="{ink}" '
                 'stroke-width="3.5" stroke-linecap="round" stroke-linejoin="round"/>')

    p.append(text("built on", "regular", 15, 88, 306, c["sub"], anchor="middle"))
    p.append(text("BMOPFTools", "semibold", 17, 88, 326, ink, anchor="middle"))
    stages = [(286, "Estimation", "states & parameters"),
              (564, "Devices", "generators, batteries & EVs"),
              (842, "Inverters", "filters, topologies & controls"),
              (1120, "Envelopes", "dynamic operating envelopes"),
              (1398, "Operability", "HELM & voltage margins")]
    for x, title, caption in stages:
        p.append(text(title, "semibold", 26, x, 362, ink, anchor="middle"))
        p.append(text(caption, "regular", 17, x, 390, c["sub"], anchor="middle"))
    return "".join(p)


def hero() -> str:
    c = LIGHT
    m, gap, size = 140, 14, 100
    _, w_wm = wordmark(size, 0, 0, c)
    x0 = (1600 - (m + gap + w_wm)) / 2
    wm, _ = wordmark(size, x0 + m + gap, 180, c)
    body = (f'<rect width="1600" height="720" rx="28" fill="{c["ground"]}"/>'
            + place(flask(c, uid="hero"), x0, 70, m) + wm
            + text(TAGLINE, "medium", 30, 800, 264, c["sub"], anchor="middle")
            + f'<g transform="translate(0 256)">{hero_band(c)}</g>')
    return svg_doc(1600, 720, body, f"{NAME}.jl: {TAGLINE}")


def social() -> str:
    c = LIGHT
    m, dy = 260, 34
    tx = 80 + m + 32
    wm, _ = wordmark(84, tx, 214 + dy, c)
    body = [f'<rect width="1280" height="640" fill="{c["ground"]}"/>',
            place(flask(c, uid="social"), 80, 140 + dy - 8, m), wm,
            text("Models and methods for uncertain", "medium", 30, tx, 268 + dy, c["sub"]),
            text("four-wire distribution networks", "medium", 30, tx, 307 + dy, c["sub"])]
    x = tx
    for label in ["state estimation", "batteries & EVs", "inverters", "operating envelopes", "Julia"]:
        w = FONTS["medium"].width(label, 18) + 32
        body.append(f'<rect x="{x:g}" y="{336 + dy}" width="{w:g}" height="38" rx="19" fill="{c["cream"]}" '
                    f'stroke="{c["ink"]}" stroke-width="2"/>')
        body.append(text(label, "medium", 18, x + 16, 361 + dy, c["ink"]))
        x += w + 10
    body.append(wire_band(c, 570, -10, 1290))
    return svg_doc(1280, 640, "".join(body), f"{NAME}.jl: {TAGLINE}")


# --------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------

def rsvg(src: Path, dst: Path, fmt: str, width: int | None = None) -> None:
    cmd = ["rsvg-convert", "-f", fmt, "-o", str(dst)]
    if width:
        cmd += ["-w", str(width), "--keep-aspect-ratio"]
    subprocess.run(cmd + [str(src)], check=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--fonts", type=Path, required=True, help="directory with the Poppins TTF files")
    args = ap.parse_args()
    if shutil.which("rsvg-convert") is None:
        raise SystemExit("rsvg-convert (librsvg) is required")
    for key, name in [("regular", "Poppins-Regular"), ("medium", "Poppins-Medium"),
                      ("semibold", "Poppins-SemiBold"), ("bold", "Poppins-Bold")]:
        FONTS[key] = Font(args.fonts / f"{name}.ttf")

    OUT.mkdir(parents=True, exist_ok=True)
    title = f"{NAME}.jl logo"
    variants = [  # (file stem, svg text, PNG width; None = natural size plus @2x)
        ("poweroptlab-mark", svg_doc(120, 120, flask(LIGHT), title), 512),
        ("poweroptlab-mark-dark", svg_doc(120, 120, flask(DARK), title), 512),
        ("poweroptlab-mark-mono", svg_doc(120, 120, flask(LIGHT, mono=True), title), 512),
        ("poweroptlab-mark-mono-reversed", svg_doc(120, 120, flask(DARK, mono=True), title), 512),
        ("poweroptlab-mark-small", svg_doc(120, 120, flask(LIGHT, small=True), title), 64),
        ("poweroptlab-mark-small-dark", svg_doc(120, 120, flask(DARK, small=True), title), 64),
        ("poweroptlab-icon", svg_doc(120, 120, tile(LIGHT), title), 512),
        ("poweroptlab-icon-dark", svg_doc(120, 120, tile(DARK), title), 512),
        ("poweroptlab-icon-small", svg_doc(120, 120, tile(LIGHT, small=True), title), 64),
        ("poweroptlab-lockup", lockup(LIGHT), None),
        ("poweroptlab-lockup-dark", lockup(DARK), None),
        ("poweroptlab-lockup-stacked", lockup_stacked(LIGHT), None),
        ("poweroptlab-lockup-stacked-dark", lockup_stacked(DARK), None),
        ("poweroptlab-hero", hero(), 1600),
        ("poweroptlab-social", social(), 1280),
    ]
    for stem, svg, width in variants:
        src = OUT / f"{stem}.svg"
        src.write_text(svg, encoding="utf-8")
        rsvg(src, OUT / f"{stem}.png", "png", width)
        if width is None:
            natural = int(float(svg.split('width="')[1].split('"')[0]))
            rsvg(src, OUT / f"{stem}@2x.png", "png", 2 * natural)
        rsvg(src, OUT / f"{stem}.pdf", "pdf")
    rsvg(OUT / "poweroptlab-hero.svg", OUT / "poweroptlab-hero@2x.png", "png", 3200)

    # Documenter picks up assets/logo.svg and assets/logo-dark.svg automatically.
    shutil.copyfile(OUT / "poweroptlab-mark.svg", DOC_ASSETS / "logo.svg")
    shutil.copyfile(OUT / "poweroptlab-mark-dark.svg", DOC_ASSETS / "logo-dark.svg")

    # Favicon: the simplified tile, rendered natively at each size.
    frames = []
    for px in (16, 32, 48):
        png = subprocess.run(["rsvg-convert", "-f", "png", "-w", str(px), "-h", str(px),
                              str(OUT / "poweroptlab-icon-small.svg")],
                             check=True, capture_output=True).stdout
        frames.append(Image.open(io.BytesIO(png)).convert("RGBA"))
    frames[-1].save(DOC_ASSETS / "favicon.ico", format="ICO",
                    sizes=[(16, 16), (32, 32), (48, 48)], append_images=frames[:-1])
    print(f"wrote {len(variants)} variants to {OUT.relative_to(REPO)}")


if __name__ == "__main__":
    main()
