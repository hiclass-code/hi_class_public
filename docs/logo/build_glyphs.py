#!/usr/bin/env python3
"""Typeset the logo's symbol set with LaTeX and emit it as canvas-ready outlines.

The hi_class logo is a Matrix rain of cosmology and Horndeski notation. Drawing
that notation on a canvas with a system font is not an option: no mono stack has
a script L, a Box, or an alpha with a subscript that sits where LaTeX puts it.
So the glyphs come from LaTeX itself.

    latex -> dvi -> dvisvgm --no-fonts -> outline paths

dvisvgm emits each page as a <defs> of reusable glyph <path>s plus a <use x y>
per placement, which is exactly the decomposition a canvas wants: the 39 symbols
here share only 60 distinct outlines, and LaTeX has already solved the hard part
(where the subscript goes). The result is written as a JS literal rather than
JSON so the page costs one request and no fetch.

Run this only when the symbol set changes; the output is committed.

    python3 docs/logo/build_glyphs.py

Requires a TeX installation and dvisvgm (TeX Live ships both).
"""

import glob
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "glyphs.js")

# The symbol set, read off the original logo by Marcos Vazquez Pingarron.
# `True` marks the alpha-functions, which the logo picks out in a lighter ink --
# they are the quantities hi_class actually solves for, so they are the ones a
# reader should catch. Keep that flag in sync with ACCENT in rain.js.
SYMBOLS = [
    (r"\alpha_{\rm B}", True),
    (r"\alpha_{\rm M}", True),
    (r"\alpha_{\rm K}", True),
    (r"\alpha_{\rm T}", True),
    (r"\alpha_{\rm H}", True),
    (r"\mathcal{L}_{\rm H}", False),
    (r"\sqrt{-g}", False),
    (r"G_2", False),
    (r"G_3", False),
    (r"G_4", False),
    (r"G_5", False),
    (r"R_{\mu\nu}", False),
    (r"R", False),
    (r"\Gamma^{\alpha}_{\mu\nu}", False),
    (r"\Box\phi", False),
    (r"\phi_{,\mu\nu}", False),
    (r"\phi", False),
    (r"M^2_*", False),
    (r"c^2_{\rm s}", False),
    (r"c^2_{\rm T}", False),
    (r"h_+", False),
    (r"h_\times", False),
    (r"V_X", False),
    (r"X", False),
    (r"k^2", False),
    (r"\Phi", False),
    (r"\Psi", False),
    (r"\delta", False),
    (r"\rho", False),
    (r"\theta", False),
    (r"\Omega", False),
    (r"\Pi", False),
    (r"\omega", False),
    (r"\mathcal{P}", False),
    (r"\mathcal{E}", False),
    (r"\mathcal{R}", False),
    (r"\mathcal{H}", False),
    (r"\Lambda", False),
    (r"H(z)", False),
]

TEX = r"""\documentclass[12pt]{article}
\usepackage{amsmath,amssymb}
\pagestyle{empty}
\begin{document}
%s
\end{document}
"""

USE = re.compile(r"<use x='([-\d.]+)' y='([-\d.]+)' xlink:href='#([^']+)'/>")
PATH = re.compile(r"<path id='([^']+)' d='([^']+)'/>")
VIEWBOX = re.compile(r"viewBox='([^']+)'")

# Outlines are in points at a 12pt design size, so the third decimal is a
# thousandth of a point -- far below a device pixel at any size the logo is
# drawn. Trimming there is a third of the payload for nothing visible.
DECIMALS = 3


def shorten(d):
    def cut(m):
        return ("%.*f" % (DECIMALS, float(m.group()))).rstrip("0").rstrip(".") or "0"
    return re.sub(r"-?\d+\.\d+", cut, d)


def main():
    for tool in ("latex", "dvisvgm"):
        if not shutil.which(tool):
            sys.exit("error: %s not found; install TeX Live (it ships both)" % tool)

    tmp = tempfile.mkdtemp(prefix="hiclass-glyphs-")
    try:
        body = "\n".join(r"\hbox{$\displaystyle %s$}\newpage" % s for s, _ in SYMBOLS)
        tex = os.path.join(tmp, "s.tex")
        with open(tex, "w") as fh:
            fh.write(TEX % body)

        subprocess.run(["latex", "-interaction=nonstopmode", "-halt-on-error", "s.tex"],
                       cwd=tmp, check=True, stdout=subprocess.DEVNULL)
        # --exact-bbox measures the ink rather than the paper, so every symbol
        # arrives tight and the rain can centre it in a cell.
        subprocess.run(["dvisvgm", "--no-fonts", "--exact-bbox", "--page=1-",
                        "--output=p%p.svg", "s.dvi"],
                       cwd=tmp, check=True, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)

        # dvisvgm zero-pads %p to the width of the highest page number, so the
        # file for page 1 is p01.svg once there are ten or more symbols.
        pages = sorted(glob.glob(os.path.join(tmp, "p*.svg")))
        if len(pages) != len(SYMBOLS):
            sys.exit("error: dvisvgm wrote %d pages for %d symbols"
                     % (len(pages), len(SYMBOLS)))

        outlines, index, symbols = [], {}, []
        for svg, (src, accent) in zip(pages, SYMBOLS):
            text = open(svg).read()
            for gid, d in PATH.findall(text):
                if gid not in index:
                    index[gid] = len(outlines)
                    outlines.append(shorten(d))
            x0, y0, w, h = (float(v) for v in VIEWBOX.search(text).group(1).split())
            placed = [[index[g], round(float(x) - x0, 2), round(float(y) - y0, 2)]
                      for x, y, g in USE.findall(text)]
            if not placed:
                sys.exit("error: no glyphs placed for %s" % src)
            symbols.append({"g": placed, "w": round(w, 2), "h": round(h, 2),
                            "a": 1 if accent else 0, "t": src})
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    payload = {"outlines": outlines, "symbols": symbols}
    blob = json.dumps(payload, separators=(",", ":"))
    with open(OUT, "w") as fh:
        fh.write("/* Generated by docs/logo/build_glyphs.py -- do not edit by hand.\n"
                 "   %d symbols sharing %d outlines, typeset by LaTeX at 12pt. */\n"
                 % (len(symbols), len(outlines)))
        fh.write("window.HI_GLYPHS=%s;\n" % blob)

    print("%d symbols, %d outlines -> %s (%.1f KB)"
          % (len(symbols), len(outlines), os.path.relpath(OUT), os.path.getsize(OUT) / 1024))


if __name__ == "__main__":
    main()
