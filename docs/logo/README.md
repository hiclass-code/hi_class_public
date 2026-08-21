# The animated logo

The hi_class logo, by [Marcos Vázquez Pingarrón](https://marcosvazquez.me/), is a
Matrix rain of cosmology and Horndeski notation with the wordmark burning through
the middle. In `imgs/hi_class.gif` that rain is a **still image** — the only thing
that ever moved was the caret, a second GIF frame blinking at 0.8 s. This puts the
rain itself in motion and leaves the design alone: same symbol set, same greens,
same clear band behind the wordmark, and the alpha-functions still picked out in a
lighter ink.

| file | what it is |
|---|---|
| `build_glyphs.py` | typesets the symbol set with LaTeX, emits `glyphs.js`. Run by hand, rarely. |
| `glyphs.js` | **generated** — 39 symbols sharing 51 outlines. Do not hand-edit. |
| `rain.js` | the animation. Hand-written. |
| `../imgs/hi_class_wordmark.png` | the wordmark, lifted from the GIF (see below) |
| `../imgs/hi_class.gif` | the original. Still the `og:image` and the no-JS fallback. |

## Why LaTeX outlines

The notation is the point — `\mathcal{L}_H`, `\Box\phi`, `\Gamma^\alpha_{\mu\nu}`,
`\sqrt{-g}`, an alpha with a subscript sitting where it belongs. No system font
stack has those, and the site loads no web fonts and no MathJax. So the glyphs come
from LaTeX itself:

    latex → dvi → dvisvgm --no-fonts → outline paths → glyphs.js

dvisvgm hands back exactly the decomposition a canvas wants: a `<defs>` of reusable
glyph outlines plus a `<use x y>` per placement, with LaTeX having already solved
where the subscript goes. The 39 symbols share only 51 outlines, so the payload is
42 KB — about 13 KB over the wire, against the 107 KB GIF it replaces, which
JS-enabled visitors no longer download at all.

To change the symbol set, edit `SYMBOLS` in `build_glyphs.py` and rerun it:

    python3 docs/logo/build_glyphs.py

That needs `latex` and `dvisvgm`; TeX Live ships both. The `True` flag in that list
marks the alpha-functions, which render in `ACCENT` — keep it in sync with `rain.js`.

## The wordmark

Lifted from the GIF's **second** frame, the one where the caret is hidden, so the
caret could become a real element again. Alpha-keyed off luminance and
un-premultiplied against black, which keeps the glow's hue instead of letting it go
grey; then masked to the eight letterform components so the stray rain glyphs caught
in the crop do not come along. It carries real alpha, so unlike the GIF it needs no
`mix-blend-mode: screen` and sits over moving rain without a black square.

The caret's slot was measured in the source artwork (x 92.54–98.38 %, y 5.84–94.16 %)
and is positioned in percentages, so it tracks the wordmark at any width. Its blink
is the 0.8 s on / 0.8 s off of the two GIF frames.

## Things worth knowing

- **The rain pauses itself.** A `<canvas>` is not a CSS animation, so the
  stylesheet's `prefers-reduced-motion` rules cannot reach it — `rain.js` checks the
  query itself, and also stops on a hidden tab and when scrolled out of view. Under
  reduced motion it paints one frame and holds it, which is the original logo.
- **Columns alternate direction**: even columns fall, odd ones rise.
- **The glow is a blit, not a blur.** Each frame is downscaled into a quarter-size
  buffer and drawn back over itself with `lighter`; the browser's own smoothing does
  the blurring. A `shadowBlur` per glyph would cost far more.
- **Holes, not bands.** The rain is cleared with soft ellipses over the wordmark and
  over the text block, rather than a full-width wash — a wash would clear the margins
  too, and the rain running either side of the text is the whole picture. The text
  block is *measured*, so this survives copy changes and reflow.
- Roughly 645 path fills per frame at 30 fps; about 5 ms of work per frame on a
  1705×511 canvas.
- If anything is missing — no `glyphs.js`, no `Path2D`, no 2D context — `rain.js`
  returns silently and the wordmark simply sits on the page unanimated.
