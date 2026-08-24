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
| `talk/` | renders a 1920×1080 clip of the same animation, for slides |
| `build_wordmark.py` | lifts the wordmark out of the GIF, emits the PNG. Run by hand, rarely. |
| `../imgs/hi_class_wordmark.png` | **generated** — the wordmark (see below) |
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

Made by `build_wordmark.py` from the GIF's **second** frame, the one where the caret
is hidden, so the caret could become a real element again. Alpha-keyed off luminance
and un-premultiplied against black, which keeps the glow's hue instead of letting it
go grey. It carries real alpha, so unlike the GIF it needs no `mix-blend-mode:
screen` and sits over moving rain without a black square.

A few rain glyphs fall inside the crop and have to be masked out. What separates them
from the wordmark is **thickness, not size** — the wordmark is a heavy face with
~20px strokes, a rain glyph is drawn in 3–4px strokes. So the mask erodes the ink and
keeps whole components that survived. The first version filtered on bounding-box size
instead, which silently dropped the dot of the "i" (16×17px, smaller than some of the
rain glyphs it was trying to reject) and left it lit only by the halo bleeding up from
the stem. If you change the crop or the palette, check the dot.

The caret's slot was measured in the source artwork (x 92.54–98.38 %, y 5.84–94.16 %)
and is positioned in percentages, so it tracks the wordmark at any width. Its blink
is the 0.8 s on / 0.8 s off of the two GIF frames.

## The talk clip

`talk/render_talk.py` produces a 1920×1080 loop for presentations. It does **not**
reimplement the rain: `talk/render.html` installs a controllable
`requestAnimationFrame`, loads this same `rain.js`, composites the wordmark and
caret over each frame and POSTs it to the server, which writes PNGs. Two attributes
on the canvas make that possible, and both are absent on the site:

- `data-standalone` — the canvas is the whole picture rather than a page header, so
  the clearing that exists only to protect the hero's text and bottom border is
  skipped and the rain reaches all four edges.
- `data-loop="12.8"` — makes the motion exactly periodic. Column speeds are
  quantised to a whole number of cycles in the period and the symbol churn is
  switched off, since random churn can never repeat. Frame 320 then comes back
  **pixel-identical** to frame 0, so the clip loops with no seam and no crossfade.
  12.8 s is also a whole number of 0.8 s caret blinks, which is why it is that and
  not 12 or 13.

```
python3 docs/logo/talk/render_talk.py serve             # open the URL it prints, wait
~/venvs/main/bin/python docs/logo/talk/render_talk.py encode   # needs imageio-ffmpeg
```

The page drives itself from the `?frames=` in that URL and counts progress in the
tab title. Encoding drops the duplicate closing frame, uses one GOP per loop, and
reuses the codec settings from `application_plots/utils/anim_encode.py`. Outputs go
to `talk/out/` (gitignored); file them into the animation store at
`~/Dropbox/Projects/animations/hi_class_public/logo_rain/`, since rendered video does
not belong in git.

No credit line is burned into the frame. The clip *is* the wordmark, and the design
is Marcos's — the attribution convention that applies to the physics animations does
not transfer here.

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
