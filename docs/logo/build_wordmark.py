#!/usr/bin/env python3
"""Lift the wordmark out of the original logo GIF as a transparent PNG.

Taken from the GIF's **second** frame, the one where the caret is hidden, so the
caret can be a real element again instead of a baked-in pixel block.

Two things make this more than a crop:

* **Alpha from luminance, un-premultiplied against black.** The source is
  green-on-black with a soft glow. Keying on luminance and dividing the colour
  back out recovers the glow's hue instead of letting it fade to grey, and gives
  a wordmark that composites over moving rain without a black square behind it.

* **Thickness, not size, separates wordmark from rain.** A few rain glyphs fall
  inside the crop and have to go. They are not distinguished by being *small* --
  the dot of the "i" is smaller than some of them -- but by being *thin*: the
  wordmark is set in a heavy face with ~20px strokes, while a rain glyph is drawn
  in 3-4px strokes. So the mask is a morphological reconstruction: erode, then
  keep whole components that survived. Filtering on bounding-box size instead
  drops the dot of the "i", which is what this script was written to fix.

    python3 docs/logo/build_wordmark.py
"""

import os
import sys
from collections import deque

import numpy as np
from PIL import Image, ImageFilter

HERE = os.path.dirname(os.path.abspath(__file__))
IMGS = os.path.abspath(os.path.join(HERE, os.pardir, "imgs"))
SRC = os.path.join(IMGS, "hi_class.gif")
OUT = os.path.join(IMGS, "hi_class_wordmark.png")

# Crop of the 720x729 source: the wordmark, the caret's slot, and room for glow.
BOX = (55, 296, 672, 450)
FRAME = 1              # the frame with the caret hidden
INK = 90               # luminance that counts as ink when finding components
ERODE = 7              # strokes thinner than this are rain, not wordmark
GROW, SOFTEN = 15, 7   # dilate the kept cores to recover their glow, then feather


def components(mask):
    """Yield (pixels, bbox) for each 4-connected component of a boolean mask."""
    h, w = mask.shape
    seen = np.zeros_like(mask)
    for sy in range(h):
        for sx in range(w):
            if not mask[sy, sx] or seen[sy, sx]:
                continue
            q = deque([(sy, sx)])
            seen[sy, sx] = True
            px = []
            y0 = y1 = sy
            x0 = x1 = sx
            while q:
                y, x = q.popleft()
                px.append((y, x))
                y0, y1 = min(y0, y), max(y1, y)
                x0, x1 = min(x0, x), max(x1, x)
                for dy, dx in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    ny, nx = y + dy, x + dx
                    if 0 <= ny < h and 0 <= nx < w and mask[ny, nx] and not seen[ny, nx]:
                        seen[ny, nx] = True
                        q.append((ny, nx))
            yield px, (x0, y0, x1, y1)


def main():
    if not os.path.exists(SRC):
        sys.exit("error: %s not found" % SRC)
    im = Image.open(SRC)
    im.seek(FRAME)
    x0, y0, x1, y1 = BOX
    rgb = np.asarray(im.convert("RGB")).astype(float)[y0:y1, x0:x1]
    lum = rgb.max(axis=2)

    core = lum > INK
    # Anything that survives the erosion is a heavy stroke; grow each survivor
    # back to its whole component so letters come through entire.
    eroded = np.asarray(
        Image.fromarray((core * 255).astype(np.uint8)).filter(ImageFilter.MinFilter(ERODE))
    ) > 127

    keep = np.zeros_like(core)
    kept = 0
    for px, _ in components(core):
        if any(eroded[y, x] for y, x in px):
            kept += 1
            for y, x in px:
                keep[y, x] = True
    if not kept:
        sys.exit("error: nothing survived the erosion; check INK/ERODE")

    grown = Image.fromarray((keep * 255).astype(np.uint8))
    grown = grown.filter(ImageFilter.MaxFilter(GROW)).filter(ImageFilter.GaussianBlur(SOFTEN))
    mask = np.asarray(grown).astype(float) / 255.0

    alpha = np.clip(lum / 255.0, 0, 1) * mask
    with np.errstate(divide="ignore", invalid="ignore"):
        hue = np.where(lum[..., None] > 0, rgb / np.maximum(lum, 1)[..., None] * 255.0, 0)
    out = np.dstack([np.clip(hue, 0, 255), np.clip(alpha, 0, 1) * 255]).astype(np.uint8)
    Image.fromarray(out, "RGBA").save(OUT, optimize=True)

    ys, xs = np.nonzero(alpha > 0.05)
    print("%d components kept" % kept)
    print("ink extent in crop: x %d..%d, y %d..%d" % (xs.min(), xs.max(), ys.min(), ys.max()))
    print("%s (%.1f KB)" % (os.path.relpath(OUT), os.path.getsize(OUT) / 1024))


if __name__ == "__main__":
    main()
