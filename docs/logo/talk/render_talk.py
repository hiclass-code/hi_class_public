#!/usr/bin/env python3
"""Render the animated logo to a 1920x1080 clip for talks.

The clip is produced by the site's own `rain.js`, driven frame by frame in a
real browser, rather than by a second implementation that would drift from it.
`render.html` sets up a controllable `requestAnimationFrame`, loads `rain.js`
with `data-standalone` and `data-loop`, composites the wordmark and caret over
each frame and POSTs it here as a PNG.

    python3 docs/logo/talk/render_talk.py serve          # then open the URL it prints
    python3 docs/logo/talk/render_talk.py encode

`data-loop="12.8"` makes the motion exactly periodic, so frame 320 comes back
pixel-identical to frame 0 and the clip loops without a seam; 12.8 s also holds
a whole number of 0.8 s caret blinks. Encoding drops that duplicate closing
frame and uses one GOP per loop.

Encoding needs ffmpeg. There is no system binary on this machine, so run the
encode step with a Python that has imageio-ffmpeg -- `~/venvs/main/bin/python`
does, the bare `python3` does not. Outputs land in `out/` next to this file;
file them into the animation store yourself, as videos do not belong in git.
"""

import base64
import http.server
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DOCS = os.path.abspath(os.path.join(HERE, os.pardir, os.pardir))   # .../docs
FRAMES = os.path.join(HERE, "frames")
OUT = os.path.join(HERE, "out")

PORT = 8814
FPS = 25
LOOP_SECONDS = 12.8
FRAMES_N = int(round(LOOP_SECONDS * FPS)) + 1        # +1 proves the loop closes
BASE = "logo_rain_talk_dark"                          # <family>_<geom>_<theme>


class Handler(http.server.SimpleHTTPRequestHandler):
    def __init__(self, *a, **kw):
        super().__init__(*a, directory=DOCS, **kw)

    def do_GET(self):
        if self.path.split("?")[0] == "/_render.html":
            body = open(os.path.join(HERE, "render.html"), "rb").read()
            self.send_response(200)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)
            return
        return super().do_GET()

    def do_POST(self):
        n = int(self.headers.get("Content-Length", 0))
        raw = self.rfile.read(n).decode("ascii")
        index = int(self.path.rsplit("/", 1)[-1])
        with open(os.path.join(FRAMES, "f%05d.png" % index), "wb") as fh:
            fh.write(base64.b64decode(raw.split(",", 1)[1]))
        self.send_response(204)
        self.end_headers()

    def log_message(self, *a):
        pass


def ffmpeg():
    exe = shutil.which("ffmpeg")
    if exe:
        return exe
    try:
        import imageio_ffmpeg
        return imageio_ffmpeg.get_ffmpeg_exe()
    except ImportError:
        sys.exit("error: no ffmpeg on PATH and imageio-ffmpeg is not installed")


def serve():
    os.makedirs(FRAMES, exist_ok=True)
    url = "http://127.0.0.1:%d/_render.html?frames=%d" % (PORT, FRAMES_N)
    print("frames -> %s" % FRAMES)
    print("open:      %s" % url)
    print("the tab title counts the frames; Ctrl-C here when it says done")
    http.server.ThreadingHTTPServer(("127.0.0.1", PORT), Handler).serve_forever()


def encode():
    got = sorted(f for f in os.listdir(FRAMES) if f.endswith(".png")) if os.path.isdir(FRAMES) else []
    if len(got) < FRAMES_N:
        sys.exit("error: found %d frames, expected %d -- run `serve` first" % (len(got), FRAMES_N))

    keep = FRAMES_N - 1                      # drop the duplicate closing frame
    os.makedirs(OUT, exist_ok=True)
    ff = ffmpeg()
    src = ["-framerate", str(FPS), "-start_number", "0",
           "-i", os.path.join(FRAMES, "f%05d.png"), "-frames:v", str(keep)]
    # Same codec settings as utils/anim_encode.py in application_plots.
    jobs = [
        (BASE + ".mp4", ["-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "20",
                         "-preset", "slow", "-movflags", "+faststart"]),
        (BASE + ".webm", ["-c:v", "libvpx-vp9", "-pix_fmt", "yuv420p", "-crf", "32",
                          "-b:v", "0", "-row-mt", "1"]),
    ]
    for name, args in jobs:
        dest = os.path.join(OUT, name)
        subprocess.run([ff, "-y", "-loglevel", "error"] + src + args
                       + ["-g", str(keep), dest], check=True)
        print("  %-28s %6.1f MB" % (name, os.path.getsize(dest) / 1e6))

    poster = int(0.5 * keep)
    shutil.copyfile(os.path.join(FRAMES, "f%05d.png" % poster),
                    os.path.join(OUT, BASE + "_still.png"))
    print("  %-28s (frame %d)" % (BASE + "_still.png", poster))
    print("out -> %s" % OUT)


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd == "serve":
        serve()
    elif cmd == "encode":
        encode()
    else:
        sys.exit(__doc__.strip())
