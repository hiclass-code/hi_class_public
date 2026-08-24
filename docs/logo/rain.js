/* Animated hi_class logo: a rain of cosmology and Horndeski notation.
 *
 * The logo by Marcos Vazquez Pingarron is a still of exactly this -- columns of
 * LaTeX-set symbols, the alpha-functions picked out in a lighter ink, the
 * wordmark burning through the middle. Only the caret ever moved. This puts the
 * rain itself in motion and leaves the design alone.
 *
 * Glyphs are real LaTeX outlines (see build_glyphs.py), so a Path2D per outline
 * is built once and every draw is a fill of a shared path. Columns alternate
 * direction: half fall, half rise.
 *
 * Falls back silently to whatever is already in the page if anything is
 * missing -- no glyph payload, no Path2D, reduced motion, a hidden tab.
 */
(function () {
	'use strict';

	var canvas = document.querySelector('.hero-rain');
	if (!canvas || !canvas.getContext || typeof Path2D === 'undefined') { return; }
	var glyphs = window.HI_GLYPHS;
	if (!glyphs || !glyphs.outlines || !glyphs.symbols.length) { return; }

	var ctx = canvas.getContext('2d', { alpha: true });
	if (!ctx) { return; }

	/* ---- palette ---------------------------------------------------------
	 * The greens are sampled from the original logo. ACCENT is the one place
	 * this departs from it: there the alpha-functions are near-white, the same
	 * ink as a trail head, so they only read as "bright". A distinct mint pulls
	 * them out of the green field without leaving it. */
	var BODY   = [0, 205, 51];
	var TAIL   = [6, 90, 13];
	var HEAD   = [222, 255, 216];
	var ACCENT = [157, 255, 208];

	var CELL       = 26;    /* css px between symbols down a column */
	var COL        = 34;    /* css px between columns */
	var TRAIL      = 15;    /* symbols behind each head */
	var SPEED      = [2.6, 6.4];   /* cells per second */
	var CHURN      = 0.006; /* chance per symbol per frame of being swapped */
	var MAX_DPR    = 2;
	var FRAME_MS   = 1000 / 30;

	var paths = glyphs.outlines.map(function (d) { return new Path2D(d); });
	var syms = glyphs.symbols;

	/* Scale LaTeX points to pixels off the median symbol height, so the set
	 * keeps its typographic proportions instead of every symbol being boxed to
	 * the same size -- R_{\mu\nu} really is taller than \delta. */
	var heights = syms.map(function (s) { return s.h; }).sort(function (a, b) { return a - b; });
	var SCALE = (CELL * 0.60) / heights[heights.length >> 1];

	var bloom = document.createElement('canvas');
	var bctx = bloom.getContext('2d');

	var cols = [], W = 0, H = 0, dpr = 1, rows = 0;
	var wordmark = document.querySelector('.hero-wordmark');

	/* Two options, both off on the site and both there so a talk render can be
	 * driven from this same file instead of a second copy of it that drifts:
	 *   data-standalone  the canvas is the whole picture, not a page header --
	 *                    keep the wordmark gap, drop the clearing that exists
	 *                    only to protect the hero's text and bottom border.
	 *   data-loop="T"    make the motion exactly periodic with period T seconds,
	 *                    so frame 0 and frame T*fps are identical and the clip
	 *                    loops without a seam. Column speeds are quantised to a
	 *                    whole number of cycles in T and the symbol churn is
	 *                    switched off, since random churn can never repeat. */
	var STANDALONE = canvas.hasAttribute('data-standalone');
	var LOOP = parseFloat(canvas.getAttribute('data-loop')) || 0;

	function rnd(a, b) { return a + Math.random() * (b - a); }
	function pick() { return (Math.random() * syms.length) | 0; }

	function build() {
		var rect = canvas.getBoundingClientRect();
		if (!rect.width || !rect.height) { return false; }
		dpr = Math.min(window.devicePixelRatio || 1, MAX_DPR);
		W = rect.width; H = rect.height;
		canvas.width = Math.round(W * dpr);
		canvas.height = Math.round(H * dpr);
		ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

		/* A quarter-size buffer upscaled back over the frame is the bloom: the
		 * browser's own smoothing does the blur, for one blit instead of a
		 * shadowBlur on every glyph. */
		bloom.width = Math.max(1, Math.round(W / 4));
		bloom.height = Math.max(1, Math.round(H / 4));

		rows = Math.ceil(H / CELL) + TRAIL + 2;
		var n = Math.ceil(W / COL) + 1;
		/* One full trip of a head, from spawning off one edge to clearing the
		 * other: the period the loop mode has to divide into. */
		var cycle = rows + TRAIL + 1;
		cols = [];
		for (var i = 0; i < n; i++) {
			var cell = [];
			for (var r = 0; r < rows; r++) { cell.push(pick()); }
			var speed = rnd(SPEED[0], SPEED[1]);
			if (LOOP) {
				/* Round to the nearest whole number of cycles in LOOP seconds.
				 * The speeds shift by a few per cent; the spread across columns
				 * survives, which is all the eye reads. */
				speed = Math.max(1, Math.round(speed * LOOP / cycle)) * cycle / LOOP;
			}
			cols.push({
				x: i * COL + COL / 2,
				dir: i % 2 ? -1 : 1,          /* alternate: half rain down, half rise */
				head: rnd(-TRAIL, rows),
				speed: speed,
				cycle: cycle,
				phase: Math.random() * cycle,
				t: 0,
				cell: cell
			});
		}
		return true;
	}

	function ink(rgb, a) {
		return 'rgba(' + rgb[0] + ',' + rgb[1] + ',' + rgb[2] + ',' + a.toFixed(3) + ')';
	}

	function symbol(s, cx, cy) {
		ctx.save();
		ctx.translate(cx - s.w * SCALE / 2, cy - s.h * SCALE / 2);
		ctx.scale(SCALE, SCALE);
		for (var i = 0; i < s.g.length; i++) {
			var g = s.g[i];
			ctx.save();
			ctx.translate(g[1], g[2]);
			ctx.fill(paths[g[0]]);
			ctx.restore();
		}
		ctx.restore();
	}

	/* Text the rain must not fight. The heading, tagline, buttons and badges are
	 * measured rather than guessed at, so this keeps working when the copy
	 * changes length or the hero reflows on a narrow screen. */
	var TEXT = ['h1', '.hero-tagline', '.hero-actions', '.badges'];

	function union(base) {
		var box = null;
		for (var i = 0; i < TEXT.length; i++) {
			var el = document.querySelector('.hero ' + TEXT[i]);
			if (!el) { continue; }
			var r = el.getBoundingClientRect();
			if (!r.width) { continue; }
			var b = { l: r.left - base.left, t: r.top - base.top,
			          r: r.right - base.left, b: r.bottom - base.top };
			box = box ? { l: Math.min(box.l, b.l), t: Math.min(box.t, b.t),
			              r: Math.max(box.r, b.r), b: Math.max(box.b, b.b) } : b;
		}
		return box;
	}

	/* Punch a soft ellipse in the rain. Deliberately elliptical and centred
	 * rather than a full-width band: a band would clear the sides too, and the
	 * rain running down the margins either side of the text is the whole
	 * picture. `soft` is where the hole starts feathering out. */
	function hole(cx, cy, rx, ry, soft, depth) {
		ctx.save();
		ctx.translate(cx, cy);
		ctx.scale(Math.max(rx, 1), Math.max(ry, 1));
		var g = ctx.createRadialGradient(0, 0, 0, 0, 0, 1);
		g.addColorStop(0, 'rgba(0,0,0,' + depth + ')');
		g.addColorStop(soft, 'rgba(0,0,0,' + (depth * 0.94).toFixed(3) + ')');
		g.addColorStop(1, 'rgba(0,0,0,0)');
		ctx.fillStyle = g;
		ctx.beginPath();
		ctx.arc(0, 0, 1, 0, Math.PI * 2);
		ctx.fill();
		ctx.restore();
	}

	function clearFor() {
		var base = canvas.getBoundingClientRect();
		ctx.save();
		ctx.globalCompositeOperation = 'destination-out';

		if (wordmark) {
			var w = wordmark.getBoundingClientRect();
			hole(w.left - base.left + w.width / 2, w.top - base.top + w.height / 2,
			     w.width * 0.66, w.height * 1.15, 0.55, 1);
		}

		/* Everything below exists to protect the page around the canvas, so a
		 * standalone render skips it and keeps rain to all four edges. */
		if (!STANDALONE) {
			var t = union(base);
			if (t) {
				hole((t.l + t.r) / 2, (t.t + t.b) / 2,
				     (t.r - t.l) / 2 + 46, (t.b - t.t) / 2 + 34, 0.72, 0.97);
			}

			/* And a light wash over the bottom edge, so columns fade out instead
			 * of being cut off by the section border below. */
			var lg = ctx.createLinearGradient(0, H * 0.72, 0, H);
			lg.addColorStop(0, 'rgba(0,0,0,0)');
			lg.addColorStop(1, 'rgba(0,0,0,0.85)');
			ctx.fillStyle = lg;
			ctx.fillRect(0, H * 0.72, W, H * 0.28);
		}
		ctx.restore();
	}

	function paint() {
		ctx.clearRect(0, 0, W, H);
		for (var i = 0; i < cols.length; i++) {
			var c = cols[i], head = Math.floor(c.head);
			for (var t = 0; t < TRAIL; t++) {
				var r = head - t * c.dir;
				if (r < 0 || r >= rows) { continue; }
				var y = r * CELL + CELL / 2;
				if (y < -CELL || y > H + CELL) { continue; }
				var s = syms[c.cell[r]];
				var f = 1 - t / TRAIL;
				var a = f * f * 0.92;
				var rgb;
				if (t === 0) { rgb = HEAD; a = 1; }
				else if (s.a) { rgb = ACCENT; a = Math.min(1, a * 1.5 + 0.16); }
				else { rgb = [
					Math.round(TAIL[0] + (BODY[0] - TAIL[0]) * f),
					Math.round(TAIL[1] + (BODY[1] - TAIL[1]) * f),
					Math.round(TAIL[2] + (BODY[2] - TAIL[2]) * f)
				]; }
				ctx.fillStyle = ink(rgb, a);
				symbol(s, c.x, y);
			}
		}
		clearFor();

		bctx.clearRect(0, 0, bloom.width, bloom.height);
		bctx.drawImage(canvas, 0, 0, bloom.width, bloom.height);
		ctx.save();
		ctx.globalCompositeOperation = 'lighter';
		ctx.globalAlpha = 0.55;
		ctx.drawImage(bloom, 0, 0, W, H);
		ctx.restore();
	}

	function advance(dt) {
		for (var i = 0; i < cols.length; i++) {
			var c = cols[i];

			if (LOOP) {
				/* Position is a pure function of elapsed time, so it repeats
				 * exactly every LOOP seconds -- no accumulated wrap error, and
				 * no churn, which could not repeat. */
				c.t += dt;
				var u = (c.phase + c.speed * c.t) % c.cycle;
				c.head = c.dir > 0 ? u - 1 : rows + 1 - u;
				continue;
			}

			c.head += c.dir * c.speed * dt;
			/* Wrap once the whole trail has left the far edge, and re-roll the
			 * column so the loop never shows the same sequence twice. */
			if (c.dir > 0 && c.head - TRAIL > rows) { c.head = -1; c.speed = rnd(SPEED[0], SPEED[1]); }
			if (c.dir < 0 && c.head + TRAIL < 0) { c.head = rows + 1; c.speed = rnd(SPEED[0], SPEED[1]); }
			for (var r = 0; r < rows; r++) {
				if (Math.random() < CHURN) { c.cell[r] = pick(); }
			}
		}
	}

	var last = 0, running = false, raf = 0;

	function tick(now) {
		raf = window.requestAnimationFrame(tick);
		if (now - last < FRAME_MS) { return; }
		var dt = Math.min((now - last) / 1000, 0.1);
		last = now;
		advance(dt);
		paint();
	}

	function start() {
		if (running) { return; }
		running = true;
		last = window.performance ? window.performance.now() : Date.now();
		raf = window.requestAnimationFrame(tick);
	}

	function stop() {
		running = false;
		if (raf) { window.cancelAnimationFrame(raf); raf = 0; }
	}

	var motionQuery = window.matchMedia ? window.matchMedia('(prefers-reduced-motion: reduce)') : null;
	function reduced() { return !!(motionQuery && motionQuery.matches); }

	var onscreen = true;

	function sync() {
		if (!W) { return; }
		/* A canvas is not a CSS animation, so the stylesheet's reduced-motion
		 * rules cannot reach it -- the loop has to opt out itself. A single
		 * painted frame is still the logo, just as the original was. */
		if (reduced() || document.hidden || !onscreen) { stop(); }
		else { start(); }
	}

	function resize() {
		var was = running;
		stop();
		if (!build()) { return; }
		paint();
		if (was || !reduced()) { sync(); }
	}

	if (!build()) { return; }
	paint();
	canvas.classList.add('is-live');

	if (window.ResizeObserver) {
		var ro = new ResizeObserver(function () { resize(); });
		ro.observe(canvas);
	} else {
		window.addEventListener('resize', resize);
	}

	document.addEventListener('visibilitychange', sync);
	if (motionQuery) {
		if (motionQuery.addEventListener) { motionQuery.addEventListener('change', sync); }
		else if (motionQuery.addListener) { motionQuery.addListener(sync); }
	}
	if (window.IntersectionObserver) {
		new IntersectionObserver(function (es) {
			onscreen = es[0].isIntersecting;
			sync();
		}, { rootMargin: '80px' }).observe(canvas);
	}

	sync();
}());
