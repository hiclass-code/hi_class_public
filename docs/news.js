/* ---- News freshness flag ------------------------------------------------
   Shared by the home page and the archive (news.html), so the rule cannot
   drift between the two lists -- which it would, since it is four lines of
   date arithmetic that nobody would think to update twice.

   The Jekyll sites (miguelzuma.github.io, the GLoW ERC site) decide this at
   build time and simply do not emit the flag for stale items. This site has no
   build step, so the same rule runs in the browser. Nothing is written into the
   markup, so with JS off no flag appears at all - better than a flag that
   blinks forever on an item from two years ago.
   ------------------------------------------------------------------------ */
(function () {
	'use strict';

	var FRESH_DAYS = 183;                  /* six months */

	var items = document.querySelectorAll('.news li');
	if (!items.length) { return; }
	var cutoff = Date.now() - FRESH_DAYS * 864e5;

	Array.prototype.forEach.call(items, function (li) {
		var t = li.querySelector('time[datetime]');
		if (!t) { return; }
		/* datetime is YYYY-MM; count from the END of that month, so an item
		   dated to the month it happened is not aged out on day one. */
		var parts = t.getAttribute('datetime').split('-');
		var stamp = Date.UTC(+parts[0], +parts[1], 0);
		if (!(stamp > cutoff)) { return; }

		var body = li.querySelector('.news-body');
		if (!body || body.querySelector('.news-flag')) { return; }
		var flag = document.createElement('span');
		flag.className = 'news-flag';
		flag.textContent = 'New';
		body.insertBefore(flag, body.firstChild);
	});
}());
