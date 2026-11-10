"""Fetch the total cross sections (SGT) of the four single-pion photoproduction channels on the proton and the neutron,
SAID SM22 and MAID2007 (solution MD07 of the same interface), from the public GO3 form of gwdac.phys.gwu.edu (no login),
144.6-1500 MeV, into data/2026/single_pion_said.dat (used by neutron2026.py).
rt = 1 gamma p -> pi0 p, 2 gamma p -> pi+ n, 3 gamma n -> pi- p, 4 gamma n -> pi0 n.
The gamma n -> pi- p pages also list a few measured points (printed with E = 0): dropped.
Usage (from the repository root): python tools/fetch_said_single_pion.py"""
import re, time, urllib.request
URL = 'https://gwdac.phys.gwu.edu/cgi-bin/go3pr2?sl=%s&rt=%d&ot=SGT&iv=E&il=%s&ii=%s&iu=%s&fv=A&fn=0&jpeg=GO99&u=0&l=0'
OUT = 'data/2026/single_pion_said.dat'
def _chunks(a, d, b, n=40):	# the form thins the grid beyond ~50 points per request: ask for <= 40 at a time
	m = int(round((b - a)/d)) + 1
	return [('%.1f' % (a + i*d), '%g' % d, '%.1f' % (a + min(i + n - 1, m - 1)*d)) for i in range(0, m, n)]
RANGES = _chunks(144.6, 0.1, 160.0) + _chunks(160.5, 0.5, 200.0) + _chunks(201., 1., 340.) + _chunks(345., 5., 1500.)
NUM = re.compile(r'^\s*(\d+\.\d+)\s+(-?\d\.\d+E[+-]\d+)\s+(-?\d\.\d+E[+-]\d+)\s*$')
COLS = [(sol, rt) for sol in ('SM22', 'MD07') for rt in (1, 2, 3, 4)]

def fetch(sol, rt):
	rows, hdr = {}, ''
	for il, ii, iu in RANGES:
		u = URL % (sol, rt, il, ii, iu)
		for k in range(3):
			try:
				txt = urllib.request.urlopen(u, timeout=120).read().decode('latin-1'); break
			except Exception:
				time.sleep(5)
		for l in re.sub('<[^>]*>', '', txt).splitlines():
			m = NUM.match(l)
			if m and float(m.group(1)) >= 144.0: rows[float(m.group(1))] = float(m.group(2))
			elif 'CHI/DP' in l and not hdr: hdr = l.strip()
		time.sleep(1)
	return rows, hdr

if __name__ == '__main__':
	data, hdrs = {}, {}
	for c in COLS:
		data[c], hdrs[c] = fetch(*c)
		print(c, len(data[c]), hdrs[c], flush=True)
	E = sorted(set.intersection(*(set(d) for d in data.values())))
	with open(OUT, 'w') as f:
		f.write('# Total cross sections [mub] of single-pion photoproduction (observable SGT) from the public SAID GO3 form\n'
			'# https://gwdac.phys.gwu.edu/cgi-bin/go3pr2 (fv=A, fn=0), retrieved %s by tools/fetch_said_single_pion.py.\n'
			'# SM22: W. J. Briscoe et al., Phys. Rev. C 108, 065205 (2023); server header: %s\n'
			'# MD07: MAID2007, D. Drechsel, S. S. Kamalov, L. Tiator, Eur. Phys. J. A 34, 69 (2007), as provided by the same interface.\n'
			'# Elab[GeV]  SM22: pi0p pi+n pi-p pi0n   MD07: pi0p pi+n pi-p pi0n\n' % (time.strftime('%Y-%m-%d'), hdrs[('SM22', 1)]))
		for e in E:
			f.write('%.4f ' % (e/1000) + ' '.join('%.6g' % data[c][e] for c in COLS) + '\n')
