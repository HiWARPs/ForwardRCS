"""Fetch SAID DSG and E (beam-target helicity asymmetry) vs angle at fixed Elab from the public GO3 form (no login), for
gamma p -> pi0 p (rt=1) and pi+ n (rt=2); integrate: sigma = 2 pi int DSG dcos, (s_3/2 - s_1/2)/2 = -2 pi int E DSG dcos
(E = (s_1/2 - s_3/2)/(s_1/2 + s_3/2): E = +1 at 0 and 180 deg, checked).  Usage: python fetch_said_pol.py SOL OUT.dat [energies MeV]"""
import sys, re, time, urllib.request, numpy as np
URL = 'https://gwdac.phys.gwu.edu/cgi-bin/go3pr2?sl=%s&rt=%d&ot=%s&iv=A&il=0&ii=5&iu=180&fv=E&fn=%s&jpeg=GO99&u=0&l=0'
num = re.compile(r'^\s*(\d+\.\d+)\s+(-?\d\.\d+E[+-]\d+)\s+(-?\d\.\d+E[+-]\d+)')
sol, out = sys.argv[1], sys.argv[2]
Es = [float(e) for e in sys.argv[3:]] or list(np.arange(145, 152.01, 0.5)) + list(range(153, 171)) + list(range(172, 221, 2)) + list(range(230, 301, 10))
hdr = ''
def get(rt, ot, E):
	global hdr
	for t in range(3):
		try:
			txt = urllib.request.urlopen(URL % (sol, rt, ot, '%.2f' % E), timeout=120).read().decode('latin-1')
			break
		except Exception as e:
			time.sleep(5)
	rows = []
	for l in re.sub('<[^>]*>', '', txt).splitlines():
		m = num.match(l)
		if m: rows.append((float(m.group(1)), float(m.group(2))))
		elif 'CHI/DP' in l and not hdr: hdr = l.strip()
	time.sleep(0.3)
	return np.array(rows)
res = []
for E in Es:
	r = [E]
	for rt in (1, 2):
		d, e = get(rt, 'DSG', E), get(rt, 'E', E)
		if len(d) != 37 or len(e) != 37:
			r += [np.nan, np.nan]; continue
		th = np.radians(d[:, 0]); w = np.sin(th)
		simp = lambda f: np.pi/36/3*np.sum(f*np.r_[1, [4, 2]*17, 4, 1])	# Simpson on 0..180 deg, 5 deg steps
		r += [2*np.pi*simp(d[:, 1]*w), -2*np.pi*simp(e[:, 1]*d[:, 1]*w)]
	res.append(r)
	print(' '.join('%.5g' % x for x in r), flush=True)
np.savetxt(out, np.array(res), fmt='%.6g', header='SAID %s via %s (DSG, E vs angle, 5 deg; Simpson), retrieved %s, no login. Server header: %s\nElab[MeV] sigma(pi0p) dsig/2(pi0p) sigma(pi+n) dsig/2(pi+n)   [mcb]; dsig/2 = (s_3/2 - s_1/2)/2' % (sol, URL.split('?')[0], time.strftime('%Y-%m-%d'), hdr))
