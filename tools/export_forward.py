#!/usr/bin/env python3
"""Re f and Re g of the 2026 update at given photon lab energies, for use as data points in other fits:
values, the full exp+fit covariance (from the covariances of unpol2026 and pol2026; f and g are independent),
and the separate systematics as signed modes (data selection of f, PWA model of g below 0.2 GeV).

	python tools/export_forward.py [E1 E2 ... in MeV]  ->  results/2026/forward_points.json
(default energies: those of the forward points of the 2026 Chiral-EFT paper, 50-150 MeV).
"""

import os, sys, io, json, contextlib, warnings
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)	# the code reads data/... with relative paths
sys.path.insert(0, ROOT)
warnings.filterwarnings('ignore')
import numpy

E_MEV = [float(a) for a in sys.argv[1:]] or [50, 60, 70, 79, 90, 100, 110, 120, 130, 135, 145, 150]

with contextlib.redirect_stdout(io.StringIO()):
	import unpol2026, pol2026
	wt, wp = unpol2026.WholeRegion2026(), pol2026.WholeRegionPol2026('bw')
	sel = [unpol2026.selection_fit(c)[0] for nm, c in unpol2026.SELECTIONS]	# fits restricted to the GHP15 fit-I and fit-II data
	pw = {a: pol2026.WholeRegionPol2026('bw', pwa=a) for a in ('maid', 'sm22', 'cm12')}	# central PWA tables below 0.2 GeV

k = numpy.array(E_MEV)/1000.
f = numpy.array([float(wt.re_f(x)) for x in k])
g = numpy.array([float(wp.re_g(x)) for x in k])
Jf = numpy.array([[wt.df(i, x).real for i in range(len(wt.get_p()))] for x in k])	# d Re f / d p (unpolarized parameters)
Jg = numpy.array([[wp.dg(i, x).real for i in range(len(wp.get_p()))] for x in k])
Cf, Cg = Jf @ wt.get_pcov() @ Jf.T, Jg @ wp.get_pcov() @ Jg.T
assert numpy.allclose(numpy.sqrt(numpy.diag(Cf)), [wt.re_f(x).err for x in k], rtol=1e-6)	# = the propagated errors of re_f, re_g
assert numpy.allclose(numpy.sqrt(numpy.diag(Cg)), [wp.re_g(x).err for x in k], rtol=1e-6)
sel_f = numpy.array([(sel[0].re_f(x) - sel[1].re_f(x))/2 for x in k])	# signed: (fit-I selection - fit-II selection)/2
alt = {a: numpy.array([float(w.integral_g_subtracted(x, err=False)) for x in k]) for a, w in pw.items()}
i, j = max(((a, b) for a in alt for b in alt if a < b), key=lambda ab: numpy.abs(alt[ab[0]] - alt[ab[1]]).max())
pwa_g = (alt[i] - alt[j])/2	# signed half difference of the most different pair of central PWA tables

def modes(C):	# eigenvalues (descending) and the fraction of the total variance in the leading mode
	w = numpy.linalg.eigvalsh(C)[::-1]
	return w.tolist(), float(w[0]/w.sum())

out = dict(
	info='Re f, Re g [mcb GeV] of the 2026 update (unpol2026.py, pol2026.py) at the photon lab energies E [MeV]. cov_f, cov_g: exp+fit '
		'covariance (f and g independent). sel_f: data-selection systematic of f, signed mode (GHP15 fit-I minus fit-II data, halved); '
		'pwa_g: PWA-model systematic of g, signed mode (half the difference of the two most different central PWA tables below 0.2 GeV: '
		'%s, %s). Made by tools/export_forward.py.' % (i, j),
	E=E_MEV, re_f=f.tolist(), re_g=g.tolist(), cov_f=Cf.tolist(), cov_g=Cg.tolist(), sel_f=sel_f.tolist(), pwa_g=pwa_g.tolist(),
	err_f=numpy.sqrt(numpy.diag(Cf)).tolist(), err_g=numpy.sqrt(numpy.diag(Cg)).tolist(),
	modes_f=modes(Cf), modes_g=modes(Cg), pwa_pair=[i, j])
os.makedirs(os.path.join('results', '2026'), exist_ok=True)
json.dump(out, open(os.path.join('results', '2026', 'forward_points.json'), 'w'), indent=1)
print('E [MeV]  Re f  err  sel   Re g  err  pwa')
for r in zip(E_MEV, f, out['err_f'], sel_f, g, out['err_g'], pwa_g):
	print('%5.0f  %9.5f %8.5f %8.5f   %9.5f %8.5f %8.5f' % r)
print('leading-mode fraction of the exp+fit variance: f %.4f, g %.4f' % (out['modes_f'][1], out['modes_g'][1]))
