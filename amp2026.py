#!/usr/bin/env python3
"""Forward observables of the 2026 update (unpol2026.py, pol2026.py):
	d sigma/d Omega (theta = 0) = |f|^2 + |g|^2,   -2 Re f g*,   Sigma_2z = -2 Re f g* / (|f|^2 + |g|^2)   (as amp.py).
All errors go through the two covariances (each contains its data, normalization, Regge and PWA-mode parameters; the
two fits share no parameter).  python amp2026.py  writes results/2026/obs.txt.
"""

import os, sys, numpy
from numpy import pi
from compton import hc

R26 = os.path.join('results', '2026')
NB = 10*float(hc)**2	# |f|^2 [mcb^2 GeV^2] / NB = nb/sr

_wt = _wp = None

def regions():
	global _wt, _wp
	if _wt is None:
		import unpol2026, pol2026
		_wt, _wp = unpol2026.WholeRegion2026(), pol2026.WholeRegionPol2026('bw')
	return _wt, _wp

def observables(k):	# -> (ds, err), (-2 Re f g*, err) in nb/sr, (Sigma_2z, err)
	wt, wp = regions()
	f, g = wt.f(k), wp.g(k)
	df = [wt.df(i, k) for i in range(len(wt.get_p()))]
	dg = [wp.dg(i, k) for i in range(len(wp.get_p()))]
	ds = abs(f)**2 + abs(g)**2
	R = (f*g.conjugate()).real
	# derivatives with respect to the parameters of the unpolarized (t) and the polarized (p) fit
	ds_t, ds_p = [2*(f*d.conjugate()).real for d in df], [2*(g*d.conjugate()).real for d in dg]
	R_t, R_p = [(d*g.conjugate()).real for d in df], [(f*d.conjugate()).real for d in dg]
	S_t = [-2*(a*ds - R*b)/ds**2 for a, b in zip(R_t, ds_t)]
	S_p = [-2*(a*ds - R*b)/ds**2 for a, b in zip(R_p, ds_p)]
	err = lambda jt, jp: (wt.dfi_err(lambda i: jt[i])**2 + wp.dfi_err(lambda i: jp[i])**2)**0.5
	return (ds/NB, err(ds_t, ds_p)/NB), (-2*R/NB, 2*err(R_t, R_p)/NB), (-2*R/ds, err(S_t, S_p))

def _row(k):
	return sum(observables(k), ())

def table(nproc=None):
	kk = [round(0.005*i, 3) for i in range(1, 501)] + [3., 4., 5., 7.5, 10., 20., 50., 100.]
	regions()
	import multiprocessing
	with multiprocessing.get_context('fork').Pool(nproc or min(12, os.cpu_count())) as p:
		rows = p.map(_row, kk, chunksize=2)
	os.makedirs(R26, exist_ok=True)
	with open(os.path.join(R26, 'obs.txt'), 'w') as fo:
		fo.write('# forward observables of the proton, 2026 update (amp2026.py): f from unpol2026.py, g from pol2026.py; errors exp+fit (both covariances)\n')
		fo.write('# separate systematics: PWA model of g below 0.2 GeV (results/2026/pol_g_pwa_model.txt); data selection of f (results/2026/obs_selection.txt)\n')
		fo.write(r' \nu    d sigma/d Omega(0)   error     -2 Re f g^*    error      Sigma_2z      error' + '\n')
		fo.write(r'(GeV)      (nb/sr)        (nb/sr)       (nb/sr)       (nb/sr)' + '\n\n')
		fo.write('\n'.join('%.3f   %.8f   %.8f   %.8f   %.8f   %.8f   %.8f' % ((k,) + tuple(r)) for k, r in zip(kk, rows)) + '\n')


def selection_table():	# data-selection systematic of the observables: f from the refits on the fit-I and fit-II data of GHP15
	# (results/2026/unpol_selection.txt, unpol2026.selection_systematic), g from results/2026/pol_g.txt; half the difference
	sel = {'%.3f' % r[0]: r for r in numpy.loadtxt(os.path.join(R26, 'unpol_selection.txt'))}
	g = {'%.3f' % r[0]: r for r in numpy.loadtxt(os.path.join(R26, 'pol_g.txt'), skiprows=5)}
	rows = []
	for key in sorted(set(sel) & set(g), key=float):
		k, (re1, re2, im1, im2), (reg, img) = float(key), sel[key][3:7], (g[key][1], g[key][3])
		o = []
		for f in (complex(re1, im1), complex(re2, im2)):
			ds, R = abs(f)**2 + reg**2 + img**2, f.real*reg + f.imag*img
			o.append((ds/NB, -2*R/NB, -2*R/ds))
		rows.append([k] + [abs(a - b)/2 for a, b in zip(*o)])
	numpy.savetxt(os.path.join(R26, 'obs_selection.txt'), numpy.array(rows), fmt='%.3f' + '   %.8f'*3,
		header='data-selection systematic of the forward observables (amp2026.selection_table): half the difference with f from the 2026 fit\n'
		'restricted to the data of GHP15 fit I and of fit II (results/2026/unpol_selection.txt), g from results/2026/pol_g.txt\n'
		'nu (GeV)   d(d sigma/d Omega(0)) (nb/sr)   d(-2 Re f g*) (nb/sr)   d Sigma_2z')


if __name__ == '__main__':
	table()
	selection_table()
