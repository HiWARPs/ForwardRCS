#!/usr/bin/env python

# 2026 update of the unpolarized forward Compton evaluation of GHP15 [O. Gryniuk, F. Hagelstein, V. Pascalutsa, PRD 92, 074031 (2015)]:
# one global fit to all sigma_tot(gamma p) data above E1 replaces fits I and II.
#  - below E1: average of MAID-2007 (data/tot/maid.dat) and SAID SM22 (data/2026/unpol_said_sm22.dat) times (1 + r(k) t), with r = 2 % (+) half the
#    MAID-SM22 spread (relative) and t ~ N(0, 1) one nuisance mode of the PWA error;
#  - above E1: the GHP parametrization (6 BW + background on [E1, E2], Regge above, Fermi junction at E2 = 2 GeV), all 27 parameters
#    fitted jointly with one normalization nuisance per data set (model-multiplied: y_i vs n_j sigma(k_i), penalty (n_j - 1)^2/delta_j^2), with t,
#    and with continuity at E1 to the same PWA table (one point, sigma_PWA(E1) (1 + r t) +- 0.1 mcb); the background is kept >= 0 (as in fits I
#    and II), and the lowest-chi2 admissible minimum of 48 starts is taken (positive couplings, 10 < Gamma < 600 MeV, no BW peak above 700 mcb,
#    no coincident terms);
#  - covariance from the Jacobian of the whitened residuals (penalty rows included), marginalized over the normalization nuisances, times S^2;
#  - WholeRegion2026 plugs this into the legacy machinery (integrate.WholeRegion) with t as the 28th parameter of the joint vector, so that
#    integral_n, re_f, im_f, f, df propagate all sources (data, normalizations, PWA mode) coherently through self.fit.pcov.
# Usage (from the repository root):
#   python unpol2026.py          results/2026/unpol_* from the stored fit data/2026/unpol_fit.json
#   python unpol2026.py fit      redo the PWA table and the fit first
# Inputs: data/2026/db_unpolarized.csv (curated database), data/tot/maid.dat, data/2026/unpol_said_sm22.dat.  Outputs: data/2026/unpol_pwa_threshold.dat,
# data/2026/unpol_fit.json; results/2026/unpol_sumrules.json, unpol_f.txt, unpol_variants.txt/.json.
# Summary of the method and results: README.md (section 2026 update).

import os, sys, json, csv, time
import numpy
from numpy import inf, pi, log, exp, sqrt
from scipy.optimize import least_squares
from scipy.special import expit
from numpy.polynomial.legendre import leggauss
from compton import *
import integrate, tot
from integrate import val, interp, Region

D = os.path.join('data', '2026')
DB = os.path.join(D, 'db_unpolarized.csv')
PWA = os.path.join(D, 'unpol_pwa_threshold.dat')
SM22 = os.path.join(D, 'unpol_said_sm22.dat')
FIT = os.path.join(D, 'unpol_fit.json')
RES = os.path.join('results', '2026')
M, W0 = float(Mp), flab2cm(E0)
NRES, NP, TAU, CONT = 6, 27, 0.1, 0.1	# CONT: error of the continuity point at E1 [mcb] (2pi starts from zero there); 6 BW (M, Gamma, coupling [mcb GeV^2]) + 5 background ([mcb GeV^(i-2)]) + Regge (C1, p1, C2, p2)
PNAMES = sum([['M%d' % i, 'G%d' % i, 'c%d' % i] for i in range(1, NRES+1)], []) + ['b%d' % i for i in range(5)] + ['C1', 'p1', 'C2', 'p2']
P_GHP1 = [1.213597887738517, 0.11759863296397906, 1.8070611453005943, 1.412837907985687, 0.08283786131491379, 0.06885671965676393,	# GHP fit I
	1.4959639761264092, 0.1365060667178851, 0.7537898615042972, 1.6975382043314127, 0.01876636147047618, 0.0016035207598593224,	# (start values: the
	1.6493591088203243, 0.13533828381678753, 0.38090249529849834, 1.8942779019348484, 0.3020056939402462, 0.7189253846531112,	# fit-I parameters as reproduced from tot.py)
	0.4408471296508886, -11.058979119513122, 74.38288407063168, 22.177863905898846, 37.69396440915787,
	124.09569004841966, 0.1838203135509914, 252.58300460332418, -0.8125910196357076]
P_GHP2 = [1.2147584661838817, 0.0989570946200765, 1.229721754110144, 1.4038985057350712, 0.11819850265151713, 0.18092129997440476,	# GHP fit II (same source),
	1.496853376671454, 0.13339353905054352, 0.7206983524267024, 1.6971722895098658, 0.021141594015444316, 0.0020943423800020133,	# resonances 4 and 5
	1.6479923045506686, 0.13516373807398008, 0.38187127335608706, 1.8937520809378103, 0.32349850288735554, 0.8292421707351484,	# swapped to the order of fit I
	0.26242904489129654, -7.968277707199332, 57.27560389673707, 54.25366836969477, 19.51395525431835,
	124.09569004841966, 0.1838203135509914, 252.58300460332418, -0.8125910196357076]


#--------------------------------------------------------------------------------------------------- PWA table below E1
def pwa_table(kind='avg'):	# columns k, sigma, err, sigma_MAID, sigma_SAID on the MAID grid closed at E0 and E1; kind: avg, MAID, SM22, avg_old (old SAID file)
	m = numpy.loadtxt(os.path.join('data', 'tot', 'maid.dat'))[:, :2]
	s = numpy.loadtxt(os.path.join('data', 'tot', 'said.dat'))[:, :2] if kind == 'avg_old' else numpy.loadtxt(SM22)[:, [0, 3]]
	x = numpy.array([E0] + [k for k in m[:, 0] if E0 < k < E1] + [E1])
	ym, ys = numpy.interp(x, *m.T), numpy.interp(x, *s.T)
	ym[0] = ys[0] = 0.
	y = {'MAID': ym, 'SM22': ys}.get(kind, (ym + ys)/2)
	e = sqrt((0.02*y)**2 + ((ym - ys)/2)**2*(kind in ('avg', 'avg_old')))
	return numpy.array([x, y, e, ym, ys]).T


def make_pwa_table(out=PWA):
	with open(out, 'w') as f:
		f.write('# gamma p -> pi N total cross section below E1 (2pi0 threshold) [mcb]: average of MAID-2007 (data/tot/maid.dat) and SAID SM22\n'
			'# (data/2026/unpol_said_sm22.dat, linearly interpolated to the MAID grid); error = 2 %% (+) |MAID - SM22|/2, fully correlated\n'
			'# k [GeV]  sigma  error  sigma_MAID  sigma_SM22   (made by unpol2026.make_pwa_table)\n')
		f.write('\n'.join('%.17g %.6f %.6f %.6f %.6f' % tuple(a) for a in pwa_table()) + '\n')
	return out


class PWAThreshold(integrate.ThresholdRegion):	# sigma = PWA average + t err below E1; its error is the mode t of the fit (WholeRegion2026)
	
	def __init__(self, k_min=E0, k_max=E1, name='threshold', path=None, filename=None, t=0., **kwargs):
		self.t = t
		integrate.ThresholdRegion.__init__(self, k_min, k_max, name, path=D.split(os.sep)[-1], filename=os.path.basename(PWA), **kwargs)
	
	def update_data(self, filename=None, k_min=None, k_max=None):
		self._cache = {}
		a = numpy.loadtxt(filename or os.path.join('data', self.path, self.filename))
		self.x, self.y, self.err = list(a[:, 0]), list(a[:, 1] + self.t*a[:, 2]), list(a[:, 2])	# err = d sigma/d t
		self.k_max = min(self.x[-1], self.k_max)
		self.x, self.y, self.err = self.cut_xye(k_min or self.k_min, k_max or self.k_max)
		return self
	
	def dsigma_dt(self, k):
		return float(interp(k, self.x, self.err))
	
	def cs_err(self, k):	# no error of its own (the legacy 2 %): it is in the fit covariance through t
		return 0.
	
	def _int_err(self, *args, **kwargs):
		return 0.


def shifted(T, t):	# the PWA table at the mode value t
	T = T.copy()
	T[:, 1] += t*T[:, 2]
	return T


def read_table(fn=PWA):
	return numpy.loadtxt(fn)	# columns k, sigma, err, MAID, SM22


#--------------------------------------------------------------------------------------------------- data
def load_sets(fn=DB, keep=None, norm=None, gap=True, rename=None):	# -> list of data sets {name, x, y, e, norm}; keep(row) -> bool overrides the 'use' column
	sets = {}
	for r in csv.DictReader(open(fn)):
		if not (keep(r) if keep else r['use'] == '1'):
			continue
		if not gap and r['dataset'] in ('Caldwell73', 'Meyer70', 'Michalowski77'):
			continue
		nm = rename(r) if rename else r['dataset']
		s = sets.setdefault(nm, dict(name=nm, x=[], y=[], e=[], norm=float(r['norm_fit'])))
		for k, c in (('x', 'Egamma_GeV'), ('y', 'sigma_ub'), ('e', 'err_fit')):
			s[k].append(float(r[c]))
	for s in sets.values():
		s.update((k, numpy.array(s[k])) for k in 'xye')
		if norm is not None:
			s['norm'] = norm(s['name'], s['norm'])
	return list(sets.values())


#--------------------------------------------------------------------------------------------------- model (vectorized copy of fit.ResonanceRegionFit, AsymptoticRegionFit, Smooth)
def model(k, p, E2=E2, a=0.1, jac=False):	# sigma(k) [mcb] for k >= E1, and d sigma/d p
	k = numpy.atleast_1d(numpy.asarray(k, float))
	W = sqrt(M*(M + 2*k))
	e = expit(-(k - E2)/a)
	d = W - W0
	J = numpy.empty((len(k), NP))
	r = 0.
	for i in range(NRES):
		m, g, c = p[3*i:3*i+3]
		Di = (W - m)**2 + g**2/4
		r = r + c/Di
		J[:, 3*i], J[:, 3*i+1], J[:, 3*i+2] = 2*c*(W - m)/Di**2, -c*g/2/Di**2, 1/Di
	for i in range(5):
		J[:, 18+i] = d**(i - 2)
		r = r + p[18+i]*J[:, 18+i]
	w1, w2 = W**p[24], W**p[26]
	h = 0.5*(p[23]*w1 + p[25]*w2)
	J[:, 23], J[:, 24], J[:, 25], J[:, 26] = 0.5*w1, 0.5*p[23]*w1*log(W), 0.5*w2, 0.5*p[25]*w2*log(W)
	J[:, :23] *= e[:, None]
	J[:, 23:] *= (1 - e)[:, None]
	s = e*r + (1 - e)*h
	return (s, J) if jac else s


#--------------------------------------------------------------------------------------------------- global fit
class GlobalFit:	# residuals: data (y - n_j sigma)/e, the PWA point at E1, penalties (n_j - 1)/delta_j, rows for background >= 0. Variants only: mode != 0
	# makes the PWA error a nuisance t ~ N(0, 1) shared by the table (sigma + t err) and the point at E1 (y + t mode, mode = err(E1)); below = data sets
	# below E1 compared with the table (columns 'tab', 'tabm': sigma and err of the table at their energies)

	def __init__(self, sets, pseudo, E2=E2, a=0.1, mode=0., pos=True, below=()):	# pseudo = (x, y, e)
		self.sets, self.pseudo, self.E2, self.a, self.mode, self.below = sets, [numpy.atleast_1d(numpy.asarray(v, float)) for v in pseudo], E2, a, mode, below
		self.nuis = [s for s in sets if s['norm'] > 0]
		W = numpy.linspace(flab2cm(E1), flab2cm(E2 + 0.5), 80)	# background >= 0 there (as in fits I and II): without it the fit runs into a flat valley
		self.pos = numpy.array([(W - W0)**(i - 2) for i in range(5)]).T if pos else None	# (Delta or broad terms with 1000-2000 mcb peaks cancelled by a negative background)
		self.names = PNAMES + ['n_' + s['name'] for s in self.nuis] + (['t_PWA'] if mode else [])

	def split(self, q):
		n = dict(zip([s['name'] for s in self.nuis], q[NP:NP+len(self.nuis)]))
		return q[:NP], n, (q[-1] if self.mode else 0.)

	def resid(self, q, jac=False):
		p, n, t = self.split(q)
		R, J = [], []
		x, y, e = self.pseudo
		blocks = [(s['x'], s['y'], s['e'], n.get(s['name'], 1.), s['name'] if s['name'] in n else None) for s in self.sets] + [(x, y + t*self.mode, e, 1., None)]
		for x, y, e, f, nm in blocks:
			s, Jm = model(x, p, self.E2, self.a, jac=True)
			R.append((y - f*s)/e)
			Ji = numpy.zeros((len(x), len(q)))
			Ji[:, :NP] = -f*Jm/e[:, None]
			if nm is not None:
				Ji[:, self.names.index('n_' + nm)] = -s/e
			J.append(Ji)
		for s in self.below:
			f = n.get(s['name'], 1.)
			R.append((s['y'] - f*(s['tab'] + t*s['tabm']))/s['e'])
			Ji = numpy.zeros((len(s['x']), len(q)))
			Ji[:, -1] = -f*s['tabm']/s['e']
			if s['name'] in n:
				Ji[:, self.names.index('n_' + s['name'])] = -(s['tab'] + t*s['tabm'])/s['e']
			J.append(Ji)
		for i, s in enumerate(self.nuis):	# penalties
			R.append(numpy.array([(q[NP+i] - 1)/s['norm']]))
			Ji = numpy.zeros((1, len(q))); Ji[0, NP+i] = 1/s['norm']; J.append(Ji)
		if self.mode:
			J[len(self.sets)][:, -1] = self.mode/self.pseudo[2]
			R.append(numpy.array([t]))
			Ji = numpy.zeros((1, len(q))); Ji[0, -1] = 1.; J.append(Ji)
		if self.pos is not None:	# stiff penalty rows (0.1 mcb) where the background is negative
			bg = self.pos @ p[18:23]
			R.append(numpy.minimum(bg, 0.)/TAU)
			Ji = numpy.zeros((len(bg), len(q))); Ji[:, 18:23] = self.pos*(bg < 0)[:, None]/TAU; J.append(Ji)
		return numpy.vstack(J) if jac else numpy.concatenate(R)

	def fit(self, p0=None, quiet=True, max_nfev=100000, cov=True):	# p0: model parameters, or the full vector
		p0 = P_GHP1 if p0 is None else p0
		q0 = numpy.asarray(p0, float) if len(p0) == len(self.names) else numpy.concatenate([p0[:NP], numpy.ones(len(self.nuis)), numpy.zeros(len(self.names) - NP - len(self.nuis))])
		sol = least_squares(self.resid, q0, jac=lambda q: self.resid(q, jac=True), method='lm', x_scale='jac', ftol=1e-14, xtol=1e-14, gtol=1e-14, max_nfev=max_nfev)
		self.sol, self.q = sol, sol.x
		self.chi2 = float(sol.fun @ sol.fun)
		self.chi2_pos = float(sol.fun[len(sol.fun) - len(self.pos):] @ sol.fun[len(sol.fun) - len(self.pos):]) if self.pos is not None else 0.
		self.ndata = sum(len(s['x']) for s in list(self.sets) + list(self.below)) + len(self.pseudo[0])
		self.dof = self.ndata - NP
		self.S = sqrt(self.chi2/self.dof)
		if not cov: return self
		J = self.resid(self.q, jac=True)
		c = sqrt((J**2).sum(0))	# column scaling before the SVD
		U, sv, Vt = numpy.linalg.svd(J/c, full_matrices=False)
		self.cov = (Vt.T/sv**2) @ Vt/numpy.outer(c, c)	# full covariance (JtJ)^-1
		self.cond = (sv[0]/sv[-1])**2	# of the scaled JtJ
		Jp = J[:, :NP]/c[:NP]	# conditional covariance (nuisances fixed), for the error budget
		self.cov_cond = numpy.linalg.inv(Jp.T @ Jp)/numpy.outer(c[:NP], c[:NP])
		if not quiet: print('chi2/dof = %.2f/%d = %.3f, S = %.3f, status %d, nfev %d' % (self.chi2, self.dof, self.chi2/self.dof, self.S, sol.status, sol.nfev))
		return self

	def chi2_sets(self):
		p, n, t = self.split(self.q)
		r = {s['name']: (float((((s['y'] - n.get(s['name'], 1.)*model(s['x'], p, self.E2, self.a))/s['e'])**2).sum()), len(s['x'])) for s in self.sets}
		x, y, e = self.pseudo
		r['PWA@E1'] = (float((((y + t*self.mode - model(x, p, self.E2, self.a))/e)**2).sum()), len(x))
		return r


#--------------------------------------------------------------------------------------------------- fast integrals (independent of the legacy quadrature: variants and checks)
def _gl(a, b, npan, m=8):
	xg, wg = leggauss(m)
	e = numpy.linspace(a, b, npan + 1)
	h = numpy.diff(e)/2
	return ((e[:-1] + h)[:, None] + h[:, None]*xg).ravel(), (h[:, None]*wg).ravel()


class Fast:	# sigma = table below E1 (trapezoid on the nodes, as integrate.ThresholdRegion), model above (Gauss-Legendre)

	def __init__(self, table, E2=E2, a=0.1):
		self.tx, self.ty, self.te = table[:, 0], table[:, 1], table[:, 2]
		self.E2, self.a = E2, a
		parts = [_gl(E1, E2, int((E2 - E1)/0.005) + 1), _gl(E2, 3., int((3. - E2)/0.005) + 1)]
		t, w = _gl(log(3.), log(1e10), 600)
		parts.append((exp(t), w*exp(t)))
		self.gx = numpy.concatenate([u for u, v in parts]); self.gw = numpy.concatenate([v for u, v in parts])

	def set_p(self, p):
		self.p = p
		self.gs, self.gJ = model(self.gx, p, self.E2, self.a, jac=True)
		return self

	def sigma(self, k):
		return float(numpy.interp(k, self.tx, self.ty)) if k <= E1 else float(model(k, self.p, self.E2, self.a)[0])

	def _lin_n(self, n, a, b, s, t):	# hc^(n-1)/(2 pi^2) int_a^b sigma/k^n: linear in s (grid values) and t (table values)
		r = 0.
		if a < E1:
			m = (self.tx >= a) & (self.tx <= min(b, E1))
			x, y = self.tx[m], t[m]/self.tx[m]**n
			r += numpy.sum(numpy.diff(x)*(y[1:] + y[:-1])/2)
		if b > E1:
			m = (self.gx > max(a, E1)) & (self.gx < b)	# a, b are panel edges
			r += numpy.sum(self.gw[m]*s[m]/self.gx[m]**n)
		return float(hc)**(n-1)*r/(2*pi**2)

	def sr(self, n, a, b, err=False):	# value; with err also the gradient in the model parameters and the derivative in t (table column 3 = d sigma/d t)
		v = self._lin_n(n, a, b, self.gs, self.ty)
		if not err: return v
		g = numpy.array([self._lin_n(n, a, b, self.gJ[:, i], 0*self.ty) for i in range(NP)])
		return v, g, self._lin_n(n, a, b, 0*self.gs, self.te)

	def _lin_f(self, k, s, t, sk):	# k^2/(2 pi^2) PV int sigma(p)/(p^2 - k^2), subtracted at k (sk = sigma(k)): linear in s, t, sk
		x = self.tx
		with numpy.errstate(divide='ignore', invalid='ignore'):
			y = (t - sk)/(x**2 - k**2)
		for i in numpy.where(~numpy.isfinite(y))[0]:	# k on a node: mean of the neighbours (as integrate.ThresholdRegion._int)
			y[i] = numpy.mean([y[j] for j in (i-1, i+1) if 0 <= j < len(x) and numpy.isfinite(y[j])])
		r = numpy.sum(numpy.diff(x)*(y[1:] + y[:-1])/2) + numpy.sum(self.gw*(s - sk)/(self.gx**2 - k**2))
		r += sk*(-log(abs((E0 - k)/(E0 + k)))/(2*k))	# sk * PV int_E0^inf dp/(p^2 - k^2)
		return k**2*r/(2*pi**2)

	def re_f(self, k, err=False):
		sk = self.sigma(k) if k >= E0 else 0.
		v = float(e_Mp) + self._lin_f(k, self.gs, self.ty, sk)
		if not err: return v
		Jk = model(k, self.p, self.E2, self.a, jac=True)[1][0] if k > E1 else numpy.zeros(NP)
		g = numpy.array([self._lin_f(k, self.gJ[:, i], 0*self.ty, Jk[i]) for i in range(NP)])
		ek = float(numpy.interp(k, self.tx, self.te)) if E0 <= k <= E1 else 0.
		return v, g, self._lin_f(k, 0*self.gs, self.te, ek)


#--------------------------------------------------------------------------------------------------- fit with several starts
def checks(p):	# problems of a solution: widths, masses, couplings, BW peaks far above any data (cancellations), overlapping pairs
	m, g, c = p[0:18:3], p[1:18:3], p[2:18:3]
	peak = c/(g**2/4)
	r = ['width %d' % (i+1) for i in range(NRES) if not 0.01 < g[i] < 0.6] + ['mass %d' % (i+1) for i in range(NRES) if not flab2cm(E1) < m[i] < 2.3]
	r += ['coupling %d' % (i+1) for i in range(NRES) if c[i] <= 0] + ['peak %d' % (i+1) for i in range(NRES) if peak[i] > 700]	# sigma_tot < 600 mcb everywhere
	r += ['pair %d%d' % (j+1, i+1) for i in range(NRES) for j in range(i) if abs(m[i] - m[j]) < min(g[i], g[j])/2]	# two terms within half a width: not separable
	return r


def _pool(n=12):
	import multiprocessing
	from concurrent.futures import ProcessPoolExecutor
	return ProcessPoolExecutor(n, mp_context=multiprocessing.get_context('fork'))	# fork: the workers inherit the objects built before


def _fit1(a):
	g = GlobalFit(*a[0], **a[1]).fit(a[2], max_nfev=20000, cov=False)
	return g.chi2, g.q, g.sol.status if numpy.isfinite(g.chi2) else -9


def multistart(sets, pseudo, nstart=24, seed=1, **kw):	# lowest chi2 among the converged starts that pass checks(); -> GlobalFit, list of the minima found
	rng = numpy.random.default_rng(seed)
	P0 = [numpy.array(P_GHP1), numpy.array(P_GHP2)] + [numpy.array((P_GHP1, P_GHP2)[i % 2])*(1 + 0.1*rng.standard_normal(NP)) for i in range(nstart)]
	with _pool() as ex:
		R = list(ex.map(_fit1, [((sets, pseudo), kw, p) for p in P0]))
	ok = [r for r in R if r[2] > 0 and not checks(r[1][:NP])]
	g = GlobalFit(sets, pseudo, **kw).fit(min(ok, key=lambda r: r[0])[1])
	T = read_table()
	F = lambda q: Fast(shifted(T, q[-1] if kw.get('mode') else 0.)).set_p(q[:NP])
	minima = sorted([(r[0], F(r[1]).sr(2, E0, inf), [F(r[1]).re_f(k) for k in (0.3, 0.5, 1.)], checks(r[1][:NP])) for r in R if r[2] > 0], key=lambda a: a[0])
	return g, minima


def fit(out=FIT, nstart=46, quiet=False):
	make_pwa_table()
	T = read_table()
	sets = load_sets()
	g, minima = multistart(sets, ([E1], [T[-1, 1]], [CONT]), nstart, mode=T[-1, 2])
	p, t = g.q[:NP], g.q[-1]
	ip = list(range(NP)) + [len(g.q) - 1]	# the model parameters and t; the normalization nuisances are integrated out
	C = g.cov[numpy.ix_(ip, ip)]
	e = sqrt(numpy.diag(g.cov))
	R = C/numpy.outer(e[ip], e[ip])
	J = g.resid(g.q, jac=True)
	J[:sum(len(s['x']) for s in sets)] /= max(g.S, 1.)	# alternative: S on the data rows only (priors and the PWA mode not scaled)
	c = sqrt((J**2).sum(0))
	U, sv, Vt = numpy.linalg.svd(J/c, full_matrices=False)
	C2 = ((Vt.T/sv**2) @ Vt/numpy.outer(c, c))[numpy.ix_(ip, ip)]
	s1, J1 = model(E1, p, jac=True)
	r = dict(info='2026 global fit of sigma_tot(gamma p) above E1 with the PWA mode t (unpol2026.fit). cov: covariance of (p, t), the normalization'
		' nuisances integrated out, not scaled by S; cov_data: of p with the nuisances and t fixed; cov_S_on_data: of (p, t) with only the data errors'
		' scaled by S', date=time.strftime('%Y-%m-%d'), db=DB, table=PWA, E2=g.E2, a=g.a, tau_positivity=TAU, continuity_err=CONT,
		names=PNAMES + ['t_PWA'], p=list(p), t=t, err=list(e[ip]), cov=C.tolist(), cov_data=g.cov_cond.tolist(), cov_S_on_data=C2.tolist(),
		chi2=g.chi2, chi2_positivity=g.chi2_pos, ndata=g.ndata, dof=g.dof, S=g.S, status=int(g.sol.status),
		nuisances={s['name']: dict(value=g.q[g.names.index('n_' + s['name'])], err=e[g.names.index('n_' + s['name'])], delta=s['norm'],
			pull=(g.q[g.names.index('n_' + s['name'])] - 1)/s['norm']) for s in g.nuis},
		chi2_sets=g.chi2_sets(), E1_point=dict(k=E1, sigma_pwa=T[-1, 1], err_pwa=T[-1, 2], t=t, t_err=sqrt(C[-1, -1]), sigma_pwa_at_t=T[-1, 1] + t*T[-1, 2],
			sigma_fit=float(s1[0]), err_fit=float(sqrt(J1[0] @ C[:NP, :NP] @ J1[0]))),
		checks=checks(p), cond_JtJ_scaled=g.cond, cond_corr=float(numpy.linalg.cond(R)), max_rho=float(numpy.max(abs(R - numpy.eye(NP + 1)))),
		minima=[dict(chi2=a, baldin=b, re_f_03_05_10=f, problems=c) for a, b, f, c in minima])
	json.dump(r, open(out, 'w'), indent=1)
	if not quiet:
		print('fit: chi2/dof = %.2f/%d = %.3f, S = %.3f, t = %.3f +- %.3f; %d converged starts, best admissible; problems: %s' % (g.chi2, g.dof, g.chi2/g.dof, g.S, t, sqrt(C[-1, -1]), len(minima), r['checks'] or 'none'))
	return r


#--------------------------------------------------------------------------------------------------- the legacy machinery with the 2026 fit
class WholeRegion2026(tot.WholeRegion):	# p = 27 model parameters + the PWA mode t; below E1 sigma depends only on t, above only on the model
	
	def __init__(self, fitfile=FIT, name='2026', quiet=True, **kwargs):
		self.result = r = json.load(open(fitfile))
		self.t = r['t']
		tot.WholeRegion.__init__(self, b=[E0, E1, r['E2'], inf], name=name, path='2026', filename=os.path.basename(DB), init=False, quiet=quiet,
			ThresholdRegion=lambda *a, **kw: PWAThreshold(*a, t=r['t'], **kw), **kwargs)
		self.smooth.set_p(r['a'])
		rf, hf = self.r_region.fit, self.h_region.fit
		self.fit.get_p = lambda: rf.get_p() + hf.get_p() + (self.t,)
		def set_p(*p):
			rf.set_p(*p[:23])
			hf.set_p(*p[23:NP])
			self.t = p[NP]	# the table keeps the fitted t: only the derivatives in t are needed
		self.fit.set_p = set_p
		self.fit.set_p(*(list(r['p']) + [r['t']]))
		self.S = max(r['S'], 1.)	# S on the full matrix: chi2/dof includes the PWA mode (its prior row and the E1 point) consistently
		self.fit.pcov = numpy.array(r['cov'])*self.S**2
		assert abs(self.fit(0.5) - model(0.5, numpy.array(r['p']), r['E2'], r['a'])[0]) < 1e-9	# legacy model = vectorized copy
	
	def dcs(self, i, k):
		if k <= self.b[1]:
			return self.data_region.dsigma_dt(k) if i == NP else 0.
		return integrate.WholeRegion.dcs(self, i, k)
	
	def cs_err(self, k):
		if k <= self.b[1]:
			return abs(self.data_region.dsigma_dt(k))*sqrt(self.fit.pcov[NP, NP])
		return integrate.WholeRegion.cs_err(self, k)
	
	def dintegral_f(self, i, k, k0=None, k1=None, subtracted=True):	# + the data region, absent from the legacy derivative
		r = integrate.WholeRegion.dintegral_f(self, i, k, k0, k1, subtracted)
		if i == NP:
			r += k**2*self.data_region._int(self.data_region.dsigma_dt, 0, k, self.b[0], self.b[1], singular=Region.singular).val/(2*pi**2)
		return r
	
	def re_f(self, k, err=True):	# the legacy re_f adds the threshold and fit errors in quadrature; here one gradient over (p, t)
		r = self.integral_f(k, err=False)
		return r + val(0, self.dfi_err(lambda i: self.dintegral_f(i, k))) if err else r
	
	def dintegral_n(self, i, n, k0, k1):
		r = 0.
		a, b = max(k0, self.b[0]), min(k1, self.b[1])
		if i == NP and a < b:
			r += self.data_region._int(self.data_region.dsigma_dt, 2-n, 0, a, b, subtracted=False).val
		a, b = max(k0, self.b[1]), min(k1, self.b[-1])
		if i < NP and a < b:
			r += self.fitted_region._int(lambda q: self.dcs(i, q), 2-n, 0, a, b, subtracted=False).val
		return float(hc)**(n-1)*r/(2*pi**2)
	
	@integrate.cached_args_method()
	def integral_n(self, n, k0, k1, err=True):
		r = self.data_region.integral_n(n, k0, k1, err=False) + self.fitted_region.integral_n(n, k0, k1, err=False)
		return r + val(0, self.dfi_err(lambda i: self.dintegral_n(i, n, k0, k1))) if err else r


#--------------------------------------------------------------------------------------------------- results
REGIONS = (('[E0,E1]', E0, E1), ('[E1,E2]', E1, E2), ('[E2,inf)', E2, inf), ('total', E0, inf))
UNITS = {2: '1e-4 fm^3', 4: '1e-4 fm^5', 6: '1e-4 fm^7'}
GHP = {	# current legacy values (tests/test_regression.py) and C3, an independent global fit made in the 2026 audit
	'GHP fit I': {2: (14.29336, 0.12725), 4: (6.082363, 0.090501), 6: (4.358922, 0.078211)},
	'GHP fit II': {2: (13.84881, 0.1257), 4: (6.008037, 0.091294), 6: (4.421464, 0.080107)},
	'referee C3': {2: (14.0222, 0.1681), 4: (5.9936, 0.0745)}}
REF_C3 = {0.05: (-2.8456, 0.0022), 0.1: (-2.2201, 0.0097), 0.15: (-0.4237, 0.0317), 0.2: (1.0051, 0.0487), 0.3: (1.9498, 0.0653), 0.5: (-4.1338, 0.0452), 1.0: (-9.5846, 0.1813)}
NU_CMP = (0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.75)
NU_TAB = [0.005*i for i in range(1, 501)] + [3., 5., 10., 20., 50., 100., 1000., 10000.]
_W = {}


def _ref(a):	# Re f, Im f with errors (worker)
	name, k = a
	w = _W[name]
	return tuple(w.re_f(k)) + tuple(w.im_f(k))


def _fast_f(F, C, k):	# Re f and its unscaled fit error from the independent integrator
	v, g, gt = F.re_f(k, err=True)
	g = numpy.append(g, gt)
	return v, float(sqrt(g @ C @ g))


def results(out=RES, quiet=False):
	w = _W['2026'] = WholeRegion2026()
	r, S = w.result, w.S
	F = Fast(shifted(read_table(), r['t']), r['E2'], r['a']).set_p(numpy.array(r['p']))
	C, Cd, C2 = numpy.array(r['cov']), numpy.array(r['cov_data']), numpy.array(r['cov_S_on_data'])
	Ct = C[:NP, :NP] - numpy.outer(C[:NP, NP], C[NP, :NP])/C[NP, NP]	# p with t fixed (nuisances integrated out)
	sr = {}
	for n in (2, 4, 6):
		v = {nm: w.integral_n(n, a, b) for nm, a, b in REGIONS}
		x, gp, gt = F.sr(n, E0, inf, err=True)
		g = numpy.append(gp, gt)
		V, Vt, Vd = g @ C @ g, gp @ Ct @ gp, gp @ Cd @ gp
		assert abs(v['total'].val - x) < 1e-3 and abs(v['total'].err - S*sqrt(V)) < 0.02*v['total'].err	# legacy quadrature = independent grid
		sr[n] = dict(units=UNITS[n], regions={nm: [a.val, a.err] for nm, a in v.items()},
			budget=dict(data=sqrt(Vd), normalizations=sqrt(Vt - Vd), pwa_mode=sqrt(V - Vt), fit_unscaled=sqrt(V), scale_factor=sqrt(V*(S**2 - 1)),
				total=v['total'].err, total_S_on_data_only=sqrt(g @ C2 @ g)), compare={k: a[n] for k, a in GHP.items() if n in a})
	_W['GHP fit I'], _W['GHP fit II'] = tot.WholeRegion_1(quiet=True), tot.WholeRegion_2(quiet=True)	# for the comparison; built before forking
	with _pool() as ex:	# Re f on the grid
		jobs = [('2026', k) for k in NU_TAB] + [(nm, k) for nm in ('GHP fit I', 'GHP fit II') for k in NU_CMP]
		f = dict(zip(jobs, ex.map(_ref, jobs, chunksize=4)))
	os.makedirs(out, exist_ok=True)
	with open(os.path.join(out, 'unpol_f.txt'), 'w') as fo:
		fo.write(r' \nu       Re f       error      Im f       error' + '\n')
		fo.write(r'(GeV)    (mcb*GeV)  (mcb*GeV)  (mcb*GeV)  (mcb*GeV)' + '\n\n')
		d = ['%.3f   %.8f   %.8f   %.8f   %.8f' % (0, e_Mp, 0, 0, 0)] + ['%.3f   %.8f   %.8f   %.8f   %.8f' % ((k,) + f['2026', k]) for k in NU_TAB]
		fo.write('\n'.join(d) + '\n')
	cmp = {'%g' % k: {'2026': f['2026', min(NU_TAB, key=lambda a: abs(a - k))][:2], 'GHP fit I': f['GHP fit I', k][:2], 'GHP fit II': f['GHP fit II', k][:2],
		'referee C3': REF_C3.get(k), 'fast integrator (2026): value, unscaled error': _fast_f(F, C, k)} for k in NU_CMP}
	res = dict(info='2026 unpolarized forward Compton (unpol2026.py): sum rules in the units given, Re f in mcb GeV. All errors from one covariance of'
		' (27 model parameters, PWA mode t), normalization nuisances integrated out, times S^2. Budget: data = nuisances and t fixed; normalizations ='
		' added by integrating them out; pwa_mode = added by integrating t out (it includes its correlation with the fit through the E1 point);'
		' fit_unscaled = all three; scale_factor = the extra from S; total = fit_unscaled*S; total_S_on_data_only = alternative with S on the data'
		' rows only.', date=time.strftime('%Y-%m-%d'),
		fit=dict((k, r[k]) for k in ('chi2', 'ndata', 'dof', 'S', 't', 'nuisances', 'chi2_sets', 'E1_point', 'checks', 'cond_corr', 'max_rho')),
		sum_rules={'%d' % n: a for n, a in sr.items()}, re_f_compare=cmp)
	json.dump(res, open(os.path.join(out, 'unpol_sumrules.json'), 'w'), indent=1)
	if not quiet:
		print('t = %.3f +- %.3f (unscaled)' % (r['t'], sqrt(C[NP, NP])))
		for n in (2, 4, 6):
			b = sr[n]['budget']
			print('SR%d [%s]: ' % (n, UNITS[n]) + '  '.join('%s %.4f(%.4f)' % (nm, *sr[n]['regions'][nm]) for nm, a, c in REGIONS) +
				'   budget: data %.4f, norm %.4f, pwa %.4f, S %.4f; total %.4f (S on data only: %.4f)' % (b['data'], b['normalizations'], b['pwa_mode'], b['scale_factor'], b['total'], b['total_S_on_data_only']))
		print('Re f [mcb GeV]   2026   GHP I   GHP II   C3')
		for k in NU_CMP:
			c = cmp['%g' % k]
			print('  %.2f  %.3f(%.3f)  %.3f(%.3f)  %.3f(%.3f)  %s' % ((k,) + tuple(c['2026']) + tuple(c['GHP fit I']) + tuple(c['GHP fit II']) + (c['referee C3'] and '%.3f' % c['referee C3'][0],)))
	return res


#--------------------------------------------------------------------------------------------------- variants (fast integrator; one fit each, started from the baseline)
C3SETS = ('LEGS96', 'MacCormick96', 'Armstrong72', 'Bartalini08', 'Caldwell78', 'H1_95', 'ZEUS02')
C3 = dict(table='avg_old', mode=0.03, pseudo_err=0.01, pos=False, gap=False,	# the choices of the independent audit fit C3 (there
	keep=lambda r: r['use'] == '1' and r['dataset'] in C3SETS and not (r['dataset'] == 'Armstrong72' and 0.69 <= float(r['Egamma_GeV']) < 1.45),	# the 3 % PWA
	norm=lambda n, d: {'LEGS96': 0.03, 'MacCormick96': 0.028, 'Armstrong72': 0.03, 'Bartalini08': 0.03, 'Caldwell78': 0.02}.get(n, d))	# scale divides the table)
ARM2 = lambda r: 'Armstrong72 >2 GeV' if r['dataset'] == 'Armstrong72' and float(r['Egamma_GeV']) > 2 else r['dataset']
VARIANTS = [('baseline (correlated PWA mode)', {}),
	('independent PWA error at the E1 point (first version)', dict(mode=None, pseudo_err=None)),
	('+ LEGS and MacCormick data below E1', dict(below=True)),
	('normalizations all 2 %', dict(norm=lambda n, d: 0.02)),
	('normalizations all 5 %', dict(norm=lambda n, d: 0.05)),
	('Armstrong: separate normalization above 2 GeV', dict(rename=ARM2)),
	('MAID table only (2 %)', dict(table='MAID')),
	('SM22 table only (2 %)', dict(table='SM22')),
	('MAID + old SAID file (said.dat)', dict(table='avg_old')),
	('without Armstrong below 0.69 GeV', dict(keep=lambda r: r['use'] == '1' and not (r['dataset'] == 'Armstrong72' and float(r['Egamma_GeV']) < 0.69))),
	('junction E2 = 1.6 GeV', dict(E2=1.6)),
	('junction E2 = 2.5 GeV', dict(E2=2.5)),
	('without the 4-18 GeV data', dict(gap=False)),
	('without the continuity point at E1', dict(pseudo_err=0)),
	('without background positivity', dict(pos=False)),
	('C3a: referee C3 choices', C3),
	('C3b: + background positivity', dict(C3, pos=True)),
	('C3c: + SM22 instead of said.dat', dict(C3, pos=True, table='avg')),
	('C3d: + PWA mode 2 % (+) half spread instead of a 3 % scale', dict(C3, pos=True, table='avg', mode='err', pseudo_err=CONT)),
	('C3e: + 2026 normalizations', dict(C3, pos=True, table='avg', mode='err', pseudo_err=CONT, norm=None)),
	('C3f: + Armstrong 0.69-1.44 GeV and 4-18 GeV data (= baseline)', dict(C3, pos=True, table='avg', mode='err', pseudo_err=CONT, norm=None, keep=None, gap=True))]


def variant(name, table='avg', keep=None, norm=None, gap=True, rename=None, E2=E2, pos=True, mode='err', pseudo_err=CONT, below=False, p0=None):	# mode: 'err'
	T = pwa_table(table)	# (the error column), a fraction (a scale mode) or None (independent PWA error); pseudo_err: None = PWA error, 0 = no point
	m = T[:, 2] if mode == 'err' else (mode*T[:, 1] if mode else 0*T[:, 1])	# d sigma/d t
	pseudo = ([E1], [T[-1, 1]], [T[-1, 2] if pseudo_err is None else pseudo_err]) if pseudo_err != 0 else ([], [], [])
	B = load_sets(keep=lambda r: r['dataset'] in ('LEGS96', 'MacCormick96') and float(r['Egamma_GeV']) < E1, norm=norm) if below else []
	for s in B:
		s['tab'], s['tabm'] = numpy.interp(s['x'], T[:, 0], T[:, 1]), numpy.interp(s['x'], T[:, 0], m)
	G = lambda: GlobalFit(load_sets(keep=keep, norm=norm, gap=gap, rename=rename), pseudo, E2=E2, pos=pos, mode=m[-1], below=B)
	R = [G().fit(q, max_nfev=20000, cov=False) for q in (p0, P_GHP1, P_GHP2)]	# the baseline solution and the GHP fits as starts: best admissible
	g = G().fit(min(R, key=lambda g: (bool(checks(g.q[:NP])) if pos else False, g.chi2)).q)
	p, n, t = g.split(g.q)
	if mode: T[:, 2] = m
	F = Fast(shifted(T, t), E2).set_p(p)
	o = dict(name=name, chi2=g.chi2, dof=g.dof, S=g.S, problems=checks(p), step_E1=float(model(E1, p, E2)[0] - T[-1, 1] - t*T[-1, 2]), t_pwa=t,
		t_err=float(sqrt(g.cov[-1, -1])) if mode else 0., nuisances={k: round(v, 4) for k, v in n.items()})
	for k in (2, 4, 6):
		v, gr, gt = F.sr(k, E0, inf, err=True)
		gr = numpy.concatenate([gr, numpy.zeros(len(g.q) - NP - 1), [gt]]) if mode else gr
		o['SR%d' % k] = [v, float(sqrt(gr @ g.cov[:len(gr), :len(gr)] @ gr))*max(g.S, 1.), 0. if mode else abs(gt)]	# value, fit (incl. PWA mode) error times S, independent PWA error
	o['B_regions'] = [F.sr(2, a, b) for nm, a, b in REGIONS[:3]]
	o['Re_f'] = {'%g' % k: F.re_f(k) for k in (0.1, 0.3, 0.5, 1.0)}
	return o


def _variant(i):
	return variant(VARIANTS[i][0], p0=numpy.array(json.load(open(FIT))['p']), **VARIANTS[i][1])


def variants(out=RES, quiet=False):
	with _pool() as ex:
		V = list(ex.map(_variant, range(len(VARIANTS))))
	json.dump(V, open(os.path.join(out, 'unpol_variants.json'), 'w'), indent=1)
	h = '%-56s %8s %6s %8s %8s %8s %8s %7s %7s %7s %7s %7s %7s %7s %7s %6s' % ('variant', 'chi2/dof', 'S', 'Baldin', '+-fit*S', '+-pwa', '[E0,E1]', '[E1,E2]', '[E2,inf)', 'SR4', 'SR6', 'Ref.1', 'Ref.3', 'Ref.5', 'step@E1', 't')
	L = [h] + ['%-56s %8.3f %6.3f %8.4f %8.4f %8.4f %8.4f %7.4f %7.4f %7.4f %7.4f %7.3f %7.3f %7.3f %7.2f %6.3f' % ((v['name'][:56], v['chi2']/v['dof'], v['S']) + tuple(v['SR2']) + tuple(v['B_regions'])
		+ (v['SR4'][0], v['SR6'][0]) + tuple(v['Re_f'][k] for k in ('0.1', '0.3', '0.5')) + (v['step_E1'], v['t_pwa'])) + ('  ' + ','.join(v['problems']) if v['problems'] else '') for v in V]
	open(os.path.join(out, 'unpol_variants.txt'), 'w').write('# unpol2026.variants: sum rules in 1e-4 fm^(n+1), Re f in mcb GeV; +-fit*S includes the PWA mode t where the variant has one, +-pwa is the\n# independent PWA error otherwise; step@E1 = sigma_fit(E1) - sigma_PWA(E1) (1 + r t) [mcb]\n' + '\n'.join(L) + '\n')
	if not quiet: print('\n'.join(L))
	return V


if __name__ == '__main__':
	t0 = time.time()
	if 'fit' in sys.argv[1:] or not os.path.exists(FIT):
		fit()
	results()
	variants()
	print('run time %.0f s' % (time.time() - t0))
