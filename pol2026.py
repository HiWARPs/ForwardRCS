#!/usr/bin/env python
# 2026 update of the polarized (helicity-difference) evaluation of GHP16 (Gryniuk, Hagelstein, Pascalutsa, PRD 94, 034043 (2016)).
#
# Convention as in pol.py: the cross-section column is (sigma_3/2 - sigma_1/2)/2 = Delta sigma/2 [mcb];
#   Im g = -nu col/(4 pi), I_GDH = 4 pi^2 integral_n(1) [mcb], gamma0 = -100 integral_n(3) [1e-6 fm^4], gamma0bar = -100 integral_n(5) [1e-6 fm^6].
# Regions (b = [E0, 0.2 GeV, E_j, inf]):
#   [E0, 0.2]   PWA table: average of the current public solutions MAID2007 (data/pol/maid.dat) and SAID SM22 (data/2026/pol_said_sm22.dat)
#               (1pi only, below the 2pi threshold), times (1 + r(E) t): one nuisance mode t (prior N(0, 1)), r = 2% (+) half the MAID-SM22
#               difference (relative), the same mode on the PWA side of the blend (data/2026/pol_threshold.dat); the SAID table of GHP16
#               (data/pol/said.dat = CM12) enters the PWA-model systematic (half the largest difference among MAID2007, SM22, CM12);
#   [0.2, E_j]  fit to data/2026/pol_db.csv (A2 2026, MAMI 2001, ELSA 2003/04) with one scale nuisance per experiment
#               (model-multiplied, penalty (n-1)^2/delta^2) and one additive nuisance for the A2 extrapolation correction;
#               joined continuously to the PWA by a cos^2 blend on [0.2, 0.2 + DB] (a step at 0.2 GeV would make Re g log-singular there);
#   [E_j, inf)  Regge: S. D. Bass, M. Skurzok, P. Moskal, PRC 98, 025209 (2018): sigma_A - sigma_P = (35.3 +- 3.6) s^(-0.69 +- 0.04) mcb
#               for sqrt(s) >= 2.5 GeV, i.e. in the AsymptoticRegionFit form 0.5*(C1 W^p1 + C2 W^p2): C1 = -35.3, p1 = -1.38, C2 = 0;
#               enters the fit as a Gaussian prior (correlation RHO assumed); joined to the fit with fit.Smooth at E_j.
# All errors are propagated through one covariance matrix (self.fit.pcov) of (model parameters, t, C1, p1, C2, p2), coherently across
# the regions: data, scale and A2-extrapolation nuisances (they act on the data, not on sigma: marginalized into the model block), Regge
# prior, PWA mode t. The PWA-model systematic is quoted separately.
# Run from the repository root:  python pol2026.py [fit] [results] [variants] [g]   (no argument: all four)

import os, sys, csv, json, time, numpy
from numpy import inf, pi
from scipy.optimize import least_squares
from scipy.interpolate import BSpline
from compton import E0, Mp, Mn, Mpi0, Mpip, IGDH, hc, fcm2lab
from fit import FunctionFit
from val import val
import integrate, pol

D26 = os.path.join('data', '2026')
R26 = os.path.join('results', '2026')
FIT_FILE = os.path.join(D26, 'pol_fit.json')
M, MPI, HC, SR_GDH = Mp.val, Mpi0.val, float(hc), float(IGDH)
B1 = 0.2	# GeV: PWA below, data above (no data below 0.2 GeV)
DB = 0.01	# GeV: width of the PWA-fit blend above B1
EJ = (2.5**2 - M**2)/(2*M)	# 2.8615 GeV: W = 2.5 GeV, lower end of the Bass et al. fit (and end of the ELSA data)
AJ = 0.1	# GeV: width of the Smooth junction (fit.Smooth, as in tot.py)
REGGE = [-35.3, -1.38, 0., 0.]	# C1 [mcb], p1, C2, p2 of AsymptoticRegionFit: 0.5*C1*W^p1 = (sigma_3/2 - sigma_1/2)/2 (Bass et al.)
RHO = 0.37	# ASSUMED corr(C1, p1) (not quoted by Bass et al.): reproduces their +-2 mcb on I_GDH(sqrt(s) > 2.5 GeV) (rho = 0 gives +-2.5)
def regge_cov(rho=RHO, e=(3.6, 0.08)):	# errors of C1 and p1 (= 2 x exponent: +-0.08); C2 = p2 = 0 fixed
	C = numpy.zeros((4, 4))
	C[:2, :2] = [[e[0]**2, rho*e[0]*e[1]], [rho*e[0]*e[1], e[1]**2]]
	return C
PWA_FRAC = 0.02	# baseline PWA normalization error (as GHP16), combined with half the MAID2007-SM22 spread
PWA_COLS = {'avg': 1, 'maid': 2, 'sm22': 3, 'cm12': 5}	# columns of data/2026/pol_threshold.dat

W = lambda k: numpy.sqrt(M*(M + 2*numpy.asarray(k, float)))
def qpi(w):	# pi N cm momentum
	return numpy.sqrt(numpy.clip((w**2 - (M + MPI)**2)*(w**2 - (M - MPI)**2), 0, None))/(2*w)


#---------------------------------------------------------------- PWA table and database

def pwa_table():	# data/2026/pol_threshold.dat on [E0, 0.22 GeV]: E, average, MAID2007, SM22, (MAID - SM22)/2, CM12 (all Delta sigma/2, mcb)
	global _T
	tab = lambda f: numpy.array(sorted([float(a) for a in l.split()[:2]] for l in open(os.path.join('data', 'pol', f)) if len(l.split()) > 1 and l.strip()[0] != '#')).T	# said.dat has wrapped 1-field lines
	m, c = tab('maid.dat'), tab('said.dat')
	s = numpy.loadtxt(os.path.join(D26, 'pol_said_sm22.dat'), usecols=(0, 5)).T
	m, s = numpy.c_[[E0, 0.], m], numpy.c_[[E0, 0.], s]	# the tables start above threshold
	x = numpy.unique(numpy.r_[E0, numpy.arange(0.1448, 0.2201, 0.0002), fcm2lab(Mn + Mpip)])	# 0.2 MeV grid plus the pi+ n threshold
	ym, ys, yc = numpy.interp(x, *m), numpy.interp(x, *s), numpy.interp(x, *c)
	numpy.savetxt(os.path.join(D26, 'pol_threshold.dat'), numpy.c_[x, (ym + ys)/2, ym, ys, (ym - ys)/2, yc], fmt='%.6f',
		header='PWA table [E0, 0.22 GeV] for pol2026.py: data/pol/maid.dat (MAID2007), data/2026/pol_said_sm22.dat (SAID SM22), data/pol/said.dat (SAID CM12: checked, see pol_said_sm22.dat)\n'
		'columns: E_gamma (GeV), Delta sigma/2 = (s_3/2 - s_1/2)/2 (mcb): average(MAID2007, SM22), MAID2007, SM22, half spread (MAID2007 - SM22)/2, CM12; linear interpolation of each table')
	_T = None

_T = None
def load_pwa():
	global _T
	if _T is None: _T = numpy.loadtxt(os.path.join(D26, 'pol_threshold.dat')).T
	return _T

def pwa_parts(k, col=1):	# PWA column a and its error mode d = r a, r = 2% (+) |MAID2007 - SM22|/(2|a|)
	t = load_pwa()
	a, h = numpy.interp(k, t[0], t[col]), numpy.interp(k, t[0], t[4])
	return a, a*numpy.sqrt(PWA_FRAC**2 + h**2/(a**2 + 1e-4))	# 1e-4 mcb^2: r regular where the column crosses 0 (|a| < 0.1 mcb, below the pi+ n threshold)

def pwa(k, t, col=1):	# PWA column times (1 + r t)
	a, d = pwa_parts(k, col)
	return a + t*d

def load_db(keep=lambda r: r['used'] == '1'):	# rows of data/2026/pol_db.csv in the code convention (column = Delta sigma/2)
	R = [r for r in csv.DictReader(open(os.path.join(D26, 'pol_db.csv'))) if keep(r)]
	a = lambda c, f=float: numpy.array([f(r[c]) for r in R])
	return dict(E=a('Egamma_GeV'), y=a('dsigma_ub')/2, s=a('stat_ub')/2, ex=a('extrap_ub')/2, g=a('norm_group', str), nf={r['norm_group']: float(r['norm_frac']) for r in R})


#---------------------------------------------------------------- parametrizations of the column above 0.2 GeV

XB = 0.165	# GeV: range of the centrifugal barrier factor (a choice; the spline variant checks the parametrization dependence)
def bw(w, A, MR, G, L):	# relativistic Breit-Wigner, width ~ q^(2L+1) with barrier factor; A = value at w = MR
	qr, q = qpi(MR), qpi(w)
	g = G*(q/qr)**(2*L + 1)*((qr**2 + XB**2)/(q**2 + XB**2))**L
	return A*MR**2*G*g/((w**2 - MR**2)**2 + MR**2*g**2)

L_BW = (1, 2, 0, 3)	# Delta(1232); second region (N(1520) + Delta pi, effective); N(1535); N(1680)
LN = numpy.log(B1), numpy.log(2.8615)	# background: Chebyshev series in ln(nu) mapped from [0.2, 2.8615] GeV to [-1, 1]
P_BW = [354.2, 1.2253, 0.1259, 182.1, 1.509, 0.323, -45.1, 1.5688, 0.088, 43.2, 1.6649, 0.070, -75.9, 49.4, 42.4, -28.8, -12.7, 12.9]	# start: the best interior minimum of a 40-start search
BOUNDS_BW = ([-inf, 1.18, 0.06, -inf, 1.44, 0.05, -inf, 1.50, 0.03, -inf, 1.62, 0.05] + [-inf]*6, [inf, 1.28, 0.2, inf, 1.58, 0.6, inf, 1.62, 0.3, inf, 1.75, 0.25] + [inf]*6)	# guard rails only: the solution must be interior
KNOTS = [0.2, 0.24, 0.27, 0.3, 0.33, 0.36, 0.4, 0.45, 0.5, 0.56, 0.62, 0.68, 0.74, 0.8, 0.86, 0.92, 0.98, 1.05, 1.15, 1.3, 1.5, 1.8, 2.3]	# + E_j

def model(kind, k, p, knots=None):	# the fitted column (mcb), vectorized in k
	if kind == 'bw':
		w = W(k)
		return sum(bw(w, *p[3*i:3*i+3], L) for i, L in enumerate(L_BW)) + numpy.polynomial.chebyshev.chebval(2*(numpy.log(k) - LN[0])/(LN[1] - LN[0]) - 1, p[12:])
	t = numpy.log(knots)	# cubic B-spline in ln(nu)
	return BSpline(numpy.r_[[t[0]]*3, t, [t[-1]]*3], numpy.asarray(p), 3)(numpy.log(k))

def blend(k, f, t, db=DB, col=1):	# column continuous at B1: cos^2 weight from PWA (k = B1) to fit (k >= B1 + db)
	w = numpy.cos(pi/2*numpy.clip((k - B1)/db, 0, 1))**2 if db > 0 else 0.
	return w*pwa(k, t, col) + (1 - w)*f

def regge(k, p):
	w = W(k)
	return 0.5*(p[0]*w**p[1] + p[2]*w**p[3])

def smooth(k, f1, f2, m, a):	# = fit.Smooth
	e = 1/(1 + numpy.exp(numpy.clip((k - m)/a, -700, 700)))
	return e*f1 + (1 - e)*f2


class Model2026(FunctionFit):	# the fitted-region column as a legacy FunctionFit (ResonanceRegion.fit); parameters: model, PWA mode t

	def __init__(self, F, col):
		self.kind, self.knots, self.db, self.col = F['model'], F.get('knots'), F['blend'], col
		FunctionFit.__init__(self, *F['p'][:-4])

	def f(self, k, *p):
		return 0. if k > 50. else float(blend(k, model(self.kind, k, p[:-1], self.knots), p[-1], self.db, self.col))	# shortcut: above E_j + 0.3 GeV the Smooth weight is < 5e-21; avoids inf*0 far out


#---------------------------------------------------------------- integrals (vectorized; the legacy machinery is the check)

def _nodes(a, b, n=48):	# Gauss-Legendre panels on [a, b]; b = inf: panels in ln(nu) up to e^19 GeV
	x, w = numpy.polynomial.legendre.leggauss(n)
	if b == inf:
		e = numpy.r_[numpy.log(a), numpy.arange(numpy.ceil(numpy.log(a)), 19.)]
		u = numpy.concatenate([(e1 - e0)/2*x + (e1 + e0)/2 for e0, e1 in zip(e[:-1], e[1:])])
		return numpy.exp(u), numpy.concatenate([(e1 - e0)/2*w for e0, e1 in zip(e[:-1], e[1:])])*numpy.exp(u)
	e = numpy.unique(numpy.clip(numpy.r_[B1 + DB, numpy.arange(0.2, 0.5, 0.025), numpy.arange(0.5, 1.4, 0.05), numpy.arange(1.4, 3.6, 0.1), a, b], a, b))
	return numpy.concatenate([(e1 - e0)/2*x + (e1 + e0)/2 for e0, e1 in zip(e[:-1], e[1:])]), numpy.concatenate([(e1 - e0)/2*w for e0, e1 in zip(e[:-1], e[1:])])

UNIT = {1: 2., 3: -100*HC**2/(2*pi**2), 5: -100*HC**4/(2*pi**2)}	# I_GDH [mcb], gamma0 [1e-6 fm^4], gamma0bar [1e-6 fm^6] from int col/nu^n
NAMES = {1: 'I_GDH', 3: 'gamma0', 5: 'gamma0bar'}
REGIONS = ('[E0,0.2]', '[0.2,Ej]', '[Ej,inf)')

def moments(f, a, b):	# {n: UNIT[n] * int_a^b f(nu)/nu^n dnu}
	k, w = _nodes(a, b)
	fk = w*f(k)
	return {n: UNIT[n]*numpy.sum(fk/k**n) for n in UNIT}

def low_moments(col=1):	# trapezoid on the table nodes of [E0, 0.2] (as integrate.ThresholdRegion): central column ('avg'), its error mode, MAID2007, SM22, CM12
	t = load_pwa()
	x = numpy.r_[t[0][t[0] < B1], B1]
	I = lambda y: {n: UNIT[n]*numpy.trapezoid(y/x**n, x) for n in UNIT}
	a, d = pwa_parts(x, col)
	return dict(avg=I(a), mode=I(d), **{k: I(numpy.interp(x, t[0], t[c])) for k, c in PWA_COLS.items() if k != 'avg'})

def pwa_model_syst(low):	# PWA-model systematic: half the largest difference among MAID2007, SM22, CM12 on [E0, 0.2]
	return {n: max(abs(low[a][n] - low[b][n]) for a in ('maid', 'sm22', 'cm12') for b in ('maid', 'sm22', 'cm12'))/2 for n in UNIT}


#---------------------------------------------------------------- fit

def fit_data(kind='bw', D=None, ej=EJ, aj=AJ, db=DB, regge0=REGGE, rcov=None, gdh=False, p0=None, groups_map=None, ex_nuisance=True, pcol=1):
	"""Fit n_d * col(E) + eta * extrap(E) to the data, col = Smooth(blend(model, PWA), Regge). Nuisances: n_d (penalty (n_d - 1)/delta_d),
	eta (A2 extrapolation, unit prior), PWA mode t (unit prior; it acts on [E0, 0.2] and the blend, but is not fitted to the data in the
	blend; only the GDH sum rule as a pseudo-measurement, gdh=True, moves it); Regge as a Gaussian prior. pcol: central PWA table.
	Returns the solution, the whitened Jacobian and its row groups."""
	D = D or load_db()
	low = low_moments(pcol)
	E, y, s, ex, g = D['E'], D['y'], D['s'], D['ex'], D['g'].copy()
	if groups_map: g = numpy.array([groups_map.get(a, a) for a in g])
	if not ex_nuisance: s, ex = numpy.hypot(s, ex), 0*ex	# variant: extrapolation error point-to-point
	G = sorted(set(g)); gi = numpy.array([G.index(a) for a in g])
	nf = numpy.array([max(f for a, f in D['nf'].items() if (groups_map or {}).get(a, a) == b) for b in G])	# merged groups: the larger scale error
	knots = [x for x in KNOTS if x < ej - 0.1] + [ej]
	m = len(P_BW) if kind == 'bw' else len(knots) + 2
	rcov = regge_cov() if rcov is None else numpy.array(rcov)
	fr = [i for i in range(4) if rcov[i, i] > 0]; nr = len(fr)	# free Regge parameters (C2 = p2 = 0 fixed for Bass et al.)
	Lr = numpy.linalg.cholesky(numpy.linalg.inv(rcov[numpy.ix_(fr, fr)])).T
	K, Kw = map(numpy.concatenate, zip(_nodes(B1, ej), _nodes(ej, inf)))
	ng = len(G)
	def phys(P):	# model parameters, PWA mode t, Regge parameters
		rp = numpy.array(regge0, float); rp[fr] = P[m:m+nr]
		return P[:m], P[-1], rp
	def col(k, P, pwa_free=True):	# pwa_free=False: t = 0 (data rows: the data do not calibrate the PWA)
		p, t, rp = phys(P)
		return smooth(k, blend(k, model(kind, k, p, knots), t if pwa_free else 0., db, pcol), regge(k, rp), ej, aj)
	def gdh_val(P):	# I_GDH total = PWA [E0, 0.2] (with its mode) + [0.2, inf)
		return low['avg'][1] + P[-1]*low['mode'][1] + 2*numpy.sum(Kw*col(K, P)/K)
	def resid(P):
		nrm, eta = P[m+nr:m+nr+ng], P[m+nr+ng]
		r = [(nrm[gi]*col(E, P, False) + eta*ex - y)/s, (nrm - 1)/nf, [eta], Lr.dot(P[m:m+nr] - numpy.array(regge0)[fr]), P[-1:]]
		if gdh: r.append([(gdh_val(P) - SR_GDH)/0.01])	# the sum-rule value is exact for this purpose
		return numpy.concatenate(r)
	lo, hi = BOUNDS_BW if kind == 'bw' else ([-inf]*m, [inf]*m)
	P0 = numpy.r_[(P_BW if kind == 'bw' else numpy.zeros(m)) if p0 is None else p0, numpy.array(regge0)[fr], numpy.ones(ng), 0., 0.]	# the spline is linear: any start
	LO, HI = numpy.r_[lo, [-inf]*(nr + ng + 2)], numpy.r_[hi, [inf]*(nr + ng + 2)]
	sol = least_squares(resid, numpy.clip(P0, LO + 1e-9, HI - 1e-9), bounds=(LO, HI), method='trf', jac='3-point', x_scale='jac', xtol=1e-12, ftol=1e-12, max_nfev=5000)
	P, J, nd = sol.x, sol.jac, len(E)
	rows = dict(data=slice(0, nd), norm=slice(nd, nd + ng + 1), regge=slice(nd + ng + 1, nd + ng + 1 + nr), pwa=slice(nd + ng + 1 + nr, nd + ng + 2 + nr))
	if gdh: rows['gdh'] = slice(nd + ng + 2 + nr, nd + ng + 3 + nr)
	chi2 = numpy.sum(sol.fun**2); dof = nd - m + (1 if gdh else 0)	# nuisance and prior rows balance their parameters
	return dict(kind=kind, P=P, J=J, rows=rows, m=m, fr=fr, nr=nr, G=G, nf=nf, knots=knots, ej=ej, aj=aj, db=db, col=col, phys=phys, low=low, gdh=gdh, pcol=pcol,
		chi2=chi2, chi2_data=numpy.sum(sol.fun[:nd]**2), chi2_sets={a: float(numpy.sum(sol.fun[:nd][g == a]**2)) for a in G}, n_sets={a: int(numpy.sum(g == a)) for a in G},
		nd=nd, dof=dof, S=max(1., numpy.sqrt(chi2/dof)), at_bound=[i for i in range(m) if min(P[i] - LO[i], HI[i] - P[i]) < 1e-6], status=sol.status)

def analyse(F):
	"""Sum rules with the error budget: variance split by row group, exact for the linearized fit (Cov = A sum_G J_G^T J_G A, A = (J^T J)^-1);
	the scale factor S = sqrt(chi2/dof), if > 1, multiplies the non-PWA part."""
	P, J, m, nr = F['P'], F['J'], F['m'], F['nr']
	A = numpy.linalg.inv(J.T.dot(J))
	ej, low = F['ej'], F['low']
	def obs(P):	# {n: [PWA region, [0.2, E_j], [E_j, inf)]}
		mo = [moments(lambda k: F['col'](k, P), a, b) for a, b in ((B1, ej), (ej, inf))]
		return {n: [low['avg'][n] + P[-1]*low['mode'][n]] + [x[n] for x in mo] for n in UNIT}
	O = obs(P)
	grad = {n: numpy.zeros((3, len(P))) for n in UNIT}
	for i in list(range(m + nr)) + [len(P) - 1]:	# the scale and extrapolation nuisances do not enter the observables
		h = 1e-5*max(1., abs(P[i]))
		up, dn = P.copy(), P.copy(); up[i] += h; dn[i] -= h
		Ou, Od = obs(up), obs(dn)
		for n in UNIT: grad[n][:, i] = (numpy.array(Ou[n]) - numpy.array(Od[n]))/(2*h)
	def budget(gv):
		b = {G: float(numpy.linalg.norm(J[F['rows'][G]].dot(A.dot(gv)))) for G in F['rows']}
		fit = numpy.sqrt(sum(b[G]**2 for G in b if G != 'pwa'))
		b['S'] = float(fit*numpy.sqrt(F['S']**2 - 1))
		b['total'] = float(numpy.sqrt(F['S']**2*fit**2 + b['pwa']**2))
		return b
	out = {}
	for n in UNIT:
		b = budget(grad[n].sum(0))
		if not F['gdh']: b.update(pwa_low=float(abs(low['mode'][n])))	# the [E0, 0.2] part of 'pwa' (prior)
		out[NAMES[n]] = dict(value=float(sum(O[n])), err=b['total'], budget=b, regions={r: [float(O[n][i]), budget(grad[n][i])['total']] for i, r in enumerate(REGIONS)})
	return out

def physics_cov(F):	# covariance of (model parameters, PWA mode t, C1, p1, C2, p2), marginalized over the data nuisances, times S^2
	A = numpy.linalg.inv(F['J'].T.dot(F['J']))
	m, n = F['m'], len(F['P'])
	idx = list(range(m)) + [n - 1] + [m + j for j in range(F['nr'])]
	pos = list(range(m + 1)) + [m + 1 + i for i in F['fr']]
	C = numpy.zeros((m + 5, m + 5))
	C[numpy.ix_(pos, pos)] = A[numpy.ix_(idx, idx)]
	return C*F['S']**2

def entry(F):	# the JSON entry read by WholeRegionPol2026
	p, t, rp = F['phys'](F['P'])
	C = physics_cov(F); sd = numpy.sqrt(numpy.diag(C)); sd[sd == 0] = 1
	R = C/numpy.outer(sd, sd)
	A = numpy.linalg.inv(F['J'].T.dot(F['J'])); i0 = F['m'] + F['nr']
	return dict(model=F['kind'], p=list(map(float, p)) + [float(t)] + list(map(float, rp)), pcov=C.tolist(), b=[E0, B1, F['ej'], 'inf'], a=F['aj'], blend=F['db'], pwa_col=F['pcol'],
		knots=F['knots'] if F['kind'] == 'spline' else None,
		chi2=float(F['chi2']), chi2_data=float(F['chi2_data']), ndata=F['nd'], dof=F['dof'], S=float(F['S']), chi2_sets=F['chi2_sets'], n_sets=F['n_sets'],
		nuisances=dict({G: [float(F['P'][i0 + i]), float(numpy.sqrt(A[i0 + i, i0 + i])), float(F['nf'][i])] for i, G in enumerate(F['G'])},
			eta_extrap=[float(F['P'][i0 + len(F['G'])]), float(numpy.sqrt(A[i0 + len(F['G']), i0 + len(F['G'])]))], pwa_mode=[float(t), float(sd[F['m']])]),
		max_corr=float(numpy.max(numpy.abs(R - numpy.eye(len(C)))[:F['m'], :F['m']])), cond_corr=float(numpy.linalg.cond(R[:F['m'], :F['m']])), at_bound=F['at_bound'], gdh_constrained=F['gdh'])

def fit(quiet=False):	# builds the PWA table, fits, writes data/2026/pol_fit.json; returns the fits
	pwa_table()
	F = {'bw': fit_data('bw'), 'spline': fit_data('spline'), 'bw_gdh': fit_data('bw', gdh=True)}
	json.dump({k: entry(v) for k, v in F.items()}, open(FIT_FILE, 'w'), indent=1)
	if not quiet:
		for k, v in F.items():
			print('%-7s chi2/dof = %.1f/%d = %.3f (data %.1f for %d points), S = %.3f, at bound: %s' % (k, v['chi2'], v['dof'], v['chi2']/v['dof'], v['chi2_data'], v['nd'], v['S'], v['at_bound']))
	return F


#---------------------------------------------------------------- legacy machinery

class ThresholdRegion2026(integrate.ThresholdRegion):	# the PWA column with its mode t at the fitted value (0 unless the GDH sum rule is imposed); no own error: t is in pcov

	def __init__(self, t, col, *args, **kwargs):
		self.t, self.col = t, col
		integrate.ThresholdRegion.__init__(self, *args, **kwargs)

	def fit(self, k):
		return 0. if k < self.k_min else float(pwa(k, self.t, self.col))

	def _int_err(self, *args, **kwargs):
		return 0.


class WholeRegionPol2026(pol.WholeRegion):	# integral_n, re_g, im_g, g, dg, cs, print_all work as for pol.WholeRegion

	def __init__(self, model='bw', pwa=None, quiet=True):	# model: 'bw' (baseline), 'spline', 'bw_gdh' (GDH sum rule imposed); pwa: other central PWA table ('maid', 'sm22', 'cm12'), same parameters
		F = json.load(open(FIT_FILE))[model]
		col = PWA_COLS[pwa] if pwa else F['pwa_col']
		integrate.WholeRegion.__init__(self, b=[E0, B1, F['b'][2], inf], name='2026', path='2026', filename='pol_threshold.dat', init=False, quiet=quiet, pre_f=False,
			ThresholdRegion=lambda *a, **kw: ThresholdRegion2026(F['p'][-5], col, *a, **kw), ResonanceRegion=lambda k_min, k_max, **kw: integrate.Region(Model2026(F, col), k_min, k_max, 'resonances', False))
		self.h_region.fit.set_p(*F['p'][-4:])
		self.smooth.set_p(F['a'])
		self.fit.pcov = numpy.array(F['pcov'])	# model parameters, PWA mode t, C1, p1, C2, p2
		self.i_pwa, self.col = len(F['p']) - 5, col
		self.baldin = self.sr_2 = self.sr_4 = self.sr_6 = None
		self.f = self.re_f = self.im_f = self.integral_f = None
		self.F = F

	# Errors: all through self.fit.pcov, coherently across b[1] (the PWA mode acts on both sides of it); the legacy sum of
	# region Values would add the PWA region in quadrature and, at k = b[1], multiply 0 by an infinite log (integrate.Region._int).

	@integrate.cached_args_method()
	def dcs(self, i, k):
		if k > self.b[1]:
			return self.fitted_region.fit.df(i, k)
		return float(pwa_parts(k, self.col)[1]) if i == self.i_pwa and k >= self.b[0] else 0.

	def _err(self, d):
		return float(numpy.sqrt(numpy.dot(d, self.fit.pcov.dot(d))))

	def _dvec(self, d):	# derivative of an integral w.r.t. each parameter: fitted region, plus the PWA region for the PWA mode
		C = self.fit.pcov
		return numpy.array([0. if C[i, i] == 0 else d(self.fitted_region, i) + (d(self.data_region, i) if i == self.i_pwa else 0.) for i in range(len(C))])

	@integrate.cached_args_method()
	def cs_err(self, k):
		return self._err(numpy.array([self.dcs(i, k) for i in range(len(self.fit.pcov))]))

	@integrate.cached_args_method()
	def integral_n(self, n, k0, k1, err=True):
		r = integrate.WholeRegion.integral_n(self, n, k0, k1, err=False)
		if not err:
			return r
		d = lambda reg, i: float(reg.integral(n, 0, k0, k1, err=False, cs=lambda a: self.dcs(i, a), subtracted=False, const=True)) if max(k0, reg.k_min) < min(k1, reg.k_max) else 0.
		return val(float(r), HC**(n - 1)*self._err(self._dvec(d)))

	@integrate.cached_args_method()
	def _dg(self, n, i, k):	# d/dp_i of the g integrals: -k^n/(2 pi^2) PV int dcs(i, p) p^(2-n)/(p^2 - k^2) dp over [E0, inf)
		if self.fit.pcov[i, i] == 0:
			return 0.
		cs = lambda a: self.dcs(i, a)
		r = float(self.fitted_region.integral(n, k, err=False, cs=cs))	# subtracted at p = k; singular parts set to 0 in both regions (WholeRegion) ...
		if i == self.i_pwa:
			r += float(integrate.Region.integral(self.data_region, n, k, err=False, cs=cs))
		return -r - integrate.Region.singular(k, self.b[0], self.b[-1])*cs(k)*k**2/(2*pi**2)	# ... and taken here over [E0, inf)

	def dintegral_g(self, i, k, k0=None, k1=None, subtracted=True):
		return self._dg(1, i, k)

	def dintegral_g_subtracted(self, i, k, k0=None, k1=None, subtracted=True):
		return self._dg(3, i, k)

	def integral_g(self, k, k0=None, k1=None, err=True):
		r = float(integrate.WholeRegion.integral_g(self, k, k0, k1, err=False))
		return val(r, self._err(numpy.array([self._dg(1, i, k) for i in range(len(self.fit.pcov))]))) if err else r

	def integral_g_subtracted(self, k, k0=None, k1=None, err=True):
		r = float(integrate.WholeRegion.integral_g_subtracted(self, k, k0, k1, err=False))
		return val(r, self._err(numpy.array([self._dg(3, i, k) for i in range(len(self.fit.pcov))]))) if err else r


#---------------------------------------------------------------- results

def _pm(v):
	return '%.2f(%.2f)' % ((v['value'], v['err']) if isinstance(v, dict) else tuple(v))

def results(F=None):
	F = F or {'bw': fit_data('bw'), 'spline': fit_data('spline'), 'bw_gdh': fit_data('bw', gdh=True)}
	t0 = time.time()
	syst = {'unconstrained': pwa_model_syst(F['bw']['low'])}	# PWA-model systematic: half the largest difference among MAID2007, SM22, CM12
	g3 = [analyse(fit_data('bw', gdh=True, pcol=PWA_COLS[a])) for a in ('maid', 'sm22', 'cm12')]	# with the sum rule imposed: refits with each table as central PWA
	syst['gdh_constrained'] = {n: max(abs(x[NAMES[n]]['value'] - y[NAMES[n]]['value']) for x in g3 for y in g3)/2 for n in UNIT}
	w = WholeRegionPol2026('bw')
	b, u = w.b, {1: 4*pi**2, 3: -100, 5: -100}	# integral_n -> I_GDH, gamma0, gamma0bar
	legacy = {NAMES[n]: {r: [u[n]*v.val, abs(u[n])*v.err] for r, v in zip(REGIONS + ('total',), [w.integral_n(n, b[i], b[j]) for i, j in ((0, 1), (1, 2), (2, 3), (0, 3))])} for n in UNIT}
	t1 = time.time() - t0
	low = F['bw']['low']
	lex = {k: [float(w.re_g(k).val), -k*SR_GDH/(4*pi**2) + k**3*legacy['gamma0']['total'][0]/(100*HC**2) + k**5*legacy['gamma0bar']['total'][0]/(100*HC**4)] for k in (0.02, 0.05)}	# check: Re g vs its low-energy expansion
	R = dict(
		conventions='column = (sigma_3/2 - sigma_1/2)/2 [mcb]; I_GDH [mcb], gamma0 [1e-6 fm^4], gamma0bar [1e-6 fm^6]; regions [E0, 0.2 GeV] (PWA), [0.2 GeV, E_j] (fit; blend with the PWA on [0.2, %.2f]), [E_j, inf) (Regge, Bass et al. 2018), E_j = %.4f GeV, Smooth width %.2f GeV' % (B1 + DB, EJ, AJ),
		error_budget='data: statistical errors of the data; norm: scale nuisances and the A2 extrapolation nuisance; regge: Bass et al. prior (corr(C1, p1) = %.2f assumed); pwa: PWA mode t (column x (1 + r t), r = 2%% (+) half the MAID2007-SM22 spread, on [E0, 0.2] and the PWA side of the blend; pwa_low: its [E0, 0.2] part); gdh: the sum-rule pseudo-measurement; S: scale factor sqrt(chi2/dof) if > 1; total = (data (+) norm (+) regge (+) gdh) x S (+) pwa = the exp+fit error; pwa_model_systematic: separate' % RHO,
		unconstrained=analyse(F['bw']), gdh_constrained=analyse(F['bw_gdh']), spline=analyse(F['spline']), legacy_machinery=legacy,
		pwa_low={k: {NAMES[n]: float(v[n]) for n in UNIT} for k, v in low.items()},
		fit={k: {a: c for a, c in entry(v).items() if a not in ('p', 'pcov')} for k, v in F.items()},
		step_at_0p2={k: float(model(k, numpy.array([B1]), v['phys'](v['P'])[0], v['knots'])[0] - pwa(B1, 0.)) for k, v in ((k, fit_data(k, db=0.)) for k in ('bw', 'spline'))},	# fit without the blend minus PWA at 0.2 GeV
		check_re_g_low_energy_expansion=lex, gdh_sum_rule=SR_GDH,
		pwa_model_systematic=dict(definition='half the largest difference among MAID2007, SAID SM22, SAID CM12 below 0.2 GeV: [E0, 0.2] contributions (unconstrained); refits with each table as central PWA (GDH-constrained)',
			**{c: {NAMES[n]: float(v) for n, v in syst[c].items()} for c in syst}),
		comparison=dict(GHP16=dict(I_GDH=[204.5, 21.4], gamma0=[-92.9, 10.5], gamma0bar=[48.4, 8.2]), A2_2026=dict(I_GDH=[210, 2, 18], below_0p2=[-30, 2], range_0p2_2p9=[255, 2, 18], above_2p9=[-15, 2]),
			audit=dict(gamma0='about -104 to -107 +- 12 (MAMI data from 0.2 GeV, PWA spread below); A2 2026 instead of MAID on 0.2-E1: gamma0 -10, I_GDH +7')))
	R['headline'] = {c: {NAMES[n]: [R[c][NAMES[n]]['value'], R[c][NAMES[n]]['err'], R['pwa_model_systematic'][c][NAMES[n]]] for n in UNIT} for c in ('unconstrained', 'gdh_constrained')}	# value, exp+fit, PWA model
	for n in UNIT:	# self-checks: legacy machinery (mpmath quadrature, pcov) = vectorized integration (row-group budget); Re g = its low-energy expansion
		a, l = R['unconstrained'][NAMES[n]], legacy[NAMES[n]]['total']
		assert abs(a['value'] - l[0]) < 1e-4*max(1, abs(l[0])) and abs(a['err'] - l[1]) < 1e-3*l[1], (NAMES[n], a['value'], a['err'], l)
	assert all(abs(a - b) < 1e-4 for a, b in lex.values()), lex
	os.makedirs(R26, exist_ok=True)
	json.dump(R, open(os.path.join(R26, 'pol_sumrules.json'), 'w'), indent=1)
	for k in ('unconstrained', 'gdh_constrained'):
		print('%-16s' % k, '  '.join('%s %.2f +- %.2f (exp+fit) +- %.2f (PWA model)' % ((a,) + tuple(R['headline'][k][a])) for a in NAMES.values()))
	print('%-16s' % 'spline', '  '.join('%s %s' % (a, _pm(R['spline'][a])) for a in NAMES.values()))
	print('%-16s' % 'legacy', '  '.join('%s %s' % (a, _pm(R['legacy_machinery'][a]['total'])) for a in NAMES.values()), '(%.0f s)' % t1)
	return R

def _re_im_g(k):
	return (0., 0., 0., 0.) if k == 0 else tuple(_W.re_g(k)) + tuple(_W.im_g(k))

def _re_im_g_val(a):	# values only, with another central PWA table below 0.2 GeV
	name, k = a
	return float(_WS[name].integral_g_subtracted(k, err=False)), -float(_WS[name].cs(k, err=False))*k/(4*pi)

def g_table(nproc=None):	# results/2026/pol_g.txt with the legacy machinery (format of the APS supplement g.txt), in parallel
	global _W, _WS
	_W = WholeRegionPol2026('bw')
	_WS = {a: WholeRegionPol2026('bw', pwa=a) for a in ('maid', 'sm22', 'cm12')}
	kk = [round(0.005*i, 3) for i in range(501)] + [3., 4., 5., 7.5, 10., 20., 50., 100.]
	ks = kk[1:]	# PWA-model systematic of g: all energies (the dispersive tail of a change below 0.2 GeV extends to high nu)
	import multiprocessing
	with multiprocessing.get_context('fork').Pool(nproc or min(12, os.cpu_count())) as p:
		rows = p.map(_re_im_g, kk, chunksize=2)
		alt = p.map(_re_im_g_val, [(a, k) for a in _WS for k in ks], chunksize=2)
	alt = numpy.array(alt).reshape(3, len(ks), 2)
	half = lambda x: (x.max(0) - x.min(0))/2
	numpy.savetxt(os.path.join(R26, 'pol_g_pwa_model.txt'), numpy.c_[ks, half(alt[:, :, 0]), half(alt[:, :, 1]), alt[:, :, 0].T], fmt='%.3f' + '   %.8f'*5,
		header='PWA-model systematic of g(nu) (pol2026.py, model "bw"): half the largest difference among MAID2007, SAID SM22, SAID CM12\n'
		'as central PWA below 0.2 GeV (and on the PWA side of the 0.2-0.21 GeV blend), fit parameters unchanged (Im g: zero above 0.21 GeV)\n'
		'nu (GeV)   dRe g   dIm g   Re g(MAID2007)   Re g(SM22)   Re g(CM12)   (mcb*GeV)')
	os.makedirs(R26, exist_ok=True)
	with open(os.path.join(R26, 'pol_g.txt'), 'w') as f:
		f.write('# forward spin-dependent amplitude g(nu) of the proton: 2026 update of GHP16 (pol2026.py, model "bw"); GDH-subtracted DR, i.e. the LET g\'(0) = -alpha kappa^2/(2 M^2)\n')
		f.write('# errors (exp+fit): one covariance (data, scale and extrapolation nuisances, Regge prior, PWA mode below 0.2 GeV); PWA-model systematic: pol_g_pwa_model.txt; nu > E_j = 2.86 GeV: Regge region (Bass et al. 2018)\n')
		f.write(r' \nu       Re g       error      Im g       error' + '\n')
		f.write(r'(GeV)    (mcb*GeV)  (mcb*GeV)  (mcb*GeV)  (mcb*GeV)' + '\n\n')
		f.write('\n'.join('%.3f   %.8f   %.8f   %.8f   %.8f' % ((k,) + r) for k, r in zip(kk, rows)) + '\n')


#---------------------------------------------------------------- variants (model dependence)

def variants():
	low = low_moments()
	base = analyse(fit_data('bw'))
	keep = lambda f: (lambda r: r['used'] == '1' and f(r))
	ghp_elsa = lambda r: r['dataset'] == 'ELSA04' or (r['used'] == '1' and (r['dataset'] != 'ELSA03' or float(r['Egamma_GeV']) < 1.25))	# GHP16: ELSA04 from 1.3 GeV
	V = {
		'spline (second parametrization)': dict(kind='spline'),
		'junction E_j = 2.5 GeV': dict(ej=2.5), 'junction E_j = 3.2 GeV': dict(ej=3.2),
		'Smooth width 0.05 GeV': dict(aj=0.05), 'Smooth width 0.2 GeV': dict(aj=0.2),
		'PWA-fit blend 5 MeV': dict(db=0.005), 'PWA-fit blend 20 MeV': dict(db=0.02), 'no blend (step at 0.2 GeV)': dict(db=0.),
		'without A2 2026': dict(D=load_db(keep(lambda r: r['dataset'] != 'A2_26'))),
		'without MAMI 2001': dict(D=load_db(keep(lambda r: r['dataset'] != 'MAMI01'))),
		'ELSA as GHP16 (ELSA04 at 1.3-1.74 GeV instead of ELSA03 at 1.29-1.80)': dict(D=load_db(ghp_elsa)),
		'common scale MAMI 2001 + ELSA (same GDH-collaboration target)': dict(groups_map={'MAMI01': 'GDH', 'ELSA': 'GDH'}),
		'common scale A2 2026 + MAMI 2001 (same tagged-photon facility)': dict(groups_map={'MAMI01': 'MAMI', 'A2_26': 'MAMI'}),
		'A2 extrapolation error point-to-point': dict(ex_nuisance=False),
		'Regge corr(C1,p1) = 0': dict(rcov=regge_cov(0.)), 'Regge corr(C1,p1) = 0.9': dict(rcov=regge_cov(0.9)),
		'Regge of GHP16 (Helbing 2006, diagonal errors)': dict(regge0=[-34.1, -1.16, 209.4, -3.32], rcov=numpy.diag([5.7, 0.46, 29., 0.44])**2),
	}
	row = lambda r: {a: [r[a]['value'], r[a]['err']] + [r[a]['regions'][x][0] for x in REGIONS] for a in NAMES.values()}
	out = {'baseline': row(base)}
	for name, kw in V.items():
		F = fit_data(**kw)
		out[name] = row(analyse(F))
		e = entry(F)
		out[name]['fit'] = dict(chi2_data=e['chi2_data'], ndata=e['ndata'], dof=e['dof'], S=e['S'], at_bound=e['at_bound'], cond_corr=e['cond_corr'], max_corr=e['max_corr'])
	for name, tab in (('PWA below 0.2 GeV: MAID2007 only', 'maid'), ('PWA below 0.2 GeV: SAID SM22 only', 'sm22'), ('PWA below 0.2 GeV: SAID CM12 only (GHP16 SAID table)', 'cm12')):	# no refit: the data start at 0.2 GeV (the 10-MeV blend is neglected here)
		out[name] = {NAMES[n]: [base[NAMES[n]]['value'] + low[tab][n] - low['avg'][n], None, low[tab][n], None, None] for n in UNIT}
	name = 'PWA below 0.2 GeV: average MAID2007 + CM12 (the GHP16 tables)'
	out[name] = {NAMES[n]: [base[NAMES[n]]['value'] + (low['maid'][n] + low['cm12'][n])/2 - low['avg'][n], None, (low['maid'][n] + low['cm12'][n])/2, None, None] for n in UNIT}
	r = 30./-low['avg'][1]	# ESTIMATE: A2 2026 quote -30 +- 2 mcb below 0.2 GeV (mean of SM22 and BnGa-2019); same shape assumed (gamma0/I_GDH of [E0, 0.2]: MAID2007 -3.26, SM22 -3.28, CM12 -3.28)
	out['PWA below 0.2 GeV scaled to I_GDH = -30 (SM22/BnGa-2019 per A2 2026; ESTIMATE)'] = {NAMES[n]: [base[NAMES[n]]['value'] + (r - 1)*low['avg'][n], None, r*low['avg'][n], None, None] for n in UNIT}
	os.makedirs(R26, exist_ok=True)
	json.dump(out, open(os.path.join(R26, 'pol_variants.json'), 'w'), indent=1)
	print('%-75s %18s %18s %18s' % ('variant', 'I_GDH', 'gamma0', 'gamma0bar'))
	for k, v in out.items():
		print('%-75s' % k[:75], ' '.join('%8.2f(%5s) %+6.2f' % (v[a][0], '%.2f' % v[a][1] if v[a][1] else '-', v[a][0] - out['baseline'][a][0]) for a in NAMES.values()), v.get('fit', ''))
	return out


if __name__ == '__main__':

	t0 = time.time()
	todo = sys.argv[1:] or ['fit', 'results', 'variants', 'g']
	F = fit() if 'fit' in todo else None
	if 'results' in todo: results(F)
	if 'variants' in todo: variants()
	if 'g' in todo: g_table()
	print('run time %.0f s' % (time.time() - t0))
