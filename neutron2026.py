"""2026 update, neutron: forward sum rules I_n = 1/(2 pi^2) int sigma/nu^n, n = 2 (Baldin, alpha+beta), 4, 6, of the
neutron, from the single-pion partial waves and the 2026 proton evaluation (unpol2026.py).

There are no total photoabsorption data on a free neutron. The cross section is built as follows:
- below E1 (2 pi threshold): single-pion production is the whole cross section; average of SM22 and MAID2007 for
  gamma n -> pi- p + pi0 n, error 2 % (+) half their difference, as for the proton (unpol2026.py);
- above E1, up to 1.5 GeV: single-pion production from the same average, error 3 % (+) half the difference;
- the rest (multi-pion production, and single pion above 1.5 GeV) is taken from the proton: R_p = (2026 proton value
  above E1) - (single-pion proton part above E1, same average), assumed equal for the neutron, with 50 % error.
Check: the Baldin sum (n = 2) agrees with the evaluation of M. I. Levchuk and A. I. L'vov, Nucl. Phys. A 674, 449 (2000)
(alpha+beta = 15.2(5) for the neutron).  The neutron - proton differences are given with the errors combined in
quadrature (the remainder is then counted in both; proton: exp+fit (+) data selection).
Data: data/2026/single_pion_said.dat (tools/fetch_said_single_pion.py).  Output: results/2026/neutron_sumrules.json.
Units: 1e-4 fm^3 (n = 2), 1e-4 fm^5 (n = 4), 1e-4 fm^7 (n = 6)."""
import json
import numpy as np
from compton import E1, hc as _hc
hc = _hc.val	# fm GeV (a Value in compton.py)

DATA, PROTON, OUT = 'data/2026/single_pion_said.dat', 'results/2026/unpol_sumrules.json', 'results/2026/neutron_sumrules.json'
EMAX = 1.5	# GeV, end of the partial-wave table

def moment(k, s, n, lo=0., hi=np.inf):	# 1/(2 pi^2) int_lo^hi s/k^n dk [s in mcb, k in GeV] -> 10^-4 fm^(n+1); trapezoid, split point interpolated
	m = (k > lo) & (k < hi)
	x, y = k[m], s[m]
	if lo > k[0]: x, y = np.r_[lo, x], np.r_[np.interp(lo, k, s), y]
	if hi < k[-1]: x, y = np.r_[x, hi], np.r_[y, np.interp(hi, k, s)]
	return np.sum((y[1:]/x[1:]**n + y[:-1]/x[:-1]**n) * np.diff(x)) / 2 * 1e-4 * hc**(n - 1) / (2*np.pi**2) * 1e4

def single_pion(nucleon):	# total single-pion cross section [mcb] of SM22, MAID2007
	d = np.loadtxt(DATA)
	c = (1, 2) if nucleon == 'p' else (3, 4)
	return d[:, 0], d[:, c[0]] + d[:, c[1]], d[:, c[0] + 4] + d[:, c[1] + 4]

def avg(k, a, b, n, lo=0., hi=np.inf):	# average of the two analyses and half their difference
	x, y = moment(k, a, n, lo, hi), moment(k, b, n, lo, hi)
	return (x + y)/2, abs(x - y)/2

def results(out=OUT, quiet=False):
	PJ = json.load(open(PROTON))
	P = PJ['sum_rules']
	kn, sn, mn = single_pion('n')
	kp, sp, mp = single_pion('p')
	res = {}
	for n in (2, 4, 6):
		reg = P[str(n)]['regions']
		lo, dlo = avg(kn, sn, mn, n, hi=E1)
		hi, dhi = avg(kn, sn, mn, n, lo=E1, hi=EMAX)
		p_hi1pi, _ = avg(kp, sp, mp, n, lo=E1, hi=EMAX)
		R = reg['[E1,E2]'][0] + reg['[E2,inf)'][0] - p_hi1pi
		err = {'below_E1': float(np.hypot(0.02*lo, dlo)), 'above_E1_single_pion': float(np.hypot(0.03*hi, dhi)), 'remainder': 0.5*R}
		val = lo + hi + R
		tot = float(np.sqrt(sum(e**2 for e in err.values())))
		p, ep = reg['total'][0], float(np.hypot(reg['total'][1], PJ['selection_systematic'][str(n)]['half_difference']))
		res[str(n)] = dict(value=val, error=tot, budget=err,
			parts={'below_E1': lo, 'below_E1_SM22_MAID_half_difference': dlo, 'above_E1_single_pion': hi, 'remainder_from_proton': R},
			proton_2026=[p, ep], n_minus_p=[val - p, float(np.hypot(tot, ep))],
			single_pion_n_minus_p=avg(kn, sn, mn, n)[0] - avg(kp, sp, mp, n)[0])
		if not quiet:
			print('I%d neutron: below E1 %.3f + single pion above E1 %.3f + remainder %.3f = %.3f +- %.3f   (n - p = %.3f +- %.3f)'
				% (n, lo, hi, R, val, tot, val - p, np.hypot(tot, ep)))
	json.dump(dict(info=__doc__.split('\n\n')[0], date='2026-11-10', units={'2': '1e-4 fm^3', '4': '1e-4 fm^5', '6': '1e-4 fm^7'},
		compare={'2': {'Levchuk-Lvov 2000': [15.2, 0.5]}}, sum_rules=res), open(out, 'w'), indent=1)
	return res

if __name__ == '__main__':
	results()
