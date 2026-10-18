#!/usr/bin/env python3
"""Regression test.  Run from the repository root:  python -m unittest discover -s tests   (2-3 minutes)

- pins the current sum rules, amplitudes and forward observables;
- checks them against the published tables of Gryniuk, Hagelstein, Pascalutsa, PRD 92, 074031 (2015) and
  PRD 94, 034043 (2016), within the quoted errors;
- pins the sum rules and two amplitudes of the 2026 fits (unpol2026.py, pol2026.py);
- guards the fixes of audit findings 1-4 (threshold precomputation near E1, table energies, Sigma_2z errors,
  error matrices built from the fitted points).
"""

import os, sys, io, contextlib, unittest, warnings

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)	# the code reads data/... with relative paths
sys.path.insert(0, ROOT)
warnings.filterwarnings('ignore')

with contextlib.redirect_stdout(io.StringIO()):
	import amp	# builds fits I, II, Ahrens and the polarized fit
import integrate
import numpy
from numpy import inf, pi
from compton import E0

FITS = {'I': amp.wt1, 'II': amp.wt2, 'pol': amp.wp}

def sum_rule(w, n):	# Baldin (n=2) and higher moments in 10^-4 fm^(n+1); I_GDH (n=1) in mcb; gamma0, gamma0bar (n=3, 5) in 10^-6 fm^(n+1)
	r = w.integral_n(n, E0, inf)
	return 4*pi**2*r if n == 1 else (-100*r if n in (3, 5) else r)

def amplitudes(k, err=True):	# Re f (fit I), Re f (fit II), Re g  in mcb GeV
	if not err:	# values only (no error propagation, much faster)
		return amp.wt1.integral_f(k, err=False), amp.wt2.integral_f(k, err=False), amp.wp.integral_g_subtracted(k, err=False)
	return amp.wt1.re_f(k), amp.wt2.re_f(k), amp.wp.re_g(k)

# current values (value, error)
SUM_RULES = {
	('I', 2): (14.29336, 0.12725), ('I', 4): (6.082363, 0.090501), ('I', 6): (4.358922, 0.078211),
	('II', 2): (13.84881, 0.1257), ('II', 4): (6.008037, 0.091294), ('II', 6): (4.421464, 0.080107),
	('pol', 1): (204.4647, 21.369), ('pol', 3): (-92.88449, 10.386), ('pol', 5): (48.42621, 8.2287),
}
AMPLITUDES = {
	0.05: ((-2.842096, 0.001689), (-2.847785, 0.00167), (-0.2622255, 0.0003514)),
	0.1: ((-2.205285, 0.007959), (-2.228372, 0.007907), (-0.5357793, 0.003391)),
	0.15: ((-0.3929199, 0.0319), (-0.4250318, 0.0323), (-0.5218651, 0.02305)),
	0.2: ((1.13035, 0.04167), (0.9211955, 0.03987), (-2.352231, 0.02737)),
	0.3: ((2.181357, 0.1404), (1.790078, 0.1532), (-3.003071, 0.1589)),
	0.44: ((-7.09996, 0.1623), (-6.326045, 0.1075), (1.806493, 0.1532)),
	1.0: ((-9.728304, 0.2021), (-9.326246, 0.1999), (1.461064, 0.3546)),
	1.05: ((-10.54067, 0.2846), (-10.13234, 0.2791), (2.802539, 0.3587)),
	1.075: ((-12.03095, 0.4006), (-11.61426, 0.387), (3.349768, 0.3838)),
	2.2: ((-10.343, 0.2041), (-9.977031, 0.2026), (2.097569, 0.5251)),
}
OBSERVABLES = {	# fit I; ds in nb/sr, refg = Re(f g*) in nb/sr, s2z dimensionless
	('ds', 0.1): (13.22705, 0.0906402), ('refg', 0.1): (3.034434, 0.022109), ('s2z', 0.1): (-0.4588224, 0.00297035),
	('ds', 0.15): (1.09597, 0.0892192), ('refg', 0.15): (0.5266031, 0.0486699), ('s2z', 0.15): (-0.9609808, 0.0245548),
	('ds', 0.2): (28.8673, 0.684016), ('refg', 0.2): (-3.60849, 0.343432), ('s2z', 0.2): (0.2500054, 0.0196919),
	('ds', 0.3): (525.7371, 15.6839), ('refg', 0.3): (-199.8044, 7.77634), ('s2z', 0.3): (0.7600924, 0.0235951),
}

# published: GHP15 Table "sumruletest" (fits I, II), GHP16 Table "sumruletest" ("This work") and Table "values_forward_amps"
PUBLISHED_SUM_RULES = {
	('I', 2): (14.29, 0.27), ('I', 4): (6.08, 0.12), ('I', 6): (4.36, 0.09),
	('II', 2): (13.85, 0.22), ('II', 4): (6.01, 0.11), ('II', 6): (4.42, 0.08),
	('pol', 1): (204.5, 21.4), ('pol', 3): (-92.9, 10.5), ('pol', 5): (48.4, 8.2),
}
PUBLISHED_AMPLITUDES = {	# nu [MeV]: Re f (fit I), Re f (fit II), Re g
	50: ((-2.84, 0.01), (-2.85, 0.01), (-0.26, 0.01)), 60: ((-2.76, 0.01), (-2.76, 0.01), (-0.32, 0.01)),
	70: ((-2.65, 0.01), (-2.66, 0.01), (-0.37, 0.01)), 79: ((-2.54, 0.01), (-2.56, 0.01), (-0.42, 0.01)),
	90: ((-2.38, 0.01), (-2.40, 0.01), (-0.48, 0.01)), 100: ((-2.21, 0.01), (-2.23, 0.01), (-0.54, 0.01)),
	110: ((-2.00, 0.01), (-2.03, 0.01), (-0.59, 0.01)), 120: ((-1.75, 0.01), (-1.78, 0.01), (-0.64, 0.01)),
	130: ((-1.45, 0.02), (-1.48, 0.02), (-0.68, 0.01)), 135: ((-1.27, 0.02), (-1.30, 0.02), (-0.69, 0.01)),
	145: ((-0.79, 0.03), (-0.83, 0.03), (-0.65, 0.02)), 150: ((-0.39, 0.03), (-0.42, 0.03), (-0.52, 0.02)),
	160: ((0.04, 0.04), (0.00, 0.04), (-0.68, 0.03)), 180: ((0.48, 0.04), (0.34, 0.04), (-1.53, 0.03)),
	200: ((1.13, 0.04), (0.92, 0.04), (-2.35, 0.03)), 225: ((2.25, 0.05), (1.94, 0.05), (-3.38, 0.03)),
	245: ((3.29, 0.05), (2.87, 0.05), (-4.14, 0.03)), 265: ((4.09, 0.04), (3.55, 0.05), (-4.59, 0.03)),
	275: ((4.14, 0.03), (3.58, 0.04), (-4.53, 0.03)), 300: ((2.15, 0.14), (1.82, 0.16), (-3.00, 0.18)),
	325: ((-2.00, 0.22), (-2.41, 0.24), (-0.44, 0.22)), 440: ((-7.12, 0.16), (-6.27, 0.11), (1.81, 0.15)),
	580: ((-2.93, 0.18), (-2.47, 0.11), (-1.18, 0.19)), 750: ((-7.12, 0.14), (-6.63, 0.13), (1.45, 0.25)),
	1000: ((-9.73, 0.19), (-9.32, 0.19), (1.46, 0.36)), 2200: ((-10.35, 0.21), (-9.97, 0.21), (2.10, 0.55)),
}


class Pinned(unittest.TestCase):	# fails on any change of the numbers: update the pins only for an intended change

	def check(self, r, ref, what):
		self.assertLess(abs(float(r) - ref[0]), 1e-4*max(1., abs(ref[0])), what + ' value')
		self.assertLess(abs(r.err - ref[1]), 1e-3*ref[1], what + ' error')

	def test_sum_rules(self):
		for (f, n), ref in SUM_RULES.items():
			self.check(sum_rule(FITS[f], n), ref, 'fit %s, n=%d' % (f, n))

	def test_amplitudes(self):
		for k, refs in AMPLITUDES.items():
			for r, ref, what in zip(amplitudes(k), refs, ('Re f I', 'Re f II', 'Re g')):
				self.check(r, ref, '%s(%g GeV)' % (what, k))

	def test_observables(self):
		for (name, k), ref in OBSERVABLES.items():
			self.check(getattr(amp, name)(k), ref, '%s(%g GeV)' % (name, k))


class Published(unittest.TestCase):	# within the quoted errors (tables rounded to the last digit)

	def test_sum_rules(self):
		for (f, n), (v, e) in PUBLISHED_SUM_RULES.items():
			self.assertLessEqual(abs(float(sum_rule(FITS[f], n)) - v), e, 'fit %s, n=%d' % (f, n))

	def test_amplitudes(self):
		for k, refs in PUBLISHED_AMPLITUDES.items():
			for r, (v, e), what in zip(amplitudes(k/1000., err=False), refs, ('Re f I', 'Re f II', 'Re g')):
				self.assertLessEqual(abs(float(r) - v), e + 0.005, '%s(%d MeV)' % (what, k))


class AuditFixes(unittest.TestCase):

	def test_threshold_precomputation_near_E1(self):	# finding 1: precomputed (interpolated) = direct integral
		for w in (amp.wt1, amp.wt2):
			dr = w.data_region
			for k, tol in ((0.300, 0.01), (0.305, 0.01), (0.3085, 0.05)):	# at 0.3085 the step of sigma at E1 enters the trapezoid end point
				pre = float(dr.integral_f(k))
				direct = float(integrate.Region.integral_f(dr, k, dr.k_min, dr.k_max))
				self.assertLess(abs(pre - direct), tol, '%s at %g GeV: %g vs %g' % (w.filename, k, pre, direct))

	def test_table_energy(self):	# finding 2: k on the threshold table grid
		k, h = 0.150, 2e-5
		self.assertIn(k, amp.wt1.data_region.x)
		mid = 0.5*(float(amp.ds(k - h)) + float(amp.ds(k + h)))
		self.assertLess(abs(float(amp.ds(k)) - mid), 1e-3*mid)

	def test_s2z_error(self):	# finding 3: |S| |rel(refg) - rel(ds)| <= err(S) <= |S| (rel(refg) + rel(ds))
		for k in (0.1, 0.2, 0.3):
			S, r, d = amp.s2z(k), amp.refg(k), amp.ds(k)
			rr, rd = r.err/abs(float(r)), d.err/float(d)
			self.assertLessEqual(S.err, abs(float(S))*(rr + rd)*(1 + 1e-9))
			self.assertGreaterEqual(S.err, abs(float(S))*abs(rr - rd)*(1 - 1e-9))

	def test_error_matrix_from_fit_points(self):	# finding 4: resonance block of the whole matrix = matrix of the resonance fit
		for w in (amp.wt1, amp.wt2):
			Cr = w.r_region.fit.pcov
			ratio = numpy.sqrt(numpy.diag(w.fit.pcov)[:15]/numpy.diag(Cr)[:15])	# first five resonances (the sixth sits at the 2 GeV junction)
			self.assertLess(abs(ratio - 1).max(), 0.02, w.filename)


class Update2026(unittest.TestCase):	# the 2026 fits as stored in data/2026/*_fit.json (no refit); f: unpol2026, g: pol2026

	SUM_RULES = {
		('f', 2): (14.0092, 0.18143), ('f', 4): (5.949925, 0.090506), ('f', 6): (4.29274, 0.07509),
		('g', 1): (207.7411, 21.406), ('g', 3): (-102.4547, 8.0119), ('g', 5): (44.85594, 5.8208),
	}
	AMPLITUDES = {	# Re f, Re g
		0.1: ((-2.221681, 0.01069), (-0.5384334, 0.002516)),
		0.3: ((2.033378, 0.06668), (-2.977288, 0.03546)),
	}

	@classmethod
	def setUpClass(cls):
		with contextlib.redirect_stdout(io.StringIO()):
			import unpol2026, pol2026
			cls.w = {'f': unpol2026.WholeRegion2026(), 'g': pol2026.WholeRegionPol2026('bw')}

	check = Pinned.check

	def test_sum_rules(self):
		for (a, n), ref in self.SUM_RULES.items():
			self.check(sum_rule(self.w[a], n), ref, '2026 %s, n=%d' % (a, n))

	def test_amplitudes(self):
		for k, (rf, rg) in self.AMPLITUDES.items():
			self.check(self.w['f'].re_f(k), rf, '2026 Re f(%g GeV)' % k)
			self.check(self.w['g'].re_g(k), rg, '2026 Re g(%g GeV)' % k)


if __name__ == '__main__':
	unittest.main()
