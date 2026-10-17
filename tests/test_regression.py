#!/usr/bin/env python3
"""Regression test.  Run from the repository root:  python -m unittest discover -s tests   (about a minute)

- pins the current sum rules, amplitudes and forward observables;
- checks them against the published tables of Gryniuk, Hagelstein, Pascalutsa, PRD 92, 074031 (2015) and
  PRD 94, 034043 (2016), within the quoted errors;
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

def amplitudes(k):	# Re f (fit I), Re f (fit II), Re g  in mcb GeV
	return amp.wt1.re_f(k), amp.wt2.re_f(k), amp.wp.re_g(k)

# current values (value, error)
SUM_RULES = {
	('I', 2): (14.29278, 0.12726), ('I', 4): (6.082344, 0.090502), ('I', 6): (4.358921, 0.078211),
	('II', 2): (13.84835, 0.12571), ('II', 4): (6.008021, 0.091294), ('II', 6): (4.421463, 0.080107),
	('pol', 1): (203.0205, 21.605), ('pol', 3): (-91.44974, 10.742), ('pol', 5): (48.92885, 8.2849),
}
AMPLITUDES = {
	0.05: ((-2.842104, 0.001689), (-2.847791, 0.00167), (-0.2621784, 0.0003628)),
	0.1: ((-2.205315, 0.007976), (-2.228396, 0.007908), (-0.5353744, 0.003478)),
	0.15: ((-0.3929873, 0.03192), (-0.4250859, 0.0323), (-0.5203053, 0.02324)),
	0.2: ((1.130228, 0.04172), (0.9210978, 0.03987), (-2.347618, 0.02876)),
	0.3: ((2.181071, 0.1404), (1.789849, 0.1532), (-2.945351, 0.1952)),
	0.44: ((-7.100639, 0.1626), (-6.32659, 0.1076), (1.798395, 0.1554)),
	1.0: ((-9.744597, 0.2032), (-9.337921, 0.2018), (1.399858, 0.3793)),
	2.2: ((-10.33842, 0.2055), (-9.973302, 0.2039), (1.998537, 0.5714)),
}
OBSERVABLES = {	# fit I; ds in nb/sr, refg = Re(f g*) in nb/sr, s2z dimensionless
	('ds', 0.1): (13.22628, 0.0906654), ('refg', 0.1): (3.032182, 0.0225357), ('s2z', 0.1): (-0.4585088, 0.0030284),
	('ds', 0.15): (1.091931, 0.0894553), ('refg', 0.15): (0.5251191, 0.0486527), ('s2z', 0.15): (-0.9618172, 0.0243762),
	('ds', 0.2): (28.81091, 0.692039), ('refg', 0.2): (-3.594366, 0.344108), ('s2z', 0.2): (0.2495142, 0.0197477),
	('ds', 0.3): (524.8521, 15.7916), ('refg', 0.3): (-199.4789, 7.80852), ('s2z', 0.3): (0.7601336, 0.0235409),
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
			for r, (v, e), what in zip(amplitudes(k/1000.), refs, ('Re f I', 'Re f II', 'Re g')):
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


if __name__ == '__main__':
	unittest.main()
