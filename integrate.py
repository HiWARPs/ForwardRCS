#!/usr/bin/env python

import os
# from scipy import pi, sqrt, inf, exp, nan, linalg
# from scipy import sqrt, nan, linalg
# from scipy.integrate import simps, trapz
from scipy.integrate import trapezoid as trapz	# py3 port: scipy.integrate.trapz/simps were removed; trapezoid is the same routine (simps is overridden below anyway)
# quad = romberg
simps = trapz

# import mpmath
# from mpmath import pi, inf, exp, log
# quad = lambda f, k0, k1, **kw: mpmath.quad(f, [k0, k1], **dict(kw, error=True))

# from scipy.misc import derivative
# import numpy

from val import *
# from compton import *
# from common import *

from fit import *
# import fit
# for a in ['log', 'exp', 'pi']:
#	setattr(fit, a, getattr(mpmath, a))

quad = lambda f, k0, k1, **kw: mpmath.quad(lambda k: f(float(k)), [k0, k1], **dict(kw, error=True, method='tanh-sinh'))


caching = True


def aargs(self, *args, **kwargs):
	try:
		p = self.get_p()
	except AttributeError:
		p = ()
	return args + (str(kwargs),) + p
	# return args + tuple(kwargs.iteritems())

def cached_args_method(size=10000000):
	def decorator(f):
		if not caching:
			return f
		else:
			name = f.__name__
			def cached(self, *args, **kwargs):	# added kwargs - not optimal solution so far (read - "shot in the foot")
				# args = tuple(list(args) + [kwargs[k] for k in sorted(kwargs)])	# hack hack hack
				cargs = aargs(self, *args, **kwargs)
				if not hasattr(self, '_cache'):
					self._cache = {}
				if name not in self._cache:
					self._cache[name] = {}
				cache = self._cache[name]
				if cargs in cache:
					return cache[cargs]
				else:
					if len(cache) == size:
						del cache[next(iter(cache), None)]	# py3 port
					r = f(self, *args, **kwargs)
					cache[cargs] = r
					return r
			return cached
	return decorator


class _oldstyle_property(object):	# py3 port: Python-2 classes here were old-style, where a property only provides a getter and an instance
	def __init__(self, fget):		# assignment (e.g. `self.baldin = None` in tot.py/pol.py) just shadows it; a non-data descriptor reproduces that
		self.fget = fget
	def __get__(self, obj, objtype=None):
		return self if obj is None else self.fget(obj)


class Data:
	def __init__(self, path, filename='data1.dat'):
		self.x, self.y, self.err = map(list, read_xye(os.path.join('data', path, filename)))


class Region:
	
	def update_data(self, filename=None, k_min=None, k_max=None):
		self._cache = {}
		y, e, n = {}, {}, {}
		for ax, ay, ae in zip(*map(list, read_xye(filename or os.path.join('data', self.path, self.filename)))):
			if ax in n:
				y[ax] = (y[ax]*n[ax] + ay)/(n[ax] + 1)
				e[ax] = (e[ax]*n[ax] + ae)/(n[ax] + 1)
				n[ax] += 1
			else:
				y[ax] = ay
				e[ax] = ae
				n[ax] = 1
		self.x, self.y, self.err = map(list, zip(*([ax, y[ax], e[ax]] for ax in sorted(n.keys()))))
		self.x, self.y, self.err = self.cut_xye(k_min or self.k_min, k_max or self.k_max)	# hack
		return self
	
	def update_fit(self, quiet=False):
		self._cache = {}
		self.fit.fit(self.x, self.y, self.err, quiet=quiet)
		return self
	
	def __init__(self, fit=None, k_min=0., k_max=inf, name='', init=True, p=None, path='tot', filename='mdata.dat', quiet=False):
		self.bw = Smooth(lambda k: 0., ResonanceFit(flab2cm(0.3), 0.1, 0.01), k_min, 0.01)
		self.filename = filename
		self.path = path
		# make_sure_path_exists(os.path.join('out', path))
		if fit is not None:
			self.fit = fit
		if p is not None:
			self.fit.set_p(*self.fit.full_p(p))
		self.k_min, self.k_max = k_min, k_max
		self.name = name
		if init:
			self.update_data()
			self.update_fit(quiet)
	
	@cached_args_method()
	def cs(self, k, err=True):
		if not (self.k_min <= k <= self.k_max):
			if err: return val(0)
			return 0.
		if err:
			return val(self.fit(k), self.cs_err(k))
		else:
			return self.fit(k)
	
	@cached_args_method()
	def dcs(self, i, k):
		return self.fit.df(i, k)
	
	def cut_xye(self, k0, k1):	# helper method (TODO: rewrite this mess some day)
		x0, x, x1 = split_sorted(self.x, k0, k1)
		l0, l1 = map(len, (x0, x1))
		y, e = self.y[l0:(-l1 or None)], self.err[l0:(-l1 or None)]
		if x == []:
			if l0>0 and l1>0:
				liny = linear(x0[-1], x1[0], self.y[l0-1], self.y[l0])
				line = linear(x0[-1], x1[0], self.err[l0-1], self.err[l0])
			elif l1>0:
				liny = linear(x1[0], 2*x1[0], self.y[l0], self.y[l0])
				line = linear(x1[0], 2*x1[0], self.err[l0], self.err[l0])
			elif l0>0:
				liny = linear(x0[-1], 2*x0[-1], self.y[l0-1], self.y[l0-1])
				line = linear(x0[-1], 2*x0[-1], self.err[l0-1], self.err[l0-1])
			else:
				return [], [], []
			x = [k0, k1]
			y = list(map(liny, x))	# py3 port
			e = list(map(line, x))	# py3 port
		else:
			if l0>0 and x[0]!=k0:
				y = [linear(x0[-1], x[0], self.y[l0-1], y[0])(k0)] + y
				e = [linear(x0[-1], x[0], self.err[l0-1], e[0])(k0)] + e
				x = [k0] + x
			if l1>0:
				y += [linear(x[-1], x1[0], y[-1], self.y[-l1])(k1)]
				e += [linear(x[-1], x1[0], e[-1], self.err[-l1])(k1)]
				x += [k1]
			if l0 == 0 and l1 == 0:
				pass	# TODO
		return x, y, e
	
	def _stat_error(self, f, k0, k1):
		x, y, e = self.cut_xye(k0, k1)
		v = [f(ax) * (self.cs(ax, err=False) - ay)**2 / ae for ax, ay, ae in zip(x, y, e)]
		return sum([(x[i+1] - x[i]) * (v[i+1] + v[i]) / 2 for i in range(len(x) - 1) if x[i] != x[i+1]])
	
	@staticmethod
	def singular(k, k0, k1):
		f = lambda p: 0. if p is inf else log(abs((p-k)/(p+k)))/(2*k)
		return f(k1) - f(k0)
	
	@cached_args_method()
	def cst(self, k):	# this is the one from bare QED (without kappa in the vertex) + artificial Breit-Weigner
		x = k/float(Mp)
		ln = log(1+2*x)
		r = (1+x)/x**3*(2*x*(1+x)/(1+2*x)-ln) + ln/(2*x) - (1+3*x)/(1+2*x)**2
		return 2*pi*float(alpha)**2/float(Mp)**2 * r * 10**4*float(hc)**2 + self.bw(k)
		return 0.
	
	@cached_args_method()
	def cst1(self, k):
		x = k/float(Mp)
		ka = float(kappa)
		r = ka**4*(2*x**6 + 4*x**5 + 3*x**4 + x**3) + ka**3*(16*x**5 + 12*x**4 + 4*x**3) + ka**2*(16*x**5 - 6*x**4 - 2*x**3) + ka*(-24*x**4 - 8*x**3) + 4*x**4 + 36*x**3 + 32*x**2 + 8*x + (-2*ka**4*x**4 - 2*ka**4*x**3 - ka**4*x**2/2 - 8*ka**3*x**4 - 8*ka**3*x**3 - 2*ka**3*x**2 + 4*ka**2*x**4 + 4*ka**2*x**3 + ka**2*x**2 + 16*ka*x**4 + 16*ka*x**3 + 4*ka*x**2 + 8*x**4 - 8*x**3 - 30*x**2 - 20*x - 4)*log(2*x + 1)
		return pi*float(alpha)**2/(2 * float(Mp)**2 * x**3*(2*x+1)**2) * r * 10**4*float(hc)**2
	
	@cached_args_method()
	def csp(self, k):	# this is the one straight from QED (without kappa in the vertex)
		x = k/float(Mp)
		ln = log(1+2*x)
		r = (1+1/x)*ln - 2*(1 + x**2/(1+2*x)**2)
		return pi*float(alpha)**2/float(Mp)**2 * r * 10**4*float(hc)**2 / x
		return 0.
	
	def get_avr_cs_err(self):
		r, n = 0, 0
		for y, e in zip(self.y, self.err):
			if y != 0 and e/y > 0.001:
				r += e/y
				n += 1
		return r/n
	
	# @cached_args_method()
	# def dcs_dp(self, k, i, fit=None):
	# 	if fit is None:
	# 		try:
	# 			fit = self.fit.f1
	# 		except AttributeError:
	# 			fit = self.r_region.fit
	# 	p0 = list(fit.get_p())
	# 	def ds_p(x):
	# 		p = p0[:]
	# 		p[i] = x
	# 		fit.set_p(*p)
	# 		r = self.fit(k)
	# 		# r = fit(k)
	# 		fit.set_p(*p0)
	# 		return r
	# 	return derivative(ds_p, p0[i], dx=10**-6)
	
	@cached_args_method()
	def cs_err(self, k):
		return self.fit.err(k)
	
	# def get_pcov(self):
	# 	r = self.r_region
	# 	fit = r.fit
	# 	m = len(fit.get_p())
	# 	x, err = zip(*[[ax, ae] for ax, ae in zip(r.x, r.err) if ax >= r.k_min and ax <= r.k_max])
	# 	J = numpy.array([[self.dcs_dp(k, i) for i in range(m)] for k in x])
	# 	return numpy.linalg.inv(J.T.dot(numpy.diag(numpy.array(err)**-2)).dot(J))
	
	#******************************************************************************************************************
	def _gen_fit_ref_method(method):
		def fit_ref(self, *args, **kwargs):
			return getattr(self.fit, method)(*args, **kwargs)
		return fit_ref
	
	for method in ['get_p', 'set_p', 'get_pcov']:
		vars()[method] = _gen_fit_ref_method(method)
	#******************************************************************************************************************
	
	def _int(self, f, m, k, k0, k1, subtracted=True, singular=None, factor=1):
		singular = singular or self.singular
		# def integrate(f):
		# 	i = quad(f, k0, k1, full_output=1)
		# 	# print i
		# 	return val(i[:2])
		# if k > fcm2lab(5e1):
		# 	cuts = [.5, 1, 2, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8]
		# 	# def integrate(f):
		# 	# 	s = map(lambda (a, b): val(quad(f, a, b)), zip([k0] + cuts, cuts + [k1]))
		# 	# 	for v in s:
		# 	# 		print ' ', v
		# 	# 	return sum(s)
		# 	# integrate = lambda f: val(quad(lambda x: exp(x)*f(exp(x)), log(k0), log(k1)))
		# 	# base = 3*exp(1)
		# 	# cut = 1e1
		# 	# integrate = lambda f: val(quad(lambda x: base**x * f(base**x), log(k0)/log(base), log(k1)/log(base)))
		# 	# integrate = lambda f: val(quad(f, k0, cut)) + val(quad(lambda x: base**x * f(base**x), log(cut)/log(base), log(k1)/log(base)))
		# 	integrate = lambda f: sum(map(lambda (a, b): val(quad(f, a, b)), zip([k0] + cuts, cuts + [k1])))
		# 	# integrate = lambda f: sum(map(lambda (a, b): val(quad(f, a, b, points=[k] if a < k < b else [])), zip([k0] + cuts, cuts + [1e10])))
		# 	# integrate = lambda f: val(quad(f, k0, 1e9, points=cuts))
		# else:
		# 	integrate = lambda f: val(quad(f, k0, k1))
		if subtracted:
			return val(quad(lambda p: factor * (p**m*f(p) - k**m*f(k))/(p**2 - k**2), k0, k1)) + factor * k**m*f(k)*singular(k, k0, k1)
			# cut = 2*k - k0
			# left = quad(lambda p: factor * p**m*f(p)/(p + k), k0, cut, weight='cauchy', wvar=k)	# THIS IS WRONG BECAUSE OF HACKS FOR THE WHOLEREGION BELOW!
			# right = quad(lambda p: factor * p**m*f(p)/(p**2 - k**2), cut, k1)
			# return val(left) + val(right)
		return val(quad(lambda p: factor * p**m*f(p)/(p**2 - k**2), k0, k1))
		# def integrate(k0, k1):
		# 	return lambda f: val(quad(f, k0, k1))
		# if subtracted:
		# 	f_k = f(k)
		# 	# if k > fcm2lab(1e2):
		# 	# 	d = .2*(k - k0)
		# 	# 	c0, c1 = k - d, k + d
		# 	# 	sub = integrate(c0, c1)(lambda p: factor * (p**m*f(p) - k**m*f_k)/(p**2 - k**2)) + factor * k**m*f_k*singular(k, c0, c1)
		# 	# 	left = integrate(k0, c0)(lambda p: factor * p**m*f(p)/(p**2 - k**2))
		# 	# 	right = integrate(c1, k1)(lambda p: factor * p**m*f(p)/(p**2 - k**2))
		# 	# 	return left + sub + right
		# 	return integrate(k0, k1)(lambda p: factor * (p**m*f(p) - k**m*f_k)/(p**2 - k**2)) + factor * k**m*f_k*singular(k, k0, k1)
		# else:
		# 	return integrate(k0, k1)(lambda p: factor * p**m*f(p)/(p**2 - k**2))
	
	def dfi_err(self, dfi):
		adf = numpy.array([dfi(i) for i in range(len(self.get_p()))])
		return adf.dot(self.get_pcov()).dot(adf)**0.5
	
	def _int_err(self, df, m, k, k0, k1, subtracted=True, factor=1):
		return self.dfi_err(lambda i: self._int(lambda p: df(i,p), m, k, k0, k1, subtracted, singular=Region.singular, factor=factor).val)	# Region.singular !!!
	
	def integral(self, n, k, k0=None, k1=None, err=True, cs=None, dcs=None, subtracted=True, const=False, factor=1, **err_kwargs):
		k0, k1 = (self.k_min if k0 is None else k0), (self.k_max if k1 is None else k1)
		k0 = min(max(k0, self.k_min), self.k_max)
		k1 = max(min(k1, self.k_max), self.k_min)
		r = self._int(cs or (lambda a: float(self.cs(a, err=False))), 2-n, k, k0, k1, subtracted, factor=factor)
		if err:
			r += val(0, self._int_err(dcs or self.dcs, 2-n, k, k0, k1, subtracted=subtracted, **err_kwargs))
			# r.err = (r.err**2 + self._int_err(dcs or self.dcs, 2-n, k, k0, k1, subtracted)**2)**0.5
		return r/(2*pi**2) if const else k**n*r/(2*pi**2)
	
	@cached_args_method()
	def integral_f(self, k, k0=None, k1=None, err=True, subtracted=True):
		# return e_Mp + self.integral(2, k, k0, k1, err, cs=(lambda a: float(self.cs(a, err=False) + self.cst(a))), subtracted=subtracted)
		return e_Mp + self.integral(2, k, k0, k1, err, cs=(lambda a: float(self.cs(a, err=False))), subtracted=subtracted)
	
	@cached_args_method()
	def dintegral_f(self, i, k, k0=None, k1=None, subtracted=True):
		return self.integral(2, k, k0, k1, err=False, cs=(lambda a: self.dcs(i, a)), subtracted=subtracted)
	
	@cached_args_method()
	def integral_g(self, k, k0=None, k1=None, err=True, subtracted=True):
		# return - self.integral(1, k, k0, k1, err, cs=(lambda a: float(self.cs(a, err=False) + self.csp(a))), subtracted=subtracted)
		return - self.integral(1, k, k0, k1, err, cs=(lambda a: float(self.cs(a, err=False))), subtracted=subtracted)
	
	@cached_args_method()
	def dintegral_g(self, i, k, k0=None, k1=None, subtracted=True):
		return - self.integral(1, k, k0, k1, err=False, cs=(lambda a: self.dcs(i, a)), subtracted=subtracted)
	
	@cached_args_method()
	def integral_g_subtracted(self, k, k0=None, k1=None, err=True, subtracted=True):
		# return - self.integral(3, k, k0, k1, err, cs=(lambda a: float(self.cs(a, err=False) + self.csp(a))), subtracted=subtracted) - k*IGDH/(4*pi**2)
		return - self.integral(3, k, k0, k1, err, cs=(lambda a: float(self.cs(a, err=False))), subtracted=subtracted) - k*IGDH/(4*pi**2)
	
	@cached_args_method()
	def dintegral_g_subtracted(self, i, k, k0=None, k1=None, subtracted=True):
		return - self.integral(3, k, k0, k1, err=False, cs=(lambda a: self.dcs(i, a)), subtracted=subtracted)
	
	@cached_args_method()
	def integral_n(self, n, k0, k1, err=True):
		k0 = min(max(k0, self.k_min), self.k_max)
		k1 = max(min(k1, self.k_max), self.k_min)
		if k0 >= k1:
			return val()
		return hc**(n-1) * self.integral(n, 0, k0, k1, err=err, subtracted=False, const=True)
		# r = val(quad(lambda p: float(self.cs(p))/p**n, k0, k1))# <<pm>> self._stat_error(lambda p: 1/p**n, k0, k1)
		# if err:
		# 	dcs = numpy.array([quad(lambda p: float(self.dcs(i, p))/p**n, k0, k1)[0] for i in range(len(self.fit.get_p()))])
		# 	r.err = (r.err**2 + dcs.dot(self.fit.get_pcov()).dot(dcs))**0.5
		# 	# r.err = quad(lambda p: self.cs_err(p)/p**n, k0, k1)[0]	# wrong!
		# return hc**(n-1) * r / (2*pi**2)
	
	@cached_args_method()
	def integral_n_log(self, n, k0, k1, err=True):
		k0 = min(max(k0, self.k_min), self.k_max)
		k1 = max(min(k1, self.k_max), self.k_min)
		if k0 >= k1:
			return val()
		return hc**(n-1) * self.integral(n, 0, k0, k1, err=err, subtracted=False, const=True, cs=(lambda a: float(self.cs(a, err=False))*log(a)))
	
	def integral_c(self, k1):
		return 4*Mp*pi**2*10**-4/hc * self.integral_n(0, self.k_min, k1)
	
	def re_f(self, k, err=True):
		return self.integral_f(k, err=err)
	
	def im_f(self, k, err=True):
		return self.cs(k, err=err) * k/(4*pi)
	
	@cached_args_method()
	def re_over_im_f(self, k, err=True, subtracted=True):
		im = self.cs(k, err=False) * k/(4*pi)
		# factor = 1./im
		factor = 1.
		re_im = val(0, e_Mp.err*factor) + float(e_Mp)*factor + self.integral(2, k, k0=None, k1=None, err=False, cs=(lambda a: float(self.cs(a, err=False))), subtracted=subtracted, factor=factor)
		re_im /= im
		if err:
			return val(re_im) |pm| self.dfi_err(lambda i: self.dre_over_im_f(i, k, subtracted=subtracted))
		return re_im
	
	@cached_args_method()
	def dre_over_im_f(self, i, k, subtracted=True):
		re_im = float(self.re_over_im_f(k, err=False, subtracted=subtracted))
		dre = self.dintegral_f(i, k, subtracted=subtracted)
		im = self.im_f(k, err=False)
		dim = self.dim_f(i, k)
		return (dre - dim*re_im)/im
	
	def f(self, k):
		return float(self.re_f(k)) + 1j * float(self.im_f(k))
	
	def dre_f(self, i, k):
		return self.dintegral_f(i, k)
	
	def dim_f(self, i, k):
		return self.dcs(i, k) * k/(4*pi)
	
	def df(self, i, k):
		return float(self.dre_f(i, k)) + 1j * float(self.dim_f(i, k))
	
	def re_g0(self, k):
		return self.integral_g(k)	# unsubtracted
	
	def dre_g0(self, i, k):
		return self.dintegral_g(i, k)	# unsubtracted
	
	def re_g(self, k):
		return self.integral_g_subtracted(k)	# using subtracted integral here
	
	def dre_g(self, i, k):
		return self.dintegral_g_subtracted(i, k)	# using subtracted integral here
	
	def im_g(self, k):
		return - self.cs(k) * k/(4*pi)	# flipping cs for polarized amplitude here
	
	def dim_g(self, i, k):
		return - self.dcs(i, k) * k/(4*pi)	# flipping cs for polarized amplitude here
	
	def g(self, k):
		return float(self.re_g(k)) + 1j * float(self.im_g(k))
	
	def dg(self, i, k):
		return float(self.dre_g(i, k)) + 1j * float(self.dim_g(i, k))
	
	def sr_1(self, k0=None, k1=None):
		return 4*pi**2 * self.integral_n(1, self.k_min if k0 is None else k0, self.k_max if k1 is None else k1)
	
	def sr_2(self, k0=None, k1=None):
		return self.integral_n(2, self.k_min if k0 is None else k0, self.k_max if k1 is None else k1)
	
	def sr_2_log(self, k0=None, k1=None):
		return self.integral_n_log(2, self.k_min if k0 is None else k0, self.k_max if k1 is None else k1)
	
	@property
	def mean_exct(self):
		return exp(float(self.sr_2_log()/self.sr_2()))
	
	def sr_3(self, k0=None, k1=None):
		return -100 * self.integral_n(3, self.k_min if k0 is None else k0, self.k_max if k1 is None else k1)
	
	def sr_4(self, k0=None, k1=None):
		return self.integral_n(4, self.k_min if k0 is None else k0, self.k_max if k1 is None else k1)
	
	def sr_5(self, k0=None, k1=None):
		return -100 * self.integral_n(5, self.k_min if k0 is None else k0, self.k_max if k1 is None else k1)
	
	def sr_6(self, k0=None, k1=None):
		return self.integral_n(6, self.k_min if k0 is None else k0, self.k_max if k1 is None else k1)
	
	gdh = _oldstyle_property(sr_1)	# py3 port: see _oldstyle_property
	baldin = _oldstyle_property(sr_2)
	
	def plot(self, f, prefix='plot', xmin=None, xmax=None, dx=None, logx=False, filename=None, raw=True, **kwargs):
		if filename is None:
			make_sure_path_exists(os.path.join('out', self.path))
			filename = os.path.join('out', self.path, prefix+(self.name and '_'+self.name)+('_logx' if logx else '')+'.out')
		xmin = self.k_min if xmin is None else xmin											#
		xmax = min(self.k_max, 1E8) if xmax is None else xmax								#  please don't "prettify" this anymore >_<
		xf = self.cut_xye(xmin, xmax)[0] if dx is None else arange(xmin, xmax, dx, logx)	#
		fx = map_verbose(f, xf, 'plotting to "'+filename+'"')
		data = [xf] + list(zip(*fx)) if is_iterable(fx[0]) else [xf, fx]
		data = list(zip(*[d for d in zip(*data) if nan not in d]))	# dirty hack	(py3 port: list)
		kwargs['plots'] = kwargs.get('plots') or ['1 '+str(i)+' label="'+self.name+'"' for i in range(2, len(data)+1)]
		plot_to_file(data, filename=filename, **({} if raw else kwargs))
		return self
	
	def plot_cs(self, xmin=0., xmax=4., dx=0.001, logx=False, **kwargs):
		kwargs.update({'title': 'Total cross section', 'xylabel': r'"$\nu$ [GeV]" "$\sigma$ [$\mu\mathrm{b}$]"'})
		return self.plot(self.cs, 'cs', xmin, xmax, dx, logx, **kwargs)
	
	def plot_f(self, xmin=0.002, xmax=4., dx=0.002, logx=False, **kwargs):
		kwargs.update({'title': 'Amplitude', 'xylabel': r'"$\nu$ [GeV]" "Re $f-f_0$ [$\mu\mathrm{b}\cdot\mathrm{GeV}$]"'})
		return self.plot(self.integral_f, 'f', xmin, xmax, dx, logx, **kwargs)
	
	def plot_g(self, xmin=0.002, xmax=4., dx=0.002, logx=False, **kwargs):
		kwargs.update({'title': 'Amplitude', 'xylabel': r'"$\nu$ [GeV]" "Re $g$ [$\mu\mathrm{b}\cdot\mathrm{GeV}$]"'})
		return self.plot(self.re_g, 'g', xmin, xmax, dx, logx, **kwargs)
	
	def plot_g0(self, xmin=0.002, xmax=4., dx=0.002, logx=False, **kwargs):
		kwargs.update({'title': 'Amplitude', 'xylabel': r'"$\nu$ [GeV]" "Re $g$ [$\mu\mathrm{b}\cdot\mathrm{GeV}$]"'})
		return self.plot(self.re_g0, 'g0', xmin, xmax, dx, logx, **kwargs)
	
	def plot_im_g(self, xmin=0.002, xmax=4., dx=0.002, logx=False, **kwargs):
		kwargs.update({'title': 'Amplitude', 'xylabel': r'"$\nu$ [GeV]" "Re $g$ [$\mu\mathrm{b}\cdot\mathrm{GeV}$]"'})
		return self.plot(self.im_g, 'im_g', xmin, xmax, dx, logx, **kwargs)
	
	# def plot_g_sub(self, xmin=0.002, xmax=4., dx=0.002, logx=False, **kwargs):
	# 	kwargs.update({'title': 'Amplitude', 'xylabel': r'"$\nu$ [GeV]" "Re $g$ [$\mu\mathrm{b}\cdot\mathrm{GeV}$]"'})
	# 	return self.plot(self.integral_g_subtracted, 'g_sub', xmin, xmax, dx, logx, **kwargs)
	
	def plot_n(self, n, xmin=None, xmax=10., dx=0.04, logx=True, **kwargs):
		kwargs.update(kw = {'title': 'SR', 'xylabel': '"cut-off [GeV]" "SR integral"'})
		return self.plot(lambda x: self.integral_n(n, 0., x), str(n), xmin, xmax, dx, logx, **kwargs)
	
	def plot_1(self, *args, **kwargs):
		kwargs.update(kw = {'title': 'GDH SR'})
		return self.plot_n(1, *args, **kwargs)
	
	def plot_2(self, *args, **kwargs):
		kwargs.update(kw = {'title': 'Baldin SR'})
		return self.plot_n(2, *args, **kwargs)
	
	def plot_3(self, *args, **kwargs):
		kwargs.update(kw = {'title': '$\gamma_0$'})
		return self.plot_n(3, *args, **kwargs)
	
	def plot_4(self, *args, **kwargs):
		kwargs.update(kw = {'title': '4th order SR'})
		return self.plot_n(4, *args, **kwargs)
	
	def plot_5(self, *args, **kwargs):
		kwargs.update(kw = {'title': '$\bar{\gamma_0}$'})
		return self.plot_n(5, *args, **kwargs)
	
	def plot_6(self, *args, **kwargs):
		kwargs.update(kw = {'title': '6th order SR'})
		return self.plot_n(6, *args, **kwargs)
	
	def plot_c(self, xmin=0.002, xmax=5., dx=0.02, logx=True, **kwargs):
		kwargs.update(kw = {'xlabel': r'"cut-off $\nu$ [GeV]"'})
		return self.plot(self.integral_c, 'c', xmin, xmax, dx, logx, **kwargs)



class ThresholdRegion(Region):
	
	def update_data(self, filename=None, k_min=None, k_max=None):
		self._cache = {}
		y, n = {}, {}
		for ax, ay in zip(*map(list, read_xy(filename or os.path.join('data', self.path, self.filename)))):
			if ax in n:
				y[ax] = (y[ax]*n[ax] + ay)/(n[ax] + 1)
				n[ax] += 1
			else:
				y[ax] = ay
				n[ax] = 1
		self.x, self.y = map(list, zip(*([ax, y[ax]] for ax in sorted(n.keys()))))
		self.err = [0.02*abs(ay) for ay in self.y]
		self.k_max = min(self.x[-1], self.k_max)	# another hack
		self.x, self.y, self.err = self.cut_xye(k_min or self.k_min, k_max or self.k_max)	# hack
		return self
	
	def update_fit(self):
		return self
	
	def _gen_method_precompute_x(method):
		def precompute_x(self, quiet=False):
			if not quiet: print('precomputing %s ...'%method)
			setattr(self, '_'+method, list(map(getattr(self, 'integral_'+method), self._x)))	# py3 port: list
		precompute_x.__name__ = 'precompute_'+method
		return precompute_x
	
	def __init__(self, k_min=E0, k_max=E1, name='threshold', path='tot', filename='maid.dat', pre_f=True, pre_g=True, quiet=False):
		Region.__init__(self, k_min=k_min, k_max=k_max, name=name, init=False, path=path, filename=filename)
		self.update_data()
		#self._interp = interp1d(self.x, self.y, fill_value=0., bounds_error=False, assume_sorted=True, kind='linear')
		self._f = self._g = self._g_subtracted = None
		self._x = [1.5*self.x[0] - 0.5*self.x[1]] + [0.5*(self.x[i] + self.x[i+1]) for i in range(len(self.x) - 1)]	# no node above E1: there sigma(k) comes from the fit, attached only after the precomputation (WholeRegion.__init__); energies above the last midpoint are integrated directly
		if pre_f:
			self.precompute_f(quiet=quiet)
		if pre_g:
			self.precompute_g(quiet=quiet)
			self.precompute_g_subtracted(quiet=quiet)
	
	def cs_err(self, k):
		return self.fit(k) * 0.02
	
	def fit(self, k):
		if k < self.k_min:
			return 0.
		else:
			return interp(k, self.x, self.y)
			#return self._interp(k)
	
	def dcs(self, i, k):
		return 0.
	
	def _int(self, f, m, k, k0, k1, subtracted=True, singular=None, factor=1):
		singular = singular or self.singular
		x = self.cut_xye(k0, k1)[0]
		if subtracted:
			f_k = f(k)
			y = [factor * (p**m*f(p) - k**m*f_k)/(p**2 - k**2) if p != k else nan for p in x]	# 0/0 when k is a table energy ...
			for i in [i for i, a in enumerate(y) if a != a]:
				y[i] = numpy.nanmean(y[max(i-1, 0):i+2])	# ... finite limit: mean of the neighbouring nodes
			return val(simps(y, x)) + factor * k**m*f_k*singular(k, k0, k1)
		else:
			return val(simps([factor * p**m*f(p)/(p**2 - k**2) for p in x], x))
	
	def _int_err(self, df, m, k, k0, k1, subtracted=True):
		return abs(self._int(lambda p: self.fit(p)*0.02, m, k, k0, k1, subtracted, singular=Region.singular).val)
	
	def _gen_method_integral_x(method):
		def integral_x(self, k, k0=None, k1=None, err=True, subtracted=True):
			if k0 is None and k1 is None and getattr(self, '_'+method) and self._x[0] <= k <= self._x[-1]:
				return val(*map(lambda fu: interp(k, self._x, fu), zip(*getattr(self, '_'+method))))
			return getattr(Region, 'integral_'+method)(self, k, k0, k1, err, subtracted=subtracted)
		integral_x.__name__ = 'integral_'+method
		return integral_x
	
	def integral_n(self, n, k0, k1, err=True):
		k0 = min(max(k0, self.k_min), self.k_max)
		k1 = max(min(k1, self.k_max), self.k_min)
		if k0 >= k1:
			return val()
		# return hc**(n-1) * self.integral(n, 0, k0, k1, err=err, subtracted=False, const=True, err_abs=True)	# doing "abs" error
		return hc**(n-1) * self.integral(n, 0, k0, k1, err=err, subtracted=False, const=True)
	
	for m in ['f', 'g', 'g_subtracted']:
		vars()['precompute_'+m] = _gen_method_precompute_x(m)
		vars()['integral_'+m] = cached_args_method()(_gen_method_integral_x(m))
		vars()['dintegral_'+m] = lambda *args, **kwargs: 0. 	# NOTE this should be rewritten, so that subtracted integral DOES NOT give zero if singular is zero (see WholeRegion.dintegral_x)



class ResonanceRegion(Region):
	
	def __init__(self, k_min=E1, k_max=E2, name='resonances', init=True, p=None, path='tot', filename='mdata.dat', quiet=False, nres=6, bgd=True):
		Region.__init__(self, ResonanceRegionFit(nres=nres, bgd=bgd), k_min, k_max, name, init, p, path, filename, quiet)



class AsymptoticRegion(Region):
	
	def __init__(self, k_min=E2, k_max=inf, name='regge', init=True, p=None, path='tot', filename='mdata.dat', quiet=False):
		Region.__init__(self, AsymptoticRegionFit(), k_min, k_max, name, init, p, path, filename, quiet)



class WholeRegion(Region):
	
	def update_data(self, filename=None, k_min=None, k_max=None):
		self._cache = {}
		self.data_region.update_data(filename, k_min, k_min and min(k_min, self.data_region.k_min))
		self.r_region.update_data(filename, k_min and max(k_min, self.r_region.k_min), k_max and min(k_max, self.r_region.k_max))
		self.h_region.update_data(filename, k_min and max(k_min, self.h_region.k_min), k_max)
		self.fitted_region.update_data(filename, k_min and max(k_min, self.fitted_region.k_min), k_max)
		#xye = sorted(set(zip(*(lambda a, b, c: (a.x + b.x + c.x, a.y + b.y + c.y, a.err + b.err + c.err))(self.data_region, self.r_region, self.h_region))), key=lambda a: a[0])
		#self.x, self.y, self.err = zip(*xye)
		return Region.update_data(self, filename, k_min, k_max)
	
	def update_fit(self, quiet=False):
		self._cache = {}
		self.r_region.update_fit(quiet=quiet)
		self.h_region.update_fit(quiet=quiet)
		return self
	
	def __init__(self, b=[E0, E1, E2, inf], name='', path='tot', filename='mdata.dat', init=True, quiet=False, nres=6, bgd=True, pre_f=True, pre_g=True, ThresholdRegion=ThresholdRegion, ResonanceRegion=ResonanceRegion, AsymptoticRegion=AsymptoticRegion):
		self.b = b
		self.data_region = ThresholdRegion(b[0], b[1], path=path, filename=filename, pre_f=False, pre_g=False)
		self.r_region = ResonanceRegion(k_min=b[1], k_max=b[2], path=path, filename=filename, init=False, nres=nres, bgd=bgd)
		self.h_region = AsymptoticRegion(k_min=b[2], k_max=b[-1], path=path, filename=filename, init=False)
		self.data_region.singular = lambda k, k0, k1: 0.	# hack (improving precision)
		if pre_f: self.data_region.precompute_f(quiet=quiet)
		if pre_g:
			self.data_region.precompute_g(quiet=quiet)
			self.data_region.precompute_g_subtracted(quiet=quiet)
		self.smooth = Smooth(self.r_region.fit, self.h_region.fit, b[2])
		# self.smooth = Smooth(Smooth(self.data_region.fit, self.r_region.fit, b[1]), self.h_region.fit, b[2])
		# self.smooth.f1.set_p(0.002)
		Region.__init__(self, FunctionFit(f=self.fit), b[0], b[-1], name, False, None, path, filename, quiet)
		self.fit.get_p = lambda : self.r_region.fit.get_p() + self.h_region.fit.get_p()
		def set_p(*p):
			l = len(self.h_region.fit.get_p())
			self.r_region.fit.set_p(*p[:-l])
			self.h_region.fit.set_p(*p[-l:])
		self.fit.set_p = set_p
		self.fitted_region = Region(self.fit, k_min=b[1], k_max=b[-1], path=path, filename=filename, init=False)
		self.fitted_region.singular = self.data_region.singular	# hack (improving precision)
		self.fitted_region.cs = self.data_region.cs = self.cs	# actually should be before precompute, but this way it's smooth (not fair though)
		self.fitted_region.dcs = self.data_region.dcs = self.dcs
		if init:
			self.update_data()
			self.update_fit(quiet)
	
	def fit(self, k, *p):
		p0 = self.fit.get_p()
		self.fit.set_p(*self.fit.full_p(p))
		# if k < self.b[1]:
		# 	r = self.data_region.fit(k)
		# else:
		# 	r = self.smooth(k)
		r = self.smooth(k)
		self.fit.set_p(*p0)
		return r
	
	@cached_args_method()
	def cs(self, k, err=True):
		if not (self.k_min <= k <= self.k_max):
			if err: return val(0)
			return 0.
		if k < self.b[1]:
			r = self.data_region.fit(k)
		else:
			r = self.fit(k)
		if err:
			r = val(r, self.cs_err(k))
		return r
	
	@cached_args_method()
	def dcs(self, i, k):
		if k < self.b[1]:
			return 0.
		else:
			return self.fitted_region.fit.df(i, k)
	
	@cached_args_method()
	def cs_err(self, k):
		if k < self.b[1]:
			return self.data_region.cs_err(k)
		else:
			return self.fit.err(k)
	
	def _gen_method_integral_x(method):
		def integral_x(self, k, k0=None, k1=None, err=True):
			ret = sum(getattr(r, 'integral_'+method)(k, k0, k1, err) for r in [self.data_region, self.fitted_region])
			if method == 'f':
				ret -= float(e_Mp)
				# cs = lambda a: float(self.cs(a, err=False) + self.cst(a))
				cs = lambda a: float(self.cs(a, err=False))
			elif 'g' in method:
				if '_subtracted' in method:
					ret += k*float(IGDH)/(4*pi**2)
				# cs = lambda a: -float(self.cs(a, err=False) + self.csp(a))
				cs = lambda a: -float(self.cs(a, err=False))
			return ret + Region.singular(k, self.b[0], self.b[-1]) * cs(k) * k**2/(2*pi**2)
		integral_x.__name__ = 'integral_'+method
		return integral_x
	
	def _gen_method_dintegral_x(method):
		def dintegral_x(self, i, k, k0=None, k1=None, subtracted=True):
			# ret = sum(getattr(r, 'dintegral_'+method)(i, k, k0, k1) for r in [self.data_region, self.fitted_region])
			ret = getattr(self.fitted_region, 'dintegral_'+method)(i, k, k0, k1, subtracted=subtracted)
			singular = Region.singular(k, self.b[1], self.b[-1]) * self.dcs(i, k) * k**2/(2*pi**2)	# note b[1] instead of b[0] here (data_region does not account for singular part right now (see above))
			if 'g' in method:
				return ret - singular
			return ret + singular
		dintegral_x.__name__ = 'dintegral_'+method
		return dintegral_x
	
	for m in ['f', 'g', 'g_subtracted']:
		vars()['integral_'+m] = _gen_method_integral_x(m)
		vars()['dintegral_'+m] = _gen_method_dintegral_x(m)
	
	def integral_n(self, n, k0, k1, err=True):
		return self.data_region.integral_n(n, k0, k1, err=err) + self.fitted_region.integral_n(n, k0, k1, err=err)
	
	def plot_regge(self, xmin=None, xmax=10**7, dx=0.02, logx=True, **kwargs):
		kwargs.update({'title': 'Regge cross section', 'xylabel': r'"$\nu$ [GeV]" "$\sigma$ [$\mu\mathrm{b}$]"'})
		return self.plot(self.h_region.fit, 'cs_regge', xmin, xmax, dx, logx, **kwargs)
	
	def print_n(self, n, b=None, short=False):
		b = b or self.b
		if b[0] is None:
			b[0] = self.k_min
		if b[-1] is None:
			b[-1] = self.k_max
		if n == 1:
			print('GDH SR:')
		elif n == 2:
			print('Baldin SR:')
		else:
			print(str(n)+'th order SR:')
		if short:
			print('  ' + str_v(self.integral_n(n, b[0], b[-1])))
		else:
			s = ['  [%.1f, %.1f]'%tuple(b[:2]) + ' threshold region:  ',
				'  [%.1f, %.1f]'%tuple(b[1:3]) + ' resonance region:  ',
				'  [%.1f, %.1f]'%(b[2], b[-1]) + ' high energies:  ',
				'  [%.1f, %.1f]'%(b[0], b[-1]) + ' TOTAL:  ']
			m = max(map(len, s))
			s = [a + ' ' * (m - len(a)) for a in s]
			print('\n'.join(s[i] + str_v(self.integral_n(n, b[i1], b[i2]), d=5) for i, i1, i2 in zip(range(4), (0,1,2,0), (1,2,3,-1))))
		return self

