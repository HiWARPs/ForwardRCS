#!/usr/bin/env python

# import numpy
# from scipy import log, exp, pi
# from mpmath import log, exp, pi
# from scipy.misc import derivative
# from mpmath import diff as derivative
# from scipy.optimize import curve_fit, leastsq

# import mpmath
# array = lambda a: numpy.array(a, dtype=mpmath.mpf)

from val import val
import compton
from compton import *
from common import *


def chisqr(f, xdata, ydata, errdata):
	return sum([((y-f(x))/err)**2 for x, y, err in zip(xdata, ydata, errdata)])


class FunctionFit:
	
	def __init__(self, *p, **kwargs):
		self.pcov = None
		self.set_p(*p)
		if 'f' in kwargs:	# py3 port
			self.f = kwargs['f']
	
	def f(self, k, *p):
		return 0.
	
	def __call__(self, k, *p):
		p = self.full_p(p)
		if type(k) == numpy.ndarray:
			# if type(k[0]) == mpmath.mpf:
			# 	return numpy.array([self.f(x, *p) for x in k], dtype=object)
			return numpy.array([self.f(x, *p) for x in k], dtype=float)
		else:
			return self.f(k, *p)
	
	def get_p(self):	# let's always return immutable tuple
		return tuple(self.p)
	
	def set_p(self, *p):
		self.p = p
		# self.pcov = None
		# self.perr = self.perr[:len(p)] + ([0.]*len(p))[len(self.perr):]
	
	def df(self, i, k, *p):
		p = self.full_p(p)
		# return derivative(lambda x: self(k, *p[:i]+(x,)+p[i+1:]), p[i], dx=1e-5, order=3)
		return derivative(lambda x: self(k, *p[:i]+(x,)+p[i+1:]), p[i])
	
	def get_pcov(self, xe=None):
		if xe is not None:
			x, err = xe
			m = len(self.get_p())
			# J = mpmath.matrix([[self.df(i, k) for i in range(m)] for k in x])
			J = numpy.array([[self.df(i, k) for i in range(m)] for k in x])
			# self.pcov = numpy.array(((J.T * mpmath.matrix(numpy.diag(numpy.array(err, dtype=float)**-2)) * J)**-1).tolist())	# hack with mpmath matrices
			self.pcov = numpy.linalg.inv(J.T.dot(numpy.diag(numpy.array(err)**-2)).dot(J))
		return self.pcov
	
	def err(self, k, *p):
		dcs = numpy.array([self.df(i, k, *p) for i in range(len(self.get_p()))])
		return (dcs.dot(self.get_pcov()).dot(dcs))**0.5
	
	def val(self, k, *p):
		return val(self(k, *p), self.err(k, *p))
	
	def get_perr(self):
		return tuple(numpy.sqrt(numpy.diag(self.get_pcov())))
	
	def full_p(self, p):
		return tuple(list(p) + list(self.get_p())[len(p):])
	
	def fit(self, xdata, ydata, errdata, quiet=False):	# let's return self
		if not quiet: print('fitting ...')
		p, pcov = curve_fit(self, xdata, ydata, p0=self.get_p(), sigma=errdata, absolute_sigma=True)
		self.set_p(*p)
		self.pcov = pcov
		if not quiet: print('chisqr/point = %.2f' % (chisqr(self, xdata, ydata, errdata)/len(xdata)), '(for %d points)' % len(xdata))
		return self
	
	def print_p(self, s=None, p=None, sep='\n  '):
		if p is None:
			p = zip(self.get_p(), self.get_perr())
		if s is not None:
			print(s)
		print('  ' + sep.join(map(str_v, p)))


class BackgroundFit(FunctionFit):
	
	def f(self, k, *p):
		W0 = flab2cm(E0)
		w = flab2cm(k)
		return sum([p[i]*(w-W0)**(i-2) for i in range(len(p))])


class ResonanceFit(FunctionFit):
	
	def f(self, k, *p):
		m, width, coupling = p
		w = flab2cm(k)
		a = width**2 / 4
		#return coupling * a / ((w - m)**2 + a)
		return coupling / ((w - m)**2 + a)
		#return coupling / ((w**2 - m**2)**2 + m**2 * width**2)
	
	def print_p(self, s=None, *args, **kwargs):
		p, perr = self.get_p(), self.get_perr()
		p1 = 1000*p[0], 1000*p[1], p[2]/(p[1]**2/4)	# [MeV], [MeV], [mcb]
		perr1 = 1000*perr[0], 1000*perr[1], perr[2]/(p[1]**2/4)#*((4/p[1]**2)**2 + (8*p[2]/p[1]**3)**2)**0.5	# is it Ok?
		FunctionFit.print_p(self, s, zip(p1, perr1), *args, **kwargs)


class ResonanceRegionFit(FunctionFit):
	
	def __init__(self, nres=6, bgd=True):
		self.resonances = [ResonanceFit() for i in range(nres)]
		self.background = BackgroundFit() if bgd else FunctionFit()
		self.bgd = bgd
		if nres == 6 and bgd:	# very dirty hack
			p = [1.2136,   0.1177,  522.6315,   1.4128,   0.0832,   40.2570,
				 1.4960,   0.1366,  161.9159,   1.6975,   0.0194,   18.3686,
				 1.6491,   0.1345,   82.7881,   1.8942,   0.3068,   31.3450,
				 0.3936, -10.2068,   69.4929,  31.5777,  32.2989]
		elif nres == 6 and not bgd:
			p = [1.2102,   0.1193,   294.3772,  1.4046,   0.4939,   -147.3924,
				1.4608,   0.2395,   136.1701,   1.5855,   0.1117,   -72.6674,
				1.6165,   0.3605,   117.5453,   1.7525,   0.1050,   -37.4879]
		else:
			p = [0.]*(3*nres + (5 if bgd else 0))
		FunctionFit.__init__(self, *p)
	
	def f(self, k, *p):
		p0 = self.get_p()
		self.set_p(*self.full_p(p))
		r = sum([r(k) for r in self.resonances]) + self.background(k)
		self.set_p(*p0)
		return r
	
	def get_p(self):
		return tuple(sum([list(r.get_p()) for r in self.resonances], []) + list(self.background.get_p()))
	
	def set_p(self, *p):
		n = len(self.resonances)
		for i in range(n):
			self.resonances[i].set_p(*p[3*i:3*i+3])
		self.background.set_p(*p[3*n:])
		# self.perr = self.perr[:len(p)] + ([0.]*len(p))[len(self.perr):]
	
	def fit_bgrd(self, xdata, ydata, errdata, quiet=False):
		if self.bgd:
			self.background.fit(xdata, [y - self(x) + self.background(x) for x, y in zip(xdata, ydata)], errdata, quiet=quiet)
			# self.perr = self.perr[:3*len(self.resonances)] + self.background.perr
		return self
	
	def fit_ress(self, xdata, ydata, errdata, quiet=False):
		n = len(self.resonances)
		ff = FunctionFit(*self.get_p()[:3*n])
		ff.f = lambda k, *p: self.f(k, *self.full_p(p))
		ff.fit(xdata, ydata, errdata, quiet=quiet)
		self.set_p(*self.full_p(ff.get_p()))
		# self.perr = ff.perr + self.perr[3*n:]
		# for i in range(n):
			# self.resonances[i].perr = self.perr[3*i:3*i+3]
		return self
	
	def fit(self, xdata, ydata, errdata, bgd_fit=True, ncycles=3, quiet=False):
		if not quiet: print('fitting ...')
		for i in range(ncycles):
			self.fit_ress(xdata, ydata, errdata, quiet=True)
			if bgd_fit: self.fit_bgrd(xdata, ydata, errdata, quiet=True)
		# ff = FunctionFit(*self.get_p())
		# ff.f = lambda k, *p: self.f(k, *self.full_p(p))
		# ff.fit(xdata, ydata, errdata, quiet=True)
		# self.set_p(*ff.get_p())
		# self.pcov = ff.pcov
		# self.perr = ff.perr
		self.pcov = self.get_pcov((xdata, errdata))
		if not quiet:
			#print self.get_p()
			print('chisqr/point = %.2f' % (chisqr(self, xdata, ydata, errdata)/len(xdata)), '(for %d points)' % len(xdata))
			self.print_p()
		return self
	
	def print_p(self):
		p = list(map(val, zip(self.get_p(), self.get_perr())))	# py3 port
		n = len(self.resonances)
		print(r'resonances parameters (M [MeV], \Gamma [MeV], A [mcb]):')
		for i in range(n):
			print(' ', '   '.join('%.4f +- %.4f'%tuple(a) for a in p[3*i:3*i+3]))
		if self.bgd:
			print('background parameters ([mcb * GeV^n]):')
			print(' ', '   '.join('%.4f +- %.4f'%tuple(a) for a in p[3*n:]))


class AsymptoticRegionFit(FunctionFit):
	
	def __init__(self):
		FunctionFit.__init__(self, 71.*2, 0.075*2, 120.*2, -0.46*2)
	
	def f(self, k, *p):
		W = flab2cm(k)
		return 0.5*(p[0]*W**p[1] + p[2]*W**p[3])
	
	def fit(self, xdata, ydata, errdata, quiet=False):
		FunctionFit.fit(self, xdata, ydata, errdata, quiet=quiet)
		if not quiet:
			self.print_p()
		return self
	
	def print_p(self, s='Regge parameters (C1, p1, C2, p2):', *args, **kwargs):
		FunctionFit.print_p(self, s, *args, **kwargs)


class Smooth(FunctionFit):
	
	def f(self, k, a):
		numpy.seterr(all='ignore')
		e = 1./(1. + exp((k - self.m)/a)) # ^ just ignore this warning about overflow in exp ...
		r = e * self.f1(k) + (1-e) * self.f2(k)
		numpy.seterr(all='warn')
		return r
	
	def __init__(self, f1, f2, m, a0=0.1):
		self.f1, self.f2, self.m = f1, f2, m
		FunctionFit.__init__(self, a0)


class BornT(FunctionFit):	# Born contribution for unpolarized CS
	
	def fp(self, k):
		M = float(Mp)
		m = float(Mpi0)
		C = float(alpha)/(4*pi) * (float(compton.gA)*M/float(fpi))**2
		s = M**2 + 2*M*k
		if s < (M+m)**2:
			return 0.
		else:
			xN = (s + M**2 - m**2) / (2*s)
			l = (s - (M + m)**2)**0.5 * (s - (M-m)**2)**0.5 / (2*s)
			r = (k**2 - (m/M)**2*xN*s) * log((xN+l)/(xN-l)) + 2*l * (k**2*(xN-2)+s*(m/M)**2)
			return pi * C * r / (M * k**3)
	
	def fn(self, k):
		M = float(Mn)
		m = float(Mpip)
		C = float(alpha)/(4*pi) * (float(compton.gA)*M/float(fpi))**2
		s = M**2 + 2*M*k
		if s < (M+m)**2:
			return 0.
		else:
			xN = (s + M**2 - m**2) / (2*s)
			xpi = (s - M**2 + m**2) / (2*s)
			l = (s - (M + m)**2)**0.5 * (s - (M-m)**2)**0.5 / (2*s)
			r = -xpi*s*(m/M)**2 * log((xpi+l)/(xpi-l)) + 2*l*(xN*k**2 +s*(m/M)**2)
			return 2 * pi * C * r / (M * k**3)
	
	def f(self, k):
		return 10**4 * float(hc)**2 * (self.fp(k) + self.fn(k)) * 0.5


class BornP(FunctionFit):	# Born contribution for polarized CS
	
	def fp(self, k):
		M = float(Mp)
		m = float(Mpi0)
		#C = 0.096127619
		C = float(alpha)/(4*pi) * (float(compton.gA)*M/float(fpi))**2
		s = M**2 + 2*M*k
		if s < (M+m)**2:
			return 0.
		else:
			xN = (s + M**2 - m**2) / (2*s)
			l = (s - (M + m)**2)**0.5 * (s - (M-m)**2)**0.5 / (2*s)
			r = -(2*xN * s/M**2 + 1 - k/M) * log((xN+l)/(xN-l)) + 2*l * (k/M*(xN-2) + s/M**2*(xN+2))
			return pi * C * r / k**2
	
	def fn(self, k):
		M = float(Mn)
		m = float(Mpip)
		fpi = 0.0924
		C = float(alpha)/(4*pi) * (float(compton.gA)*M/float(fpi))**2
		s = M**2 + 2*M*k
		if s < (M+m)**2:
			return 0.
		else:
			xN = (s + M**2 - m**2) / (2*s)
			xpi = (s - M**2 + m**2) / (2*s)
			l = (s - (M + m)**2)**0.5 * (s - (M-m)**2)**0.5 / (2*s)
			r = (m/M)**2 * log((xpi+l)/(xpi-l)) - 2*l*(s/M**2*xpi - k/M*xN)
			return 2 * pi * C * r / k**2
	
	def f(self, k):
		return - 10**4 * float(hc)**2 * (self.fp(k) + self.fn(k)) * 0.5


if __name__ == '__main__':
	
	xdata, ydata, errdata = read_xye('data/tot/mdata.dat', xmin=2.)
	AsymptoticRegionFit().fit(xdata, ydata, errdata)
	
	print('\n------------------ fit I -----------------')
	xdata, ydata, errdata = read_xye('data/tot/mdata.dat', xmin=E1, xmax=2.)
	f = ResonanceRegionFit().fit(xdata, ydata, errdata)
	
	print('\n------------------ fit II ----------------')
	xdata, ydata, errdata = read_xye('data/tot/sdata.dat', xmin=E1, xmax=2.)
	f = ResonanceRegionFit().fit(xdata, ydata, errdata)
	
	bgdata = [f.background(x) for x in xdata]
	resdata = [[r(x)+f.background(x) for x in xdata] for r in f.resonances]
	fdata = [f(x) for x in xdata]
	data = [xdata, ydata, errdata, fdata, bgdata] + resdata
	plots = ['1 2 yerr=3 label=data marker=. linestyle=None', '1 4 label=fit linewidth=2', '1 5 label=background'] + ['1 '+str(i+6) for i in range(len(resdata))]
	plot_to_file(data, plots, 'outfit')
