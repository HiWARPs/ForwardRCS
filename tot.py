#!/usr/bin/env python

from compton import *
import integrate


class WholeRegion(integrate.WholeRegion):
	
	# #..........................................................................
	# def update_data(self, filename=None, k_min=None, k_max=None):
	# 	self.fitted_region.update_data(filename, k_min and max(k_min, self.fitted_region.k_min), k_max)
	# 	return integrate.WholeRegion.update_data(self, filename, k_min, k_max)
	# #..........................................................................
	
	def __init__(self, *args, **kwargs):
		integrate.WholeRegion.__init__(self, *args, pre_g=False, **kwargs)
		# integrate.WholeRegion.__init__(self, *args, **kwargs)
		# #-----------------------------------------------------------------------------------------------------------------
		# self.b = b
		# self.data_region = ThresholdRegion(b[0], b[1], path=path, filename=filename, pre_f=False, pre_g=False)
		# self.r_region = ResonanceRegion(k_min=b[1], k_max=b[2], path=path, filename=filename, init=False, nres=nres, bgd=bgd)
		# self.h_region = AsymptoticRegion(k_min=b[2], k_max=b[-1], path=path, filename=filename, init=False)
		# self.data_region.singular = lambda k, k0, k1: 0.	# hack (improving precision)
		# if pre_f: self.data_region.precompute_f(quiet=quiet)
		# if pre_g:
		# 	self.data_region.precompute_g(quiet=quiet)
		# 	self.data_region.precompute_g_subtracted(quiet=quiet)
		# self.smooth = Smooth(self.r_region.fit, self.h_region.fit, b[2])
		# Region.__init__(self, FunctionFit(f=self.fit), b[0], b[-1], name, False, None, path, filename, quiet)
		# self.fit.get_p = lambda : self.r_region.fit.get_p() + self.h_region.fit.get_p()
		# def set_p(*p):
		# 	l = len(self.h_region.fit.get_p())
		# 	self.r_region.fit.set_p(*p[:-l])
		# 	self.h_region.fit.set_p(*p[-l:])
		# self.fit.set_p = set_p
		# self.fitted_region = Region(self.fit, k_min=b[1], k_max=b[-1], path=path, filename=filename, init=False)
		# self.fitted_region.singular = self.data_region.singular	# hack (improving precision)
		# self.fitted_region.cs = self.data_region.cs = self.cs	# actually should be before precompute, but this way it's smooth (not fair though)
		# if init:
		# 	self.update_data()
		# 	self.update_fit(quiet)
		# #-----------------------------------------------------------------------------------------------------------------
		self.gdh = self.sr_1 = self.sr_3 = self.sr_5 = None
		self.g = self.re_g = self.im_g = self.integral_g = None
	
	# #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	# def fit(self, k, *p):
	# 	p0 = self.fit.get_p()
	# 	self.fit.set_p(*self.fit.full_p(p))
	# 	if k < self.b[1]:
	# 		r = self.data_region.fit(k)
	# 	else:
	# 		r = self.smooth(k)
	# 	self.fit.set_p(*p0)
	# 	return r
	
	# cs = integrate.Region.cs
	
	# def cs_err(self, k):
	# 	if k < self.b[1]:
	# 		return self.data_region.cs_err(k)
	# 	else:
	# 		return self.fit.err(k)
	
	# def integral_f(self, k, k0=None, k1=None):
	# 	cs = lambda a: float(self.cs(a) + self.cst(a))
	# 	singular = Region.singular(k, self.b[0], self.b[-1]) * cs(k) * k**2/(2*pi**2)	#
	# 	return singular + self.data_region.integral_f(k, k0, k1) + self.fitted_region.integral_f(k, k0, k1) - e_Mp	#
	
	# def integral_g(self, k, k0=None, k1=None):
	# 	cs = lambda a: float(self.cs(a) + self.csp(a))
	# 	singular = Region.singular(k, self.b[0], self.b[-1]) * (-cs(k)) * k**2/(2*pi**2)	#
	# 	return singular + self.data_region.integral_g(k, k0, k1) + self.fitted_region.integral_g(k, k0, k1)		#
	
	# def integral_g_subtracted(self, k, k0=None, k1=None):
	# 	cs = lambda a: float(self.cs(a) + self.csp(a))
	# 	singular = Region.singular(k, self.b[0], self.b[-1]) * (-cs(k)) * k**2/(2*pi**2)	#
	# 	return singular + self.data_region.integral_g_subtracted(k, k0, k1) + self.fitted_region.integral_g_subtracted(k, k0, k1) + k*IGDH/(4*pi**2)		#
	
	# def integral_n(self, n, k0, k1):
	# 	return self.data_region.integral_n(n, k0, k1) + self.fitted_region.integral_n(n, k0, k1)
	# #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
	def integral_c(self, k1):
		r = Region(self.h_region.fit, k_min=self.k_min, k_max=self.k_max, filename=self.filename, init=False)
		r.err = [1.]*len(r.err)	# temporary hack
		return self.data_region.integral_c(k1) + self.fitted_region.integral_c(k1) - r.integral_c(k1)[0]
	
	def print_2(self, b=None, short=False):
		return self.print_n(2, b, short)
	
	print_baldin = print_2
	
	def print_4(self, b=None, short=False):
		return self.print_n(4, b, short)
	
	def print_6(self, b=None, short=False):
		return self.print_n(6, b, short)
	
	def print_all(self):
		return self.print_2().print_4().print_6()
	
	def plot_all(self):
		return self.plot_cs().plot_f().plot_2().plot_4().plot_6().plot_cs(xmin=self.k_min, xmax=None, logx=True)


class AhrensWholeRegion(integrate.ThresholdRegion, WholeRegion):
	
	# py3 port: Python 2 old-style classes resolve attributes depth-first (ThresholdRegion -> Region -> WholeRegion), new-style (C3)
	# would pick WholeRegion.cs and tot.WholeRegion.integral_c instead of Region's; restore the Python 2 resolution
	cs = integrate.Region.cs
	integral_c = integrate.Region.integral_c
	
	def __init__(self, k_min=E0, k_max=1E10, name='ahrens', path='tot', filename='ahrens_fit.dat', quiet=False):
		integrate.ThresholdRegion.__init__(self, k_min=k_min, k_max=k_max, name=name, path=path, filename=filename, pre_g=False, quiet=quiet)
		self.b = [self.k_min, E1, E2, self.k_max]


class WholeRegion_1(WholeRegion):
	
	def __init__(self, name='1', init=True, quiet=False, **kwargs):
		WholeRegion.__init__(self, name=name, path='tot', filename='mdata.dat', init=False, quiet=quiet, **kwargs)
		p = [1.2136,   0.1177,  522.6315,   1.4128,   0.0832,   40.2570,
			 1.4960,   0.1366,  161.9159,   1.6975,   0.0194,   18.3686,
			 1.6491,   0.1345,   82.7881,   1.8942,   0.3068,   31.3450,
			 0.3936, -10.2068,   69.4929,  31.5777,  32.2989]
		self.r_region.fit.set_p(*p)
		if init:
			self.update_data()
			self.update_fit(quiet=quiet)
		self.fit.pcov = self.fit.get_pcov((self.r_region.x + self.h_region.x, self.r_region.err + self.h_region.err))	# the points and errors of the fit


class WholeRegion_2(WholeRegion):
	
	def __init__(self, name='2', init=True, quiet=False, **kwargs):
		WholeRegion.__init__(self, name=name, path='tot', filename='sdata.dat', init=False, quiet=quiet, **kwargs)
		p = [1.21471, 0.10218, 0.10218**2/4*495.9433, 1.41770, 0.13359, 0.13359**2/4*78.28695,
			 1.49898, 0.09562, 0.09562**2/4*138.9695, 1.58756, 0.11560, 0.11560**2/4*51.76039,
			 1.67504, 0.10624, 0.10624**2/4*78.73006, 1.87538, 0.37208, 0.37208**2/4*33.18966,
			 0.19899, -6.7412, 49.62000, 70.8614, 8.85919]
		self.r_region.fit.set_p(*p)
		if init:
			self.update_data()
			self.update_fit(quiet=quiet)
		self.fit.pcov = self.fit.get_pcov((self.r_region.x + self.h_region.x, self.r_region.err + self.h_region.err))	# the points and errors of the fit


if __name__ == '__main__':
	
	print ("\n------------------ fit I ------------------")
	WholeRegion_1().print_all().plot_all()
	exit()
	print ("\n----------------- fit II ------------------")
	WholeRegion_2().print_all().plot_all()
	print ("\n-------------- Ahrens' fit  ---------------")
	AhrensWholeRegion().print_all().plot_f()
