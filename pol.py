#!/usr/bin/env python

from common import *
from compton import *
import integrate, numpy, scipy
from scipy.integrate import quad


data = integrate.Data('pol')


class ThresholdRegion(integrate.ThresholdRegion):
	
	def __init__(self, *args, **kwargs):
		integrate.ThresholdRegion.__init__(self, *args, **kwargs)
		kwargs.update(path='tot', filename='maid.dat', pre_f=False, pre_g=False, quiet=True)
		self.fit_tot = integrate.ThresholdRegion(*args, **kwargs).fit 	# for error estimation
	
	def cs_err(self, k):
		return self.fit_tot(k) * 0.02
	
	def _int_err(self, df, m, k, k0, k1, subtracted=True):
		return abs(self._int(lambda p: self.fit_tot(p)*0.02, m, k, k0, k1, subtracted, singular=integrate.Region.singular).val)


class WholeRegion(integrate.WholeRegion):
	
	def update_fit(self, quiet=False):
		self._cache = {}
		self.r_region.update_fit(quiet=quiet)
		#self.fitted_region.update_fit(quiet=quiet)
		return self
	
	def __init__(self, init=True, quiet=False):
		integrate.WholeRegion.__init__(self, name='', path='pol', filename='data.dat', init=init, quiet=quiet, bgd=False, pre_f=False, ThresholdRegion=ThresholdRegion)
		self.h_region.fit.set_p(-34.1, -0.58*2, 209.4, -1.66*2)
		self.smooth.m = self.r_region.k_max = self.h_region.k_min = self.b[2] = 2.03005
		self.smooth.set_p(0.15)
		# self.r_region.singular = self.h_region.singular = self.data_region.singular
		# self.r_region.fit = self.h_region.fit = self.fit
		# self.r_region.fit.get_pcov(zip(*[(x, e) for x, e in zip(data.x, data.err) if self.r_region.k_min <= x]))	# NOTE we use here all the points (not only up to 2 GeV)
		# self.r_region.fit.get_pcov(zip(*[(x, e) for x, e in zip(data.x, data.err) if self.r_region.k_min <= x <= self.r_region.k_max]))
		self.r_region.fit.get_pcov(zip(*[(x, e) for x, e in zip(data.x, data.err) if x <= self.r_region.k_max]))
		self.h_region.fit.pcov = numpy.diag([5.7, 0.23*2, 29., 0.22*2])**2
		#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		# self.h_region.fit.f = self.h_f
		# self.h_region.fit.pcov = scipy.linalg.block_diag(self.r_region.fit.get_pcov(), self.h_region.fit.get_pcov())
		# self.h_region.fit.set_p(*(self.r_region.fit.get_p() + self.h_region.fit.get_p()))
		#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		# self.r_region.update_data(filename='data/pol/data1.dat')
		# self.h_region.fit.fit(*zip(*[a for a in zip(self.x, self.y, self.err) if a[0]>=1.45]), quiet=quiet)
		# self.h_region.fit.perr = [5.7, 0.23*2, 29., 0.22*2]
		# self.h_region.fit.perr = [5.7, 0, 29., 0]
		# self.h_region.fit.pcov = numpy.diag(self.h_region.fit.perr)**2
		# self.fit.pcov = self.get_pcov()
		# self.fit.pcov = self.fit.get_pcov(zip(*[(x, e) for x, e in zip(data.x, data.err) if self.r_region.k_min <= x <= self.h_region.k_max]))
		# self.fit.pcov = scipy.linalg.block_diag(self.fit.pcov[:-4,:-4], numpy.diag([5.7, 0, 29., 0])**2)
		# self.fit.pcov = scipy.linalg.block_diag(self.fit.pcov[:-4,:-4], numpy.diag([5.7, 0.23*2, 29., 0.22*2])**2)
		# self.fit.pcov = scipy.linalg.block_diag(self.r_region.fit.pcov, numpy.diag([5.7, 0.23*2, 29., 0.22*2])**2)
		self.fit.pcov = scipy.linalg.block_diag(self.r_region.fit.pcov, self.h_region.fit.pcov)
#		I1 = quad(lambda p: float(self.dcs(-4, p))/p, self.k_min, self.k_max)[0]
#		I2 = quad(lambda p: float(self.dcs(-2, p))/p, self.k_min, self.k_max)[0]
#		self.fit.pcov[-4,-2] = self.fit.pcov[-2,-4] = (IGDH.err**2/4 - (I1*5.7)**2 - (I2*29)**2)/(2*I1*I2)
		# self.fit.pcov = scipy.linalg.block_diag(self.fit.pcov[:-4,:-4], numpy.diag([0, 0, 0, 0])**2)
		# self.fit.pcov[-4:,-4:] = numpy.diag([5.7, 0, 29., 0])**2
		# self.fit.pcov[-1,:] = numpy.zeros(self.fit.pcov.shape[0])
		# self.fit.pcov[-3,:] = numpy.zeros(self.fit.pcov.shape[0])
		# self.fit.pcov[:,-1] = numpy.zeros(self.fit.pcov.shape[0])
		# self.fit.pcov[:,-3] = numpy.zeros(self.fit.pcov.shape[0])
		# print self.fit.pcov
		self.baldin = self.sr_2 = self.sr_4 = self.sr_6 = None
		self.f = self.re_f = self.im_f = self.integral_f = None
	
	# def get_pcov(self):
	# 	x, err = zip(*[(x, e) for x, e in zip(data.x, data.err) if self.r_region.k_min <= x <= self.h_region.k_max])
	# 	m = len(self.fit.get_p())
	# 	J = numpy.array([[self.fit.df(i, k) for i in range(m)] for k in x])
	# 	# J = 
	# 	err += [IGDH.err]
	# 	return numpy.linalg.inv(J.T.dot(numpy.diag(numpy.array(err)**-2)).dot(J))
	
	def print_1(self, b=None, short=False):
		return self.print_n(1, b, short)
	
	print_gdh = print_1
	
	def print_3(self, b=None, short=False):
		return self.print_n(3, b, short)
	
	def print_5(self, b=None, short=False):
		return self.print_n(5, b, short)
	
	def print_all(self):
		return self.print_1().print_3().print_5()
	
	def plot_all(self):
		return self.plot_cs().plot_g().plot_1().plot_3().plot_5().plot_cs(xmin=self.k_min, xmax=None, logx=True)


if __name__ == '__main__':
	
	WholeRegion().print_all().plot_all()
