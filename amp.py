#!/usr/bin/env python

import os
from scipy import pi, nan

from val import val
from common import *
from compton import *
from fit import BornP
import integrate, tot, pol


# print '\n----------------------------- unpolarized: fit 1 ----------------------------'
#wt1 = tot.WholeRegion_1(quiet=True, pre_g=True)#.print_all()	# pre_g for wp error estimation...
wt1 = tot.WholeRegion_1(quiet=True)#.print_all()
# print '\n----------------------------- unpolarized: fit 2 ----------------------------'
wt2 = tot.WholeRegion_2(quiet=True)#.print_all()
# print "\n-------------------------- unpolarized: Ahrens' fit -------------------------"
wa = tot.AhrensWholeRegion(quiet=True)#.print_all()
# print '\n--------------------------------- polarized ---------------------------------'
wp = pol.WholeRegion(quiet=True)#.print_all()
# print '\n-----------------------------------------------------------------------------\n'


# def f(k, wt=wt1):	# mcb * GeV
# 	return wt.f(k)

# def g(k, wp=wp):	# mcb * GeV
# 	return wp.g(k)

def singular(cs, k, k_min=E0, k_max=E1):
	return integrate.Region.singular(k, k_min, k_max) * cs(k) * k**2/(2*pi**2)

def s2z(k, wt=wt1, wp=wp):
	fk, dfk, gk, dgk = wt.f(k), (lambda i: wt.df(i, k)), wp.g(k), (lambda i: wp.dg(i, k))
	refgk = fk.real*gk.real + fk.imag*gk.imag
	dsk = fk.real**2 + fk.imag**2 + gk.real**2 + gk.imag**2
	r = val(refgk/dsk)
	e1 = wt.dfi_err(lambda i: ((dfk(i).real*gk.real + dfk(i).imag*gk.imag)*dsk-(2*fk.real*dfk(i).real + 2*fk.imag*dfk(i).imag)*refgk)/(refgk*dsk))
	e2 = wp.dfi_err(lambda i: ((fk.real*dgk(i).real + fk.imag*dgk(i).imag)*dsk-(2*gk.real*dgk(i).real + 2*gk.imag*dgk(i).imag)*refgk)/(refgk*dsk))
	efk = (wt.data_region.f(k) + singular(lambda a: wt.data_region.cs(a, err=False), k, wt.data_region.k_min, wt.data_region.k_max) - float(e_Mp)) * 0.02 		#
	if k > wt.data_region.k_max: efk = efk.real + 0j		###
	# egk = (wp.data_region.g(k) - singular(lambda a: wp.data_region.cs(a, err=False), k, wp.data_region.k_min, wp.data_region.k_max) + k*float(IGDH)/(4*pi**2)) * 0.02 		#
	egk = (wt1.data_region.g(k) - singular(lambda a: wt1.data_region.cs(a, err=False), k, wp.data_region.k_min, wp.data_region.k_max) + k*float(IGDH)/(4*pi**2)) * 0.02 		#
	if k > wp.data_region.k_max: egk = egk.real + 0j		###
	e1t = ((efk.real*gk.real + efk.imag*gk.imag)*dsk - 2*(fk.real*efk.real + fk.imag*efk.imag)*refgk)/(refgk*dsk)	# threshold region contrib
	e2t = ((fk.real*egk.real + fk.imag*egk.imag)*dsk - 2*(gk.real*egk.real + gk.imag*egk.imag)*refgk)/(refgk*dsk)	# threshold region contrib
	r.err = (e1**2 + e2**2 + e1t**2 + e2t**2)**0.5
	return -2 * r

def ds(k, wt=wt1, wp=wp):	# nb
	fk, dfk, gk, dgk = wt.f(k), (lambda i: wt.df(i, k)), wp.g(k), (lambda i: wp.dg(i, k))
	r = val(fk.real**2 + fk.imag**2 + gk.real**2 + gk.imag**2)
	e1 = wt.dfi_err(lambda i: 2*fk.real*dfk(i).real + 2*fk.imag*dfk(i).imag)
	e2 = wp.dfi_err(lambda i: 2*gk.real*dgk(i).real + 2*gk.imag*dgk(i).imag)
	# efk = wt.data_region.re_f(k).err + 1j*wt.data_region.im_f(k).err 		#
	# egk = wp.data_region.re_g(k).err + 1j*wp.data_region.im_g(k).err 		#
	# efk = (wt.data_region.f(k) - float(e_Mp)) * 0.02 		#
	# egk = (wp.data_region.g(k) + k*float(IGDH)/(4*pi**2)) * 0.02 		#
	efk = (wt.data_region.f(k) + singular(lambda a: wt.data_region.cs(a, err=False), k, wt.data_region.k_min, wt.data_region.k_max) - float(e_Mp)) * 0.02 		#
	if k > wt.data_region.k_max: efk = efk.real + 0j		###
	# egk = (wp.data_region.g(k) - singular(lambda a: wp.data_region.cs(a, err=False), k, wp.data_region.k_min, wp.data_region.k_max) + k*float(IGDH)/(4*pi**2)) * 0.02 		#
	egk = (wt1.data_region.g(k) - singular(lambda a: wt1.data_region.cs(a, err=False), k, wp.data_region.k_min, wp.data_region.k_max) + k*float(IGDH)/(4*pi**2)) * 0.02 		#
	if k > wp.data_region.k_max: egk = egk.real + 0j		###
	e1t = 2*fk.real*efk.real + 2*fk.imag*efk.imag		# threshold region contrib
	e2t = 2*gk.real*egk.real + 2*gk.imag*egk.imag		# threshold region contrib
	r.err = (e1**2 + e2**2 + e1t**2 + e2t**2)**0.5
	return r / (10*float(hc)**2)

def refg(k, wt=wt1, wp=wp):	# nb
	fk, dfk, gk, dgk = wt.f(k), (lambda i: wt.df(i, k)), wp.g(k), (lambda i: wp.dg(i, k))
	r = val(fk.real*gk.real + fk.imag*gk.imag)
	e1 = wt.dfi_err(lambda i: dfk(i).real*gk.real + dfk(i).imag*gk.imag)
	e2 = wp.dfi_err(lambda i: fk.real*dgk(i).real + fk.imag*dgk(i).imag)
	efk = (wt.data_region.f(k) + singular(lambda a: wt.data_region.cs(a, err=False), k, wt.data_region.k_min, wt.data_region.k_max) - float(e_Mp)) * 0.02 		#
	if k > wt.data_region.k_max: efk = efk.real + 0j		###
	# egk = (wp.data_region.g(k) - singular(lambda a: wp.data_region.cs(a, err=False), k, wp.data_region.k_min, wp.data_region.k_max) + k*float(IGDH)/(4*pi**2)) * 0.02 		#
	egk = (wt1.data_region.g(k) - singular(lambda a: wt1.data_region.cs(a, err=False), k, wp.data_region.k_min, wp.data_region.k_max) + k*float(IGDH)/(4*pi**2)) * 0.02 		#
	if k > wp.data_region.k_max: egk = egk.real + 0j		###
	e1t = efk.real*gk.real + efk.imag*gk.imag		# threshold region contrib
	e2t = fk.real*egk.real + fk.imag*egk.imag		# threshold region contrib
	r.err = (e1**2 + e2**2 + e1t**2 + e2t**2)**0.5
	return r / (10*float(hc)**2)

def plot(f, prefix='plot', xmin=0., xmax=2.5, dx=0.002, logx=False, filename=None, raw=True, **kwargs):
	path = os.path.join('out', 'amp')
	make_sure_path_exists(path)
	filename = os.path.join(path, prefix+('_logx' if logx else '')+'.out') if filename is None else filename
	xf = arange(xmin, xmax, dx, logx)
	fx = map_verbose(f, xf, 'plotting to "'+filename+'"')
	data = [xf] + list(zip(*fx)) if is_iterable(fx[0]) else [xf, fx]
	data = zip(*[d for d in zip(*data) if nan not in d])	# dirty hack
	kwargs['plots'] = kwargs.get('plots') or ['1 '+str(i) for i in range(2, len(data)+1)]
	plot_to_file(data, filename=filename, **({} if raw else kwargs))

def plot_s2z(xmin=0.002, xmax=2.5, dx=0.002, postfix='', **kwargs):
	plot(lambda k: s2z(k, wt=wt1), 's2z_1'+(postfix and '_'+postfix), xmin, xmax, dx, **kwargs)
	plot(lambda k: s2z(k, wt=wt2), 's2z_2'+(postfix and '_'+postfix), xmin, xmax, dx, **kwargs)

def plot_ds(xmin=0.002, xmax=2.5, dx=0.002, postfix='', **kwargs):
	plot(lambda k: ds(k, wt=wt1), 'ds_1'+(postfix and '_'+postfix), xmin, xmax, dx, **kwargs)
	plot(lambda k: ds(k, wt=wt2), 'ds_2'+(postfix and '_'+postfix), xmin, xmax, dx, **kwargs)

def plot_2refg(xmin=0.002, xmax=2.5, dx=0.002, postfix='', **kwargs):
	plot(lambda k: -2*refg(k, wt=wt1), '2refg_1'+(postfix and '_'+postfix), xmin, xmax, dx, **kwargs)
	plot(lambda k: -2*refg(k, wt=wt2), '2refg_2'+(postfix and '_'+postfix), xmin, xmax, dx, **kwargs)

def plot_f(xmin=0.002, xmax=2.5, dx=0.002, postfix='', **kwargs):
	l = lambda c: (c.real, c.imag)
	plot(lambda k: l(wt1.f(k)), 'f_1'+(postfix and '_'+postfix), xmin, xmax, dx, **kwargs)
	plot(lambda k: l(wt2.f(k)), 'f_2'+(postfix and '_'+postfix), xmin, xmax, dx, **kwargs)

def plot_g(xmin=0.002, xmax=2.5, dx=0.002, postfix='', **kwargs):
	plot(lambda k: (lambda c: (c.real, c.imag))(wp.g(k)), 'g'+(postfix and '_'+postfix), xmin, xmax, dx, **kwargs)

def plot_cs_born(xmin=0., xmax=4., dx=0.001, **kwargs):	# Born contribution for polarized CS
	lo = BornP()
	plot(lo, '', xmin, xmax, dx, filename=os.path.join('out', 'pol', 'born.out'))


if __name__ == '__main__':
	
	# plot_ds()
	# plot_2refg()
	# plot_s2z()
	# exit()
	wt1.plot_cs().plot_f().plot_cs(xmin=wt1.k_min, xmax=None, logx=True)
	wt2.plot_cs().plot_f().plot_cs(xmin=wt2.k_min, xmax=None, logx=True)
	# wt1.plot_all()
	# wt2.plot_all()
	wa.plot_f()
	wp.plot_cs().plot_g().plot_cs(xmin=wp.k_min, xmax=None, logx=True)
	wp.plot_1(xmax=1000.).plot_3(xmax=1000.).plot_5(xmax=1000.)
	wp.plot_im_g()
	# wp.plot_all()
	plot_cs_born()
	plot_f()
	plot_f(xmin=0.1, xmax=0.2, dx=0.0004, postfix='small')
	plot_g()
	plot_g(xmin=0.1, xmax=0.2, dx=0.0004, postfix='small')
	plot_ds()
	plot_ds(xmin=0.0002, xmax=0.4, dx=0.0002, postfix='small')
	# plot_s2z()
	# plot_s2z(xmin=0.0002, xmax=0.4, dx=0.0002, postfix='small')
	plot_2refg()
	plot_2refg(xmin=0.0002, xmax=0.4, dx=0.0002, postfix='small')
	# plot_s2z(postfix='0')
