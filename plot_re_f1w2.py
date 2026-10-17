#!/usr/bin/env python3

import tot

w = tot.WholeRegion_1(quiet=True)

def re_f1(k2, **kw):
	# return ref1.re_f1(w, k2**0.5, **kw)
	if k2 == 0:
		return w.baldin
	return tot.hc * (w.re_f(k2**0.5)/k2 - tot.e_Mp/k2)

def f_to_f1(dataset):
	x, y, err = dataset
	x = [float(k**2) for k in x]
	y = [float(tot.hc * (f - tot.e_Mp) / k2) for k2, f in zip(x, y)]
	err = [float(tot.hc * (f - tot.e_Mp) / k2) for k2, f in zip(x, err)]
	return [x, y, err]

# print re_f1(0, subtracted=False)
# print w.integral(2, 0, err=True, subtracted=False, const=True) * tot.hc
# print w.integral_n(2, w.k_min, w.k_max)

# print re_f1(0)
# exit()

from common import *

xmin = 0.005
x = arange(xmin, 10, 0.01, logx=True)
y, err = zip(*map(re_f1, x))

damashek = f_to_f1(read_xye('data/tot/f_damashek.dat', xmin=0.01))[:-1]
armstrong = f_to_f1(read_xye('data/tot/f_armstrong.dat', xmin=xmin))[:-1]
ahrens = f_to_f1(read_xye('out/tot/f_ahrens.out', xmin=xmin))[:-1]

data = list(map(list, [x, y, err])) + damashek + armstrong + ahrens

plot_to_file(data, filename='out/re_f1.out')
