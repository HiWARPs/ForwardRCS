#!/usr/bin/env python

from numpy import pi	# py3 port
from compton import *
from common import *
from tot import *
from custom import *

x = arange(0.1, 1000, 0.01, logx=True)
# w1, w2 = BHWholeRegion(), DLWholeRegion()
w = BHWholeRegion()

r, n = 0, 1
for y, e in zip(w.y, w.err):
	if y != 0 and e/y > 0.001:
		r += e/y
		n += 1

r /= n
print('average cs relative error:', r)

plot_to_file([x, [(1+r)*w.cs(k) for k in x]], ['1 2'], filename='tmp/cs_p', xlogscale=True)
plot_to_file([x, [(1-r)*w.cs(k) for k in x]], ['1 2'], filename='tmp/cs_m', xlogscale=True)

# plot_to_file([x, map(lambda a: float(w1.integral_f(6, w1.k_min, fcm2lab(a))/w2.integral_f(6, w1.k_min, fcm2lab(a))), x)], ['1 2'], 'tmp/dfcutoff_6', xlogscale=True)
# print 'amp. relative error:', 
print(w.integral_f(6))

p, perr = w.r_region.fit.get_p() + w.h_region.fit.get_p(), w.r_region.fit.get_perr() + w.h_region.fit.get_perr()
pp = p
for i in range(len(p)):
	pp_ = []
	for a in pp:
		a1 = list(p)
		a1[i] += perr[i]
		a2 = list(p)
		a2[i] -= perr[i]
		pp_ += [a1, a2]
	pp = pp_

print(len(p), len(pp))

# for a in pp:
	# w._cache = {}
	# w.fit.set_p()