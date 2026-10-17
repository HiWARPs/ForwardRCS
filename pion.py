#!/usr/bin/env python

from compton import *
from common import *
from fit import *
from integrate import *

t = Region(fit=BornT(), k_min=fcm2lab(Mp+Mpi0), init=False)
p = Region(fit=BornP(), k_min=fcm2lab(Mp+Mpi0), init=False)

def s2z(k):
	fk = t.f(k)
	#gk = p.re_g0(k) + 1j*p.im_g(k)
	gk = p.g(k)
	return -2*(fk*gk.conjugate()).real / (fk.real**2 + fk.imag**2 + gk.real**2 + gk.imag**2)

# x = arange(0.1, 2.1, 0.1)
# thr = fcm2lab(Mp+Mpi0)
# y = [s2z(thr-0.001)-s2z(thr-0.002) for compton.gA in map(lambda c: gA*c**0.5, x)]
# plot_to_file([x,y],['1 2'],filename='tmp/s2ztest')

# exit()

c = 1
compton.gA = gA*(c)**0.5
d = 0.002
x = arange(d, 1, d)

plot_to_file([x, list(map(t.im_f, x))],['1 2 label="Im f"'],filename='tmp/pi_imf')	# py3 port: list(map(...)) here and below
plot_to_file([x, list(map(p.im_g, x))],['1 2 label="Im g"'],filename='tmp/pi_img')
plot_to_file([x, list(map(t.re_f, x))],['1 2 label="Re f"'],filename='tmp/pi_ref')
plot_to_file([x, list(map(p.re_g, x))],['1 2 label="Re g"'],filename='tmp/pi_reg')

plot_to_file([x, list(map(s2z, x))],[r'1 2 label=$\\Sigma_{2z}$'],filename='tmp/pi_s2z')