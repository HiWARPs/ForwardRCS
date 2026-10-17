#!/usr/bin/env python
from custom import *
from common import *

# # quad = lambda f, k0, k1, **kw: mpmath.quad(f, [k0, k1], **dict(kw, error=True))
quad = lambda f, k0, k1, **kw: mpmath.quad(lambda k: f(float(k)), [k0, k1], **dict(kw, error=True))

def singular(k, k0, k1):
	f = lambda p: 0. if p is inf else log(abs((p-k)/(p+k)))/(2*k)
	return f(k1) - f(k0)

def re_f(k, w=None):
	cs = lambda k: float(w.cs(k, err=False))
	cs_k = cs(k)
	k0, k1 = w.k_min, w.k_max
	r = val(quad(lambda p: (cs(p) - cs_k)/(p**2 - k**2), k0, k1)) + cs_k*singular(k, k0, k1)
	return e_Mp + k**2*r/(2*pi**2)

def im_f(k, w=None):
	return w.im_f(k, err=False)

def re_over_im_f(k, w=None):
	return re_f(k, w)/im_f(k, w)

w1 = BHWholeRegion()
w1.smooth.set_p(0.02)
w2 = DLWholeRegion()
print "plotting..."

x = arange(flab2cm(w1.k_min+2e-3), 1e5, dx=1e-2, logx=True)
# x = arange(1.08, 1e8, dx=0.01, logx=True)
# x = [10]
y_e = zip(*(sum((tuple(w.re_over_im_f(k=fcm2lab(W), err=True)) for w in (w1, w2)), ()) for W in x))
# y_e = zip(*(sum((tuple(re_over_im_f(k=fcm2lab(W), w=w)) for w in (w1, w2)), ()) for W in x))
# print y_e
plot_to_file([x] + y_e, filename='paper_llpg/re_im.dat')
