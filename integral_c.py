#!/usr/bin/env python
import tot
from compton import *
from numpy import pi	# py3 port
from integrate import Region, Smooth

def change_regge(w, f):
	w.h_region.fit = f
	w.fitted_region.fit.f2 = f
	w.update_data()

for w, n in (tot.WholeRegion_1(quiet=True), '1'), (tot.WholeRegion_2(quiet=True), '2'):
	
	print('fit '+n+' ' + '-'*20)
	
	p = w.h_region.fit.get_p()
	g1, a1, g2, a2 = p[0], p[1]*0.5, p[2], p[3]*0.5
	s0 = Mp**2
	f = lambda k: g1*(lab2cm(k)**2 - s0)**a1 + (a1 + 1)*g1*s0*(lab2cm(k)**2 - s0)**(a1-1) + g2*(lab2cm(k)**2 - s0)**a2

	i = 200

	r1 = Region(Smooth(w.h_region.fit, lambda x: 0., w.b[2]), k_min=w.b[0], k_max=w.b[-1], filename=w.filename, init=False)
	r2 = Region(Smooth(f, lambda x: 0., w.b[2]), k_min=w.b[0], k_max=w.b[-1], filename=w.filename, init=False)
	r1.err = [1.]*len(r1.err)
	r2.err = [1.]*len(r2.err)
	#for i in [1, 2, 20, 2000]:
	#	print (r'\nu from \nu_0 to %.1f GeV :' + '\t %.4f') % (i, r1.integral_c(i)[0] - r2.integral_c(i)[0])

	print(r1.integral_c(i)[0] - r2.integral_c(i)[0])

	change_regge(w, f)
	s0 = 0
	print(w.integral_c(i)[0])
#	w.plot_regge(filename='out/regge_mine_'+n, logx=True, xmax=10**7)#, raw=False, plots=['1 2 label=mine'])
#	w.plot_c(filename='out/c_mine_'+n, logx=False, xmax=4)

	w.update_data()	# cleaning wrong cache
	s0 = Mp**2
	print(w.integral_c(i)[0])
	#w.plot_regge(filename='out/regge_mine_s0_%.2f'%s0, logx=True, xmax=10**7)#, raw=False, plots=['1 2 label=mine_s0_%.2f'%s0])

	#---------------------------
	X0, X1, X2 = 0.001461, 0.5954, 1.154
	e0, e1, e2 = 0.4372, 0.0808, -0.4525
	f = lambda k: sum(x * 4*pi**2*alpha * 10**4*hc**2 * (lab2cm(k)**2 - s0) ** e for x, e in zip([X0, X1, X2], [e0, e1, e2]))
	change_regge(w, f)

	w.update_data()
	s0 = 0
	print(w.integral_c(i)[0])
#	w.plot_regge(filename='out/regge_dl_'+n, logx=True, xmax=10**7)#, raw=False, plots=['1 2 label=D-L'])
#	w.plot_c(filename='out/c_dl_'+n, logx=False, xmax=4)

	w.update_data()
	s0 = Mp**2
	print(w.integral_c(i)[0])
	#w.plot_regge(filename='out/regge_dl_s0_%.2f'%s0, logx=True, xmax=10**7)#, raw=False, plots=['1 2 label=D-L_s0_%.2f'%s0])
