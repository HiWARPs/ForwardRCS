#!/usr/bin/env python

from common import *
from compton import e_Mp
from amp import *
from os.path import join

outdir = 'suppl'

kk = arange(0, 2, 0.005)

with open(join(outdir, 'cs_pol.txt'), 'w') as f:
	f.write(r' \nu   s_3/2 - s_1/2   error' + '\n')
	f.write(r'(GeV)      (mcb)       (mcb)' + '\n\n')
	d = ['%.3f   %.8f   %.8f' % ((k,) + tuple(2 * wp.cs(k))) for k in kk]
	f.write('\n'.join(d))

with open(join(outdir, 'cs_tot_1.txt'), 'w') as f:
	f.write(r' \nu  (s_3/2 + s_1/2)/2  error' + '\n')
	f.write(r'(GeV)      (mcb)         (mcb)' + '\n\n')
	d = ['%.3f   %.8f   %.8f' % ((k,) + tuple(wt1.cs(k))) for k in kk]
	f.write('\n'.join(d))

with open(join(outdir, 'cs_tot_2.txt'), 'w') as f:
	f.write(r' \nu  (s_3/2 + s_1/2)/2  error' + '\n')
	f.write(r'(GeV)      (mcb)         (mcb)' + '\n\n')
	d = ['%.3f   %.8f   %.8f' % ((k,) + tuple(wt2.cs(k))) for k in kk]
	f.write('\n'.join(d))

with open(join(outdir, 'g.txt'), 'w') as f:
	f.write(r' \nu       Re g       error      Im g       error' + '\n')
	f.write(r'(GeV)    (mcb*GeV)  (mcb*GeV)  (mcb*GeV)  (mcb*GeV)' + '\n\n')
	d = ['%.3f   %.8f   %.8f   %.8f   %.8f' % (0, 0, 0, 0, 0)]
	d += ['%.3f   %.8f   %.8f   %.8f   %.8f' % ((k,) + tuple(wp.re_g(k)) + tuple(wp.im_g(k))) for k in kk[1:]]
	f.write('\n'.join(d))

with open(join(outdir, 'f_1.txt'), 'w') as f:
	f.write(r' \nu       Re f       error      Im f       error' + '\n')
	f.write(r'(GeV)    (mcb*GeV)  (mcb*GeV)  (mcb*GeV)  (mcb*GeV)' + '\n\n')
	d = ['%.3f   %.8f   %.8f   %.8f   %.8f' % (0, e_Mp, 0, 0, 0)]
	d += ['%.3f   %.8f   %.8f   %.8f   %.8f' % ((k,) + tuple(wt1.re_f(k)) + tuple(wt1.im_f(k))) for k in kk[1:]]
	f.write('\n'.join(d))

with open(join(outdir, 'f_2.txt'), 'w') as f:
	f.write(r' \nu       Re f       error      Im f       error' + '\n')
	f.write(r'(GeV)    (mcb*GeV)  (mcb*GeV)  (mcb*GeV)  (mcb*GeV)' + '\n\n')
	d = ['%.3f   %.8f   %.8f   %.8f   %.8f' % (0, e_Mp, 0, 0, 0)]
	d += ['%.3f   %.8f   %.8f   %.8f   %.8f' % ((k,) + tuple(wt2.re_f(k)) + tuple(wt2.im_f(k))) for k in kk[1:]]
	f.write('\n'.join(d))

with open(join(outdir, 'ds_1.txt'), 'w') as f:
	f.write(r' \nu    dsi/dO   error' + '\n')
	f.write(r'(GeV)    (nb)     (nb)' + '\n\n')
	d = ['%.3f   %.8f   %.8f' % (0, e_Mp**2/hc**2/10, 0)]
	d += ['%.3f   %.8f   %.8f' % ((k,) + tuple(ds(k, wt=wt1))) for k in kk[1:]]
	f.write('\n'.join(d))

with open(join(outdir, 'ds_2.txt'), 'w') as f:
	f.write(r' \nu    dsi/dO   error' + '\n')
	f.write(r'(GeV)    (nb)     (nb)' + '\n\n')
	d = ['%.3f   %.8f   %.8f' % (0, e_Mp**2/hc**2/10, 0)]
	d += ['%.3f   %.8f   %.8f' % ((k,) + tuple(ds(k, wt=wt2))) for k in kk[1:]]
	f.write('\n'.join(d))

with open(join(outdir, 'm2refg_1.txt'), 'w') as f:
	f.write(r' \nu   -2 Re f g^*   error' + '\n')
	f.write(r'(GeV)     (nb)        (nb)' + '\n\n')
	d = ['%.3f   %.8f   %.8f' % (0, 0, 0)]
	d += ['%.3f   %.8f   %.8f' % ((k,) + tuple(-2 * refg(k, wt=wt1))) for k in kk[1:]]
	f.write('\n'.join(d))

with open(join(outdir, 'm2refg_2.txt'), 'w') as f:
	f.write(r' \nu   -2 Re f g^*   error' + '\n')
	f.write(r'(GeV)     (nb)        (nb)' + '\n\n')
	d = ['%.3f   %.8f   %.8f' % (0, 0, 0)]
	d += ['%.3f   %.8f   %.8f' % ((k,) + tuple(-2 * refg(k, wt=wt2))) for k in kk[1:]]
	f.write('\n'.join(d))
