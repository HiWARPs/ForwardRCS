#!/usr/bin/env python

import matplotlib.pyplot as plt
from matplotlib import rc

from common import *

# using LaTeX rendering:
rc('text', usetex=True)
rc('text.latex', preamble = ','.join(r'''
\usepackage{upgreek}
\usepackage{mathptmx}
'''.split()))
# ---------------------
rc('font', family='serif')
rc('mathtext', fontset='custom', rm='serif')

fig, (plt1, plt2) = plt.subplots(2)

#fig = plt.gcf()
fig.patch.set_facecolor('white')

#plt1.set_xlabel(r'$\upnu$ [GeV]')
#plt1.set_ylabel(r'$(\mathrm{d}\sigma/\mathrm{d}\Omega_{\mathrm{lab}})\,(\theta=0) \;\; [\mathrm{nb}]$')
fig.text(0.5, 0.02, r'$\upnu$ [GeV]', ha='center')
fig.text(0.02, 0.5, r'$\Sigma_{2z}\,(\theta=0)$', va='center', rotation='vertical')

d1 = read_xy('out/amp/s2z_1.out')
d2 = read_xy('out/amp/s2z_2.out')
x, y1, y2 = d1[0], d1[1], d2[1]
plt1.axis([0, 2.1, -1.2, 1.0])
plt1.fill_between(x, y1, y2, alpha=0.65, color='black')

d1 = read_xy('out/amp/s2z_1_small.out')
d2 = read_xy('out/amp/s2z_2_small.out')
d3 = read_xy('out/cpt/s2z')
x, y1, y2 = d1[0], d1[1], d2[1]
#plt2.set_yscale('log')
plt2.axis([0, 0.4, -1.2, 1.0])
plt2.fill_between(x, y1, y2, alpha=0.65, color='black', label='dispertion relation')
p = plt2.plot(d3[0], d3[1], color='red', label='$\chi$PT')
p[0].set_dashes([6,3])
#plt2.legend(loc=0, numpoints=1)

fig.subplots_adjust(hspace=0.2)
#fig.tight_layout()

#plt.show()
fig.set_size_inches(6,5)
fig.savefig('fig/s2z.svg')
