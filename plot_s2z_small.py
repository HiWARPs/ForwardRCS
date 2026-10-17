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

#fig, (plt1, plt2) = plt.subplots(2)
fig = plt.gcf()

fig.patch.set_facecolor('white')

#plt1.set_xlabel(r'$\upnu$ [GeV]')
#plt1.set_ylabel(r'$(\mathrm{d}\sigma/\mathrm{d}\Omega_{\mathrm{lab}})\,(\theta=0) \;\; [\mathrm{nb}]$')
fig.text(0.5, 0.02, r'$\upnu$ [GeV]', ha='center')
fig.text(0.02, 0.5, r'$\Sigma_{2z}\,(\theta=0)$', va='center', rotation='vertical')

d1 = read_xy('out/amp/s2z_1_small.out')
d2 = read_xy('out/amp/s2z_2_small.out')
x, y1, y2 = d1[0], d1[1], d2[1]
plt.axis([0.12, 0.18, -1.2, 0.])
plt.fill_between(x, y1, y2, alpha=0.65, color='black')

#plt.show()
fig.set_size_inches(3,2)
fig.savefig('fig/s2z_small.svg')
