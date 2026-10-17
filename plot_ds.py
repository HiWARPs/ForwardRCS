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
fig.text(0.02, 0.5, r'$(\mathrm{d}\sigma/\mathrm{d}\Omega_{\mathrm{lab}})\,(\theta=0) \;\; [\mathrm{nb}]$', va='center', rotation='vertical')

d1 = read_xye('out/amp/ds_1.out')
d2 = read_xye('out/amp/ds_2.out')

x, y1, e1, y2, e2 = d1[0], d1[1], d1[2], d2[1], d2[2]
plt1.axis([0, 2.1, 0, 1700])
plt1.fill_between(x, [ay-ae for ay,ae in zip(y1,e1)], [ay+ae for ay,ae in zip(y1,e1)], alpha=1, linewidth=0, color='0.68')
plt1.fill_between(x, [ay-ae for ay,ae in zip(y2,e2)], [ay+ae for ay,ae in zip(y2,e2)], alpha=1, linewidth=0, color='0.68')
plt1.plot(x, y1, color='black')
p = plt1.plot(x, y2, color='black')
p[0].set_dashes([1,1])

d1 = read_xye('out/amp/ds_1_small.out')
d2 = read_xye('out/amp/ds_2_small.out')
d3 = read_xy('out/cpt/ds')

x, y1, e1, y2, e2 = d1[0], d1[1], d1[2], d2[1], d2[2]
plt2.set_yscale('log')
plt2.axis([0, 0.4, 0.1, 2000])
plt2.fill_between(x, [ay-ae for ay,ae in zip(y1,e1)], [ay+ae for ay,ae in zip(y1,e1)], alpha=1, linewidth=0, color='0.68')
plt2.fill_between(x, [ay-ae for ay,ae in zip(y2,e2)], [ay+ae for ay,ae in zip(y2,e2)], alpha=1, linewidth=0, color='0.68')
plt2.plot(x, y1, color='black', label=r'fit I of $\sigma_\mathrm{abs}$')
p = plt2.plot(x, y2, color='black', label=r'fit II of $\sigma_\mathrm{abs}$', linestyle=':')
p[0].set_dashes([1,1])
#plt2.fill_between(x, y1, y2, alpha=0.65, color='black', label='dispertion relation')
p = plt2.plot(d3[0], d3[1], color='red', label='$\chi$PT')
p[0].set_dashes([6,3])

plt2.legend(loc=0, numpoints=1)

fig.subplots_adjust(hspace=0.2)
#fig.tight_layout()

#plt.show()
fig.set_size_inches(6,5)
fig.savefig('fig/ds.svg')
fig.savefig('fig/ds.pdf')
