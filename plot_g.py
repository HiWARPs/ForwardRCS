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
fig.text(0.02, 0.5, r'$g \;\; [\upmu \mathrm{b} \cdot \mathrm{GeV}]$', va='center', rotation='vertical')

plt1.axis([0, 2.1, -7.2, 4.5])
x, y1, e1 = read_xye('out/pol/g.out')
_, y2, e2 = read_xye('out/pol/im_g.out')

plt1.fill_between(x, [ay-ae for ay,ae in zip(y1,e1)], [ay+ae for ay,ae in zip(y1,e1)], alpha=0.36, linewidth=0, color='purple')
plt1.plot(x, y1, color='purple', label='$\mathrm{Re}\,g$')
plt1.fill_between(x, [ay-ae for ay,ae in zip(y2,e2)], [ay+ae for ay,ae in zip(y2,e2)], alpha=0.36, linewidth=0, color='blue')
p = plt1.plot(x, y2, color='blue', label='$\mathrm{Im}\,g$')
p[0].set_dashes([5,1.5,1.5,1.5])
xd, yd, ed = read_xye('data/pol/data_img.dat')
plt1.errorbar(xd, yd, yerr=ed, label="MAMI + ELSA", marker='o', markersize=5, color='black', linestyle='')

plt2.axis([0, 0.4, -7, 3])
#plt2.set_xscale('log')
#x, y1, y2 = read_xye('out/amp/g.out')
plt2.fill_between(x, [ay-ae for ay,ae in zip(y1,e1)], [ay+ae for ay,ae in zip(y1,e1)], alpha=0.36, linewidth=0, color='purple')
plt2.plot(x, y1, color='purple')#, label='$\mathrm{Re}\,g$')
plt2.fill_between(x, [ay-ae for ay,ae in zip(y2,e2)], [ay+ae for ay,ae in zip(y2,e2)], alpha=0.36, linewidth=0, color='blue')
p = plt2.plot(x, y2, color='blue')#, label='$\mathrm{Im}\,g$')
p[0].set_dashes([5,1.5,1.5,1.5])
x, y1, y2 = read_xye('out/cpt/g')
p = plt2.plot(x, y1, color='red', label='$\chi$PT $\mathrm{Re}\,g$')
p[0].set_dashes([5.5,2])
p = plt2.plot(x, y2, color='green', label='$\chi$PT $\mathrm{Im}\,g$')
p[0].set_dashes([1.5,1.5])

plt1.legend(loc=0, numpoints=1)
plt2.legend(loc=0, numpoints=1)

fig.subplots_adjust(hspace=0.2)
#fig.tight_layout()

#plt.show()
fig.set_size_inches(6,6)
fig.savefig('fig/g.svg')
fig.savefig('fig/g.pdf')
