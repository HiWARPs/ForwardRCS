#!/usr/bin/env python

import matplotlib.pyplot as plt
from matplotlib import rc
from numpy import array, array_split	# py3 port

from common import *
from compton import *

# using LaTeX rendering:
rc('text', usetex=True)
rc('text.latex', preamble = ','.join(r'''
\usepackage{upgreek}
\usepackage{mathptmx}
'''.split()))
# ---------------------
rc('font', family='serif')
rc('mathtext', fontset='custom', rm='serif')

fig = plt.gcf()

fig.patch.set_facecolor('white')

fig.text(0.45, 0.01, r'$\Lambda$ [GeV]', ha='center', fontsize=14)
# fig.text(0.01, 0.5, r'$\mathrm{Re}f_\Lambda\;/\;\mathrm{Re}f$', va='center', rotation='vertical', fontsize=14)
# plt.xlabel(r'$\Lambda$ [GeV]', fontsize=14)
plt.ylabel(r'$\mathrm{Re}f_\Lambda\;/\;\mathrm{Re}f$', fontsize=14)

path = 'paper_llpg/'
d1 = list(map(array, read_xye(path+'f_cutoffn_bh.dat')))	# py3 port: list(map(...)) here and below
d2 = list(map(array, read_xye(path+'f_cutoffn_dl.dat')))
plt.axis([5, 25, 0.9, 1.05])
p = plt.plot(d1[0], d1[1], color='black', label='Block - Halzen $\,$fit')
p[0].set_dashes([6,3])
plt.fill_between(d1[0], d1[1]-d1[2], d1[1]+d1[2], alpha=0.4, color='black', linewidth=0)
p = plt.plot(d2[0], d2[1], color='orange', label='Donnachie - Landshoff $\,$fit')
n, x0, x1 = 10, 8, 10
xchunks = array_split([x for x in d2[0] if x0 <= x <= x1], n)
xmax = min(xchunks[0])
for i, xchunk in zip(range(1,n+1), xchunks):
	if len(xchunk):
		xmin, xmax = xmax, max(xchunk)
		if xmax != xmin:
			d = list(map(array, zip(*[l for l in zip(*d2) if xmin <= l[0] <= xmax])))
			plt.fill_between(d[0], d[1]-d[2], d[1]+d[2], alpha=0.4*i/(n+1), color='orange', linewidth=0)
d = list(map(array, zip(*[l for l in zip(*d2) if l[0] >= x1])))
plt.fill_between(d[0], d[1]-d[2], d[1]+d[2], alpha=0.4, color='orange', linewidth=0)
plt.legend(loc=4, numpoints=1)

# fig.tight_layout(rect = [0.024, 0.036, 0.96, 1])
# fig.tight_layout(rect = [0.07, 0.06, 0.96, 0.95])
fig.set_tight_layout(True)
# fig.set_size_inches(7.3,4.9)
# fig.set_size_inches(5.1,3.8)
# fig.set_size_inches(5.8,4.2)
fig.set_size_inches(5.8,3.8)
fig.savefig('paper_llpg/cutoff.pdf')
