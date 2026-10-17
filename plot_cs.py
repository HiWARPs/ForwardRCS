#!/usr/bin/env python

import matplotlib.pyplot as plt
from matplotlib import rc
from scipy import array

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

# fig, (plt1, plt2) = plt.subplots(2)
plt1 = plt.subplot2grid((1,2), (0, 0))
plt2 = plt.subplot2grid((1,2), (0, 1))
fig = plt.gcf()

fig.patch.set_facecolor('white')

fig.text(0.5, 0.01, r'$W$ [GeV]', ha='center', fontsize=14)
# # Haaaaaack...
# plt1.set_xlabel(r'$W$', fontsize=14, x=0.97)
# fig.text(0.520, 0.053, r'[GeV]', ha='center', fontsize=14)
# fig.text(0.01, 0.5, r'$\sigma$ [$\mu$b]', va='center', rotation='vertical', fontsize=14)
# plt1.set_xlabel(r'$W$ [GeV]', fontsize=14, x=1)
plt1.set_ylabel(r'$\sigma$ [$\mu$b]', fontsize=14)

data_path = 'data/wtot/'
datasets = [('armstrong.dat', r'Armstrong $et\,al.$', 'o', 4, (152/255., 223/255., 138/255.))]
			# ('maccormick.dat', r'MacCormick $et$ $al.$', '.', 6, 'black'),
			# ('bartalini.dat', r'Bartalini $et$ $al.$', '.', 6, 'green')]
for s,l,m,ms,c in datasets:
	x, y, e = map(array, read_xye(data_path+s))
	plt1.errorbar(x, y, yerr=e, color=c, marker=m, markersize=ms, linestyle='None', label=l)

path = 'paper_llpg/'
d1 = map(array, read_xye(path+'cs_1.dat'))
d2 = map(array, read_xye(path+'cs_2.dat'))
plt1.axis([1.08, 2.4, 0, 600])
p = plt1.plot(d1[0], d1[1], color='red', label='fit I')
plt1.fill_between(d1[0], d1[1]-d1[2], d1[1]+d1[2], alpha=0.4, color='red', linewidth=0)
p = plt1.plot(d2[0], d2[1], color='blue', label='fit II')
p[0].set_dashes([6,3])
plt1.fill_between(d1[0], d2[1]-d2[2], d2[1]+d2[2], alpha=0.4, color='blue', linewidth=0)
plt1.legend(loc=1, numpoints=1)

datasets = [('armstrong.dat', None, 'o', 4, (152/255., 223/255., 138/255.)),
			('pdg/caldwell.dat', r'Caldwell $et\,al.$', 'v', 6, 'cyan'),
			('pdg/aid.dat', 'H1 Collaboration', '>', 6, 'red'),
			('pdg/chekanov.dat', 'ZEUS Collaboration', '<', 6, 'purple')]
for s,l,m,ms,c in datasets:
	x, y, e = map(array, read_xye(data_path+s))
	plt2.errorbar(x, y, yerr=e, color=c, marker=m, markersize=ms, linestyle='None', label=l)

d1 = map(array, read_xye(path+'cs_bh.dat'))
d2 = map(array, read_xye(path+'cs_dl.dat'))
plt2.axis([2.4, 250, 100, 200])
plt2.yaxis.tick_right()
plt2.yaxis.set_ticks_position('both')
p = plt2.plot(d1[0], d1[1], color='black', label='Block - Halzen fit')
plt2.fill_between(d1[0], d1[1]-d1[2], d1[1]+d1[2], alpha=0.4, color='black', linewidth=0)
p[0].set_dashes([6,3])
p = plt2.plot(d2[0], d2[1], color='orange', label='Donnachie - Landshoff fit')
plt2.fill_between(d1[0], d2[1]-d2[2], d2[1]+d2[2], alpha=0.4, color='orange', linewidth=0)
plt2.set_xscale('log')
plt2.legend(loc=2, bbox_to_anchor=(0.08, 1), numpoints=1)

# fig.subplots_adjust(wspace=0.15)
# fig.subplots_adjust(wspace=0.0)
# fig.tight_layout()	# the same as bbox_inches='tight'
# plt.tight_layout(rect = [-0.01, 0.032, 1.01, 1])
# plt.tight_layout(rect = [0.01, 0.06, 1.01, 1])
fig.set_tight_layout(True)
# fig.subplots_adjust(wspace=0.04)
fig.subplots_adjust(wspace=0.0)

#plt.show()
# fig.set_size_inches(7,3.8)
# fig.set_size_inches(7.3,4.1)
# fig.set_size_inches(14.6,4.9)
# fig.set_size_inches(10.2,3.8)
# fig.set_size_inches(11.6,4.2)
fig.set_size_inches(11.6,3.8)
# fig.savefig('paper_llpg/cs.pdf', bbox_inches='tight')
fig.savefig('paper_llpg/cs.pdf')
