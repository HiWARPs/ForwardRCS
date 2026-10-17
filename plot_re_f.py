#!/usr/bin/env python

import matplotlib.pyplot as plt
from matplotlib import rc
from numpy import array	# py3 port

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

fig.text(0.50, 0.01, r'$W$ [GeV]', ha='center', fontsize=14)
# fig.text(0.01, 0.5, r'$\mathrm{Re}f$ [$\mu$b$\cdot$GeV]', va='center', rotation='vertical', fontsize=14)
# plt.xlabel(r'$W$ [GeV]', fontsize=14)
# # Haaaaaack...
# plt1.set_xlabel(r'$W$', fontsize=14, x=0.97)
# fig.text(0.513, 0.053, r'[GeV]', fontsize=14)
plt1.set_ylabel(r'$\mathrm{Re}f$ [$\mu$b$\cdot$GeV]', fontsize=14)

lns = []

path = 'paper_llpg/'
d1 = list(map(array, read_xye(path+'re_f_low_1.dat')))	# py3 port: list(map(...)) here and below
d2 = list(map(array, read_xye(path+'re_f_low_2.dat')))
plt1.axis([1.08, 2.4, -16, 5])
p = plt1.plot(d1[0], d1[1], color='red', label='based on fit I')
lns += p
plt1.fill_between(d1[0], d1[1]-d1[2], d1[1]+d1[2], alpha=0.4, color='red', linewidth=0)
p = plt1.plot(d2[0], d2[1], color='blue', label='based on fit II')
lns += p
p[0].set_dashes([6,3])
plt1.fill_between(d1[0], d2[1]-d2[2], d2[1]+d2[2], alpha=0.4, color='blue', linewidth=0)
p = plt1.errorbar([flab2cm(2.2)], [-12.3], yerr=[2.4], color='black', marker='s', linestyle='None', label=r'Alvensleben $et$ $al.$')
lns += [p]
# print p.get_label()
# exit()
# print list(p)
# exit()
# lns += p[2]
# plt1.legend(loc=3, numpoints=1)

d1 = list(map(array, read_xye(path+'re_f_high_bh.dat')))
d2 = list(map(array, read_xye(path+'re_f_high_dl.dat')))
plt2.axis([2.4, 15, -20, 20])
plt2.yaxis.tick_right()
plt2.yaxis.set_ticks_position('both')
p = plt2.plot((flab2cm(11), flab2cm(11)), (-20, 20), linewidth=1.5, color='grey')
p[0].set_dashes([2,2])
p = plt2.plot(d1[0], d1[1], color='black', label='Block - Halzen')
lns += p
plt2.fill_between(d1[0], d1[1]-d1[2], d1[1]+d1[2], alpha=0.4, color='black', linewidth=0)
p[0].set_dashes([6,3])
p = plt2.plot(d2[0], d2[1], color='orange', label='Donnachie - Landshoff')
lns += p
plt2.fill_between(d1[0], d2[1]-d2[2], d2[1]+d2[2], alpha=0.4, color='orange', linewidth=0)
plt2.set_xscale('log')
# plt2.legend(loc=1, numpoints=1)
plt.legend(lns, [l.get_label() for l in lns], loc='upper center', bbox_to_anchor=(0.0, 1.), ncol=1, numpoints=1)

# plt2.set_xticks(list(plt2.get_xticks()) + [flab2cm(11)])
# plt2.axis([2.4, 15, -20, 20])
# ll = plt2.get_xticks().tolist()
# ll[-1] = 'E=11.0'
# plt2.set_xticklabels(ll)


# fig.subplots_adjust(wspace=0.15)
fig.subplots_adjust(wspace=0.0)
# fig.subplots_adjust(wspace=0.0)
# fig.tight_layout()	# the same as bbox_inches='tight'
# plt.tight_layout(rect = [0.024, 0.036, 1, 1])
# plt.tight_layout(rect = [0.03, 0.05, 0.96, 1])
# plt.tight_layout(rect = [0.06, 0.06, 0.98, 0.95])
fig.set_tight_layout(True)

#plt.show()
# fig.set_size_inches(7.3,4.9)
# fig.set_size_inches(7.0,4.2)
# fig.set_size_inches(5.1,3.8)
# fig.set_size_inches(5.8,4.2)
fig.set_size_inches(5.8,3.8)

# fig.savefig('paper_llpg/re_f.pdf', bbox_inches='tight')
fig.savefig('paper_llpg/re_f.pdf')
