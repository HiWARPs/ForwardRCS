#!/usr/bin/env python

from common import *
from tot import *

w = WholeRegion_1(quiet=True)
dx = 0.001
x = arange(dx, 10, dx)
plot_to_file([x,map(w.cst,x)],['1 2'],'tmp/cst')
plot_to_file([x,map(w.cst1,x)],['1 2'],'tmp/cst1')

dx = 0.01
x = arange(dx, 1000, dx, logx=True)
plot_to_file([x,map(w.cst1,x)],['1 2'],'tmp/cst1_logx',xlogscale=True)