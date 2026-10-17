#!/usr/bin/env python2

from tot import *
from common import *

w = WholeRegion_2(init=True)

xdata = arange(0.002, 2, 0.002)
resdata = [[r(x) for x in xdata] for r in w.r_region.fit.resonances]
fdata = [float(w.cs(x)) for x in xdata]
edata = [w.cs(x).err for x in xdata]

data = [xdata, fdata, edata] + resdata

plots = ['1 2 yband=3 label=fit linewidth=1 color=black'] + ['1 '+str(i+4) for i in range(len(resdata))]
plot_to_file(data, plots, 'outfit')