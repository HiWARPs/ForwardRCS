#!/usr/bin/env python

from compton import *
import integrate

class GDHDataRegion(integrate.ThresholdRegion, integrate.WholeRegion):
	
	cs = integrate.Region.cs	# py3 port: keep the Python 2 old-style (depth-first) resolution ThresholdRegion -> Region before WholeRegion
	
	def __init__(self, k_min=0.2026, k_max=2.820, name='data', path='pol', filename='data1.dat'):
		integrate.ThresholdRegion.__init__(self, k_min=k_min, k_max=k_max, name=name, path=path, filename=filename, pre_f=False)
		self.b = [self.k_min, E1, E2, self.k_max]

r = GDHDataRegion()
