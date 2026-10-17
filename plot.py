#!/usr/bin/env python3

import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.ticker import MultipleLocator
import argparse
import shlex
import os


class Plot:
	
	def __init__(self, iy, ix=None, iex=None, iey=None, iby=None, color=None, mrkr=None, mrkrsz=None, lnwd=None, lnstyle=None, label=None, fiy=lambda y: y):
		self.kwargs = {'linewidth':lnwd or 1, 'linestyle':lnstyle or '-', 'marker':mrkr, 'label':(label.replace('\eq', '=') if label is not None else None)}
		if mrkrsz is not None:
			self.kwargs['markersize'] = mrkrsz
		if color is not None:
			self.kwargs['color'] = color
		self.ix, self.iy, self.iex, self.iey, self.iby = ix, iy, iex, iey, iby
		self.fiy = fiy
	
	def shift(self, shift):
		self.ix, self.iy, self.iex, self.iey, self.iby = map(lambda i: i if i is None else i+shift, (self.ix, self.iy, self.iex, self.iey, self.iby))
	
	def plot(self, data, plt=plt):
		y0 = [d and self.fiy(d) for d in data[self.iy-1]]
		if self.ix is None:
			x0 = range(len(y0))
		else:
			x0 = data[self.ix-1]
		if self.iex is None:
			ex0 = [None]*len(x0)
		else:
			ex0 = map(lambda d: 0. if d is None else d, data[self.iex-1])
		if self.iey is None:
			ey0 = [None]*len(y0)
		else:
			ey0 = map(lambda d: 0. if d is None else self.fiy(d), data[self.iey-1])
		if self.iby is None:
			by0 = [None]*len(y0)
		else:
			by0 = map(lambda d: 0. if d is None else self.fiy(d), data[self.iby-1])
		x, y, ex, ey, by = zip(*[[dx, dy, dex, dey, dby] for dx, dy, dex, dey, dby in zip(x0,y0,ex0,ey0,by0) if dx is not None and dy is not None])
		if self.iex is not None:
			self.kwargs['xerr'] = ex
		if self.iey is not None:
			self.kwargs['yerr'] = ey
		kwargs = {k:v for k,v in self.kwargs.items() if not (k == 'linestyle' and v[0] == '[')}	# py3 port
		self.p = plt.errorbar(x, y, **kwargs)
		if self.iby:
			plt.fill_between(x, [ay - ab for ay, ab in zip(y, by)], [ay + ab for ay, ab in zip(y, by)], alpha=0.4, linewidth=0, facecolor=kwargs.get('color'))
		if self.kwargs['linestyle'][0] == '[':
			self.p[0].set_dashes(list(map(float, self.kwargs['linestyle'][1:-1].split(','))))
		if self.kwargs['label'] is not None:
			plt.legend(numpoints=1)


def isNumber(a):
	try:
		float(a)
	except ValueError:
		return False
	else:
		return True


_default_axis = [None]*4
_opened_files = set()

def parseFile(filename):
	
	global _default_axis
	default_style = {'color':None, 'mrkr':None, 'mrkrsz':None, 'lnstyle':None, 'lnwd':None}
	included_files, data, inc_data, plots, inc_plots = [], [], [], [], []
	
	global _opened_files
	_opened_files.add(filename)
	
	with open(filename) as inputfile:
		for eachline in inputfile:
			splitted = eachline.split(':')
			a, val = splitted[0].strip().lower(), ':'.join(splitted[1:]).strip()
			if a == '' or a[0] == '#':
				pass
			elif a == 'title':
				plt.title(val)
			elif a == 'xlogscale':
				plt.xscale('log')
			elif a == 'ylogscale':
				plt.yscale('log')
			elif a == 'xlabel':
				# plt.xlabel(shlex.split(val)[0], fontsize=14)
				plt.gcf().text(0.50, 0.01, xlabel, ha='center', fontsize=14)
			elif a == 'ylabel':
				plt.ylabel(shlex.split(val)[0], fontsize=14)
			elif a == 'xylabel':
				xlabel, ylabel = shlex.split(val)
				# plt.xlabel(xlabel, fontsize=14)
				plt.gcf().text(0.50, 0.01, xlabel, ha='center', fontsize=14)
				plt.ylabel(ylabel, fontsize=14)
			elif a == 'size':
				plt.gcf().set_size_inches(*map(float, val.split()))
			elif a == 'minx':
				_default_axis[0] = float(val)
			elif a == 'maxx':
				_default_axis[1] = float(val)
			elif a == 'miny':
				_default_axis[2] = float(val)
			elif a == 'maxy':
				_default_axis[3] = float(val)
			elif a == 'xminorticks':
				plt.axes().xaxis.set_minor_locator(MultipleLocator(float(val)))
			elif a == 'yminorticks':
				plt.axes().yaxis.set_minor_locator(MultipleLocator(float(val)))
			elif a == 'xticks':
				axis = plt.gca()
				vals = shlex.split(val)
				axis.xaxis.set_ticks(list(map(float, vals)))	# py3 port
				axis.xaxis.set_ticklabels(vals)
			elif a == 'yticks':
				axis = plt.gca()
				vals = shlex.split(val)
				axis.yaxis.set_ticks(list(map(float, vals)))	# py3 port
				axis.yaxis.set_ticklabels(vals)
			elif a == 'color':
				default_style['color'] = val
			elif a == 'marker':
				default_style['mrkr'] = val
			elif a == 'markersize':
				default_style['mrkrsz'] = val
			elif a == 'linewidth':
				default_style['lnwd'] = val
			elif a == 'linestyle':
				default_style['lnstyle'] = val
			elif a == 'plot':
				l = shlex.split(val)
				if all(map(isNumber, l[1].lower().split('x'))):
					ly = l[1].lower().split('x')
					ix = int(l[0])
					l = l[2:]
				else:
					ly = l[0].lower().split('x')
					ix = None
					l = l[1:]
				iy = int(ly[0])
				fiy = len(ly) == 2 and (lambda n: lambda y: y*n)(float(ly[1])) or (lambda y: y)
				color = default_style['color']
				mrkr = default_style['mrkr']
				mrkrsz = default_style['mrkrsz']
				lnwd = default_style['lnwd']
				lnstyle = default_style['lnstyle']
				iex = iey = iby = label = source = None
				for s in l:
					a = s.split('=')
					pn = a[0].strip()
					if 'xerr' in pn:
						iex = int(a[1])
					elif 'yerr' in pn:
						iey = int(a[1])
					elif 'yband' in pn:
						iby = int(a[1])
					elif pn == 'color':
						color = a[1].strip()
					elif pn == 'marker':
						mrkr = a[1].strip()
					elif pn == 'markersize':
						mrkrsz = a[1].strip()
					elif pn == 'linewidth':
						lnwd = float(a[1])
					elif pn == 'linestyle':
						lnstyle = a[1]
					elif pn == 'label':
						label = a[1]
					elif pn == 'source':
						source = a[1]
				shift = 0
				#label = label or source
				if source is not None:
					source = os.path.join(os.path.dirname(filename), source)
					if source in _opened_files:
						raise Exception('Recursive inclusion!')
					if source in included_files:
						shift = sum(map(len, inc_data[:included_files.index(source)]))
					else:
						nd, np = parseFile(source)
						shift = sum(map(len, inc_data))
						for p in np:
							p.shift(shift)
						inc_data.append(nd)
						inc_plots += np
						included_files.append(source)
				plots.append(Plot(iy,ix,iex,iey,iby,color,mrkr,mrkrsz,lnwd,lnstyle,label,fiy))
				plots[-1].shift(shift)
			elif a == 'include':
				for f in shlex.split(val):
					f = os.path.join(os.path.dirname(filename), f)
					if f in _opened_files:
						raise Exception('Recursive inclusion!')
					if f not in included_files:
						nd, np = parseFile(f)
						shift = sum(map(len, inc_data))
						for p in np:
							p.shift(shift)
						inc_data.append(nd)
						inc_plots += np
						included_files.append(f)
			else:
				data.append(a.split())
	
	_opened_files.remove(filename)
	
	if data != []:
		dlen = max(map(len, data))
		data = list(zip(*list(map(lambda s: list(map(lambda a: float(a) if isNumber(a) else None, s)) + [None]*(dlen-len(s)), data))))
	
	shift = len(data)
	for p in inc_plots:
		p.shift(shift)
	
	return data + sum(inc_data,[]), plots + inc_plots


def parseArgs():
	parser = argparse.ArgumentParser()
	parser.add_argument('--subplots', action='store_true')
	parser.add_argument('-o')
	parser.add_argument('-s', '--size')
	parser.add_argument('filenames', nargs='+')
	return parser.parse_args()


if __name__ == '__main__':
	
	args = parseArgs()
	
	filenames = args.filenames or ['plot.txt']
	
	data, plots = [], []
	for f in filenames:
		ndata, nplots = parseFile(f)
		for plot in nplots:
			plot.shift(len(data))
		data += ndata
		plots += nplots
	
	if plots == []:
		plots = [Plot(i+1) for i in range(len(data))]
	
	x1, x2, y1, y2 = _default_axis
	all_x = sum([[d for d in data[plot.ix-1] if d is not None] if plot.ix is not None else list(range(len(data[plot.iy-1]))) for plot in plots],[])	# py3 port: list(range)
	all_y = sum([[d for d in data[plot.iy-1] if d is not None] for plot in plots],[])
	x2 = x2 or max(all_x)
	x1 = x1 or min(all_x)
	y2 = y2 or max(all_y)
	y1 = y1 or min(all_y)
	xs = x2 - x1
	ys = y2 - y1
	#plt.axis([x1-xs/20., x2+xs/20., y1-ys/20., y2+ys/20.])
	plt.axis([x1, x2, y1, y2])
	
	if args.subplots:
		fig, subplots = plt.subplots(len(plots))
		try:
			subplots[0]
		except IndexError:
			subplots = (subplots,)
		for plot, subplot in zip(plots, subplots):
			plot.plot(data, subplot)
		fig.subplots_adjust(hspace=0)
		plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
	else:
		fig = plt.gcf()
		for plot in plots:
			plot.plot(data)
	
	# using LaTeX rendering:
	rc('text', usetex=True)
	rc('text.latex', preamble = ','.join(r'''
	\usepackage{upgreek}
	\usepackage{mathptmx}
	'''.split()))
	# ---------------------
	rc('font', family='serif')
	#rc('text', dvipnghack=False, usetex=True)
	rc('mathtext', fontset='custom', rm='serif')
	fig.patch.set_facecolor('white')
	#plt.grid(True, which="both", ls="-", color='gray')
	#plt.minorticks_on()
	#plt.gca().tick_params(which='minor', bottom='on', labelbottom='on')
	#plt.setp(plt.gca().get_xticklabels(minor=True), visible=True)
	fig.set_size_inches(10,6)
	plt.legend(loc=0, numpoints=1, fontsize=14)
	
	if args.o:
		if args.size:
			fig.set_size_inches(*map(float, args.size.split('x')))
		fig.set_tight_layout(True)
		fig.savefig(args.o)
	else:
		plt.show()
