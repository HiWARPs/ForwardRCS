def str_v(v_err, d=4, td=8):	# py3 port: no tuple parameter unpacking in py3
	v, err = v_err
	d, td = map(str, (d, td))
	return ('%'+td+'.'+d+r'f \pm %'+td+'.'+d+'f') % (v, err)


from functools import cmp_to_key	# py3 port


def read_xy(filename, xmin = None, xmax = None):
	data = []
	with open(filename) as fi:
		for l in fi:
			if l.strip()!='' and l.strip()[0]!='#':
				try:
					data.append(list(map(float, l.split()[:2])))
				except:
					pass
	compare = lambda a, b: 1 if a[0]>b[0] else -1
	data = sorted(data, key=cmp_to_key(compare))	# py3 port (same timsort, same '<' comparisons)
	r = list(zip(*[(d[0], d[1]) for d in data if (xmin or d[0]) <= d[0] <= (xmax or d[0])]))
	return [(),()] if r==[] else r


def read_xye(filename, xmin = None, xmax = None):
	data = []
	with open(filename) as fi:
		for l in fi:
			if l.strip()!='' and l.strip()[0]!='#':
				try:
					data.append(list(map(float, l.split()[:3])))
				except:
					pass
	compare = lambda a, b: 1 if a[0]>b[0] else -1
	data = sorted(data, key=cmp_to_key(compare))	# py3 port (same timsort, same '<' comparisons)
	r = list(zip(*[(d[0], d[1], 0. if len(d) < 3 else d[2]) for d in data if (xmin or d[0]) <= d[0] <= (xmax or d[0])]))
	return [(),(),()] if r==[] else r


def lfind_sorted(a, x):
	if len(a) == 1:
		return 0
	l2 = len(a)//2	# py3 port: integer division
	return lfind_sorted(a[:l2], x) if a[l2]>x else l2 + lfind_sorted(a[l2:], x)

def rfind_sorted(a, x):
	if len(a) == 1:
		return 0
	l2 = len(a)//2	# py3 port: integer division
	return rfind_sorted(a[:l2], x) if a[l2-1]>=x else l2 + rfind_sorted(a[l2:], x)

def split_sorted(a, *args):
	r = []
	for x in sorted(args):
		if len(a) == 0:
			r += [[]]
			continue
		elif x > a[-1]:
			i = len(a)
		else:
			i = rfind_sorted(a, x)
		r += [a[:i]]
		a = a[i:]
	return r + [a]


import math

def arange(a, b, dx, logx=False):
	if a > b:
		return []
	elif logx:
		return [math.exp(math.log(a) + dx*i) for i in range(int((math.log(b)-math.log(a))/dx))] + [b]
	else:
		return [a + dx*i for i in range(int((b-a)/dx))] + [b]


from scipy.interpolate import interp1d

_interp_cache = {}
def interp(x, xd, yd):
# def interp1(x, xd, yd):
	key = tuple(map(tuple, (xd, yd)))
	_interp_cache[key] = f = _interp_cache.get(key) or interp1d(key[0], key[1], kind='cubic', copy=False, bounds_error=False, fill_value=(0, yd[-1]), assume_sorted=True)
	return f(x)


def interp0(x, xd, yd):
# def interp(x, xd, yd):
	i = lfind_sorted(xd, x)
	if i == len(xd) - 1:
		i -= 1
	x1, x2, y1, y2 = xd[i], xd[i+1], yd[i], yd[i+1]
	return (x-x1)*(y2-y1)/(x2-x1) + y1


def linear(x1, x2, y1, y2):
	return lambda x: (x-x1)*(y2-y1)/(x2-x1) + y1


from itertools import zip_longest as izip_longest	# py3 port

def plot_to_file(data, plots=[], filename='out/out', title=None, xylabel=None, xlogscale=False, **kwargs):
	with open(filename, 'w') as fout:
		for k, v in kwargs.items():	# py3 port
			fout.write('%s: %s\n' % (k, v))
		if xlogscale:
			fout.write('xlogscale:\n')
		if title is not None:
			fout.write('title: ' + title + '\n')
		if xylabel is not None:
			fout.write('xylabel: ' + xylabel + '\n')
		for p in plots:
			fout.write('plot: ' + p + '\n')
		# fout.write('\n'.join([' '.join(map(str, map(float, l))) for l in izip_longest(*data)]))
		fout.write('\n'.join([' '.join(map(str, l)) for l in izip_longest(*data, fillvalue='None')]))


from sys import stdout
import numpy

def map_verbose(f, l, s=''):
	r = []
	stdout.write('\r[0%  ] ' + s)
	for i in range(len(l)):
		stdout.write('\r[%d%%' % (100*(i+1)/len(l)))
		stdout.flush()
		numpy.seterr(all='ignore')	# ignoring warnings
		r.append(f(l[i]))
		numpy.seterr(all='warn')
	stdout.write('\n')
	return r


import mpmath
# from mpmath import pi, nan, inf

# # vectorize mpmath functions
# for a in ('sin', 'cos', 'exp', 'log', 'sqrt'):
# 	vars()[a] = numpy.vectorize(getattr(mpmath, a))

# from scipy import pi, nan, inf, sin, cos, exp, log, sqrt, tan
from numpy import pi, nan, inf, sin, cos, exp, tan	# py3 port: these scipy aliases of numpy objects were removed from scipy; same objects
from numpy.lib.scimath import log, sqrt	# py3 port: scipy.log/scipy.sqrt were the numpy.lib.scimath versions (complex result for negative input)


import os

def make_sure_path_exists(path):
	try: 
		os.makedirs(path)
	except OSError:
		if not os.path.isdir(path):
			raise


def is_iterable(a):	# duck typing
	try:
		iter(a)
	except TypeError:
		return False
	else:
		return True


# from mpmath import diff as derivative
# from scipy.misc import derivative as diff
def diff(func, x0, dx=1.0, n=1, args=(), order=3):	# py3 port: scipy.misc.derivative was removed from SciPy; this is its exact algorithm for n=1, order=3 (the only case used)
	if n != 1 or order != 3:
		raise NotImplementedError('py3 port: only n=1, order=3 of scipy.misc.derivative is re-implemented')
	weights = numpy.array([-1,0,1])/2.0
	val = 0.0
	ho = order >> 1
	for k in range(order):
		val += weights[k]*func(x0+(k-ho)*dx,*args)
	return val / numpy.prod((dx,)*n,axis=0)
derivative = lambda *a, **kw: diff(*a, **dict(kw, dx=1e-6))

from scipy.optimize import curve_fit, leastsq