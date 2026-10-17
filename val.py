from functools import wraps

def _binop(method):
	@wraps(method)
	def wrapped(self, v):
		try:
			return method(self, v)
		except:
			return NotImplemented
	return wrapped

pm_str = '+-'


class Value:
	
	def __init__(self, val=None, err=None):
		try:
			self.val, self.err = map(float, val.replace('+-', r'\pm').split(r'\pm'))
		except:
			try:
				self.val, self.err = map(float, val)
			except (IndexError, TypeError, ValueError):
				self.val = 0. if val is None else float(val)
				self.err = 0. if err is None else float(err)
			else:
				try:
					self.err = (self.err**2 + float(err)**2)**0.5	# e.g. stat. + syst.
				except (TypeError, ValueError):
					pass
		self.err = abs(self.err)
	
	def __getitem__(self, i):
		if i in self.__dict__:
			return self.__dict__[i]
		else:
			return (self.val, self.err)[i]
	
	def __str__(self):
		return ('%s ' + pm_str + ' %s') % tuple(map(repr, (self.val, self.err)))
	
	def __float__(self):
		return float(self.val)
	
	def __int__(self):
		return int(self.__float__())
	
	def binop_v(method):
		return wraps(method)(_binop(lambda self, v: method(self, Value(v))))
	
	def binop_f(method):
		return wraps(method)(_binop(lambda self, v: method(self, float(v))))
	
	@binop_v
	def __eq__(self, v):
		return self.val == v.val and self.err == v.err
	
	def __ne__(self, v):
		return not self == v
	
	@binop_f
	def __lt__(self, v):
		return self.val < v
	
	@binop_f
	def __gt__(self, v):
		return self.val > v
	
	@binop_f
	def __le__(self, v):
		return self.val <= v
	
	@binop_f
	def __ge__(self, v):
		return self.val >= v
	
	def __nonzero__(self):
		return self.val != 0 or self.err != 0
	
	__bool__ = __nonzero__	# py3 port
	
	def __neg__(self):
		return Value(-self.val, self.err)
	
	def __pos__(self):
		return Value(self)
	
	def __abs__(self):
		return Value(abs(self.val), self.err)
	
	@binop_v
	def __add__(self, v):
		return Value(self.val + v.val, (self.err**2 + v.err**2)**0.5)
		# return Value(self.val + v.val, self.err + v.err)	# assuming here 2 values are "fully" correlated (other operations assume uncorrelated values)
	
	def __radd__(self, v):
		return self + v
	
	def __sub__(self, v):
		return self + (-v)
	
	def __rsub__(self, v):
		return v + (-self)
	
	@binop_v
	def __mul__(self, v):
		return Value(self.val * v.val, ((self.val*v.err)**2 + (self.err*v.val)**2)**0.5)
	
	def __rmul__(self, v):
		return self * v
	
	@binop_v
	def __div__(self, v):
		return Value(self.val / v.val, ((self.val*v.err/v.val**2)**2 + (self.err/v.val)**2)**0.5)
	
	def __rdiv__(self, v):
		return Value(v)/self
	
	__truediv__, __rtruediv__ = __div__, __rdiv__	# py3 port: '/' calls __truediv__
	
	def __pow__(self, v):
		return Value(self.val**v, abs(v*self.val**(v-1)*self.err))


val = Value


class _PmInfix:
	
	def __ror__(self, x):
		self.x = x
		return self
	
	def __or__(self, y):
		if 'x' in self.__dict__:
			x = self.__dict__['x']
			del self.__dict__['x']
		else:
			x = None
		return Value(x, y)
	
	__rlshift__ = __rand__ = __rxor__ = __ror__
	__rshift__ = __and__ = __xor__ = __or__


pm = _PmInfix()