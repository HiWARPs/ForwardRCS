from scipy import pi
from val import pm

h = 4.135667662*10**-15 |pm| 0.000000025*10**-15	# eV * s (CODATA)
c = 299792458	# m/s (exact) (CODATA)
hc = h*c/(2*pi)*10**6	# fm * GeV (CODATA)

alpha = 0.0072973525664 |pm| 0.0000000000017	# (CODATA)

Mp = 0.9382720813 |pm| 0.0000000058	# GeV (CODATA)
Mn = 0.9395654133 |pm| 0.0000000058	# GeV (CODATA)
Mpip = 0.13957018 |pm| 0.00000035	# GeV (PDG)
Mpi0 = 0.1349766 |pm| 0.0000006	# GeV (PDG)

kappa = 1.7928473508 |pm| 0.0000000085	# (CODATA)
IGDH = 4*pi**2 * alpha/(2*Mp**2) * 10**4*hc**2 * kappa**2

def cm2lab(w):
	return 0.5 * (w**2/Mp - Mp)

def lab2cm(k):
	return (Mp*(Mp + 2*k))**0.5

def fcm2lab(w):
	return 0.5 * (float(w)**2/Mp.val - Mp.val)

def flab2cm(k):
	return (Mp.val*(Mp.val + 2*float(k)))**0.5

gA = 1.2695 |pm| 0.0029
fpi = 0.09221 |pm| 0.00014

e_Mp = - alpha/Mp * hc**2 * 10**4

E0 = fcm2lab(Mp + Mpi0)	# ~0.14469 GeV
E1 = fcm2lab(Mp + 2*Mpi0)	# ~0.30879 GeV
E2 = 2.	# GeV
