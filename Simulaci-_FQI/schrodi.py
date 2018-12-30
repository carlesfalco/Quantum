from numpy import sin, cos, sinh, cosh, tanh, sqrt, exp # Vectorized
from numpy import asarray, zeros, linspace
from scipy.integrate import quad

# Valors numèrics, variables globals.

k = 5.12665			  	# Factor d'unitats (  eV^(-1/2)*nm^(-1)  )
hbar = 6.582e-4			# hbar (eV · ps)
L = 6.0				   	# Mitja longitud de la caixa (nm)
mu = -L/2			   	# <x> inicial (nm)
l = 0.6				   	# Mitja amplada de la barrera (nm)
sigma = 0.6				# Stdev inicial (nm)
T = 0.5				   	# Energia cinetica (eV)
V0 = 1.0				# Barrera de potencial (eV)
N = 50				# Numero d'estats propis que utilitzarem
dE = 0.0001		  		# Precisio en energies (eV)

# Equacions transcedentals de l'energia
# a -> 0 < E < V0
# b -> V0 < E, E>0
# c -> V0 < E < 0

def even_a(E):
	return (sqrt(V0-E))*tanh(k*(sqrt(V0-E))*l)*sin(k*(sqrt(E))*(L-l)) + \
	(sqrt(E))*cos(k*(sqrt(E))*(L-l))
	
def odd_a(E):
	return (sqrt(V0-E))*sin(k*(sqrt(E))*(L-l)) + \
		(sqrt(E))*tanh(k*(sqrt(V0-E))*l)*cos(k*(sqrt(E))*(L-l))

def even_b(E):
	return (sqrt(E-V0))*sin(k*(sqrt(E-V0))*l)*sin(k*(sqrt(E))*(L-l)) - \
		(sqrt(E))*cos(k*(sqrt(E-V0))*l)*cos(k*(sqrt(E))*(L-l))

def odd_b(E):
	return (sqrt(E-V0))*cos(k*(sqrt(E-V0))*l)*sin(k*(sqrt(E))*(L-l)) + \
		(sqrt(E))*sin(k*(sqrt(E-V0))*l)*cos(k*(sqrt(E))*(L-l))

def even_c(E):
    return sqrt(E-V0)*tanh(k*sqrt(-E)*(l-L))*sin(k*sqrt(E-V0)*l) + \
        sqrt(-E)*cos(k*sqrt(E-V0)*l)

def odd_c(E):
    return sqrt(E-V0)*tanh(k*sqrt(-E)*(l-L))*cos(k*sqrt(E-V0)*l) - \
        sqrt(-E)*sin(k*sqrt(E-V0)*l)
	
# Allowed energy values

def find_bound_states():
	if V0 > 0:
		even = even_a
		odd = odd_a
	elif V0 < 0:
		even = even_c
		odd = odd_c
	else:
		return []
	E_spectra = []
	j = 0
	Eb = [0,V0]
	E = min(Eb)
	while ( j < N ):
		e, o = even(E), odd(E)
		E += dE
		if E > max(Eb):
			break
		if e*even(E) < 0:
			E_spectra.append(E)
			j += 1
		if o*odd(E) < 0:
			E_spectra.append(E)
			j += 1
	return E_spectra
	
def find_states():
	E_spectra = find_bound_states()
	j = len(E_spectra)
	Eb = [0,V0]
	E = max(Eb)
	while( j < N ): # Finding scattering states
		e, o = even_b(E), odd_b(E)
		E += dE
		if e*even_b(E) < 0:
			E_spectra.append(E)
			j += 1
		if o*odd_b(E) < 0:
			E_spectra.append(E)
			j += 1
	return E_spectra
		
# Defining wavefunctions, vectorized functions.
	
def phi_even_a(x,E):
	x = asarray(x)
	phi = zeros(x.shape)
	phi += ( (-L < x) & (x < -l) ) * sin(k*(sqrt(E))*(x+L))
	phi += ( abs(x) < l ) * sin(k*(sqrt(E))*(L-l))*cosh(k*(sqrt(V0-E))*x)/(cosh(k*(sqrt(V0-E))*l))
	phi += ( (x > l) & (x < L) ) * (-sin(k*(sqrt(E))*(x-L)))
	return phi

def phi_even_b(x,E):
	x = asarray(x)
	phi = zeros(x.shape)
	phi += ( (-L < x) & (x < -l) ) * sin(k*(sqrt(E))*(x+L))
	phi += ( abs(x) < l ) * sin(k*(sqrt(E))*(L-l))*cos(k*(sqrt(E-V0))*x)/(cos(k*(sqrt(E-V0))*l))
	phi += ( (x > l) & (x < L) ) * (-sin(k*(sqrt(E))*(x-L)))
	return phi
    
def phi_even_c(x,E):
	x = asarray(x)
	phi = zeros(x.shape)
	phi += ( (-L < x) & (x < -l) ) * sinh(k*sqrt(-E)*(x+L))
	phi += ( abs(x) < l ) * sinh(k*(sqrt(-E))*(L-l))*cos(k*(sqrt(E-V0))*x)/(cos(k*(sqrt(E-V0))*l))
	phi += ( (x > l) & (x < L) ) * (-sinh(k*(sqrt(-E))*(x-L)))
	return phi

def phi_odd_a(x,E):
	x = asarray(x)
	phi = zeros(x.shape)
	phi += ( (-L < x) & (x < -l) ) * sin(k*(sqrt(E))*(x+L))
	phi += ( abs(x) < l ) * (-sin(k*(sqrt(E))*(L-l))*sinh(k*(sqrt(V0-E))*x)/(sinh(k*(sqrt(V0-E))*l)))
	phi += ( (x > l) & (x < L) ) * sin(k*(sqrt(E))*(x-L))
	return phi

def phi_odd_b(x,E):
	x = asarray(x)
	phi = zeros(x.shape)
	phi += ( (-L < x) & (x < -l) ) * sin(k*(sqrt(E))*(x+L))
	phi += ( abs(x) < l ) * (-sin(k*(sqrt(E))*(L-l))*sin(k*(sqrt(E-V0))*x)/(sin(k*(sqrt(E-V0))*l)))
	phi += ( (x > l) & (x < L) ) * sin(k*(sqrt(E))*(x-L))
	return phi
    
def phi_odd_c(x,E):
	x = asarray(x)
	phi = zeros(x.shape)
	phi += ( (-L < x) & (x < -l) ) * sinh(k*(sqrt(-E))*(x+L))
	phi += ( abs(x) < l ) * (-sinh(k*(sqrt(-E))*(L-l))*sin(k*(sqrt(E-V0))*x)/(sin(k*(sqrt(E-V0))*l)))
	phi += ( (x > l) & (x < L) ) * sinh(k*(sqrt(-E))*(x-L))
	return phi

def phi_odd(x,E): 
	if E < V0:
		return phi_odd_a(x,E)
	elif E > 0:
		return phi_odd_b(x,E)
	else:
		return phi_odd_c(x,E)

def phi_even(x,E): 
	if E < V0:
		return phi_even_a(x,E)
	elif E > 0:
		return phi_even_b(x,E)
	else:
		return phi_even_c(x,E)
		
def evaluate_wavefunction(x,E):
	phi = zeros((len(E),len(x)))
	for i in range(0,len(E)):
		phi[i] = ( i%2 == 0 ) * phi_even(x,E[i]) + ( i%2 == 1 ) * phi_odd(x,E[i])
		phi[i] /= sqrt( 2*L/len(x) * sum(phi[i]*phi[i]) ) 
	return phi

# Defining gaussian 

def gaussian(x):
	return exp( - (x - mu)**2/4/sigma**2 )
	
# Computing basis coefficients

def basis_coefficients(phi_basis,phi_init,x):
	cn = []
	for phi_n in phi_basis:
		cn.append( 2*L/len(x) * sum(phi_n*phi_init) )
	return cn

# Time evolution operator

def time_evolve(x, phi_basis, cn, E, t):
	return sum( [ phi_basis[i]*cn[i]*exp(-1j*E[i]*t/hbar) for i in range(0,len(E)) ] )

