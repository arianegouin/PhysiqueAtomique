import numpy as np
from matplotlib import pyplot as plt
import scipy.special as ss
import scipy


def f(x):
    return x
    
def f0(x, Omega2):
    return - (1 - x**2)**0.5 * (np.tan((Omega2 * (1 - x**2))**0.5))**-1

def F(L, x, Omega2):
     a = scipy.sqrt((Omega2 * (1 - x**2)))
     b = 1j * Omega2**0.5 * x
     coeff = -1j * scipy.sqrt((1 - x**2))
     firstTerm = ss.spherical_jn(L - 1,a) / ss.spherical_jn(L,a)
     secdTerm = (ss.spherical_jn(L,b) + 1j*ss.spherical_yn(L,b)) / (ss.spherical_jn(L-1,b) + 1j * ss.spherical_yn(L-1,b))
     return coeff * firstTerm * secdTerm

# Omega^2 values
singulet = 2.23 * 4
triplet = 3.57 * 4

# w values knowing w range is [0,1]
w = np.array([(i + 1)/10000 for i in range(9999)])

# left term of the transcendant equation
f = np.array([f(i) for i in w])

# right term of the transcendant equation
f_00 = np.array([f0(i, singulet) for i in w])       # when (L,S) = (0,0)
f_01 = np.array([f0(i, triplet) for i in w])        # when (L,S) = (0,1)
F_10 = np.array([F(1,i, singulet) for i in w])      # when (L,S) = (1,0)
F_11 = np.array([F(1,i, triplet) for i in w])       # when (L,S) = (1,1)

# when (L,S) = (0,0)
plt.figure()
plt.plot(w,f)
plt.plot(w,f_00)
idx = np.argwhere(np.diff(np.sign(f_00 - f))).flatten()     # index where intersection is
print('intersection 00', w[idx])                               # value of w for this index
plt.plot(w[idx], f[idx], 'ro')
plt.ylim([-1,10])
plt.savefig('00', bbox_inches = 'tight')

# when (L,S) = (0,1)
plt.figure()
plt.plot(w,f)
plt.plot(w,f_01)
idx = np.argwhere(np.diff(np.sign(f_01 - f))).flatten()     # index where intersection is
print('intersection 01', w[idx])                               # value of w for this index
plt.plot(w[idx], f[idx], 'ro')
plt.ylim([-1,10])
plt.savefig('01', bbox_inches = 'tight')

# when (L,S) = (1,0)
plt.figure()
plt.plot(w,f)
plt.plot(w,np.real(F_10))
idx = np.argwhere(np.diff(np.sign(F_10 - f))).flatten()     # index where intersection is
print('intersection 10', w[idx])                               # value of w for this index
plt.plot(w[idx], f[idx], 'ro')
plt.ylim([-10,10])
plt.savefig('10', bbox_inches = 'tight')

# when (L,S) = (1,1)
plt.figure()
plt.plot(w,f)
plt.plot(w,np.real(F_11))
idx = np.argwhere(np.diff(np.sign(F_11 - f))).flatten()     # index where intersection is
print('intersection 11', w[idx])                               # value of w for this index
plt.plot(w[idx], f[idx], 'ro')
plt.ylim([-1,10])
plt.savefig('11', bbox_inches = 'tight')




ww = np.array([(i + 1)/100 for i in range(99)])
ff = np.array([f(i) for i in ww])
F_LS = np.array([F(6,i, triplet) for i in ww])       # when (L,S) = (1,1)
plt.figure()
plt.plot(ww,ff)
plt.plot(ww,np.real(F_LS))
idx = np.argwhere(np.diff(np.sign(F_LS - ff))).flatten()     # index where intersection is
print('intersection', ww[idx])                               # value of w for this index
plt.plot(ww[idx], ff[idx], 'ro')
plt.savefig('try', bbox_inches = 'tight')