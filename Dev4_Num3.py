import numpy as np
from matplotlib import pyplot as plt
import scipy.special as ss
import scipy

# Omega = 2.23
Omega = 3.57

def f(x):
    return x

def f0(x):
    return - (1 - x**2)**0.5 * (np.tan((Omega * (1 - x**2))**0.5))**-1

def f1(x):
    # return ((np.tan((Omega * (1 - x**2))**0.5))**-1 * (1 - x**2)**-0.5 - (Omega * x**2 * (1 - x**2))**-1)**-1
    a = (1 - x**2)**0.5
    b = np.tan(Omega**0.5 * a)**-1
    c = Omega * x**2 * a**2
    return (b / a - 1 / c)**-1


# def F(L,x):
#     a = scipy.sqrt((Omega * (1 - x**2)))
#     b = 1j * Omega**0.5 * x
#
#     coeff = -1j * scipy.sqrt((1 - x**2))
#     firstTerm = ss.spherical_jn(L - 1,a) / ss.spherical_jn(L,a)
#     secdTerm = (ss.spherical_jn(L,b) + 1j*ss.spherical_yn(L,b)) / (ss.spherical_jn(L-1,b) + 1j * ss.spherical_yn(L-1,b))
#     return coeff * firstTerm * secdTerm

def omegafunc(omega, big_omega, L):
    root = scipy.sqrt(1 - omega**2)

    term1 = -1j * root
    term2 = ss.spherical_jn(L-1, big_omega * root) / ss.spherical_jn(L, big_omega * root)
    term3 = ss.spherical_jn(L, 1j * big_omega * omega) + 1j * ss.spherical_yn(L, 1j * big_omega * omega)
    term4 = ss.spherical_jn(L-1, 1j * big_omega * omega) + 1j * ss.spherical_yn(L-1, 1j * big_omega * omega)

    return term1 * term2 * term3 / term4

print(omegafunc(0.1,2.23,1))

# def BesselJn(L,x):
#     return ss.spherical_jn(L, x)

def BesselYn(L,x):
    return ss.spherical_yn(L,x)

w = [(i + 1)/100 for i in range(99)]

f = [f(i) for i in w]
f0 = [f0(i) for i in w]
# f1 = [f1(i) for i in w]

# Bessel = [BesselYn(0,i) for i in w]
# print(F(1,0.3))

# plt.plot(w,f)
# plt.plot(w,f0)
# plt.plot(w,f1)
# plt.plot(w,F1)

# plt.plot(w,Bessel)
# plt.plot(w,F0)

# plt.show()

