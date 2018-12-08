"ATOMIC AND NUCLEAR PHYSICS - HOMEWORK 5"
import numpy as np

def E_shellModel(A, n, l, s, j):
    V0 = 50    # Mev
    hw = 41 * A ** (-1 / 3)    # Mev
    D = 0.0225 * hw
    C = 0.1 * hw
    return - V0 + hw * (2 * (n - 1) + l + 3 / 2) - D * l * (l + 1) - C * (j * (j + 1) - l * (l + 1) - s * (s + 1)) / 2

def E_SEMF(Z, A):
    aV = 15.75    # Mev
    aS = 17.80    # Mev
    aC = 0.711    # Mev
    aA = 94.780    # Mev
    aP = 11.18    # Mev
    return ( aV * A
           - aS * A ** (2/3)
           - aC * Z**2 / A**(1/3)
           - aA * (Z - A/2)**2 / A
           + aP * ((-1)**Z + (-1)**(A - Z)) / (2 * A**(1/2)) )

def R_alpha(Z, A, E_alpha):
    c = 299792458    # m/s
    r0 = 1
    # m_alpha = 2 * 938.27205 + 2 * 939.56541    # Mev/c^2
    # E_alpha = m_alpha * c**2
    m_alpha = E_alpha / c**2
    Z_alpha = 4
    alpha = 1 / 137.035999
    return ( c / (r0 * A**(1/3))
             * ( E_alpha / (2 * m_alpha * c**2))**0.5
             * np.exp(- np.pi * (Z - Z_alpha) * Z_alpha * alpha * (2 * m_alpha * c**2 / E_alpha)**0.5))


"NUMBER 1"
Q_gamma = E_shellModel(A=7, n=1, l=1, s=0.5, j=2.5) - E_shellModel(A=7, n=1, l=2, s=0.5, j=0.5)
print('1b)', 'Q_gamma: ', Q_gamma, 'Mev')


"NUMBER 2"
Z = 298
A = 92
Q_alpha = E_SEMF(Z - 2, A - 4) + E_SEMF(2, 4) - E_SEMF(Z, A)
print('2)', 'Q_alpha: ', Q_alpha, 'MeV')
R_alpha = R_alpha(Z, A, 4.27)
print('  ', 'R_alpha', R_alpha, 'Mev')

"NUMBER 3"
me = 0.51099891    # Mev/c^2
mp = 938.27205    # Mev/c^2
mn = 939.56541    # Mev/c^2
c = 299792458    # m/s
Z = 29
A = 11
Q_betaPlus = E_SEMF(Z - 1, A) - E_SEMF(Z,A) + (mp - mn - me) * c**2
print('3)', 'Q_betaPlus', Q_betaPlus, 'MeV')


