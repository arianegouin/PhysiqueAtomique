import matplotlib.pyplot as plt
import math
import numpy as np
from scipy import special

Omega_sing = math.sqrt(2.23*4)
Omega_trip = math.sqrt(3.57*4)
Vos_sing = 4*11.8
Vos_trip = 4*35.1

w = np.linspace(1e-10, 0.9999999999, 100000)
L_0_sing = []
L_0_trip = []
L_1_sing = []
L_1_trip = []
L_i_sing = []
L_i_trip = []

def spherique(L, arg1, arg2):
    # singulet
    j_Ls = special.spherical_jn(L, arg1)
    y_Ls = special.spherical_yn(L, arg1)
    # triplet
    j_Lt = special.spherical_jn(L, arg2)
    y_Lt = special.spherical_yn(L, arg2)

    return [j_Ls, y_Ls, j_Lt, y_Lt]

for i in range(0, len(w)):
    #### L = 0 ####
    L_0_sing.append(-math.sqrt(1-w[i]**2)*
                    (math.tan(Omega_sing*math.sqrt(1-(w[i])**2))**-1))
    L_0_trip.append(-math.sqrt(1-w[i]**2)*
                    (math.tan(Omega_trip*math.sqrt(1-(w[i])**2))**-1))

    #### L = 1 ####
    L_1_sing.append(((math.tan(Omega_sing*math.sqrt(1-(w[i])**2))**-1)/
                     math.sqrt(1-w[i]**2) -
                     1/(Omega_sing*(w[i]**2)*(1-w[i]**2)))**-1)
    L_1_trip.append(((math.tan(Omega_trip*math.sqrt(1-(w[i])**2))**-1)/
                     math.sqrt(1-w[i]**2) -
                     1/(Omega_trip*(w[i]**2)*(1-w[i]**2)))**-1)
    #### L > 1 ####
    # singulet
    L = 2
    [j_Ls, _, j_Lt, _] = spherique(L, Omega_sing*math.sqrt(1-(w[i])**2),
                                   Omega_trip*math.sqrt(1-(w[i])**2))
    [j_Lsm, _, j_Ltm, _] = spherique(L-1, Omega_sing*math.sqrt(1-(w[i])**2),
                                     Omega_trip*math.sqrt(1-(w[i])**2))
    [j_LsC, y_LsC,j_LtC, y_LtC] = spherique(L, 1j*Omega_sing*w[i],
                                            1j*Omega_trip*w[i])
    [j_LsCm, y_LsCm, j_LtCm, y_LtCm] = spherique(L-1, 1j*Omega_sing*w[i],
                                                 1j*Omega_trip*w[i])

    L_i_sing.append(-1j*math.sqrt(1-w[i]**2) * j_Lsm/j_Ls *
                    (j_LsC+1j*y_LsC)/(j_LsCm+1j*y_LsCm))
    L_i_trip.append(-1j*math.sqrt(1-w[i]**2) * j_Ltm/j_Lt *
                    (j_LtC+1j*y_LtC)/(j_LtCm+1j*y_LtCm))

# Intersections
idx1 = list(np.argwhere(np.diff(np.sign(w - L_0_sing))).flatten())
idx2 = list(np.argwhere(np.diff(np.sign(w - L_0_trip))).flatten())
idx3 = list(np.argwhere(np.diff(np.sign(w - L_1_sing))).flatten())
idx4 = list(np.argwhere(np.diff(np.sign(w - L_1_trip))).flatten())
idx5 = list(np.argwhere(np.diff(np.sign(w - L_i_sing))).flatten())
idx6 = list(np.argwhere(np.diff(np.sign(w - L_i_trip))).flatten())
print("Les solutions 'w' sont (singulet, triplet):\n", w[idx1], w[idx2], "\n",
      w[idx3], w[idx4], '\n', w[idx5], w[idx6])

#########################################
# Graphe L = 0 singulet
plt.plot(w, L_0_sing)
plt.plot(w, w)
plt.plot(w[idx1], w[idx1], 'ro')
plt.title("L = 0, singulet")
#plt.show()

# Graphe L = 0 triplet
plt.plot(w, L_0_trip)
plt.plot(w, w)
plt.plot(w[idx2], w[idx2], 'ro')
plt.title("L = 0, triplet")
#plt.show()
########################################
# Graphe L = 1 singulet
plt.plot(w, L_1_sing)
plt.plot(w, w)
plt.plot(w[idx3], w[idx3], 'ro')
plt.title("L = 1, singulet")
#plt.show()

# Graphe L = 1 triplet
plt.plot(w, L_1_trip)
plt.plot(w, w)
plt.plot(w[idx4], w[idx4], 'ro')
plt.title("L = 1, triplet")
#plt.show()
#########################################
# Graphe L = i singulet
plt.plot(w, np.real(L_i_sing))
plt.plot(w, w)
plt.plot(w[idx5], w[idx5], 'ro')
plt.title("L = i, singulet")
plt.ylim((0, 1))
#plt.show()

# Graphe L = i triplet
plt.plot(w, np.real(L_i_trip))
plt.plot(w, w)
plt.plot(w[idx6], w[idx6], 'ro')
plt.title("L = i, triplet")
plt.ylim((0, 1))
#plt.show()


############################################################

Vo = 50
hw = 6.9198
n = [3, 2, 3, 1, 1, 1, 2, 3, 1, 2, 1, 2, 1, 2, 1, 1, 2, 1, 1, 1, 1]
l = [1, 3, 1, 6, 5, 5, 2, 0, 4, 2, 4, 1, 3, 1, 3, 2, 0, 2, 1, 1, 0]
j = [0.5, 2.5, 1.5, 6.5, 4.5, 5.5, 1.5, 0.5, 3.5, 2.5, 4.5, 0.5, 2.5,
     1.5, 3.5, 1.5, 0.5, 2.5, 0.5, 1.5, 0.5]
s = 0.5

E_neutron = []
E_proton = []
for i in range(0, len(n)):
    E_neutron.append((-Vo + hw*(2*(n[i]-1) + l[i] + 1.5) - 0.0225*hw*l[i]*(l[i] + 1) -
                0.05*hw*(j[i]*(j[i] + 1) - l[i]*(l[i] + 1) - s*(s+1)))*(2*j[i] + 1))
    if i >= 5:
        E_proton.append((-Vo + hw*(2*(n[i]-1) + l[i] + 1.5) - 0.0225*hw*l[i]*(l[i] + 1) -
                0.05*hw*((-1)**(j[i] - l[i] - 0.5))*(2*l[i] - j[i] + 0.5))*(2*j[i] + 1))

E_tot = sum(E_proton) + sum(E_neutron)

print("Énergie liaison totale = ", E_tot, "MeV")
print(E_proton[0]/12, E_neutron[0]/2)


