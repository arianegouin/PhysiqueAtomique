

class PB:
    def __init__(self):
        self.A = 208    # Nombre de masse
        self.Z = 82     # Nombre atomique

        self.Vo = 50    # MeV
        self.hw = 41 * self.A**(-1/3)    # MeV
        self.D = 0.0225 * self.hw
        self.C = 0.1 * self.hw

        self.protons = [                 # couches complètées pour 82 protons
                     [1, 0, 0.5],
                     [1, 1, 1.5],
                     [1, 1, 0.5],
                     [1, 2, 2.5],
                     [2, 0, 0.5],
                     [1, 2, 1.5],
                     [1, 3, 3.5],
                     [2, 1, 1.5],
                     [1, 3, 2.5],
                     [2, 1, 0.5],
                     [1, 4, 4.5],
                     [2, 2, 2.5],
                     [1, 4, 3.5],
                     [3, 0, 0.5],
                     [2, 2, 1.5],
                     [1, 5, 5.5],         # 82 protons
                     ]

        self.neutrons = [                 # couches complétées pour 126 neutrons
                     [1, 0, 0.5],
                     [1, 1, 1.5],
                     [1, 1, 0.5],
                     [1, 2, 2.5],
                     [2, 0, 0.5],
                     [1, 2, 1.5],
                     [1, 3, 3.5],
                     [2, 1, 1.5],
                     [1, 3, 2.5],
                     [2, 1, 0.5],
                     [1, 4, 4.5],
                     [2, 2, 2.5],
                     [1, 4, 3.5],
                     [3, 0, 0.5],
                     [2, 2, 1.5],
                     [1, 5, 5.5],         # 82 neutrons

                     [2, 3, 3.5],
                     [1, 5, 4.5],
                     [1, 6, 6.5],
                     [3, 1, 1.5],
                     [2, 3, 2.5],
                     [3, 1, 0.5]          # 126 neutrons
                     ]

    def Enlsj(self,n,l,s,j):
        return - self.Vo + self.hw * (2 * (n - 1) + l + 3/2) - self.D * l * (l + 1) - self.C * (j * (j + 1) - l * (l + 1) - s * (s + 1)) / 2

    def lastProton(self):
        shell = self.protons[-1]
        n, l, s, j = shell[0], shell[1], 0.5, shell[2]
        return - self.Enlsj(n, l, s, j)

    def lastNeutron(self):
        shell = self.neutrons[-1]
        n, l, s, j = shell[0], shell[1], 0.5, shell[2]
        return - self.Enlsj(n, l, s, j)

    def bindingEnergy(self):
        protonsEnergy = 0
        for shell in self.protons:
            n, l, s, j = shell[0], shell[1], 0.5, shell[2]
            protonsEnergy += self.Enlsj(n, l, s, j) * (2 * j + 1)
        neutronsEnergy = 0
        for shell in self.neutrons:
            n, l, s, j = shell[0], shell[1], 0.5, shell[2]
            neutronsEnergy += self.Enlsj(n, l, s, j)  * (2 * j + 1)
        return - (protonsEnergy + neutronsEnergy)


pb = PB()
print('The last proton energy is', pb.lastProton(), 'MeV.')
print('The last neutron energy is', pb.lastNeutron(), 'MeV.')
print('The binding energy is', pb.bindingEnergy(), 'MeV.')


