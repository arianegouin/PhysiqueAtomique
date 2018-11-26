

class PB:
    def __init__(self):
        self.A = 208
        self.Z = 82

        self.Vo = 50    # MeV
        self.hw = 41 * self.A**(-1/3)
        self.D = 0.0225 * self.hw
        self.C = 0.1*self.hw


    def Enlsj(self,n,l,s,j):
        return -self.Vo + self.hw * (2 * (n - 1) + l + 3/2) - self.D * l * (l + 1) - C * (j * (j + 1) - l * (l + 1) - s * (s + 1)) / 2
