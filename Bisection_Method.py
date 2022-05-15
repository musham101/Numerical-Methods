class BiSection:
    a = 0
    b = 0
    tolerance = 0
    max_itrations = 0
    eq = "0"
    
    def setA(self, a):
        self.a = a

    def setA(self, b):
        self.b = b

    def setTolerance(self, tolerance):
        self.tolerance = tolerance

    def setMaxItrations(self, max_itrations):
        self.max_itrations = max_itrations

    def setEq(self, eq):
        self.eq = eq

    def calculate(self):
        i = 1
        print("{:>0} {:>15} {:>15} {:>15} {:>15} {:>15} {:>15}\n" .format("#", "a", "b", "c", "f(a)", "f(b)", "f(c)"))
        while i <= self.max_itrations:
            X = self.a
            FA = eval(self.eq)
            X = self.b
            FB = eval(self.eq)
            c = (self.a + self.b) / 2
            X = c
            FC = eval(self.eq)
            print("{:>0} {:>15} {:>15} {:>15} {:>15} {:>15} {:>15}\n" .format(int(i), self.a, self.b, c, FA, FB, FC))
            if FC == 0 or abs(FC) < self.tolerance:
                print("X = " + str(c) + '\n')
                return
            if FA * FC < 0:
                self.b = c
            else:
                self.a = c
            i = i + 1
        print("\t-----Failed to get root for the given number of itrations\n")
        return

while True:
    print("\tNM Project\n")
    print("Select the method you want to use\n")
    print("- Press 1 for biseection Method\n")
    print("- Press 0 to exit\n")
    op = int(print("INPUT: "))


