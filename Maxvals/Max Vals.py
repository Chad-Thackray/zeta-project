import sympy as sym
import numpy as np
import matplotlib.pyplot as plt

p = sym.symbols("p")

psi0 = sym.cos(2*sym.pi *(p**2 -p -0.0625))/sym.cos(2*sym.pi*p)
psi0 = sym.series(psi0, p, x0 = 0.5, n = 20)
psi0 = psi0.removeO()
psi0 = sym.collect(psi0, p)


c1 = -1/(2**5 * 3 * sym.pi**2)* sym.diff( psi0, p, 3)
c2 = (1/(2**11 * 3**2 * sym.pi**4)) * sym.diff(psi0,p,6) + (1/(2**6 * sym.pi**2 ))* sym.diff(psi0, p, 2)
c3 = (-1)/(2**16 * 3**4 * sym.pi**6)*sym.diff(psi0,p,9) - 1/(2**8 * 3 * 5 * sym.pi**4)*sym.diff(psi0,p,5) - 1/(2**6 * sym.pi**2)*sym.diff(psi0,p,1)
c4 = 1/(2**23 * 3**5 * sym.pi**8)*sym.diff(psi0,p,12) + 11/(2**17 * 3**2 * 5 * sym.pi**6)*sym.diff(psi0,p,8) + 19/(2**13 * 3* sym.pi**4)*sym.diff(psi0,p,4) + 1/(2**7 * sym.pi**2)*psi0


interval = np.linspace(0,1,1000)
lists = []
function = c4
function = sym.series(function, x0=0.5, n=6)
function = function.removeO()

for x in interval:
    lists.append(sym.N(function.subs(p,x), 3))

print(max(lists))





