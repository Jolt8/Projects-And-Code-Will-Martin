import sympy
from sympy import symbols

x, y, z = symbols("x y z")

integer = 10

f = 2*x + integer + z
eq = sympy.Eq(y,f)
print(sympy.solve(eq,x)[0])
 