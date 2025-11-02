import preos 

from CoolProp.CoolProp import PropsSI
import CoolProp.CoolProp as cp
"""
def peng_robinson(fluid, temperature, pressure):
    critical_tempearture = cp.PropsSI("Tcrit", fluid)
    critical_pressure = cp.PropsSI("Pcrit", fluid)
    acentric_factor = cp.PropsSI("acentric", fluid)
    fugacity = fugacity("Isobutane")

    fluid = preos.Molecule(fluid, critical_tempearture, critical_pressure, acentric_factor)
    
    print(critical_tempearture, critical_pressure, acentric_factor)

    props = preos.preos(fluid, temperature, (pressure / 101325), plotcubic=False, printresults=True)

peng_robinson("Isobutane", 400, 300000)"""
"""
E = cp.PropsSI("H", "T", 305, "P", 300000, "butane")
print(
    cp.PropsSI("H", "P", 100000, "Q", 0, "water"),
    cp.PropsSI("Q", "P", 1000000, "H", 2600000, "water"),
    cp.PropsSI("T", "P", 300000, "Q", 1, "butane"),
    cp.PropsSI("P", "T", 300, "Q", 1, "butane"),
    cp.PropsSI("P", "T", 360, "Q", 1, "butane"),
    cp.PropsSI("P", "T", 360, "Q", 1, "butane"),
    cp.PropsSI("H", "T", 305, "P", 300000, "butane")   
      
      
      
      )

A = cp.PropsSI("T", "P", 300000, "Q", 0, "butane")
B = cp.PropsSI("T", "P", 300000, "Q", 1, "butane")
C = cp.PropsSI("P", "T", 300, "Q", 1, "butane")
D = cp.PropsSI("P", "T", 360, "Q", 1, "butane")


print(A,B,C,D)
"""
print (list([1, 2, 3]))
print (([1, 2, 3]))