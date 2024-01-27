# Input and output of ODE

import ODEfunctions as odefx

ode_input = input("Enter the ODE (x as indep. var., y as dep. var): ")
ode = sp.Eq(sp.sympify(ode_input), 0)

