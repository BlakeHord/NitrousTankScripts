import CoolProp.CoolProp as CP
from CoolProp.CoolProp import PropsSI
import numpy as np
import matplotlib
import matplotlib.pyplot as plt 

def N2O_density_liquid(T):
	b1 = 1.72328
	b2 = -0.83950
	b3 = 0.51060
	b4 = -0.10412
	Tr = T/309.57
	rho_c = 452
	rho = rho_c*np.exp(b1*(1-Tr)**(1/3)+b2*(1-Tr)**(2/3)+b3*(1-Tr)+b4*(1-Tr)**(4/3))
	return rho


print()

#constants
Pa_over_psi = 6895
K_over_R = 5/9
C_to_K = 273.15
m_over_in = 0.0254
m3_s_over_gpm = 0.000063

fluid = "NitrousOxide"

fill_cv = 0.06
check_cv = 1.1

total_cv = np.sqrt(1/(1/fill_cv**2 + 1/check_cv**2)) #~0.0599

N2O_T_bottle = 25 + C_to_K #K
N2O_T_inlet = 5 + C_to_K #K
N2O_P_inlet = PropsSI('P','T',N2O_T_bottle,'Q',0,fluid)/Pa_over_psi #psi

N2O_rho = N2O_density_liquid(N2O_T_inlet) #kg/m^3
print("Density of Nitrous =", round(N2O_rho), "kg/m3")
SG = N2O_rho/1000

N2O_P_outlet = 14 #psi - atmospheric
dP = N2O_P_inlet - N2O_P_outlet

Q = total_cv*np.sqrt(dP/SG)*m3_s_over_gpm #m3/s
mdot = Q*N2O_rho #kg/s
print("Max mdot =", round(mdot,3), "kg/s")


print(N2O_density_liquid(5 + C_to_K))
print(N2O_density_liquid(30 + C_to_K))