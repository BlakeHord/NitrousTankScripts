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


#Constants:
Pa_over_psi = 6895
K_over_R = 5/9
C_to_K = 273.15
m_over_in = 0.0254

print()

fluid = 'NitrousOxide'
N2O_Tcrit = PropsSI(fluid,'Tcrit')
print(round(N2O_Tcrit - C_to_K, 2))

P_N2O_inlet = PropsSI('P','T',N2O_Tcrit,'Q',0,fluid)-10 #750 * Pa_over_psi #Pa - average bottle pressure
print(CP.PhaseSI('P',P_N2O_inlet,'Q',0,fluid))

print(P_N2O_inlet/Pa_over_psi)

pipe_id = 0.25 * m_over_in #m - 1/4" ID
pipe_A = pipe_id**2/4*np.pi #m2

N2O_cp = PropsSI('CPMASS','T',N2O_Tcrit-3,'P|liquid',P_N2O_inlet,fluid) #J/kg/K - mass specific constant pressure specific heat
N2O_cv = PropsSI('CVMASS','T',N2O_Tcrit,'P|liquid',P_N2O_inlet,fluid) #J/kg/K - mass specific constant volume specific heat
gamma = N2O_cp/N2O_cv #[-] - ratio of specific heats
N2O_R = 8.314472/PropsSI('MOLARMASS','T',N2O_Tcrit,'P|liquid',P_N2O_inlet,fluid) #J/kg-K - mass specific constant volume specific heat

N2O_cp_20 = PropsSI('CPMASS','T',N2O_Tcrit-15,'P|liquid',P_N2O_inlet,fluid) #J/kg/K - mass specific constant pressure specific heat

print("N2O cp at Tcrit(~35C) =", N2O_cp)
print("N2O cp at T=20C =", N2O_cp_20)

v_sound = np.sqrt(gamma*N2O_R*N2O_Tcrit) #m/s
print("Speed of sound =", v_sound, "m/s")

N2O_rho = PropsSI('D','T',N2O_Tcrit,'P|liquid',P_N2O_inlet,fluid) #kg/m^3 - density
N2O_rho_2 = N2O_density_liquid(N2O_Tcrit)

#print(N2O_rho)
#print(N2O_rho_2)

max_mdot = v_sound*N2O_rho_2*pipe_A

print("Max mdot =", max_mdot, "kg/s")

N2O_specific_gravity = N2O_rho_2/1000
print("Specific Gravity =", N2O_specific_gravity)