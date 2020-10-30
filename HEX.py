import CoolProp.CoolProp as CP
from CoolProp.CoolProp import PropsSI
import numpy as np
import matplotlib
import matplotlib.pyplot as plt 

'''
Script for modeling length of copper pipe in ice water bath required to cool nitrous oxide
Must set input temperature and desired output temperature

Assumptions:
- Constant ice water temperature (O C)
- Constant nitrous oxide mass flow rate
- Constant thermal conductivity of copper
- Constant pipee outside temperature
'''

#Function for calculating N2O saturated liquid viscosity
#from "Thermophysical properties of nitrous oxide" - http://edge.rit.edu/edge/P07106/public/Nox.pdf
def N2O_viscosity(T):
	b1 = 1.6089
	b2 = 2.0439
	b3 = 5.24
	b4 = 0.0293423
	Tc = 309.57 #K - critical temperature

	theta = (Tc-b3)/(T-b3)
	visc = b4*np.exp(b1*(theta-1)**(1/3)+b2*(theta-1)**(4/3))
	return visc

#Function for calculating N2O saturated liquid thermal conductivity
#from "Thermophysical properties of nitrous oxide" - http://edge.rit.edu/edge/P07106/public/Nox.pdf
def N2O_thermal_conductivity(T):
	b1 = 72.35
	b2 = 1.5
	b3 = -3.5
	b4 = 4.5
	Tc = 309.57 #K - critical temperature
	Tr = T/Tc
	k = b1*(1+b2*(1-Tr)**(1/3)+b3*(1-Tr)**(2/3)+b4*(1-Tr))
	return k


#Constants:
Pa_over_psi = 6895
K_over_R = 5/9
C_to_K = 273.15
m_over_in = 0.0254

print()

#Nitrous Oxide Parameters
fluid = 'NitrousOxide'
N2O_Tcrit = PropsSI(fluid,'Tcrit')

print("Initial Temperature =", round(N2O_Tcrit - C_to_K, 2), "C")

P_N2O_inlet = PropsSI('P','T',N2O_Tcrit,'Q',0,fluid) #Pa - find source
print("Initial Pressure =", round(P_N2O_inlet/Pa_over_psi), "psi")
T_N2O_inlet = N2O_Tcrit #K - worst case critical in bottle
T_N2O_goal = 5 + C_to_K #K  - arbitrary goal around temperature of ice
P_N2O_outlet = PropsSI('P','T',T_N2O_goal,'Q',0,fluid)

N2O_total_mass = 8.48 #kg - from propsim 18kNs 
N2O_fill_time = 5*60 #s - fill time, estimated as 5 minutes
N2O_mdot = N2O_total_mass/N2O_fill_time #kg/s - mass flow rate, must assume some kind of average constant, unless we want to simulate time-dependent flow

N2O_cp = PropsSI('CPMASS','T',(T_N2O_inlet+T_N2O_goal)/2,'P|liquid',(P_N2O_inlet+P_N2O_outlet)/2,fluid) #J/kg-K - mass specific constant pressure specific heat
Power_req = N2O_cp*N2O_mdot*(T_N2O_inlet-T_N2O_goal) #W
print("Cooling power required =", round(Power_req), "W")

#Ice Water Parameters
T_water = 0 + C_to_K #K
h_water = 1000 #W/m^2-K - convective heat transfer coefficient (limited natural convection)

#Pipe Parameters
Copper_k = 385 #W-K/m - thermal conductivity of copper
pipe_id = 0.25 * m_over_in #m - 1/4" ID
pipe_od = 0.375 * m_over_in #m - 3/8" OD
pipe_A = pipe_id**2/4*np.pi #m2
pipe_t = (pipe_od-pipe_id)/2 #m - wall thickness
pipe_e = 10**-4 * m_over_in #m - conservative roughness estimate for pipe

#Step parameters
T_N2O_i = np.array([T_N2O_inlet])
P_N2O_i = np.array([P_N2O_inlet])

#Get density and speed
N2O_rho = PropsSI('D','T',T_N2O_i[0],'P|liquid',P_N2O_i[0],fluid) #kg/m^3 - density
N2O_v = N2O_mdot/N2O_rho/pipe_A #m/s - velocity

i = 0
dx = 0.001 #m - step size (1 cm)

while (T_N2O_i[i] > T_N2O_goal):
	#1st, calculate h (convective heat transfer coefficient) for N2O
	mu_N2O = N2O_viscosity(T_N2O_i[i]) #Pa-s - dynamic viscosity
	Re_N2O = N2O_rho*N2O_v*pipe_id/mu_N2O
	k_N2O = N2O_thermal_conductivity(T_N2O_i[i]) #W/m-K - Thermal Conductivity of N2O

	if (Re_N2O > 2600): #flow is turbulent
		#from Lienhard & Lienhard "A Heat Transfer Textbook" eq. 7.41 & 7.42 for turbulent flow
		f  = (1.82*np.log10(Re_N2O) - 1.64)**-2
		N2O_cp = PropsSI('CPMASS','T',T_N2O_i[i],'P|liquid',P_N2O_i[i],fluid) #J/kg/K - mass specific constant pressure specific heat
		Pr_N2O = N2O_cp*mu_N2O/k_N2O #[-] - Prandtl Number
		Nu = (f/8)*(Re_N2O-1000)*Pr_N2O/(1+12.7*(f/8)**0.5*(Pr_N2O**(2/3)-1)) #[-] - Nusselt Number from Gnielinski Correlation
	else: #flow is laminar
		#from Lienhard & Lienhard "A Heat Transfer Textbook" eq. 7.23 for isothermal walls laminar flow
		Nu = 3.657 
	
	h_N2O = Nu*k_N2O/pipe_id
	
	#2nd, calculate the rate of energy addition into the N2O in the pipe
	#U_tot = 1/(1/h_water + pipe_t/Copper_k + 1/h_N2O) #W/m^2-K - overall heat transfer coefficient, convection from water
	U_tot = 1/(pipe_t/Copper_k + 1/h_N2O) #W/m^2-K - overall heat transfer coefficient, isothermal wall
	#print(U_tot)
	Qdot = U_tot*(T_water - T_N2O_i[i]) #W/m^2 - heat transfer per unit area
	Q = Qdot*dx*np.pi*pipe_od #W - energy per unit time, multiply by length of step

	#3rd, get N2O cp and find dT from addition of energy
	N2O_cp = PropsSI('CPMASS','T',T_N2O_i[i],'P|liquid',P_N2O_i[i]+10,fluid) #J/kg/K - mass specific constant pressure specific heat
	dT = Q/N2O_mdot/N2O_cp #K - change in temperature from heat addition
	T_N2O_i = np.append(T_N2O_i,T_N2O_i[i]+dT)

	#4th, re-calculate density and velocity from temperature change
	N2O_rho = PropsSI('D','T',T_N2O_i[i+1],'P|liquid',P_N2O_i[i],fluid) #kg/m^3 - density
	N2O_v = N2O_mdot/N2O_rho/pipe_A #m/s - velocity

	#5th, calculate pressure drop with Darcy Weisbach friction factor - may not hold because nitrous is self-pressurizing
	'''
	mu_N2O = N2O_viscosity(T_N2O_i[i+1]) #Pa-s - dynamic viscosity
	Re_N2O = N2O_rho*N2O_v*pipe_id/mu_N2O
	dwff = (-1.8*np.log10((pipe_e/pipe_id/3.7)**1.11 + 6.9/Re_N2O))**-2
	dP = dwff*dx/pipe_id*N2O_rho*N2O_v**2/2
	P_N2O_i = np.append(P_N2O_i,P_N2O_i[i]-dP)
	'''
	P_N2O_i = np.append(P_N2O_i,PropsSI('P','T',T_N2O_i[i+1],'Q',0,fluid)) #assumes saturated liquid

	#print(P_N2O_i[i+1])

	#6th, increment the counter
	i += 1

print("End Pressure =", round(P_N2O_i[i]/Pa_over_psi), "psi")

print("Length required =", dx*i, "m")
#print("Total pressure drop =", (P_N2O_i[0]-P_N2O_i[i])/Pa_over_psi, "psi")

fig, ax1 = plt.subplots()
x = np.arange(0,(i+1)*dx, dx)
color = 'tab:red'
ax1.plot(x,T_N2O_i-C_to_K, color=color)
ax1.set_xlabel('Distance in pipe [m]')
ax1.set_ylabel('Temperature [C]', color=color)
ax1.grid()

ax2 = ax1.twinx()
color = 'tab:blue'
ax2.plot(x,P_N2O_i/Pa_over_psi, color=color)
ax2.set_ylabel('Pressure [psi]', color=color)
ax2.set_ylim([0,P_N2O_i[0]/Pa_over_psi*1.1])
#plt.show()
