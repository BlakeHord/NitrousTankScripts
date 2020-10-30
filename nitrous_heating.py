import CoolProp.CoolProp as CP
from CoolProp.CoolProp import PropsSI
import numpy as np
import matplotlib
import matplotlib.pyplot as plt 

print()

#constants
Pa_over_psi = 6895
K_over_R = 5/9
C_to_K = 273.15
m_over_in = 0.0254
m3_s_over_gpm = 0.000063

al_t6_k = 167 #W/m-K


#options
nitrous_starting_temperature = 5 #C
ambient_temperature = 40 #C
include_solar_radiation = False
include_diffuse_radiation = False


#nitrous tank parameters
fluid = 'NitrousOxide'
T_ambient = ambient_temperature + C_to_K #K - temperature outside at launch site
T_N2O_0 = nitrous_starting_temperature + C_to_K #K - initial nitrous temperature
T_crit = PropsSI(fluid,'Tcrit')-3 #K - critical temperature of nitrous (plus several degree safety factor)

print("Nitrous starting temperature =", T_N2O_0 - C_to_K, "C")
print("Ambient temperature =", T_ambient - C_to_K, "C")
print("Nitrous critical temperature =", round(T_crit - C_to_K, 1), "C")

N2O_tank_od = 5 * m_over_in #m - od of nitrous tank
N2O_tank_wall_t = 0.125 * m_over_in #m - nitrous tank wall thickness (1/8")
N2O_tank_L = 1 #m - nitrous tank length (approximate)
m_N2O = 8.48 #kg - mass of nitrous


#calculate h for natural convection of air around cylinder  - eq 8.27 from Lienhard & Lienhard
Pr = 0.71 #[-] Prandtl Number for air
g = 9.81 #m/s2
beta = 1/T_ambient #1/K
dT = T_ambient - T_N2O_0 #K
L = 2 #m - vertical length of rocket
v = 1.81 * 10**-5 #Pa-s - kinematic viscosity of air
Gl = g*beta*dT*L**3/v**2 #[-] Grashof Number

Ra = Pr * Gl #[-] Rayleigh Number
Nu = 0.68 + 0.67*Ra**(1/4) * (1 + (0.492/Pr)**(9/16))**(-4/9)

k_air = 30 * 10**-3#W/K
h_air = Nu/L*k_air #W/m-K - Thermal Conductivity of air

#Not using internal h of nitrous (conservative estimate) because it's very tricky to calculate natural internal convection


#solar radiation (without shroud)
al_a = 0.65 #aluminum absorptivity (dull)
sun_flux = 1000 #W/m^2 - solar irradiance on clear day at sea level

sun_angle = np.radians(45) #rad - angle of solar radiation in sky (above horizon)
I_dir = sun_flux * np.cos(sun_angle) #W/m^2 - direct solar beam irradiance
tank_x_area = N2O_tank_L*N2O_tank_od #m^2 - tank cross sectional area
I_cylinder_dir = tank_x_area * I_dir #W

beta = np.radians(90) #rad - angle of surface from horizontal
Id = 60 #W/m^2 - approximate diffuse irradiance on horizontal plate
I_diffuse = 0.5 * Id * (1 + np.cos(beta)) #W/m^2 - irradiance from sky/diffused sunlights
tank_s_area = N2O_tank_L*N2O_tank_od*np.pi #m^2
I_cylinder_diffuse = tank_s_area * I_diffuse #W


#time-dependent calculation
U_tot = 1/(N2O_tank_wall_t/al_t6_k + 1/h_air) #W/m-K - total heat transfer coefficient (air + aluminum wall)
time_step = 5 #s
i = 0
T_N2O_i = np.array([T_N2O_0])

while (T_N2O_i[i] <= T_crit):
	Qdot = U_tot*(T_ambient - T_N2O_i[i]) #W/m^2 - heat transfer per unit area
	Q = Qdot*N2O_tank_L*np.pi*N2O_tank_od #W - energy per unit time
	N2O_cp = PropsSI('CPMASS','T|liquid',T_N2O_i[i],'Q',0,fluid) #J/kg-K - mass specific constant pressure specific heat
	
	if include_solar_radiation and include_diffuse_radiation:
		Q_tot = Q + I_cylinder_dir + I_cylinder_diffuse #W
	elif include_solar_radiation:
		Q_tot = Q + I_cylinder_dir #W
	else:
		Q_tot = Q #W

	energy = Q_tot * time_step #J - energy added in unit time
	dT = energy/(N2O_cp*m_N2O) #K - temperature change
	T_N2O_i = np.append(T_N2O_i,T_N2O_i[i]+dT)
	#print(T_N2O_i[i+1])

	i += 1

print("It took", i*time_step/60, "minutes to reach critical temperature")

fig, ax1 = plt.subplots()
t = np.arange(0,(i+1)*time_step, time_step)
color = 'tab:red'
ax1.plot(t/60,T_N2O_i-C_to_K, color=color)
ax1.set_xlabel('Time [min]')
ax1.set_ylabel('Temperature [C]', color=color)
if include_solar_radiation:
	title = 'Including Solar Radiation with T_i =' + str(nitrous_starting_temperature) + 'C'
else:
	title = 'Not Including Solar Radiation with T_i =' + str(nitrous_starting_temperature) + 'C'

ax1.set_title(title)
ax1.axhline(T_crit - C_to_K, color='tab:blue')
ax1.text(10,T_crit - C_to_K+0.05,'Critical Temperature Region', color='tab:blue')
ax1.grid()

plt.show()
