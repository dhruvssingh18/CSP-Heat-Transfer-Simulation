import numpy as np

from scipy.integrate import quad

#Heat Transfer Fluid - Calcium Chloride(CaCL_2)
htf_vis = 2.03
htf_den = 1398 # Density of the material (kg/m^3)  
htf_velocity = 3 #m/s
htf_thermal_cond = 1.583 #(W/(m⋅K))

#Storage Fluid - Sodium Chloride(NaCL)
sf_vis = 2.03
sf_den = 1.51 # Density of the material (kg/m^3)  
sf_velocity = 1 #m/s
sf_thermal_cond = 1.004 #(W/(m⋅K))
#Simulation Parameters: 
absorbed = 836000 * 53 #Joules - 53(time of one cycle)
pipevolpos1 = 3.14
specific_heat_cap = 31.485 #J/kg°C

rey_num = (htf_den*htf_velocity*0.1)/(htf_vis) #reynolds number
darcy =  64/rey_num 
pr_num = (specific_heat_cap * htf_vis)/(htf_thermal_cond)
Nu_num = 3.66/(rey_num**1/3)*(pr_num**1/2)  
heat_cof = ((htf_thermal_cond)/(0.1)) * Nu_num

def heatlosspos1(t, T):
    T_surroundings = 700 #Celsius
    
    return -heat_cof * (T - T_surroundings)

T0 = 900  # Initial temp
timepos1 = 100/3 #Time of Fluid in Pos1 = Distance/Velocity 
time_intervalpos1 = (0, timepos1)  #Interval Step 

# Integrate the rate of heat loss function over the time interval
heat_losspos1, _ = quad(heatlosspos1, *time_intervalpos1, args=(T0,))

def fricheatlosspos1(): 
    pipe_length = 100 # Pipeline length of the system (m)
    pipe_rad = 0.1 #Radius of heat transfer pipe (m)
    pressurepos1 = darcy *((pipe_length*htf_den*htf_velocity**2)/(2*0.2))
    v = 0.0314
    
    return timepos1*(-1*pressurepos1*v)

def conductionhtrpos1():
    k = htf_thermal_cond
    A =  62.89
    L =  0.3
    hot = 900
    cold = 700
    return ((k * A * (hot - cold)) / L) * timepos1

def thermalradpos1():
    sigma = 5.67e-8  #(W/(m²·K⁴))
    A = 62.89 # Surface area of the object (m²)
    Ts = 1173.15 # Surface temperature of the object (K)
    Tsurr = 283.15  # Temperature of the surroundings (K)

    return 0.001*(sigma * A * ((Ts ** 4) - (Tsurr ** 4))) * timepos1
totalheatlosspos1 = ((-1*heat_losspos1) + (-1 * fricheatlosspos1()) + conductionhtrpos1() + thermalradpos1())
finalheatlosspos1 = (absorbed - totalheatlosspos1)
print("Absorbed:", absorbed)
print("Total Heat Loss:", totalheatlosspos1)
print("Heat Remaining in Position 1:", finalheatlosspos1)
print("     Convection Loss:", -1* heat_losspos1)
print("     Fluid Friction Loss:", -1 * fricheatlosspos1())
print("     Conduction Loss:", conductionhtrpos1())
print("     Thermal Radiation Loss:", thermalradpos1())
print("--------------------")
#########################################################################################

#Position 2: Heat Exchange
def heatlosspos2(t, T):
    T_surroundings = 700
    return -heat_cof * (T - T_surroundings)

T0 = 900  # Initial temp
timepos2 = 10/3 #Time of Fluid in Pos1 = Distance/Velocity 
time_intervalpos2 = (0, timepos2)  #Interval Step 

heat_losspos2, _ = quad(heatlosspos2, *time_intervalpos2, args=(T0,))

def conductionhtrpos2():
    k = htf_thermal_cond
    A = 39.96
    L =  0.3
    hot = 900
    cold = 700
    return ((k * A * (hot - cold)) / L) * timepos2

def thermalradpos2():
    sigma = 5.67e-8  # Stefan-Boltzmann constant (W/(m²·K⁴))
    A = 39.96  # Surface area of the object (m²)
    Ts = 1173.15 # Surface temperature of the object (K)
    Tsurr = 973.15  # Temperature of the surroundings (K)

    return 0.001*(sigma * A * ((Ts ** 4) - (Tsurr ** 4))) * timepos2

def fricheatlosspos2(): 
    pipe_length = 10 # Pipeline length of the system (m)
    pressurepos2 = darcy *((pipe_length*htf_den*htf_velocity**2)/(2*0.2))
    v = 0.0314
    
    return timepos2*(-1*pressurepos2*v)

fricheatlosspos2()

totalheatlosspos2 = ((-1*heat_losspos2) + (-1 * fricheatlosspos2()) + conductionhtrpos2() + thermalradpos2())
finalheatlosspos2 = (finalheatlosspos1 - totalheatlosspos2)
print("Total Heat Loss:", totalheatlosspos2)
print("Heat Remaining in Position 2:", finalheatlosspos2)
print("     Convection Loss:", -1* heat_losspos2)
print("     Fluid Friction Loss:", -1 * fricheatlosspos2())
print("     Conduction Loss:", conductionhtrpos2())
print("     Thermal Radiation Loss:", thermalradpos2())
print("----------------------")

###################################################################3

#Position 3: Final 
def heatlosspos3(t, T):
    T_surroundings = 700
    A = 31.48
    return -heat_cof * (T - T_surroundings)

T0 = 900  # Initial temp
timepos3 = 50/3 #Time of Fluid in Pos1 = Distance/Velocity 
time_intervalpos3 = (0, timepos3)  #Interval Step 

heat_losspos3, _ = quad(heatlosspos3, *time_intervalpos3, args=(T0,))

def conductionhtrpos3():
    k = htf_thermal_cond
    A = 31.48
    L =  0.3
    hot = 900
    cold = 700
    return ((k * A * (hot - cold)) / L) * timepos3

def thermalradpos3():
    sigma = 5.67e-8  # Stefan-Boltzmann constant (W/(m²·K⁴))
    A = 31.48  # Surface area of the object (m²)
    Ts = 1173.15 # Surface temperature of the object (K)
    Tsurr = 973.15  # Temperature of the surroundings (K)

    return 0.001*(sigma * A * ((Ts ** 4) - (Tsurr ** 4))) * timepos3

def fricheatlosspos3(): 
    pipe_length = 50 # Pipeline length of the system (m)
    pipe_rad = 0.1 #Radius of heat transfer pipe (m)
    pressurepos3 = darcy *((pipe_length*htf_den*htf_velocity**2)/(2*0.2))
    v = 0.0314

    
    return timepos3*(-1*pressurepos3*v)

fricheatlosspos3()

totalheatlosspos3 = ((-1*heat_losspos3) + (-1 * fricheatlosspos3()) + conductionhtrpos3() + thermalradpos3())
finalheatlosspos3 = (finalheatlosspos2 - totalheatlosspos3)
print("Total Heat Loss:", totalheatlosspos3)
print("Heat Remaining in Position 3:", finalheatlosspos3)
print("     Convection Loss:", -1* heat_losspos3)
print("     Fluid Friction Loss:", -1 * fricheatlosspos3())
print("     Conduction Loss:", conductionhtrpos3())
print("     Thermal Radiation Loss:", thermalradpos3())
print("-----------------")

####################################################################

#Heat Efficiency Output
finaleff = ((finalheatlosspos3)/absorbed)* 100
print('Total heat lost in simulation:', totalheatlosspos1+totalheatlosspos2+totalheatlosspos3)
print("Total heat transferred in simulation:", absorbed - ((totalheatlosspos1+totalheatlosspos2+totalheatlosspos3)))
print("The final efficiency of this system:", finaleff, )
print("Time", timepos1 + timepos2 + timepos3 )

import seaborn as sns
from matplotlib import pyplot as plt
import numpy as np

'''
# Define data using lists and variables


# Create DataFrame
df = pd.DataFrame({
    #'Simulation Positions': Positions,
    'Heat Loss': HeatL,
    'Frictional Heat Loss': FrictionHeatL, 
    'Conduction Heat Loss': ConductionHeatL,
    'Convection Heat Loss': Convection,
    'Thermal Radiation Heat Loss': ThermRad,
    
#Positions = ['Position1', 'Position 2', 'Position 3']
HeatL = [totalheatlosspos1, totalheatlosspos2, totalheatlosspos3]
FrictionHeatL = [-1*fricheatlosspos1(), -1*fricheatlosspos2(), -1*fricheatlosspos3()]
ConductionHeatL = [conductionhtrpos1(), conductionhtrpos2(), conductionhtrpos3()]
Convection = [-1*heat_losspos1, -1*heat_losspos2, -1*heat_losspos3]
ThermRad = [thermalradpos1(), thermalradpos2(), thermalradpos3()]
})  
'''

plt.rc('xtick', labelsize=14)
plt.rc('ytick', labelsize=14)

ylabels = ['Position 1', 'Position 2', 'Position 3']
xlabels = ['Total Heat Loss', 'Friction', 'Conduction ', 'Convection', 'ThermalRad']

orders = np.array([[totalheatlosspos1,-1*fricheatlosspos1(),conductionhtrpos1(), -1*heat_losspos1, thermalradpos1()], 
                   [totalheatlosspos2, -1*fricheatlosspos2(), conductionhtrpos2(),-1*heat_losspos2, thermalradpos2()],
                   [totalheatlosspos3, -1*fricheatlosspos3() , conductionhtrpos3(), -1*heat_losspos3, thermalradpos3()]]
                   
                 )

plt.figure(figsize=(8,5))
sns.heatmap(orders, 
            cmap='YlOrBr',
            vmin=0,
            xticklabels=xlabels,
            yticklabels=ylabels,
            annot=True,
            square=True,
            annot_kws={'fontsize':14, 'fontweight': 'bold', 'color': 'black'}
           )
plt.yticks(rotation=0)
plt.tick_params(
    which='both',      
    bottom=False,      
    left=False,      
    labelbottom=False,
    labeltop=True) 
plt.tight_layout();
plt.show() 