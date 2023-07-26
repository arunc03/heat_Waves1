# -*- coding: utf-8 -*-
"""
Created on Sun Jul 16 20:18:01 2023

@author: chonk
"""

import pickle
import numpy as np
import matplotlib.pyplot as plt

from Simulation import ThermalSimulation
from Input_Heat_Waves import HeatingFunction
from Fourier import FourierAnalysis
  
#%%
# 10x5x5 block, three drills on top
with open('block_three_holes.pkl', 'rb') as sim_file_a:
    sim1a = pickle.load(sim_file_a)

sim1a.plot_temperature_harmonics()

x = sim1a.plot_times
y1 = sim1a.T_start
y2 = sim1a.T_quarter
y3 = sim1a.T_half
y4 = sim1a.T_three_quarter

plt.figure(figsize=(12, 4))

plt.subplot(2,2,1)
plt.plot(x[5813:-1], y1[5813:-1])
plt.xlabel('Time (s)')
plt.ylabel('Temperature (°C)')

plt.subplot(2,2,2)
plt.plot(x[5813:-1], y2[5813:-1])
plt.xlabel('Time (s)')
plt.ylabel('Temperature (°C)')

plt.subplot(2,2,3)
plt.plot(x[5813:-1], y3[5813:-1])
plt.xlabel('Time (s)')
plt.ylabel('Temperature (°C)')

plt.subplot(2,2,4)
plt.plot(x[5813:-1], y4[5813:-1])
plt.xlabel('Time (s)')
plt.ylabel('Temperature (°C)')

plt.tight_layout()
plt.show()

print("Fourier Analysis - 10x5x5 block, three drills at top")
fourier1a_quarter = FourierAnalysis(np.array(sim1a.plot_times[5813:-1]), np.array(sim1a.T_quarter[5813:-1]), 5)
print(fourier1a_quarter.calculate(0.025, 1, 75))

fourier1a_half = FourierAnalysis(np.array(sim1a.plot_times[5813:-1]), np.array(sim1a.T_half[5813:-1]), 5)
print(fourier1a_half.calculate(0.05, 1, 75))

fourier1a_three_quarter = FourierAnalysis(np.array(sim1a.plot_times[5813:-1]), np.array(sim1a.T_three_quarter[5813:-1]), 5)
print(fourier1a_three_quarter.calculate(0.075, 1, 75))

#%%
# 10x5x5 block, one drill along centre
with open('block_one_hole.pkl', 'rb') as sim_file_b:
    sim2a = pickle.load(sim_file_b)
    
sim2a.plot_temperature_harmonics()

x = sim2a.plot_times
y1 = sim2a.T_start
y2 = sim2a.T_quarter
y3 = sim2a.T_half
y4 = sim2a.T_three_quarter

plt.figure(figsize=(12, 4))

plt.subplot(2,2,1)
plt.plot(x[5813:-1], y1[5813:-1])
plt.xlabel('Time (s)')
plt.ylabel('Temperature (°C)')

plt.subplot(2,2,2)
plt.plot(x[5813:-1], y2[5813:-1])
plt.xlabel('Time (s)')
plt.ylabel('Temperature (°C)')

plt.subplot(2,2,3)
plt.plot(x[5813:-1], y3[5813:-1])
plt.xlabel('Time (s)')
plt.ylabel('Temperature (°C)')

plt.subplot(2,2,4)
plt.plot(x[5813:-1], y4[5813:-1])
plt.xlabel('Time (s)')
plt.ylabel('Temperature (°C)')

plt.tight_layout()
plt.show()

print("Fourier Analysis - 10x5x5 block, one drill along centre")
fourier2a_quarter = FourierAnalysis(np.array(sim2a.plot_times[5813:-1]), np.array(sim2a.T_quarter[5813:-1]), 5)
print(fourier2a_quarter.calculate(0.025, 1, 75))

fourier2a_half = FourierAnalysis(np.array(sim2a.plot_times[5813:-1]), np.array(sim2a.T_half[5813:-1]), 5)
print(fourier2a_half.calculate(0.05, 1, 75))

fourier2a_three_quarter = FourierAnalysis(np.array(sim2a.plot_times[5813:-1]), np.array(sim2a.T_three_quarter[5813:-1]), 5)
print(fourier2a_three_quarter.calculate(0.075, 1, 75))

#%%
# 10x3x3 block, three drills at top
with open('block_three_hole_3.pkl', 'rb') as sim_file_c:
    sim1b = pickle.load(sim_file_c)
    
sim1b.plot_temperature_harmonics()

x = sim1b.plot_times
y1 = sim1b.T_start
y2 = sim1b.T_quarter
y3 = sim1b.T_half
y4 = sim1b.T_three_quarter

plt.figure(figsize=(12, 4))

plt.subplot(2,2,1)
plt.plot(x[7912:-1], y1[7912:-1])
plt.xlabel('Time (s)')
plt.ylabel('Temperature (°C)')

plt.subplot(2,2,2)
plt.plot(x[7912:-1], y2[7912:-1])
plt.xlabel('Time (s)')
plt.ylabel('Temperature (°C)')

plt.subplot(2,2,3)
plt.plot(x[7912:-1], y3[7912:-1])
plt.xlabel('Time (s)')
plt.ylabel('Temperature (°C)')

plt.subplot(2,2,4)
plt.plot(x[7912:-1], y4[7912:-1])
plt.xlabel('Time (s)')
plt.ylabel('Temperature (°C)')

plt.tight_layout()
plt.show()

print("Fourier Analysis - 10x5x5 block, three drills at top")
fourier1a_quarter = FourierAnalysis(np.array(sim1b.plot_times[7912:-1]), np.array(sim1b.T_quarter[7912:-1]), 5)
print(fourier1a_quarter.calculate(0.025, 1, 75))

fourier1a_half = FourierAnalysis(np.array(sim1b.plot_times[7912:-1]), np.array(sim1b.T_half[7912:-1]), 5)
print(fourier1a_half.calculate(0.05, 1, 75))

fourier1a_three_quarter = FourierAnalysis(np.array(sim1b.plot_times[7912:-1]), np.array(sim1b.T_three_quarter[7912:-1]), 5)
print(fourier1a_three_quarter.calculate(0.075, 1, 75))

#%%
# 10x3x3 block, one drill along centre
with open('block_one_hole_3.pkl', 'rb') as sim_file_d:
    sim2b = pickle.load(sim_file_d)

sim2b.plot_temperature_harmonics()

x = sim2b.plot_times
y1 = sim2b.T_start
y2 = sim2b.T_quarter
y3 = sim2b.T_half
y4 = sim2b.T_three_quarter

plt.figure(figsize=(12, 4))

plt.subplot(2,2,1)
plt.plot(x[5813:-1], y1[5813:-1])
plt.xlabel('Time (s)')
plt.ylabel('Temperature (°C)')

plt.subplot(2,2,2)
plt.plot(x[5813:-1], y2[5813:-1])
plt.xlabel('Time (s)')
plt.ylabel('Temperature (°C)')

plt.subplot(2,2,3)
plt.plot(x[5813:-1], y3[5813:-1])
plt.xlabel('Time (s)')
plt.ylabel('Temperature (°C)')

plt.subplot(2,2,4)
plt.plot(x[5813:-1], y4[5813:-1])
plt.xlabel('Time (s)')
plt.ylabel('Temperature (°C)')

plt.tight_layout()
plt.show()

print("Fourier Analysis - 10x3x3 block, one drill along centre")
fourier2a_quarter = FourierAnalysis(np.array(sim2b.plot_times[5813:-1]), np.array(sim2b.T_quarter[5813:-1]), 5)
print(fourier2a_quarter.calculate(0.025, 1, 75))

fourier2a_half = FourierAnalysis(np.array(sim2b.plot_times[5813:-1]), np.array(sim2b.T_half[5813:-1]), 5)
print(fourier2a_half.calculate(0.05, 1, 75))

fourier2a_three_quarter = FourierAnalysis(np.array(sim2b.plot_times[5813:-1]), np.array(sim2b.T_three_quarter[5813:-1]), 5)
print(fourier2a_three_quarter.calculate(0.075, 1, 75))

#%%
heating = HeatingFunction(peak=75, room_temperature=25, frequency=1/5)

#%%
simulation1 = ThermalSimulation(block_length=0.10, block_width=0.03, block_height=0.03, total_time=500.0,
                                grid_size=35, density=2700, specific_heat=897, initial_conductivity=237,
                                heating_function=heating.sawtooth_wave_temperature, boundary_conductivity=0.024,
                                boundary_specific_heat=1000, boundary_density=15)
simulation1.drill_hole('z', (0.025, 0.015, 0.015), 0.015, 0.003)
simulation1.drill_hole('z', (0.05, 0.015, 0.015), 0.015, 0.003)
simulation1.drill_hole('z', (0.075, 0.015, 0.015), 0.015, 0.003)
simulation1.simulate()

#%%
# simulation2 = ThermalSimulation(block_length=0.10, block_width=0.03, block_height=0.03, total_time=500.0,
#                                 grid_size=35, density=2700, specific_heat=897, initial_conductivity=237,
#                                 heating_function=heating.sawtooth_wave_temperature)
# simulation2.drill_hole('x', (0.025, 0.025, 0.025), 0.075, 0.0025)
# simulation2.simulate()