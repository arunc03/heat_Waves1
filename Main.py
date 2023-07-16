# -*- coding: utf-8 -*-
"""
Created on Sun Jul 16 20:18:01 2023

@author: chonk
"""
from Simulation import ThermalSimulation
from Input_Heat_Waves import HeatingFunction

# %%
heating = HeatingFunction(peak=75, room_temperature=25, frequency=1/5)

# %%
# simulation = ThermalSimulation(block_length=0.10, block_width=0.05, block_height=0.05, total_time=500.0,
#                                grid_size=30, density=2700, specific_heat=897, initial_conductivity=237,
#                                heating_function=heating.square_wave_temperature)

# hole_positions = [0.025]  
# hole_width = 0.0001  
# hole_depth = 0.05  
# hole_face = '+y'
# simulation.drill_holes(hole_positions, hole_depth, hole_width, hole_face)
# simulation.simulate()

# %%
simulation2 = ThermalSimulation(block_length=0.10, block_width=0.05, block_height=0.05, total_time=1000.0,
                                grid_size=30, density=2700, specific_heat=897, initial_conductivity=237,
                                heating_function=heating.square_wave_temperature)
hole_positions = [0.025, 0.05]  
hole_width = 0.002 
hole_depth = 0.05
hole_face = '+y'
# simulation2.drill_holes(hole_positions, hole_depth, hole_width, hole_face)
simulation2.visualize_holes()

simulation2.simulate()

#%%
simulation.plot_temperature_harmonics()
simulation2.plot_temperature_harmonics()
# %%