# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 10:10:37 2023

@author: chonk
"""
from Simulation import ThermalSimulation
from Input_Heat_Waves import HeatingFunction

#%%
# Create a heating function
heating = HeatingFunction(peak=75, room_temperature=25, frequency = 1/10)

#%%
# Create an instance of the ThermalSimulation class

simulation = ThermalSimulation(block_length=0.10, block_width=0.05, block_height=0.05, total_time=500.0,
                               grid_size=30, density=2700, specific_heat=897, initial_conductivity=237,
                               heating_function=heating.square_wave_temperature)

# Drill holes at 2.5cm, 5cm, and 7.5cm with a depth of 2.5cm in the 'x' face of the block
hole_positions = [0.025, 0.05, 0.075]  # Modify the hole positions as needed
hole_depth = 0.025  # Modify the hole depth as needed
hole_face = 'x'  # Specify the face of the block where the holes are drilled ('x', 'y', or 'z')
simulation.drill_holes(hole_positions, hole_depth, hole_face)

# Run the thermal simulation
simulation.simulate()