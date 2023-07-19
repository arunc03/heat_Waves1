# -*- coding: utf-8 -*-
"""
Created on Sun Jul 16 20:18:01 2023

@author: chonk
"""
import pickle
from Simulation import ThermalSimulation
from Input_Heat_Waves import HeatingFunction
from Fourier import FourierAnalysis

#%%
heating = HeatingFunction(peak=75, room_temperature=25, frequency=1/5)

#%%
simulation1 = ThermalSimulation(block_length=0.10, block_width=0.05, block_height=0.05, total_time=100.0,
                                grid_size=50, density=2700, specific_heat=897, initial_conductivity=237,
                                heating_function=heating.square_wave_temperature)
simulation1.drill_hole('z', (0.025, 0.025, 0.025), 0.025, 0.003)
simulation1.drill_hole('z', (0.05, 0.025, 0.025), 0.025, 0.003)
simulation1.drill_hole('z', (0.075, 0.025, 0.025), 0.025, 0.003)
simulation1.simulate()

#%%
simulation1.plot_temperature_harmonics()