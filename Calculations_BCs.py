# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 10:31:34 2023

@author: chonk
"""

import numpy as np
from Calculations import HeatEquationSolver


class HeatEquationSolverWithBoundary(HeatEquationSolver):
    def __init__(self, dx, dy, dz, dt, num_steps, temperature, heating_function, specific_heat, initial_conductivity, initial_density, boundary_conductivity, boundary_specific_heat, boundary_density):
        super().__init__(dx, dy, dz, dt, num_steps, temperature, heating_function, specific_heat, initial_conductivity, initial_density)
        self.boundary_conductivity = boundary_conductivity
        self.boundary_specific_heat = boundary_specific_heat
        self.boundary_density = boundary_density

    def solve(self, time):
        inside = (slice(1, -1),) * 3
        
        new_temperature = np.copy(self.temperature)
        self.thermal_conductivity = self.thermal_conductivity_field(self.temperature, self.initial_conductivity)
        self.density = self.density_field(self.temperature)

        self.alpha = self.thermal_conductivity / (self.density * self.specific_heat)
        self.dt = min(self.dx**2, self.dy**2, self.dz**2) / (6 * self.alpha.max())

        laplacian = ((self.temperature[:-2, 1:-1, 1:-1] - 2*self.temperature[1:-1, 1:-1, 1:-1] + self.temperature[2:, 1:-1, 1:-1]) / self.dx**2 +
                     (self.temperature[1:-1, :-2, 1:-1] - 2*self.temperature[1:-1, 1:-1, 1:-1] + self.temperature[1:-1, 2:, 1:-1]) / self.dy**2 +
                     (self.temperature[1:-1, 1:-1, :-2] - 2*self.temperature[1:-1, 1:-1, 1:-1] + self.temperature[1:-1, 1:-1, 2:]) / self.dz**2)

        new_temperature[inside] += self.dt * self.alpha[inside] * laplacian

        boundary_alpha = self.boundary_conductivity / (self.boundary_density * self.boundary_specific_heat)
        insulation_factor = boundary_alpha / self.alpha[inside].max()
        
        new_temperature[-1, :, :] = new_temperature[-2, :, :] * insulation_factor  # x = block_length
        new_temperature[:, 0, :] = new_temperature[:, 1, :] * insulation_factor  # y = 0
        new_temperature[:, -1, :] = new_temperature[:, -2, :] * insulation_factor  # y = block_width
        new_temperature[:, :, 0] = new_temperature[:, :, 1] * insulation_factor  # z = 0
        new_temperature[:, :, -1] = new_temperature[:, :, -2] * insulation_factor  # z = block_height

        heater_temp = self.heating_function(time)
        new_temperature[0, :, :] = heater_temp  # x = 0 (heater)
        
        self.temperature, new_temperature = new_temperature, self.temperature

        return self.temperature