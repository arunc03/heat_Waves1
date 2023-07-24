# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 10:31:34 2023

@author: chonk
"""

import numpy as np

class HeatEquationSolverWithBoundary:
    def __init__(self, dx, dy, dz, dt, num_steps, temperature, heating_function, specific_heat, initial_conductivity, initial_density, boundary_conductivity, boundary_specific_heat, boundary_density):
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.dt = dt
        self.alpha = 1.0
        self.num_steps = num_steps
        self.temperature = temperature
        self.heating_function = heating_function
        self.specific_heat = specific_heat
        self.initial_conductivity = initial_conductivity
        self.initial_density = initial_density
        self.boundary_conductivity = boundary_conductivity
        self.boundary_specific_heat = boundary_specific_heat
        self.boundary_density = boundary_density

    def thermal_conductivity_field(self, temperature, initial_conductivity):
        """
        Function to calculate the thermal conductivity field.

        Parameters
        ----------
        temperature (float): The current temperature.
        initial_conductivity (float): The initial thermal conductivity.

        Returns
        ----------
        conductivity (float): The calculated thermal conductivity.
        """
        
        conductivity  = initial_conductivity - 0.001 * (temperature - 25) 

        return conductivity

    def density_field(self, temperature):
        """
        Function to calculate the density field.

        Parameters
        ----------
        temperature (float): The current temperature.

        Returns
        ----------
        density (float): The calculated density.
        """
        
        density_0 = self.initial_density
        alpha_rho = 0.001
        density = density_0 * (1 - alpha_rho * (temperature - 25))
        
        return density

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