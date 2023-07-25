# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 10:57:49 2023

@author: chonk
"""
import numpy as np

class HeatEquationSolver:
    def __init__(self, dx, dy, dz, dt, num_steps, temperature, heating_function, specific_heat, initial_conductivity, initial_density):
        """
        Initialize the HeatEquationSolver object.

        Parameters
        ----------
        dx (float): Spatial step size in the x-direction.
        dy (float): Spatial step size in the y-direction.
        dz (float): Spatial step size in the z-direction.
        dt (float): Temporal step size.
        num_steps (int): Number of time steps to be solved.
        temperature (ndarray): Initial temperature distribution as a 3D numpy array.
        heating_function (callable): A function that takes the current time as input and returns the temperature at the heater position (x=0).
        specific_heat (float): Specific heat capacity of the material.
        initial_conductivity (float): Initial thermal conductivity of the material.
        initial_density (float): Initial density of the material.
        """
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

    def thermal_conductivity_field(self, temperature, initial_conductivity):
        """
        Calculate the thermal conductivity field based on the current temperature.

        Parameters
        ----------
        temperature (float): The current temperature.
        initial_conductivity (float): The initial thermal conductivity.

        Returns
        ----------
        conductivity (float): The calculated thermal conductivity.
        """
        conductivity = initial_conductivity - 0.001 * (temperature - 25)
        return conductivity

    def density_field(self, temperature):
        """
        Calculate the density field based on the current temperature.

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
        """
        Solve the heat equation for the given time.

        Parameters
        ----------
        time (float): The current time.

        Returns
        ----------
        temperature (ndarray): The updated temperature distribution as a 3D numpy array.
        """
        inside = (slice(1, -1),) * 3

        new_temperature = np.copy(self.temperature)
        self.thermal_conductivity = self.thermal_conductivity_field(self.temperature, self.initial_conductivity)
        self.density = self.density_field(self.temperature)

        self.alpha = self.thermal_conductivity / (self.density * self.specific_heat)
        self.dt = min(self.dx ** 2, self.dy ** 2, self.dz ** 2) / (6 * self.alpha.max())

        laplacian = ((self.temperature[:-2, 1:-1, 1:-1] - 2 * self.temperature[1:-1, 1:-1, 1:-1] + self.temperature[2:, 1:-1, 1:-1]) / self.dx ** 2 +
                     (self.temperature[1:-1, :-2, 1:-1] - 2 * self.temperature[1:-1, 1:-1, 1:-1] + self.temperature[1:-1, 2:, 1:-1]) / self.dy ** 2 +
                     (self.temperature[1:-1, 1:-1, :-2] - 2 * self.temperature[1:-1, 1:-1, 1:-1] + self.temperature[1:-1, 1:-1, 2:]) / self.dz ** 2)

        new_temperature[inside] += self.dt * self.alpha[inside] * laplacian

        new_temperature[-1, :, :] = new_temperature[-2, :, :]  # x = block_length
        new_temperature[:, 0, :] = new_temperature[:, 1, :]  # y = 0
        new_temperature[:, -1, :] = new_temperature[:, -2, :]  # y = block_width
        new_temperature[:, :, 0] = new_temperature[:, :, 1]  # z = 0
        new_temperature[:, :, -1] = new_temperature[:, :, -2]  # z = block_height

        heater_temp = self.heating_function(time)
        new_temperature[0, :, :] = heater_temp  # x = 0 (heater)

        self.temperature, new_temperature = new_temperature, self.temperature

        return self.temperature