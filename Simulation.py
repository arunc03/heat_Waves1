# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 10:09:36 2023

@author: chonk
"""
import numpy as np
import matplotlib.pyplot as plt
from Calculations import HeatEquationSolver

class ThermalSimulation:
    """
    Class for thermal simulation of a block.
    """

    def __init__(self, block_length, block_width, block_height, total_time, grid_size,
                 density, specific_heat, initial_conductivity, heating_function,
                 boundary_conductivity, boundary_specific_heat, boundary_density):
        """
        Initialize the ThermalSimulation class.

        Parameters
        ----------
        block_length (float): The length of the block in meters.
        block_width (float): The width of the block in meters.
        block_height (float): The height of the block in meters.
        total_time (float): The total time of simulation in seconds.
        grid_size (int): The size of the grid.
        density (float): The initial density of the block in kg/m^3.
        specific_heat (float): The specific heat of the block in J/(kg*K).
        initial_conductivity (float): The initial thermal conductivity of the block in W/(m*K).
        """

        self.block_length = block_length
        self.block_width = block_width
        self.block_height = block_height
        self.total_time = total_time
        self.grid_size = grid_size
        self.initial_density = np.full(
            (grid_size, grid_size, grid_size), density)
        self.specific_heat = specific_heat
        self.initial_conductivity = np.full(
            (grid_size, grid_size, grid_size), initial_conductivity)
        self.boundary_conductivity = boundary_conductivity
        self.boundary_specific_heat = boundary_specific_heat
        self.boundary_density = boundary_density
        self.heating_function = heating_function
        self.T_start = []
        self.T_quarter = []
        self.T_half = []
        self.T_three_quarter = []
        self.T_end = []
        self.plot_times = []

        self.dx = self.block_length / self.grid_size
        self.dy = self.block_width / self.grid_size
        self.dz = self.block_height / self.grid_size

        self.temperature = np.ones(
            (self.grid_size, self.grid_size, self.grid_size)) * 25

        self.thermal_conductivity = self.initial_conductivity - \
            0.001 * (self.temperature - 25)
        self.density = self.initial_density * \
            (1 - 0.001 * (self.temperature - 25))
        self.alpha = self.thermal_conductivity / \
            (self.density * self.specific_heat)

        self.dt = min(self.dx**2, self.dy**2, self.dz**2) / \
            (6 * self.alpha.max())
        self.num_steps = int(self.total_time / self.dt)
        self.time = np.linspace(0, self.total_time, self.num_steps)

        # Determine the range of cells to be drilled
    def drill_hole(self, face, start_point, depth, width):
        """
        Function to drill a hole filled with air on any face.

        Parameters
        ----------
        face (str): The face to be drilled on. This can be 'x', 'y', or 'z'.
        start_point (tuple): The start of the hole. This is a tuple of x, y, z coordinates.
        depth (float): The depth of the hole in meters.
        width (float): The width of the hole in meters.
        """
        start_x, start_y, start_z = start_point
        start_x_idx = int(start_x / self.dx)
        start_y_idx = int(start_y / self.dy)
        start_z_idx = int(start_z / self.dz)

        # Determine the range of cells to be drilled
        if face == 'x':
            depth_cells = int(depth / self.dx)
            width_cells = int(width / self.dy)
            end_cell = start_x_idx + depth_cells
            self.initial_conductivity[start_x_idx:end_cell, start_y_idx:start_y_idx+width_cells, start_z_idx:start_z_idx+width_cells] = 0.024
            self.initial_density[start_x_idx:end_cell, start_y_idx:start_y_idx+width_cells, start_z_idx:start_z_idx+width_cells] = 1.224
        
        elif face == 'y':
            depth_cells = int(depth / self.dy)
            width_cells = int(width / self.dx)
            end_cell = start_y_idx + depth_cells
            self.initial_conductivity[start_x_idx:start_x_idx+width_cells, start_y_idx:end_cell, start_z_idx:start_z_idx+width_cells] = 0.024
            self.initial_density[start_x_idx:start_x_idx+width_cells, start_y_idx:end_cell, start_z_idx:start_z_idx+width_cells] = 1.224
       
        elif face == 'z':
            depth_cells = int(depth / self.dz)
            width_cells = int(width / self.dx)
            end_cell = start_z_idx + depth_cells
            self.initial_conductivity[start_x_idx:start_x_idx+width_cells, start_y_idx:start_y_idx+width_cells, start_z_idx:end_cell] = 0.024
            self.initial_density[start_x_idx:start_x_idx+width_cells, start_y_idx:start_y_idx+width_cells, start_z_idx:end_cell] = 1.224


    def simulate(self):
        """
        Function to run the thermal simulation.

        Functionality
        ----------
        Inside the loop, it performs the following steps:
        - Copy the current temperature grid to a new variable.
        - Update the thermal conductivity, density and alpha.
        - Calculate the Laplacian of the temperature.
        - Calculate the new temperature grid.
        - Enforce boundary conditions.
        - Swap the new and old temperature grids.
        - Plot the temperature at certain intervals.
        """
        solver = HeatEquationSolver(self.dx, self.dy, self.dz, self.dt, self.num_steps, self.temperature,
                                    self.heating_function, self.specific_heat, self.initial_conductivity, 
                                    self.initial_density)

        for step in range(self.num_steps):
            self.temperature = solver.solve(self.time[step])

            if step % 100 == 0:
                self.plot_temperature(step)
                self.T_start.append(self.temperature[0, int(self.grid_size / 2), int(self.grid_size / 2)])
                self.T_quarter.append(self.temperature[int(0.25*self.grid_size), int(self.grid_size / 2)+3, int(self.grid_size / 2)+3])
                self.T_half.append(self.temperature[int(0.5*self.grid_size), int(self.grid_size / 2)+3, int(self.grid_size / 2)+3])
                self.T_three_quarter.append(self.temperature[int(0.75*self.grid_size), int(self.grid_size / 2)+3, int(self.grid_size / 2)+3])
                self.T_end.append(self.temperature[-1, int(self.grid_size / 2)+3, int(self.grid_size / 2)+3])
                self.plot_times.append(self.time[step])
        
        plt.show()

        self.plot_temperature_harmonics()

    def plot_temperature_harmonics(self):
        """
        Function to plot the temperature.

        Functionality
        ----------
        The function creates a 5-part plot:
        - A line plot showing temperature with time at the start.
        - A line plot showing temperature with time at a quarter of the block length.
        - A line plot showing temperature with time at a half of the block length.
        - A line plot showing temperature with time at three quarters of the block length.
        - A line plot showing temperature with time at the block length.        
        """

        times_start = np.linspace(0, int(self.total_time), 100000)
        T_start = self.heating_function(times_start)
        plt.figure(figsize=(12, 4))

        plt.subplot(151)
        plt.plot(times_start, T_start)
        plt.xlabel("Time (s)")
        plt.ylabel('Temperature (°C)')
        title1 = "Distance: 0 m"
        plt.title(title1)

        plt.subplot(152)
        plt.plot(self.plot_times, self.T_quarter)
        plt.xlabel("Time (s)")
        plt.ylabel('Temperature (°C)')
        title2 = "Distance: " + str(0.25*self.block_length) + " m"
        plt.title(title2)

        plt.subplot(153)
        plt.plot(self.plot_times, self.T_half)
        plt.xlabel("Time (s)")
        plt.ylabel('Temperature (°C)')
        title3 = "Distance: " + str(0.5*self.block_length) + " m"
        plt.title(title3)

        plt.subplot(154)
        plt.plot(self.plot_times, self.T_three_quarter)
        plt.xlabel("Time (s)")
        plt.ylabel('Temperature (°C)')
        title4 = "Distance: " + str(0.75*self.block_length) + " m"
        plt.title(title4)

        plt.subplot(155)
        plt.plot(self.plot_times, self.T_end)
        plt.xlabel("Time (s)")
        plt.ylabel('Temperature (°C)')
        title5 = "Distance: " + str(self.block_length) + " m"
        plt.title(title5)

        plt.tight_layout()

    def plot_temperature(self, step):
        """
        Function to plot the temperature.

        Parameters
        ----------
        step (int): The current time step.

        Functionality
        ----------
        The function creates a 3-part plot:
        - A 2D color plot showing temperature on the x-y plane.
        - A 2D color plot showing temperature on the x-z plane.
        - A line plot showing temperature along the x axis.
        """
        title = "Time (s): " + str(self.time[step])
        plt.figure(figsize=(12, 4))
        plt.suptitle(title, fontsize=12)
        
        plt.subplot(131)
        plt.imshow(self.temperature[:, :, int(self.grid_size / 2)].T, cmap='hot', origin='lower',
                   extent=[0, self.block_length, 0, self.block_width], vmin=25, vmax=80)
        plt.colorbar()
        plt.xlabel('x (m)')
        plt.ylabel('y (m)')
        plt.title('x-y Plane')

        plt.subplot(132)
        plt.imshow(self.temperature[:, int(self.grid_size / 2), :].T, cmap='hot', origin='lower',
                   extent=[0, self.block_length, 0, self.block_height], vmin=25, vmax=80)
        plt.colorbar()
        plt.xlabel('x (m)')
        plt.ylabel('z (m)')
        plt.title('x-z Plane')
        
        plt.subplot(133)
        plt.plot(np.linspace(0, self.block_length, self.grid_size),
                 self.temperature[:, int(self.grid_size / 2), int(self.grid_size / 2)])
        plt.xlabel('x (m)')
        plt.ylabel('Temperature (°C)')
        plt.ylim(20, 80)
        
        plt.tight_layout()
        plt.show()