# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 10:09:36 2023

@author: chonk
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from Calculations import HeatEquationSolver
from Input_Heat_Waves import HeatingFunction


class ThermalSimulation:
    """
    Class for thermal simulation of a block.
    """

    def __init__(self, block_length, block_width, block_height, total_time, grid_size,
                 density, specific_heat, initial_conductivity, heating_function):
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
        self.heating_function = heating_function
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

        self.is_hole = np.zeros(
            (self.grid_size, self.grid_size, self.grid_size), dtype=bool)

    def drill_holes(self, hole_positions, hole_depth, hole_width, face):
        """
        Method to represent square holes in the material, considering the face from which the holes are drilled.
    
        Parameters
        ----------
        hole_positions (list of floats): The positions of the hole centers along the drilling direction (x, y, z) in meters.
        hole_depth (float): The depth of the holes in meters.
        hole_width (float): The width of the holes in meters (considered in the plane orthogonal to drilling direction).
        face (str): The face from which the drilling starts. Could be '+x', '-x', '+y', '-y', '+z', '-z'.
        """
    
        # Mapping from face to drill direction
        face_to_direction = {
            "+x": (1, 0, 0),
            "-x": (-1, 0, 0),
            "+y": (0, 1, 0),
            "-y": (0, -1, 0),
            "+z": (0, 0, 1),
            "-z": (0, 0, -1),
        }
    
        drill_direction = face_to_direction[face]
    
        self.is_hole = np.full_like(self.temperature, fill_value=False)
    
        for hole_position in hole_positions:
            hole_center_indices = np.round(
                np.array(hole_position) / np.array([self.dx, self.dy, self.dz])).astype(int)
            hole_depth_indices = int(hole_depth / max(self.dx, self.dy, self.dz))
            hole_width_indices = int(hole_width / max(self.dx, self.dy, self.dz))
    
            x_indices = np.arange(self.grid_size)
            y_indices = np.arange(self.grid_size)
            z_indices = np.arange(self.grid_size)
    
            X, Y, Z = np.meshgrid(x_indices, y_indices, z_indices, indexing='ij')
    
            rel_X = X - hole_center_indices[0]
            rel_Y = Y - hole_center_indices[1]
            rel_Z = Z - hole_center_indices[2]
    
            distance_from_center_along_drill_direction = rel_X * \
                drill_direction[0] + rel_Y * \
                drill_direction[1] + rel_Z*drill_direction[2]
    
            depth_mask = (0 <= distance_from_center_along_drill_direction) & (
                distance_from_center_along_drill_direction <= hole_depth_indices)
    
            distance_from_drill_line_X = rel_X - \
                distance_from_center_along_drill_direction*drill_direction[0]
            distance_from_drill_line_Y = rel_Y - \
                distance_from_center_along_drill_direction*drill_direction[1]
            distance_from_drill_line_Z = rel_Z - \
                distance_from_center_along_drill_direction*drill_direction[2]
    
            width_mask = (np.abs(distance_from_drill_line_X) <= hole_width_indices/2) & (np.abs(distance_from_drill_line_Y)
                                                                                         <= hole_width_indices/2) & (np.abs(distance_from_drill_line_Z) <= hole_width_indices/2)
    
            self.is_hole = np.logical_or(self.is_hole, depth_mask & width_mask)
    
        self.initial_conductivity[self.is_hole] = 0.026  # thermal conductivity of air in W/m/K
        self.initial_density[self.is_hole] = 1.2041  # density of air in kg/m^3 at 20 degrees Celsius


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
                                    self.heating_function, self.specific_heat, self.initial_conductivity, self.initial_density)

        for step in range(self.num_steps):
            self.temperature = solver.solve(self.time[step])

            if step % 100 == 0:
                self.plot_temperature(step)
                plt.pause(0.05)
                self.T_quarter.append(self.temperature[int(
                    0.25*self.grid_size), int(self.grid_size / 2), int(self.grid_size / 2)])
                self.T_half.append(self.temperature[int(
                    0.5*self.grid_size), int(self.grid_size / 2), int(self.grid_size / 2)])
                self.T_three_quarter.append(self.temperature[int(
                    0.75*self.grid_size), int(self.grid_size / 2), int(self.grid_size / 2)])
                self.T_end.append(
                    self.temperature[-1, int(self.grid_size / 2), int(self.grid_size / 2)])
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
        plt.imshow(self.temperature[:, :, int(self.grid_size / 2)+1].T, cmap='hot', origin='lower',
                   extent=[0, self.block_length, 0, self.block_width], vmin=25, vmax=80)
        plt.colorbar()
        plt.xlabel('x (m)')
        plt.ylabel('y (m)')
        plt.title('x-y Plane')

        plt.subplot(132)
        plt.imshow(self.temperature[:, int(self.grid_size / 2)-1, :].T, cmap='hot', origin='lower',
                   extent=[0, self.block_length, 0, self.block_height], vmin=25, vmax=80)
        plt.colorbar()
        plt.xlabel('x (m)')
        plt.ylabel('z (m)')
        plt.title('x-z Plane')
        
        plt.subplot(133)
        plt.plot(np.linspace(0, self.block_length, self.grid_size),
                 self.temperature[:, int(self.grid_size / 2)-1, int(self.grid_size / 2)]+1)
        plt.xlabel('x (m)')
        plt.ylabel('Temperature (°C)')
        plt.ylim(20, 80)
        
        plt.tight_layout()
        plt.show()
        
    def visualize_holes(self):
        """
        Function to visualize the holes in the block.
        """
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        x, y, z = np.where(self.is_hole == True)
        ax.scatter(x, y, z)
        plt.show()

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
