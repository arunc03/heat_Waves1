# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 10:08:44 2023

@author: chonk
"""

import numpy as np
from scipy.signal import square

class HeatingFunction:
    def __init__(self, peak, room_temperature, frequency):
        self.peak = peak
        self.room_temperature = room_temperature
        self.frequency = frequency

    def constant_temperature(self, time, amplitude=100):
        """
        Function to simulate constant temperature.

        Parameters
        ----------
        time (float): The current time in seconds.
        peak (float): The constant temperature in Celsius.

        Returns
        -------
        float
            The constant temperature.
        """
        return amplitude
    
    def sinusoidal_temperature(self, time, frequency=1, peak=100, offset=0):
        """
        Function to simulate sinusoidal temperature.

        Parameters
        ----------
        time (float): The current time in seconds.
        frequency (float): The frequency of the sinusoidal wave.
        peak (float): The peak of the sinusoidal wave in Celsius.
        offset (float): The offset of the sinusoidal wave in Celsius.

        Returns
        -------
        float
            The sinusoidal wave temperature.
        """
        
        return np.sin(2 * np.pi * self.frequency * time) * ((self.peak - self.room_temperature) / 2) +  self.room_temperature + ((self.peak - self.room_temperature) / 2)

    def square_wave_temperature(self, time, frequency=1, amplitude=100, offset=0):
        """
        Function to simulate square wave temperature.

        Parameters
        ----------
        time (float): The current time in seconds.
        frequency (float): The frequency of the square wave.
        peak (float): The peak of the square wave in Celsius.
        offset (float): The offset of the square wave in Celsius.

        Returns
        -------
        float
            The square wave temperature.
        """
        return (square(2 * np.pi * self.frequency * time) + 1) / 2 * (self.peak - self.room_temperature) + self.room_temperature    