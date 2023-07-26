# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 09:56:49 2023

@author: chonk
"""

import numpy as np

class FourierAnalysis:
    def __init__(self, x, y, period):
        """
        Initialize the FourierAnalysis object.

        Parameters
        ----------
        x : array-like
            The x-values of the data.
        y : array-like
            The y-values of the data.
        period : float
            The period of the periodic function.

        Attributes
        ----------
        x_values : array-like
            The x-values of the data.
        y_values : array-like
            The y-values of the data.
        period : float
            The period of the periodic function.
        y_FS : float
            The mean value of the y-values.
        beta : list
            List of Fourier coefficients (amplitudes).
        phi : list
            List of Fourier coefficients (phases).
        mode_number : list
            List of mode numbers corresponding to the Fourier components.
        diffusivity_transmission : list
            List of diffusivity values for transmission.
        diffusivity_phase : list
            List of diffusivity values for phase.

        """
        self.x_values = x
        self.y_values = y
        self.period = period
        self.y_FS = 0
        self.beta = []
        self.phi = []
        self.mode_number = []
        self.diffusivity_transmission = []
        self.diffusivity_phase = []

    def trapezium(self, new_x, new_y):
        """
        Calculate the area under the curve using the trapezium rule.

        Parameters
        ----------
        new_x : array-like
            The x-values of the data.
        new_y : array-like
            The y-values of the data.

        Returns
        -------
        area : float
            The calculated area under the curve.

        """
        area = 0
        for i in range(0, self.period):
            area += 0.5 * (new_x[i + 1] - new_x[i]) * (new_y[i] + new_y[i + 1])
        return area

    def fourier_components(self, modes):
        """
        Calculate the Fourier components for the given number of modes.

        Parameters
        ----------
        modes : int
            Number of Fourier modes to calculate.

        Returns
        -------
        y_FS : float
            The mean value of the y-values.
        beta : list
            List of Fourier coefficients (amplitudes).
        phi : list
            List of Fourier coefficients (phases).

        """
        self.y_FS = np.mean(self.y_values)
        self.beta = []
        self.phi = []
        for n in range(1, modes + 1):
            y_sin = np.sin((2 * np.pi * n * self.x_values) / self.period)
            y_cos = np.cos((2 * np.pi * n * self.x_values) / self.period)

            a_n = (2 / self.period) * self.trapezium(self.x_values, self.y_values * y_cos)
            b_n = (2 / self.period) * self.trapezium(self.x_values, self.y_values * y_sin)

            a_i = np.sqrt(a_n ** 2 + b_n ** 2)
            self.beta.append(a_i)
            b_i = -np.arctan2(b_n, a_n)
            self.phi.append(b_i)

            self.y_FS += a_i * np.sin((2 * np.pi * n * self.x_values) / self.period - b_i)

        return self.y_FS, self.beta, self.phi

    def diffusivity(self, l, gamma, phase):
        """
        Calculate the diffusivity values for transmission and phase.

        Parameters
        ----------
        l : float
            Length parameter.
        gamma : float
            Gamma parameter.
        phase : float
            Phase parameter.

        Returns
        -------
        D_TF : float
            Diffusivity value for transmission.
        D_PL : float
            Diffusivity value for phase.

        """
        w = 2 * np.pi / self.period

        D_TF = (w * l ** 2) / (2 * (np.log(gamma)) ** 2)
        D_PL = (w * l ** 2) / (2 * phase ** 2)

        return D_TF, D_PL

    def calculate(self, length, modes, amplitude):
        """
        Calculate diffusivity values for the given length and number of modes.

        Parameters
        ----------
        length : float
            Length parameter.
        modes : int
            Number of Fourier modes to calculate.

        Returns
        -------
        mode_number : list
            List of mode numbers corresponding to the Fourier components.
        diffusivity_transmission : list
            List of diffusivity values for transmission.
        diffusivity_phase : list
            List of diffusivity values for phase.

        """
        self.mode_number = []
        self.diffusivity_transmission = []
        self.diffusivity_phase = []
        self.amplitude = amplitude
        y_fs, beta_i, phi_i = self.fourier_components(modes)
        for i in range(1, modes + 1, 2):
            self.mode_number.append(i)
            if phi_i[i - 1] <= np.pi:
                phi = phi_i[i - 1] + 2 * np.pi
            else:
                phi = phi_i[i - 1]

            DTFi, DPLi = self.diffusivity(length, beta_i[i - 1] / (2*self.amplitude / ((i) * np.pi)), phi)
            self.diffusivity_transmission.append(DTFi)
            self.diffusivity_phase.append(DPLi)

        return self.mode_number, self.diffusivity_transmission, self.diffusivity_phase
