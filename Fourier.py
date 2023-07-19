# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 09:56:49 2023

@author: chonk
"""

import numpy as np

class FourierAnalysis:
    def __init__(self, x, y, period):
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
        area = 0
        for i in range(0, self.period):
            area += 0.5*(new_x[i+1]-new_x[i])*(new_y[i]+new_y[i+1])
        return area
    
    def fourier_components(self, modes):
        self.y_FS = np.mean(self.y)
        self.beta = []
        self.phi = []
        for n in range(1, modes+1):
            y_sin = np.sin((2*np.pi*n*self.x)/self.period)
            y_cos = np.cos((2*np.pi*n*self.x)/self.period)
        
            a_n = (2/self.period) * self.trapezium(self.x, self.y*y_cos, self.period) 
            b_n = (2/self.period) * self.trapezium(self.x, self.y*y_sin, self.period)
           
            a_i = np.sqrt(a_n**2 + b_n**2)
            self.beta.append(a_i)
            b_i = -np.arctan2(a_n,b_n)
            self.phi.append(b_i)
            
            self.y_FS += a_i * np.sin((2*np.pi*n*self.x)/self.period - b_i)
            
        return self.y_FS, self.beta, self.phi
    
    def diffusivity(self, l, gamma, phase):
        w = 2*np.pi/self.period
        
        D_TF = (w*l**2)/(2*(np.log(gamma))**2)
        D_PL = (w*l**2)/(2*phase**2)
        
        return D_TF, D_PL
    
    def calculate(self, length, modes):
        self.mode_number = []        
        self.diffusivity_transmission = []
        self.diffusivity_phase = []
        y_fs, beta_i, phi_i = self.fourier_components(modes)
        for i in range(1,modes+1, 2):
            self.mode_number.append(i)
            if phi_i[i-1] <= np.pi:
                phi = phi_i[i-1]+4*np.pi
            else:
                phi = phi_i[i-1]
                
            DTFi, DPLi = self.diffusivity(length, beta_i[i-1]/(200/((i)*np.pi)), phi)
            self.diffusivity_transmission.append(DTFi)
            self.diffusivity_phase.append(DPLi)
            
        return self.mode_number, self.diffusivity_transmission, self.diffusivity_phase