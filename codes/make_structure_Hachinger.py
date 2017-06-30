#!/usr/bin/env python

"""
THIS CODE MAKES DENSITY AND ABUNDANCE FILES ACCORDING TO THE PARAMETERS
FROM THE INPUT FILES.
"""

import os                                                               
import sys
import time
import copy
import numpy as np
import itertools                                                        
import random
import astropy.units as u
import matplotlib.pyplot as plt
  
class Structure_Hachinger(object):

    def __init__(self, mass_scaling, energy_scaling):
        
        self.ms = mass_scaling
        self.es = energy_scaling

        self.model = 'models-05bl-w7'
        if self.ms != '1.0':
            self.model = self.model + 'm' + self.ms 
        if self.es != '1.0':
            self.model = self.model + 'e' + self.es

        self.rise_time = 17. # days
        self.time = None
        self.velocity_array = None  

        for epoch in ['m6', 'm5', 'm3', 'p48', 'p129']:
            self.phase = epoch
            self.convert_phase()
            self.get_data()
            self.make_density_file()
            self.make_abundance_file()
            #self.check_decay()
        plt.legend()
        #plt.show()

    def convert_phase(self):
        converter = {'m6': -6., 'm5': -5., 'm3': -3., 'p48': 4.8, 'p129': 12.9}
        self.time = self.rise_time + converter[self.phase]
        self.time = self.time * u.day

    def get_data(self):
        
        path = ('./../INPUT_FILES/Hachinger_2005bl/' + self.model + '/' 
                + 'SN2005bl_' + self.phase + '/' + 'abuplot.dat')
                
        self.rows = []
        with open(path, 'r') as inp:
            for line in itertools.islice(inp, 2, None, 2):
                column = line.rstrip('\n').split(' ') 
                column = filter(None, column)
                self.rows.append(column)

    def make_density_file(self):
        
        self.velocity_array = np.asarray([float(row[1]) for row in self.rows])
        density_array = np.asarray([10.**float(row[3]) for row in self.rows])
                
        density_file = ('./../INPUT_FILES/DENSITY_FILES/density_05bl_es-'
                        + self.es + '_ms-' + self.ms + '_' 
                        + str(self.time.to(u.day).value) + '_day.dat')
                        
        #Density in Hachinger files are given at each time, whereas we want
        #the density at 100s, for standardization purposes.
        #density_ratio = (self.time.to(u.s).value / 100.) ** 3.
        #density_array = density_array * density_ratio
        
        plt.plot(self.velocity_array, density_array, label=self.phase)

        with open(density_file, 'w') as out_density:    
            out_density.write(str(int(self.time.to(u.s).value)) + ' s  \n')      
            out_density.write('# index velocity (km/s) density (g/cm^3)')
            
            for i, (velocity, density) in enumerate(zip(
              self.velocity_array, density_array)):
               
                out_density.write('\n' + str(i) + '    ' + str(velocity)
                                  + '    ' + str(density))                 

    def make_abundance_file(self):

        abundance_file = ('./../INPUT_FILES/STRATIFIED_COMPOSITION_FILES/'
                          + 'abundance_05bl_' + str(self.time.to(u.day).value)
                          + '_day.dat')
        
        with open(abundance_file, 'w') as out_abundance:
            out_abundance.write('# index Z=1 - Z=30')   

            for i, velocity in enumerate(self.velocity_array):
                sum_elem = 0.    
                out_abundance.write('\n' + str(i))
                
                for j in range(30):
                    out_abundance.write(' ' + str(self.rows[i][4 + j]))

                    sum_elem += float(float(self.rows[i][4 + j]))

                if abs(sum_elem - 1.) > 1.e-5:
                    raise ValueError("Error: In index %s , velocity \
                    = %s , the abundance does not add to 1. Needs %s"
                    % (i, velocity, 1-sum_elem))

    def check_decay(self):

        def decay_Ni_Co_Fe(time,Ni_formed):
            """Calculates the decay chain Ni->Co->Fe.
            
            Source: Nadyozhin+ 1994
            http://adsabs.harvard.edu/full/1994ApJS...92..527N
            
            Remarks: Difference between halflife and lifetime: the 
            former is the intuitive (1/2)^(t/Tau), while lifetime uses
            the exponential form. The correlation between constants
            is T_exp = T_1/2 *(1/ln(2))
            
            Values: Half-lifes are taken from Jonghwa Chang's
            (Korea Atomic Energy Research Institute)
            http://atom.kaeri.re.kr/nuchart/
            """
           
            Ni_half_life = 6.07/np.log(2.)
            Co_half_life = 77.24/np.log(2.)
           
            Ni_remain = Ni_formed*np.exp(-time/Ni_half_life) 
             
            Co_remain = Ni_formed*(Co_half_life/(Co_half_life
            - Ni_half_life))*(np.exp(-time/Co_half_life)
            - np.exp(-time/Ni_half_life))
            
            Fe_remain = Ni_formed*(1. + Ni_half_life/(Co_half_life
            - Ni_half_life)*np.exp(-time/Ni_half_life)
            - Co_half_life/(Co_half_life - Ni_half_life)
            *np.exp(-time/Co_half_life)
                                  )
            
            Ni_change = Ni_remain - Ni_formed
            Co_change = Co_remain
            Fe_change = Fe_remain
                        
            """Round numbers to write the abundance file"""
            Ni_change = float(format(Ni_change, '.4f'))
            Co_change = float(format(Co_change, '.4f'))
            Fe_change = float(format(Fe_change, '.4f'))           
            
            if abs(Ni_change + Co_change + Fe_change) > 1.e-5:
                #Test if rounded decay products sum to initial Ni_formed
                Ni_change -= Ni_change + Co_change + Fe_change
            
            return Ni_change, Co_change, Fe_change

        for i, velocity in enumerate(self.velocity_array):                
            print 'row ', i
       
            Ni_initial = float(self.rows[i][-1])            
            Ni_change, Co_change, Fe_change = (
              decay_Ni_Co_Fe(self.time.to(u.day).value, Ni_initial))
                            
            #print Ni_initial, Ni_change, Co_change, Fe_change
            
            Ni = str(format(float(self.rows[i][-1]) + Ni_change, '.4f'))  
            Co = str(format(float(0.) + Co_change, '.4f'))                  
            Fe = str(format(float(self.rows[i][-2]) + Fe_change, '.4f'))
            
            print '    My Ni = ', Ni, '| Hac Ni = ', self.rows[i][-5]
            print '    My Co = ', Co, '| Hac Co = ', self.rows[i][-6]
            print '    My Fe = ', Fe, '| Hac Fe = ', self.rows[i][-7]    


Structure_Hachinger(mass_scaling='1.0', energy_scaling='0.7')


