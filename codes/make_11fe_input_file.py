#!/usr/bin/env python

import os                                                               
import sys
import numpy as np
from scipy.interpolate import interp1d
from astropy import constants as const

M_sun = const.M_sun.to('g').value

class Make_11fe_file(object):
    """
    Code Description
    ----------

    TBW
    """
    
    def __init__(self):
        
        self.velocity_requested = (
          list(np.logspace(np.log10(3550.), np.log10(40000.), 100).astype(str)))   

        self.D = {}
        self.out_abun = {}
        self.density_fine = []
        self.velocity_fine = []
        self.mass_fine = None

        self.density_coarse = []
        self.velocity_coarse = []
        self.mass_coarse = []
        
        self.N_coarse_shells = None
        
        self.top_dir = './../INPUT_FILES/Mazzali_2011fe/'
        
        self.el = ['C', 'O', 'Mg', 'Si', 'S', 'Ca', 'TiCr', 'Fe0', 'Ni0']
        
        self.el_all = [
          'H', 'He', 'Li', 'B', 'Be', 'C', 'N', 'O', 'F', 'Ne', 'Na',
          'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc',
          'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Fe0', 'Ni0'] 
        
        self.run_make()

    def initialize_dict(self):
        for el in self.el:
            self.D[el] = []
        
    def load_density(self):
        fpath = self.top_dir + 'density.csv'
        
        with open(fpath, 'r') as inp:
            for line in inp:
                column = line.rstrip('\n').split(',') 
                column = filter(None, column)
                self.velocity_fine.append(float(column[0]))                        
                self.density_fine.append(float(column[1]))                        

        self.velocity_fine = np.asarray(self.velocity_fine)
        self.density_fine = np.asarray(self.density_fine)
        
    def compute_mass(self):
        time = 100.
        km2cm = 1.e5

        volume_array = (4./3. * np.pi * time**3. * km2cm**3. *
                        (np.power(self.velocity_fine[1:], 3)
                        - np.power(self.velocity_fine[0:-1], 3)))
        
        avg_density = (self.density_fine[0:-1] + self.density_fine[1:]) / 2.
        
        self.mass_fine = np.multiply(avg_density, volume_array)
        self.mass_fine = np.cumsum(self.mass_fine) / M_sun
        
    def interpolate_abun(self, Mass, Abun):
        
        #If the last value is smaller than 0.02% (0.0002), then it should 
        #actually be zero, see original plot in Mazzali's paper.
        if Abun[-1] < 0.02:
            Abun[-1] = 0.
            
        #Check whether first value should be zero. This is important for C and Mg.
        if Mass[0] > 0.4:
            Abun[0] = 0.
            
        #For the interpolation to always work, extrapolate the first and last
        #values to masses 0 and 1.5,m by conserving the extreme abundances.
        Mass = np.asarray([0.] + Mass + [1.5])
        Abun = np.asarray([Abun[0]] + Abun + [Abun[-1]])
        
        mass2abun = interp1d(Mass, Abun)
        return mass2abun 

    def load_abundances(self):
        """ For the Ni curve, the points were there seemed to be a change in
        slope of the abundance track were recorded by hand. These points will
        provide the velocity (and consequently mass at were one needs to
        retrieve the abundances of the other elements.)
        """
        fpath = self.top_dir + '/abundance_Ni0.csv'

        with open(fpath, 'r') as inp:
            for line in inp:
                column = line.rstrip('\n').split(',') 
                column = filter(None, column)
                self.mass_coarse.append(float(column[0]))            
        
        self.mass_coarse = np.asarray(self.mass_coarse)
               
        for el in self.el:
            fpath = self.top_dir + '/abundance_' + el + '.csv'
            mass = []
            abun = []

            with open(fpath, 'r') as inp:
                for line in inp:
                    column = line.rstrip('\n').split(',') 
                    column = filter(None, column)
                    mass.append(float(column[0]))
                    abun.append(float(column[1]))                        

            mass2abun = self.interpolate_abun(mass, abun)

            self.D[el] = mass2abun(self.mass_coarse) / 100.
        self.N_coarse_shells = len(self.mass_coarse)    

    def check_abun_normalization(self):
        
        total_per_shell = []
                
        for i in range(self.N_coarse_shells):
            
            total = 0.
            for el in self.el:
                total += self.D[el][i]
            
            if total < 0.8:
                self.D['O'][i] += 1. - total
                
            else:
                for el in self.el:
                    self.D[el][i] *= 1. / total                            

    def get_velocity_coarse(self):
        avg_velocity = (self.velocity_fine[0:-1] + self.velocity_fine[1:]) / 2.
        mass2velocity = interp1d(self.mass_fine, avg_velocity)
        self.velocity_coarse = mass2velocity(self.mass_coarse)       

    def get_density_coarse(self):
        avg_density = (self.density_fine[0:-1] + self.density_fine[1:]) / 2.
        mass2density = interp1d(self.mass_fine, avg_density)
        self.density_coarse = mass2density(self.mass_coarse)       

    def make_out_dict(self):
        
        zero_array = np.zeros(self.N_coarse_shells)
        nan_array = np.zeros(self.N_coarse_shells)
        nan_array[:] = np.nan

        other_el = ['H', 'He', 'Li', 'B', 'Be', 'N', 'F', 'Ne', 'Na', 'Al',
                    'P', 'Cl', 'Ar', 'K', 'Sc', 'V', 'Mn', 'Co', 'Cu', 'Zn'] 
        
        for el in other_el:
            self.out_abun[el] = zero_array

        for el in ['C', 'O', 'Mg', 'Si', 'S', 'Ca', 'Fe0', 'Ni0']:
            self.out_abun[el] = self.D[el]
            
        self.out_abun['Ti'] = self.D['TiCr'] / 2.    
        self.out_abun['Cr'] = self.D['TiCr'] / 2. 
        self.out_abun['Fe'] = nan_array
        self.out_abun['Ni'] = nan_array
        
    def make_output(self):
        
        fpath = self.top_dir + '/ejecta_layers.dat'
        
        #Make header
        line1 = 'shno velocity cum_M logdens'
        line2 = 'shno velocity cum_M logdens'
        for i, el in enumerate(self.el_all):
            el_num = i + 1
            if el_num == 31:
                el_num = 98
            if el_num == 32:
                el_num = 99
                            
            line1 += ' ' + str(el_num) 
            line2 += ' ' + el 

        line1 += ' \n'
        line2 += ' \n'

        #Write outfile.      
        with open(fpath, 'w') as out:
            out.write(line1) 
            out.write(line2) 
            for i in range(self.N_coarse_shells):            
                shno = str(i + 1)
                vel = str(format(self.velocity_coarse[i], '.0f'))
                cum_m = str(format(self.mass_coarse[i], '.3f'))
                logdens = str(format(np.log10(self.density_coarse[i]), '.5f'))
                out.write(' ' + shno + ' ' + vel + ' ' + cum_m + ' ' + logdens)
                for el in self.el_all:
                    out.write(' ' + str(format(self.out_abun[el][i], '.6f')))
                out.write(' \n')    


    def run_make(self):
        self.initialize_dict()
        self.load_density()
        self.compute_mass()
        self.load_abundances()
        self.check_abun_normalization()
        self.get_velocity_coarse()
        self.get_density_coarse()
        self.make_out_dict()
        self.make_output()
        
        
if __name__ == '__main__':
    Make_11fe_file()

