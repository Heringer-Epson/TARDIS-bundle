#!/usr/bin/env python

import numpy as np
from astropy import constants as const

class Input_Parameters(object):

    """
    Input parameters for quick run (<10s). Adopts homogenous ejecta.
    Used to check functionality of supporting codes.
    """
        
    def __init__(self):     

        self.input_file = __file__.split('/')[-1]
        self.subdir = 'test_fast/'
        self.filename_structure = 'test_fast'
        self.mode = 'sequential'
        self.luminosity = '9.44'
        self.time_explosion = '13.0' 
        self.velocity_start = '11000'
        self.velocity_stop = '20000'
        self.luminosity_units = 'logsolar'
        self.structure_type = 'specific'
        self.abundance_type = 'uniform' 
        self.pass_density_as = 'by_hand'
        self.ionization = 'nebular'
        self.excitation = 'dilute-lte'
        self.rad_rates_type = 'dilute-blackbody'
        self.line_interaction = 'macroatom'

        self.seeds = ['23111963']#list(np.arange(1.e4,1.0001e7,1.e4).astype(int).astype(str))
        #self.seeds = list(np.arange(1.e4,1.0001e7,1.e4).astype(int).astype(str))
        self.num_packs = '4.0e+4'
        self.iterations = '1'
        self.last_num_packs = '5.0e+4'
        self.num_virtual_packs = '5'

        self.run_uncertainties = True
        self.smoothing_window = 21
        self.N_MC_runs = 300#3000

        self.make_kromer = True
    
        self.convert_luminosity_to_logsolar()

        self.structure_num = '20'
        self.density_type = 'branch85_w7'
        self.time_0, self.rho_0, self.v_0, self.density_value   = 'blank', 'blank', 'blank', 'blank'            

        self.abun_Ni = '0.0'
        self.abun_Si = '0.52'
        self.abun_Fe = '0.0'
        self.abun_Co = '0.0'
        self.abun_Ca = '0.03' 
        self.abun_S  = '0.19'
        self.abun_Mg = '0.03'
        self.abun_Na = '0.0'
        self.abun_C  = '0.0'
        self.abun_O  = '0.18'
        self.abun_Ti = '0.01'
        self.abun_Ar = '0.04'

    def convert_luminosity_to_logsolar(self):
        """
        Ensure that luminosity is passed in consistent units to TARDIS.
        """
        L_sun = const.L_sun.cgs.value
      
        if self.luminosity_units == 'logsolar':
            pass        
        elif self.luminosity_units == 'logcgs':            
            if isinstance(self.luminosity, str):
                self.luminosity = str(format(
                  np.log10(10.**float(self.luminosity)/L_sun), '.2f')
                  )           
            elif isinstance(self.luminosity, list):
                self.luminosity = [str(format(np.log10(10.**float(lum)/L_sun), '.2f'))
                                   for lum in self.luminosity]
        elif self.luminosity_units == 'cgs':
            if isinstance(self.luminosity, str):
                self.luminosity = str(format(
                np.log10(float(self.luminosity)/L_sun), '.2f'))
            elif isinstance(self.luminosity, list):
                self.luminosity = [str(format(np.log10(float(lum)/L_sun), '.2f'))
                                   for lum in self.luminosity]
        elif self.luminosity_units == 'solar':
            if isinstance(self.luminosity, str):
                self.luminosity = str(format(
                np.log10(float(self.luminosity)), '.2f'))
            elif isinstance(self.luminosity, list):
                self.luminosity = [str(format(np.log10(float(lum)), '.2f'))
                                   for lum in self.luminosity]             
        return None
