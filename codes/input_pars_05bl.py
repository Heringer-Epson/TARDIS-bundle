#!/usr/bin/env python

import numpy as np
from astropy import constants as const
from density_interpolator import W7_velocity_to_density

class Input_Parameters(object):

    """
    Input parameters to reproduce the tomography of 2005bl according
    to my interpretation of the data in Hachinger+ 2009
    [[http://adsabs.harvard.edu/abs/2009MNRAS.399.1238H]].
    """

    def __init__(self):     

        self.input_file = __file__.split('/')[-1]
        self.subdir = '05bl_standard_downbranch/'
        self.filename_structure = '05bl_standard_downbranch'
        self.mode = 'sequential'

        """
        Below, several pre-defined options for simulation 05bl ejecta.
        Uncomment preferred set. 
        """
        
        """
        Default
        """
        #self.luminosity = ['3.35e8', '4.24e8', '5.68e8', '7.73e8', '3.93e8']
        #self.time_explosion = ['11.0', '12.0', '14.0', '21.8', '29.9']
        #self.velocity_start = ['8400', '8100', '7500', '6600', '3300']
        
        """
        Used to luminosity variations at all epochs
        """
        self.luminosity = [
          '3.35e8', '4.24e8', '5.68e8', '7.73e8', '3.93e8',
          '6.70e8', '8.48e8', '1.14e9', '1.55e9', '7.86e8',
          '1.01e9', '1.27e9', '1.70e9', '2.32e9', '1.02e9',
          '1.34e9', '1.70e9', '2.24e9', '3.01e9', '1.57e9'
          ]
        self.time_explosion = ['11.0', '12.0', '14.0', '21.8', '29.9']*4
        self.velocity_start = ['8400', '8100', '7500', '6600', '3300']*4

        """
        Used to compute L-grid
        """
        #self.luminosity = list(np.logspace(8.544, 9.72, 20).astype(str))
        #self.time_explosion = '21.8'
        #self.velocity_start = '6600'

        """
        Below, default parameters that usually need not be changed.
        """

        self.velocity_stop = '17500'
        self.luminosity_units = 'solar'
        self.structure_type = 'file' #specific #file
        self.abundance_type = 'file' #file,uniform,branch_w7
        self.pass_density_as = 'by_hand'

        self.ionization = 'nebular'
        self.excitation = 'dilute-lte'
        self.rad_rates_type = 'dilute-blackbody'
        self.line_interaction = 'downbranch'

        #For high S/N runs
        self.seeds = '23111963'
        self.num_packs = '2.0e+5'
        self.iterations = '15'
        self.last_num_packs = '5.0e+5'
        self.num_virtual_packs = '5'

        #For faster runs
        #self.num_packs = '1.0e+5'
        #self.iterations = '10'
        #self.last_num_packs = '1.0e+5'
        #self.num_virtual_packs = '5'

        self.run_uncertainties = True
        self.smoothing_window = 21
        self.N_MC_runs = 3000

        self.make_kromer = False
        
        self.convert_luminosity_to_logsolar()
        
        self.mass_ratio_scaling = '1.0'
        self.energy_ratio_scaling = '1.0'

        #Velocity grid used for the simulation.
        #Note, this grid can be finer than the self.velocity_zones defined 
        #below. In which case, the corresponding abundance will be that of the
        #closest velocity zone **interior** to the requested velocity step.
        self.velocity_array = (list(np.logspace(np.log10(3500.),
                               np.log10(17000.), 100).astype(str)))
                               
        self.density_array = W7_velocity_to_density(self.velocity_array)
        self.time_0 = '100 s'
        self.rho_0, self.v_0 = None, None

        #-Do not modify the following four lines-#
        self.abun_raw, self.abun = {}, {}
        self.elements = [
          'H', 'He', 'Li', 'B', 'Be', 'C', 'N', 'O', 'F', 'Ne', 'Na',
          'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc',
          'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn'
          ] 

        """
        Note:
        Hachinger+ 2009
        [[http://adsabs.harvard.edu/abs/2009MNRAS.399.1238H]]
        paper does not especify how high high-velocity is.
        Here assumed to be 17,500 km/s. Values for other layers come
        from their table A1 for the W7 model. Needs adding finer layers.
        """            
        self.velocity_zones = ['3500',     '4000',   '4500',   '5000',   '5500',   '6000',   '6500',   '7000',   '7500',   '8000',   '8500',   '9000', '10000',  '11000',  '12000',  '13000',  '14000',  '15000',  '15500',  '16000',  '17000']
        
        self.abun_raw['H']  = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['He'] = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Li'] = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Be'] = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['B']  = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['C']  = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0500', '0.0500', '0.0500', '0.0500', '0.0500', '0.0500', '0.0500', '0.4184', '0.4184', '0.4184']
        self.abun_raw['N']  = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['O']  = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0100', '0.1700', '0.8431', '0.8431', '0.8431', '0.8431', '0.8431', '0.8431', '0.8431', '0.5726', '0.5726', '0.5726']
        self.abun_raw['F']  = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Ne'] = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Na'] = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0013', '0.0025', '0.0025', '0.0025', '0.0025', '0.0025', '0.0025', '0.0025', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Mg'] = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0400', '0.0300', '0.0300', '0.0300', '0.0300', '0.0300', '0.0300', '0.0300', '0.0030', '0.0030', '0.0030']
        self.abun_raw['Al'] = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0011', '0.0032', '0.0032', '0.0032', '0.0032', '0.0032', '0.0032', '0.0032', '0.0032', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Si'] = ['0.6447', '0.6447', '0.6447', '0.6447', '0.6447', '0.6447', '0.6447', '0.6947', '0.6947', '0.6886', '0.6477', '0.0500', '0.0500', '0.0500', '0.0500', '0.0500', '0.0500', '0.0500', '0.0050', '0.0050', '0.0050']
        self.abun_raw['P']  = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['S']  = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0600', '0.0600', '0.1200', '0.1000', '0.0200', '0.0200', '0.0200', '0.0200', '0.0200', '0.0200', '0.0200', '0.0010', '0.0010', '0.0010']
        self.abun_raw['Cl'] = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Ar'] = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['K']  = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Ca'] = ['0.0003', '0.0003', '0.0003', '0.0003', '0.0003', '0.0003', '0.0003', '0.0003', '0.0003', '0.0003', '0.0003', '0.0003', '0.0003', '0.0003', '0.0003', '0.0003', '0.0003', '0.0003', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Sc'] = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Ti'] = ['0.0350', '0.0350', '0.0350', '0.0350', '0.0350', '0.0350', '0.0350', '0.0675', '0.0675', '0.0500', '0.0100', '0.0004', '0.0004', '0.0004', '0.0004', '0.0004', '0.0004', '0.0004', '0.0000', '0.0000', '0.0000']
        self.abun_raw['V']  = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Cr'] = ['0.0350', '0.0350', '0.0350', '0.0350', '0.0350', '0.0350', '0.0350', '0.0675', '0.0675', '0.0500', '0.0100', '0.0004', '0.0004', '0.0004', '0.0004', '0.0004', '0.0004', '0.0004', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Mn'] = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Fe'] = ['0.0950', '0.0950', '0.0950', '0.0950', '0.0950', '0.0950', '0.0950', '0.1000', '0.1000', '0.0800', '0.0175', '0.0001', '0.0001', '0.0001', '0.0001', '0.0001', '0.0001', '0.0001', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Co'] = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Ni'] = ['0.1900', '0.1900', '0.1900', '0.1900', '0.1900', '0.1900', '0.1900', '0.0100', '0.0100', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Cu'] = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Zn'] = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']

        """
        Format abundances.
        """
        #The velocity array can be finer than the zones. To each velocity in
        #the array, find the closest zone below and attribute the mass fraction
        #of that zone. Use this procedure for all elements.
        v_zones = np.asarray(self.velocity_zones).astype(np.float)
        for element in self.elements:
            self.abun[element] = []
            for v in self.velocity_array:
                v = float(v)
                condition = (v_zones <= v)
                idx_zone = v_zones[condition].argmax()
                self.abun[element].append(self.abun_raw[element][idx_zone])

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
