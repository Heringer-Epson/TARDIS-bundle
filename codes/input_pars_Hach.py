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
        #Standard model is 'w7e0.7' from Hachinger files and
        #filename_structure should not be changed from 05bl.
        self.subdir = '05bl_Lgrid_macroatom/'
        self.filename_structure = '05bl'
        self.extinction = -0.20

        """
        Below, several pre-defined options for simulation 05bl ejecta.
        Uncomment preferred set. 
        """
        
        """
        Default
        """
        #self.luminosity = [str(l) for l in [8.520, 8.617, 8.745, 8.861, 8.594]]
        #self.time_explosion = ['11.0', '12.0', '14.0', '21.8', '29.9']
        #self.velocity_start = ['8350', '8100', '7600', '6800', '3350']
        
        """
        Used to luminosity variations at all epochs
        """
        #self.luminosity = []
        #lum = [8.617, 8.861, 8.594]
        #for scale in [1., 2., 3., 4.]:
        #    self.luminosity += [str(format(np.log10(10.**l * scale), '.3f')) for l in list(lum)]
        #self.time_explosion = ['12.0', '21.8', '29.9']*4
        #self.velocity_start = ['8100', '6800', '3350']*4

        """
        Used to compute L-grid
        """
        self.luminosity = [str(format(l, '.3f')) for l in np.log10(np.logspace(8.544, 9.72, 20))]
        self.time_explosion = '21.8'
        self.velocity_start = '6800'

        """
        Below, default parameters that usually need not be changed.
        """

        self.velocity_stop = '48000'
        self.luminosity_units = 'logsolar'
        self.structure_type = 'file' #specific #file
        self.abundance_type = 'file' #file,uniform,branch_w7
        self.pass_density_as = 'by_hand'

        self.ionization = 'nebular'
        self.excitation = 'dilute-lte'
        self.rad_rates_type = 'dilute-blackbody'
        self.line_interaction = 'macroatom'

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
        
        self.velocity_array = None                
        self.mass_ratio_scaling = '1.0'
        self.energy_ratio_scaling = '0.7'

          
