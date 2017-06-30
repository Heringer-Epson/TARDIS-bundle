#!/usr/bin/env python

import numpy as np
from astropy import constants as const
from density_interpolator import rho11fe_velocity_to_density

class Input_Parameters(object):

    """
    Input parameters to reproduce the tomography of 2011fe according
    to my interpretation of the data in Mazzali+ 2014
    [[http://adsabs.harvard.edu/abs/2014MNRAS.439.1959M]].
    """
        
    def __init__(self):     

        self.input_file = __file__.split('/')[-1]
        self.subdir = '11fe_testiii/'
        self.filename_structure = '11fe_test'      
        self.extinction = -0.014

        """
        Temp - for a quick test
        """
        self.luminosity = [str(format(np.log10(3.5e7), '.3f'))]
        self.time_explosion = ['19.1'] 
        self.velocity_start = ['18000']
        
        """
        Default
        """
        #lum_solar = [0.08e9, 0.32e9, 1.1e9, 2.3e9, 3.2e9, 3.5e9, 3.2e9, 2.3e9]
        #self.luminosity  = [str(np.log10(l)) for l in lum_solar]        
        #self.time_explosion = ['3.7', '5.9', '9.0', '12.1', '16.1', '19.1', '22.4', '28.3'] 
        #self.velocity_start = ['13300', '12400', '11300', '10700', '9000', '7850', '6700', '4550']

        """
        Kromer
        """
        #self.luminosity = ['2.3e9', '1.15e9']
        #self.time_explosion = ['12.1', '12.1'] 
        #self.velocity_start = ['10700', '10700']

        """
        Used to luminosity variations at all epochs
        """
        #self.luminosity = []
        #lum = [np.log10(l) for l in [0.08e9, 0.32e9, 1.1e9, 2.3e9, 3.2e9, 3.5e9, 3.2e9, 2.3e9]]
        #for scale in [1., 0.5, 0.33, 0.25]:
        #    self.luminosity += [str(format(np.log10(10.**l * scale), '.3f')) for l in list(lum)]        
        #self.time_explosion = ['3.7', '5.9', '9.0', '12.1', '16.1', '19.1', '22.4', '28.3']*4
        #self.velocity_start = ['13300', '12400', '11300', '10700', '9000', '7850', '6700', '4550']*4      

        """
        Used to compute L-grid
        """
        #self.luminosity = [str(format(l, '.3f')) for l in np.log10(np.logspace(8.544, 9.72, 20))]
        #self.time_explosion = '19.1'    
        #self.velocity_start = '7850'

        """
        Change photosphere position at pre-maximum for 05bl Lum.
        """
        #self.luminosity =    str(format(np.log10(4.13e8), '.3f')) #'0.77e9' #'2.3e9' #'2.3e9' #
        #self.time_explosion = '12.1'  
        #self.velocity_start = ['5000', '6000', '8100', '10700', '12200']

        """
        Change time explosion at pre-maximum for 05bl Lum.
        """
        #self.luminosity =    '8.617' #'0.77e9' #'2.3e9' #'2.3e9' #
        #self.time_explosion = ['12.1']#['6', '7', '8', '9', '12.1', '13', '14']   
        #self.velocity_start = '10700'

        """
        Change time explosion at post-maximum for 05bl Lum.
        """
        #self.luminosity =    '8.594' #'0.77e9' #'2.3e9' #'2.3e9' #
        #self.time_explosion = ['23', '25', '28.3', '30', '32', '34']   
        #self.velocity_start = '4550'

        """
        Change photosphere position at maximum for default and quarter L.
        """
        #self.luminosity = '0.88e9' #'3.5e9' 
        #self.time_explosion = '19.1'   
        #self.velocity_start = ['6700', '6900', '7100', '7300', '7500', '7700',  '7850', '8000', '8200',  '8400',  '8600',  '8800']

        """
        Change photosphere position at post-maximum for default, 0.33 and L.
        """
        #self.luminosity =   '2.3e9' #'3.93e8' #'0.77e9'
        #self.time_explosion = '28.3'   
        #self.velocity_start = ['3000', '3300', '3600', '3900', '4200', '4550', '4800', '5300', '5800', '6500']


        """
        Used to ejecta with scaled Titanium and Chromium
        """
        #self.luminosity = []
        #lum = [np.log10(l) for l in [3.5e9, 3.5e9 / 4.]]
        #self.luminosity += [str(format(np.log10(10.**l), '.3f')) for l in list(lum)]        
        #self.time_explosion = '19.1'#['3.7', '5.9', '9.0', '12.1', '16.1', '19.1', '22.4', '28.3']*4    
        #self.velocity_start = '7850'#['13300', '12400', '11300', '10700', '9000', '7850', '6700', '4550']*4
        """"""

        #Velocity stop does not matter, as long as it is above the highest
        #velocity in the velocity_array.
        self.velocity_stop = '40000'
        
        self.luminosity_units = 'logsolar'
        self.structure_type = 'file' #specific #file
        self.abundance_type = 'file' #file,uniform,branch_w7
        self.pass_density_as = 'by_hand'

        self.ionization = 'nebular'
        self.excitation = 'dilute-lte'
        self.rad_rates_type = 'dilute-blackbody'
        self.line_interaction = 'downbranch'

        #For high S/N runs
        #self.seeds = '23111963'
        #self.num_packs = '2.0e+5'
        #self.iterations = '20'
        #self.last_num_packs = '5.0e+5'
        #self.num_virtual_packs = '5'

        #For faster runs
        #self.seeds = '23111963'
        #self.num_packs = '1.0e+5'
        #self.iterations = '15'
        #self.last_num_packs = '1.0e+5'
        #self.num_virtual_packs = '5'

        #For quite faster runs
        self.seeds = '23111963'
        self.num_packs = '1.0e+4'
        self.iterations = '15'
        self.last_num_packs = '5.0e+4'
        self.num_virtual_packs = '5'

        self.run_uncertainties = False
        self.smoothing_window = 21
        self.N_MC_runs = 3000

        self.make_kromer = False
        
        self.TiCr_scaling = 1.00
        self.Fe_scaling = 1.00
        self.convert_luminosity_to_logsolar()
        
        self.mass_ratio_scaling = '1.0'
        self.energy_ratio_scaling = '1.0'
        
        #Velocity grid used for the simulation.
        #Note, this grid can be finer than the self.velocity_zones defined 
        #below. In which case, the corresponding abundance will be that of the
        #closest velocity zone **interior** to the requested velocity step.
        self.velocity_array = (list(np.logspace(np.log10(3550.),
                               np.log10(24000.), 100).astype(str)))       
        self.pass_density_as = 'by_hand'#'from_exponential'# or 'by_hand'
        self.density_array = rho11fe_velocity_to_density(self.velocity_array)        
        
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
        Trying to approach Mazzali paper
        self.velocity_zones = ['3500',   '6700',   '7500',   '7850',   '8500',   '9000',   '10700',  '11300', '13300',  '16000',  '19500']

        self.abun_raw['H']  = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['He'] = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Li'] = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Be'] = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['B']  = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['C']  = ['0.0000', '0.0000', '0.0000', '0.0080', '0.0080', '0.0080', '0.0070', '0.000', '0.0310', '0.0250', '0.9804']
        self.abun_raw['N']  = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['O']  = ['0.0000', '0.0000', '0.0000', '0.0200', '0.0900', '0.0900', '0.0800', '0.241', '0.6990', '0.8600', '0.0120']
        self.abun_raw['F']  = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Ne'] = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Na'] = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Mg'] = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.020', '0.0300', '0.0300', '0.0030']
        self.abun_raw['Al'] = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Si'] = ['0.0900', '0.2000', '0.2165', '0.2285', '0.1985', '0.4785', '0.4800', '0.480', '0.2000', '0.0600', '0.0040']
        self.abun_raw['P']  = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['S']  = ['0.0300', '0.0600', '0.0700', '0.0700', '0.0700', '0.1500', '0.1500', '0.150', '0.0300', '0.0200', '0.0005']
        self.abun_raw['Cl'] = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Ar'] = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['K']  = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Ca'] = ['0.0001', '0.0030', '0.0030', '0.0030', '0.0030', '0.0030', '0.0030', '0.003', '0.0030', '0.0020', '0.0000']
        self.abun_raw['Sc'] = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Ti'] = ['0.0010', '0.0050', '0.0050', '0.0050', '0.0050', '0.0050', '0.0050', '0.003', '0.0005', '0.0000', '0.0000']
        self.abun_raw['V']  = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Cr'] = ['0.0010', '0.0050', '0.0050', '0.0050', '0.0050', '0.0050', '0.0050', '0.003', '0.0005', '0.0000', '0.0000']
        self.abun_raw['Mn'] = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Fe'] = ['0.1500', '0.0005', '0.0005', '0.0005', '0.0005', '0.0005', '0.0600', '0.060', '0.0040', '0.0020', '0.0001']
        self.abun_raw['Co'] = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Ni'] = ['0.7279', '0.7265', '0.7000', '0.6600', '0.6200', '0.2600', '0.2100', '0.040', '0.0020', '0.0010', '0.0000']
        self.abun_raw['Cu'] = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Zn'] = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.0000', '0.0000', '0.0000']
        """

        """Former attempt to model MAzzali's paper that seems to work better to
        reproduce 05bl, but the boundary of some zone are slightly different
        than in the text of the paper."""
        #self.velocity_start = ['13300', '12400', '11300', '10700', '9000', '7850', '6700', '4550']
        self.velocity_zones = ['3500',   '7000',   '7500',   '8000',   '8500',   '9000',   '11000', '12000', '13500',  '16000',  '19500']
        
        self.abun_raw['H']  = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['He'] = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Li'] = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Be'] = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['B']  = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['C']  = ['0.0000', '0.0000', '0.0000', '0.0080', '0.0080', '0.0080', '0.008', '0.000', '0.0310', '0.0260', '0.9804']
        self.abun_raw['N']  = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['O']  = ['0.0000', '0.0000', '0.0000', '0.0200', '0.0900', '0.0900', '0.110', '0.351', '0.7030', '0.8605', '0.0120']
        self.abun_raw['F']  = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Ne'] = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Na'] = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Mg'] = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.020', '0.0300', '0.0300', '0.0030']
        self.abun_raw['Al'] = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Si'] = ['0.0900', '0.2000', '0.2165', '0.2285', '0.1985', '0.4785', '0.563', '0.440', '0.2000', '0.0600', '0.0040']
        self.abun_raw['P']  = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['S']  = ['0.0300', '0.0600', '0.0700', '0.0700', '0.0700', '0.1500', '0.150', '0.080', '0.0300', '0.0200', '0.0005']
        self.abun_raw['Cl'] = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Ar'] = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['K']  = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Ca'] = ['0.0001', '0.0030', '0.0030', '0.0030', '0.0030', '0.0030', '0.003', '0.003', '0.0030', '0.0020', '0.0000']
        self.abun_raw['Sc'] = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Ti'] = ['0.0010', '0.0050', '0.0050', '0.0050', '0.0050', '0.0050', '0.005', '0.003', '0.0005', '0.0000', '0.0000']
        self.abun_raw['V']  = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Cr'] = ['0.0010', '0.0050', '0.0050', '0.0050', '0.0050', '0.0050', '0.005', '0.003', '0.0005', '0.0000', '0.0000']
        self.abun_raw['Mn'] = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Fe'] = ['0.1500', '0.0005', '0.0005', '0.0005', '0.0005', '0.0005', '0.006', '0.060', '0.0010', '0.0005', '0.0001']
        self.abun_raw['Co'] = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Ni'] = ['0.7279', '0.7265', '0.7000', '0.6600', '0.6200', '0.2600', '0.150', '0.040', '0.0010', '0.0010', '0.0000']
        self.abun_raw['Cu'] = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.000', '0.0000', '0.0000', '0.0000']
        self.abun_raw['Zn'] = ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.000', '0.0000', '0.0000', '0.0000']
  
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
                
        """
        Scale Ti/Cr at the expense of most abundant element.
        """
        for i in range(len(self.abun['H'])):
            abundance_all = np.asarray([self.abun[el][i] for el in self.elements])
            el_most = self.elements[abundance_all.argmax()]
            orig_Ti = float(self.abun['Ti'][i])    
            orig_Cr = float(self.abun['Cr'][i])
            self.abun['Ti'][i] = str(self.TiCr_scaling*float(self.abun['Ti'][i]))       
            self.abun['Cr'][i] = str(self.TiCr_scaling*float(self.abun['Cr'][i]))             
            self.abun[el_most][i] = str( float(self.abun[el_most][i])
              + (orig_Ti - float(self.abun['Ti'][i]))
              + (orig_Cr - float(self.abun['Cr'][i])))

        """
        Scale Fe at the expense of most abundant element.
        """
        for i in range(len(self.abun['H'])):
            abundance_all = np.asarray([self.abun[el][i] for el in self.elements])
            el_most = self.elements[abundance_all.argmax()]
            orig_Fe = float(self.abun['Fe'][i])
            self.abun['Fe'][i] = str(self.Fe_scaling*float(self.abun['Fe'][i]))            
            self.abun[el_most][i] = str( float(self.abun[el_most][i])
            + (orig_Fe - float(self.abun['Fe'][i])))

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
                  np.log10(10.**float(self.luminosity)/L_sun), '.2f'))           
            elif isinstance(self.luminosity, list):
                self.luminosity = [str(format(np.log10(10.**float(lum)/L_sun), '.2f'))
                                   for lum in self.luminosity]
        elif self.luminosity_units == 'cgs':
            if isinstance(self.luminosity, str):
                self.luminosity = str(format(
                np.log10(float(self.luminosity)/L_sun), '.2f')
                )
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
