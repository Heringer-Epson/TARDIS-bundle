#!/usr/bin/env python

import numpy as np
import tardis
from density_interpolator import rho11fe_velocity_to_density

tardis_path = tardis.__path__[0]

def L2logL(_L):
    return str(format(np.log10(_L), '.3f'))

t_11fe = ['3.7', '5.9', '9.0', '12.1', '16.1', '19.1', '22.4', '28.3'] 
t2L_11fe = {
  '3.7': L2logL(0.08e9), '5.9': L2logL(0.32e9), '9.0': L2logL(1.1e9),
  '12.1': L2logL(2.3e9), '16.1': L2logL(3.2e9), '19.1': L2logL(3.5e9),
  '22.4': L2logL(3.2e9), '28.3': L2logL(2.3e9)}
t2vphot_11fe = {
  '3.7': '13300', '5.9': '12400', '9.0': '11300', '12.1': '10700',
  '16.1': '9000', '19.1': '7850', '22.4': '6700', '28.3': '4550'}

class Input_Parameters(object):
    """Set input parameters for running TARDIS.
    """

    def __init__(self, event, case=None, StoN=None, custom_par=None,
                 run_uncertainties=False):
        
        self.event = event
        self.custom_par = custom_par
        self.run_uncertainties = run_uncertainties
        self.smoothing_window = 21
        self.N_MC_runs = 3000
        
        #Variables that are common amonng all pre-defined events. 
        self.atomic_data = tardis_path + '/data/kurucz_cd23_chianti_H_He.h5'
        self.spec_start = '500'
        self.spec_stop = '20000'
        self.spec_num = '10000'       
        #self.seeds = '23111963'
        self.seeds = '23111970'
        self.delta = 'Default'

        #Set nlte species.
        self.ntle = ''

        self.el_setting = {'el': 'None', 'v_start': '0', 'v_stop': '100000',
                          'set_value': '0.0'}
        self.el_adding = {'el': 'None', 'v_start': '0', 'v_stop': '100000',
                          'add': '0.0'}
        self.el1_scaling = {'el': 'None', 'v_start': '0', 'v_stop': '100000',
                            'factors': '1.0'}
        self.el2_scaling = {'el': 'None', 'v_start': '0', 'v_stop': '100000',
                            'factors': '1.0'}

        self.elements = [
          'H', 'He', 'Li', 'B', 'Be', 'C', 'N', 'O', 'F', 'Ne', 'Na',
          'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc',
          'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Fe0', 'Ni0'] 

        #Publication quality.
        if StoN == 'super-high':
            self.num_packs = '2.0e+5'
            self.iterations = '20'
            self.last_num_packs = '2.0e+6'
            self.num_virtual_packs = '5'
          
        #Publication quality.
        elif StoN == 'high':
            self.num_packs = '2.0e+5'
            self.iterations = '20'
            self.last_num_packs = '5.0e+5'
            self.num_virtual_packs = '5'

        #For high SN in last iter.
        elif StoN == 'medium-high':
            self.num_packs = '1.0e+4'
            self.iterations = '15'
            #self.last_num_packs = '1.0e+5'
            self.last_num_packs = '5.0e+5'
            self.num_virtual_packs = '5'

        #For Kromer plots.
        elif StoN == 'medium':
            self.num_packs = '1.0e+4'
            self.iterations = '15'
            #self.last_num_packs = '1.0e+5'
            self.last_num_packs = '2.0e+5'
            self.num_virtual_packs = '5'

        #For fast runs.
        elif StoN == 'low':
            self.num_packs = '1.0e+4'
            #self.iterations = '15'
            self.iterations = '9'
            self.last_num_packs = '1.0e+4'
            self.num_virtual_packs = '5'

        elif StoN == 'very-low':
            self.num_packs = '4.0e+4'
            self.iterations = '1'
            self.last_num_packs = '5.0e+4'
            self.num_virtual_packs = '5'
        else:
            raise ValueError('Signal to noise (StoN) "%s" is not implemented'\
                             '.\n\n'%(StoN))
         
        #=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#    
        #=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#    
        #=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#    

        if event == 'fast':     
            #Input parameters for quick run (<10s). Adopts homogenous ejecta.
            #Used to check functionality of supporting codes.
            
            self.subdir = 'test-kromer/'
            #self.luminosity = '9.44'
            self.luminosity = '9.05'
            self.time_explosion = '13.0' 
            self.velocity_start = '11000'
            self.velocity_stop = '20000'
            self.luminosity_units = 'logsolar'
            self.structure_type = 'specific'
            self.abundance_type = 'uniform' 
            self.ionization = 'nebular'
            self.excitation = 'dilute-lte'
            self.rad_rates_type = 'dilute-blackbody'
            self.line_interaction = 'downbranch'
            self.extinction = 0.

            self.structure_num = '20'
            self.density_type = 'branch85_w7'

            self.abun = {}
            for el in (self.elements + ['Fe0', 'Ni0']):
                self.abun[el] = ['0.0']
          
            '''
            self.abun['Ni0'] = ['0.0']
            self.abun['Si'] = ['0.43']
            self.abun['Fe0'] = ['0.01']
            self.abun['Fe'] = ['0.01']
            self.abun['Co'] = ['0.0']
            self.abun['Ca'] = ['0.03']
            self.abun['S'] = ['0.19']
            self.abun['Mg'] = ['0.03']
            self.abun['Na'] = ['0.0']
            self.abun['C'] = ['0.09']
            self.abun['O'] = ['0.18']
            self.abun['Ti'] = ['0.00']
            self.abun['Ar'] = ['0.04']
            '''
            
            self.abun['Ni0'] = ['0.0']
            self.abun['Si'] = ['0.22']
            self.abun['Fe0'] = ['0.0']
            self.abun['Co'] = ['0.0']
            self.abun['Ca'] = ['0.03']
            self.abun['S'] = ['0.19']
            self.abun['Mg'] = ['0.03']
            self.abun['Na'] = ['0.0']
            self.abun['He'] = ['0.1']
            self.abun['C'] = ['0.2']
            self.abun['O'] = ['0.18']
            self.abun['Ti'] = ['0.01']
            self.abun['Ar'] = ['0.04']

            #=-=-=-=-= Pre-defined cases with simulation parameters.=-=-=-=-=#
            
            if case == 'single':
                self.subdir = 'fast_single_blank_up/'
                self.seeds = [self.seeds]
                self.delta = '1'
            elif case == 'multiple':
                self.subdir = 'fast_multiple/'
                self.seeds = list(np.arange(1,201,1).astype(int).astype(str))
            elif case == 'test_features_sine':
                self.subdir = 'fast_test-sine/'
                self.seeds = [self.seeds]                
            else:
                self.case_error(case, event)
                
            #Unsued variables in this setting (event):
            self.time_0 = 'blank'
            self.rho_0 = 'blank'
            self.v_0 = 'blank'
            self.density_value = 'blank'             
            self.velocity_array = 'blank'
            self.density_array = 'blank'
            self.distance = 'blank'
            #self.distance = '24.2 Mpc'
            self.exponent = 'blank'
            self.temperature_requested = 'blank'

        #=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#    
        #=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#    
        #=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#    

        elif event == '05bl':
            #Input parameters to reproduce the tomography of 2005bl according
            #to the data provided by Hachinger, which was used in
            #[[http://adsabs.harvard.edu/abs/2009MNRAS.399.1238H]].
                
            self.subdir = 'test-05bl/'
            self.extinction = -0.20

            #Velocity stop does not matter, as long as it is above the highest
            #velocity in the velocity_array.
            self.velocity_stop = '48000'
            
            self.luminosity_units = 'logsolar'
            self.structure_type = 'file' #specific #file
            self.abundance_type = 'file' #file,uniform,branch_w7
            self.time_0 = '100 s'
            
            self.ionization = 'nebular'
            self.excitation = 'dilute-lte'
            self.rad_rates_type = 'dilute-blackbody'

            #Below are default cases, which might be overwritten by some
            #specific cases.
            self.line_interaction = 'macroatom'
                
            #=-=-=-=-= Pre-defined cases with simulation parameters.=-=-=-=-=#
            
            #Used to compute the default tomography analysis.
            if case == 'default':    
                self.subdir = '05bl_default/'
                self.luminosity = ['8.520', '8.617', '8.745', '8.861', '8.594']
                self.time_explosion = ['11.0', '12.0', '14.0', '21.8', '29.9']
                self.velocity_start = ['8350', '8100', '7600', '6800', '3350']
                
            #Compute spectra with standard settings, where L is scaled by
            #1, 1/2, 1/3 and 1/4.
            elif case == 'default_L-scaled':    
                self.subdir = '05bl_default_L-scaled/'
                self.luminosity = []
                lum = [8.520, 8.617, 8.745, 8.861, 8.594]
                for scale in [1., 2., 3., 4.]:
                    self.luminosity += [str(format(np.log10(10.**l * scale), '.3f')) for l in list(lum)]                
                self.luminosity = self.luminosity * 2
                self.time_explosion = ['11.0', '12.0', '14.0', '21.8', '29.9'] * 8
                self.velocity_start = ['8350', '8100', '7600', '6800', '3350'] * 8
                self.line_interaction = ['downbranch'] * 20 + ['macroatom'] * 20

            #Used to compute L-grid.1
            elif case == 'L-grid':    
                self.subdir = '05bl_L-grid/'
                self.luminosity = [str(format(l, '.3f')) for l in np.log10(np.logspace(8.544, 9.72, 20))]
                self.luminosity = self.luminosity * 2 #One for each line_int mode.
                self.time_explosion = '21.8'
                self.velocity_start = '6800'
                self.line_interaction = ['downbranch'] * 20 + ['macroatom'] * 20

            #Compute spectra with standard settings, where L is scaled by
            #1, 1/2, 1/3 and 1/4.
            elif case == 'default_L-scaled_extra0.01Fe':    
                self.subdir = '05bl_default_L-scaled_extra0.01Fe/'
                self.el_adding = {'el': 'Fe0', 'v_start': 8000.,
                                  'v_stop': 1.e6, 'add': ['+0.01'] * 6}
                self.luminosity = []
                lum = [8.617, 8.745, 8.861]
                for scale in [0.8, 1.]:
                    self.luminosity += [str(format(np.log10(10.**l * scale), '.3f')) for l in list(lum)]                
                self.time_explosion = ['12.0', '14.0', '21.8'] * 2
                self.velocity_start = ['8100', '7600', '6800'] * 2
                self.line_interaction = ['macroatom'] * 6

            #Compute spectra with standard settings, where L is scaled by
            #1, 1/2, 1/3 and 1/4.
            elif case == 'default_L-scaled_extra0.01Fe-outer':    
                self.subdir = '05bl_default_L-scaled_extra0.01Fe-outer/'
                self.el_adding = {'el': 'Fe0', 'v_start': 15000.,
                                  'v_stop': 1.e6, 'add': ['+0.01'] * 3}
                self.luminosity = []
                lum = [8.617, 8.745, 8.861]
                for scale in [1.]:
                    self.luminosity += [str(format(np.log10(10.**l * scale), '.3f')) for l in list(lum)]                
                self.time_explosion = ['12.0', '14.0', '21.8'] * 1
                self.velocity_start = ['8100', '7600', '6800'] * 1
                self.line_interaction = ['macroatom'] * 3

            elif case == '12d_C-scaled_v0':
                self.subdir = '05bl_12d_C-scaled_v0/'    

                self.luminosity = '8.617'
                self.time_explosion = '12.0'
                self.velocity_start = '8100'

                scales = ['0.00', '0.20', '0.50', '1.00', '2.00', '5.00']
                                
                self.line_interaction = 'macroatom'
                self.el1_scaling = {'el': 'C', 'v_start': 0., 'v_stop': 1.e6,
                                    'factors': scales}             
            
            #Test case - Useful to check relevant line trnasitions for Si, C and Fe.
            elif case == 'test':       
                self.subdir = '05bl_test-case/'
                self.el_adding = {'el': 'Fe0', 'v_start': '8000',
                                  'v_stop': '100000', 'add': ['+0.01']}

                self.luminosity = [str(format(np.log10(10.**8.745 * 0.8), '.3f'))]                
                self.time_explosion = ['14.0']    
                self.velocity_start = ['7600']
                self.line_interaction = 'macroatom'

            else:
                self.case_error(case, event) 
                
            #Unsued variable in this setting (event):
            self.velocity_array = 'blank'                
            self.density_array = 'blank'
            self.time_0 = 'blank'
            self.abun = 'blank'
            self.rho_0 = 'blank'
            self.v_0 = 'blank'
            self.density_value = 'blank'             
            self.distance = 'blank' 
            self.density_type = 'blank'
            self.structure_num = 'blank'
            self.exponent = 'blank'
            self.abun_Ni = 'blank'
            self.abun_Si = 'blank'
            self.abun_Fe = 'blank'
            self.abun_Co = 'blank'
            self.abun_Ca = 'blank' 
            self.abun_S  = 'blank'
            self.abun_Mg = 'blank'
            self.abun_Na = 'blank'
            self.abun_C  = 'blank'
            self.abun_O  = 'blank'
            self.abun_Ti = 'blank'
            self.abun_Ar = 'blank'
            self.temperature_requested = 'blank'

        #=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#    
        #=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#    
        #=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#    
                             
        elif event == '11fe':     

            self.subdir = 'q-11fe/'
            self.extinction = -0.014

            #Velocity stop does not matter, as long as it is above the highest
            #velocity in the velocity_array.
            self.velocity_stop = '40000'
            
            self.luminosity_units = 'logsolar'
            self.structure_type = 'file' #specific #file
            self.abundance_type = 'file' #file,uniform,branch_w7
            self.time_0 = '100 s'
            
            self.ionization = 'nebular'
            self.excitation = 'dilute-lte'
            self.rad_rates_type = 'dilute-blackbody'

            #Below are default cases, which might be overwritten by some
            #specific cases.
            self.line_interaction = 'macroatom'

            #=-=-=-=-= Pre-defined cases with simulation parameters.=-=-=-=-=#

            #Used to compute the default tomography analysis.
            if case == 'default':    
                self.subdir = '11fe_default/'
                self.luminosity  = [str(format(np.log10(l), '.3f')) for l in [0.08e9, 0.32e9, 1.1e9, 2.3e9, 3.2e9, 3.5e9, 3.2e9, 2.3e9]]
                self.time_explosion = ['3.7', '5.9', '9.0', '12.1', '16.1', '19.1', '22.4', '28.3'] 
                self.velocity_start = ['13300', '12400', '11300', '10700', '9000', '7850', '6700', '4550']
                self.line_interaction = ['downbranch'] * 8 #Similar to original work.

                #self.subdir = '11fe_t37_kromer/'
                #self.luminosity = [self.luminosity[0]]
                #self.time_explosion = [self.time_explosion[0]]
                #self.velocity_start = [self.velocity_start[0]]
                #self.line_interaction = [self.line_interaction[0]]

            #Used to compute the default tomography analysis.
            elif case == 'best_delta':    
                self.subdir = '11fe_best_delta/'
                self.luminosity  = [str(format(np.log10(l), '.3f')) for l in [0.08e9, 0.32e9, 1.1e9, 2.3e9, 3.2e9, 3.5e9, 3.2e9, 2.3e9]]
                self.time_explosion = ['3.7', '5.9', '9.0', '12.1', '16.1', '19.1', '22.4', '28.3'] 
                self.velocity_start = ['13300', '12400', '11300', '10700', '9000', '7850', '6700', '4550']
                self.line_interaction = ['downbranch'] * 8 #Similar to original work.
                self.delta = '1'

                self.el_setting = {'el': 'C', 'v_start': '7850', 'v_stop': '16000',
                                    'set_value': '0.01'} 
                self.el1_scaling = {'el': 'C', 'v_start': '7850', 'v_stop': '13300',
                                    'factors': ['0.2'] * 8} 
                self.el2_scaling = {'el': 'C', 'v_start': '13300', 'v_stop': '16000',
                                    'factors': ['1.00'] * 8}   

            #Used to compute 11fe models that are slightly hotter to check
            #what happens to the carbon trough.
            elif case == 'best_hot':    
                self.subdir = '11fe_best_hot/'
                self.luminosity  = [str(format(np.log10(l * 1.3), '.3f')) for l in [0.32e9, 1.1e9, 2.3e9, 3.2e9, 3.5e9]]
                self.time_explosion = ['5.9', '9.0', '12.1', '16.1', '19.1'] 
                self.velocity_start = ['12400', '11300', '10700', '9000', '7850']
                self.line_interaction = ['downbranch'] * 5 #Similar to original work.

                self.el_setting = {'el': 'C', 'v_start': '7850', 'v_stop': '16000',
                                    'set_value': '0.01'} 
                self.el1_scaling = {'el': 'C', 'v_start': '7850', 'v_stop': '13300',
                                    'factors': ['0.2'] * 5} 
                self.el2_scaling = {'el': 'C', 'v_start': '13300', 'v_stop': '16000',
                                    'factors': ['1.00'] * 5}   
                                    
                                    
            #Compute spectra with standard settings, where L is scaled by
            #1, 1/2, 1/3 and 1/4.
            elif case == 'default_L-scaled':    
                self.subdir = '11fe_default_L-scaled_UP/'
                self.luminosity = []
                lum = [np.log10(l) for l in [0.08e9, 0.32e9, 1.1e9, 2.3e9, 3.2e9, 3.5e9, 3.2e9, 2.3e9]]
                for scale in [1., 0.5, 0.33, 0.25]:
                    self.luminosity += [str(format(np.log10(10.**l * scale), '.3f')) for l in list(lum)]        
                self.time_explosion = ['3.7', '5.9', '9.0', '12.1', '16.1', '19.1', '22.4', '28.3']*4
                self.velocity_start = ['13300', '12400', '11300', '10700', '9000', '7850', '6700', '4550']*4     
                self.line_interaction = ['downbranch'] * 32

            #Change photosphere position at pre-maximum using 05bl luminosity.
            elif case == 'vphoto_premax':
                self.subdir = '11fe_vphoto_premax/'                     
                self.luminosity = ['8.617'] * 5 + [str(format(np.log10(2.3e9), '.3f'))] * 4
                self.time_explosion = '12.1'  
                self.velocity_start = ['5000', '6000', '8100', '10700', '12200'] * 2
                self.line_interaction = 'downbranch' #Similar to original work.

            #Change photosphere position at maximum for L_11fe and L_05bl.
            elif case == 'vphoto_max':
                self.subdir = '11fe_vphoto_max/'       
                self.luminosity = ['8.861'] * 4 + [str(format(np.log10(3.5e9), '.3f'))] * 4 
                self.time_explosion = '19.1'   
                self.velocity_start = ['6800', '7850', '8600',  '9000'] * 2
                self.line_interaction = 'downbranch' #Similar to original work.

            #Change photosphere position at post-maximum for L_11fe and L_05bl.
            elif case == 'vphoto_postmax':
                self.subdir = '11fe_vphoto_postmax/'       
                self.luminosity = ['8.594'] * 5 + [str(format(np.log10(2.3e9), '.3f'))] * 5 
                self.time_explosion = '28.3'   
                self.velocity_start = ['3050', '3600', '4550', '5500', '6500'] * 2
                self.line_interaction = 'downbranch' #Similar to original work.

            #Makes Fig. 1 of the carbon paper.
            elif case == 'C-scan':
                t2vstop = {'3.7': 13000., '5.9': 12000., '9.0': 11000.,
                           '12.1': 10500, '16.1': 8500., '19.1': 7500.}
           
                _t = self.custom_par
                self.subdir = '11fe_' + str(int(float(_t))) + 'd_C-scan/'     
                self.luminosity = t2L_11fe[_t]
                self.time_explosion = _t    
                self.velocity_start = t2vphot_11fe[_t]
                v_stop = np.arange(t2vstop[_t], 15000., 500.)
                v_stop = list(v_stop.astype('int').astype('str'))
                self.line_interaction = 'downbranch' 
                self.el1_scaling = {'el': 'C', 'v_start': '0',
                                    'v_stop': v_stop, 'factors': '0'}
            
            elif case == '6d_C-scan':
                self.subdir = '11fe_6d_C-scan/'     

                v_stop_list = list(np.arange(12000., 15000., 500.).astype('int').astype('str'))
                self.luminosity = [self.lum2loglum(0.32e9)]                
                self.luminosity = self.luminosity * len(v_stop_list)

                self.time_explosion = '5.9'    
                self.velocity_start = '12400'
                self.line_interaction = 'downbranch' 
                self.el1_scaling = {'el': 'C', 'v_start': '0',
                                    'v_stop': v_stop_list, 'factors': '0'}

            elif case == '9d_C-scan':
                self.subdir = '11fe_9d_C-scan/'     

                v_stop_list = list(np.arange(11000., 15000., 500.).astype('int').astype('str'))
                self.luminosity = [self.lum2loglum(1.1e9)]                
                self.luminosity = self.luminosity * len(v_stop_list)

                self.time_explosion = '9.0'    
                self.velocity_start = '11300'
                self.line_interaction = 'downbranch' 
                self.el1_scaling = {'el': 'C', 'v_start': '0',
                                    'v_stop': v_stop_list, 'factors': '0'}

           
            elif case == '12d_C-scan':
                self.subdir = '11fe_12d_C-scan/'     

                v_stop_list = list(np.arange(10500., 15000., 500.).astype('int').astype('str'))
                self.luminosity = [self.lum2loglum(2.3e9 )]                
                self.luminosity = self.luminosity * len(v_stop_list)

                self.time_explosion = '12.1'    
                self.velocity_start = '10700'
                self.line_interaction = 'downbranch' 
                self.el1_scaling = {'el': 'C', 'v_start': '0',
                                    'v_stop': v_stop_list, 'factors': '0'}

            elif case == '16d_C-scan':
                self.subdir = '11fe_16d_C-scan/'     

                v_stop_list = list(np.arange(8500., 15000., 500.).astype('int').astype('str'))
                self.luminosity = [self.lum2loglum(3.2e9)]                
                self.luminosity = self.luminosity * len(v_stop_list)

                self.time_explosion = '16.1'    
                self.velocity_start = '9000'
                self.line_interaction = 'downbranch' 
                self.el1_scaling = {'el': 'C', 'v_start': '0',
                                    'v_stop': v_stop_list, 'factors': '0'}

            elif case == '19d_C-scan':
                self.subdir = '11fe_19d_C-scan/'     

                v_stop_list = list(np.arange(7500., 15000., 500.).astype('int').astype('str'))
                self.luminosity = [self.lum2loglum(3.5e9) ]                
                self.luminosity = self.luminosity * len(v_stop_list)

                self.time_explosion = '19.1'    
                self.velocity_start = '7850'
                self.line_interaction = 'downbranch' 
                self.el1_scaling = {'el': 'C', 'v_start': '0',
                                    'v_stop': v_stop_list, 'factors': '0'}

            elif case == '6d_C-plateaus_scaling':
                self.subdir = '11fe_6d_C-plateaus_scaling_SN/'     
                self.smoothing_window = 7

                self.luminosity = self.lum2loglum(0.32e9)                
                
                C1_scaling, C2_scaling = [], []
                scales = ['0.00', '0.05', '0.1', '0.2', '0.5',
                          '1.00', '2.00', '5.00', '10.00']
                for s1 in scales:
                    for s2 in scales:
                        if float(s2) >= float(s1):
                           C1_scaling.append(s1)                             
                           C2_scaling.append(s2)                             
                
                self.time_explosion = '5.9'    
                self.velocity_start = '12400'
                self.line_interaction = 'downbranch' 
                self.el_setting = {'el': 'C', 'v_start': '7850', 'v_stop': '16000',
                                    'set_value': '0.01'} 
                self.el1_scaling = {'el': 'C', 'v_start': '7850', 'v_stop': '13300',
                                    'factors': C1_scaling} 
                self.el2_scaling = {'el': 'C', 'v_start': '13300', 'v_stop': '16000',
                                    'factors': C2_scaling}  

            elif case == '9d_C-plateaus_scaling':
                self.subdir = '11fe_9d_C-plateaus_scaling_SN/'     
                self.smoothing_window = 7

                self.luminosity = self.lum2loglum(1.1e9)                
                
                C1_scaling, C2_scaling = [], []
                scales = ['0.00', '0.05', '0.1', '0.2', '0.5',
                          '1.00', '2.00', '5.00', '10.00']
                for s1 in scales:
                    for s2 in scales:
                        if float(s2) >= float(s1):
                           C1_scaling.append(s1)                             
                           C2_scaling.append(s2)                             
                
                self.time_explosion = '9.0'    
                self.velocity_start = '11300'
                self.line_interaction = 'downbranch' 
                self.el_setting = {'el': 'C', 'v_start': '7850', 'v_stop': '16000',
                                    'set_value': '0.01'} 
                self.el1_scaling = {'el': 'C', 'v_start': '7850', 'v_stop': '13300',
                                    'factors': C1_scaling} 
                self.el2_scaling = {'el': 'C', 'v_start': '13300', 'v_stop': '16000',
                                    'factors': C2_scaling}                         
            
            elif case == '12d_C-plateaus_scaling':
                self.subdir = '11fe_12d_C-plateaus_scaling_SN/'     
                self.smoothing_window = 7

                self.luminosity = self.lum2loglum(2.3e9)                
                
                C1_scaling, C2_scaling = [], []
                scales = ['0.00', '0.05', '0.1', '0.2', '0.5',
                          '1.00', '2.00', '5.00', '10.00']
                for s1 in scales:
                    for s2 in scales:
                        if float(s2) >= float(s1):
                           C1_scaling.append(s1)                             
                           C2_scaling.append(s2)                             
                
                self.time_explosion = '12.1'    
                self.velocity_start = '10700'
                self.line_interaction = 'downbranch' 
                self.el_setting = {'el': 'C', 'v_start': '7850', 'v_stop': '16000',
                                    'set_value': '0.01'} 
                self.el1_scaling = {'el': 'C', 'v_start': '7850', 'v_stop': '13300',
                                    'factors': C1_scaling} 
                self.el2_scaling = {'el': 'C', 'v_start': '13300', 'v_stop': '16000',
                                    'factors': C2_scaling}   

            elif case == '16d_C-plateaus_scaling':
                self.subdir = '11fe_16d_C-plateaus_scaling_SN/'     
                self.smoothing_window = 7

                self.luminosity = self.lum2loglum(3.2e9)                
                
                C1_scaling, C2_scaling = [], []
                scales = ['0.00', '0.05', '0.1', '0.2', '0.5',
                          '1.00', '2.00', '5.00', '10.00']
                for s1 in scales:
                    for s2 in scales:
                        if float(s2) >= float(s1):
                           C1_scaling.append(s1)                             
                           C2_scaling.append(s2)                             
                
                self.time_explosion = '16.1'    
                self.velocity_start = '9000'
                self.line_interaction = 'downbranch' 
                self.el_setting = {'el': 'C', 'v_start': '7850', 'v_stop': '16000',
                                    'set_value': '0.01'} 
                self.el1_scaling = {'el': 'C', 'v_start': '7850', 'v_stop': '13300',
                                    'factors': C1_scaling} 
                self.el2_scaling = {'el': 'C', 'v_start': '13300', 'v_stop': '16000',
                                    'factors': C2_scaling}   
                                    
                                                                        
            elif case == '19d_C-plateaus_scaling':
                self.subdir = '11fe_19d_C-plateaus_scaling_SN/'     
                self.smoothing_window = 7

                self.luminosity = self.lum2loglum(3.5e9)                
                
                C1_scaling, C2_scaling = [], []
                scales = ['0.00', '0.05', '0.1', '0.2', '0.5',
                          '1.00', '2.00', '5.00', '10.00']
                for s1 in scales:
                    for s2 in scales:
                        if float(s2) >= float(s1):
                           C1_scaling.append(s1)                             
                           C2_scaling.append(s2)                             
                
                self.time_explosion = '19.1'    
                self.velocity_start = '7850'
                self.line_interaction = 'downbranch' 
                self.el_setting = {'el': 'C', 'v_start': '7850', 'v_stop': '16000',
                                    'set_value': '0.01'} 
                self.el1_scaling = {'el': 'C', 'v_start': '7850', 'v_stop': '13300',
                                    'factors': C1_scaling} 
                self.el2_scaling = {'el': 'C', 'v_start': '13300', 'v_stop': '16000',
                                    'factors': C2_scaling}                         

            elif case == '6d_C-best':
                self.subdir = '11fe_6d_C-best/'     
                self.smoothing_window = 7

                self.luminosity = self.lum2loglum(0.32e9)                
                
                C1_scaling = ['0.00', '0.2'] * 2
                C2_scaling = ['2.00', '1.00'] * 2
                self.line_interaction = ['downbranch'] * 2 + ['macroatom'] * 2
                self.excitation = ['dilute-lte'] * 4
                
                self.time_explosion = '5.9'    
                self.velocity_start = '12400'
                self.el_setting = {'el': 'C', 'v_start': '7850', 'v_stop': '16000',
                                    'set_value': '0.01'} 
                self.el1_scaling = {'el': 'C', 'v_start': '7850', 'v_stop': '13300',
                                    'factors': C1_scaling} 
                self.el2_scaling = {'el': 'C', 'v_start': '13300', 'v_stop': '16000',
                                    'factors': C2_scaling}  

            elif case == '9d_C-best':
                self.subdir = '11fe_9d_C-best/'     
                self.smoothing_window = 7

                self.luminosity = self.lum2loglum(1.1e9)                
                
                C1_scaling = ['0.00', '0.2'] * 2
                C2_scaling = ['2.00', '1.00'] * 2
                self.line_interaction = ['downbranch'] * 2 + ['macroatom'] * 2
                self.excitation = ['dilute-lte'] * 4                          
                
                self.time_explosion = '9.0'    
                self.velocity_start = '11300'
                self.el_setting = {'el': 'C', 'v_start': '7850', 'v_stop': '16000',
                                    'set_value': '0.01'} 
                self.el1_scaling = {'el': 'C', 'v_start': '7850', 'v_stop': '13300',
                                    'factors': C1_scaling} 
                self.el2_scaling = {'el': 'C', 'v_start': '13300', 'v_stop': '16000',
                                    'factors': C2_scaling}                         
            
            elif case == '12d_C-best':
                self.subdir = '11fe_12d_C-best/'     
                self.smoothing_window = 7

                self.luminosity = self.lum2loglum(2.3e9)                
                
                C1_scaling = ['0.00', '0.2'] * 2
                C2_scaling = ['2.00', '1.00'] * 2
                self.line_interaction = ['downbranch'] * 2 + ['macroatom'] * 2
                self.excitation = ['dilute-lte'] * 4                    
                self.seeds = '23111963'
                
                self.time_explosion = '12.1'    
                self.velocity_start = '10700'
                self.el_setting = {'el': 'C', 'v_start': '7850', 'v_stop': '16000',
                                    'set_value': '0.01'} 
                self.el1_scaling = {'el': 'C', 'v_start': '7850', 'v_stop': '13300',
                                    'factors': C1_scaling} 
                self.el2_scaling = {'el': 'C', 'v_start': '13300', 'v_stop': '16000',
                                    'factors': C2_scaling}   

            elif case == '16d_C-best':
                self.subdir = '11fe_16d_C-best/'     
                self.smoothing_window = 7

                self.luminosity = self.lum2loglum(3.2e9)                
                
                C1_scaling = ['0.00', '0.2'] * 2
                C2_scaling = ['2.00', '1.00'] * 2
                self.line_interaction = ['downbranch'] * 2 + ['macroatom'] * 2
                self.excitation = ['dilute-lte'] * 4                    
                
                self.time_explosion = '16.1'    
                self.velocity_start = '9000'
                self.el_setting = {'el': 'C', 'v_start': '7850', 'v_stop': '16000',
                                    'set_value': '0.01'} 
                self.el1_scaling = {'el': 'C', 'v_start': '7850', 'v_stop': '13300',
                                    'factors': C1_scaling} 
                self.el2_scaling = {'el': 'C', 'v_start': '13300', 'v_stop': '16000',
                                    'factors': C2_scaling}   
                                    
                                                                        
            elif case == '19d_C-best':
                self.subdir = '11fe_19d_C-best/'     
                self.smoothing_window = 7

                self.luminosity = self.lum2loglum(3.5e9)                
                
                C1_scaling = ['0.00', '0.2'] * 2
                C2_scaling = ['2.00', '1.00'] * 2
                self.line_interaction = ['downbranch'] * 2 + ['macroatom'] * 2
                self.excitation = ['dilute-lte'] * 4                          
                
                self.time_explosion = '19.1'    
                self.velocity_start = '7850'
                self.el_setting = {'el': 'C', 'v_start': '7850', 'v_stop': '16000',
                                    'set_value': '0.01'} 
                self.el1_scaling = {'el': 'C', 'v_start': '7850', 'v_stop': '13300',
                                    'factors': C1_scaling} 
                self.el2_scaling = {'el': 'C', 'v_start': '13300', 'v_stop': '16000',
                                    'factors': C2_scaling}  

            elif case == '12d_C-neutral':
                self.subdir = '11fe_12d_C-neutral_test/'     

                self.luminosity = [self.lum2loglum(2.3e9)]           
                C_scaling = ['10.00']
                
                self.luminosity = ['8.5'] * 2 + [self.lum2loglum(2.3e9)] * 2           
                C_scaling = ['0.00', '10.00'] * 2
                                            
                self.time_explosion = '12.1'    
                self.velocity_start = '10700'
                self.line_interaction = 'downbranch' 
                self.el_setting = {'el': 'C', 'v_start': '7850', 'v_stop': '16000',
                                    'set_value': '0.01'} 
                self.el1_scaling = {'el': 'C', 'v_start': '7850', 'v_stop': '16000',
                                    'factors': C_scaling}  


            elif case == 'W7_C-prof':    
                self.subdir = '11fe_W7_C-prof/'
                self.luminosity  = [str(format(np.log10(l), '.3f')) for l in [0.08e9, 0.32e9, 1.1e9, 2.3e9, 3.2e9, 3.5e9, 3.2e9, 2.3e9]]
                self.time_explosion = ['3.7', '5.9', '9.0', '12.1', '16.1', '19.1', '22.4', '28.3'] 
                self.velocity_start = ['13300', '12400', '11300', '10700', '9000', '7850', '6700', '4550']
                self.line_interaction = 'downbranch'

                self.el_setting = {'el': 'C', 'v_start': '7850', 'v_stop': '20000',
                                    'set_value': '0.01'} 
                self.el1_scaling = {'el': 'C', 'v_start': '7850', 'v_stop': '14200',
                                    'factors': '0.0001'} 
                self.el2_scaling = {'el': 'C', 'v_start': '14200', 'v_stop': '20000',
                                    'factors': '50.00'}  

            elif case == 'ddet_C-prof':    
                self.subdir = '11fe_ddet_C-prof/'
                self.luminosity  = [str(format(np.log10(l), '.3f')) for l in [0.08e9, 0.32e9, 1.1e9, 2.3e9, 3.2e9, 3.5e9, 3.2e9, 2.3e9]]
                self.time_explosion = ['3.7', '5.9', '9.0', '12.1', '16.1', '19.1', '22.4', '28.3'] 
                self.velocity_start = ['13300', '12400', '11300', '10700', '9000', '7850', '6700', '4550']
                self.line_interaction = 'downbranch'

                self.el_setting = {'el': 'C', 'v_start': '7850', 'v_stop': '18000',
                                    'set_value': '0.00001'} 
                self.el1_scaling = {'el': 'C', 'v_start': '13500', 'v_stop': '14500',
                                    'factors': '15000.00'} 
                self.el2_scaling = {'el': 'C', 'v_start': '14500', 'v_stop': '18000',
                                    'factors': '60.00'}  
            
            
            else:
                self.case_error(case, event) 
            
            #Unsued variable in this setting (event):
            self.velocity_array = 'blank'                
            self.density_array = 'blank'
            self.abun = 'blank'                       
            self.rho_0 = 'blank'
            self.v_0 = 'blank'
            self.density_value = 'blank'             
            self.distance = 'blank' 
            self.density_type = 'blank'
            self.structure_num = 'blank'
            self.exponent = 'blank'
            self.abun_Ni = 'blank'
            self.abun_Si = 'blank'
            self.abun_Fe = 'blank'
            self.abun_Co = 'blank'
            self.abun_Ca = 'blank' 
            self.abun_S  = 'blank'
            self.abun_Mg = 'blank'
            self.abun_Na = 'blank'
            self.abun_C  = 'blank'
            self.abun_O  = 'blank'
            self.abun_Ti = 'blank'
            self.abun_Ar = 'blank'
            self.temperature_requested = 'blank'

        #=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#    
        #=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#    
        #=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#    


        #If event is not defined in this input file, raise ValueError.
        else:
            raise ValueError('Event "%s" is currently not implemented.\n\n'
                             %event)    

        #=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#    
        #=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#    
        #=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#    

    def case_error(self, case, event):
        raise ValueError('Case "%s" is not allowed in event "%s".\n\n'
                         %(case, event))

    def lum2loglum(self, lum):
        return str(format(np.log10(lum), '.3f'))
