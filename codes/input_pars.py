#!/usr/bin/env python

import numpy as np
from density_interpolator import rho11fe_velocity_to_density

class Input_Parameters(object):
    """Set input parameters for running TARDIS.
    """

    def __init__(self, event, case=None, StoN=None, run_uncertainties=False):
        
        self.event = event
        self.run_uncertainties = run_uncertainties
        self.smoothing_window = 21
        self.N_MC_runs = 3000
        
        #Variables that are common amonng all pre-defined events. 
        self.atomic_data = 'kurucz_cd23_chianti_H_He.h5'
        self.spec_start = '500'
        self.spec_stop = '20000'
        self.spec_num = '10000'       
        self.seeds = '23111963'

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
        if StoN == 'high':
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
            self.iterations = '15'
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
            self.luminosity = '9.44'
            #self.luminosity = '9.05'
            self.time_explosion = '13.0' 
            self.velocity_start = '11000'
            self.velocity_stop = '20000'
            self.luminosity_units = 'logsolar'
            self.structure_type = 'specific'
            self.abundance_type = 'uniform' 
            self.ionization = 'nebular'
            self.excitation = 'dilute-lte'
            self.rad_rates_type = 'dilute-blackbody'
            self.line_interaction = 'macroatom'
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
            self.abun['Si'] = ['0.52']
            self.abun['Fe0'] = ['0.0']
            self.abun['Co'] = ['0.0']
            self.abun['Ca'] = ['0.03']
            self.abun['S'] = ['0.19']
            self.abun['Mg'] = ['0.03']
            self.abun['Na'] = ['0.0']
            self.abun['C'] = ['0.0']
            self.abun['O'] = ['0.18']
            self.abun['Ti'] = ['0.01']
            self.abun['Ar'] = ['0.04']

            #=-=-=-=-= Pre-defined cases with simulation parameters.=-=-=-=-=#
            
            if case == 'single':
                self.subdir = 'fast_single_blank/'
                self.seeds = [self.seeds]
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
                self.luminosity  = ['0.08e9', '0.32e9', '1.1e9', '2.3e9', '3.2e9', '3.5e9', '3.2e9', '2.3e9']       
                self.time_explosion = ['3.7', '5.9', '9.0', '12.1', '16.1', '19.1', '22.4', '28.3'] 
                self.velocity_start = ['13300', '12400', '11300', '10700', '9000', '7850', '6700', '4550']
                self.line_interaction = 'downbranch' #Similar to original work.

            #Compute spectra with standard settings, where L is scaled by
            #1, 1/2, 1/3 and 1/4.
            elif case == 'default_L-scaled':    
                self.subdir = '11fe_default_L-scaled_UP/'
                self.luminosity = []
                lum = [np.log10(l) for l in [0.08e9, 0.32e9, 1.1e9, 2.3e9, 3.2e9, 3.5e9, 3.2e9, 2.3e9]]
                for scale in [1., 0.5, 0.33, 0.25]:
                    self.luminosity += [str(format(np.log10(10.**l * scale), '.3f')) for l in list(lum)]        
                self.luminosity = self.luminosity * 2 #One for each line_int mode.                
                self.time_explosion = ['3.7', '5.9', '9.0', '12.1', '16.1', '19.1', '22.4', '28.3']*8
                self.velocity_start = ['13300', '12400', '11300', '10700', '9000', '7850', '6700', '4550']*8      
                self.line_interaction = ['downbranch'] * 32 + ['macroatom'] * 32
                
            #Used to compute L-grid.
            elif case == 'L-grid':    
                self.subdir = '11fe_L-grid_UP/'   
                self.luminosity = [str(format(l, '.3f')) for l in np.log10(np.logspace(8.544, 9.72, 20))]
                self.luminosity = self.luminosity * 2 #One for each line_int mode.
                self.time_explosion = '19.1'    
                self.velocity_start = '7850'
                self.line_interaction = ['downbranch'] * 20 + ['macroatom'] * 20

            #Used to compute Ti-grid.
            elif case == 'Ti-grid':    
                self.subdir = '11fe_Ti-grid/'   
                self.luminosity = ([str(format(np.log10(3.5e9), '.3f'))] * 10        
                                   + [str(format(np.log10(3.5e9 / 4.), '.3f'))] * 10)        
                self.luminosity = self.luminosity * 2 #One for each line_int mode.
                self.time_explosion = '19.1'    
                self.velocity_start = '7850'
                self.TiCr_scaling = ['0.00', '0.05', '0.1', '0.2', '0.5', '1.00',
                                     '2.00', '5.00', '10.00', '20.00'] * 4
                self.line_interaction = ['downbranch'] * 20 + ['macroatom'] * 20

            #Used to compute Fe-grid.
            elif case == 'Fe-grid':    
                self.subdir = '11fe_Fe-grid/'   
                self.luminosity = ([str(format(np.log10(3.5e9), '.3f'))] * 10        
                                   + [str(format(np.log10(3.5e9 / 4.), '.3f'))] * 10)        
                self.luminosity = self.luminosity * 2 #One for each line_int mode.
                self.time_explosion = '19.1'    
                self.velocity_start = '7850'
                self.Fe_scaling = ['0.00', '0.05', '0.1', '0.2', '0.5', '1.00',
                                     '2.00', '5.00', '10.00', '20.00'] * 4
                self.line_interaction = ['downbranch'] * 20 + ['macroatom'] * 20                

            #Used to compute a 2D grid in luminosity and Fe.
            elif case == 'L_Fe-grid_19d':    
                self.subdir = '11fe_L_Fe-grid/'   
                self.luminosity = [str(format(l, '.3f')) for l in
                                   np.log10(np.logspace(
                                   np.log10(3.5e9 / 4), np.log10(3.5e9 * 1.5), 15))]      
                self.luminosity = self.luminosity * 10 #One for each Fe value.
                Fe_list = ['0.00', '0.05', '0.1', '0.2', '0.5', '1.00',
                           '2.00', '5.00', '10.00', '20.00']
                self.Fe_scaling = []
                for Fe in Fe_list:
                    self.Fe_scaling += [Fe] * 15 #One for each L value.
                self.time_explosion = '19.1'    
                self.velocity_start = '7850'
                self.line_interaction = ['downbranch'] * 150

            #Used to compute a 2D grid in luminosity and Fe.
            elif case == 'L_Fe-grid_12d':    
                self.subdir = '11fe_L_Fe-grid_12d/'   
                self.luminosity = [str(format(l, '.3f')) for l in
                                   np.log10(np.logspace(
                                   np.log10(2.3e9 / 4), np.log10(2.3e9 * 1.5), 15))]      
                self.luminosity = self.luminosity * 10 #One for each Fe value.
                Fe_list = ['0.00', '0.05', '0.1', '0.2', '0.5', '1.00',
                           '2.00', '5.00', '10.00', '20.00']
                self.Fe_scaling = []
                for Fe in Fe_list:
                    self.Fe_scaling += [Fe] * 15 #One for each L value.
                self.time_explosion = '12.1'    
                self.velocity_start = '10700'
                self.line_interaction = ['downbranch'] * 150
    
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


            #Used to test the presence of the carbon trough -- run with kromer plots
            #texp= 12, 16 and 19 days. L scaled by 1, 0.5 and 0.25.
            elif case == 'kromer_L-scaled':
                self.subdir = '11fe_kromer_L-scaled/'                     
                self.luminosity = []
                lum = [np.log10(l) for l in [2.3e9, 3.2e9, 3.5e9]]
                for scale in [1., 0.5, 0.25]:
                    self.luminosity += [str(format(np.log10(10.**l * scale), '.3f')) for l in list(lum)]        
                self.time_explosion = ['12.1', '16.1', '19.1'] * 3
                self.velocity_start = ['10700', '9000', '7850'] * 3
                #self.line_interaction = 'downbranch' #Similar to original work.
                self.line_interaction = ['macroatom'] * 9 #Similar to original work.

            #Change time explosion at pre-maximum for 05bl Lum.
            elif case == 'texp_premax':
                self.subdir = '11fe_texp_premax/' 
                self.luminosity =    '8.617'
                self.time_explosion = ['6', '7', '8', '9', '12.1', '13', '14']   
                self.velocity_start = '10700'
                self.line_interaction = 'downbranch' #Similar to original work.

            #Note that at maximum the spectrum already agrees.
            #Change time explosion at post-maximum for 05bl Lum.
            elif case == 'texp_postmax':
                self.luminosity =    '8.594'
                self.time_explosion = ['23', '25', '28.3', '30', '32', '34']   
                self.velocity_start = '4550'
                self.line_interaction = 'downbranch' #Similar to original work.

            elif case == 'early_2D-grid':
                self.subdir = 'early-carbon-grid_kromer_new2/'     

                #L_scal = np.arange(0.6, 2.01, 0.1)
                L_scal = np.arange(0.2, 1.61, 0.1)
                self.luminosity = [self.lum2loglum(0.32e9 * l) for l in L_scal]                
                self.luminosity = self.luminosity * 10
                
                #scales = ['0.00', '0.05', '0.10', '0.20', '0.50',
                #          '1.00', '2.00', '5.00', '10.00', '20.00']

                scales = ['0.00', '0.20', '0.50', '1.00', '2.00', '5.00',
                          '10.00', '20.00', '50.00', '100.00']
                                
                #By default, Fe will be scaled at the expense of the most 
                #abundant element, which is Si in this case.
                Fe_scaling = []
                for s in scales:                    
                    Fe_scaling += [s] * 15
                                            
                self.time_explosion = '5.9'    
                self.velocity_start = '12400'
                self.line_interaction = ['macroatom'] * 150  
                self.el1_scaling = {'el': 'Fe0', 'v_start': 17000.,
                                    'v_stop': 1.e6, 'factors': Fe_scaling}                

            elif case == 'early_2D-grid_no-C':
                self.subdir = 'early-carbon-grid_no-C_new/'     
                
                #Exclude the layers that are nearly 100% carbon, as one
                #cannot scale another element to be 100% without significantly
                #changing the ejecta.
                #self.velocity_stop = '16800'
                L_scal = np.arange(0.6, 2.01, 0.1)
                self.luminosity = [self.lum2loglum(0.32e9 * l) for l in L_scal]                
                self.luminosity = self.luminosity * 10
                
                scales = ['0.00', '0.05', '0.10', '0.20', '0.50',
                          '1.00', '2.00', '5.00', '10.00', '20.00']
                                
                #By default, Fe will be scaled at the expense of the most 
                #abundant element, which is Si in this case.
                Fe_scaling = []
                for s in scales:                    
                    Fe_scaling += [s] * 15
                                            
                self.time_explosion = '5.9'    
                self.velocity_start = '12400'
                self.line_interaction = ['macroatom'] * 150  
                self.el1_scaling = {'el': 'Fe0', 'v_start': 13100.,
                                    'v_stop': 1.e6, 'factors': Fe_scaling}               
                self.el2_scaling = {'el': 'C', 'v_start': 13100.,
                                    'v_stop': 1.e6, 'factors': '0.0'}   

            elif case == '6d_2D-grid_v19590':
                self.subdir = '11fe_2D-grid_6d_v19590_UP/'     

                L_scal = np.arange(0.2, 1.61, 0.1)
                self.luminosity = [self.lum2loglum(0.32e9 * l) for l in L_scal]                
                self.luminosity = self.luminosity * 10

                scales = ['0.00', '0.20', '0.50', '1.00', '2.00', '5.00',
                          '10.00', '20.00', '50.00', '100.00']
                                
                #By default, Fe will be scaled at the expense of the most 
                #abundant element, which is Si in this case.
                Fe_scaling = []
                for s in scales:                    
                    Fe_scaling += [s] * 15
                                            
                self.time_explosion = '5.9'    
                self.velocity_start = '12400'
                self.line_interaction = ['macroatom'] * 150  
                self.el1_scaling = {'el': 'Fe0', 'v_start': 19590.,
                                    'v_stop': 1.e6, 'factors': Fe_scaling} 
            
            
            elif case == '12d_2D-grid_v19590':
                self.subdir = '11fe_2D-grid_12d_v19590_UP/'     

                L_scal = np.arange(0.2, 1.61, 0.1)
                self.luminosity = [self.lum2loglum(2.3e9 * l) for l in L_scal]                
                self.luminosity = self.luminosity * 10

                scales = ['0.00', '0.20', '0.50', '1.00', '2.00', '5.00',
                          '10.00', '20.00', '50.00', '100.00']
                                
                #By default, Fe will be scaled at the expense of the most 
                #abundant element, which is Si in this case.
                Fe_scaling = []
                for s in scales:                    
                    Fe_scaling += [s] * 15
                                            
                self.time_explosion = '12.1'    
                self.velocity_start = '10700'
                self.line_interaction = ['macroatom'] * 150  
                self.el1_scaling = {'el': 'Fe0', 'v_start': 19590.,
                                    'v_stop': 1.e6, 'factors': Fe_scaling} 

            elif case == '19d_2D-grid_v19590':
                self.subdir = '11fe_2D-grid_19d_v19590_UP/'     

                L_scal = np.arange(0.2, 1.61, 0.1)
                self.luminosity = [self.lum2loglum(3.5e9 * l) for l in L_scal]                
                self.luminosity = self.luminosity * 10

                scales = ['0.00', '0.20', '0.50', '1.00', '2.00', '5.00',
                          '10.00', '20.00', '50.00', '100.00']
                                
                #By default, Fe will be scaled at the expense of the most 
                #abundant element, which is Si in this case.
                Fe_scaling = []
                for s in scales:                    
                    Fe_scaling += [s] * 15
                                            
                self.time_explosion = '19.1'    
                self.velocity_start = '7850'
                self.line_interaction = ['macroatom'] * 150  
                self.el1_scaling = {'el': 'Fe0', 'v_start': 19780.,
                                    'v_stop': 1.e6, 'factors': Fe_scaling} 
            
            


            elif case == '6d_2D-grid_v13400':
                self.subdir = '11fe_2D-grid_6d_v13400_UP/'     

                L_scal = np.arange(0.2, 1.61, 0.1)
                self.luminosity = [self.lum2loglum(0.32e9 * l) for l in L_scal]                
                self.luminosity = self.luminosity * 10

                scales = ['0.00', '0.20', '0.50', '1.00', '2.00', '5.00',
                          '10.00', '20.00', '50.00', '100.00']
                                
                #By default, Fe will be scaled at the expense of the most 
                #abundant element, which is Si in this case.
                Fe_scaling = []
                for s in scales:                    
                    Fe_scaling += [s] * 15
                                            
                self.time_explosion = '5.9'    
                self.velocity_start = '12400'
                self.line_interaction = ['macroatom'] * 150  
                self.el1_scaling = {'el': 'Fe0', 'v_start': 13400.,
                                    'v_stop': 1.e6, 'factors': Fe_scaling} 
            
            
            elif case == '12d_2D-grid_v13400':
                self.subdir = '11fe_2D-grid_12d_v13400_UP/'     

                L_scal = np.arange(0.2, 1.61, 0.1)
                self.luminosity = [self.lum2loglum(2.3e9 * l) for l in L_scal]                
                self.luminosity = self.luminosity * 10

                scales = ['0.00', '0.20', '0.50', '1.00', '2.00', '5.00',
                          '10.00', '20.00', '50.00', '100.00']
                                
                #By default, Fe will be scaled at the expense of the most 
                #abundant element, which is Si in this case.
                Fe_scaling = []
                for s in scales:                    
                    Fe_scaling += [s] * 15
                                            
                self.time_explosion = '12.1'    
                self.velocity_start = '10700'
                self.line_interaction = ['macroatom'] * 150  
                self.el1_scaling = {'el': 'Fe0', 'v_start': 13400.,
                                    'v_stop': 1.e6, 'factors': Fe_scaling} 

            elif case == '19d_2D-grid_v13400':
                self.subdir = '11fe_2D-grid_19d_v13400_UP/'     

                L_scal = np.arange(0.2, 1.61, 0.1)
                self.luminosity = [self.lum2loglum(3.5e9 * l) for l in L_scal]                
                self.luminosity = self.luminosity * 10

                scales = ['0.00', '0.20', '0.50', '1.00', '2.00', '5.00',
                          '10.00', '20.00', '50.00', '100.00']
                                
                #By default, Fe will be scaled at the expense of the most 
                #abundant element, which is Si in this case.
                Fe_scaling = []
                for s in scales:                    
                    Fe_scaling += [s] * 15
                                            
                self.time_explosion = '19.1'    
                self.velocity_start = '7850'
                self.line_interaction = ['macroatom'] * 150  
                self.el1_scaling = {'el': 'Fe0', 'v_start': 13400.,
                                    'v_stop': 1.e6, 'factors': Fe_scaling}             
            
            
            elif case == '12d_C-scaled_v0':
                self.subdir = '11fe_12d_C-scaled/'     

                L_scal = np.arange(0.2, 1.61, 0.1)
                self.luminosity = [self.lum2loglum(2.3e9 * l) for l in L_scal]                
                self.luminosity = self.luminosity * 12

                scales = ['0.00', '0.01', '0.02', '0.05', '0.10', '0.20',
                          '0.50', '1.00', '2.00', '5.00', '10.00', '20.00']

                C_scaling = []
                for s in scales:                    
                    C_scaling += [s] * 15
                                            
                self.time_explosion = '12.1'    
                self.velocity_start = '10700'
                self.line_interaction = 'downbranch' 
                self.el1_scaling = {'el': 'C', 'v_start': '0', 'v_stop': '100000',
                                    'factors': C_scaling} 

            elif case == '19d_C-scaled':
                self.subdir = '11fe_19d_C-scaled/'     

                L_scal = np.arange(0.2, 1.61, 0.1)
                self.luminosity = [self.lum2loglum(3.5e9 * l) for l in L_scal]                
                self.luminosity = self.luminosity * 12

                scales = ['0.00', '0.01', '0.02', '0.05', '0.10', '0.20',
                          '0.50', '1.00', '2.00', '5.00', '10.00', '20.00']

                C_scaling = []
                for s in scales:                    
                    C_scaling += [s] * 15
                                            
                self.time_explosion = '19.1'    
                self.velocity_start = '7850'
                self.line_interaction = 'downbranch' 
                self.el1_scaling = {'el': 'C', 'v_start': '0', 'v_stop': '100000',
                                    'factors': C_scaling} 
            
            #Test case.
            elif case == 'test_5.9d':       
                self.subdir = '11fe_test_5.9d/'     

                self.luminosity = [str(format(np.log10(0.32e9), '.3f'))] * 2
                self.time_explosion = '5.9'    
                self.velocity_start = '12400'
                self.line_interaction = 'downbranch'                
                
                #self.luminosity = ['9.544']
                #self.time_explosion = '19.1'    
                #self.velocity_start = '7850'
                #self.line_interaction = 'macroatom'
                self.el1_scaling = {'el': 'Z', 'v_start': ['15000', '15000'],
                                    'v_stop': '100000', 'factors': ['0.1', '3.0']}   
                #self.el1_scaling = {'el': 'Mg', 'v_start': ['12400', '14000'],
                #                    'v_stop': '100000', 'factors': ['0.1', '0.1']}            
                #self.el2_scaling = {'el': 'Ca', 'v_start': 12400.,
                #                   'v_stop': 1.e6, 'factors': ['5.0']} 
           
            elif case == '12d_C-scan':
                self.subdir = '11fe_12d_C-scan/'     

                self.luminosity = self.lum2loglum(2.3e9)                

                v_stop = list(np.arange(10500., 18000., 500.).astype('int').astype('str'))
                self.time_explosion = '12.1'    
                self.velocity_start = '10700'
                self.line_interaction = 'downbranch' 
                self.el1_scaling = {'el': 'C', 'v_start': '0',
                                    'v_stop': v_stop, 'factors': '0'}
                                    
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
        
                                              
