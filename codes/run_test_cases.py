#!/usr/bin/env python

import os   
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import tardis.tardistools.compute_features as cp
from master_run import Master
from append_features import Analyse_Features
from input_pars import Input_Parameters as class_input



#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= INPUTS =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=                          
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#

class Input_Test(object):
    """THIS CODE USES THE MASTER CODE TO CREATE DEFAULT FILES WHICH ARE THEN
    COMPARED AGAINST PREVIOUSLY WRITTEN FILES (KNOWN TO WORK.)
    """
    
    def __init__(self, test_case='11fe'):
        
        self.test_case = test_case
        
        self._created_ymlfiles_list = None
        self._yml_orig, self._yml_new = None, None
        self._abun_orig, self._abun_new = None, None
        self._dens_orig, self._dens_new = None, None

        os.system('clear')
        print '\n\n\n'
        print '****************************************************'
        print '************ RUNNING TEST CASE FOR ' + self.test_case +' ************'         
        print '****************************************************'
        print '\n'
    
        self.run_test()
    
    def make_files(self):
        
        self._created_ymlfiles_list = Master(
          event=self.test_case, case='test', StoN='high',
          flag_run_simulation=False, run_uncertainties=False,
          flag_compute_features=False, make_kromer=False,
          flag_display_interface=False, verbose=False).run_master() 
        
    def get_dirs(self):
        
        if self.test_case == '11fe':
            top_dir_orig = os.path.abspath('./../test_cases/11fe') + '/'
            self._yml_orig = top_dir_orig + 'loglum-9.544.yml'
            self._abun_orig = (top_dir_orig + 'abundance_19.1_day.dat')
            self._dens_orig = (top_dir_orig + 'density_es-1.0_ms-1.0.dat')
                              
        elif self.test_case == '05bl':
            top_dir_orig = os.path.abspath('./../test_cases/05bl') + '/'
            self._yml_orig = (top_dir_orig + 'velocity_start-8100_loglum-8.617'
                              + '_time_explosion-12.0.yml')
            self._abun_orig = (top_dir_orig + 'abundance_es-0.7_ms-1.0_12.0'
                               + '_day.dat')
            self._dens_orig = (top_dir_orig + 'density_es-0.7_ms-1.0_12.0'
                               + '_day.dat')                              
                          
        self._yml_new = self._created_ymlfiles_list[0]
        path_new = os.path.dirname(self._yml_new)
        for fname in os.listdir(path_new):
            if fname[0:9] == 'abundance':
                self._abun_new = path_new + '/' + fname
            elif fname[0:7] == 'density':
                self._dens_new = path_new + '/' + fname        
    
    def compare_files(self):
        
        for fname, f1, f2 in zip(
                            ['ABUNDANCE', 'DENSITY', 'YML'],
                            [self._abun_orig, self._dens_orig, self._yml_orig],
                            [self._abun_new, self._dens_new, self._yml_new]
                            ):
        
            print '\n\n---->COMPARING ' + fname + ' FILES:\n\n'
            with open(f1, 'r') as inp1, open(f2, 'r') as inp2:
                for i, (line1, line2) in enumerate(zip(inp1, inp2)):
                    if line1 != line2:
                        print '-------->LINE ' + str(i + 1) +':\n'
                        print '    |---->OLD', line1.strip('\n') 
                        print '    |---->NEW', line2.strip('\n')
                        print '\n' 
    
        
    def run_test(self):
        self.make_files()
        self.get_dirs()
        self.compare_files()

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= FEATURES =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=                          
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#

class Feature_Test(object):
    """THIS CODE CHECKS WHETHER THE CALCULATION OF FEATURES SUCH AS PEW AND
    DEPTH ARE BEING DONE PROPERLY. CURRENT OPTIONS ALLOW TO COMPARE VALUES
    AGAINST THE ANALYTICAL EXPECTATION OF A SIMPLE SINOIDAL CASE OR AGAINST
    A PREVIOUSLY MEASURED SPECTRA.

    Parameters
    ----------
    test_case : ~str
        Defines which test case to run. Implemented: 'sine' or 'tardis'.

    noise : ~float
       Sigma value to be used to generate noisy by np.random. Only used under
       the 'sine' test case.

    spectra_grid : ~float
       Separation (in angs) of the wavelength bins. Only used under
       the 'sine' test case. In TARDIS, the default value is ~2angs.
    """
    
    def __init__(self, test_case='sine', noise=0.05, spectra_grid=2):     
        
        self.test_case = test_case
        self.noise = noise
        self.spectra_grid = spectra_grid
        
        self.run_test()        
        
    def make_files(self):
        
        if self.test_case == 'sine':
            #In the following, the chosen analytical function is sin(x). For
            #the feature to be computed, the signal has to be in a wavelength
            #space near 6000angs.
            
            #f(x) = sin(x), where the integration range is [pi/2, 5*pi/2].
            #w(x) := 6000 + 100*x/pi => x(w) = (w - 6000)*pi/100
            #Therefore f(w)  = f(x(w)) = sin((w-6000)*pi/100)
            #pEW = int(f_cont(w) - f(w))/f_cont(w) dw
            #f_cont(w) = 1 for all w, since the peak of sin(w) is at 1.
            #pEW = int(1 - f(w))dw = int(1 - f(w(x)))dw/dx dx
            #pEW = int(1 - sin((w-6000)*pi/100) )100/pi dx; x from pi/2 to 5pi/2
            #sin part gives 0, so pEW = 100/pi * (5pi/2 - pi/2) = 200
            #The depth is 2, given that sin oscilates from -1 to 1.
            
            #Therefore, expected pEW is 200angs and depth = 2 [flux units].
            #Also note that the max occur at w=6050 and w=6250angs.

            #class_input is called only to grab the default smoothing window
            #and N_MC_runs. Event and case variables are unimportant.
            inputs = class_input(event='fast', case='single', StoN='low',
                                 run_uncertainties=False, make_kromer=False)
                
            #Steps of 2 (angs) are similar to the default TARDIS config, where
            #the spectra starts at 500angs, ends at 20000angs, with 10000 steps.            
            #wavelength = np.arange(6000., 6300., 0.01)
            wavelength = np.arange(6000., 6300., self.spectra_grid)
            x = np.linspace(0., 3 * np.pi, len(wavelength))           
            flux = np.sin(x)
            flux_noise = np.random.normal(np.sin(x), self.noise)
            
            #Compute pEW when flux has no noise.
            D = {}
            D = pd.DataFrame({'wavelength_raw': [wavelength],
                             'flux_raw': [flux], 'host_redshift': [0.]})
            D = cp.Analyse_Spectra(
                  D, extinction = 0., smoothing_mode='savgol',
                  smoothing_window=inputs.smoothing_window,
                  verbose=True).run_analysis()

            #Compute pEW and its uncertainty when flux is noisy.
            D_noise = {}
            D_noise = pd.DataFrame({'wavelength_raw': [wavelength],
                             'flux_raw': [flux_noise], 'host_redshift': [0.]})            
            D_noise = cp.Analyse_Spectra(
                  D_noise, extinction = 0., smoothing_mode='savgol',
                  smoothing_window=inputs.smoothing_window,
                  verbose=True).run_analysis()
            D_noise = cp.Compute_Uncertainty(
                  D_noise, smoothing_mode='savgol',
                  smoothing_window=inputs.smoothing_window,
                  N_MC_runs=inputs.N_MC_runs,
                  verbose=True).run_uncertainties()             
    
            #Create mock spectra to get the "true" uncertainty of the noisy
            #flux. To do so, create another 1000 mock noisy spectra
            mock_pEW = []
            mock_spectra = [np.random.normal(flux, self.noise) for i in range(100)]
            
            for flux_mock in mock_spectra:
                D_mock = {}
                D_mock = pd.DataFrame({'wavelength_raw': [wavelength],
                             'flux_raw': [flux_mock], 'host_redshift': [0.]})
                D_mock = cp.Analyse_Spectra(
                  D_mock, extinction = 0., smoothing_mode='savgol',
                  smoothing_window=inputs.smoothing_window, verbose=False).run_analysis()             
                mock_pEW.append(D_mock['pEW_f7'].tolist()[0])                         
            
            pEW_nonoise = str(format(D['pEW_f7'].tolist()[0], '.2f'))
            pEW_noise = str(format(D_noise['pEW_f7'].tolist()[0], '.2f'))
            pEW_unc_noise = str(format(D_noise['pEW_unc_f7'].tolist()[0], '.2f'))
            pEW_unc_true = str(format(np.std(mock_pEW), '.2f'))
            
            print '\n\n************ RESULTS ************\n\n'
            
            print 'Using:' 
            print '  Noise = ' + str(self.noise)
            print '  Wavelength bin = ' + str(self.spectra_grid) + '\n'

            print '-->pEW:'
            print '---->Analytical: 200'
            print '---->Computed (w/o noise): ' + pEW_nonoise
            print '---->Computed (w noise): ' + pEW_noise + '\n'
            print '-->pEW uncertainty:'
            print '---->From mock data: ' + pEW_unc_true
            print '---->Computed: ' + pEW_unc_noise

            print '\n\n************* NOTES *************\n\n'

            print ('In order for the computed pEW to match the analytical one,\n'
              + 'it is necessary to use a fine spectra_grid (~0.01). This \n'
              + 'means that the TARDIS grid already introduces some\n'
              + 'uncertainty to the pEW. Also, the uncertainty computed\nwhen '
              + 'noise=0.05 agrees within 2% of the expected\none from the '
              + 'mock spectra.\n\n\n')
                
    def run_test(self):
        self.make_files()       
            

#Input_Test(test_case='11fe')      
#Input_Test(test_case='05bl')
Feature_Test(test_case='sine', noise=0.05, spectra_grid=2)

      
