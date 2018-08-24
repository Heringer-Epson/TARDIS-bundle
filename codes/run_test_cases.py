#!/usr/bin/env python

import os                                                               
import sys
import time

path_tardis_output = os.environ['path_tardis_output']

import math
import cPickle
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import tardis.tardistools.compute_features as cp
from master_run import Master
from append_features import Analyse_Features
from input_pars import Input_Parameters as class_input
from scipy.signal import savgol_filter
from matplotlib.ticker import MultipleLocator
from astropy import units as u

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= QUICK TEST =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=                          
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#

class Quick_Test(object):
    """TBW.
    """
    
    def __init__(self):     
        yml_list = Master(event='fast', case='single', StoN='very-low',
                          flag_run_simulation=True, flag_compute_features=True,
                          run_uncertainties=False, make_kromer=False,
                          plot_spectra=False, show_figs=True,
                          verbose=True).run_master()   


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

    quantity : ~str
        Which quantity to use in 'run_seeds'. Can be 'pEW', 'depth' or
        'velocity'. Note that in 'run_sine', the analytical function is placed
        in the wavelength region of feature 7.
        
    feature_number : ~str
        Which feature to use in 'run_seeds'.  Currently only 'f6' and 'f7' are
        implemented.

    noise : ~float
       Sigma value to be used to generate noisy by np.random. Only used under
       the 'sine' test case.

    spectra_grid : ~float
       Separation (in angs) of the wavelength bins. Only used under
       the 'sine' test case. In TARDIS, the default value is ~2angs.
    """
    
    def __init__(self, test_case='sine', quantity='pEW', feature_number='f7',
                 noise=0.05, spectra_grid=2):     
        
        self.test_case = test_case
        self.q = quantity
        self.n = feature_number
        self.noise = noise
        self.spectra_grid = spectra_grid
        
        self.run_test()        
        
    def run_single(self):

        yml_list = Master(event='fast', case='single', StoN='low',
                          flag_run_simulation=True, flag_compute_features=True,
                          run_uncertainties=True, make_kromer=False,
                          plot_spectra=False, show_figs=False,
                          verbose=False).run_master()   

        file_path = yml_list[0].split('.yml')[0]
        with open(file_path + '.pkl', 'r+') as inp:
            D = cPickle.load(inp)
            print '\n\nComparison of a few quantities. The default values are'\
            ' from previous calculations made in stable versions of the code,'\
            ' when the test case - event: fast, case: single, is run with a '\
            'StoN: low.'
            
            print 'pEW_f7 - current:', D['pEW_f7'], '+/-', D['pEW_unc_f7'],\
            'flag: ', D['pEW_flag_f7']
            print 'pEW_f7 - default: 78.516112254 +/- 5.47670163015',\
            'flag: False\n'
            
            print 'depth_f7 - current:', D['depth_f7'], '+/-', D['depth_unc_f7'],\
            'flag:', D['depth_flag_f7']
            print 'depth_f7 - default: 0.506081861312 +/- 0.0457031793996',\
            'flag: False\n'
            
            print 'velocity_f7 - current:', D['velocity_f7'], '+/-',\
            D['velocity_unc_f7'], 'flag: ', D['velocity_flag_f7']
            
            print 'velocity_f7 - default: -11.9338201259 +/- 0.432073982855',\
            'flag: False\n'
            
    def run_sine(self):
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
                             run_uncertainties=False)
            
        #Steps of 2 (angs) are similar to the default TARDIS config, where
        #the spectra starts at 500angs, ends at 20000angs, with 10000 steps.            
        #wavelength = np.arange(6000., 6300., 0.01)
        wavelength = np.arange(6000., 6300., self.spectra_grid)
        x = np.linspace(0., 3. * np.pi, len(wavelength))           
        flux = np.sin(x)
        flux_noise = np.random.normal(np.sin(x), self.noise)
                
        #Compute pEW when flux has no noise.
        D = cp.Analyse_Spectra(
          wavelength=wavelength, flux=flux, redshift=0.,
          extinction=0., smoothing_window=inputs.smoothing_window).run_analysis()
        
        #Compute pEW and its uncertainty when flux is noisy.
        D_noise = cp.Analyse_Spectra(
          wavelength=wavelength, flux=flux_noise, redshift=0.,
          extinction=0., smoothing_window=inputs.smoothing_window).run_analysis()
       
        D_noise = cp.Compute_Uncertainty(
              D=D_noise, smoothing_window=inputs.smoothing_window,
              N_MC_runs=inputs.N_MC_runs).run_uncertainties()             

        #Create mock spectra to get the "true" uncertainty of the noisy
        #flux. To do so, create another 1000 mock noisy spectra
        mock_pEW = []
        mock_spectra = [np.random.normal(flux, self.noise) for i in range(100)]
        
        for flux_mock in mock_spectra:
              
            D_mock = cp.Analyse_Spectra(
              wavelength=wavelength, flux=flux_mock, redshift=0., extinction=0.,            
              smoothing_window=inputs.smoothing_window).run_analysis()

            mock_pEW.append(D_mock['pEW_f7'])                         
        mock_pEW = np.asarray(mock_pEW)
        mock_pEW = mock_pEW[~np.isnan(mock_pEW)]

        pEW_nonoise = str(format(D['pEW_f7'], '.2f'))
        pEW_noise = str(format(D_noise['pEW_f7'], '.2f'))
        pEW_unc_mock = str(format(np.std(np.asarray(mock_pEW)), '.2f'))
        pEW_unc_noise = str(format(D_noise['pEW_unc_f7'], '.2f'))
                
        print '\n\n************ RESULTS ************\n\n'
        
        print 'Using:' 
        print '  Noise = ' + str(self.noise)
        print '  Wavelength bin = ' + str(self.spectra_grid) + '\n'

        print '-->pEW:'
        print '---->Analytical: 200'
        print '---->Computed (w/o noise): ' + pEW_nonoise
        print '---->Computed (w noise): ' + pEW_noise + '\n'
        print '-->pEW uncertainty:'
        print '---->From mock data: ' + pEW_unc_mock
        print '---->Mean of MC computed: ' + pEW_unc_noise

        print '\n\n************* NOTES *************\n\n'

        print ('In order for the computed pEW to match the analytical one,\n'
          + 'it is necessary to use a fine spectra_grid (~0.01). This \n'
          + 'means that the TARDIS grid already introduces some\n'
          + 'uncertainty to the pEW. Also, the uncertainty computed\nwhen '
          + 'noise=0.05 agrees within 2% of the expected\none from the '
          + 'mock spectra.\n\n\n')

    def run_seeds(self):
        
        quantity_list = []
        unc_list = []
        
        #This is used to grab the correct directory, and to run the simulations
        #if they do not yet exist.
        inputs = class_input(event='fast', case='multiple', StoN='low',
                             run_uncertainties=False)        

        inp_dir = os.path.abspath('./../OUTPUT_FILES/' + inputs.subdir) + '/'

        #Check if directory with 'fast' runs exists; if not, run simulations.
        if not os.path.exists(inp_dir):
            Master(event='fast', case='multiple', StoN='high',
              flag_run_simulation=False, flag_compute_features=True,
              run_uncertainties=True, make_kromer=False,
              flag_display_interface=False, verbose=False).run_master()


        #Grab wanted quantities by looping over files in the input directory.
        for subdirs, dirs, files in os.walk(inp_dir):
            for subdir in dirs:
                subdir_fullpath = inp_dir + subdir + '/'
                for fname in os.listdir(subdir_fullpath): 
                    if fname.endswith('.pkl'):
                        with open(subdir_fullpath + fname) as inp:
                            pkl = cPickle.load(inp) 
                            quantity_list.append(pkl[self.q + '_' + self.n])
                            unc_list.append(pkl[self.q + '_unc_' + self.n])
                            
        #Mean of the measured quantity.
        quantity_mean = str(format(np.mean(quantity_list), '.2f'))
        
        #Uncertainty derived by computing the standard deviation of the same
        #quantity given multiple identical runs.
        quantity_unc_seeds = str(format(np.std(quantity_list), '.2f'))
        
        #Uncertainty computed as the mean (or meadian) of the uncertainty 
        #estiamted for each of the runs using the MC approach.
        quantity_unc_MC = str(format(np.mean(unc_list), '.2f'))
        
        #Standard deviation of the MC uncertainties.
        quantity_unc_std = str(format(np.std(unc_list), '.2f'))        
        
        print '\n\n************ RESULTS ************\n\n'
        
        print 'Using:'
        print '  Quantity is ' + self.q + '_' + self.n 
        print '  Number of simulations = ' + str(len(quantity_list)) + '\n'

        print '-->Quantity value:'
        print '---->Quantity"s mean value from multiple runs: ' + quantity_mean + '\n'     

        print '-->Uncertainty from:'
        print '---->Standard deviation of the quantity values from multiple runs: ' + quantity_unc_seeds        
        print '---->Mean of the uncertainties from multiple MC run: ' + quantity_unc_MC + '\n'        
        
        print '---->Standard deviation of the uncertainties from MC runs: ' + quantity_unc_std        
 
        print '\n\n************* NOTES *************\n\n'               

        print ('Analysis based on runs adopting the "high" S/N option.\n'
          + 'The uncertainty in the pEW features are consistent between\n'
          + 'the calcultion using MC runs (which exhibit no more than a\n'
          + 'reasonable spread) and different seeds. Depth and velocity\n'
          + 'calculation show worse agreement. No need to correct\n'
          + 'for 20% in the pEW uncertainty, as claimed in paper I.\n\n\n')
          
    def run_kromer(self):

        #Master(event='fast', case='single', StoN='very-low',
        #       flag_run_simulation=True, flag_compute_features=False,
        #       run_uncertainties=False, make_kromer=True,
        #       plot_spectra=False, show_figs=True,
        #       verbose=False).run_master()   

        Master(event='11fe', case='test', StoN='very-low',
               flag_run_simulation=True, flag_compute_features=False,
               run_uncertainties=False, make_kromer=True,
               plot_spectra=False, show_figs=True, verbose=False).run_master()
                       
    def run_test(self):
        if self.test_case == 'single':
            self.run_single()            
        if self.test_case == 'sine':
            self.run_sine()       
        if self.test_case == 'seeds':
            self.run_seeds()
        if self.test_case == 'kromer':
            self.run_kromer()    

class Compute_Smoothing_Factor(object):
    
    def __init__(self):
        
        #self.smoothing_windows = [17, 19, 21, 51]
        self.smoothing_windows = [51]
        self.noises = [0.001, 0.01, 0.05, 0.1]
        self.wavelength = None
        self.flux = None
        
        self.run_smoothing_factor()
        
    def make_mock_spectrum(self):
        self.wavelength = np.arange(6000., 6300., 2.)
        x = np.linspace(0., 3. * np.pi, len(self.wavelength))           
        self.flux = np.sin(x)

    def compute_factors(self):

        def rms(y_data, y_smot):
            #Given a noisy and a smoothed data, compute an array of the
            #squared differences and take the square-root of its mean.
            #Used as a proxy of the noise.
            rms_out = np.sqrt(((y_data - y_smot)**2.).mean())
            if rms_out < 1.e-10: rms_out = 1.e-5     
            return rms_out
        
        w_window = ((self.wavelength >= 6100.) & (self.wavelength <= 6200.))
        
        for noise in self.noises:
            flux_noise = np.random.normal(self.flux, noise)                        
            for sw in self.smoothing_windows:
                
                flux_smoothed = savgol_filter(flux_noise, sw, 3)
                estimated_noise = rms(flux_noise[w_window], flux_smoothed[w_window])
                
                
                correction_factor = noise / estimated_noise 
                
                print sw, noise, estimated_noise, correction_factor
                
                
        

    def run_smoothing_factor(self):
        self.make_mock_spectrum()
        self.compute_factors()
   
class Plot_11fe_Sequence(object):
    """Whenever TARDIS is updated to a new version, one can call this class
    to run and plot the spectra of 11fe at maximum, at the correct luminosity
    and a cooler version. Visual inspection of the plot serves a basis to 
    check whether the update is behaving accordingly.
    """
    
    def __init__(self, show_fig=True, save_fig=False):

        self.show_fig = show_fig
        self.save_fig = save_fig 
        self.fig, self.ax = plt.subplots(figsize=(10.,18))

        self.obs_list = [
          '2011_08_25', '2011_08_28', '2011_08_31', '2011_09_03',
          '2011_09_07', '2011_09_10', '2011_09_13', '2011_09_19']
        
        self.syn_list = [
          ['7.903', '13300', '3.7'], ['8.505', '12400', '5.9'],
          ['9.041', '11300', '9.0'], ['9.362', '10700', '12.1'],
          ['9.505', '9000', '16.1'], ['9.544', '7850', '19.1'],
          ['9.505', '6700', '22.4'], ['9.362', '4550', '28.3']]  

        self.make_plot()     
        
    def set_fig_frame(self):
        """Define the configuration of the figure axes."""
        x_label = r'$\mathrm{rest \ wavelength} \ \mathrm{[\AA}]}$'
        y_label = r'$\mathrm{Relative \ f}_{\lambda}$'
        
        fs = 26.
        self.ax.set_xlabel(x_label,fontsize=fs)
        self.ax.set_ylabel(y_label,fontsize=fs)
        self.ax.set_xlim(1500., 10000.)
        self.ax.set_ylim(-5.5, 12.)      
        self.ax.tick_params(axis='y', which='major', labelsize=fs, pad=8)      
        self.ax.tick_params(axis='x', which='major', labelsize=fs, pad=8)
        self.ax.minorticks_on()
        self.ax.tick_params('both', length=8, width=1, which='major')
        self.ax.tick_params('both', length=4, width=1, which='minor')
        self.ax.xaxis.set_minor_locator(MultipleLocator(500.))
        self.ax.xaxis.set_major_locator(MultipleLocator(1000.))
        self.ax.yaxis.set_minor_locator(MultipleLocator(0.5))
        self.ax.yaxis.set_major_locator(MultipleLocator(2.))        
        self.ax.tick_params(labelleft='off')  

    def load_and_plot_observational_spectrum(self):
        
        directory = '/home/heringer/Research/routines_11fe-05bl/INPUT_FILES/'\
                    + 'observational_spectra/2011fe/'
        
        offset = 9.
        for date in self.obs_list:
            with open(directory + date + '.pkl', 'r') as inp:
                pkl = cPickle.load(inp)
                self.ax.plot(
                  pkl['wavelength_corr'], pkl['flux_smoothed'] + offset,
                  color='k', ls='-', lw=3., zorder=2.)
            offset -= 2.  
                    
    def load_and_plot_standard_synthetic_spectrum(self):
        
        def make_fpath(L, v, t):
            fname = ('velocity_start-' + v + '_loglum-' + L
                     + '_line_interaction-downbranch_time_explosion-' + t)
            return (path_tardis_output + '11fe_default_L-scaled/' + fname
                    + '/' + fname + '.pkl')
        
        syn_fpath_list = [make_fpath(s[0], s[1], s[2]) for s in self.syn_list]    
                    
        offset = 9.
        for fpath in syn_fpath_list:
            with open(fpath, 'r') as inp:
                pkl = cPickle.load(inp)
                self.ax.plot(
                  pkl['wavelength_corr'], pkl['flux_smoothed'] + offset,
                  color='b', ls='-', lw=3., zorder=2.)
            offset -= 2.                      

    def load_and_plot_updated_synthetic_spectrum(self):
        
        def make_fpath(L, v, t):
            fname = ('velocity_start-' + v + '_loglum-' + L
                     + '_line_interaction-downbranch_time_explosion-' + t)
            return (path_tardis_output + '11fe_default/' + fname
                    + '/' + fname + '.pkl')
        
        syn_fpath_list = [make_fpath(s[0], s[1], s[2]) for s in self.syn_list]    
                    
        offset = 9.
        for fpath in syn_fpath_list:
            with open(fpath, 'r') as inp:
                pkl = cPickle.load(inp)
                self.ax.plot(
                  pkl['wavelength_corr'], pkl['flux_smoothed'] + offset,
                  color='r', ls='-', lw=3., zorder=2.)
            offset -= 2.                      


    def add_text(self):
        offset = 9.
        for s in self.syn_list:
            self.ax.text(1800., offset + 0.4, s[2] + ' d',
                         fontsize=20., horizontalalignment='left')        
            offset -= 2.                      
                    
    def save_figure(self):        
        if self.save_fig:
            directory = './../OUTPUT_FILES/FIGURES/'
            plt.savefig(directory + 'Fig_11fe_sequence.png',
                        format='png', dpi=360)

    def make_plot(self):
        self.set_fig_frame()
        self.load_and_plot_observational_spectrum()
        self.load_and_plot_standard_synthetic_spectrum()
        #self.load_and_plot_updated_synthetic_spectrum()
        self.add_text()
        plt.grid(False)
        plt.tight_layout()
        self.save_figure()
        if self.show_fig:
            plt.show()
        plt.close()    

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= S/N TEST =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=                          
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#

class SN_Test(object):
    """THIS CODE READS AND PLOT THE SPECTRA FROM 3 SIMULATIONS RAN WITH
    DIFFERENT OPTIONS OF S/N.
    """
    
    def __init__(self, show_fig=True, save_fig=True):     

        self.show_fig = show_fig
        self.save_fig = save_fig 
                
        self.cases = ['high', 'sh1', 'sh2']
        self.labels = ['lp = 5.e5, vp = 5',
                       'lp = 1.e6, vp = 5',
                       'lp = 5.e5, vp = 10']
        self.colors = ['r', 'g', 'b']

        self.fig, self.ax = plt.subplots(figsize=(12.,8.))
        self.fs = 26.
        self.offset = -1.e-14
        self.d = (6.4e6 * u.pc).to(u.cm)
        self.run_test()        

    def set_fig_frame(self):
        """Define the configuration of the figure axes."""
        x_label = r'$\mathrm{rest \ wavelength} \ \mathrm{[\AA}]}$'
        y_label = r'$\mathrm{Relative \ f}_{\lambda}$'
        
        fs = 26.
        self.ax.set_xlabel(x_label,fontsize=fs)
        self.ax.set_ylabel(y_label,fontsize=fs)
        self.ax.set_xlim(1500., 10000.)
        #self.ax.set_ylim(-5.5, 12.)      
        self.ax.tick_params(axis='y', which='major', labelsize=fs, pad=8)      
        self.ax.tick_params(axis='x', which='major', labelsize=fs, pad=8)
        self.ax.minorticks_on()
        self.ax.tick_params('both', length=8, width=1, which='major')
        self.ax.tick_params('both', length=4, width=1, which='minor')
        self.ax.xaxis.set_minor_locator(MultipleLocator(500.))
        self.ax.xaxis.set_major_locator(MultipleLocator(1000.))
        #self.ax.yaxis.set_minor_locator(MultipleLocator(0.5))
        #self.ax.yaxis.set_major_locator(MultipleLocator(2.))        
        self.ax.tick_params(labelleft='off')  

    def load_and_plot_observational_spectrum(self):
        
        fpath = '/home/heringer/Research/routines_11fe-05bl/INPUT_FILES/'\
                    + 'observational_spectra/2011fe/2011_08_25.pkl'
        
        with open(fpath, 'r') as inp:
            pkl = cPickle.load(inp)
            for i in range(len(self.cases)):
                self.ax.plot(
                  pkl['wavelength_corr'], pkl['flux_raw'] + self.offset * i,
                  color='k', ls='-', lw=3., zorder=2.)
        
    def load_and_plot_synthetic_spectrum(self):

        fname = ('velocity_start-13300_loglum-7.903_line_interaction-'
                 + 'downbranch_time_explosion-3.7')
        
        for i, case in enumerate(self.cases):
            fpath =  (path_tardis_output + '11fe_default_' + case + '/'
                      + fname + '/' + fname + '.pkl')
            with open(fpath, 'r') as inp:
                pkl = cPickle.load(inp)
                
                distance_factor = 4. * np.pi * self.d**2.
                flux = pkl['flux_raw'] / distance_factor.value  + self.offset * i
                
                self.ax.plot(
                  pkl['wavelength_corr'], flux,
                  color=self.colors[i], ls='-', lw=2., zorder=1.,
                  label=self.labels[i])  

        self.ax.legend(
          frameon=False, fontsize=self.fs, numpoints=1, ncol=1,
          labelspacing=0.05, handletextpad=0., loc='best')  

    def save_figure(self):        
        if self.save_fig:
            directory = './../OUTPUT_FILES/FIGURES/'
            plt.savefig(directory + 'Fig_11fe_S-to-N_comparison.png',
                        format='png', dpi=360)

    def run_test(self):
        self.set_fig_frame()
        self.load_and_plot_observational_spectrum()
        self.load_and_plot_synthetic_spectrum()
        plt.grid(False)
        plt.tight_layout()
        self.save_figure()
        if self.show_fig:
            plt.show()
        plt.close()    

    
if __name__ == '__main__':

    Quick_Test()

    #SN_Test()
    
    #Plot_11fe_Sequence()
    
    #Feature_Test(test_case='single')
    #Feature_Test(test_case='kromer')
    #Feature_Test(test_case='sine', quantity='pEW', feature_number='f7',
    #              noise=0.20, spectra_grid=2)

    #Feature_Test(test_case='seeds', quantity='pEW', feature_number='f7')
    #Feature_Test(test_case='seeds', quantity='velocity', feature_number='f7')
    #Feature_Test(test_case='seeds', quantity='depth', feature_number='f7')

    #Feature_Test(test_case='seeds', quantity='pEW', feature_number='f6')
    #Feature_Test(test_case='seeds', quantity='velocity', feature_number='f6')
    #Feature_Test(test_case='seeds', quantity='depth', feature_number='f6')

    #Compute_Smoothing_Factor()
