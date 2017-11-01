#!/usr/bin/env python

import os
import time   

from make_tardis_input import Make_Inputs
from run_simulations import Simulate_Spectra
from append_features import Analyse_Features
from input_pars import Input_Parameters as class_input

class Master(object):
   
    """THIS CODE MAKES THE YML FILES, RUN THE SIMULATIONS AND COMPARE
    THE OUTPUT SPECTRA.
    """
    
    def __init__(self, event, case, StoN,
                 flag_run_simulation=True, flag_compute_features=True, 
                 run_uncertainties=False, make_kromer=False,
                 plot_spectra=False, plot_abun=True, show_figs=False,
                 parallel=False, verbose=True):

        self.inputs = class_input(event=event, case=case, StoN=StoN,
                                  run_uncertainties=run_uncertainties)

        self.flag_run_simulation = flag_run_simulation
        self.flag_compute_features = flag_compute_features
        self.run_uncertainties = run_uncertainties
        self.make_kromer = make_kromer
        self.plot_spectra = plot_spectra
        self.plot_abun = plot_abun
        self.show_figs = show_figs
        self.parallel = parallel
        self.verbose = verbose
        self.created_ymlfiles_list = None

        if self.verbose:        
            os.system('clear')
            print '****************************************************'
            print '***************** TARDIS SIMULATION ****************'         
            print '****************************************************'
            print '\n'
            
            print 'RUN SIMULATIONS------->', self.flag_run_simulation
            print 'COMPUTE FEATURES------>', self.flag_compute_features
            print 'MAKE KROMER PLOT------>', self.make_kromer
            print 'PLOT SPECTRA----->', self.plot_spectra
            print 'SHOW FIGURES----->', self.show_figs
            print '\n\n'

        #self.run_master()

    def run_master(self):
        
        #Create top level output directory.
        output_dir = './../OUTPUT_FILES/' + self.inputs.subdir
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)                        

        #Make a sub-dir and a .yml file for each simulation.
        object_make_yml = Make_Inputs(self.inputs, self.plot_abun, self.verbose)                 
        self.created_ymlfiles_list = object_make_yml.run()

        if self.flag_run_simulation:
            """Run sequential TARDIS simulations."""
            module_tardis_sim = Simulate_Spectra(
              created_ymlfiles_list=self.created_ymlfiles_list,
              make_kromer=self.make_kromer,
              run_uncertainties=self.inputs.run_uncertainties,
              smoothing_window=self.inputs.smoothing_window,
              N_MC_runs=self.inputs.N_MC_runs,
              extinction = self.inputs.extinction,
              parallel = self.parallel,
              show_figs=self.show_figs)
              
            module_tardis_sim.run_SIM()

        if self.flag_compute_features:
			Analyse_Features(
			  created_ymlfiles_list=self.created_ymlfiles_list,
              run_uncertainties=self.run_uncertainties,
			  smoothing_window=self.inputs.smoothing_window,
              N_MC_runs=self.inputs.N_MC_runs,
              plot_spectra = self.plot_spectra,
              show_fig=self.show_figs)
              
        #This return is used by the code 'run_test_cases.py' code.
        return self.created_ymlfiles_list      				

if __name__ == '__main__':

    """Tests"""

    #Master(event='11fe', case='12d_C-scan', StoN='medium-high',
    #       flag_run_simulation=True, flag_compute_features=True,
    #       run_uncertainties=False, make_kromer=False,
    #       plot_spectra=False, show_figs=False, parallel=False,
    #       verbose=True).run_master()  

    #Master(event='fast', case='single', StoN='low',
    #       flag_run_simulation=True, flag_compute_features=True,
    #       run_uncertainties=False, make_kromer=False,
    #       plot_spectra=False, show_figs=False, parallel=False,
    #       verbose=True).run_master()  

    Master(event='11fe', case='test_5.9d', StoN='low',
           flag_run_simulation=False, flag_compute_features=False,
           run_uncertainties=False, make_kromer=False,
           plot_spectra=False, plot_abun=True,
           show_figs=True, parallel=False,
           verbose=True).run_master()  

    #Master(event='11fe', case='19d_C-scaled', StoN='medium-high',
    #       flag_run_simulation=True, flag_compute_features=True,
    #       run_uncertainties=True, make_kromer=False,
    #       plot_spectra=True, show_figs=False, parallel=False,
    #       verbose=True).run_master()  

    #Master(event='05bl', case='12d_C-scaled_v0', StoN='low',
    #       flag_run_simulation=True, flag_compute_features=True,
    #       run_uncertainties=True, make_kromer=True,
    #       plot_spectra=True, show_figs=False, parallel=False,
    #       verbose=True).run_master() 

    #Master(event='05bl', case='12d_2D-grid_v0', StoN='low',
    #       flag_run_simulation=True, flag_compute_features=True,
    #       run_uncertainties=True, make_kromer=True,
    #       plot_spectra=True, show_figs=False, parallel=False,
    #       verbose=True).run_master()  





    #Master(event='11fe', case='12d_2D-grid_v13400', StoN='low',
    #       flag_run_simulation=True, flag_compute_features=True,
    #       run_uncertainties=True, make_kromer=False,
    #       plot_spectra=True, show_figs=False, parallel=False,
    #       verbose=True).run_master()  





    #Master(event='05bl', case='test', StoN='low',
    #       flag_run_simulation=True, flag_compute_features=False,
    #       run_uncertainties=False, make_kromer=False,
    #       plot_spectra=False, show_figs=True, parallel=False,
    #       verbose=True).run_master() 


    
    """Runs used in the 11fe-05bl paper"""
    #Master(event='11fe', case='default_L-scaled', StoN='high', flag_run_simulation=False,
    #        flag_compute_features=True, run_uncertainties=True, make_kromer=False,
    #       plot_spectra=False, show_figs=False, parallel=True,
    #       verbose=False).run_master()         

    #Master(event='11fe', case='L-grid', StoN='high', flag_run_simulation=False,
    #       flag_compute_features=True, run_uncertainties=True, make_kromer=False,
    #       plot_spectra=False, show_figs=False, parallel=True,
    #       verbose=False).run_master()         

    #Master(event='11fe', case='Ti-grid', StoN='high', flag_run_simulation=False,
    #       flag_compute_features=True, run_uncertainties=True, make_kromer=False,
    #       plot_spectra=False, show_figs=False, parallel=True,
    #       verbose=False).run_master()         

    #Master(event='11fe', case='L_Fe-grid_12d', StoN='medium', flag_run_simulation=True,
    #       flag_compute_features=True, run_uncertainties=True, make_kromer=True,
    #       plot_spectra=True, show_figs=False, parallel=False,
    #       verbose=False).run_master()  

    #Master(event='05bl', case='default_L-scaled', StoN='high', flag_run_simulation=False,
    #       flag_compute_features=True, run_uncertainties=True, make_kromer=False,
    #       plot_spectra=False, show_figs=False, parallel=True,
    #       verbose=False).run_master()                                    
    
    #Master(event='05bl', case='L-grid', StoN='high', flag_run_simulation=False,
    #       flag_compute_features=True, run_uncertainties=True, make_kromer=False,
    #       plot_spectra=False, show_figs=False, parallel=True,
    #       verbose=False).run_master()               

    """To be used to run multiple seeds for testing uncertainty calculation."""
    #Master(event='fast', case='multiple', StoN='high', flag_run_simulation=False,
    #       flag_compute_features=True, run_uncertainties=True, make_kromer=False,
    #       plot_spectra=False, show_figs=False, parallel=True,
    #       verbose=False).run_master()            

       
