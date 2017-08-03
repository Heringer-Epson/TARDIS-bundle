#!/usr/bin/env python

import os   
import inspect                                                            
from shutil import copyfile

from make_tardis_input import Make_Inputs
from run_simulations import Simulate_Spectra
from append_features import Analyse_Features
from input_pars import Input_Parameters as class_input

class Master(object):
   
    """THIS CODE MAKES THE YML FILES, RUN THE SIMULATIONS AND COMPARE
    THE OUTPUT SPECTRA.
    """
    
    def __init__(self, event, case, StoN,
                 flag_run_simulation=True, run_uncertainties=False, 
                 flag_compute_features=False, make_kromer=False,
                 flag_display_interface=False, verbose=True):

        self.inputs = class_input(event=event, case=case, StoN=StoN,
                                  run_uncertainties=run_uncertainties,
                                  make_kromer=make_kromer)

        self.flag_run_simulation = flag_run_simulation
        self.flag_display_interface = flag_display_interface
        self.flag_compute_features = flag_compute_features
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
            print 'MAKE KROMER PLOT------>', self.inputs.make_kromer
            print 'DISPLAY INTERFACE----->', self.flag_display_interface
            print '\n\n'

        #self.run_master()

    def run_master(self):
        
        #Create top level output directory.
        output_dir = './../OUTPUT_FILES/' + self.inputs.subdir
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)                        

        #Make a sub-dir and a .yml file for each simulation.
        object_make_yml = Make_Inputs(self.inputs, self.verbose)                 
        self.created_ymlfiles_list = object_make_yml.run()

        if self.flag_run_simulation:
            """Run sequential TARDIS simulations."""
            module_tardis_sim = Simulate_Spectra(
              created_ymlfiles_list=self.created_ymlfiles_list,
              make_kromer=self.inputs.make_kromer,
              run_uncertainties=self.inputs.run_uncertainties,
              smoothing_window=self.inputs.smoothing_window,
              N_MC_runs=self.inputs.N_MC_runs,
              display_interface=self.flag_display_interface)
              
            module_tardis_sim.run_SIM()

        if self.flag_compute_features:
			Analyse_Features(
			  self.subdir, created_ymlfiles_list=self.created_ymlfiles_list,
			  extinction = self.extinction,
              run_uncertainties=self.run_uncertainties,
			  smoothing_window=self.smoothing_window, N_MC_runs=self.N_MC_runs,
              show_fig=self.flag_display_interface)
              
        return self.created_ymlfiles_list      				

if __name__ == '__main__':

    Master(event='fast', case='single', StoN='low', flag_run_simulation=True,
           run_uncertainties=False, flag_compute_features=False, make_kromer=False,
           flag_display_interface=False, verbose=True).run_master()      
           
    #Master(event='11fe', case='test', StoN='high', flag_run_simulation=False,
    #       run_uncertainties=False, flag_compute_features=False, make_kromer=False,
    #       flag_display_interface=False, verbose=True).run_master()      
           

       
