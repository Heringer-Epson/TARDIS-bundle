#!/usr/bin/env python

import os   
import inspect                                                            
from shutil import copyfile

import make_yml as module_make_yml
import run_simulations as module_run_simulations
import make_structure_files as module_make_structure_files

from make_structure_files import make_structure as class_make_structure
from make_yml import make_yml as class_make_yml
from run_simulations import simulate_spectra as class_simulate_spectra

"""Uncomment input file to be imported"""
#from input_pars_11fe import input_parameters as class_input
from input_pars_05bl import input_parameters as class_input
#from input_pars_04eo import input_parameters as class_input
#from input_pars_fast import input_parameters as class_input

##########################  CALLING CLASSES  ###########################

class master(class_input):
    """
    THIS CODE MAKES THE YML FILES, RUN THE SIMULATIONS AND COMPARE THE
    OUTPUT SPECTRA.
    """
    def __init__(self, flag_make_structure=True, flag_make_yml=True,
                 flag_run_simulation=True, flag_display_interface=False,
                 verbose=True):

        class_input.__init__(self)

        self.flag_make_structure = flag_make_structure
        self.flag_make_yml = flag_make_yml
        self.flag_run_simulation = flag_run_simulation
        self.flag_display_interface = flag_display_interface
        self.verbose = verbose
        #Never create structure file if a fast run is required, as the run adopts a homogenous model.
        #This prevents the run from crashing if flag_make_structure is passed as 'True' by accident.
        if self.input_file == 'input_pars_fast':
            self.flag_make_structure = False

        self.created_ymlfiles_list = None
        self.tardis_sim = None

        if self.verbose:        
            os.system('clear')
            print '*******************************************************************************************************************************************'
            print '************************************************************ TARDIS SIMULATION ************************************************************'         
            print '*******************************************************************************************************************************************'
            print '\n\n'
            
            print 'MAKE STRUCTURE FILES---->', self.flag_make_structure
            print 'MAKE YML FILES---------->', self.flag_make_yml
            print 'RUN SIMULATIONS--------->', self.flag_run_simulation
            print 'MAKE KROMER PLOT-------->', self.make_kromer
            print 'DISPLAY INTERFACE------->', self.flag_display_interface
            print '\n\n'

        self.run_master()

    def run_master(self):

        #Only make structure files if requested. Structure files can be common for multiple simulations.
        if self.flag_make_structure:
            class_make_structure(abundance_dict=self.abun,velocity_array=self.velocity_array, filename=self.filename_structure, es=self.energy_ratio_scaling, ms=self.mass_ratio_scaling, t_exp=self.time_explosion,
                                 pass_density_as=self.pass_density_as, time_0=self.time_0, rho_0=self.rho_0, v_0=self.v_0, density_array_given=self.density_array)                  

        #Run yml code. This will create the yml files if self.flag_make_yml is True, or simply try to use pre-existing
        #files if False. In either case, the list with the names of the yml files built from the input parameters is returned.                                          
        if self.structure_type == 'file':
            object_make_yml = class_make_yml(mode=self.mode, write_files=self.flag_make_yml, log_lum=self.luminosity, time_explosion=self.time_explosion, velocity_start=self.velocity_start, velocity_stop=self.velocity_stop, velocity_array=self.velocity_array,
                                             structure_type=self.structure_type, filename_structure=self.filename_structure,
                                             energy_scaling=self.energy_ratio_scaling, mass_scaling=self.mass_ratio_scaling,
                                             abundance_type=self.abundance_type,
                                             ionization=self.ionization, excitation=self.excitation, rad_rates_type=self.rad_rates_type, line_interaction=self.line_interaction,
                                             seed=self.seeds, num_packs=self.num_packs, iterations=self.iterations, last_num_packs=self.last_num_packs, num_virtual_packs=self.num_virtual_packs,
                                             folder_name=self.subdir)

        elif self.structure_type == 'specific':
            object_make_yml = class_make_yml(mode=self.mode, write_files=self.flag_make_yml, log_lum=self.luminosity, time_explosion=self.time_explosion,
                                             structure_type=self.structure_type, velocity_start=self.velocity_start, velocity_stop=self.velocity_stop, structure_num=self.structure_num,
                                             density_type=self.density_type, time_0=self.time_0, rho_0=self.rho_0, v_0=self.v_0,
                                             abundance_type=self.abundance_type, density_value=self.density_value,
                                             abun_Ni=self.abun_Ni, abun_Si=self.abun_Si,    abun_Fe=self.abun_Fe, abun_Co=self.abun_Co, 
                                             abun_Ca=self.abun_Ca, abun_S=self.abun_S, abun_Mg=self.abun_Mg, abun_Na=self.abun_Na, 
                                             abun_C=self.abun_C, abun_O=self.abun_O, abun_Ti=self.abun_Ti, abun_Ar=self.abun_Ar,            
                                             ionization=self.ionization, excitation=self.excitation, rad_rates_type=self.rad_rates_type, line_interaction=self.line_interaction,
                                             seed=self.seeds, num_packs=self.num_packs, iterations=self.iterations, last_num_packs=self.last_num_packs, num_virtual_packs=self.num_virtual_packs,
                                             folder_name=self.subdir)               

        self.created_ymlfiles_list = object_make_yml.run()

        if self.flag_run_simulation:
            module_tardis_sim = class_simulate_spectra(self.subdir, created_ymlfiles_list=self.created_ymlfiles_list, make_kromer=self.make_kromer, run_uncertainties=self.run_uncertainties, smoothing_window=self.smoothing_window, N_MC_runs=self.N_MC_runs, display_interface=self.flag_display_interface)
            output_folder = module_tardis_sim.run_SIM()
            copyfile(module_input.__file__.split('/')[-1].split('.py')[0]+'.py', output_folder+module_input.__file__.split('/')[-1].split('.py')[0]+'.py')
        return None 

master_obj                                                              = master()      
