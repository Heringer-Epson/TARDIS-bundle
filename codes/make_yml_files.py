#!/usr/bin/env python

import os                                                               
import sys
import time
import shutil

import numpy as np

                                  
class Make_Yml(object):

    """
    Code Description
    ----------

    THIS CODE MAKES A SERIES OF YML FILES ACCORDING TO THE INPUT
    PARAMETERS PROVIDED HERE.

    THE NON-DEFAULT INPUT PARAMETERS *MUST* BE PASSED AS LISTS, EVEN IF
    THE LIST CONTAINS A SINGLE ELEMENT. EXCEPTIONS TO THIS ARE THE
    PARAMETERS THAT DO *NOT* BELONG TO THE MASTER DICTIONARY.

    THE CODE WILL LOOP THROUGH ALL THE INPUT LISTS.
    ONE YML FILE IS CREATED PER INPUT COMBINATION.

    Output
    ------

    run_SIM:
        One .yml file per combination of parameters.

    Returns
    -------

    created_ymlfiles_list:
        List contaning the string of the created .yml files.

    Notes
    -----
    1)  The atomic data used should be 'kurucz_cd23_chianti_H_He.h5'.
        It is not available with the tardis developer's version in
        GitHub and should be taken from the tardis_example.
    """


    def __init__(self,
      write_files=True, log_lum='9.4472', time_explosion='19',
      distance='0', atomic_data='kurucz_cd23_chianti_H_He.h5',
      structure_type='specific', filename_structure='None',
      velocity_start='1.1e+4', velocity_stop='2.0e+4',
      velocity_array=None, structure_num='20',
      density_type='branch85_w7', time_0='2.', rho_0='6.e-10',
      v_0='3000.', exponent='1.0', density_value='1.e-40',
      abundance_type='uniform', abun_O='0.00', abun_C='0.00',
      abun_Na='0.00', abun_Mg='0.00', abun_Si='0.00', abun_S='0.00',
      abun_Ar='0.00', abun_Ca='0.00', abun_Ti='0.00', abun_Fe='0.00',
      abun_Co='0.00', abun_Ni='0.00', energy_scaling='1.',
      mass_scaling='1', ionization='nebular', excitation='dilute-lte',
      rad_rates_type='dilute-blackbody', line_interaction='macroatom',
      seed='23111963', num_packs='2.0e+5', iterations='20',
      last_num_packs='5.0e+5', num_virtual_packs='5', spec_start='500',
      spec_stop='20000', spec_num='10000', folder_name=False,
      copy_this_code=True, subdir_safety=False, clean_subdir=True,
      verbose=True):
        
        self.MASTER = {}
        
        self.write_files = write_files
        
        self.MASTER['loglum'] = log_lum
        self.MASTER['time_explosion'] = time_explosion
        self.MASTER['distance'] = distance

        self.MASTER['atomic_data'] = atomic_data
                
        self.MASTER['structure_type'] = structure_type
        self.MASTER['velocity_start'] = velocity_start
        self.MASTER['velocity_stop'] = velocity_stop
        self.MASTER['structure_num'] = structure_num
        self.MASTER['density_type'] = density_type

        self.MASTER['es'] = energy_scaling
        self.MASTER['ms'] = mass_scaling

        self.MASTER['time_0'] = time_0
        self.MASTER['rho_0'] = rho_0
        self.MASTER['v_0'] = v_0
        self.MASTER['exponent'] = exponent
        self.MASTER['density_value'] = density_value
        
        
        self.MASTER['abundance_type'] = abundance_type
        self.MASTER['abun_O'] = abun_O
        self.MASTER['abun_C'] = abun_C
        self.MASTER['abun_Na'] = abun_Na
        self.MASTER['abun_Mg'] = abun_Mg
        self.MASTER['abun_Si'] = abun_Si
        self.MASTER['abun_S'] = abun_S
        self.MASTER['abun_Ar'] = abun_Ar
        self.MASTER['abun_Ca'] = abun_Ca
        self.MASTER['abun_Ti'] = abun_Ti
        self.MASTER['abun_Fe'] = abun_Fe
        self.MASTER['abun_Co'] = abun_Co
        self.MASTER['abun_Ni'] = abun_Ni
        
        self.MASTER['ionization'] = ionization
        self.MASTER['excitation'] = excitation
        self.MASTER['rad_rates_type'] = rad_rates_type
        self.MASTER['line_interaction'] = line_interaction
        
        self.MASTER['seed'] = seed
        self.MASTER['num_packs'] = num_packs
        self.MASTER['iterations'] = iterations
        self.MASTER['last_num_packs'] = last_num_packs
        self.MASTER['num_virtual_packs'] = num_virtual_packs
                
        self.MASTER['spec_start'] = spec_start
        self.MASTER['spec_stop'] = spec_stop
        self.MASTER['spec_num'] = spec_num      
        
        self.filename_structure = filename_structure
        self.velocity_array = velocity_array
        self.folder_name = folder_name
        self.subdir_safety = subdir_safety
        self.clean_subdir = clean_subdir
        self.verbose = verbose
                                
        self.default_pars = []
        self.num_files = 1
        self.non_default_pars = []    
        self.subdir_fullpath = ('./../INPUT_FILES/YML_FILES/'
                                +self.folder_name)
        self.simulation_list = []
        self.created_ymlfiles_list = []
        self.start_time = time.time()

    def print_run_start(self):
        if self.verbose:
            print '----------------------------------------------------'
            print '--------------- Making .yml files ------------------'
            print '----------------------------------------------------'
            print '\n'
        return None
        
    def collect_non_default_pars(self):
        """ Organise variables so that .yml files are created according
        to non-default parameters (passed as lists.) 
        """
        for key in self.MASTER:         
            if isinstance(self.MASTER[key], str):
                self.MASTER[key] = [(self.MASTER[key])]
            elif isinstance(self.MASTER[key], list):
                self.non_default_pars.append(key)
            else:
                sys.exit('Error: Input variable '+key+' is neither a '
                          +'list nor a string.')   
        self.num_files = len(self.MASTER[self.non_default_pars[0]])
                          
        if self.verbose:
            print 'INPUT PARAMETERS (non-default):'
            for key in self.non_default_pars:
                print '    '+key+':', self.MASTER[key]
            print '\nSTATUS:'   
            print ('    '+str(self.num_files)+' files will be created at'
                   +self.subdir_fullpath+'\n') 
        return None  

    def make_outfolder(self):
        """Set outfolder name.
        """
        if not isinstance(self.folder_name, str):
            self.folder_name = 'Default/'
            for par in self.non_default_pars:
                self.folder_name = self.folder_name[:-1]+'_AND_'+par
        
        """Create the directory if it doesn't already exist.
        """
        outfolder = self.subdir_fullpath
        if os.path.exists(outfolder):
            if self.subdir_safety:
                sys.exit('Error: Output directory already exists.')
            else:
                if self.clean_subdir:
                    shutil.rmtree(outfolder)
                    os.makedirs(outfolder)
        else:
            os.makedirs(outfolder)  
        return None                     

    def make_combination_of_parameters(self):
        """ 'self.simulation_list' is a master list where each entry is
        another list containing all the paramters that defined a
        similation. Prefixes are added (such as loglum=X) because
        the loop through the keys is not ordenate.
        """
        list_of_all_parameters = []
        for par in self.MASTER.keys():
            list_of_all_parameters.append([par+'|'+value
                                          for value in self.MASTER[par]])
            
        for k in range(self.num_files):
            self.simulation_list.append([])
            for par in list_of_all_parameters:
                if len(par) > 1:
                    self.simulation_list[k].append(par[k])
                else:   
                    self.simulation_list[k].append(par[0])
        return None

    def filename_of_simulation(self, list_of_pars):
        """ Returns a 'standardised' file name according to the
        input parameters of the simulation. 
        """
        filename = ''
        for entry in list_of_pars:
            variable, value = entry.split('|')
            if variable in self.non_default_pars:
                if variable[0:5] != 'abun_':
                    filename += variable+'-'+value+'_'
                else:
                    if abs(float(value)) > 0.0001:  
                        filename += variable+'-'+value+'_'
        filename = filename[:-1]
        filename  += '.yml' 
        return filename

    def write_yml(self):
        """ Main function to write the output .yml files.
        """
        if self.verbose:
            print 'MAKING FILES:'
        
        for k, simulation in enumerate(self.simulation_list):
            filename = self.filename_of_simulation(simulation)
            ymlfile_fullpath = self.subdir_fullpath+filename         
            yml_file = open(ymlfile_fullpath  ,'w')
            self.created_ymlfiles_list.append(filename)
            
            if self.write_files:
                if self.verbose:
                    print '    CREATED: '+ymlfile_fullpath
                
                """Organize each simulation as a dictionary.
                i.e. unfold the simulation parameters (with prefix)
                into a dictionary.
                """
                PARS = {}    
                for entry in simulation:
                    variable, value = entry.split('|')
                    PARS[variable] = value
                
                """Write the .yml file.
                """
                yml_file.write('tardis_config_version: v1.0\n')
                yml_file.write('supernova:\n')
                yml_file.write('    luminosity_requested: '
                                    +PARS['loglum']+' log_lsun\n')
                yml_file.write('    time_explosion: '
                                    +PARS['time_explosion']+' day\n')
                if PARS['distance'] != '0':
                    yml_file.write('distance: '
                                   +PARS['distance']+' Mpc\n')  
                yml_file.write('\n')
                yml_file.write('atom_data: '+PARS['atomic_data']+'\n')
                yml_file.write('\n')
                yml_file.write('model:\n')
                yml_file.write('\n')
                yml_file.write('    structure:\n')
                yml_file.write('        type: '
                               +PARS['structure_type']+'\n')
                if PARS['structure_type'] == 'specific':
                    yml_file.write('        velocity:\n')
                    yml_file.write('            start: '
                                   +PARS['velocity_start']+' km/s\n')
                    yml_file.write('            stop: '
                                   +PARS['velocity_stop']+' km/s\n')
                    yml_file.write('            num: '
                                   +PARS['structure_num']+'\n')
                    yml_file.write('\n')
                    yml_file.write('        density:\n')
                    yml_file.write('            type: '
                                   +PARS['density_type']+'\n')
                    if PARS['density_type'] == 'exponential':
                        yml_file.write('            time_0: '
                                       +PARS['time_0']+' day\n')
                        yml_file.write('            rho_0: '
                                       +PARS['rho_0']+' g/cm^3\n')
                        yml_file.write('            v_0: '
                                       +PARS['v_0']+' km/s\n')
                        yml_file.write('            exponent: '
                                       +PARS['exponent']+'\n')
                    elif PARS['density_type'] == 'uniform':
                        yml_file.write('            value: '
                                       +PARS['density_value']+' g/cm^3\n') 
                elif PARS['structure_type'] == 'file':
                    yml_file.write('        filename: '+os.path.abspath(
                                   './../INPUT_FILES/DENSITY_FILES/')
                                   +'/density_'+self.filename_structure
                                   +'_es-'+PARS['es']+'_ms-'+PARS['ms']
                                   +'_'+PARS['time_explosion']+'_day.dat\n')                                    
                    yml_file.write('        filetype: simple_ascii\n')
                    yml_file.write('        v_inner_boundary: '
                                   +PARS['velocity_start']+' km/s\n')
                    yml_file.write('        v_outer_boundary: '
                                   +PARS['velocity_stop']+' km/s\n')                        
                yml_file.write('\n')
                yml_file.write('    abundances:\n')
                yml_file.write('        type: '
                               +PARS['abundance_type']+'\n')
                if PARS['abundance_type'] == 'uniform':
                    for entry in PARS.keys():
                        if entry[0:5] == 'abun_':
                            if float(PARS[entry]) > 1.e-6:
                                element = entry.split('abun_')[1]
                                yml_file.write('        '+element
                                               +': '+PARS[entry]+'\n')
                                
                elif PARS['abundance_type'] == 'file':
                    yml_file.write('        filename: '+os.path.abspath(
                                   './../INPUT_FILES/STRATIFIED_COMPOSITION_FILES/')
                                   +'/abundance_'+self.filename_structure
                                   +'_'+PARS['time_explosion']+'_day.dat\n') 
                    yml_file.write('        filetype: simple_ascii\n')              
                yml_file.write('\n')
                yml_file.write('plasma:\n')
                yml_file.write('    ionization: '+PARS['ionization']+'\n')
                yml_file.write('    excitation: '+PARS['excitation']+'\n')
                yml_file.write('    radiative_rates_type: '
                               +PARS['rad_rates_type']+'\n')
                yml_file.write('    line_interaction_type: '
                               +PARS['line_interaction']+'\n')
                yml_file.write('\n')
                yml_file.write('montecarlo:\n')
                yml_file.write('    seed: '+PARS['seed']+'\n')
                yml_file.write('    no_of_packets: '+PARS['num_packs']+'\n')
                yml_file.write('    iterations: '+PARS['iterations']+'\n')
                yml_file.write('    last_no_of_packets: '
                               +PARS['last_num_packs']+'\n')
                yml_file.write('    no_of_virtual_packets: '
                               +PARS['num_virtual_packs']+'\n')
                yml_file.write('\n')
                yml_file.write('spectrum:\n')
                yml_file.write('    start: '+PARS['spec_start']+' angstrom\n')
                yml_file.write('    stop: '+PARS['spec_stop']+' angstrom\n')
                yml_file.write('    num: '+PARS['spec_num'])

                yml_file.close()
        
            else:
                if self.verbose:
                    print '    USING: '+ymlfile_fullpath                    
            
        if self.verbose:
            print '\n*** DONE - SUCCESSFUL RUN\n\n'
        
        return None
        
    def run(self):
        self.print_run_start()
        self.collect_non_default_pars()
        self.make_outfolder()
        self.make_combination_of_parameters()
        self.write_yml()
        return self.created_ymlfiles_list
        
