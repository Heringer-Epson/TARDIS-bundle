#!/usr/bin/env python

import os                                                               
import sys
import shutil
import copy
import itertools                                                        
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator

class Make_Inputs(object):
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

    '''
    def __init__(self,
      write_files=True, pars_fname=None, log_lum='9.4472',
      time_explosion='19',
      distance='0', atomic_data='kurucz_cd23_chianti_H_He.h5',
      structure_type='specific', filename_structure='None',
      velocity_start='1.1e+4', velocity_stop='2.0e+4',
      velocity_array=None, density_array = None, structure_num='20',
      density_type='branch85_w7', time_0=None, rho_0='6.e-10',
      abundance_dict = None,
      v_0='3000.', exponent='1.0', density_value='1.e-40',
      abundance_type='uniform', abun_O='0.00', abun_C='0.00',
      abun_Na='0.00', abun_Mg='0.00', abun_Si='0.00', abun_S='0.00',
      abun_Ar='0.00', abun_Ca='0.00', abun_Ti='0.00', abun_Fe='0.00',
      abun_Co='0.00', abun_Ni='0.00', energy_scaling='1.',
      mass_scaling='1', TiCr_scaling='1.00', Fe_scaling='1.00',
      ionization='nebular', excitation='dilute-lte',
      rad_rates_type='dilute-blackbody', line_interaction='macroatom',
      seed='23111963', num_packs='2.0e+5', iterations='20',
      last_num_packs='5.0e+5', num_virtual_packs='5', spec_start='500',
      spec_stop='20000', spec_num='10000', folder_name=False,
      copy_this_code=True, subdir_safety=False, clean_subdir=True,
      verbose=True):
      '''
    def __init__(self, inputs, plot_abun, verbose):
          
        self.plot_abun = plot_abun
        self.verbose = verbose

        self.event = inputs.event
        self.atomic_data = inputs.atomic_data
        self.structure_type = inputs.structure_type
        self.density_type = inputs.density_type        
        self.abundance_type = inputs.abundance_type
        self.spec_start = inputs.spec_start
        self.spec_stop = inputs.spec_stop
        self.spec_num = inputs.spec_num          
        self.velocity_array = inputs.velocity_array
        self.density_array = inputs.density_array
        self.abun = inputs.abun
        self.time_0 = inputs.time_0
        self.subdir = inputs.subdir

        self.ionization = inputs.ionization
        self.excitation = inputs.excitation
        self.rad_rates_type = inputs.rad_rates_type
        
        self.elements = inputs.elements
                
        self.MASTER = {}
        
        self.MASTER['loglum'] = inputs.luminosity
        self.MASTER['temperature_requested'] = inputs.temperature_requested
        self.MASTER['time_explosion'] = inputs.time_explosion
        self.distance = inputs.distance
                
        self.MASTER['velocity_start'] = inputs.velocity_start
        self.MASTER['velocity_stop'] = inputs.velocity_stop
        self.MASTER['structure_num'] = inputs.structure_num

        #Up to two elements whose mass fraction will be be scaled/added in
        #a region of the ejecta.
        self.el1_scaling = inputs.el1_scaling
        self.el2_scaling = inputs.el2_scaling
        self.el_adding = inputs.el_adding
        self.MASTER[self.el1_scaling['el'] + '-F1'] = inputs.el1_scaling['factors']
        self.MASTER['v_start_F1'] = inputs.el1_scaling['v_start']
        self.MASTER['v_stop_F1'] = inputs.el1_scaling['v_stop']
        self.MASTER[self.el2_scaling['el'] + '-F2'] = inputs.el2_scaling['factors']
        self.MASTER['v_start_F2'] = inputs.el2_scaling['v_start']
        self.MASTER['v_stop_F2'] = inputs.el2_scaling['v_stop']
        self.MASTER[self.el_adding['el'] + '-A'] = inputs.el_adding['add']
        self.MASTER['v_start_A'] = inputs.el_adding['v_start']
        self.MASTER['v_stop_A'] = inputs.el_adding['v_stop']
                
        self.MASTER['rho_0'] = inputs.rho_0
        self.MASTER['v_0'] = inputs.v_0
        self.MASTER['exponent'] = inputs.exponent
        self.MASTER['density_value'] = inputs.density_value
        
        self.MASTER['seed'] = inputs.seeds
        self.MASTER['num_packs'] = inputs.num_packs
        self.MASTER['iterations'] = inputs.iterations
        self.MASTER['last_num_packs'] = inputs.last_num_packs
        self.MASTER['num_virtual_packs'] = inputs.num_virtual_packs
        
        self.MASTER['line_interaction'] = inputs.line_interaction
        
        self.default_pars = []
        self.num_files = 1
        self.non_default_pars = []    
        self.subdir_fullpath = ('./../OUTPUT_FILES/' + self.subdir)
        self.simulation_list = []
        self.created_ymlfiles_list = []
        self.dens_fpath = None
        self.abun_fpath = None
        
    def print_run_start(self):
        if self.verbose:
            print '----------------------------------------------------'
            print '--------------- Making .yml files ------------------'
            print '----------------------------------------------------'
            print '\n'
        
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
            print ('    '+str(self.num_files)+' .yml files will be created at '
                   + self.subdir_fullpath + '\n') 

    def make_outfolder(self):
        """Set outfolder name."""
        if not isinstance(self.subdir, str):
            self.subdir = 'Default/'
            for par in self.non_default_pars:
                self.subdir = self.subdir[:-1]+'_AND_'+par
        
        """Create the directory if it doesn't already exist."""
        outfolder = self.subdir_fullpath
        if not os.path.exists(outfolder):
            os.makedirs(outfolder)

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

    def filename_of_simulation(self, list_of_pars):
        """ Returns a 'standardised' file name according to the
        input parameters of the simulation. 
        """
        filename = ''
        for entry in list_of_pars:
            variable, value = entry.split('|')
            if variable in self.non_default_pars:
                    filename += variable+'-'+value+'_'
        filename = filename[:-1]
        return filename

    def read_structure(self, event, time):
        #For 05bl and 11fe, the published data is stored in input files.
        #For 05bl, the files are provided by Hachinger in private comm.
        #For 11fe, the files were created by the code 'make_11fe_input_file'
        #based on digitalizing the published figures 3 and 10 in
        #http://adsabs.harvard.edu/abs/2014MNRAS.439.1959M        
        if event == '05bl':
            sys.exit('05bl is currently disabled. Need uptade for auto Ni56.')
            #read_every = 2
            #t_exp2str = {'11.0': 'm6', '12.0': 'm5', '14.0': 'm3', '21.8': 'p48',
            #             '29.9': 'p129'}                
            #phase = t_exp2str[time]
            #fpath = ('./../INPUT_FILES/Hachinger_2005bl/models-05bl-w7e0.7/' 
            #         + 'SN2005bl_' + phase + '/' + 'abuplot.dat')
            #self.time_0 = str(float(time) * 24. * 3600.) + ' s'
                
        if event == '11fe':
            read_every = 1 
            fpath = ('./../INPUT_FILES/Mazzali_2011fe/ejecta_layers.dat')                   
        
        abun = {} 
        
        velocity_array, density_array = np.loadtxt(
          fpath, skiprows=2, usecols=(1, 3), unpack=True)
        density_array = 10. ** density_array        
        
        for i, el in enumerate(self.elements):
            abun[el] = np.loadtxt(fpath, skiprows=2, usecols=(4 + i,),
                                  unpack=True)
        
        return velocity_array, density_array, abun

    def scale_abun(self, inp_abun, el_scaling, scale, v_start, v_stop):
        
        element = el_scaling['el']
        #Make independent copy of the abun dictionary to prevent differentt_exp
        #from modifying an already modified (decayed, or scaled) abundundaces. 
        abun = copy.deepcopy(inp_abun)
    
        if element != 'None' and element != 'Z':
            shells = np.arange(0, len(self.velocity_array), 1)
            condition = ((self.velocity_array >= v_start)
                         & (self.velocity_array <= v_stop))
            shell_window = shells[condition]
            
            for i in shell_window:
                abundance_all = np.asarray([abun[el][i] for el in self.elements])
                abundance_all = np.nan_to_num(abundance_all)
                el_most = self.elements[abundance_all.argmax()]
                orig_abun = float(abun[element][i])
                new_abun = scale * orig_abun
                
                #Prevent layers rich in a given element (say carbon) to be
                #scaled and receive a larger than 1 value.
                if new_abun < 1.:                                
					#Update abundances.
					abun[element][i] = str(new_abun)            
					abun[el_most][i] = str(float(abun[el_most][i])
										   + (orig_abun - new_abun))        

        if element == 'Z':

            #Modified the abundance of elements with Z>20.
            el_mod = ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu',
                      'Zn', 'Fe0', 'Ni0'] 

            shells = np.arange(0, len(self.velocity_array), 1)
            condition = ((self.velocity_array >= v_start)
                         & (self.velocity_array <= v_stop))
            shell_window = shells[condition]
            
            for i in shell_window:
                abundance_all = np.asarray([abun[el][i] for el in self.elements])
                abundance_all = np.nan_to_num(abundance_all)
                el_most = self.elements[abundance_all.argmax()]
                
                changed = 0.
                for _element in el_mod:
                    orig_abun = float(abun[_element][i])
                    new_abun = scale * orig_abun
                    abun[_element][i] = str(new_abun)            
                    changed += (orig_abun - new_abun)
                    
                abun[el_most][i] = str(float(abun[el_most][i]) + changed)        

        return abun

    def add_abun(self, inp_abun, el_adding, parcel, v_start, v_stop):

        element = el_adding['el']
        #Make independent copy of the abun dictionary to prevent differentt_exp
        #from modifying an already modified (decayed, or scaled) abundundaces. 
        abun = copy.deepcopy(inp_abun)
    
        if element != 'None':
            shells = np.arange(0, len(self.velocity_array), 1)
            condition = ((self.velocity_array >= v_start)
                         & (self.velocity_array <= v_stop))
            shell_window = shells[condition]
            
            for i in shell_window:
                abundance_all = np.asarray([abun[el][i] for el in self.elements])
                abundance_all = np.nan_to_num(abundance_all)
                el_most = self.elements[abundance_all.argmax()]
                orig_abun = float(abun[element][i])
                new_abun = orig_abun + parcel
                                                
                #Update abundances.
                abun[element][i] = str(new_abun)            
                abun[el_most][i] = str(float(abun[el_most][i]) + (orig_abun - new_abun))        

        return abun            

    def make_structure_file(self, fpath, velocity_array, density_array, abun):

        N_shells = len(abun['H'])
        header1 = 't0: ' + self.time_0
        header2 = 'index velocity density'
        header3 = '1 km/s g/cm^3'
        for el in self.elements:
            header2 += ' ' + el  
            header3 += ' 1'  

        with open(fpath, 'w') as out_structure:         
            out_structure.write(header1 + ' \n')
            out_structure.write(header2 + ' \n')
            out_structure.write(header3)
            for i, (v, d) in enumerate(zip(velocity_array, density_array)):
                out_structure.write('\n' + str(i) + ' ' + str(v) + ' ' + str(d))                 
                for el in self.elements:
                    out_structure.write(' ' + str(abun[el][i]))

    def make_density_file(self, fpath, velocity_array, density_array):
        with open(fpath, 'w') as out_density:
            out_density.write(self.time_0 +'  \n')      
            out_density.write('# index velocity (km/s) density (g/cm^3)')
            for i, (v, d) in enumerate(zip(velocity_array, density_array)):
                out_density.write('\n' + str(i) + '    ' + str(v) + '    ' + str(d))                 

    def make_abundance_file(self, fpath, abun):
        with open(fpath, 'w') as out_abundance:
            out_abundance.write('# index Z=1 - Z=30')   
            N_shells = len(abun['H'])
            for i in range(N_shells):
                out_abundance.write('\n' + str(i))
                for el in self.elements:
                    out_abundance.write(' ' + str(abun[el][i]))
                    
    def make_abundance_plot(self, fpath_plot, vel, inp_abun):
        
        abun = copy.deepcopy(inp_abun)
        #Arrange figure frame.
        fs = 26.
        fig, ax = plt.subplots(figsize=(14,8))
        x_label = r'$v\ \ \rm{[km\ \ s^{-1}]}$'
        y_label = r'$\rm{mass\ \ fraction}$'
        ax.set_xlabel(x_label, fontsize=fs)
        ax.set_ylabel(y_label, fontsize=fs)
        ax.set_yscale('log')
        ax.set_xlim(2500., 23000.)
        ax.set_ylim(1.e-4, 1.1)
        ax.tick_params(axis='y', which='major', labelsize=fs, pad=8)       
        ax.tick_params(axis='x', which='major', labelsize=fs, pad=8)
        #ax.minorticks_off()
        ax.tick_params('both', length=8, width=1, which='major')
        ax.tick_params('both', length=4, width=1, which='minor')
        ax.xaxis.set_minor_locator(MultipleLocator(1000.))
        ax.xaxis.set_major_locator(MultipleLocator(5000.))  
        
        #Plot abundances.
        list_el_plot = ['C', 'O', 'Mg', 'Si', 'S', 'Ca', 'Ti', 'Fe0', 'Ni0']
        label = ['C', 'O', 'Mg', 'Si', 'S', 'Ca', 'Ti+Cr', 'Fe0', 'Ni0']
        color= ['y', 'r', 'b', 'lightgreen', 'peru', 'grey', 'g', 'purple', 'k']
        for i, el in enumerate(list_el_plot):
            x = vel
            y = abun[el]
            if el == 'Ti':
                y = np.asarray(y) * 2.
            
            ax.step(x, y, color=color[i], label=label[i], lw=3., where='post')

        ax.legend(frameon=False, fontsize=20., numpoints=1, ncol=1,
                  labelspacing=0.05, loc='best')          
        plt.tight_layout()        
        plt.savefig(fpath_plot, format='png', dpi=360)
        plt.close()
                                    
    def control_structure_files(self, spawn_dir, t_exp,
                                scale_F1, v_start_F1, v_stop_F1,
                                scale_F2, v_start_F2, v_stop_F2,
                                parcel_A, v_start_A, v_stop_A):
        """This function uses the name of the input_pars_X.py file to determine
        whether to make a density and an abundance files. Note that
        input_pars_Hach.py has its own routine to read the data provided by
        Hachinger and create structure files readable by TARDIS.
        """

        self.dens_fpath = spawn_dir + 'density_' + t_exp + '_day.dat'
        self.abun_fpath = spawn_dir + 'abundance_' + t_exp + '_day.dat'
        #self.model_fpath = spawn_dir + 'model_' + t_exp + '_day.dat'
        abun_plot_fpath = spawn_dir + 'abundance.png'
            
        if self.structure_type == 'file' and self.abundance_type == 'file':
            if self.event == '05bl' or '11fe':                
                self.velocity_array, self.density_array, self.abun = \
                  self.read_structure(self.event, t_exp)

        else:
            flag_make_structure = False
                                    
        #Call routine to scale the mass fraction of elements.
        abun_up1 = self.scale_abun(self.abun, self.el1_scaling,
                                   scale_F1, v_start_F1, v_stop_F1)
        abun_up2 = self.scale_abun(abun_up1, self.el2_scaling,
                                   scale_F2, v_start_F2, v_stop_F2)

        #Call routine to add the mass fraction of elements.
        abun_up3 = self.add_abun(abun_up2, self.el_adding,
                                 parcel_A, v_start_A, v_stop_A)
            
        #If necessary, write abundance and density files and make abun plot.
        if self.structure_type == 'file' and self.abundance_type == 'file':                
            self.make_density_file(self.dens_fpath, self.velocity_array,
                                   self.density_array)
            self.make_abundance_file(self.abun_fpath, abun_up3)            
            #self.make_structure_file(self.model_fpath, self.velocity_array,
            #                         self.density_array, abun_up3)
            if self.plot_abun:
                self.make_abundance_plot(abun_plot_fpath, self.velocity_array,
                                         abun_up3)    
    
    def write_yml(self):
        """ Main function to write the output .yml files."""
        
        for k, simulation in enumerate(self.simulation_list):

            #Organize each simulation as a dictionary. i.e. unfold the
            #simulation parameters (with prefix) into a dictionary.
            PARS = {}    
            for entry in simulation:
                variable, value = entry.split('|')
                PARS[variable] = value

            #Set name variables.
            filename = self.filename_of_simulation(simulation)
            spawn_dir = os.path.abspath(self.subdir_fullpath + filename) + '/'
            ymlfile_fullpath = spawn_dir + filename + '.yml'         
            self.created_ymlfiles_list.append(ymlfile_fullpath)

            #If necessary, make density and abundance files.            
            if not os.path.exists(spawn_dir):
                os.mkdir(spawn_dir)
            self.control_structure_files(
              spawn_dir, PARS['time_explosion'],
              float(PARS[self.el1_scaling['el'] + '-F1']),
              float(PARS['v_start_F1']), float(PARS['v_stop_F1']),
              float(PARS[self.el2_scaling['el'] + '-F2']),
              float(PARS['v_start_F2']), float(PARS['v_stop_F2']),
              float(PARS[self.el_adding['el'] + '-A']),
              float(PARS['v_start_A']), float(PARS['v_stop_A']))

            if self.verbose:
                print '    CREATED: ' + ymlfile_fullpath

            with open(ymlfile_fullpath  ,'w') as yml_file:
                yml_file.write('tardis_config_version: v1.0\n')
                yml_file.write('supernova:\n')
                yml_file.write('    luminosity_requested: '
                                    +PARS['loglum']+' log_lsun\n')
                if PARS['temperature_requested'] != 'blank': 
                    yml_file.write('    temperature_requested: '
                                   +PARS['temperature_requested']+' K\n')
                yml_file.write('    time_explosion: '
                                    +PARS['time_explosion']+' day\n')
                if self.distance != 'blank':
                    yml_file.write('    distance: ' + self.distance + '\n')  
                yml_file.write('\n')
                yml_file.write('atom_data: ' + self.atomic_data + '\n')
                yml_file.write('\n')
                yml_file.write('model:\n')
                yml_file.write('\n')
                yml_file.write('    structure:\n')
                yml_file.write('        type: ' + self.structure_type + '\n')
                if self.structure_type == 'specific':
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
                                   + self.density_type + '\n')
                    
                    if self.density_type == 'exponential':
                        yml_file.write('            time_0:  '
                                       +PARS['time_0']+' day\n')
                        yml_file.write('            rho_0: '
                                       +PARS['rho_0']+' g/cm^3\n')
                        yml_file.write('            v_0: '
                                       +PARS['v_0']+' km/s\n')
                        yml_file.write('            exponent: '
                                       +PARS['exponent']+'\n')
                    elif self.density_type == 'uniform':
                        yml_file.write('            value: '
                                       +PARS['density_value']+' g/cm^3\n') 
                elif self.structure_type == 'file':   
                    yml_file.write('        filename: ' + self.dens_fpath + '\n')                                    
                    yml_file.write('        filetype: simple_ascii\n')
                    yml_file.write('        v_inner_boundary: '
                                   +PARS['velocity_start']+' km/s\n')
                    yml_file.write('        v_outer_boundary: '
                                   +PARS['velocity_stop']+' km/s\n')                        
                yml_file.write('\n')
                yml_file.write('    abundances:\n')
                yml_file.write('        type: ' + self.abundance_type + '\n')
                if self.abundance_type == 'uniform':
                    for el in self.elements[:-2]:
                        yml_file.write('        ' + el
                                       + ': ' + self.abun[el][0] + '\n')                        
                elif self.abundance_type == 'file':
                    yml_file.write('        filename: ' + self.abun_fpath + '\n') 
                    yml_file.write('        filetype: simple_ascii\n')              
                yml_file.write('\n')
                yml_file.write('plasma:\n')
                yml_file.write('    ionization: ' + self.ionization + '\n')
                yml_file.write('    excitation: ' + self.excitation + '\n')
                yml_file.write('    radiative_rates_type: '
                               + self.rad_rates_type + '\n')
                yml_file.write('    line_interaction_type: '
                               + PARS['line_interaction'] + '\n')
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
                yml_file.write('    start: ' + self.spec_start + ' angstrom\n')
                yml_file.write('    stop: ' + self.spec_stop + ' angstrom\n')
                yml_file.write('    num: ' + self.spec_num)                          
            
        if self.verbose:
            print '\n*** DONE - SUCCESSFUL RUN\n\n'
                
    def run(self):
        self.print_run_start()
        self.collect_non_default_pars()
        self.make_outfolder()
        self.make_combination_of_parameters()
        self.write_yml()
        return self.created_ymlfiles_list
        
