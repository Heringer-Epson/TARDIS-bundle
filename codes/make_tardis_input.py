#!/usr/bin/env python

import os                                                               
import sys
import shutil
import copy
import itertools                                                        
import numpy as np
import astropy.units as u

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
    def __init__(self, inputs, verbose):
          
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

        #Up to two elements whose mass fraction will be be scaled in
        #a region of the ejecta.
        self.el1_scaling = inputs.el1_scaling
        self.el2_scaling = inputs.el2_scaling
        self.MASTER[self.el1_scaling['el']] = inputs.el1_scaling['factors']
        self.MASTER[self.el2_scaling['el']] = inputs.el2_scaling['factors']
                
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
            read_every = 2
            t_exp2str = {'11.0': 'm6', '12.0': 'm5', '14.0': 'm3', '21.8': 'p48',
                         '29.9': 'p129'}                
            phase = t_exp2str[t_exp]
            fpath = ('./../INPUT_FILES/Hachinger_2005bl/models-05bl-w7e0.7/' 
                            + 'SN2005bl_' + phase + '/' + 'abuplot.dat')
                
        elif event == '11fe':
            read_every = 1 
            fpath = ('./../INPUT_FILES/Mazzali_2011fe/ejecta_layers.dat')                   
        
        rows = []
        with open(fpath, 'r') as inp:
            for line in itertools.islice(inp, 2, None, read_every):
                column = line.rstrip('\n').split(' ') 
                column = filter(None, column)
                rows.append(column)               

        #Extract velocity, density and abundances.
        velocity_array = np.asarray([float(row[1]) for row in rows])
        density_array = np.asarray([10.**float(row[3]) for row in rows])
                
        #Get abundances.
        abun = {} 
        el_loop = self.elements + ['Fe0', 'Ni0']
        
        for el in el_loop:
            abun[el] = np.zeros(len(velocity_array))
        
        for i, velocity in enumerate(velocity_array):
            sum_elem = 0.               
            for j, el in enumerate(el_loop):
                
                abun_value = float(rows[i][4 + j])
                abun[el][i] = abun_value
                sum_elem += abun_value

            if abs(sum_elem - 1.) > 1.e-5:
                raise ValueError("Error: In index %s , velocity \
                = %s , the abundance does not add to 1. Needs %s"
                % (i, velocity, 1. - sum_elem))        

        return velocity_array, density_array, abun

    def scale_abun(self, inp_abun, el_scaling, scale):
        
        element = el_scaling['el']
        v_start = el_scaling['v_start']
        #Make independent copy of the abun dictionary to prevent differentt_exp
        #from modifying an already modified (decayed, or scaled) abundundaces. 
        abun = copy.deepcopy(inp_abun)
    
        if element != 'None':
            shells = np.arange(0, len(self.velocity_array), 1)
            condition = (self.velocity_array >= v_start)
            shell_window = shells[condition]
            
            for i in shell_window:
                abundance_all = np.asarray([abun[el][i] for el in self.elements])
                abundance_all = np.nan_to_num(abundance_all)
                el_most = self.elements[abundance_all.argmax()]
                orig_abun = float(abun[element][i])
                new_abun = scale * orig_abun
                                
                #Update abundances.
                abun[element][i] = str(new_abun)            
                abun[el_most][i] = str(float(abun[el_most][i])
                                       + (orig_abun - new_abun))        
        
        return abun
        
    def compute_Ni_decay(self, inp_abun, time):
        
        #Make independent copy of the abun dictionary to prevent differentt_exp
        #from modifying an already modified (decayed, or scaled) abundundaces. 
        abun = copy.deepcopy(inp_abun)
                
        def decay_Ni_Co_Fe(time, Ni_formed):
            """Calculates the decay chain Ni->Co->Fe.
            
            Source: Nadyozhin+ 1994
            http://adsabs.harvard.edu/full/1994ApJS...92..527N
            
            Remarks: Difference between halflife and lifetime: the 
            former is the intuitive (1/2)^(t/Tau), while lifetime uses
            the exponential form. The correlation between constants
            is T_exp = T_1/2 *(1/ln(2))
            
            Values: Half-lifes are taken from Jonghwa Chang's
            (Korea Atomic Energy Research Institute)
            http://atom.kaeri.re.kr/nuchart/
            """
           
            Ni_half_life = 6.07/np.log(2.)
            Co_half_life = 77.24/np.log(2.)
           
            Ni_remain = Ni_formed*np.exp(-time/Ni_half_life) 
             
            Co_remain = Ni_formed*(Co_half_life/(Co_half_life
            - Ni_half_life))*(np.exp(-time/Co_half_life)
            - np.exp(-time/Ni_half_life))
            
            Fe_remain = Ni_formed*(1. + Ni_half_life/(Co_half_life
            - Ni_half_life)*np.exp(-time/Ni_half_life)
            - Co_half_life/(Co_half_life - Ni_half_life)
            *np.exp(-time/Co_half_life))
            
            Ni_change = Ni_remain - Ni_formed
            Co_change = Co_remain
            Fe_change = Fe_remain
                        
            """Round numbers to write the abundance file"""
            Ni_change = float(format(Ni_change, '.6f'))
            Co_change = float(format(Co_change, '.6f'))
            Fe_change = float(format(Fe_change, '.6f'))           
            
            if abs(Ni_change + Co_change + Fe_change) > 1.e-5:
                #Test if rounded decay products sum to initial Ni_formed
                Ni_change -= Ni_change + Co_change + Fe_change
            
            return Ni_change, Co_change, Fe_change

        #for i, velocity in enumerate(self.velocity_array):
        for i in range(len(abun['Si'])):
            
            """Compute changes due to decay of Ni and Co."""     
            Ni_initial = float(abun['Ni0'][i])
            Ni_change, Co_change, Fe_change = decay_Ni_Co_Fe(
                                            float(time),Ni_initial)

            abun['Ni'][i] = str(format(float(abun['Ni0'][i])
            + Ni_change, '.6f'))  
            
            abun['Co'][i] = str(format(float(abun['Co'][i])
            + Co_change, '.6f'))  
            
            abun['Fe'][i] = str(format(float(abun['Fe0'][i])
            + Fe_change, '.6f'))  

            """Test if abundances add up to 1 in each layer."""
            sum_elem = 0.    
            for element in self.elements:
                sum_elem += float(abun[element][i])
            if abs(sum_elem - 1.) > 1.e-5:
                raise ValueError("Error: In index %s,\
                the abundance does not add to 1. Needs %s"
                % (i, 1-sum_elem))      

        return abun

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
                                    
    def control_structure_files(self, spawn_dir, t_exp, scale1, scale2):
        """This function uses the name of the input_pars_X.py file to determine
        whether to make a density and an abundance files. Note that
        input_pars_Hach.py has its own routine to read the data provided by
        Hachinger and create structure files readable by TARDIS.
        """

        self.dens_fpath = spawn_dir + 'density_' + t_exp + '_day.dat'
        self.abun_fpath = spawn_dir + 'abundance_' + t_exp + '_day.dat'
            
        if self.structure_type == 'file' and self.abundance_type == 'file':
            if self.event == '05bl' or '11fe':                
                self.velocity_array, self.density_array, self.abun = \
                  self.read_structure(self.event, t_exp)

        else:
            flag_make_structure = False
                                    
        #Call routine to scale the mass fraction of elements.
        abun_up1 = self.scale_abun(self.abun, self.el1_scaling, scale1)
        abun_up2 = self.scale_abun(abun_up1, self.el2_scaling, scale2)
        
        #Call routine to compute the decay of 56Ni.
        abun_decayed = self.compute_Ni_decay(abun_up2, t_exp)
            
        #If necessary, write abundance and density files.
        if self.structure_type == 'file' and self.abundance_type == 'file':                
            self.make_density_file(self.dens_fpath, self.velocity_array,
                                   self.density_array)
            self.make_abundance_file(self.abun_fpath, abun_decayed)
    
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
            self.control_structure_files(spawn_dir, PARS['time_explosion'],
                                         float(PARS[self.el1_scaling['el']]),
                                         float(PARS[self.el2_scaling['el']]))
        
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
                        yml_file.write('            time_0: '
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
                    for el in self.elements:
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
        
