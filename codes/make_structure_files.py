#!/usr/bin/env python

"""
THIS CODE MAKES DENSITY AND ABUNDANCE FILES ACCORDING TO THE PARAMETERS
FROM THE INPUT FILES.
"""

import os                                                               
import sys
import time
import copy
import numpy as np
import itertools                                                        
import random

class make_structure(object):

    def __init__(self,  abundance_dict, velocity_array, pass_density_as,
                 filename, es, ms, t_exp, time_0=None, rho_0=None,
                 v_0=None, density_array_given=None, verbose=True):

        self.abun_dict = abundance_dict
        self.filename_structure = filename
        self.pass_density_as = pass_density_as
        self.density_array_given = density_array_given
        self.velocity_array = velocity_array
        self.es = es
        self.ms = ms
        self.t_exp = t_exp
        self.time_0 = time_0
        self.rho_0 = rho_0
        self.v_0 = v_0
        self.verbose = verbose
        self.elements = ['H', 'He', 'Li', 'B', 'Be', 'C', 'N', 'O', 'F',
                         'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl',
                         'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn',
                         'Fe', 'Co', 'Ni', 'Cu', 'Zn']

        if self.verbose:
            print '----------------------------------------------------'
            print '----------- RUNNING make_structure_files -----------'
            print '----------------------------------------------------'
            print '\n' 
            
            print 'DESCRIPTION:'
            print '    ONE abundance and one velocity/density .dat \
                   FILES WILL BE CREATED ACCORDING TO THE PASSED INPUT \
                   PARAMETERS.\n'
            print 'MAKING FILES:'           

        self.make_density_file()
        self.make_abundance_file()

    def make_density_file(self):                        
        if self.pass_density_as == 'from_exponential':              
            def rho_of_v(v):
                return self.rho_0*np.exp(-v/self.v_0)           
            density_array = rho_of_v(np.asarray(self.velocity_array)
                            .astype(float))

            density_file = './../INPUT_FILES/DENSITY_FILES/'+'density_'
                           +self.filename_structure+'.dat'   
            out_density = open(density_file, 'w')
            
            out_density.write(self.time_0 +'  \n')      
            out_density.write('# index velocity (km/s) density (g/cm^3)'
                              )
            
            for i, (velocity, density) in enumerate(zip(
                                   self.velocity_array, density_array)):
                out_density.write('\n'+str(i)+'    '+str(velocity)
                                  +'    '+str(density))                 
        
            out_density.close()

            if self.verbose:
                print '    DENSITY FILE: '+density_file
    
        elif self.pass_density_as == 'by_hand':
            density_array = np.asarray(self.density_array_given)
                                       .astype(float)    
    
            if isinstance(self.es, str):
                self.es = [self.es]
            if isinstance(self.ms, str):
                self.ms = [self.ms]     
        
            for e_s in self.es:
                for m_s in self.ms:

                    velocity_scaling, density_scaling =
                                     float(e_s)**(-1.5)*float(m_s)**2.5,
                                     float(e_s)**0.5*float(m_s)**(-1.5)
                                     
                    velocity_array = velocity_scaling*np.array(
                                   self.velocity_array).astype(np.float)
                                   
                    density_array = density_scaling*density_array

                    density_file = './../INPUT_FILES/DENSITY_FILES/'
                                   +'density_'+self.filename_structure
                                   +'_es:'+e_s+'_ms:'+m_s+'.dat' 
                                   
                    out_density = open(density_file, 'w')
                    
                    out_density.write(self.time_0 +'  \n')      
                    out_density.write('# index velocity (km/s)  \
                                       density (g/cm^3)')
                    
                    for i, (velocity, density) in enumerate(zip(
                                        velocity_array, density_array)):
                                            
                        out_density.write('\n'+str(i)+'    '
                                     +str(velocity)+'    '+str(density))                 

                    out_density.close()

                    if self.verbose:
                        print '    DENSITY FILE: '+density_file

        else:   
            raise ValueError("density type 'exponential' is allowed \
                              (supplied %s) " % (self.density_type))
        
    def make_abundance_file(self):

        def decay_Ni_Co_Fe(time,Ni_formed):
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
            *np.exp(-time/Co_half_life)
                                  )
            
            Ni_change = Ni_remain - Ni_formed
            Co_change = Co_remain
            Fe_change = Fe_remain
                        
            """Round numbers to write the abundance file"""
            Ni_change = float(format(Ni_change, '.4f'))
            Co_change = float(format(Co_change, '.4f'))
            Fe_change = float(format(Fe_change, '.4f'))           
            
            if abs(Ni_change + Co_change + Fe_change) > 1.e-5:
                #Test if rounded decay products sum to initial Ni_formed
                Ni_change -= Ni_change + Co_change + Fe_change
            
            return Ni_change, Co_change, Fe_change

        if isinstance(self.t_exp, list):
            aux_t_exp = self.t_exp
        else:
            aux_t_exp = [self.t_exp]
                
        for t_exp in aux_t_exp:
            abun = copy.deepcopy(self.abun_dict)
            
            abundance_file = './../INPUT_FILES/STRATIFIED_COMPOSITION\
            _FILES/'+'abundance_'+self.filename_structure+'_'+t_exp
            +'_day.dat'
            
            out_abundance = open(abundance_file, 'w')
            out_abundance.write('# index Z=1 - Z=30')   

            for i, velocity in enumerate(self.velocity_array):
                out_abundance.write('\n'+str(i))
                
                """Compute changes due to decay of Ni and Co."""     
                Ni_initial = float(abun['Ni'][i])
                            
                Ni_change, Co_change, Fe_change = 
                                decay_Ni_Co_Fe(float(t_exp),Ni_initial)
                                
                abun['Ni'][i] = str(format(float(abun['Ni'][i])
                + Ni_change, '.4f'))  
                
                abun['Co'][i] = str(format(float(abun['Co'][i])
                + Co_change, '.4f'))  
                
                abun['Fe'][i] = str(format(float(abun['Fe'][i])
                + Fe_change, '.4f'))  

                """Test if abundances add up to 1 in each layer."""
                sum_elem = 0.    
                for element in self.elements:
                    sum_elem += float(abun[element][i])
                if abs(sum_elem - 1.) > 1.e-5:
                    raise ValueError("Error: In index %s , velocity \
                    = %s , the abundance does not add to 1. Needs %s"
                    % (i, velocity, 1-sum_elem))      
            
                """Write in the output file."""
                for element in self.elements:
                    out_abundance.write(' '+str(abun[element][i]))
        
            out_abundance.close()

            if self.verbose:
                print '    ABUNDANCE FILE: '+abundance_file
                
        if self.verbose:
            print '\n*** DONE - SUCCESSFUL RUN\n\n'

        return None
        
