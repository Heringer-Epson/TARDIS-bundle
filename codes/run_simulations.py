#!/usr/bin/env python

import os                                                               
import sys
import shutil

version = os.path.basename(__file__).split('_v')[-1].split('.py')[0]

import tardis
import matplotlib.pyplot as plt

from tardis.gui import interface
from tardis.tardistools.compute_features import analyse_spectra as analyse
from tardis.tardistools.compute_features import uncertainty as uncertainty_routine
import tardis.tardistools.tardis_kromer_plot as tkp
import tardis.tardistools.tardis_minimal_model as tmm

import pandas as pd
import cPickle

import numpy as np


class Simulate_Spectra(object):

    """Runs one TARDIS model.
    
    Parameters
    ----------
    
    subdir : ~str
        Subdirectory within 'INPUT_FILES/YML_FILES/' that contains
        the input parameter files to be simulated by TARDIS.
        
    created_ymlfiles_list : ~list
        Each entry in the list corresponds to the a created .yml file.
        
    run_uncertainties : ~boolean
        If True, computes the uncertainty of some spectra features.
                     
    smoothing_window : ~float
        Window to be used by the Savitzky-Golay filter to smooth the spectra.
        
    N_MC_runs : ~int
        Number of MC runs to be used to compute uncertainties.
        Only applicable if run_uncertainties is True.
        
    make_kromer : ~boolean
        Flag to whether or not make a 'Kromer style' plot.
        Requires the tardis_kromer_plot and tardis_minimal_model modules
        available under the tardistools, which is to be installed concomitantly
        to TARDIS.
        
    display_interface : ~boolean
        Flag to whether or not display the default TARDIS graphical interface.
        
    verbose : ~boolean
        Flag to whether or not print extra information.                                
                
    Output
    -------
    1) A .pkl file containing the synthetic spectra and other info about
    the simulation (such as temperature at the photosphere.
    2) A 'Kromer' plot (if make_kromer == True).
    
    Files are acreated at
    './../INPUT_FILES/YML_FILES/'+subdir
    Where 'subdir' defined in the input file.

    Returns
    -------
    self.output_dir: Directory where the synthetic spectra is saved.
                     Useful for copying the input file to that directory
                     in the master code.                            
    """

    def __init__(self, subdir, created_ymlfiles_list,
                 run_uncertainties=True, smoothing_window=21,
                 N_MC_runs=3000, make_kromer=False,
                 display_interface=False, verbose=True):

        self.subdir = subdir
        self.input_dir = './../INPUT_FILES/YML_FILES/'+subdir
        self.output_dir = './../OUTPUT_FILES/'+subdir
        self.atom_data_dir = os.path.join(os.path.split(tardis.__file__)[0],
          'data', 'kurucz_cd23_chianti_H_He.h5')
        self.created_ymlfiles_list = created_ymlfiles_list
        self.make_kromer = make_kromer
        self.display_interface = display_interface
        self.verbose = verbose

        self.run_uncertainties = run_uncertainties
        self.smoothing_window = smoothing_window
        self.N_MC_runs = N_MC_runs 

        self.simulation = None
        self.kromer_figures = []
        self.hdf_files = []

        if self.verbose:
            print '----------------------------------------------------'
            print ('------------ RUNNING TARDIS v'+str(tardis.__version__)
            +' -----------')
            print '----------------------------------------------------'
            print '\n' 
            
            print 'DESCRIPTION:'
            print '    A TARDIS SIMULATION WILL BE CREATED FOR EACH .yml FILE.\n'
        
    def check_subdir(self): 
        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)
        return None 

    def show_interface(self):
        if self.display_interface:
            interface.show(self.simulation)
        return None

    def make_kromer_plot(self, output, ylim=None):
        minmodel = tmm.minimal_model(mode="virtual")
        minmodel.from_interactive(self.simulation)
        plotter = tkp.tardis_kromer_plotter(minmodel, mode="virtual")
        plotter.generate_plot(xlim=(2500,11000), twinx=False)      
        plt.tight_layout()
        plt.savefig(output, format='png', dpi=360)
        #plt.show()
        self.kromer_figures.append(output)
        return None

    def analyse_and_add_quantities(self):
        """ Add some extra information about simulation to the .pkl
        output file. This output file contains not only the synthetic
        spectra but also info such as the temperature at the photosphere.
        If more infor is needed, check dir(self.simulation),
        dir(self.simulation.model), etc...
        """        
        wavelength = self.simulation.runner.spectrum_virtual.wavelength[::-1]
        flux = (self.simulation.runner.spectrum_virtual
          .luminosity_density_lambda[::-1])
        D = {}
        D = pd.DataFrame({'wavelength_raw': [wavelength],
                         'flux_raw': [flux], 'host_redshift': [0.]})
                    
        D['t_rad'] = [self.simulation.model.t_rad.cgs]
        D['luminosity_requested'] = [self.simulation.luminosity_requested.cgs]
        D['seed'] = [self.simulation.runner.seed]
        D['t_inner'] = [self.simulation.model.t_inner.cgs.value]
        D['v_inner'] = [self.simulation.model.v_inner.cgs]
        D['v_outer'] = [self.simulation.model.v_outer.cgs]
        D['w'] = [self.simulation.model.w]
        D['time_explosion'] = [self.simulation.model.time_explosion.cgs]
        D['density'] = [self.simulation.model.density.cgs]
        D['r_outer'] = [self.simulation.model.r_outer.cgs]
        D['volume'] = [self.simulation.model.volume.cgs]

        #Add the number fraction of the ions for some of the most
        #important elements in the ejecta. The fraction is the integrated
        #fraction across all the shells. 
        for element, el in zip([6, 8, 14, 20, 22, 24], ['C', 'O', 'Si', 'Ca', 'Ti', 'Cr']): 
            for ion, num in zip([0, 1, 2], ['I', 'II', 'III']):
                """Runs through neutral, singly and doubly ionized states"""
                total_number_density = 0.
                total_ion_density = 0.
                try:
                    for i in range(len(self.simulation.model.density.cgs.value)):
                        """Runs through the shells"""
                        total_number_density += (self.simulation.plasma.
                          number_density[i].ix[element])            
                        total_ion_density += (self.simulation.plasma.
                          ion_number_density[i].ix[element].tolist()[ion])
                        
                    total_ion_fraction = total_ion_density / total_number_density
                    D['Integrated_number_density_'+str(el)] = total_ion_density
                    D['Integrated_ion_density_'+el+'_'+num] = total_number_density
                    D['Fraction_'+el+'_'+num] = total_ion_fraction
                except:
                    D['Integrated_number_density_'+str(el)] = 'Failed'
                    D['Integrated_ion_density_'+el+'_'+num] = 'Failed'
                    D['Fraction_'+el+'_'+num] = 'Failed'
                
        #Delete the simulation to make sure the memory is being freed.
        del self.simulation
        
        #Retrieve the pEW, depth, and velocity of some feature and if
        #self.run_uncertainties == True, also compute the corresponding
        #uncertainties via a MC method.
        D = analyse(D, smoothing_mode='savgol',
                    smoothing_window=self.smoothing_window, 
                    verbose=True).run_analysis()       
        if self.run_uncertainties:
            D = uncertainty_routine(
              D, smoothing_mode = 'savgol',
              smoothing_window = self.smoothing_window,
              N_MC_runs = self.N_MC_runs, verbose=True).run_uncertainties()                     
        return D

    def run_SIM(self):      
        self.check_subdir()
        for i, inpfile in enumerate(self.created_ymlfiles_list):
            outfile = self.output_dir+inpfile[:-4]+'.pkl'
                
            self.simulation = tardis.run_tardis(
              self.input_dir+inpfile, self.atom_data_dir)            
        
            if self.display_interface:
                interface.show(self.simulation)
            
            if self.make_kromer:
                kromer_output = self.output_dir+inpfile[:-4]+'.png'
                self.make_kromer_plot(kromer_output)
        
            D = self.analyse_and_add_quantities()     
            D.to_pickle(outfile)
            with open(outfile, 'r') as inp:
                pkl = cPickle.load(inp)     
    
        if self.verbose:            
            print '\n*** DONE - SUCCESSFUL RUN.'
            print '           TARDIS version used: ' \
                  +str(tardis.__version__)+'\n\n'
            
        return self.output_dir

