#!/usr/bin/env python

import os                                                               
import sys
import time
import shutil
import cPickle
import tardis
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import tardis.tardistools.tardis_kromer_plot as tkp
import tardis.tardistools.tardis_minimal_model as tmm
import multiprocessing as mp
from tardis.gui import interface

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
        
    show_figs : ~boolean
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

    def __init__(self, created_ymlfiles_list,
                 run_uncertainties=True, smoothing_window=21,
                 N_MC_runs=3000, make_kromer=False,
                 extinction=0., show_figs=False):

        self.created_ymlfiles_list = created_ymlfiles_list
        self.make_kromer = make_kromer
        self.show_figs = show_figs
        self.extinction = extinction

        self.run_uncertainties = run_uncertainties
        self.smoothing_window = smoothing_window
        self.N_MC_runs = N_MC_runs 

        print '----------------------------------------------------'
        print ('------------ RUNNING TARDIS v'+str(tardis.__version__)
        +' -----------')
        print '----------------------------------------------------'
        print '\n' 
        
        print 'DESCRIPTION:'
        print '    A TARDIS SIMULATION WILL BE CREATED FOR EACH .yml FILE.\n'

    def make_kromer_plot(self, simulation, output, ylim=None):
        minmodel = tmm.minimal_model(mode="virtual")
        minmodel.from_interactive(simulation)
        plotter = tkp.tardis_kromer_plotter(minmodel, mode="virtual")
        plotter.generate_plot(xlim=(1500,11000), ylim=(0.,5.), twinx=False)      
        plt.savefig(output, format='png', dpi=360)
        if self.show_figs:
            plt.show()

    def analyse_and_add_quantities(self, sim):
        """ Add some extra information about simulation to the .pkl
        output file. This output file contains not only the synthetic
        spectra but also info such as the temperature at the photosphere.
        If more infor is needed, check dir(self.simulation),
        dir(self.simulation.model), etc...
        """        
        D = {}

        wavelength = sim.runner.spectrum_virtual.wavelength[::-1]
        flux = sim.runner.spectrum_virtual.luminosity_density_lambda[::-1]
        
        #Note that the synthetic spectra are not corrected for redshift.
        #Instead, the observed spectra are. Only the wavelength and flux are
        #stored with no units to allow re-computing quantities without
        #re-running the simulation.
        D['wavelength_raw'] = wavelength.value
        D['flux_raw'] = flux.value
        D['host_redshift'] = 0.
        D['extinction'] = self.extinction
        D['t_rad'] = sim.model.t_rad.cgs
        D['luminosity_requested'] = sim.luminosity_requested.cgs
        D['seed'] = sim.runner.seed
        D['t_inner'] = sim.model.t_inner.cgs
        D['v_inner'] = sim.model.v_inner.cgs
        D['v_outer'] = sim.model.v_outer.cgs
        D['w'] = sim.model.w
        D['time_explosion'] = sim.model.time_explosion.cgs
        D['density'] = sim.model.density.cgs
        D['r_inner'] = sim.model.r_inner.cgs
        D['r_outer'] = sim.model.r_outer.cgs
        D['volume'] = sim.model.volume.cgs
                
        #Add the number fraction of the ions for some of the most
        #important elements in the ejecta. The fraction is the integrated
        #fraction across all the shells. 
        for element, el in zip([6, 8, 14, 20, 22, 24], ['C', 'O', 'Si', 'Ca', 'Ti', 'Cr']): 
            for ion, num in zip([0, 1, 2], ['I', 'II', 'III']):
                """Runs through neutral, singly and doubly ionized states"""
                total_number_density = 0.
                total_ion_density = 0.
                try:
                    for i in range(len(sim.model.density.cgs.value)):
                        """Runs through the shells"""
                        total_number_density += (sim.plasma.
                          number_density[i].ix[element])            
                        total_ion_density += (sim.plasma.
                          ion_number_density[i].ix[element].tolist()[ion])
                    
                    total_ion_fraction = total_ion_density / total_number_density
                    D['Integrated_number_density_'+str(el)] = total_number_density
                    D['Integrated_ion_density_'+el+'_'+num] = total_ion_density
                    D['Fraction_'+el+'_'+num] = total_ion_fraction
                except:
                    D['Integrated_number_density_'+str(el)] = np.nan
                    D['Integrated_ion_density_'+el+'_'+num] = np.nan
                    D['Fraction_'+el+'_'+num] = np.nan

        #for i in range(len(sim.model.density.cgs.value)):
            
            #Fix this.
            #C_II_10 = sim.plasma.level_number_density[i].ix[6].ix[1].ix[10] 
            
            #print ('T = ', sim.model.t_rad.cgs.value[i], 'K'\
            #       'C_II at level 10 # density = ', (C_II_10 / D['Integrated_number_density_C']))
             
                             
        ###Tests with plasma
        #print dir(sim)
        #print dir(sim.plasma)
        #print sim.plasma.ion_number_density
        #print sim.plasma.number_density
        #print dir(sim.plasma.atomic_data)
        #print sim.plasma.atomic_data.levels
        #print sim.plasma.atomic_data.ionization_data
        #print sim.plasma.atomic_data.lines
        
        #print sim.plasma.level_number_density[shell].ix[El#].ix[ion#].ix[level#]
        #print sim.plasma.level_number_density[0].ix[6].ix[1].ix[10]
        #print sim.plasma.level_number_density[0].ix[6].ix[1]
                
        return D, wavelength, flux

    def run_SIM(self):
        
        #Function to actually run the simulation.
        def perform_run(yml):
            spawn_dir = os.path.dirname(yml) + '/'
            file_prefix = yml.split('/')[-1].split('.yml')[0]            
            outfile = spawn_dir + file_prefix
            
            simulation = tardis.run_tardis(yml)  

            if self.show_figs:
                interface.show(simulation)

            if self.make_kromer:
                kromer_output = spawn_dir + file_prefix +'_kromer.png'
                self.make_kromer_plot(simulation, kromer_output)

            D, w, f = self.analyse_and_add_quantities(simulation) 
            
            #Save simulation as hdf.
            simulation.to_hdf(outfile + '.hdf')
            
            #Create .pkl containg the spectrum and derived qquantities.
            with open(outfile + '.pkl', 'w') as out_pkl:
                cPickle.dump(D, out_pkl, protocol=cPickle.HIGHEST_PROTOCOL)            
            
            #Create .dat file containg spectra only for uploading into WISEREP.
            with open(outfile + '.dat'  ,'w') as out_spec:
                for x, y in zip(w[:-1], f[:-1]):
                    out_spec.write(str(x.value) + '    ' + str(y.value) + '\n')
                out_spec.write(str(w[-1].value) + '    ' + str(f[-1].value))

            #Delete the simulation to make sure the memory is being freed.
            del simulation
            
        time_start = time.time()
        for i, yml in enumerate(self.created_ymlfiles_list):            
            print '\n\nRunning yml #', str(i + 1) + '/'\
                  + str(len(self.created_ymlfiles_list))
            
            if i != 0:
                completed_sims = float(i)
                avg_sim_time = elapsed_time / completed_sims / 3600.
                #estimated_time = (len(yml) - completed_sims) * avg_sim_time                    
                #estimated_time = str(format(estimated_time, '.1f'))
                #print '..Estimated remaining time: ' + estimated_time + 'h\n\n'    
            
            perform_run(yml)
            elapsed_time = time.time() - time_start 
        
        print '\n*** DONE - SUCCESSFUL RUN.'
        print '           TARDIS version used: ' \
              + str(tardis.__version__) + '\n\n'
               
