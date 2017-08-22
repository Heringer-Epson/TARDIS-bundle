#!/usr/bin/env python

import os                                                               
import sys
import shutil
import pickle
import tardis
import numpy as np
import matplotlib.pyplot as plt
import tardis.tardistools.tardis_kromer_plot as tkp
import tardis.tardistools.tardis_minimal_model as tmm
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
                 extinction=0., show_figs=False, verbose=True):

        self.atom_data_dir = os.path.join(os.path.split(tardis.__file__)[0],
          'data', 'kurucz_cd23_chianti_H_He.h5')
        self.created_ymlfiles_list = created_ymlfiles_list
        self.make_kromer = make_kromer
        self.show_figs = show_figs
        self.extinction = extinction
        self.verbose = verbose

        self.run_uncertainties = run_uncertainties
        self.smoothing_window = smoothing_window
        self.N_MC_runs = N_MC_runs 

        self.simulation = None

        if self.verbose:
            print '----------------------------------------------------'
            print ('------------ RUNNING TARDIS v'+str(tardis.__version__)
            +' -----------')
            print '----------------------------------------------------'
            print '\n' 
            
            print 'DESCRIPTION:'
            print '    A TARDIS SIMULATION WILL BE CREATED FOR EACH .yml FILE.\n'

    def show_interface(self):
        if self.show_figs:
            interface.show(self.simulation)

    def make_kromer_plot(self, output, ylim=None):
        minmodel = tmm.minimal_model(mode="virtual")
        minmodel.from_interactive(self.simulation)
        plotter = tkp.tardis_kromer_plotter(minmodel, mode="virtual")
        plotter.generate_plot(xlim=(1500,11000), ylim=(0.,5.), twinx=False)      
        plt.savefig(output, format='png', dpi=360)
        if self.show_figs:
            plt.show()

    def analyse_and_add_quantities(self):
        """ Add some extra information about simulation to the .pkl
        output file. This output file contains not only the synthetic
        spectra but also info such as the temperature at the photosphere.
        If more infor is needed, check dir(self.simulation),
        dir(self.simulation.model), etc...
        """        
        D = {}

        wavelength = self.simulation.runner.spectrum_virtual.wavelength[::-1]
        flux = self.simulation.runner.spectrum_virtual.luminosity_density_lambda[::-1]                       
        
        #Note that the synthetic spectra are not corrected for redshift.
        #Instead, the observed spectra are. Only the wavelength and flux are
        #stored with no units to allow re-computing quantities without
        #re-running the simulation.
        D['wavelength_raw'] = wavelength.value
        D['flux_raw'] = flux.value
        D['host_redshift'] = 0.
        D['extinction'] = self.extinction
        D['t_rad'] = self.simulation.model.t_rad.cgs
        D['luminosity_requested'] = self.simulation.luminosity_requested.cgs
        D['seed'] = self.simulation.runner.seed
        D['t_inner'] = self.simulation.model.t_inner.cgs
        D['v_inner'] = self.simulation.model.v_inner.cgs
        D['v_outer'] = self.simulation.model.v_outer.cgs
        D['w'] = self.simulation.model.w
        D['time_explosion'] = self.simulation.model.time_explosion.cgs
        D['density'] = self.simulation.model.density.cgs
        D['r_outer'] = self.simulation.model.r_outer.cgs
        D['volume'] = self.simulation.model.volume.cgs
                
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
                    D['Integrated_number_density_'+str(el)] = np.nan
                    D['Integrated_ion_density_'+el+'_'+num] = np.nan
                    D['Fraction_'+el+'_'+num] = np.nan
                
        #Delete the simulation to make sure the memory is being freed.
        del self.simulation
                             
        return D, wavelength, flux

    def run_SIM(self):      
        for i, ymlfile_fullpath in enumerate(self.created_ymlfiles_list):
            spawn_dir = os.path.dirname(ymlfile_fullpath) + '/'
            file_prefix = ymlfile_fullpath.split('/')[-1].split('.yml')[0]
            
            outfile = spawn_dir + file_prefix
                
            self.simulation = tardis.run_tardis(
              ymlfile_fullpath, self.atom_data_dir)            
        
            #if self.show_figs:
            #    interface.show(self.simulation)
            
            if self.make_kromer:
                kromer_output = spawn_dir + file_prefix +'_kromer.png'
                self.make_kromer_plot(kromer_output)
        
            D, w, f = self.analyse_and_add_quantities() 
            
            #Create .pkl containg the spectrum and derived qquantities.
            with open(outfile + '.pkl', 'w') as out_pkl:
                pickle.dump(D, out_pkl, protocol=pickle.HIGHEST_PROTOCOL)

            #Create .dat file containg spectra only for uploading into WISEREP.
            with open(outfile + '.dat'  ,'w') as out_spec:
                for x, y in zip(w[:-1], f[:-1]):
                    out_spec.write(str(x.value) + '    ' + str(y.value) + '\n')
                out_spec.write(str(w[-1].value) + '    ' + str(f[-1].value))
    
        if self.verbose:            
            print '\n*** DONE - SUCCESSFUL RUN.'
            print '           TARDIS version used: ' \
                  + str(tardis.__version__) + '\n\n'
            
