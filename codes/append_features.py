#!/usr/bin/env python

import numpy as np
import pandas as pd
import pickle

import tardis.tardistools.compute_features as cp

class Analyse_Features(object):

    """This code appends spectral feature values to a pkl file.
    The pkl file might be newly generated by 'run_simulations' code, or pre-
    existing from a previous run. If the file is pre-existing, the features
    are re-computed.
    
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
        
    verbose : ~boolean
        Flag to whether or not print extra information.                                
                
    Output
    -------
    1) A .pkl file where spectral features have been added in addition to the
    spectrum.
    
    Files are acreated at
    './../INPUT_FILES/YML_FILES/'+subdir
    Where 'subdir' defined in the input file.                           
    """

    def __init__(self, created_ymlfiles_list,
                 extinction=0., run_uncertainties=True, smoothing_window=21,
                 N_MC_runs=3000, show_fig=False, verbose=True):

        self.created_ymlfiles_list = created_ymlfiles_list
        self.show_fig = show_fig
        self.verbose = verbose

        self.extinction = extinction
        self.run_uncertainties = run_uncertainties
        self.smoothing_window = smoothing_window
        self.N_MC_runs = N_MC_runs 

        if self.verbose:
            print '----------------------------------------------------'
            print '---------------- COMPUTING FEATURES ----------------'
            print '----------------------------------------------------'
            print '\n' 
            
            print 'DESCRIPTION:'
            print '    Appends pEW, depth and velocity to the simulations.\n'
        
        self.get_features()
    
    def get_features(self):
        """Reads each .yml file generated by the master_run code and append
        a pEW, velocity and depth value to some spectral features. The
        spectral features are defined in the imported 'Analyse_Spectra' class.
        """ 
        for inpfile in self.created_ymlfiles_list:
            file_path = inpfile.split('.yml')[0]+'.pkl'
            with open(file_path, 'r+') as inp:
                D = pickle.load(inp)

                #Check if features had been computed in a previous run.
                #If so, then wipe feature values to prevent the error where
                #the column already exists when adding a new column.
                if 'pEW_f7' in D.keys():
                    for key in D.keys():
                       
                        keys_minimal = [
                          'wavelength_raw', 'flux_raw','host_redshift',
                          'extinction', 't_rad', 'luminosity_requested',
                          'seed', 't_inner', 'v_inner', 'v_outer', 'w',
                          'time_explosion', 'density', 'r_outer', 'volume']

                        if key not in keys_minimal:
                            del D[key]
                
                #Perform feature analysis.              
                #pkl = cp.Analyse_Spectra(
                #  pkl, smoothing_mode='savgol',
                #  smoothing_window=self.smoothing_window,
                #  verbose=True).run_analysis()

                w = D['wavelength_raw']
                f = D['flux_raw']
                
                D = cp.Analyse_Spectra(
                  wavelength=D['wavelength_raw'].value,
                  flux=D['flux_raw'].value,
                  redshift=D['host_redshift'], extinction=D['extinction'],
                  D=D, smoothing_window=self.smoothing_window,
                  verbose=True).run_analysis()
                                            
                #Perfomer calclulation of uncertainties.
                if self.run_uncertainties:
                    D = cp.Compute_Uncertainty(
                      D=D, smoothing_window=self.smoothing_window,
                      N_MC_runs=self.N_MC_runs, verbose=True).run_uncertainties()       

                if self.show_fig:
                    cp.Plot_Spectra(D, show_fig=self.show_fig,
                                    save_fig=False)
                            
                #Create .pkl containg the spectrum and derived qquantities.
                with open(file_path + '.pkl', 'w') as handle:
                    pickle.dump(D, handle, protocol=pickle.HIGHEST_PROTOCOL)                
                                                           

