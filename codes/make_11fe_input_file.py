#!/usr/bin/env python

import os                                                               
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.interpolate import interp1d
from scipy.integrate import cumtrapz
from astropy import constants as const
from matplotlib.ticker import MultipleLocator

M_sun = const.M_sun.to('g').value

class Make_11fe_file(object):
    """
    Code Description
    ----------

    TBW
    """
    
    def __init__(self):
        
        #Quantities computed in a number of requested shells. This number can
        #be higher (finer) structure than the number of layers where the
        #composition changes. A fined grid helps to increase the 'precision'
        #of the simulation.
        self.velocity_requested = (
          np.logspace(np.log10(3550.), np.log10(24000.), 100)) 
        self.density_requested = None
        self.mass_requested = None
        self.abun_requested = {}
        self.N_requested_shells = len(self.velocity_requested)

        #The coarse shells demarcate where the composition of the ejecta changes.
        #These shells were chosen by hand using the changes of slope in the Ni56
        #curve of the original plot.
        self.abun_coarse = {}
        self.format_abun_format = {}
        self.density_coarse = []
        self.velocity_coarse = []
        self.mass_coarse = []
               
        #Finer shells are the shells read from the digitalized plots. These
        #shells need not be the same for each element and are used to compute
        #the abundance **at** the position of the coarse shells.
        self.density_fine = []
        self.velocity_fine = []
        self.mass_fine = None
     
        self.N_coarse_shells = None
        
        self.top_dir = './../INPUT_FILES/Mazzali_2011fe/'
        
        self.el = ['C', 'O', 'Mg', 'Si', 'S', 'Ca', 'TiCr', 'Fe58', 'Ni56']
        
        self.el_all = [
          'H', 'He', 'Li', 'B', 'Be', 'C', 'N', 'O', 'F', 'Ne', 'Na',
          'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V',
          'Cr', 'Mn', 'Fe56', 'Fe58', 'Co', 'Ni56', 'Ni58', 'Cu', 'Zn'] 
        
        self.run_make()

    def initialize_dict(self):
        for el in self.el:
            self.abun_coarse[el] = []
        
    def load_density(self):
        fpath = self.top_dir + 'density.csv'
        
        with open(fpath, 'r') as inp:
            for line in inp:
                column = line.rstrip('\n').split(',') 
                column = filter(None, column)
                self.velocity_fine.append(float(column[0]))                        
                self.density_fine.append(float(column[1]))                        

        self.velocity_fine = np.asarray(self.velocity_fine)
        self.density_fine = np.asarray(self.density_fine)
        
    def compute_mass(self):
        time = 100.
        km2cm = 1.e5
        
        r_array = self.velocity_fine * time * km2cm
        int_array = 4. * np.pi * r_array**2. * self.density_fine
        self.mass_fine = cumtrapz(int_array, r_array) / M_sun
        
    def interpolate_abun(self, Mass, Abun):
        
        #If the last value is smaller than 0.02% (0.0002), then it should 
        #actually be zero, see original plot in Mazzali's paper.
        if Abun[-1] < 0.02:
            Abun[-1] = 0.
            
        #Check whether first value should be zero. This is important for C and Mg.
        if Mass[0] > 0.4:
            Abun[0] = 0.
            
        #For the interpolation to always work, extrapolate the first and last
        #values to masses 0 and 1.5m by conserving the extreme abundances.
        Mass = np.asarray([0.] + Mass + [1.5])
        Abun = np.asarray([Abun[0]] + Abun + [Abun[-1]])
        
        mass2abun = interp1d(Mass, Abun)
        return mass2abun 

    def load_abundances(self):
        """ For the Ni curve, the points were there seemed to be a change in
        slope of the abundance track were recorded by hand. These points will
        provide the velocity (and consequently mass at were one needs to
        retrieve the abundances of the other elements.)
        """
        
        #Read the mass at specified points (where slopes change).
        #Note that the Ni56 abundance file was digitalized differently, so that
        #the slope changing points were determined by hand.
        fpath = self.top_dir + '/abundance_Ni56.csv'
        with open(fpath, 'r') as inp:
            for line in inp:
                column = line.rstrip('\n').split(',') 
                column = filter(None, column)
                self.mass_coarse.append(float(column[0]))            
        self.mass_coarse = np.asarray(self.mass_coarse)
               
        #Get abundances at those specified points by interpolation.
        for el in self.el:
            fpath = self.top_dir + '/abundance_' + el + '.csv'
            mass = []
            abun = []

            with open(fpath, 'r') as inp:
                for line in inp:
                    column = line.rstrip('\n').split(',') 
                    column = filter(None, column)
                    mass.append(float(column[0]))
                    abun.append(float(column[1]))                        

            mass2abun = self.interpolate_abun(mass, abun)
            self.abun_coarse[el] = mass2abun(self.mass_coarse) / 100.
        
        self.N_coarse_shells = len(self.mass_coarse)    

    def enforce_abun_normalization(self):
        
        total_per_shell = []
                
        for i in range(self.N_coarse_shells):
            
            total = 0.
            for el in self.el:
                total += self.abun_coarse[el][i]
            
            #In one of the coarse shells the value of the oxygen abundance
            #happend to fall in the middle of a steep slope. So the total value
            #is corrected by adding oxygen.
            if total < 0.8:
                self.abun_coarse['O'][i] += 1. - total
                
            else:
                for el in self.el:
                    self.abun_coarse[el][i] *= 1. / total                            

    def format_abun_dict(self):
        
        zero_array = np.zeros(self.N_coarse_shells)

        self.el = ['C', 'O', 'Mg', 'Si', 'S', 'Ca', 'TiCr', 'Fe58', 'Ni56']
        
        self.el_all = [
          'H', 'He', 'Li', 'B', 'Be', 'C', 'N', 'O', 'F', 'Ne', 'Na',
          'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V',
          'Cr', 'Mn', 'Fe56', 'Fe58', 'Co', 'Ni56', 'Ni58', 'Cu', 'Zn'] 

        other_el = ['H', 'He', 'Li', 'B', 'Be', 'N', 'F', 'Ne', 'Na', 'Al',
                    'P', 'Cl', 'Ar', 'K', 'Sc', 'V', 'Mn', 'Co', 'Fe56',
                    'Ni58', 'Cu', 'Zn'] 
        
        for el in other_el:
            self.format_abun_format[el] = zero_array

        for el in ['C', 'O', 'Mg', 'Si', 'S', 'Ca', 'Fe58', 'Ni56']:
            self.format_abun_format[el] = self.abun_coarse[el]
            
        self.format_abun_format['Ti'] = self.abun_coarse['TiCr'] / 2.    
        self.format_abun_format['Cr'] = self.abun_coarse['TiCr'] / 2. 
        self.format_abun_format['Fe56'] = zero_array
        self.format_abun_format['Ni58'] = zero_array

    def get_velocity_coarse(self):
        avg_velocity = (self.velocity_fine[0:-1] + self.velocity_fine[1:]) / 2.
        mass2velocity = interp1d(self.mass_fine, avg_velocity)
        self.velocity_coarse = mass2velocity(self.mass_coarse)
        
    def get_density_requested(self):
        vel2dens = interp1d(self.velocity_fine, np.log10(self.density_fine))
        self.density_requested = 10.**vel2dens(self.velocity_requested)

    def get_mass_requested(self):
        avg_velocity = (self.velocity_fine[0:-1] + self.velocity_fine[1:]) / 2.
        vel2mass = interp1d(avg_velocity, self.mass_fine)
        self.mass_requested = vel2mass(self.velocity_requested)

    def get_abun_requested(self):
        """Format abundances. Note that the velocity requested can be finer than
        the coarse zones. To each velocity in the requested array, find the
        closest zone below
        and attribute the mass fraction of that zone. Use this procedure for
        all elements.
        """
        self.abun_requested = {}
        for element in self.el_all:
            self.abun_requested[element] = []
            for v in self.velocity_requested:
                condition = (self.velocity_coarse <= v)
                if not True in condition:
                    idx_zone = 0
                else:    
                    idx_zone = self.velocity_coarse[condition].argmax()

                self.abun_requested[element].append(
                  self.format_abun_format[element][idx_zone])
                        
    def make_output(self):
        
        fpath = self.top_dir + '/ejecta_layers.dat'
        
        #Make header
        header = 'shno velocity cum_M logdens'
        for i, el in enumerate(self.el_all):
            header += ' ' + el 
        header += ' \n'

        #Write outfile.             
        with open(fpath, 'w') as out:
            out.write('#SN2011fe\n') 
            out.write(header) 
            for i in range(self.N_requested_shells):            
                shno = str(i)
                vel = str(format(self.velocity_requested[i], '.0f'))
                cum_m = str(format(self.mass_requested[i], '.3f'))
                logdens = str(format(np.log10(self.density_requested[i]), '.5f'))
                out.write(' ' + shno + ' ' + vel + ' ' + cum_m + ' ' + logdens)
                for el in self.el_all:
                    out.write(' ' + str(format(self.abun_requested[el][i], '.6f')))
                out.write(' \n')    

    def plot_abundances_xvel(self):

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
        list_el_plot = ['C', 'O', 'Mg', 'Si', 'S', 'Ca', 'Ti', 'Fe58', 'Ni56']
        label = ['C', 'O', 'Mg', 'Si', 'S', 'Ca', 'Ti+Cr', 'Fe58', 'Ni56']
        color= ['y', 'r', 'b', 'lightgreen', 'peru', 'grey', 'g', 'purple', 'k']
        for i, el in enumerate(list_el_plot):
            x = self.velocity_requested
            y = self.abun_requested[el]
            if el == 'Ti':
                y = np.asarray(y) * 2.
            
            ax.step(x, y, color=color[i], label=label[i], lw=3., where='post')

        #plot and print velocity at the coarse shell transitions (i.e. where
        #abundances change.)
        for i, v in enumerate(self.velocity_coarse):
            plt.axvline(x=v, color='k', alpha=0.5, lw=1., ls=':')
            print 'Shell ' + str(i) + ' starts at v = ' + str(v) + ' [km/s]'

        #Print requested velocities layers used in the simulations.
        print '\n\n'
        for i, v in  enumerate(self.velocity_requested):
            print 'layer ' + str(i) + ' starts at v = ' + str(v) + ' [km/s]'
        
        #Add legend
        ax.legend(frameon=False, fontsize=20., numpoints=1, ncol=1,
                  labelspacing=0.05, loc='best')          

        plt.tight_layout()        

        #Save and show figure.
        directory = './../INPUT_FILES/Mazzali_2011fe/'
        plt.savefig(directory + 'Fig_11fe_abundance_xvel.png', format='png')
        plt.show()    

    def plot_abundances_xmass(self):

        #Arrange figure frame.
        fs = 26.
        fig, ax = plt.subplots(figsize=(14,8))
        x_label = r'enclosed mass $\rm{[M_\odot]}$'
        y_label = r'$\rm{mass\ \ fraction}$'
        ax.set_xlabel(x_label, fontsize=fs)
        ax.set_ylabel(y_label, fontsize=fs)
        ax.set_yscale('log')
        ax.set_xlim(0., 1.4)
        ax.set_ylim(1.e-4, 1.1)
        ax.tick_params(axis='y', which='major', labelsize=fs, pad=8)       
        ax.tick_params(axis='x', which='major', labelsize=fs, pad=8)
        #ax.minorticks_off()
        ax.tick_params('both', length=8, width=1, which='major')
        ax.tick_params('both', length=4, width=1, which='minor')
        ax.xaxis.set_minor_locator(MultipleLocator(0.1))
        ax.xaxis.set_major_locator(MultipleLocator(0.2))  

        ax_top = ax.twiny()
        x_label = r'$v\ \ \rm{[1000\ \ km\ \ s^{-1}]}$'
        ax_top.set_xlabel(x_label, fontsize=fs)
        ax_top.set_xlim(0., 1.4)
        ticks_top = np.arange(2500, 17501, 2500)
        avg_velocity = (self.velocity_fine[0:-1] + self.velocity_fine[1:]) / 2.
        vel2mass = interp1d(avg_velocity, self.mass_fine)
        tick_pos = vel2mass(ticks_top)
        ax_top.set_xticks(tick_pos)
        ax_top.set_xticklabels(ticks_top / 1000.)
        ax_top.tick_params(axis='x', which='major', labelsize=fs, pad=8)       
        ax_top.tick_params('x', length=8, width=1, which='major')
        

        
        #Plot abundances.
        list_el_plot = ['C', 'O', 'Mg', 'Si', 'S', 'Ca', 'Ti', 'Fe58', 'Ni56']
        label = ['C', 'O', 'Mg', 'Si', 'S', 'Ca', 'Ti+Cr', 'Fe58', 'Ni56']
        color= ['y', 'r', 'b', 'lightgreen', 'peru', 'grey', 'g', 'purple', 'k']
        for i, el in enumerate(list_el_plot):
            x = self.mass_requested
            y = self.abun_requested[el]
            if el == 'Ti':
                y = np.asarray(y) * 2.
            
            ax.step(x, y, color=color[i], label=label[i], lw=3., where='post')

        #plot mass at the coarse shell transitions (i.e. where
        #abundances change.)
        for i, m in enumerate(self.mass_coarse):
            plt.axvline(x=m, color='k', alpha=0.5, lw=1., ls=':')

        #Add legend
        ax.legend(frameon=False, fontsize=20., numpoints=1, ncol=1,
                  labelspacing=0.05, loc='best')          

        plt.tight_layout()        

        #Save and show figure.
        directory = './../INPUT_FILES/Mazzali_2011fe/'
        plt.savefig(directory + 'Fig_11fe_abundance_xmass.png', format='png')
        plt.show()    
    
    def run_make(self):
        self.initialize_dict()
        self.load_density()
        self.compute_mass()
        self.load_abundances()
        self.enforce_abun_normalization()
        self.format_abun_dict()
        self.get_velocity_coarse()
        self.get_density_requested()
        self.get_mass_requested()
        self.get_abun_requested()
        self.make_output()
        self.plot_abundances_xvel()
        self.plot_abundances_xmass()
        
if __name__ == '__main__':
    Make_11fe_file()

