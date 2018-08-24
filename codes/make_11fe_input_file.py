#!/usr/bin/env python

import os                                                               
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.interpolate import interp1d
from scipy.integrate import cumtrapz
from astropy import constants as const
from astropy import units as u
from matplotlib.ticker import MultipleLocator

M_sun = const.M_sun.to('g').value

class Make_11fe_file(object):
    """
    Code Description
    ----------
    
    Abundance files have datapoints all located at the same mass coordinate,
    read by digitalizing the published plot and only moving the 'cursor'
    vertically for each element.
    
    The original code used by Mazzali does not interpolate abundances between
    layers and such the slopes plotted are artificial - in reality the plot
    should be step wise.
    
    Comparing the plotted curves and the info in the paper, one can conclude
    that the slopes represent the data to the left. For instance X(C) = 98%
    for v>=19400km/s, therefore a point placed on the steep slope preceding 98%
    would actually correspond to the value to the left (which is ~13% as there
    is an actual data point in the slope region.)
    
    Velocities are converted from the mass coordinate read from the
    digitalized abundance plot using the velocity-density plot. Precise
    velocity number are then adjusted to match the description in the paper
    where necessary. This is important as some photospheric velocities are
    placed in a layer transition and one does not want to start at the layer
    below the correct one.
    
    """
    
    def __init__(self):
        
        #Quantities computed in a number of requested shells. This number can
        #be higher (finer) structure than the number of layers where the
        #composition changes. A fined grid helps to increase the 'precision'
        #of the simulation.
        self.velocity_requested = (
          #np.logspace(np.log10(3550.), np.log10(35000.), 140)) 
          np.logspace(np.log10(3550.), np.log10(24000.), 120)) 
        self.density_requested = None
        self.mass_requested = None
        self.abun_requested = {}
        self.N_requested_shells = len(self.velocity_requested)

        #The coarse shells demarcate where the composition of the ejecta changes.
        #These shells were chosen by hand using the changes of slope in the Ni0
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
        self.mass2vel = None
     
        self.N_coarse_shells = None
        
        self.top_dir = './../INPUT_FILES/Mazzali_2011fe/'
        
        self.el = ['C', 'O', 'Mg', 'Si', 'S', 'Ca', 'TiCr', 'Fe0', 'Ni0']
        
        self.el_all = [
          'H', 'He', 'Li', 'B', 'Be', 'C', 'N', 'O', 'F', 'Ne', 'Na',
          'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc',
          'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Fe0', 'Ni0'] 
        
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
                
        time = 100. * u.s
        dens = self.density_fine * u.g / u.cm**3

        #For each shell volumne, use the average density between the shell edges.
        avg_dens = (dens.value[0:-1] + dens.value[1:]) / 2.
        avg_dens = np.array(avg_dens) * u.g / u.cm**3 
        
        v = self.velocity_fine * u.km / u.s
        r = v.to(u.cm / u.s) * time    
        vol = 4. / 3. * np.pi * r**3.
        vol_step = np.diff(vol)
        mass_step = np.multiply(vol_step, avg_dens) / const.M_sun.to('g')
        
        #Add zero value to array so that the length is preserved.
        #mass_step = np.array([0] + list(mass_step))
        
        mass_cord = np.cumsum(mass_step)
        
        #In the original publication, the position of the first point in the
        #density-velocity relationship is not obvious, possibly creating an
        #offset. Empirically, subtracting 0.005 M_sun seems to provide very
        #good agreement between the velocity and mass scales in the abundance
        #plots. Also note that it is okay for the mass_fine array to have
        #len = len(velocity_fine) - 1.
        self.mass_fine = mass_cord.value - 0.01

    def get_mass2vel_conversor(self):
        avg_velocity = (self.velocity_fine[0:-1] + self.velocity_fine[1:]) / 2.
        self.mass2vel = interp1d(self.mass_fine, avg_velocity)

    def get_velocity_array(self):  
        #Use any of the abundance files to retrieve the mass_coarse array,
        #as all abundance files use the exact same mass coordinates.
        fpath = self.top_dir + '/abundance_Ni0.csv'

        self.mass_coarse = np.loadtxt(fpath, dtype=float, delimiter=',', 
                                      usecols=(0,), unpack=True)         
        
        self.velocity_coarse = self.mass2vel(self.mass_coarse)
        self.N_coarse_shells = len(self.mass_coarse)
        
        #Make small corrections necessary to match the paper text. Question
        #marks are modifications based on the photospheric speed used.
        
        #velocity_coarse[-3] is ~16500 but should be 16000km/s.
        self.velocity_coarse[-3] = 15990.

        #velocity_coarse[-6] is ~13600 but should be 13300km/s.
        self.velocity_coarse[-6] = 13290.        

        #velocity_coarse[-6] is ~11400 but should be 11300km/s?
        self.velocity_coarse[-9] = 11290.        

        #velocity_coarse[8] is ~9300 but should be 9000km/s.
        self.velocity_coarse[8] = 8990.        

        #velocity_coarse[4] is ~7900 but should be 7850km/s?
        self.velocity_coarse[4] = 7840.        

        #velocity_coarse[2] is ~6800 but should be 6700km/s?
        self.velocity_coarse[2] = 6690. 

        #velocity_coarse[0] is ~4600 but should be 4550km/s?
        self.velocity_coarse[0] = 4540. 
                        
    def load_abundances(self):
        for el in self.el:
            fpath = self.top_dir + '/abundance_' + el + '.csv'
            abun = np.loadtxt(fpath, dtype=float, delimiter=',', usecols=(1,),
                              unpack=True)   
            
            self.abun_coarse[el] = abun / 100.
            
            #There are no zeros in the log plot. So when digitalizing the data,
            #The 'zero values' were attributed a small number, but they should
            #actually be zero.
            self.abun_coarse[el][self.abun_coarse[el] < 1.e-4] = 0.
        
    def enforce_abun_normalization(self):
        for i in range(self.N_coarse_shells):            
            total = 0.           
            for el in self.el:
                total += self.abun_coarse[el][i]
            for el in self.el:
                self.abun_coarse[el][i] *= 1. / total                            

    def format_abun_dict(self):
        
        zero_array = np.zeros(self.N_coarse_shells)
        nan_array = np.zeros(self.N_coarse_shells)
        nan_array[:] = np.nan

        other_el = ['H', 'He', 'Li', 'B', 'Be', 'N', 'F', 'Ne', 'Na', 'Al',
                    'P', 'Cl', 'Ar', 'K', 'Sc', 'V', 'Mn', 'Co', 'Cu', 'Zn'] 
        
        for el in other_el:
            self.format_abun_format[el] = zero_array

        for el in ['C', 'O', 'Mg', 'Si', 'S', 'Ca', 'Fe0', 'Ni0']:
            self.format_abun_format[el] = self.abun_coarse[el]
            
        self.format_abun_format['Ti'] = self.abun_coarse['TiCr'] / 2.    
        self.format_abun_format['Cr'] = self.abun_coarse['TiCr'] / 2. 
        self.format_abun_format['Fe'] = nan_array
        self.format_abun_format['Ni'] = nan_array

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
        closest zone below and attribute the mass fraction of that zone - this
        is how the Mazzali code operates, despite the confusing way it is
        plotted (not as steps). Use this procedure for all elements.
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
        line1 = 'shno velocity cum_M logdens'
        line2 = 'shno velocity cum_M logdens'
        for i, el in enumerate(self.el_all):
            el_num = i + 1
            if el_num == 31:
                el_num = 98
            if el_num == 32:
                el_num = 99
                            
            line1 += ' ' + str(el_num) 
            line2 += ' ' + el 

        line1 += ' \n'
        line2 += ' \n'

        #Write outfile.             
        with open(fpath, 'w') as out:
            out.write(line1) 
            out.write(line2) 
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
        ax.set_xlim(2500., 24000.)
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
            x = self.velocity_requested
            y = self.abun_requested[el]
            if el == 'Ti':
                y = np.asarray(y) * 2.
            
            ax.step(x, y, color=color[i], label=label[i], lw=3., where='post')

        #plot and print velocity at the coarse shell transitions (i.e. where
        #abundances change.)
        for i, v in enumerate(self.velocity_coarse):
            plt.axvline(x=v, color='k', alpha=0.5, lw=1., ls=':')
            #print 'Shell ' + str(i) + ' starts at v = ' + str(v) + ' [km/s]'

        #Print requested velocities layers used in the simulations.
        print '\n\n'
        for i, v in  enumerate(self.velocity_requested):
            #print 'layer ' + str(i) + ' starts at v = ' + str(v) + ' [km/s]'
            pass
        
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
        list_el_plot = ['C', 'O', 'Mg', 'Si', 'S', 'Ca', 'Ti', 'Fe0', 'Ni0']
        label = ['C', 'O', 'Mg', 'Si', 'S', 'Ca', 'Ti+Cr', 'Fe0', 'Ni0']
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
        self.get_mass2vel_conversor()
        self.get_velocity_array()
        self.load_abundances()
        self.enforce_abun_normalization()
        self.format_abun_dict()
        self.get_density_requested()
        self.get_mass_requested()
        self.get_abun_requested()
        self.make_output()
        self.plot_abundances_xvel()
        self.plot_abundances_xmass()
        
if __name__ == '__main__':
    Make_11fe_file()

