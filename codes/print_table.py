#!/usr/bin/env python

import os   

"""Uncomment input file to be imported"""
from input_pars_11fe import Input_Parameters as class_input

class Make_Table(class_input):
   
    """ Prints a table in .tex format (to be copied and pasted) containing the
    element stratification of 11fe. 05bl is not included as the data from
    Hachinger is in a different format.
    """
    
    def __init__(self):

        class_input.__init__(self)
        self.SN = self.input_file.split('.py')[0][-4::]
        zones = self.velocity_zones[::-1]           
        self.elements = ['C', 'O', 'Na', 'Mg', 'Al', 'Si', 'S', 'Ca', 'Ti',
                         'Cr', 'Fe', 'Ni']
        
        if self.SN == '11fe':
            n_layers = '11'
        elif self.SN == '05bl':
            n_layers = '6'
                                                       
        for i, v in enumerate(zones):
            
            if i == 0:
                print ('\multirow{' + n_layers + '}{*}{' + self.SN + '} & '
                        + v + ' $-$ ' + self.velocity_stop),
            else:
                print (' & ' + v + ' $-$ ' + zones[i - 1]),
                
            for el in self.elements:
                print (' & ' + self.abun_raw[el][::-1][i]),
            print ' \\\\'               

Make_Table()      
