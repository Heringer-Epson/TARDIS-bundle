#!/usr/bin/env python

import numpy as np
from astropy															import constants as const

#######################  CODE DESCRIPTION  #############################

"""
TBW
"""

############################  UPDATE LOG  ##############################

"""
VERSION1:
	-Not stable.
VERSION2:
	-Still building the structure of the program. More detailed comments will be given after stable versions.	
"""

###############################  MAIN  #################################

class input_parameters(object):
	
	def __init__(self):		



        self.input_file = __file__.split('/')[-1].split('.py')[0]
		self.subdir														= 'test_kromer/'#'seeds_with_unc_noisier/'
		self.filename_structure											= '2011fe'
		self.mode														= 'sequential'
		self.luminosity													= '9.44'#= ['0.08e9', '0.32e9', '1.1e9', '2.3e9', '3.2e9', '3.5e9', '3.2e9', '2.3e9']
		self.time_explosion												= '13.0'#['3.7', '5.9', '9.0', '12.1', '16.1', '19.1', '22.4', '28.3']	
		self.velocity_start												= '11000'#['13300', '12400', '11300', '10700', '9000', '7850', '6700', '4550']
		self.velocity_stop												= '20000'
		self.luminosity_units											= 'logsolar'
		self.structure_type												= 'specific' #specific #file
		self.abundance_type												= 'uniform' #file,uniform,branch_w7

		self.ionization													= 'lte'#'nebular'
		self.excitation													= 'lte'#'dilute-lte'
		self.rad_rates_type												= 'dilute-blackbody'
		self.line_interaction											= 'macroatom'

		self.seeds														= ['23111963']#list(np.arange(1.e4,1.0001e7,1.e4).astype(int).astype(str))#['23111963']
		self.num_packs													= '4.0e+4'
		self.iterations													= '1'#'5'
		self.last_num_packs												= '5.0e+4'#'1.0e+5'
		self.num_virtual_packs											= '5'#'10'

		self.run_uncertainties											= False
		self.smoothing_window											= 21
		self.N_MC_runs													= 3000

		self.make_kromer												= True
	
		self.convert_luminosity_to_logsolar()

		if self.structure_type == 'specific':
			self.structure_num											= '20'
			self.density_type											= 'branch85_w7' #or_exponential
			
			if self.density_type == 'exponential': #don't use
				self.time_0												= 1.
				self.rho_0												= 9.e-9
				self.v_0												= 2670.
				self.density_value										= 'blank'
			elif self.density_type == 'branch85_w7':
				self.time_0, self.rho_0, self.v_0, self.density_value	= 'blank', 'blank', 'blank', 'blank'			
			elif self.density_type == 'uniform':
				self.density_value										= '1.e-11'
				self.time_0, self.rho_0, self.v_0						= 'blank', 'blank', 'blank'			
				
					
		elif self.structure_type == 'file':
			"""
			To make an stratified abundance file first write the desired velocities in 'self.velocity_array'.
			Then in the matrix below below modify then corresponding columns only.
			The program will ensure that the 'velocity_start' and 'velocity_stop' parameters match the given array.  
			"""

			self.mass_ratio_scaling										= '1.0'
			self.energy_ratio_scaling									= '1.0'

			#----------------------------------------------------------#
			#Prototype - 2011fe rho11 based
			self.velocity_array											= ['4500', '5000', '5500', '6000', '6500', '7000', '7500', '8000', '8500', '9000', '9500', '10000', '10500', '11000', '11500', '12000', '12500', '13000', '13500', '14000', '14500', '15000', '15500', '16000', '17000', '18000', '19000', '19500', '21000', '24000']
			#----------------------------------------------------------#

			
			#----------------------------------------------------------#
			#2011fe rho11 velocity layers
			#self.velocity_array											= ['4500', '5000', '5500', '6000', '6500', '7000', '7500', '8000', '8500', '9000', '9500', '10000', '10500', '11000', '11500', '12000', '12500', '13000', '13500', '14000', '14500', '15000', '15500', '16000', '17000', '18000', '19000', '19500', '21000', '24000']
			#----------------------------------------------------------#

			#----------------------------------------------------------#
			#2005bl w7 velocity layers
			#self.velocity_array											= ['3500', '4000', '4500', '5000', '5500', '6000', '6500', '7000', '7500', '8000', '8500', '9000', '10000', '11000', '12000',  '13000', '14000', '15000', '15500', '16000', '17000']
			#----------------------------------------------------------#

			self.pass_density_as										= 'by_hand'#'from_exponential'# or 'by_hand'
		
			if self.pass_density_as == 'from_exponential':
				self.time_0												= '100 s' #check this.
				self.rho_0												= 1.8
				self.v_0												= 2320.
				self.density_array										= None
				
			elif self.pass_density_as == 'by_hand':

				#-------------------------------------------------------#
				#Prototype - 2011fe rho11 based
				self.density_array										= ['0.65', '0.58', '0.5', '0.45', '0.35', '0.32', '0.3', '0.28', '0.23', '0.21', '0.205', '0.205', '0.2', '0.2', '0.17', '0.12', '0.08', '0.06', '0.04', '0.03', '0.025', '0.02', '0.015', '0.01', '0.006', '0.004', '0.003', '0.002', '0.001', '0.0002']				
				
				#self.density_array										= ['0.65', '0.58', '0.5', '0.45', '0.35', '0.32', '0.3', '0.28', '0.23', '0.21', '0.205', '0.205', '0.2', '0.2', '0.17', '0.12', '0.08', '0.06', '0.04', '0.03', '0.025', '0.02', '0.015', '0.01', '0.006', '0.008', '0.006', '0.004', '0.002', '0.0004']				
				#self.density_array										= ['0.65', '0.58', '0.5', '0.45', '0.35', '0.32', '0.3', '0.28', '0.23', '0.21', '0.205', '0.205', '0.2', '0.2', '0.17', '0.12', '0.08', '0.06', '0.04', '0.03', '0.025', '0.02', '0.015', '0.01', '0.006', '0.020', '0.015', '0.010', '0.005', '0.0010']				
				#-------------------------------------------------------#

				#-------------------------------------------------------#
				#2011fe rho11 density profile#
				#self.density_array										= ['0.65', '0.58', '0.5', '0.45', '0.35', '0.32', '0.3', '0.28', '0.23', '0.21', '0.205', '0.205', '0.2', '0.2', '0.17', '0.12', '0.08', '0.06', '0.04', '0.03', '0.025', '0.02', '0.015', '0.01', '0.006', '0.004', '0.003', '0.002', '0.001', '0.0002']				
				#-------------------------------------------------------#
			
				#-------------------------------------------------------#
				#2005bl w7 density profile 
				#self.density_array										= ['1.278', '1.150','0.859', '0.733', '0.659', '0.493', '0.443', '0.359', '0.298', '0.268', '0.223', '0.19', '0.17', '0.16', '0.10', '0.079', '0.057', '0.062', '0.048', '0.04', '0.0058'] #w7 8k to 23k
				#-------------------------------------------------------#

				self.time_0												= '100 s'
				self.rho_0, self.v_0									= None, None
			
		if self.abundance_type == 'file':		
			#-Do not modify the following four lines-#
			self.abun_raw, self.abun									= {}, {}
			self.elements												= ['H', 'He', 'Li', 'B', 'Be', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']
			self.velocity_zones											= ['3500',   '4000',   '4500',   '5000',   '5500',   '6000',   '6500',   '7000',   '7500',   '8000',   '8500',   '9000',   '9500',  '10000',  '10500', '11000', '11500', '12000', '12500', '13000',  '13500',  '14000',  '14500',  '15000',  '15500',  '16000', '16500',  '17000', '17500',  '18000', '18500', '19000',   '19500',  '20000', '21000',  '24000']
			##########################################


			#Prototype - 2011fe based.
			self.abun_raw['H']											= ['0.000', '0.000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.000', '0.000', '0.000', '0.000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.0000', '0.000', '0.0000', '0.000', '0.0000', '0.0000', '0.000', '0.0000', '0.0000']
			self.abun_raw['He']											= ['0.000', '0.000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.000', '0.000', '0.000', '0.000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.0000', '0.000', '0.0000', '0.000', '0.0000', '0.0000', '0.000', '0.0000', '0.0000']
			self.abun_raw['Li']											= ['0.000', '0.000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.000', '0.000', '0.000', '0.000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.0000', '0.000', '0.0000', '0.000', '0.0000', '0.0000', '0.000', '0.0000', '0.0000']
			self.abun_raw['Be']											= ['0.000', '0.000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.000', '0.000', '0.000', '0.000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.0000', '0.000', '0.0000', '0.000', '0.0000', '0.0000', '0.000', '0.0000', '0.0000']
			self.abun_raw['B']											= ['0.000', '0.000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.000', '0.000', '0.000', '0.000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.0000', '0.000', '0.0000', '0.000', '0.0000', '0.0000', '0.000', '0.0000', '0.0000']
			self.abun_raw['C']											= ['0.000', '0.000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0080', '0.0080', '0.0080', '0.0080', '0.0080', '0.0080', '0.008', '0.008', '0.000', '0.000', '0.000', '0.0310', '0.0310', '0.0310', '0.0310', '0.0310', '0.0260', '0.026', '0.0260', '0.060', '0.0260', '0.250', '0.0260', '0.9804', '0.000', '0.9804', '0.9804']
			self.abun_raw['N']											= ['0.000', '0.000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.000', '0.000', '0.000', '0.000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.0000', '0.000', '0.0000', '0.000', '0.0000', '0.0000', '0.000', '0.0000', '0.0000']
			self.abun_raw['O']											= ['0.020', '0.100', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0200', '0.0900', '0.0900', '0.0900', '0.0900', '0.0900', '0.110', '0.110', '0.351', '0.351', '0.351', '0.7030', '0.7030', '0.7030', '0.7030', '0.7030', '0.8605', '0.916', '0.8605', '0.831', '0.8605', '0.860', '0.8605', '0.0120', '0.000', '0.0120', '0.0120']
			self.abun_raw['F']											= ['0.000', '0.000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.000', '0.000', '0.000', '0.000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.0000', '0.000', '0.0000', '0.000', '0.0000', '0.0000', '0.000', '0.0000', '0.0000']
			self.abun_raw['Ne']											= ['0.000', '0.000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.000', '0.000', '0.000', '0.000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.0000', '0.000', '0.0000', '0.000', '0.0000', '0.0000', '0.000', '0.0000', '0.0000']
			self.abun_raw['Na']											= ['0.001', '0.000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.000', '0.000', '0.000', '0.000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.0000', '0.000', '0.0000', '0.000', '0.0000', '0.0000', '0.000', '0.0000', '0.0000']
			self.abun_raw['Mg']											= ['0.000', '0.020', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.000', '0.020', '0.020', '0.020', '0.0300', '0.0300', '0.0300', '0.0300', '0.0300', '0.0300', '0.010', '0.0300', '0.040', '0.0300', '0.030', '0.0300', '0.0030', '0.000', '0.0030', '0.0030']
			self.abun_raw['Al']											= ['0.000', '0.000', '0.0900', '0.0900', '0.0900', '0.0900', '0.0900', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.000', '0.000', '0.000', '0.000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.0000', '0.000', '0.0000', '0.000', '0.0000', '0.0000', '0.000', '0.0000', '0.0000']
			self.abun_raw['Si']											= ['0.120', '0.000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.2000', '0.2165', '0.2285', '0.1985', '0.4785', '0.4785', '0.4785', '0.4740', '0.563', '0.563', '0.440', '0.440', '0.440', '0.2000', '0.2000', '0.2000', '0.2000', '0.2000', '0.0600', '0.001', '0.0600', '0.060', '0.0600', '0.000', '0.0600', '0.0040', '0.000', '0.0040', '0.0040']
			self.abun_raw['P']											= ['0.000', '0.000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.000', '0.000', '0.000', '0.000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.0000', '0.000', '0.0000', '0.000', '0.0000', '0.0000', '0.000', '0.0000', '0.0000']
			self.abun_raw['S']											= ['0.030', '0.040', '0.0300', '0.0300', '0.0300', '0.0300', '0.0300', '0.0600', '0.0700', '0.0700', '0.0700', '0.1500', '0.1500', '0.1500', '0.1500', '0.150', '0.150', '0.080', '0.080', '0.080', '0.0300', '0.0300', '0.0300', '0.0300', '0.0300', '0.0200', '0.001', '0.0200', '0.007', '0.0200', '0.000', '0.0200', '0.0005', '0.000', '0.0005', '0.0005']
			self.abun_raw['Cl']											= ['0.000', '0.000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.000', '0.000', '0.000', '0.000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.0000', '0.000', '0.0000', '0.000', '0.0000', '0.0000', '0.000', '0.0000', '0.0000']
			self.abun_raw['Ar']											= ['0.000', '0.000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.000', '0.000', '0.000', '0.000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.0000', '0.000', '0.0000', '0.000', '0.0000', '0.0000', '0.000', '0.0000', '0.0000']
			self.abun_raw['K']											= ['0.000', '0.000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.000', '0.000', '0.000', '0.000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.0000', '0.000', '0.0000', '0.000', '0.0000', '0.0000', '0.000', '0.0000', '0.0000']
			self.abun_raw['Ca']											= ['0.010', '0.020', '0.0001', '0.0001', '0.0001', '0.0001', '0.0001', '0.0030', '0.0030', '0.0030', '0.0030', '0.0030', '0.0030', '0.0030', '0.0030', '0.003', '0.003', '0.003', '0.003', '0.003', '0.0030', '0.0030', '0.0030', '0.0030', '0.0030', '0.0020', '0.020', '0.0020', '0.001', '0.0020', '0.000', '0.0020', '0.0000', '0.000', '0.0000', '0.0000']
			self.abun_raw['Sc']											= ['0.000', '0.000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.000', '0.000', '0.000', '0.000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.0000', '0.000', '0.0000', '0.000', '0.0000', '0.0000', '0.000', '0.0000', '0.0000']
			self.abun_raw['Ti']											= ['0.000', '0.000', '0.0010', '0.0010', '0.0010', '0.0010', '0.0010', '0.0050', '0.0050', '0.0050', '0.0050', '0.0050', '0.0050', '0.0050', '0.0050', '0.005', '0.005', '0.003', '0.003', '0.003', '0.0005', '0.0005', '0.0005', '0.0005', '0.0005', '0.0000', '0.000', '0.0000', '0.000', '0.0000', '0.000', '0.0000', '0.0000', '0.000', '0.0000', '0.0000']
			self.abun_raw['V']											= ['0.000', '0.000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.000', '0.000', '0.000', '0.000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.0000', '0.000', '0.0000', '0.000', '0.0000', '0.0000', '0.000', '0.0000', '0.0000']
			self.abun_raw['Cr']											= ['0.000', '0.000', '0.0010', '0.0010', '0.0010', '0.0010', '0.0010', '0.0050', '0.0050', '0.0050', '0.0050', '0.0050', '0.0050', '0.0050', '0.0050', '0.005', '0.005', '0.003', '0.003', '0.003', '0.0005', '0.0005', '0.0005', '0.0005', '0.0005', '0.0000', '0.000', '0.0000', '0.000', '0.0000', '0.000', '0.0000', '0.0000', '0.000', '0.0000', '0.0000']
			self.abun_raw['Mn']											= ['0.000', '0.000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.000', '0.000', '0.000', '0.000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.0000', '0.000', '0.0000', '0.000', '0.0000', '0.0000', '0.000', '0.0000', '0.0000']
			self.abun_raw['Fe']											= ['0.200', '0.120', '0.1500', '0.1500', '0.1500', '0.1500', '0.1500', '0.0005', '0.0005', '0.0005', '0.0005', '0.0005', '0.0005', '0.0005', '0.0050', '0.006', '0.006', '0.060', '0.060', '0.060', '0.0010', '0.0010', '0.0010', '0.0010', '0.0010', '0.0005', '0.001', '0.0005', '0.001', '0.0005', '0.000', '0.0005', '0.0001', '0.000', '0.0001', '0.0001']
			self.abun_raw['Co']											= ['0.000', '0.000', '0.5000', '0.5000', '0.5000', '0.5000', '0.5000', '0.5000', '0.3500', '0.3300', '0.3100', '0.1300', '0.1300', '0.1300', '0.1300', '0.075', '0.075', '0.020', '0.020', '0.020', '0.0005', '0.0005', '0.0005', '0.0005', '0.0005', '0.0005', '0.000', '0.0005', '0.000', '0.0005', '0.000', '0.0005', '0.0000', '0.000', '0.0000', '0.0000']
			self.abun_raw['Ni']											= ['0.619', '0.580', '0.2279', '0.2279', '0.2279', '0.2279', '0.2279', '0.2265', '0.3500', '0.3300', '0.3100', '0.1300', '0.1300', '0.1300', '0.1300', '0.075', '0.075', '0.020', '0.020', '0.020', '0.0005', '0.0005', '0.0005', '0.0005', '0.0005', '0.0005', '0.000', '0.0005', '0.000', '0.0005', '0.000', '0.0005', '0.0000', '0.000', '0.0000', '0.0000']
			self.abun_raw['Cu']											= ['0.000', '0.000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.000', '0.000', '0.000', '0.000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.0000', '0.000', '0.0000', '0.000', '0.0000', '0.0000', '0.000', '0.0000', '0.0000']
			self.abun_raw['Zn']											= ['0.000', '0.000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.000', '0.000', '0.000', '0.000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.000', '0.0000', '0.000', '0.0000', '0.000', '0.0000', '0.0000', '0.000', '0.0000', '0.0000']

			
			'''
			canvas
			self.abun_raw['H']											= ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
			self.abun_raw['He']											= ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
			self.abun_raw['Li']											= ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
			self.abun_raw['Be']											= ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
			self.abun_raw['B']											= ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
			self.abun_raw['C']											= ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
			self.abun_raw['N']											= ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
			self.abun_raw['O']											= ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
			self.abun_raw['F']											= ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
			self.abun_raw['Ne']											= ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
			self.abun_raw['Na']											= ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
			self.abun_raw['Mg']											= ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
			self.abun_raw['Al']											= ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
			self.abun_raw['Si']											= ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
			self.abun_raw['P']											= ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
			self.abun_raw['S']											= ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
			self.abun_raw['Cl']											= ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
			self.abun_raw['Ar']											= ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
			self.abun_raw['K']											= ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
			self.abun_raw['Ca']											= ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
			self.abun_raw['Sc']											= ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
			self.abun_raw['Ti']											= ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
			self.abun_raw['V']											= ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
			self.abun_raw['Cr']											= ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
			self.abun_raw['Mn']											= ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
			self.abun_raw['Fe']											= ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
			self.abun_raw['Co']											= ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
			self.abun_raw['Ni']											= ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
			self.abun_raw['Cu']											= ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
			self.abun_raw['Zn']											= ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
			'''

			self.velocity_zones											= np.asarray(self.velocity_zones)
			for element in self.elements:
				self.abun[element]										= [value for value,velocity in zip(self.abun_raw[element],self.velocity_zones) if velocity in self.velocity_array]

		elif self.abundance_type == 'uniform':
			self.abun_Ni												= '0.0'
			self.abun_Si												= '0.52'
			self.abun_Fe												= '0.0'
			self.abun_Co												= '0.0'
			self.abun_Ca												= '0.03' 
			self.abun_S													= '0.19'
			self.abun_Mg												= '0.03'
			self.abun_Na												= '0.0'
			self.abun_C													= '0.0'
			self.abun_O													= '0.18'
			self.abun_Ti												= '0.01'
			self.abun_Ar												= '0.04'

	def convert_luminosity_to_logsolar(self):
		L_sun															= const.L_sun.cgs.value
		if self.luminosity_units == 'logsolar':
			pass		
		elif self.luminosity_units == 'logcgs':
			if isinstance(self.luminosity, str):
				self.luminosity											= str(format(np.log10(10.**float(self.luminosity)/L_sun), '.2f'))
			elif isinstance(self.luminosity, list):
				self.luminosity											= [str(format(np.log10(10.**float(lum)/L_sun), '.2f')) for lum in self.luminosity]
		elif self.luminosity_units == 'cgs':
			if isinstance(self.luminosity, str):
				self.luminosity											= str(format(np.log10(float(self.luminosity)/L_sun), '.2f'))
			elif isinstance(self.luminosity, list):
				self.luminosity											= [str(format(np.log10(float(lum)/L_sun), '.2f')) for lum in self.luminosity]
		elif self.luminosity_units == 'solar':
			if isinstance(self.luminosity, str):
				self.luminosity											= str(format(np.log10(float(self.luminosity)), '.2f'))
			elif isinstance(self.luminosity, list):
				self.luminosity											= [str(format(np.log10(float(lum)), '.2f')) for lum in self.luminosity]				
		return None
