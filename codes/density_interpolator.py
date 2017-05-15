#!/usr/bin/env python

import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

def W7_velocity_to_density(velocity_list):
	"""W7 profile read from Mazzali+ 2014
	[[http://adsabs.harvard.edu/abs/2014MNRAS.439.1959M]].
	The y-axis is interpolated in logspace and returned in linear space.
	"""	
	input_velocity = map(float, velocity_list)
	
	"""Interpolates in logspace."""
	W7_velocity_data = np.array([
	  1000., 5000., 9000., 10000., 13360., 14120., 15170., 15644., 16000.,
	  16680., 18020., 19790., 21560., 22720., 23400., 23900., 24600.
	  ])
	
	W7_density_data = np.array([
	  np.log10(3.068), np.log10(0.7367), np.log10(0.1962),
	  np.log10(0.1596), np.log10(0.0580), np.log10(0.06307),
	  np.log10(0.0482), np.log10(0.05578), np.log10(0.03842),
	  np.log10(0.00762), np.log10(0.00333), np.log10(0.00194),
	  np.log10(0.000983), np.log10(0.0004294), np.log10(0.0001163),
	  np.log10(0.00001955), np.log10(0.00000226)
	  ])
	
	W7_function = interp1d(W7_velocity_data, W7_density_data)	
	return list(10.**np.asarray(W7_function(np.asarray(input_velocity))))

def rho11fe_velocity_to_density(velocity_list):
	"""rho11fe profile read from Mazzali+ 2014
	[[http://adsabs.harvard.edu/abs/2014MNRAS.439.1959M]].
	The y-axis is interpolated in logspace and returned in linear space.
	"""
	input_velocity = map(float, velocity_list)
	
	"""Interpolates in logspace."""
	rho11_velocity_data = np.array([
	  280., 780., 1910., 3880., 4640., 5330., 5870., 6815., 8010., 9210.,
	  11130., 11750., 15135., 16735., 18042., 22440., 24220, 29080., 37110.])
	
	rho11_density_data = np.array([
	  np.log10(2.65), np.log10(4.094), np.log10(2.385), np.log10(0.671),
	  np.log10(0.6303), np.log10(0.512), np.log10(0.4806),
	  np.log10(0.3743), np.log10(0.2915), np.log10(0.2088),
	  np.log10(0.1802), np.log10(0.1347), np.log10(0.0165),
	  np.log10(0.00648), np.log10(0.00410), np.log10(0.000462),
	  np.log10(0.000156), np.log10(0.0000413), np.log10(0.00000206)
	  ])
	
	rho11_function = interp1d(rho11_velocity_data, rho11_density_data)	
	
	return list(10.**np.asarray(rho11_function(np.asarray(input_velocity))))
