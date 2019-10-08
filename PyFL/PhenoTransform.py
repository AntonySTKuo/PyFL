from numba import jit
import numpy as np
from numba import vectorize, float64


@jit(nopython = True) # Set "nopython" mode for best performance, equivalent to @njit	
def hill_func(input_nparray, Vmax = 1.0, Michaelis_const = 0.5, n = 2):
	output = Vmax * input_nparray**n / (Michaelis_const**n + input_nparray**n)
	return output

@jit(nopython = True)
def boltzman(input_nparray, xmid = 0.0, slope = 1.0, maxi = 1.0, mini = 0.0):
	output = mini + (maxi - mini)/(1.0 + np.exp(-(input_nparray-xmid)/slope))
	return output


@jit(nopython = True)
def gaussian_selection(input_nparray, mu = 1.0, sigma = 0.05):
	output = 1.0/(sigma * np.sqrt(2.0 * np.pi)) * np.exp(-np.power(input_nparray - mu, 2.0)/(2.0 * np.power(sigma, 2.0)))
	return output


@jit(nopython = True)
def arctan(input_nparray):
	output = np.arctan(input_nparray/np.max(input_nparray)*100.0) / np.pi * 2.0
	return output

@vectorize([float64(float64, float64)])
def sum(x, y):
    return x + y
