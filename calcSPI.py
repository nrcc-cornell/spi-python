'''Calculate SPI.

	- available functions use incomplete gamma to fit data and calculate probabilities
	- main section provides example of SPI calculation using these functions.
	- 20220207: updated to work with both python2 and python3

	bnb2
'''
import sys
import os
import numpy as np
import scipy.stats as ss
import scipy.special as ssp

def gammainc_parameters(l):
	'''Calculate estimates of parameters for incomplete gamma.

		Incomplete gamma used to deal with zeros in the data we are fitting.
		Method of using the zero percentage was acquired from FORTRAN code
		obtained from Colorado Climate Center. The incomplete gamma was
		the method used by McKee when calculating SPI.

		We do a few things:
		1. calculate the percentage of the data that are zero.
		2. for data that are not zero maximum likelihood estimation is used
			to determine the shape and scale parameters.
			(Thom, 1958; Wilks, pg 89)
		3. return the estimated parameters and percentage of zeros.

		l = 2-D numpy array of values to fit
		D = sample statistic used in maximum likelihood estimation

		returned parameters:
		pzero : probability of zero (1-D numpy array)
		shape: shape parameter (1-D numpy array)
		scale: scale parameter (1-D numpy array)

		bnb2
	'''
	### create masked array
	### here the zeros are masked and ignored in subsequent calcs
	l_ma = np.ma.masked_less_equal(l,0.00)

	### calculate probability of zero from the mask.
	### in the mask, values of True correspond to zeros.
	l_mask = l<=0.00
	pzero = l_mask.sum(axis=1)/float(l_mask.shape[1])

	### calculate sample statistic (D)
	### NOTE: np.log continues to calculate the log of each element in an array,
	### even when some values are masked. This really shouldn't happen, the masked
	### values should be completely ignored. If the data array contains zeros, this
	### would produce a runtime warning, but does not affect results as the mask
	### continues to be honored and applied in np.log output. However, to suppress
	### this warning, a conditional is added to np.log. This affects performance
	### on large datasets, but does not affect results. If speed is an issue, remove
	### the 'where' conditional and tolerate the runtime warning.
	term1 = np.log(np.mean(l_ma, axis=1))
	term2 = np.sum( np.log( l_ma, where=np.logical_not(l_mask) ), axis=1 ) / ((1.-pzero)*l_ma.shape[1])
	D = term1 - term2

	### calculate shape and scale using Thom's method
	num = 1 + np.sqrt( 1+(4.*D/3.) )
	den = 4.*D
	shape = num/den
	scale = np.mean(l_ma, axis=1)/shape

	### if all of the values in a row are zero, or
	### if all but one of the values are zero,
	### then McKee "assigns something reasonable" for shape/scale.
	### Since we are using masked arrays, not necessary - the
	### masked arrays are returned.

	return shape,scale,pzero

def calcGammaincProbs(d,vals):
	'''Calculate probabilities from incomplete gamma distribution fitted to list of data.
		From a 2-D numpy array of given values (d):
		- fit incomplete gamma distribution to each row of values
		- calculate probabilities for fitted distributions using cdf
		
		variables:
		d = 2-D numpy array of data values (masked arrays)
		vals = 2-D numpy array of values to be evaluated by cdf
		result = 2-D list of probabilities

		bnb2
	'''
	### calculate incomplete gamma parameters
	shape,scale,pzero = gammainc_parameters(d)

	### calculate probs of interest for distributions that have these parameters.
	### use the incomplete gamma function built-in to scipy
	if isinstance(vals,list): vals = np.array(vals)
	result = [
		np.where( x<=0.0, v3, v3+(1.0-v3)*ssp.gammainc(v1,x/v2) )
		for v1,v2,v3,x in zip(shape,scale,pzero,vals)
		]

	return result

if __name__ == "__main__":

	# sample SPI calculations using these functions

	# some sample precipitation data to perform test on
	sample_data = [
		1.20, 0.49, 1.13, 1.13, 1.22, 0.93, 0.90, 0.50, 3.62, 0.51, 0.22, 0.03, 0.13, 0.50,
		2.13, 0.41, 1.24, 1.04, 1.15, 1.37, 0.82, 0.87, 1.51, 2.21, 0.40, 0.51, 0.005, 0.17,
		0.75, 0.60, 0.95, 0.61, 0.31, 1.56, 0.19, 0.58, 0.47, 1.32, 0.22, 1.77, 1.03, 0.67,
		0.24, 0.87, 2.71, 0.61, 1.77, 0.53, 0.37, 0.71, 2.17, 0.34, 1.62, 0.61, 0.48, 1.69,
		2.43, 0.62, 3.29, 1.15, 1.61, 2.61, 0.59, 0.02, 1.97, 0.74, 0.09, 0.53, 1.02, 0.76,
		0.07, 3.07, 2.26, 1.42, 1.93, 0.50, 0.64, 0.54, 0.42, 0.00, 0.62, 0.73, 0.96, 0.55,
		0.34, 2.81, 0.71, 0.03, 1.65, 1.21, 1.38, 0.51, 0.30, 1.71, 0.34, 0.66, 2.09, 1.24,
		0.31, 0.07]

	# 1) Data list converted to a 2-D numpy array. This is expected by calcGammaincProbs, and
	# allows for calculation of multiple locations (stations or grid points), simultaneously
	# when multiple rows are present. Each row of the 2-D array is evaluated independently.
	# In this simple example, there is only one row of data to evaluate in the 2-D array.
	sample_data = np.array(sample_data).reshape(1,len(sample_data))

	# 2) Define the data values to calculate SPI for. In this test, we'll find SPI for all values.
	data_to_evaluate = sample_data.tolist()

	# 3) Fit distribution to sample_data, then calculate probabilities from CDF for data_to_evaluate.
	probs = calcGammaincProbs(sample_data, data_to_evaluate)

	# 4) Calculate SPI from probabilities. Use inverse CDF of standard normal distribution to do this.
	spi = ss.norm.ppf(probs)

	# print to check values
	for values in zip(data_to_evaluate[0],probs[0],spi[0]): print('%8.3f %8.3f %8.3f' % tuple(values))

