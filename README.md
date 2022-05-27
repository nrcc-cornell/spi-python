# Standardized Precipitation Index (SPI) using Python

## Calculating SPI using the python script: calcSPI.py

#### Docstring of calcSPI.py
```shell
'''Calculate SPI.

	- available functions use incomplete gamma to fit data and calculate probabilities
	- main section provides example of SPI calculation using these functions.
	- 20220207: updated to work with both python2 and python3

'''
```

#### Sample SPI calculations using calcSPI.py functions (available in main section of calcSPI.py)
```shell
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
```

## Results comparison with original FORTRAN code (compare-with-fortran-src)

We successfully compared results of this python script with FORTRAN code obtained from the Colorado Climate Center.

1. edit calcSPI.f, specifying which precipitation data file to use

2. compile FORTRAN source code to create executable
```shell
$ gfortran -o calcSPI.e calcSPI.f spiRoutines.f
```

3. run executable to calculate and write out SPI calculations.
```shell
$ ./calcSPI.e
```

Results of SPI calculation tests are in 'fortran' and 'python' *.out files.

