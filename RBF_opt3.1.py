# -*- coding: utf-8 -*-
"""
Program to estimate 7 orbital parameters ofa binary system containing an 
RR Lyrae Star


@author HildboldCF

Last modified by Caleb Hildbold - 15:35 5/4/21

"""

import numpy as np
from numpy import sin, cos, pi, sqrt
import matplotlib.pylab as plt
import rbfopt
from scipy.optimize import least_squares
import pickle
import pandas as pd

from numba import njit

#Importing Data 
#csv must be pre-processed and have the form
# HJD     Uncert.     O-C     E
# and contain no quotes
#the easiest way to process this is through the 'Data' tab in excel
folder_location = (r'C:\Users\hildboldcf17\OneDrive - Grove City College\Senior' 
                   ' Year\The Last Semester\RR Lyrae Research\GEOS_data')
file_name = r'\XX_And-2021-03-09T16 09 41.csv'
loc = folder_location + file_name
df = pd.read_csv(loc, error_bad_lines=False,header=0)   # Creates a dataframe to work with
xi = df['E'].to_numpy()             # Extract epochs to use as indep. variable
yi0 = df['HJD'].to_numpy()          # Extract HJD to get times of maxima for dep. variable
yi = yi0-2422500                    # Scaling to make times of maximum less ridiculous


@njit
def make_data(t, params):
    '''
    This function can be used to generate fake data based on given parameters

    '''
    M = 2*pi*sin(2*pi*xi/int(params[3]*365/params[1]))
    E = M + params[4]*sin(M + params[4]*sin(M))
    return (t*params[1] + params[0] + params[2]/2*t**2 + 
            params[6]*(sqrt(1 - params[4]**2)*sin(E)*cos(params[5]) + 
            sin(params[5]*cos(E))))

@njit
def oc1_func(t, params):
    '''
    Represents the linear portion of the O-C diagram
    
    '''
    return params[0] + t*params[1]

@njit 
def oc2_func(t, params):
    ''' 
    Represents the quadratic portion of the O-C diagram
    
    '''
    return params[0]*((t+params[1])**2)/2 + params[2]

@njit
def oc_func(t, pul_p, params):
    '''
    Represents the period portion of the O-C diagram that is caused
    by the orbit in a binary system

    Parameters: [T, e, w, A]
    
    '''
    M = 2*pi*sin(2*pi*xi/int(params[0]*365/pul_p))
    E = M + params[1]*sin(M + params[1]*sin(M))
    return (params[3]*(sqrt(1 - params[1]**2)*sin(E)*cos(params[2]) + 
            sin(params[2]*cos(E))))



starname = 'XX And'

'''
Parameters
Initial Epoch for graph
init_Epoch = 1

Pulsation period
pul_p = 0.59208 days

Parabolic decay term
B = -0.025E-6 days/year

Orbital period
T = 75.3 years

Eccentricity of the orbit 
e = .77

Longitude of periastron
w = 161 degrees
w = w*pi/180 radians

A = 0.064 days
'''


@njit
def optimizing_function(par2, pul_p, y):
    '''
    Returns the sum of squared errors between the data and the periodic variation
    Used in : rbfopt
    Parameters
    ----------
    par2 : periodic term parameters
    
    pulp : pulsation period
    
    y : times of maxima data points (only periodic portion)

    Returns
    -------
    sum of squared errors
    
    '''
    
    resid = (y - oc_func(xi, pul_p, par2))**2
    return np.sum(resid)

@njit 
def opt_lsq_fxn1(params):
    '''
    Least squares uses an array of residuals to optimize parameters
    Used in : calculating y1 (subtracting linear portion)

    Parameters
    ----------
    params : list
            linear portion parameters

    Returns
    -------
    resid : array of residuals between initial data and linear portion of fit

    '''
    resid = yi - oc1_func(xi, params)
    return resid

@njit 
def opt_lsq_fxn2(params, y1, x):
    '''
    Least squares uses an array of residuals to optimize parameters
    Used in : calculating y2 (subtracting quadratic portion)

    Parameters
    ----------
    params : list
            quadratic portion parameters

    Returns
    -------
    resid : array of residuals between O-C 1 and quadratic portion of fit

    '''
    resid = y1 - oc2_func(x, params)
    return resid

@njit
def opt_lsq_fxn(params, pulp, y):
    '''
    Least squares uses an array of residuals to optimize parameters
    Used in : calculating y3 (subtracting periodic portion)

    Parameters
    ----------
    params : list
            periodic portion parameters

    Returns
    -------
    resid : array of residuals between O-C 2 and periodic portion of fit

    '''
    resid = y - oc_func(xi, pulp, params)
    return resid

def main():
    '''
    main function evaluates O-C 1 and O-C 2 using least squares,
    a final O-C using rbfopt and least squares

    Outputs: 
        graphs of each O-C diagram
        status and setting or each optimization algorithm
        guesses value of the parameters

    '''
    
    output_parameters = np.zeros([7])       # Initial values array of parameters 
                                            # that we care about
    
    # Finding the linear parameters
    two_guesses = np.array([0,1])
    two_lower_bounds = np.array([-1e6,1e-3])
    two_upper_bounds = np.array([1e6, 10])
    
    two_params = least_squares(opt_lsq_fxn1, two_guesses, bounds=(
        two_lower_bounds, two_upper_bounds))
    
    output_parameters[0:2] = two_params.x   # Save the linear parameters
    
    # Finding the quadratic parameter
    y1 = yi - oc1_func(xi, two_params.x)    # Creating O-C 1
    pulsation = two_params.x[1]
    B_guess = np.array([.5e-6,-45000,-40])
    B_lower_bound = np.array([-.1E-5,-1e6,-100])
    B_upper_bound = np.array([.1E-5, 1e6, 100])
    
    B_param = least_squares(opt_lsq_fxn2, B_guess, bounds=(
        B_lower_bound, B_upper_bound), args=(y1, xi))
    output_parameters[2] = B_param.x[0]
    
    y2 = y1 - oc2_func(xi, B_param.x)       # Creating O-C 2
    
    @njit
    # This wrap allows rbf to use the optimizing function in the form we have
    # Errors with parameter sizes may have something to do with this
    def opt_wrap(par1):
        temp = optimizing_function(par1, pulsation, y2)
        return temp
    
    # rbfopt algorithm to find the four periodic parameters 
    number_of_paramters = 4
    list_of_parameter_lower_bounds = np.array([58,.6,0.015,5e-3])
    list_of_parameter_upper_bounds = np.array([63,.75,6.28319,.1])
    numIterations = 200         # iterations rbfopt will run for
    
    # Creates the boundaries of the rbf function.
    bb = rbfopt.RbfoptUserBlackBox(
                number_of_paramters, list_of_parameter_lower_bounds, 
                list_of_parameter_upper_bounds, 
                
                # This tells rbfopt that all the parameters are real (it can 
                # also handle integer ones which would be [I]).
                np.array(['R'] * number_of_paramters), 
                # your function that you want optimized. Note: this function must 
                # return only 1 number, usually the sum of residuals.
                opt_wrap
                ) 
    
    # Creates the settings of the rbf function.
    settings = rbfopt.RbfoptSettings(
            # This sets the maximum number of parameter guesses before it quits. If
            # it cannot improve at some point, it will quit before. This is 
            # just to keep it from running forever.
            
            # The number of iterations should statistically be at least 200*(N+1)
            # for N parameters.
            # I've found it works best to allow slightly more evaluations.
            # max_evaluations = int(1.08*numIterations),
            # RBFopt has built-in parallelization. For a GCC laptop, use 2 or 4. 
            num_cpus=2,
            max_iterations = numIterations,
            max_evaluations = numIterations
            )
    
    # Creates an algorithm from the settings and boundaries.
    alg = rbfopt.RbfoptAlgorithm(settings, bb) 
    # optimizes the algorithm.
    (value_of_optimizing_function_at_best_parameters, optimized_parameters, 
                itercount, evalcount, fast_evalcount) = alg.optimize()
    
    print('___________RBF OUTPUT____________')
    print(value_of_optimizing_function_at_best_parameters, optimized_parameters)
    
    
    import sys
    print()
    print("___________STD SETTINGS__________")
    settings.print(sys.stdout)
 
    
    par = optimized_parameters
    lower_bounds = np.array([55, .6, 1e-3, 1e-3])
    upper_bounds = np.array([63, .75, 6.28319, .1])
    res = least_squares(opt_lsq_fxn, par, args=(pulsation, y2),
                        bounds = (lower_bounds, upper_bounds))
    output_parameters[3:7] = res.x          # Save the four period parameters
    
    y3 = y2 - oc_func(xi, pulsation, res.x) # Create the final O-C
    
    # This line should fit O-C 2 perfectly (currently it does not)
    generated_line = oc_func(xi, pulsation, output_parameters[3:7])
    
    print(res)
    print()
    print('Output parameters: {}'.format(output_parameters))
    
    # From here, this is all plotting stuff
    fig, axs = plt.subplots(2,2)
    fig.suptitle("O-C Plots for {}".format(starname))
    #Inital data
    axs[0,0].plot(xi, yi, c='firebrick', marker='.',linestyle='None')
    axs[0,0].set_title('Inital Data',weight='bold')
    
    #Linear subtracted O-C 1
    axs[0,1].plot(xi, y1, c='firebrick', marker='.',linestyle='None')
    axs[0,1].set_title('O-C 1',weight='bold')
    
    #Parabolic subtracted O-C 2
    axs[1,0].plot(xi, y2, c='firebrick', marker='.',linestyle='None')
    axs[1,0].set_title('O-C 2',weight='bold')
    axs[1,0].plot(xi, generated_line)
    
    #Periodic subtracted Final O-C
    axs[1,1].plot(xi, y3, c='firebrick', marker='.',ms=4,linestyle='None')
    axs[1,1].set_title('Final O-C',weight='bold')
    
    for ax in axs.flat:
        ax.set(xlabel='Epochs', ylabel='Days')
        
    fig.tight_layout()
    
    # These lines save the figures as pickle objects (interactable later)
    # and as PNGs
    fig_location = r'C:\Users\hildboldcf17\OneDrive - Grove City College\Senior Year\The Last Semester\RR Lyrae Research\OC Plots\{}'.format(starname)
    fig.savefig('{} OC'.format(fig_location)) # This line saves the fig as a png
    pickle.dump(fig, open('{}.fig.pickle'.format(fig_location), 'wb')) # This line saves the fig as a binary file; interactable and extractable

    
if __name__ == '__main__':
    main()