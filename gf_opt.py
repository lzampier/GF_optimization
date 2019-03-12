""" Optimize parameters with Green's Functions
 necessary data:
 y - vector of observations
 sigma - vector of observational errors
 x - vector of model data of baseline run
 G_j - vector of model data of perturbed run
"""
from __future__ import print_function

import numpy as np
import example_gf_opt as ex

def calcDataKernel(x, x_pert, pert_vals): 
    """ Calculate Data Kernel Matrix G
    each column is the approximation of Green's Function for one perturbation experiment:
    the difference (perturbed - baseline), divided by the magnitude of the perturbation

    :param x: original model values
    :param x_pert: list of arrays with model values from perturbation experiments
    :param pert_vals: list of magnitudes of parameter perturbations
    """ 
    n_y = max(x.shape)
    n_pert = len(pert_vals)

    G = np.zeros((n_y,n_pert))
    for i in range(n_pert):
        G[:,i] = (x_pert[i] - x) / pert_vals[i]
    
    return G


def optimiseGF(x, y, sigma, G): 
    """ Calculates set of optimal parameter perturbations.
    Attention: eta_new must be added to the original parameter values!
    
    Solves:
    $ /eta_{opt} = (G^T R^{-1} G)^{-1} G^T R^{-1} /Delta y $,
    where R is diagonal with measurement uncertainties sigma on the main diagonal
    and $\Delta y$ is the model error.

    :param G: Data Kernel Matrix (array)

    :returns: optimal parameter perturbation
    """ 
    npert = G.shape[1]
    y_diff = (y - x)
    R_diag = (1 / np.square(sigma))

    # P = (G^T R^{-1} G)
    P = np.zeros((npert, npert))
    for i in range(npert):
        for j in range(npert):
            P[i, j] = np.dot((R_diag * G[:, i]), G[:, j] )

    aux_vec = np.dot(G.transpose(), R_diag * y_diff)
    eta_new = np.dot(np.linalg.inv(P), aux_vec)
    
    return eta_new


def costFunction(x, y, sigma): 
    """ Calculate the quadratic cost function F that is minimized.
    F = \sum_i ( (y_i - x_i) / sigma_i)**2

    :param x: model values
    :param y: observed values
    :param sigma: uncertainties for the observations

    :returns: Cost function value F(x, y, sigma)
    """
    y_diff = (y - x)
    R_diag = (1 / np.square(sigma))
    res = np.dot(np.square(y_diff), R_diag)

    return res


def printParms(parm_names, orig_parms, pert, eta_new): 
    """ Print the optimised parameters.

    :param parm_names: iterable of names of tested parameters
    :param orig_parms: iterable of original values
    :param pert: iterable of perturbations used in the sensitivity experiments
    :param eta_new: iterable of optimized perturbations
    """ 
    print(" *** parameter estimations *********************")
    print('{0:21}'.format('Parameters:'), end=',')
    for s in parm_names:
        print('{0:>10}'.format(s), end=',')
    print()
    print('{0:21}'.format('original values:'), end=',')
    for s in orig_parms:
        print('{0:10.5}'.format(s), end=',')
    print()
    print('{0:21}'.format('perturbed values:'), end=',')
    for s in pert+orig_parms:
        print('{0:10.5}'.format(s), end=',')
    print()
    print('{0:21}'.format('optimised values:'), end=',')
    for s in (orig_parms+eta_new):
        print('{0:10.5}'.format(s), end=',')
    print()
    
if __name__ == '__main__':
    ## set up choice of validation data and simulations
    dataflags = ex.Dataflags()
    dataflags.flag_strgf_Pdist = 1
    dataflags.use_osiconc=0 # produces MemoryError (19/03/12)

    # read data
    (basel_path, pert_path, optim_path, toolpath, orig_parms, pert, parm_names,
    ) = ex.fill_info(dataflags)
    x, y, x_pert, sigma, pertvals, nr_sigma = ex.getData(dataflags)

    # caculate optimal parameters
    G = calcDataKernel(x, x_pert, pertvals)
    eta_new = optimiseGF(x, y, sigma, G)
    new_parms = orig_parms + eta_new

    # output
    printParms(parm_names,orig_parms,pert,eta_new)
