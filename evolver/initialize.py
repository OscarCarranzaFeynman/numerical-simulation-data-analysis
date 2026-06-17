"""
initialize.py

Initializes parameters for a run
"""

import numpy as np
from evolver.utilities import pack, unpack
from evolver.eoms import (compute_YeqDM, compute_Hubble, EOMParameters, 
                       compute_Rate3to2, compute_RateXiXitoh2h2, compute_RateXih2toXjXk, compute_Rate3toh2, compute_RateHtoInv, compute_RateDMtoSM, compute_neqDM, compute_neqh1, compute_neqh2, compute_sent, compute_Corr, compute_boltzmann)
                        
def create_package(infmodel,
                    start_time,
                    end_time,
                    basefilename,
                    fulloutput=True,
                    **kwargs):
                    
    """
    Packages all the settings that need to be set for a run into a dictionary.
    This package can then be initialized using create_parameters.
                    
    A package can be used multiple times, and can be easily modified to
    provide different initial conditions.
                    
    Some arguments are mandatory, while others are optional. Any extra
    package settings can be specified by keyword arguments.
                    
    Arguments:
        *Ydm0: Starting value for DM abundance
        *infmodel: Instance of an DarkMattter class
        *start_time: start time to begin evolution
        *end_time: Maximum time to which to evolve
        *basefilename: The base name of the desired output file
                        (different extensions will be added)
        *fulloutput: Do you want the full output from a run, or only the initial/final conditions?
    """
    package = {
        'infmodel': infmodel,
        'start_time': start_time,
        'end_time': end_time,
        'basefilename': basefilename,
        'fulloutput': fulloutput,
        **kwargs
    }
    return package
    
def create_parameters(package):
    """
    Takes in a package and generates the parameters dictionary required by
    a Model. This involves constructing grids and initial conditions. Once this
    dictionary has been generated, it should not be modified, as doing so may
    cause inconsistencies.
    """
    # Pass through a consistency of initial conditions check
    packagecopy = {**package}
    try:
        parameters = _create_parameters(packagecopy)
        return parameters
    except BadK:
        #Exception raised if initial conditions not consistent
        return None
        
class BadK(Exception):
    """An exception that is raised when initial conditions are inconsistent"""
    
def _create_parameters(package):
    #Start by copying everything in the package
    parameters = {**package}
    
    infmodel = parameters['infmodel']
    
    ####Calculate needed quantities to establish initial conditions
    #### and check their consistency here
    
    ###########################
    # Construct EOMParameters #
    ###########################
    parameters['eomparams'] = EOMParameters(infmodel)
    
    ################################
    # Construct initial conditions #
    ################################
    
    # The starting value of x
    xstart = parameters['start_time']
    
    #Starting value of Ydm0
    Ydm = compute_YeqDM(xstart, infmodel)
    
    #Compute the rates of number-changing processes
    #Change for your rates
    Rate3to20 = compute_Rate3to2(xstart, infmodel, Ydm)
    RateXiXitoh2h20 = compute_RateXiXitoh2h2(xstart, infmodel, Ydm)
    RateXih2toXjXk0 = compute_RateXih2toXjXk(xstart, infmodel, Ydm)
    Rate3toh20 = compute_Rate3toh2(xstart, infmodel, Ydm)
    RateHtoInv0 = compute_RateHtoInv(xstart, infmodel, Ydm)
    RateDMtoSM0 = compute_RateDMtoSM(xstart, infmodel, Ydm)
    
    
    #Compute the Hubble rate
    Hubble0 = compute_Hubble(xstart, infmodel)
    
    #Make sure that dark and SM sectors are in equilibrium on the initial conditions
    #Debugging: Initial rates -- equilibrium conditions
   # print("RateHtoInv0/Hubble0:", RateHtoInv0/Hubble0)
   # print("RateHtoInv0:", RateHtoInv0)
   # print("RateDMtoSM0/Hubble0:", RateDMtoSM0/Hubble0)
   # print("Rate3to20/Hubble0:", Rate3to20/Hubble0)
   # print("RateXiXitoh2h20/Hubble0:", RateXiXitoh2h20/Hubble0)
   # print("RateXih2toXjXk0/Hubble0:", RateXih2toXjXk0/Hubble0)
   # print("Rate3toh20/Hubble0:", Rate3toh20/Hubble0)
   # print("Hubble0:", Hubble0)

    if RateHtoInv0/Hubble0 < 1 and RateDMtoSM0/Hubble0 < 1:
        #Initial conditions fail to produce equilibrium between DM and SM
        raise BadK()
        
    #Pack all the initial data together
    parameters['initial_data'] = pack(Ydm)
    
    #Debugging: initial step size
    initial_data = parameters['initial_data']
    unpacked_data = unpack(initial_data)
    neqDM0 = compute_neqDM(xstart, infmodel)
    neqh10 = compute_neqh1(xstart, infmodel)
    neqh20 = compute_neqh2(xstart, infmodel)
    sent0 = compute_sent(xstart, infmodel)
    Corr0 = compute_Corr(xstart, infmodel)
    dYdmdx0 = compute_boltzmann(unpacked_data, neqDM0, neqh10, neqh20, Hubble0, sent0, Corr0, parameters['eomparams'], xstart)
    
    #print("dYdmdx0:",dYdmdx0)
    #print("Ydm0:", Ydm)
    #print("sent:", sent0)

    # 2 to 2
    #print("SigmavXiXitoh2h2:", infmodel.SigmavXiXitoh2h2(Ydm))
    #print("SigmavXih2toXjXk:", infmodel.SigmavXih2toXjXk(Ydm))

    # 3 to 2 
   # print("SigmavXiXiXitoXih2:", infmodel.SigmavXiXiXitoXih2(Ydm))
   # print("SigmavXiXiXitoXjXk:", infmodel.SigmavXiXiXitoXjXk(Ydm))
   # print("SigmavXiXiXjtoXiXk:", infmodel.SigmavXiXiXjtoXiXk(Ydm))
   # print("SigmavXiXjXjtoXih2:", infmodel.SigmavXiXjXjtoXih2(Ydm))
   # print("SigmavXiXjXktoXiXi:", infmodel.SigmavXiXjXktoXiXi(Ydm))
    # Higgs invisible (h1 → Xi Xi)
  #  print("GammaHtoInv:", infmodel.GammaHtoInv(Ydm))
    #print("GammaHtoInvh1XiXi", infmodel.SigmavXiXjXjtoXih2(Ydm))
    #print("GammaHtoInvh1h2h2", infmodel.GammaHtoInvh1h2h2(Ydm))
 #   print("SigmavDMtoSM:", infmodel.SigmavDMtoSM(Ydm))
    
    #Return everything
    return parameters
