"""
utilities.py

Helpful utilities
"""
import numpy as np

def pack(Ydm):
    """
    Pack all field values into a data array for integration.
    
    Arguments: 
        * Ydm: Respective values to pack
    
    Returns:
        * data: A numpy array containing all the data
    """
    background = np.array([Ydm])
    return np.concatenate([background])
    
def unpack(data):
    """
    Unpack field values from a data array into a meaningful data structure.
    This reverses the operations performed in pack.
    
    Arguments:
        * data: The full data array of all fields, their derivatives, and auxiliary values
    
    Returns: 
        * (Ydm)
          where these quantities are initialized as in make_initial_data
    """
    Ydm = data[0]
    
    #Return the results
    return Ydm
    
def load_data(filename):
    """Loads data from a run"""
    with open(filename + ".dat") as f: #With this notation the files is closed when block witch ends
        data = f.readlines()
    with open(filename + ".dat2") as f:
        data2 = f.readlines()
    with open(filename + ".info") as f:
        data3 = f.readlines()
        
    #Process the data
    results = np.array([list(map(float, line.split(", "))) for line in data]).transpose()
    results2 = np.array([list(map(float, line.split(", "))) for line in data2]).transpose()
    #unpack into nicer variables
    t = results[0]
    Ydm = unpack(results[1:])
    #Change for the variables, rates, etc which you want to track
    (Rate3to2, RateXiXitoh2h2, RateXih2toXjXk, Rate3toh2, RateHtoInv, RateDMtoSM, Hubble, sent, neqh1, neqh2, neqDM) = results2[1:]
    
    # Put the data into a container
    #Change here too
    fulldata = {
        "t": t,
        "Ydm": Ydm,
        "Rate3to2": Rate3to2,
        "RateXiXitoh2h2": RateXiXitoh2h2,
        "RateXih2toXjXk": RateXih2toXjXk,
        "Rate3toh2": Rate3toh2,
        "RateHtoInv": RateHtoInv,
        "RateDMtoSM": RateDMtoSM,
        "Hubble": Hubble,
        "sent": sent,
        "neqh1": neqh1,
        "neqh2": neqh2,
        "neqDM": neqDM
    }
    
    #Return __everything__
    return fulldata
