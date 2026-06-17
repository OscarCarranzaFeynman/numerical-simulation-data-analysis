"""
eoms.py

Computes the equations of motion associated with a model
"""
from math import pi, sqrt
import numpy as np
from scipy.special import kn
from evolver.constants import Mpl


class EOMParameters(object):
    """
    A class to store all the parameters that are needed for computing the EOMs, and also for plotting
    """
    #add more parameters later as needed
    def __init__(self, infmodel):
        """Save parameters for later computations"""
        self.model = infmodel
    
    
def eoms(unpacked_data, params, time):
    """Workhorse to compute all equations of motion.
       This gets used by derivatives() in model.py to engine the integrator
    
    Arguments:
        * unpacked_data: Tuple containing all data solved by ODE integrator:
                            (Ydm, DTemp)
        * params: EOMParameters object
        * time: The current x (needed for EOMs)
    
    Returns a tuple of derivatives:
        * (dYdmdx, dDTempdx, Others)
    """
    #Initialization
    Ydm = unpacked_data
    
    #compute rates and thermodynamic quantities, you can look file 7 in boltzmann equations notes
    (Rate3to2, RateXiXitoh2h2, RateXih2toXjXk, Rate3toh2, RateHtoInv, RateDMtoSM, neqDM, neqh1, neqh2, Hubble, sent, Corr) = compute_all(unpacked_data, params, time)
    
    #derivatives
    dYdmdx = compute_boltzmann(unpacked_data, neqDM, neqh1, neqh2, Hubble, sent, Corr, params, time)
    
    #Return results
    return (dYdmdx, Rate3to2, RateXiXitoh2h2, RateXih2toXjXk, Rate3toh2, RateHtoInv, RateDMtoSM, Hubble)
    
    
def compute_all(unpacked_data, params, time):
    """
    Routine to compute all quantities of interest
    All quantities passed to compute_boltzmann, rates used for conditions and writing to file
    
    Arguments:
        *unpacked_data: Tuple containing all data solved by ODE integrator:
                        (Ydm, DTemp)
        *params: EOMParameters object
    """
    #Initialization
    Ydm = unpacked_data
    
    #Compute thermodynamic quantities
    neqDM = compute_neqDM(time, params.model)
    neqh1 = compute_neqh1(time, params.model)
    neqh2 = compute_neqh2(time, params.model)
    Hubble = compute_Hubble(time, params.model)
    sent = compute_sent(time, params.model)
    Corr = compute_Corr(time, params.model)
    
    Rate3to2 = compute_Rate3to2(time, params.model, Ydm)
    RateXiXitoh2h2 = compute_RateXiXitoh2h2(time, params.model, Ydm)
    RateXih2toXjXk = compute_RateXih2toXjXk(time, params.model, Ydm)
    Rate3toh2 = compute_Rate3toh2(time, params.model, Ydm)
    RateHtoInv = compute_RateHtoInv(time, params.model, Ydm)
    RateDMtoSM = compute_RateDMtoSM(time, params.model, Ydm)
    
    #Return results
    return Rate3to2, RateXiXitoh2h2, RateXih2toXjXk, Rate3toh2, RateHtoInv, RateDMtoSM, neqDM, neqh1, neqh2, Hubble, sent, Corr
    
    

def compute_boltzmann(unpacked_data, neqDM, neqh1, neqh2, Hubble, sent, Corr, params, time):
    """
    Compute the RHS of the Boltzmann equation for evolution of the DM number density, to be solved by the ODE integrator later down the pipeline
    """
    #Initialization
    Ydm = unpacked_data
    model = params.model

    s = params.model.XiXitoh2h2_mode
    if s == "forward":
        term_XiXitoh2h2 = (2.0/3.0) * model.SigmavXiXitoh2h2(Ydm) * sent * (-Ydm**2 + (neqDM/sent)**2)
    elif s == "reverse":
     
        term_XiXitoh2h2 = 6.0 * model.SigmavXiXitoh2h2(Ydm) * (neqh2**2) * (1.0/sent) * (1.0 - (Ydm/(neqDM/sent))**2)
    else:
        raise ValueError("s debe ser 'forward' o 'reverse'")
    
    return ( (1/(Hubble*time*Corr))*( term_XiXitoh2h2 + neqh2*model.SigmavXih2toXjXk(Ydm)*Ydm*(1-Ydm/(neqDM/sent)) + (sent**2)*( -Ydm**3+(Ydm**2)*(neqDM/sent) )* ( 1/9*model.SigmavXiXiXitoXjXk(Ydm) + 2/9*model.SigmavXiXiXjtoXiXk(Ydm)+ 1/9*model.SigmavXiXjXktoXiXi(Ydm)) +(sent**2)*(-Ydm**3+Ydm*(neqDM/sent)**2)*(2/9*model.SigmavXiXiXitoXih2(Ydm)+4/9*model.SigmavXiXjXjtoXih2(Ydm))+2/3*model.SigmavDMtoSM(Ydm)*(sent)*(-Ydm**2+(neqDM/sent)**2)  ))
             
                #+2*model.GammaHtoInv(Ydm)*(neqh1)/(sent*time)*(1.0 - (Ydm / (neqDM/sent))**2)
def eoms_terms(unpacked_data, params, time):

    Ydm = unpacked_data
    model = params.model

    neqDM = compute_neqDM(time, model)
    neqh1 = compute_neqh1(time, model)
    neqh2 = compute_neqh2(time, model)
    Hubble = compute_Hubble(time, model)
    sent   = compute_sent(time, model)
    Corr   = compute_Corr(time, model)

    # XiXi <->h2 h2 
    s = model.XiXitoh2h2_mode
    if s == "forward":
        t_XiXitoh2h2_raw = (2.0/3.0) * model.SigmavXiXitoh2h2(Ydm) * sent * (-Ydm**2 + (neqDM/sent)**2)
    elif s == "reverse":
      
        t_XiXitoh2h2_raw = 6.0 * model.SigmavXiXitoh2h2(Ydm) * (neqh2**2) * (1.0/sent) * (1.0 - (Ydm/(neqDM/sent))**2)
    else:
        raise ValueError("s debe ser 'forward' o 'reverse'")


    # Xi h2 -> Xj Xk
    t_Xih2toXjXk_raw = neqh2 * model.SigmavXih2toXjXk(Ydm) * Ydm * (1.0 - Ydm/(neqDM/sent))

    #  3->2: XiXiXi -> XjXk  
    combo_3to2 = ((1.0/9.0)*model.SigmavXiXiXitoXjXk(Ydm) + (2.0/9.0)*model.SigmavXiXiXjtoXiXk(Ydm) + (1.0/9.0)*model.SigmavXiXjXktoXiXi(Ydm))
    t_3to2_XX_raw = (sent**2) * (-Ydm**3 + (Ydm**2)*(neqDM/sent)) * combo_3to2

    # 3->2  XiXiXi -> Xi h2  (dos sigmas)
    combo_3toh2 = ((2.0/9.0)*model.SigmavXiXiXitoXih2(Ydm) + (4.0/9.0)*model.SigmavXiXjXjtoXih2(Ydm) )
    t_3to2_Xh2_raw = (sent**2) * (-Ydm**3 + Ydm*(neqDM/sent)**2) * combo_3toh2

    # DM -> SM
    t_DMtoSM_raw = (2.0/3.0) * model.SigmavDMtoSM(Ydm) * sent * (-Ydm**2 + (neqDM/sent)**2)

    terms_raw = {
        "XiXi_to_h2h2_raw": t_XiXitoh2h2_raw,
        "Xih2_to_XjXk_raw": t_Xih2toXjXk_raw,
        "3to2_to_XX_raw":   t_3to2_XX_raw,
        "3to2_to_Xh2_raw":  t_3to2_Xh2_raw,
        "DM_to_SM_raw":     t_DMtoSM_raw
    }

    return terms_raw


def compute_Rate3to2(time, model, Ydm):
    """
    Compute the per-DM particle number-changing rate of 3-->2 self-annihilations
    """
    #Here we write simmetry factors as they appear in the equation
    return ( 1/9*model.SigmavXiXjXktoXiXi(Ydm) + 2/9*model.SigmavXiXiXjtoXiXk(Ydm) + 1/9*model.SigmavXiXiXitoXjXk(Ydm))*(compute_sent(time,model)**2)*(Ydm**2)
    
def compute_RateXiXitoh2h2(time, model, Ydm):
    """
    Compute the per-DM particle number-changing rate of 3-->2 self-annihilations
    """
    #Here we write simmetry factors as they appear in the equation, 2/3*model.SigmavXiXitoh2h2(Ydm)*(compute_sent(time,model)*Ydm
    s = model.XiXitoh2h2_mode  

    if s == "forward":
        # Xi Xi → h₂ h₂
        return (2.0/3.0) * model.SigmavXiXitoh2h2(Ydm) * compute_sent(time, model) * Ydm

    elif s == "reverse":
        # h₂ h₂ → Xi Xi
        return 6.0 * model.SigmavXiXitoh2h2(Ydm)* (compute_neqh2(time, model)**2) / (Ydm * compute_sent(time, model))
    
def compute_RateXih2toXjXk(time, model, Ydm):
    """
    Compute the per-DM particle number-changing rate of 3-->2 self-annihilations
    """
    return model.SigmavXih2toXjXk(Ydm)*(compute_neqh2(time,model))      
    
def compute_Rate3toh2(time, model, Ydm):
    """
    Compute the per-DM particle number-changing rate of 3-->2 self-annihilations
    """
       #Here we write simmetry factors as they appear in the equation
    return ( 2/9*model.SigmavXiXiXitoXih2(Ydm) + 4/9*model.SigmavXiXjXjtoXih2(Ydm) )*(compute_sent(time,model)**2)*(Ydm**2)   

#def compute_RateHtoInv_h1XiXi(time, model, Ydm):
#    return model.GammaHtoInvh1XiXi(Ydm) * compute_neqh2(time,model)/( compute_sent(time,model)*Ydm )

#def compute_RateHtoInv_h1h2h2(time, model, Ydm):
#    return model.GammaHtoInvh1h2h2(Ydm) * compute_neqh2(time,model)/( compute_sent(time,model)*Ydm )
    
    
def compute_RateHtoInv(time, model, Ydm):
    """
    Compute the per-DM particle number-changing rate of h-->phi phi decays
    """
    return model.GammaHtoInv(Ydm) * compute_neqh1(time,model)/( compute_sent(time,model)*Ydm )
    
def compute_RateDMtoSM(time, model, Ydm):
    """
    Compute the per-DM particle number-changing rate of phi phi --> SM SM annihilations
    """
       #Here we write simmetry factors as they appear in the equation
    return 2/3*model.SigmavDMtoSM(Ydm)*( compute_sent(time,model)**2 )*( Ydm**2 )/( compute_sent(time,model)*Ydm )
    
def compute_YeqDM(time, model):
    """
    Compute the equilibrium abundance Y
    """
    return compute_neqDM(time,model)/compute_sent(time,model)
    
    
def compute_neqDM(time, model):
    """
    Compute equilibrium number densities
    """
    factor = 1/(2*pi**2)
    return factor*3* model.mx**2 * (model.mx/time) * kn(2, time) #We should set it 3
    
def compute_neqh1(time, model):
    """
    Compute equilibrium number densities
    """
    factor = 1/(2*pi**2)
    return factor * model.mh1**2 * (model.mx/time) * kn(2, (model.mh1/model.mx)*time)    #Compare with page 106 of Patrick's Thesis
def compute_neqh2(time, model):
    """
    Compute equilibrium number densities
    """
    factor = 1/(2*pi**2)
    return factor * model.mh2**2 * (model.mx/time) * kn(2, (model.mh2/model.mx)*time)
    
def compute_Hubble(time, model):
    """
    Compute the Hubble parameter for a given x
    """
    #Mpl is sqrt(1/(8*Pi*G))
    factor = 1/(3*Mpl**2)
    rho = compute_rho(model.mx/time)
    
    return sqrt(factor * rho)
        
def compute_sent(time, model):
    """
    Compute the entropy density of the Universe
    """
    factor = 2*pi**2/45
    return factor * heffBWInt(model.mx/time) * (model.mx/time)**3
    
    
def compute_Corr(time, model):
    """
    Compute the correction factor from thermodynamic degrees of freedom
    """
    #for now setting it to 1
    return 1
        
def compute_rho(T):
    """
    Compute the energy density of the Universe for a given T
    """
    factor = pi**2/30
    return factor*geff(T)*T**4
    
def geff(T):
    """
    Effective number of degrees of freedom contributing to the energy density
    """
    #for now setting it to its late-time value
    return 100
    
def heffBWInt(T):
    """
    Effective number of degrees of freedom contributing to the entropy density
    """
    #for now setting it to its late-time value
    return 100
