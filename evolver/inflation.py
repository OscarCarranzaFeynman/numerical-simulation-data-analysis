"""
inflation.py, this is just a name you can use this for other models

Describes a dark matter model
"""
from math import pi, sqrt, sin, cos
import numpy as np
from evolver.constants import Mpl
from evolver.interpolation import Gamma_toth2 
from evolver.parameters import get_portal_params


def SigmavDMtoSMlepton(mx, mf, beta, mh1, mh2, vEW, gX):
    """
    The part of <sigma v>_DMtoSM that depends on the particular final state fermion
    """
    if mf >= mx:
        return 0.0
    return ( sqrt(1-(mf/mx)**2) * (mf**2) * (sin(2.0*beta)*sin(2.0*beta)) * (mh1**2-mh2**2)**2 * (mx**2-mf**2)* (gX**2) *1/ (8.0*pi*(vEW**2)*3.0 * (mh1**2-4.0*mx**2)**2 * (mh2**2-4.0*mx**2)**2 ) )

def SigmavDMtoSMquark(mx, mf, beta, mh1, mh2, vEW, gX):
    """
    The part of <sigma v>_DMtoSM that depends on the particular final state fermion
    """
    if mf >= mx:
        return 0.0
    return sqrt(1-(mf/mx)**2) * (mf**2) * (sin(2.0*beta)*sin(2.0*beta)) * (mh1**2-mh2**2)**2 * (mx**2-mf**2)* (3.0*gX**2) * 1/ (8.0*pi*(vEW**2) * (mh1**2-4.0*mx**2)**2 * (mh2**2-4.0*mx**2)**2 )    


class DarkMatterModel(object):#In this case DarkMatterModel means anymodel, this is just an abstract class
    def __init__(self, **kwargs): #This dictionary is never used here but imposes a condition to declare it in subclasses, the same is for the rest of attributes
        """
        Initialize the model
        """
        raise NotImplementedError()
        
    def SigmavXiXitoh2h2(self, Ydm):
        """
        <sigma v> thermal averaged cross section of 2-2 self-annihilations Xi Xi -> h2 h2
        """
        raise NotImplementedError()
        
    def SigmavXih2toXjXk(self, Ydm):
        """
        <sigma v> thermal averaged cross section of 2-2 self-annihilations Xi h2 -> Xj Xk
        """
        raise NotImplementedError()
         
    def SigmavXiXiXitoXih2(self, Ydm):
        """
        <sigma v> thermal averaged cross section of 3-2 self-annihilations Xi Xi Xi -> Xi h2
        """
        raise NotImplementedError()
         
    def SigmavXiXiXitoXjXk(self, Ydm):
        """
        <sigma v> thermal averaged cross section of 3-2 self-annihilations Xi Xi Xi -> Xj Xk
        """
        raise NotImplementedError()
         
    def SigmavXiXiXjtoXiXk(self, Ydm):
        """
        <sigma v> thermal averaged cross section of 3-2 self-annihilations Xi Xi Xj -> Xi Xk
        """
        raise NotImplementedError()  
   
    def SigmavXiXjXjtoXih2(self, Ydm):
        """
        <sigma v> thermal averaged cross section of 3-2 self-annihilations Xi Xj Xj -> Xi h2
        """
        raise NotImplementedError()  
        
    def SigmavXiXjXktoXiXi(self, Ydm):
        """
        <sigma v> thermal averaged cross section of 3-2 self-annihilations Xi Xj Xk -> Xi Xi
        """
        raise NotImplementedError()                    
            
    def GammaHtoInv(self, Ydm):
        """
        Gamma_h1--> Xi Xi (Higgs to invisible) decays
        """
        raise NotImplementedError()
        
  #  def GammaHtoInvh1h2h2(self, Ydm):
  #      """
   #     Gamma_h1--> h2 h2 (Higgs to invisible) decays
   #     """
   #     raise NotImplementedError()    
        
    def SigmavDMtoSM(self, Ydm):
        """
        <sigma v> thermal averaged cross section of phi phi rightarrow SM SM annihilations
        """
        raise NotImplementedError()
        
    def info(self):
        """
        Returns a string to write to file describing the model
        """
        raise NotImplementedError()
        
        
#class HiddenVector(DarkMatterModel):, you can change HiddenVector for another model
class HiddenVector(DarkMatterModel):
    def __init__(self, **kwargs):
        self.mx = kwargs['mx']
        self.mh2 = kwargs['mh2'] #this is the dark scalar of the model 
        self.gX = kwargs['gX']
        self.lambdaH = kwargs['lambdaH'] #Lambada Higss SM
        self.lambdaPhi = kwargs['lambdaPhi']
        self.lambdaM = kwargs['lambdaM']
        self.muPhi= kwargs['muPhi'] #Mass of mh2
        self.beta= kwargs['beta'] #Mixing angle
        self.mh1 = kwargs['mh1']
        self.vEW = kwargs['vEW']
        self.vPhi = kwargs['vPhi']
        self.me = kwargs['me']
        self.mmu = kwargs['mmu']
        self.mtau = kwargs['mtau']
        self.mup = kwargs['mup']
        self.mdown = kwargs['mdown']
        self.mcharm = kwargs['mcharm']
        self.mstrange = kwargs['mstrange']
        self.mtop = kwargs['mtop']
        self.mbottom = kwargs['mbottom']
        self.XiXitoh2h2_mode = kwargs['XiXitoh2h2_mode']

    def Gammah2toee(self,Ydm,m,me):
        if m<= 2*me:
            return 0.0
        
        return me**2*(np.sin(self.beta)**2)*m*(1-4*(me/m)**2)**1.5/(8*np.pi*self.vEW**2)
            
    
    def Gammah2_total(self, Ydm,m):
        """
        Total decay width of h2 into SM (hadrons + leptons),
        obtained from the digitized plot (plot-data.csv).
        Units: GeV.
        """
        return float(Gamma_toth2(m))       
            
    def SigmavXiXitoh2h2(self, Ydm):
        """
        <σv> thermal averaged cross section of 2-2  Xi Xi -> h2 h2 or h2 h2 -> Xi Xi 
        """
        s= self.XiXitoh2h2_mode

        if s == "forward":
            # Kinematic condition
            if self.mx <= self.mh2:
                return 0.0
    
            num_sqrt = sqrt(1.0 - (self.mh2**2)/(self.mx**2))
            num_inner = (self.gX**2 * self.mh2**2 * (self.mh2**2 - 16.0*self.mx**2) 
                     - 24.0*self.mx**2 * self.lambdaPhi * (self.mh2**2 - 4.0*self.mx**2))**2
            denom = 1536.0 * pi * self.mx**2 * (self.mh2**4 - 20.0*self.mh2**2*self.mx**2 + 64.0*self.mx**4)**2
    
            return num_sqrt * num_inner / denom
            
        elif s == "reverse":
            # Kinematic condition
            if self.mx >= self.mh2:  
                return 0.0
        
            num_sqrt = sqrt(1.0 - (self.mx**2)/(self.mh2**2))
            poly = (4*self.mh2**4 - 4.0*self.mh2**2*self.mx**2 + 3.0*self.mx**4)
            inner = (self.gX**2 * (self.mh2**2 - 4.0*self.mx**2) + 8.0*self.mx**2*self.lambdaPhi)**2
            denom = 512.0 * pi * (self.mh2**6) * (self.mx**4)

            return num_sqrt * poly * inner / denom

        else: 
            raise ValueError(f"Invalid option for s: {s}")
        
    def SigmavXih2toXjXk(self, Ydm):
       """
       <sigma v> thermal averaged cross section of 2-2 self-annihilations Xi h2 -> Xj Xk
       """
       # Kinematic condition
       if self.mh2 <=self.mx:
           return 0.0
       
       root = 1.0 - (4.0*self.mx**2) / ((self.mh2 + self.mx)**2)
       num_poly = (self.mh2**6 + 4.0*(self.mh2**5)*self.mx - 4.0*(self.mh2**4)*self.mx**2 - 10.0*(self.mh2**3)*(self.mx**3)
                + 144.0*(self.mh2**2)*(self.mx**4) + 396.0*self.mh2*self.mx**5 + 297.0*self.mx**6)
       num = self.gX**4 * (self.mh2 - self.mx) * (self.mh2 + 3.0*self.mx) * num_poly * sqrt(root)
       den = 384.0 * pi * self.mh2**3 * self.mx**5 * (self.mh2 + 2.0*self.mx)**2

       return num / den
       
    def SigmavXiXiXitoXih2(self, Ydm):
        """
        <sigma v> thermal averaged cross section of 3-2 self-annihilations Xi Xi Xi -> Xi h2
        """
        root = (self.mh2**4)/(self.mx**4) - 20.0*(self.mh2**2)/(self.mx**2) + 64.0
        if root <= 0.0:
            return 0.0

        return 5.0 * sqrt(root) * (self.mh2**4 - 20.0*self.mh2**2*self.mx**2 + 172.0*self.mx**4) * (
            self.gX**3 * (-11.0*self.mh2**4 + 10.0*self.mh2**2*self.mx**2 + 64.0*self.mx**4)
            + 144.0 * self.gX * self.lambdaPhi * self.mx**2 * (self.mh2**2 - 4.0*self.mx**2))**2 / (71663616.0 * pi * self.mx**5 * (self.mh2**2 -               4.0*self.mx**2)**4 * (self.mh2**2 + 2.0*self.mx**2)**2)

    def SigmavXiXiXitoXjXk(self, Ydm):
        """
        <sigma v> thermal averaged cross section of 3-2 self-annihilations Xi Xi Xi -> Xj Xk
        """
        return  (25.0 * sqrt(5.0) * self.gX**6 *
            (8375.0*self.mh2**4 - 66638.0*self.mh2**2*self.mx**2 + 134265.0*self.mx**4)
           ) / (2654208.0 * self.mx**5 * (self.mh2**2 - 4.0*self.mx**2)**2 * pi)
        
    def SigmavXiXiXjtoXiXk(self, Ydm):
        """
        <sigma v> thermal averaged cross section of 3-2 self-annihilations Xi Xi Xj -> Xi Xk
        """
        return (5.0 * sqrt(5.0) * self.gX**6 *
            (14377.0*self.mh2**8 - 85320.0*self.mh2**6*self.mx**2 + 22816.0*self.mh2**4*self.mx**4
             + 299340.0*self.mh2**2*self.mx**6 + 346027.0*self.mx**8)
           ) / (2654208.0 * self.mx**5 * (self.mh2**4 - 3.0*self.mh2**2*self.mx**2 - 4.0*self.mx**4)**2 * pi)
        
    def SigmavXiXjXjtoXih2(self, Ydm):
        """
        <sigma v> thermal averaged cross section of 3-2 self-annihilations Xi Xj Xj -> Xi h2
        """
        root_arg = 64.0 + (self.mh2**4)/(self.mx**4) - 20.0*(self.mh2**2)/(self.mx**2)
        if root_arg <= 0.0:
            return 0.0
        return (
            self.gX**2
            * sqrt(64.0 + (self.mh2**4)/(self.mx**4) - 20.0*(self.mh2**2)/(self.mx**2))
            * (
            self.gX**4 * (self.mh2**2 + 2.0*self.mx**2)**2
            * (
                12.0*self.mh2**16 - 688.0*self.mh2**14*self.mx**2
                + 41335.0*self.mh2**12*self.mx**4 - 1072062.0*self.mh2**10*self.mx**6
                + 13436087.0*self.mh2**8*self.mx**8 - 87386228.0*self.mh2**6*self.mx**10
                + 312481572.0*self.mh2**4*self.mx**12 - 611437312.0*self.mh2**2*self.mx**14
                + 538513408.0*self.mx**16
            )
            + 288.0 * self.gX**2 * (self.mh2 - 2.0*self.mx) * self.mx**4 * (self.mh2 + 2.0*self.mx)
            * (self.mh2**2 - 7.0*self.mx**2) * (self.mh2**2 + 2.0*self.mx**2)
            * (
                4.0*self.mh2**10 - 173.0*self.mh2**8*self.mx**2 + 4459.0*self.mh2**6*self.mx**4
                - 44968.0*self.mh2**4*self.mx**6 + 143732.0*self.mh2**2*self.mx**8 - 88960.0*self.mx**10
            ) * self.lambdaPhi
            + 62208.0 * self.mx**8
            * (self.mh2**4 - 11.0*self.mh2**2*self.mx**2 + 28.0*self.mx**4)**2
            * (self.mh2**4 - 20.0*self.mh2**2*self.mx**2 + 172.0*self.mx**4)
            * (self.lambdaPhi**2)
            )
            ) / (
            214990848.0 * self.mx**9 * (self.mh2**2 - 4.0*self.mx**2)**4
            * (self.mh2**4 - 5.0*self.mh2**2*self.mx**2 - 14.0*self.mx**4)**2 * pi
            )   
        
    def SigmavXiXjXktoXiXi(self, Ydm):
        """
        <sigma v> thermal averaged cross section of 3-2 self-annihilations Xi Xj Xk -> Xi Xi
        """
        num = 5.0 * sqrt(5.0) * self.gX**6 * (
        347.0 * self.mh2**4 +
        586.0 * self.mh2**2 * self.mx**2 +
        707.0 * self.mx**4
        )
        den = 331776.0 * self.mx**5 * (self.mh2**2 + self.mx**2)**2 * pi

        return num / den     
    
        
    def GammaHtoInv(self, Ydm):
        """
        Gamma_h1--> Xi Xi (Higgs to invisible) decays
        """
        if self.mh1<=2.0*self.mx:
            return 0.0

        num_sqrt = sqrt(self.mh1**2 - 4.0*self.mx**2)
        num_poly = (self.mh1**4 - 4.0*self.mh1**2*self.mx**2 + 12.0*self.mx**4)
        num = (np.sin(self.beta))**2 * (self.gX**2) * num_sqrt * num_poly
        den = 128.0 * pi * self.mh1**2 * self.mx**2

        return num / den
        
        
    def SigmavDMtoSM(self, Ydm):
        """
        <sigma v> thermal averaged cross section of phi phi rightarrow SM SM annihilations
        """
        den=( ((self.mx+self.mx)**2-self.mh2**2)**2+(self.mh2**2)*(self.Gammah2_total(Ydm,self.mh2)*np.sin(self.beta)**2+self.Gammah2toee(Ydm,self.mh2,self.me))**2 )
        return ( 4/np.sqrt((self.mx+self.mx)**2)*(self.gX**2/2*2*self.mx/self.gX*np.cos(self.beta))**2*(1/den)*(self.Gammah2_total(Ydm,2*self.mx)*np.sin(self.beta)**2+self.Gammah2toee(Ydm,2*self.mx,self.me)))                   
                                                                                       
    
    def info(self):
        """Returns a multi-line string describing the HiddenVector model."""
        lines = []
    
        # Identificación del modelo
        lines.append("Model: HiddenVector")
    
        # Masas (asume GeV en tu input)
        lines.append(
            "Masses [GeV]: "
            f"mx={self.mx}, mh1={self.mh1}, mh2={self.mh2}"
        )
    
        # Acoplamientos del sector oscuro / portal
        lines.append(
            "Couplings: "
            f"gX={self.gX}, lambdaH={self.lambdaH}, lambdaPhi={self.lambdaPhi}, lambdaM={self.lambdaM}"
        )
    
        # Escalas y mezcla
        sB = np.sin(self.beta)
        cB = np.cos(self.beta)
        lines.append(
            "Scales [GeV]: "
            f"vEW={self.vEW}, vPhi={self.vPhi}, muPhi={self.muPhi}"
        )
        lines.append(
            "Mixing: "
            f"beta={self.beta}, sin(beta)={sB}, cos(beta)={cB}"
        )
    
        # Masas de fermiones SM (para canales a SM)
        lines.append(
            "SM fermion masses [GeV]: "
            f"me={self.me}, mmu={self.mmu}, mtau={self.mtau}, "
            f"mu={self.mup}, md={self.mdown}, ms={self.mstrange}, "
            f"mc={self.mcharm}, mb={self.mbottom}, mt={self.mtop}"
        )
    
        return "\n".join(lines)