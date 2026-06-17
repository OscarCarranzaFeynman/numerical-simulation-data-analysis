"""model.py

Defines the model for evolution, built on top of AbstractModel
"""
import pickle
import time as systime
from math import pi, sqrt, log
from evolver.eoms import eoms, compute_all
from evolver.integrator import AbstractModel
from evolver.utilities import pack, unpack

class Model(AbstractModel):
    """Set up the model for evolution
        Model is a subclass of AbstractModel -- it inherits all the methods and properties of AbstractModel and can override or extend them
        Model implements the interface defined generically by AbstractModel for a specific physical system
    """
    
    def __init__(self, run_params):
        """Set defaults and save parameters for rapid access"""
        #Defaults
        defaults = {
            'timestepinfo': [200,10]
        }
        defaults.update(run_params)
        if not run_params['perform_run']:
            defaults['end_time']=defaults['start_time']
            
        #Call the AbstractModel initializer (call the __init__method of the parent class, passing it the defaults dictionary -- give acess to all methods and attributes of superclass -- initializes Model as an AbstractModel)
        super().__init__(defaults)
        
        #Save parameters for future use (faster dictionary access)
        self.basefilename = self.parameters['basefilename']
        self.timestepinfo = self.parameters['timestepinfo']
        self.infmodel = self.parameters['infmodel']
        self.eomparams = self.parameters['eomparams']
        self.fulloutput = self.parameters['fulloutput']
        self.stiff_solver = self.parameters['stiff_solver']
        
        #Set internal flags
        self.inflationstart = None
        self.slowroll = False
        self.slowrollend = None
        
    def begin(self, time):
        """Open file handles for evolution and write initial description"""
        #Set up file handles for writing data
        if self.fulloutput:
            self.datafile = open(self.basefilename + ".dat", "w")
            self.extrafile = open(self.basefilename + ".dat2", "w")
            
        #Write the initial description to file
        params = self.eomparams
        unpacked_data = unpack(self.initial_data)
        Ydm = unpacked_data
        
        #Change for your rates and quantities
        (Rate3to2, RateXiXitoh2h2, RateXih2toXjXk, Rate3toh2, RateHtoInv, RateDMtoSM, neqDM, neqh1, neqh2, Hubble, sent, Corr) =      compute_all(unpacked_data, params, time)
        
        #Change here too
        with open(self.basefilename + ".info", "w") as f:
            f.write(f"""Evolution Parameters and Initial Conditions
            Model: {type(params.model).__name__}
            {params.model.info()}
            Initial Ydm: {Ydm}
            Initial Hubble: {Hubble}
            Initial Rate3to2: {Rate3to2}
            Initial RateXiXitoh2h2: {RateXiXitoh2h2}
            Initial RateXih2toXjXk: {RateXih2toXjXk}
            Initial Rate3toh2: {Rate3toh2}
            Initial RateHtoInv: {RateHtoInv}
            Initial RateDMtoSM: {RateDMtoSM}
            """)
            f.write("\n")
        
        #Save initial data for later dumping
        #Change here too
        self.quickdata = {
            "Ydm": Ydm,
            "Rate3to2": Rate3to2,
            "RateHtoInv": RateHtoInv,
            "RateDMtoSM": RateDMtoSM,
            "runtime": systime.time(),
            "slowroll": False,
            "inflationended": False,
            "xfreezeout": None,
            "inflationstart": None
        }
        
    def cleanup(self, time, data):
        """A function that is called after evolution finishes"""
        if self.fulloutput:
            self.datafile.close()
            self.extrafile.close()
            
        #Finish up the quickdata
        self.quickdata["runtime"] = systime.time() - self.quickdata["runtime"]
        self.quickdata["slowroll"] = self.slowroll
        if self.halt:
            self.quickdata["inflationended"] = True
            self.quickdata["xfreezeout"] = time
        self.quickdata["inflationstart"] = self.inflationstart

        #Relic Density

        Yinf = unpack(data)                         # Last value of Ydm
        mx = self.eomparams.model.mx               # GeV
        RELIC = 2.74372479e8                       # s0 / (rho_c/h^2)
        Omega_h2 = RELIC * mx * Yinf

        # Guárdalo en quickdata y muéstralo en consola
        self.quickdata["Yinf"] = Yinf
        self.quickdata["Omega_h2"] = Omega_h2
        print(f"Y∞ (final Ydm) = {Yinf:.12e}")
        print(f"Ωχ h² = {Omega_h2:.6e}")
        
        with open(self.basefilename + ".quick", 'wb') as f:
            pickle.dump(self.quickdata, f)
            
        with open(self.basefilename + ".info", "a") as f:
            f.write(f"Runtime (seconds): {self.quickdata['runtime']}\n")
            
        
    def derivatives(self, time, data):
        """Core function used by the numerical integrator to compute time derivatives of the system's variables"""
        unpacked_data = unpack(data)
        #data contains all dynamical variables being solved by integrator
        
        #Change for your rates
        (dYdmdx,Rate3to2, RateXiXitoh2h2, RateXih2toXjXk, Rate3toh2, RateHtoInv, RateDMtoSM, Hubble) = eoms(unpacked_data, self.eomparams, time)
        
        #Check for Boltzmann number-changing evolution
        #Change here too
        if (RateHtoInv/Hubble > 1.0 or RateDMtoSM/Hubble > 1.0) and self.inflationstart is None:
            self.inflationstart = time
            self.slowroll = True
        #elif self.slowroll and Rate4to2/Hubble < 1.0 and RateHtoInv/Hubble < 1.0 and RateDMtoSM/Hubble < 1.0:
        #    self.halt = True
        #    self.haltmsg = "Number-changing evolution has ended"
            
        #recombine the derivatives into a flat arry so the ODE solver can integrate them
        return pack(dYdmdx)
        
        
    def compute_timestep(self, time, data):
        """Compute the desired time step at this point in the evolution"""
        factor1, factor2 = self.timestepinfo
        #later perhaps come back here and adjust timesteps, for now just take a fixed timestep
        timestep = 0.01
        return timestep
        
        
    def write_data(self, time, data):
        """
        A function that is called between timesteps (and also at the start/end of integration), when data is ready for writing.
        """
        if not self.fulloutput:
            return
            
        #write the raw data
        text = self.format_data(time, data)
        self.datafile.write(text)
        
        #compute derived data
        unpacked_data = unpack(data)
        Ydm = unpacked_data
        
        #Change for your rates and desired quantities
        (Rate3to2, RateXiXitoh2h2, RateXih2toXjXk, Rate3toh2, RateHtoInv, RateDMtoSM, neqDM, neqh1, neqh2, Hubble, sent, Corr) =  compute_all(unpacked_data, self.eomparams, time)

                # --- Detectar primer cruce Rate/H < 1 ---
        if Rate3to2/Hubble < 1 and not hasattr(self, "x_cross_3to2"):
            self.x_cross_3to2 = time
            print(f"Rate3to2 cross at x = {time:.3f}")

        if RateXiXitoh2h2/Hubble < 1 and not hasattr(self, "x_cross_XiXitoh2h2"):
            self.x_cross_XiXitoh2h2 = time
            print(f"RateXiXitoh2h2 cross at x = {time:.3f}")

        if RateXih2toXjXk/Hubble < 1 and not hasattr(self, "x_cross_Xih2toXjXk"):
            self.x_cross_Xih2toXjXk = time
            print(f"RateXih2toXjXk cross at x = {time:.3f}")

        if Rate3toh2/Hubble < 1 and not hasattr(self, "x_cross_3toh2"):
            self.x_cross_3toh2 = time
            print(f"Rate3toh2 cross at x = {time:.3f}")

        if RateHtoInv/Hubble < 1 and not hasattr(self, "x_cross_HtoInv"):
            self.x_cross_HtoInv = time
            print(f"RateHtoInv cross at x = {time:.3f}")

        if RateDMtoSM/Hubble < 1 and not hasattr(self, "x_cross_DMtoSM"):
            self.x_cross_DMtoSM = time
            print(f"RateDMtoSM cross at x = {time:.3f}")

        
        #Change here too
        #In accordance with utilities.py (writing to and loading from .dat2)
        extradata = [Rate3to2, RateXiXitoh2h2, RateXih2toXjXk, Rate3toh2, RateHtoInv, RateDMtoSM, Hubble, sent, neqh1, neqh2, neqDM]
        
        #write the derived data
        text = self.format_data(time, extradata)
        self.extrafile.write(text)
