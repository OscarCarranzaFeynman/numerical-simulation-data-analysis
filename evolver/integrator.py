"""
integrator.py

Integration library built on top of scipy.integrate.ode and scipy.integrate.solve_ivp

This library is agnostic as to the model, which must be implemented as a subclass of AbstractModel
"""
from enum import Enum
import pickle
from scipy.integrate import ode
from scipy.integrate import solve_ivp


class IntegrationError(Exception):
    """Raised when an integration error arises"""
    
class Status(Enum):
    """Status of Driver object"""
    Initializing = 0
    OK = 1
    Finished = 2
    IntegrationError = -1
    Terminated = -2
    
class AbstractModel(object):
    """Abstract class that handles all model and run-dependent information"""
    
    def __init__(self, run_params):
        """Sets up all data for an evolution
        
        Arguments:
        *run_params: Dictionary of all initialization information
        
        run_params can include whatever keys you like. The following keys have meaning for AbstractModel:
            * initial_data: The initial data for the integration (required)
            * start_time: The initial value for x \equiv mDM/T (default 1)
            * end_time: The maximum value of x to which to integrate (default 100)
            * atol: Absolute tolerance in integration (default 1e-10
            * rtol: Relative tolerance in integration (default 1e-10)
            * separator: Separator to use between fields when outputting data (default ", ")
        """
        #Defaults (initialize the attribute self.parameters as a dictionary with settings for integration run)
        self.parameters = {
            'start_time': 1.0,
            'end_time': 100.0,
            'atol': 1e-10,
            'rtol':1e-9,
            'separator': ', '
        }
        
        #Add in new information/update default information (update self.parameters dictionary with contents of run_params dictionary -- the one passed by the user when initializing the module)
        self.parameters.update(run_params)
        
        #store values that will be used repeatedly (initialize instance attributes of AbstractModel class using values from the self.parameters dictionary)
        self.initial_data = self.parameters['initial_data']
        self.start_time = self.parameters['start_time']
        self.end_time = self.parameters['end_time']
        self.atol = self.parameters['atol']
        self.rtol = self.parameters['rtol']
        self.separator = self.parameters['separator']
        
        #Initialize internal flags
        self.halt = False
        self.haltmsg = None
        
    def begin(self, time):
        """A function (method) that is called before evolution begins"""
        raise NotImplementedError()
        
    def cleanup(self, time, data):
        """A function that is called after evolution finishes"""
        raise NotImplementedError()
        
    def derivatives(self, time, data):
        """Computes derivatives for evolution"""
        raise NotImplementedError()
        
    def compute_timestep(self, time, data):
        """Compute the desired timestep at this point in the evolution"""
        raise NotImplementedError()
        
    def solout(self, time, data):
        """A function that is called after evey internal timestep.
        If this returns -1, it terminates the evolution.
        """
        if self.halt:
            return -1
        return 0
        
    def write_data(self, time, data):
        """ A function that is called between timesteps (and also at the start/end of integration), when data is ready for writing
        """
        raise NotImplementedError()
        
    def format_data(self, time, data):
        """Format the current data for output
        This method returns a string-formatted version of the time and current state of the system contained in data to be written to an output file"""
        return (str(time) + self.separator + self.separator.join(map(str,data)) + "\n")
        
    def save(self, filename):
        """Save this model to file (using pickle)
        Allows you to save model setups and reload them later without redefining anything manually"""
        with open(filename, 'wb') as f:
            pickle.dump(self.parameters, f)
            
    @classmethod
    def load(cls, filename):
        """Instantiate a model from a previously-saved file"""
        with open(filename, 'rb') as f:
            data = pickle.load(f)
        return cls(data)
        
        
class Driver(object):
    """Sets up and runs the evolution"""
    
    def __init__(self, model):
        """Set parameters for driving the evolution
        
        Arguments:
            * model: an AbstractModel suclass instance
        """
        self.status = Status.Initializing
        self.error_msg = ""
        self.model = model
        
        #Initializing the integrator
        self.data = model.initial_data
        if not self.model.stiff_solver:
            self.integrator = ode(model.derivatives)
            self.integrator.set_integrator('dopri5',
                                            nsteps=10000,
                                            rtol=model.rtol,
                                            atol=model.atol)
            self.integrator.set_initial_value(model.initial_data, model.start_time)
            self.integrator.set_solout(model.solout)
        
        self.status = Status.OK
        
    @property
    def time(self):
        """Returns the current time that the integrator is at"""
        """Only works for the Old run code"""
        return self.integrator.t
        
    def run_stiff(self):
        """
        Runs the evolution for stiff system
        """
        #Make sure we're ready to roll
        if self.status != Status.OK:
            raise ValueError("Cannot begin evolution as status is not OK.")
        
        #Initialization
        self.model.begin(self.model.start_time)
        
        #write initial data
        self.model.write_data(self.model.start_time, self.data)
            
        #Build the timeline for evaluation
        eval_times=[self.model.start_time]
        act_time=self.model.start_time
        while act_time < self.model.end_time:
            act_time += self.model.compute_timestep(act_time, self.data)
            eval_times.append(act_time)
            
        try:
            #result = solve_ivp(self.model.derivatives, (self.model.start_time, self.model.end_time), self.data, method='BDF', t_eval=eval_times, 
            #                    atol=self.model.atol, rtol=self.model.rtol, vectorized=False, dense_output=False)
            result = solve_ivp(self.model.derivatives, (self.model.start_time, self.model.end_time), self.data, method='Radau', 
                                            atol=self.model.atol, rtol=self.model.rtol, vectorized=False, dense_output=True)
                                
            for i, t in enumerate(result.t):
                self.data = result.y[:,i]
                self.model.write_data(t, self.data)
                
            
            if not result.success:
                self.status = Status.IntegrationError
                self.error_msg = result.message
            
            elif result.status == 1:
                self.status = Status.Terminated
                
            else:
                self.status = Status.Finished
                
                    
            self.model.cleanup(result.t[-1], result.y[:,-1])
            
        except Exception as e:
            self.status = Status.IntegrationError
            self.error_msg = str(e)
        
        
    def run(self):
        """
        Runs the evolution
        
        The resulting status is stored in self.status.
        If it stores Status.IntegrationError, the error message is available through self.error_msg.
        """
        #Make sure we're ready to roll
        if self.status != Status.OK:
            raise ValueError("Cannot begin evolution as status is not OK.")
            
        #Initialization
        self.model.begin(self.time)
        
        #Write initial data
        self.model.write_data(self.time, self.data)
        
        #Integration Loop
        newtime = self.time
        while True:
            #Check to see if we're finished
            if newtime >= self.model.end_time:
                self.status = Status.Finished
                break
            elif self.model.halt:
                self.status = Status.Terminated
                self.error_msg = self.model.haltmsg
                break
                
            #Construct the time to integrate to
            newtime = self.time + self.model.compute_timestep(self.time, self.data)
            
            #Take a step to newtime
            try:
                results = self.integrator.integrate(newtime, relax=True)
                if not self.integrator.successful():
                    raise IntegrationError("DOPRI Error Code {}"
                                            .format(self.integrator.get_return_code()))
                self.data = results
            except IntegrationError as e:
                self.status = Status.IntegrationError
                self.error_msg = e.args[0]
                break
                
            #write the data
            self.model.write_data(self.time, self.data)
            
        #Clean up
        self.model.cleanup(self.time, self.data)
