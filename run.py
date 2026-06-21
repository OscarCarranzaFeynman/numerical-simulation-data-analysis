"""
run.py

Performs a single evolution
"""
import sys
import time
from math import sqrt
from evolver.integrator import Driver, Status
from evolver.initialize import create_package, create_parameters
from evolver.inflation import HiddenVector
from evolver.model import Model
from evolver.parameters import get_portal_params


#Initialize all settings
#Masses in GeV
mx      = 0.5   # masa del vector oscuro X (GeV)
mh2     = 0.55    # escalar oscuro (GeV)
gX      = 0.0043
beta    = 0.1       # ángulo de mezcla (rad)

mh1     = 125.0       # Higgs SM (GeV) 
vEW = 246           #Electroweak vev

# Compute parameters
portal_params = get_portal_params(
    mx=mx, mh2=mh2, gX=gX, beta=beta, mh1=mh1, vEW=vEW
)
print("Portal params:", portal_params)

me = 0.511e-3       #Electron mass
mmu = 105.7e-3      #Muon mass
mtau = 1.77         #Tauon mass
mup = 2e-3          #Up quark mass
mdown = 5e-3        #Down quark mass
mcharm = 1.27       #Charm quark mass
mstrange = 96.4e-3  #Strange quark mass
mtop = 173          #Top quark mass
mbottom = 4.18      #Bottom quark mass
filename = "data/output"
if mh2 > mx:
    XiXitoh2h2_mode = "reverse"
elif mh2 < mx:
    XiXitoh2h2_mode = "forward"

# Built the dictionary
hv_kwargs = {
    **portal_params,   # contiene: mx, mh2, gX, beta, mh1, vEW, vPhi, lambdaH, lambdaPhi, lambdaM, muPhi, ...
    "me": me, "mmu": mmu, "mtau": mtau,
    "mup": mup, "mdown": mdown, "mcharm": mcharm, "mstrange": mstrange,
    "mtop": mtop, "mbottom": mbottom, "XiXitoh2h2_mode": XiXitoh2h2_mode
}

#This is a subclass of DarkMatterModel, thing that we have are print(infmodel.mx) -> 10.0,   print(infmodel.info())  ->  descripción completa del modelo,  print(infmodel.SigmavXiXitoh2h2(Ydm=None))-> <σv> for a certain process

infmodel = HiddenVector(**hv_kwargs)


package = create_package(
    infmodel=infmodel,
    start_time=1.0, end_time=300.0,
    basefilename=filename,
    timestepinfo=[200, 10],
    perform_run=True, stiff_solver=True
)

parameters = create_parameters(package)
if parameters is None:
    print("Unable to construct initial conditions")
    sys.exit(1)
    
#Create the model
model = Model(parameters)
model.save(filename + ".params")

#Construct the driver
driver = Driver(model)

#Perform the evolution
print("Beginning evolution")
start = time.time()
if model.stiff_solver:
    driver.run_stiff()
else:
    driver.run()
end = time.time()
print("Finished in {} s".format(round(end - start, 4)))

#Check to see what our status is
if driver.status == Status.IntegrationError:
    print("Unable to integrate further: {}".format(driver.error_msg))
elif driver.status == Status.Terminated:
    print("Evolution completed with message: {}".format(driver.error_msg))
elif driver.status == Status.Finished:
    print("Evolution completed!")
