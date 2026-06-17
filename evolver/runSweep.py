"""
runSweep.py

Performs a sweep of evolutions

use python3 -m evolver.runSweep to run this code
"""
import multiprocessing as mp
import numpy as np
import pickle
#from typing import Dict, Iterable, Tuple, Any, List, Optional
import sys
import time
from math import sqrt, log10
from evolver.integrator import Driver, Status
from evolver.initialize import create_package, create_parameters
from evolver.inflation import HiddenVector
from evolver.model import Model
from evolver.parameters import get_portal_params
import signal
import os


def run_once(*, mx, mh2, gX, beta, basefilename="data/output"):
    # Defaults for fixed SM massses
    ################### GeV
    me = 0.511e-3       #Electron mass
    mmu = 105.7e-3      #Muon mass
    mtau = 1.77         #Tauon mass
    mup = 2e-3          #Up quark mass
    mdown = 5e-3        #Down quark mass
    mcharm = 1.27       #Charm quark mass
    mstrange = 96.4e-3  #Strange quark mass
    mtop = 173          #Top quark mass
    mbottom = 4.18      #Bottom quark mass
    mh1 = 125            #Higgs boson mass
    vEW = 246           #Electroweak vev
    #####################
    filename = basefilename
    if mh2 > mx:
        XiXitoh2h2_mode = "reverse"
    elif mh2 < mx:
        XiXitoh2h2_mode = "forward"
    
    # Portal-derived params (vPhi, lambdas, muPhi, etc.)
    portal_params = get_portal_params(mx=mx,mh2=mh2, gX=gX, beta=beta, mh1=mh1, vEW=vEW)
   
    hv_kwargs = {
        **portal_params,
        "me": me, "mmu": mmu, "mtau": mtau,
        "mup": mup, "mdown": mdown, "mcharm": mcharm, "mstrange": mstrange,
        "mtop": mtop, "mbottom": mbottom,
        "XiXitoh2h2_mode": XiXitoh2h2_mode,
    }

    infmodel = HiddenVector(**hv_kwargs)
    
    package = create_package(
        infmodel=infmodel,
        start_time=1.0,
        end_time=300.0,             
        basefilename=basefilename,
        timestepinfo=[200, 10],
        perform_run=True,
        stiff_solver=True,
        fulloutput=False,            
    )
    
    parameters = create_parameters(package)

    if parameters is None:
        return (None, None, 1)
    #Create the model
    model = Model(parameters)
    model.save(filename + ".params")

    #Construct the driver
    driver = Driver(model)

    #Perform the evolution
    #print("Beginning evolution")
    start = time.time()
    if model.stiff_solver:
        driver.run_stiff()
    else:
        driver.run()
    end = time.time()
    #print("Finished in {} s".format(round(end - start, 4)))

    #Check to see what our status is
    if driver.status == Status.IntegrationError:
        #print("Unable to integrate further: {}".format(driver.error_msg))
        Switch = 2
    elif driver.status == Status.Terminated:
        #print("Evolution completed with message: {}".format(driver.error_msg))
        Switch = 0
    elif driver.status == Status.Finished:
        #print("Evolution completed!")
        Switch = 0
    
    #File 3
    with open(filename + ".quick", 'rb') as f:
        quickdata = pickle.load(f)
    
    Ydm = quickdata.get("Yinf", quickdata.get("Ydm", None)) 
    Omega_h2 = quickdata.get("Omega_h2", None)

    return (Ydm, Omega_h2, Switch)

def worker_task(args):
    (job_id, mx, beta, r, gX, basefilename, timeout_s) = args

    # filename único para evitar conflictos
    run_base = f"{basefilename}_tmp_{job_id}"

    # tu r es mh2/mx  =>  mh2 = r*mx
    mh2 = r * mx

    def alarm_handler(signum, frame):
        raise TimeoutError()

    old = signal.signal(signal.SIGALRM, alarm_handler)
    signal.alarm(int(timeout_s))

    try:
        # correr y usar lo que regresa run_once
        Ydm, Omega_h2, Switch = run_once(
            mx=mx,
            mh2=mh2,
            gX=gX,
            beta=beta,
            basefilename=run_base,
        )

    except TimeoutError:
        Ydm, Omega_h2, Switch = None, None, 9  # timeout

    except Exception:
        Ydm, Omega_h2, Switch = None, None, 2  # otro error

    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)

        # borrar archivos temporales
        for ext in [".quick", ".params"]:
            try:
                os.remove(run_base + ext)
            except FileNotFoundError:
                pass

    relic_min = 1e-5
    relic_max = 1.0

    if (
        Switch != 0
        or Ydm is None
        or Omega_h2 is None
        or Omega_h2 <= relic_min
        or Omega_h2 >= relic_max
    ):
        Ydm_out = -999
        Omega_out = -999
    else:
        Ydm_out = Ydm
        Omega_out = Omega_h2

    separator = ", "
    text = (
        str(mx) + separator
        + str(beta) + separator
        + str(r) + separator
        + str(gX) + separator
        + str(Ydm_out) + separator
        + str(Omega_out) + separator
        + str(Switch) + "\n"
    )

    return text

mx_start   = 10.0**-3
mx_stop    = 10.0**1
mx_steps   = 25

beta_start = 10.0**-5
beta_stop  = np.pi
beta_steps = 15

r_start=1.01
r_stop=4
r_steps=15

gX_start   = 0.01
gX_stop    = 10
gX_steps   = 20

mxs   = np.logspace(log10(mx_start),   log10(mx_stop),   mx_steps)
rs  = np.linspace(r_start, r_stop,  r_steps)
betas = np.linspace(beta_start, beta_stop, beta_steps)
gXs   = np.logspace(log10(gX_start), log10(gX_stop), gX_steps)

print(mxs)
print(betas)
print(rs)
print(gXs)

basefilename = "data/SweepHiddenVector1"




if __name__ == "__main__":
    os.makedirs("data", exist_ok=True)
    datafile = open(basefilename + ".dat", "w")
    total_start_time = time.time()

    timeout_s = 30
    nproc = max(1, mp.cpu_count()-1)

    jobs = []
    counter = 0

    for mx in mxs:
        for beta in betas:
            for r in rs:
                for gX in gXs:
                    counter += 1
                    jobs.append((counter, mx, beta, r, gX, basefilename, timeout_s))

    print("Total runs:", len(jobs))
    print("Using processes:", nproc)

    pool = mp.Pool(nproc)

    done = 0

    for line in pool.imap_unordered(worker_task, jobs, chunksize=1):
        if line is None:
            continue
        datafile.write(line)

        done += 1
        print(done, flush=True)

        if done % 100 == 0:
            datafile.flush()

    pool.close()
    pool.join()

    datafile.close()
    
    total_end_time = time.time()
    total_seconds = total_end_time - total_start_time
    print("\n======================================")
    print(f"Total sweep time: {total_seconds:.2f} seconds")
    print(f"Total sweep time: {total_seconds/60:.2f} minutes")
    print("======================================\n")




