import subprocess
import sys
import numpy as np
from scipy.optimize import curve_fit

#### Initialize ####

experiment = sys.argv[1]
mode = sys.argv[2]
start_itt = int(sys.argv[3]) if len(sys.argv) > 3 else 0

if experiment not in ('incoherent', 'powder', 'magnon', 'single_crystal'):
    print('Error: Invalid experiment. Choose: incoherent, powder, magnon or single_crystal')
    sys.exit(1)
if mode not in ('high_res', 'wide_dyn_range'):
    print('Error: Invalid mode. Choose: high_res or wide_dyn_range')
    sys.exit(1)

from KVASIR import *
if mode == 'high_res':
    from ToF_conversion_high_res import *
if mode == 'wide_dyn_range':
    from ToF_conversion_wide_dyn_range import *


## Redction functions
def lin(x,a,b):
    return a*x + b

def FWHM(data_x, data_y):
    if np.sum(data_y) == 0:
        data_y = np.ones(len(data_y))*1e-5
        mean = np.sum(data_x*data_y)/np.sum(data_y)
        std = 0
        FWHM = 0
    else:
        mean = np.sum(data_x*data_y)/np.sum(data_y)
        std = np.sqrt( np.sum(data_y*(data_x - mean)**2)/np.sum(data_y) )
        FWHM = 2.355*std
    return FWHM, mean, std

def reduce_raw(pos, exp, mode):
    data = np.zeros((330,500,80))
    for j in range(330):
        for i in range(24):
            if j+i-12 < 0 or j+i-12 >=330:
                continue
            data[j+i-12]+=(np.genfromtxt(f'Data_raw/{exp}_{mode}_n{j}/Detector_{i}_ToF_{pos}.dat')[:500,:])

    np.savetxt(f"Data_reduced/data_{exp}_{mode}_{pos}.txt", np.vstack(data))
    print(f"{exp}_{pos} reduced")

def prismatic_Ef_file(pos, mode):
    if mode == 'high_res':
        E_scale = {'upper':np.linspace(1.84,1.96,500),'cent':np.linspace(1.93,2.1,500),'lower':np.linspace(1.84,1.96,500)}
    if mode == 'wide_dyn_range':
        E_scale = {'upper':np.linspace(7.4,7.9,500),'cent':np.linspace(7.7,8.4,500),'lower':np.linspace(7.4,7.9,500)}
    h = np.linspace(-0.2,0.2,80)

    data = np.genfromtxt(f'Data_raw/signal_on_det_{mode}/Detector_{pos}.dat')[:500,:]

    Ef_model = []
    Ef_model_err = []
    for j in range(80):
        _, Ef_mean, Ef_err = FWHM(E_scale[pos], data[:,j])
        Ef_model.append(Ef_mean)
        Ef_model_err.append(Ef_err)
    params = [(Ef_model[0]-Ef_model[-1])/len(h), Ef_model[40]]
    par, cov = curve_fit(lin, h, Ef_model, sigma=Ef_model_err, p0=params)

    np.savetxt(f'Storage/Ef_model_{pos}_{mode}.txt',lin(h,*par))

#### Pipeline ####
## Run Simulations ##
if start_itt == 0:
    subprocess.run(['rm', '-r', f'Data_raw/{experiment}_{mode}_n{start_itt}'], check=True)
    subprocess.run(['mkdir', 'Data_raw'], check=True)

if subprocess.run(['test', '-e', f'Storage/Ef_model_cent_{mode}.txt']).returncode != 0:
    print(f"RUNNING: generating prismatic Ef files in {mode} mode")
    KVASIR_instrument('gen_prismatic_Ef', mode)
    for p in ['upper', 'cent','lower']:
        prismatic_Ef_file(p, mode)

print(f"RUNNING: {experiment} in {mode} mode")
KVASIR_instrument(experiment, mode, start_itt)
for p in ['upper', 'cent','lower']:
    reduce_raw(p, experiment, mode)
subprocess.run(['rm', '-r', 'Data_raw'], check=True)

## ToF conversion ##
for p in ['upper', 'cent','lower']:
    o = conv(p, experiment, mode)
    o.dataset()

## plotting is done in plotexperiment.ipynb

