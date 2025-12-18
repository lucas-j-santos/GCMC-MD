import numpy as np
import pandas as pd
import os
import subprocess
import mdtraj
import MDAnalysis as mda
from utils import *

temperature = np.arange(50.0,650.0,50.)
# temperature = np.array([87.0,195.0])
# temperature = np.array([169.0, 258.0])
pressure = 1e5
Lz = np.empty_like(temperature)
rc = 14.0
current_dir = os.getcwd()

for k in range(len(temperature)):

    new_dirn = f'temperature_{temperature[k]:.0f}_K'
    os.mkdir(new_dirn)
    os.chdir(new_dirn)

    write_mdp(temperature[k], pressure, rc, 1000000, 2.0, 'yes')
    
    subprocess.run('gmx grompp -f npt.mdp -c ../box.gro -p ../topol.top -o npt -maxwarn 2', 
                shell=True, check=True, text=True)
    subprocess.run('gmx mdrun -v -deffnm npt', shell=True, check=True, text=True)
    subprocess.run('rm mdout.mdp npt.log npt.tpr npt.gro npt.edr', shell=True, check=True, text=True)

    traj = mdtraj.load('npt.trr', top='../box.gro')
    traj[-1].save_gro('npt.gro')
    
    write_mdp(temperature[k], pressure, rc, 2000000, 1.0, 'no')
    subprocess.run('gmx grompp -f npt.mdp -c npt.gro -t npt.cpt -p ../topol.top -o npt -maxwarn 2', 
                shell=True, check=True, text=True)
    subprocess.run('gmx mdrun -v -deffnm npt', shell=True, check=True, text=True)
    
    aux = mda.auxiliary.EDR.EDRReader('npt.edr')
    data_dictionary = aux.data_dict
    df_npt = pd.DataFrame(data_dictionary)
    df_npt.to_csv(f'npt_output.csv', index=False)

    subprocess.run('rm *# mdout.mdp npt.cpt npt.log npt.tpr npt.gro npt.edr npt_prev.cpt', shell=True, check=True, text=True)
    traj = mdtraj.load('npt.trr', top='../box.gro')
    traj[-1].save_gro('npt.gro')
    Lz[k] = traj.unitcell_lengths.mean()
    os.remove('npt.trr')
    
    os.chdir(current_dir)
    
df = pd.DataFrame()
df['Temperature (K)'] = temperature
df['Unit-cell length (nm)'] = Lz
df.to_csv('IRMOF-1_MD.csv')
