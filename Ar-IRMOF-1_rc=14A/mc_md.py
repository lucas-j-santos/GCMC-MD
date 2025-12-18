import torch
import numpy as np
import pandas as pd
import os
import re
import subprocess
import mdtraj
import MDAnalysis as mda
from cdft.lj_eos import lj_eos
from utils import *

framework_name = 'IRMOF-1'
temperature = 78.0
# pressure = torch.hstack((torch.tensor([5.0]), torch.arange(1e1, 6e1, 1e1), torch.tensor([60.,65.0,66.]), torch.arange(7e1,2.1e2,1e1))) 
pressure = np.hstack((np.array([5.0]), np.arange(1e1, 6e1, 1e1), np.array([60.,70.0,80.,84.,85.]), np.arange(9e1,2.1e2,1e1)))
pressure = torch.tensor(pressure)

sigma = 3.40
epsilon = 119.8
parameters = {'sigma':sigma, 'epsilon':epsilon}
eos = lj_eos(parameters, temperature)
bulk_density = torch.empty_like(pressure)
fugacity_coefficient = torch.empty(len(bulk_density))
bulk_density[0] = eos.density(pressure[0],'vap')
fugacity_coefficient[0] = eos.fugacity_coefficient(bulk_density[0])
for i in range(1,len(pressure)):
    bulk_density[i] = eos.density(pressure[i],bulk_density[i-1])
    fugacity_coefficient[i] = eos.fugacity_coefficient(bulk_density[i])

rc = 14.0
molecule_name = 'argon'
iterations = 100
mol = 377

raspa_dir = '/home/lucas/Programs/gRASPA'

Nads = np.empty(iterations)
unitcell_lengths = np.empty((iterations,3))
volumes_mean = np.empty(iterations)
volumes_var = np.empty(iterations)

current_dir = os.getcwd()

for k in range(9,10):

    new_dirn = f'pressure_{pressure[k].numpy():.2f}_Pa'
    os.mkdir(new_dirn)
    os.chdir(new_dirn)
    subprocess.run(f'cp ../force_field_mixing_rules.def ../force_field.def ../{molecule_name}.def ../pseudo_atoms.def .', 
                   shell=True, check=True, text=True)

    write_mdp(temperature, pressure[k], rc)
    
    subprocess.run('gmx grompp -f npt.mdp -c ../box.gro -p ../topol.top -o npt.tpr -maxwarn 2', 
                shell=True, check=True, text=True)
    subprocess.run('gmx mdrun -v -deffnm npt', shell=True, check=True, text=True)
    subprocess.run('rm mdout.mdp npt.cpt npt.log npt.tpr npt.gro npt.edr', shell=True, check=True, text=True)

    traj = mdtraj.load('npt.trr', top='../box.gro')
    traj[-1].save_gro('npt.gro')
    os.remove('npt.trr')

    for i in range(iterations):

        write_inp(framework_name,temperature,pressure[k],rc,molecule_name,fugacity_coefficient[k],mol)

        if i == 0:
            framework_mc = mdtraj.load('npt.gro', top='npt.gro')
            write_cif(framework_mc, f'{framework_name}.cif')    
            unitcell_lengths[i] = framework_mc.unitcell_lengths
            volumes_mean[i] = traj.unitcell_volumes.mean()
            volumes_var[i] = traj.unitcell_volumes.var()

            framework_mc.save_gro('npt_mod.gro')
            os.remove('npt.gro')

            del traj, framework_mc

        else:
            traj = mdtraj.load('npt.trr',top='minim.gro')
            traj[-1].save_gro('npt.gro')
            framework_md = mdtraj.load('npt.gro', top='npt.gro')
            argon_residues = [res for res in framework_md.topology.residues if any(atom.name == 'Ar' and atom.element.symbol == 'Ar' for atom in res.atoms)]
            argon_atoms = set(atom.index for res in argon_residues for atom in res.atoms)
            keep_atoms = [atom.index for atom in framework_md.topology.atoms if atom.index not in argon_atoms]

            framework_mc = framework_md.atom_slice(keep_atoms)
            framework_mc.save_gro('npt_mod.gro')

            subprocess.run('rm minim.gro npt.gro npt.trr', shell=True, check=True, text=True)

            write_cif(framework_mc, f'{framework_name}.cif')
            
            unitcell_lengths[i] = framework_mc.unitcell_lengths
            volumes_mean[i] = traj.unitcell_volumes.mean()
            volumes_var[i] = traj.unitcell_volumes.var()

            del traj, framework_md, framework_mc
        
        subprocess.run(f'{raspa_dir}/src_clean/nvc_main.x', shell=True, check=True, text=True)
        subprocess.run(f'mv Output/*.data gcmc_output_{i}.data', shell=True, check=True, text=True)
        subprocess.run('rm -r AllData FirstBead Lambda Movies Output Restart TMMC', shell=True, check=True, text=True)

        with open(f'gcmc_output_{i}.data', 'r') as file:
            content = file.read()

        loading_pattern = re.compile(r'# MOLECULES')
        average_pattern = re.compile(r'Overall: Average:\s+([-+]?[0-9]*\.?[0-9]+)')
        match = loading_pattern.search(content)
        matches = list(average_pattern.finditer(content, match.end()))
        Nads[i] = float(matches[1].group(1)) 
        mol = int(np.round(Nads[i]))
        write_topol(framework_name,molecule_name,mol)
        
        if i < iterations-1:
        
            subprocess.run(f'gmx insert-molecules -f npt_mod.gro -ci ../{molecule_name}.gro -nmol {mol} -o box.gro', shell=True, check=True, text=True)
            subprocess.run(f'gmx grompp -f ../minim.mdp -c box.gro -p topol.top -o minim.tpr -maxwarn 2', shell=True, check=True, text=True)
            subprocess.run(f'gmx mdrun -v -deffnm minim', shell=True, check=True, text=True)
            subprocess.run(f'gmx grompp -f npt.mdp -c minim.gro -p topol.top -o npt.tpr -maxwarn 2', shell=True, check=True, text=True)
            subprocess.run(f'gmx mdrun -v -deffnm npt', shell=True, check=True, text=True)
            
            aux = mda.auxiliary.EDR.EDRReader('npt.edr')
            data_dictionary = aux.data_dict
            df_npt = pd.DataFrame(data_dictionary)
            df_npt.to_csv(f'npt_output_{i}.csv', index=False)
            
            subprocess.run(f'rm minim.edr box.gro minim.log minim.tpr minim.trr mdout.mdp npt.cpt npt.edr npt.log npt.tpr npt.gro npt_mod.gro', shell=True, check=True, text=True)
    
    df = pd.DataFrame()
    df['molecules'] = Nads
    df['Lx'] = unitcell_lengths[:,0]
    df['Ly'] = unitcell_lengths[:,1]
    df['Lz'] = unitcell_lengths[:,2]
    df['volumes_mean'] = volumes_mean
    df['volumes_var'] = volumes_var
    df.to_csv(f'output.csv')
    os.chdir(current_dir)
