import torch
import numpy as np
import pandas as pd
import os
import re
import subprocess
import mdtraj
from cdft.pcsaft_eos import pcsaft

torch.set_default_dtype(torch.float64)

def write_inp(framework_name,temperature,pressure,rc,molecule_name,fugacity_coefficient,create_molecules):

    with open('simulation.input', 'w+') as f:

        f.write(f"""NumberOfInitializationCycles 15000     
NumberOfEquilibrationCycles  0
NumberOfProductionCycles     15000

UseMaxStep  no

UseChargesFromCIFFile no

RestartFile no
RandomSeed  0

NumberOfTrialPositions 10
NumberOfTrialOrientations 10

AdsorbateAllocateSpace 10240
NumberOfSimulations 1
SingleSimulation yes

InputFileType cif
FrameworkName {framework_name}
UnitCells 0 1 1 1

ChargeMethod Ewald
Temperature {temperature:.1f}
Pressure    {pressure:.1f}

OverlapCriteria 1e5
CutOffVDW {rc:.2f}
CutOffCoulomb {rc:.2f}
EwaldPrecision 1e-6

SaveOutputToFile yes

Component 0 MoleculeName             {molecule_name}
            IdealGasRosenbluthWeight 1.0
            FugacityCoefficient      {fugacity_coefficient:.4f}
            TranslationProbability   1.0
            RotationProbability      1.0
            ReinsertionProbability   1.0
            SwapProbability          1.0
            CreateNumberOfMolecules  {create_molecules}
""")
        
def write_cif(traj, cif_file):

    with open(cif_file, 'w') as f:

        box_lengths = traj.unitcell_lengths[0]
        a, b, c = box_lengths

        f.write("data_IRMOF-1\n")
        f.write("_symmetry_space_group_name_H-M    'P 1'\n")
        f.write("_symmetry_Int_Tables_number       1\n")
        f.write(f"_cell_length_a                    {a*10:.5f}\n")
        f.write(f"_cell_length_b                    {b*10:.5f}\n")
        f.write(f"_cell_length_c                    {c*10:.5f}\n")
        f.write("_cell_angle_alpha                 90.0\n")
        f.write("_cell_angle_beta                  90.0\n")
        f.write("_cell_angle_gamma                 90.0\n")
        f.write("\n")
        f.write("loop_\n")
        f.write("_atom_site_label\n")
        f.write("_atom_site_type_symbol\n")
        f.write("_atom_site_fract_x\n")
        f.write("_atom_site_fract_y\n")
        f.write("_atom_site_fract_z\n")

        for i, atom in enumerate(traj.topology.atoms):
            pos = traj.xyz[0, i]  # posição em nanômetros
            # Converte posição cartesiana em fração da célula unitária
            fx = pos[0] / a
            fy = pos[1] / b
            fz = pos[2] / c

            f.write(f"{atom.name} {atom.name} {fx:.5f} {fy:.5f} {fz:.5f}\n")

m = torch.tensor([1.5131], dtype=torch.float64)
sigma = torch.tensor([3.1869], dtype=torch.float64)
epsilon = torch.tensor([163.33], dtype=torch.float64)
q = torch.tensor([4.4], dtype=torch.float64) # DA
parameters = {'m':m, 'sigma':sigma, 'epsilon':epsilon, 'q':q}


framework_name = 'IRMOF-1'
temperature = 195.0

pressure = torch.tensor([  1000.0000,   1599.8684,   5332.8947,   9865.8553, 
                         11920.0, 12946., 13459.0921, 13972.1842, 15492.0592, 16705.2928,  
                         18078.5132,  19331.7434,  19998.3553, 20664.9671,  21998.1908,  
                         23331.4145,  25331.2500,  28664.3092, 33330.5921,  40130.0329])

bulk_density = torch.empty_like(pressure)
fugacity_coefficient = torch.empty_like(pressure)
composition = torch.tensor([1.0])

eos = pcsaft(parameters, temperature)
bulk_density[0] = eos.density(pressure[0],composition,'vap')
fugacity_coefficient[0] = eos.fugacity_coefficient(bulk_density[0],composition)
for i in range(1,len(pressure)):
    bulk_density[i] = eos.density(pressure[i],composition,bulk_density[i-1])
    fugacity_coefficient[i] = eos.fugacity_coefficient(bulk_density[i],composition)

rc = 14.0
molecule_name = 'CO2'

framework_mc = mdtraj.load('../box.gro', top='../box.gro')
write_cif(framework_mc, f'{framework_name}.cif')
Nads = np.zeros_like(pressure.numpy())
mol = 1515

current_dir = os.getcwd()

for k in range(12,len(pressure)):

    new_dirn = f'pressure_{1e-3*pressure[k].numpy():.2f}_kPa'
    os.mkdir(new_dirn)
    os.chdir(new_dirn)

    subprocess.run(f'cp ../force_field_mixing_rules.def ../force_field.def ../{framework_name}.cif ../{molecule_name}.def ../pseudo_atoms.def .', 
                   shell=True, check=True, text=True)
    
    write_inp(framework_name,temperature,pressure[k],rc,molecule_name,fugacity_coefficient[k],mol)
    subprocess.run('/home/lucas/Programs/gRASPA/src_clean/nvc_main.x', shell=True, check=True, text=True)
    subprocess.run('mv Output/*.data output.data', shell=True, check=True, text=True)
    subprocess.run('rm -r AllData FirstBead Lambda Movies Output Restart TMMC', shell=True, check=True, text=True)

    with open(f'output.data', 'r') as file:
        content = file.read()

    loading_pattern = re.compile(r'# MOLECULES')
    average_pattern = re.compile(r'Overall: Average:\s+([-+]?[0-9]*\.?[0-9]+)')
    match = loading_pattern.search(content)
    matches = list(average_pattern.finditer(content, match.end()))
    Nads[k] = float(matches[1].group(1)) 
    mol = int(np.round(Nads[k]))

    os.chdir(current_dir)

df = pd.DataFrame()
df['Pressure (Pa)'] = pressure.numpy()
df['Absolute adsorption (molecules/uc)'] = Nads
# df.to_csv(f'{molecule_name}_{framework_name}_{temperature:.0f}K.csv')
