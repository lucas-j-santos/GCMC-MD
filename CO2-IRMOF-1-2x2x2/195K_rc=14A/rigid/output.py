import torch
import numpy as np 
import pandas as pd
import re
from pymatgen.io.cif import CifParser
import matplotlib.pyplot as plt

torch.set_default_dtype(torch.float64)

kB = 1.380649e-23
NA = 6.02214076e23

plt.rcParams.update({'text.usetex':True, 
'font.family':'serif', 
'font.size':18, 
'axes.linewidth':1.0, 
'lines.linewidth':1.5,
'legend.fontsize': 18,
'legend.frameon':False,
# 'figure.figsize':(6.4, 4.8)
})

class Framework:

    def __init__(self, cif_path):  

        # Define all geometric information of the framework
        self.cif_path = cif_path
        self.cif_parser = CifParser(self.cif_path, occupancy_tolerance=100)
        self.structure_adsorbent = self.cif_parser.parse_structures(primitive=False)[0]
        self.volume_cell = self.structure_adsorbent.volume
        self.get_lattice_and_coordinates()
        self.get_system_size()
        self.angles = np.array([self.alpha_deg, self.beta_deg, self.gamma_deg])
        self.num_atoms = self.coordinates.shape[0]

    def get_lattice_and_coordinates(self):
        # Get the coordinates of each atom
        # dic = self.cif_parser.as_dict()
        # dic = dic[list(dic.keys())[0]]
        self.a = self.structure_adsorbent.lattice.a
        self.b = self.structure_adsorbent.lattice.b
        self.c = self.structure_adsorbent.lattice.c
        self.alpha_rad = self.structure_adsorbent.lattice.angles[0]*np.pi/180.
        self.beta_rad = self.structure_adsorbent.lattice.angles[1]*np.pi/180.
        self.gamma_rad = self.structure_adsorbent.lattice.angles[2]*np.pi/180.
        self.alpha_deg = self.structure_adsorbent.lattice.angles[0]
        self.beta_deg = self.structure_adsorbent.lattice.angles[1]
        self.gamma_deg = self.structure_adsorbent.lattice.angles[2]
        self.lattice = self.get_lattice_RASPA()
        self.frac_coordinates = (self.structure_adsorbent.frac_coords)
        self.coordinates = (self.structure_adsorbent.frac_coords)@self.lattice 
        
    def get_lattice_RASPA(self):
        # This function calculates the lattice like RASPA because it is 
        # different from the one in pymatgen
        zeta = (np.cos(self.alpha_rad)-np.cos(self.gamma_rad)*np.cos(self.beta_rad))/np.sin(self.gamma_rad)
        lat_test = np.array([
            [self.a, self.b*np.cos(self.gamma_rad), self.c*np.cos(self.beta_rad)],
            [0, self.b*np.sin(self.gamma_rad), self.c*zeta],
            [0, 0, self.c*np.sqrt(1 - np.cos(self.beta_rad)**2 - zeta**2)]])
        return lat_test.T

    def get_system_size(self):
        # Get the size of the cell.
        self.system_size = [l for l in self.structure_adsorbent.lattice.lengths]
        self.system_size = np.array(self.system_size)

def pore(structure):

    forcefield = pd.DataFrame()
    forcefield['type'] = ['C','H','O','Zn']
    forcefield['mass'] = np.array([12.000, 1.0000, 15.999, 65.37]) 
    
    mass = 0.0

    for k, site in enumerate(structure.structure_adsorbent):

        atom_idx = forcefield["type"]==site.species_string
        mass += float(forcefield["mass"][atom_idx].iloc[0])

    return mass

MM = 15.999*2+12.0
framework_name = 'IRMOF-1'
structure = Framework(f'{framework_name}.cif')
T = 195.0
exp_data = pd.read_csv('../CO2_IRMOF-1_195K_exp_ads.csv', skiprows=11)
exp_data['pressure'] = exp_data['pressure']*133.3223684211
pressure = torch.tensor([  1000.0000,   1599.8684,   5332.8947,   9865.8553, 
                         11920.0, 12946., 13459.0921, 13972.1842, 15492.0592, 16705.2928,  
                         18078.5132,  19331.7434,  19998.3553, 20664.9671,  21998.1908,  
                         23331.4145,  25331.2500,  28664.3092, 33330.5921,  40130.0329])
framework_mass = pore(structure)
Nads = np.empty_like(pressure.numpy())

for i in range(len(pressure)):

    with open(f'pressure_{1e-3*pressure[i].numpy():.2f}_kPa/output.data', 'r') as file:
        content = file.read()

    production_pattern = re.compile(r'PRODUCTION Cycle:\s+([-+]?[0-9]*\.?[0-9]+),\s+([-+]?[0-9]*\.?[0-9]+)')
    matches = list(production_pattern.finditer(content))

    step = np.zeros(len(matches))
    molecules = np.zeros(len(matches)) 
    k = 0
    for match in matches:
        step[k] = match.group(1) 
        molecules[k] = match.group(2)
        k = k+1

    loading_pattern = re.compile(r'# MOLECULES')
    average_pattern = re.compile(r'Overall: Average:\s+([-+]?[0-9]*\.?[0-9]+)')
    match = loading_pattern.search(content)
    matches = list(average_pattern.finditer(content, match.end()))
    Nads[i] = float(matches[1].group(1)) 

    plt.plot(step,molecules)
    plt.axhline(y=np.mean(molecules), color="k", linestyle="--", linewidth=1.5)
    plt.xlabel("Steps")
    plt.ylabel("Molecules")

plt.show()

exp_data['total_adsorption'] = exp_data['total_adsorption']/MM

with plt.style.context('seaborn-v0_8'):
    plt.plot(1e-3*pressure, Nads*(1e3/framework_mass),'C0-o', markersize=7, markeredgewidth=1.5, mfc='none')
    # plt.plot(exp_data['pressure'], exp_data['total_adsorption'],'C0o', markersize=7)
    plt.xlabel('Pressure (kPa)')
    plt.ylabel('Absolute adsorption (mol/kg)')
plt.show()

df = pd.DataFrame()
df['Pressure (Pa)'] = pressure.numpy()
df['Absolute adsorption (molecules/uc)'] = Nads/8
df['Absolute adsorption (mol/kg)'] = Nads*(1e3/framework_mass)
df.to_csv(f'CO2_{framework_name}_{T:.2f}K.csv')
