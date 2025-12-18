import numpy as np 
import pandas as pd
import seaborn as sns
from pymatgen.io.cif import CifParser
import matplotlib.pyplot as plt

kB = 1.380649e-23
NA = 6.02214076e23

plt.rcParams.update({'text.usetex':True, 
'font.family':'serif', 
'font.size':18, 
})

colors = sns.color_palette("mako")

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

exp_data = pd.read_csv('CO2_IRMOF-1_195K_exp_ads.csv', skiprows=11)
exp_data['pressure'] = exp_data['pressure']*0.1333223684
exp_data['total_adsorption'] = exp_data['total_adsorption']/MM

pressure = np.array([  1000.0000,   1599.8684,   5332.8947,   9865.8553,  13972.1842,
         15492.0592, 16000.,  16705.2928,  18078.5132,  19331.7434,  19998.3553,
         20664.9671,  21998.1908,  23331.4145,  25331.2500,  28664.3092,
         33330.5921,  40130.0329])

Nads = np.zeros_like(pressure)  
L0 = np.zeros_like(pressure) 
Lz = np.zeros_like(pressure)  
strain = np.zeros_like(pressure)  
volume_mean = np.zeros_like(pressure)  
volume_var = np.zeros_like(pressure)  

framework_mass = pore(structure)
V0 = 2.58312011**3

points = len(pressure)

for i in range(points):

    df = pd.read_csv(f'pressure_{1e-3*pressure[i]:.2f}_kPa/output.csv')

    # plt.plot(df.index.values,df.molecules.values,'-')
    # plt.axhline(y=df.molecules.mean(), linestyle="--", linewidth=1.5)
    # plt.xlabel("Steps")
    # plt.ylabel("Molecules")

    plt.plot(df.index.values,df.molecules.values,'-')
    plt.axhline(y=df.Lz.mean(), linestyle="--", linewidth=1.5)
    plt.xlabel("Steps")
    plt.ylabel("Molecules")

    Nads[i] = df.molecules.values.mean() 
    L0[i] = df.Lz.values[0] 
    Lz[i] = df.Lz.values.mean() 
    volume_mean[i] = 1e-27*df.volumes_mean.values.mean() 
    volume_var[i] = 1e-54*df.volumes_var.values.mean()  
    strain[i] = ((df.Lz.values**3-V0)/V0).mean() 

plt.show()

plt.plot(1e-3*pressure[:points], Nads[:points]*(1e3/framework_mass),'-o', color=colors[2], markersize=9, label='GCMC/MD')
# plt.plot(1e-3*gcmc['Pressure (Pa)'], gcmc['Absolute adsorption (mol/kg)'],'-o',markersize=8,markeredgewidth=1.5,mfc='none',color=colors[4], label='GCMC-rigid')
plt.plot(exp_data['pressure'][0:16], exp_data['total_adsorption'][0:16],'o', color=colors[2], markersize=8,markeredgewidth=1.5,mfc='none',label='Experimental')
# plt.plot(des['pressure'], des['total_adsorption'],'o', color=colors[2], markersize=8,label='Experimental')
plt.xlabel('Pressure (kPa)', fontsize=18)
plt.ylabel('Adsorption (mol/kg)',fontsize=18)
plt.minorticks_on()
# plt.xscale('log')
plt.tick_params(direction='in',right=True, top=True)
plt.tick_params(labelsize=18)
plt.tick_params(labelbottom=True, labeltop=False, labelright=False, labelleft=True)
plt.tick_params(direction='in',which='minor', length=4, bottom=True, top=True, left=True, right=True)
plt.tick_params(direction='in',which='major', length=8, bottom=True, top=True, left=True, right=True)
plt.legend(fontsize=16, frameon=False, edgecolor='k')
# plt.text(0.05, 0.95, r'CO$_{2}$, 195 K', transform=plt.gca().transAxes, 
#             fontsize=18, verticalalignment='top')
# plt.savefig(f'figures/adsorption_isotherm.svg', bbox_inches='tight')
plt.show()

plt.plot(1e-3*pressure[:points], 1e2*strain[:points],'-o', color=colors[2], markersize=9, label='GCMC/MD')
plt.xlabel('Pressure (kPa)', fontsize=18)
plt.ylabel('Strain (\%)',fontsize=18)
plt.minorticks_on()
# plt.xscale('log')
plt.tick_params(direction='in',right=True, top=True)
plt.tick_params(labelsize=18)
plt.tick_params(labelbottom=True, labeltop=False, labelright=False, labelleft=True)
plt.tick_params(direction='in',which='minor', length=4, bottom=True, top=True, left=True, right=True)
plt.tick_params(direction='in',which='major', length=8, bottom=True, top=True, left=True, right=True)
plt.legend(fontsize=16, frameon=False, edgecolor='k')
plt.text(0.05, 0.95, r'CO$_{2}$, 195 K', transform=plt.gca().transAxes, 
            fontsize=18, verticalalignment='top')
# plt.savefig(f'figures/strain_isotherm.svg', bbox_inches='tight')

plt.show()

plt.plot(1e-3*pressure[:points], 1e9*volume_var[:points]/(kB*T*volume_mean[:points]),'-o', color=colors[2], markersize=9, label='GCMC/MD')
plt.xlabel('Pressure (kPa)', fontsize=18)
plt.ylabel('Compressibility (GPa)',fontsize=18)
plt.minorticks_on()
# plt.xscale('log')
plt.tick_params(direction='in',right=True, top=True)
plt.tick_params(labelsize=18)
plt.tick_params(labelbottom=True, labeltop=False, labelright=False, labelleft=True)
plt.tick_params(direction='in',which='minor', length=4, bottom=True, top=True, left=True, right=True)
plt.tick_params(direction='in',which='major', length=8, bottom=True, top=True, left=True, right=True)
plt.legend(fontsize=16, frameon=False, edgecolor='k')
# plt.text(0.05, 0.95, r'CO$_{2}$, 195 K', transform=plt.gca().transAxes, 
#             fontsize=18, verticalalignment='top')
# plt.savefig(f'figures/adsorption_isotherm.svg', bbox_inches='tight')
plt.show()

compressibility = volume_var[:points]/(kB*T*volume_mean[:points])  
# strain = (Lz_mean[:points]**3-V0)/V0
adsorption_stress = pressure[:points]+strain/compressibility

plt.plot(1e-3*pressure[:points], 1e-6*adsorption_stress,'-o', color=colors[2], markersize=9, label='GCMC/MD')
plt.xlabel('Pressure (kPa)', fontsize=18)
plt.ylabel('Adsoprtion stress',fontsize=18)
plt.minorticks_on()
# plt.xscale('log')
plt.tick_params(direction='in',right=True, top=True)
plt.tick_params(labelsize=18)
plt.tick_params(labelbottom=True, labeltop=False, labelright=False, labelleft=True)
plt.tick_params(direction='in',which='minor', length=4, bottom=True, top=True, left=True, right=True)
plt.tick_params(direction='in',which='major', length=8, bottom=True, top=True, left=True, right=True)
plt.legend(fontsize=16, frameon=False, edgecolor='k')
# plt.text(0.05, 0.95, r'CO$_{2}$, 195 K', transform=plt.gca().transAxes, 
#             fontsize=18, verticalalignment='top')
# plt.savefig(f'figures/adsorption_isotherm.svg', bbox_inches='tight')
plt.show()

df = pd.DataFrame()
df['Pressure (Pa)'] = pressure[:points]
df['Absolute adsorption (molecues/uc)'] = Nads[:points]
df['Absolute adsorption (mol/kg)'] = Nads[:points]*(1e3/framework_mass)
df['Unit cell (nm)'] = Lz
df['Strain (%)'] = 1e2*strain
df['Compressibility (GPa)'] = compressibility*1e9 
df['Adsorption stress (GPa)'] = adsorption_stress*1e-9 
df.to_csv(f'CO2_{framework_name}_{T:.2f}K_GCMC-MD.csv')
