import numpy as np
import pandas as pd
import seaborn as sns
import re
import matplotlib.pyplot as plt

kB = 1.380649e-23  # in J/K
Na = 6.02214076e23  # Avogadro's number

plt.rcParams.update({'text.usetex':True, 
'font.family':'sans-serif', 
'font.size':18, 
'axes.linewidth':1.1,
'text.latex.preamble': r'\usepackage{sfmath}'
})

colors = sns.color_palette(palette='mako', n_colors=7)

framework_name = 'IRMOF-1'
T = 78.0

pressure = np.hstack((np.array([5.0]), np.arange(1e1, 6e1, 1e1), np.array([60.,65.0,66.]), np.arange(7e1,2.1e2,1e1))) 

points = len(pressure)

volume = np.zeros(points)
Nads = np.zeros(points)
hads = np.zeros(points)
vdw_sf = np.zeros(points)
vdw_ff = np.zeros(points)
coulomb_sf = np.zeros(points)
coulomb_ff = np.zeros(points)
ewald_sf = np.zeros(points)
ewald_ff = np.zeros(points)

for k in range(points):

    with open(f'pressure_{pressure[k]:.2f}_Pa/output.data', 'r') as file:
        content = file.read()

    box_pattern = re.compile(r'SIMULATION BOX PARAMETERS')
    volume_pattern = re.compile(r'Box Volume:\s+([-+]?[0-9]*\.?[0-9]+)')
    match = box_pattern.search(content)
    matches = list(volume_pattern.finditer(content, match.end()))
    volume[k] = float(matches[0].group(1)) 

    loading_pattern = re.compile(r'# MOLECULES')
    average_pattern = re.compile(r'Overall: Average:\s+([-+]?[0-9]*\.?[0-9]+)')
    match = loading_pattern.search(content)
    matches = list(average_pattern.finditer(content, match.end()))
    Nads[k] = float(matches[1].group(1))  

    heat_pattern = re.compile(r'HEAT OF ADSORPTION')
    average_pattern = re.compile(r'Overall: Average:\s+([-+]?[0-9]*\.?[0-9]+)')
    match = heat_pattern.search(content)
    matches = list(average_pattern.finditer(content, match.end()))
    hads[k] = float(matches[0].group(1))

    energy_pattern = re.compile(r'PRODUCTION PHASE AVERAGE ENERGY')

    vdw_sf_pattern = re.compile(r'VDW \[Host-Guest\]:\s+([-+]?[0-9]*\.?[0-9]+)')
    match = energy_pattern.search(content)
    matches = list(vdw_sf_pattern.finditer(content, match.end()))
    vdw_sf[k] = float(matches[0].group(1)) 

    vdw_ff_pattern = re.compile(r'VDW \[Guest-Guest\]:\s+([-+]?[0-9]*\.?[0-9]+)')
    match = energy_pattern.search(content)
    matches = list(vdw_ff_pattern.finditer(content, match.end()))
    vdw_ff[k] = float(matches[0].group(1)) 

    coulomb_sf_pattern = re.compile(r'Real Coulomb \[Host-Guest\]:\s+([-+]?[0-9]*\.?[0-9]+)')
    match = energy_pattern.search(content)
    matches = list(coulomb_sf_pattern.finditer(content, match.end()))
    coulomb_sf[k] = float(matches[0].group(1)) 

    coulomb_ff_pattern = re.compile(r'Real Coulomb \[Guest-Guest\]:\s+([-+]?[0-9]*\.?[0-9]+)')
    match = energy_pattern.search(content)
    matches = list(coulomb_ff_pattern.finditer(content, match.end()))
    coulomb_ff[k] = float(matches[0].group(1))

    ewald_sf_pattern = re.compile(r'Ewald \[Host-Guest\]:\s+([-+]?[0-9]*\.?[0-9]+)')
    match = energy_pattern.search(content)
    matches = list(ewald_sf_pattern.finditer(content, match.end()))
    ewald_sf[k] = float(matches[0].group(1))

    ewald_ff_pattern = re.compile(r'Ewald \[Guest-Guest\]:\s+([-+]?[0-9]*\.?[0-9]+)')
    match = energy_pattern.search(content)
    matches = list(ewald_ff_pattern.finditer(content, match.end()))
    ewald_ff[k] = float(matches[0].group(1))

plt.plot(pressure[:points], Nads/8,'o-', color=colors[3], 
         lw=1.75, ms=8, mec='k')
plt.xlabel('Pressure (Pa)', fontsize=18)
plt.ylabel('Adsorption (molecules/uc)',fontsize=18)
plt.minorticks_on()
# plt.xscale('log')
plt.tick_params(direction='in',right=True, top=True)
plt.tick_params(labelsize=18)
plt.tick_params(labelbottom=True, labeltop=False, labelright=False, labelleft=True)
plt.tick_params(direction='in',which='minor', length=4, bottom=True, top=True, left=True, right=True)
plt.tick_params(direction='in',which='major', length=8, bottom=True, top=True, left=True, right=True)
plt.legend(fontsize=16, frameon=False, edgecolor='k')
plt.show()

plt.plot(pressure[:points], 0.5*volume**(1/3),'-o', color=colors[3], 
         lw=1.75, ms=8, mec='k')
plt.xlabel('Pressure (Pa)', fontsize=18)
plt.ylabel('Unit cell lenght',fontsize=18)
plt.minorticks_on()
# plt.xscale('log')
plt.tick_params(direction='in',right=True, top=True)
plt.tick_params(labelsize=18)
plt.tick_params(labelbottom=True, labeltop=False, labelright=False, labelleft=True)
plt.tick_params(direction='in',which='minor', length=4, bottom=True, top=True, left=True, right=True)
plt.tick_params(direction='in',which='major', length=8, bottom=True, top=True, left=True, right=True)
plt.legend(fontsize=16, frameon=False, edgecolor='k')
plt.show()

plt.plot(pressure[:points], hads,'-o', color=colors[3], 
         lw=1.75, ms=8, mec='k')
plt.xlabel('Pressure (Pa)', fontsize=18)
plt.ylabel('Heat of adsorption (kJ/mol)',fontsize=18)
plt.minorticks_on()
# plt.xscale('log')
plt.tick_params(direction='in',right=True, top=True)
plt.tick_params(labelsize=18)
plt.tick_params(labelbottom=True, labeltop=False, labelright=False, labelleft=True)
plt.tick_params(direction='in',which='minor', length=4, bottom=True, top=True, left=True, right=True)
plt.tick_params(direction='in',which='major', length=8, bottom=True, top=True, left=True, right=True)
plt.legend(fontsize=16, frameon=False, edgecolor='k')
plt.show()

plt.plot(pressure[:points], (vdw_sf+coulomb_sf+ewald_sf)/(Na*kB*1e-1),'-o', color=colors[2], 
         lw=1.75, ms=8, mec='k', label='Solid-fluid')
plt.plot(pressure[:points], (vdw_ff+coulomb_ff+ewald_ff)/(Na*kB*1e-1),'-o', color=colors[4], 
         lw=1.75, ms=8, mec='k', label='Fluid-fluid')
plt.xlabel('Pressure (Pa)', fontsize=18)
plt.ylabel('Energy (K)',fontsize=18)
plt.minorticks_on()
# plt.xscale('log')
plt.tick_params(direction='in',right=True, top=True)
plt.tick_params(labelsize=18)
plt.tick_params(labelbottom=True, labeltop=False, labelright=False, labelleft=True)
plt.tick_params(direction='in',which='minor', length=4, bottom=True, top=True, left=True, right=True)
plt.tick_params(direction='in',which='major', length=8, bottom=True, top=True, left=True, right=True)
plt.legend(fontsize=16, frameon=False, edgecolor='k')
plt.show()

df = pd.DataFrame({
    'Pressure (Pa)': pressure[:points],
    'Adsorption (molecules/uc)': Nads/8,
    'Unit cell length (Angstrom)': 0.5*volume**(1/3),
    'Heat of adsorption (kJ/mol)': hads,
    'Solid-fluid energy (K)': (vdw_sf + coulomb_sf + ewald_sf)/(Na*kB*1e-1),
    'Fluid-fluid energy (K)': (vdw_ff + coulomb_ff + ewald_ff)/(Na*kB*1e-1)
})

df.to_csv(f'Ar_{framework_name}_{T:.2f}K.csv')
