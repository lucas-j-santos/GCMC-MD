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
T = 195.0

pressure = np.array([  1000.0000,   1599.8684,   5332.8947,   9865.8553,  13972.1842,
         15492.0592, 16000., 16705.2928,  18078.5132,  19331.7434,  19998.3553,
         20664.9671,  21998.1908,  23331.4145,  25331.2500,  28664.3092,
         33330.5921,  40130.0329])

points = len(pressure)
iterations = 100

Lz_it = np.zeros(iterations-1)
Lz = np.zeros(points)

volume_it = np.zeros(iterations)
volume = np.zeros(points)

Nads_it = np.zeros(iterations)
Nads = np.zeros(points)

hads_it = np.zeros(iterations)
hads = np.zeros(points)

vdw_sf_it = np.zeros(iterations)
vdw_sf = np.zeros(points)

vdw_ff_it = np.zeros(iterations)
vdw_ff = np.zeros(points)

coulomb_sf_it = np.zeros(iterations)
coulomb_sf = np.zeros(points)

coulomb_ff_it = np.zeros(iterations)
coulomb_ff = np.zeros(points)

ewald_sf_it = np.zeros(iterations)
ewald_sf = np.zeros(points)

ewald_ff_it = np.zeros(iterations)
ewald_ff = np.zeros(points)

for k in range(points):

    for i in range(iterations-1):
        
        df = pd.read_csv(f'pressure_{1e-3*pressure[k]:.2f}_kPa/npt_output_{i}.csv')
        Lz_it[i] = df['Box-Z'].mean()

    for i in range(iterations):

        with open(f'pressure_{1e-3*pressure[k]:.2f}_kPa/gcmc_output_{i}.data', 'r') as file:
            content = file.read()

        box_pattern = re.compile(r'SIMULATION BOX PARAMETERS')
        volume_pattern = re.compile(r'Box Volume:\s+([-+]?[0-9]*\.?[0-9]+)')
        match = box_pattern.search(content)
        matches = list(volume_pattern.finditer(content, match.end()))
        volume_it[i] = float(matches[0].group(1)) 

        loading_pattern = re.compile(r'# MOLECULES')
        average_pattern = re.compile(r'Overall: Average:\s+([-+]?[0-9]*\.?[0-9]+)')
        match = loading_pattern.search(content)
        matches = list(average_pattern.finditer(content, match.end()))
        Nads_it[i] = float(matches[1].group(1))  

        heat_pattern = re.compile(r'HEAT OF ADSORPTION')
        average_pattern = re.compile(r'Overall: Average:\s+([-+]?[0-9]*\.?[0-9]+)')
        match = heat_pattern.search(content)
        matches = list(average_pattern.finditer(content, match.end()))
        hads_it[i] = float(matches[0].group(1))

        energy_pattern = re.compile(r'PRODUCTION PHASE AVERAGE ENERGY')

        vdw_sf_pattern = re.compile(r'VDW \[Host-Guest\]:\s+([-+]?[0-9]*\.?[0-9]+)')
        match = energy_pattern.search(content)
        matches = list(vdw_sf_pattern.finditer(content, match.end()))
        vdw_sf_it[i] = float(matches[0].group(1)) 

        vdw_ff_pattern = re.compile(r'VDW \[Guest-Guest\]:\s+([-+]?[0-9]*\.?[0-9]+)')
        match = energy_pattern.search(content)
        matches = list(vdw_ff_pattern.finditer(content, match.end()))
        vdw_ff_it[i] = float(matches[0].group(1)) 

        coulomb_sf_pattern = re.compile(r'Real Coulomb \[Host-Guest\]:\s+([-+]?[0-9]*\.?[0-9]+)')
        match = energy_pattern.search(content)
        matches = list(coulomb_sf_pattern.finditer(content, match.end()))
        coulomb_sf_it[i] = float(matches[0].group(1)) 

        coulomb_ff_pattern = re.compile(r'Real Coulomb \[Guest-Guest\]:\s+([-+]?[0-9]*\.?[0-9]+)')
        match = energy_pattern.search(content)
        matches = list(coulomb_ff_pattern.finditer(content, match.end()))
        coulomb_ff_it[i] = float(matches[0].group(1))

        ewald_sf_pattern = re.compile(r'Ewald \[Host-Guest\]:\s+([-+]?[0-9]*\.?[0-9]+)')
        match = energy_pattern.search(content)
        matches = list(ewald_sf_pattern.finditer(content, match.end()))
        ewald_sf_it[i] = float(matches[0].group(1))

        ewald_ff_pattern = re.compile(r'Ewald \[Guest-Guest\]:\s+([-+]?[0-9]*\.?[0-9]+)')
        match = energy_pattern.search(content)
        matches = list(ewald_ff_pattern.finditer(content, match.end()))
        ewald_ff_it[i] = float(matches[0].group(1))

    volume[k] = np.mean(volume_it[20:])
    Lz[k] = np.mean(Lz_it[20:])
    Nads[k] = np.mean(Nads_it[20:])
    hads[k] = np.mean(hads_it[20:])
    vdw_sf[k] = np.mean(vdw_sf_it[20:])
    vdw_ff[k] = np.mean(vdw_ff_it[20:])
    coulomb_sf[k] = np.mean(coulomb_sf_it[20:])
    coulomb_ff[k] = np.mean(coulomb_ff_it[20:])
    ewald_sf[k] = np.mean(ewald_sf_it[20:])
    ewald_ff[k] = np.mean(ewald_ff_it[20:])

plt.plot(1e-3*pressure[:points], Nads,'-o', color=colors[3], 
         lw=1.75, ms=8, mec='k')
plt.xlabel('Pressure (kPa)', fontsize=18)
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

plt.plot(1e-3*pressure[:points], Lz,'-o', color=colors[3], 
         lw=1.75, ms=8, mec='k')
plt.xlabel('Pressure (kPa)', fontsize=18)
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

plt.plot(1e-3*pressure[:points], hads,'-o', color=colors[3], 
         lw=1.75, ms=8, mec='k')
plt.xlabel('Pressure (kPa)', fontsize=18)
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

plt.plot(1e-3*pressure[:points], (vdw_sf+coulomb_sf+ewald_sf)/(Na*kB*1e-1),'-o', color=colors[2], 
         lw=1.75, ms=8, mec='k', label='Solid-fluid')
plt.plot(1e-3*pressure[:points], (vdw_ff+coulomb_ff+ewald_ff)/(Na*kB*1e-1),'-o', color=colors[4], 
         lw=1.75, ms=8, mec='k', label='Fluid-fluid')
plt.xlabel('Pressure (kPa)', fontsize=18)
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
    'Pressure (kPa)': 1e-3*pressure[:points],
    'Adsorption (molecules/uc)': Nads,
    'Unit cell length (Angstrom)': Lz,
    'Heat of adsorption (kJ/mol)': hads,
    'Solid-fluid energy (K)': (vdw_sf + coulomb_sf + ewald_sf)/(Na*kB*1e-1),
    'Fluid-fluid energy (K)': (vdw_ff + coulomb_ff + ewald_ff)/(Na*kB*1e-1)
})

df.to_csv(f'CO2_{framework_name}_{T:.2f}K_GCMC-MD.csv')
