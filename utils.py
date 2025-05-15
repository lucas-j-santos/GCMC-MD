def write_mdp(temperature,pressure,rc):
    with open('npt.mdp', 'w+') as f:
        f.write(f"""title                   = NPT      
; Run control
integrator              = md        ; leap-frog integrator
nsteps                  = 100000   ; 2 * 50000 = 100 ps
dt                      = 0.001     ; 1 fs
; Output control
nstxout                 = 500       ; save coordinates every 1.0 ps
nstvout                 = 500       ; save velocities every 1.0 ps
nstenergy               = 500       ; save energies every 1.0 ps
nstlog                  = 500       ; update log file every 1.0 ps
; Neighbor searching
cutoff-scheme           = Verlet    
nstlist                 = 10       
pbc                     = xyz      
; Electrostatics
coulombtype             = PME
coulomb-modifier        = Potential-shift
rcoulomb                = {1e-1*rc:.2f}     
; Van der Waals
vdwtype                 = Cut-off 
vdw-modifier            = Potential-shift
rvdw                    = {1e-1*rc:.2f}
DispCorr                = no
; Ewald
fourierspacing          = 0.16
pme_order               = 4
; Temperature coupling
tcoupl                  = v-rescale
tc-grps                 = System
tau_t                   = 0.1	
ref_t                   = {temperature:.1f}
; Pressure coupling 
pcoupl                  = c-rescale
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = {1e-5*pressure:.1e}
compressibility         = 4.5e-5
refcoord-scaling        = com
; Velocity generation
gen_vel                 = no
""")
        
def write_topol(framework_name,molecule_name,N):
    with open('topol.top', 'w+') as f:
        f.write(f""";      
[ defaults ]
; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
1               2               no              0            0           

[ atomtypes ]
; name    at.num    mass    charge ptype  sigma      epsilon
OC             8  15.999430  0.00000000  A          0.305     0.65684255
CO             6  12.010780  0.00000000  A           0.28     0.22449049
Zn_a           0  65.380000  0.00000000  A           0.27   0.0034920731
H_a            0   1.007900  0.00000000  A          0.285    0.063605618
C_a            0  12.011000  0.00000000  A          0.374     0.39077961
C_b            0  12.011000  0.00000000  A          0.347     0.39793005
C_c            0  12.011000  0.00000000  A          0.347     0.39793005
O_a            0  15.999400  0.00000000  A          0.298      5.8201219
O_b            0  15.999400  0.00000000  A          0.311     0.58616942


#include "../{framework_name}.itp"
#include "../{molecule_name}_trappe.itp"

[ system ]
; Name
Flexible {framework_name}

[ molecules ]
; Compound       #mols
{framework_name} 1
{molecule_name}  {N:}
""")

def write_inp(framework_name,temperature,pressure,rc,molecule_name,fugacity_coefficient,create_molecules):

    with open('simulation.input', 'w+') as f:

        f.write(f"""NumberOfInitializationCycles 100      
NumberOfEquilibrationCycles  0
NumberOfProductionCycles     900

UseMaxStep  no
MaxStepPerCycle 1

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
