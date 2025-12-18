def write_mdp(temperature,pressure,rc,nsteps,tau_p,gen_vel):
    with open('npt.mdp', 'w+') as f:
        f.write(f"""; Run control
integrator              = md        ; leap-frog      
nsteps                  = {nsteps}   ; 2 * 50000 = 100 ps
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
ewald-rtol              = 1e-6

; Temperature coupling
tcoupl                  = v-rescale
tc-grps                 = System
tau_t                   = 0.1
ref_t                   = {temperature:.2f}

; Pressure coupling 
pcoupl                  = c-rescale
pcoupltype              = isotropic
tau_p                   = {tau_p:.1f}
ref_p                   = {1e-5*pressure:.2f}
compressibility         = 4.5e-5
refcoord-scaling        = all

; Velocity generation
gen_vel                 = {gen_vel}
gen-temp                = {temperature:.2f}
gen-seed                = -1 
""")
