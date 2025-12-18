gmx grompp -f npt.mdp -c box.gro -p topol.top -o npt.tpr -maxwarn 2
gmx mdrun -v -deffnm npt
gmx energy -f npt.edr -o temperature.xvg
gmx energy -f npt.edr -o pressure.xvg
