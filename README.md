# Ensemble-Average-MSD
This Fortran code calculates the ensemble average mean square displacement (MSD) of atoms from LAMMPS trajectory file.

The CONTROL file contains MD parameters such as
LAMMPS trajectory file location (InFile)
Location of MSD output file (OutFile)
no. of molecules in simulation box (MoleculePerBox)
no. of atoms per molecule (AtomsPerMolecule)
first time step in LAMMPS trajectory file (FirstStep)
last time step in LAMMPS trajectory file (LastStep)
how often is the trajectory sampled (TimeStep)
no. of equilibration steps that should be ignored for MSD calculation (Step_Equilibration)
last time step at which MSD should be calculated (LastStep_MSD)
Note: while calculating ensemble average, it is a good idea to have LastStep_MSD ≈ 1/2(LastStep).

The ensemble average MSD is calculated as:

![image](https://github.com/grvkr0807/Ensemble-Average-MSD/assets/35864049/85103494-224f-440d-8c9f-c0e4eef9804d)
