#!/bin/bash -l
#SBATCH --job-name=espic2d
#SBATCH --time=12:00:00
#SBATCH --hint=nomultithread
#SBATCH -c 4
#SBATCH -n 2
#srun amplxe-cl -c hotspots -r results_8 -- ../src/espic2d < stable8.in
#export I_MPI_HYDRA_BOOTSTRAP=ssh
export FI_PROVIDER=tcp
export I_MPI_DEBUG=1000
cat >job.in  << EOM
A Skeleton for MPI Time-Dependent program
=========================================
T.M. Tran   SPC/EPFL
-
&BASIC
  job_time=144000.0, extra_time=360.0,
  nrun=1000,  !# of steps
  nlres=f,
  femorder=3,3,
  ngauss=4,4,
  nlppform=.TRUE.
  nlxg=f,       ! Display graphical interface
  nlclassical=t, ! Solve classical equations of motion
  dt=5E-11,     !timestep [s]
  nz=128,        !# of intervals in z
  nnr=10,40,30  !# of intervals between radii(1) and radii(2), between radii(2) and radii(3) and between radii(3) and radii(4)
  it0d=10,
  it2d=100,
  ittext=100,
  itparts=1,
  itgraph=100,
  resfile='stable_dt_11.h5'

  radii=0.001,6e-3,9e-3,0.01  !radius of inner metallic wall, plasma boundaries and outer metallic wall
  plasmadim=-0.10, 0.10, 0.0070, 0.0074,     !zmin zmax rmin rmax of initial plasma loading
  distribtype=7         ! 1: uniform RZ gaussian in V, 2: stable eq 4.85 from Davidson
  lz=-0.32,0.32,    		       !Cylinder length
  nbspecies=1,
  partfile='ion_tracers.in',
  nplasma=529200 !132300,   		! 100000# of superparticles
  npartsalloc=600000
  nblock=1000,
  nlPhis=t,           ! Deactivate calculation of self electric field
!  nplasma=100
  potinn=-20000,       		!potential at inner wall [V]
  potout=0,     !30000,      		!potential at outer wall [V]
  B0=0.21,     			!Davidson chapter 4 formula 4.89 [T]
  Rcurv=1,!1.001,           !Davidson chapter 4 formula 4.89 
  width=0.64,
!  n0=-5e13, dt=8E-13                    !Initial plasma density in plasmadim
!  n0=-4E16, dt=8E-12
!  n0=-5.e15, dt=8E-12,
!  n0=-3.e19, dt=8e-13,
!  n0=-5.e8
  temp=10000,                   !Initial temperature of plasma [K] (thermal velocity)
  !H0=3.2e-14,
  !P0=8.66e-25,
/

&geomparams
r_a=0.001,
r_b=0.01
nlweb=t,
walltype=0,
above1=1,
above2=-1,
/


EOM
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
mpirun -np 2  ../../src/espic2d job.in


