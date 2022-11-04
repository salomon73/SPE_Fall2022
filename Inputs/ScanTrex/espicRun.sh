#!/bin/bash -l
#SBATCH --job-name=espic2d
#SBATCH --time=12:00:00
#SBATCH -n 4
#SBATCH -c 4
#srun amplxe-cl -c hotspots -r results_$nbprocs -- ../src/espic2d < stable.in


export I_MPI_PLATFORM=auto
export espicsrc='/home/sguincha/espic2d/src'
export espicwk='/home/sguincha/espic2d/wk'
res_folder='/scratch/sguincha/ScanTRex/Results/'
mkdir -p $res_folder
cat > job.in << EOM 
A Skeleton for MPI Time-Dependent program
=========================================
T.M. Tran   SPC/EPFL
-
&BASIC
  job_time=21400.0, extra_time=200.0,
  nrun=500,  !# of steps
  nlres=t,
  newres=f,
  femorder=3,3,
  ngauss=6,6,
  nlppform=.TRUE.
  partperiodic=f
  nlxg=f,       ! Display graphical interface
  nlclassical=t, ! Solve classical equations of motion
  nz=480,        !# of intervals in z
  nnr=50,130,35  !# of intervals between radii(1) and radii(2), between radii(2) and radii(3) and between radii(3) and radii(4)
  it0d=10,
  it2d=100,
  ittext=100,
  itparts=50000,
  ittracer=1, !10
  itgraph=100, 
  nbcelldiag=1,
  itcelldiag=20,
  resfile='${res_folder}resultrestart_5e-12.h5'
  rstfile='${res_folder}restartrestart_5e-12.h5'
  nbspecies=3,
  partfile='electrons_gauss.in','electron_tracers.in','ion_tracers.in',
  nbaddtestspecies = 2,
  addedtestspecfile = 'electrons.in', 'electronV0.in',
  radii=0.009, 0.019, 0.033,0.04  !radius of inner metallic wall, plasma boundaries and outer metallic wall
  distribtype=7          ! 1: uniform RZ gaussian in V, 2: stable eq 4.85 from Davidson
  lz=0.25,0.445,    		       !Cylinder length
  nlfreezephi=t, !pas de calcul à chaque t de E
  nplasma=2116800, !132300,   		! 100000# of superparticles
  nlmaxwellsource=t,
  nlPhis=t,           ! if false deactivate calculation of self electric fieldq
  potinn=-20000,       		!potential at inner wall [V]
  potout=0,     !30000,      		!potential at outer wall [V]
  B0=0.28,     			!Davidson chapter 4 formula 4.89 [T]
  n0=-3e17, dt=1E-12                    !Initial plasma density in plasmadim
  magnetfile='${res_folder}10T_DNPW.h5'
  bscaling=0,
/

&maxwellsourceparams
frequency=5E10, 
temperature=22000, 
rlimits=0.01,0.027, 
zlimits=0.35,0.4
time_start=-1.0,
radialtype=2
time_end=-1.0
/

&celldiagparams
specieid=1,
rindex=112,
zindex=117,
/

&geomparams
!x=0.028/2*cos(t)+0.012;y=0.002*sin(t)+0.081
r_a=0.01
r_bLeft=0.02358
alpha=0.1745
z_0=0.375
r_0=0.028
z_r=0.025
r_r=0.005
r_b=0.028
r_bRight=0.0375
above2=-1
above1=1
Interior=-1
walltype=6 ! ellipse plus flat bottom and flat top
nlweb=t
testkr=5,
testkz=5,
/

&neutcolparams
neutdens=2.5e21,
!Neon parameters
Eion=15.418,
scatter_fac=8.3,
nlcol=f, !f no coll neutral electrons
io_cross_sec_file='${espicwk}/H2_io_cross_sec.in',
ela_cross_sec_file='${espicwk}/H2_ela_cross_sec.in',
/

&psupplyparams
expneutdens=2.5e18,
PsResistor=1e6
geomcapacitor=337e-12
targetbias=20000
nbhdt=5
active=false
bdpos=-1,1
/
EOM
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#srun -n $nbprocs ../src/espic2d "$VAR"
mpirun -np 4  /home/sguincha/espic2d/src/espic2d job.in
