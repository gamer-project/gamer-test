


# ===========================================================================
# INSTRUCTIONS
#
# To compile GAMER, please set the following configurations properly:
#
# (1) simulation options
# (2) compiler and flags
# (3) library paths
# ===========================================================================



# executable
#######################################################################################################
EXECUTABLE := gamer



# output detailed compilation commands (0/1 = off/on)
#######################################################################################################
COMPILE_VERBOSE := 0



# simulation options
#######################################################################################################

# (a) physical modules
# -------------------------------------------------------------------------------
# physical model: HYDRO/ELBDM/PAR_ONLY (ELBDM: wave dark matter; PAR_ONLY: particle-only)
# --> must be set in any cases; PAR_ONLY is not supported yet;
#     ELBDM is not publicly available yet but will be released soon
SIMU_OPTION += -DMODEL=HYDRO

# enable gravity
#SIMU_OPTION += -DGRAVITY

# enable particles
# --> must enable GRAVITY
#SIMU_OPTION += -DPARTICLE

# support Grackle, a chemistry and radiative cooling library
# --> must set NCOMP_PASSIVE_USER according to the primordial chemistry network set by GRACKLE_PRIMORDIAL
# --> please enable OpenMP when compiling Grackle (by "make omp-on")
#SIMU_OPTION += -DSUPPORT_GRACKLE


# (b) options of different physical modules
# (b-1) hydro options
# ------------------------------------------------------------------------------------
ifeq "$(filter -DMODEL=HYDRO, $(SIMU_OPTION))" "-DMODEL=HYDRO"
# hydrodynamic scheme: RTVD/MHM/MHM_RP/CTU
# --> must be set when MODEL=HYDRO
# --> MHD only supports MHM_RP and CTU
SIMU_OPTION += -DFLU_SCHEME=CTU

# scheme of spatial data reconstruction: PLM/PPM (piecewise-linear/piecewise-parabolic)
# --> useless for RTVD
SIMU_OPTION += -DLR_SCHEME=PPM

# Riemann solver: EXACT/ROE/HLLE/HLLC/HLLD
# --> pure hydro: EXACT/ROE/HLLE/HLLC(*); MHD: ROE/HLLE/HLLD(*);
#                 (recommended solvers are highlighted by *)
# --> useless for RTVD
SIMU_OPTION += -DRSOLVER=ROE

# dual energy formalism: DE_ENPY/DE_EINT (evolve entropy or internal energy)
# --> DE_EINT is not supported yet; useless for RTVD
#SIMU_OPTION += -DDUAL_ENERGY=DE_ENPY

# number of user-defined passively advected scalars
# --> set it to 0 or comment it out if none is required
# --> useless for RTVD
SIMU_OPTION += -DNCOMP_PASSIVE_USER=0

# magnetohydrodynamics
#SIMU_OPTION += -DMHD

# cosmic rays (not supported yet)
#SIMU_OPTION += -DCOSMIC_RAY

# equation of state: EOS_GAMMA/EOS_ISOTHERMAL/EOS_NUCLEAR/EOS_TABULAR/EOS_USER
# --> must be set when MODEL=HYDRO; must also enable BAROTROPIC_EOS for EOS_ISOTHERMAL
SIMU_OPTION += -DEOS=EOS_GAMMA

# whether or not the EOS set above is barotropic
# --> mandatory for EOS_ISOTHERMAL; optional for EOS_TABULAR/EOS_USER
#SIMU_OPTION += -DBAROTROPIC_EOS


# (b-2) ELBDM options
# ------------------------------------------------------------------------------------
else ifeq "$(filter -DMODEL=ELBDM, $(SIMU_OPTION))" "-DMODEL=ELBDM"
# enforce the mass conservation
SIMU_OPTION += -DCONSERVE_MASS

# 4-th order Laplacian
SIMU_OPTION += -DLAPLACIAN_4TH

# include the quartic self-interaction potential
# --> must turn on GRAVITY; does not support COMOVING
#SIMU_OPTION += -DQUARTIC_SELF_INTERACTION

# number of passively advected scalars
# --> not supported yet and must be set to 0 (or comment it out)
SIMU_OPTION += -DNCOMP_PASSIVE_USER=0

endif # MODEL


# (b-3) gravity options
# -------------------------------------------------------------------------------
ifeq "$(filter -DGRAVITY, $(SIMU_OPTION))" "-DGRAVITY"
# Poisson solver: SOR/MG (successive-overrelaxation (recommended)/multigrid)
# --> must be set when GRAVITY is enabled
SIMU_OPTION += -DPOT_SCHEME=SOR

# store GRA_GHOST_SIZE ghost-zone potential for each patch on each side
# --> recommended when PARTICLE is enabled for improving accuracy for particles around the patch boundaries
# --> must be enabled for STAR_FORMATION + STORE_PAR_ACC
SIMU_OPTION += -DSTORE_POT_GHOST

# use unsplitting method to couple gravity to the target model (recommended)
# --> for HYDRO only; ELBDM is not supported yet
SIMU_OPTION += -DUNSPLIT_GRAVITY

# comoving frame for cosmological simulations
# --> not supported for PARTICLE
#SIMU_OPTION += -DCOMOVING
endif # GRAVITY


# (b-4) particle options
# ------------------------------------------------------------------------------------
ifeq "$(filter -DPARTICLE, $(SIMU_OPTION))" "-DPARTICLE"
# store particle acceleration (recommended)
SIMU_OPTION += -DSTORE_PAR_ACC

# allow for creating new particles after initialization
# --> must turn on STORE_POT_GHOST when STORE_PAR_ACC is adopted
#SIMU_OPTION += -DSTAR_FORMATION

# feedback from particles to grids
#SIMU_OPTION += -DFEEDBACK

# number of user-defined particle attributes
# --> set it to 0 or comment it out if none is required
SIMU_OPTION += -DPAR_NATT_USER=0
endif # PARTICLE


# (c) miscellaneous options (AMR, GPU, parallelization, optimizations, ...)
# ------------------------------------------------------------------------------------
# maximum number of AMR levels including the root level (level id = 0 ... NLEVEL-1)
# --> must be set in any cases
SIMU_OPTION += -DNLEVEL=10

# maximum number of patches on each AMR level
# --> must be set in any cases
SIMU_OPTION += -DMAX_PATCH=1000000

# GPU acceleration
# --> must set GPU_ARCH as well
#SIMU_OPTION += -DGPU

# GPU architecture: FERMI/KEPLER/MAXWELL/PASCAL/VOLTA/TURING
SIMU_OPTION += -DGPU_ARCH=KEPLER

# debug mode
#SIMU_OPTION += -DGAMER_DEBUG

# bitwise reproducibility
#SIMU_OPTION += -DBITWISE_REPRODUCIBILITY

# measure the wall-clock time of GAMER
SIMU_OPTION += -DTIMING

# measure the wall-clock time of GPU solvers (will disable CPU/GPU overlapping)
# --> must enable TIMING
#SIMU_OPTION += -DTIMING_SOLVER

# double precision
#SIMU_OPTION += -DFLOAT8

# serial mode (in which no MPI libraries are required)
# --> must disable LOAD_BALANCE
SIMU_OPTION += -DSERIAL

# load-balance parallelization: HILBERT
# --> must disable SERIAL
#SIMU_OPTION += -DLOAD_BALANCE=HILBERT

# overlap MPI communication with computation
# --> NOT supported yet; must enable LOAD_BALANCE
#SIMU_OPTION += -DOVERLAP_MPI

# enable OpenMP parallelization
SIMU_OPTION += -DOPENMP

# work on the NAOC Laohu GPU cluster
#SIMU_OPTION += -DLAOHU

# support HDF5 format
#SIMU_OPTION += -DSUPPORT_HDF5

# support GNU scientific library
#SIMU_OPTION += -DSUPPORT_GSL

# support yt inline analysis
#SIMU_OPTION += -DSUPPORT_LIBYT

# random number implementation: RNG_GNU_EXT/RNG_CPP11 (GNU extension drand48_r/c++11 <random>)
# --> use RNG_GNU_EXT for compilers supporting GNU extensions (**not supported on some macOS**)
#     use RNG_CPP11   for compilers supporting c++11 (**may need to add -std=c++11 to CXXFLAG**)
SIMU_OPTION += -DRANDOM_NUMBER=RNG_GNU_EXT





# compiler and flags
#######################################################################################################

# intel
# -------------------------------------------------------------------------------
# CXX         = icpc                                      # serial compiler
##CXX         = $(MPI_PATH)/bin/mpicxx                    # MPI compiler
##CXX         = CC                                        # cray wrapper script
# CXXFLAG     = -g -O3                                    # general flags
##CXXFLAG     = -g -O3 -std=c++11 #-gxx-name=YOUR_G++
##CXXFLAG     = -g -fast
# CXXFLAG    += -w1                                       # warning flags
# OPENMPFLAG  = -fopenmp                                  # openmp flag
# LIB         = -limf                                     # libraries and linker flags
#
## for debug only
#ifeq "$(filter -DGAMER_DEBUG, $(SIMU_OPTION))" "-DGAMER_DEBUG"
##CXXFLAG    += -fstack-protector-all
##LIB        += -lssp
#endif
#
## suppress warning when OpenMP is disabled
#ifeq "$(filter -DOPENMP, $(SIMU_OPTION))" ""
# CXXFLAG    += -Wno-unknown-pragmas -diag-disable 3180
#endif


# gnu
# -------------------------------------------------------------------------------
 CXX         = g++                                       # serial compiler
#CXX         = $(MPI_PATH)/bin/mpicxx                    # MPI compiler
#CXX         = CC                                        # cray wrapper script
 CXXFLAG     = -g -O3                                    # general flags
#CXXFLAG     = -g -O3 -std=c++11
#CXXFLAG     = -g -Ofast
 CXXFLAG    += -Wall -Wextra                             # warning flags
 CXXFLAG    += -Wno-unused-variable -Wno-unused-parameter \
               -Wno-maybe-uninitialized -Wno-unused-but-set-variable \
               -Wno-unused-result -Wno-unused-function
 OPENMPFLAG  = -fopenmp                                  # openmp flag
 LIB         =                                           # libraries and linker flags

# for debug only
ifeq "$(filter -DGAMER_DEBUG, $(SIMU_OPTION))" "-DGAMER_DEBUG"
#CXXFLAG    += -fstack-protector-all
endif

# suppress warning when OpenMP is disabled
ifeq "$(filter -DOPENMP, $(SIMU_OPTION))" ""
 CXXFLAG    += -Wno-unknown-pragmas
endif


# CUDA
# -------------------------------------------------------------------------------
 NVCC = $(CUDA_PATH)/bin/nvcc                            # CUDA compiler
#NVCC = $(CUDA_PATH)/bin/nvcc -ccbin YOUR_G++            # with a specified host compiler
#NVCC = nvcc -ccbin CC                                   # for the Cray environment





# library paths
#######################################################################################################

# template
CUDA_PATH    :=
FFTW_PATH    :=
MPI_PATH     :=
HDF5_PATH    :=
GRACKLE_PATH :=
GSL_PATH     :=
LIBYT_PATH   :=

# NCSA-bluewaters (only need to specify the path of user-installed packages)
#GRACKLE_PATH :=

# UIUC-campus (intel-15.0 + openmpi-1.8.4)
#CUDA_PATH    := /usr/local/cuda/7.0
#FFTW_PATH    := /projects/ncsa/grav/softwares/fftw/2.1.5
#MPI_PATH     := /usr/local/mpi/openmpi-1.8.4-intel-15.0
#HDF5_PATH    := /projects/ncsa/grav/softwares/miniconda2
#GRACKLE_PATH := /home/hyschive/local
#GSL_PATH     := /usr/local/gsl/gsl-1.16

# UIUC-campus (gcc-4.9.2 + openmpi-1.8.4)
#CUDA_PATH    := /usr/local/cuda/7.0
#FFTW_PATH    := /projects/ncsa/grav/softwares/fftw/2.1.5-openmpi-1.8.4-gcc-4.9.2
#MPI_PATH     := /usr/local/mpi/openmpi-1.8.4-gcc-4.9.2
#HDF5_PATH    := /projects/ncsa/grav/softwares/miniconda2
#GRACKLE_PATH := /home/hyschive/local
#GSL_PATH     := /usr/local/gsl/gsl-1.16

# NCSA-hybrid (gcc)
#CUDA_PATH    := /usr/local/cuda-8.0
#MPI_PATH     := /home/hyschive/project/software/openmpi/2.0.2
#FFTW_PATH    := /home/hyschive/project/software/fftw/2.1.5
#HDF5_PATH    := /home/hyschive/project/software/miniconda2
#GSL_PATH     :=

# NCSA-nano (intel 2016)
#CUDA_PATH    := /usr/local/cuda
#MPI_PATH     := /mnt/lustre/allen/hyschive/software/openmpi/2.0.2
#FFTW_PATH    := /mnt/lustre/allen/hyschive/software/fftw/2.1.5
#HDF5_PATH    := /mnt/lustre/allen/hyschive/software/miniconda2
#GSL_PATH     :=

# NCSA-dxl (gcc, serial)
#CUDA_PATH    :=
#FFTW_PATH    := /home/hyschive/software/fftw/2.1.5
#MPI_PATH     := /usr
#HDF5_PATH    := /home/hyschive/software/miniconda2
#GSL_PATH     := /usr/lib
#LIBYT_PATH   := /home/hyschive/project/libyt/libyt.dev/libyt

# NTU-eureka (default: openmpi-intel)
#CUDA_PATH    := /software/cuda/default
#FFTW_PATH    := /software/fftw/default
#MPI_PATH     := /software/openmpi/default
#HDF5_PATH    := /software/hdf5/default
#GRACKLE_PATH :=
#GSL_PATH     := /software/gsl/default
#LIBYT_PATH   :=

# NTU-hulk (openmpi-intel-qlc)
#CUDA_PATH    := /opt/gpu/cuda/4.2
#FFTW_PATH    := /opt/math/fftw/2.1.5-intel-openmpi-1.4.3-qlc
#MPI_PATH     := /usr/mpi/intel/openmpi-1.4.3-qlc
#HDF5_PATH    := /opt/hdf5/1.8.20
#GRACKLE_PATH := /work1/fish/Software/grackle
#GSL_PATH     := /opt/math/gsl/2.5
#LIBYT_PATH   := /project/fish/libyt/libyt
#FFTW_PATH    := /opt/math/fftw/2.1.5-intel-openmpi-1.6.0-qlc
#MPI_PATH     := /opt/mpi/openmpi/1.6.0-intel-qlc
#FFTW_PATH    := /opt/math/fftw/2.1.5-intel-openmpi-2.1.5-qlc
#MPI_PATH     := /opt/mpi/openmpi/2.1.5-intel-qlc
#CUDA_PATH    := /usr/local/cuda-7.5

# NTU-hulk (openmpi-gcc-qlc)
#CUDA_PATH    := /opt/gpu/cuda/4.2
#FFTW_PATH    := /opt/math/fftw/2.1.5-gcc-qlc
#MPI_PATH     := /usr/mpi/gcc/openmpi-1.4.3-qlc
#HDF5_PATH    := /opt/hdf5/1.8.9
#GRACKLE_PATH := /work1/fish/Software/grackle
#GSL_PATH     := /opt/math/gsl
#LIBYT_PATH   := /project/fish/libyt/libyt
#MPI_PATH     := /opt/mpi/openmpi/1.8.1-gcc-qlc
#CUDA_PATH    := /usr/local/cuda-7.5

# NTHU-fomalhaut (openmpi-gnu)
#CUDA_PATH    := /usr/local/cuda-10.0
#FFTW_PATH    := /home/hyschive/software/fftw/2.1.5-openmpi-gnu
#MPI_PATH     := /storage/app/gnu/openmpi-3.1.3-gcc8
#HDF5_PATH    := /storage/app/gnu/hdf5-1.8.21-gcc8
#GRACKLE_PATH :=
#GSL_PATH     :=

# NAOC-laohu (OpenMPI-1.3.2-intel)
#CUDA_PATH    := $(CUDA_INSTALL_PATH)
#FFTW_PATH    := /ifs/data/hsiyu/Software/fftw/2.1.5-openmpi-1.3.2-intel
#MPI_PATH     := /usr/mpi/intel/openmpi-1.3.2-qlc
#GPUID_PATH   := /ifs/data/hsiyu/Software/gpudevmgr

# A "typical" macOS installation without GPU
#CUDA_PATH    :=
#FFTW_PATH    := /usr/local/fftw-2.1.5
#MPI_PATH     := /usr/local/mpich-3.2
#HDF5_PATH    := ${HOME}/miniconda3
#GRACKLE_PATH :=
#GSL_PATH     := /usr/local/gsl-1.16

# NASA Pleiades (intel compilers and SGI MPI)
#CUDA_PATH    := /nasa/cuda/9.0/
#FFTW_PATH    := /u/jzuhone/opt/fftw-2.1.5
#MPI_PATH     := /nasa/sgi/mpt/2.14r19
#HDF5_PATH    := /nasa/hdf5/1.8.18_serial
#GRACKLE_PATH :=
#GSL_PATH     :=




# source files
#######################################################################################################

# common source files
# -------------------------------------------------------------------------------
# CUDA source files (compiled with nvcc)
GPU_FILE    := CUAPI_Asyn_FluidSolver.cu  CUAPI_DiagnoseDevice.cu  CUAPI_MemAllocate_Fluid.cu \
               CUAPI_MemFree_Fluid.cu  CUAPI_Set_Default_GPU_Parameter.cu  CUAPI_SetDevice.cu \
               CUAPI_Synchronize.cu  CUAPI_Asyn_dtSolver.cu  CUAPI_SetConstMemory.cu  CUAPI_SetConstMemory_EoS.cu


# C/C++ source files (compiled with c++ compiler)
CPU_FILE    := Main.cpp  EvolveLevel.cpp  InvokeSolver.cpp  Prepare_PatchData.cpp \
               InterpolateGhostZone.cpp

CPU_FILE    += Aux_Check_Parameter.cpp  Aux_Check_Conservation.cpp  Aux_Check.cpp  Aux_Check_Finite.cpp \
               Aux_Check_FluxAllocate.cpp  Aux_Check_PatchAllocate.cpp  Aux_Check_ProperNesting.cpp \
               Aux_Check_Refinement.cpp  Aux_Check_Restrict.cpp  Aux_Error.cpp  Aux_GetCPUInfo.cpp \
               Aux_GetMemInfo.cpp  Aux_Message.cpp  Aux_Record_PatchCount.cpp  Aux_TakeNote.cpp  Aux_Timing.cpp \
               Aux_Check_MemFree.cpp  Aux_Record_Performance.cpp  Aux_CheckFileExist.cpp  Aux_Array.cpp \
               Aux_Record_User.cpp  Aux_Record_CorrUnphy.cpp  Aux_SwapPointer.cpp  Aux_Check_NormalizePassive.cpp \
               Aux_LoadTable.cpp  Aux_IsFinite.cpp  Aux_ComputeProfile.cpp

CPU_FILE    += CPU_FluidSolver.cpp  Flu_AdvanceDt.cpp  Flu_Prepare.cpp  Flu_Close.cpp  Flu_FixUp_Flux.cpp \
               Flu_FixUp_Restrict.cpp  Flu_AllocateFluxArray.cpp  Flu_BoundaryCondition_User.cpp  Flu_ResetByUser.cpp \
               Flu_CorrAfterAllSync.cpp  Flu_ManageFixUpTempArray.cpp  Flu_DerivedField_BuiltIn.cpp \
               Flu_DerivedField_User.cpp

CPU_FILE    += End_GAMER.cpp  End_MemFree.cpp  End_MemFree_Fluid.cpp  End_StopManually.cpp  End_User.cpp \
               Init_BaseLevel.cpp  Init_GAMER.cpp  Init_Load_DumpTable.cpp \
               Init_Load_FlagCriteria.cpp  Init_Load_Parameter.cpp  Init_MemAllocate.cpp \
               Init_MemAllocate_Fluid.cpp  Init_Parallelization.cpp  Init_RecordBasePatch.cpp  Init_Refine.cpp \
               Init_ByRestart_v1.cpp  Init_ByFunction.cpp  Init_TestProb.cpp  Init_ByFile.cpp  Init_OpenMP.cpp \
               Init_ByRestart_HDF5.cpp  Init_ResetParameter.cpp  Init_ByRestart_v2.cpp  Init_MemoryPool.cpp \
               Init_Unit.cpp  Init_UniformGrid.cpp  Init_Field.cpp  Init_User.cpp

CPU_FILE    += Interpolate.cpp  Int_CQuadratic.cpp  Int_MinMod1D.cpp  Int_MinMod3D.cpp  Int_vanLeer.cpp \
               Int_Quadratic.cpp  Int_Table.cpp  Int_CQuartic.cpp  Int_Quartic.cpp

CPU_FILE    += Mis_CompareRealValue.cpp  Mis_GetTotalPatchNumber.cpp  Mis_GetTimeStep.cpp  Mis_Heapsort.cpp \
               Mis_BinarySearch.cpp  Mis_1D3DIdx.cpp  Mis_Matching.cpp  Mis_GetTimeStep_User.cpp \
               Mis_dTime2dt.cpp  Mis_CoordinateTransform.cpp  Mis_BinarySearch_Real.cpp  Mis_InterpolateFromTable.cpp \
               CPU_dtSolver.cpp  dt_Prepare_Flu.cpp  dt_Prepare_Pot.cpp  dt_Close.cpp  dt_InvokeSolver.cpp

CPU_FILE    += Output_DumpData_Total.cpp  Output_DumpData.cpp  Output_DumpManually.cpp  Output_PatchMap.cpp \
               Output_DumpData_Part.cpp  Output_FlagMap.cpp  Output_Patch.cpp  Output_PreparedPatch_Fluid.cpp \
               Output_PatchCorner.cpp  Output_Flux.cpp  Output_User.cpp  Output_BasePowerSpectrum.cpp \
               Output_DumpData_Total_HDF5.cpp  Output_L1Error.cpp

CPU_FILE    += Flag_Real.cpp  Refine.cpp   SiblingSearch.cpp  SiblingSearch_Base.cpp  FindFather.cpp \
               Flag_User.cpp  Flag_Check.cpp  Flag_Lohner.cpp  Flag_Region.cpp

CPU_FILE    += Table_01.cpp  Table_02.cpp  Table_03.cpp  Table_04.cpp  Table_05.cpp  Table_06.cpp \
               Table_07.cpp  Table_SiblingSharingSameEdge.cpp  Table_SiblingPatch.cpp

vpath %.cu     GPU_API
vpath %.cpp    Main  Init  Refine  Fluid  Interpolation  Tables  Output  Miscellaneous  Auxiliary


# hydrodynamic source files (included only if "MODEL=HYDRO")
# ------------------------------------------------------------------------------------
ifeq "$(filter -DMODEL=HYDRO, $(SIMU_OPTION))" "-DMODEL=HYDRO"
GPU_FILE    += CUFLU_dtSolver_HydroCFL.cu  CUFLU_FluidSolver_RTVD.cu  CUFLU_FluidSolver_MHM.cu  CUFLU_FluidSolver_CTU.cu \
               GPU_EoS_Gamma.cu  GPU_EoS_User_Template.cu  GPU_EoS_Isothermal.cu

CPU_FILE    += CPU_FluidSolver_RTVD.cpp  CPU_FluidSolver_MHM.cpp  CPU_FluidSolver_CTU.cpp \
               CPU_Shared_DataReconstruction.cpp  CPU_Shared_FluUtility.cpp  CPU_Shared_ComputeFlux.cpp \
               CPU_Shared_FullStepUpdate.cpp  CPU_Shared_RiemannSolver_Exact.cpp  CPU_Shared_RiemannSolver_Roe.cpp \
               CPU_Shared_RiemannSolver_HLLE.cpp  CPU_Shared_RiemannSolver_HLLC.cpp  CPU_Shared_DualEnergy.cpp \
               CPU_dtSolver_HydroCFL.cpp  CPU_EoS_Gamma.cpp  CPU_EoS_User_Template.cpp  CPU_EoS_Isothermal.cpp

CPU_FILE    += Hydro_Init_ByFunction_AssignData.cpp  Hydro_Aux_Check_Negative.cpp \
               Hydro_BoundaryCondition_Reflecting.cpp  Hydro_BoundaryCondition_Outflow.cpp \
               EoS_Init.cpp  EoS_End.cpp

vpath %.cu     Model_Hydro/GPU_Hydro  EoS  EoS/Gamma  EoS/User_Template  EoS/Isothermal
vpath %.cpp    Model_Hydro/CPU_Hydro  Model_Hydro  EoS  EoS/Gamma  EoS/User_Template  EoS/Isothermal

ifeq "$(filter -DGRAVITY, $(SIMU_OPTION))" "-DGRAVITY"
GPU_FILE    += CUPOT_HydroGravitySolver.cu  CUPOT_dtSolver_HydroGravity.cu

CPU_FILE    += CPU_HydroGravitySolver.cpp  CPU_dtSolver_HydroGravity.cpp

vpath %.cu     Model_Hydro/GPU_HydroGravity
vpath %.cpp    Model_Hydro/CPU_HydroGravity
endif

ifeq "$(filter -DMHD, $(SIMU_OPTION))" "-DMHD"
CPU_FILE    += MHD_GetCellCenteredBInPatch.cpp  MHD_InterpolateBField.cpp  MHD_AllocateElectricArray.cpp \
               MHD_Aux_Check_InterfaceB.cpp  MHD_FixUp_Electric.cpp  MHD_Aux_Check_DivergenceB.cpp \
               MHD_BoundaryCondition_Outflow.cpp  MHD_BoundaryCondition_Reflecting.cpp  MHD_BoundaryCondition_User.cpp \
               MHD_CopyPatchInterfaceBField.cpp MHD_Init_BField_ByFile.cpp

CPU_FILE    += CPU_Shared_ConstrainedTransport.cpp  CPU_Shared_RiemannSolver_HLLD.cpp

ifeq "$(findstring -DLOAD_BALANCE, $(SIMU_OPTION))" "-DLOAD_BALANCE"
CPU_FILE    += MHD_LB_EnsureBFieldConsistencyAfterRestrict.cpp  MHD_LB_AllocateElectricArray.cpp \
               MHD_LB_ResetBufferElectric.cpp  MHD_LB_Refine_GetCoarseFineInterfaceBField.cpp
endif # LOAD_BALANCE
endif # MHD


# ELBDM source files (included only inf "MODEL=ELBDM")
# -------------------------------------------------------------------------------
else ifeq "$(filter -DMODEL=ELBDM, $(SIMU_OPTION))" "-DMODEL=ELBDM"
GPU_FILE    += CUFLU_ELBDMSolver.cu

CPU_FILE    += CPU_ELBDMSolver.cpp  ELBDM_Init_ByFunction_AssignData.cpp  ELBDM_Init_ByFile_AssignData.cpp \
               ELBDM_GetTimeStep_Fluid.cpp  ELBDM_Flag_EngyDensity.cpp  ELBDM_UnwrapPhase.cpp \
               ELBDM_GetTimeStep_Phase.cpp  ELBDM_SetTaylor3Coeff.cpp

vpath %.cu     Model_ELBDM/GPU_ELBDM
vpath %.cpp    Model_ELBDM/CPU_ELBDM  Model_ELBDM

ifeq "$(filter -DGRAVITY, $(SIMU_OPTION))" "-DGRAVITY"
GPU_FILE    += CUPOT_ELBDMGravitySolver.cu

CPU_FILE    += CPU_ELBDMGravitySolver.cpp  ELBDM_GetTimeStep_Gravity.cpp

vpath %.cu     Model_ELBDM/GPU_ELBDMGravity
vpath %.cpp    Model_ELBDM/CPU_ELBDMGravity
endif

endif # MODEL


# self-gravity source files (included only if "GRAVITY" is turned on)
# ------------------------------------------------------------------------------------
ifeq "$(filter -DGRAVITY, $(SIMU_OPTION))" "-DGRAVITY"
GPU_FILE    += CUAPI_MemAllocate_PoissonGravity.cu  CUAPI_MemFree_PoissonGravity.cu \
               CUAPI_Asyn_PoissonGravitySolver.cu  CUAPI_SetConstMemory_ExtAccPot.cu \
               CUAPI_SendExtPotTable2GPU.cu

GPU_FILE    += CUPOT_PoissonSolver_SOR_10to14cube.cu  CUPOT_PoissonSolver_SOR_16to18cube.cu \
               CUPOT_PoissonSolver_MG.cu  CUPOT_ExtAcc_PointMass.cu  CUPOT_ExtPot_PointMass.cu \
               CUPOT_ExtPotSolver.cu  CUPOT_ExtPot_Tabular.cu

CPU_FILE    += CPU_PoissonGravitySolver.cpp  CPU_PoissonSolver_SOR.cpp  CPU_PoissonSolver_FFT.cpp \
               CPU_PoissonSolver_MG.cpp  CPU_ExtPotSolver.cpp  CPU_ExtPotSolver_BaseLevel.cpp

CPU_FILE    += Init_FFTW.cpp  Gra_Close.cpp  Gra_Prepare_Flu.cpp  Gra_Prepare_Pot.cpp  Gra_Prepare_Corner.cpp \
               Gra_AdvanceDt.cpp  Poi_Close.cpp  Poi_Prepare_Pot.cpp  Poi_Prepare_Rho.cpp \
               Output_PreparedPatch_Poisson.cpp  Init_MemAllocate_PoissonGravity.cpp \
               End_MemFree_PoissonGravity.cpp  Init_Set_Default_SOR_Parameter.cpp  Init_GreenFuncK.cpp \
               Init_Set_Default_MG_Parameter.cpp  Poi_GetAverageDensity.cpp  Poi_AddExtraMassForGravity.cpp \
               Poi_BoundaryCondition_Extrapolation.cpp  Gra_Prepare_USG.cpp  Poi_StorePotWithGhostZone.cpp \
               Init_ExtAccPot.cpp  End_ExtAccPot.cpp  CPU_ExtAcc_PointMass.cpp  CPU_ExtPot_PointMass.cpp \
               Poi_UserWorkBeforePoisson.cpp  Init_LoadExtPotTable.cpp  CPU_ExtPot_Tabular.cpp

vpath %.cu     SelfGravity/GPU_Poisson  SelfGravity/GPU_Gravity
vpath %.cpp    SelfGravity/CPU_Poisson  SelfGravity/CPU_Gravity  SelfGravity
endif # GRAVITY


# particle source files (included only if "PARTICLE" is turned on)
# ------------------------------------------------------------------------------------
ifeq "$(filter -DPARTICLE, $(SIMU_OPTION))" "-DPARTICLE"
GPU_FILE    +=

CPU_FILE    += Par_Init_ByFunction.cpp  Par_Output_TextFile.cpp  Par_FindHomePatch_UniformGrid.cpp \
               Par_Aux_Check_Particle.cpp  Par_PassParticle2Father.cpp  Par_CollectParticle2OneLevel.cpp \
               Par_MassAssignment.cpp  Par_UpdateParticle.cpp  Par_GetTimeStep_VelAcc.cpp \
               Par_PassParticle2Sibling.cpp  Par_CountParticleInDescendant.cpp  Par_Aux_GetConservedQuantity.cpp \
               Par_Aux_InitCheck.cpp  Par_Aux_Record_ParticleCount.cpp  Par_PassParticle2Son_MultiPatch.cpp \
               Par_Synchronize.cpp  Par_PredictPos.cpp  Par_Init_ByFile.cpp  Par_Init_Attribute.cpp \
               Par_AddParticleAfterInit.cpp  Par_PassParticle2Son_SinglePatch.cpp  Par_Sort.cpp

vpath %.cu     Particle/GPU
vpath %.cpp    Particle/CPU  Particle

ifeq "$(findstring -DLOAD_BALANCE, $(SIMU_OPTION))" "-DLOAD_BALANCE"
CPU_FILE    += Par_LB_SendParticleData.cpp  Par_LB_CollectParticle2OneLevel.cpp \
               Par_LB_CollectParticleFromRealPatch.cpp  Par_LB_RecordExchangeParticlePatchID.cpp \
               Par_LB_MapBuffer2RealPatch.cpp  Par_LB_ExchangeParticleBetweenPatch.cpp \
               Par_LB_Refine_SendParticle2Father.cpp

vpath %.cpp    Particle/LoadBalance
endif # LOAD_BALANCE

endif # PARTICLE


# parallelization source files (included only if "SERIAL" is turned off)
# ------------------------------------------------------------------------------------
ifeq "$(filter -DSERIAL, $(SIMU_OPTION))" ""
CPU_FILE    += Flu_AllocateFluxArray_Buffer.cpp

CPU_FILE    += Flag_Buffer.cpp  Refine_Buffer.cpp

CPU_FILE    += Buf_AllocateBufferPatch.cpp  Buf_AllocateBufferPatch_Base.cpp  Buf_GetBufferData.cpp \
               Buf_RecordExchangeDataPatchID.cpp  Buf_RecordExchangeFluxPatchID.cpp Buf_SortBoundaryPatch.cpp \
               Buf_RecordBoundaryFlag.cpp  Buf_RecordBoundaryPatch.cpp  Buf_RecordBoundaryPatch_Base.cpp \
               Buf_ResetBufferFlux.cpp

CPU_FILE    += MPI_ExchangeBoundaryFlag.cpp  MPI_ExchangeBufferPosition.cpp  MPI_ExchangeData.cpp \
               Init_MPI.cpp  MPI_Exit.cpp

CPU_FILE    += Output_BoundaryFlagList.cpp  Output_ExchangeDataPatchList.cpp  Output_ExchangeFluxPatchList.cpp \
               Output_ExchangePatchMap.cpp

CPU_FILE    += Aux_Record_BoundaryPatch.cpp

vpath %.cpp    Buffer  MPI
endif # !SERIAL


# load-balance source files (included only if "LOAD_BALANCE" is turned on)
# ------------------------------------------------------------------------------------
CPU_FILE    += LB_HilbertCurve.cpp  LB_Utility.cpp

ifeq "$(findstring -DLOAD_BALANCE, $(SIMU_OPTION))" "-DLOAD_BALANCE"
CPU_FILE    += LB_Init_LoadBalance.cpp  LB_AllocateBufferPatch_Sibling.cpp  LB_RecordOvelapMPIPatchID.cpp \
               LB_Output_LBIdx.cpp  LB_AllocateBufferPatch_Father.cpp  LB_FindFather.cpp  LB_SiblingSearch.cpp \
               LB_RecordExchangeDataPatchID.cpp  LB_GetBufferData.cpp  LB_AllocateFluxArray.cpp \
               LB_RecordExchangeRestrictDataPatchID.cpp  LB_GrandsonCheck.cpp  LB_ExchangeFlaggedBuffer.cpp \
               LB_Refine.cpp  LB_Refine_GetNewRealPatchList.cpp  LB_Refine_AllocateNewPatch.cpp \
               LB_FindSonNotHome.cpp  LB_Refine_AllocateBufferPatch_Sibling.cpp \
               LB_AllocateBufferPatch_Sibling_Base.cpp  LB_RecordExchangeFixUpDataPatchID.cpp \
               LB_EstimateWorkload_AllPatchGroup.cpp  LB_EstimateLoadImbalance.cpp  LB_SetCutPoint.cpp \
               LB_Init_ByFunction.cpp  LB_Init_Refine.cpp

endif # LOAD_BALANCE

vpath %.cpp    LoadBalance


# yt inline analysis source files (included only if "SUPPORT_LIBYT" is turned on)
# ------------------------------------------------------------------------------------
ifeq "$(filter -DSUPPORT_LIBYT, $(SIMU_OPTION))" "-DSUPPORT_LIBYT"
CPU_FILE    += YT_Init.cpp  YT_End.cpp  YT_SetParameter.cpp  YT_AddAllGrid.cpp  YT_Inline.cpp

vpath %.cpp    YT
endif # SUPPORT_LIBYT


# local source terms source files
# ------------------------------------------------------------------------------------
GPU_FILE    += CUAPI_Asyn_SrcSolver.cu  CUSRC_SrcSolver_IterateAllCells.cu  CUSRC_Src_Deleptonization.cu \
               CUSRC_Src_User_Template.cu

CPU_FILE    += CPU_SrcSolver.cpp  CPU_SrcSolver_IterateAllCells.cpp  CPU_Src_Deleptonization.cpp \
               CPU_Src_User_Template.cpp

CPU_FILE    += Src_AdvanceDt.cpp  Src_Prepare.cpp  Src_Close.cpp  Src_Init.cpp  Src_End.cpp \
               Src_WorkBeforeMajorFunc.cpp

vpath %.cu     SourceTerms  SourceTerms/User_Template  SourceTerms/Deleptonization
vpath %.cpp    SourceTerms  SourceTerms/User_Template  SourceTerms/Deleptonization


# Grackle source files (included only if "SUPPORT_GRACKLE" is turned on)
# ------------------------------------------------------------------------------------
ifeq "$(filter -DSUPPORT_GRACKLE, $(SIMU_OPTION))" "-DSUPPORT_GRACKLE"
CPU_FILE    += CPU_GrackleSolver.cpp

CPU_FILE    += Grackle_Init.cpp  Grackle_End.cpp  Init_MemAllocate_Grackle.cpp  End_MemFree_Grackle.cpp \
               Grackle_Prepare.cpp  Grackle_Close.cpp  Grackle_Init_FieldData.cpp  Grackle_AdvanceDt.cpp

vpath %.cpp    Grackle  Grackle/CPU_Grackle
endif # SUPPORT_GRACKLE


# star formation source files (included only if "STAR_FORMATION" is turned on)
# ------------------------------------------------------------------------------------
ifeq "$(filter -DSTAR_FORMATION, $(SIMU_OPTION))" "-DSTAR_FORMATION"
CPU_FILE    += SF_CreateStar.cpp  SF_CreateStar_AGORA.cpp

vpath %.cpp    StarFormation
endif # STAR_FORMATION


# feedback source files (included only if "FEEDBACK" is turned on)
# ------------------------------------------------------------------------------------
ifeq "$(filter -DFEEDBACK, $(SIMU_OPTION))" "-DFEEDBACK"
CPU_FILE    += FB_AdvanceDt.cpp  FB_Init.cpp  FB_End.cpp  FB_User_Template.cpp  FB_SNe.cpp

vpath %.cpp    Feedback  Feedback/User_Template  Feedback/SNe
endif # FEEDBACK


# test problem source files
# --> just compile all .cpp and .cu files under TestProblem/*/*/
# ------------------------------------------------------------------------------------
CPU_FILE    += $(notdir $(wildcard TestProblem/*/*/*.cpp))
GPU_FILE    += $(notdir $(wildcard TestProblem/*/*/*.cu))

VPATH := $(dir $(wildcard TestProblem/*/*/))



# rules and targets
#######################################################################################################

# object files
# -------------------------------------------------------------------------------
# add filename prefixes to distinguish CPU and GPU object files
PREFIX_CPU   := __cpu__
PREFIX_GPU   := __gpu__
OBJ_PATH     := Object
OBJ_CPU      := $(patsubst %.cpp, $(OBJ_PATH)/$(PREFIX_CPU)%.o, $(CPU_FILE))
ifeq "$(filter -DGPU, $(SIMU_OPTION))" "-DGPU"
OBJ_GPU      := $(patsubst %.cu,  $(OBJ_PATH)/$(PREFIX_GPU)%.o, $(GPU_FILE))
OBJ_GPU_LINK := $(OBJ_PATH)/gpu_link.o
endif


# libraries
# -------------------------------------------------------------------------------
ifeq "$(filter -DGPU, $(SIMU_OPTION))" "-DGPU"
LIB += -L$(CUDA_PATH)/lib64
#LIB += -L$(CUDA_PATH)/lib
LIB += -lcudart
endif

ifeq "$(filter -DLAOHU, $(SIMU_OPTION))" "-DLAOHU"
LIB += -L$(GPUID_PATH) -lgpudevmgr
endif

ifeq "$(filter -DGRAVITY, $(SIMU_OPTION))" "-DGRAVITY"
   LIB += -L$(FFTW_PATH)/lib
   ifeq "$(filter -DFLOAT8, $(SIMU_OPTION))" "-DFLOAT8"
      ifeq "$(filter -DSERIAL, $(SIMU_OPTION))" "-DSERIAL"
         LIB += -ldrfftw -ldfftw
      else
         LIB += -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
      endif
   else
      ifeq "$(filter -DSERIAL, $(SIMU_OPTION))" "-DSERIAL"
         LIB += -lsrfftw -lsfftw
      else
         LIB += -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
      endif
   endif
endif

ifeq "$(filter -DSUPPORT_GRACKLE, $(SIMU_OPTION))" "-DSUPPORT_GRACKLE"
LIB += -L$(GRACKLE_PATH)/lib -lgrackle
endif

ifeq "$(filter -DSUPPORT_HDF5, $(SIMU_OPTION))" "-DSUPPORT_HDF5"
LIB += -L$(HDF5_PATH)/lib -lhdf5
endif

ifeq "$(filter -DSUPPORT_GSL, $(SIMU_OPTION))" "-DSUPPORT_GSL"
LIB += -L$(GSL_PATH)/lib -lgsl -lgslcblas
endif

ifeq "$(filter -DSUPPORT_LIBYT, $(SIMU_OPTION))" "-DSUPPORT_LIBYT"
LIB += -L$(LIBYT_PATH)/lib -lyt
endif


# headers
# -------------------------------------------------------------------------------
INCLUDE := -I../include

ifeq "$(filter -DMODEL=HYDRO, $(SIMU_OPTION))" "-DMODEL=HYDRO"
INCLUDE += -IModel_Hydro/GPU_Hydro
endif

ifeq "$(filter -DSERIAL, $(SIMU_OPTION))" ""
INCLUDE += -I$(MPI_PATH)/include
endif

ifeq "$(filter -DGRAVITY, $(SIMU_OPTION))" "-DGRAVITY"
INCLUDE += -I$(FFTW_PATH)/include
endif

ifeq "$(filter -DSUPPORT_GRACKLE, $(SIMU_OPTION))" "-DSUPPORT_GRACKLE"
INCLUDE += -I$(GRACKLE_PATH)/include
endif

ifeq "$(filter -DSUPPORT_HDF5, $(SIMU_OPTION))" "-DSUPPORT_HDF5"
INCLUDE += -I$(HDF5_PATH)/include
endif

ifeq "$(filter -DSUPPORT_GSL, $(SIMU_OPTION))" "-DSUPPORT_GSL"
INCLUDE += -I$(GSL_PATH)/include
endif

ifeq "$(filter -DSUPPORT_LIBYT, $(SIMU_OPTION))" "-DSUPPORT_LIBYT"
INCLUDE += -I$(LIBYT_PATH)/include
endif


# CXX flags
# -------------------------------------------------------------------------------
# remove the OpenMP flag if OPENMP is disabled
ifeq "$(filter -DOPENMP, $(SIMU_OPTION))" ""
   OPENMPFLAG =
endif

# fixes compilation issues on Intel MPI
ifeq "$(filter -DSERIAL, $(SIMU_OPTION))" ""
   CXXFLAG += -DMPICH_IGNORE_CXX_SEEK
endif

COMMONFLAG := $(INCLUDE) $(SIMU_OPTION)
CXXFLAG    += $(COMMONFLAG) $(OPENMPFLAG)


# NVCC flags
# -------------------------------------------------------------------------------
# common flags
NVCCFLAG_COM := $(COMMONFLAG) -O3 # -Xcompiler #-Xopencc -OPT:Olimit=0

ifeq      "$(filter -DGPU_ARCH=FERMI,   $(SIMU_OPTION))" "-DGPU_ARCH=FERMI"
   NVCCFLAG_ARCH += -gencode arch=compute_20,code=\"compute_20,sm_20\"
#  NVCCFLAG_ARCH += -gencode arch=compute_20,code=\"compute_20,sm_21\"
else ifeq "$(filter -DGPU_ARCH=KEPLER,  $(SIMU_OPTION))" "-DGPU_ARCH=KEPLER"
   NVCCFLAG_ARCH += -gencode arch=compute_30,code=\"compute_30,sm_30\"
   NVCCFLAG_ARCH += -gencode arch=compute_35,code=\"compute_35,sm_35\"
   NVCCFLAG_ARCH += -gencode arch=compute_37,code=\"compute_37,sm_37\"
else ifeq "$(filter -DGPU_ARCH=MAXWELL, $(SIMU_OPTION))" "-DGPU_ARCH=MAXWELL"
   NVCCFLAG_ARCH += -gencode arch=compute_50,code=\"compute_50,sm_50\"
   NVCCFLAG_ARCH += -gencode arch=compute_52,code=\"compute_52,sm_52\"
else ifeq "$(filter -DGPU_ARCH=PASCAL,  $(SIMU_OPTION))" "-DGPU_ARCH=PASCAL"
   NVCCFLAG_ARCH += -gencode arch=compute_60,code=\"compute_60,sm_60\"
   NVCCFLAG_ARCH += -gencode arch=compute_61,code=\"compute_61,sm_61\"
else ifeq "$(filter -DGPU_ARCH=VOLTA,  $(SIMU_OPTION))" "-DGPU_ARCH=VOLTA"
   NVCCFLAG_ARCH += -gencode arch=compute_70,code=\"compute_70,sm_70\"
else ifeq "$(filter -DGPU_ARCH=TURING,  $(SIMU_OPTION))" "-DGPU_ARCH=TURING"
   NVCCFLAG_ARCH += -gencode arch=compute_75,code=\"compute_75,sm_75\"
else ifeq "$(filter -DGPU, $(SIMU_OPTION))" "-DGPU"
   $(error unknown GPU_ARCH (please set it in the Makefile))
endif

NVCCFLAG_COM += $(NVCCFLAG_ARCH)

ifeq "$(filter -DGAMER_DEBUG, $(SIMU_OPTION))" "-DGAMER_DEBUG"
   NVCCFLAG_COM += -g -Xptxas -v   # "-G" may cause the GPU Poisson solver to fail
endif

#NVCCFLAG_COM += -use_fast_math

# fluid solver flags
NVCCFLAG_FLU += -Xptxas -dlcm=ca -prec-div=false -ftz=true

ifeq      "$(filter -DGPU_ARCH=FERMI, $(SIMU_OPTION))" "-DGPU_ARCH=FERMI"
   ifeq "$(filter -DFLOAT8, $(SIMU_OPTION))" "-DFLOAT8"
#    NVCCFLAG_FLU += --maxrregcount=XX
   else
#    NVCCFLAG_FLU += --maxrregcount=XX
   endif
else ifeq "$(filter -DGPU_ARCH=KEPLER, $(SIMU_OPTION))" "-DGPU_ARCH=KEPLER"
   ifeq "$(filter -DFLOAT8, $(SIMU_OPTION))" "-DFLOAT8"
     NVCCFLAG_FLU += --maxrregcount=128
   else
     NVCCFLAG_FLU += --maxrregcount=64
   endif
else ifeq "$(filter -DGPU_ARCH=MAXWELL, $(SIMU_OPTION))" "-DGPU_ARCH=MAXWELL"
   ifeq "$(filter -DFLOAT8, $(SIMU_OPTION))" "-DFLOAT8"
     NVCCFLAG_FLU += --maxrregcount=192
   else
     NVCCFLAG_FLU += --maxrregcount=128
   endif
else ifeq "$(filter -DGPU_ARCH=PASCAL, $(SIMU_OPTION))" "-DGPU_ARCH=PASCAL"
   ifeq "$(filter -DFLOAT8, $(SIMU_OPTION))" "-DFLOAT8"
     NVCCFLAG_FLU += --maxrregcount=192
   else
     NVCCFLAG_FLU += --maxrregcount=128
   endif
else ifeq "$(filter -DGPU_ARCH=VOLTA, $(SIMU_OPTION))" "-DGPU_ARCH=VOLTA"
   ifeq "$(filter -DFLOAT8, $(SIMU_OPTION))" "-DFLOAT8"
     NVCCFLAG_FLU += --maxrregcount=192
   else
     NVCCFLAG_FLU += --maxrregcount=128
   endif
else ifeq "$(filter -DGPU_ARCH=TURING, $(SIMU_OPTION))" "-DGPU_ARCH=TURING"
   ifeq "$(filter -DFLOAT8, $(SIMU_OPTION))" "-DFLOAT8"
     NVCCFLAG_FLU += --maxrregcount=192
   else
     NVCCFLAG_FLU += --maxrregcount=128
   endif
endif

# Poisson/gravity solvers flags
NVCCFLAG_POT += -Xptxas -dlcm=ca


# remove extra whitespaces
# -------------------------------------------------------------------------------
CXX          := $(strip $(CXX))
CXXFLAG      := $(strip $(CXXFLAG))
OPENMPFLAG   := $(strip $(OPENMPFLAG))
LIB          := $(strip $(LIB))
NVCC         := $(strip $(NVCC))
NVCCFLAG_COM := $(strip $(NVCCFLAG_COM))
NVCCFLAG_FLU := $(strip $(NVCCFLAG_FLU))
NVCCFLAG_POT := $(strip $(NVCCFLAG_POT))


# implicit rules (do NOT modify the order of the following rules)
# -------------------------------------------------------------------------------
# output detailed compilation commands or not
ifeq "$(COMPILE_VERBOSE)" "0"
ECHO = @
else
ECHO =
endif

# GPU codes
ifeq "$(filter -DGPU, $(SIMU_OPTION))" "-DGPU"
$(OBJ_PATH)/$(PREFIX_GPU)CUFLU_%.o : CUFLU_%.cu
	@echo "Compiling $<"
	$(ECHO)$(NVCC) $(NVCCFLAG_COM) $(NVCCFLAG_FLU) -o $@ -dc $<

$(OBJ_PATH)/$(PREFIX_GPU)CUPOT_%.o : CUPOT_%.cu
	@echo "Compiling $<"
	$(ECHO)$(NVCC) $(NVCCFLAG_COM) $(NVCCFLAG_POT) -o $@ -dc $<

$(OBJ_PATH)/$(PREFIX_GPU)CUSRC_%.o : CUSRC_%.cu
	@echo "Compiling $<"
	$(ECHO)$(NVCC) $(NVCCFLAG_COM) -o $@ -dc $<

$(OBJ_PATH)/$(PREFIX_GPU)CUAPI_%.o : CUAPI_%.cu
	@echo "Compiling $<"
	$(ECHO)$(NVCC) $(NVCCFLAG_COM) -o $@ -dc $<

$(OBJ_PATH)/$(PREFIX_GPU)%.o : %.cu
	@echo "Compiling $<"
	$(ECHO)$(NVCC) $(NVCCFLAG_COM) -o $@ -dc $<
endif # GPU

# CPU codes
$(OBJ_PATH)/$(PREFIX_CPU)%.o : %.cpp
	@echo "Compiling $<"
	$(ECHO)$(CXX) $(CXXFLAG) -o $@ -c $<


# linking
# -------------------------------------------------------------------------------
$(EXECUTABLE) : $(OBJ_CPU) $(OBJ_GPU)
# GPU linker
ifeq "$(filter -DGPU, $(SIMU_OPTION))" "-DGPU"
	@echo "Linking GPU codes"
	$(ECHO)$(NVCC) -o $(OBJ_GPU_LINK) $(OBJ_GPU) $(NVCCFLAG_ARCH) -dlink
endif

# CPU linker
	@echo "Linking CPU codes"
ifeq "$(COMPILE_VERBOSE)" "1"
	$(CXX) -o $@ $^ $(OBJ_GPU_LINK) $(LIB) $(OPENMPFLAG)
	@printf "\nCompiling GAMER --> Successful!\n\n"; \
	cp $(EXECUTABLE) ../bin/
else
	@$(CXX) -o $@ $^ $(OBJ_GPU_LINK) $(LIB) $(OPENMPFLAG); \
	(if [ -e $@ ]; then \
		printf "\nCompiling GAMER --> Successful!\n\n"; \
		cp $(EXECUTABLE) ../bin/; \
	else \
		printf "\nCompiling GAMER --> Failed!\n\n"; \
	fi)
endif
	@rm -f ./*.linkinfo


# clean
# -------------------------------------------------------------------------------
.PHONY: clean
clean :
	@rm -f $(OBJ_PATH)/*
	@rm -f $(EXECUTABLE)
	@rm -f ./*.linkinfo



