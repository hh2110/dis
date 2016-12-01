!-----------------------------------------------------------------------
! MODULE DEFINING COMMON VARIABLES
!-----------------------------------------------------------------------
      MODULE VAR
      USE, INTRINSIC :: ISO_C_BINDING
      SAVE
! DEFINE FFTW3 VARIABLES
      TYPE(C_PTR) :: FORWARD_PLAN, BACKWARD_PLAN, data 
      COMPLEX(C_DOUBLE_COMPLEX), POINTER :: FFT_DATA(:,:,:)     
      INTEGER :: FFT_PLANS_CREATED      
! DEFINE COMMON PARAMETERS
      REAL :: TSTART,TFINISH,TRUN            !! SIMULATION START, FINISH AND RUN TIMES
      REAL :: FFTSTART,FFTFINISH,FFTRUN      !! FFT PLANNING START, FINISH AND RUN TIMES
      REAL :: BPQSTART,BPQFINISH,BPQRUN      !! BPQ START, FINISH AND RUN TIMES
      REAL :: EVOSTART,EVOFINISH,EVOAVG      !! SIMULATION START, FINISH AND AVERAGE TIMES
      INTEGER(C_INT) :: FFTP                 !! FFT PLAN TYPE
      INTEGER :: NOMP                        !! NO. OF OPENMP THREADS
      INTEGER :: OUTGE                       !! OUTPUT GRADIENT ENERGY (YES/NO)
      INTEGER :: OUTFE                       !! OUTPUT FAULT ENERGY (YES/NO)
      INTEGER :: OUTETAS                     !! OUTPUT ETAS (YES/NO)
      INTEGER :: INITBENCH                   !! BENCHMARK INITIALISATION SPEED
      INTEGER :: NX,NY,NZ                    !! SYSTEM CELL DIMENSIONS
      INTEGER :: NXYZ                        !! SYSTEM VOLUME
      INTEGER :: NP                          !! NO OF SLIP PLANES
      INTEGER :: NQ                          !! NO OF IN-PLANE SLIP DIR. 
      REAL, PARAMETER :: SQRT3= 1.732050808  !! SQUARE ROOT OF 3
      REAL, PARAMETER :: PI=3.14159265359    !! PI
      REAL, PARAMETER :: TWOPI=2*PI          !! 2xPI
! DEFINE COMMON VARIABLES FOR 
      REAL :: AXX(3),AXY(3),AXZ(3)           !! ORTHOGONAL CRYSTAL AXES X,Y,Z
      REAL :: ROT(3,3)                       !! CRYSTAL ROTATION MATRIX
      INTEGER :: ROTFLAG                     !! CRYSTAL ROTATION FLAG 1 OR 0       
      REAL :: NV(4,3)                        !! NORMAL VECTORS TO SLIP PLANES
      REAL :: BV(4,3,3)                      !! BURGERS VECTORS
      REAL,ALLOCATABLE :: G2(:,:,:)          !! K^2 CONSTANTS, LAPLACE OPERATOR IN K-SPACE
      REAL,ALLOCATABLE :: G(:,:,:,:)         !! K (IJK)COMPONENTS
      REAL,ALLOCATABLE :: BPQ(:,:,:,:,:,:,:) !! BPQ MATRIX
      REAL :: CIJKL(3,3,3,3)                 !! COMPLIANCE MATRIX
      REAL,ALLOCATABLE :: E0(:,:,:,:)        !! TRANSFORMATION STRAIN
      REAL,ALLOCATABLE :: S0(:,:,:,:)        !! TRANSFORMATION STRESS
      REAL :: S_APP(3,3)                     !! EXTERNALlY APPLIED STRESS
      COMPLEX,ALLOCATABLE :: ETAK(:,:,:,:,:) !! FIELD VARIABLES, K-SPACE
      REAL,ALLOCATABLE :: ETA(:,:,:,:,:)     !! FIELD VARIABLES, R-SPACE
      INTEGER :: OUT1                        !! NUMBER OF OUTPUT TIME-SLICES
      INTEGER :: OUTDIM                      !! OUTPUT 2D OR 3D DATA ?     
      INTEGER :: SLICE_AXIS                    !! AXIS ALONG WHICH 2D TIME SLICES AR EOUTPUT X,Y,Z=1,2,3      
      INTEGER :: SLICE_POS                   !! POSITION ALONG AXIS FOR TIME SLICE OUTPUT         
      REAL,ALLOCATABLE :: TPR(:,:,:)         !! TEMPORARY VARIABLES, K-SPACE
      COMPLEX,ALLOCATABLE :: TPK(:,:,:)      !! TEMPORARY VARIABLES, R-SPACE
      REAL,ALLOCATABLE :: PHI(:,:,:)         !! MICROSTRUCTURE VARIABLE
! DEFINE COMMON VARIABLES FOR INO.F
      REAL :: AA,A0                          !! LATTICE PARAMETER, REDUCED LATTICE PARAMETER
      REAL :: C11,C12,C44,C11X,C12X,C44X     !! ELASTIC CONSTANTS
      REAL :: S11,S12,S13,S21,S22,S23,S31,S32,S33  !! STRESS TENSOR COMPONENTS
      REAL :: BETA                           !! GRADIENT ENERGY COEFFICEINT BETA
      REAL :: T0, DT                         !! TOTAL TIME AND TIME STEP
      REAL :: S_APP_MAG                      !! APPLIED STRESS MAGNITUDE (REDUCED)
      INTEGER :: NIN                         !! INITIAL CONFIGURATION IDENTIFIER
      INTEGER :: NPHIS                       !! NUMBER OF PHASES - PHI      
      INTEGER :: PHID                        !! PHASE ID NUMBER
      INTEGER :: GSURF(9)                    !! GAMMA SURFACE DATA IDENTIFIER
      REAL,ALLOCATABLE :: GC(:,:)            !! GAMMA SURFACE FOURIER COEFFICIENTS
      INTEGER,ALLOCATABLE :: GSF_EQ(:)       !! GAMMA SURFACE EQUATION
      REAL :: DXD                            !! ELASTIC SCALING FACTOR
      INTEGER :: BORDER                      !! SIZE OF FROZEN EVOLUTION FRAME
      INTEGER :: NSS                         !! SLIP SYSTEM IDENTIFIER
      CHARACTER(LEN=25) :: FILENAME      
!      
      REAL,ALLOCATABLE :: DER(:,:,:,:,:,:)   !! FOR CALCULATING GRADIENT ENERGY
      END MODULE VAR