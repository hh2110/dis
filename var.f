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
      REAL(KIND=8) :: TSTART,TFINISH,TRUN            !! SIMULATION START, FINISH AND RUN TIMES
      REAL(KIND=8) :: FFTSTART,FFTFINISH,FFTRUN      !! FFT PLANNING START, FINISH AND RUN TIMES
      REAL(KIND=8) :: BPQSTART,BPQFINISH,BPQRUN      !! BPQ START, FINISH AND RUN TIMES
      REAL(KIND=8) :: EVOSTART,EVOFINISH,EVOAVG      !! EVOLUTION START, FINISH AND AVERAGE TIMES
      REAL(KIND=8) :: DIFSTART,DIFFINISH,DIFAVG      !! CONCEVOLUTION START, FINISH AND AVERAGE TIMES
      INTEGER(C_INT) :: FFTP                 !! FFT PLAN TYPE
      INTEGER :: NOMP                        !! NO. OF OPENMP THREADS
      INTEGER :: OUTGE                       !! OUTPUT GRADIENT ENERGY (YES/NO)
      INTEGER :: OUTFE                       !! OUTPUT FAULT ENERGY (YES/NO)
      INTEGER :: OUTETAS                     !! OUTPUT ETAS (YES/NO)
      INTEGER :: OUTCONC                     !! OUTPUT CONC (YES/NO)
      INTEGER :: OUTSTRA                     !! OUTPUT STRA (YES/NO)
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
      REAL,ALLOCATABLE :: SE0(:,:,:,:,:)     !! SYSTEM TRANSFORMATION STRAIN
      REAL,ALLOCATABLE :: S0(:,:,:,:)        !! TRANSFORMATION STRESS
      REAL :: S_APP(3,3)                     !! EXTERNALlY APPLIED STRESS
      COMPLEX,ALLOCATABLE :: ETAK(:,:,:,:,:) !! FIELD VARIABLES, K-SPACE
      REAL,ALLOCATABLE :: ETA(:,:,:,:,:)     !! FIELD VARIABLES, R-SPACE
      INTEGER :: OUT1                        !! NUMBER OF OUTPUT TIME-SLICES
      INTEGER :: OUTDIM                      !! OUTPUT 2D OR 3D DATA ?     
      INTEGER :: SLICE_AXIS                  !! AXIS ALONG WHICH 2D TIME SLICES AR EOUTPUT X,Y,Z=1,2,3      
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
      INTEGER :: BORDERSTARTX, BORDERENDX     !! SIZE OF FROZEN EVOLUTION FRAME IN X
      INTEGER :: BORDERSTARTY, BORDERENDY     !! SIZE OF FROZEN EVOLUTION FRAME IN Y
      INTEGER :: BORDERSTARTZ, BORDERENDZ     !! SIZE OF FROZEN EVOLUTION FRAME IN Z
      INTEGER :: NSS                         !! SLIP SYSTEM IDENTIFIER
      CHARACTER(LEN=25) :: FILENAME      
!      
      REAL,ALLOCATABLE :: DER(:,:,:,:,:,:)   !! FOR CALCULATING GRADIENT ENERGY
!
!
! DEFINE COMMON VARIABLES FOR DIF.F
      INTEGER :: ONDIFF                      !! INLCUDE DIFFUSION (1/0?)
! VARIABLES
      REAL,ALLOCATABLE :: CONC(:,:,:)        !! CONC VARIABLE, R-SPACE
      REAL,ALLOCATABLE :: ONE(:,:,:), TWO(:,:,:), THREE(:,:,:) !!SEE NOTES
      REAL,ALLOCATABLE :: GDIFF(:,:,:)
      REAL,ALLOCATABLE :: C0(:,:,:)           !! position of chemical BFE curve
      REAL,ALLOCATABLE :: CONCFE(:,:,:)
      REAL,ALLOCATABLE :: ALPDG(:,:,:),FBENERGY(:,:,:), CONCGE(:,:,:)
      REAL :: DELTAGSISF          !!difference in SISFE
! PARAMETERS
      REAL, PARAMETER :: CG = 0.02    !!equilibrium concentration Gamma
      REAL, PARAMETER :: CGP = 0.3    !!equilibrium concentration GPrime
      REAL, PARAMETER :: CS = 0.4    !!equilibrium concentration SISF
      REAL, PARAMETER :: A1=0.5*(CS+CGP)     !!parameter for alpha...
      REAL, PARAMETER :: A2=0.02            !!...function of conc
      REAL, PARAMETER :: BETA1=1.0          !! BFE param
      REAL, PARAMETER :: MG = 0.05      !!concTIME EVOLUTION CONSTANT...
      REAL, PARAMETER :: EG = 1.0           !!gradient CONSTANT for conc
      REAL, PARAMETER :: DDX=1          !!DELTA_X FOR CONC EVO
      REAL, PARAMETER :: DELTA = 0.001     !!tolerance in conc evol
      REAL, PARAMETER :: IRAND = 0.005     !!initial randomisation
!      
      END MODULE VAR
