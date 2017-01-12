      PROGRAM DIS
      USE VAR 
      IMPLICIT NONE
! DEFINE VARIABLES     
      INTEGER :: STEPNO             !! TIMESTEP NUMBER 
      INTEGER :: NSTEPS             !! TOTAL NUMBER OF ITERATIONS 
      INTEGER :: SFREQ              !! OUTPUT SAMPLING FREQUENCY
      INTEGER :: ICOUNT             !! ITERATION COUNTER
      INTEGER :: IX,IY,IZ           !! POSITION VECTOR INDICES
! SET UP OPENMP STUFF     
      REAL(KIND=8) :: OMP_GET_WTIME !! OMP TIMER
! LOG THE SIMULATION START TIME      
      TSTART=OMP_GET_WTIME()   
! CALL SUBROUTINES AND PRINT KEY SIMULATION PARAMETERS
      PRINT *, '_|_|_|    _|_|_|    _|_|_|'
      PRINT *, '_|    _|    _|    _|' 
      PRINT *, '_|    _|    _|      _|_|'
      PRINT *, '_|    _|    _|          _|' 
      PRINT *, '_|_|_|    _|_|_|  _|_|_|'
      PRINT *, ''
      PRINT *, 'DISLOCATION INTERACTION SIMILATOR'
      PRINT *, 'VERSION 1.0 BETA'
      PRINT *, ''      
      PRINT *, 'BASED ON THE ORIGINAL PHASE FIELD MODEL BY DR. CHEN SHEN and PROF. YUNZHI WANG'
      PRINT *, 'IMPLEMENTING THE KHACHATURYAN-SHATALOV ELASTIC INCLUSION THEORY'
      PRINT *, ''      
      PRINT *, 'CREATED BY: VASSILI VORONTSOV and HIKMATYAR HASAN'
      PRINT *, 'EMAIL QUESTIONS, BUG REPORTS AND SUGGESTIONS TO: vassili.vorontsov@gmail.com'      
      PRINT *, ''      
      PRINT *, '<<BEGINNING RUN>>'
      PRINT *, ''        
      PRINT *, 'INITIALISING - CELL PARAMETERS FROM "PAR.INP"'
      CALL PAR
      PRINT *, ''         
      PRINT *, 'INITIALISING - SLIP SYSTEMS'
      CALL SLP      
      PRINT *, ''        
      PRINT *, 'INITIALISING - ALLOCATING ARRAYS'
      CALL ALO      
      PRINT *, ''        
      PRINT *, 'INITIALISING - FAST FOURIER TRANSFORM FFTW3'
      CALL FFT_CREATE_PLANS(FFTP)
! INITIALISE PARALLELISATION IN REMAINING CODE      
      CALL OMP_SET_NUM_THREADS(NOMP)
      PRINT *, ''        
      PRINT *, 'INITIALISING - GSF DATA'
      CALL GSF
      PRINT *, ''         
      PRINT *, 'INITIALISING - SLIP SYSTEMS'
      CALL SLP
      PRINT *, ''         
      PRINT *, 'INITIALISING - K^2 TABLE'     
      CALL GRD
      PRINT *, ''       
      PRINT *, 'INITIALISING - BPQ MATRIX'
      CALL BMA
      PRINT *, ''       
      PRINT *, 'INITIALISING - FIELD VARIABLES & MICROSTRUCTURE'
      CALL MIC
      PRINT *, ''       
      PRINT *, 'INITIALISING - APPLIED STRESS'
      CALL STR
      PRINT *, ''       
      PRINT *, 'INITIALISING - CALCULATING FILE SIZES'
      CALL BIG      
!
      NSTEPS=T0/DT
!
      IF(INITBENCH == 0) THEN
        PRINT *, ''       
        PRINT *, '<<BEGINNING SYSTEM EVOLUTION>>'
!        
        IF(OUT1 == 0) THEN
          SFREQ =0
        ELSE
          SFREQ =T0/DT/OUT1
        END IF
!        
        CALL OPN
!        
        ICOUNT =0
!        
        EVOSTART=OMP_GET_WTIME()         
        DO STEPNO=0,NSTEPS      
          CALL EVO(STEPNO)
!
          IF(ONDIFF == 0) THEN
          ELSE IF (ONDIFF == 1) THEN
            CALL DIF(STEPNO)
          ELSE
            STOP '>>ERROR - INVALID ONDIFF IN PAR.INP - ENTER 0-NO OR 1-YES'
          END IF
!          
          IF(ICOUNT >= SFREQ .AND. OUT1 /= 0) THEN
            PRINT '(" TIME=",F8.2)',STEPNO*DT        
!
            IF(OUTCONC == 0) THEN
            ELSE IF (OUTCONC == 1) THEN
              CALL OCE(STEPNO)
            ELSE
              STOP '>>ERROR - INVALID OUTCONC IN PAR.INP - ENTER 0-NO or 1-YES'
            END IF
!
            IF(OUTGE == 0) THEN
            ELSE IF (OUTGE == 1) THEN
              CALL OGE(STEPNO)
            ELSE
              STOP '>>ERROR - INVALID OUTGE IN PAR.INP - ENTER 0-NO or 1-YES'
            END IF
!        
            IF(OUTFE == 0) THEN
            ELSE IF (OUTFE == 1) THEN
              CALL OFE(STEPNO)
            ELSE
              STOP '>>ERROR - INVALID OUTFE IN PAR.INP - ENTER 0-NO or 1-YES'
            END IF
!        
            IF(OUTETAS == 0) THEN
            ELSE IF (OUTETAS == 1) THEN
              CALL OET(STEPNO)
            ELSE
              STOP '>>ERROR - INVALID OUTETAS IN PAR.INP - ENTER 0-NO or 1-YES'
            END IF         
            ICOUNT =1
          ELSE
            ICOUNT =ICOUNT+1
          END IF
        END DO
!        
        EVOFINISH=OMP_GET_WTIME()
        EVOAVG=(EVOFINISH-EVOSTART)/(NSTEPS+1)
        IF(EVOAVG < 60.0) THEN
          PRINT '(" >>AVERAGE TIME PER ITERATION = ",F8.5," SECONDS.")',EVOAVG
        ELSE IF(EVOAVG >= 60.0 .AND. EVOAVG<3600.0) THEN
          PRINT '(" >>AVERAGE TIME PER ITERATION = ",F8.5," MINUTES.")',EVOAVG/60.0
        ELSE IF(EVOAVG>=3600.0) THEN
          PRINT '(" >>AVERAGE TIME PER ITERATION = ",ES14.5," HOURS.")',EVOAVG/3600.0
        END IF
!          
        PRINT *, '<<GENERATING OUTPUT FILES>>'
        CALL CLS
        CALL FET
        ELSE IF(INITBENCH==1) THEN
        PRINT *, '<<INITIALISATION BENCHMARK ONLY>>'
        PRINT *, '<<NO DATA WILL BE COMPUTED>>' 
        END IF
        CALL FFT_DESTROY_PLANS
!
        PRINT *, '<<PROGRAM FINISHED>>'
! CALCULATE SIMULATION TIME   
      TFINISH=OMP_GET_WTIME()
      TRUN=TFINISH-TSTART
      IF(TRUN < 60.0) THEN
      PRINT '(" >>CALULATION TIME = ",F8.5," SECONDS.")',TFINISH-TSTART
      ELSE IF(TRUN >= 60.0 .AND. TRUN<3600.0) THEN
      PRINT '(" >>CALULATION TIME = ",F8.5," MINUTES.")',(TFINISH-TSTART)/60.0
      ELSE IF(TRUN>=3600.0) THEN
      PRINT '(" >>CALULATION TIME = ",ES14.5," HOURS.")',(TFINISH-TSTART)/3600.0
      END IF
      END PROGRAM DIS
