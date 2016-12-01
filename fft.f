!-----------------------------------------------------------------------
! SUBROUTINE TO INITIALISE FFTW3 LIBRARY
!-----------------------------------------------------------------------
      SUBROUTINE FFT_CREATE_PLANS(PLAN)
      USE VAR
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE
      INCLUDE 'fftw3.f03'
      INTEGER :: PLAN, VOID
      INTEGER(C_INT) :: RET
      REAL(KIND=8) :: OMP_GET_WTIME
! LOG THE TIME      
      FFTSTART=OMP_GET_WTIME()       
! CHECK IF FFT IS ALREADY INITIALISED
      IF(FFT_PLANS_CREATED == 1) THEN
        PRINT *, '>>WARNING - FFTW PLANS ARE ALREADY INITIALISED!'
      END IF
! ALLOCATE SYSTEM SIZE     
      data = fftw_alloc_complex(int(NX*NY*NZ,C_SIZE_T))
      CALL c_f_pointer(data,FFT_DATA,[NX,NY,NZ])
! INITIALISE THREADING
      VOID = fftw_init_threads()
      IF (VOID==1) THEN
        PRINT *,'>>NOTE - FFTW THREADS INITIALISED SUCCESSFULLY'     
      ELSE
        PRINT *,''
        PRINT *,'>>WARNING - FFTW THREADS FAILED TO INITIALISE!'     
        PRINT *,''
      END IF
      CALL OMP_SET_NUM_THREADS(NOMP)
      CALL fftw_plan_with_nthreads(NOMP)     
! INITIALISE FFT ACCORDING TO FFT PLAN'
      IF(PLAN == 1) THEN
         FORWARD_PLAN = fftw_plan_dft_3d(NZ, NY, NX, FFT_DATA, FFT_DATA, FFTW_FORWARD, FFTW_ESTIMATE)
         BACKWARD_PLAN = fftw_plan_dft_3d(NZ, NY, NX, FFT_DATA, FFT_DATA, FFTW_BACKWARD, FFTW_ESTIMATE)     
         PRINT *,'>>NOTE - USING ESTIMATE FFT PLAN'       
      ELSE IF(PLAN == 2) THEN
         FORWARD_PLAN = fftw_plan_dft_3d(NZ, NY, NX, FFT_DATA, FFT_DATA, FFTW_FORWARD, FFTW_MEASURE)
         BACKWARD_PLAN = fftw_plan_dft_3d(NZ, NY, NX, FFT_DATA, FFT_DATA, FFTW_BACKWARD, FFTW_MEASURE)     
         PRINT *,'>>NOTE - USING MEASURE FFT PLAN'         
      ELSE IF(PLAN == 3) THEN
         FORWARD_PLAN = fftw_plan_dft_3d(NZ, NY, NX, FFT_DATA, FFT_DATA, FFTW_FORWARD, FFTW_PATIENT)    
         RET = fftw_export_wisdom_to_filename(C_CHAR_'wisdom_fwd.dis' // C_NULL_CHAR)
         IF (RET .EQ. 0) STOP '>>ERROR - FAILED EXPORTING FORWARD PLAN WISDOM TO FILE' 
!         
         BACKWARD_PLAN = fftw_plan_dft_3d(NZ, NY, NX, FFT_DATA, FFT_DATA, FFTW_BACKWARD, FFTW_PATIENT)      
         RET = fftw_export_wisdom_to_filename(C_CHAR_'wisdom_bkwd.dis' // C_NULL_CHAR)
         IF (RET .EQ. 0) STOP '>>ERROR - FAILED EXPORTING BACKWARD PLAN WISDOM TO FILE'
!              
         PRINT *,'>>NOTE - USING PATIENT FFT PLAN'       
      ELSE IF(PLAN == 4) THEN
         FORWARD_PLAN = fftw_plan_dft_3d(NZ, NY, NX, FFT_DATA, FFT_DATA, FFTW_FORWARD, FFTW_EXHAUSTIVE)
         RET = fftw_export_wisdom_to_filename(C_CHAR_'wisdom_fwd.dis' // C_NULL_CHAR)
         IF (RET .EQ. 0) STOP '>>ERROR - FAILED EXPORTING FORWARD PLAN WISDOM TO FILE' 
!          
         BACKWARD_PLAN = fftw_plan_dft_3d(NZ, NY, NX, FFT_DATA, FFT_DATA, FFTW_BACKWARD, FFTW_EXHAUSTIVE)
         RET = fftw_export_wisdom_to_filename(C_CHAR_'wisdom_bkwd.dis' // C_NULL_CHAR)
         IF (RET .EQ. 0) STOP '>>ERROR - FAILED EXPORTING BACKWARD PLAN WISDOM TO FILE'  
!               
         PRINT *,'>>NOTE - USING EXHAUSTIVE FFT PLAN'
      ELSE IF(PLAN == 5) THEN
         RET = fftw_import_wisdom_from_filename(C_CHAR_'wisdom_fwd.dis' // C_NULL_CHAR)
         IF (RET .EQ. 0) STOP '>>ERROR - FAILED IMPORTING FORWARD PLAN WISDOM TO FILE' 
         FORWARD_PLAN = fftw_plan_dft_3d(NZ, NY, NX, FFT_DATA, FFT_DATA, FFTW_FORWARD, FFTW_WISDOM_ONLY)
!         
         RET = fftw_import_wisdom_from_filename(C_CHAR_'wisdom_bkwd.dis' // C_NULL_CHAR)
         IF (RET .EQ. 0) STOP '>>ERROR - FAILED IMPORTING BACKWARD PLAN WISDOM TO FILE' 
         BACKWARD_PLAN = fftw_plan_dft_3d(NZ, NY, NX, FFT_DATA, FFT_DATA, FFTW_BACKWARD, FFTW_WISDOM_ONLY)
!              
         PRINT *,'>>NOTE - USING FFTW WISDOM'     
      ELSE
      STOP '>>ERROR - VALID OPTION (1-ESTIMATE,2-MEASURE,3-PATIENT,4-EXHAUSTIVE,5-WISDOM)!'
      END IF
!      
      FFT_PLANS_CREATED = 1
      PRINT *,'>>NOTE - FFTW PLANS CREATED SUCCESSFULLY'
! CALCULATE FFT PLANNING TIME   
      FFTFINISH=OMP_GET_WTIME()
      FFTRUN=FFTFINISH-FFTSTART
      IF(FFTRUN < 60.0) THEN
      PRINT '(" >>NOTE - FFT PLANNING TIME = ",F8.5," SECONDS.")',FFTRUN
      ELSE IF(FFTRUN >= 60.0 .AND. FFTRUN<3600.0) THEN
      PRINT '(" >>NOTE - FFT PLANNING TIME = ",F8.5," MINUTES.")',FFTRUN/60.0
      ELSE IF(FFTRUN>=3600.0) THEN
      PRINT '(" >>NOTE - FFT PLANNING TIME = ",ES14.5," HOURS.")',FFTRUN/3600.0
      END IF
!       
      END SUBROUTINE FFT_CREATE_PLANS
!-----------------------------------------------------------------------
! SUBROUTINE TO DESTROY FFTW3 PLANS
!-----------------------------------------------------------------------
      SUBROUTINE FFT_DESTROY_PLANS
      USE VAR
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE
      INCLUDE 'fftw3.f03'
! CHECK IF FFT IS INITIALISED
      IF(FFT_PLANS_CREATED /= 1) THEN
         STOP '>>ERROR - FFTW PLAN CREATION HAS FAILED!'
      END IF
      CALL fftw_destroy_plan(FORWARD_PLAN)
      CALL fftw_destroy_plan(BACKWARD_PLAN)
      FFT_PLANS_CREATED = 0
      END SUBROUTINE FFT_DESTROY_PLANS
!-----------------------------------------------------------------------
!-  FORWARD FFT , REAL TO COMPLEX
!-----------------------------------------------------------------------
      SUBROUTINE FFT_R2C(IN_DATA, OUT_DATA)
      USE VAR
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE
      INCLUDE 'fftw3.f03'
! DEFINE VARIABLES
      REAL :: IN_DATA(NX,NY,NZ)
      COMPLEX :: OUT_DATA(NX,NY,NZ)
      INTEGER :: I,J,K
! CHECK FFT IS INITIALISED
      IF(FFT_PLANS_CREATED /= 1) THEN
         STOP '>>ERROR - FFTW PLAN CREATION HAS FAILED!'
      END IF
! CONVERT INPUT ARRAY TO DOUBLE COMPLEX TYPE
      DO K=1,NZ
      DO J=1,NY
      DO I=1,NX
        FFT_DATA(I,J,K) =IN_DATA(I,J,K)
      END DO
      END DO 
      END DO
      CALL fftw_execute_dft(FORWARD_PLAN,FFT_DATA,FFT_DATA)
! COPY TO OUT_DATA
      DO K=1,NZ
      DO J=1,NY
      DO I=1,NX
        OUT_DATA(I,J,K) =CMPLX(FFT_DATA(I,J,K))
      END DO
      END DO
      END DO
      END SUBROUTINE FFT_R2C
!-----------------------------------------------------------------------
!-   INVERSE FFT , COMPLEX TO REAL
!-----------------------------------------------------------------------
      SUBROUTINE INVFFT_C2R(IN_DATA, OUT_DATA)
      USE VAR
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE
      INCLUDE 'fftw3.f03'
! DEFINE VARIABLES
      COMPLEX :: IN_DATA(NX,NY,NZ)
      REAL :: OUT_DATA(NX,NY,NZ)
      INTEGER :: I,J,K
! CHECK FFT IS INITIALISED
      IF(FFT_PLANS_CREATED /= 1) THEN
         STOP '>>ERROR - FFTW PLAN CREATION HAS FAILED!'
      END IF
! CONVERT INPUT ARRAY TO DOUBLE COMPLEX TYPE
      DO K=1,NZ
      DO J=1,NY
      DO I=1,NX
        FFT_DATA(I,J,K) =IN_DATA(I,J,K)
      END DO
      END DO 
      END DO
      CALL fftw_execute_dft(BACKWARD_PLAN,FFT_DATA,FFT_DATA)
! COPY TO OUT_DATA
      DO K=1,NZ
      DO J=1,NY
      DO I=1,NX
        OUT_DATA(I,J,K) =FFT_DATA(I,J,K)/NXYZ
      END DO
      END DO
      END DO
      END SUBROUTINE INVFFT_C2R