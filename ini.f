!-----------------------------------------------------------------------
! SUBROUTINE DEFINING INITIAL MICROSTRUCTURE (VARIABLE FIELDS)
!-----------------------------------------------------------------------
      SUBROUTINE MIC
      USE VAR 
      IMPLICIT NONE
! DEFINE VARIABLES
      INTEGER :: SA, SMA
      INTEGER :: IX,IY,IZ
      INTEGER :: OX,OY,OZ
      INTEGER :: CFG(NX,NZ)
      INTEGER :: XPOS
      INTEGER :: UNITCOUNT
! GAMMA-GAMMA PRIME SYSTEM DEFINITION
!-----------------------------------------------------------------------
! CONDITION-1 EDGE CONFIGURATION a/2<110> SPHERICAL PRECIPITATE
!-----------------------------------------------------------------------
      IF(NIN  ==  1) THEN
      PRINT *,'>>NOTE - EDGE a/2[10-1] vs SPHERICAL PRECIPITATE'      
! DEFINE INITIAL VARIABLE FIELDS - ETA(1,2)
!$OMP PARALLEL DO PRIVATE(IY,IX)
        DO IZ=1,NZ
        DO IY=NY/2,NY/2
        DO IX=1,NX
          IF(IX >= 0.4*NX .AND. IX <= 0.6*NX) THEN
            ETA(1,2,IX,IY,IZ)=1.0
          END IF
        END DO
        END DO 
        END DO
!$OMP END PARALLEL DO        
!        
! DETERMINE PRECIPITATE AREA
      XPOS=0.4*NX-NZ/4
!$OMP PARALLEL DO PRIVATE(IY,IX)
        DO IZ=1,NZ
        DO IY=1,NY
        DO IX=XPOS,NX
          IF(IX <= ((NZ/4)**2-(IZ-NZ/2)**2)**0.5+(0.4*NX-NZ/4-3))THEN
          PHI(IX,IY,IZ)=1.0
          END IF        
        END DO
        END DO
        END DO
!$OMP END PARALLEL DO         
!
!$OMP PARALLEL DO PRIVATE(IY,IX)
        DO IZ=1,NZ
        DO IY=1,NY
        DO IX=1,XPOS-1
          IF(IX >= (-((NZ/4)**2-(IZ-NZ/2)**2)**0.5)+(0.4*NX-NZ/4-3))THEN
          PHI(IX,IY,IZ)=1.0
          END IF        
        END DO
        END DO
        END DO
!$OMP END PARALLEL DO         
!        
!-----------------------------------------------------------------------
! CONDITION-1 EDGE CONFIGURATION a/2<110> IN SINGLE PHASE
!-----------------------------------------------------------------------
      ELSE IF(NIN  ==  2) THEN
      PRINT *,'>>NOTE - EDGE a/2[10-1] vs SPHERICAL PRECIPITATE'      
! DEFINE INITIAL VARIABLE FIELDS - ETA(1,2)
        DO IZ=1,NZ
        DO IY=NY/2,NY/2
        DO IX=1,NX
          IF(IX >= 0.4*NX .AND. IX <= 0.6*NX) THEN
            ETA(1,2,IX,IY,IZ)=2.0
          END IF
        END DO
        END DO 
        END DO
!        
! DETERMINE MATRIX AREA
        DO IZ=1,NZ
        DO IY=1,NY
        DO IX=1,NX
          PHI(IX,IY,IZ)=0.0
        END DO
        END DO
        END DO         
!-----------------------------------------------------------------------
! CONDITION-2 SCREW CONFIGURATION a/2<110> 180х660
!-----------------------------------------------------------------------
!       ELSE IF(NIN  ==  2) THEN
!       PRINT *,'>>NOTE - NME THIS CONFIG'
! C DEFINE INITIAL VARIABLE FIELDS - ETA(1,2)
!         DO IZ=1,NZ
!         DO IY=NY/2,NY/2
!         DO IX=1,NX
!           IF(IZ >= 270 .AND. IZ <= 390) THEN
!             ETA(1,2,IX,IY,IZ)=1.0
!           END IF
!         END DO
!         END DO
!         END DO
!     DETERMINE PRECIPITATE AREA
!         DO IX=1,NX
!         DO IY=1,NY
!         DO IZ=212,NZ
!           IF(IZ <= (2496-(IX-90)**2)**0.5+212)THEN
!           PHI(IX,IY,IZ)=1.0
!           END IF        
!         END DO
!         END DO
!         END DO
!         DO IX=1,NX
!         DO IY=1,NY
!         DO IZ=1,211
!           IF(IZ >= (-(2496-(IX-90)**2)**0.5)+212)THEN
!           PHI(IX,IY,IZ)=1.0
!           END IF        
!         END DO
!         END DO
!         END DO    
!-----------------------------------------------------------------------
! CONDITION-11 READ FIELD VARIABLES FROM FILES NXхNZ
!-----------------------------------------------------------------------  
      ELSE IF(NIN  ==  99) THEN
       PRINT *,'>>NOTE - READING FIELD VARIABLES FROM FILES'
!
! LOAD 2D MICROSTRUCTURE AND ETAS FROM FILE
        IF (OUTDIM == 2) THEN 
! SLICE ALONG X-AXIS AT SLICE_POS
          IF (SLICE_AXIS == 1) THEN
!          
            UNITCOUNT=91
!            
            DO IX=1,NX
              WRITE (FILENAME,'(A7)')'phi.dis'              
              OPEN (UNIT = UNITCOUNT, FILE=FILENAME, FORM = 'FORMATTED')
              READ (UNITCOUNT,*) ((PHI(IX,IY,IZ),IZ=1,NZ),IY=1,NY)
              CLOSE(UNITCOUNT)              
            END DO
!              
            UNITCOUNT=91
!
            DO SA=1,NP
            DO SMA=1,NQ
              WRITE (FILENAME,'(A3,I1,I1,A10)')'eta',SA,SMA,'_final.dis'
              OPEN (UNIT=UNITCOUNT, FILE=FILENAME, FORM = 'FORMATTED')        
              READ(UNITCOUNT,*) ((ETA(SA,SMA,SLICE_POS,IY,IZ),IZ=1,NZ),IY=1,NY)
              CLOSE(UNITCOUNT)
              UNITCOUNT=UNITCOUNT+1
            END DO
            END DO
! SLICE ALONG Y-AXIS AT SLICE_POS            
          ELSE IF (SLICE_AXIS == 2) THEN
!          
            UNITCOUNT=91
!            
            DO IY=1,NY
              WRITE (FILENAME,'(A7)')'phi.dis'                
              OPEN (UNIT = UNITCOUNT, FILE=FILENAME, FORM = 'FORMATTED')
              READ (UNITCOUNT,*) ((PHI(IX,IY,IZ),IZ=1,NZ),IX=1,NX)
              CLOSE(UNITCOUNT)              
            END DO
!              
            UNITCOUNT=91
!            
            DO SA=1,NP
            DO SMA=1,NQ
              WRITE (FILENAME,'(A3,I1,I1,A10)')'eta',SA,SMA,'_final.dis'
              OPEN (UNIT=UNITCOUNT, FILE=FILENAME, FORM = 'FORMATTED')        
              READ(UNITCOUNT,*) ((ETA(SA,SMA,IX,SLICE_POS,IZ),IZ=1,NZ),IX=1,NX)
              CLOSE(UNITCOUNT)
              UNITCOUNT=UNITCOUNT+1
            END DO
            END DO
! SLICE ALONG Z-AXIS AT SLICE_POS            
!          
            UNITCOUNT=91
!            
            DO IZ=1,NZ
              WRITE (FILENAME,'(A7)')'phi.dis'                 
              OPEN (UNIT = UNITCOUNT, FILE=FILENAME, FORM = 'FORMATTED')
              READ (UNITCOUNT,*) ((PHI(IX,IY,IZ),IY=1,NY),IX=1,NX)
              CLOSE(UNITCOUNT)              
            END DO
!              
            UNITCOUNT=91
!
            DO SA=1,NP
            DO SMA=1,NQ
              WRITE (FILENAME,'(A3,I1,I1,A10)')'eta',SA,SMA,'_final.dis'
              OPEN (UNIT=UNITCOUNT, FILE=FILENAME, FORM = 'FORMATTED')        
              READ(UNITCOUNT,*) ((ETA(SA,SMA,IX,IY,SLICE_POS),IY=1,NY),IX=1,NX)
              CLOSE(UNITCOUNT)
              UNITCOUNT=UNITCOUNT+1
            END DO
            END DO
          END IF   
!
! LOAD 3D MICROSTRUCTURE FROM FILE 
        ELSE IF (OUTDIM == 3) THEN 
!
          UNITCOUNT=91
!          
            WRITE (FILENAME,'(A7)')'phi.dis'                 
            OPEN (UNIT = UNITCOUNT, FILE=FILENAME, FORM = 'FORMATTED')
            READ (UNITCOUNT,*) (((PHI(IX,IY,IZ),IZ=1,NZ),IY=1,NY),IX=1,NX)
            CLOSE(UNITCOUNT)  
! LOAD 3D ETAS FROM FILE 
          UNITCOUNT=91
!
          DO SA=1,NP
          DO SMA=1,NQ
            WRITE (FILENAME,'(A3,I1,I1,A10)')'eta',SA,SMA,'_final.dis'
            OPEN (UNIT=UNITCOUNT, FILE=FILENAME, FORM = 'FORMATTED')        
            READ(UNITCOUNT,*) (((ETA(SA,SMA,IX,IY,IZ),IZ=1,NZ),IY=1,NY),IX=1,NX)
            CLOSE(UNITCOUNT)
            UNITCOUNT=UNITCOUNT+1
          END DO
          END DO
        END IF  
!
!-----------------------------------------------------------------------
! CONDITION-ERROR MESSAGE
!----------------------------------------------------------------------- 
      ELSE
        STOP '>>ERROR - INVALID INITIAL CONDITION !'
      END IF
!
!-----------------------------------------------------------------------
! OUTPUT MICROSTRUCTURE GEOMETRY TO A FILE
!-----------------------------------------------------------------------       
!
! WRITE 2D MICROSTRUCTURE TO FILE
      IF (OUTDIM == 2) THEN
        IF (SLICE_AXIS == 1) THEN
! SLICE ALONG X-AXIS AT SLICE_POS                  
          UNITCOUNT=91            
          WRITE (FILENAME,'(A7)')'phi.dis'                 
          OPEN (UNIT = UNITCOUNT, FILE=FILENAME, FORM = 'FORMATTED')
          WRITE (UNITCOUNT,*) ((PHI(IX,IY,IZ),IZ=1,NZ),IY=1,NY)
          CLOSE(UNITCOUNT)              
        ELSE IF (SLICE_AXIS == 2) THEN
! SLICE ALONG Y-AXIS AT SLICE_POS        
          UNITCOUNT=91            
          WRITE (FILENAME,'(A7)')'phi.dis'                 
          OPEN (UNIT = UNITCOUNT, FILE=FILENAME, FORM = 'FORMATTED')
          WRITE (UNITCOUNT,*) ((PHI(IX,IY,IZ),IZ=1,NZ),IX=1,NX)
          CLOSE(UNITCOUNT)              
        ELSE IF (SLICE_AXIS == 3) THEN
! SLICE ALONG Z-AXIS AT SLICE_POS        
          UNITCOUNT=91            
          WRITE (FILENAME,'(A7)')'phi.dis'                 
          OPEN (UNIT = UNITCOUNT, FILE=FILENAME, FORM = 'FORMATTED')
          WRITE (UNITCOUNT,*) ((PHI(IX,IY,IZ),IY=1,NY),IX=1,NX)
          CLOSE(UNITCOUNT)                       
        END IF
! WRITE 3D MICROSTRUCTURE TO FILE        
      ELSE IF (OUTDIM == 3) THEN
        ELSE IF (OUTDIM == 3) THEN 
!
          UNITCOUNT=91
!          
          WRITE (FILENAME,'(A7)')'phi.dis'                 
          OPEN (UNIT = UNITCOUNT, FILE=FILENAME, FORM = 'FORMATTED')
          WRITE(UNITCOUNT,*) (((PHI(IX,IY,IZ),IZ=1,NZ),IY=1,NY),IX=1,NX)
          CLOSE(UNITCOUNT)    
      END IF
      END SUBROUTINE MIC