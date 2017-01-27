!-----------------------------------------------------------------------
! SUBROUTINE DEFINING INITIAL MICROSTRUCTURE (VARIABLE FIELDS)
!-----------------------------------------------------------------------
      SUBROUTINE MIC
      USE VAR 
      IMPLICIT NONE
! DEFINE VARIABLES
      INTEGER :: SA, SMA
      INTEGER :: IX,IY,IZ,P,P1,P2,II,IG,IGP
      INTEGER :: OX,OY,OZ
      INTEGER :: CFG(NX,NZ)
      INTEGER :: XPOS
      INTEGER :: UNITCOUNT
      INTEGER :: POSI(NX*NY*NZ,3)
      INTEGER :: LIST(NX*NY*NZ)
      REAL :: RAN
      INTEGER :: RANI, STORE, HALF, GCOUNT, GPCOUNT
      INTEGER :: POSIG(NXYZ,3), POSIGP(NXYZ,3)
      INTEGER :: LISTG(NXYZ), LISTGP(NXYZ)
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
! CONDITION-2 EDGE CONFIGURATION a/2<110> IN SINGLE PHASE
!-----------------------------------------------------------------------
      ELSE IF(NIN  ==  2) THEN
      PRINT *,'>>NOTE - EDGE a/2[10-1] vs single phase'      
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
          PHI(IX,IY,IZ)=1.0
          CONC(IX,IY,IZ)=CGP
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
! CONDITION-5 READ FIELD VARIABLES FROM FILES init_micro.inp  and also
! iniitalise conc config according to g and gp phases
!-----------------------------------------------------------------------  
      ELSE IF(NIN  ==  5) THEN 
! INITIALISING ETA
        DO IZ=1,NZ
        DO IY=NY/2,NY/2
        DO IX=1,NX
!          IF(IZ >= 15 .AND. IZ < 60) THEN
          IF(IZ >= 5 .AND. IZ < 15) THEN
            ETA(1,3,IX,IY,IZ)=-1.0
            ETA(1,2,IX,IY,IZ)= 1.0
          END IF
        END DO
        END DO
        END DO
! INITIALISING PHI
        DO IY=1,NY
!          OPEN (UNIT = 91, FILE = 'i_micro_238.inp',STATUS = 'UNKNOWN')
          OPEN (UNIT = 91, FILE = '58_120_VF55.inp',STATUS = 'UNKNOWN')
          READ (91,*) ((PHI(IX,IY,IZ),IZ=1,NZ),IX=1,NX)
          CLOSE(91)
        END DO
! INITIALISING CONC 
! COUNTING G AND G' PHASE POINTS
        GCOUNT=0
        DO IZ=1,NZ
        DO IX=1,NX
          IF(PHI(IX,1,IZ) <= 0.5) THEN
            GCOUNT=GCOUNT+1
          END IF
        END DO
        END DO
        PRINT *, 'the 2D GCOUNT ', GCOUNT
        GPCOUNT=(NX*NZ)-GCOUNT
        PRINT *, 'the 2D GPCOUNT ', GPCOUNT
        GCOUNT=NY*GCOUNT
        GPCOUNT=NY*GPCOUNT
! DETERMINE BASIC INITIAL POINTS OF G AND GP
! STORE THEM IN POSIG AND POSIGP
        CALL RANDOM_SEED()
        IG=1
        IGP=1
        DO IZ=1,NZ
        DO IY=1,NY
        DO IX=1,NX
          CONC(IX,IY,IZ)= ( PHI(IX,IY,IZ)*CGP ) + ((1-PHI(IX,IY,IZ))*CG)
          IF(PHI(IX,IY,IZ)<0.5) THEN
            POSIG(IG,1)=IX
            POSIG(IG,2)=IY
            POSIG(IG,3)=IZ
            LISTG(IG)=IG
            IG=IG+1
          ELSE 
            POSIGP(IGP,1)=IX
            POSIGP(IGP,2)=IY
            POSIGP(IGP,3)=IZ
            LISTGP(IGP)=IGP
            IGP=IGP+1
          END IF
        END DO
        END DO
        END DO
        PRINT *, 'the initial sum conc is ', SUM(CONC)
! FOR LOOP TO RANDOMISE POSITIONS IN POSIG AND POSIGP SEPARATELY
        DO P=IG-1,1,-1
          CALL RANDOM_NUMBER(RAN)
          RAN = (RAN*GCOUNT)+1
          RANI= INT(RAN)
          STORE = LISTG(P)
          LISTG(P) = LISTG(RANI)
          LISTG(RANI)= STORE
        END DO
        DO P=IGP-1,1,-1
          CALL RANDOM_NUMBER(RAN)
          RAN = (RAN*GPCOUNT)+1
          RANI= INT(RAN)
          STORE = LISTGP(P)
          LISTGP(P) = LISTGP(RANI)
          LISTGP(RANI)= STORE
        END DO
! PRINT INITIAL CONC BEFORE AND THEN LATER AFTER RANDOMISATION
        OPEN(UNIT=40,FILE='conc_init.dis',STATUS='UNKNOWN')
        WRITE(40,*) (((CONC(IX,IY,IZ),IZ=1,NZ),IY=1,NY),IX=1,NX)
! USING THE RANDOMISED LISTS WE CAN NOW PLUS/MINUS SMALL CONC TO SYSTEM
! HALF IS 1/2 THE NUMBER OF POSITIONS IN G OR G' - SINCE WE PLUS SMALL
! CONC FOR ONE POSIITION AND THEN MINUS FOR ANOTHER ONE (NEXT ONE) DOWN
! THE LIST
        HALF = (IG-1)/2
        DO P=1,HALF
          CALL RANDOM_NUMBER(RAN)
          P1 = 2*P
          CONC(POSIG(LISTG(P1),1),POSIG(LISTG(P1),2), POSIG(LISTG(P1),3)) = CG + (RAN*IRAND)
          P2 = P1 - 1
          CONC(POSIG(LISTG(P2),1),POSIG(LISTG(P2),2), POSIG(LISTG(P2),3)) = CG - (RAN*IRAND)
        END DO
        HALF = (IGP-1)/2
        PRINT *, 'HALFIGP', HALF
        DO P=1,HALF
          CALL RANDOM_NUMBER(RAN)
          P1 = 2*P
          CONC(POSIGP(LISTGP(P1),1),POSIGP(LISTGP(P1),2), POSIGP(LISTGP(P1),3)) = CGP + (RAN*IRAND)
          P2 = P1 - 1
          CONC(POSIGP(LISTGP(P2),1),POSIGP(LISTGP(P2),2), POSIGP(LISTGP(P2),3)) = CGP - (RAN*IRAND)
        END DO
! CHECK TO MAKE SURE CONC HAS BEEN CONSERVED
        PRINT *, 'the new sum conc is ', SUM(CONC)
        WRITE(40,*) (((CONC(IX,IY,IZ),IZ=1,NZ),IY=1,NY),IX=1,NX)
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
          WRITE (UNITCOUNT,*) ((PHI(SLICE_POS,IY,IZ),IZ=1,NZ),IY=1,NY)
          CLOSE(UNITCOUNT)              
        ELSE IF (SLICE_AXIS == 2) THEN
! SLICE ALONG Y-AXIS AT SLICE_POS        
          UNITCOUNT=91            
          WRITE (FILENAME,'(A7)')'phi.dis'                 
          OPEN (UNIT = UNITCOUNT, FILE=FILENAME, FORM = 'FORMATTED')
          WRITE (UNITCOUNT,*) ((PHI(IX,SLICE_POS,IZ),IZ=1,NZ),IX=1,NX)
          CLOSE(UNITCOUNT)              
        ELSE IF (SLICE_AXIS == 3) THEN
! SLICE ALONG Z-AXIS AT SLICE_POS        
          UNITCOUNT=91            
          WRITE (FILENAME,'(A7)')'phi.dis'                 
          OPEN (UNIT = UNITCOUNT, FILE=FILENAME, FORM = 'FORMATTED')
          WRITE (UNITCOUNT,*) ((PHI(IX,IY,SLICE_POS),IY=1,NY),IX=1,NX)
          CLOSE(UNITCOUNT)                       
        END IF
! WRITE 3D MICROSTRUCTURE TO FILE        
      ELSE IF (OUTDIM == 3) THEN
!        ELSE IF (OUTDIM == 3) THEN 
!
          UNITCOUNT=91
!          
          WRITE (FILENAME,'(A7)')'phi.dis'                 
          OPEN (UNIT = UNITCOUNT, FILE=FILENAME, FORM = 'FORMATTED')
          WRITE(UNITCOUNT,*) (((PHI(IX,IY,IZ),IZ=1,NZ),IY=1,NY),IX=1,NX)
          CLOSE(UNITCOUNT)    
      END IF
      END SUBROUTINE MIC
