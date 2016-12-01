!-----------------------------------------------------------------------
! SUBROUTINE DEFINING ACTIVE SLIP SYSTEMS
!-----------------------------------------------------------------------
      SUBROUTINE SLP
      USE VAR 
      IMPLICIT NONE
      INTEGER :: I,SA,SMA        !! INDICES AND IDENTIFIERS
      REAL :: TV1(3),TV2(3)      !! TEMPORARY VECTORS
! DEFINE REDUCED LATTICE PARAMETERS
      A0=AA/AA
! DEFINE BURGERS VECTORS AND SLIP PLANE NORMALS
!-----------------------------------------------------------------------
! SLIP SYSTEM #1 - (111) plane in FCC, 
!-----------------------------------------------------------------------
! (111) plane in FCC standard crystal axes x,y,z=[100],[010],[001]
! N.B. Useful to rotate the crystal axes to x,y,z=[10-1],[111],[1-21]
      IF(NSS  ==  1)  THEN
!
      NP=1
      NQ=3
!
      PRINT *,'>>NOTE - (111) PLANE IN FCC CRYSTAL'
      PRINT *,'>>NOTE - CRYSTAL AXES X,Y,Z=[100],[010],[001]'   
!        
          NV(1,1)=1.0/SQRT(3.0)      
          NV(1,2)=1.0/SQRT(3.0)
          NV(1,3)=1.0/SQRT(3.0)
!        
          BV(1,1,1)= 0.0     
          BV(1,1,2)=-A0
          BV(1,1,3)= A0
!        
          BV(1,2,1)= A0
          BV(1,2,2)= 0.0
          BV(1,2,3)=-A0
!        
          BV(1,3,1)=-A0
          BV(1,3,2)= A0
          BV(1,3,3)= 0.0
!        
!-----------------------------------------------------------------------
! SLIP SYSTEM #2 - (111) FULL OCTAHEDRAL SLIP IN FCC 
!-----------------------------------------------------------------------
! Full octahedral FCC with standard crystal axes x,y,z=[100],[010],[001]
      ELSE IF (NSS == 2) THEN
!
      NP=4
      NQ=3
!
      PRINT *,'>>NOTE - FULL OCTAHEDRAL SLIP IN FCC CRYSTAL'
      PRINT *,'>>NOTE - CRYSTAL AXES X,Y,Z=[100],[010],[001]'
!        
          NV(1,1)= 1.0/SQRT(3.0)      
          NV(1,2)= 1.0/SQRT(3.0)
          NV(1,3)= 1.0/SQRT(3.0)
!
          BV(1,1,1)= 0.0     
          BV(1,1,2)=-A0
          BV(1,1,3)= A0
!        
          BV(1,2,1)= A0
          BV(1,2,2)= 0.0
          BV(1,2,3)=-A0
!        
          BV(1,3,1)=-A0
          BV(1,3,2)= A0
          BV(1,3,3)= 0.0
!        
          NV(2,1)=-1.0/SQRT(3.0)      
          NV(2,2)= 1.0/SQRT(3.0)
          NV(2,3)= 1.0/SQRT(3.0)
!
          BV(2,1,1)= 0.0      
          BV(2,1,2)= A0
          BV(2,1,3)=-A0
!        
          BV(2,2,1)= A0
          BV(2,2,2)= 0.0
          BV(2,2,3)= A0
!        
          BV(2,3,1)=-A0
          BV(2,3,2)=-A0
          BV(2,3,3)= 0.0                  
!        
          NV(3,1)= 1.0/SQRT(3.0)      
          NV(3,2)=-1.0/SQRT(3.0)
          NV(3,3)= 1.0/SQRT(3.0)
!
          BV(3,1,1)= 0.0     
          BV(3,1,2)= A0
          BV(3,1,3)= A0
!        
          BV(3,2,1)= A0
          BV(3,2,2)= 0.0
          BV(3,2,3)=-A0
!        
          BV(3,3,1)=-A0
          BV(3,3,2)=-A0
          BV(3,3,3)= 0.0
!        
          NV(4,1)= 1.0/SQRT(3.0)      
          NV(4,2)= 1.0/SQRT(3.0)
          NV(4,3)=-1.0/SQRT(3.0)
!
          BV(4,1,1)= 0.0     
          BV(4,1,2)=-A0
          BV(4,1,3)=-A0
!        
          BV(4,2,1)= A0
          BV(4,2,2)= 0.0
          BV(4,2,3)= A0
!        
          BV(4,3,1)=-A0 
          BV(4,3,2)= A0
          BV(4,3,3)= 0.0
!              
      ELSE
        STOP '>>ERROR - INVALID SLIP SYSTEM!'
      END IF
!-----------------------------------------------------------------------
! ROTATE SLIP SYSTEM TO NEW AXIS SYSTEM IF ROTFLAG=1 
!-----------------------------------------------------------------------      
      IF (ROTFLAG == 0) THEN          
      PRINT *,'>>NOTE - NO ROTATION APPLIED TO CRYSTAL'            
      ELSE IF (ROTFLAG == 1) THEN
      PRINT *,'>>NOTE - APPLYING ROTATION MATRIX TO CRYSTAL AXES'           
!        
      CALL ORT
!
      CALL ROM
!        
        DO SA=1,NP
          DO I=1,3
            TV1(I)=NV(SA,I)
          END DO
          CALL VMM(TV2,TV1,ROT)
          DO I=1,3
            NV(SA,I)=TV2(I)
          END DO
        END DO  
!         
        DO SA=1,NP
        DO SMA=1,NQ
          DO I=1,3
            TV1(I)=BV(SA,SMA,I)
          END DO
          CALL VMM(TV2,TV1,ROT)
          DO I=1,3
            BV(SA,SMA,I)=TV2(I)
          END DO
        END DO
        END DO
!          
      END IF      
      END SUBROUTINE SLP
!-----------------------------------------------------------------------
! SUBROUTINE TO CALCULATE CRYSTAL ROTATION MATRIX
!-----------------------------------------------------------------------
      SUBROUTINE ROM
      USE VAR 
      IMPLICIT NONE
!      
        ROT(1,1)=AXX(1)/(AXX(1)**2.0+AXX(2)**2.0+AXX(3)**2.0)**0.5
        ROT(1,2)=AXX(2)/(AXX(1)**2.0+AXX(2)**2.0+AXX(3)**2.0)**0.5
        ROT(1,3)=AXX(3)/(AXX(1)**2.0+AXX(2)**2.0+AXX(3)**2.0)**0.5
        ROT(2,1)=AXY(1)/(AXY(1)**2.0+AXY(2)**2.0+AXY(3)**2.0)**0.5
        ROT(2,2)=AXY(2)/(AXY(1)**2.0+AXY(2)**2.0+AXY(3)**2.0)**0.5
        ROT(2,3)=AXY(3)/(AXY(1)**2.0+AXY(2)**2.0+AXY(3)**2.0)**0.5
        ROT(3,1)=AXZ(1)/(AXZ(1)**2.0+AXZ(2)**2.0+AXZ(3)**2.0)**0.5
        ROT(3,2)=AXZ(2)/(AXZ(1)**2.0+AXZ(2)**2.0+AXZ(3)**2.0)**0.5
        ROT(3,3)=AXZ(3)/(AXZ(1)**2.0+AXZ(2)**2.0+AXZ(3)**2.0)**0.5
!
      PRINT *,'>>NOTE - ORIGINAL AXES ROTATED BY THE FOLLOWING MATRIX'
      PRINT '(F9.5,F9.5,F9.5)',ROT(1,1),ROT(1,2),ROT(1,3)
      PRINT '(F9.5,F9.5,F9.5)',ROT(2,1),ROT(2,2),ROT(2,3)
      PRINT '(F9.5,F9.5,F9.5)',ROT(3,1),ROT(3,2),ROT(3,3)
      END SUBROUTINE ROM
!-----------------------------------------------------------------------
! SUBROUTINE TO CHECK ORTHOGONALITY OF AXES
!-----------------------------------------------------------------------
      SUBROUTINE ORT
      USE VAR 
      IMPLICIT NONE
      REAL :: CHECK
!      
      CALL DOT(CHECK,AXX,AXY)
      IF (CHECK /= 0) THEN
        STOP '>>ERROR - CRYSTAL AXES NOT ORTHOGONAL!'
      END IF
!      
      CALL DOT(CHECK,AXX,AXZ)
      IF (CHECK /= 0) THEN
        STOP '>>ERROR - CRYSTAL AXES NOT ORTHOGONAL!'
      END IF
!      
      CALL DOT(CHECK,AXZ,AXY)
      IF (CHECK /= 0) THEN
        STOP '>>ERROR - CRYSTAL AXES NOT ORTHOGONAL!'
      END IF
!      
      END SUBROUTINE ORT   