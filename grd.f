!-----------------------------------------------------------------------
! SUBROUTINE TO CALCULATE K^2 COMPONENTS IN K SPACE
!-----------------------------------------------------------------------
      SUBROUTINE GRD
      USE VAR 
      IMPLICIT NONE
! DEFINE VARIABLES
      INTEGER :: IX,IY,IZ
! BEGIN LOOP
!$OMP PARALLEL DO PRIVATE(IY,IX)
      DO IZ =1, NZ
      DO IY =1, NY
      DO IX =1, NX
        IF(IX <= NX/2+1) THEN
          G(1,IX,IY,IZ) =TWOPI/NX *(IX-1)
        ELSE
          G(1,IX,IY,IZ) =TWOPI/NX *(IX-1-NX)
        END IF
        IF(IY <= NY/2+1) THEN
          G(2,IX,IY,IZ) =TWOPI/NY *(IY-1)
        ELSE
          G(2,IX,IY,IZ) =TWOPI/NY *(IY-1-NY)
        END IF
        IF(IZ <= NZ/2+1) THEN
          G(3,IX,IY,IZ) =TWOPI/NZ *(IZ-1)
        ELSE
          G(3,IX,IY,IZ) =TWOPI/NZ *(IZ-1-NZ)
        END IF
        G2(IX,IY,IZ) =G(1,IX,IY,IZ)**2 + G(2,IX,IY,IZ)**2 +
     &    G(3,IX,IY,IZ)**2
      END DO
      END DO
      END DO
!$OMP END PARALLEL DO      
      END SUBROUTINE GRD
!-----------------------------------------------------------------------
! SUBROUTINE TO CALCULATE DETA(SA,SMA)/DXK
!-----------------------------------------------------------------------
      SUBROUTINE GRE
      USE VAR 
      IMPLICIT NONE
! DEFINE VARIABLES
      INTEGER :: IX,IY,IZ
      INTEGER :: SA,SMA
      INTEGER :: K
      COMPLEX :: COMPLEX_I=CMPLX(0,1)
!     CALCULATE DETA(SA,SMA)/DXK
      DO SA =1,NP
      DO SMA=1,NQ
        DO K=1,3
!$OMP PARALLEL DO PRIVATE(IY,IX)        
          DO IZ=1,NZ
          DO IY=1,NY
          DO IX=1,NX
            TPK(IX,IY,IZ) =COMPLEX_I*G(K,IX,IY,IZ)*
     &        ETAK(SA,SMA,IX,IY,IZ)
          END DO
          END DO
          END DO
!$OMP END PARALLEL DO           
            CALL INVFFT_C2R(TPK,TPR)
!$OMP PARALLEL DO PRIVATE(IY,IX)             
          DO IZ=1,NZ
          DO IY=1,NY
          DO IX=1,NX
            DER(SA,SMA,K,IX,IY,IZ) =TPR(IX,IY,IZ)
          END DO
          END DO
          END DO
!$OMP END PARALLEL DO           
        END DO
      END DO
      END DO
      END SUBROUTINE GRE