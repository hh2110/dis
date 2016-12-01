!-----------------------------------------------------------------------
! SUBROUTINE FOR 3X3 MATRIX CREATOR
!-----------------------------------------------------------------------
      SUBROUTINE MAT(M,M11,M12,M13,M21,M22,M23,M31,M32,M33)
      IMPLICIT NONE
! DEFINE VARIABLES
      REAL :: M(3,3) !! 3X3 MATRIX
      REAL :: M11,M12,M13,M21,M22,M23,M31,M32,M33 !! MATRIX COMPONENTS
! ASSIGN INDECES TO MATRIX COMPONENTS
      M(1,1)=M11
      M(1,2)=M12
      M(1,3)=M13
      M(2,1)=M21
      M(2,2)=M22
      M(2,3)=M23
      M(3,1)=M31
      M(3,2)=M32
      M(3,3)=M33
      END SUBROUTINE MAT
!-----------------------------------------------------------------------
! SUBROUTINE FOR CALCULATING INVERSE OF A 3X3 MATRIX
!-----------------------------------------------------------------------
      SUBROUTINE INV(A,RA) 
! DEFINE VARIABLES (MATRICES)
      REAL, INTENT(IN)  :: A(3,3)   !! INPUT MATRIX FOR COMPUTING INVERSE OF
      REAL, INTENT(OUT) :: RA(3,3)  !! INVERSE OF INPUT MATRIX
      REAL :: B(3,3)                !! DUMMY VARIABLE MATRIX
! DEFINE DETERMINANT OF A
      DET_A = A(1,1)*A(2,2)*A(3,3)
     &  + A(2,1)*A(3,2)*A(1,3)+A(1,2)*A(2,3)*A(3,1)
     &  - A(3,1)*A(2,2)*A(1,3) 
     &  - A(2,1)*A(1,2)*A(3,3)-A(3,2)*A(2,3)*A(1,1)
! A^-1=B(I,J)/DET_A 
      B(1,1)=A(2,2)*A(3,3)-A(3,2)*A(2,3)
      B(1,2)=A(1,2)*A(3,3)-A(3,2)*A(1,3)
      B(1,3)=A(1,2)*A(2,3)-A(2,2)*A(1,3)
      B(2,1)=A(2,1)*A(3,3)-A(3,1)*A(2,3)
      B(2,2)=A(1,1)*A(3,3)-A(3,1)*A(1,3)
      B(2,3)=A(1,1)*A(2,3)-A(2,1)*A(1,3)
      B(3,1)=A(2,1)*A(3,2)-A(3,1)*A(2,2)
      B(3,2)=A(1,1)*A(3,2)-A(3,1)*A(1,2)
      B(3,3)=A(1,1)*A(2,2)-A(2,1)*A(1,2)
! CALCULATE INVERSE
      DO I=1,3
      DO J=1,3
        IF(MOD(I+J,2) == 0)THEN
          RA(I,J)=B(I,J)/DET_A 
        ELSE
          RA(I,J)=-B(I,J)/DET_A
        END IF
      END DO
      END DO
      END SUBROUTINE INV
!-----------------------------------------------------------------------
! SUBROUTINE TO FIND DOT PRODUCT OF TWO VECTORS
!-----------------------------------------------------------------------
      SUBROUTINE DOT(PROD,V1,V2)
      IMPLICIT NONE
! DEFINE VARIABLES (MATRICES)
      REAL, INTENT(IN) :: V1(3)   !! VECTOR 1 
      REAL, INTENT(IN) :: V2(3)   !! VECTOR 2
      REAL, INTENT(OUT):: PROD    !! DOT PRODUCT OF V1 AND V2 
! CALCULATE THE TERMS
      PROD=V1(1)*V2(1)+V1(2)*V2(2)+V1(3)*V2(3)
      END SUBROUTINE DOT
!-----------------------------------------------------------------------
! SUBROUTINE FOR MULTIPLYING 3X3 MATRICES
!-----------------------------------------------------------------------
      SUBROUTINE MUL (PROD,M1,M2)
      IMPLICIT NONE
! DEFINE VARIABLES (MATRICES)
      REAL, INTENT(IN) :: M1(3,3)     !! FIRST MATRIX 
      REAL, INTENT(IN) :: M2(3,3)     !! SECOND MATRIX
      REAL, INTENT(OUT):: PROD(3,3)   !! PRODUCT MATRIX
! CALCULATE THE TERMS
      PROD(1,1)=M1(1,1)*M2(1,1)+M1(1,2)*M2(2,1)+M1(1,3)*M2(3,1)
      PROD(1,2)=M1(1,1)*M2(1,2)+M1(1,2)*M2(2,2)+M1(1,3)*M2(3,2)
      PROD(1,3)=M1(1,1)*M2(1,3)+M1(1,2)*M2(2,3)+M1(1,3)*M2(3,3)
      PROD(2,1)=M1(2,1)*M2(1,1)+M1(2,2)*M2(2,1)+M1(2,3)*M2(3,1)
      PROD(2,2)=M1(2,1)*M2(1,2)+M1(2,2)*M2(2,2)+M1(2,3)*M2(3,2)
      PROD(2,3)=M1(2,1)*M2(1,3)+M1(2,2)*M2(2,3)+M1(2,3)*M2(3,3)
      PROD(3,1)=M1(3,1)*M2(1,1)+M1(3,2)*M2(2,1)+M1(3,3)*M2(3,1)
      PROD(3,2)=M1(3,1)*M2(1,2)+M1(3,2)*M2(2,2)+M1(3,3)*M2(3,2)
      PROD(3,3)=M1(3,1)*M2(1,3)+M1(3,2)*M2(2,3)+M1(3,3)*M2(3,3)
      END SUBROUTINE MUL
!-----------------------------------------------------------------------
! SUBROUTINE FOR MULTIPLYING A VECTOR BY A MATRIX
!-----------------------------------------------------------------------
      SUBROUTINE VMM (PROD,V1,M1)
      IMPLICIT NONE
! DEFINE VARIABLES (MATRICES)
      REAL, INTENT(IN) :: V1(3)
      REAL, INTENT(IN) :: M1(3,3)
      REAL, INTENT(OUT):: PROD(3)
! CALCULATE THE TERMS
      PROD(1)=M1(1,1)*V1(1)+M1(1,2)*V1(2)+M1(1,3)*V1(3)
      PROD(2)=M1(2,1)*V1(1)+M1(2,2)*V1(2)+M1(2,3)*V1(3)
      PROD(3)=M1(3,1)*V1(1)+M1(3,2)*V1(2)+M1(3,3)*V1(3)
      END SUBROUTINE VMM