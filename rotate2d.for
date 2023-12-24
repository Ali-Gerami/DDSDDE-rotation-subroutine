c ===========================================================================
c
c       **   DDSDDE or Stiffness Matrix Rotation Subroutine    **
c      ***  for 2d Simulations (Plane Stress or Plane Strain)  ***
c      
c      -Inputs
c        ALPHA : Rotation angle (in rad, CCW is positive) 
c           D0 : Initial DDSDDE or stiffness matrix
c        NTENS : Size of the stress or strain component array
c      
c      -Output
c           D1 : New DDSDDE or stiffness matrix
c          
c ===========================================================================
       
      SUBROUTINE ROTATE(ALPHA,D0,NTENS,D1)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION ROT(3,3),D0(NTENS,NTENS),D1(NTENS,NTENS)
      DIMENSION D(6,6),RD(6,6),T(6,6),TT(6,6),M(6)
      
      
      IF ((NTENS==3).OR.(NTENS==4)) THEN       ! *Acceptable values of NTENS*
      
      ROT = 0.0D0
      ROT(1,1) = COS(ALPHA)
      ROT(1,2) = SIN(ALPHA)
      ROT(2,1) = -SIN(ALPHA)
      ROT(2,2) = COS(ALPHA)
      ROT(3,3) = 1.0D0
      
      D=0.0D0
      DO i=1,NTENS-1
      DO j=1,NTENS-1
        D(i,j) = D0(i,j)
      END DO
      END DO
      D(1,6) = D0(1,NTENS)
      D(2,6) = D0(2,NTENS)
      D(6,6) = D0(NTENS,NTENS)
      D(6,1) = D(1,6)
      D(6,2) = D(2,6)
      
      M(1) = 1
      M(2) = 2
      M(3) = 3
      M(4) = 1
      M(5) = 2
      
      
C     CALCULATE T-MATRIX:
      
      DO i=1,3
       DO j=1,3
          T(i,j) = ROT(i,j)*ROT(i,j)
       END DO
       DO j=4,6
          T(i,j) = 2.0D0*(ROT(i,M(j-2))*ROT(i,M(j-1)))
       END DO
      END DO
      
      DO i=4,6
       DO j=1,3
          T(i,j) = ROT(M(i-2),j)*ROT(M(i-1),j)
       END DO
       DO j=4,6
          T(i,j) = (ROT(M(i-2),M(j-2))*ROT(M(i-1),M(j-1))) +
     +            (ROT(M(i-2),M(j-1))*ROT(M(i-1),M(j-2)))
       END DO
      END DO
      
      
C     CALCULATE TRANSPOSE OF T-MATRIX:
      
      DO i=1,6
      DO j=1,6
        TT(i,j)=T(j,i)
      END DO
      END DO
      
      
C     ROTATE D-MATRIX:
      
      RD = MATMUL(MATMUL(T,D),TT)
      
      
C     OUTPUT PREPARATION:
      
      D1 = 0.0D0
      DO i=1,NTENS-1
      DO j=1,NTENS-1
        D1(i,j) = RD(i,j)
      END DO
      END DO
      D1(1,NTENS) = RD(1,6)
      D1(2,NTENS) = RD(2,6)
      D1(NTENS,NTENS) = RD(6,6)
      D1(NTENS,1) = D1(1,NTENS)
      D1(NTENS,2) = D1(2,NTENS)
      
      
      ELSE            ! *Unacceptable values of NTENS*
          
                      !This massage will print in .log file:
       print *,"* ERROR! NTENS IS NOT ACCEPTABLE IN ROTATE SUBROUTINE *"
       
      END IF
      
      RETURN
      END
          

