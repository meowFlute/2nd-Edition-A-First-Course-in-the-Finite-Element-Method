! TRUSS Finite Element Analysis Program
! Solves 3D truss structures using the direct stiffness method
! Fortran 90/95 implementation

PROGRAM TRUSS
    IMPLICIT NONE
    
    ! Maximum dimensions
    INTEGER, PARAMETER :: MAXNODE = 100
    INTEGER, PARAMETER :: MAXELE = 100
    INTEGER, PARAMETER :: MAXDOF = 300
    
    ! Problem variables
    INTEGER :: NELE, NNODE, NDOF, MUD
    INTEGER :: IFIX(3, MAXNODE)
    INTEGER :: NODE(2, MAXELE)
    REAL(8) :: XC(MAXNODE), YC(MAXNODE), ZC(MAXNODE)
    REAL(8) :: FORCE(3, MAXNODE)
    REAL(8) :: E(MAXELE), A(MAXELE)
    REAL(8) :: K(MAXDOF, MAXDOF), F(MAXDOF), D(MAXDOF)
    REAL(8) :: STRESS(MAXELE)
    CHARACTER(LEN=80) :: TITLE
    CHARACTER(LEN=256) :: FILENAME
    
    INTEGER :: I, J, KE, IERR
    
    ! Get input filename from command line
    IF (COMMAND_ARGUMENT_COUNT() .NE. 1) THEN
        WRITE(*,*) 'Usage: truss <input_file>'
        STOP
    END IF
    
    CALL GET_COMMAND_ARGUMENT(1, FILENAME)
    
    ! Read input data
    CALL READ_DATA(FILENAME, TITLE, NELE, NNODE, IFIX, XC, YC, ZC, &
                   FORCE, NODE, E, A, MAXNODE, MAXELE)
    
    NDOF = 3 * NNODE
    
    ! Initialize stiffness matrix and force vector
    K = 0.0D0
    F = 0.0D0
    
    ! Assemble global stiffness matrix and force vector
    CALL ASSEMBLE(NELE, NNODE, NODE, XC, YC, ZC, E, A, K, F, FORCE, &
                  MAXNODE, MAXELE, MAXDOF, MUD)
    
    ! Apply boundary conditions and solve
    CALL SOLVE_SYSTEM(K, F, D, IFIX, NNODE, NDOF, MAXNODE, MAXDOF)
    
    ! Compute element stresses
    CALL COMPUTE_STRESS(NELE, NODE, XC, YC, ZC, E, D, STRESS, &
                        MAXNODE, MAXELE, MAXDOF)
    
    ! Output results
    CALL WRITE_OUTPUT(TITLE, NELE, NNODE, IFIX, XC, YC, ZC, FORCE, &
                      NODE, E, A, D, STRESS, MUD, MAXNODE, MAXELE)
    
END PROGRAM TRUSS


!***********************************************************************
SUBROUTINE READ_DATA(FILENAME, TITLE, NELE, NNODE, IFIX, XC, YC, ZC, &
                     FORCE, NODE, E, A, MAXNODE, MAXELE)
!***********************************************************************
    IMPLICIT NONE
    
    CHARACTER(LEN=*), INTENT(IN) :: FILENAME
    CHARACTER(LEN=80), INTENT(OUT) :: TITLE
    INTEGER, INTENT(IN) :: MAXNODE, MAXELE
    INTEGER, INTENT(OUT) :: NELE, NNODE
    INTEGER, INTENT(OUT) :: IFIX(3, MAXNODE)
    INTEGER, INTENT(OUT) :: NODE(2, MAXELE)
    REAL(8), INTENT(OUT) :: XC(MAXNODE), YC(MAXNODE), ZC(MAXNODE)
    REAL(8), INTENT(OUT) :: FORCE(3, MAXNODE)
    REAL(8), INTENT(OUT) :: E(MAXELE), A(MAXELE)
    
    INTEGER :: I, J, K, IUNIT
    
    IUNIT = 10
    OPEN(UNIT=IUNIT, FILE=TRIM(FILENAME), STATUS='OLD', ACTION='READ')
    
    ! Read title
    READ(IUNIT, '(A)') TITLE
    
    ! Read number of elements and nodes
    READ(IUNIT, *) NELE, NNODE
    
    ! Initialize arrays
    IFIX = 0
    XC = 0.0D0
    YC = 0.0D0
    ZC = 0.0D0
    FORCE = 0.0D0
    
    ! Read node data
    DO I = 1, NNODE
        READ(IUNIT, *) J, IFIX(1,J), IFIX(2,J), IFIX(3,J), &
                       XC(J), YC(J), ZC(J), &
                       FORCE(1,J), FORCE(2,J), FORCE(3,J)
    END DO
    
    ! Read element data
    DO I = 1, NELE
        READ(IUNIT, *) K, NODE(1,K), NODE(2,K), E(K), A(K)
    END DO
    
    CLOSE(IUNIT)
    
END SUBROUTINE READ_DATA


!***********************************************************************
SUBROUTINE ASSEMBLE(NELE, NNODE, NODE, XC, YC, ZC, E, A, K, F, FORCE, &
                    MAXNODE, MAXELE, MAXDOF, MUD)
!***********************************************************************
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: NELE, NNODE, MAXNODE, MAXELE, MAXDOF
    INTEGER, INTENT(IN) :: NODE(2, MAXELE)
    REAL(8), INTENT(IN) :: XC(MAXNODE), YC(MAXNODE), ZC(MAXNODE)
    REAL(8), INTENT(IN) :: E(MAXELE), A(MAXELE)
    REAL(8), INTENT(IN) :: FORCE(3, MAXNODE)
    REAL(8), INTENT(INOUT) :: K(MAXDOF, MAXDOF), F(MAXDOF)
    INTEGER, INTENT(OUT) :: MUD
    
    INTEGER :: KE, NODE1, NODE2, I, J, II, JJ
    INTEGER :: DOF(6)
    REAL(8) :: X1, Y1, Z1, X2, Y2, Z2
    REAL(8) :: DX, DY, DZ, EL, CX, CY, CZ
    REAL(8) :: KELEM(6,6), FACTOR
    
    ! Loop over all elements
    DO KE = 1, NELE
        NODE1 = NODE(1, KE)
        NODE2 = NODE(2, KE)
        
        X1 = XC(NODE1)
        Y1 = YC(NODE1)
        Z1 = ZC(NODE1)
        X2 = XC(NODE2)
        Y2 = YC(NODE2)
        Z2 = ZC(NODE2)
        
        ! Calculate element length
        DX = X2 - X1
        DY = Y2 - Y1
        DZ = Z2 - Z1
        EL = DSQRT(DX*DX + DY*DY + DZ*DZ)
        
        ! Direction cosines
        CX = DX / EL
        CY = DY / EL
        CZ = DZ / EL
        
        ! Element stiffness matrix
        FACTOR = E(KE) * A(KE) / EL
        KELEM = 0.0D0
        
        ! Build 6x6 element stiffness matrix
        DO I = 1, 3
            DO J = 1, 3
                IF (I == 1) THEN
                    IF (J == 1) THEN
                        KELEM(I,J) = FACTOR * CX * CX
                        KELEM(I,J+3) = -FACTOR * CX * CX
                        KELEM(I+3,J) = -FACTOR * CX * CX
                        KELEM(I+3,J+3) = FACTOR * CX * CX
                    ELSE IF (J == 2) THEN
                        KELEM(I,J) = FACTOR * CX * CY
                        KELEM(I,J+3) = -FACTOR * CX * CY
                        KELEM(I+3,J) = -FACTOR * CX * CY
                        KELEM(I+3,J+3) = FACTOR * CX * CY
                    ELSE IF (J == 3) THEN
                        KELEM(I,J) = FACTOR * CX * CZ
                        KELEM(I,J+3) = -FACTOR * CX * CZ
                        KELEM(I+3,J) = -FACTOR * CX * CZ
                        KELEM(I+3,J+3) = FACTOR * CX * CZ
                    END IF
                ELSE IF (I == 2) THEN
                    IF (J == 1) THEN
                        KELEM(I,J) = FACTOR * CY * CX
                        KELEM(I,J+3) = -FACTOR * CY * CX
                        KELEM(I+3,J) = -FACTOR * CY * CX
                        KELEM(I+3,J+3) = FACTOR * CY * CX
                    ELSE IF (J == 2) THEN
                        KELEM(I,J) = FACTOR * CY * CY
                        KELEM(I,J+3) = -FACTOR * CY * CY
                        KELEM(I+3,J) = -FACTOR * CY * CY
                        KELEM(I+3,J+3) = FACTOR * CY * CY
                    ELSE IF (J == 3) THEN
                        KELEM(I,J) = FACTOR * CY * CZ
                        KELEM(I,J+3) = -FACTOR * CY * CZ
                        KELEM(I+3,J) = -FACTOR * CY * CZ
                        KELEM(I+3,J+3) = FACTOR * CY * CZ
                    END IF
                ELSE IF (I == 3) THEN
                    IF (J == 1) THEN
                        KELEM(I,J) = FACTOR * CZ * CX
                        KELEM(I,J+3) = -FACTOR * CZ * CX
                        KELEM(I+3,J) = -FACTOR * CZ * CX
                        KELEM(I+3,J+3) = FACTOR * CZ * CX
                    ELSE IF (J == 2) THEN
                        KELEM(I,J) = FACTOR * CZ * CY
                        KELEM(I,J+3) = -FACTOR * CZ * CY
                        KELEM(I+3,J) = -FACTOR * CZ * CY
                        KELEM(I+3,J+3) = FACTOR * CZ * CY
                    ELSE IF (J == 3) THEN
                        KELEM(I,J) = FACTOR * CZ * CZ
                        KELEM(I,J+3) = -FACTOR * CZ * CZ
                        KELEM(I+3,J) = -FACTOR * CZ * CZ
                        KELEM(I+3,J+3) = FACTOR * CZ * CZ
                    END IF
                END IF
            END DO
        END DO
        
        ! Global DOF indices for this element
        DOF(1) = (NODE1 - 1) * 3 + 1
        DOF(2) = (NODE1 - 1) * 3 + 2
        DOF(3) = (NODE1 - 1) * 3 + 3
        DOF(4) = (NODE2 - 1) * 3 + 1
        DOF(5) = (NODE2 - 1) * 3 + 2
        DOF(6) = (NODE2 - 1) * 3 + 3
        
        ! Add element stiffness to global stiffness
        DO I = 1, 6
            DO J = 1, 6
                K(DOF(I), DOF(J)) = K(DOF(I), DOF(J)) + KELEM(I, J)
            END DO
        END DO
    END DO
    
    ! Assemble force vector
    DO J = 1, NNODE
        F((J-1)*3 + 1) = FORCE(1, J)
        F((J-1)*3 + 2) = FORCE(2, J)
        F((J-1)*3 + 3) = FORCE(3, J)
    END DO
    
    ! Calculate MUD (maximum upper codiagonal)
    MUD = 0
    DO I = 1, 3 * NNODE
        DO J = I + 1, 3 * NNODE
            IF (DABS(K(I,J)) .GT. 1.0D-10) THEN
                MUD = MAX(MUD, J - I)
            END IF
        END DO
    END DO
    
END SUBROUTINE ASSEMBLE


!***********************************************************************
SUBROUTINE SOLVE_SYSTEM(K, F, D, IFIX, NNODE, NDOF, MAXNODE, MAXDOF)
!***********************************************************************
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: NNODE, NDOF, MAXNODE, MAXDOF
    INTEGER, INTENT(IN) :: IFIX(3, MAXNODE)
    REAL(8), INTENT(INOUT) :: K(MAXDOF, MAXDOF), F(MAXDOF)
    REAL(8), INTENT(OUT) :: D(MAXDOF)
    
    INTEGER :: I, J, IDOF
    REAL(8) :: BIG
    
    ! Apply boundary conditions using penalty method
    BIG = 1.0D20
    
    DO J = 1, NNODE
        DO I = 1, 3
            IF (IFIX(I, J) .EQ. 1) THEN
                IDOF = (J - 1) * 3 + I
                K(IDOF, IDOF) = K(IDOF, IDOF) + BIG
                F(IDOF) = 0.0D0
            END IF
        END DO
    END DO
    
    ! Solve using Gaussian elimination
    CALL GAUSS(K, F, D, NDOF, MAXDOF)
    
    ! Explicitly zero out fixed displacements to avoid numerical noise
    DO J = 1, NNODE
        DO I = 1, 3
            IF (IFIX(I, J) .EQ. 1) THEN
                IDOF = (J - 1) * 3 + I
                D(IDOF) = 0.0D0
            END IF
        END DO
    END DO
    
END SUBROUTINE SOLVE_SYSTEM


!***********************************************************************
SUBROUTINE GAUSS(A, B, X, N, MAXDOF)
!***********************************************************************
! Gaussian elimination with partial pivoting
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: N, MAXDOF
    REAL(8), INTENT(INOUT) :: A(MAXDOF, MAXDOF), B(MAXDOF)
    REAL(8), INTENT(OUT) :: X(MAXDOF)
    
    INTEGER :: I, J, K, IMAX
    REAL(8) :: AMAX, TEMP, FACTOR
    
    ! Forward elimination
    DO K = 1, N - 1
        ! Partial pivoting
        AMAX = DABS(A(K, K))
        IMAX = K
        DO I = K + 1, N
            IF (DABS(A(I, K)) .GT. AMAX) THEN
                AMAX = DABS(A(I, K))
                IMAX = I
            END IF
        END DO
        
        ! Swap rows if necessary
        IF (IMAX .NE. K) THEN
            DO J = K, N
                TEMP = A(K, J)
                A(K, J) = A(IMAX, J)
                A(IMAX, J) = TEMP
            END DO
            TEMP = B(K)
            B(K) = B(IMAX)
            B(IMAX) = TEMP
        END IF
        
        ! Eliminate column
        DO I = K + 1, N
            IF (DABS(A(K, K)) .GT. 1.0D-10) THEN
                FACTOR = A(I, K) / A(K, K)
                DO J = K + 1, N
                    A(I, J) = A(I, J) - FACTOR * A(K, J)
                END DO
                B(I) = B(I) - FACTOR * B(K)
            END IF
        END DO
    END DO
    
    ! Back substitution
    X(N) = B(N) / A(N, N)
    DO I = N - 1, 1, -1
        X(I) = B(I)
        DO J = I + 1, N
            X(I) = X(I) - A(I, J) * X(J)
        END DO
        X(I) = X(I) / A(I, I)
    END DO
    
END SUBROUTINE GAUSS


!***********************************************************************
SUBROUTINE COMPUTE_STRESS(NELE, NODE, XC, YC, ZC, E, D, STRESS, &
                          MAXNODE, MAXELE, MAXDOF)
!***********************************************************************
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: NELE, MAXNODE, MAXELE, MAXDOF
    INTEGER, INTENT(IN) :: NODE(2, MAXELE)
    REAL(8), INTENT(IN) :: XC(MAXNODE), YC(MAXNODE), ZC(MAXNODE)
    REAL(8), INTENT(IN) :: E(MAXELE), D(MAXDOF)
    REAL(8), INTENT(OUT) :: STRESS(MAXELE)
    
    INTEGER :: K, NODE1, NODE2
    REAL(8) :: X1, Y1, Z1, X2, Y2, Z2
    REAL(8) :: DX, DY, DZ, EL, CX, CY, CZ
    REAL(8) :: D1(3), D2(3), DELTA_L, STRAIN
    
    DO K = 1, NELE
        NODE1 = NODE(1, K)
        NODE2 = NODE(2, K)
        
        X1 = XC(NODE1)
        Y1 = YC(NODE1)
        Z1 = ZC(NODE1)
        X2 = XC(NODE2)
        Y2 = YC(NODE2)
        Z2 = ZC(NODE2)
        
        ! Element length and direction cosines
        DX = X2 - X1
        DY = Y2 - Y1
        DZ = Z2 - Z1
        EL = DSQRT(DX*DX + DY*DY + DZ*DZ)
        
        CX = DX / EL
        CY = DY / EL
        CZ = DZ / EL
        
        ! Displacements at nodes
        D1(1) = D((NODE1-1)*3 + 1)
        D1(2) = D((NODE1-1)*3 + 2)
        D1(3) = D((NODE1-1)*3 + 3)
        D2(1) = D((NODE2-1)*3 + 1)
        D2(2) = D((NODE2-1)*3 + 2)
        D2(3) = D((NODE2-1)*3 + 3)
        
        ! Axial elongation
        DELTA_L = CX * (D2(1) - D1(1)) + CY * (D2(2) - D1(2)) + &
                  CZ * (D2(3) - D1(3))
        
        ! Strain and stress
        STRAIN = DELTA_L / EL
        STRESS(K) = E(K) * STRAIN
    END DO
    
END SUBROUTINE COMPUTE_STRESS


!***********************************************************************
SUBROUTINE WRITE_OUTPUT(TITLE, NELE, NNODE, IFIX, XC, YC, ZC, FORCE, &
                        NODE, E, A, D, STRESS, MUD, MAXNODE, MAXELE)
!***********************************************************************
    IMPLICIT NONE
    
    CHARACTER(LEN=80), INTENT(IN) :: TITLE
    INTEGER, INTENT(IN) :: NELE, NNODE, MUD, MAXNODE, MAXELE
    INTEGER, INTENT(IN) :: IFIX(3, MAXNODE), NODE(2, MAXELE)
    REAL(8), INTENT(IN) :: XC(MAXNODE), YC(MAXNODE), ZC(MAXNODE)
    REAL(8), INTENT(IN) :: FORCE(3, MAXNODE), E(MAXELE), A(MAXELE)
    REAL(8), INTENT(IN) :: D(3*MAXNODE), STRESS(MAXELE)
    
    INTEGER :: J, K
    
    WRITE(*, '(A)') TRIM(TITLE)
    WRITE(*, *)
    WRITE(*, '(A,I3)') 'NUMBER OF ELEMENTS(NELE) = ', NELE
    WRITE(*, '(A,I3)') 'NUMBER OF NODES(NNODE)   = ', NNODE
    WRITE(*, *)
    WRITE(*, '(A)') 'NODE POINTS'
    WRITE(*, '(A)') 'K     IFIX        XC(K)           YC(K)           ZC(K)'
    
    DO J = 1, NNODE
        WRITE(*, '(I1,4X,I1,1X,I1,1X,I1,4X,E12.6,4X,E12.6,4X,E12.6)') &
              J, IFIX(1,J), IFIX(2,J), IFIX(3,J), XC(J), YC(J), ZC(J)
    END DO
    
    WRITE(*, *)
    WRITE(*, '(A)') '               FORCE(1,K)      FORCE(2,K)      FORCE(3,K)'
    DO J = 1, NNODE
        WRITE(*, '(14X,E12.6,4X,E12.6,4X,E12.6)') &
              FORCE(1,J), FORCE(2,J), FORCE(3,J)
    END DO
    
    WRITE(*, *)
    WRITE(*, '(A)') '   ELEMENTS'
    WRITE(*, '(A)') 'K      NODE(1,K)          E(K)               A(K)'
    DO K = 1, NELE
        WRITE(*, '(I1,8X,I1,3X,I1,9X,E12.4,8X,E12.4)') &
              K, NODE(1,K), NODE(2,K), E(K), A(K)
    END DO
    
    WRITE(*, *)
    WRITE(*, '(A,I3)') 'NUMBER OF NONZERO UPPER CODIAGONALS(MUD) = ', MUD
    WRITE(*, *)
    WRITE(*, '(A)') ' DISPLACEMENTS       X            Y            Z'
    
    DO J = 1, NNODE
        WRITE(*, '(A,I1,3X,E12.4,3X,E12.4,2X,E12.4)') &
              'NODE NUMBER ', J, D((J-1)*3+1), D((J-1)*3+2), D((J-1)*3+3)
    END DO
    
    WRITE(*, *)
    WRITE(*, '(A)') 'STRESSES IN ELEMENTS (IN CURRENT UNITS)'
    WRITE(*, *)
    WRITE(*, '(A)') 'ELEMENT NUMBER     STRESS'
    DO K = 1, NELE
        WRITE(*, '(10X,I1,A,4X,E12.5)') K, ' =', STRESS(K)
    END DO
    
END SUBROUTINE WRITE_OUTPUT
