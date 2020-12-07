!Модуль для решения задачи обтекания пластины
!При решении применяется численное решение уравнений Прандтля
!Входные параметры задаются в resource/Input.txt
!Решение выводиться в файл resource/data_pt.plt

module  PrandtlSolver_Plate
    use SweepMethods
    IMPLICIT NONE
    contains

    subroutine PrandtlSolve_Plate()
        INTEGER, PARAMETER:: IO = 12 ! input-output unit
        REAL, PARAMETER :: Eps = 3e-5
        INTEGER NI, NJ
        INTEGER I,J, S
        REAL L,H,dx,dy, visk, U0
        REAL,ALLOCATABLE :: X_Node(:,:),Y_Node(:,:)
        REAL,ALLOCATABLE :: X_Cell(:,:),Y_Cell(:,:)
        REAL,ALLOCATABLE :: U_c(:,:),V_c(:,:),P_c(:,:)
        REAL,ALLOCATABLE :: U_n(:,:),V_n(:,:),P_n(:,:)
        real, allocatable :: A(:), B(:), C(:), D(:)

        write(*,*) 'Read input file'
        open(IO,FILE='source\resource\Input.txt')
        read(IO,*) L
        read(IO,*) H
        read(IO,*) NI
        read(IO,*) NJ
        read(IO,*) U0
        read(IO,*) visk
        close(IO)

        allocate(X_Node(NI,NJ)) ! mesh nodes X-coordinates
        allocate(Y_Node(NI,NJ)) ! mesh nodes Y-coordinates
        allocate(X_Cell(0:NI,0:NJ)) ! Cell Centers
        allocate(Y_Cell(0:NI,0:NJ)) ! Cell Centers

        allocate(A(1:NJ))
        allocate(B(1:NJ))
        allocate(C(1:NJ))
        allocate(D(1:NJ))

    !*******************  Cell-centered variables **********
        allocate(U_c(0:NI,0:NJ))   ! Velocity U
        allocate(V_c(0:NI,0:NJ))   ! Velocity V
        allocate(P_c(0:NI,0:NJ))   ! Pressure

    !*******************  Node variables ******************
        allocate(U_n(NI,NJ))   ! Velocity U
        allocate(V_n(NI,NJ))   ! Velocity V
        allocate(P_n(NI,NJ))   ! Pressure

        dx=L/(NI-1)
        dy=H/(NJ-1)

        do I=1,NI
          do J=1,NJ
            X_Node(I,J)=(I-1)*dx
            Y_Node(I,J)=(J-1)*dy
          end do
        end do

        X_Cell(0,1:NJ)=-dx/2
        Y_Cell(0,1:NJ)=Y_Node(1,1:NJ)+dy/2
        X_Cell(1:NI,0)=X_Node(1:NI,1)+dx/2
        Y_Cell(1:NI,0)=-dy/2
        do I=1,NI
          do J=1,NJ
            X_Cell(I,J)=X_Node(I,J)+dx/2
            Y_Cell(I,J)=Y_Node(I,J)+dy/2
          end do
        end do

    !************************* INITIAL FIELD *********************

        U_c = 0.
        V_c = 0.
        U_n = 0.
        V_n = 0.
        call InitValue(U_c,NI,NJ,U0)
        call BoundValue(U_c,V_c,NI,NJ,U0)
        call InitValue(U_n,NI,NJ,U0)
        call BoundValue(U_n,V_c,NI,NJ,U0)

    !****************** Solve equation ********************
        do I = 2, NI
            U_c(I,1:NJ) = U_c(I-1,1:NJ)
            V_c(I,1:NJ) = V_c(I-1,1:NJ)
            S = 0
            do
                S = S + 1
                A(1) = 0
                B(1) = 1
                C(1) = 0
                D(1) = 0

                A(NJ) = 0
                B(NJ) = 1
                C(NJ) = 0
                D(NJ) = U0

                do J=2,NJ-1
                        A(J) = -V_c(I,J-1)/(2*dy) - visk/(dy**2)
                        B(J) = U_c(I,J)/dx + 2*visk/(dy**2)
                        C(J) = V_c(I,J+1)/(2*dy) - visk/(dy**2)
                        D(J) = U_n(I-1,J)**2 / dx
                enddo

                call progonka(A,B,C,D,NJ,U_n(I,1:NJ))

                V_n(I,1) = 0.
                do J=2,NJ
                        V_n(I,J) = V_n(I,J-1) - (dy/(2.0 * dx)) * (U_n(I,J) - U_n(I-1,J) + U_n(I,J-1) - U_n(I-1,J-1))
                enddo

                If (((maxval(abs(U_n(I,1:NJ)-U_c(I,1:NJ)))/maxval(abs(U_n(I,1:NJ)))).LE.Eps).and.&
                    &((maxval(abs(V_n(I,1:NJ)-V_c(I,1:NJ)))/maxval(abs(V_n(I,1:NJ)))).LE.Eps)) then
                        Write(*,*) "s = ", s,  "U_n(", I, "," , NJ, ")=", U_n(I,NJ)
                        exit
                endif

                If (S>1000) then
                    write(*,*) "error s ","I=",I
                    Write(*,*) "errotU", maxval(abs(U_n(I,1:NJ)-U_c(I,1:NJ)))/maxval(abs(U_n(I,1:NJ)))
                    Write(*,*) "errotV", maxval(abs(V_n(I,1:NJ)-V_c(I,1:NJ)))/maxval(abs(V_n(I,1:NJ)))
                    stop
                endif

                U_c=U_n
                V_c=V_n

            enddo
        enddo

    !****************** Output Results ********************

        write(*,*) 'Output data node (Prandtl)'
        Open(IO,FILE='source/resource/data_pr.tec')
        Call OutputFields_Node(IO,NI,NJ,X_Node,Y_Node,U_n,V_n,P_n)
        Close(IO)
        return

    end subroutine

    SUBROUTINE OutputFields_Cell(IO,NI,NJ,X,Y,U,V,P)
        IMPLICIT NONE

        INTEGER NI,NJ,IO
        REAL, DIMENSION(NI,NJ):: X,Y
        REAL, DIMENSION(0:NI,0:NJ)::U,V,P

        Write(IO,*) 'VARIABLES = "X", "Y", "U", "V", "P"'
        Write(IO,*) 'ZONE I=',NI,', J=',NJ,', DATAPACKING=BLOCK, VARLOCATION=([3-5]=CELLCENTERED)'
        Write(IO,'(100E25.16)') X(1:NI,1:NJ)
        Write(IO,'(100E25.16)') Y(1:NI,1:NJ)
        Write(IO,'(100E25.16)') U(1:NI-1,1:NJ-1)
        Write(IO,'(100E25.16)') V(1:NI-1,1:NJ-1)
        Write(IO,'(100E25.16)') P(1:NI-1,1:NJ-1)

    END SUBROUTINE

    SUBROUTINE OutputFields_Node(IO,NI,NJ,X,Y,U,V,P)
        IMPLICIT NONE

        INTEGER NI,NJ,IO
        REAL, DIMENSION(NI,NJ):: X,Y
        REAL, DIMENSION(NI,NJ):: U,V,P

        Write(IO,*) 'VARIABLES = "X", "Y", "U", "V", "P"'
        Write(IO,*) 'ZONE I=',NI,', J=',NJ, ', DATAPACKING=BLOCK'
        Write(IO,'(100E25.16)') X(1:NI,1:NJ)
        Write(IO,'(100E25.16)') Y(1:NI,1:NJ)
        Write(IO,'(100E25.16)') U(1:NI,1:NJ)
        Write(IO,'(100E25.16)') V(1:NI,1:NJ)
        Write(IO,'(100E25.16)') P(1:NI,1:NJ)

    END  SUBROUTINE

    !----------------------- Set Boundary Condition ------------
    subroutine BoundValue(U,V,NI,NJ,U0)
        IMPLICIT NONE

        INTEGER NI,NJ
        REAL :: U(1:NI,1:NJ), V(1:NI,1:NJ)
        REAL U0

        U(1:NI,1) = 0.
        V(1:NI,1) = 0.
        U(1:NI,NJ) = U0

    end subroutine

    !---------------Set Initial Values---------------------
    subroutine InitValue(U,NI,NJ,U0)
        IMPLICIT NONE

        Real :: U(1:NI,1:NJ)
        Real U0
        INTEGER NI,NJ

        U(1,1:NJ) = U0

    end subroutine
End module
