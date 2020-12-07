!Модуль для решения задачи обтекания пластины
!При решении применяется численное решение уравнений Навье-Стокса методом установления
!Входные параметры задаются в resource/Input.txt
!Решение выводиться в файл resource/data_ns.plt
MODULE MethodOfEstablishing_Plate
    implicit none

    CONTAINS

    SUBROUTINE MethodOfEstablishinglSolve_Plate()
        integer, parameter:: IO = 1, IO_Residuals = 2 ! input-output unit
        REAL, parameter :: Eps = 3e-5
        integer NI, NJ, NITER
        integer I,J, N
        real L,H,dx,dy, visk, U0, CFL
        real dt, A, U_Residuals, V_Residuals, P_Residuals
        real,allocatable :: X_Node(:,:),Y_Node(:,:)
        real,allocatable :: X_Cell(:,:),Y_Cell(:,:)
        real,allocatable :: U(:,:),V(:,:),P(:,:)
        real,allocatable :: U_n(:,:),V_n(:,:),P_n(:,:)
        real,allocatable :: U_half_i(:,:),V_half_i(:,:),P_half_i(:,:), U_half_j(:,:),V_half_j(:,:),P_half_j(:,:)

        write(*,*) 'Read input file'
        open(IO,FILE='source\resource\Input.txt')
        read(IO,*) L
        read(IO,*) H
        read(IO,*) NI
        read(IO,*) NJ
        read(IO,*) U0
        read(IO,*) visk
        read(IO,*) NITER
        read(IO,*) CFL
        close(IO)

        allocate(X_Node(NI,NJ)) ! mesh nodes X-coordinates
        allocate(Y_Node(NI,NJ)) ! mesh nodes Y-coordinates
        allocate(X_Cell(0:NI,0:NJ)) ! Cell Centers
        allocate(Y_Cell(0:NI,0:NJ)) ! Cell Centers

        !Cell-centered variables
        allocate(U(0:NI,0:NJ))   ! Velocity U
        allocate(V(0:NI,0:NJ))   ! Velocity V
        allocate(P(0:NI,0:NJ))   ! Pressure

        !Node variables
        allocate(U_n(0:NI,0:NJ))   ! Velocity U
        allocate(V_n(0:NI,0:NJ))   ! Velocity V
        allocate(P_n(0:NI,0:NJ))   ! Pressure

        !1/2 variables
        allocate(U_half_i(0:NI-1,0:NJ-1))   ! Velocity U
        allocate(V_half_i(0:NI-1,0:NJ-1))   ! Velocity V
        allocate(P_half_i(0:NI-1,0:NJ-1))   ! Pressure
        allocate(U_half_j(0:NI-1,0:NJ-1))   ! Velocity U
        allocate(V_half_j(0:NI-1,0:NJ-1))   ! Velocity V
        allocate(P_half_j(0:NI-1,0:NJ-1))   ! Pressure

        dx=L/(NI-1)
        dy=H/(NJ-1)
        A = 1/(U0**2)
        dt = CFL * dx / U0

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

    !Initial field
    call InitValue(U, V, P, NI, NJ, U0)
    call BoundValue(U, V, P, NI, NJ, U0)

    !Solve equation
    open(IO_Residuals,FILE='source/resource/residuals.dat', status = "replace")
    do N = 1, NITER
        !Заполнение половинных ячеек
        do J = 0, NJ-1
            do I = 0, NI-1
                U_half_i(I,J) = UH(U(I,J), U(I+1,J), I, NI)
                V_half_i(I,J) = UH(V(I,J), V(I+1,J), I, NI)
                U_half_j(I,J) = UH(U(I,J), U(I,J+1), J, NJ)
                V_half_j(I,J) = UH(V(I,J), V(I,J+1), J, NJ)
                P_half_i(I,J) = UP(P(I+1,J), P(I,J), U(I,J), U(I+1,J), I, NI)
                P_half_j(I,J) = UP(P(I+1,J), P(I,J), V(I,J), V(I,J+1), J, NJ)
            end do
        end do
        do J = 1, NJ-1                                                          !Мэйби коректирование индексов или разбивание циклов
            do I = 1, NI-1
                !Вычисление давления
                P_n(I,J) = P(I,J) - (dt/A)*((U_half_i(I,J) - U_half_i(I-1,J))/dx + (V_half_j(I,J) - V_half_j(I,J-1))/dy)

                !Высление продольной компоненты скорости
                U_n(I,J) = U(I,J) - dt*( ( (U(I,J) + U(I+1,J))*U_half_i(I,J)/2 - (U(I,J) + U(I-1,J))*U_half_i(I-1,J)/2 )/dx &
                & + ( (V(I,J) + V(I,J+1))*U_half_j(I,J)/2 - (V(I,J) + V(I,J-1))*U_half_j(I,J-1)/2 )/dy &
                & + (P_half_i(I,J) - P_half_i(I-1,J))/dx &
                & -(1/dx)*(visk * (U(I+1,J) - U(I,J))/dx - visk*(U(I,J) - U(I-1,J))/dx) &
                & -(1/dy)*(visk * (U(I,J+1) - U(I,J))/dy - visk*(U(I,J) - U(I,J-1))/dy))

                !Вычисление поперечной компоненты скорости
                V_n(I,J) = V(I,J) - dt*( ( (V(I,J) + V(I,J+1))*V_half_j(I,J)/2 - (V(I,J) + V(I,J-1))*V_half_j(I,J-1)/2 )/dy &
                & + ( (U(I,J) + U(I+1,J))*V_half_i(I,J)/2 - (U(I,J) + U(I-1,J))*V_half_i(I-1,J)/2 )/dx &
                & + (P_half_j(I,J) - P_half_j(I,J-1))/dy &
                & - (1/dx)*(visk * (V(I+1,J) - V(I,J))/dx - visk*(V(I,J) - V(I-1,J))/dx) &
                & - (1/dy)*(visk * (V(I,J+1) - V(I,J))/dy - visk*(V(I,J) - V(I,J-1))/dy))

                !Пересчет в граничных ячейках
                call BoundValue(U_n, V_n, P_n, NI, NJ, U0)

            end do
        end do

        !Проверяем сходимость и выводим невязки
        U_Residuals = maxval(abs(U_n-U))/maxval(abs(U_n))
        V_Residuals = maxval(abs(V_n-V))/maxval(abs(V_n))
        P_Residuals = maxval(abs(P_n-P))/maxval(abs(P_n))
        if( (U_Residuals.le.Eps ) .and. (V_Residuals.le.Eps ) .and. (P_Residuals.le.Eps ) ) then
            write(*,*) "MethodOfEstablishinglSolve_Plate:Complete"
            exit
        endif
        write(*,*) "N=", N
        write(IO_Residuals, *) dt*N, U_Residuals, V_Residuals, P_Residuals

        !Переопределяем для следующего шага
        U = U_n
        V = V_n
        P = P_n
    end do
    close(IO_Residuals)

    call writeAnswer(IO,NI,NJ,X_Cell,Y_Cell,U_n,V_n,P_n)
    return
    END SUBROUTINE MethodOfEstablishinglSolve_Plate

    !Функция вычисления давления в 1/2
    REAL FUNCTION UP(P1, P2, U1, U2, J, NJ)
        real :: U1, U2, UCap, P1, P2
        integer:: J, NJ
        Ucap = (U1 + U2)/2
        if( (J .eq. NJ-1) .or. (J .eq. 1)) then
            UP = Ucap
            return
        end if
        if(Ucap .ge. 0) then
            UP = P1
            return
        end if
        UP = P2
    END FUNCTION UP

    !Функция вычисления скорости в 1/2
    REAL FUNCTION UH(U1, U2, J, NJ)
        real :: U1, U2, UCap
        integer:: J, NJ
        Ucap = (U1 + U2)/2
        if( (J .eq. NJ-1) .or. (J .eq. 1)) then
            UH = Ucap
            return
        end if
        if(Ucap .ge. 0) then
            UH = U1
            return
        end if
        UH = U2
    END FUNCTION UH

    !Функция пишет ответ в файл
    SUBROUTINE writeAnswer(IO,NI,NJ,X,Y,U,V,P)
       implicit none

       integer NI,NJ,IO
       real, dimension(NI,NJ):: X,Y
       real, dimension(0:NI,0:NJ)::U,V,P

       write(*,*) 'Output data cell (Navier - Stokes)   '
       open(IO,FILE='source/resource/data_ns.plt', status = "replace")
       call OutputFields_Cell(IO,NI,NJ,X,Y,U,V,P)
       close(IO)
    END SUBROUTINE writeAnswer

    !Функция для вывода в формате техплот
    SUBROUTINE OutputFields_Cell(IO,NI,NJ,X,Y,U,V,P)
        implicit none

        integer NI,NJ,IO
        real, dimension(NI,NJ):: X,Y
        real, dimension(0:NI,0:NJ)::U,V,P

        write(IO,*) 'VARIABLES = "X", "Y", "U", "V", "P"'
        write(IO,*) 'ZONE I=',NI,', J=',NJ,', DATAPACKING=BLOCK, VARLOCATION=([3-5]=CELLCENTERED)'
        write(IO,'(100E25.16)') X(1:NI,1:NJ)
        write(IO,'(100E25.16)') Y(1:NI,1:NJ)
        write(IO,'(100E25.16)') U(1:NI-1,1:NJ-1)
        write(IO,'(100E25.16)') V(1:NI-1,1:NJ-1)
        write(IO,'(100E25.16)') P(1:NI-1,1:NJ-1)

    END SUBROUTINE OutputFields_Cell

    SUBROUTINE BoundValue(U,V,P,NI,NJ,U0)
        implicit none

        integer NI,NJ
        real :: U(0:NI,0:NJ), V(0:NI,0:NJ), P(0:NI,0:NJ)
        real U0
        !нижняя непроницаемая граница
        U(1:NI-1,0) = - U(1:NI-1,1)
        V(1:NI-1,0) = - V(1:NI-1,1)
        P(1:NI-1,0) = P(1:NI-1,1)

        !левая граница, вход
        U(0,1:NJ-1) =  U0
        V(0,1:NJ-1) = 0
        P(0,1:NJ-1) = P(1,1:NJ-1)

        !правая граница, выход
        U(NI,1:NJ-1) = U(NI-1,1:NJ-1)
        V(NI,1:NJ-1) = V(NI-1,1:NJ-1)
        P(NI,1:NJ-1) = 0

        !верхняя граница, выход
        U(1:NI-1,NJ) = U(1:NI-1,NJ-1)
        V(1:NI-1,NJ) = V(1:NI-1,NJ-1)
        P(1:NI-1,NJ) = 0
    END SUBROUTINE BoundValue

    SUBROUTINE InitValue(U,V,P,NI,NJ,U0)
        implicit none

        real :: U(1:NI-1,1:NJ-1), V(1:NI-1,1:NJ-1),P(1:NI-1,1:NJ-1)
        real U0
        integer NI,NJ

        U(1:NI-1,1:NJ-1) = U0
        V(1:NI-1,1:NJ-1) = 0
        P(1:NI-1,1:NJ-1) = 0

    END SUBROUTINE InitValue

END MODULE MethodOfEstablishing_Plate
