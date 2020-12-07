MODULE MethodOfEstablishing_Plate_Test
    use MethodOfEstablishing_Plate
    implicit none
    real, parameter :: eps = 0.000001
    integer, parameter :: logs = 1
    CONTAINS

    SUBROUTINE UH_Test_GU()
        open(logs, file = "test/testlogs.txt", status = "old")
        if((UH(2.,10.,1,10) .eq. 6.) .and. (UH(2.,10.,9,10).eq.6.) ) then
            write(*,*) ' MethodOfEstablishing_Plate_Test:UH_Test_GU::complete'
        else
            write(logs,*) ' MethodOfEstablishing_Plate_Test:UH_Test_GU::Ошибка в граничных условиях'
            write(*,*) ' MethodOfEstablishing_Plate_Test:UH_Test_GU::failed'
        end if
        close(logs)
    END SUBROUTINE UH_Test_GU

    SUBROUTINE UH_Test_NoGU()
        open(logs, file = "test/testlogs.txt", status = "old")
        if((UH(2.,10.,4,10) .eq. 2.) .and. (UH(2.,-10.,4,10).eq. -10.) .and. (UH(10.,-10.,4,10).eq. 10.)) then
            write(*,*) ' MethodOfEstablishing_Plate_Test:UH_Test_NoGU::complete'
        else
            write(logs,*) ' MethodOfEstablishing_Plate_Test:UH_Test_NoGU::Ошибка в условиях на середине'
            write(*,*) ' MethodOfEstablishing_Plate_Test:UH_Test_NoGU::failed'
        end if
        close(logs)
    END SUBROUTINE UH_Test_NoGU

END MODULE MethodOfEstablishing_Plate_Test
