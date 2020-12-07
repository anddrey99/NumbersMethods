module TestAllModules
    use SweepMethods_Test
    use MethodOfEstablishing_Plate_Test
    implicit none

    private :: refresh_file
    public :: test_all

    contains

    subroutine test_all()
        call refresh_file()
        write(*,*) 'Тесты начаты'
        call progonka_test()
        call UH_Test_GU()
        call UH_Test_NoGU()
        Write(*,*) 'Тесты завершены'
        return
    end subroutine

    subroutine refresh_file()
        open(1, file = "test/testlogs.txt", status = "replace")
        close(1)
    end subroutine


end module TestAllModules
