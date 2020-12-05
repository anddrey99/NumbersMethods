module TestAllModules
    use SweepMethods_Test
    implicit none

    contains

    subroutine test_all()
        write(*,*) 'Тесты начаты'
        call progonka_test()
        Write(*,*) 'Тесты завершены'
        return
    end subroutine

end module TestAllModules
