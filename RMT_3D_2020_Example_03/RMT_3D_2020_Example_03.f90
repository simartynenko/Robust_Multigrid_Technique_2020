!
! Multigrid structure generation in unit cube for the uniform grids 
!
! Chapter II, Section 3
!
! File: RMT_3D_2020_Example_03.f90
!    
program RMT_3D_2020_Example_03
USE     RMT_3D_2020                                ! 001
USE     RMT_3D_UFG                                 ! 002

        NXX_FG = 10; N_coarsest_X = 2; LevelXd = 0 ! 003
        NYY_FG = 10; N_coarsest_Y = 2; LevelYd = 0 ! 004
        NZZ_FG = 10; N_coarsest_Z = 2; LevelZd = 0 ! 005
        
        call MSSize                                ! 006
        call U_Finest_Grid(1,1,1)                  ! 007
        call Coarse_Grids('-')                     ! 008
        call Copy_Grid                             ! 009

end program RMT_3D_2020_Example_03