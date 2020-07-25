!
! Multigrid structure generation in unit cube for the nonuniform grids 
!
! Chapter II, Section 6
!
! File: RMT_3D_2020_Example_08.f90
!    
program RMT_3D_2020_Example_08
USE     RMT_3D_2020
USE     RMT_3D_NFG2                    ! 001  

        NXX_FG = 50; N_coarsest_X = 2; LevelXd = 0
        NYY_FG = 50; N_coarsest_Y = 2; LevelYd = 0
        NZZ_FG = 50; N_coarsest_Z = 2; LevelZd = 0

        call MSSize                    
        call Finest_Grid               ! 002
        call Coarse_Grids('+')
        call Copy_Grid
        
end program RMT_3D_2020_Example_08