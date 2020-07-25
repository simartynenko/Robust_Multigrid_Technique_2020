!
! First multigrid iteration in the dynamic cycle
!
! Chapter IV, Section 4
!
! File: RMT_3D_2020_Example_16_M.f90    
!
program RMT_3D_2020_Example_16_M
USE     RMT_3D_2020_Example_16_1
USE     RMT_3D_2020_Example_16_2
USE     RMT_3D_2020 
USE     RMT_3D_UFG 

    NXX_FG = 270; N_coarsest_X = 5; LevelXd = 1; a_ = 1.d0 ! 001 
    NYY_FG = 270; N_coarsest_Y = 5; LevelYd = 1; b_ = 1.d0 ! 002
    NZZ_FG = 270; N_coarsest_Z = 5; LevelZd = 1; c_ = 1.d0 ! 003
           
    The_number_of_smoothing_iterations = 6                 ! 004
    The_number_of_multigrid_iterations = 1                 ! 005
                                MGI_DC = 5                 ! 006
             call MSSize                                   ! 007
             call U_Finest_Grid(1,1,1)                     ! 008    
             call Coarse_Grids('+')                        ! 009
             call Copy_Grid                                ! 010
open(1,file='RMT_3D_2020_Example_16.res')                  ! 011 
             call Starting_Guess_and_Boundary_Conditions   ! 012
    do     GridsXd = 1,3**LevelXd; call Dynamic_Grids(1)   ! 013
     do    GridsYd = 1,3**LevelYd; call Dynamic_Grids(2)   ! 014
      do   GridsZd = 1,3**LevelZd; call Dynamic_Grids(3)   ! 015
                                   call Mini_Cycle         ! 016 
      end do                                               ! 017
     end do                                                ! 018
    end do                                                 ! 019
                                   call Convergence_Test   ! 020
close(1)                                                   ! 021
end program RMT_3D_2020_Example_16_M