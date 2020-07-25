!
! Test of the index mapping in the dynamic cycle
!
! Chapter IV Section 2
!
! File: RMT_3D_2020_Example_14.f90
!    
program RMT_3D_2020_Example_14
USE     RMT_3D_2020                                                   ! 001    
USE     RMT_3D_UFG                                                    ! 002

        NXX_FG =  29; N_coarsest_X = 5                                ! 003
        NYY_FG = 187; N_coarsest_Y = 2                                ! 004
        NZZ_FG = 335; N_coarsest_Z = 2                                ! 005

        call MSSize                                                   ! 006
        call U_Finest_Grid(1,2,1)                                     ! 007
        call Coarse_Grids('+')                                        ! 008
        call Copy_Grid                                                ! 009
            
        open(1,file='RMT_3D_2020_Example_14.res')                     ! 010
                                     write(1,1)                       ! 011
           LevelXd = 1; GridsXd = 1; call Dynamic_Grids(1)            ! 012        
           LevelX  = 0; GridsX  = 1; call Index_Mapping('-',1)        ! 013
            do   i = 1,Nx+1                                           ! 014
                                     write(1,2) i,TxV(i),TxV_(TxV(i)) ! 015
            end do                                                    ! 016 
        close(1)                                                      ! 017 
1 format('     CG    DG     FG')        
2 format(5x,i2,' -> ',i2,' ->  ',i2)        
end program RMT_3D_2020_Example_14
