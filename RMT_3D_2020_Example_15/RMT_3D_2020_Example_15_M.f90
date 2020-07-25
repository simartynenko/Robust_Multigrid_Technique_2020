!
! Evaluation of the integrals in the dynamic cycle
!
! Chapter IV, Section 3
!
! File: RMT_3D_2020_Example_15_M.f90    
!
program RMT_3D_2020_Example_15_M
USE     RMT_3D_2020_Example_15_1
USE     RMT_3D_2020_Example_15_2
USE     RMT_3D_UFG                                   
USE     RMT_3D_2020 
character*2 strX,strY,strZ

           NXX_FG = 100; N_coarsest_X = 5; LevelXd = 1; a_ = 1.d0  ! 001 
           NYY_FG = 100; N_coarsest_Y = 5; LevelYd = 1; b_ = 1.d0  ! 002
           NZZ_FG = 100; N_coarsest_Z = 5; LevelZd = 1; c_ = 1.d0  ! 003 

            call MSSize                                            ! 004
            call U_Finest_Grid(1,1,1)                              ! 005    
            call Coarse_Grids('+')                                 ! 006
            call Copy_Grid                                         ! 007
 
            call Starting_Guess_and_Boundary_Conditions            ! 008

   do     GridsXd = 1,3**LevelXd;      call Dynamic_Grids(1)       ! 009
    do    GridsYd = 1,3**LevelYd;      call Dynamic_Grids(2)       ! 010
     do   GridsZd = 1,3**LevelZd;      call Dynamic_Grids(3)       ! 011
         
                                    write(strX,'(i2)'   ) GridsXd         ! 012
                                       if(strX(1:1)==' ')   strX(1:1)='0' ! 012
                                       if(strX(2:2)==' ')   strX(2:2)='0' ! 012
                                    write(strY,'(i2)'   ) GridsYd         ! 012
                                       if(strY(1:1)==' ')   strY(1:1)='0' ! 012
                                       if(strY(2:2)==' ')   strY(2:2)='0' ! 012
                                    write(strZ,'(i2)'   ) GridsZd         ! 012
                                       if(strZ(1:1)==' ')   strZ(1:1)='0' ! 012
                                       if(strZ(2:2)==' ')   strZ(2:2)='0' ! 012
         
open(1,file='_GridsXd_'//strX(1:2)//' GridsYd_'//strY(1:2)//' GridsZd_'//strZ(1:2)//'.txt')   
                                       call Mini_Cycle             ! 013 
close(1)         
         
     end do                                                        ! 014
    end do                                                         ! 015
   end do                                                          ! 016

deallocate(MyF)   

end program RMT_3D_2020_Example_15_M