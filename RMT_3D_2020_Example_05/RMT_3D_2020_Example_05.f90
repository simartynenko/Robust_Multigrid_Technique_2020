!
! Test of the index mapping
!
! Chapter II, Section 5
!
! File: RMT_3D_2020_Example_05.f90
!    
program RMT_3D_2020_Example_05
USE     RMT_3D_2020                                                            ! 001    
USE     RMT_3D_UFG                                                             ! 002

        NXX_FG =  10; N_coarsest_X = 7; LevelXd = 0                            ! 003
        NYY_FG = 187; N_coarsest_Y = 2; LevelYd = 0                            ! 004
        NZZ_FG = 335; N_coarsest_Z = 2; LevelZd = 0                            ! 005

        call MSSize                                                            ! 006
        call U_Finest_Grid(1,2,1)                                              ! 007
        call Coarse_Grids('+')                                                 ! 008
        call Copy_Grid                                                         ! 009
            
        open(1,file='RMT_3D_2020_Example_05.res')
            LevelX = 1                                                         ! 010
         do GridsX = 1,3**LevelX;  call Index_Mapping('-',1)                   ! 011
                                   write(1,1) GridsX, Nx
                                   write(1,2) 
           do    i = 1,Nx+1;       write(1,3) i,TxV(i);        end do          ! 012
                                   write(1,4)           
           if(TxF(0)>0)            write(1,3) 0,TxF(0)                         ! 013
           do    i = 1,Nx  ;       write(1,3) i,TxF(i);        end do          ! 014
           if(TxF(Nx+1)<NXX_FG+1)  write(1,3) Nx+1,TxF(Nx+1)                   ! 015
         end do                                                                ! 016    
        close(1)
        
1 format(//5x,'GridsX = ',i1,'    Nx =',i2,'  = = = =')            
2 format( /5x,'Mapping of the vertex indices:   Xv -> Xv'/)      
3 format(  5x,i2,' -> ',i2)                        
4 format( /5x,'Mapping of the face indices:     Xf -> Xf'/)             
    
end program RMT_3D_2020_Example_05
