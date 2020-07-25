!
! Approximation of the second derivative on uniform grid
!
! Chapter II, Section 5 
!
! File: RMT_3D_2020_Example_06.f90
!    
program RMT_3D_2020_Example_06
USE     RMT_3D_2020                                                            ! 001
USE     RMT_3D_UFG                                                             ! 002
real*8  DFun(1000),DF,Error,ErrMAX

            NXX_FG = 123; N_coarsest_X = 4; LevelXd = 0                        ! 003
            NYY_FG = 187; N_coarsest_Y = 2; LevelYd = 0                        ! 004
            NZZ_FG = 335; N_coarsest_Z = 2; LevelZd = 0                        ! 005 

            call MSSize
            call U_Finest_Grid(1,2,1)
            call Coarse_Grids('+')
            
              do i = 1,NXX_FG+1; Dfun(i) = dexp(Xv_FG(i)); end do              ! 006    
            
      open(1,file='RMT_3D_2020_Example_06.res')            
      do    LevelX = LeXmaxFG,0,-1                                             ! 007
                                   write(1,1) LevelX
       do   GridsX = 1,3**LevelX;  call Index_Mapping('-',1)                   ! 008 
         ErrMAX    = 0.d0                                                      ! 009   
         do      i = 2,Nx                                                      ! 010  
                DF =(Dfun(TxV(i-1))-2.d0*Dfun(TxV(i))+Dfun(TxV(i+1))) &        ! 011
                   /(h_x*dfloat(3**LevelX))**2                                 ! 011      
             Error = dabs(DF - dexp(Xv_FG(TxV(i))))                            ! 012
          if(Error > ErrMAX) ErrMAX = Error                                    ! 013
         end do                                                                ! 014 
                                   write(1,2) GridsX,Nx,ErrMAX                 ! 015  
       end do                                                                  ! 016 
      end do                                                                   ! 017 
      close(1)
        
1 format(//5x,'LevelX = ',i1/)
2 format(  5x,'GridsX = ',i2,' Nx =',i3,'  ErrMAX = ',D12.5) 
  
end program RMT_3D_2020_Example_06
