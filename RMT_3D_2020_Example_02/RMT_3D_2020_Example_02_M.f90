!
! Block Gauss-Seidel iterations for solving Poisson equation (1.6)
!
! Chapter I, Section 3
!
! File: RMT_3D_2020_Example_02_M.f90
!    
program RMT_3D_2020_Example_02_M
USE     RMT_3D_2020_Example_02_1    
USE     RMT_3D_2020_Example_02_2

real*8  A_W,A_E,A_S,A_N,A_D,A_U,A_P
integer bGS_IterMAX

            NXX_FG = 50; NX_block = 3; a_ = 1.d0                                    ! 001 
            NYY_FG = 50; NY_block = 3; b_ = 1.d0                                    ! 002
            NZZ_FG = 50; NZ_block = 3; c_ = 1.d0                                    ! 003
       bGS_IterMAX = 1000000;     ItCheck = 400                                     ! 004
            ResLIM = 1.d-7                                                          ! 005
            
                Nx = NXX_FG; h_x = 1.d0/dfloat(Nx)                                  ! 006  
                Ny = NYY_FG; h_y = 1.d0/dfloat(Ny)                                  ! 007
                Nz = NZZ_FG; h_z = 1.d0/dfloat(Nz)                                  ! 008                

open(1,file='RMT_3D_2020_Example_02.res')
call MSSize                                                                         ! 009  
call Starting_Guess_and_Boundary_Conditions                                         ! 010

do       Iteration = 1,bGS_IterMAX                                                  ! 011

 do       I_block  = 1,X_block                                                      ! 012
                                I_top =(I_block-1)*NX_block+1                       ! 013                   
                                I_bot = I_block   *NX_block                         ! 014  
       if(I_block ==   X_block) I_top = Nx+2      -NX_block                         ! 015 
       if(I_block ==   X_block) I_bot = Nx+1                                        ! 016 

 do       J_block  = 1,Y_block                                                      ! 017
                                J_top =(J_block-1)*NY_block+1                       ! 018
                                J_bot = J_block   *NY_block                         ! 019
       if(J_block ==   Y_block) J_top = Ny+2      -NY_block                         ! 020
       if(J_block ==   Y_block) J_bot = Ny+1                                        ! 021

 do       K_block  = 1,Z_block                                                      ! 022
                                K_top =(K_block-1)*NZ_block+1                       ! 023
                                K_bot = K_block   *NZ_block                         ! 024
       if(K_block ==   Z_block) K_top = Nz+2      -NZ_block                         ! 025
       if(K_block ==   Z_block) K_bot = Nz+1                                        ! 026
       
                AG = 0.D0;         xG = 0.D0                                        ! 027
   allocate(M_A_P(I_top:I_bot,J_top:J_bot,K_top:K_bot))                             ! 028
                                   m  = 0                                           ! 029
   do       i      = I_top,I_bot                                                    ! 030  
    do        j    = J_top,J_bot                                                    ! 031  
     do         k  = K_top,K_bot;  m  = m+1                                         ! 032
      M_A_P(i,j,k) =               m                                                ! 033    
     end do                                                                         ! 034
    end do                                                                          ! 035 
   end do                                                                           ! 036   
                                   N_ = m                                           ! 037
                                   m  = 0                                           ! 038
   do       i      = I_top,I_bot                                                    ! 039
    do        j    = J_top,J_bot                                                    ! 040
     do          k = K_top,K_bot                                                    ! 041
                                   m  = m+1                                         ! 042   

               A_W = 0.D0; A_E = 0.D0; A_S = 0.D0; A_N = 0.D0                       ! 043
               A_D = 0.D0; A_U = 0.D0; A_P = 0.D0;                    way = 0       ! 044

! Plane x = 0 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(i==1)                                                                 then ! 045
                             AG(m,M_A_P(i,j,k)) = 1.D0;               way = 1       ! 046 
                             AG(m,N_+1)         = BC_0YZ(j,k)                       ! 047
      end if                                                                        ! 048
! Plane x = 1 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(i==Nx+1)                                                              then ! 049 
                             AG(m,M_A_P(i,j,k)) = 1.D0;               way = 1       ! 050 
                             AG(m,N_+1)         = BC_1YZ(j,k)                       ! 051 
      end if                                                                        ! 052 
! Plane y = 0 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
      if(j==1)                                                                 then ! 053 
                             AG(m,M_A_P(i,j,k)) = 1.D0;               way = 1       ! 054 
                             AG(m,N_+1)         = BC_X0Z(i,k)                       ! 055 
      end if                                                                        ! 056 
! Plane y = 1 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
      if(j==Ny+1)                                                              then ! 057 
                             AG(m,M_A_P(i,j,k)) = 1.D0;               way = 1       ! 058 
                             AG(m,N_+1)         = BC_X1Z(i,k)                       ! 059 
      end if                                                                        ! 060 
! Plane z = 0 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
      if(k==1)                                                                 then ! 061 
                             AG(m,M_A_P(i,j,k)) = 1.D0;               way = 1       ! 062
                             AG(m,N_+1)         = BC_XY0(i,j)                       ! 063
      end if                                                                        ! 064
! Plane z = 1 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(k==Nz+1)                                                              then ! 065
                             AG(m,M_A_P(i,j,k)) = 1.D0;               way = 1       ! 066
                             AG(m,N_+1)         = BC_XY1(i,j)                       ! 067
      end if                                                                        ! 068

      if(way == 0)                                                             then ! 069    
               A_W =       1.d0/h_x**2                                              ! 070
               A_E =       1.d0/h_x**2                                              ! 071
               A_S =                     1.d0/h_y**2                                ! 072
               A_N =                     1.d0/h_y**2                                ! 073
               A_U =                                   1.d0/h_z**2                  ! 074 
               A_D =                                   1.d0/h_z**2                  ! 075
               A_P =-2.d0*(1.d0/h_x**2 + 1.d0/h_y**2 + 1.d0/h_z**2)                 ! 076

      if(i-1 >= I_top)  AG(m,M_A_P(i-1,j  ,k  )) = A_W                              ! 077   
      if(i-1 >= I_top)                             A_W = 0.d0                       ! 078
	     
      if(i+1 <= I_bot)  AG(m,M_A_P(i+1,j  ,k  )) = A_E                              ! 079
      if(i+1 <= I_bot)                             A_E = 0.d0                       ! 080
	                                                                                   
      if(j-1 >= J_top)  AG(m,M_A_P(i  ,j-1,k  )) = A_S                              ! 081  
      if(j-1 >= J_top)                             A_S = 0.d0                       ! 082  
	                                                                                   
      if(j+1 <= J_bot)  AG(m,M_A_P(i  ,j+1,k  )) = A_N                              ! 083  
      if(j+1 <= J_bot)                             A_N = 0.d0                       ! 084  
	                                                                                   
      if(k-1 >= K_top)  AG(m,M_A_P(i  ,j  ,k-1)) = A_D                              ! 085  
      if(k-1 >= K_top)                             A_D = 0.d0                       ! 086  
	                                                                                   
      if(k+1 <= K_bot)  AG(m,M_A_P(i  ,j  ,k+1)) = A_U                              ! 087  
      if(k+1 <= K_bot)                             A_U = 0.d0                       ! 088  
	   
                        AG(m,M_A_P(i  ,j  ,k  )) = A_P                              ! 089
                        AG(m,N_+1)               =-f(i,j,k) &                       ! 090        
                  - A_W*U(i-1,j  ,k  ) - A_E*U(i+1,j  ,k  ) &                       ! 090
                  - A_S*U(i  ,j-1,k  ) - A_N*U(i  ,j+1,k  ) &                       ! 090
                  - A_D*U(i  ,j  ,k-1) - A_U*U(i  ,j  ,k+1)                         ! 090

      end if                                                                        ! 091
     end do                                                                         ! 092 
    end do                                                                          ! 093
   end do                                                                           ! 094
                        call GAUSSIAN_ELIMINATION(N_)                               ! 095 
                                   m  = 0                                           ! 096 
   do       i      = I_top,I_bot                                                    ! 097   
    do        j    = J_top,J_bot                                                    ! 098
     do         k  = K_top,K_bot;  m  = m+1                                         ! 099 
          U(i,j,k) = xG(m)                                                          ! 100
     end do                                                                         ! 101
    end do                                                                          ! 102
   end do                                                                           ! 103  
                     deallocate(M_A_P)                                              ! 104 
  end do                                                                            ! 105 
  end do                                                                            ! 106 
 end do                                                                             ! 107  
 
 if((Iteration / ItCheck)*ItCheck==Iteration)                                  then ! 108  
   call Convergence_Test(way)                                                       ! 109
   if(way==1) exit                                                                  ! 110
 end if                                                                             ! 111 
 
 end do                                                                             ! 112

close(1)    
end program RMT_3D_2020_Example_02_M    