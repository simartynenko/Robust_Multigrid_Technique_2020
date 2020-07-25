MODULE  RMT_3D_2020_Example_13_2
USE     RMT_3D_2020_Example_13_1
USE     RMT_3D_2020

PRIVATE

PUBLIC :: Starting_Guess_and_Boundary_Conditions, Matrix_and_Vector
PUBLIC :: Convergence_Test, Vanka_type_Smoother
    
CONTAINS

subroutine  Starting_Guess_and_Boundary_Conditions                                       
       corU = 0.d0;   hatU = 0.d0;     RHSF = 0.d0;  J3 = 0.d0                           ! 001
     BC_0YZ = 0.d0; BC_1YZ = 0.d0; Lambda_X = 1.d0                                       ! 002
     BC_X0Z = 0.d0; BC_X1Z = 0.d0; Lambda_Y = 1.d0                                       ! 003   
     BC_XY0 = 0.d0; BC_XY1 = 0.d0; Lambda_Z = 1.d0                                       ! 004
   
                       hx0 = 1.d0/dfloat(NXX_FG)                                         ! 005
                       hy0 = 1.d0/dfloat(NYY_FG)                                         ! 006
                       hz0 = 1.d0/dfloat(NZZ_FG)                                         ! 007
! 
! Right Hand Side Function
                  Res_000  = 0.D0                                                        ! 008
                  Err_000  = 0.D0                                                        ! 009
   do              i       = 1,NXX_FG+1                                                  ! 010
    do               j     = 1,NYY_FG+1                                                  ! 011
     do                k   = 1,NZZ_FG+1                                                  ! 012
              RHSF(i,j,k)  =(1.D0-a_)*dexp((1.D0-a_)*(Xv_FG(i) + Yv_FG(j) + Zv_FG(k))) & ! 013
                           +(1.D0-b_)*dexp((1.D0-b_)*(Xv_FG(i) + Yv_FG(j) + Zv_FG(k))) & ! 013
                           +(1.D0-c_)*dexp((1.D0-c_)*(Xv_FG(i) + Yv_FG(j) + Zv_FG(k)))   ! 013
      if(dabs(RHSF(i,j,k)) > Res_000) Res_000 = dabs(RHSF(i,j,k))                        ! 014    
              Error        = dexp(Xv_FG(i) + Yv_FG(j) + Zv_FG(k))                        ! 015 
           if(Error        > Err_000) Err_000 = Error                                    ! 016
     end do                                                                              ! 017
    end do                                                                               ! 018
   end do                                                                                ! 019
!
! Boundary Conditions
!
    do        j    = 1,NYY_FG+1                                                          ! 020  
     do         k  = 1,NZZ_FG+1                                                          ! 021
       BC_0YZ(j,k) = dexp(           Yv_FG(j) + Zv_FG(k)) ! Plane x = 0                  ! 022 
       BC_1YZ(j,k) = dexp(1.D0     + Yv_FG(j) + Zv_FG(k)) ! Plane x = 1	                 ! 023
     end do                                                                              ! 024
    end do                                                                               ! 025

   do         i    = 1,NXX_FG+1                                                          ! 026   
     do         k  = 1,NZZ_FG+1                                                          ! 027
       BC_X0Z(i,k) = dexp(Xv_FG(i)            + Zv_FG(k)) ! Plane y = 0                  ! 028
       BC_X1Z(i,k) = dexp(Xv_FG(i) + 1.D0     + Zv_FG(k)) ! Plane y = 1                  ! 029   
     end do                                                                              ! 030
   end do                                                                                ! 031

   do         i    = 1,NXX_FG+1                                                          ! 032 
    do          j  = 1,NYY_FG+1                                                          ! 033  
       BC_XY0(i,j) = dexp(Xv_FG(i) + Yv_FG(j)           ) ! Plane z = 0                  ! 034
       BC_XY1(i,j) = dexp(Xv_FG(i) + Yv_FG(j) + 1.D0    ) ! Plane z = 1                  ! 035     
    end do                                                                               ! 036  
   end do                                                                                ! 037 
                     write(*,1)                                                          ! 038    
                     write(1,1)                                                          ! 039     
                     write(*,2) 0,1.d0,Err_000,0                                         ! 040
                     write(1,2) 0,1.d0,Err_000,0                                         ! 041
                     call TimeS(Time0)                                                   ! 042 
1 format(/'       q     ResMAX        ErrMAX          Time')  
2 format(i8,2(2x,D12.5),2x,i8)
end subroutine  Starting_Guess_and_Boundary_Conditions


subroutine  Averaged_Coefficients
 do   i      = 1,NXX_FG+1                           ! 001 
  do    j    = 1,NYY_FG+1                           ! 002 
   do     k  = 1,NZZ_FG+1                           ! 003 
   J3(i,j,k) = Lambda_X(i,j,k)       *F3Y(j)*F3Z(k) ! 004   
   end do                                           ! 005 
  end do                                            ! 006
 end do                                             ! 007
               call LHS_yz; L_A_X = J3              ! 008                            
               
 do   i      = 1,NXX_FG+1                           ! 009  
  do    j    = 1,NYY_FG+1                           ! 010 
   do     k  = 1,NZZ_FG+1                           ! 011  
   J3(i,j,k) = Lambda_Y(i,j,k)*F3X(i)       *F3Z(k) ! 012  
   end do                                           ! 013  
  end do                                            ! 014  
 end do                                             ! 015  
               call LHS_xz; L_A_Y = J3              ! 016                              
                         
  do  i      = 1,NXX_FG+1                           ! 017                           
  do    j    = 1,NYY_FG+1                           ! 018
   do     k  = 1,NZZ_FG+1                           ! 019
   J3(i,j,k) = Lambda_Z(i,j,k)*F3X(i)*F3Y(j)        ! 020
   end do                                           ! 021
  end do                                            ! 022
 end do                                             ! 023
               call LHS_xy; L_A_Z = J3              ! 024                                          
end subroutine  Averaged_Coefficients        


subroutine  Matrix_and_Vector                                                       
real*8      La_w,La_e,La_s,La_n,La_d,La_u,h,D2uDx2,D2uDy2,D2uDz2  

                           F3X    = 0.d0                                            ! 001
        do i = 1,NXX_FG+1; F3X(i) = Xf_FG(i) - Xf_FG(i-1); end do                   ! 002

                           F3Y    = 0.d0                                            ! 003
        do j = 1,NYY_FG+1; F3Y(j) = Yf_FG(j) - Yf_FG(j-1); end do                   ! 004
 
                           F3Z    = 0.d0                                            ! 005
        do k = 1,NZZ_FG+1; F3Z(k) = Zf_FG(k) - Zf_FG(k-1); end do                   ! 006

if(q==1)                                                                       then ! 007
 if(LevelC  == LevelM)                                                         then ! 008
               call Averaged_Coefficients                                           ! 009
 else                                                                 
    Lambda_X =(hatU+corU)**(-a_)                                                    ! 010
    Lambda_Y =(hatU+corU)**(-b_)                                                    ! 011
    Lambda_Z =(hatU+corU)**(-c_)                                                    ! 012
               call Averaged_Coefficients                                           ! 013
 end if                                                                             ! 014  
         J3  = 0.d0                                                                 ! 015  
 do   i      = 2,NXX_FG                                                             ! 016
  do    j    = 2,NYY_FG                                                             ! 017  
   do     k  = 2,NZZ_FG                                                             ! 018
   J3(i,j,k) = RHSF(i,j,k)                            *F3X(i)*F3Y(j)*F3Z(k)         ! 019
   end do                                                                           ! 020   
  end do                                                                            ! 021
 end do                                                                             ! 022  

else
    Lambda_X =(hatU+corU)**(-a_)                                                    ! 023
    Lambda_Y =(hatU+corU)**(-b_)                                                    ! 024
    Lambda_Z =(hatU+corU)**(-c_)                                                    ! 025
               call Averaged_Coefficients                                           ! 026  
		 
         J3  = 0.d0                                                                 ! 027  
 do   i      = 2,NXX_FG                                                             ! 028 
  do    j    = 2,NYY_FG                                                             ! 029   
   do     k  = 2,NZZ_FG                                                             ! 030                                                              
         h   =      Lambda_X(i  ,j  ,k  )                                           ! 031
      La_w   = 2.D0*Lambda_X(i-1,j  ,k  )*h/(Lambda_X(i-1,j  ,k  )+h)               ! 032
      La_e   = 2.D0*Lambda_X(i+1,j  ,k  )*h/(Lambda_X(i+1,j  ,k  )+h)               ! 033

         h   =      Lambda_Y(i  ,j  ,k  )                                           ! 034
      La_s   = 2.D0*Lambda_Y(i  ,j-1,k  )*h/(Lambda_Y(i  ,j-1,k  )+h)               ! 035 
      La_n   = 2.D0*Lambda_Y(i  ,j+1,k  )*h/(Lambda_Y(i  ,j+1,k  )+h)               ! 036 

         h   =      Lambda_Z(i  ,j  ,k  )                                           ! 037
      La_d   = 2.D0*Lambda_Z(i  ,j  ,k-1)*h/(Lambda_Z(i  ,j  ,k-1)+h)               ! 038
      La_u   = 2.D0*Lambda_Z(i  ,j  ,k+1)*h/(Lambda_Z(i  ,j  ,k+1)+h)               ! 039

         h   =       hatU(i  ,j  ,k  )                                              ! 040
      D2uDx2 =(La_e*(hatU(i+1,j  ,k  )-h) - La_w*(h-hatU(i-1,j  ,k  )))/hx0**2      ! 041 
      D2uDy2 =(La_n*(hatU(i  ,j+1,k  )-h) - La_s*(h-hatU(i  ,j-1,k  )))/hy0**2      ! 042
      D2uDz2 =(La_u*(hatU(i  ,j  ,k+1)-h) - La_d*(h-hatU(i  ,j  ,k-1)))/hz0**2      ! 043
   J3(i,j,k) =(RHSF(i,j,k) - D2uDx2 - D2uDy2 - D2uDz2)*F3X(i)*F3Y(j)*F3Z(k)         ! 044
   end do                                                                           ! 045 
  end do                                                                            ! 046
 end do                                                                             ! 047
end if                                                                              ! 048  
               call RHSI3D                                                          ! 049  
end subroutine Matrix_and_Vector


subroutine  Convergence_Test
real*8      La_w,La_e,La_s,La_n,La_d,La_u,h 
real*8      D2uDx2,D2uDy2,D2uDz2

                      hatU = corU + hatU                                       ! 001
                      corU = 0.d0                                              ! 002
                    BC_0YZ = 0.d0; BC_1YZ = 0.d0                               ! 003 
                    BC_X0Z = 0.d0; BC_X1Z = 0.d0                               ! 004
                    BC_XY0 = 0.d0; BC_XY1 = 0.d0                               ! 005

                  Lambda_X = hatU**(-a_)                                       ! 006
                  Lambda_Y = hatU**(-b_)                                       ! 007
                  Lambda_Z = hatU**(-c_)                                       ! 008  

                   ResMAX  = 0.d0                                              ! 009
                   ErrMAX  = 0.d0                                              ! 010  
   do              i       = 2,NXX_FG                                          ! 011
    do               j     = 2,NYY_FG                                          ! 012
     do                k   = 2,NZZ_FG                                          ! 013 
         h   =      Lambda_X(i  ,j  ,k  )                                      ! 014
      La_w   = 2.D0*Lambda_X(i-1,j  ,k  )*h/(Lambda_X(i-1,j  ,k  )+h)          ! 015
      La_e   = 2.D0*Lambda_X(i+1,j  ,k  )*h/(Lambda_X(i+1,j  ,k  )+h)          ! 016

         h   =      Lambda_Y(i  ,j  ,k  )                                      ! 017
      La_s   = 2.D0*Lambda_Y(i  ,j-1,k  )*h/(Lambda_Y(i  ,j-1,k  )+h)          ! 018
      La_n   = 2.D0*Lambda_Y(i  ,j+1,k  )*h/(Lambda_Y(i  ,j+1,k  )+h)          ! 019

         h   =      Lambda_Z(i  ,j  ,k  )                                      ! 020
      La_d   = 2.D0*Lambda_Z(i  ,j  ,k-1)*h/(Lambda_Z(i  ,j  ,k-1)+h)          ! 021
      La_u   = 2.D0*Lambda_Z(i  ,j  ,k+1)*h/(Lambda_Z(i  ,j  ,k+1)+h)          ! 022 

         h   =       hatU(i  ,j  ,k  )                                         ! 023
      D2uDx2 =(La_e*(hatU(i+1,j  ,k  )-h) - La_w*(h-hatU(i-1,j  ,k  )))/hx0**2 ! 024  
      D2uDy2 =(La_n*(hatU(i  ,j+1,k  )-h) - La_s*(h-hatU(i  ,j-1,k  )))/hy0**2 ! 025    
      D2uDz2 =(La_u*(hatU(i  ,j  ,k+1)-h) - La_d*(h-hatU(i  ,j  ,k-1)))/hz0**2 ! 026 
         Res = RHSF(i,j,k) - D2uDx2 - D2uDy2 - D2uDz2                          ! 027
 if(dabs(Res)> ResMAX) ResMAX = dabs(Res)                                      ! 028 
    Error    = dabs(dexp(Xv_FG(i) + Yv_FG(j) + Zv_FG(k)) - hatU(i,j,k))        ! 029
 if(Error    > ErrMAX) ErrMAX = Error                                          ! 030
     end do                                                                    ! 031
    end do                                                                     ! 032
   end do                                                                      ! 033
call TimeS(Time1)                                                              ! 034
write(*,1) q,ResMAX/Res_000,ErrMAX,Time1-Time0
write(1,1) q,ResMAX/Res_000,ErrMAX,Time1-Time0
1 format(i8,2(2x,D12.5),2x,i8)
end subroutine Convergence_Test


subroutine  Vanka_type_Smoother
implicit    real*8 (a-h,o-z)
integer     way,X_block,Y_block,Z_block,M_A_P_
real*8      La_w,La_e,La_s,La_n,La_d,La_u
real*8      h,hx_,hx2,hy_,hy2,hz_,hz2,H_W,H_E,H_S,H_N,H_U,H_D
real*8      A_W,A_E,A_S,A_N,A_U,A_D,A_Px,A_Py,A_Pz,S_x,S_y,S_z   

          hx_ = hx0*dfloat(3**LevelX); hx2 = 1.D0/hx_**2                              ! 001
          hy_ = hy0*dfloat(3**LevelY); hy2 = 1.D0/hy_**2                              ! 002 
          hz_ = hz0*dfloat(3**LevelZ); hz2 = 1.D0/hz_**2                              ! 003
          
    NX_block_ =                   NX_block                                            ! 004  
    NY_block_ =          NY_block                                                     ! 005 
    NZ_block_ = NZ_block                                                              ! 006
    NSIL_     = NSIL                                                                  ! 007
    
            i = NX_block*NY_block*NZ_block; allocate( AG(i,i+1),xG(i))                ! 008
    
 if(NX_block  > Nx+1) NX_block = Nx+1                                                 ! 009  
 if(NY_block  > Ny+1) NY_block = Ny+1                                                 ! 010
 if(NZ_block  > Nz+1) NZ_block = Nz+1                                                 ! 011  

if(((Nx+1)/NX_block)*NX_block== Nx+1)                                          then   ! 012  
                      X_block =(Nx+1)/NX_block                                        ! 013
else                                                                                  ! 014
                      X_block =(Nx+1)/NX_block + 1                                    ! 015
end if                                                                                ! 016
if(((Ny+1)/NY_block)*NY_block== Ny+1)                                          then   ! 017
                      Y_block =(Ny+1)/NY_block                                        ! 018     
else                                                                                  ! 019
                      Y_block =(Ny+1)/NY_block + 1                                    ! 020   
end if                                                                                ! 021
if(((Nz+1)/NZ_block)*NZ_block== Nz+1)                                          then   ! 022
                      Z_block =(Nz+1)/NZ_block                                        ! 023  
else                                                                                  ! 024
                      Z_block =(Nz+1)/NZ_block + 1                                    ! 025 
end if                                                                                ! 026 

if((levelX==LeXmax).and.(levelY==LeYmax).and.(levelZ==LeZmax))       NSIL_ = 2*NSIL   ! 027
if(X_block*Y_block*Z_block==1)                                       NSIL_ = 1        ! 028

 do Iteration = 1,NSIL_                                                               ! 029

  do I_block  = 1,X_block                                                             ! 030
                           I_top =(I_block-1)*NX_block+1                              ! 031 
                           I_bot = I_block   *NX_block                                ! 032  
  if(I_block ==   X_block) I_top = Nx+2      -NX_block                                ! 033
  if(I_block ==   X_block) I_bot = Nx+1                                               ! 034  

  do J_block  = 1,Y_block                                                             ! 035
                           J_top =(J_block-1)*NY_block+1                              ! 036 
                           J_bot = J_block   *NY_block                                ! 037
  if(J_block ==   Y_block) J_top = Ny+2      -NY_block                                ! 038
  if(J_block ==   Y_block) J_bot = Ny+1                                               ! 039  

  do K_block  = 1,Z_block                                                             ! 040
                           K_top =(K_block-1)*NZ_block+1                              ! 041
                           K_bot = K_block   *NZ_block                                ! 042   
  if(K_block ==   Z_block) K_top = Nz+2      -NZ_block                                ! 043
  if(K_block ==   Z_block) K_bot = Nz+1                                               ! 044   

       AG = 0.D0; xG = 0.D0;  allocate(M_A_P(I_top:I_bot,J_top:J_bot,K_top:K_bot))    ! 045

        m = 0                                                                         ! 046
   do   i = I_top,I_bot                                                               ! 047
    do  j = J_top,J_bot                                                               ! 048
     do k = K_top,K_bot                                                               ! 049 
        m = m+1                                                                       ! 050 
                                   M_A_P(i,j,k) = m                                   ! 051
     end do                                                                           ! 052  
    end do                                                                            ! 053
   end do                                                                             ! 054
       N_ = m                                                                         ! 055

        m = 0                                                                         ! 056
   do   i = I_top,I_bot; II = TxV(i)                                                  ! 057
    do  j = J_top,J_bot; JJ = TyV(j)                                                  ! 058  
     do k = K_top,K_bot; KK = TzV(k)                                                  ! 059 
        m = m+1                                                                       ! 060  
	  
      way = 0;                           M_A_P_ = M_A_P(i,j,k)                        ! 061 

! Plane x = 0 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if((i==   1).and.(II==       1))                                         then   ! 062
                                   AG(m,M_A_P_) = 1.D0;          way = 1              ! 063
                                   AG(m,N_+1)   = BC_0YZ(JJ,KK)                       ! 064   
      end if                                                                          ! 065
! Plane x = 1 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if((i==Nx+1).and.(II==NXX_FG+1))                                         then   ! 066
                                   AG(m,M_A_P_) = 1.D0;          way = 1              ! 067
                                   AG(m,N_+1)   = BC_1YZ(JJ,KK)                       ! 068
      end if                                                                          ! 069
! Plane y = 0 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if((j==   1).and.(JJ==       1))                                         then   ! 070
                                   AG(m,M_A_P_) = 1.D0;          way = 1              ! 071
                                   AG(m,N_+1)   = BC_X0Z(II,KK)                       ! 072
      end if                                                                          ! 073
! Plane y = 1 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if((j==Ny+1).and.(JJ==NYY_FG+1))                                         then   ! 074
                                   AG(m,M_A_P_) = 1.D0;          way = 1              ! 075
                                   AG(m,N_+1)   = BC_X1Z(II,KK)                       ! 076
      end if                                                                          ! 077
! Plane z = 0 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if((k==   1).and.(KK==       1))                                         then   ! 078
                                   AG(m,M_A_P_) = 1.D0;          way = 1              ! 079
                                   AG(m,N_+1)   = BC_XY0(II,JJ)                       ! 080
      end if                                                                          ! 081
! Plane z = 1 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if((k==Nz+1).and.(KK==NZZ_FG+1))                                         then   ! 082
                                   AG(m,M_A_P_) = 1.D0;          way = 1              ! 083
                                   AG(m,N_+1)   = BC_XY1(II,JJ)                       ! 084
      end if                                                                          ! 085

      if(way == 0)                                                             then   ! 086  
       if(i==1)                                                                then   ! 087
         h   =      Xv_FG(II)  / hx_                                                  ! 088   
         call  Lambda_X_Faces(i,II,JJ,KK,La_w,La_e,h)                                 ! 089
        A_W  = 0.D0                                 ;  H_W = 0.D0                     ! 090
        A_Px =(-La_e + (h-2.D0)/ h      *La_w) * hx2                                  ! 091  
        A_E  =( La_e - (h-1.D0)/(h+1.D0)*La_w) * hx2;  H_E = A_E*corU(TxV(i+1),JJ,KK) ! 092
        S_x  = 2.D0  / (h+1.D0)/ h      *La_w  * hx2*          BC_0YZ(         JJ,KK) ! 093
       end if                                                                         ! 094   
       if((1<i).and.(i<Nx+1))                                                  then   ! 095
         call  Lambda_X_Faces(i,II,JJ,KK,La_w,La_e,h)                                 ! 096  
        A_W  =  La_w                           * hx2;  H_W = A_W*corU(TxV(i-1),JJ,KK) ! 097
        A_Px =-(La_w + La_e)                   * hx2                                  ! 098 
        A_E  =         La_e                    * hx2;  H_E = A_E*corU(TxV(i+1),JJ,KK) ! 099
        S_x  = 0.D0                                                                   ! 100 
       end if                                                                         ! 101 
       if(i==Nx+1)                                                             then   ! 102
         h   =(1.D0-Xv_FG(II)) / hx_                                                  ! 103  
         call  Lambda_X_Faces(i,II,JJ,KK,La_w,La_e,h)                                 ! 104
        A_W  =( La_w - (h-1.D0)/(h+1.D0)*La_e) * hx2;  H_W = A_W*corU(TxV(i-1),JJ,KK) ! 105
        A_Px =(-La_w + (h-2.D0)/ h      *La_e) * hx2                                  ! 106   
        A_E  = 0.D0                                 ;  H_E = 0.D0                     ! 107  
        S_x  = 2.D0  / (h+1.D0)/ h      *La_e  * hx2*          BC_1YZ(         JJ,KK) ! 108
       end if                                                                         ! 109

       if(j==1)                                                                then   ! 110
         h   =      Yv_FG(JJ)  / hy_                                                  ! 111
         call  Lambda_Y_Faces(j,II,JJ,KK,La_s,La_n,h)                                 ! 112
        A_S  = 0.D0                                 ;  H_S = 0.D0                     ! 113
        A_Py =(-La_n + (h-2.D0)/ h      *La_s) * hy2                                  ! 114
        A_N  =( La_n - (h-1.D0)/(h+1.D0)*La_s) * hy2;  H_N = A_N*corU(II,TyV(j+1),KK) ! 115 
        S_y  = 2.D0  / (h+1.D0)/ h      *La_s  * hy2*          BC_X0Z(II,         KK) ! 116
       end if                                                                         ! 117
       if((1<j).and.(j<Ny+1))                                                  then   ! 118
         call  Lambda_Y_Faces(j,II,JJ,KK,La_s,La_n,h)                                 ! 119
        A_S  =  La_s                           * hy2;  H_S = A_S*corU(II,TyV(j-1),KK) ! 120 
        A_Py =-(La_s + La_n)                   * hy2                                  ! 121 
        A_N  =         La_n                    * hy2;  H_N = A_N*corU(II,TyV(j+1),KK) ! 122 
        S_y  = 0.D0                                                                   ! 123 
       end if                                                                         ! 124 
       if(j==Ny+1)                                                             then   ! 125 
         h   =(1.D0-Yv_FG(JJ)) / hy_                                                  ! 126 
         call  Lambda_Y_Faces(j,II,JJ,KK,La_s,La_n,h)                                 ! 127 
        A_S  =( La_s - (h-1.D0)/(h+1.D0)*La_n) * hy2;  H_S = A_S*corU(II,TyV(j-1),KK) ! 128 
        A_Py =(-La_s + (h-2.D0)/ h      *La_n) * hy2                                  ! 129 
        A_N  = 0.D0                                 ;  H_N = 0.D0                     ! 130  
        S_y  = 2.D0  / (h+1.D0)/ h      *La_n  * hy2*          BC_X1Z(II,         KK) ! 131 
       end if                                                                         ! 132 
                                                                                      
       if(k==1)                                                                then   ! 133  
         h   =      Zv_FG(KK)  / hz_                                                  ! 134  
         call  Lambda_Z_Faces(k,II,JJ,KK,La_d,La_u,h)                                 ! 135  
        A_D  = 0.D0                                 ;  H_D = 0.D0                     ! 136  
        A_Pz =(-La_u + (h-2.D0)/ h      *La_d) * hz2;                                 ! 137  
        A_U  =( La_u - (h-1.D0)/(h+1.D0)*La_d) * hz2;  H_U = A_U*corU(II,JJ,TzV(k+1)) ! 138  
        S_z  = 2.D0  / (h+1.D0)/ h      *La_d  * hz2*          BC_XY0(II,JJ         ) ! 149  
       end if                                                                         ! 140  
       if((1<k).and.(k<Nz+1))                                                  then   ! 141  
         call  Lambda_Z_Faces(k,II,JJ,KK,La_d,La_u,h)                                 ! 142  
        A_D  =  La_d                           * hz2;  H_D = A_D*corU(II,JJ,TzV(k-1)) ! 143  
        A_Pz =-(La_d + La_u)                   * hz2                                  ! 144  
        A_U  =         La_u                    * hz2;  H_U = A_U*corU(II,JJ,TzV(k+1)) ! 145  
        S_z  = 0.D0                                                                   ! 146  
       end if                                                                         ! 147  
       if(k==Nz+1)                                                             then   ! 148  
         h   =(1.D0-Zv_FG(KK)) / hz_                                                  ! 149       
         call  Lambda_Z_Faces(k,II,JJ,KK,La_d,La_u,h)                                 ! 150      
        A_D  =( La_d - (h-1.D0)/(h+1.D0)*La_u) * hz2;  H_D = A_D*corU(II,JJ,TzV(k-1)) ! 151      
        A_Pz =(-La_d + (h-2.D0)/ h      *La_u) * hz2                                  ! 152       
        A_U  = 0.D0                                 ;  H_U = 0.D0                     ! 153       
        S_z  = 2.D0  / (h+1.D0)/ h      *La_u  * hz2*          BC_XY1(II,JJ         ) ! 154       
       end if                                                                         ! 155       
                                                                                        
       if(i-1 >= I_top)     AG(m,M_A_P(i-1,j  ,k  )) = A_W                            ! 156       
       if(i-1 >= I_top)                                H_W = 0.D0                     ! 157       

       if(i+1 <= I_bot)     AG(m,M_A_P(i+1,j  ,k  )) = A_E                            ! 158       
       if(i+1 <= I_bot)                                H_E = 0.D0                     ! 159       

       if(j-1 >= J_top)     AG(m,M_A_P(i  ,j-1,k  )) = A_S                            ! 160       
       if(j-1 >= J_top)                                H_S = 0.D0                     ! 161       

       if(j+1 <= J_bot)     AG(m,M_A_P(i  ,j+1,k  )) = A_N                            ! 162       
       if(j+1 <= J_bot)                                H_N = 0.D0                     ! 163       

       if(k-1 >= K_top)     AG(m,M_A_P(i  ,j  ,k-1)) = A_D                            ! 164       
       if(k-1 >= K_top)                                H_D = 0.D0                     ! 165       

       if(k+1 <= K_bot)     AG(m,M_A_P(i  ,j  ,k+1)) = A_U                            ! 166       
       if(k+1 <= K_bot)                                H_U = 0.D0                     ! 167       

                            AG(m,M_A_P(i  ,j  ,k  )) = A_Px + A_Py + A_Pz             ! 168       
                            AG(m,N_+1)               = J3(II,JJ,KK) &                 ! 169       
                       - S_x - S_y - S_z - H_W - H_E - H_S - H_N - H_D - H_U          ! 169       
      end if                                                                          ! 170       
     end do                                                                           ! 171       
    end do                                                                            ! 172       
   end do                                                                             ! 173       
                              call GAUSSIAN_ELIMINATION(N_)                           ! 174       
        m = 0                                                                         ! 175       
   do   i = I_top,I_bot; II = TxV(i)                                                  ! 176       
    do  j = J_top,J_bot; JJ = TyV(j)                                                  ! 177       
     do k = K_top,K_bot; KK = TzV(k)                                                  ! 178       
        m = m+1                                                                       ! 179       
                              corU(II,JJ,KK) = xG(m)                                  ! 180       
     end do                                                                           ! 181     
    end do                                                                            ! 182     
   end do                                                                             ! 183     
                              deallocate(M_A_P)                                       ! 184     
  end do                                                                              ! 185     
  end do                                                                              ! 186     
  end do                                                                              ! 187     

 end do                                                                               ! 188     
     NX_block = NX_block_                                                             ! 189     
     NY_block = NY_block_                                                             ! 190     
     NZ_block = NZ_block_                                                             ! 191
     deallocate(AG,xG)                                                                ! 192
end subroutine  Vanka_type_Smoother


subroutine  Lambda_X_Faces(i,II,JJ,KK,La_w,La_e,h)
real*8      La_w,La_e,h,hL_m,hL_0,hL_p,hL_b0,hL_b1

       hL_0  = L_A_X( TxV(i  ),JJ,KK)                                               ! 001 

       if(i==1)                                                                then ! 002
       hL_p  = L_A_X( TxV(i+1),JJ,KK)                                               ! 003 
       hL_b0 = L_A_X(       1 ,JJ,KK)                                               ! 004 
       hL_m  = 2.D0/ h /         (h+1.D0) * hL_b0 &                                 ! 005
             + 2.D0/ h *(h-1.D0)          * hL_0  &                                 ! 005
             -          (h-1.D0)/(h+1.D0) * hL_p                                    ! 005 
       end if                                                                       ! 006 
       if((1<i).and.(i<Nx+1))                                                  then ! 007
       hL_m  = L_A_X( TxV(i-1),JJ,KK)                                               ! 008 
       hL_p  = L_A_X( TxV(i+1),JJ,KK)                                               ! 009
       end if                                                                       ! 010  
       if(i==Nx+1)                                                             then ! 011
       hL_m  = L_A_X( TxV(i-1),JJ,KK)                                               ! 012  
       hL_b1 = L_A_X(NXX_FG+1 ,JJ,KK)                                               ! 013
       hL_p  = 2.D0/ h /         (h+1.D0) * hL_b1 &                                 ! 014  
             + 2.D0/ h *(h-1.D0)          * hL_0  &                                 ! 014     
             -          (h-1.D0)/(h+1.D0) * hL_m                                    ! 014 
       end if                                                                       ! 015
       La_w  = 2.D0 * hL_0 * hL_m / (hL_0 + hL_m)                                   ! 016
       La_e  = 2.D0 * hL_0 * hL_p / (hL_0 + hL_p)                                   ! 017
end subroutine  Lambda_X_Faces

subroutine  Lambda_Y_Faces(j,II,JJ,KK,La_s,La_n,h)
real*8      La_s,La_n,h,hL_m,hL_0,hL_p,hL_b0,hL_b1

       hL_0  = L_A_Y(II, TyV(j  ),KK) 

       if(j==1)                                                                then
       hL_p  = L_A_Y(II, TyV(j+1),KK) 
       hL_b0 = L_A_Y(II,       1 ,KK)  
       hL_m  = 2.D0/ h /         (h+1.D0) * hL_b0 &
             + 2.D0/ h *(h-1.D0)          * hL_0  &
             -          (h-1.D0)/(h+1.D0) * hL_p
       end if
       if((1<j).and.(j<Ny+1))                                                  then
       hL_m  = L_A_Y(II, TyV(j-1),KK) 
       hL_p  = L_A_Y(II, TyV(j+1),KK) 
       end if
       if(j==Ny+1)                                                             then
       hL_m  = L_A_Y(II, TyV(j-1),KK) 
       hL_b1 = L_A_Y(II,NYY_FG+1 ,KK) 
       hL_p  = 2.D0/ h /         (h+1.D0) * hL_b1 &
             + 2.D0/ h *(h-1.D0)          * hL_0  &
             -          (h-1.D0)/(h+1.D0) * hL_m
       end if
       La_s  = 2.D0 * hL_0 * hL_m / (hL_0 + hL_m)
       La_n  = 2.D0 * hL_0 * hL_p / (hL_0 + hL_p)
end subroutine  Lambda_Y_Faces

subroutine  Lambda_Z_Faces(k,II,JJ,KK,La_d,La_u,h)
real*8      La_d,La_u,h,hL_m,hL_0,hL_p,hL_b0,hL_b1

       hL_0  = L_A_Z(II,JJ, TzV(k  )) 

       if(k==1)                                                                then
       hL_p  = L_A_Z(II,JJ, TzV(k+1)) 
       hL_b0 = L_A_Z(II,JJ,       1 )  
       hL_m  = 2.D0/ h /         (h+1.D0) * hL_b0 &
             + 2.D0/ h *(h-1.D0)          * hL_0  &
             -          (h-1.D0)/(h+1.D0) * hL_p
       end if
       if((1<k).and.(k<Nz+1))                                                  then
       hL_m  = L_A_Z(II,JJ, TzV(k-1)) 
       hL_p  = L_A_Z(II,JJ, TzV(k+1)) 
       end if
       if(k==Nz+1)                                                             then
       hL_m  = L_A_Z(II,JJ, TzV(k-1)) 
       hL_b1 = L_A_Z(II,JJ,NZZ_FG+1 ) 
       hL_p  = 2.D0/ h /         (h+1.D0) * hL_b1 &
             + 2.D0/ h *(h-1.D0)          * hL_0  &
             -          (h-1.D0)/(h+1.D0) * hL_m
       end if
       La_d  = 2.D0 * hL_0 * hL_m / (hL_0 + hL_m)
       La_u  = 2.D0 * hL_0 * hL_p / (hL_0 + hL_p)
end subroutine  Lambda_Z_Faces


SUBROUTINE   GAUSSIAN_ELIMINATION(NG)
Implicit     Real*8(A-H,O-Z)
N1         = NG+1
DO      K  = 1,NG
 K1        = K+1
 S         = AG(K,K)
 J         = K
 DO     I  = K1,NG
  R        = AG(I,K)
  IF(DABS(R).Gt.DABS(S)) then
    S      = R
    J      = I
  end if
 end do
 IF(J.Ne.K)              then
  DO    I  = K,N1
   R       = AG(K,I)
   AG(K,I) = AG(J,I)
   AG(J,I) = R
  end do
 end if
 DO     J  = K1,N1
   AG(K,J) = AG(K,J)/S
 end do
 DO   I    = K1,NG
 R         = AG(I,K)
  DO    J  = K1,N1
   AG(I,J) = AG(I,J)-AG(K,J)*R
  end do
 end do
end do
xG(NG)     = AG(NG,N1)
DO    I    = NG-1,1,-1
S          = AG(I,N1)
 DO     J  = I+1,NG
  S        = S-AG(I,J)*xG(J)
 end do
xG(I)      = S
end do
END SUBROUTINE GAUSSIAN_ELIMINATION


subroutine  TimeS(The_Time)
integer     IHR,IMIN,ISEC,I100TH,The_Time
call GETTIM(IHR,IMIN,ISEC,I100TH)
The_Time  = IHR*3600 + IMIN*60 + ISEC
end subroutine TimeS

END MODULE RMT_3D_2020_Example_13_2


