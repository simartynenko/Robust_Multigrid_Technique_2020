MODULE  RMT_3D_2020_Example_12_2
USE     RMT_3D_2020_Example_12_1
USE     RMT_3D_2020
USE     RMT_3D_TIMES

PRIVATE

PUBLIC :: Starting_Guess_and_Boundary_Conditions, Matrix_and_Vector
PUBLIC :: Convergence_Test, point_Seidel_Smoother
    
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


subroutine  point_Seidel_Smoother
implicit    real*8 (a-h,o-z)
integer     way
real*8      La_w,La_e,La_s,La_n,La_d,La_u
real*8      h,hx_,hx2,hy_,hy2,hz_,hz2,H_W,H_E,H_S,H_N,H_U,H_D
real*8      A_W,A_E,A_S,A_N,A_U,A_D,A_Px,A_Py,A_Pz,S_x,S_y,S_z,A_P,B_   

          hx_ = hx0*dfloat(3**LevelX); hx2 = 1.D0/hx_**2                              ! 001
          hy_ = hy0*dfloat(3**LevelY); hy2 = 1.D0/hy_**2                              ! 002 
          hz_ = hz0*dfloat(3**LevelZ); hz2 = 1.D0/hz_**2                              ! 003

                                                                 NSIL_ =   NSIL       ! 004
if((levelX==LeXmax).and.(levelY==LeYmax).and.(levelZ==LeZmax))   NSIL_ = 2*NSIL       ! 005

 do Iteration = 1,NSIL_                                                               ! 006

   do       i = 1,Nx+1; II = TxV(i)                                                   ! 007
    do      j = 1,Ny+1; JJ = TyV(j)                                                   ! 008  
     do     k = 1,Nz+1; KK = TzV(k)                                                   ! 009 

          A_W = 0.d0; A_Px = 0.d0; A_E = 0.d0; S_x = 0.d0
          A_S = 0.d0; A_Py = 0.d0; A_N = 0.d0; S_y = 0.d0;  A_P = 1.D0
          A_U = 0.d0; A_Pz = 0.d0; A_D = 0.d0; S_z = 0.d0;  way = 0    

! Plane x = 0 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if((i==   1).and.(II==       1))                                         then   ! 010
                                   B_  = BC_0YZ(JJ,KK);     way = 1                   ! 011   
      end if                                                                          ! 012
! Plane x = 1 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if((i==Nx+1).and.(II==NXX_FG+1))                                         then   ! 013
                                   B_  = BC_1YZ(JJ,KK);     way = 1                   ! 014
      end if                                                                          ! 015
! Plane y = 0 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if((j==   1).and.(JJ==       1))                                         then   ! 016
                                   B_  = BC_X0Z(II,KK);     way = 1                   ! 017
      end if                                                                          ! 018
! Plane y = 1 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if((j==Ny+1).and.(JJ==NYY_FG+1))                                         then   ! 019
                                   B_  = BC_X1Z(II,KK);     way = 1                   ! 020
      end if                                                                          ! 021
! Plane z = 0 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if((k==   1).and.(KK==       1))                                         then   ! 022
                                   B_  = BC_XY0(II,JJ);     way = 1                   ! 023
      end if                                                                          ! 024
! Plane z = 1 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if((k==Nz+1).and.(KK==NZZ_FG+1))                                         then   ! 025
                                   B_  = BC_XY1(II,JJ);     way = 1                   ! 026
      end if                                                                          ! 027

      if(way == 0)                                                             then   ! 028  
       if(i==1)                                                                then   ! 029
         h   =      Xv_FG(II)  / hx_                                                  ! 030   
         call  Lambda_X_Faces(i,II,JJ,KK,La_w,La_e,h)                                 ! 031
        A_W  = 0.D0                                 ;  H_W = 0.D0                     ! 032
        A_Px =(-La_e + (h-2.D0)/ h      *La_w) * hx2                                  ! 033  
        A_E  =( La_e - (h-1.D0)/(h+1.D0)*La_w) * hx2;  H_E = A_E*corU(TxV(i+1),JJ,KK) ! 034
        S_x  = 2.D0  / (h+1.D0)/ h      *La_w  * hx2*          BC_0YZ(         JJ,KK) ! 035
       end if                                                                         ! 036   
       if((1<i).and.(i<Nx+1))                                                  then   ! 037
         call  Lambda_X_Faces(i,II,JJ,KK,La_w,La_e,h)                                 ! 038  
        A_W  =  La_w                           * hx2;  H_W = A_W*corU(TxV(i-1),JJ,KK) ! 039
        A_Px =-(La_w + La_e)                   * hx2                                  ! 040 
        A_E  =         La_e                    * hx2;  H_E = A_E*corU(TxV(i+1),JJ,KK) ! 041
        S_x  = 0.D0                                                                   ! 042 
       end if                                                                         ! 043 
       if(i==Nx+1)                                                             then   ! 044
         h   =(1.D0-Xv_FG(II)) / hx_                                                  ! 045  
         call  Lambda_X_Faces(i,II,JJ,KK,La_w,La_e,h)                                 ! 046
        A_W  =( La_w - (h-1.D0)/(h+1.D0)*La_e) * hx2;  H_W = A_W*corU(TxV(i-1),JJ,KK) ! 047
        A_Px =(-La_w + (h-2.D0)/ h      *La_e) * hx2                                  ! 048   
        A_E  = 0.D0                                 ;  H_E = 0.D0                     ! 049  
        S_x  = 2.D0  / (h+1.D0)/ h      *La_e  * hx2*          BC_1YZ(         JJ,KK) ! 050
       end if                                                                         ! 051

       if(j==1)                                                                then   ! 052
         h   =      Yv_FG(JJ)  / hy_                                                  ! 053
         call  Lambda_Y_Faces(j,II,JJ,KK,La_s,La_n,h)                                 ! 054
        A_S  = 0.D0                                 ;  H_S = 0.D0                     ! 055
        A_Py =(-La_n + (h-2.D0)/ h      *La_s) * hy2                                  ! 056
        A_N  =( La_n - (h-1.D0)/(h+1.D0)*La_s) * hy2;  H_N = A_N*corU(II,TyV(j+1),KK) ! 057
        S_y  = 2.D0  / (h+1.D0)/ h      *La_s  * hy2*          BC_X0Z(II,         KK) ! 058
       end if                                                                         ! 059
       if((1<j).and.(j<Ny+1))                                                  then   ! 060
         call  Lambda_Y_Faces(j,II,JJ,KK,La_s,La_n,h)                                 ! 061
        A_S  =  La_s                           * hy2;  H_S = A_S*corU(II,TyV(j-1),KK) ! 062
        A_Py =-(La_s + La_n)                   * hy2                                  ! 063
        A_N  =         La_n                    * hy2;  H_N = A_N*corU(II,TyV(j+1),KK) ! 064
        S_y  = 0.D0                                                                   ! 065
       end if                                                                         ! 066
       if(j==Ny+1)                                                             then   ! 067
         h   =(1.D0-Yv_FG(JJ)) / hy_                                                  ! 068
         call  Lambda_Y_Faces(j,II,JJ,KK,La_s,La_n,h)                                 ! 069
        A_S  =( La_s - (h-1.D0)/(h+1.D0)*La_n) * hy2;  H_S = A_S*corU(II,TyV(j-1),KK) ! 070
        A_Py =(-La_s + (h-2.D0)/ h      *La_n) * hy2                                  ! 071
        A_N  = 0.D0                                 ;  H_N = 0.D0                     ! 072 
        S_y  = 2.D0  / (h+1.D0)/ h      *La_n  * hy2*          BC_X1Z(II,         KK) ! 073
       end if                                                                         ! 074
                                                                                        
       if(k==1)                                                                then   ! 075
         h   =      Zv_FG(KK)  / hz_                                                  ! 076
         call  Lambda_Z_Faces(k,II,JJ,KK,La_d,La_u,h)                                 ! 077
        A_D  = 0.D0                                 ;  H_D = 0.D0                     ! 078
        A_Pz =(-La_u + (h-2.D0)/ h      *La_d) * hz2;                                 ! 079
        A_U  =( La_u - (h-1.D0)/(h+1.D0)*La_d) * hz2;  H_U = A_U*corU(II,JJ,TzV(k+1)) ! 080
        S_z  = 2.D0  / (h+1.D0)/ h      *La_d  * hz2*          BC_XY0(II,JJ         ) ! 081
       end if                                                                         ! 082
       if((1<k).and.(k<Nz+1))                                                  then   ! 083
         call  Lambda_Z_Faces(k,II,JJ,KK,La_d,La_u,h)                                 ! 084
        A_D  =  La_d                           * hz2;  H_D = A_D*corU(II,JJ,TzV(k-1)) ! 085
        A_Pz =-(La_d + La_u)                   * hz2                                  ! 086
        A_U  =         La_u                    * hz2;  H_U = A_U*corU(II,JJ,TzV(k+1)) ! 087
        S_z  = 0.D0                                                                   ! 088
       end if                                                                         ! 089
       if(k==Nz+1)                                                             then   ! 090
         h   =(1.D0-Zv_FG(KK)) / hz_                                                  ! 091
         call  Lambda_Z_Faces(k,II,JJ,KK,La_d,La_u,h)                                 ! 092
        A_D  =( La_d - (h-1.D0)/(h+1.D0)*La_u) * hz2;  H_D = A_D*corU(II,JJ,TzV(k-1)) ! 093
        A_Pz =(-La_d + (h-2.D0)/ h      *La_u) * hz2                                  ! 094
        A_U  = 0.D0                                 ;  H_U = 0.D0                     ! 095
        S_z  = 2.D0  / (h+1.D0)/ h      *La_u  * hz2*          BC_XY1(II,JJ         ) ! 096
       end if                                                                         ! 097

         A_P = A_Px + A_Py + A_Pz                                                     ! 098      
         B_  = J3(II,JJ,KK)- S_x - S_y - S_z - H_W - H_E - H_S - H_N - H_D - H_U      ! 099       
      end if                                                                          ! 100       
           
                                            corU(II,JJ,KK) = B_/A_P                   ! 101
     end do                                                                           ! 102       
    end do                                                                            ! 103       
   end do                                                                             ! 104       
  end do                                                                              ! 105 
end subroutine  point_Seidel_Smoother


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

END MODULE RMT_3D_2020_Example_12_2


