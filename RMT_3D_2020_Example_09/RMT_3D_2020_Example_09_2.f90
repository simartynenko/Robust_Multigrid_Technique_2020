MODULE  RMT_3D_2020_Example_09_2
USE     RMT_3D_2020_Example_09_1
USE     RMT_3D_2020
    
PRIVATE    
PUBLIC :: Matrix_and_Vector, Vanka_type_Smoother

CONTAINS     

subroutine  Matrix_and_Vector                           
                                 F3X    = 0.d0                                     ! 001
              do i = 1,NXX_FG+1; F3X(i) = Xf_FG(i) - Xf_FG(i-1); end do            ! 002

                                 F3Y    = 0.d0                                     ! 003
              do j = 1,NYY_FG+1; F3Y(j) = Yf_FG(j) - Yf_FG(j-1); end do            ! 004
 
                                 F3Z    = 0.d0                                     ! 005
              do k = 1,NZZ_FG+1; F3Z(k) = Zf_FG(k) - Zf_FG(k-1); end do            ! 006
 
           J3 = 0.d0                                                               ! 007  
  do    i     = 1,NXX_FG+1                                                         ! 008
   do     j   = 1,NYY_FG+1                                                         ! 009
    do      k = 1,NZZ_FG+1                                                         ! 010
     J3(i,j,k)= dexp(a_*Xv_FG(i) + b_*Yv_FG(j) + c_*Zv_FG(k))*F3X(i)*F3Y(j)*F3Z(k) ! 011
    end do                                                                         ! 012
   end do                                                                          ! 013
  end do                                                                           ! 014
  call RHSI3D                                                                      ! 015
  
end subroutine Matrix_and_Vector


subroutine  Vanka_type_Smoother
real*8 ErrMAX,Error,Exact_In
real*8 x_a,x_b,y_a,y_b,z_a,z_b
       ErrMAX = 0.D0                                        ! 001   
  do        i = 1,Nx+1                                      ! 002
          x_a = Xf_FG(max(       0, TxF(i-1)))              ! 003     
          x_b = Xf_FG(min(NXX_FG+1, TxF(i  )))              ! 004     
   do       j = 1,Ny+1                                      ! 005
          y_a = Yf_FG(max(       0, TyF(j-1)))              ! 006        
          y_b = Yf_FG(min(NYY_FG+1, TyF(j  )))              ! 007        
    do      k = 1,Nz+1                                      ! 008  
          z_a = Zf_FG(max(       0, TzF(k-1)))              ! 009
          z_b = Zf_FG(min(NZZ_FG+1, TzF(k  )))              ! 010
     Exact_In =(dexp(a_*x_b)-dexp(a_*x_a))/(x_b-x_a)/a_ &   ! 011
              *(dexp(b_*y_b)-dexp(b_*y_a))/(y_b-y_a)/b_ &   ! 011
              *(dexp(c_*z_b)-dexp(c_*z_a))/(z_b-z_a)/c_     ! 011 
        Error = dabs(Exact_In - J3(TxV(i),TyV(j),TzV(k)))   ! 012
     if(Error > ErrMAX) ErrMAX = Error                      ! 013
    end do                                                  ! 014 
   end do                                                   ! 015 
  end do                                                    ! 016
write(1,1) LevelX,LevelY,LevelZ,GridsX,GridsY,GridsZ,ErrMAX ! 017            
1 format(3i2,3i4,3x,D12.5)        
end subroutine  Vanka_type_Smoother

END MODULE RMT_3D_2020_Example_09_2


