MODULE  RMT_3D_2020_Example_15_2
USE     RMT_3D_2020_Example_15_1
USE     RMT_3D_2020
    
PRIVATE    
PUBLIC :: Starting_Guess_and_Boundary_Conditions
PUBLIC :: Matrix_and_Vector_in_Mini_Cycle, Mini_Cycle

CONTAINS     

subroutine Starting_Guess_and_Boundary_Conditions

allocate(MyF(0:NXX_FG+2,0:NYY_FG+2,0:NZZ_FG+2)); MyF = 0.d0    ! 001

do     i      = 1,NXX_FG+1                                     ! 002
 do      j    = 1,NYY_FG+1                                     ! 003
  do       k  = 1,NZZ_FG+1                                     ! 004 
   MyF(i,j,k) = dexp(a_*Xv_FG(i) + b_*Yv_FG(j) + c_*Zv_FG(k))  ! 005   
  end do                                                       ! 006
 end do                                                        ! 007
end do                                                         ! 008 
end subroutine Starting_Guess_and_Boundary_Conditions


subroutine     Matrix_and_Vector_in_Mini_Cycle
real*8         h_x_,h_y_,h_z_
!
!- - - - -    X-Integral        - - - - -
!
              Ih   =(3**LevelXd-1)/2                                                  ! 001   
   do         i    =-Ih,Ih+NxxD+2                                                     ! 002 
          F3X(i)   = 0.d0                                                             ! 003 
    do        i_   = TxV_(i)-Ih,TxV_(i)+Ih                                            ! 004 
     if((1.le.i_).and.(i_.le.NXX_FG+1))     F3X(i) = F3X(i) + Xf_FG(i_) - Xf_FG(i_-1) ! 005 
    end do                                                                            ! 006 
   end do                                                                             ! 007 
!                                                                                   
!- - - - -    Y-Integral        - - - - -
!
              Jh   =(3**LevelYd-1)/2                                                  ! 008        
   do         j    =-Jh,Jh+NyyD+2                                                     ! 009      
          F3Y(j)   = 0.d0                                                             ! 010      
    do        j_   = TyV_(j)-Jh,TyV_(j)+Jh                                            ! 011      
     if((1.le.j_).and.(j_.le.NYY_FG+1))     F3Y(j) = F3Y(j) + Yf_FG(j_) - Yf_FG(j_-1) ! 012      
    end do                                                                            ! 013      
   end do                                                                             ! 014      
!                                                                                           
!- - - - -    Z-Integral        - - - - -
!
              Kh   =(3**LevelZd-1)/2                                                  ! 015       
   do         k    =-Kh,Kh+NzzD+2                                                     ! 016      
          F3Z(k)   = 0.d0                                                             ! 017      
    do        k_   = TzV_(k)-Kh,TzV_(k)+Kh                                            ! 018      
     if((1.le.k_).and.(k_.le.NZZ_FG+1))     F3Z(k) = F3Z(k) + Zf_FG(k_) - Zf_FG(k_-1) ! 019      
    end do                                                                            ! 020      
   end do                                                                             ! 021      

   do         i    =-Ih,Ih+NxxD+2                                                     ! 022     
    do        j    =-Jh,Jh+NyyD+2                                                     ! 023   
     do       k    =-Kh,Kh+NzzD+2                                                     ! 024   
                                                                                      ! 025   
         J3(i,j,k) = 0.d0                                                             ! 026   
      do      i_   = TxV_(i)-Ih,TxV_(i)+Ih                                            ! 027   
      if((1.le.i_).and.(i_.le.NXX_FG+1))                                       then   ! 028   
                                              h_x_ = Xf_FG(i_) - Xf_FG(i_-1)          ! 029   
                                                                                      ! 030   
       do     j_   = TyV_(j)-Jh,TyV_(j)+Jh                                            ! 031   
       if((1.le.j_).and.(j_.le.NYY_FG+1))                                      then   ! 032   
                                              h_y_ = Yf_FG(j_) - Yf_FG(j_-1)          ! 033   
                                                                                      ! 034   
        do    k_   = TzV_(k)-Kh,TzV_(k)+Kh                                            ! 035   
        if((1.le.k_).and.(k_.le.NZZ_FG+1))                                     then   ! 036   
                                              h_z_ = Zf_FG(k_) - Zf_FG(k_-1)          ! 037   
                                                                                      ! 038   
         J3(i,j,k) = J3(i,j,k) + h_x_*h_y_*h_z_*MyF(i_,j_,k_)                         ! 039   
                                                                                      ! 040   
        end if                                                                        ! 041  
        end do                                                                        ! 042   
       end if                                                                         ! 043   
       end do                                                                         ! 044   
      end if                                                                          ! 045   
      end do                                                                          ! 046   
                                                                                      ! 047   
     end do                                                                           ! 048   
    end do                                                                            ! 049   
   end do                                                                             ! 050   
                     call RHSI3D                                                      ! 051     
end subroutine Matrix_and_Vector_in_Mini_Cycle  


subroutine    Mini_Cycle
real*8        x_a,x_b,y_a,y_b,z_a,z_b,Exact_In,Error,ErrMAXG    
do    LevelCC = max(LeXmax,LeYmax,LeZmax        ),0,-1              ! 001
      LevelX  = min(LeXmax              ,LevelCC)                   ! 002   
      LevelY  = min(       LeYmax       ,LevelCC)                   ! 003   
      LevelZ  = min(              LeZmax,LevelCC)                   ! 004   

      call Matrix_and_Vector_in_Mini_Cycle                          ! 005   

 do   GridsX  = 1,3**LevelX; call Index_Mapping('-',1)              ! 006   
  do  GridsY  = 1,3**LevelY; call Index_Mapping('-',2)              ! 007   
   do GridsZ  = 1,3**LevelZ; call Index_Mapping('-',3)              ! 008   

      ErrMAXG = 0.D0                                                ! 009   
    do      i = 1,Nx+1                                              ! 010   
          x_a = Xf_FG(max(       0, TxF_(TxF(i-1))))                ! 011   
          x_b = Xf_FG(min(NXX_FG+1, TxF_(TxF(i  ))))                ! 012   
     do     j = 1,Ny+1                                              ! 013   
          y_a = Yf_FG(max(       0, TyF_(TyF(j-1))))                ! 014   
          y_b = Yf_FG(min(NYY_FG+1, TyF_(TyF(j  ))))                ! 015   
      do    k = 1,Nz+1                                              ! 016   
          z_a = Zf_FG(max(       0, TzF_(TzF(k-1))))                ! 017   
          z_b = Zf_FG(min(NZZ_FG+1, TzF_(TzF(k  ))))                ! 018   
     Exact_In =(dexp(a_*x_b)-dexp(a_*x_a))/(x_b-x_a)/a_ &           ! 019   
              *(dexp(b_*y_b)-dexp(b_*y_a))/(y_b-y_a)/b_ &           ! 019
              *(dexp(c_*z_b)-dexp(c_*z_a))/(z_b-z_a)/c_             ! 019
        Error = dabs(Exact_In - J3(TxV(i),TyV(j),TzV(k)))           ! 020
     if(Error > ErrMAXG) ErrMAXG = Error                            ! 021
      end do                                                        ! 022
     end do                                                         ! 023
    end do                                                          ! 024
write(1,1) LevelXd,LevelYd,LevelZd,GridsXd,GridsYd,GridsZd, &       ! 025
           LevelX ,LevelY ,LevelZ ,GridsX ,GridsY ,GridsZ ,ErrMAXG  ! 025
   end do                                                           ! 026
  end do                                                            ! 027
 end do                                                             ! 028
end do                                                              ! 029
1 format(6(2x,i1),6(2x,i3),3x,D13.6)                                   
end subroutine  Mini_Cycle                                            
                                                                       

END MODULE RMT_3D_2020_Example_15_2






















