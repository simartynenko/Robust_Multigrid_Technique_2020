!
! Point Gauss-Seidel iterations for solving Poisson equation (1.6)
!
! Chapter I, Section 2
!
! File: RMT_3D_2020_Example_01.f90
!    
program RMT_3D_2020_Example_01
integer NXX_FG,NYY_FG,NZZ_FG,pGS_IterMAX,Iteration,way
real*8  a_,b_,c_,h_x,h_y,h_z,psi,hf,hx,hy,hz,ResLIM,Fnorm

real*8, allocatable :: u(:,:,:), f(:,:,:)                                 ! 001

NXX_FG = 50;  h_x = 1.d0/dfloat(NXX_FG);  a_ = 1.d0                       ! 002
NYY_FG = 50;  h_y = 1.d0/dfloat(NYY_FG);  b_ = 1.d0                       ! 003
NZZ_FG = 50;  h_z = 1.d0/dfloat(NZZ_FG);  c_ = 1.d0                       ! 004
      pGS_IterMAX = 1000000;         ItCheck = 1000                       ! 005
           ResLIM = 1.d-7;         Iteration = 0                          ! 006

allocate(u(0:NXX_FG+2,0:NYY_FG+2,0:NZZ_FG+2), &                           ! 007 
         f(0:NXX_FG+2,0:NYY_FG+2,0:NZZ_FG+2))                             ! 007

open(1,file='RMT_3D_2020_Example_01.res')                                 ! 008
call Starting_Guess_and_Boundary_Conditions                               ! 009
               psi = 2.d0*(1.d0/h_x**2+1.d0/h_y**2+1.d0/h_z**2)           ! 010
               hf  = 1.d0                                     /psi        ! 011
               hx  =       1.d0/h_x**2                        /psi        ! 012
               hy  =                   1.d0/h_y**2            /psi        ! 013
               hz  =                               1.d0/h_z**2/psi        ! 014                  
 do      Iteration = 1,pGS_IterMAX                                        ! 015 
  do         i     = 2,NXX_FG                                             ! 016
   do          j   = 2,NYY_FG                                             ! 017
    do           k = 2,NZZ_FG                                             ! 018
           u(i,j,k)= hf*f(i,j,k) + hx*(u(i-1,j  ,k  ) + u(i+1,j  ,k  )) & ! 019
                                 + hy*(u(i  ,j-1,k  ) + u(i  ,j+1,k  )) & ! 019
                                 + hz*(u(i  ,j  ,k-1) + u(i  ,j  ,k+1))   ! 019
    end do                                                                ! 020
   end do                                                                 ! 021
  end do                                                                  ! 022
     if((Iteration / ItCheck)*ItCheck==Iteration) &                       ! 023     
     call Convergence_Test(way); if(way==1) exit                          ! 023
 end do                                                                   ! 024 
close(1)                                                                  ! 025
deallocate(u,f)                                                           ! 026  
        
CONTAINS
    
Subroutine Starting_Guess_and_Boundary_Conditions
real       x,y,z
                 u = 0.d0;       f = 0.d0                        ! 001             
             Fnorm =(a_**2+b_**2+c_**2)*dexp(a_+b_+c_)           ! 002
      do     i     = 1,NXX_FG+1;  x = dfloat(i-1)/dfloat(NXX_FG) ! 003
       do      j   = 1,NYY_FG+1;  y = dfloat(j-1)/dfloat(NYY_FG) ! 004
        do       k = 1,NZZ_FG+1;  z = dfloat(k-1)/dfloat(NZZ_FG) ! 005
           f(i,j,k)=-(a_**2+b_**2+c_**2)*dexp(a_*x+b_*y+c_*z)    ! 006
        end do                                                   ! 007
       end do                                                    ! 008
      end do                                                     ! 009
!= = = = = = = = = = = = =  X planes  = = = = = = = = = = = = =
       do      j   = 1,NYY_FG+1;  y = dfloat(j-1)/dfloat(NYY_FG) ! 010 
        do       k = 1,NZZ_FG+1;  z = dfloat(k-1)/dfloat(NZZ_FG) ! 011 
                     u(       1,j,k)= dexp(     b_*y+c_*z)       ! 012 
                     u(NXX_FG+1,j,k)= dexp(a_  +b_*y+c_*z)       ! 013  
        end do                                                   ! 014 
       end do                                                    ! 015 
!= = = = = = = = = = = = =  Y planes  = = = = = = = = = = = = =   
      do     i     = 1,NXX_FG+1;  x = dfloat(i-1)/dfloat(NXX_FG) ! 016 
        do       k = 1,NZZ_FG+1;  z = dfloat(k-1)/dfloat(NZZ_FG) ! 017    
                     u(i,       1,k)= dexp(a_*x     +c_*z)       ! 018    
                     u(i,NYY_FG+1,k)= dexp(a_*x+b_  +c_*z)       ! 019               
        end do                                                   ! 020    
      end do                                                     ! 021    
!= = = = = = = = = = = = =  Z planes  = = = = = = = = = = = = =      
      do     i     = 1,NXX_FG+1;  x = dfloat(i-1)/dfloat(NXX_FG) ! 022    
       do      j   = 1,NYY_FG+1;  y = dfloat(j-1)/dfloat(NYY_FG) ! 023    
                     u(i,j,       1)= dexp(a_*x+b_*y     )       ! 024
                     u(i,j,NZZ_FG+1)= dexp(a_*x+b_*y+c_  )       ! 025           
       end do                                                    ! 026
      end do                                                     ! 027
                     write(*,1); write(1,1)                      ! 028
                     call Convergence_Test(way)                  ! 029
1 format('    Iter   ResMAX/Fnorm      ErrMAX')                     
end Subroutine Starting_Guess_and_Boundary_Conditions

Subroutine Convergence_Test(way)
integer    way
real*8     Res,ResMAX,Err,ErrMAX,D2uDx2,D2uDy2,D2uDz2,x,y,z
            ResMAX = 0.d0;   ErrMAX = 0.d0;  way = 0                         ! 001 
  do         i     = 2,NXX_FG;    x = dfloat(i-1)/dfloat(NXX_FG)             ! 002 
   do          j   = 2,NYY_FG;    y = dfloat(j-1)/dfloat(NYY_FG)             ! 003 
    do           k = 2,NZZ_FG;    z = dfloat(k-1)/dfloat(NZZ_FG)             ! 004 
            D2uDx2 =(u(i-1,j  ,k  ) - 2.d0*u(i,j,k) + u(i+1,j  ,k  ))/h_x**2 ! 005 
            D2uDy2 =(u(i  ,j-1,k  ) - 2.d0*u(i,j,k) + u(i  ,j+1,k  ))/h_y**2 ! 006 
            D2uDz2 =(u(i  ,j  ,k-1) - 2.d0*u(i,j,k) + u(i  ,j  ,k+1))/h_z**2 ! 007            
               Res = dabs(f(i,j,k) + D2uDx2 + D2uDy2 + D2uDz2)               ! 008 
            if(Res > ResMAX) ResMAX = Res                                    ! 009
               Err = dabs(u(i,j,k)-dexp(a_*x+b_*y+c_*z))                     ! 010
            if(Err > ErrMAX) ErrMAX = Err                                    ! 011      
    end do                                                                   ! 012
   end do                                                                    ! 013
  end do                                                                     ! 014       
  if(ResMAX/Fnorm<ResLIM) way = 1                                            ! 015       
                     write(*,1) Iteration,ResMAX/Fnorm,ErrMAX                ! 016       
                     write(1,1) Iteration,ResMAX/Fnorm,ErrMAX                ! 017       
1 format(i8,2x,2(D13.6,2x))                                          
end subroutine Convergence_Test                                            

end program RMT_3D_2020_Example_01
