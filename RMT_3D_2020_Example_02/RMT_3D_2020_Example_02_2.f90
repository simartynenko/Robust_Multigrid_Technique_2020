MODULE RMT_3D_2020_Example_02_2
USE    RMT_3D_2020_Example_02_1

PRIVATE

PUBLIC :: MSSize, Starting_Guess_and_Boundary_Conditions, GAUSSIAN_ELIMINATION
PUBLIC :: TimeS,Convergence_Test

CONTAINS
    
SUBROUTINE   MSSize    
if(((Nx+1)/NX_block)*NX_block== Nx+1)                                          then  ! 001
                      X_block =(Nx+1)/NX_block                                       ! 002
else                                                   
                      X_block =(Nx+1)/NX_block + 1                                   ! 003  
end if               
if(((Ny+1)/NY_block)*NY_block== Ny+1)                                          then
                      Y_block =(Ny+1)/NY_block 
else
                      Y_block =(Ny+1)/NY_block + 1 
end if               
if(((Nz+1)/NZ_block)*NZ_block== Nz+1)                                          then
                      Z_block =(Nz+1)/NZ_block 
else
                      Z_block =(Nz+1)/NZ_block + 1 
end if      
                           N_ = NX_block*NY_block*NZ_block                           ! 004                                 
 allocate( AG(N_,N_+1), xG(N_))                                                      ! 005            
 
 allocate(      u(0:NXX_FG+2,0:NYY_FG+2,0:NZZ_FG+2))                                 ! 006
 allocate(      f(0:NXX_FG+2,0:NYY_FG+2,0:NZZ_FG+2))                                 ! 007  
 allocate( BC_0YZ(           0:NYY_FG+2,0:NZZ_FG+2))                                 ! 008 
 allocate( BC_1YZ(           0:NYY_FG+2,0:NZZ_FG+2))                                 ! 009
 allocate( BC_X0Z(0:NXX_FG+2,           0:NZZ_FG+2))                                 ! 010
 allocate( BC_X1Z(0:NXX_FG+2,           0:NZZ_FG+2))                                 ! 011
 allocate( BC_XY0(0:NXX_FG+2,0:NYY_FG+2           ))                                 ! 012
 allocate( BC_XY1(0:NXX_FG+2,0:NYY_FG+2           ))                                 ! 013
END SUBROUTINE   MSSize  
    

Subroutine Starting_Guess_and_Boundary_Conditions
real       x,y,z
                 u = 0.d0;       f = 0.d0                        ! 001             
             Fnorm =(a_**2+b_**2+c_**2)*dexp(a_+b_+c_)           ! 002
      do     i     = 2,NXX_FG;    x = dfloat(i-1)/dfloat(NXX_FG) ! 003
       do      j   = 2,NYY_FG;    y = dfloat(j-1)/dfloat(NYY_FG) ! 004
        do       k = 2,NZZ_FG;    z = dfloat(k-1)/dfloat(NZZ_FG) ! 005
           F(i,j,k)=-(a_**2+b_**2+c_**2)*dexp(a_*x+b_*y+c_*z)    ! 006
        end do                                                   ! 007
       end do                                                    ! 008
      end do                                                     ! 009
!= = = = = = = = = = = = =  X planes  = = = = = = = = = = = = =
       do      j   = 2,NYY_FG;    y = dfloat(j-1)/dfloat(NYY_FG) ! 010 
        do       k = 2,NZZ_FG;    z = dfloat(k-1)/dfloat(NZZ_FG) ! 011 
      BC_0YZ(  j,k)= dexp(     b_*y+c_*z)                        ! 012 
      BC_1YZ(  j,k)= dexp(a_  +b_*y+c_*z)                        ! 013  
        end do                                                   ! 014 
       end do                                                    ! 015 
!= = = = = = = = = = = = =  Y planes  = = = = = = = = = = = = =   
      do     i     = 2,NXX_FG;    x = dfloat(i-1)/dfloat(NXX_FG) ! 016 
        do       k = 2,NZZ_FG;    z = dfloat(k-1)/dfloat(NZZ_FG) ! 017    
      BC_X0Z(i,  k)= dexp(a_*x     +c_*z)                        ! 018    
      BC_X1Z(i,  k)= dexp(a_*x+b_  +c_*z)                        ! 019               
        end do                                                   ! 020    
      end do                                                     ! 021    
!= = = = = = = = = = = = =  Z planes  = = = = = = = = = = = = =      
      do     i     = 2,NXX_FG;    x = dfloat(i-1)/dfloat(NXX_FG) ! 022    
       do      j   = 2,NYY_FG;    y = dfloat(j-1)/dfloat(NYY_FG) ! 023    
      BC_XY0(i,j  )= dexp(a_*x+b_*y     )                        ! 024
      BC_XY1(i,j  )= dexp(a_*x+b_*y+c_  )                        ! 025           
       end do                                                    ! 026
      end do                                                     ! 027
                     call TimeS(Time0)                           ! 028      
                     write(*,1); write(1,1)                      ! 029
                     call Convergence_Test(way)                  ! 030
1 format('    Iter   ResMAX/Fnorm      ErrMAX         Time')                  
end Subroutine Starting_Guess_and_Boundary_Conditions


Subroutine Convergence_Test(way)
integer    way,TimeC
real*8     Res,ResMAX,Err,ErrMAX,D2uDx2,D2uDy2,D2uDz2,x,y,z
            ResMAX = 0.d0;   ErrMAX = 0.d0;  way = 0                         ! 001 
  do         i     = 2,NXX_FG;    x = dfloat(i-1)/dfloat(NXX_FG)             ! 002 
   do          j   = 2,NYY_FG;    y = dfloat(j-1)/dfloat(NYY_FG)             ! 003 
    do           k = 2,NZZ_FG;    z = dfloat(k-1)/dfloat(NZZ_FG)             ! 004 
            D2uDx2 =(U(i-1,j  ,k  ) - 2.d0*U(i,j,k) + U(i+1,j  ,k  ))/h_x**2 ! 005 
            D2uDy2 =(U(i  ,j-1,k  ) - 2.d0*U(i,j,k) + U(i  ,j+1,k  ))/h_y**2 ! 006 
            D2uDz2 =(U(i  ,j  ,k-1) - 2.d0*U(i,j,k) + U(i  ,j  ,k+1))/h_z**2 ! 007            
               Res = dabs(F(i,j,k) + D2uDx2 + D2uDy2 + D2uDz2)               ! 008 
            if(Res > ResMAX) ResMAX = Res                                    ! 009
               Err = dabs(u(i,j,k)-dexp(a_*x+b_*y+c_*z))                     ! 010
            if(Err > ErrMAX) ErrMAX = Err                                    ! 012      
    end do                                                                   ! 013
   end do                                                                    ! 014
  end do                                                                     ! 015   
                     call TimeS(TimeC)                                       ! 016        
  if(ResMAX/Fnorm<ResLIM) way = 1                                            ! 017       
                     write(*,1) Iteration,ResMAX/Fnorm,ErrMAX,TimeC-Time0    ! 018       
                     write(1,1) Iteration,ResMAX/Fnorm,ErrMAX,TimeC-Time0    ! 019       
1 format(i8,2x,2(D13.6,2x),2x,i6)                                          
end subroutine Convergence_Test           


SUBROUTINE   GAUSSIAN_ELIMINATION(NG)
Implicit     Real*8(A-H,O-Z)
N1         = NG+1
DO      K  = 1,NG

 K1        = K+1
 S         = AG(K,K)
 J         = K
 DO     I  = K1,NG
  R        = AG(I,K)
  IF(DABS(R).Gt.DABS(S))                                                       then
    S      = R
    J      = I
  end if
 end do
 IF(J.Ne.K)                                                                    then
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
DO   I     = NG-1,1,-1
S          = AG(I,N1)
 DO     J  = I+1,NG
  S        = S-AG(I,J)*xG(J)
 end do
xG(I)      = S
end do
END SUBROUTINE GAUSSIAN_ELIMINATION

    
SUBROUTINE  TimeS(The_Time)
integer     IHR,IMIN,ISEC,I100TH,The_Time
!call GETTIM(IHR,IMIN,ISEC,I100TH)
The_Time  = IHR*3600 + IMIN*60 + ISEC
END SUBROUTINE TimeS
                             

END MODULE RMT_3D_2020_Example_02_2
