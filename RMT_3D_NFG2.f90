MODULE   RMT_3D_NFG2
USE      RMT_3D_2020                             
         PRIVATE                                               ! 001 
         PUBLIC :: Finest_Grid                                 ! 002
    
CONTAINS

subroutine   Finest_Grid
real*8       STEPA,STEPB,mesh,W

             STEPA = 0.004d0                                   ! 003                                           
             STEPB = 0.090d0                                   ! 004             
             call CELLS(0.0d0,0.5d0,STEPA,STEPB,   0,N,W,mesh) ! 005
                N1 = N-1                                       ! 006
        Xf_FG      = 0.d0                                      ! 007
       do     i    = 0,N1                                      ! 008
        Xf_FG(i+1) = Xf_FG(i) + mesh*W**i                      ! 009
       end do                                                  ! 010
  
             call CELLS(0.5d0,1.0d0,STEPB,STEPA,N1+1,N,W,mesh) ! 011    
       do     i    = N1+1,N-1                                  ! 012    
        Xf_FG(i+1) = Xf_FG(i) + mesh*W**(i-N1-1)               ! 013    
       end do                                                  ! 014    
            NXX_FG = N-1                                       ! 015    
                                                                   
       do     i    = 1,NXX_FG+1                                ! 016    
        Xv_FG(i  ) =(Xf_FG(i)+Xf_FG(i-1))/2.d0                 ! 017    
       end do                                                  ! 018
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - Y-Direction - - - - - - - - - 
         NYY_FG    = NXX_FG                                    ! 019
          Yv_FG    = Xv_FG                                     ! 020
          Yf_FG    = Xf_FG                                     ! 021
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - Z-Direction - - - - - - - - - 
         NZZ_FG    = NXX_FG                                    ! 022          
          Zv_FG    = Xv_FG                                     ! 023
          Zf_FG    = Xf_FG                                     ! 024
end subroutine  Finest_Grid

    
Subroutine     CELLS(A0,B0,STEPA,STEPB,N1,N2,W,STEP)
implicit       REAL*8(A-H,O-Z)
AL           = B0-A0
W            =(AL-STEPA)/(AL-STEPB)
N2           = 1+N1+Int(Dlog(STEPB/STEPA)/Dlog(W))
do     loop  = 1,1000
 STNEW       = AL*(W-1.d0)/(W**(N2-N1)-1.d0)
 DSTEP       = Dabs(1.d0-STEPA/STNEW)*100.d0
 if(DSTEP.LE.1.D-12) exit
 STEPA       = STNEW
end do
if(loop>999)   write(*,*) '   Error in the grid generation!!!'
STEP         = STEPA
end subroutine CELLS

END MODULE RMT_3D_NFG2
