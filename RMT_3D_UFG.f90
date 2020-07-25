MODULE   RMT_3D_UFG
USE      RMT_3D_2020                                         ! 000        
REAL*8 , PUBLIC :: h_x,h_y,h_z
INTEGER, PUBLIC :: KindGX,KindGY,KindGZ
    
CONTAINS

subroutine  U_Finest_Grid(KindGX,KindGY,KindGZ)
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - X-Direction - - - - - - - - - 
if(KindGX==1)                                                                  then
               h_x = 1.d0         /dfloat(  NXX_FG  )        ! 001
       do       i  = 1,                     NXX_FG+1         ! 002
          Xv_FG(i) = dfloat(  i-1)/dfloat(  NXX_FG  )        ! 003  (2.1a)
       end do                                                ! 004
       do       i  = 1,                     NXX_FG           ! 005
          Xf_FG(i) = dfloat(2*i-1)/dfloat(2*NXX_FG  )        ! 006  (2.1b)
       end do                                                ! 007 
end if
if(KindGX==2)                                                                  then
               h_x = 2.d0         /dfloat(2*NXX_FG+1)        ! 008
       do       i  = 1,                     NXX_FG+1         ! 009
          Xv_FG(i) = dfloat(2*i-1)/dfloat(2*NXX_FG+1)        ! 010  (2.2a)
       end do                                                ! 011
       do       i  = 0,                     NXX_FG           ! 012
          Xf_FG(i) = dfloat(2*i  )/dfloat(2*NXX_FG+1)        ! 013  (2.2b)
       end do                                                ! 014
end if
if(KindGX==3)                                                                  then
               h_x = 2.d0         /dfloat(2*NXX_FG+1)        ! 015
       do       i  = 1,                     NXX_FG+1         ! 016
          Xv_FG(i) = dfloat(2*i-2)/dfloat(2*NXX_FG+1)        ! 017  (2.3a)
          Xf_FG(i) = dfloat(2*i-1)/dfloat(2*NXX_FG+1)        ! 018  (2.3b)
       end do    
end if    
if(KindGX==4)                                                                  then
               h_x = 1.d0         /dfloat(  NXX_FG+1)        ! 019
       do       i  = 1,                     NXX_FG+1         ! 020
          Xv_FG(i) = dfloat(2*i-1)/dfloat(2*NXX_FG+2)        ! 021  (2.4a)
       end do                                                ! 022
       do       i  = 0,                     NXX_FG+1         ! 023
          Xf_FG(i) = dfloat(  i  )/dfloat(  NXX_FG+1)        ! 024  (2.4b)
       end do                                                ! 025
end if

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - Y-Direction - - - - - - - - - 
if(KindGY==1)                                                                  then
               h_y = 1.d0         /dfloat(  NYY_FG  )
       do       j  = 1,                     NYY_FG+1
          Yv_FG(j) = dfloat(  j-1)/dfloat(  NYY_FG  )
       end do
       do       j  = 1,                     NYY_FG
          Yf_FG(j) = dfloat(2*j-1)/dfloat(2*NYY_FG  )
       end do
end if
if(KindGY==2)                                                                  then
               h_y = 2.d0         /dfloat(2*NYY_FG+1)
       do       j  = 1,                     NYY_FG+1
          Yv_FG(j) = dfloat(2*j-1)/dfloat(2*NYY_FG+1)
       end do
       do       j  = 0,                     NYY_FG
          Yf_FG(j) = dfloat(2*j  )/dfloat(2*NXX_FG+1)
       end do
end if
if(KindGY==3)                                                                  then
               h_y = 2.d0         /dfloat(2*NYY_FG+1)
       do       j  = 1,                     NYY_FG+1
          Yv_FG(j) = dfloat(2*j-2)/dfloat(2*NYY_FG+1)
          Yf_FG(j) = dfloat(2*j-1)/dfloat(2*NYY_FG+1)
       end do    
end if    
if(KindGY==4)                                                                  then
               h_y = 1.d0         /dfloat(  NYY_FG+1)
       do       j  = 1,                     NYY_FG+1
          Yv_FG(j) = dfloat(2*j-1)/dfloat(2*NYY_FG+2)
       end do
       do       j  = 0,                     NYY_FG+1
          Yf_FG(j) = dfloat(  j  )/dfloat(  NYY_FG+1)
       end do 
end if

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - Z-Direction - - - - - - - - - 
if(KindGZ==1)                                                                  then
               h_z = 1.d0         /dfloat(  NZZ_FG  )
       do       k  = 1,                     NZZ_FG+1
          Zv_FG(k) = dfloat(  k-1)/dfloat(  NZZ_FG  )
       end do
       do       k  = 1,                     NZZ_FG
          Zf_FG(k) = dfloat(2*k-1)/dfloat(2*NZZ_FG  )
       end do
end if
if(KindGZ==2)                                                                  then
               h_z = 2.d0         /dfloat(2*NZZ_FG+1)
       do       k  = 1,                     NZZ_FG+1
          Zv_FG(k) = dfloat(2*k-1)/dfloat(2*NZZ_FG+1)
       end do
       do       k  = 0,                     NZZ_FG
          Zf_FG(k) = dfloat(2*k  )/dfloat(2*NZZ_FG+1)
       end do
end if
if(KindGZ==3)                                                                  then
               h_z = 2.d0         /dfloat(2*NZZ_FG+1)
       do       k  = 1,                     NZZ_FG+1
          Zv_FG(k) = dfloat(2*k-2)/dfloat(2*NZZ_FG+1)
          Zf_FG(k) = dfloat(2*k-1)/dfloat(2*NZZ_FG+1)
       end do    
end if    
if(KindGZ==4)                                                                  then
               h_z = 1.d0         /dfloat(  NZZ_FG+1)
       do       k  = 1,                     NZZ_FG+1
          Zv_FG(k) = dfloat(2*k-1)/dfloat(2*NZZ_FG+2)
       end do
       do       k  = 0,                     NZZ_FG+1
          Zf_FG(k) = dfloat(  k  )/dfloat(  NZZ_FG+1)
       end do 
end if
end subroutine  U_Finest_Grid

END MODULE RMT_3D_UFG