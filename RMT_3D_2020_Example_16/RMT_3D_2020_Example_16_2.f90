MODULE  RMT_3D_2020_Example_16_2
USE     RMT_3D_2020_Example_16_1
USE     RMT_3D_2020
USE     RMT_3D_UFG 
    
PRIVATE    
PUBLIC :: Starting_Guess_and_Boundary_Conditions,Mini_Cycle,Convergence_Test

CONTAINS     

subroutine Starting_Guess_and_Boundary_Conditions
real*8     h

allocate(    cu(0:NXX_FG+2,0:NYY_FG+2,0:NZZ_FG+2),  J3_(-10:NXX_FG+12,-10:NYY_FG+12,-10:NZZ_FG+12) ) 
allocate(  hatU(0:NXX_FG+2,0:NYY_FG+2,0:NZZ_FG+2),  RHSF(0:NXX_FG+2,0:NYY_FG+2,0:NZZ_FG+2) )
allocate(  BCX0(           0:NYY_FG+2,0:NZZ_FG+2),  BCX1(           0:NYY_FG+2,0:NZZ_FG+2) ) 
allocate(  BCY0(0:NXX_FG+2,           0:NZZ_FG+2),  BCY1(0:NXX_FG+2,           0:NZZ_FG+2) )
allocate(  BCZ0(0:NXX_FG+2,0:NYY_FG+2           ),  BCZ1(0:NXX_FG+2,0:NYY_FG+2           ) ) 
allocate(  BUX0(           0:NYY_FG+2,0:NZZ_FG+2),  BUX1(           0:NYY_FG+2,0:NZZ_FG+2) ) 
allocate(  BUY0(0:NXX_FG+2,           0:NZZ_FG+2),  BUY1(0:NXX_FG+2,           0:NZZ_FG+2) )
allocate(  BUZ0(0:NXX_FG+2,0:NYY_FG+2           ),  BUZ1(0:NXX_FG+2,0:NYY_FG+2           ) ) 
                                                  
          corU = 0.d0;   hatU = 0.d0;    RHSF = 0.d0;      J3 = 0.d0 ! 001
          BCX0 = 0.d0;   BCX1 = 0.d0;    BUX0 = 0.d0;    BUX1 = 0.d0 ! 002     
          BCY0 = 0.d0;   BCY1 = 0.d0;    BUY0 = 0.d0;    BUY1 = 0.d0 ! 003     
          BCZ0 = 0.d0;   BCZ1 = 0.d0;    BUZ0 = 0.d0;    BUZ1 = 0.d0 ! 004     
! 
! Right Hand Side Function
!            
                h  = a_**2+b_**2+c_**2                               ! 005 
          Res_000  = dexp(a_+b_+c_)*h                                ! 006
          Err_000  = dexp(a_+b_+c_)                                  ! 007
   do       i      = 1,NXX_FG+1                                      ! 008
    do        j    = 1,NYY_FG+1                                      ! 009
     do         k  = 1,NZZ_FG+1                                      ! 010
       RHSF(i,j,k) = dexp(a_*Xv_FG(i) + b_*Yv_FG(j) + c_*Zv_FG(k))*h ! 011
     end do                                                          ! 012
    end do                                                           ! 013
   end do                                                            ! 014
!
! Boundary Conditions
!
    do        j    = 1,NYY_FG+1                                      ! 015  
     do         k  = 1,NZZ_FG+1                                      ! 016
         BCX0(j,k) = dexp(              b_*Yv_FG(j) + c_*Zv_FG(k))   ! 017 
         BCX1(j,k) = dexp(a_          + b_*Yv_FG(j) + c_*Zv_FG(k))   ! 018
     end do                                                          ! 019
    end do                                                           ! 020

   do         i    = 1,NXX_FG+1                                      ! 021   
     do         k  = 1,NZZ_FG+1                                      ! 022
         BCY0(i,k) = dexp(a_*Xv_FG(i)               + c_*Zv_FG(k))   ! 023
         BCY1(i,k) = dexp(a_*Xv_FG(i) + b_          + c_*Zv_FG(k))   ! 024   
     end do                                                          ! 025
   end do                                                            ! 026

   do         i    = 1,NXX_FG+1                                      ! 027 
    do          j  = 1,NYY_FG+1                                      ! 028  
         BCZ0(i,j) = dexp(a_*Xv_FG(i) + b_*Yv_FG(j)              )   ! 029
         BCZ1(i,j) = dexp(a_*Xv_FG(i) + b_*Yv_FG(j) + c_         )   ! 030     
    end do                                                           ! 031  
   end do                                                            ! 032 
                     write(*,1)                                      ! 033    
                     write(1,1)                                      ! 034     
                     write(*,2) 0,1.d0,Err_000,0                     ! 035
                     write(1,2) 0,1.d0,Err_000,0                     ! 036
                     call TimeS(Time0)                               ! 037 
1 format(/'       q     ResMAX        ErrMAX          Time')  
2 format(i8,2(2x,D12.5),2x,i8)
end subroutine Starting_Guess_and_Boundary_Conditions


subroutine    Mini_Cycle
                     call TimeS(MC_Time0)                                 ! 000
    do         MGI = 1,MGI_DC                                             ! 001
     do    LevelCC = max(LeXmax,LeYmax,LeZmax        ),0,-1               ! 002
           LevelX  = min(LeXmax              ,LevelCC)                    ! 003   
           LevelY  = min(       LeYmax       ,LevelCC)                    ! 004   
           LevelZ  = min(              LeZmax,LevelCC)                    ! 005   
                     call Matrix_and_Vector_in_Mini_Cycle                 ! 006   
      do   GridsX  = 1,3**LevelX; call Index_Mapping('-',1)               ! 007   
       do  GridsY  = 1,3**LevelY; call Index_Mapping('-',2)               ! 008   
        do GridsZ  = 1,3**LevelZ; call Index_Mapping('-',3)               ! 009   
                     call Point_Seidel_Smoother_in_Mini_Cycle(LevelCC)    ! 010
        end do                                                            ! 011
       end do                                                             ! 012
      end do                                                              ! 013
                     call Convergence_Test_in_MiniCycle(LevelCC,MC_Time0) ! 014
     end do                                                               ! 015
                     call Update_in_MiniCycle                             ! 016
    end do                                                                ! 017 
                     call BC_Restore                                      ! 018
end subroutine  Mini_Cycle  


subroutine    Convergence_Test_in_MiniCycle(LevelCC,MC_Time0)

            ErrMAX = 0.                                                      ! 001   
 do       i        = 1,NxxD+1;   II = TxV_(i)                                ! 002 
  do         j     = 1,NyyD+1;   JJ = TyV_(j)                                ! 003 
   do           k  = 1,NzzD+1;   KK = TzV_(k)                                ! 004 
             Error = dabs(dexp(a_*Xv_FG(II) + b_*Yv_FG(JJ) + c_*Zv_FG(KK)) & ! 005 
                   - hatU(II,JJ,KK) - cu(II,JJ,KK))                          ! 006             
          if(Error > ErrMAX) ErrMAX = Error                                  ! 007 
   end do                                                                    ! 008 
  end do                                                                     ! 009 
 end do                                                                      ! 010 
                     call TimeS(MC_Time1)                                    ! 011 
 if((MGI==    1).and.(LevelCC==max(LeXmax,LeYmax,LeZmax))) write(1,*) 
 if( MGI==MGI_DC.and.(LevelCC==                       0))  &
 write(*,1)  LevelXd,LevelYd,LevelZd,GridsXd,GridsYd,GridsZd,            ErrMAX,MC_Time1-MC_Time0 
 write(1,2)  LevelXd,LevelYd,LevelZd,GridsXd,GridsYd,GridsZd,LevelCC,MGI,ErrMAX,MC_Time1-MC_Time0

1 format('   Dynamic Level (',i1,',',i1,',',i1,')    Dynamic Grid(',i1,',',i1,',',i1, & 
')    ErrMAX =',E10.3,'   Time = ',i2,' s')
2 format('   Dynamic Level (',i1,',',i1,',',i1,')    Dynamic Grid(',i1,',',i1,',',i1, &
')    LevelCC = ',i1,'    MGI = ',i1,'   ErrMAX =',E10.3,'   Time = ',i2,' s')
end subroutine  Convergence_Test_in_MiniCycle  


subroutine    Update_in_MiniCycle

 do       i        = 1,NxxD+1;           II        = TxV_(i)          ! 001   
  do         j     = 1,NyyD+1;              JJ     = TyV_(j)          ! 002   
   do           k  = 1,NzzD+1;                 KK  = TzV_(k)          ! 003   
    hatU(II,JJ,KK) = hatU(II,JJ,KK) + cu(II,JJ,KK)                    ! 004   
                                      cu(II,JJ,KK) = 0.d0             ! 005   
   end do                                                             ! 006   
  end do                                                              ! 007   
 end do                                                               ! 008   
                                                                   
    do       j     = 1,NyyD+1;                       JJ     = TyV_(j) ! 009   
     do         k  = 1,NzzD+1;                          KK  = TzV_(k) ! 010   
       BUX0(JJ,KK) = BUX0(JJ,KK) + BCX0(JJ,KK); BCX0(JJ,KK) = 0.d0    ! 011   
       BUX1(JJ,KK) = BUX1(JJ,KK) + BCX1(JJ,KK); BCX1(JJ,KK) = 0.d0    ! 012   
     end do                                                           ! 013   
    end do                                                            ! 014   
                                                                      ! 015   
   do        i     = 1,NxxD+1;                       II     = TxV_(i) ! 016   
     do         k  = 1,NzzD+1;                          KK  = TzV_(k) ! 017   
       BUY0(II,KK) = BUY0(II,KK) + BCY0(II,KK); BCY0(II,KK) = 0.d0    ! 018    
       BUY1(II,KK) = BUY1(II,KK) + BCY1(II,KK); BCY1(II,KK) = 0.d0    ! 019    
     end do                                                           ! 020 
   end do                                                             ! 021

   do        i     = 1,NxxD+1;                       II     = TxV_(i) ! 022 
    do          j  = 1,NyyD+1;                          JJ  = TyV_(j) ! 023 
       BUZ0(II,JJ) = BUZ0(II,JJ) + BCZ0(II,JJ); BCZ0(II,JJ) = 0.d0    ! 024 
       BUZ1(II,JJ) = BUZ1(II,JJ) + BCZ1(II,JJ); BCZ1(II,JJ) = 0.d0    ! 025 
    end do                                                            ! 026 
   end do                                                             ! 027      

end subroutine    Update_in_MiniCycle                                  
                                                                       
                                                                       
subroutine    BC_Restore
    do       j     = 1,NyyD+1;                       JJ     = TyV_(j) ! 001    
     do         k  = 1,NzzD+1;                          KK  = TzV_(k) ! 002    
       BCX0(JJ,KK) = BUX0(JJ,KK);               BUX0(JJ,KK) = 0.d0    ! 003    
       BCX1(JJ,KK) = BUX1(JJ,KK);               BUX1(JJ,KK) = 0.d0    ! 004    
     end do                                                           ! 005    
    end do                                                            ! 006    

   do        i     = 1,NxxD+1;                       II     = TxV_(i) ! 007    
     do         k  = 1,NzzD+1;                          KK  = TzV_(k) ! 008    
       BCY0(II,KK) = BUY0(II,KK);               BUY0(II,KK) = 0.d0    ! 009    
       BCY1(II,KK) = BUY1(II,KK);               BUY1(II,KK) = 0.d0    ! 010    
     end do                                                           ! 011 
   end do                                                             ! 012 

   do        i     = 1,NxxD+1;                       II     = TxV_(i) ! 013 
    do          j  = 1,NyyD+1;                          JJ  = TyV_(j) ! 014 
       BCZ0(II,JJ) = BUZ0(II,JJ);               BUZ0(II,JJ) = 0.d0    ! 015 
       BCZ1(II,JJ) = BUZ1(II,JJ);               BUZ1(II,JJ) = 0.d0    ! 016 
    end do                                                            ! 017 
   end do                                                             ! 018 
end subroutine    BC_Restore                                           


subroutine     Matrix_and_Vector_in_Mini_Cycle
integer        IV_m,IV_0,IV_p,JV_m,JV_0,JV_p,KV_m,KV_0,KV_p
real*8         AX_m,AX_0,AX_p,Sx,AY_m,AY_0,AY_p,Sy,AZ_m,AZ_0,AZ_p,Sz
real*8         h_x_,h_y_,h_z_,h_Lx,h_Ly,h_Lz,d2udx2,d2udy2,d2udz2,xi,h___
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

   do       i      =-Ih,Ih+NxxD+2                                                     ! 022     
    do        j    =-Jh,Jh+NyyD+2                                                     ! 023   
     do         k  =-Kh,Kh+NzzD+2                                                     ! 024   
                                                                                      ! 025   
        J3_(i,j,k) = 0.d0                                                             ! 026   
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
         J3_(i,j,k)= J3_(i,j,k) + h_x_*h_y_*h_z_*RHSF(i_,j_,k_)                       ! 039  

        end if                                                                        ! 040  
        end do                                                                        ! 041   
       end if                                                                         ! 042   
       end do                                                                         ! 043   
      end if                                                                          ! 044   
      end do                                                                          ! 045   
                                                                                      ! 046   
     end do                                                                           ! 047   
    end do                                                                            ! 048   
   end do                                                                             ! 049   

              h_Lx = 1.d0/(h_x*dfloat(3**LevelXd))**2                                     ! 050        
              h_Ly = 1.d0/(h_y*dfloat(3**LevelYd))**2                                     ! 051 
              h_Lz = 1.d0/(h_z*dfloat(3**LevelZd))**2                                     ! 052 
              h___ = h_x*dfloat(3**LevelXd)*h_y*dfloat(3**LevelYd)*h_z*dfloat(3**LevelZd) ! 053 

                J3 = 0.d0                                                                 ! 054 

      do    i      = 1,NxxD+1; IV_m = TxV_(i-1); IV_0 = TxV_(i); IV_p = TxV_(i+1)         ! 055 
       do     j    = 1,NyyD+1; JV_m = TyV_(j-1); JV_0 = TyV_(j); JV_p = TyV_(j+1)         ! 056 
        do      k  = 1,NzzD+1; KV_m = TzV_(k-1); KV_0 = TzV_(k); KV_p = TzV_(k+1)         ! 057       

        if((IV_0==       1).or.(JV_0==       1).or.(KV_0==       1).or. &                 ! 058 
           (IV_0==NXX_FG+1).or.(JV_0==NYY_FG+1).or.(KV_0==NZZ_FG+1)) then                 ! 059 

        else                                                                              ! 060 
!- - - - - - - - - - - - - - - - - - - - - - - - - - -  X-Direction                   
         if( i==1)                                                   then                 ! 061 
                xi =      Xv_FG(IV_0) *dsqrt(h_Lx)                                        ! 062 
              AX_m = 0.d0                                                                 ! 063 
              AX_o = 2.d0          /xi                                                    ! 064 
              AX_p = 2.d0/(xi+1.d0)                                                       ! 065 
              Sx   = 2.d0/(xi+1.d0)/xi                  * h_Lx &                          ! 066 
                   *             BUX0(     JV_0,KV_0)                                     ! 067 
            d2udx2 =(      -AX_o*hatU(IV_0,JV_0,KV_0)          &                          ! 068     
                           +AX_p*hatU(IV_p,JV_0,KV_0) ) * h_Lx + Sx                       ! 069      
         end if                                                                           ! 070
         if((1<i).and.(i<NxxD+1))                                    then                 ! 071 
            d2udx2 = 1.d0*(      hatU(IV_m,JV_0,KV_0)          &                          ! 072 
                           -2.d0*hatU(IV_0,JV_0,KV_0)          &                          ! 073 
                           +     hatU(IV_p,JV_0,KV_0) ) * h_Lx                            ! 074 
         end if                                                                           ! 075 
         if(           i==NxxD+1)                                    then                 ! 076 
                xi =(1.d0-Xv_FG(IV_0))*dsqrt(h_Lx)                                        ! 077 
              AX_p = 0.d0                                                                 ! 078 
              AX_o = 2.d0          /xi                                                    ! 079 
              AX_m = 2.d0/(xi+1.d0)                                                       ! 080  
              Sx   = 2.d0/(xi+1.d0)/xi                  * h_Lx &                          ! 081  
                   *             BUX1(     JV_0,KV_0)                                     ! 082  
            d2udx2 =(       AX_m*hatU(IV_m,JV_0,KV_0)          &                          ! 083  
                           -AX_o*hatU(IV_0,JV_0,KV_0) ) * h_Lx + Sx                       ! 084  
         end if                                                                           ! 085  
!- - - - - - - - - - - - - - - - - - - - - - - - - - -  Y-Direction                       
                                                                                          ! 086  
         if( j==1)                                                   then                 ! 087  
                xi =      Yv_FG(JV_0) *dsqrt(h_Ly)                                        ! 088  
              AY_m = 0.d0                                                                 ! 089 
              AY_o = 2.d0          /xi                                                    ! 090 
              AY_p = 2.d0/(xi+1.d0)                                                       ! 091 
              Sy   = 2.d0/(xi+1.d0)/xi                  * h_Ly &                          ! 092 
                   *             BUY0(IV_0,     KV_0)                                     ! 093 
            d2udy2 =(      -AY_o*hatU(IV_0,JV_0,KV_0)          &                          ! 094     
                           +AY_p*hatU(IV_0,JV_p,KV_0) ) * h_Ly + Sy                       ! 095 
         end if                                                                           ! 096 
         if((1<j).and.(j<NyyD+1))                                    then                 ! 097 
            d2udy2 =(            hatU(IV_0,JV_m,KV_0)          &                          ! 098 
                           -2.d0*hatU(IV_0,JV_0,KV_0)          &                          ! 099 
                           +     hatU(IV_0,JV_p,KV_0) ) * h_Ly                            ! 100 
         end if                                                                           ! 101 
         if(           j==NyyD+1)                                    then                 ! 102 
                xi =(1.d0-Yv_FG(JV_0))*dsqrt(h_Ly)                                        ! 103 
              AY_p = 0.d0                                                                 ! 104 
              AY_o = 2.d0          /xi                                                    ! 105 
              AY_m = 2.d0/(xi+1.d0)                                                       ! 106 
              Sy   = 2.d0/(xi+1.d0)/xi                  * h_Ly &                          ! 107 
                   *             BUY1(IV_0,     KV_0)                                     ! 108 
            d2udy2 =(       AY_m*hatU(IV_0,JV_m,KV_0)          &                          ! 109 
                           -AY_o*hatU(IV_0,JV_0,KV_0) ) * h_Ly + Sy                       ! 110 
         end if                                                                           ! 111 
!- - - - - - - - - - - - - - - - - - - - - - - - - - -  Z-Direction
         if( k==1)                                                   then                 ! 112 
                xi =      Zv_FG(KV_0) *dsqrt(h_Lz)                                        ! 113 
              AZ_m = 0.d0                                                                 ! 114 
              AZ_o = 2.d0          /xi                                                    ! 115 
              AZ_p = 2.d0/(xi+1.d0)                                                       ! 116 
              Sz   = 2.d0/(xi+1.d0)/xi                  * h_Lz &                          ! 117 
                   *             BUZ0(IV_0,JV_0     )                                     ! 118 
            d2udz2 =(      -AZ_o*hatU(IV_0,JV_0,KV_0)          &                          ! 119     
                           +AZ_p*hatU(IV_0,JV_0,KV_p) ) * h_Lz + Sz                       ! 120 
         end if                                                                           ! 121 
         if((1<k).and.(k<NzzD+1))                                    then                 ! 122 
            d2udz2 =(            hatU(IV_0,JV_0,KV_m)          &                          ! 123
                           -2.d0*hatU(IV_0,JV_0,KV_0)          &                          ! 124
                           +     hatU(IV_0,JV_0,KV_p) ) * h_Lz                            ! 125
         end if                                                                           ! 126
         if(           k==NzzD+1)                                    then                 ! 127
                xi =(1.d0-Zv_FG(KV_0))*dsqrt(h_Lz)                                        ! 128
              AZ_p = 0.d0                                                                 ! 129
              AZ_o = 2.d0          /xi                                                    ! 130
              AZ_m = 2.d0/(xi+1.d0)                                                       ! 131
              Sz   = 2.d0/(xi+1.d0)/ xi                 * h_Lz &                          ! 132
                   *             BUZ1(IV_0,JV_0     )                                     ! 133
            d2udz2 =(       AZ_m*hatU(IV_0,JV_0,KV_m)          &                          ! 134   
                           -AZ_o*hatU(IV_0,JV_0,KV_0) ) * h_Lz + Sz                       ! 135   
         end if                                                                           ! 136   
        J3(i,j,k)  = J3_(i,j,k) - (d2udx2 + d2udy2 + d2udz2)*h___                         ! 137   
        end if                                                                            ! 138   
        end do                                                                            ! 139   
       end do                                                                             ! 140   
      end do                                                                              ! 141
                            call RHSI3D                                                   ! 142
 end subroutine Matrix_and_Vector_in_Mini_Cycle  


subroutine Point_Seidel_Smoother_in_Mini_Cycle(LevelCC)
real*8     h_Lx,h_Ly,h_Lz,xi,Sx,Sy,Sz
       h_Lx  = 1.d0/(h_x*dfloat(3**(LevelX + LevelXd)))**2            ! 001  
       h_Ly  = 1.d0/(h_y*dfloat(3**(LevelY + LevelYd)))**2            ! 002 
       h_Lz  = 1.d0/(h_z*dfloat(3**(LevelZ + LevelZd)))**2            ! 003 

do        i  = 1,Nx+1; IV_0 = TxV_(TxV(i))                            ! 004 
 do       j  = 1,Ny+1; JV_0 = TyV_(TyV(j))                            ! 005 
  do      k  = 1,Nz+1; KV_0 = TzV_(TzV(k))                            ! 006 

  if((IV_0==1).or.(IV_0==NXX_FG+1).or. &                              ! 007 
     (JV_0==1).or.(JV_0==NYY_FG+1).or. &                              ! 008 
     (KV_0==1).or.(KV_0==NZZ_FG+1))                             then  ! 009 
                                                                      
  else                                                                ! 010 
                                                                      
!- - - - - - - - - - - - - - - - - - - - - - - - - - -   X-Direction   
          Sx = 0.d0                                                   ! 011 
   if(i==   1)                                                  then  ! 012 
          xi =      Xv_FG(IV_0) *dsqrt(h_Lx)                          ! 013 
          Sx = 2.d0*h_Lx/(xi+1.d0)/xi * BCX0(     JV_0,KV_0)          ! 014
   end if                                                             ! 015
   if(i==Nx+1)                                                  then  ! 016
          xi =(1.d0-Xv_FG(IV_0))*dsqrt(h_Lx)                          ! 017
          Sx = 2.d0*h_Lx/(xi+1.d0)/xi * BCX1(     JV_0,KV_0)          ! 018
   end if                                                             ! 019 
                                                                      
!- - - - - - - - - - - - - - - - - - - - - - - - - - -   Y-Direction  
          Sy = 0.d0                                                   ! 020   
   if(j==   1)                                                  then  ! 021   
          xi =      Yv_FG(JV_0) *dsqrt(h_Ly)                          ! 022   
          Sy = 2.d0*h_Ly/(xi+1.d0)/xi * BCY0(IV_0,     KV_0)          ! 023   	
   end if                                                             ! 024   
   if(j==Ny+1)                                                  then  ! 025   
          xi =(1.d0-Yv_FG(JV_0))*dsqrt(h_Ly)                          ! 026   
          Sy = 2.d0*h_Ly/(xi+1.d0)/xi * BCY1(IV_0,     KV_0)          ! 027   
   end if                                                             ! 028   
                                                                      
!- - - - - - - - - - - - - - - - - - - - - - - - - - -   Z-Direction  
          Sz = 0.d0                                                   ! 029   
   if(k==   1)                                                  then  ! 030   
          xi =      Zv_FG(KV_0) *dsqrt(h_Lz)                          ! 031   
          Sz = 2.d0*h_Lz/(xi+1.d0)/xi * BCZ0(IV_0,JV_0     )          ! 032   	
   end if                                                             ! 033   
   if(k==Nz+1)                                                  then  ! 034   
          xi =(1.d0-Zv_FG(KV_0))*dsqrt(h_Lz)                          ! 035   
          Sz = 2.d0*h_Lz/(xi+1.d0)/xi * BCZ1(IV_0,JV_0     )          ! 036   	
   end if                                                             ! 037    
                                                                     
   J3(TxV(i),TyV(j),TzV(k)) = J3(TxV(i),TyV(j),TzV(k)) - Sx - Sy - Sz ! 038   
 
  end if                                                              ! 039   
  end do                                                              ! 040   
 end do                                                               ! 041   
end do                                                                ! 042   

     MMM_ = The_number_of_smoothing_iterations                        ! 043   
if(LevelCC== max(LeXmax,LeYmax,LeZmax)) MMM_ = MMM_*2                 ! 044   
do    MMM = 1,MMM_                                                    ! 045  
 do    i  = 1,Nx+1                                                    ! 046  
  do   j  = 1,Ny+1                                                    ! 047  
   do  k  = 1,Nz+1;          call Point_Gauss_Seidel_Iteration(i,j,k) ! 048  
   end do                                                             ! 049       
  end do                                                              ! 050   
 end do                                                               ! 051        
end do                                                                ! 052  
end subroutine Point_Seidel_Smoother_in_Mini_Cycle                     
                                                                       
subroutine Point_Gauss_Seidel_Iteration(i,j,k)
real*8     h1,h2,h3,h_Lx,h_Ly,h_Lz
real*8      AX_m,  AX_o,  AX_p,  AY_m,  AY_o,  AY_p,  AZ_m,  AZ_o,  AZ_p
real*8     HAX_m,        HAX_p, HAY_m,        HAY_p, HAZ_m,        HAZ_p

       h_Lx  = 1.d0/(h_x*dfloat(3**(LevelX + LevelXd)))**2 
       h_Ly  = 1.d0/(h_y*dfloat(3**(LevelY + LevelYd)))**2	
       h_Lz  = 1.d0/(h_z*dfloat(3**(LevelZ + LevelZd)))**2 
                                                                                 
IV_m=max(TxV_(TxV(i-1)),1); IV_0=TxV_(TxV(i)); IV_p=min(TxV_(TxV(i+1)),NXX_FG+1) 
JV_m=max(TyV_(TyV(j-1)),1); JV_0=TyV_(TyV(j)); JV_p=min(TyV_(TyV(j+1)),NYY_FG+1) 
KV_m=max(TzV_(TzV(k-1)),1); KV_0=TzV_(TzV(k)); KV_p=min(TzV_(TzV(k+1)),NZZ_FG+1) 

  if((IV_0==1).or.(IV_0==NXX_FG+1).or. &
     (JV_0==1).or.(JV_0==NYY_FG+1).or. &
     (KV_0==1).or.(KV_0==NZZ_FG+1))                             then

   if(IV_0==       1) cu(IV_0,JV_0,KV_0) = BCX0(     JV_0,KV_0)
   if(IV_0==NXX_FG+1) cu(IV_0,JV_0,KV_0) = BCX1(     JV_0,KV_0)
   if(JV_0==       1) cu(IV_0,JV_0,KV_0) = BCY0(IV_0,     KV_0) 
   if(JV_0==NYY_FG+1) cu(IV_0,JV_0,KV_0) = BCY1(IV_0,     KV_0) 
   if(KV_0==       1) cu(IV_0,JV_0,KV_0) = BCZ0(IV_0,JV_0     )
   if(KV_0==NZZ_FG+1) cu(IV_0,JV_0,KV_0) = BCZ1(IV_0,JV_0     )   
   else
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - -   X-Direction
!
   HAX_m = 0.d0; HAX_p = 0.d0       
   if(i==1)                                                     then
      xi =          Xv_FG(IV_0)        *dsqrt(h_Lx)
    AX_m =  0.d0
    AX_o =-(1.d0 - (xi-1.d0)/ xi     ) * 2.d0*h_Lx 
    AX_p = (1.d0 - (xi-1.d0)/(xi+1.d0))      *h_Lx 
   HAX_p =  AX_p*cu(IV_p,JV_0,KV_0)
   end if 
   if((1<i).and.(i<Nx+1))                                       then
    AX_m =                                    h_Lx 
    AX_p =                                    h_Lx 
    AX_o =-AX_m-AX_p
   HAX_m = AX_m*cu(IV_m,JV_0,KV_0)
   HAX_p = AX_p*cu(IV_p,JV_0,KV_0)
   end if
   if(i==Nx+1)                                                  then
      xi = (1.d0 -  Xv_FG(IV_0))       *dsqrt(h_Lx)
    AX_p =  0.d0
    AX_o =-(1.d0 - (xi-1.d0)/ xi     ) * 2.d0*h_Lx 
    AX_m = (1.d0 - (xi-1.d0)/(xi+1.d0))      *h_Lx 
   HAX_m =  AX_m*cu(IV_m,JV_0,KV_0)
   end if 
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - -   Y-Direction
!
   HAY_m = 0.d0; HAY_p = 0.d0   
   if(j==1)                                                     then
      xi =          Yv_FG(JV_0)        *dsqrt(h_Ly)
    AY_m =  0.d0
    AY_o =-(1.d0 - (xi-1.d0)/ xi     ) * 2.d0*h_Ly 
    AY_p = (1.d0 - (xi-1.d0)/(xi+1.d0))      *h_Ly 
   HAY_p =  AY_p*cu(IV_0,JV_p,KV_0)
   end if 
   if((1<j).and.(j<Ny+1))                                       then
    AY_m =                                    h_Ly
    AY_p =                                    h_Ly 
    AY_o =-AY_m-AY_p
   HAY_m = AY_m*cu(IV_0,JV_m,KV_0)
   HAY_p = AY_p*cu(IV_0,JV_p,KV_0)   
   end if
   if(j==Ny+1)                                                  then
      xi = (1.d0 -  Yv_FG(JV_0)      ) *dsqrt(h_Ly)
    AY_p =  0.d0
    AY_o =-(1.d0 - (xi-1.d0)/ xi     ) * 2.d0*h_Ly 
    AY_m = (1.d0 - (xi-1.d0)/(xi+1.d0))      *h_Ly 
   HAY_m =  AY_m*cu(IV_0,JV_m,KV_0)    
   end if 
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - -   Z-Direction
!
   HAZ_m = 0.d0; HAZ_p = 0.d0
   if(k==1)                                                     then
      xi =          Zv_FG(KV_0)        *dsqrt(h_Lz)
    AZ_m =  0.d0
    AZ_o =-(1.d0 - (xi-1.d0)/ xi     ) * 2.d0*h_Lz 
    AZ_p = (1.d0 - (xi-1.d0)/(xi+1.d0))      *h_Lz 
   HAZ_p =  AZ_p*cu(IV_0,JV_0,KV_p)   
   end if 
   if((1<k).and.(k<Nz+1))                                       then
    AZ_m =                                    h_Lz
    AZ_p =                                    h_Lz
    AZ_o =-AZ_m-AZ_p
   HAZ_m = AZ_m*cu(IV_0,JV_0,KV_m)
   HAZ_p = AZ_p*cu(IV_0,JV_0,KV_p)   
    end if
   if(k==Nz+1)                                                  then
      xi = (1.d0 -  Zv_FG(KV_0)      ) *dsqrt(h_Lz)
    AZ_p =  0.d0
    AZ_o =-(1.d0 - (xi-1.d0)/ xi     ) * 2.d0*h_Lz 
    AZ_m = (1.d0 - (xi-1.d0)/(xi+1.d0))      *h_Lz 
   HAZ_m =  AZ_m*cu(IV_0,JV_0,KV_m)    
   end if 
      h1 = HAX_m + HAX_p
      h2 = HAY_m + HAY_p
      h3 = HAZ_m + HAZ_p
    cu(IV_0,JV_0,KV_0) =(J3(TxV(i),TyV(j),TzV(k))-h1-h2-h3)/(AX_o+AY_o+AZ_o)
end if
end subroutine Point_Gauss_Seidel_Iteration


subroutine  Convergence_Test
real*8      u_a,Error,ErrMAX,h,D2UDx2,D2UDy2,D2UDz2,Res,ResMAX,h_x,h_y_h_z
               h_x = 1.d0/dfloat(NXX_FG)
               h_y = 1.d0/dfloat(NYY_FG)
               h_z = 1.d0/dfloat(NZZ_FG)
            ResMAX = 0.d0
            ErrMAX = 0.d0
 do         i      = 2,NXX_FG
  do          j    = 2,NYY_FG
   do           k  = 2,NZZ_FG
               u_a = dexp(a_*Xv_FG(i) + b_*Yv_FG(j) + c_*Zv_FG(k))
             Error = dabs(u_a - hatU(i,j,k))
          if(Error > ErrMAX)    ErrMAX = Error
                 h = hatU(i  ,j  ,k  )                         *2.d0
            D2uDx2 =(hatU(i-1,j  ,k  ) - h + hatU(i+1,j  ,k  ))/h_x**2  
            D2uDy2 =(hatU(i  ,j-1,k  ) - h + hatU(i  ,j+1,k  ))/h_y**2
            D2uDz2 =(hatU(i  ,j  ,k-1) - h + hatU(i  ,j  ,k+1))/h_z**2
               Res = dabs(RHSF(i,j,k) - D2uDx2 - D2uDy2 - D2uDz2)
            if(Res > ResMAX)    ResMAX = Res
   end do
  end do
 end do
                     call TimeS(ITime)
write(*,1) The_number_of_multigrid_iterations,ResMAX/Res_000,ErrMAX,ITime-Time0
write(1,1) The_number_of_multigrid_iterations,ResMAX/Res_000,ErrMAX,ITime-Time0
1 format(i8,2(2x,D12.5),2x,i8)
end subroutine Convergence_Test


subroutine  TimeS(The_Time)
integer     IHR,IMIN,ISEC,I100TH,The_Time
call GETTIM(IHR,IMIN,ISEC,I100TH)
The_Time  = IHR*3600 + IMIN*60 + ISEC
end subroutine TimeS

END MODULE RMT_3D_2020_Example_16_2





















