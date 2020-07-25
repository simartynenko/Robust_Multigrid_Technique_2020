MODULE   RMT_3D_2020_Example_13_1

real*8   hx0,hy0,hz0,Res,ResMAX,Res_000,Error,ErrMAX,Err_000,a_,b_,c_
integer  NX_block,NY_block,NZ_block,q,q_max,LevelC,Time0,Time1,NSIL  

real*8 , allocatable, save ::       AG(:,:),         xG(:)
real*8 , allocatable, save ::   BC_0YZ(:,:),     BC_1YZ(:,:)
real*8 , allocatable, save ::   BC_X0Z(:,:),     BC_X1Z(:,:)
real*8 , allocatable, save ::   BC_XY0(:,:),     BC_XY1(:,:)
real*8 , allocatable, save ::     corU(:,:,:),     hatU(:,:,:),     RHSF(:,:,:)
real*8 , allocatable, save :: Lambda_X(:,:,:), Lambda_Y(:,:,:), Lambda_Z(:,:,:)
real*8 , allocatable, save ::    L_A_X(:,:,:),    L_A_Y(:,:,:),    L_A_Z(:,:,:)
integer, allocatable, save ::    M_A_P(:,:,:)

END MODULE RMT_3D_2020_Example_13_1
