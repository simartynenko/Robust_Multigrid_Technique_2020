MODULE RMT_3D_2020_Example_02_1

REAL*8   a_,b_,c_,Fnorm, h_x,h_y,h_z,ResLIM                               
INTEGER  NXX_FG,NYY_FG,NZZ_FG,Nx,Ny,Nz,Time0,way              
INTEGER  NX_block,NY_block,NZ_block,X_block,Y_block,Z_block,Iteration     

INTEGER, ALLOCATABLE ::  M_A_P(:,:,:)
REAL*8 , ALLOCATABLE ::      U(:,:,:),      F(:,:,:)
REAL*8 , ALLOCATABLE ::     AG(:,:)  ,     xG(:)          
REAL*8 , ALLOCATABLE :: BC_0YZ(:,:)  , BC_1YZ(:,:)
REAL*8 , ALLOCATABLE :: BC_X0Z(:,:)  , BC_X1Z(:,:)
REAL*8 , ALLOCATABLE :: BC_XY0(:,:)  , BC_XY1(:,:)

END MODULE RMT_3D_2020_Example_02_1    
