MODULE  RMT_3D_2020_Example_16_1
integer The_number_of_smoothing_iterations, Time0 
integer The_number_of_multigrid_iterations, MGI_DC, MGI
real*8  a_,b_,c_,Res_000,Err_000

real*8, allocatable, save ::   cu(:,:,:),   J3_(:,:,:)
real*8, allocatable, save :: hatU(:,:,:),  RHSF(:,:,:)
real*8, allocatable, save :: BCX0(:,:),    BCX1(:,:)
real*8, allocatable, save :: BCY0(:,:),    BCY1(:,:)
real*8, allocatable, save :: BCZ0(:,:),    BCZ1(:,:)
real*8, allocatable, save :: BUX0(:,:),    BUX1(:,:)
real*8, allocatable, save :: BUY0(:,:),    BUY1(:,:)
real*8, allocatable, save :: BUZ0(:,:),    BUZ1(:,:)

END MODULE RMT_3D_2020_Example_16_1
