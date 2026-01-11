! module with declarations of global variables, arrays and structures

! Fortran code created by J.Sochor   ( https://github.com/JNSresearcher )
! e-mail: JNSresearcher@gmail.com

MODULE m_vxc2data
implicit none

! structures for calculating functions
TYPE tFun
    CHARACTER (LEN=1)               :: ex           ! source direction identifier: X, Y, Z or 0
    REAL(8)                         :: vely         ! function value
    CHARACTER (LEN=50)              :: eqn          ! function expression
    INTEGER                         :: args,nomsch  ! number of arguments and domain number
    CHARACTER (LEN=8),ALLOCATABLE   :: namex(:)     ! argument names
    REAL(8),ALLOCATABLE             :: velx(:)      ! argument values
END TYPE tFun
    TYPE (tFun), ALLOCATABLE :: Fun(:), Vmech(:), Venv(:)    !  structures of functions for sources and functions for sources movement

! structure for cells in which functions are located
TYPE tfun_nod
!  number of cells for vector and scalar functions
    INTEGER             :: numnod_Fx,  numnod_Fy,  numnod_Fz
! arrays with cells numbers for vector and scalar functions
    INTEGER,ALLOCATABLE :: nods_Fx(:), nods_Fy(:), nods_Fz(:)
    REAL(8),ALLOCATABLE :: nods_Fx_cos(:), nods_Fy_cos(:), nods_Fz_cos(:)
! speed, distance and relative shift of the source movement along the X Y Z axes
    REAL(8)             :: vel_Vmech(3),  Distance(3),  shift(3)  
    INTEGER             :: num_Vmech(3),  &    ! X Y Z axis motion function number
                           move(3)             ! sign of motion in input data
    INTEGER             :: length(3)           ! distance in integers: length = nint(Distance)
END TYPE tfun_nod
    TYPE (tfun_nod), ALLOCATABLE :: fun_nod(:)  

TYPE t_dom_C
    INTEGER, ALLOCATABLE  :: nods(:)
    INTEGER numdom_C, siznod_C
    real(8) valdom_C
END TYPE t_dom_C
TYPE (t_dom_C), ALLOCATABLE :: dom_C(:)  
INTEGER sizdom_C

!  for domain C - conduction region
INTEGER ::  siznod_C_sum, siznod_Fe, mode
INTEGER(1),ALLOCATABLE:: geoPHYS(:,:,:)     ! array with physical domain numbers
!array with number local cell of  conduction and ferrum  regions
INTEGER,   ALLOCATABLE:: geoPHYS_C(:,:,:), geoPHYS_Cd(:,:,:),geoPHYS_Fe(:,:,:)

! valPHYS - array of domain parameters
! valPHYS(:,1) - "diffusion" parameter D
! valPHYS(:,2) - "inertial"  parameter C
! valPHYS(:,3) - domain speed along the X axis
! valPHYS(:,4) - domain speed along the Y axis
! valPHYS(:,5) - domain speed along the Z axis
REAL(8), ALLOCATABLE  :: valPHYS(:,:) 

CHARACTER (LEN=6), ALLOCATABLE  :: typPHYS(:), &  ! array of symbolic domain types
                                   namePHYS(:)    ! array of domain names

END MODULE m_vxc2data
