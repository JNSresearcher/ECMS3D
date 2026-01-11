! 3D-Eddy Current Calculation Program
! Fortran codes created by J.Sochor   ( https://github.com/JNSresearcher )
! e-mail: JNSresearcher@gmail.com

PROGRAM EC3D

USE m_fparser
USE m_vxc2data
IMPLICIT NONE

! --------------internal working variables--------------!
! variables to measure calculation time
INTEGER:: k0,k1,k2
REAL(8)   T, Tcalc, Tsavedata

INTEGER nCells,     &              ! number of cells with vector field sources
        nCellsGlob, &              ! total number of cells
        NcellsX, NcellsY, NcellsZ  ! number of cells with vector field sources along the axes "x", "y", "z"
 
INTEGER  i,j,k, L,LL, m, n, nl, nn,ios, mm ! m1, m2, 
INTEGER  kim ,kjm, kkm, kip,kjp,kkp,  kdz 
INTEGER  kim2 ,kjm2, kkm2, kip2,kjp2,kkp2
INTEGER  nim,njm,nkm,  nip,njp,nkp,  nijm
INTEGER  nim2,njm2,nkm2,  nip2,njp2,nkp2
INTEGER  kipjp, kipjm, kimjp, kimjm, kipkp, kipkm, kimkp, kimkm, kjpkp, kjpkm,  kjmkp, kjmkm

REAL(8) Vex,   Ve_dx,sx,dsx,   Vey,Ve_dy,sy,dsy,    Vez,Ve_dz,sz,dsz, a, s
REAL(8), ALLOCATABLE:: Vexsum(:),Ve_dxsum(:),  Veysum(:),Ve_dysum(:),  Vezsum(:),Ve_dzsum(:)

REAL(8) sxy, sxz, syz
REAL(8) Z0, Zim, Zip, Zjm, Zjp, Zkm, Zkp

! arrays for new cells of field sources as they move 
INTEGER, ALLOCATABLE ::  new_nodesX(:), new_nodesY(:), new_nodesZ(:) 
INTEGER  nij, nX, nY, nZ, Lnew, jnew, inew, flag_move,  movestop(3)

!                               base  
INTEGER, ALLOCATABLE::  irow(:), jcol(:)
REAL(8), ALLOCATABLE :: valA(:)  
INTEGER                 num_nz, num_nzU  ! for vect and scalar potencials

TYPE t_div_C 
    INTEGER, ALLOCATABLE:: irow(:),jcol(:),bndx(:),bndy(:),bndz(:)
    INTEGER siz_nods, nnz
    real(8), ALLOCATABLE::valC(:),Uf(:),divJ(:),Jxsum(:),Jysum(:),Jzsum(:)
END TYPE t_div_C
TYPE (t_div_C), ALLOCATABLE :: div_C(:)  

! arrays for storing numbers of boundary cells
! for vector potential 
integer num_bndX, num_bndY, num_bndZ
integer, allocatable:: cel_bndX(:),  cel_bndY(:),  cel_bndZ(:)
! for scalar potential
integer num_bndUx, num_bndUy, num_bndUz
integer, allocatable::  cel_bndUx(:),  cel_bndUy(:),  cel_bndUz(:)

integer num_nzX, num_nzY, num_nzZ

! variables to control the output of calculation results
INTEGER  Ntime, Npoint, Nout, Nprint, Nprint_display, Nout_display

! Uaf - solution vector, Jaf(:) - right hand side vector
REAL(8),ALLOCATABLE ::  Uaf(:),   Jaf(:), Bfex(:,:,:),Bfey(:,:,:),Bfez(:,:,:)   ! arrays for vector fields (dimension nCells)

! Jsum - the sum of the speed, transformer and scalar components of the eddy current
REAL(8),ALLOCATABLE :: Jxsum(:), Jysum(:), Jzsum(:), Uf(:),   divJ(:)
                       
!---------------------------------------!
! input data from SUBROUTINE vxc2data   !
!---------------------------------------!
REAL(8) :: delta(3)         ! array grid spacing along X, Y and Z   
REAL(8) :: dt,   &          ! time step
           Time, &          ! stop time
           dtt              ! jump step

! for update rhs
REAL(8) :: valrhs(13), valU(40)
INTEGER :: colrhs(13), colU(40), Lfi,  Lx !, d, Ly, Lz 

INTEGER :: sdx,sdy,sdz      ! number of cells along X,  Y and  Z   
INTEGER :: nsub,     &      ! number of physical domains
           nsub_air, &      ! number of environment domains
           nsub_glob        ! total number of domains
INTEGER :: numfun,   &      ! number of functions for calculating external sources
           numMech,  &      ! number of functions for calculating movements of external sources
           numVenv

! solv - character string corresponding to the name of methods: 'BCG' 'GMR' or 'PRD'
CHARACTER(LEN=3) :: solv 
REAL(8)      :: BND(3,2)  ! values  boundary conditions

CHARACTER(16):: files      ! name for output files

!//=============== declare variables BCG for Divergence
real(8) tol_bcg 
integer itmax_bcg

!// =====================declare the variables of pardiso====================
integer(kind=8)           :: pt(64)   !// kind=8:x64; kind=4:win32
integer                   :: maxfct, mnum, mtype, phase, nrhs, error, msglvl
integer, allocatable      :: perm(:)
integer                   :: iparm(64)

!//=======================declare the variables of GMRES
real(8), allocatable ::  au_lu(:), vv_lu(:,:), w_lu(:), rhs(:)
integer, allocatable ::  ju_lu(:), jau_lu(:),  iw_lu(:,:), jw_lu(:),    iperm(:)
integer ierr, im_lu,  iwk, nscale,  iout_lu,   mbloc
REAL(8) ::  tol, eps,  permtol      ! convergence criterion
INTEGER ::  iter, itmax, elem       ! actual number of iterations,  maximum number of iterations,  maximum of elements in row

!// =====================the variables of LU====================
im_lu = 10
iout_lu = 6
!//===============initialization variables BCGSTAB for Divergence====
tol_bcg = 5.0d-3 
itmax_bcg = 10000

!-------------------------------!
! call to data entry subroutine !
!-------------------------------!
CALL vxc2data ( delta,    dt,        Time,     dtt,                   &
                sdx,    sdy,      sdz,       nsub,     nsub_air, nsub_glob, &
                numfun, numMech,  numVenv,     &
                solv,   files,    tol, eps, elem, itmax,   BND ) 
!----------------------!
! START of calculation !
!----------------------!
nCells = sdx*sdy*sdz; kdz=sdx*sdy
sx = 1.d0/(delta(1)*delta(1)) 
sy = 1.d0/(delta(2)*delta(2))
sz = 1.d0/(delta(3)*delta(3))

dsx = 0.5d0/delta(1) ! 1/2dx
dsy = 0.5d0/delta(2) ! 1/2dy
dsz = 0.5d0/delta(3) ! 1/2dz

sxy = 1.d0/(delta(1)*delta(2)) 
sxz = 1.d0/(delta(1)*delta(3)) 
syz = 1.d0/(delta(2)*delta(3)) 

! calculated of number of cells with scalar field in conduc domains
nCellsGlob = 3*nCells + siznod_C_sum

!----------------------------------------------------!
! formation of a sparse matrix                       !
!----------------------------------------------------!
allocate (irow(nCellsGlob + 1), source=1)

IF ( siznod_C_sum /= 0 ) THEN
    allocate (Vexsum(sizdom_C),Ve_dxsum(sizdom_C),Veysum(sizdom_C),Ve_dysum(sizdom_C),Vezsum(sizdom_C),Ve_dzsum(sizdom_C), source=0.0_8 )
    allocate (Jxsum(siznod_C_sum),Jysum(siznod_C_sum),Jzsum(siznod_C_sum),  source=0.0_8)
    allocate ( divJ(siznod_C_sum), Uf(siznod_C_sum),  source=0.0_8)

    allocate (div_C(sizdom_C) )
    do i=1,sizdom_C
        j = dom_C(i)%siznod_C
        div_C(i)%siz_nods = j
        allocate(div_C(i)%Uf(j), div_C(i)%divJ(j), div_C(i)%Jxsum(j), div_C(i)%Jysum(j), div_C(i)%Jzsum(j), source=0.d0 )
        allocate(div_C(i)%irow(j+1),source=1 )
    enddo
    
    call gen_sparse_C_matrix
ENDIF

print'(*(a,i8,1x))', 'nCells=',nCells, ' nCellsGlob=',nCellsGlob, ' siznod_C_sum=',siznod_C_sum, ' numfun=' ,numfun

!detect speed
nX=0;nY=0;nZ=0
IF (numMech /= 0) THEN
    Ve_dx=0.d0; Ve_dy=0.d0; Ve_dz=0.d0
    Vex=0.d0;  Vey=0.d0;  Vez=0.d0
    do n=1,numMech
        if     (Vmech(n)%ex == 'X' .and. nX==0 ) then 
            nX = n
        elseif (Vmech(n)%ex == 'Y' .and. nY==0 ) then 
            nY = n
        elseif (Vmech(n)%ex == 'Z' .and. nZ==0 ) then 
            nZ = n
        endif
    enddo
elseIF (numVenv /= 0) THEN
    DO n=1,numVenv
        if     (Venv(n)%ex == 'X' .and. nX==0 ) then 
            nX = n
        elseif (Venv(n)%ex == 'Y' .and. nY==0 ) then 
            nY = n
        elseif (Venv(n)%ex == 'Z' .and. nZ==0 ) then 
            nZ = n
        endif
    enddo
endif

if (mode == 0 .and. ( nX /= 0 .or. nY /= 0 .or. nZ /= 0 ) ) mode = 1

CALL gen_sparse_matrix         ! generation of sparse matrix

! magnetic density flow
allocate (  Bfex(sdx,sdy,sdz),Bfey(sdx,sdy,sdz),Bfez(sdx,sdy,sdz), source=0.0_8)

!-------------------------------!
! preparation of time intervals !
!-------------------------------!
T=0.
Ntime=0                              ! counter of each point with step dt
Nout_display=nint(Time/(100.0*DT))   ! jump size for display of symbol  ">"
Nprint_display=Ntime + Nout_display  ! counter of points for display symbol ">" with jump Nout_display
Npoint=0                             ! counter of calculated points with step dtt
Nout=nint(DTT/DT)                    ! jump size for output to files
Nprint = Ntime + Nout                ! point counter with step Nout for output to files with step Nout

!-----------------------------------!
! preparing arrays of vector fields !
!-----------------------------------!
allocate (Jaf(nCellsGlob), Uaf(nCellsGlob),  source=0.0_8)

!---------------------------------------------------------------!
! preparation of arrays for calculating the movement of sources !
!---------------------------------------------------------------!
flag_move = 0 ! global flag of motion in input data

IF (numfun /=0 ) THEN
    DO i=1,numfun
    ! motion search processing
        fun_nod(i)%Distance = 0.d0 
    
        ! unlock
        IF (fun_nod(i)% num_Vmech(1) == 0  .and. fun_nod(i)%move(1) /=0) THEN ! it is not function
        !  distance traveled in time dt with speed Vx in fractions of cell size along dX
            fun_nod(i)%shift(1) = fun_nod(i)%vel_Vmech(1)*dt/delta(1)   !Vx*dt/dX
            flag_move = 1
        ELSEIF (  fun_nod(i)%move(1) /=0) THEN 
            flag_move = 1
        ENDIF
        
        IF (fun_nod(i)% num_Vmech(2) == 0  .and. fun_nod(i)%move(2) /=0) THEN ! it is not function
            fun_nod(i)%shift(2) = fun_nod(i)%vel_Vmech(2)*dt/delta(2)
            flag_move = 1; 
        ELSEIF (  fun_nod(i)%move(2) /=0) THEN 
            flag_move = 1
        ENDIF
    
        IF (fun_nod(i)% num_Vmech(3) == 0 .and. fun_nod(i)%move(3) /=0 ) THEN ! it is not function, but not a shift either
            fun_nod(i)%shift(3) = fun_nod(i)%vel_Vmech(3)*dt/delta(3)
            flag_move = 1; 
        ELSEIF (  fun_nod(i)%move(3) /=0) THEN 
            flag_move = 1
        ENDIF
    ENDDO
! initial allocation of array memory for source nodes, regardless of whether there is movement or not
    NcellsX = 0; NcellsY = 0; NcellsZ = 0; 
    DO n=1,numfun
        NcellsX = NcellsX + fun_nod(n)%numnod_Fx
        NcellsY = NcellsY + fun_nod(n)%numnod_Fy
        NcellsZ = NcellsZ + fun_nod(n)%numnod_Fz
    ENDDO
    IF (NcellsX /= 0) ALLOCATE (new_nodesX(NcellsX), source=0)
    IF (NcellsY /= 0) ALLOCATE (new_nodesY(NcellsY), source=0)
    IF (NcellsZ /= 0) ALLOCATE (new_nodesZ(NcellsZ), source=0)
ENDIF !numfun

! filling source node arrays if there is no movement
IF (flag_move == 0 .and. numfun /=0 ) THEN
    NcellsX = 0; NcellsY = 0; NcellsZ = 0; ! counters
    DO n=1,numfun
        DO k = 1,fun_nod(n)%numnod_Fx
            m = fun_nod(n)%nods_Fx(k)
            NcellsX = NcellsX + 1
            new_nodesX(NcellsX) = m
        ENDDO
        DO k = 1,fun_nod(n)%numnod_Fy
            m = fun_nod(n)%nods_Fy(k) - nCells
            NcellsY = NcellsY + 1
            new_nodesY(NcellsY) = m !+ 
        ENDDO
        DO k = 1,fun_nod(n)%numnod_Fz
            m = fun_nod(n)%nods_Fz(k) - 2*nCells
            NcellsZ = NcellsZ + 1
            new_nodesZ(NcellsZ) = m !+ 
        ENDDO 
    ENDDO
ENDIF

call execute_command_line ('mkdir '//trim(files), exitstat=i)
PRINT *, 'created dir '//trim(files)//': status=', i

!!!!!!  START time counter
CALL system_clock(k1,k0)
Tsavedata = 0.0
movestop = 1  ! flag of local start-stop movement in all coordinates
write(*,'( "start solver: ", a)') solv

if (solv == 'GMR') then
    nscale=4 
    allocate(jau_lu(nscale*num_nz),  ju_lu(nCellsGlob),  iw_lu(nCellsGlob, 3), jw_lu(2*nCellsGlob),  source=0)
    allocate( au_lu(nscale*num_nz),  vv_lu(nCellsGlob,20), w_lu(nCellsGlob+1), rhs(nCellsGlob), source=0.d0)

    iwk=nscale*num_nz
    permtol = 0.1d0  !  0 for never permute; or values 0.1 to 0.01
    mbloc = nCellsGlob 
    allocate( iperm(2*nCellsGlob), source=0)
! ju_lu diagonal indices
    do i = 1, nCellsGlob
        do j = irow(i), irow(i+1)-1
            if (jcol(j) == i) then
                ju_lu(i) =  j
                exit
            endif
        enddo
    enddo
    
    call ilutp(nCellsGlob,valA,jcol,irow, elem, tol, permtol,mbloc, au_lu, jau_lu, ju_lu, iwk, w_lu, jw_lu,iperm,ierr)

    if (ierr == 0) then
        print*, 'done'
    elseif (ierr > 0) then
            ! Zero pivot
        print*, 'ierr > 0, call lusol,  zero pivot ierr=', ierr
    elseif (ierr == -1) then
    ! The input matrix is not formatted correctly
        print*, 'ierr == -1, The input matrix is not formatted correctly'
        stop
    else if (ierr == -2 .or. ierr == -3) then
        print*,'au_lu(iwk) and jau_lu(iwk) are too small  new iwk='
        stop
    else if (ierr == -4) then
    ! Illegal value for elem - reset and try again
        elem = n
        print*,'Illegal value for ELEM elem - reset and try again ELEM=',ELEM
        stop
    else if (ierr == -5) then
    ! Zero row encountered
        print*,'ierr == -5, Zero row encountered vsretil'
        stop
    else
    ! We should never get here, but just in case
        print*, 'error ILUT encountered, an unknown error.  The error code from the ILUT routine is' 
        stop
    end if
endif

if (solv == 'PRD') then
    allocate( perm( nCellsGlob ), source=0  ) 
    pt = 0  !// pointer initialization
    maxfct = 1; mnum = 1; 
    ! !// mtype = -2: symmetric nonpositive definite matrix, mtype = 2: symmetric positive definite matrix
    mtype = 11  !real and unsymmetric matrix
    nrhs = 1
    iparm(1) = 0  !// iparm use default values
    iparm(5) = 1  
    error = 0  
    msglvl = 0
    phase = 12  !// LU decompose
    call pardiso( pt, maxfct, mnum, mtype, phase, nCellsGlob, valA, irow, jcol, perm, nrhs, iparm, msglvl, Jaf, Uaf, error )
    write(*, '(a, i9, a, i3)' ) 'pardiso for base: iparm15+16+63', iparm(15) + iparm(16) + iparm(63), ' error=',error
    phase = 33  !// solve equations
endif

!------------!
! Main loop  !
!------------!
2000 CONTINUE   

!-------------------------------------!
! Calculation of independent sources  !
!-------------------------------------!
IF (numfun /= 0) THEN
    CALL initf (numfun) 
    DO i=1,numfun
        DO m=1,Fun(i)%args
            IF ( trim(Fun(i)%namex(m)) == 'T') THEN
                Fun(i)%velx(m)=T
            ENDIF
        ENDDO
        CALL parsef (i, trim(fun(i)%eqn), fun(i)%namex(1:fun(i)%args) )
        fun(i)%vely = evalf (i, fun(i)%velx(1:fun(i)%args)) * 0.12566370964050292d-5  ! * mu0
    ENDDO   
ENDIF
!-------------------------------------------------!
! Calculation of the speed of movement of sources !
!-------------------------------------------------!
IF (numMech /= 0) THEN
    CALL initf (numMech) 
    DO i=1,numMech
        DO m=1,Vmech(i)%args
            IF ( trim(Vmech(i)%namex(m)) == 'T') THEN
                Vmech(i)%velx(m)=T
            ENDIF
        ENDDO
        CALL parsef (i, trim(Vmech(i)%eqn), Vmech(i)%namex(1:Vmech(i)%args) )
        Vmech(i)%vely = evalf (i, Vmech(i)%velx(1:Vmech(i)%args))
    ENDDO
ENDIF
!-----------------------------------------------------!
! Calculation of the speed of movement of environment !
!-----------------------------------------------------!
IF (numVenv /= 0) THEN
    CALL initf (numVenv) 
    DO i=1,numVenv
        DO m=1,Venv(i)%args
            IF ( trim(Venv(i)%namex(m)) == 'T') THEN
                Venv(i)%velx(m)=T
            ENDIF
        ENDDO
        CALL parsef (i, trim(Venv(i)%eqn), Venv(i)%namex(1:Venv(i)%args) )
        Venv(i)%vely = evalf (i, Venv(i)%velx(1:Venv(i)%args))
        
        if (Venv(i)%ex == 'X') then     ! Vex 
            valPHYS( Venv(i)%nomsch,3)  = Venv(i)%vely * valPHYS( Venv(i)%nomsch,2)
        elseif (Venv(i)%ex == 'Y') then ! Vey 
            valPHYS( Venv(i)%nomsch,4)  = Venv(i)%vely * valPHYS( Venv(i)%nomsch,2)
        elseif (Venv(i)%ex == 'Z') then ! Vez
            valPHYS( Venv(i)%nomsch,5)  = Venv(i)%vely * valPHYS( Venv(i)%nomsch,2)
        endif
    ENDDO
ENDIF

!--------------------------------------------------!
! Distribution of independent sources into domains !
!--------------------------------------------------!
IF (flag_move == 1) THEN 
    ! there are moving sources ! 
    IF (siznod_C_sum /= 0  ) THEN      ! only inertial sources 
        Jaf = 0.0d0
    ENDIF
    
    ! calculation of new coordinates of nodes of independent sources 
    IF (numfun /= 0) THEN
        NcellsX=0; NcellsY=0; NcellsZ=0
        DO n=1,numfun
            CALL motion_calc    ! motion calculation
            a = Fun(n)%vely
            DO k = 1,fun_nod(n)%numnod_Fx
                m = fun_nod(n)%nods_Fx(k)
                CALL new_m !(m)
                Jaf(m) = a*  fun_nod(n)%nods_Fx_cos(k) 
                NcellsX = NcellsX + 1
                new_nodesX(NcellsX) = m
            ENDDO
            DO k = 1,fun_nod(n)%numnod_Fy
                m = fun_nod(n)%nods_Fy(k) - nCells
                CALL new_m 
                NcellsY = NcellsY + 1
                new_nodesY(NcellsY) = m 
                m = m + nCells
                Jaf(m) = a*  fun_nod(n)%nods_Fy_cos(k) 
            ENDDO
            DO k = 1,fun_nod(n)%numnod_Fz
                m = fun_nod(n)%nods_Fz(k) - 2*nCells
                CALL new_m 
                NcellsZ = NcellsZ + 1
                new_nodesZ(NcellsZ) = m 
                m = m + 2*nCells
                Jaf(m) = a*  fun_nod(n)%nods_Fz_cos(k) 
            ENDDO                  
        ENDDO
    ENDIF
ELSE
    !--------------------!
    ! no moving sources  !
    !--------------------!
    DO n=1,numfun
        a = Fun(n)%vely
        DO k = 1,fun_nod(n)%numnod_Fx
            m = fun_nod(n)%nods_Fx(k)                 ! k - local number m - global number
            Jaf(m) = a *  fun_nod(n)%nods_Fx_cos(k) 
            m = fun_nod(n)%nods_Fy(k)
            Jaf(m) = a * fun_nod(n)%nods_Fy_cos(k) 
            m = fun_nod(n)%nods_Fz(k)
            Jaf(m) = a * fun_nod(n)%nods_Fz_cos(k) 
        ENDDO
    ENDDO
ENDIF

IF (sizdom_C /=0) THEN
    m = 0
    do LL = 1, sizdom_C
        do n = 1, dom_C(LL)%siznod_C
            m = m + 1        ! counter nodes in domain of conductivity
            L = dom_C(LL)%nods(n)
            nn  = L + nCells
            nl  = L + 2*nCells
            Jaf(L)  = dom_C(LL)%valdom_C * Uaf(L)  + div_C(LL)%Jxsum(n)  
            Jaf(nn) = dom_C(LL)%valdom_C * Uaf(nn) + div_C(LL)%Jysum(n)   ! Jb = 2C/dt * U_ + Ib  
            Jaf(nl) = dom_C(LL)%valdom_C * Uaf(nl) + div_C(LL)%Jzsum(n) 
            s=0.d0
            do i = irow(3*nCells + m), irow(3*nCells + m +1)-1   
                k = jcol(i)
                if ( k < 3*nCells+1) then
                    s = s + valA(i) * Uaf(k) 
                endif
            enddo
            Jaf( 3*nCells + m ) =  s 
        enddo
    enddo
    
    Jaf(cel_bndX) = 0.0d0 
    Jaf(cel_bndY) = 0.0d0 
    Jaf(cel_bndZ) = 0.0d0  
 
!======================= update right parts of vector field=== and scalar field =========================
    IF (numMech /= 0) THEN
    ! Automatically sets Venv environment speed if field source speed is set
        do i = 1,sizdom_C
            j = dom_C(i)%numdom_C 
            if (nX /= 0 )  valPHYS(j,3) = -Vmech(nX)%vely*valPHYS(j,2) !  valPHYS(numdom_C,2) = sigma*mu0
            if (nY /= 0 )  valPHYS(j,4) = -Vmech(nY)%vely*valPHYS(j,2)
            if (nZ /= 0 )  valPHYS(j,5) = -Vmech(nZ)%vely*valPHYS(j,2)
        enddo
    endif

    IF (numVenv /= 0  .or. numMech /= 0) THEN
!  constants for the velocity shift of the vector valrhs(1:7) and scalar valU field   
        do i = 1,sizdom_C
            j = dom_C(i)%numdom_C
            if ( nx /=0) then
                Ve_dxsum(i)=valPHYS(j,3)/(2.d0*delta(1)) ! valPHYS(m,3) = Vex * mu0*sigma
                Vexsum(i) = valPHYS(j,3)/valPHYS(j,2)    ! Vex = (Vex * mu0*sigma)/mu0*sigma
            endif
            if ( ny /=0) then
                Ve_dysum(i)=valPHYS(j,4)/(2.d0*delta(2))
                Veysum(i) = valPHYS(j,4)/valPHYS(j,2)
            endif
            if ( nz /=0) then
                Ve_dzsum(i)=valPHYS(j,5)/(2.d0*delta(3))
                Vezsum(i) = valPHYS(j,5)/valPHYS(j,2)
            endif
        enddo

        m = 0
        do LL = 1, sizdom_C
            do n = 1, dom_C(LL)%siznod_C
                m = m + 1
                L = dom_C(LL)%nods(n)
        
                k= ceiling(REAL(L)/( REAL(sdx*sdy) ) )
                IF (k == 1) THEN
                    nij = L
                ELSE
                    nij = L - (k-1)*sdx*sdy
                ENDIF
                j = ceiling( REAL(nij) / REAL(sdx) )  
                i = nij - (j - 1) * sdx  
        
                nim  = geoPHYS_C(i-1,j,k);   nip  = geoPHYS_C(i+1,j,k)
                njm  = geoPHYS_C(i,j-1,k);   njp  = geoPHYS_C(i,j+1,k)
                nkm  = geoPHYS_C(i,j,k-1);   nkp  = geoPHYS_C(i,j,k+1)
                kim = L-1;  kjm = L-sdx;  kkm = L-kdz;
                kip = L+1;  kjp = L+sdx;  kkp = L+kdz;
                    ! vector field stretching due to velocity
                if (  nim /= 0 .and. nip /= 0 .and. njm /= 0 .and. njp /= 0 .and. nkm /= 0 .and. nkp /= 0 ) then 
                    Lx=0;  colrhs = 0;  valrhs = 0.d0 ! 
                    if ( nx /=0 ) then
                        Ve_dx = Ve_dxsum(LL)
                        colrhs(Lx+1:Lx+2) =       [kip,  kim ]
                        valrhs(Lx+1:Lx+2) = Ve_dx*[-1d0, 1d0]
                        Lx = Lx + 2  
                    endif
                    if ( ny /=0 ) then
                        Ve_dy = Ve_dysum(LL)
                        colrhs(Lx+1:Lx+2) =       [kjp, kjm  ]
                        valrhs(Lx+1:Lx+2) = Ve_dy*[-1d0, 1d0]
                        Lx = Lx + 2
                    endif
                    if ( nz /=0 ) then
                        Ve_dz = Ve_dzsum(LL)
                        colrhs(Lx+1:Lx+2) =      [kkp, kkm  ]
                        valrhs(Lx+1:Lx+2) = Ve_dz*[-1d0, 1d0]
                        Lx = Lx + 2
                    endif
                    Jaf( colrhs(1:Lx)) = Jaf( colrhs(1:Lx) ) + valrhs(1:Lx) * Uaf( colrhs(1:Lx) )
                    colrhs(1:Lx) = colrhs(1:Lx) + nCells
                    Jaf( colrhs(1:Lx)) = Jaf( colrhs(1:Lx) ) + valrhs(1:Lx) * Uaf( colrhs(1:Lx) )
                    colrhs(1:Lx) = colrhs(1:Lx) + nCells
                    Jaf( colrhs(1:Lx)) = Jaf( colrhs(1:Lx) ) + valrhs(1:Lx) * Uaf( colrhs(1:Lx) )
                endif
            
                kim2= L-2;  kjm2 = L-2*sdx; kkm2 = L-2*kdz;
                kip2= L+2;  kjp2 = L+2*sdx; kkp2 = L+2*kdz; 
                kimjm = L - 1 - sdx; kimjp = L - 1 + sdx; kimkm = L - 1 - kdz;  kimkp = L - 1 + kdz;
                kipjm = L + 1 - sdx; kipjp = L + 1 + sdx; kipkm = L + 1 - kdz;  kipkp = L + 1 + kdz; 
                kjmkm = L - sdx - kdz;  kjmkp = L - sdx + kdz;
                kjpkm = L + sdx - kdz;  kjpkp = L + sdx + kdz;

                Lfi = 0; colU = 0; valU = 0.d0
                call update_rhs
                !adding velocity components of the vector field to the central node Vx*ddA/ddx of the scalar field
                Jaf( 3*nCells + m ) = Jaf( 3*nCells + m ) + dot_product(valU(1:Lfi), Uaf(colU(1:Lfi)) ) 
            ENDDO
        enddo
    
    ENDIF  ! numVenv /=0  
ENDIF  !  siznod_C /=0

        !---------------------------!
        ! SOLVERS for vector fields !
        !---------------------------!
if     (solv == 'BCG') then
    CALL sprsBCGstabwr (valA, irow, jcol, nCellsGlob, Jaf, Uaf,  tol, itmax, iter)
elseif (solv == 'PRD') then
    call pardiso( pt, maxfct, mnum, mtype, phase, nCellsGlob, valA, irow, jcol, perm, nrhs, iparm, msglvl, Jaf, Uaf, error )
elseif (solv == 'GMR') then
    rhs=Jaf
    call pgmres (nCellsGlob,   im_lu, rhs, Uaf,  vv_lu, eps, itmax, iout_lu, valA, jcol, irow, au_lu, jau_lu, ju_lu, ierr)
else
    print*,'not solver ', trim(solv)
    stop
endif
        ! calculation induction flow  
m=0
do k=1,sdz;  do j=1,sdy;  do i=1,sdx            
    m = m+1
    nim = m - 1;  njm = m - sdx;  nkm = m - kdz;
    nip = m + 1;  njp = m + sdx;  nkp = m + kdz;
    sx=0.5d0;  sy=0.5d0;  sz=0.5d0; 
if ( i>1 .and. i<sdx .and. j>1 .and. j<sdy .and.k>1 .and. k<sdz ) then 
    IF (siznod_Fe /=0) THEN
        n   = geoPHYS(i,j,k)
        kim = geoPHYS(i-1,j,k);  kip = geoPHYS(i+1,j,k);
        kjm = geoPHYS(i,j-1,k);  kjp = geoPHYS(i,j+1,k);
        kkm = geoPHYS(i,j,k-1);  kkp = geoPHYS(i,j,k+1);
        Z0  = valPHYS(n,1); 
        Zim = valPHYS(kim,1); Zip = valPHYS(kip,1); 
        Zjm = valPHYS(kjm,1); Zjp = valPHYS(kjp,1); 
        Zkm = valPHYS(kkm,1); Zkp = valPHYS(kkp,1); 
        BOUNDARY_B: BLOCK
            real(8) Zbase(6) 
            integer i_loc
            logical mask(6)
            Zbase = [Zim, Zip, Zjm, Zjp, Zkm, Zkp]
            do i_loc = 1, 6
                mask(i_loc) = les(Z0, Zbase(i_loc) )
            enddo
            if  (mask(1))  then; sx=1.d0; nim=m; endif
            if  (mask(2))  then; sx=1.d0; nip=m; endif
            if  (mask(3))  then; sy=1.d0; njm=m; endif
            if  (mask(4))  then; sy=1.d0; njp=m; endif
            if  (mask(5))  then; sz=1.d0; nkm=m; endif
            if  (mask(6))  then; sz=1.d0; nkp=m; endif
        END BLOCK BOUNDARY_B
        Bfex(i,j,k) = sy*(Uaf(2*nCells+njp)- Uaf(2*nCells+njm))/delta(2) - sz*( Uaf(nCells+nkp)   - Uaf(nCells+nkm))  /delta(3)  
        Bfey(i,j,k) = sz*(Uaf(nkp)         - Uaf(nkm))         /delta(3) - sx*( Uaf(2*nCells+nip) - Uaf(2*nCells+nim))/delta(1) 
        Bfez(i,j,k) = sx*(Uaf(nCells+nip)  - Uaf(nCells+nim))  /delta(1) - sy*( Uaf(njp)          - Uaf(njm))         /delta(2)
    else  ! siznod_Fe ==0
        Bfex(i,j,k) = 0.5d0*(Uaf(2*nCells+njp)- Uaf(2*nCells+njm))/delta(2) - 0.5d0*( Uaf(nCells+nkp)   - Uaf(nCells+nkm))  /delta(3)  
        Bfey(i,j,k) = 0.5d0*(Uaf(nkp)         - Uaf(nkm))         /delta(3) - 0.5d0*( Uaf(2*nCells+nip) - Uaf(2*nCells+nim))/delta(1) 
        Bfez(i,j,k) = 0.5d0*(Uaf(nCells+nip)  - Uaf(nCells+nim))  /delta(1) - 0.5d0*( Uaf(njp)          - Uaf(njm))         /delta(2)
    endif 
else 
    if (i==1) nim = m; if (i==sdx) nip  =m;      
    if (j==1) njm = m; if (j==sdy) njp = m;   
    if (k==1) nkm = m; if (k==sdz) nkp = m;  
    Bfex(i,j,k) = 0.5d0*(Uaf(2*nCells+njp)- Uaf(2*nCells+njm))/delta(2) - 0.5d0*( Uaf(nCells+nkp)   - Uaf(nCells+nkm))  /delta(3)  
    Bfey(i,j,k) = 0.5d0*(Uaf(nkp)         - Uaf(nkm))         /delta(3) - 0.5d0*( Uaf(2*nCells+nip) - Uaf(2*nCells+nim))/delta(1) 
    Bfez(i,j,k) = 0.5d0*(Uaf(nCells+nip)  - Uaf(nCells+nim))  /delta(1) - 0.5d0*( Uaf(njp)          - Uaf(njm))         /delta(2)
endif

enddo;  enddo;  enddo  ! end calculation induction flow 

        !---------------------------------------------!
        ! calculation of sources of inertial branches !
        !---------------------------------------------!
IF (sizdom_C /=0) THEN 
    m = 0
    do LL = 1, sizdom_C
        do n = 1, dom_C(LL)%siznod_C
            m = m + 1               
            L = dom_C(LL)%nods(n)   
            nn  = L + nCells
            nl  = L + 2*nCells

            k= ceiling(REAL(L)/( REAL(sdx*sdy) ) )
            IF (k == 1) THEN
                nij = L
            ELSE
                nij = L - (k-1)*sdx*sdy
            ENDIF
            j = ceiling( REAL(nij) / REAL(sdx) )  
            i = nij - (j - 1) * sdx  
            Jaf(L)  = dom_C(LL)%valdom_C * Uaf(L)  - Jaf(L)  
            Jaf(nn) = dom_C(LL)%valdom_C * Uaf(nn) - Jaf(nn)  ! Jb = 2C/dt * U_ - Ib  
            Jaf(nl) = dom_C(LL)%valdom_C * Uaf(nl) - Jaf(nl)

            mm   = geoPHYS_C(i,j,k)
            kim  = geoPHYS_C(i-1,j,k);   kip  = geoPHYS_C(i+1,j,k)
            kjm  = geoPHYS_C(i,j-1,k);   kjp  = geoPHYS_C(i,j+1,k)
            kkm  = geoPHYS_C(i,j,k-1);   kkp  = geoPHYS_C(i,j,k+1)
            kim2 = geoPHYS_C(i-2,j,k);   kip2 = geoPHYS_C(i+2,j,k)
            kjm2 = geoPHYS_C(i,j-2,k);   kjp2 = geoPHYS_C(i,j+2,k)
            kkm2 = geoPHYS_C(i,j,k-2);   kkp2 = geoPHYS_C(i,j,k+2)
            if     (kim == 0) then
                sx = 0.5d0*( -Uaf(kip2) + 4.d0*Uaf(kip) - 3.d0*Uaf(mm) )/delta(1)
            elseif (kip == 0) then
                sx = 0.5d0*(  Uaf(kim2) - 4.d0*Uaf(kim) + 3.d0*Uaf(mm) )/delta(1)
            else
                sx = 0.5d0*( Uaf(kip) - Uaf(kim) )/delta(1)
            endif
            if (kjm == 0) then
                sy = 0.5d0*( -Uaf(kjp2) + 4.d0*Uaf(kjp) - 3.d0*Uaf(mm) )/delta(2)
            elseif (kjp == 0) then
                sy = 0.5d0*(  Uaf(kjm2) - 4.d0*Uaf(kjm) + 3.d0*Uaf(mm) )/delta(2)
            else
                sy = 0.5d0*( Uaf(kjp) - Uaf(kjm) )/delta(2)
            endif
            if (kkm == 0) then
                sz = 0.5d0*( -Uaf(kkp2) + 4.d0*Uaf(kkp) - 3.d0*Uaf(mm) )/delta(3)
            elseif (kkp == 0) then
                sz = 0.5d0*(  Uaf(kkm2) - 4.d0*Uaf(kkm) + 3.d0*Uaf(mm) )/delta(3)
            else
                sz = 0.5d0*( Uaf(kkp) - Uaf(kkm) )/delta(3)
            endif
            ! current from scalar potential
            nn = dom_C(LL)%numdom_C
            
            div_C(LL)%Jxsum(n) = valPHYS(nn,2)*sx 
            div_C(LL)%Jysum(n) = valPHYS(nn,2)*sy
            div_C(LL)%Jzsum(n) = valPHYS(nn,2)*sz 

            IF (numVenv /= 0  .or. numMech /= 0) THEN
                sx=0d0; sy=0d0; sz=0d0
                if ( nx /=0) then
                    nim = L - 1;  nip = L + 1; 
                    nim2= L - 2;  nip2= L + 2; 
                    if     (kim == 0) then
                        sx = 0.5d0*( -Uaf(nip2) + 4.d0*Uaf(nip) - 3.d0*Uaf(L) )/delta(1)
                        sy = 0.5d0*( -Uaf(nip2+nCells) + 4.d0*Uaf(nip+nCells) - 3.d0*Uaf(L+nCells) )/delta(1)
                        sz = 0.5d0*( -Uaf(nip2+2*nCells) + 4.d0*Uaf(nip+2*nCells) - 3.d0*Uaf(L+2*nCells) )/delta(1)
                    elseif (kip == 0) then
                        sx = 0.5d0*(  Uaf(nim2) - 4.d0*Uaf(nim) + 3.d0*Uaf(L) )/delta(1)
                        sy = 0.5d0*(  Uaf(nim2+nCells) - 4.d0*Uaf(nim+nCells) + 3.d0*Uaf(L+nCells) )/delta(1)
                        sz = 0.5d0*(  Uaf(nim2+2*nCells) - 4.d0*Uaf(nim+2*nCells) + 3.d0*Uaf(L+2*nCells) )/delta(1)
                    else
                        sx = 0.5d0*( Uaf(nip)          - Uaf(nim) )        /delta(1)
                        sy = 0.5d0*( Uaf(nCells+nip)   - Uaf(nCells+nim))  /delta(1)
                        sz = 0.5d0*( Uaf(2*nCells+nip) - Uaf(2*nCells+nim))/delta(1)
                    endif
                    ! velocity component of current for Vx 
                    div_C(LL)%Jxsum(n) = div_C(LL)%Jxsum(n) +  sx*valPHYS(nn,3) 
                    div_C(LL)%Jysum(n) = div_C(LL)%Jysum(n) +  sy*valPHYS(nn,3)
                    div_C(LL)%Jzsum(n) = div_C(LL)%Jzsum(n) +  sz*valPHYS(nn,3) 
                endif
                sx=0d0; sy=0d0; sz=0d0
                if ( ny /=0) then
                    njm  = L - sdx;  njp  = L + sdx;    
                    njm2 = L-2*sdx;  njp2 = L+2*sdx;    
                    if     (kjm == 0) then
                        sx = 0.5d0*( -Uaf(njp2) + 4.d0*Uaf(njp) - 3.d0*Uaf(L) )/delta(2)
                        sy = 0.5d0*( -Uaf(njp2+nCells) + 4.d0*Uaf(njp+nCells) - 3.d0*Uaf(L+nCells) )/delta(2)
                        sz = 0.5d0*( -Uaf(njp2+2*nCells) + 4.d0*Uaf(njp+2*nCells) - 3.d0*Uaf(L+2*nCells) )/delta(2)
                    elseif (kjp == 0) then
                        sx = 0.5d0*(  Uaf(njm2) - 4.d0*Uaf(njm) + 3.d0*Uaf(L) )/delta(2)
                        sy = 0.5d0*(  Uaf(njm2+nCells) - 4.d0*Uaf(njm+nCells) + 3.d0*Uaf(L+nCells) )/delta(2)
                        sz = 0.5d0*(  Uaf(njm2+2*nCells) - 4.d0*Uaf(njm+2*nCells) + 3.d0*Uaf(L+2*nCells) )/delta(2)
                    else
                        sx = 0.5d0*( Uaf(njp)          - Uaf(njm) )        /delta(2)
                        sy = 0.5d0*( Uaf(nCells+njp)   - Uaf(nCells+njm))  /delta(2)
                        sz = 0.5d0*( Uaf(2*nCells+njp) - Uaf(2*nCells+njm))/delta(2)
                    endif
                    ! velocity component of current for Vy
                    div_C(LL)%Jxsum(n) = div_C(LL)%Jxsum(n) +  sx*valPHYS(nn,4) 
                    div_C(LL)%Jysum(n) = div_C(LL)%Jysum(n) +  sy*valPHYS(nn,4)
                    div_C(LL)%Jzsum(n) = div_C(LL)%Jzsum(n) +  sz*valPHYS(nn,4) 
                endif
                sx=0d0; sy=0d0; sz=0d0
                if ( nz /=0) then
                    nkm  = L - kdz; nkp  = L + kdz;   
                    nkm2 = L-2*kdz; nkp2 = L+2*kdz;   
                    if     (kkm == 0) then
                        sx = 0.5d0*( -Uaf(nkp2) + 4.d0*Uaf(nkp) - 3.d0*Uaf(L) )/delta(3)
                        sy = 0.5d0*( -Uaf(nkp2+nCells) + 4.d0*Uaf(nkp+nCells) - 3.d0*Uaf(L+nCells) )/delta(3)
                        sz = 0.5d0*( -Uaf(nkp2+2*nCells) + 4.d0*Uaf(nkp+2*nCells) - 3.d0*Uaf(L+2*nCells) )/delta(3)
                    elseif (kkp == 0) then
                        sx = 0.5d0*(  Uaf(nkm2) - 4.d0*Uaf(nkm) + 3.d0*Uaf(L) )/delta(3)
                        sy = 0.5d0*(  Uaf(nkm2+nCells) - 4.d0*Uaf(nkm+nCells) + 3.d0*Uaf(L+nCells) )/delta(3)
                        sz = 0.5d0*(  Uaf(nkm2+2*nCells) - 4.d0*Uaf(nkm+2*nCells) + 3.d0*Uaf(L+2*nCells) )/delta(3)
                    else
                        sx = 0.5d0*( Uaf(nkp)          - Uaf(nkm) )        /delta(3)
                        sy = 0.5d0*( Uaf(nCells+nkp)   - Uaf(nCells+nkm))  /delta(3)
                        sz = 0.5d0*( Uaf(2*nCells+nkp) - Uaf(2*nCells+nkm))/delta(3)
                    endif
                    ! velocity component of current for Vz
                    div_C(LL)%Jxsum(n) = div_C(LL)%Jxsum(n) +  sx*valPHYS(nn,5) 
                    div_C(LL)%Jysum(n) = div_C(LL)%Jysum(n) +  sy*valPHYS(nn,5)
                    div_C(LL)%Jzsum(n) = div_C(LL)%Jzsum(n) +  sz*valPHYS(nn,5)
                endif
            endif
        enddo
    enddo
    Jaf(cel_bndX) = 0.0d0 
    Jaf(cel_bndY) = 0.0d0 
    Jaf(cel_bndZ) = 0.0d0  

!=====================================================================
!== cycle for Jaf  Calculating Divergence Jaf
!=====================================================================
    do LL = 1, sizdom_C
        div_C(LL)%Uf = 0.d0
        m =  dom_C(LL)%siznod_C
        
        do nijm = 1,10
            do n = 1, m
                L = dom_C(LL)%nods(n)

                k= ceiling(REAL(L)/( REAL(sdx*sdy) ) )
                IF (k == 1) THEN
                    nij = L
                ELSE
                    nij = L - (k-1)*sdx*sdy
                ENDIF
                j = ceiling( REAL(nij) / REAL(sdx) )  
                i = nij - (j - 1) * sdx  
            
                nim = L - 1;  njm  = L - sdx;   nkm  = L - kdz;
                nip = L + 1;  njp  = L + sdx;   nkp  = L + kdz; 
                kim = geoPHYS_C(i-1,j,k);   kip = geoPHYS_C(i+1,j,k)
                kjm = geoPHYS_C(i,j-1,k);   kjp = geoPHYS_C(i,j+1,k)
                kkm = geoPHYS_C(i,j,k-1);   kkp = geoPHYS_C(i,j,k+1) 
                sx=0.5d0;sy=0.5d0;sz=0.5d0;
                if (kim == 0) then
                    nim = L; sx=1.d0
                endif
                if (kip == 0)  then
                    nip = L; sx=1.d0
                endif
                if (kjm == 0) then 
                    njm = L; sy=1.d0
                endif
                if (kjp == 0) then 
                    njp = L; sy=1.d0
                endif
                if (kkm == 0) then 
                    nkm = L; sz=1.d0
                endif
                if (kkp == 0) then 
                    nkp = L;  sz=1.d0
                endif
                div_C(LL)%divJ(n) = -0.5d0*( sx*( Jaf(nip) - Jaf(nim) )/delta(1) + &  ! dJx/dx + dJy/dy + dJz/dz
                          sy*( Jaf(njp+  nCells) - Jaf(njm+   nCells) )/delta(2) + &
                          sz*( Jaf(nkp + 2*nCells) - Jaf(nkm+ 2*nCells) )/delta(3) ) 
            enddo  
            CALL sprsBCGstabwr (div_C(LL)%valC, div_C(LL)%irow, div_C(LL)%jcol, div_C(LL)%siz_nods, &
                                div_C(LL)%divJ, div_C(LL)%Uf,  tol_bcg, itmax_bcg, iter)
            div_C(LL)%Uf(div_C(LL)%bndx)=0.d0
            div_C(LL)%Uf(div_C(LL)%bndy)=0.d0
            div_C(LL)%Uf(div_C(LL)%bndz)=0.d0 
            
            do n = 1, dom_C(LL)%siznod_C 
                L = dom_C(LL)%nods(n)
                nn  = L + nCells
                nl  = L + 2*nCells
                k= ceiling(REAL(L)/( REAL(sdx*sdy) ) )
                IF (k == 1) THEN
                    nij = L
                ELSE
                    nij = L - (k-1)*sdx*sdy
                ENDIF
                j = ceiling( REAL(nij) / REAL(sdx) )  
                i = nij - (j - 1) * sdx 

                kim = geoPHYS_Cd(i-1,j,k);   kip = geoPHYS_Cd(i+1,j,k)
                kjm = geoPHYS_Cd(i,j-1,k);   kjp = geoPHYS_Cd(i,j+1,k)
                kkm = geoPHYS_Cd(i,j,k-1);   kkp = geoPHYS_Cd(i,j,k+1) 
                sx=0.5d0;sy=0.5d0;sz=0.5d0;
                if (kim <= 0) then
                    kim = n;sx=1.d0
                endif
                if (kip <= 0) then
                    kip = n;sx=1.d0
                endif
                if (kjm <= 0) then
                    kjm = n;sy=1.d0
                endif
                if (kjp <= 0) then
                    kjp = n;sy=1.d0
                endif
                if (kkm <= 0) then
                    kkm = n;sz=1.d0
                endif
                if (kkp <= 0) then
                    kkp = n;sz=1.d0
                endif
                ! Subtracting the gradient of a scalar field
                Jaf(L)  =  Jaf(L)  - sx*(div_C(LL)%Uf(kip) - div_C(LL)%Uf(kim) )/delta(1)
                Jaf(nn) =  Jaf(nn) - sy*(div_C(LL)%Uf(kjp) - div_C(LL)%Uf(kjm) )/delta(2)
                Jaf(nl) =  Jaf(nl) - sz*(div_C(LL)%Uf(kkp) - div_C(LL)%Uf(kkm) )/delta(3)
            enddo ! end cycle n
        enddo ! end cycle nijm
    enddo ! end cycle LL 
!=====================================================================
!== end cycle for calc divJaf
!=====================================================================
    
    ! adding current  eddy from scalar potential + eddy from speed component
    do LL = 1, sizdom_C
        do n = 1, dom_C(LL)%siznod_C 
            L = dom_C(LL)%nods(n)
            nn  = L + nCells
            nl  = L + 2*nCells
            div_C(LL)%Jxsum(n) = div_C(LL)%Jxsum(n) + Jaf(L)  
            div_C(LL)%Jysum(n) = div_C(LL)%Jysum(n) + Jaf(nn)
            div_C(LL)%Jzsum(n) = div_C(LL)%Jzsum(n) + Jaf(nl) 
        enddo
        div_C(LL)%Jxsum(div_C(LL)%bndx)=0.d0
        div_C(LL)%Jysum(div_C(LL)%bndy)=0.d0
        div_C(LL)%Jzsum(div_C(LL)%bndz)=0.d0
    enddo

    !=====================================================================
    !== cycle for calc divJsum
    !=====================================================================
    do LL = 1, sizdom_C
        do nijm = 1, 20 
            do n = 1, dom_C(LL)%siznod_C 
                L = dom_C(LL)%nods(n)
                k= ceiling(REAL(L)/( REAL(sdx*sdy) ) )
                IF (k == 1) THEN
                    nij = L
                ELSE
                    nij = L - (k-1)*sdx*sdy
                ENDIF
                j = ceiling( REAL(nij) / REAL(sdx) )  
                i = nij - (j - 1) * sdx  
                kim = geoPHYS_Cd(i-1,j,k);   kip = geoPHYS_Cd(i+1,j,k)
                kjm = geoPHYS_Cd(i,j-1,k);   kjp = geoPHYS_Cd(i,j+1,k)
                kkm = geoPHYS_Cd(i,j,k-1);   kkp = geoPHYS_Cd(i,j,k+1) 
                sx=0.5d0;sy=0.5d0;sz=0.5d0;
                if (kim <= 0) then
                    kim = n;sx=1.d0
                endif
                if (kip <= 0) then
                    kip = n;sx=1.d0
                endif
                if (kjm <= 0) then
                    kjm = n;sy=1.d0
                endif
                if (kjp <= 0) then
                    kjp = n;sy=1.d0
                endif
                if (kkm <= 0) then
                    kkm = n;sz=1.d0
                endif
                if (kkp <= 0) then
                    kkp = n;sz=1.d0
                endif
                    ! 0.5   0.25    ! dJx/dx + dJy/dy + dJz/dz
                div_C(LL)%divJ(n) =  -0.5d0*( sx*(div_C(LL)%Jxsum(kip) - div_C(LL)%Jxsum(kim) )/delta(1) + & 
                               sy*( div_C(LL)%Jysum(kjp) - div_C(LL)%Jysum(kjm) )/delta(2) + &
                               sz*( div_C(LL)%Jzsum(kkp) - div_C(LL)%Jzsum(kkm) )/delta(3) ) 
            enddo
            div_C(LL)%Uf(div_C(LL)%bndx)=0.d0
            div_C(LL)%Uf(div_C(LL)%bndy)=0.d0
            div_C(LL)%Uf(div_C(LL)%bndz)=0.d0
            
            CALL sprsBCGstabwr (div_C(LL)%valC, div_C(LL)%irow, div_C(LL)%jcol, div_C(LL)%siz_nods, &
                                div_C(LL)%divJ, div_C(LL)%Uf,  tol_bcg, itmax_bcg, iter)
                                
            do n = 1, dom_C(LL)%siznod_C 
                L = dom_C(LL)%nods(n)
                k= ceiling(REAL(L)/( REAL(sdx*sdy) ) )
                IF (k == 1) THEN
                    nij = L
                ELSE
                    nij = L - (k-1)*sdx*sdy
                ENDIF
                j = ceiling( REAL(nij) / REAL(sdx) )  
                i = nij - (j - 1) * sdx 

                kim = geoPHYS_Cd(i-1,j,k);   kip = geoPHYS_Cd(i+1,j,k)
                kjm = geoPHYS_Cd(i,j-1,k);   kjp = geoPHYS_Cd(i,j+1,k)
                kkm = geoPHYS_Cd(i,j,k-1);   kkp = geoPHYS_Cd(i,j,k+1) 
                sx=0.5d0;sy=0.5d0;sz=0.5d0;
                if (kim <= 0) then
                    kim = n;sx=1.d0
                endif
                if (kip <= 0) then
                    kip = n;sx=1.d0
                endif
                if (kjm <= 0) then
                    kjm = n;sy=1.d0
                endif
                if (kjp <= 0) then
                    kjp = n;sy=1.d0
                endif
                if (kkm <= 0) then
                    kkm = n;sz=1.d0
                endif
                if (kkp <= 0) then
                    kkp = n;sz=1.d0
                endif
                ! Subtracting the gradient of a scalar field
                div_C(LL)%Jxsum(n) = div_C(LL)%Jxsum(n) - sx*(div_C(LL)%Uf(kip) - div_C(LL)%Uf(kim) )/delta(1)
                div_C(LL)%Jysum(n) = div_C(LL)%Jysum(n) - sy*(div_C(LL)%Uf(kjp) - div_C(LL)%Uf(kjm) )/delta(2) 
                div_C(LL)%Jzsum(n) = div_C(LL)%Jzsum(n) - sz*(div_C(LL)%Uf(kkp) - div_C(LL)%Uf(kkm) )/delta(3) 
            enddo 
            div_C(LL)%Jxsum(div_C(LL)%bndx)=0.d0
            div_C(LL)%Jysum(div_C(LL)%bndy)=0.d0
            div_C(LL)%Jzsum(div_C(LL)%bndz)=0.d0
        enddo ! fin nijm for LL
    enddo ! fin LL

    !=====================================================================
    !== end cycle for calc divJsum
    !=====================================================================
    Jaf(cel_bndX) = 0.0d0 
    Jaf(cel_bndY) = 0.0d0 
    Jaf(cel_bndZ) = 0.0d0  

    Uaf(cel_bndx) = 0.0d0    
    Uaf(cel_bndy) = 0.0d0 
    Uaf(cel_bndz) = 0.0d0 

    m = 0
    do LL = 1, sizdom_C
        do n = 1, dom_C(LL)%siznod_C  
            m = m + 1   
            Jxsum(m) =  div_C(LL)%Jxsum(n)
            Jysum(m) =  div_C(LL)%Jysum(n)
            Jzsum(m) =  div_C(LL)%Jzsum(n)
        enddo
    enddo
    
endif ! end domains conductivity

!===========================================================================================
IF (Ntime >= Nprint .and. Ntime /=0 ) THEN  ! output with step dtt skip 1st point
    Nprint = Ntime  + Nout
    Npoint= Npoint + 1                      ! point counter with step dtt
    ios=0
    
    call writeVtk_field (Npoint,  sdx, sdy, sdz, nCells, delta, Uaf, Jaf,Jxsum,Jysum,Jzsum, Bfex, Bfey, Bfez, files)
    IF (numMech /= 0) THEN
        CALL writeVtk_src ( Npoint, numfun,  NcellsX, NcellsY, NcellsZ,  new_nodesX, new_nodesY,  new_nodesZ,  &
                            sdx, sdy, delta, files)  
    endif
ENDIF

IF (Ntime >= Nprint_display) THEN
    Nprint_display = Ntime + Nout_display
    if (mod(Npoint,5) == 0) then
        write(*,'( i0,$ )') Npoint
    else 
        write(*,'( a,$ )') '>'
    endif
ENDIF

Ntime = Ntime + 1   ! this is every point
T = T + DT

IF (T < Time) GOTO 2000
PRINT*,'|'

CALL system_clock(k2);
Tcalc =  REAL(k2-k1)/REAL(k0)
   
write(*,'( "solve complet. Tcalc= ", g10.3)') Tcalc

CONTAINS

SUBROUTINE gen_sparse_C_matrix

real(8) v_row(7), valtmp(7)
real(8) s, sz, sy, sx
logical mask(7)
integer i_row(7),coltmp(7), countC, Ltmp, nim, nip, njm, njp, nkm, nkp, nc, nFix, nFiy, nFiz, num_nztmp
integer, allocatable :: temp(:)
real(8), allocatable :: v_row_add(:)

type espm
    integer im
    DOUBLE PRECISION em
    type (espm),pointer ::prec
end type espm
type(espm),pointer ::  sp_coltmp,next_tmp

sz = 1.d0/(delta(3)*delta(3))
sy = 1.d0/(delta(2)*delta(2))
sx = 1.d0/(delta(1)*delta(1)) 
s = 2.d0*(sx + sy + sz);

nullify(sp_coltmp,next_tmp)

do LL = 1, sizdom_C
    div_C(LL)%nnz = 0;
    countC = 0
    div_C(LL)%bndx = [integer:: ];
    div_C(LL)%bndy = [integer:: ];
    div_C(LL)%bndz = [integer:: ];
    num_nztmp=0
    
    do n = 1, dom_C(LL)%siznod_C    
        L = dom_C(LL)%nods(n)
        k= ceiling(REAL(L)/( REAL(sdx*sdy) ) )
        IF (k == 1) THEN
            nij = L
        ELSE
            nij = L - (k-1)*sdx*sdy
        ENDIF
        j = ceiling( REAL(nij) / REAL(sdx) )  
        i = nij - (j - 1) * sdx
        valtmp = 0.0d0;  coltmp = 0; Ltmp=0; 
    
        nFix=0; nFiy=0; nFiz=0  
        countC = countC + 1 
        nc = geoPHYS_Cd(i,j,k)
        nim = geoPHYS_Cd(i-1,j,k);   nip = geoPHYS_Cd(i+1,j,k)
        njm = geoPHYS_Cd(i,j-1,k);   njp = geoPHYS_Cd(i,j+1,k)
        nkm = geoPHYS_Cd(i,j,k-1);   nkp = geoPHYS_Cd(i,j,k+1)
        i_row(1:7) =   [nim, nip, njm, njp, nkm, nkp, nc ] 
        v_row(1:7) =  [-sx, -sx, -sy, -sy, -sz, -sz,  s ]
        mask = .false.
        mask = i_row > 0
        temp = pack([(i, i=1, 7)], mask)
        
        Ltmp = size(temp)
        coltmp(1:Ltmp) = i_row(temp)
        valtmp(1:Ltmp) = v_row(temp)
        if ( XOR(mask(1), mask(2))  ) nFix = 1
        if ( XOR(mask(3), mask(4))  ) nFiy = 1
        if ( XOR(mask(5), mask(6))  ) nFiz = 1
        k0=0 ! for div
        mm2: do k1=1,Ltmp-1  
            do k2=k1+1, Ltmp
                if  (coltmp(k1) == coltmp(k2) ) then 
                    k0 = coltmp(k2)
                    exit mm2
                endif
            enddo
        enddo mm2
        if (k0 /=0 ) then
            print*,'node Fi2 double', k0, 'i=',i, 'j=',j, 'k=',k
            stop
        endif
        if ( nFix == 1 ) div_C(LL)%bndx = [div_C(LL)%bndx, nc]
        if ( nFiy == 1 ) div_C(LL)%bndy = [div_C(LL)%bndy, nc]
        if ( nFiz == 1 ) div_C(LL)%bndz = [div_C(LL)%bndz, nc]
        call full_sort(coltmp, valtmp, Ltmp, 1,1)

        do m = 1, Ltmp                           ! for div
            if (coltmp(m) <= 0) then
                print*, 'cell=',nn, 'm',m,'coltmp(m)',coltmp(m)
                stop
            endif
            allocate(sp_coltmp)
            sp_coltmp=espm(coltmp(m),valtmp(m),next_tmp) 
            num_nztmp = num_nztmp + 1 
            next_tmp => sp_coltmp
        enddo

        div_C(LL)%irow(countC + 1) = div_C(LL)%irow( countC) + Ltmp  ! for div
        
        deallocate(temp)
        if (allocated (v_row_add)) deallocate ( v_row_add)
    enddo
    
    if (num_nztmp /= 0) then
        div_C(LL)%nnz = num_nztmp
        if (countC /= dom_C(LL)%siznod_C) stop 'countC /= siznod_C'
        if (countC /= div_C(LL)%siz_nods) stop 'countC /= div_C(LL)%siz_nods'
        
        allocate(div_C(LL)%jcol(num_nztmp),source=0)
        allocate(div_C(LL)%valC(num_nztmp),source=0.0d0)
    endif
            
    i= num_nztmp
    DO WHILE (associated(sp_coltmp))
        div_C(LL)%jcol(i) = sp_coltmp%im  
        div_C(LL)%valC(i) = sp_coltmp%em   !  valAu
        i=i-1
        sp_coltmp => sp_coltmp%prec  
    END DO
    nullify(sp_coltmp, next_tmp)
 
PRINT '(a,    a,i9, a,g12.5,a)', 'sparse matrix for div Je:',  &
                 ' non zero=',div_C(LL)%nnz, ' density=', 100.0* REAL(num_nztmp)/REAL(countC)/REAL(countC),'%'

enddo;

end SUBROUTINE gen_sparse_C_matrix

!==================================================================
SUBROUTINE gen_sparse_matrix

integer colX(16), colY(16), colZ(16), colU(40) 
real(8) valX(16), valY(16), valZ(16), valU(40) 

real(8) s, a, b, c, d, e, sz, sy, sx, sxy, sxz, syz, dsx, dsy, dsz
integer countU, nAx, nAy, nAz, Lx, Ly, Lz, Lfi,Lfi2, nFix, nFiy, nFiz,  kFi, n

 !  indices for vect potential
INTEGER  kim, kip,  kjm, kjp, kkm, kkp 
INTEGER  kipkp,kipjp,kjpkp,kimkm,kimjm,kjmkm,kipkm,kipjm,kjpkm, kimkp,kimjp,kjmkp

 !  indices for geoPHYS  (magnetic permeability)
INTEGER       mim, mip, mjm, mjp, mkm, mkp 
real(8)  Z0,  Zim, Zip, Zjm, Zjp, Zkm, Zkp
INTEGER mipkp,mipjp,mjpkp,mimkm,mimjm,mjmkm,mipkm,mipjm,mjpkm, mimkp,mimjp,mjmkp
real(8) Zipkp,Zipjp,Zjpkp,Zimkm,Zimjm,Zjmkm,Zipkm,Zipjm,Zjpkm, Zimkp,Zimjp,Zjmkp

! indices for conductivity
INTEGER  nc, nim, nip,  njm, njp,  nkm, nkp, nim2, nip2,  njm2, njp2,  nkm2, nkp2

integer face_short, model
type espm
    integer im
    DOUBLE PRECISION em
    type (espm),pointer ::prec
end type espm
type(espm),pointer :: sp_colX,next_X, sp_colY,next_Y, sp_colZ,next_Z, sp_colU, next_U 
nullify(sp_colX,next_X,  sp_colY,next_Y,  sp_colZ,next_Z,  sp_colU,next_U) 
if (sizdom_C /= 0) then
    num_bndX=0; num_bndY=0; num_bndZ=0
    cel_bndX = [integer:: ];  cel_bndY = [integer:: ];  cel_bndZ = [integer:: ] 

    num_bndUx=0; num_bndUy=0; num_bndUz=0
    cel_bndUx = [integer:: ];  cel_bndUy = [integer:: ];  cel_bndUz = [integer:: ]
endif

num_nzX=0;  num_nzY=0;  num_nzZ=0;  num_nz=0; num_nzU=0

sz = 1.d0/(delta(3)*delta(3))
sy = 1.d0/(delta(2)*delta(2))
sx = 1.d0/(delta(1)*delta(1)) 
s = 2.d0*(sx + sy + sz);

dsx = 0.5d0/delta(1) ! 1/2dx
dsy = 0.5d0/delta(2) ! 1/2dy
dsz = 0.5d0/delta(3) ! 1/2dz

a = 2.0d0/(dt*delta(1))  
b = 2.0d0/(dt*delta(2))
c = 2.0d0/(dt*delta(3))

IF (siznod_Fe /= 0  ) THEN 
    IF (sizdom_C /= 0  ) THEN  
        if      ( mode == 0) then ! if not contact Fe-AL and not move environ  AL
            face_short = 1 
            model = 1 
        elseif ( mode == 1 ) then ! if not contact Fe-AL and move environ  AL
            face_short = 1 
            model = 0 
        elseif ( mode == 2 ) then ! if contact Fe-AL and move environ  AL
            face_short = 0 
            model = 1 
        endif
    else
        face_short = 1 
    ENDIF
endif

!  X-1:nCells  Y-nCells+1: 2*nCells  Fi-2*nCells+1 : nCellsGlob (2*nCells + nCellsFi)
nn = 0; 
countU = 0 
do k=1,sdz; do j=1,sdy; do i=1,sdx

    n = geoPHYS(i,j,k)   
    nn =  nn + 1    
    valX = 0.0d0;   colX = 0;  Lx=0; 
    valY = 0.0d0;   colY = 0;  Ly=0; 
    valZ = 0.0d0;   colZ = 0;  Lz=0;
    nAx=0; nAy=0; nAz=0  
    valU = 0.0d0;   colU = 0;  Lfi=0; Lfi2=0; 
    nFix=0; nFiy=0; nFiz=0  
    kFi=0;
    IF (geoPHYS_C(i,j,k) /=0) THEN
        nc = geoPHYS_C(i,j,k)
        nim = geoPHYS_C(i-1,j,k);   nip = geoPHYS_C(i+1,j,k)
        njm = geoPHYS_C(i,j-1,k);   njp = geoPHYS_C(i,j+1,k)
        nkm = geoPHYS_C(i,j,k-1);   nkp = geoPHYS_C(i,j,k+1) 
        nim2 = geoPHYS_C(i-2,j,k);  nip2 = geoPHYS_C(i+2,j,k);
        njm2 = geoPHYS_C(i,j-2,k);  njp2 = geoPHYS_C(i,j+2,k);
        nkm2 = geoPHYS_C(i,j,k-2);  nkp2 = geoPHYS_C(i,j,k+2);
        kFi=1;                           
        countU = countU + 1 
    endif
    kim = nn-1;  kjm = nn-sdx;  kkm = nn-kdz;
    kip = nn+1;  kjp = nn+sdx;  kkp = nn+kdz;
    ! 8 corners
    IF (i==1 .or.j==1 .or.k==1 .or.i==sdx .or.j==sdy .or.k==sdz ) THEN
        BOUNDARY_A: BLOCK
            real(8) v_row1(6), v_row2(6) 
            integer i_row(6), base_row(6), fact_row(6), i_loc, j_loc
            logical mask(6)
            integer, allocatable :: temp(:)
            real(8), allocatable :: v_row_add(:)
                    
            base_row = [1, sdx, 1, sdy, 1, sdz] 
            i_row (1:6) =   [kim,     kip,      kjm,      kjp,     kkm,       kkp] 
            v_row1(1:6) =  [sx, sx, sy, sy, sz, sz]
            v_row2(1:6) =  [BND(1,1), BND(1,2), BND(2,1), BND(2,2), BND(3,1), BND(3,2)]
                    
            fact_row = [i, i, j, j, k, k]
            mask = .false.
            mask = (fact_row - base_row ) /= 0
            temp = pack([(i_loc, i_loc=1, 6)], mask)
            Lx = size(temp)
            colX(1:Lx) = i_row(temp)
            valX(1:Lx) = v_row1(temp)
            valX(Lx+1) = sum(v_row1(temp))
            colX(Lx+1) = nn
            allocate( v_row_add(Lx) ,source=0.d0)
            j_loc=0
            do i_loc = 1,6,2 
                if    (mask(i_loc) .and. mask(i_loc+1)) then
                    j_loc=j_loc+1; v_row_add(j_loc) = -1.d0 
                    j_loc=j_loc+1; v_row_add(j_loc) = -1.d0 
                ELSEIF (mask(i_loc) ) then
                    j_loc=j_loc+1; v_row_add(j_loc) = v_row2(j_loc) 
                ELSEIF (mask(i_loc+1)) then
                    j_loc=j_loc+1; v_row_add(j_loc) = v_row2(j_loc) 
                endif
            enddo
            valX(1:j_loc) = valX(1:j_loc)*v_row_add(1:j_loc);
            Lx = Lx + 1
            deallocate(temp)
            if (allocated (v_row_add)) deallocate ( v_row_add)
        END BLOCK BOUNDARY_A

        Ly = Lx; Lz = Lx
        colY = nCells+colX; colZ = 2*nCells+colX; valY = valX; valZ = valX
    else
!----------------------------------- 
        mim = geoPHYS(i-1,j,k);  mip = geoPHYS(i+1,j,k);
        mjm = geoPHYS(i,j-1,k);  mjp = geoPHYS(i,j+1,k);
        mkm = geoPHYS(i,j,k-1);  mkp = geoPHYS(i,j,k+1);
        Z0  = valPHYS(n,1);
        Zim = valPHYS(mim,1); Zip = valPHYS(mip,1); 
        Zjm = valPHYS(mjm,1); Zjp = valPHYS(mjp,1); 
        Zkm = valPHYS(mkm,1); Zkp = valPHYS(mkp,1); 
        if (siznod_Fe /=0 )   then
            mipkp = geoPHYS(i+1,j,k+1);  Zipkp = valPHYS(mipkp,1);
            mipjp = geoPHYS(i+1,j+1,k);  Zipjp = valPHYS(mipjp,1);
            mjpkp = geoPHYS(i,j+1,k+1);  Zjpkp = valPHYS(mjpkp,1);
            mimkm = geoPHYS(i-1,j,k-1);  Zimkm = valPHYS(mimkm,1);
            mimjm = geoPHYS(i-1,j-1,k);  Zimjm = valPHYS(mimjm,1);
            mjmkm = geoPHYS(i,j-1,k-1);  Zjmkm = valPHYS(mjmkm,1);
            mipkm = geoPHYS(i+1,j,k-1);  Zipkm = valPHYS(mipkm,1);
            mipjm = geoPHYS(i+1,j-1,k);  Zipjm = valPHYS(mipjm,1);
            mjpkm = geoPHYS(i,j+1,k-1);  Zjpkm = valPHYS(mjpkm,1);
            mimkp = geoPHYS(i-1,j,k+1);  Zimkp = valPHYS(mimkp,1);
            mimjp = geoPHYS(i-1,j+1,k);  Zimjp = valPHYS(mimjp,1);
            mjmkp = geoPHYS(i,j-1,k+1);  Zjmkp = valPHYS(mjmkp,1);
            Lx=7; 
            colX(1:Lx) = [kim, kip, kjm, kjp, kkm, kkp,     nn]
    valX(1:Lx) = [-0.5d0*(Zim+Zip)*sx, -0.5d0*(Zim+Zip)*sx, -0.5d0*(Zim+Zip)*sy, -0.5d0*(Zim+Zip)*sy, -Zim*sz, -Zip*sz, 3.0d0*(Zim+Zip)*s ]  
            Ly=7; 
            colY(1:Ly) = [kim, kip, kjm, kjp, kkm, kkp,     nn]+ nCells 
    valY(1:Ly) = [-0.5d0*(Zjm+Zjp)*sx, -0.5d0*(Zjm+Zjp)*sx, -0.5d0*(Zjm+Zjp)*sy, -0.5d0*(Zjm+Zjp)*sy, -Zjm*sz, -Zjp*sz, 3.0d0*(Zjm+Zjp)*s ] 
            Lz=7; 
            colZ(1:Lz) = [kim, kip, kjm, kjp, kkm, kkp,     nn] + 2*nCells
    valZ(1:Lz) = [-0.5d0*(Zkm+Zkp)*sx, -0.5d0*(Zkm+Zkp)*sx, -0.5d0*(Zkm+Zkp)*sy, -0.5d0*(Zkm+Zkp)*sy, -Zkm*sz, -Zkp*sz, 3.0d0*(Zkm+Zkp)*s ] 
               ! edges ****************
            if ( lar(Z0, Zipkp) .and. equ(Z0, Zip) .and. equ(Z0, Zkp)  .and. equ(Zjp, Zjm) ) then  ! i-1; k-1 (1) 
                Ly=7
                e = 0.5d0*(Z0 + Zipkp)
                d = (1.5d0*Z0 + 0.5d0*Zipkp)*(sx+sz) + 2.d0*e*sy
                colY(1:Ly) = nCells+[kim,    kkm,     kip,   kkp,    kjm,   kjp,  nn]
                valY(1:Ly) =        [-Z0*sx, -Z0*sz, -e*sx, -e*sz,  -e*sy, -e*sy,  d]
            elseif ( lar(Z0, Zjpkp) .and. equ(Z0, Zjp) .and. equ(Z0, Zkp) .and.  equ(Zim, Zip)  ) then  ! j-1; k-1 (2) 
                Lx=7
                e = 0.5d0*(Z0 + Zjpkp)
                d = (1.5d0*Z0 + 0.5d0*Zjpkp)*(sy+sz) + 2.d0*e*sx
                colX(1:Lx) = [ kjm,   kkm,     kjp,   kkp,    kim,   kip,  nn]
                valX(1:Lx) = [-Z0*sy, -Z0*sz, -e*sy, -e*sz,  -e*sx, -e*sx,  d]
            elseif ( lar(Z0, Zimkp) .and. equ(Z0, Zim)  .and.  equ(Z0, Zkp) .and. equ(Zjm, Zjp)   ) then  ! i+1; k-1 (3) 
                Ly=7
                e = 0.5d0*(Zipkm + Z0)    
                d = (1.5d0*Zipkm + 0.5d0*Z0)*(sx+sz) + 2.d0*e*sy
                colY(1:Ly) = nCells+[ kip,       kkm,      kim,   kkp,    kjm,   kjp,  nn]
                valY(1:Ly) =        [-Zipkm*sx, -Zipkm*sz, -e*sx, -e*sz,  -e*sy, -e*sy,  d]
            elseif ( lar(Z0, Zjmkp) .and. equ(Z0, Zjm) .and. equ(Z0, Zkp) .and. equ(Zim, Zip) ) then  ! j+1; k-1 (4)
                Lx=7
                e = 0.5d0*(Z0 + Zjmkp)
                d = (1.5d0*Z0 + 0.5d0*Zjmkp)*(sy+sz) + 2.d0*e*sx
                colX(1:Lx) = [ kjp,    kkm,     kjm,   kkp,    kim,   kip,  nn]
                valX(1:Lx) = [-Z0*sy, -Z0*sz, -e*sy, -e*sz,  -e*sx, -e*sx,  d]
            elseif (  lar(Z0, Zipjp) .and. equ(Z0, Zip) .and. equ(Z0, Zjp) .and.  equ(Zkp, Zkm)    ) then  ! i-1; j-1 (5)  
                Lz=7
                e = 0.5d0*(Z0 + Zipjp)
                d = (1.5d0*Z0 + 0.5d0*Zipjp)*(sx+sy) + 2.d0*e*sz
                colZ(1:Lz)=2*nCells+[kim,     kjm,     kip,   kjp,    kkm,   kkp,  nn]
                valZ(1:Lz)=         [-Z0*sx, -Z0*sy, -e*sx, -e*sy,  -e*sz, -e*sz,  d]/Zipjp
            elseif (  lar(Z0, Zimjp) .and. equ(Z0, Zim) .and. equ(Z0, Zjp) .and.  equ(Zkp, Zkm)    ) then  ! i+1; j-1 (6)  
                Lz=7
                e = 0.5d0*(Z0 + Zimjp)
                d = (1.5d0*Z0 + 0.5d0*Zimjp)*(sx+sy) + 2.d0*e*sz
                colZ(1:Lz)=2*nCells+[kip,     kjm,     kim,   kjp,    kkm,   kkp,  nn]
                valZ(1:Lz)=         [-Z0*sx, -Z0*sy, -e*sx, -e*sy,  -e*sz, -e*sz,  d]/Zimjp
            elseif (  lar(Z0, Zipjm) .and. equ(Z0, Zip) .and. equ(Z0, Zjm)  .and.  equ(Zkp, Zkm)    ) then ! i-1; j+1 (7) 
                Lz=7
                e = 0.5d0*(Z0 + Zipjm)
                d = (1.5d0*Z0 + 0.5d0*Zipjm)*(sx+sy) + 2.d0*e*sz
                colZ(1:Lz)=2*nCells+[kim,     kjp,     kip,   kjm,    kkm,   kkp,  nn]
                valZ(1:Lz)=         [-Z0*sx, -Z0*sy, -e*sx, -e*sy,  -e*sz, -e*sz,  d]/Zipjm
            elseif ( lar(Z0, Zimjm) .and.  equ(Z0, Zim) .and. equ(Z0, Zjm) .and. equ(Zkp, Zkm)  ) then ! i+1; j+1 (8)  
                Lz=7
                e = 0.5d0*(Z0 + Zimjm)
                d = (1.5d0*Z0 + 0.5d0*Zimjm)*(sx+sy) + 2.d0*e*sz
                colZ(1:Lz)=2*nCells+[kip,     kjp,     kim,   kjm,    kkm,   kkp,  nn]
                valZ(1:Lz)=         [-Z0*sx, -Z0*sy, -e*sx, -e*sy,  -e*sz, -e*sz,  d]/Zimjm
            elseif ( lar(Z0, Zipkm) .and. equ(Z0, Zip)  .and. equ(Z0, Zkm)  .and. equ(Zjp, Zjm) ) then ! i-1; k+1 (9)
                Ly=7
                e = 0.5d0*(Z0 + Zipkm)
                d = (1.5d0*Z0 + 0.5d0*Zipkm)*(sx+sz) + 2.d0*e*sy
                colY(1:Ly) = nCells+[ kim,    kkp,     kip,   kkm,    kjm,   kjp,  nn]
                valY(1:Ly) =        [-Z0*sx, -Z0*sz, -e*sx, -e*sz,  -e*sy, -e*sy,  d]
            elseif ( lar(Z0, Zjpkm) .and. equ(Z0, Zjp) .and. equ(Z0, Zkm) .and. equ(Zim, Zip)    ) then ! j-1; k+1 (10)
                Lx=7
                e = 0.5d0*(Z0 + Zjpkm)
                d = (1.5d0*Z0 + 0.5d0*Zjpkm)*(sy+sz) + 2.d0*e*sx
                colX(1:Lx) = [ kjm,    kkp,     kjp,   kkm,    kim,   kip,  nn]
                valX(1:Lx) = [-Z0*sy, -Z0*sz, -e*sy, -e*sz,  -e*sx, -e*sx,  d]
            elseif (  lar(Z0,Zimkm) .and. equ(Z0, Zim)  .and.  equ(Z0, Zkm) .and. equ(Zjp, Zjm) ) then  ! i+1; k+1 (11)
                Ly=7
                e = 0.5d0*(Z0 + Zimkm)
                d = (1.5d0*Z0 + 0.5d0*Zimkm)*(sx+sz) + 2.d0*e*sy
                colY(1:Ly) = nCells+[ kip,    kkp,     kim,   kkm,    kjm,   kjp,  nn]
                valY(1:Ly) =        [-Z0*sx, -Z0*sz, -e*sx, -e*sz,  -e*sy, -e*sy,  d]
            elseif ( lar(Z0,Zjmkm) .and. equ(Z0, Zjm).and. equ(Z0, Zkm) .and.   equ(Zim, Zip)  ) then  ! j+1; k+1 (12)
                Lx=7
                e = 0.5d0*(Z0 + Zjmkm)
                d = (1.5d0*Z0 + 0.5d0*Zjmkm)*(sy+sz) + 2.d0*e*sx
                colX(1:Lx) = [ kjp,    kkp,     kjm,   kkm,    kim,   kip,  nn]
                valX(1:Lx) = [-Z0*sy, -Z0*sz, -e*sy, -e*sz,  -e*sx, -e*sx,  d]
        
                ! 6 faces  =======================
            elseif ( neq(Z0, Zip) .and. equ(Zjp, Zjm) .and. equ(Zkp, Zkm) ) then ! not   i-1  
                if (face_short == 1) then
                    colY(1:3) =         [kim,   kip,     nn    ] + nCells
                    valY(1:3) = 2d0*dsx*[-Zim, -Zip,  (Zim+Zip)] 
                    Ly = 3
                else
                    Ly = 7
                    colY(1:Ly) =              [kim,                 kip,                 kjm, kjp, kkm, kkp,  nn] + nCells
                    valY(1:Ly)=(Zim+Zip)*[-2.d0*Zim/(Zim+Zip)*sx, -2.d0*Zip/(Zim+Zip)*sx, -sy, -sy, -sz, -sz,  s]
                endif
                Lz = Ly
                colZ = colY + nCells
                valZ = valY
            elseif ( neq(Z0, Zim) .and. equ(Zjp, Zjm) .and. equ(Zkp, Zkm) ) then ! not   i+1 
                if (face_short == 1) then
                    colY(1:3) =         [kim,   kip,     nn    ] + nCells
                    valY(1:3) = 2d0*dsx*[-Zim, -Zip,  (Zim+Zip)]
                    Ly = 3
                else
                    Ly = 7
                    colY(1:Ly) =         [         kim,                   kip,            kjm, kjp, kkm, kkp, nn] + nCells
                    valY(1:Ly)=(Zim+Zip)*[-2.d0*Zim/(Zim+Zip)*sx, -2.d0*Zip/(Zim+Zip)*sx, -sy,-sy,  -sz, -sz,  s]
                endif
                Lz = Ly
                colZ = colY + nCells
                valZ = valY
            elseif ( neq(Z0, Zjp) .and. equ(Zip, Zim) .and. equ(Zkp, Zkm) ) then ! not   j-1  
                if (face_short == 1) then
                    colX(1:3) =         [kjm,   kjp,     nn    ]
                    valX(1:3) = 2d0*dsy*[-Zjm, -Zjp,  (Zjm+Zjp)] 
                    Lx = 3
                else
                    Lx=7
                    colX(1:Lx) =          [kim, kip,  kjm,                 kjp,             kkm,  kkp, nn]
                    valX(1:Lx) =(Zjm+Zjp)*[-sx, -sx, -2.d0*Zjm/(Zjm+Zjp)*sy,-2.d0*Zjp/(Zjm+Zjp)*sy,-sz, -sz,  s]
                endif
                Lz=Lx
                colZ(1:Lz) = colX(1:Lx) + 2*nCells
                valZ(1:Lz) = valX(1:Lx)
            elseif ( neq(Z0, Zjm) .and. equ(Zip, Zim) .and. equ(Zkp, Zkm) ) then ! not   j+1  
                if (face_short == 1) then
                    colX(1:3) =         [kjm,   kjp,     nn    ]
                    valX(1:3) = 2d0*dsy*[-Zjm, -Zjp,  (Zjm+Zjp)]
                    Lx = 3
                else
                    Lx=7
                    colX(1:Lx) =          [kim, kip,  kjm,                 kjp,             kkm,  kkp, nn]
                    valX(1:Lx) =(Zjm+Zjp)*[-sx, -sx, -2.d0*Zjm/(Zjm+Zjp)*sy,-2.d0*Zjp/(Zjm+Zjp)*sy,-sz, -sz,  s]
                endif
                Lz=Lx
                colZ(1:Lz) = colX(1:Lx) + 2*nCells
                valZ(1:Lz) = valX(1:Lx)
            elseif ( neq(Z0, Zkp) .and. equ(Zip, Zim) .and. equ(Zjp, Zjm)) then ! not   k-1 
                if (face_short == 1) then
                    colX(1:3) =         [kkm,   kkp,     nn    ] 
                    valX(1:3) = 2d0*dsz*[-Zkm, -Zkp,  (Zkm+Zkp)] 
                    Lx = 3
                else
                    Lx=7
                    colX(1:Lx) =          [kim, kip, kjm,kjp,       kkm,                  kkp,              nn]
                    valX(1:Lx)=(Zkm+Zkp)*[-sx,-sx,  -sy,- sy,-2.d0*Zkm/(Zkm+Zkp)*sz,-2.d0*Zkp/(Zkm+Zkp)*sz, s]
                endif
                Ly=Lx
                colY(1:Ly) = colX(1:Lx) + nCells
                valY(1:Ly) = valX(1:Lx)
            elseif ( neq(Z0, Zkm) .and. equ(Zip, Zim) .and. equ(Zjp, Zjm)) then ! not   k+1 
                if (face_short == 1) then
                    colX(1:3) =         [kkm,   kkp,     nn    ]
                    valX(1:3) = 2d0*dsz*[-Zkm, -Zkp,  (Zkm+Zkp)] 
                    Lx = 3
                else
                    Lx=7
                    colX(1:Lx) =          [kim, kip, kjm,kjp,       kkm,                  kkp,              nn]
                    valX(1:Lx)=(Zkm+Zkp)*[-sx,-sx,  -sy,- sy,-2.d0*Zkm/(Zkm+Zkp)*sz,-2.d0*Zkp/(Zkm+Zkp)*sz, s]
                endif
                Ly=Lx
                colY(1:Ly) = colX(1:Lx) + nCells
                valY(1:Ly) = valX(1:Lx)
            else       ! end boundary
                IF (sizdom_C == 0  ) THEN ! =========== only for mag domains 
                    sxy = 0.25d0/(delta(1)*delta(2))  
                    sxz = 0.25d0/(delta(1)*delta(3)) 
                    syz = 0.25d0/(delta(2)*delta(3))  
                    kipkp = nn + 1 + kdz;   kimkm = nn - 1 - kdz;   kipkm = nn + 1 - kdz;   kimkp = nn - 1 + kdz;
                    kipjp = nn + 1 + sdx;   kimjm = nn - 1 - sdx;   kipjm = nn + 1 - sdx;   kimjp = nn - 1 + sdx;  
                    kjpkp = nn + sdx + kdz; kjmkm = nn - sdx - kdz; kjpkm = nn + sdx - kdz; kjmkp = nn - sdx + kdz;
            
                    Lx=13
                    d= 2.d0*sy + 2.d0*sz;
                    valX(1:Lx)   = [d,         &   !1
                    -sy,   -sy,   -sz,    -sz,       &   !2
                    +sxy,  -sxy,  -sxy,   +sxy,      &   !3
                    +sxz,  -sxz,  -sxz,   +sxz    ]       !
                    colX(1:Lx)= [nn,                                                & !1 dAx(1/dy+1/dz)
                    kjp,           kjm,           kkp,           kkm,               & !2
                    kipjp+nCells,   kimjp+nCells,   kipjm+nCells,   kimjm+nCells,   & !3 Ay/dxdy
                    kipkp+2*nCells, kjmkp+2*nCells, kipkm+2*nCells, kimkm+2*nCells]   !4 Az/dxdz
                    
                    Ly=13
                    colY(1:Ly) = [nn +nCells,                                       & !1 (dAy/dz+dAy/dx) 
                    kip +nCells,    kim+nCells,     kkp+nCells,   kkm+nCells,       & !2 
                    kipjp,          kipjm,          kimjp,          kimjm,          & !3 dAx/dxdy
                    kjpkp+2*nCells, kjmkp+2*nCells, kjpkm+2*NCells, kjmkm+2*NCells]   !4 dAz/dzdy
                    d=2.d0*sx + 2.d0*sz 
                    valY(1:Ly)   = [d,           &  !1
                    -sx,  -sx,  -sz,   -sz,     &  !2
                    +sxy, -sxy, -sxy,  +sxy,    &  !3 1/dxdy
                    +syz, -syz, -syz,  +syz ]      !4 1/dzdy
                    
                    Lz=13
                    colZ(1:Lz) = [nn +2*nCells,                                &  !1 dAz/dx+dAz/dy) 
                    kip+2*NCells, kim+2*nCells, kjp+2*nCells, kjm+2*nCells,    &  !2
                    kipkp,        kipkm,        kimkp,          kimkm,         &  !3 dAx/dxdz
                    kjpkp+nCells, kjpkm+nCells, kjmkp+nCells,   kjmkm+nCells ]    !4 dAy/dydz 
                    d= 2.d0*sx + 2.d0*sy 
                    valZ(1:Lz) = [d,      &  !1
                    -sx,   -sx,   -sy,  -sy,    &  !2
                    +sxz,  -sxz,  -sxz, +sxz,   &  !3
                    +syz,  -syz,  -syz, +syz  ]    !4
        
                elseif (sizdom_C /= 0  ) then ! mag and conductivity
                    if ( model == 1 ) then
                        kimkm = nn - 1 - kdz;  
                        kimjm = nn - 1 - sdx;  
                        kjmkm = nn - sdx - kdz; 
                        sxy = 0.250d0/(delta(1)*delta(2))
                        sxz = 0.250d0/(delta(1)*delta(3))
                        syz = 0.250d0/(delta(2)*delta(3))
                        Lx=13;  
                        colX(1:Lx) = [  nn,             nn   +nCells, nn+ 2*nCells,       &  !1
                                        kjp,            kjm,          kkp,         kkm,   &  !2
                                        kimjm  +nCells, kjm  +nCells, kim  +nCells,       &  !3
                                        kimkm+2*nCells, kim+2*nCells, kkm+2*nCells]          !4
                        d= 2.d0*sy + 2.d0*sz; 
                        valX(1:Lx)   = [d,      sxy,          sxz,                    &   !1
                                        -sy,         -sy,           -sz,          -sz,     &   !2
                                        +sxy,        -sxy,          -sxy,                  &   !3
                                        +sxz,        -sxz,          -sxz            ]          !4
                        Ly=13;
                        colY(1:Ly) = [  nn   +nCells,  nn,            nn+2*nCells,            &  !1
                                        kip   +nCells, kim  +nCells, kkp+nCells, kkm+nCells,  &  !2
                                        kimjm,         kjm,          kim,                      &  !3
                                        kjmkm+2*nCells, kjm+2*nCells, kkm+2*nCells]               !4
                        d= 2.d0*sx + 2.d0*sz;
                        valY(1:Ly)   = [d,      sxy,          syz,                  &  !1
                                        -sx,         -sx,           -sz,       -sz,      &  !2
                                        +sxy,        -sxy,           -sxy,               &  !3
                                        +syz,        -syz,          -syz            ]       !4
                        Lz=13; 
                        colZ(1:Lz) = [  nn +2*nCells,  nn,           nn+nCells,                 &  !1
                                        kip+2*nCells,  kim+2*nCells, kjp+2*nCells, kjm+2*nCells, &  !2
                                        kjmkm +nCells, kkm+nCells,   kjm+nCells,                &  !3
                                        kimkm,         kkm,          kim]                         !4                                     
                        d= 2.d0*sx + 2.d0*sy;
                        valZ(1:Lz)   = [d,      sxz,           syz,                   &  !1
                                        -sx,         -sx,           -sy,       -sy,         &  !2
                                        +syz,        -syz,          -syz,                    &  !3
                                        +sxz,        -sxz,          -sxz            ]          !4
                    elseif ( model == 0 ) then 
                        Lx=7; 
                        colX(1:Lx) = [kim, kip, kjm, kjp, kkm, kkp,     nn]
                        valX(1:Lx) =  0.5d0*[-sx,-sx, -sy, -sy, -sz, -sz, 2.d0*(sx + sy + sz)  ]
                        Ly=Lx; Lz=Lx 
                        colY = nCells+colX; colZ = 2*nCells+colX; valY = valX; valZ = valX
                    else
                        stop 'not model'
                    endif 
                endif
            endif 
        
        else     !  if siznod_Fe ==0  
            Lx=7; ! for Ax
            colX(1:Lx) = [kim, kip, kjm, kjp, kkm, kkp,     nn]
            valX(1:Lx) =  [-sx,-sx, -sy, -sy, -sz, -sz, 2.d0*(sx + sy + sz)  ]
            Ly=Lx; Lz=Lx 
            colY = nCells+colX; colZ = 2*nCells+colX; valY = valX; valZ = valX
        endif 
    
        if (kFi /= 0) then 
            if ( model == 0 ) then
                valX(7) = valX(7) + 2.d0*valPHYS(n,2)/dt    ! 2*sigms*mu0/dt
                valY = valX; valZ = valX
            else
                valX(1) = valX(1) + 2.d0*valPHYS(n,2)/dt    ! 2*sigms*mu0/dt
                valY(1) = valY(1) + 2.d0*valPHYS(n,2)/dt
                valZ(1) = valZ(1) + 2.d0*valPHYS(n,2)/dt
            endif
            if     (nip == 0) then
                colX(Lx+1:Lx+3 ) = [nc, nim, nim2];  valX(Lx+1:Lx+3 ) = [-3.d0, +4.d0, -1.d0]*valPHYS(n,2)*dsx
                Lx = Lx + 3;  nAx  = 1
            elseif (nim == 0) then
                colX(Lx+1:Lx+3 ) = [nc, nip, nip2];  valX(Lx+1:Lx+3 ) = [3.d0, -4.d0, +1.d0]*valPHYS(n,2)*dsx
                Lx = Lx + 3;  nAx  = 1
            else
                colX(Lx+1:Lx+2 ) = [ nip, nim];  valX(Lx+1:Lx+2 ) = [-1.d0, 1.d0]*valPHYS(n,2)*dsx
                Lx = Lx + 2
            endif
            if     (njp == 0) then
                colY(Ly+1:Ly+3 ) = [nc, njm, njm2];  valY(Ly+1:Ly+3 ) = [-3.d0, +4.d0, -1.d0]*valPHYS(n,2)*dsy
                Ly = Ly + 3;  nAy = 1
            elseif (njm == 0) then
                colY(Ly+1:Ly+3 ) = [nc, njp, njp2];  valY(Ly+1:Ly+3 ) = [3.d0, -4.d0, 1.d0]*valPHYS(n,2)*dsy
                Ly = Ly + 3;  nAy = 1
            else 
                colY(Ly+1:Ly+2 ) = [ njp, njm];  valY(Ly+1:Ly+2 ) = [-1.d0, 1.d0]*valPHYS(n,2)*dsy
                Ly = Ly + 2
            endif
            if     (nkp == 0) then
                colZ(Lz+1:Lz+3 ) = [nc, nkm, nkm2];  valZ(Lz+1:Lz+3 ) = [-3.d0, +4.d0, -1.d0]*valPHYS(n,2)*dsz
                Lz = Lz + 3;  nAz = 1
            elseif (nkm == 0) then
                colZ(Lz+1:Lz+3 ) = [nc, nkp, nkp2];  valZ(Lz+1:Lz+3 ) = [3.d0, -4.d0, 1.d0]*valPHYS(n,2)*dsz
                Lz = Lz + 3;  nAz = 1
            else
                colZ(Lz+1:Lz+2 ) = [ nkp, nkm];  valZ(Lz+1:Lz+2 ) = [-1.d0, 1.d0]*valPHYS(n,2)*dsz
                Lz = Lz + 2
            endif
        endif
    endif   ! end form row for A

    !--------------------------X------------------------------------
    call full_sort(colX, valX, Lx, 1,1)
    do m = 1, Lx
        if (colX(m) <= 0) then
            print*, 'maybe the conductivity region by X is too narrow'
            print*, 'cell=',nn, 'm',m,'colX(m)',colX(m),' coord=',i,j,k
            stop
        endif
        allocate(sp_colX)
        sp_colX=espm(colX(m),valX(m),next_X) 
        num_nzX = num_nzX + 1
        next_X => sp_colX
    enddo
    irow(nn+1) = irow(nn) + Lx  !   
    !---------------------------------Y----------------------
    call full_sort(colY, valY, Ly, 1,1)
    do m = 1, Ly
        if (colY(m) <= 0) then
            print*, 'maybe the conductivity region by Y is too narrow'
            print*, 'cell=',nn, 'm',m,'colY(m)',colY(m),' coord=',i,j,k
            stop
        endif
        allocate(sp_colY)
        sp_colY=espm(colY(m),valY(m),next_Y) 
        num_nzY = num_nzY + 1
        next_Y => sp_colY
    enddo
    irow( nCells + nn+1) = irow( nCells + nn) + Ly  
    !-------------------------------Z---------------------------------------
    call full_sort(colZ, valZ, Lz, 1,1)
    do m = 1, Lz
        if (colZ(m) <= 0) then
            print*, 'maybe the conductivity region by Z is too narrow'
            print*, 'cell=',nn, 'm',m,'colZ(m)',colZ(m),' coord=',i,j,k
            stop
        endif
        allocate(sp_colZ)
        sp_colZ=espm(colZ(m),valZ(m),next_Z) 
        num_nzZ = num_nzZ + 1
        next_Z => sp_colZ
    enddo
    irow( 2*nCells + nn+1) = irow( 2*nCells + nn) + Lz  !  
    
    if ( nAx==1 )  cel_bndX = [cel_bndX, nn]
    if ( nAy==1 )  cel_bndY = [cel_bndY, (nn + nCells) ]
    if ( nAz==1 )  cel_bndZ = [cel_bndZ, (nn + 2*nCells) ]
!================================================
!===========END vect cells A
!=================================================
    if (kFi /=0) then  
        nc = geoPHYS_C(i,j,k)
        if (nc /= 0 ) then
            BCGSTAB:BLOCK
                real(8) v_row(7), v_row2(3) 
                integer i_row(7), i_row2(3),i,j, k
                logical mask(7)
                integer, allocatable :: temp(:)
                real(8), allocatable :: v_row_add(:)
                
                i_row(1:7) =   [nim, nip, njm, njp, nkm, nkp, nc ] 
                v_row(1:7) =  [-sx, -sx, -sy, -sy, -sz, -sz,  s ]
                i_row2(1:3) =   [nn, nCells+nn, 2*nCells + nn ] 
                v_row2(1:3) =   [a,  b,          c ]
                mask = .false.
                mask = i_row > 0
                temp = pack([(i, i=1, 7)], mask)
                Lfi2 = size(temp)
                ! 8 corners, 12 edges, 6 faces   
                if ( Lfi2 /= 7) then 
                    Lfi = 7
                    colU(1:Lfi2) = i_row(temp) 
                    valU(1:Lfi2) = v_row(temp)
                    allocate( v_row_add(Lfi2) ,source=0.d0)
                    j=0
                    do i = 1,6,2
                        if    (mask(i) .and. mask(i+1)) then
                            j=j+1; v_row_add(j) = 1.d0
                            j=j+1; v_row_add(j) = 1.d0
                        ELSEIF (mask(i) ) then
                            j=j+1; v_row_add(j) = 2.d0
                        ELSEIF (mask(i+1)) then
                            j=j+1; v_row_add(j) = 2.d0
                        endif
                    enddo
                    valU(1:Lfi2-1) = v_row(temp(1:Lfi2-1))*v_row_add(1:Lfi2-1);
                    j=0; k=0
                    do i = 1,6,2
                        if    (mask(i) .and. mask(i+1)) then
                            j = j+1
                        ELSEIF (mask(i) ) then
                            j=j+1; k = k+1
                            colU(Lfi2+k) = i_row2(j)
                            valU(Lfi2+k) = v_row2(j) 
                        ELSEIF (mask(i+1)) then
                            j=j+1;  k = k+1
                            colU(Lfi2+k) = i_row2(j)
                            valU(Lfi2+k) =-v_row2(j) 
                        endif
                    enddo
                    if ( XOR(mask(1), mask(2))  ) nFix = 1
                    if ( XOR(mask(3), mask(4))  ) nFiy = 1
                    if ( XOR(mask(5), mask(6))  ) nFiz = 1
                elseif ( Lfi2 == 7) then  
                    colU(1:7) = [nim,  nip, njm, njp, nkm, nkp, nc]
                    valU(1:7) = [-sx, -sx, -sy, -sy, -sz, -sz,  s ]
                    colU(8:13) =        [ kip,  kim,  nCells+kjp,  nCells+kjm,  2*nCells+kkp,  2*nCells+kkm  ]
                    valU(8:13) = 0.25d0*[-a,    a,     -b,                   b,      -c,              c]
                    Lfi = 13
                endif
                deallocate(temp)
                if (allocated (v_row_add)) deallocate ( v_row_add)
            END BLOCK BCGSTAB
              
            k0=0 ! for base
            mm: do k1=1,Lfi-1  
                do k2=k1+1, Lfi
                    if  (colU(k1) == colU (k2) ) then 
                        k0 = colU(k2)
                        exit mm
                    endif
                enddo
            enddo mm
            if (k0 /=0 ) then
                print*,'node Fi double', k0, 'i=',i, 'j=',j, 'k=',k
                stop
            endif
            if ( nFix == 1 ) cel_bndUx = [cel_bndUx, nc]
            if ( nFiy == 1 ) cel_bndUy = [cel_bndUy, nc]
            if ( nFiz == 1 ) cel_bndUz = [cel_bndUz, nc]
            call full_sort(colU, valU, Lfi, 1,1)
            do m = 1, Lfi                         ! for base U
                if (colU(m) <= 0) then
                    print*, 'cell=',nn, 'm',m,'colU(m)',colU(m)
                    stop
                endif
                allocate(sp_colU)
                sp_colU=espm(colU(m),valU(m),next_U) 
                num_nzU = num_nzU + 1 
                next_U => sp_colU
            enddo
            irow( 3*nCells + countU + 1) = irow( 3*nCells + countU) + Lfi ! base U
            
        endif   !geoPHYS_C /=0
    endif       ! kFi/=0  
enddo;enddo;enddo   ! i,j,k 

num_nz = num_nzX + num_nzY + num_nzZ + num_nzU
num_bndX = size(cel_bndX );  num_bndY = size(cel_bndY );  num_bndZ = size(cel_bndZ );  

irow(nCells+1) = num_nzX + 1
do i=nCells+2, 2*nCells
    irow(i) = irow(i) + num_nzX
enddo
m = num_nzX + num_nzY 
irow(2*nCells+1) = m + 1
do i=2*nCells+2, 3*nCells
    irow(i) = irow(i) + m 
enddo
m =  num_nzX + num_nzY + num_nzZ
irow(3*nCells+1) = m + 1
do i=3*nCells+2, nCellsGlob+1
    irow(i) = irow(i) + m
enddo

allocate(jcol(num_nz),source=0)
allocate(valA(num_nz),source=0.0d0)

if (size(jcol) < irow( 3*nCells + countU + 1)-1 ) then
    print*,'base ', i,j,k, countU, ' irow=', irow( 3*nCells + countU + 1) , ' size jcol=', size(jcol)
    stop
endif

i=num_nzX 
DO WHILE (associated(sp_colX))
    jcol(i) = sp_colX%im  !
    valA(i) = sp_colX%em   !
    i=i-1
    sp_colX => sp_colX%prec
END DO
nullify(sp_colX, next_X)

i=  num_nzX + num_nzY 
DO WHILE (associated(sp_colY))
    jcol(i) = sp_colY%im  
    valA(i) = sp_colY%em   !
    i=i-1
    sp_colY => sp_colY%prec
END DO
nullify(sp_colY, next_Y)

i=  num_nzX + num_nzY + num_nzZ
DO WHILE (associated(sp_colZ))
    jcol(i) = sp_colZ%im  
    valA(i) = sp_colZ%em   !
    i=i-1
    sp_colZ => sp_colZ%prec
END DO
nullify(sp_colZ, next_Z)
! tau = 18*0.005 = 0.09   2*tau*50 = 9 m/c
i=num_nz
DO WHILE (associated(sp_colU))
    jcol(i) = sp_colU%im  
    valA(i) = sp_colU%em   !
    i=i-1
    sp_colU => sp_colU%prec
END DO
nullify(sp_colU, next_U)

PRINT '(a,i3,a,i3,a,i3,     a,i9, a,g12.5,a)', 'base sparse matrix on grid (', sdx,' x ',sdy,' x ',sdz, &
                 ' ), Non zero= ',num_nz, ' Density=', 100.0* REAL(num_nz)/REAL(nCells)/REAL(nCells),'%'
                 

END SUBROUTINE gen_sparse_matrix
!==============================================
!================================================
function equ (a, b)
logical equ
real(8) a, b
    equ = abs(a - b) < 1.0d-7
end function

function neq (a, b)
logical neq
real(8) a, b
    neq = abs(a - b) > 1.0d-7
end function

function lar (a, b)
logical lar
real(8) a, b
    lar = a > b + 1.d-7
end function

function les (a, b)
logical les
real(8) a, b
    les = a < b - 1.d-7
end function


SUBROUTINE update_rhs
    ! 8 corners
    if     ( nim == 0 .and. njm == 0 .and. nkm == 0  ) then ! not i-1 j-1 k-1 
        if (Nx /=0 .or. Ny /= 0 .or. Nz /= 0 ) then
            colU(Lfi+1:Lfi+3) =      [    L,                      L+nCells,                     L+2*nCells ]
            valU(Lfi+1:Lfi+3) = -3d0*[ Vex*sx+Vey*sxy+Vez*sxz,  Vex*sxy+Vey*sy+Vez*syz,  Vex*sxz+Vey*syz+Vez*sz] 
            Lfi = Lfi + 3
        endif
        if (Nx /=0) then
            colU(Lfi+1:Lfi+6) =     [kip,       kip2, nCells+kip, nCells+kip2, 2*nCells+kip, 2*nCells+kip2 ] 
            valU(Lfi+1:Lfi+6) = -Vex*[ -4d0*sx,   sx,   -4d0*sxy,    sxy,        -4d0*sxz,      sxz    ]
            Lfi = Lfi + 6
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+6) =     [kjp,       kjp2, nCells+kjp, nCells+kjp2, 2*nCells+kjp, 2*nCells+kjp2 ] 
            valU(Lfi+1:Lfi+6) = -Vey*[ -4d0*sxy,   sxy,   -4d0*sy,    sy,        -4d0*syz,      syz    ]
            Lfi = Lfi + 6
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+6) =     [kkp,       kkp2, nCells+kkp, nCells+kkp2, 2*nCells+kkp, 2*nCells+kkp2 ] 
            valU(Lfi+1:Lfi+6) = -Vez*[ -4d0*sxz,   sxz,   -4d0*syz,    syz,         -4d0*sz,      sz    ]
            Lfi = Lfi + 6
        endif
    elseif ( nip == 0 .and. njm == 0 .and. nkm == 0  ) then ! not i+1 j-1 k-1 
        if (Nx /=0 .or. Ny /= 0 .or. Nz /= 0 ) then
            colU(Lfi+1:Lfi+3) =      [    L,                      L+nCells,                     L+2*nCells ]
            valU(Lfi+1:Lfi+3) = -3d0*[ Vex*sx-Vey*sxy-Vez*sxz,  -Vex*sxy+Vey*sy+Vez*syz,  -Vex*sxz+Vey*syz+Vez*sz] 
            Lfi = Lfi + 3
        endif
        if (Nx /=0) then
            colU(Lfi+1:Lfi+6) =     [kim,       kim2, nCells+kim, nCells+kim2, 2*nCells+kim, 2*nCells+kim2 ] 
            valU(Lfi+1:Lfi+6) = -Vex*[ -4d0*sx,   sx,   4d0*sxy,    -sxy,        4d0*sxz,      -sxz    ]
            Lfi = Lfi + 6
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+6) =     [kjp,       kjp2, nCells+kjp, nCells+kjp2, 2*nCells+kjp, 2*nCells+kjp2 ] 
            valU(Lfi+1:Lfi+6) = -Vey*[ 4d0*sxy,   -sxy,   -4d0*sy,    sy,        -4d0*syz,      syz    ]
            Lfi = Lfi + 6
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+6) =     [kkp,       kkp2, nCells+kkp, nCells+kkp2, 2*nCells+kkp, 2*nCells+kkp2 ] 
            valU(Lfi+1:Lfi+6) = -Vez*[ 4d0*sxz,   -sxz,   -4d0*syz,    syz,         -4d0*sz,      sz    ]
            Lfi = Lfi + 6
        endif
    elseif ( nim == 0 .and. njp == 0 .and. nkm == 0 ) then ! not i-1 j+1 k-1 
        if (Nx /=0 .or. Ny /= 0 .or. Nz /= 0 ) then
            colU(Lfi+1:Lfi+3) =      [    L,                      L+nCells,                     L+2*nCells ]
            valU(Lfi+1:Lfi+3) = -3d0*[ Vex*sx-Vey*sxy+Vez*sxz,  -Vex*sxy+Vey*sy-Vez*syz,  Vex*sxz-Vey*syz+Vez*sz] 
            Lfi = Lfi + 3
        endif
        if (Nx /=0) then
            colU(Lfi+1:Lfi+6) =     [kip,       kip2, nCells+kip, nCells+kip2, 2*nCells+kip, 2*nCells+kip2 ] 
            valU(Lfi+1:Lfi+6) = -Vex*[ -4d0*sx,   sx,   4d0*sxy,    -sxy,        -4d0*sxz,      sxz    ]
            Lfi = Lfi + 6
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+6) =     [kjm,       kjm2, nCells+kjm, nCells+kjm2, 2*nCells+kjm, 2*nCells+kjm2 ] 
            valU(Lfi+1:Lfi+6) = -Vey*[ 4d0*sxy,   -sxy,   -4d0*sy,    sy,        4d0*syz,      -syz    ]
            Lfi = Lfi + 6
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+6) =     [kkp,       kkp2, nCells+kkp, nCells+kkp2, 2*nCells+kkp, 2*nCells+kkp2 ] 
            valU(Lfi+1:Lfi+6) = -Vez*[ -4d0*sxz,   sxz,   4d0*syz,    -syz,         -4d0*sz,      sz    ]
            Lfi = Lfi + 6
        endif
    elseif ( nip == 0 .and. njp == 0 .and. nkm == 0 ) then ! not i+1 j+1 k-1 
        if (Nx /=0 .or. Ny /= 0 .or. Nz /= 0 ) then
            colU(Lfi+1:Lfi+3) =      [    L,                      L+nCells,                     L+2*nCells ]
            valU(Lfi+1:Lfi+3) = -3d0*[ Vex*sx+Vey*sxy-Vez*sxz,  Vex*sxy+Vey*sy-Vez*syz,  -Vex*sxz-Vey*syz+Vez*sz] 
            Lfi = Lfi + 3
        endif
        if (Nx /=0) then
            colU(Lfi+1:Lfi+6) =     [kim,       kim2, nCells+kim, nCells+kim2, 2*nCells+kim, 2*nCells+kim2 ] 
            valU(Lfi+1:Lfi+6) = -Vex*[ -4d0*sx,   sx,   -4d0*sxy,    sxy,        4d0*sxz,      -sxz    ]
            Lfi = Lfi + 6
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+6) =     [kjm,       kjm2, nCells+kjm, nCells+kjm2, 2*nCells+kjm, 2*nCells+kjm2 ] 
            valU(Lfi+1:Lfi+6) = -Vey*[ -4d0*sxy,   sxy,   -4d0*sy,    sy,        4d0*syz,      -syz    ]
            Lfi = Lfi + 6
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+6) =     [kkp,       kkp2, nCells+kkp, nCells+kkp2, 2*nCells+kkp, 2*nCells+kkp2 ] 
            valU(Lfi+1:Lfi+6) = -Vez*[ 4d0*sxz,   -sxz,   4d0*syz,    -syz,         -4d0*sz,      sz    ]
            Lfi = Lfi + 6
        endif
    elseif ( nim == 0 .and. njm == 0 .and. nkp == 0 ) then ! not i-1 j-1 k+1 
        if (Nx /=0 .or. Ny /= 0 .or. Nz /= 0 ) then
            colU(Lfi+1:Lfi+3) =      [    L,                      L+nCells,                     L+2*nCells ]
            valU(Lfi+1:Lfi+3) = -3d0*[ Vex*sx+Vey*sxy-Vez*sxz,  Vex*sxy+Vey*sy-Vez*syz,  -Vex*sxz-Vey*syz+Vez*sz] 
            Lfi = Lfi + 3
        endif
        if (Nx /=0) then
            colU(Lfi+1:Lfi+6) =     [kip,       kip2, nCells+kip, nCells+kip2, 2*nCells+kip, 2*nCells+kip2 ] 
            valU(Lfi+1:Lfi+6) = -Vex*[ -4d0*sx,   sx,   -4d0*sxy,    sxy,        4d0*sxz,      -sxz    ]
            Lfi = Lfi + 6
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+6) =     [kjp,       kjp2, nCells+kjp, nCells+kjp2, 2*nCells+kjp, 2*nCells+kjp2 ] 
            valU(Lfi+1:Lfi+6) = -Vey*[ -4d0*sxy,   sxy,   -4d0*sy,    sy,        4d0*syz,      -syz    ]
            Lfi = Lfi + 6
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+6) =     [kkm,       kkm2, nCells+kkm, nCells+kkm2, 2*nCells+kkm, 2*nCells+kkm2 ] 
            valU(Lfi+1:Lfi+6) = -Vez*[ 4d0*sxz,   -sxz,   4d0*syz,    -syz,         -4d0*sz,      sz    ]
            Lfi = Lfi + 6
        endif
    elseif ( nip == 0 .and. njm == 0 .and. nkp == 0 ) then ! not i+1 j-1 k+1 
        if (Nx /=0 .or. Ny /= 0 .or. Nz /= 0 ) then
            colU(Lfi+1:Lfi+3) =      [    L,                      L+nCells,                     L+2*nCells ]
            valU(Lfi+1:Lfi+3) = -3d0*[ Vex*sx-Vey*sxy+Vez*sxz,  -Vex*sxy+Vey*sy-Vez*syz,  Vex*sxz-Vey*syz+Vez*sz] 
            Lfi = Lfi + 3
        endif
        if (Nx /=0) then
            colU(Lfi+1:Lfi+6) =     [kim,       kim2, nCells+kim, nCells+kim2, 2*nCells+kim, 2*nCells+kim2 ] 
            valU(Lfi+1:Lfi+6) = -Vex*[ -4d0*sx,   sx,   4d0*sxy,    -sxy,        -4d0*sxz,      sxz    ]
            Lfi = Lfi + 6
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+6) =     [kjp,       kjm2, nCells+kjp, nCells+kjp2, 2*nCells+kjp, 2*nCells+kjp2 ] 
            valU(Lfi+1:Lfi+6) = -Vey*[ 4d0*sxy,   -sxy,   -4d0*sy,    sy,        4d0*syz,      -syz    ]
            Lfi = Lfi + 6
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+6) =     [kkm,       kkm2, nCells+kkm, nCells+kkm2, 2*nCells+kkm, 2*nCells+kkm2 ] 
            valU(Lfi+1:Lfi+6) = -Vez*[ -4d0*sxz,   sxz,   4d0*syz,    -syz,         -4d0*sz,      sz    ]
            Lfi = Lfi + 6
        endif
    elseif ( nim == 0 .and. njp == 0 .and. nkp == 0 ) then ! not i-1 j+1 k+1 
        if (Nx /=0 .or. Ny /= 0 .or. Nz /= 0 ) then
            colU(Lfi+1:Lfi+3) =      [    L,                      L+nCells,                     L+2*nCells ]
            valU(Lfi+1:Lfi+3) = -3d0*[ Vex*sx-Vey*sxy-Vez*sxz,  -Vex*sxy+Vey*sy+Vez*syz,  -Vex*sxz+Vey*syz+Vez*sz] 
            Lfi = Lfi + 3
        endif
        if (Nx /=0) then
            colU(Lfi+1:Lfi+6) =      [kip,       kip2, nCells+kip, nCells+kip2, 2*nCells+kip, 2*nCells+kip2 ] 
            valU(Lfi+1:Lfi+6) = -Vex*[ -4d0*sx,   sx,   4d0*sxy,    -sxy,        4d0*sxz,      -sxz    ]
            Lfi = Lfi + 6
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+6) =      [kjm,       kjm2, nCells+kjm, nCells+kjm2, 2*nCells+kjm, 2*nCells+kjm2 ] 
            valU(Lfi+1:Lfi+6) = -Vey*[ 4d0*sxy,   -sxy,   -4d0*sy,    sy,        -4d0*syz,      syz    ]
            Lfi = Lfi + 6
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+6) =      [kkm,       kkm2, nCells+kkm, nCells+kkm2, 2*nCells+kkm, 2*nCells+kkm2 ] 
            valU(Lfi+1:Lfi+6) = -Vez*[ 4d0*sxz,   -sxz,   -4d0*syz,    syz,         -4d0*sz,      sz    ]
            Lfi = Lfi + 6
        endif
    elseif ( nip == 0 .and. njp == 0 .and. nkp == 0 ) then  ! not  i+1 j+1 k+1     
        if (Nx /=0 .or. Ny /= 0 .or. Nz /= 0 ) then
            colU(Lfi+1:Lfi+3) =      [    L,                      L+nCells,                     L+2*nCells ]
            valU(Lfi+1:Lfi+3) = -3d0*[ Vex*sx+Vey*sxy+Vez*sxz,  Vex*sxy+Vey*sy+Vez*syz,  Vex*sxz+Vey*syz+Vez*sz] 
            Lfi = Lfi + 3
        endif
        if (Nx /=0) then
            colU(Lfi+1:Lfi+6) =      [kim,       kim2, nCells+kim, nCells+kim2, 2*nCells+kim, 2*nCells+kim2 ] 
            valU(Lfi+1:Lfi+6) = -Vex*[ -4d0*sx,   sx,   -4d0*sxy,    sxy,        -4d0*sxz,      sxz    ]
            Lfi = Lfi + 6
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+6) =      [kjm,       kjm2, nCells+kjm, nCells+kjm2, 2*nCells+kjm, 2*nCells+kjm2 ] 
            valU(Lfi+1:Lfi+6) = -Vey*[ -4d0*sxy,   sxy,   -4d0*sy,    sy,        -4d0*syz,      syz    ]
            Lfi = Lfi + 6
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+6) =      [kkm,       kkm2, nCells+kkm, nCells+kkm2, 2*nCells+kkm, 2*nCells+kkm2 ] 
            valU(Lfi+1:Lfi+6) = -Vez*[ -4d0*sxz,   sxz,   -4d0*syz,    syz,         -4d0*sz,      sz    ]
            Lfi = Lfi + 6
        endif
        !                        12 edges
        ! edges  along X
    elseif ( njp == 0  .and. nkm == 0 ) then  ! not  j+1  k-1  
        if (Ny /= 0 .or. Nz /= 0 ) then
            colU(Lfi+1:Lfi+2) =      [        L+nCells,    L+2*nCells ]
            valU(Lfi+1:Lfi+2) = -3d0*[ Vey*sy - Vez*sxz,  -Vey*syz + Vez*sz] 
            Lfi = Lfi + 2
        endif
        if (Nx /=0) then
            colU(Lfi+1:Lfi+4) =     [nCells+kip, nCells+kim, 2*nCells+kip, 2*nCells+kim ] 
            valU(Lfi+1:Lfi+4) = -Vex*[ sxy,        -sxy,        -sxz,          sxz ]
            Lfi = Lfi + 4
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+4) =     [nCells+kjm,  nCells+kjm2, 2*nCells+kjm, 2*nCells+kjm2 ] 
            valU(Lfi+1:Lfi+4) = -Vey*[ -4d0*sy,       sy,         4d0*syz,     -syz ]
            Lfi = Lfi + 4
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+4) =     [nCells+kkp,   nCells+kkp2, 2*nCells+kkp, 2*nCells+kkp2  ]
            valU(Lfi+1:Lfi+4) = -Vez*[ 4d0*syz,      -syz,        -4d0*sz,       sz ]
            Lfi = Lfi + 4
        endif
    elseif ( njm == 0  .and. nkm == 0 ) then   ! not  j-1  k-1   
        if (Ny /= 0 .or. Nz /= 0 ) then
            colU(Lfi+1:Lfi+2) =      [        L+nCells,    L+2*nCells ]
            valU(Lfi+1:Lfi+2) = -3d0*[ Vey*sy + Vez*sxz,  Vey*syz + Vez*sz] 
            Lfi = Lfi + 2
        endif
        if (Nx /=0) then
            colU(Lfi+1:Lfi+4) =     [nCells+kip, nCells+kim, 2*nCells+kip, 2*nCells+kim ] 
            valU(Lfi+1:Lfi+4) = -Vex*[ -sxy,        sxy,        -sxz,          sxz ]
            Lfi = Lfi + 4
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+4) =     [nCells+kjp,  nCells+kjp2, 2*nCells+kjp, 2*nCells+kjp2 ] 
            valU(Lfi+1:Lfi+4) = -Vey*[ -4d0*sy,       sy,         -4d0*syz,     syz ]
            Lfi = Lfi + 4
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+4) =     [nCells+kkp,   nCells+kkp2, 2*nCells+kkp, 2*nCells+kkp2  ]
            valU(Lfi+1:Lfi+4) = -Vez*[ -4d0*syz,      syz,        -4d0*sz,       sz ]
            Lfi = Lfi + 4
        endif
    elseif ( njp == 0  .and. nkp == 0 ) then   ! not  j+1  k+1
        if (Ny /= 0 .or. Nz /= 0 ) then
            colU(Lfi+1:Lfi+2) =      [        L+nCells,    L+2*nCells ]
            valU(Lfi+1:Lfi+2) = -3d0*[ Vey*sy + Vez*sxz,  Vey*syz + Vez*sz] 
            Lfi = Lfi + 2
        endif
        if (Nx /=0) then
            colU(Lfi+1:Lfi+4) =     [nCells+kip, nCells+kim, 2*nCells+kip, 2*nCells+kim ] 
            valU(Lfi+1:Lfi+4) = -Vex*[ sxy,        -sxy,        sxz,          -sxz ]
            Lfi = Lfi + 4
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+4) =     [nCells+kjm,  nCells+kjm2, 2*nCells+kjm, 2*nCells+kjm2 ] 
            valU(Lfi+1:Lfi+4) = -Vey*[ -4d0*sy,       sy,         -4d0*syz,      syz ]
            Lfi = Lfi + 4
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+4) =     [nCells+kkm,   nCells+kkm2, 2*nCells+kkm, 2*nCells+kkm2  ]
            valU(Lfi+1:Lfi+4) = -Vez*[ -4d0*syz,      syz,        -4d0*sz,        sz ]
            Lfi = Lfi + 4
        endif
    elseif ( njm == 0  .and. nkp == 0 ) then   ! not  j-1  k+1 
        if (Ny /= 0 .or. Nz /= 0 ) then
            colU(Lfi+1:Lfi+2) =      [        L+nCells,    L+2*nCells ]
            valU(Lfi+1:Lfi+2) = -3d0*[ Vey*sy - Vez*sxz,  -Vey*syz + Vez*sz] 
            Lfi = Lfi + 2
        endif
        if (Nx /=0) then
            colU(Lfi+1:Lfi+4) =      [nCells+kip, nCells+kim, 2*nCells+kip, 2*nCells+kim ] 
            valU(Lfi+1:Lfi+4) = -Vex*[ -sxy,        sxy,        sxz,          -sxz ]
            Lfi = Lfi + 4
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+4) =      [nCells+kjp,  nCells+kjp2, 2*nCells+kjp, 2*nCells+kjp2 ] 
            valU(Lfi+1:Lfi+4) = -Vey*[ -4d0*sy,       sy,         4d0*syz,     -syz ]
            Lfi = Lfi + 4
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+4) =      [nCells+kkm,   nCells+kkm2, 2*nCells+kkm, 2*nCells+kkm2  ]
            valU(Lfi+1:Lfi+4) = -Vez*[ 4d0*syz,      -syz,        -4d0*sz,       sz ]
            Lfi = Lfi + 4
        endif
        ! edges  along  Y
    elseif ( nip == 0  .and. nkm == 0 ) then  ! not  i+1  k-1  
        if (Nx /=0 .or. Nz /= 0 ) then
            colU(Lfi+1:Lfi+2) =      [        L,             L+2*nCells ]
            valU(Lfi+1:Lfi+2) = -3d0*[ Vex*sx - Vez*sxz,  -Vex*sxz + Vez*sz] 
            Lfi = Lfi + 2
        endif
        if (Nx /=0) then
            colU(Lfi+1:Lfi+4) =     [ kim,     kim2, 2*nCells+kim, 2*nCells+kim2 ] 
            valU(Lfi+1:Lfi+4) = -Vex*[ -4d0*sx,  sx,   4d0*sxz,      -sxz ]
            Lfi = Lfi + 4
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+4) =     [ kjp,  kjm, 2*nCells+kjp, 2*nCells+kjm ] 
            valU(Lfi+1:Lfi+4) = -Vey*[ sxy,  -sxy,   -syz,         syz ]
            Lfi = Lfi + 4
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+4) =     [ kkp,       kkp2, 2*nCells+kkp, 2*nCells+kkp2  ]
            valU(Lfi+1:Lfi+4) = -Vez*[ 4d0*sxz,  -sxz,    -4d0*sz,            sz ]
            Lfi = Lfi + 4
        endif
    elseif ( nim == 0  .and. nkm == 0 ) then   ! not  i-1  k-1  
        if (Nx /=0 .or. Nz /= 0 ) then
            colU(Lfi+1:Lfi+2) =      [        L,             L+2*nCells ]
            valU(Lfi+1:Lfi+2) = -3d0*[ Vex*sx + Vez*sxz,  +Vex*sxz + Vez*sz] 
            Lfi = Lfi + 2
        endif
        if (Nx /=0) then
            colU(Lfi+1:Lfi+4) =     [ kip,     kip2, 2*nCells+kip, 2*nCells+kip2 ] 
            valU(Lfi+1:Lfi+4) = -Vex*[ -4d0*sx,  sx,   -4d0*sxz,      sxz ]
            Lfi = Lfi + 4
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+4) =     [ kjp,  kjm, 2*nCells+kjp, 2*nCells+kjm ] 
            valU(Lfi+1:Lfi+4) = -Vey*[ -sxy,  sxy,   -syz,         syz ]
            Lfi = Lfi + 4
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+4) =     [ kkp,       kkp2, 2*nCells+kkp, 2*nCells+kkp2  ]
            valU(Lfi+1:Lfi+4) = -Vez*[ -4d0*sxz,  sxz,    -4d0*sz,        sz ]
            Lfi = Lfi + 4
        endif
    elseif ( nip == 0  .and. nkp == 0 ) then  ! not  i+1  k+1  
        if (Nx /=0 .or. Ny /= 0 ) then
            colU(Lfi+1:Lfi+2) =     [        L,             L+nCells    ]
            valU(Lfi+1:Lfi+2) = -3d0*[ Vex*sx + Vey*sxy,  +Vex*sxy + Vey*sy] 
            Lfi = Lfi + 2
        endif
        if (Nx /=0) then
            colU(Lfi+1:Lfi+4) =     [ kim,     kim2, 2*nCells+kim, 2*nCells+kim2 ] 
            valU(Lfi+1:Lfi+4) = -Vex*[ -4d0*sx,  sx,   -4d0*sxy,    sxy ]
            Lfi = Lfi + 4
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+4) =     [ kjp,  kjm, 2*nCells+kjp, 2*nCells+kjm ] 
            valU(Lfi+1:Lfi+4) = -Vey*[ sxy,  -sxy,   syz,         -syz ]
            Lfi = Lfi + 4
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+4) =     [ kkm,     kkm2, 2*nCells+kkm, 2*nCells+kkm2  ]
            valU(Lfi+1:Lfi+4) = -Vez*[ -4d0*sxz,   sxz,    -4d0*sz,            sz ]
            Lfi = Lfi + 4
        endif
    elseif ( nim == 0  .and. nkp == 0 ) then   ! not  i-1  k+1   
        if (Nx /=0 .or. Ny /= 0 ) then
            colU(Lfi+1:Lfi+2) =     [        L,             L+nCells    ]
            valU(Lfi+1:Lfi+2) = -3d0*[ Vex*sx - Vey*sxy,  -Vex*sxy + Vey*sy] 
            Lfi = Lfi + 2
        endif
        if (Nx /=0) then
            colU(Lfi+1:Lfi+4) =     [ kip,     kip2, 2*nCells+kip, 2*nCells+kip2 ] 
            valU(Lfi+1:Lfi+4) = -Vex*[ -4d0*sx,  sx,   4d0*sxz,      -sxz ]
            Lfi = Lfi + 4
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+4) =     [ kjp,  kjm, 2*nCells+kjp, 2*nCells+kjm ] 
            valU(Lfi+1:Lfi+4) = -Vey*[ -sxy,  sxy,   syz,         -syz ]
            Lfi = Lfi + 4
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+4) =     [ kkm,       kkm2, 2*nCells+kkm, 2*nCells+kkm2  ]
            valU(Lfi+1:Lfi+4) = -Vez*[ 4d0*sxz,  -sxz,    -4d0*sz,        sz ]
            Lfi = Lfi + 4
        endif
        ! edges  along Z 
    elseif ( nim == 0  .and. njm == 0 ) then   ! not  i-1 j-1 
        if (Nx /=0 .or. Ny /= 0 ) then
            colU(Lfi+1:Lfi+2) =     [        L,             L+nCells    ]
            valU(Lfi+1:Lfi+2) = -3d0*[ Vex*sx + Vey*sxy,  +Vex*sxy + Vey*sy] 
            Lfi = Lfi + 2
        endif
        if (Nx /=0) then
            colU(Lfi+1:Lfi+4) =     [ kip,     kip2, nCells+kip, nCells+kip2 ] 
            valU(Lfi+1:Lfi+4) = -Vex*[ -4d0*sx,  sx,   -4d0*sxy,    sxy ]
            Lfi = Lfi + 4
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+4) =     [ kjm,     kjm2, nCells+kjm, nCells+kjm2 ] 
            valU(Lfi+1:Lfi+4) = -Vey*[ -4d0*sxy,  sxy,   -4d0*sy,    sy ]
            Lfi = Lfi + 4
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+4) =     [ kkp,     kkm, nCells+kkp, nCells+kkm  ]
            valU(Lfi+1:Lfi+4) = -Vez*[ -sxz,   sxz,    -syz,       syz ]
            Lfi = Lfi + 4
        endif
    elseif ( nip == 0 .and. njm == 0 ) then    ! not  i+1 j-1 ! 
        if (Nx /=0 .or. Ny /= 0 ) then
            colU(Lfi+1:Lfi+2) =     [        L,             L+nCells    ]
            valU(Lfi+1:Lfi+2) = -3d0*[ Vex*sx - Vey*sxy,  -Vex*sxy - Vey*sy] 
            Lfi = Lfi + 2
        endif
        if (Nx /=0) then
            colU(Lfi+1:Lfi+4) =      [ kim,     kim2, nCells+kim, nCells+kim2 ] 
            valU(Lfi+1:Lfi+4) = -Vex*[ -4d0*sx,  sx,   4d0*sxy,    -sxy ]
            Lfi = Lfi + 4
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+4) =      [ kjp,     kjp2, nCells+kjp, nCells+kjp2 ] 
            valU(Lfi+1:Lfi+4) = -Vey*[ 4d0*sxy, -sxy,   -4d0*sy,    sy ]
            Lfi = Lfi + 4
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+4) =      [ kkp,   kkm, nCells+kkp, nCells+kkm  ]
            valU(Lfi+1:Lfi+4) = -Vez*[ sxz, -sxz,    -syz,       syz ]
            Lfi = Lfi + 4
        endif
    elseif ( nim == 0  .and. njp == 0 ) then   ! not  i-1 j+1 
        if (Nx /=0 .or. Ny /= 0 ) then
            colU(Lfi+1:Lfi+2) =     [        L,             L+nCells    ]
            valU(Lfi+1:Lfi+2) = -3d0*[Vex*sx - Vey*sxy,  -Vex*sxy + Vey*sy] 
            Lfi = Lfi + 2
        endif
        if (Nx /=0) then
            colU(Lfi+1:Lfi+4) =     [ kip,     kip2, nCells+kip, nCells+kip2 ] 
            valU(Lfi+1:Lfi+4) = -Vex*[ -4d0*sx,  sx,   4d0*sxy,    -sxy ]
            Lfi = Lfi + 4
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+4) =     [ kjm,     kjm2, nCells+kjm, nCells+kjm2 ] 
            valU(Lfi+1:Lfi+4) = -Vey*[  4d0*sxy, -sxy,   -4d0*sy,    sy ]
            Lfi = Lfi + 4
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+4) =     [ kkp,   kkm, nCells+kkp, nCells+kkm  ]
            valU(Lfi+1:Lfi+4) = -Vez*[ -sxz, sxz,    syz,       -syz ]
            Lfi = Lfi + 4
        endif
    elseif ( nip == 0  .and. njp == 0 ) then   ! not   i+1 j+1  
        if (Nx /=0 .or. Ny /= 0 ) then
            colU(Lfi+1:Lfi+2) =      [        L,              L+nCells   ]
            valU(Lfi+1:Lfi+2) = -3d0*[ (Vex*sx+Vey*sxy),  (Vex*sxy+Vey*sy) ]
            Lfi = Lfi + 2
        endif
        if (Nx /=0) then
            colU(Lfi+1:Lfi+4) =      [ kim,     kim2, nCells+kim, nCells+kim2 ] 
            valU(Lfi+1:Lfi+4) = -Vex*[ -4d0*sx,  sx,   -4d0*sxy,    sxy ]
            Lfi = Lfi + 4
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+4) =      [ kjm,     kjm2, nCells+kjm, nCells+kjm2 ] 
            valU(Lfi+1:Lfi+4) = -Vey*[ -4d0*sxy,  sxy,   -4d0*sy,    sy ]
            Lfi = Lfi + 4
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+4) =      [ kkp,   kkm, nCells+kkp, nCells+kkm  ]
            valU(Lfi+1:Lfi+4) = -Vez*[ sxz, -sxz,    syz,       -syz ]
            Lfi = Lfi + 4
        endif
        ! 6 faces  
    elseif (  nip == 0 .and. njp /= 0 .and. njm /= 0 .and. nkp /= 0 .and. nkm /= 0  ) then
         ! not   i+1
        if (Nx /=0) then
            colU(Lfi+1:Lfi+3) =          [ L,  kim,  kim2 ] 
            valU(Lfi+1:Lfi+3) = Vex*sx * [ -3d0, 4d0, -1d0 ]
            Lfi = Lfi + 3
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+2) =           [  kjp,  kjm]
            valU(Lfi+1:Lfi+2) = Vey*sxy * [ -1d0,  1d0]
            Lfi = Lfi + 2
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+2) =           [ kkp,  kkm ]
            valU(Lfi+1:Lfi+2) = Vez*sxz * [ -1d0, 1d0 ]
            Lfi = Lfi + 2
        endif
    elseif (  nim == 0 .and. njp /= 0 .and. njm /= 0 .and. nkp /= 0 .and. nkm /= 0  ) then
        ! not   i-1 x-
        if (Nx /=0) then
            colU(Lfi+1:Lfi+3) =          [ L,  kip,  kip2 ] 
            valU(Lfi+1:Lfi+3) = Vex*sx * [ -3d0, 4d0, -1d0 ]
            Lfi = Lfi + 3
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+2) =           [  kjp,  kjm]
            valU(Lfi+1:Lfi+2) = Vey*sxy * [ 1d0,  -1d0]
            Lfi = Lfi + 2
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+2) =          [ kkp,  kkm ]
            valU(Lfi+1:Lfi+2) =Vez*sxz * [ 1d0,  -1d0 ]
            Lfi = Lfi + 2
        endif
    elseif (  njp == 0 .and. nip /= 0 .and. nim /= 0 .and. nkp /= 0 .and. nkm /= 0  ) then
         ! not   j+1 
        if (Nx /=0) then
            colU(Lfi+1:Lfi+2) =  nCells + [ kip,  kim ] 
            valU(Lfi+1:Lfi+2) = Vex*sxy * [ -1d0, 1d0 ]
            Lfi = Lfi + 2
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+3) = nCells + [ L,  kjm,  kjm2 ]
            valU(Lfi+1:Lfi+3) = Vey*sy * [ -3d0, 4d0,  -1d0  ]
            Lfi = Lfi + 3
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+2) =  nCells + [ kkp,  kkm ]
            valU(Lfi+1:Lfi+2) = Vez*syz * [ -1d0,  1d0 ]
            Lfi = Lfi + 2
        endif
    elseif (  njm == 0 .and. nip /= 0 .and. nim /= 0 .and. nkp /= 0 .and. nkm /= 0  ) then
        ! not   j-1
        if (Nx /=0) then
            colU(Lfi+1:Lfi+2) =  nCells + [ kip,  kim ] 
            valU(Lfi+1:Lfi+2) = Vex*sxy * [ 1d0, -1d0 ]
            Lfi = Lfi + 2
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+3) = nCells + [ L,  kjp,  kjp2 ]
            valU(Lfi+1:Lfi+3) = Vey*sy * [ -3d0, 4d0,  -1d0  ]
            Lfi = Lfi + 3
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+2) =  nCells + [ kkp,  kkm ]
            valU(Lfi+1:Lfi+2) = Vez*syz * [ 1d0,  -1d0 ]
            Lfi = Lfi + 2
        endif
    elseif (  nkp == 0 .and. nip/= 0 .and. nim /= 0 .and. njp /= 0 .and. njm /= 0  ) then
        ! not   k+1 
        if (Nx /=0) then
            colU(Lfi+1:Lfi+2) = 2*nCells + [ kip,  kim ] 
            valU(Lfi+1:Lfi+2) =  Vex*sxz * [ -1d0, 1d0 ]
            Lfi = Lfi + 2
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+2) =2*nCells + [  kjp,  kjm ]
            valU(Lfi+1:Lfi+2) = Vey*syz * [ -1d0,  1d0 ]
            Lfi = Lfi + 2
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+3) =2*nCells + [L,  kkm,  kkm2 ]
            valU(Lfi+1:Lfi+3) =  Vez*sz * [-3d0, 4d0,  -1d0 ]
            Lfi = Lfi + 3
        endif
            
     elseif (  nkm == 0 .and. nip/= 0 .and. nim /= 0 .and. njp /= 0 .and. njm /= 0  ) then
        ! not   k-1
        if (Nx /=0) then
            colU(Lfi+1:Lfi+2) = 2*nCells + [ kip,  kim ] 
            valU(Lfi+1:Lfi+2) =  Vex*sxz * [ 1d0, -1d0 ]
            Lfi = Lfi + 2
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+2) =2*nCells + [  kjp,  kjm ]
            valU(Lfi+1:Lfi+2) = Vey*syz * [ 1d0,  -1d0 ]
            Lfi = Lfi + 2
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+3) =2*nCells + [ L,  kkp,  kkp2 ]
            valU(Lfi+1:Lfi+3) =  Vez*sz * [ -3d0, 4d0,  -1d0 ]
            Lfi = Lfi + 3
        endif
    else 
        Lfi = 0
        if (Nx /=0) then
            colU(Lfi+1:Lfi+3)=         [  L,  kip, kim ]
            valU(Lfi+1:Lfi+3)= Vex*sx* [ -2d0, 1d0,  1d0] 
            Lfi= Lfi + 3
            
            colU(Lfi+1:Lfi+4) =         nCells + [ kipjp, kipjm, kimjp, kimjm ]
            valU(Lfi+1:Lfi+4) = 0.25d0*Vex*sxy * [ 1d0,   -1d0,   -1d0,   1d0 ]
            Lfi = Lfi + 4 
            
            colU(Lfi+1:Lfi+4) =       2*nCells + [kipkp, kipkm, kimkp, kimkm ]
            valU(Lfi+1:Lfi+4) = 0.25d0*Vex*sxz * [ 1d0,   -1d0,   -1d0,   1d0]
            Lfi = Lfi + 4     
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+4) =                  [ kipjp, kipjm, kimjp, kimjm ]
            valU(Lfi+1:Lfi+4) = 0.25d0*Vey*sxy * [ 1d0,   -1d0,   -1d0,   +1d0 ]
            Lfi = Lfi + 4
            
            colU(Lfi+1:Lfi+3)= nCells + [ L,    kjp, kjm ]
            valU(Lfi+1:Lfi+3)=  Vey*sy* [ -2.d0, 1d0,  1d0] 
            Lfi= Lfi + 3

            colU(Lfi+1:Lfi+4) =       2*nCells + [kjpkp, kjpkm, kjmkp, kjmkm ]
            valU(Lfi+1:Lfi+4) = 0.25d0*Vey*syz * [ 1d0,   -1d0,   -1d0,   1d0]
            Lfi = Lfi + 4
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+4) =                  [ kipkp, kimkp, kipkm, kimkm ]
            valU(Lfi+1:Lfi+4) = 0.25d0*Vez*sxz * [ 1d0,   -1d0,   -1d0,   1d0 ]
            Lfi = Lfi + 4
            
            colU(Lfi+1:Lfi+4) =         nCells + [kjpkp, kjpkm, kjmkp, kjmkm ]
            valU(Lfi+1:Lfi+4) = 0.25d0*Vez*syz * [ 1d0,   -1d0,   -1d0,   1d0]
            Lfi = Lfi + 4
            
            colU(Lfi+1:Lfi+3)= 2*nCells + [ L,    kjp, kjm ]
            valU(Lfi+1:Lfi+3)=   Vez*sz * [ -2.d0, 1d0,  1d0] 
            Lfi= Lfi + 3
        endif
    endif

END SUBROUTINE update_rhs

SUBROUTINE motion_calc 
    DO i=1,3
        IF (fun_nod(n)% num_Vmech(i) == 0 ) THEN     !
            fun_nod(n)%Distance(i) = fun_nod(n)%Distance(i) + movestop(1)*fun_nod(n)%shift(i)
            fun_nod(n)%length(i) = nint(fun_nod(n)%Distance(i))
        ELSE                                         ! 
            fun_nod(n)%Distance(i) = fun_nod(n)%Distance(i) +  Vmech(fun_nod(n)%num_Vmech(i))%vely*dt/delta(i)
            fun_nod(n)%length(i) = nint(fun_nod(n)%Distance(i))
        ENDIF
    ENDDO
END SUBROUTINE motion_calc 

SUBROUTINE new_m 
    L = ceiling(REAL(m)/( REAL(sdx*sdy) ) ) 
    Lnew = L + fun_nod(n)%length(3)
    
    IF  ( Lnew > sdz-2  ) THEN 
        movestop(3) =0; Lnew = sdz-2
    ELSEIF(Lnew < 2    ) THEN
        movestop(3) =0;  Lnew = 2
    ELSEIF ( movestop(3) == 0 .and. (Lnew < sdz-2 .or. Lnew > 2)  ) THEN
             movestop(3) = 1 
    ENDIF
    IF (L == 1) THEN
        nij = m
    ELSE
        nij = m - (L-1)*sdx*sdy
    ENDIF
    j = ceiling( REAL(nij) / REAL(sdx) )  
    jnew = j + fun_nod(n)%length(2)
!-------------------------------------------------------------------------------------------
! comment/uncomment to check for out of bounds along the y-axis
!  1 variant
    IF  ( jnew > sdy-2  ) THEN 
        movestop(2) =0; jnew = sdy-2
    ELSEIF(jnew < 2    ) THEN
        movestop(2) =0;  jnew = 2
    ELSEIF ( movestop(2) == 0 .and. (jnew < sdy-2 .or. jnew > 2)  ) THEN
             movestop(2) = 1 
    ENDIF
!======================================================
!  2 variant
    ! IF  ( jnew > sdy  ) THEN 
        ! movestop(2) =0; jnew = sdy
    ! ELSEIF(jnew < 0    ) THEN
        ! movestop(2) =0;  jnew = 1
    ! ELSEIF ( movestop(2) == 0 .and. (jnew < sdy .or. jnew > 0 )  ) THEN
             ! movestop(2) = 1 
    ! ENDIF
!-----------------------------------------------------------------------------------------
    i = nij - (j - 1) * sdx  
    inew = i + fun_nod(n)%length(1)

    IF  ( inew > sdx-2  ) THEN  !
        movestop(1) =0;  inew = sdx-2
    ELSEIF(inew < 2    ) THEN !
        movestop(1) =0; inew = 2
    ELSEIF ( movestop(1) == 0 .and. (inew < sdx-2 .or. inew > 2)  ) THEN
             movestop(1) = 1 
    ENDIF

    m  = inew + sdx*(jnew-1) + sdx*sdy*(Lnew-1)
END SUBROUTINE new_m

END program EC3D

