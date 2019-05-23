! The seismicFDTD2D_gpu module is designed to be integrated with python via 
! f2py and CUDA gpu computing
!
! Compile using
!     f2py3 -c --fcompiler=gnu95 -m seismicfdtd2d_dp seismicFDTD2D_dp.f95
!
! Created by Steven Bernsen with T-minus one week to AGU
! University of Maine
! Department of Earth and Environmental Sciences 
! 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


module seismicFDTD2d_gpu
 
implicit none

contains

! =============================================================================
! ------------------------------- Sigma Updates--------------------------------
attributes (global) subroutine sigmaupdates(nz, nx, dx, dz, dt, &
            pml_x, pml_z, memory_vxx, memory_vzz, memory_vzx, memory_vxz, &
            vx, vz, sxx, szz, sxz, C, rho)

! pml_? is:     a, b, K, a_half, b_half, K_half
! C is C(:,:,1) = C11, C(:,:,2) = C22, C(:,:,3) = C12, C(:,:,4) = C66

    implicit none 

    integer, parameter :: dp = kind(0.d0)

    ! Define Inputs
    integer :: nz, nx 
    real(kind=dp) :: dx, dz, dt 
    real(kind-dp),dimension(nx,6) :: pml_x 
    real(kind=dp),dimension(nz,6) :: pml_z
    real(kind=dp),dimension(nx,nz) :: memory_vxx, memory_vzz, memory_vzx, memory_vxz
    real(kind=dp),dimension(nx,nz) :: vx, vz, sxx, szz, sxz, rho
    real(kind=dp),dimension(nx, nz, 4) :: C

    ! Define Disposable Variables 
    real(kind=dp) :: value_vxx, value_vzz, value_vzx, value_vxz

    do j = 2,nz
        do i = 1,NX-1

            value_vxx = (vx(i+1,j) - vx(i,j)) / dx
            value_vzz = (vz(i,j) - vz(i,j-1)) / dz

            memory_vxx(i,j) = pml_x(i,5) * memory_vxx(i,j) + pml_x(i,4) * value_vxx
            memory_vzz(i,j) = pml_z(j,2) * memory_vzz(i,j) + pml_z(j,1) * value_vzz

            value_vxx = value_vxx / pml_x(i,6) + memory_vxx(i,j)
            value_vzz = value_vzz / pml_z(j,3) + memory_vz(i,j)

            sxx(i,j) = sxx(i,j) + ( C(i,j,1) * value_vxx + C(i,j,3) * value_vzz) * dt
            szz(i,j) = szz(i,j) + ( C(i,j,3) * value_vxx + C(i,j,2) * value_vzz) * dt


        enddo
    enddo

    do j = 1,nz-1
        do i = 2,NX

            value_vzx = (vz(i,j) - vz(i-1,j)) / dx
            value_vxz = (vx(i,j+1) - vx(i,j)) / dz

            memory_vzx(i,j) = pml_x(i,2) * memory_vzx(i,j) + pml_x(i,1) * value_vzx
            memory_vxz(i,j) = pml_z(j,5) * memory_vxz(i,j) + pml_z(j,4) * value_vxz

            value_vzx = value_vzx / pml_x(i, 3) + memory_vzx(i,j)
            value_vxz = value_vxz / pml_z(j, 6) + memory_vxz(i,j)

            sxz(i,j) = sxz(i,j) + C(i,j,4) * (value_vzx + value_vxz) * dt

        enddo
    enddo

end subroutine sigmaupdates


! ----------------------------- Velocity Updates ------------------------------
attributes (global) subroutine velocityupdates(nz, nx, dx, dz, dt &
            pml_x, pml_z, memory_sxxx, memory_sxzz, memory_sxzx, memory_szzz, &
            vx, vz, sxx, szz, sxz)

    implicit none 

    integer, parameter :: dp = kind(0.d0)

    ! Define Inputs
    integer :: nz, nx 
    real(kind=dp) :: dx, dz, dt 
    real(kind-dp),dimension(nx,6) :: pml_x 
    real(kind=dp),dimension(nz,6) :: pml_z
    real(kind=dp),dimension(nx,nz) :: memory_sxxx, memory_sxzz, memory_sxzx, memory_szzz
    real(kind=dp),dimension(nx,nz) :: vx, vz, sxx, szz, sxz, rho


    do j = 2,nz
        do i = 2,NX

            value_sxxx = (sxx(i,j) - sxx(i-1,j)) / dx
            value_sxzz = (sxz(i,j) - sxz(i,j-1)) / dz

            memory_sxxx(i,j) = pml_x(i,2) * memory_sxxx(i,j) + pml_x(i,1) * value_sxxx
            memory_sxzz(i,j) = pml_z(j,2) * memory_sxzz(i,j) + pml_z(j,1) * value_sxzz

            value_sxxx = value_sxxx / pml_x(i,3) + memory_sxxx(i,j)
            value_sxzz = value_sxzz / pml_z(j,3) + memory_sxzz(i,j)

            vx(i,j) = vx(i,j) + (value_sxxx + value_sxzz) * dt / rho(i,j)

        enddo
    enddo

    do j = 1,nz-1
        do i = 1,NX-1

            value_sxzx = (sxz(i+1,j) - sxz(i,j)) / DX
            value_szzz = (szz(i,j+1) - szz(i,j)) / dz

            memory_sxzx(i,j) = pml_x(i,5) * memory_sxzx(i,j) + pml_x(i,4) * value_sxzx
            memory_szzz(i,j) = pml_z(j,5) * memory_szzz(i,j) + pml_z(j,4) * value_szzz

            value_sxzx = value_sxzx / pml_x(i,6) + memory_sxyx(i,j)
            value_szzz = value_szzz / pml_z(j,6) + memory_szzz(i,j)

            vz(i,j) = vz(i,j) + (value_sxzx + value_szzz) * dt / rho(i,j)

        enddo
    enddo

end subroutine velocityupdates



end module seismicFDTD2d_gpu




! ================================ Main Program ===============================
program seismicFDTD2D

use cudafor 
use seismicFDTD2d_gpu

implicit none

integer, parameter :: dp=kind(0.d0)

! total number of grid points in each direction of the grid
integer :: nx
integer :: nz

! thickness of the PML layer in grid points
integer :: npoints_pml
real(kind=dp), dimension(nx,nz, 4) :: C
real(kind=dp) :: f0

! total number of time steps
integer :: nstep

! time step in seconds. decreasing the time step improves the pml attenuation
! but it should be inversely proportional to the center frequency of the 
! source frequency 
real(kind=dp) :: DT
real(kind=dp) :: dx, dz 
! parameters for the source
real(kind=dp) :: t0
real(kind=dp), parameter :: factor = 1.d7

! source
integer,dimension(:) :: src
integer :: isource, jsource
! angle of source force clockwise with respect to vertical (Y) axis
real(kind=dp) :: ANGLE_FORCE


! display information on the screen from time to time
integer :: IT_DISPLAY

! value of PI
real(kind=dp), parameter :: PI = 3.141592653589793238462643d0

! conversion from degrees to radians
real(kind=dp), parameter :: DEGREES_TO_RADIANS = PI / 180.d0

! zeros
real(kind=dp), parameter :: ZERO = 0.d0

! large value for maximum
real(kind=dp), parameter :: HUGEVAL = 1.d+30

! velocity threshold above which we consider that the code became unstable
real(kind=dp), parameter :: STABILITY_THRESHOLD = 1.d+25

! main arrays
real(kind=dp), dimension(nx,nz) :: vx,vz,sxx,szz,sxz

! power to compute d0 profile. Increasing this value allows for a larger dampening gradient in the PML
real(kind=dp), parameter :: NPOWER = 2.d0

! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-11
real(kind=dp), parameter :: K_MAX_PML = 1.d1
real(kind=dp) :: ALPHA_MAX_PML

! arrays for the memory variables
! could declare these arrays in PML only to save a lot of memory, but proof of concept only here
real(kind=dp), dimension(nx,nz) :: &
    memory_vxx, &
    memory_vxz, &
    memory_vzx, &
    memory_vzz, &
    memory_sxxx, &
    memory_szzz, &
    memory_sxzx, &
    memory_sxzz

! 1D arrays for the damping profiles
real(kind=dp), dimension(NX) :: d_x,K_x,alpha_x,a_x,b_x,d_x_half,K_x_half,alpha_x_half,a_x_half,b_x_half
real(kind=dp), dimension(nz) :: d_z,K_z,alpha_z,a_z,b_z,d_z_half,K_z_half,alpha_z_half,a_z_half,b_z_half

real(kind=dp) :: thickness_PML_x,thickness_PML_z,xoriginleft,xoriginright,zoriginbottom,zorigintop
real(kind=dp) :: Rcoef,d0_x,d0_z,xval,yval,abscissa_in_PML,abscissa_normalized

! for the source
real(kind=dp) :: a,t,force_x,force_z,source_term

integer :: i,j,it

real(kind=dp) :: velocnorm

! for stability estimate
real(kind=dp) :: quasi_cp_max


! --------------------------- Assign some constants ---------------------------

isource = src(1)+npoints_pml
jsource = src(2)+npoints_pml

t0 = 1.0d0/f0
DT = minval( (/dx,dz/) )/ ( 2.0* sqrt( ( max( c/rho) ) ) )
ALPHA_MAX_PML = PI*f0 ! from Festa and Vilotte


! ------------------------------------ GPU ------------------------------------
! GPU variables
type(dim2) :: BlockH, GridH 


! Assign GPU Threads 
BlockH = dim2(Tx,Tz)
GridH = dim2( (Nx+Tx-1)/Tx, (Nz+Tz-1)/Tz )

! -----------------------------------------------------------------------------
!---
!--- program starts here
!---

!--- define profile of absorption in PML region

! thickness of the PML layer in meters
thickness_PML_x = NPOINTS_PML * dx
thickness_PML_z = NPOINTS_PML * dz

! reflection coefficient (INRIA report section 6.1) http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
Rcoef = 0.001d0

! check that NPOWER is okay
  if (NPOWER < 1) stop 'NPOWER must be greater than 1'


  d_x(:) = ZERO
  d_x_half(:) = ZERO
  K_x(:) = 1.d0
  K_x_half(:) = 1.d0
  alpha_x(:) = ZERO
  alpha_x_half(:) = ZERO
  a_x(:) = ZERO
  a_x_half(:) = ZERO

  d_z(:) = ZERO
  d_z_half(:) = ZERO
  K_z(:) = 1.d0
  K_z_half(:) = 1.d0
  alpha_z(:) = ZERO
  alpha_z_half(:) = ZERO
  a_z(:) = ZERO
  a_z_half(:) = ZERO

! damping in the X direction

! origin of the PML layer (position of right edge minus thickness, in meters)
  xoriginleft = dble(thickness_PML_x)
  xoriginright = dx * dble(NX-1) - thickness_PML_x

do i = 1,NX
  ! to compute d0 below, and for stability estimate
  quasi_cp_max = sqrt( max(c/rho) )
  ! compute d0 from INRIA report section 6.1 http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
  d0_x = - dble(NPOWER + 1) * quasi_cp_max * log(Rcoef) / (2.d0 * thickness_PML_x)
  d0_z = - dble(NPOWER + 1) * quasi_cp_max * log(Rcoef) / (2.d0 * thickness_PML_z)


    ! abscissa of current grid point along the damping profile
    xval = DX * dble(i-1)

    !---------- left edge
    ! define damping profile at the grid points
    abscissa_in_PML = xoriginleft - xval
    if (abscissa_in_PML >= ZERO) then
      abscissa_normalized = abscissa_in_PML / thickness_PML_x
      d_x(i) = d0_x * abscissa_normalized**NPOWER
      ! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
      K_x(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
      alpha_x(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
    endif

  ! define damping profile at half the grid points
    abscissa_in_PML = xoriginleft - (xval + DX/2.d0)
    if (abscissa_in_PML >= ZERO) then
      abscissa_normalized = abscissa_in_PML / thickness_PML_x
      d_x_half(i) = d0_x * abscissa_normalized**NPOWER
      ! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
      K_x_half(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
      alpha_x_half(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
    endif

!---------- right edge
      ! define damping profile at the grid points
      abscissa_in_PML = xval - xoriginright
      if (abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_x
        d_x(i) = d0_x * abscissa_normalized**NPOWER
        K_x(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_x(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
      endif

      ! define damping profile at half the grid points
      abscissa_in_PML = xval + DX/2.d0 - xoriginright
      if (abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_x
        d_x_half(i) = d0_x * abscissa_normalized**NPOWER
        K_x_half(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_x_half(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
      endif

    ! just in case, for -5 at the end
    if (alpha_x(i) < ZERO) alpha_x(i) = ZERO
    if (alpha_x_half(i) < ZERO) alpha_x_half(i) = ZERO

    b_x(i) = exp(- (d_x(i) / K_x(i) + alpha_x(i)) * DT)
    b_x_half(i) = exp(- (d_x_half(i) / K_x_half(i) + alpha_x_half(i)) * DT)

    ! this to avoid division by zero outside the PML
    if (abs(d_x(i)) > 1.d-6) a_x(i) = d_x(i) * (b_x(i) - 1.d0) / (K_x(i) * (d_x(i) + K_x(i) * alpha_x(i)))
    if (abs(d_x_half(i)) > 1.d-6) a_x_half(i) = d_x_half(i) * &
      (b_x_half(i) - 1.d0) / (K_x_half(i) * (d_x_half(i) + K_x_half(i) * alpha_x_half(i)))

  enddo

! damping in the Z direction

! origin of the PML layer (position of right edge minus thickness, in meters)
zoriginbottom = dble(thickness_PML_z)
zorigintop = dz * dble(nz-1) - thickness_PML_z

do j = 1,nz
  ! to compute d0 below, and for stability estimate
  quasi_cp_max = max(sqrt(c22(npoints_pml,j)/rho(npoints_pml,j)),&
    sqrt(c11(npoints_pml,j)/rho(npoints_pml,j)))
  ! compute d0 from INRIA report section 6.1 http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
  d0_z = - dble(NPOWER + 1) * quasi_cp_max * log(Rcoef) / (2.d0 * thickness_PML_z)

  ! abscissa of current grid point along the damping profile
  zval = dz * dble(j-1)

  !---------- bottom edge
    ! define damping profile at the grid points
    abscissa_in_PML = zoriginbottom - zval
    if (abscissa_in_PML >= ZERO) then
      abscissa_normalized = abscissa_in_PML / thickness_PML_z
      d_z(j) = d0_z * abscissa_normalized**NPOWER
      K_z(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
      alpha_z(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
    endif

    ! define damping profile at half the grid points
    abscissa_in_PML = zoriginbottom - (zval + dz/2.d0)
    if (abscissa_in_PML >= ZERO) then
      abscissa_normalized = abscissa_in_PML / thickness_PML_z
      d_z_half(j) = d0_z * abscissa_normalized**NPOWER
      K_z_half(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
      alpha_z_half(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
    endif


  !---------- top edge
    ! define damping profile at the grid points
    abscissa_in_PML = zval - zorigintop
    if (abscissa_in_PML >= ZERO) then
      abscissa_normalized = abscissa_in_PML / thickness_PML_z
      d_z(j) = d0_z * abscissa_normalized**NPOWER
      K_z(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
      alpha_z(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
    endif

    ! define damping profile at half the grid points
    abscissa_in_PML = zval + dz/2.d0 - zorigintop
    if (abscissa_in_PML >= ZERO) then
      abscissa_normalized = abscissa_in_PML / thickness_PML_z
      d_z_half(j) = d0_z * abscissa_normalized**NPOWER
      K_z_half(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
      alpha_z_half(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
    endif

  b_z(j) = exp(- (d_z(j) / K_z(j) + alpha_z(j)) * DT)
  b_z_half(j) = exp(- (d_z_half(j) / K_z_half(j) + alpha_z_half(j)) * DT)

  if (abs(d_z(j)) > 1.d-6) a_z(j) = d_z(j) * (b_z(j) - 1.d0) / (K_z(j) * (d_z(j) + K_z(j) * alpha_z(j)))
  if (abs(d_z_half(j)) > 1.d-6) a_z_half(j) = d_z_half(j) * &
    (b_z_half(j) - 1.d0) / (K_z_half(j) * (d_z_half(j) + K_z_half(j) * alpha_z_half(j)))

enddo

! ------------ Define Variables

! pml_x = a_x, b_x, K_x, a_x_half, b_x_half, K_x_half 
! pml_y = a_y, b_y, K_y, a_y_half, b_y_half, K_y_half 

! memory_vel = memory_vxx, memory_vyy, memory_vxy, memory_vyx
! memory_sig = memory_sxxx, memory_syyy, memory_sxyx, memory_sxyy 


! initialize arrays
vx(:,:) = ZERO
vz(:,:) = ZERO
sxx(:,:) = ZERO
szz(:,:) = ZERO
sxz(:,:) = ZERO

! PML
memory_vxx(:,:) = ZERO
memory_vxz(:,:) = ZERO
memory_vzx(:,:) = ZERO
memory_vzz(:,:) = ZERO
memory_sxxx(:,:) = ZERO
memory_szzz(:,:) = ZERO
memory_sxzx(:,:) = ZERO
memory_sxzz(:,:) = ZERO


!---
!---  beginning of time loop
!---

do it = 1,NSTEP

    istat=cudaThreadSynchronize()
    call sigmaupdates<<<GridH,BlockH>>>(nz, nx, dx, dz, dt, &
        pml_x, pml_z, memory_vxx, memory_vzz, memory_vzx, memory_vxz, &
        vx, vz, sxx, szz, sxz, C, rho)

        istat=cudaThreadSynchronize()
    call velocityupdates<<<GridH,BlockH>>>(nz, nx, dx, dz, dt &
        pml_x, pml_z, memory_sxxx, memory_sxzz, memory_sxzx, memory_szzz, &
        vx, vz, sxx, szz, sxz)

enddo


! add the source (force vector located at a given grid point)
a = pi*pi*f0*f0
t = dble(it-1)*DT

! Gaussian
source_term = factor * 2.d0*exp(-a*(t-t0)**2)

force_x = sin(ANGLE_FORCE * DEGREES_TO_RADIANS) * source_term
force_y = cos(ANGLE_FORCE * DEGREES_TO_RADIANS) * source_term

vx(isource,jsource) = vx(isource,jsource) + force_x * DT / rho(i,j)
vz(isource,jsource) = vz(isource,jsource) + force_z * DT / rho(i,j)

! Dirichlet conditions (rigid boundaries) on the edges or at the bottom of the PML layers
vx(1,:) = ZERO
vx(NX,:) = ZERO

vx(:,1) = ZERO
vx(:,nz) = ZERO

vz(1,:) = ZERO
vz(NX,:) = ZERO

vz(:,1) = ZERO
vz(:,nz) = ZERO

if (mod(it,IT_DISPLAY) == 0 .or. it == 1) then

    ! print maximum of norm of velocity
    velocnorm = maxval(sqrt(vx**2 + vz**2))
    print *,'Time step # ',it,' out of ',NSTEP
    print *,'Time: ',sngl((it-1)*DT),' seconds'
  
    if (velocnorm > STABILITY_THRESHOLD) stop 'code became unstable and blew up'
  
    call write_image(vx, nx, nz, it, 'Vx')
    call write_image(vz, nx, nz, it, 'Vz')
    
endif



end program seismicFDTD2D