! The seismicFDTD2D module is designed to be integrated with python via f2py. 
!
! Compile using
!     f2py3 -c --fcompiler=gnu95 -m seismicfdtd2d_dp seismicFDTD2D_dp.f95
!
! Created by Steven Bernsen with T-minus one week to AGU
! University of Maine
! Department of Earth and Environmental Sciences 
! 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


module seismicFDTD2d

implicit none

contains

subroutine doall(im, mlist, nx, ny, dx, dy, npoints_pml, & 
                  src, f0, nstep, it_display, angle_force)
! DOALL This is kind of a wrapper function for the subsequent subroutines 
! because this will be implemented via Python or some other dynamic front end 
! language. Of course I would name this in the fashion of the Computer Programs in
! Seismology naming. 
!
! INPUT
!   im (INTEGER) - m-by-n array of integer values corresponding to different
!         materials.
!   mlist (REAL) - the p-by-13 array containing in each column:
!
!   MATERIAL_ID,TEMPERATURE,PRESSURE,C11,C12,C13,C22,C23,C33,C44,C55,C66,RHO
!   
!   nx,ny (INTEGER) - the shape variables of the input arrays in order to 
!         allocate space and other static language headaches 
!   dx,dy (REAL) - the inteval length values in the x and y directions
!   npoints_pml (INTEGER) - the number of points for the CPML layer. This is 
!         a constant value for all four sides. 
!   rcx,src (INTEGER) - the indices of the locations of the receivers 
!   f0 (REAL) - the center frequency of the source function. The time step is
!         inversely proportional to the center frequency 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

implicit none

integer,parameter :: dp = kind(0.d0)
integer :: nx, ny, nstep, npoints_pml, it_display
integer,dimension(nx,ny) :: im
! integer,dimension(:,:) :: rcx
integer,dimension(:) :: src 
real(kind=dp), dimension(:,:) :: mlist
real(kind=dp) :: f0, angle_force
real(kind=dp) :: dx, dy
real(kind=dp), dimension(nx+2*npoints_pml,ny+2*npoints_pml) :: c11, c12, c22, c66, rho
! character(len=6) :: src_type

!f2py3 intent(in) :: im, mlist, nx, ny, dx, dy, npoints_pml, src
!f2py3 intent(in) :: f0, nstep, it_display, angle_force
!f2py3 intent(hide), depend(im) :: nx = shape(im, 0), ny = shape(im,1)

! Preallocate arrays
c11(:,:) = 0.0
c12(:,:) = 0.0
c66(:,:) = 0.0
rho(:,:) = 0.0


! Setup arrays
call stiffness_arrays(im, mlist, c11, c12, c22, c66, rho, npoints_pml)

call seismic_cpml_2d(nx+2*npoints_pml, ny+2*npoints_pml, c11, c12, c22, c66, rho, dx, dy, &
                      npoints_pml, src, f0, nstep, it_display, angle_force)


end subroutine doall


!==============================================================================
subroutine stiffness_arrays(im, mlist, c11, c12, c22, c66, rho, npoints_pml) 
! STIFFNESS_ARRAYS takes a matrix containing the material integer identifiers 
! and creates the same size array for each independent coefficient of the 
! stiffness matrix along with a density matrix. Since we ae using PML
! boundaries, we will extend the the boundary values through the PML region.
! 
! INPUT 
!   im (INTEGER)  
!   mlist (REAL)
!   c11(i,j), c12(i,j), c22(i,j), c66, rho(i,j) (REAL) -
!   npoints_pml (INTEGER) - the 
!   
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

implicit none 

integer,parameter :: dp = kind(0.d0)
integer,dimension(:,:) :: im
integer :: m, n, i, j, npoints_pml
real(kind=dp), dimension(:,:) :: mlist
real(kind=dp), dimension(:, :) :: c11, c12, c22, c66, rho

m = size(im, 1)
n = size(im, 2)

do i=npoints_pml+1, m+npoints_pml
  do j=npoints_pml+1, n+npoints_pml
    c11(i,j) = mlist(im(i-npoints_pml,j-npoints_pml), 2)
    c12(i,j) = mlist(im(i-npoints_pml,j-npoints_pml), 3) 
    c22(i,j) = mlist(im(i-npoints_pml,j-npoints_pml), 5)
    c66(i,j) = mlist(im(i-npoints_pml,j-npoints_pml), 10)
    rho(i,j) = mlist(im(i-npoints_pml,j-npoints_pml), 11) 
  end do
end do


! Extend the boundary values of the stiffnesses into the PML region
do i = 1,npoints_pml+1
  ! top and bottom
  c11( i, : ) = c11(npoints_pml+1,:)
  c11( m+npoints_pml-1+i, : ) = c11(m+npoints_pml,:)

  c12( i, : ) = c12(npoints_pml+1,:)
  c12( m+npoints_pml-1+i, : ) = c12(m+npoints_pml,:)
  
  c22( i, : ) = c22(npoints_pml+1,:)
  c22( m+npoints_pml-1+i, : ) = c22(m+npoints_pml,:)
  
  c66( i, : ) = c66(npoints_pml+1,:)
  c66( m+npoints_pml-1+i, : ) = c66(m+npoints_pml,:)
  
  rho( i, : ) = rho(npoints_pml+1,:)
  rho( m+npoints_pml-1+i, : ) = rho(m+npoints_pml,:)

  ! ! left and right
  c11( :, i ) = c11(:, npoints_pml+1)
  c11( :, n+npoints_pml-1+i ) = c11(:,n+npoints_pml)
  
  c12( :, i ) = c12(:, npoints_pml+1)
  c12( :, n+npoints_pml-1+i ) = c12(:,n+npoints_pml)
  
  c22( :, i ) = c22(:, npoints_pml+1)
  c22( :, n+npoints_pml-1+i ) = c22(:,n+npoints_pml)
  
  c66( :, i ) = c66(:, npoints_pml+1)
  c66( :, n+npoints_pml-1+i ) = c66(:,n+npoints_pml)

  rho( :, i ) = rho(:, npoints_pml+1)
  rho( :, n+npoints_pml-1+i ) = rho(:,n+npoints_pml)

end do 


end subroutine stiffness_arrays

!==============================================================================


subroutine seismic_cpml_2d(nx, ny, c11, c12, c22, c66, rho, dx, dy, &
                      npoints_pml, src, f0, nstep, it_display, angle_force)

! 2D elastic finite-difference code in velocity and stress formulation
! with Convolutional-PML (C-PML) absorbing conditions for an anisotropic medium

! Dimitri Komatitsch, University of Pau, France, April 2007.
! Anisotropic implementation by Roland Martin and Dimitri Komatitsch, University of Pau, France, April 2007.

! The second-order staggered-grid formulation of Madariaga (1976) and Virieux (1986) is used:
!
!            ^ y
!            |
!            |
!
!            +-------------------+
!            |                   |
!            |                   |
!            |                   |
!            |                   |
!            |        v_y        |
!   sigma_xy +---------+         |
!            |         |         |
!            |         |         |
!            |         |         |
!            |         |         |
!            |         |         |
!            +---------+---------+  ---> x
!           v_x    sigma_xx
!                  sigma_yy
!
! IMPORTANT : all our CPML codes work fine in single precision as well (which is significantly faster).
!             If you want you can thus force automatic conversion to single precision at compile time
!             or change all the declarations and constants in the code from real(kind=dp) to single.
!
! INPUT
!   im (INTEGER)
!   nx, ny (INTEGER)
!   c11, c12, c22, c66, rho (REAL)
!   dx, dy (REAL)
!   npoints_pml (INTEGER) - the thickness of the pml
!


implicit none

integer, parameter :: dp=kind(0.d0)

! total number of grid points in each direction of the grid
integer :: nx
integer :: ny

! thickness of the PML layer in grid points
integer :: npoints_pml
! integer, dimension(nx,ny)
real(kind=dp), dimension(nx,ny) :: c11, c12, c22, c66, rho
real(kind=dp) :: f0

! total number of time steps
integer :: nstep

! time step in seconds. decreasing the time step improves the pml attenuation
! but it should be inversely proportional to the center frequency of the 
! source frequency 
real(kind=dp) :: DT
real(kind=dp) :: dx, dy 
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
real(kind=dp), dimension(nx,ny) :: vx,vy,sigmaxx,sigmayy,sigmaxy

! power to compute d0 profile. Increasing this value allows for a larger dampening gradient in the PML
real(kind=dp), parameter :: NPOWER = 2.d0

! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-11
real(kind=dp), parameter :: K_MAX_PML = 1.d1
real(kind=dp) :: ALPHA_MAX_PML

! arrays for the memory variables
! could declare these arrays in PML only to save a lot of memory, but proof of concept only here
real(kind=dp), dimension(NX,NY) :: &
    memory_dvx_dx, &
    memory_dvx_dy, &
    memory_dvy_dx, &
    memory_dvy_dy, &
    memory_dsigmaxx_dx, &
    memory_dsigmayy_dy, &
    memory_dsigmaxy_dx, &
    memory_dsigmaxy_dy

real(kind=dp) :: &
    value_dvx_dx, &
    value_dvx_dy, &
    value_dvy_dx, &
    value_dvy_dy, &
    value_dsigmaxx_dx, &
    value_dsigmayy_dy, &
    value_dsigmaxy_dx, &
    value_dsigmaxy_dy

! 1D arrays for the damping profiles
real(kind=dp), dimension(NX) :: d_x,K_x,alpha_x,a_x,b_x,d_x_half,K_x_half,alpha_x_half,a_x_half,b_x_half
real(kind=dp), dimension(NY) :: d_y,K_y,alpha_y,a_y,b_y,d_y_half,K_y_half,alpha_y_half,a_y_half,b_y_half

real(kind=dp) :: thickness_PML_x,thickness_PML_y,xoriginleft,xoriginright,yoriginbottom,yorigintop
real(kind=dp) :: Rcoef,d0_x,d0_y,xval,yval,abscissa_in_PML,abscissa_normalized

! for the source
real(kind=dp) :: a,t,force_x,force_y,source_term

integer :: i,j,it

real(kind=dp) :: velocnorm

! for stability estimate
real(kind=dp) :: quasi_cp_max


! ------------------------ Assign some constants -----------------------

isource = src(1)+npoints_pml
jsource = src(2)+npoints_pml


! nrec = size(rcx,1) 
! allocate( sisvx(nstep, nrec ), &
  ! sisvy(nstep, nrec ) )

 ! dx/minval( (/ c11/rho, c22/rho, c12/rho, c66/rho /) )

t0 = 1.0d0/f0
DT = dx/ ( 2.0* sqrt( ( maxval( (/ c11/rho, c22/rho, c12/rho, c66/rho /) ) ) ) )!dx/(256*f0)!dt)
ALPHA_MAX_PML = PI*f0 ! from Festa and Vilotte
! ----------------------------------------------------------------------
!---
!--- program starts here
!---

!--- define profile of absorption in PML region

! thickness of the PML layer in meters
thickness_PML_x = NPOINTS_PML * DX
thickness_PML_y = NPOINTS_PML * DY

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

  d_y(:) = ZERO
  d_y_half(:) = ZERO
  K_y(:) = 1.d0
  K_y_half(:) = 1.d0
  alpha_y(:) = ZERO
  alpha_y_half(:) = ZERO
  a_y(:) = ZERO
  a_y_half(:) = ZERO

! damping in the X direction

! origin of the PML layer (position of right edge minus thickness, in meters)
  xoriginleft = dble(thickness_PML_x)
  xoriginright = dx * dble(NX-1) - thickness_PML_x

do i = 1,NX
  ! to compute d0 below, and for stability estimate
  quasi_cp_max = max(sqrt(c22(i,npoints_pml)/rho(i,npoints_pml)),&
    sqrt(c11(i,npoints_pml)/rho(i,npoints_pml)))
  ! compute d0 from INRIA report section 6.1 http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
  d0_x = - dble(NPOWER + 1) * quasi_cp_max * log(Rcoef) / (2.d0 * thickness_PML_x)
  d0_y = - dble(NPOWER + 1) * quasi_cp_max * log(Rcoef) / (2.d0 * thickness_PML_y)


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

! damping in the Y direction

! origin of the PML layer (position of right edge minus thickness, in meters)
yoriginbottom = dble(thickness_PML_y)
yorigintop = dy * dble(NY-1) - thickness_PML_y

do j = 1,NY
  ! to compute d0 below, and for stability estimate
  quasi_cp_max = max(sqrt(c22(npoints_pml,j)/rho(npoints_pml,j)),&
    sqrt(c11(npoints_pml,j)/rho(npoints_pml,j)))
  ! compute d0 from INRIA report section 6.1 http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
  d0_x = - dble(NPOWER + 1) * quasi_cp_max * log(Rcoef) / (2.d0 * thickness_PML_x)
  d0_y = - dble(NPOWER + 1) * quasi_cp_max * log(Rcoef) / (2.d0 * thickness_PML_y)

  ! abscissa of current grid point along the damping profile
  yval = DY * dble(j-1)

  !---------- bottom edge
    ! define damping profile at the grid points
    abscissa_in_PML = yoriginbottom - yval
    if (abscissa_in_PML >= ZERO) then
      abscissa_normalized = abscissa_in_PML / thickness_PML_y
      d_y(j) = d0_y * abscissa_normalized**NPOWER
      K_y(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
      alpha_y(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
    endif

    ! define damping profile at half the grid points
    abscissa_in_PML = yoriginbottom - (yval + DY/2.d0)
    if (abscissa_in_PML >= ZERO) then
      abscissa_normalized = abscissa_in_PML / thickness_PML_y
      d_y_half(j) = d0_y * abscissa_normalized**NPOWER
      K_y_half(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
      alpha_y_half(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
    endif


  !---------- top edge
    ! define damping profile at the grid points
    abscissa_in_PML = yval - yorigintop
    if (abscissa_in_PML >= ZERO) then
      abscissa_normalized = abscissa_in_PML / thickness_PML_y
      d_y(j) = d0_y * abscissa_normalized**NPOWER
      K_y(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
      alpha_y(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
    endif

    ! define damping profile at half the grid points
    abscissa_in_PML = yval + DY/2.d0 - yorigintop
    if (abscissa_in_PML >= ZERO) then
      abscissa_normalized = abscissa_in_PML / thickness_PML_y
      d_y_half(j) = d0_y * abscissa_normalized**NPOWER
      K_y_half(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
      alpha_y_half(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
    endif

  b_y(j) = exp(- (d_y(j) / K_y(j) + alpha_y(j)) * DT)
  b_y_half(j) = exp(- (d_y_half(j) / K_y_half(j) + alpha_y_half(j)) * DT)

  if (abs(d_y(j)) > 1.d-6) a_y(j) = d_y(j) * (b_y(j) - 1.d0) / (K_y(j) * (d_y(j) + K_y(j) * alpha_y(j)))
  if (abs(d_y_half(j)) > 1.d-6) a_y_half(j) = d_y_half(j) * &
    (b_y_half(j) - 1.d0) / (K_y_half(j) * (d_y_half(j) + K_y_half(j) * alpha_y_half(j)))

enddo

! initialize arrays
vx(:,:) = ZERO
vy(:,:) = ZERO
sigmaxx(:,:) = ZERO
sigmayy(:,:) = ZERO
sigmaxy(:,:) = ZERO

! PML
memory_dvx_dx(:,:) = ZERO
memory_dvx_dy(:,:) = ZERO
memory_dvy_dx(:,:) = ZERO
memory_dvy_dy(:,:) = ZERO
memory_dsigmaxx_dx(:,:) = ZERO
memory_dsigmayy_dy(:,:) = ZERO
memory_dsigmaxy_dx(:,:) = ZERO
memory_dsigmaxy_dy(:,:) = ZERO

!---
!---  beginning of time loop
!---

do it = 1,NSTEP
  !------------------------------------------------------------
  ! compute stress sigma and update memory variables for C-PML
  !------------------------------------------------------------
  do j = 2,NY
    do i = 1,NX-1

      value_dvx_dx = (vx(i+1,j) - vx(i,j)) / DX
      value_dvy_dy = (vy(i,j) - vy(i,j-1)) / DY

      memory_dvx_dx(i,j) = b_x_half(i) * memory_dvx_dx(i,j) + a_x_half(i) * value_dvx_dx
      memory_dvy_dy(i,j) = b_y(j) * memory_dvy_dy(i,j) + a_y(j) * value_dvy_dy

      value_dvx_dx = value_dvx_dx / K_x_half(i) + memory_dvx_dx(i,j)
      value_dvy_dy = value_dvy_dy / K_y(j) + memory_dvy_dy(i,j)

      sigmaxx(i,j) = sigmaxx(i,j) + ( c11(i,j) * value_dvx_dx + c12(i,j) * value_dvy_dy) * DT
      sigmayy(i,j) = sigmayy(i,j) + ( c12(i,j) * value_dvx_dx + c22(i,j) * value_dvy_dy) * DT


    enddo
  enddo

  do j = 1,NY-1
    do i = 2,NX

      value_dvy_dx = (vy(i,j) - vy(i-1,j)) / DX
      value_dvx_dy = (vx(i,j+1) - vx(i,j)) / DY

      memory_dvy_dx(i,j) = b_x(i) * memory_dvy_dx(i,j) + a_x(i) * value_dvy_dx
      memory_dvx_dy(i,j) = b_y_half(j) * memory_dvx_dy(i,j) + a_y_half(j) * value_dvx_dy

      value_dvy_dx = value_dvy_dx / K_x(i) + memory_dvy_dx(i,j)
      value_dvx_dy = value_dvx_dy / K_y_half(j) + memory_dvx_dy(i,j)

      sigmaxy(i,j) = sigmaxy(i,j) + c66(i,j) * (value_dvy_dx + value_dvx_dy) * DT

    enddo
  enddo

!--------------------------------------------------------
! compute velocity and update memory variables for C-PML
!--------------------------------------------------------

  do j = 2,NY
    do i = 2,NX

      value_dsigmaxx_dx = (sigmaxx(i,j) - sigmaxx(i-1,j)) / DX
      value_dsigmaxy_dy = (sigmaxy(i,j) - sigmaxy(i,j-1)) / DY

      memory_dsigmaxx_dx(i,j) = b_x(i) * memory_dsigmaxx_dx(i,j) + a_x(i) * value_dsigmaxx_dx
      memory_dsigmaxy_dy(i,j) = b_y(j) * memory_dsigmaxy_dy(i,j) + a_y(j) * value_dsigmaxy_dy

      value_dsigmaxx_dx = value_dsigmaxx_dx / K_x(i) + memory_dsigmaxx_dx(i,j)
      value_dsigmaxy_dy = value_dsigmaxy_dy / K_y(j) + memory_dsigmaxy_dy(i,j)

      vx(i,j) = vx(i,j) + (value_dsigmaxx_dx + value_dsigmaxy_dy) * DT / rho(i,j)

    enddo
  enddo

  do j = 1,NY-1
    do i = 1,NX-1

      value_dsigmaxy_dx = (sigmaxy(i+1,j) - sigmaxy(i,j)) / DX
      value_dsigmayy_dy = (sigmayy(i,j+1) - sigmayy(i,j)) / DY

      memory_dsigmaxy_dx(i,j) = b_x_half(i) * memory_dsigmaxy_dx(i,j) + a_x_half(i) * value_dsigmaxy_dx
      memory_dsigmayy_dy(i,j) = b_y_half(j) * memory_dsigmayy_dy(i,j) + a_y_half(j) * value_dsigmayy_dy

      value_dsigmaxy_dx = value_dsigmaxy_dx / K_x_half(i) + memory_dsigmaxy_dx(i,j)
      value_dsigmayy_dy = value_dsigmayy_dy / K_y_half(j) + memory_dsigmayy_dy(i,j)

      vy(i,j) = vy(i,j) + (value_dsigmaxy_dx + value_dsigmayy_dy) * DT / rho(i,j)

    enddo
  enddo

  ! add the source (force vector located at a given grid point)
  a = pi*pi*f0*f0
  t = dble(it-1)*DT

  ! Gaussian
  source_term = factor * 2.d0*exp(-a*(t-t0)**2)

  ! first derivative of a Gaussian
  ! source_term = - factor * 2.d0*a*(t-t0)*exp(-a*(t-t0)**2)

  ! Ricker source time function (second derivative of a Gaussian)
  ! source_term = factor * (1.d0 - 2.d0*a*(t-t0)**2)*exp(-a*(t-t0)**2)

  force_x = sin(ANGLE_FORCE * DEGREES_TO_RADIANS) * source_term
  force_y = cos(ANGLE_FORCE * DEGREES_TO_RADIANS) * source_term

  ! define location of the source
  ! i = ISOURCE
  ! j = JSOURCE

  vx(isource,jsource) = vx(isource,jsource) + force_x * DT / rho(i,j)
  vy(isource,jsource) = vy(isource,jsource) + force_y * DT / rho(i,j)

  ! Dirichlet conditions (rigid boundaries) on the edges or at the bottom of the PML layers
  vx(1,:) = ZERO
  vx(NX,:) = ZERO

  vx(:,1) = ZERO
  vx(:,NY) = ZERO

  vy(1,:) = ZERO
  vy(NX,:) = ZERO

  vy(:,1) = ZERO
  vy(:,NY) = ZERO

  ! store seismograms
  ! do irec = 1,NREC
  !   sisvx(it,irec) = vx(rcx(irec,1),rcx(irec,2))
  !   sisvy(it,irec) = vy(rcx(irec,1),rcx(irec,2))
  ! enddo

  ! output information
  if (mod(it,IT_DISPLAY) == 0 .or. it == 1) then

  ! print maximum of norm of velocity
  velocnorm = maxval(sqrt(vx**2 + vy**2))
  print *,'Time step # ',it,' out of ',NSTEP
  print *,'Time: ',sngl((it-1)*DT),' seconds'
  print *,'Max norm velocity vector V (m/s) = ',velocnorm
  print *
  ! check stability of the code, exit if unstable
  if (velocnorm > STABILITY_THRESHOLD) stop 'code became unstable and blew up'

  call create_color_image(vx,NX,NY,it,ISOURCE,JSOURCE, &
                       NPOINTS_PML,1)
  call create_color_image(vy,NX,NY,it,ISOURCE,JSOURCE, &
                       NPOINTS_PML,2)

  call write_image(vx, nx, ny, it, 'Vx')
  call write_image(vy, nx, ny, it, 'Vy')
  
  endif

enddo   ! end of time loop



! save seismograms
! call write_seismograms(sisvx,sisvy,nstep,nrec,DT)

end subroutine seismic_cpml_2d


!==============================================================================
subroutine write_image(image_data, nx, ny, it, channel)

implicit none

integer, parameter :: dp = kind(0.d0)
integer :: nx, ny, it
real(kind=dp) :: image_data(nx, ny)
character(len=2) :: channel
character(len=100) :: filename

WRITE (filename, "(a2, i6.6, '.dat')" ) channel, it

open(unit = 10, form = 'unformatted', file = trim(filename) )
write(10) image_data

close(unit = 10)

end subroutine write_image

!==============================================================================
! subroutine write_seismograms(sisvx,sisvy,nt,nrec,DELTAT)

! implicit none

! integer, parameter :: dp=kind(0.d0)
! integer nt,nrec
! real(kind=dp) DELTAT

! real(kind=dp) sisvx(nt,nrec)
! real(kind=dp) sisvy(nt,nrec)

! integer irec,it

! character(len=100) file_name

! ! X component
! do irec=1,nrec
!   write(file_name,"('Vx_file_',i3.3,'.dat')") irec
!   open(unit=11,file=file_name,status='unknown')
!   do it=1,nt
!     write(11,*) sngl(dble(it-1)*DELTAT),',',sngl(sisvx(it,irec))
!   end do
!     close(11)
! enddo

! ! Y component
! do irec=1,nrec
!   write(file_name,"('Vy_file_',i3.3,'.dat')") irec
!   open(unit=11,file=file_name,status='unknown')
!   do it=1,nt
!     write(11,*) sngl(dble(it-1)*DELTAT),',',sngl(sisvy(it,irec))
!   end do
!   close(11)
! end do

! end subroutine write_seismograms

!==============================================================================
!----
!----  routine to create a color image of a given vector component
!----  the image is created in PNM format and then converted to GIF
!----

subroutine create_color_image(image_data_2D,NX,NY,it,ISOURCE,JSOURCE, &
            NPOINTS_PML,field_number)

implicit none
integer, parameter :: dp = kind(0.d0)

! non linear display to enhance small amplitudes for graphics
real(kind=dp), parameter :: POWER_DISPLAY = 0.30d0

! amplitude threshold above which we draw the color point
real(kind=dp), parameter :: cutvect = 1.0d-1

! use black or white background for points that are below the threshold
logical, parameter :: WHITE_BACKGROUND = .true.

! size of cross and square in pixels drawn to represent the source and the receivers
integer, parameter :: width_cross = 5, thickness_cross = 1, size_square = 3

integer NX,NY,it,field_number,ISOURCE,JSOURCE,NPOINTS_PML

real(kind=dp), dimension(NX,NY) :: image_data_2D

integer :: ix,iy

character(len=100) :: file_name,system_command

integer :: R, G, B

real(kind=dp) :: normalized_value,max_amplitude

! open image file and create system command to convert image to more convenient format
! use the "convert" command from ImageMagick http://www.imagemagick.org
if (field_number == 1) then
  write(file_name,"('image',i6.6,'_Vx.pnm')") it
  write(system_command,"('convert image',i6.6,'_Vx.pnm image',i6.6,'_Vx.gif ; rm image',i6.6,'_Vx.pnm')") it,it,it
else if (field_number == 2) then
  write(file_name,"('image',i6.6,'_Vz.pnm')") it
  write(system_command,"('convert image',i6.6,'_Vz.pnm image',i6.6,'_Vz.gif ; rm image',i6.6,'_Vz.pnm')") it,it,it
endif

open(unit=27, file=file_name, status='unknown')

write(27,"('P3')") ! write image in PNM P3 format

write(27,*) NX,NY ! write image size
write(27,*) '255' ! maximum value of each pixel color

! compute maximum amplitude
max_amplitude = maxval(abs(image_data_2D))

! image starts in upper-left corner in PNM format
do iy=NY,1,-1
  do ix=1,NX

! define data as vector component normalized to [-1:1] and rounded to nearest integer
! keeping in mind that amplitude can be negative
  normalized_value = image_data_2D(ix,iy) / max_amplitude

! suppress values that are outside [-1:+1] to avoid small edge effects
  if (normalized_value < -1.d0) normalized_value = -1.d0
  if (normalized_value > 1.d0) normalized_value = 1.d0

! draw an orange cross to represent the source
  if ((ix >= ISOURCE - width_cross .and. ix <= ISOURCE + width_cross .and. &
      iy >= JSOURCE - thickness_cross .and. iy <= JSOURCE + thickness_cross) .or. &
     (ix >= ISOURCE - thickness_cross .and. ix <= ISOURCE + thickness_cross .and. &
      iy >= JSOURCE - width_cross .and. iy <= JSOURCE + width_cross)) then
    R = 255
    G = 157
    B = 0

! display two-pixel-thick black frame around the image
else if (ix <= 2 .or. ix >= NX-1 .or. iy <= 2 .or. iy >= NY-1) then
    R = 0
    G = 0
    B = 0

! display edges of the PML layers
else if ((ix == NPOINTS_PML) .or. &
        (ix == NX - NPOINTS_PML) .or. &
        (iy == NPOINTS_PML) .or. &
        (iy == NY - NPOINTS_PML)) then
    R = 255
    G = 150
    B = 0

! suppress all the values that are below the threshold
  else if (abs(image_data_2D(ix,iy)) <= max_amplitude * cutvect) then

! use a black or white background for points that are below the threshold
    if (WHITE_BACKGROUND) then
      R = 255
      G = 255
      B = 255
    else
      R = 0
      G = 0
      B = 0
    endif

! represent regular image points using red if value is positive, blue if negative
  else if (normalized_value >= 0.d0) then
    R = nint(255.d0*normalized_value**POWER_DISPLAY)
    G = 0
    B = 0
  else
    R = 0
    G = 0
    B = nint(255.d0*abs(normalized_value)**POWER_DISPLAY)
  endif

! write color pixel
  write(27,"(i3,' ',i3,' ',i3)") R,G,B

  enddo
enddo

! close file
close(27)

! call the system to convert image to Gif (can be commented out if "call system" is missing in your compiler)
! call system(system_command)

end subroutine create_color_image


end module seismicFDTD2d