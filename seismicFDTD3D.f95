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

subroutine doall(im, mlist, nx, ny, nz, dx, dy, dz, npoints_pml, & 
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
integer :: nx, ny, nz, nstep, npoints_pml, it_display
integer,dimension(:,:) :: im
! integer,dimension(:,:) :: rcx
integer,dimension(:) :: src 
real(kind=dp),dimension(:,:) :: mlist
real(kind=dp) :: f0
real(kind=dp),dimension(2) :: angle_force
real(kind=dp) :: dx, dy, dz
real(kind=dp),dimension(nx+2*npoints_pml,ny+2*npoints_pml,nz+2*npoints_pml) :: &
  c11, c12, c13, c22, c23, c33, c44, c55, c66, rho
! character(len=6) :: src_type

!f2py3 intent(in) :: im, mlist, nx, ny, dx, dy, npoints_pml, src
!f2py3 intent(in) :: f0, nstep, it_display, angle_force
!f2py3 intent(hide), depend(im) :: nx = shape(im, 0), ny = shape(im,1)

! Preallocate arrays
c11(:,:,:) = 0.d0
c12(:,:,:) = 0.d0
c13(:,:,:) = 0.d0
c22(:,:,:) = 0.d0
c23(:,:,:) = 0.d0
c33(:,:,:) = 0.d0
c44(:,:,:) = 0.d0
c55(:,:,:) = 0.d0
c66(:,:,:) = 0.d0
rho(:,:,:) = 0.d0


! Setup arrays
call stiffness_arrays(im, mlist, c11, c12, c13, c22, c23, c33, c44, c55, c66, rho, npoints_pml)

call seismic_cpml_2d(nx+2*npoints_pml, ny+2*npoints_pml, nz+2*npoints_pml, &
                      c11, c12, c13, c22, c23, c33, c44, c55, c66, &
                      rho, dx, dy, dz, npoints_pml, src, f0, nstep, &
                      it_display, angle_force)


end subroutine doall


!==============================================================================
subroutine stiffness_arrays(im, mlist, c11, c12, c13, c22, c23, c33, c44, &
                            c55, c66, rho, npoints_pml) 
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
integer :: m, n, p, i, j, k, npoints_pml
real(kind=dp), dimension(:,:) :: mlist
real(kind=dp), dimension(:,:,:) :: c11, c12, c13, c22, c23, c33, c44, c55, c66, rho

m = size(rho, 1)
n = size(rho, 2)
p = size(rho, 3)

do i=npoints_pml+1, m+npoints_pml
  do j=npoints_pml+1, n+npoints_pml
    do k=npoints_pml+1, p+npoints_pml 
      c11(i,j,k) = mlist(im(i-npoints_pml,j-npoints_pml), 2)
      c12(i,j,k) = mlist(im(i-npoints_pml,j-npoints_pml), 3) 
      c22(i,j,k) = mlist(im(i-npoints_pml,j-npoints_pml), 5)
      c66(i,j,k) = mlist(im(i-npoints_pml,j-npoints_pml), 10)
      rho(i,j,k) = mlist(im(i-npoints_pml,j-npoints_pml), 11) 
    enddo
  end do
end do


! Extend the boundary values of the stiffnesses into the PML region
do i = 1,npoints_pml+1
  ! top and bottom
  c11( i, :, : ) = c11(npoints_pml+1,:, :)
  c11( m+npoints_pml-1+i, :, : ) = c11(m+npoints_pml,:, :)

  c12( i, :, : ) = c12(npoints_pml+1,:, :)
  c12( m+npoints_pml-1+i, : , :) = c12(m+npoints_pml,:, :)
  
  c22( i, :, : ) = c22(npoints_pml+1,:, :)
  c22( m+npoints_pml-1+i, : , :) = c22(m+npoints_pml,:, :)
  
  c66( i, :, : ) = c66(npoints_pml+1,:, :)
  c66( m+npoints_pml-1+i, : , :) = c66(m+npoints_pml,:, :)
  
  rho( i, :, : ) = rho(npoints_pml+1,:, :)
  rho( m+npoints_pml-1+i, : , :) = rho(m+npoints_pml,:, :)

  ! ! left and right
  c11( :, i, : ) = c11(:, npoints_pml+1, :)
  c11( :, n+npoints_pml-1+i , :) = c11(:,n+npoints_pml, :)
  
  c12( :, i, : ) = c12(:, npoints_pml+1, :)
  c12( :, n+npoints_pml-1+i , :) = c12(:,n+npoints_pml, :)
  
  c22( :, i, : ) = c22(:, npoints_pml+1, :)
  c22( :, n+npoints_pml-1+i , :) = c22(:,n+npoints_pml, :)
  
  c66( :, i, : ) = c66(:, npoints_pml+1, :)
  c66( :, n+npoints_pml-1+i , :) = c66(:,n+npoints_pml, :)

  rho( :, i, : ) = rho(:, npoints_pml+1, :)
  rho( :, n+npoints_pml-1+i , :) = rho(:,n+npoints_pml, :)

end do 


end subroutine stiffness_arrays

!==============================================================================


subroutine seismic_cpml_2d(nx, ny, nz, &
                      c11, c12, c13, c22, c23, c33, c44, c55, c66,  &
                      rho, dx, dy, dz, npoints_pml, src, f0, nstep, &
                      it_display, angle_force)

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
integer :: nx, ny, nz

! thickness of the PML layer in grid points
integer :: npoints_pml
real(kind=dp), dimension(nx,ny,nz) :: c11, c12, c13, c22, c23, c33, c44, c55, c66, rho
real(kind=dp) :: f0

! total number of time steps
integer :: nstep

! time step in seconds. decreasing the time step improves the pml attenuation
! but it should be inversely proportional to the center frequency of the 
! source frequency 
real(kind=dp) :: DT
real(kind=dp) :: dx, dy, dz
! parameters for the source
real(kind=dp) :: t0
real(kind=dp), parameter :: factor = 1.d7

! source
integer,dimension(:) :: src
integer :: isource, jsource, ksource
! angle of source force clockwise with respect to vertical (Y) axis
real(kind=dp), dimension(2) :: ANGLE_FORCE !angle_force = (/ theta, phi /)


! display information on the screen from time to time
integer :: IT_DISPLAY

! value of PI
real(kind=dp), parameter :: PI = 3.141592653589793238462643d0

! conversion from degrees to radians
real(kind=dp), parameter :: DEGREES_TO_RADIANS = PI / 180.d0

! large value for maximum
real(kind=dp), parameter :: HUGEVAL = 1.d+30

! velocity threshold above which we consider that the code became unstable
real(kind=dp), parameter :: STABILITY_THRESHOLD = 1.d+25

! main arrays
real(kind=dp), dimension(nx,ny,nz) :: vx,vy,vz,sigmaxx,sigmayy,sigmazz,sigmaxy,sigmaxz,sigmayz

! power to compute d0 profile. Increasing this value allows for a larger dampening gradient in the PML
real(kind=dp), parameter :: NPOWER = 2.d0

! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-11
real(kind=dp), parameter :: K_MAX_PML = 1.d1
real(kind=dp) :: ALPHA_MAX_PML

! arrays for the memory variables
! could declare these arrays in PML only to save a lot of memory, but proof of concept only here
real(kind=dp), dimension(nx,ny,nz) :: &
    memory_dvx_dx, memory_dvx_dy, memory_dvx_dz, &
    memory_dvy_dx, memory_dvy_dy, memory_dvy_dz, &
    memory_dvz_dx, memory_dvz_dy, memory_dvz_dz, &
    memory_dsigmaxx_dx, memory_dsigmayy_dy, memory_dsigmazz_dz, &
    memory_dsigmaxy_dx, memory_dsigmaxy_dy, &
    memory_dsigmaxz_dx, memory_dsigmaxz_dz, &
    memory_dsigmayz_dy, memory_dsigmayz_dz


real(kind=dp) :: &
    value_dvx_dx, value_dvx_dy, value_dvx_dz, &
    value_dvy_dx, value_dvy_dy, value_dvy_dz, &
    value_dvz_dx, value_dvz_dy, value_dvz_dz, &
    value_dsigmaxx_dx, value_dsigmayy_dy, value_dsigmazz_dz, &
    value_dsigmaxy_dx, value_dsigmaxy_dy, &
    value_dsigmaxz_dx, value_dsigmaxz_dz, &
    value_dsigmayz_dy, value_dsigmayz_dz

! 1D arrays for the damping profiles
real(kind=dp), dimension(nx) :: d_x,K_x,alpha_x,a_x,b_x,d_x_half,K_x_half,alpha_x_half,a_x_half,b_x_half
real(kind=dp), dimension(ny) :: d_y,K_y,alpha_y,a_y,b_y,d_y_half,K_y_half,alpha_y_half,a_y_half,b_y_half
real(kind=dp), dimension(nz) :: d_z,K_z,alpha_z,a_z,b_z,d_z_half,K_z_half,alpha_z_half,a_z_half,b_z_half

real(kind=dp) :: thickness_PML_x,thickness_PML_y,thickness_PML_z
real(kind=dp) :: xoriginleft,xoriginright,yoriginout, yoriginin, zoriginbottom,zorigintop
real(kind=dp) :: Rcoef,d0_x,d0_y,d0_z,xval,yval,zval,abscissa_in_PML,abscissa_normalized

! for the source
real(kind=dp) :: a,t,force_x,force_y,force_z,source_term

integer :: i,j,k,it

real(kind=dp) :: velocnorm

! for stability estimate
real(kind=dp) :: quasi_cp_max


! ------------------------ Assign some constants -----------------------

isource = src(1)+npoints_pml
jsource = src(2)+npoints_pml

t0 = 1.0d0/f0
DT = minval( (/dx,dy,dz/) )/ ( 2.0* sqrt( ( maxval( (/ c11/rho, c22/rho, c66/rho /) ) ) ) )
ALPHA_MAX_PML = PI*f0 ! from Festa and Vilotte
a = pi*pi*f0*f0
angle_force = angle_force * degrees_to_radians
! ----------------------------------------------------------------------
!---
!--- program starts here
!---

! reflection coefficient (INRIA report section 6.1) http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
Rcoef = 0.001d0

! check that NPOWER is okay
  if (NPOWER < 1) stop 'NPOWER must be greater than 1'

! ==================================== PML ====================================
! ---------------- define profile of absorption in PML region -----------------

! thickness of the PML layer in meters
thickness_PML_x = NPOINTS_PML * dx
thickness_PML_y = NPOINTS_PML * dy
thickness_PML_z = npoints_pml * dz


! Initialize PML 
  d_x(:) = 0.d0
  d_x_half(:) = 0.d0
  K_x(:) = 1.d0
  K_x_half(:) = 1.d0
  alpha_x(:) = 0.d0
  alpha_x_half(:) = 0.d0
  a_x(:) = 0.d0
  a_x_half(:) = 0.d0

  d_y(:) = 0.d0
  d_y_half(:) = 0.d0
  K_y(:) = 1.d0
  K_y_half(:) = 1.d0
  alpha_y(:) = 0.d0
  alpha_y_half(:) = 0.d0
  a_y(:) = 0.d0
  a_y_half(:) = 0.d0

  d_z(:) = 0.d0
  d_z_half(:) = 0.d0
  K_z(:) = 1.d0
  K_z_half(:) = 1.d0 
  alpha_z(:) = 0.d0
  alpha_z_half(:) = 0.d0
  a_x(:) = 0.d0
  a_x_half(:) = 0.d0

! ------------------------- damping in the X direction-------------------------

! origin of the PML layer (position of right edge minus thickness, in meters)
  xoriginleft = dble(thickness_PML_x)
  xoriginright = dx * dble(NX-1) - thickness_PML_x

do i = 1,NX
  ! to compute d0 below, and for stability estimate
  quasi_cp_max = max( &
    sqrt(c22(i,npoints_pml,npoints_pml)/rho(i,npoints_pml,npoints_pml) ),&
    sqrt(c11(i,npoints_pml,npoints_pml)/rho(i,npoints_pml,npoints_pml) ), &
    sqrt(c33(i,npoints_pml,npoints_pml)/rho(i,npoints_pml,npoints_pml) ) )

  ! compute d0 from INRIA report section 6.1 http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
  d0_x = - dble(NPOWER + 1) * quasi_cp_max * log(Rcoef) / (2.d0 * thickness_PML_x)


    ! abscissa of current grid point along the damping profile
    xval = DX * dble(i-1)

    !---------- left edge
    ! define damping profile at the grid points
    abscissa_in_PML = xoriginleft - xval
    if (abscissa_in_PML >= 0.d0) then
      abscissa_normalized = abscissa_in_PML / thickness_PML_x
      d_x(i) = d0_x * abscissa_normalized**NPOWER
      ! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
      K_x(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
      alpha_x(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
    endif

  ! define damping profile at half the grid points
    abscissa_in_PML = xoriginleft - (xval + DX/2.d0)
    if (abscissa_in_PML >= 0.d0) then
      abscissa_normalized = abscissa_in_PML / thickness_PML_x
      d_x_half(i) = d0_x * abscissa_normalized**NPOWER
      ! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
      K_x_half(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
      alpha_x_half(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
    endif

!---------- right edge
      ! define damping profile at the grid points
      abscissa_in_PML = xval - xoriginright
      if (abscissa_in_PML >= 0.d0) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_x
        d_x(i) = d0_x * abscissa_normalized**NPOWER
        K_x(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_x(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
      endif

      ! define damping profile at half the grid points
      abscissa_in_PML = xval + DX/2.d0 - xoriginright
      if (abscissa_in_PML >= 0.d0) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_x
        d_x_half(i) = d0_x * abscissa_normalized**NPOWER
        K_x_half(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_x_half(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
      endif

    ! just in case, for -5 at the end
    if (alpha_x(i) < 0.d0) alpha_x(i) = 0.d0
    if (alpha_x_half(i) < 0.d0) alpha_x_half(i) = 0.d0

    b_x(i) = exp(- (d_x(i) / K_x(i) + alpha_x(i)) * DT)
    b_x_half(i) = exp(- (d_x_half(i) / K_x_half(i) + alpha_x_half(i)) * DT)

    ! this to avoid division by z'ero outside the PML
    if (abs(d_x(i)) > 1.d-6) a_x(i) = d_x(i) * (b_x(i) - 1.d0) / (K_x(i) * (d_x(i) + K_x(i) * alpha_x(i)))
    if (abs(d_x_half(i)) > 1.d-6) a_x_half(i) = d_x_half(i) * &
      (b_x_half(i) - 1.d0) / (K_x_half(i) * (d_x_half(i) + K_x_half(i) * alpha_x_half(i)))

  enddo

! ------------------------ damping in the Y direction -------------------------

! origin of the PML layer (position of right edge minus thickness, in meters)
yoriginout = dble(thickness_PML_y)
yoriginin = dy * dble(ny-1) - thickness_PML_y

do j = 1,NY
  ! to compute d0 below, and for stability estimate
  quasi_cp_max = max( &
    sqrt(c22(npoints_pml,j,npoints_pml)/rho(npoints_pml,j,npoints_pml) ),&
    sqrt(c11(npoints_pml,j,npoints_pml)/rho(npoints_pml,j,npoints_pml) ), &
    sqrt(c33(npoints_pml,j,npoints_pml)/rho(npoints_pml,j,npoints_pml) ) )
  ! compute d0 from INRIA report section 6.1 http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
  d0_y = - dble(NPOWER + 1) * quasi_cp_max * log(Rcoef) / (2.d0 * thickness_PML_y)

  ! abscissa of current grid point along the damping profile
  yval = DY * dble(j-1)

  !---------- bottom edge
    ! define damping profile at the grid points
    abscissa_in_PML = yoriginout - yval
    if (abscissa_in_PML >= 0.d0) then
      abscissa_normalized = abscissa_in_PML / thickness_PML_y
      d_y(j) = d0_y * abscissa_normalized**NPOWER
      K_y(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
      alpha_y(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
    endif

    ! define damping profile at half the grid points
    abscissa_in_PML = yoriginout - (yval + DY/2.d0)
    if (abscissa_in_PML >= 0.d0) then
      abscissa_normalized = abscissa_in_PML / thickness_PML_y
      d_y_half(j) = d0_y * abscissa_normalized**NPOWER
      K_y_half(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
      alpha_y_half(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
    endif


  !---------- top edge
    ! define damping profile at the grid points
    abscissa_in_PML = yval - yoriginin
    if (abscissa_in_PML >= 0.d0) then
      abscissa_normalized = abscissa_in_PML / thickness_PML_y
      d_y(j) = d0_y * abscissa_normalized**NPOWER
      K_y(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
      alpha_y(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
    endif

    ! define damping profile at half the grid points
    abscissa_in_PML = yval + DY/2.d0 - yoriginin
    if (abscissa_in_PML >= 0.d0) then
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

! ------------------------- damping in the z direction ------------------------

! origin of the PML layer (position of right edge minus thickness, in meters)
zoriginbottom = dble(thickness_PML_z)
zorigintop = dz * dble(nz-1) - thickness_PML_z

do j = 1,nz
  ! to compute d0 below, and for stability estimate
  quasi_cp_max = max( &
    sqrt( c22(npoints_pml,npoints_pml,j) / rho(npoints_pml,npoints_pml,j) ), &
    sqrt( c11(npoints_pml,npoints_pml,j) / rho(npoints_pml,npoints_pml,j) ), &
    sqrt( c33(npoints_pml,npoints_pml,j) / rho(npoints_pml,npoints_pml,j) ) )
  ! compute d0 from INRIA report section 6.1 http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
  d0_x = - dble(NPOWER + 1) * quasi_cp_max * log(Rcoef) / (2.d0 * thickness_PML_y)

  ! abscissa of current grid point along the damping profile
  zval = dz * dble(j-1)

  !---------- bottom edge
    ! define damping profile at the grid points
    abscissa_in_PML = zoriginbottom - zval
    if (abscissa_in_PML >= 0.d0) then
      abscissa_normalized = abscissa_in_PML / thickness_PML_z
      d_z(j) = d0_z * abscissa_normalized**NPOWER
      K_z(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
      alpha_z(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
    endif

    ! define damping profile at half the grid points
    abscissa_in_PML = zoriginbottom - (zval + dz/2.d0)
    if (abscissa_in_PML >= 0.d0) then
      abscissa_normalized = abscissa_in_PML / thickness_PML_z
      d_z_half(j) = d0_z * abscissa_normalized**NPOWER
      K_z_half(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
      alpha_z_half(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
    endif


  !---------- top edge
    ! define damping profile at the grid points
    abscissa_in_PML = zval - zorigintop
    if (abscissa_in_PML >= 0.d0) then
      abscissa_normalized = abscissa_in_PML / thickness_PML_z
      d_z(j) = d0_z * abscissa_normalized**NPOWER
      K_z(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
      alpha_z(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
    endif

    ! define damping profile at half the grid points
    abscissa_in_PML = zval + dz/2.d0 - zorigintop
    if (abscissa_in_PML >= 0.d0) then
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

! =============================== Forward Model ===============================
! initialize arrays
vx(:,:,:) = 0.d0
vy(:,:,:) = 0.d0
vz(:,:,:) = 0.d0

sigmaxx(:,:,:) = 0.d0
sigmayy(:,:,:) = 0.d0
sigmazz(:,:,:) = 0.d0
sigmaxy(:,:,:) = 0.d0
sigmaxz(:,:,:) = 0.d0
sigmayz(:,:,:) = 0.d0

! PML
memory_dvx_dx(:,:,:) = 0.d0
memory_dvx_dy(:,:,:) = 0.d0
memory_dvx_dz(:,:,:) = 0.d0

memory_dvy_dx(:,:,:) = 0.d0
memory_dvy_dy(:,:,:) = 0.d0
memory_dvy_dz(:,:,:) = 0.d0

memory_dvz_dx(:,:,:) = 0.d0
memory_dvz_dy(:,:,:) = 0.d0 
memory_dvz_dz(:,:,:) = 0.d0

memory_dsigmaxx_dx(:,:,:) = 0.d0
memory_dsigmayy_dy(:,:,:) = 0.d0
memory_dsigmazz_dz(:,:,:) = 0.d0

memory_dsigmaxy_dx(:,:,:) = 0.d0
memory_dsigmaxy_dy(:,:,:) = 0.d0
memory_dsigmaxz_dx(:,:,:) = 0.d0
memory_dsigmaxz_dz(:,:,:) = 0.d0
memory_dsigmayz_dy(:,:,:) = 0.d0
memory_dsigmayz_dz(:,:,:) = 0.d0


!---
!---  beginning of time loop
!---

do it = 1,NSTEP
  !------------------------------------------------------------
  ! compute stress sigma and update memory variables for C-PML
  !------------------------------------------------------------
  ! Update in the x direction
  do k = 2,nz
    do j = 2,NY
      do i = 1,NX-1

        value_dvx_dx = (vx(i+1,j,k) - vx(i,j,k) ) / DX
        value_dvy_dy = (vy(i,j,k) -vy(i,j-1,k) ) / dy
        value_dvz_dz = (vz(i,j,k) - vz(i,j-1,k) ) / dz

        memory_dvx_dx(i,j,k) = b_x_half(i) * memory_dvx_dx(i,j,k) + a_x_half(i) * value_dvx_dx
        memory_dvy_dy(i,j,k) = b_y(j) * memory_dvy_dy(i,j,k) + a_y(j) * value_dvy_dy
        memory_dvz_dz(i,j,k) = b_z(k) * memory_dvz_dz(i,j,k) + a_z(k) * value_dvz_dz

        value_dvx_dx = value_dvx_dx / K_x_half(i) + memory_dvx_dx(i,j,k)
        value_dvy_dy = value_dvy_dy / K_y(j) + memory_dvy_dy(i,j,k)
        value_dvz_dz = value_dvz_dz / K_z(k) + memory_dvz_dz(i,j,k)

        sigmaxx(i,j,k) = sigmaxx(i,j,k) + ( c11(i,j,k) * value_dvx_dx + c12(i,j,k) * value_dvy_dy + c13(i,j,k) * value_dvz_dz) * dt
        sigmayy(i,j,k) = sigmayy(i,j,k) + ( c12(i,j,k) * value_dvx_dx + c22(i,j,k) * value_dvy_dy + c23(i,j,k) * value_dvz_dz) * dt
        sigmazz(i,j,k) = sigmazz(i,j,k) + ( c13(i,j,k) * value_dvx_dx + c23(i,j,k) * value_dvy_dy + c33(i,j,k) * value_dvz_dz) * dt

      enddo
    enddo
  enddo

  ! Update sigmaxy, x-direction is full nodes
  do k = 1,nz
    do j = 1,ny-1
      do i = 2,NX

        value_dvy_dx = (vy(i,j,k) - vy(i-1,j,k)) / DX
        value_dvx_dy = (vx(i,j+1,k) - vx(i,j,k)) / DY

        memory_dvy_dx(i,j,k) = b_x(i) * memory_dvy_dx(i,j,k) + a_x(i) * value_dvy_dx
        memory_dvx_dy(i,j,k) = b_y_half(j) * memory_dvx_dy(i,j,k) + a_y_half(j) * value_dvx_dy

        value_dvy_dx = value_dvy_dx / K_x(i) + memory_dvy_dx(i,j,k)
        value_dvx_dy = value_dvx_dy / K_y_half(j) + memory_dvx_dy(i,j,k)

        sigmaxy(i,j,k) = sigmaxy(i,j,k) + c66(i,j,k) * (value_dvy_dx + value_dvx_dy) * DT

      enddo
    enddo
  enddo

  ! Update sigmaxz, z-direction is full nodes
  do k = 1,nz-1
    do j = 1,ny
      do i = 2,nx

        value_dvz_dx = (vz(i,j,k) - vz(i-1,j,k) ) / dx
        value_dvx_dz = (vx(i,j,k+1) - vx(i,j,k) ) / dz

        memory_dvz_dx(i,j,k) = b_x_half(i) * memory_dvz_dx(i,j,k) + a_x_half(i) * value_dvz_dx
        memory_dvx_dz(i,j,k) = b_z(k) * memory_dvx_dz(i,j,k) + a_x(k) * value_dvx_dz

        value_dvz_dx = value_dvz_dx / K_x_half(i) + memory_dvz_dx(i,j,k) 
        value_dvx_dz = value_dvx_dz / K_z(k) + memory_dvx_dz(i,j,k)

        sigmaxz(i,j,k) = sigmaxz(i,j,k) + c55(i,j,k) * ( value_dvx_dz + value_dvz_dx) * dt 

      enddo
    enddo
  enddo

  ! update sigmayz, y-direction is full nodes
  do k = 1,nz-1
    do j = 1,ny-1
      do i = 1,nx

        value_dvy_dz = (vy(i,j+1,k) - vy(i,j,k) ) / dz
        value_dvz_dy = (vz(i,j,k+1) - vz(i,j,k) ) / dy

        memory_dvy_dz(i,j,k) = b_z_half(k) * memory_dvy_dz(i,j,k) + a_z_half(k) * value_dvy_dz 
        memory_dvz_dy(i,j,k) = b_y(j) * memory_dvz_dy(i,j,k) + a_y(j) * value_dvz_dy

        value_dvy_dz = value_dvy_dz / K_z_half(k) + memory_dvy_dz(i,j,k)
        value_dvz_dy = value_dvz_dy / K_y(j) + memory_dvz_dy(i,j,k)

        sigmayz(i,j,k) = sigmayz(i,j,k)  + c44(i,j,k) * ( value_dvy_dz + value_dvz_dy) * dt 

      enddo
    enddo
  enddo


!--------------------------------------------------------
! compute velocity and update memory variables for C-PML
!--------------------------------------------------------
  do k = 2,nz
    do j = 2,NY
      do i = 2,NX
        ! ds1/dx, ds6/dy, ds5,dz
        value_dsigmaxx_dx = (sigmaxx(i,j,k) - sigmaxx(i-1,j,k) ) / dx
        value_dsigmaxy_dy = (sigmaxy(i,j,k) - sigmaxy(i,j-1,k) ) / dy
        value_dsigmaxz_dz = (sigmaxz(i,j,k) - sigmaxz(i,j,k-1) ) /dz

        memory_dsigmaxx_dx(i,j,k) = b_x(i) * memory_dsigmaxx_dx(i,j,k) + a_x(i) * value_dsigmaxx_dx
        memory_dsigmaxy_dy(i,j,k) = b_y(j) * memory_dsigmaxy_dy(i,j,k) + a_y(j) * value_dsigmaxy_dy
        memory_dsigmaxz_dz(i,j,k) = b_z(k) * memory_dsigmaxz_dz(i,j,k) + a_z(k) * value_dsigmaxz_dz

        value_dsigmaxx_dx = value_dsigmaxx_dx / K_x(i) + memory_dsigmaxx_dx(i,j,k)
        value_dsigmaxy_dy = value_dsigmaxy_dy / K_y(j) + memory_dsigmaxy_dy(i,j,k)
        value_dsigmaxz_dz = value_dsigmaxz_dz / K_z(k) + memory_dsigmaxz_dz(i,j,k) 

        vx(i,j,k) = vx(i,j,k) + (value_dsigmaxx_dx + value_dsigmaxy_dy + value_dsigmaxz_dz) * dt / rho(i,j,k)

      enddo
    enddo
  enddo

  do k = 2,nz
    do j = 1,ny-1
      do i = 1,nx-1
        ! ds6/dx, ds2/dy, ds4/dz
        value_dsigmaxy_dx = ( sigmaxy(i+1,j,k) - sigmaxy(i,j,k) ) / dx
        value_dsigmayy_dy = ( sigmayy(i,j+1,k) - sigmayy(i,j,k) ) / dy
        value_dsigmayz_dz = ( sigmayz(i,j,k) - sigmayz(i,j,k-1) ) / dz

        memory_dsigmaxy_dx(i,j,k) = b_x_half(i) * memory_dsigmaxy_dx(i,j,k) + a_x_half(i) * value_dsigmaxy_dx
        memory_dsigmayy_dy(i,j,k) = b_y_half(j) * memory_dsigmayy_dy(i,j,k) + a_y_half(j) * value_dsigmayy_dy
        memory_dsigmayz_dz(i,j,k) = b_z(k) * memory_dsigmayz_dz(i,j,k) + a_z(k) * value_dsigmayz_dz

        value_dsigmaxy_dx = value_dsigmaxy_dx / K_x_half(i) + memory_dsigmaxy_dx(i,j,k)
        value_dsigmayy_dy = value_dsigmayy_dy / K_y_half(j) + memory_dsigmayy_dy(i,j,k)
        value_dsigmayz_dz = value_dsigmayz_dz / K_z(k) + memory_dsigmayz_dz(i,j,k)

        vy(i,j,k) = vy(i,j,k) + (value_dsigmaxy_dx + value_dsigmayy_dy + value_dsigmayz_dz) * dt / rho(i,j,k)

      enddo
    enddo
  enddo

  do k = 1,nz-1
    do j = 2,ny
      do i = 1,nx-1

        ! ds5/dx, ds4/dy, ds3/dz
        value_dsigmaxz_dx = ( sigmaxz(i+1,j,k) - sigmaxz(i,j,k) ) / dx
        value_dsigmayz_dy = ( sigmayz(i,j,k) - sigmayz(i,j-1,k) ) / dy
        value_dsigmazz_dz = ( sigmazz(i,j,k+1) - sigmazz(i,j,k) ) / dz

        memory_dsigmaxz_dx(i,j,k) = b_x_half(i) * memory_dsigmaxz_dx(i,j,k) + a_x_half(i) * value_dsigmaxz_dx
        memory_dsigmayz_dy(i,j,k) = b_y(j) * memory_dsigmayz_dy(i,j,k) + a_y(j) * value_dsigmayz_dy
        memory_dsigmazz_dz(i,j,k) = b_z_half(k) * memory_dsigmazz_dz(i,j,k) + a_z_half(k) * value_dsigmazz_dz

        value_dsigmaxz_dx = value_dsigmaxz_dx / K_x_half(i) + memory_dsigmaxz_dx(i,j,k)
        value_dsigmayz_dy = value_dsigmayz_dy / K_y(j) + memory_dsigmayz_dy(i,j,k)
        value_dsigmazz_dz = value_dsigmazz_dz / K_z_half(k) + memory_dsigmazz_dz(i,j,k)

        vz(i,j,k) = vz(i,j,k) + (value_dsigmaxz_dx + value_dsigmayz_dy + value_dsigmazz_dz) * dt / rho(i,j,k)

      enddo
    enddo
  enddo

  ! add the source (force vector located at a given grid point)
  t = dble(it-1)*DT

  ! Gaussian
  source_term = factor * 2.d0*exp(-a*(t-t0)**2)

  ! Use spherical coordinates for the source rotation
  force_x = sin( angle_force(1) ) * cos( angle_force(2) ) * source_term
  force_z = sin( angle_force(1) ) * sin( angle_force(2) ) * source_term
  force_y = cos( angle_force(1) ) * source_term

  vx(isource,jsource,ksource) = vx(isource,jsource,ksource) + force_x * DT / rho(i,j,k)
  vy(isource,jsource,ksource) = vy(isource,jsource,ksource) + force_y * DT / rho(i,j,k)

  ! Dirichlet conditions (rigid boundaries) on the edges or at the bottom of the PML layers
  vx(1,:,:) = 0.d0
  vx(:,1,:) = 0.d0
  vx(:,:,1) = 0.d0
  vx(NX,:,:) = 0.d0
  vx(:,NY,:) = 0.d0
  vx(:,:,NZ) = 0.d0

  vy(1,:,:) = 0.d0
  vy(:,1,:) = 0.d0
  vy(:,:,1) = 0.d0
  vy(NX,:,:) = 0.d0
  vy(:,NY,:) = 0.d0
  vy(:,:,NZ) = 0.d0
  
  vz(1,:,:) = 0.d0
  vz(:,1,:) = 0.d0
  vz(:,:,1) = 0.d0
  vz(NX,:,:) = 0.d0
  vz(:,NY,:) = 0.d0
  vz(:,:,NZ) = 0.d0

  
  ! output information
  if (mod(it,IT_DISPLAY) == 0 .or. it == 1) then

  ! print maximum of norm of velocity
  velocnorm = maxval( sqrt(vx**2 + vy**2 + vy**2) )
  print *,'Time step # ',it,' out of ',NSTEP
  print *,'Time: ',sngl((it-1)*DT),' seconds'

  if (velocnorm > STABILITY_THRESHOLD) stop 'code became unstable and blew up'
    call write_image(vx, nx, ny, nz, it, 'Vx')
    call write_image(vy, nx, ny, nz, it, 'Vy')
    call write_image(vz, nx, ny, nz, it, 'Vz')
  endif

enddo   ! end of time loop

end subroutine seismic_cpml_2d


!==============================================================================
subroutine write_image(image_data, nx, ny, nz, it, channel)

implicit none

integer, parameter :: dp = kind(0.d0)
integer :: nx, ny, nz, it
real(kind=dp) :: image_data(nx, ny, nz)
character(len=2) :: channel
character(len=100) :: filename

WRITE (filename, "(a2, i6.6, '.dat')" ) channel, it

open(unit = 10, form = 'unformatted', file = trim(filename) )
write(10) sngl(image_data)

close(unit = 10)

end subroutine write_image



end module seismicFDTD2d