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


module seismicFDTD25d

implicit none

contains


  !==============================================================================
  subroutine stiffness_write(im, mlist, npoints_pml, nx, nz) 
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
    integer :: nx, nz
    integer,dimension(nx,nz) :: im
    integer :: i, j, npoints_pml
    real(kind=dp), dimension(:,:) :: mlist
    real(kind=dp), dimension(2*npoints_pml+nx,2*npoints_pml+nz) :: c11, c12, c13, &
                                                                c22, c23, c33, &
                                                                c44, c55, c66, rho

    !f2py3 intent(in):: im, mlist, npoints_pml, nx, nz

    c11(:,:) = 0.d0 
    c12(:,:) = 0.d0 
    c13(:,:) = 0.d0
    c22(:,:) = 0.d0 
    c23(:,:) = 0.d0 
    c33(:,:) = 0.d0 
    c44(:,:) = 0.d0 
    c55(:,:) = 0.d0 
    c66(:,:) = 0.d0 
    rho(:,:) = 0.d0 

    !Assign between the PML regions
    do i=npoints_pml+1, nx+npoints_pml
      do j=npoints_pml+1, nz+npoints_pml
        c11(i,j) = mlist(im(i-npoints_pml,j-npoints_pml), 2)
        c12(i,j) = mlist(im(i-npoints_pml,j-npoints_pml), 3)
        c13(i,j) = mlist(im(i-npoints_pml,j-npoints_pml), 4) 
        c22(i,j) = mlist(im(i-npoints_pml,j-npoints_pml), 5)
        c23(i,j) = mlist(im(i-npoints_pml,j-npoints_pml), 6)
        c33(i,j) = mlist(im(i-npoints_pml,j-npoints_pml), 7)
        c44(i,j) = mlist(im(i-npoints_pml,j-npoints_pml), 8)
        c55(i,j) = mlist(im(i-npoints_pml,j-npoints_pml), 9)
        c66(i,j) = mlist(im(i-npoints_pml,j-npoints_pml), 10)
        rho(i,j) = mlist(im(i-npoints_pml,j-npoints_pml), 11) 
      enddo
    enddo

    ! Extend the boundary values of the stiffnesses into the PML region
    do i = 1,npoints_pml+1
      ! top 
      c11( i, :) = c11(npoints_pml+1,:)
      c12( i, :) = c12(npoints_pml+1,:)
      c13( i, :) = c13(npoints_pml+1,:)
      c22( i, :) = c22(npoints_pml+1,:)
      c23( i, :) = c23(npoints_pml+1,:)
      c33( i, :) = c33(npoints_pml+1,:)
      c44( i, :) = c44(npoints_pml+1,:)
      c55( i, :) = c55(npoints_pml+1,:)
      c66( i, :) = c66(npoints_pml+1,:)
      rho( i, :) = rho(npoints_pml+1,:)

      ! bottom
      c11( nx+npoints_pml-1+i, :) = c11(nx+npoints_pml,:)
      c12( nx+npoints_pml-1+i, :) = c12(nx+npoints_pml,:)
      c13( nx+npoints_pml-1+i, :) = c13(nx+npoints_pml,:)
      c22( nx+npoints_pml-1+i, :) = c22(nx+npoints_pml,:)
      c23( nx+npoints_pml-1+i, :) = c23(nx+npoints_pml,:)
      c33( nx+npoints_pml-1+i, :) = c33(nx+npoints_pml,:)
      c44( nx+npoints_pml-1+i, :) = c44(nx+npoints_pml,:)
      c55( nx+npoints_pml-1+i, :) = c55(nx+npoints_pml,:)
      c66( nx+npoints_pml-1+i, :) = c66(nx+npoints_pml,:)
      rho( nx+npoints_pml-1+i, :) = rho(nx+npoints_pml,:)

      ! left 
      c11( :, i) = c11(:, npoints_pml+1)
      c12( :, i) = c12(:, npoints_pml+1)
      c13( :, i) = c13(:, npoints_pml+1)
      c22( :, i) = c22(:, npoints_pml+1)
      c23( :, i) = c23(:, npoints_pml+1)
      c33( :, i) = c33(:, npoints_pml+1)
      c44( :, i) = c44(:, npoints_pml+1)
      c55( :, i) = c55(:, npoints_pml+1)
      c66( :, i) = c66(:, npoints_pml+1)
      rho( :, i) = rho(:, npoints_pml+1)

      ! right
      c11( :, nz+npoints_pml-1+i) = c11(:,nz+npoints_pml)
      c12( :, nz+npoints_pml-1+i) = c12(:,nz+npoints_pml)
      c13( :, nz+npoints_pml-1+i) = c12(:,nz+npoints_pml)      
      c22( :, nz+npoints_pml-1+i) = c22(:,nz+npoints_pml)
      c23( :, nz+npoints_pml-1+i) = c23(:,nz+npoints_pml)
      c33( :, nz+npoints_pml-1+i) = c33(:,nz+npoints_pml)
      c44( :, nz+npoints_pml-1+i) = c44(:,nz+npoints_pml)
      c55( :, nz+npoints_pml-1+i) = c55(:,nz+npoints_pml)      
      c66( :, nz+npoints_pml-1+i) = c66(:,nz+npoints_pml)
      rho( :, nz+npoints_pml-1+i) = rho(:,nz+npoints_pml)

    end do 

    ! Write each of the matrices to file
    call material_rw('c11.dat', c11, .FALSE.)
    call material_rw('c12.dat', c12, .FALSE.)
    call material_rw('c13.dat', c13, .FALSE.)
    call material_rw('c22.dat', c22, .FALSE.)
    call material_rw('c23.dat', c23, .FALSE.)
    call material_rw('c33.dat', c33, .FALSE.)
    call material_rw('c44.dat', c44, .FALSE.)
    call material_rw('c55.dat', c55, .FALSE.)
    call material_rw('c66.dat', c66, .FALSE.)
    call material_rw('rho.dat', rho, .FALSE. )

  end subroutine stiffness_write

  ! ---------------------------------------------------------------------------
  subroutine material_rw(filename, image_data, readfile)

  implicit none

  integer,parameter :: dp = kind(0.d0)
  character(len=7) :: filename
  real(kind=dp),dimension(:,:) :: image_data
  logical :: readfile

  
  open(unit = 13, form="unformatted", file = trim(filename))

  if ( readfile ) then
    read(13) image_data
  else
    write(13) image_data
  endif

  close(unit = 13)
  
  end subroutine material_rw

! -----------------------------------------------------------------------------

subroutine cpml_coeffs(nx, dx, dt, npml, sig_max, k_max, alpha_max, &
            kappa, alpha, acoeff, bcoeff, HALF)

implicit none

integer,parameter :: dp=kind(0.d0)
integer :: i

! Define real inputs 
real(kind=dp) :: dx, dt, sig_max, k_max, alpha_max 
integer :: nx, npml
logical :: HALF

! define the output arrays
real(kind=dp),dimension(nx) :: kappa, alpha, acoeff, bcoeff

! Define all other variables needed in the program
real(kind=dp) :: xoriginleft, xoriginright
real(kind=dp),dimension(nx) :: xval, sigma
integer,parameter :: NP = 2, NPA = 2

real(kind=dp) :: abscissa_in_pml, abscissa_normalized

! -----------------------------------------------------------------------------
sigma(:) = 0.d0

do i=1,nx 
  xval(i) = dx * dble(i - 1)
enddo

if (HALF) then 
    xval = xval + dx/2.0
endif

xoriginleft = dx * dble( npml )
xoriginright = dx * dble( (NX-1) - npml )

do i=1,nx
    !---------- left edge
    abscissa_in_PML = xoriginleft - xval(i)
    if (abscissa_in_PML >= 0.d0) then
        abscissa_normalized = abscissa_in_PML / dble(dx * npml)
        sigma(i) = sig_max * abscissa_normalized**NP
        ! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
        kappa(i) = 1.d0 + (K_MAX - 1.d0) * abscissa_normalized**NP
        alpha(i) = ALPHA_MAX * (1.d0 - abscissa_normalized)**NPA
    endif

    !---------- right edge
    ! define damping profile at the grid points
    abscissa_in_PML = xval(i) - xoriginright
    if (abscissa_in_PML >= 0.d0) then
      abscissa_normalized = abscissa_in_PML / dble(dx * npml)
      sigma(i) = sig_max * abscissa_normalized**NP
      kappa(i) = 1.d0 + (k_max - 1.d0) * abscissa_normalized**NP
      alpha(i) = alpha_max * (1.d0 - abscissa_normalized)**NPA
    endif

    ! just in case, for -5 at the end
    if (alpha(i) < 0.d0) alpha(i) = 0.d0
    ! Compute the b_i coefficents
    bcoeff(i) = exp( - (sigma(i) / kappa(i) + alpha(i)) * DT )
    
    ! Compute the a_i coefficients
    ! this to avoid division by zero outside the PML
    if (abs(sigma(i)) > 1.d-6) then 
      acoeff(i) = sigma(i) * (bcoeff(i) - 1.d0) / ( (sigma(i) + kappa(i) * alpha(i)) ) / kappa(i)
    endif

enddo 

end subroutine cpml_coeffs

!==============================================================================


!==============================================================================
subroutine seismic_cpml_25d(nx, ny, nz, dx, dy, dz, &
                      npoints_pml, src, f0, nstep, angle_force)

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
real(kind=dp), dimension(nx,nz) :: c11, c12, c13, c22, c23, c33, c44, c55, c66, rho
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
real(kind=dp), dimension(2) :: ANGLE_FORCE 


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
real(kind=dp), parameter :: K_MAX = 1.d1
real(kind=dp) :: ALPHA_MAX

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

! Name the f2py inputs 
!f2py3 intent(in) :: nx, ny, nz, dx, dy, dz,
!f2py3 intent(in) :: noints_pml, src, f0, nstep, angle_force

! ------------------------ Load Stiffness Coefficients ------------------------

call material_rw('c11.dat', c11, .TRUE.)
call material_rw('c12.dat', c12, .TRUE.)
call material_rw('c13.dat', c13, .TRUE.)
call material_rw('c22.dat', c22, .TRUE.)
call material_rw('c23.dat', c23, .TRUE.)
call material_rw('c33.dat', c33, .TRUE.)
call material_rw('c44.dat', c44, .TRUE.)
call material_rw('c55.dat', c55, .TRUE.)
call material_rw('c66.dat', c66, .TRUE.)
call material_rw('rho.dat', rho, .TRUE.)

! ------------------------ Assign some constants -----------------------

isource = src(1)+npoints_pml
jsource = src(2)+npoints_pml
ksource = src(3)+npoints_pml


t0 = 1.0d0/f0
DT = minval( (/dx,dy,dz/) )/ ( 2.0* sqrt( ( maxval( (/ c11/rho, c22/rho, c33/rho /) ) ) ) )


ALPHA_MAX = PI*f0 ! from Festa and Vilotte
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
thickness_PML_x = npoints_pml * dx
thickness_PML_y = npoints_pml * dy
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
  a_z(:) = 0.d0
  a_z_half(:) = 0.d0


  ! to compute d0 below, and for stability estimate
  ! quasi_cp_max = sqrt( maxval( (/ c11/rho, c22/rho, c33/rho /) ) )
  quasi_cp_max = ( minval( (/dx,dy,dz/) )/ ( 2.0 * dt) )

  ! compute d0 from INRIA report section 6.1 http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
  d0_x = - dble(NPOWER + 1) * quasi_cp_max * log(Rcoef) / (2.d0 * thickness_PML_x)
  d0_y = - dble(NPOWER + 1) * quasi_cp_max * log(Rcoef) / (2.d0 * thickness_PML_y)
  d0_z = - dble(NPOWER + 1) * quasi_cp_max * log(Rcoef) / (2.d0 * thickness_PML_z)


! ------------------------- damping in the X direction-------------------------

call cpml_coeffs(nx, dx, dt, npoints_pml, d0_x, k_max, alpha_max, &
            K_x, alpha_x, a_x, b_x, .FALSE.)

call cpml_coeffs(nx, dx, dt, npoints_pml, d0_x, k_max, alpha_max, &
            K_x_half, alpha_x_half, a_x_half, b_x_half, .TRUE.)

! ------------------------ damping in the Y direction -------------------------
call cpml_coeffs(ny, dy, dt, npoints_pml, d0_y, k_max, alpha_max, &
            K_y, alpha_y, a_y, b_y, .FALSE.)

call cpml_coeffs(ny, dy, dt, npoints_pml, d0_y, k_max, alpha_max, &
            K_y_half, alpha_y_half, a_y_half, b_y_half, .TRUE.)

! ------------------------- damping in the Z direction-------------------------

call cpml_coeffs(nz, dz, dt, npoints_pml, d0_z, k_max, alpha_max, &
            K_z, alpha_z, a_z, b_z, .FALSE.)

call cpml_coeffs(nz, dz, dt, npoints_pml, d0_z, k_max, alpha_max, &
            K_z_half, alpha_z_half, a_z_half, b_z_half, .TRUE.)




! Print the PML values to a file to check the values
  open(unit = 15, file = "x_values_sub.txt")
  do i=1,nx
    write(15,"(E10.3,E10.3,E10.3,E10.3,E10.3,E10.3,E10.3,E10.3)") &
          a_x(i), a_x_half(i), b_x(i), b_x_half(i), alpha_x(i), alpha_x_half(i), K_x(i), K_x_half(i)
  enddo
  close(15)

  open(unit = 16, file = "y_values_sub.txt")
  do i = 1,ny
    write(16,"(E10.3,E10.3,E10.3,E10.3,E10.3,E10.3,E10.3,E10.3)") &
          a_y(i), a_y_half(i), b_y(i), b_y_half(i), alpha_y(i), alpha_y_half(i), K_y(i), K_y_half(i)
  enddo
  close(16)

  open(unit = 17, file = "z_values_sub.txt")
  do i = 1,nz
    write(17, "(E10.3,E10.3,E10.3,E10.3,E10.3,E10.3,E10.3,E10.3)")  &
          a_z(i), a_z_half(i), b_z(i), b_z_half(i), alpha_z(i), alpha_z_half(i), K_z(i), K_z_half(i)
  enddo
  close(17)

stop

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


! !---
! !---  beginning of time loop
! !---

do it = 1,NSTEP
  !------------------------------------------------------------
  ! compute stress sigma and update memory variables for C-PML
  !------------------------------------------------------------
  ! Update in the x direction
  do k = 2,nz
    do j = 2,NY
      do i = 1,NX-1

        value_dvx_dx = (vx(i+1,j,k) - vx(i,j,k) ) / dx
        value_dvy_dy = (vy(i,j,k) - vy(i,j-1,k) ) / dy
        value_dvz_dz = (vz(i,j,k) - vz(i,j,k-1) ) / dz

        memory_dvx_dx(i,j,k) = b_x_half(i) * memory_dvx_dx(i,j,k) + a_x_half(i) * value_dvx_dx
        memory_dvy_dy(i,j,k) = b_y(j) * memory_dvy_dy(i,j,k) + a_y(j) * value_dvy_dy
        memory_dvz_dz(i,j,k) = b_z(k) * memory_dvz_dz(i,j,k) + a_z(k) * value_dvz_dz

        value_dvx_dx = value_dvx_dx / K_x_half(i) + memory_dvx_dx(i,j,k)
        value_dvy_dy = value_dvy_dy / K_y(j) + memory_dvy_dy(i,j,k)
        value_dvz_dz = value_dvz_dz / K_z(k) + memory_dvz_dz(i,j,k)

        sigmaxx(i,j,k) = sigmaxx(i,j,k) + ( c11(i,k) * value_dvx_dx + c12(i,k) * value_dvy_dy + c13(i,k) * value_dvz_dz) * dt
        sigmayy(i,j,k) = sigmayy(i,j,k) + ( c12(i,k) * value_dvx_dx + c22(i,k) * value_dvy_dy + c23(i,k) * value_dvz_dz) * dt
        sigmazz(i,j,k) = sigmazz(i,j,k) + ( c13(i,k) * value_dvx_dx + c23(i,k) * value_dvy_dy + c33(i,k) * value_dvz_dz) * dt

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

        sigmaxy(i,j,k) = sigmaxy(i,j,k) + c66(i,k) * (value_dvy_dx + value_dvx_dy) * DT

      enddo
    enddo
  enddo

  ! Update sigmaxz, z-direction is full nodes
  do k = 1,nz-1
    do j = 1,ny
      do i = 2,nx

        value_dvz_dx = (vz(i,j,k) - vz(i-1,j,k) ) / dx
        value_dvx_dz = (vx(i,j,k+1) - vx(i,j,k) ) / dz

        memory_dvz_dx(i,j,k) = b_x(i) * memory_dvz_dx(i,j,k) + a_x(i) * value_dvz_dx
        memory_dvx_dz(i,j,k) = b_z_half(k) * memory_dvx_dz(i,j,k) + a_z_half(k) * value_dvx_dz

        value_dvz_dx = value_dvz_dx / K_x(i) + memory_dvz_dx(i,j,k) 
        value_dvx_dz = value_dvx_dz / K_z_half(k) + memory_dvx_dz(i,j,k)

        sigmaxz(i,j,k) = sigmaxz(i,j,k) + c55(i,k) * ( value_dvx_dz + value_dvz_dx) * dt 

      enddo
    enddo
  ! enddo

  !   ! update sigmayz, y-direction is full nodes
  ! do k = 1,nz-1
    do j = 1,ny-1
      do i = 1,nx

        value_dvy_dz = (vy(i,j+1,k) - vy(i,j,k) ) / dz
        value_dvz_dy = (vz(i,j,k+1) - vz(i,j,k) ) / dy

        memory_dvz_dy(i,j,k) = b_y_half(j) * memory_dvz_dy(i,j,k) + a_y_half(j) * value_dvz_dy
        memory_dvy_dz(i,j,k) = b_z_half(k) * memory_dvy_dz(i,j,k) + a_z_half(k) * value_dvy_dz 

        value_dvy_dz = value_dvy_dz / K_z_half(k) + memory_dvy_dz(i,j,k)
        value_dvz_dy = value_dvz_dy / K_y_half(j) + memory_dvz_dy(i,j,k)

        sigmayz(i,j,k) = sigmayz(i,j,k)  + c44(i,k) * ( value_dvy_dz + value_dvz_dy) * dt 

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
        value_dsigmaxz_dz = (sigmaxz(i,j,k) - sigmaxz(i,j,k-1) ) / dz

        memory_dsigmaxx_dx(i,j,k) = b_x(i) * memory_dsigmaxx_dx(i,j,k) + a_x(i) * value_dsigmaxx_dx
        memory_dsigmaxy_dy(i,j,k) = b_y(j) * memory_dsigmaxy_dy(i,j,k) + a_y(j) * value_dsigmaxy_dy
        memory_dsigmaxz_dz(i,j,k) = b_z(k) * memory_dsigmaxz_dz(i,j,k) + a_z(k) * value_dsigmaxz_dz

        value_dsigmaxx_dx = value_dsigmaxx_dx / K_x(i) + memory_dsigmaxx_dx(i,j,k)
        value_dsigmaxy_dy = value_dsigmaxy_dy / K_y(j) + memory_dsigmaxy_dy(i,j,k)
        value_dsigmaxz_dz = value_dsigmaxz_dz / K_z(k) + memory_dsigmaxz_dz(i,j,k) 

        vx(i,j,k) = vx(i,j,k) + (value_dsigmaxx_dx + value_dsigmaxy_dy + value_dsigmaxz_dz) * dt / rho(i,k)

      enddo
    enddo
  ! enddo

  ! do k = 2,nz
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

        vy(i,j,k) = vy(i,j,k) + (value_dsigmaxy_dx + value_dsigmayy_dy + value_dsigmayz_dz) * dt / rho(i,k)

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

        vz(i,j,k) = vz(i,j,k) + (value_dsigmaxz_dx + value_dsigmayz_dy + value_dsigmazz_dz) * dt / rho(i,k)

      enddo
    enddo
  enddo

  ! add the source (force vector located at a given grid point)
  t = dble(it-1)*DT

  ! Gaussian
  source_term = factor * 2.d0*exp(-a*(t-t0)**2)

!   ! first derivative of a Gaussian
!   source_term = - factor * 2.d0*a*(t-t0)*exp(-a*(t-t0)**2)


  ! Use spherical coordinates for the source rotation
  force_x = sin( angle_force(1) ) * cos( angle_force(2) ) * source_term
  force_y = sin( angle_force(1) ) * sin( angle_force(2) ) * source_term
  force_z = cos( angle_force(1) ) * source_term

  vx(isource,jsource,ksource) = vx(isource,jsource,ksource) + force_x * DT / rho(i,k)
  vy(isource,jsource,ksource) = vy(isource,jsource,ksource) + force_y * DT / rho(i,k)
  vz(isource,jsource,ksource) = vz(isource,jsource,ksource) + force_z * DT / rho(i,k)

  ! Dirichlet conditions (rigid boundaries) on the edges or at the bottom of the PML layers
  vx(1,:,:) = 0.d0
  vy(1,:,:) = 0.d0
  vz(1,:,:) = 0.d0

  vx(:,1,:) = 0.d0
  vy(:,1,:) = 0.d0
  vz(:,1,:) = 0.d0
  
  vx(:,:,1) = 0.d0
  vy(:,:,1) = 0.d0
  vz(:,:,1) = 0.d0
  
  vx(NX,:,:) = 0.d0
  vy(NX,:,:) = 0.d0
  vz(NX,:,:) = 0.d0
  
  vx(:,NY,:) = 0.d0
  vy(:,NY,:) = 0.d0
  vz(:,NY,:) = 0.d0
  
  vx(:,:,NZ) = 0.d0
  vy(:,:,NZ) = 0.d0
  vz(:,:,NZ) = 0.d0


  ! output information
  if (mod(it,40) == 0 .or. it == 1) then

  ! print maximum of norm of velocity
  velocnorm = maxval( sqrt(vx**2 + vy**2 + vz**2) )
  print *,'Time step # ',it,' out of ',NSTEP
  print *,'Time: ',(it-1)*DT,' seconds'
  print *,'Max vals for vx, vy, vz: ', maxval(vx), maxval(vy), maxval(vz)

  if (velocnorm > STABILITY_THRESHOLD) stop 'code became unstable and blew up'

  call write_image(vx, nx, ny, nz, it, 'Vx')
  call write_image(vy, nx, ny, nz, it, 'Vy')
  call write_image(vz, nx, ny, nz, it, 'Vz')

  endif
  
enddo   ! end of time loop

end subroutine seismic_cpml_25d


!==============================================================================

subroutine write_image(image_data, nx, ny, nz, it, channel)
! Write the 3D array out as single precision

implicit none

integer, parameter :: dp = kind(0.d0)
integer :: nx, ny, nz, it
real(kind=dp) :: image_data(nx, ny, nz)
character(len=2) :: channel
character(len=100) :: filename

WRITE (filename, "(a2, i6.6, '.dat')" ) channel, it

open(unit = 10, form = "unformatted", file = trim(filename) )
write(10) sngl(image_data)

close(unit = 10)

end subroutine write_image



end module seismicFDTD25d