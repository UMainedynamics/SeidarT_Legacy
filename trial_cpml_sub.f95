program trial_cpml 

implicit none

contains



subroutine cpml_coeffs(nx, dx, dt, npml, sig_max, k_max, alpha_max, &
						sigma, kappa, alpha, acoeff, bcoeff, HALF)

implicit none

integer,parameter :: dp=kind(0.d0)

! Define real inputs 
real(kind=dp) :: dx, dt 
integer :: nx, npml
logical :: HALF

! define the output arrays
real(kind=dp),dimension(nx) :: sigma, kappa, alpha, acoeff, bcoeff

! Define all other variables needed in the program
real(kind=dp) :: xoriginleft, xoriginright 
real(kind=dp),dimension(nx) :: xval = dble( (/ (i, i=1,nx, 1) /) ) 
real(kind=dp),parameter :: eps0 = 8.85418782d-12
integer,parameter :: NP = 2, NPA = 3


! ===========================================================

if (HALF) then 
    xval = xval + dx/2.0
endif

xoriginleft = dx * dble( npml )
xoriginright = dx * dble( (NX-1) - npml )

do i=1,nx
    !---------- left edge
    abscissa_in_PML = xoriginleft - xval
    if (abscissa_in_PML >= 0.d0) then
        abscissa_normalized = abscissa_in_PML / dble(dx * npml)
        sigma(i) = sig_max * abscissa_normalized**NP
        ! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
        kappa(i) = 1.d0 + (K_MAX - 1.d0) * abscissa_normalized**NP
        alpha(i) = ALPHA_MAX * (1.d0 - abscissa_normalized)**NPA
    endif

    !---------- right edge
    ! define damping profile at the grid points
    abscissa_in_PML = xval - xoriginright
    if (abscissa_in_PML >= 0.d0) then
      abscissa_normalized = abscissa_in_PML / dble(dx * npml)
      sigma(i) = sig_max * abscissa_normalized**NP
      kappa(i) = 1.d0 + (k_max - 1.d0) * abscissa_normalized**NP
      alpha(i) = alpha_max * (1.d0 - abscissa_normalized)**NPA
    endif

    ! just in case, for -5 at the end
    if (alpha(i) < 0.d0) alpha(i) = 0.d0
    ! Compute the b_i coefficents
    bcoeff(i) = exp( - (sigma(i) / kappa(i) + alpha(i)) * DT/eps0 )
    
    ! Compute the a_i coefficients
    ! this to avoid division by zero outside the PML
    if (abs(sigh_x(i)) > 1.d-6) then 
      acoeff(i) = sigma(i) * (bcoeff(i) - 1.d0) / ( (sigma(i) + kappa(i) * alpha(i)) ) / kappa(i)
    endif

enddo 


end subroutine cpml_coeffs


end program trial_cpml