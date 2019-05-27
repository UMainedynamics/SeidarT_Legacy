

module material_matrices_rw_debug

  implicit none

  !==============================================================================
  subroutine stiffness_write(im, mlist, c11, c12, c13, c22, c23, c33, c44, &
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
    integer :: m, n, i, j, k, npoints_pml
    real(kind=dp), dimension(:,:) :: mlist
    real(kind=dp), dimension(:,:) :: c11, c12, c13, c22, c23, c33, c44, c55, c66, rho

    m = size(rho, 1)
    n = size(rho, 2)

    !Assign between the PML regions
    do i=npoints_pml+1, m+npoints_pml
      do j=npoints_pml+1, n+npoints_pml
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
      c11( m+npoints_pml-1+i, :) = c11(m+npoints_pml,:)
      c12( m+npoints_pml-1+i, :) = c12(m+npoints_pml,:)
      c13( m+npoints_pml-1+i, :) = c13(m+npoints_pml,:)
      c22( m+npoints_pml-1+i, :) = c22(m+npoints_pml,:)
      c23( m+npoints_pml-1+i, :) = c23(m+npoints_pml,:)
      c33( m+npoints_pml-1+i, :) = c33(m+npoints_pml,:)
      c44( m+npoints_pml-1+i, :) = c44(m+npoints_pml,:)
      c55( m+npoints_pml-1+i, :) = c55(m+npoints_pml,:)
      c66( m+npoints_pml-1+i, :) = c66(m+npoints_pml,:)
      rho( m+npoints_pml-1+i, :) = rho(m+npoints_pml,:)

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
      c11( :, n+npoints_pml-1+i) = c11(:,n+npoints_pml)
      c12( :, n+npoints_pml-1+i) = c12(:,n+npoints_pml)
      c13( :, n+npoints_pml-1+i) = c12(:,n+npoints_pml)      
      c22( :, n+npoints_pml-1+i) = c22(:,n+npoints_pml)
      c23( :, n+npoints_pml-1+i) = c23(:,n+npoints_pml)
      c33( :, n+npoints_pml-1+i) = c33(:,n+npoints_pml)
      c44( :, n+npoints_pml-1+i) = c44(:,n+npoints_pml)
      c55( :, n+npoints_pml-1+i) = c55(:,n+npoints_pml)      
      c66( :, n+npoints_pml-1+i) = c66(:,n+npoints_pml)
      rho( :, n+npoints_pml-1+i) = rho(:,n+npoints_pml)

    end do 

    ! Write each of the matrices to file
    call material_write('c11.dat', c11, )
    call material_write('c12.dat',)
    call material_write('c13.dat',)
    call material_write('c22.dat',)
    call material_write('c23.dat',)
    call material_write('c33.dat',)
    call material_write('c44.dat',)
    call material_write('c55.dat',)
    call material_write('c66.dat',)
    call material_write('rho.dat', )

  end subroutine material_assign

  ! ---------------------------------------------------------------------------
  subroutine material_rw(filename, image_data, readfile)

  implicit none

  character(len=6) :: filename
  real(kind=dp) :: image_data
  logical :: readfile

  open(unit = 10, form = 'unformatted', file = trim(filename) )

  if (readfile) then  
    read(10) image_data
  else
    write(10) image_data
  endif

  close(unit = 10)
  
  end subroutine material_rw

  ! ---------------------------------------------------------------------------


end module material_matrices_rw_debug