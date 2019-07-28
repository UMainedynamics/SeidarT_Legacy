! Python Compute the orientation tensor of a fabric
!
! compile using 
!	
!	f2py3 -c --fcompiler=gnu95 -m orientsynth orientsynth.f95  
!
! 
! =====================================================================



module orientsynth

implicit NONE

contains


!------------------------------------------------------------------------------
subroutine fabric_ortosynth(trend, plunge, anglemin, anglemax, npts, &
    euler_list, orientation_tensor)
! FABRIC_ORTOSYNTH computes the euler angles of a fabric from its statistical
! description trend, plunge, anglemin and anglemax. Angles are input as 
! degrees. The number of points (npts) corresponds to the number of equal sized 
! grains. The outputs are written to the specified file as comma delimited text.
!
! INPUTS
!	trend (float)
!	plunge (float)
!	anglemin (float)
!	anglemax (float)
!	npts (integer)
!	filename (character)
!
! Created by Steven Bernsen 5/16/2019
! University of Maine
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

implicit none

! Define the precision
integer,parameter :: sp = kind(1.0)

! Declare the inputs 
real(kind=sp),dimension(1) :: trend, plunge, anglemax, anglemin 
integer :: npts


! Declare necessary values to compute the 
integer :: rownum
real(kind=sp),parameter :: pi = 3.141592653589793
real(kind=sp),dimension(1) :: RND,RND1,RND2,angle

! Define arrays
real(kind=sp),dimension(3) :: euler_avg, euler_dev
real(kind=sp),dimension(3,3) :: euler_axis_rot, euler_dev_rot, composite_rot
real(kind=sp), dimension(npts,3) :: end_orientation
real(kind=sp),intent(out) :: orientation_tensor(3,3)
real(kind=sp),intent(out) :: euler_list(npts,3) 

! Make this a python module
!f2py3 intent(in) :: trend, plunge, anglemax, anglemin, filename
!f2py3 intent(out) :: euler_list, orientation_tensor

! -----------------------------------------------------------------------------

! -----
! ----- Get Started 
! -----

euler_avg = (/ 0.0, 90.0-plunge, -trend /)

! Convert everything to radians
euler_avg = euler_avg*pi/180.0
anglemin = anglemin*pi/180.0
anglemax = anglemax*pi/180.0

! Get the rotation matrix from the euler_avg
euler_axis_rot = transpose( rotator(euler_avg(1), euler_avg(2), euler_avg(3) ) )

! start the counter
rownum=1
do while( rownum .LE. npts)
    ! Get a random number 
    call random_number(RND) 

    ! Just for asthetics lets do this calculation
    angle = acos(2.0*RND-1.0)

    if  ( (angle(1)<=anglemax(1) ) .AND. (abs( angle(1) )>=anglemin(1) ) ) then
        call random_number(RND1)
        call random_number(RND2)

        ! Assign the angles then create the rotation matrix 
        euler_dev(:) = (/ RND1*2.0*pi, angle, RND2*2.0*pi /)
        euler_dev_rot = transpose( rotator( euler_dev(1), euler_dev(2), euler_dev(3) ) )

        ! Matrix multiplication of the two rotation matrices
        composite_rot = matmul(euler_axis_rot, euler_dev_rot)

        if ( composite_rot(3,3) .GT. 0.0) then 
            end_orientation(rownum,:) = -composite_rot(:,3)
        else
            end_orientation(rownum,:) = composite_rot(:,3)
        end if

        ! Assign the euler_list values
        euler_list(rownum,2) = acos( composite_rot(3,3) )

        if( composite_rot(3,3) .LT. 0.0 ) then 
            euler_list(rownum,1) = pi + atan2(-composite_rot(1,3),&
                -composite_rot(2,3) )
            euler_list(rownum,3) = atan2(-composite_rot(3,1), &
                composite_rot(3,2) )
        else 
            euler_list(rownum,1) = atan2(composite_rot(1,3), &
                composite_rot(2,3) )
            euler_list(rownum,3) = atan2(composite_rot(3,1), &
                -composite_rot(3,2) )
        end if 

        rownum = rownum + 1
    end if
end do 


euler_list(:,1) = -( euler_list(:,1) + pi/2)

! Get the orientation tensor
orientation_tensor = orten(end_orientation) 

end subroutine fabric_ortosynth


! -----------------------------------------------------------------------------

function rotator(phi,theta,psi)
! get the z-x-z rotation matrix from three angles

real :: phi, theta, psi
real,dimension(3,3) :: D,C,B,rotator 

D(1,:) = (/ cos(phi), -sin(phi), 0.0 /)
D(2,:) = (/sin(phi), cos(phi), 0.0 /)
D(3,:) = (/0.0, 0.0, 1.0 /)

C(1,:) = (/ 1.0, 0.0, 0.0 /)
C(2,:) = (/ 0.0, cos(theta), -sin(theta) /)
C(3,:) = (/ 0.0, sin(theta), cos(theta) /)

B(1,:) = (/ cos(psi), -sin(psi), 0.0 /)
B(2,:) = (/sin(psi), cos(psi), 0.0 /)
B(3,:) = (/0.0, 0.0, 1.0 /)

B = matmul(C, B)
rotator = matmul(D, B)


end function rotator


! -----------------------------------------------------------------------------

function orten(end_orientation)

implicit none

! compute the rotation tensor from the orientation matrix
integer :: npts
real :: end_orientation(:,:) 
real,allocatable :: ortho_comp(:,:)
real,dimension(3,3) :: orten 

npts = size(end_orientation, 1)

allocate( ortho_comp(npts, 9) )

ortho_comp(:,1:3) = end_orientation(:,1:3)
ortho_comp(:,4) = end_orientation(:,1)*end_orientation(:,1)

ortho_comp(:,5) = end_orientation(:,1)*end_orientation(:,2)
ortho_comp(:,6) = end_orientation(:,1)*end_orientation(:,3)
ortho_comp(:,7) = end_orientation(:,2)*end_orientation(:,2)
ortho_comp(:,8) = end_orientation(:,2)*end_orientation(:,3)
ortho_comp(:,9) = end_orientation(:,3)*end_orientation(:,3)


orten(1,:) = (/ sum(ortho_comp(:,4) ), sum(ortho_comp(:,5) ), sum(ortho_comp(:,6) ) /)
orten(2,:) = (/ sum(ortho_comp(:,5) ), sum(ortho_comp(:,7) ), sum(ortho_comp(:,8) ) /)
orten(3,:) = (/ sum(ortho_comp(:,6) ), sum(ortho_comp(:,8) ), sum(ortho_comp(:,9) ) /)

orten = orten/npts 

end function orten


! -----------------------------------------------------------------------


end module orientsynth
