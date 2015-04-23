!============
! Author: Elmo Tempel
! Date (last changed): 10.02.2015
!============
!
! rk defines real_kind precision... use double precision everywhere
!
module constants
	implicit none
	!integer,parameter :: int_kind = selected_int_kind(9)
	integer,parameter  :: rk = kind(1.0D0)  !realkind: double
    integer,parameter  :: rsp = kind(1.0)  !RealKind: single
    integer,parameter  :: rdp = kind(1.0D0)  !RealKind: double
	!
	real(rk),parameter:: pi = 4*atan(1.0_rk)
	real(rk),parameter:: twopi = 8*atan(1.0_rk)
	real(rk),parameter:: pihalf = 2*atan(1.0_rk)
end module constants