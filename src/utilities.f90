!============
! Author: Elmo Tempel
! Date (last changed): 08.04.2015
!============
!
module utilities
	use constants
	implicit none
	logical,private:: get_time_first_time=.true.
	integer,dimension(1:8),private,save:: t0
	integer,private,save:: tdays0
contains
	!
	! incomplete beta function
	function betainc(a,b,p) result(res)
		use nr, only: betai
		implicit none
		real(rk):: res
		real(rk),intent(in):: a,b,p
		if (a<=0.0) then
			res=1.0
		else
			res = betai(real(a),real(b),real(p))
		end if
	end function betainc
	!
	! d is input vecotr and ex,ey are two output vectors
	! ... ex,ey define the coordinate system perpendicular to d
	subroutine get_perpendicular_vectors(d,ex,ey)
		implicit none
		real(rk),dimension(1:3),intent(in):: d ! input orientation
		real(rk),dimension(1:3),intent(out):: ex,ey ! two output vectors
		real(rk),dimension(1:3):: dd
		integer:: ii(1)
		! normalise the input vector if not normalised
		dd=d/sqrt(sum(d**2))
		! find largest axis to avoid division with zero
		ii=maxloc(dd)
		! find the ex vector using scalar product
		ex=1.0_rk
		ex(ii(1))=-(sum(dd)-dd(ii(1)))/dd(ii(1))
		ex=ex/sqrt(sum(ex**2)) ! normalise
		! find the ey vecotr using vector product
		ey(1)=dd(2)*ex(3)-dd(3)*ex(2)
		ey(2)=dd(3)*ex(1)-dd(1)*ex(3)
		ey(3)=dd(1)*ex(2)-dd(2)*ex(1)
		ey=ey/sqrt(sum(ey**2)) ! normalise
	end subroutine get_perpendicular_vectors
	!
	! returns the time since first call
	subroutine get_time_min(tmin)
		implicit none
		real(rk),intent(out):: tmin
		integer,dimension(1:8):: t
		character(len = 12):: rc(3)
		integer:: i,iday,tdays
		if (get_time_first_time) then
			call date_and_time(rc(1),rc(2),rc(3),t0)
			tdays0=0
			do i=1,t0(2)-1,1
				select case (i)
					case (1); iday=31
					case (2); iday=28
					if (t0(1)==2016 .or. t0(1)==2020) iday=29
					case (3); iday=31
					case (4); iday=30
					case (5); iday=31
					case (6); iday=30
					case (7); iday=31
					case (8); iday=31
					case (9); iday=30
					case (10); iday=31
					case (11); iday=30
					case (12); iday=31
				end select
				tdays0=tdays0+iday
			end do
			tdays0=tdays0+t0(3)
			get_time_first_time=.false.
		end if
		call date_and_time(rc(1),rc(2),rc(3),t)
		!
		tdays=0
		do i=1,t(2)-1,1
			select case (i)
				case (1); iday=31
				case (2); iday=28
				if (t(1)==2016 .or. t(1)==2020) iday=29
				case (3); iday=31
				case (4); iday=30
				case (5); iday=31
				case (6); iday=30
				case (7); iday=31
				case (8); iday=31
				case (9); iday=30
				case (10); iday=31
				case (11); iday=30
				case (12); iday=31
			end select
			tdays=tdays+iday
		end do
		tdays=tdays+t(3)
		if (t(1)>t0(1)) then
			tdays=tdays+365
			if (t0(1)==2016 .or. t0(1)==2020) tdays=tdays+1
		end if
		!
		tmin=(tdays-tdays0)*24*60
		tmin=tmin+(t(5)-t0(5))*60
		tmin=tmin+(t(6)-t0(6))
		tmin=tmin+(t(7)-t0(7))/60.0_rk
		tmin=tmin+(t(8)-t0(8))/60000.0_rk
	end subroutine get_time_min
	!
    function get_file_counter(nr, mitu) result(res)
    	implicit none
    	integer, intent(in) :: nr, mitu
    	integer :: i, N
    	character(len=5) :: fm, fm2
    	character(len=mitu) :: res, res1
    	N=floor(log10(real(nr)))+1
    	if (mitu<N) then
    		print *,  "Liiga pikk arv"
    		forall (i=1:mitu) res(i:i)= "_"
    	else if (nr<0) then
    		print *,  "Kirjutatav nr on neg"
    		forall (i=1:mitu) res(i:i)= " "
    	else if (mitu<1) then
    		print *,  "Liiga suured noudmised"
    		forall (i=1:mitu) res(i:i)= " "
    	else
    		write(fm, "(A2,I2,A1)") "(I",mitu,")"
    		write(fm2, "(A2,I2,A1)") "(I",N,")"
    		forall (i=1:mitu) res(i:i)= "0"
    		write(res(mitu-N+1:mitu), fmt=fm2) nr
    	end if 
    end function get_file_counter
end module utilities