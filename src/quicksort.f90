!============
! Author: Elmo Tempel
! Date (last changed): 24.03.2015
!============
!
module quicksort
	use constants
	implicit none
	private
    interface qsort
        module procedure sort_rrri, sort_ii,sort_r, sort_rri, sort_riii
    end interface
	public:: qsort
contains
	subroutine sort_rrri(arr,slave1,slave2,slaveint)
		implicit none
		real(rk), dimension(:), intent(inout) :: arr,slave1,slave2
		integer,dimension(:),intent(inout):: slaveint
		integer, dimension(size(arr)) :: index
		call myindexx(arr,index)
		arr=arr(index)
		slave1=slave1(index)
		slave2=slave2(index)
		slaveint=slaveint(index)
	end subroutine sort_rrri
	!
	subroutine sort_riii(arr,slave1,slave2,slave3)
		implicit none
		real(rk), dimension(:), intent(inout) :: arr
		integer,dimension(:),intent(inout):: slave1,slave2,slave3
		integer, dimension(size(arr)) :: index
		call myindexx(arr,index)
		arr=arr(index)
		slave1=slave1(index)
		slave2=slave2(index)
		slave3=slave3(index)
	end subroutine sort_riii
	!
	subroutine sort_rri(arr,slave1,slaveint)
		implicit none
		real(rk), dimension(:), intent(inout) :: arr,slave1
		integer,dimension(:),intent(inout):: slaveint
		integer, dimension(size(arr)) :: index
		call myindexx(arr,index)
		arr=arr(index)
		slave1=slave1(index)
		slaveint=slaveint(index)
	end subroutine sort_rri
	!
	subroutine sort_ii(arr,slave1)
		implicit none
		integer, dimension(:), intent(inout) :: arr
		integer,dimension(:),intent(inout):: slave1
		real(rk),dimension(1:size(arr)):: arr2
		integer, dimension(size(arr)) :: index
		arr2=arr
		call myindexx(arr2,index)
		arr=arr(index)
		slave1=slave1(index)
	end subroutine sort_ii
	!
	subroutine sort_r(arr)
		implicit none
		real(rk), dimension(:), intent(inout) :: arr
		integer, dimension(size(arr)) :: index
		call myindexx(arr,index)
		arr=arr(index)
	end subroutine sort_r
	!
	! performs a quicksort for array arr
	! arr is not modified, index vector gives an order for sorted array
	! ... can be used to sort several arrays simultaneously
	subroutine myindexx(arr,index)
		implicit none
		real(rk), dimension(:), intent(in) :: arr
		integer, dimension(:), intent(out) :: index
		integer, parameter :: nn=16, nstack=128
		real(rk) :: a
		integer :: n,k,i,j,indext,jstack,l,r
		integer, dimension(nstack) :: istack
		n=size(arr)
		forall(i=1:n) index(i)=i
		jstack=0
		l=1
		r=n
		do
			if (r-l < nn) then
				do j=l+1,r
					indext=index(j)
					a=arr(indext)
					do i=j-1,l,-1
						if (arr(index(i)) <= a) exit
						index(i+1)=index(i)
					end do
					index(i+1)=indext
				end do
				if (jstack == 0) return
				r=istack(jstack)
				l=istack(jstack-1)
				jstack=jstack-2
			else
				k=(l+r)/2
				call swap_i(index(k),index(l+1))
				call icomp_xchg(index(l),index(r))
				call icomp_xchg(index(l+1),index(r))
				call icomp_xchg(index(l),index(l+1))
				i=l+1
				j=r
				indext=index(l+1)
				a=arr(indext)
				do
					do
						i=i+1
						if (arr(index(i)) >= a) exit
					end do
					do
						j=j-1
						if (arr(index(j)) <= a) exit
					end do
					if (j < i) exit
					call swap_i(index(i),index(j))
				end do
				index(l+1)=index(j)
				index(j)=indext
				jstack=jstack+2
				if (jstack > nstack) stop 'indexx: nstack too small'
				if (r-i+1 >= j-l) then
					istack(jstack)=r
					istack(jstack-1)=i
					r=j-1
				else
					istack(jstack)=j-1
					istack(jstack-1)=l
					l=i
				end if
			end if
		end do
	contains
		subroutine icomp_xchg(i,j)
			integer, intent(inout) :: i,j
			integer :: swp
			if (arr(j) < arr(i)) then
				swp=i
				i=j
				j=swp
			end if
		end subroutine icomp_xchg
	end subroutine myindexx
	subroutine swap_i(a,b)
		integer, intent(inout) :: a,b
		integer :: dum
		dum=a
		a=b
		b=dum
	end subroutine swap_i
end module quicksort
