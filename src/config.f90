!========================
! Author: Elmo Tempel
! Date (last changed): 07.11.2014
!========================
!
! Defines functions to handle configuration files
! .. all input parameters are that control the program defined in config.ini file
!
! It is not very intelligent module, but it does it job
! 
! Configuration parameters are divided into section: [section_name]
! The order of sections and parameters is arbitrary
! comments can be everywhere, and start with #
! ... see example config.ini file for more details.
!
module config
	use constants
	implicit none
	private
    interface get_conf
        module procedure get_conf_real,get_conf_int,get_conf_char,get_conf_log
    end interface
	character(len=:),allocatable,private:: cfile ! global config file
	logical:: warn_key_not_present
	public:: get_conf,set_conf_file
contains
	subroutine set_conf_file(conf_file,warn)
		implicit none
		character(len=*),intent(in):: conf_file
		logical,intent(in),optional:: warn
		cfile=conf_file
		warn_key_not_present=.true.
		if (present(warn)) warn_key_not_present=warn
	end subroutine set_conf_file
	!
	!==========================================================
	!
	subroutine get_conf_real(sect,key,value,default)
		implicit none
		character(len=*),intent(in):: sect,key
		real(rk),intent(in),optional:: default
		real(rk),intent(out):: value
		character(len=:),allocatable:: cval
		integer:: ierr
		call get_cval(sect=sect,key=key,cval=cval)
		read(cval, iostat=ierr, fmt=*) value
		if (ierr/=0) then
			if (present(default)) then
				value=default
				if (warn_key_not_present) then
					print*, "WARNING: use default value: ", sect,' ', key,' ',default
				end if
			else
				print*, "ERROR default not present and value not in file: ", sect,' ',key
				stop "Program terminated by conf_handle>get_conf_real"
			end if
		end if
	end subroutine get_conf_real
	subroutine get_conf_log(sect,key,value,default)
		implicit none
		character(len=*),intent(in):: sect,key
		logical,intent(in),optional:: default
		logical,intent(out):: value
		character(len=:),allocatable:: cval
		integer:: ierr
		call get_cval(sect=sect,key=key,cval=cval)
		read(cval, iostat=ierr, fmt=*) value
		if (ierr/=0) then
			if (present(default)) then
				value=default
				if (warn_key_not_present) then
					print*, "WARNING: use default value: ", sect,' ', key,' ',default
				end if
			else
				print*, "ERROR default not present and value not in file: ", sect,' ',key
				stop "Program terminated by conf_handle>get_conf_log"
			end if
		end if
	end subroutine get_conf_log
	subroutine get_conf_int(sect,key,value,default)
		implicit none
		character(len=*),intent(in):: sect,key
		integer,intent(in),optional:: default
		integer,intent(out):: value
		character(len=:),allocatable:: cval
		integer:: ierr
		call get_cval(sect=sect,key=key,cval=cval)
		read(cval, iostat=ierr, fmt=*) value
		if (ierr/=0) then
			if (present(default)) then
				value=default
				if (warn_key_not_present) then
					print*, "WARNING: use default value: ", sect,' ', key,' ',default
				end if
			else
				print*, "ERROR default not present and value not in file: ", sect,' ',key
				stop "Program terminated by conf_handle>get_conf_int"
			end if
		end if
	end subroutine get_conf_int
	subroutine get_conf_char(sect,key,value,default)
		implicit none
		character(len=*),intent(in):: sect,key
		character(len=*),intent(in),optional:: default
		character(len=:),allocatable,intent(out):: value
		call get_cval(sect=sect,key=key,cval=value)
		if (.not.allocated(value)) then
			if (present(default)) then
				value=default
				if (warn_key_not_present) then
					print*, "WARNING: use default value: ", sect,' ', key,' ',default
				end if
			else
				print*, "ERROR default not present and value not in file: ", sect,' ',key
				stop "Program terminated by conf_handle>get_conf_char"
			end if
		end if
	end subroutine get_conf_char
	!
	!==========================================================
	!
	subroutine get_cval(sect,key,cval)
		implicit none
		character(len=*),intent(in):: key
		character(len=*),intent(in),optional:: sect
		character(len=:),allocatable,intent(out):: cval
		integer:: iunit,ierr,n
		character(len=300):: cline
		open(newunit=iunit,file=cfile,iostat=ierr,status='old',readonly)
		if (ierr/=0) then
			print*, "ERROR reading configuration file: ", cfile
			print*, " .. file do not exist/corrupted!"
			stop "Program terminated by conf_handle>get_cval"
		end if
		! go to the beginning of a section
		if (present(sect)) then
			dosect: do
				read(iunit,fmt='(A)',iostat=ierr) cline
				if (ierr/=0) then
					print*, "ERROR reading configuration file: ", cfile
					print*, " .. requested section do not exist: ", sect,' ',key
					stop "Program terminated by conf_handle>get_cval"
				end if
				cline=adjustl(cline)
				if (cline(1:1)=='[') then
					n=index(cline,']')
					if (sect==cline(2:n-1)) exit dosect
				end if
			end do dosect
		end if
		! start scanning key
		dokey: do
			read(iunit,fmt='(A)',iostat=ierr) cline
			if (ierr/=0) then
				exit dokey
			end if
			cline=adjustl(cline)
			! test if starting next section
			if (present(sect)) then
				if (cline(1:1)=='[') then
					exit dokey
				end if
			end if
			! find the requested key
			n=index(cline,'=')
			if (key==cline(1:n-1)) then
				cline=adjustl(cline(n+1:))
				n=index(cline,'#')
				if (n>0) then
					cval=cline(1:n-1)
					cval=cval(1:len_trim(cval))
				else
					cval=cline(1:len_trim(cline))
				end if
				exit dokey
			end if
		end do dokey
		close(iunit,status='keep',iostat=ierr)
	end subroutine get_cval
	!
end module config
