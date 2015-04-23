!=========================
! rwmode:: 0-readonly, 1-readwrite
! lines: 572
!
! Datatype DATACODE value 
! bit,X 1
! byte, B 11
! logical, L 14
! ASCII character, A 16 
! short integer, I 21 
! integer, J 41 
! real, E 42 
! double precision, D 82 
! complex 83 
! double complex 163
!
!=========================
module fits_op
    use constants
    implicit none
    !
    interface read_fits_header
      module procedure read_fits_header_real, read_fits_header_integer
    end interface read_fits_header
    !
    interface write_fits_header
      module procedure write_fits_header_string, write_fits_header_integer, write_fits_header_real
    end interface write_fits_header
    !
    interface fits_read_table_col
        module procedure fits_read_table_col_real_sp, fits_read_table_col_real_dp, &
        fits_read_table_col_int, fits_read_table_col_intlong
    end interface fits_read_table_col
    !
contains
    !================================================
    ! siin on tabelite jaoks vajalikud subroutined
    subroutine fits_read_table_col_real_sp(file,colname,vec,iunit)
        character(len=*),intent(in):: file ! extension on sisendis [n] failinime taga. fail.fits[1]
        character(len=*),intent(in):: colname
        real(rsp),dimension(:),intent(inout):: vec
        real(rdp),dimension(1:size(vec)):: vec2
        integer:: status,unit,nrows,ncol,colnum,repeat,width,datacode
        logical:: anyf
        integer,intent(in),optional:: iunit
        if (present(iunit)) then
            unit=iunit
            nrows=size(vec)
        else
            call ftgiou(unit,status) ! get unit number
            call ftnopn(unit,file,0,status) ! avame faili rwmode=0
            call ftgnrw(unit,nrows,status) ! loeme kokku read...kontrollime, kas on sisendiga kooskolas
            call ftgncl(unit,ncol,status) ! loeme kokku tulbad
        end if
        call ftgcno(unit,.false.,colname, colnum,status) ! leiame tulba vastavalt nimele
        if (status==219) print*, "fitsio:table: sellist tulpa ei ole: ", colname
        call ftgtcl(unit,colnum, datacode,repeat,width,status) ! leiame tulba formaadi... andmete lugemiseks
        if (repeat>1) print*, "fitsio:table: tegemist on vektoriga... ", colname, repeat
        select case (datacode)
            case default
            print*, "fitsio:table: tulp ei ole real formaadis ", colname,datacode
            case (42)
            call ftgcve(unit,colnum,1,1,nrows,1, vec,anyf,status)
            if (anyf) print*, "fitsio:table: some colvalues are undefined! ", colname
            case (82)
            call ftgcvd(unit,colnum,1,1,nrows,1, vec2,anyf,status)
            if (anyf) print*, "fitsio:table: some colvalues are undefined! ", colname
            vec=vec2
        end select
        if (.not.present(iunit)) then
            call ftclos(unit,status)
            call ftfiou(unit,status)
        end if
    end subroutine fits_read_table_col_real_sp
    subroutine fits_read_table_col_real_dp(file,colname,vec,iunit)
        character(len=*),intent(in):: file ! extension on sisendis [n] failinime taga. fail.fits[1]
        character(len=*),intent(in):: colname
        real(rdp),dimension(:),intent(inout):: vec
        integer:: status,unit,nrows,ncol,colnum,repeat,width,datacode
        logical:: anyf
        integer,intent(in),optional:: iunit
        if (present(iunit)) then
            unit=iunit
            nrows=size(vec)
        else
            call ftgiou(unit,status) ! get unit number
            call ftnopn(unit,file,0,status) ! avame faili rwmode=0
            call ftgnrw(unit,nrows,status) ! loeme kokku read...kontrollime, kas on sisendiga kooskolas
            call ftgncl(unit,ncol,status) ! loeme kokku tulbad
        end if
        call ftgcno(unit,.false.,colname, colnum,status) ! leiame tulba vastavalt nimele
        if (status/=0) call fits_print_error(status)
        if (status==219) print*, "fitsio:table: sellist tulpa ei ole: ", colname
        call ftgtcl(unit,colnum, datacode,repeat,width,status) ! leiame tulba formaadi... andmete lugemiseks
        if (repeat>1) print*, "fitsio:table: tegemist on vektoriga... ", colname, repeat
        select case (datacode)
            case default
            print*, "fitsio:table: tulp ei ole real formaadis ", colname,datacode
            case (42,82)
            call ftgcvd(unit,colnum,1,1,nrows,1, vec,anyf,status)
            if (anyf) print*, "fitsio:table: some colvalues are undefined! ", colname
        end select
        if (.not.present(iunit)) then
            call ftclos(unit,status)
            call ftfiou(unit,status)
        end if
    end subroutine fits_read_table_col_real_dp
    subroutine fits_read_table_col_int(file,colname,vec,iunit)
        character(len=*),intent(in):: file ! extension on sisendis [n] failinime taga. fail.fits[1]
        character(len=*),intent(in):: colname
        integer,dimension(:),intent(inout):: vec
        integer:: status,unit,nrows,ncol,colnum,repeat,width,datacode
        logical:: anyf
        integer,intent(in),optional:: iunit
        if (present(iunit)) then
            unit=iunit
            nrows=size(vec)
        else
            call ftgiou(unit,status) ! get unit number
            call ftnopn(unit,file,0,status) ! avame faili rwmode=0
            call ftgnrw(unit,nrows,status) ! loeme kokku read...kontrollime, kas on sisendiga kooskolas
            call ftgncl(unit,ncol,status) ! loeme kokku tulbad
        end if
        !print*, unit,nrows,colname
        call ftgcno(unit,.false.,colname, colnum,status) ! leiame tulba vastavalt nimele
        if (status/=0) call fits_print_error(status)
        if (status==219) print*, "fitsio:table: sellist tulpa ei ole: ", colname
        call ftgtcl(unit,colnum, datacode,repeat,width,status) ! leiame tulba formaadi... andmete lugemiseks
        if (repeat>1) print*, "fitsio:table: tegemist on vektoriga... ", colname, repeat
        select case (datacode)
            case default
            print*, "fitsio:table: tulp ei ole real formaadis ", colname,datacode
            case (21,41)
            call ftgcvj(unit,colnum,1,1,nrows,1, vec,anyf,status)
            if (anyf) print*, "fitsio:table: some colvalues are undefined! ", colname
        end select
        if (.not.present(iunit)) then
            call ftclos(unit,status)
            call ftfiou(unit,status)
        end if
    end subroutine fits_read_table_col_int
    subroutine fits_read_table_col_intlong(file,colname,vec,iunit)
        character(len=*),intent(in):: file ! extension on sisendis [n] failinime taga. fail.fits[1]
        character(len=*),intent(in):: colname
        integer*8,dimension(:),intent(inout):: vec
        integer:: status,unit,nrows,ncol,colnum,repeat,width,datacode
        logical:: anyf
        integer,intent(in),optional:: iunit
        if (present(iunit)) then
            unit=iunit
            nrows=size(vec)
        else
            call ftgiou(unit,status) ! get unit number
            call ftnopn(unit,file,0,status) ! avame faili rwmode=0
            call ftgnrw(unit,nrows,status) ! loeme kokku read...kontrollime, kas on sisendiga kooskolas
            call ftgncl(unit,ncol,status) ! loeme kokku tulbad
        end if
        call ftgcno(unit,.false.,colname, colnum,status) ! leiame tulba vastavalt nimele
        if (status/=0) call fits_print_error(status)
        if (status==219) print*, "fitsio:table: sellist tulpa ei ole: ", colname
        call ftgtcl(unit,colnum, datacode,repeat,width,status) ! leiame tulba formaadi... andmete lugemiseks
        if (repeat>1) print*, "fitsio:table: tegemist on vektoriga... ", colname, repeat
        select case (datacode)
            case default
            print*, "fitsio:table: tulp ei ole real formaadis ", colname,datacode
            case (21,41,81)
            call ftgcvk(unit,colnum,1,1,nrows,1, vec,anyf,status)
            if (anyf) print*, "fitsio:table: some colvalues are undefined! ", colname
        end select
        if (.not.present(iunit)) then
            call ftclos(unit,status)
            call ftfiou(unit,status)
        end if
    end subroutine fits_read_table_col_intlong
    !================================================
    ! allpool on abisubroutined
    subroutine fits_print_error(status)
        integer,intent(in):: status
        character(len=30):: text
        real:: ver
        call ftvers(ver)
        print*, "cfitsio: ", ver
        call ftgerr(status,text)
        print*, "Error status: ", status
        print*, text
    end subroutine fits_print_error
    subroutine fits_print_errmsg()
        character(len=1):: bye
        character(len=80):: text
        print*, "printing messages. Insert nonspace for exit:"
        do
            call ftgmsg(text)
            read*, bye
            if (bye/=' ') exit
        end do
    end subroutine fits_print_errmsg


















! ================================
!   SIIT EDASI ON VANAD SUBROUTINED
! ================================
  SUBROUTINE read_fits_to_matrix(filename,pix)
    implicit none
    INTEGER                                     :: status, unit, blocksize, rw, nfound, group, fpixel
    CHARACTER(len=*),intent(in):: filename ! fits file name
    REAL(rk),DIMENSION(:,:),ALLOCATABLE,intent(out):: pix
    !REAL(rk),DIMENSION(:,:),ALLOCATABLE:: dpix
    INTEGER,dimension(2):: pix_dim
    real(rk)                                    :: nullval
    logical                                     :: anynull
    status = 0
    call ftgiou(unit, status) ! get io number for fits file....seda ei pea tegelikult kasutama...voib ka ise maarata

    rw = 0 ! 0-readonly, 1-readwrite
    ! blocksize -- ylearune asi...sellele ei tasu tahelepanu poorata.
    call ftopen(unit, filename, rw, blocksize, status)
    ! Get a sequence of numbered keyword values:
    call ftgknj(unit,'NAXIS', 1, 2, pix_dim, nfound, status)
    ! allocate memory
    ALLOCATE(pix(1:pix_dim(1),1:pix_dim(2)))
    !ALLOCATE( dpix( 0:pix_dim(1)-1,0:pix_dim(2)-1 ) )
    ! read values
    group=1
    fpixel=1
    nullval = -999.0
    IF (rk==rdp) THEN
       CALL ftgpvd(unit, group, fpixel, SIZE(pix), nullval, pix, anynull, status)
    else if (rk==rsp) then
       call ftgpve(unit, group, fpixel, size(pix), nullval, pix, anynull, status)
    else
       PRINT*, "ERROR: no subroutine fitsop: ftgpvx (rk do not mach!)"
    END IF
    !pix=dpix
    ! call ftgkyd(unit, 'delta', delta, comment, status)
    ! call ftgkyd(unit, 'scale', scale, comment, status)
    ! call ftgkyd(unit, 'x0', z_corner(1), comment, status)
    ! call ftgkyd(unit, 'y0', z_corner(2), comment, status)
    ! call ftgkyd(unit, 'z0', z_corner(3), comment, status)
    ! call ftgkyd(unit, 'margin', margin, comment, status)
    ! if (present(rank)) call ftgkyj(unit, 'rank', rank, comment, status)
    ! if (present(in_file)) call ftgkls(unit, 'dat file', in_file, comment, status)
    ! if (present(wt_file)) call ftgkls(unit, 'wts file', wt_file, comment, status)
    call ftclos(unit, status)
    call ftfiou(unit, status)
  END SUBROUTINE read_fits_to_matrix

  SUBROUTINE read_fits_to_matrix_int(filename,pix2)
    implicit none
    INTEGER                                     :: status, unit, blocksize, rw, nfound, group, fpixel
    CHARACTER(len=*),intent(in):: filename ! fits file name
    REAL(rk),DIMENSION(:,:),ALLOCATABLE, intent(out):: pix2
    integer(rk),DIMENSION(:,:),ALLOCATABLE :: pix
    !REAL(rk),DIMENSION(:,:),ALLOCATABLE:: dpix
    INTEGER,dimension(2):: pix_dim
    !real(rk)                                    :: nullval
    integer 					:: nullval
    logical                                     :: anynull
    status = 0
    call ftgiou(unit, status) ! get io number for fits file....seda ei pea tegelikult kasutama...voib ka ise maarata

    rw = 0 ! 0-readonly, 1-readwrite
    ! blocksize -- ylearune asi...sellele ei tasu tahelepanu poorata.
    call ftopen(unit, filename, rw, blocksize, status)
    ! Get a sequence of numbered keyword values:
    call ftgknj(unit,'NAXIS', 1, 2, pix_dim, nfound, status)
    ! allocate memory
    print *,  "seda ikka saab allokeerida"
    print *,  ""
    print *,  "pikslite dimensioonid"
    print *,  pix_dim(1), pix_dim(2)
    print *,  ""
    ALLOCATE(pix(1:pix_dim(1),1:pix_dim(2)))
    !ALLOCATE( dpix( 0:pix_dim(1)-1,0:pix_dim(2)-1 ) )
    ! read values
    group=1
    fpixel=1
    !nullval = -999.0
    nullval = -1000
!    IF (rk==rdp) THEN
!       CALL ftgpvd(unit, group, fpixel, SIZE(pix), nullval, pix, anynull, status)
!    else if (rk==rsp) then
!       call ftgpve(unit, group, fpixel, size(pix), nullval, pix, anynull, status)
!    else
!       PRINT*, "ERROR: no subroutine fitsop: ftgpvx (rk do not mach!)"
!    END IF
call ftgpvi(unit, group, fpixel, size(pix), nullval, pix, anynull, status)
    !pix=dpix
    ! call ftgkyd(unit, 'delta', delta, comment, status)
    ! call ftgkyd(unit, 'scale', scale, comment, status)
    ! call ftgkyd(unit, 'x0', z_corner(1), comment, status)
    ! call ftgkyd(unit, 'y0', z_corner(2), comment, status)
    ! call ftgkyd(unit, 'z0', z_corner(3), comment, status)
    ! call ftgkyd(unit, 'margin', margin, comment, status)
    ! if (present(rank)) call ftgkyj(unit, 'rank', rank, comment, status)
    ! if (present(in_file)) call ftgkls(unit, 'dat file', in_file, comment, status)
    ! if (present(wt_file)) call ftgkls(unit, 'wts file', wt_file, comment, status)
    call ftclos(unit, status)
    call ftfiou(unit, status)
    allocate(pix2(1:size(pix),1:size(pix)))
    pix2 = real(pix, kind=rk)
  END SUBROUTINE read_fits_to_matrix_int

  SUBROUTINE read_fits_to_cube(filename,pix)
    implicit none
    INTEGER                                     :: status, unit, blocksize, rw, nfound, group, fpixel
    CHARACTER(len=*),intent(in):: filename ! fits file name
    REAL(rk),DIMENSION(:,:,:),ALLOCATABLE,intent(out):: pix
    !REAL(rk),DIMENSION(:,:),ALLOCATABLE:: dpix
    INTEGER,dimension(3):: pix_dim
    real(rk)                                    :: nullval
    logical                                     :: anynull
    status = 0
    call ftgiou(unit, status) ! get io number for fits file....seda ei pea tegelikult kasutama...voib ka ise maarata
!print*, filename
    rw = 0 ! 0-readonly, 1-readwrite
    ! blocksize -- ylearune asi...sellele ei tasu tahelepanu poorata.
    call ftopen(unit, filename, rw, blocksize, status)
    ! Get a sequence of numbered keyword values:
    call ftgknj(unit,'NAXIS', 1, 3, pix_dim, nfound, status)
    ! allocate memory
    ALLOCATE(pix(1:pix_dim(1),1:pix_dim(2),1:pix_dim(3)))
!print*, pix_dim
    group=1
    fpixel=1
    nullval = -999.0
    IF (rk==rdp) THEN
       CALL ftgpvd(unit, group, fpixel, SIZE(pix), nullval, pix, anynull, status)
    else if (rk==rsp) then
       call ftgpve(unit, group, fpixel, size(pix), nullval, pix, anynull, status)
    else
       PRINT*, "ERROR: no subroutine fitsop: ftgpvx (rk do not mach!)"
    END IF
    call ftclos(unit, status)
    call ftfiou(unit, status)
  END SUBROUTINE read_fits_to_cube

  SUBROUTINE write_matrix_to_fits(pix,filename)
    use ifport
    implicit none
    CHARACTER(len=*),intent(in):: filename ! fits file name
    REAL(rk),DIMENSION(:,:),intent(in):: pix
    REAL(rsp),DIMENSION(:,:),allocatable:: pix_s
    INTEGER                                     :: status, unit, blocksize, bitpix, group, fpixel,rw, errnum
    INTEGER,dimension(2):: pix_dim
    logical                                     :: simple, extend
    ALLOCATE(pix_s(1:SIZE(pix,dim=1),1:SIZE(pix,dim=2)))
    pix_s=REAL(pix,kind=rsp)
    status = 0
    call ftgiou(unit, status)
    blocksize = 1
    status = SYSTEM("rm -f "//filename)
    If (status .eq. -1) then
       errnum = ierrno( )
       print *, 'Error ', errnum
    end if
    status=0

    call ftinit(unit, filename, blocksize, status)
    if (status==105) then ! fail eksisteerib
       rw=1
       status=0
       call ftopen(unit, filename, rw, blocksize, status)
       !call ftdopn(unit, filename, rw, status)
    end if
    simple = .true.
    !        bitpix = -64                                           ! double
    bitpix = -32                                                    ! single
    extend = .true.
    pix_dim(1)=SIZE(pix,dim=1)
    pix_dim(2)=SIZE(pix,dim=2)
    call ftphpr(unit, simple, bitpix, size(pix_dim), pix_dim, 0, 1, extend, status)
    group = 1
    fpixel = 1
    IF (rk==rdp) THEN
       !call ftpprd(unit, group, fpixel, size(pix), pix, status)        ! double (Peab vastama pix tyybile!!!)
       call ftppre(unit, group, fpixel, size(pix_s), pix_s, status)        ! single
    else if (rk==rsp) then
       call ftppre(unit, group, fpixel, size(pix), pix, status)        ! single
       PRINT*, "ERROR: no subroutine fitsop: ftpprx (rk do not mach!)"
    END IF

    ! call ftpkyd(unit, 'delta', delta, 5, 'Grid cell size', status)
    ! call ftpkyd(unit, 'scale', scale, 5, 'B3 kernel radius', status)
    ! call ftpkyd(unit, 'margin', margin, 5, 'Margin', status)
    ! call ftpkyd(unit, 'x0', z_corner(1), 10, 'Grid corner location', status)
    ! call ftpkyd(unit, 'y0', z_corner(2), 10, 'Grid corner location', status)
    ! call ftpkyd(unit, 'z0', z_corner(3), 10, 'Grid corner location', status)
    ! if (present(rank)) call ftpkyj(unit, 'rank', rank, 'A trous rank', status)
    ! if (present(in_file)) call ftpkls(unit, 'dat file', in_file, ' ', status)
    ! if (present(wt_file)) call ftpkls(unit, 'wts file', wt_file, ' ', status)

    call ftclos(unit, status)
    call ftfiou(unit, status)
    deallocate(pix_s)
  end subroutine write_matrix_to_fits

  SUBROUTINE write_cube_to_fits(pix,filename)
    use ifport
    implicit none
    CHARACTER(len=*),intent(in):: filename ! fits file name
    REAL(rk),DIMENSION(:,:,:),intent(in):: pix
    INTEGER                                     :: status, unit, blocksize, bitpix, group, fpixel,rw, errnum
    INTEGER,dimension(3):: pix_dim
    logical                                     :: simple, extend
    status = 0
    call ftgiou(unit, status)
    blocksize = 1
    status = SYSTEM("rm -f "//filename)
    If (status .eq. -1) then
       errnum = ierrno( )
       print *, 'Error ', errnum
    end if
    status=0

    call ftinit(unit, filename, blocksize, status)
    simple = .true.
    !        bitpix = -64                                           ! double
    bitpix = -32                                                    ! single
    extend = .true.
    pix_dim(1)=SIZE(pix,dim=1)
    pix_dim(2)=SIZE(pix,dim=2)
    pix_dim(3)=SIZE(pix,dim=3)
    call ftphpr(unit, simple, bitpix, size(pix_dim), pix_dim, 0, 1, extend, status)
    group = 1
    fpixel = 1
    IF (rk==rdp) THEN
       call ftpprd(unit, group, fpixel, size(pix), pix, status)        ! double (Peab vastama pix tyybile!!!)
       !call ftppre(unit, group, fpixel, size(pix), real(pix,kind=rsp), status)        ! single
    else if (rk==rsp) then
       call ftppre(unit, group, fpixel, size(pix), pix, status)        ! single
       PRINT*, "ERROR: no subroutine fitsop: ftpprx (rk do not mach!)"
    END IF

    call ftclos(unit, status)
    call ftfiou(unit, status)
  end subroutine write_cube_to_fits

  SUBROUTINE write_cube_to_fits_sp(pix,filename)
    use ifport
    implicit none
    CHARACTER(len=*),intent(in):: filename ! fits file name
    REAL(rsp),DIMENSION(:,:,:),intent(in):: pix
    INTEGER                                     :: status, unit, blocksize, bitpix, group, fpixel,rw, errnum
    INTEGER,dimension(3):: pix_dim
    logical                                     :: simple, extend
    status = 0
    call ftgiou(unit, status)
    blocksize = 1
    status = SYSTEM("rm -f "//filename)
    If (status .eq. -1) then
       errnum = ierrno( )
       print *, 'Error ', errnum
    end if
    status=0

    call ftinit(unit, filename, blocksize, status)
    simple = .true.
    !        bitpix = -64                                           ! double
    bitpix = -32                                                    ! single
    extend = .true.
    pix_dim(1)=SIZE(pix,dim=1)
    pix_dim(2)=SIZE(pix,dim=2)
    pix_dim(3)=SIZE(pix,dim=3)
    call ftphpr(unit, simple, bitpix, size(pix_dim), pix_dim, 0, 1, extend, status)
    group = 1
    fpixel = 1
    IF (rk==rdp) THEN
       !call ftpprd(unit, group, fpixel, size(pix), pix, status)        ! double (Peab vastama pix tyybile!!!)
       call ftppre(unit, group, fpixel, size(pix), real(pix,kind=rsp), status)        ! single
    else if (rk==rsp) then
       call ftppre(unit, group, fpixel, size(pix), pix, status)        ! single
       !PRINT*, "ERROR: no subroutine fitsop: ftpprx (rk do not mach!)"
    END IF

    call ftclos(unit, status)
    call ftfiou(unit, status)
  end subroutine write_cube_to_fits_sp
  !
  subroutine read_fits_header_real(filename, keyword, vastus)
    	use ifport
  	implicit none
      	CHARACTER(len=*),intent(in):: filename, keyword ! fits file name
  	character(len=80)          :: comment
    	INTEGER                    :: status=0, unit , blocksize, rw=0, card=0!, bitpix, group, fpixel,rw, errnum
     	logical                    :: simple=.true., extend
  	real(rk) :: vastus
      	call ftgiou(unit, status)				!vaatab, kas unit vaba
  	call FTOPEN(unit, filename,rw, blocksize, status)	!avab read-only
  !	call FTGKY[EDJKLS](unit,keyword, > keyval,comment,status)
  	call FTGKYD(unit,keyword, vastus,comment,status)
  	call FTCLOS(unit, status)				!kinni panek
    end subroutine read_fits_header_real

    subroutine read_fits_header_integer(filename, keyword, vastus)
    	use ifport
  	implicit none
      	CHARACTER(len=*),intent(in):: filename, keyword ! fits file name
  	character(len=80)          :: comment
    	INTEGER                    :: status=0, unit , blocksize, rw=0, card=0, vastus!, bitpix, group, fpixel,rw, errnum
     	logical                    :: simple=.true., extend
      	call ftgiou(unit, status)				!vaatab, kas unit vaba
  	call FTOPEN(unit, filename,rw, blocksize, status)	!avab read-only
  !	call FTGKY[EDJKLS](unit,keyword, > keyval,comment,status)
  	call FTGKYK(unit,keyword, vastus,comment,status)
  	call FTCLOS(unit, status)				!kinni panek
    end subroutine read_fits_header_integer
    !
    subroutine write_fits_header_string(filename, keyword, sisu)
	use ifport
  	implicit none
      	CHARACTER(len=*),intent(in):: filename, keyword ! fits file name
	character(len=*), intent(in) 	   :: sisu !tekst, mida kirjtatakse 
  	character(len=80)          :: comment
    	INTEGER                    :: status=0, unit , blocksize, rw=1, card=0 !, bitpix, group, fpixel,rw, errnum
     	logical                    :: simple=.true., extend
      	call ftgiou(unit, status)				!vaatab, kas unit vaba
  	call FTOPEN(unit, filename,rw, blocksize, status)	!avab 
    	call FTPKYS(unit, keyword, sisu, comment, status) 	
	call FTCLOS(unit, status)
    end subroutine write_fits_header_string
	!
    subroutine write_fits_header_real(filename, keyword, sisu, decimals)
	use ifport
  	implicit none
      	CHARACTER(len=*),intent(in):: filename, keyword ! fits file name
	real(rk), intent(in) 	   :: sisu !tekst, mida kirjtatakse 
	integer, optional, intent(in) :: decimals !komakohtade arv, mida kirjutatakse faili
	integer 			:: dec
  	character(len=80)          :: comment
    	INTEGER                    :: status=0, unit , blocksize, rw=1, card=0 !, bitpix, group, fpixel,rw, errnum
     	logical                    :: simple=.true., extend
	if(.not.present(decimals)) then
		dec = 4
	else 
		dec = decimals
	end if
      	call ftgiou(unit, status)				!vaatab, kas unit vaba
  	call FTOPEN(unit, filename,rw, blocksize, status)	!avab 
	comment=''
    	call FTPKYD(unit, keyword, sisu, dec,  comment, status) 	
	call FTCLOS(unit, status)
    end subroutine write_fits_header_real
	!
    subroutine write_fits_header_integer(filename, keyword, sisu)
	use ifport
  	implicit none
      	CHARACTER(len=*),intent(in):: filename, keyword ! fits file name
	integer, intent(in) 	   :: sisu !tekst, mida kirjtatakse 
  	character(len=80)          :: comment
    	INTEGER                    :: status=0, unit , blocksize, rw=1, card=0 !, bitpix, group, fpixel,rw, errnum
     	logical                    :: simple=.true., extend
      	call ftgiou(unit, status)				!vaatab, kas unit vaba
  	call FTOPEN(unit, filename,rw, blocksize, status)	!avab 
    	call FTPKYK(unit, keyword, sisu, comment, status) 
	call FTCLOS(unit, status)
    end subroutine write_fits_header_integer
end module fits_op
