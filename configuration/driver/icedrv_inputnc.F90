MODULE icedrv_inputnc   

    USE netcdf
    USE icedrv_constants, only: nu_diag
    USE icedrv_system
    USE icedrv_calendar, only: time_after_reference

    IMPLICIT NONE
    PRIVATE

    PUBLIC nc_read_var
    PUBLIC nc_read_dim
    PUBLIC nc_open_read
    PUBLIC string2date_nc
    PUBLIC nc_findtimes
    PUBLIC nc_getscalar
    PUBLIC nc_getscalar1d


    INTERFACE nc_read_var
        MODULE PROCEDURE   nc_read_var_1d_i        , nc_read_var_1d_d       &
                         , nc_read_var_2d_d        , nc_read_var_0d_d
    END INTERFACE

    CONTAINS

    SUBROUTINE nc_findtimes(kdate,kdate0,hreal_units,ptimes,ktimesid,pweights)
       !!----------------------------------------------------------------------
       !!                 ***  nc_findtimes  ***
       !!
       !! **Purpose  :   read current model date, dataset reference date, times
       !!                and units. Returns time levels that contain model date
       !!                and interpolation weights
       !!
       !!
       !! ** Action  :
       !!
       !! References :   G. De Cillis, first write, 2022-10
       !!----------------------------------------------------------------------
       IMPLICIT NONE
 
       !arguments
       INTEGER(KIND=8)   , INTENT(IN   )         ::      kdate(7)     ! date to search
       INTEGER(KIND=8)   , INTENT(IN   )         ::      kdate0(7)       ! reference date in dataset
       CHARACTER(len=10) , INTENT(IN   )         ::      hreal_units    ! dataset time units
       DOUBLE PRECISION  , INTENT(IN   )         ::      ptimes(:)       ! dataset times
       INTEGER           , INTENT(  OUT)         ::      ktimesid(2)     ! time level
       DOUBLE PRECISION  , INTENT(  OUT)         ::      pweights(2)     !weights for time interpolation
 
       !local
       DOUBLE PRECISION                  ::      zmodel_time
       DOUBLE PRECISION                  ::      ztime(2)
       INTEGER                           ::      ji
 
 
       !! convert in seconds after 0001/01/01
       zmodel_time       = time_after_reference(kdate)
 
       ktimesid = 0
       !DO ji = 1,SIZE(ptimes)
       !    print*,zmodel_time,time_after_reference(kdate0,hreal_units,ptimes(ji))
       !ENDDO
 
       !! scan dataset dates
       DO ji = 1,SIZE(ptimes)
          IF ( zmodel_time .LE. time_after_reference(kdate0,hreal_units,ptimes(ji)) ) THEN
             ktimesid(2) = ji
             ztime(2)  = time_after_reference(kdate0,hreal_units,ptimes(ji))
             EXIT
          END IF
       END DO
       ktimesid(1) = ktimesid(2) -1
 
       IF (ktimesid(1) .EQ. 0) THEN
           ktimesid(1) = ktimesid(2)
           ktimesid(2) = ktimesid(1) + 1
           ztime(2)  = time_after_reference(kdate0,hreal_units,ptimes(ktimesid(2)))
       ELSEIF (ktimesid(1) .LT. 0 )  THEN
           CALL handle_err(0,'model time not in dataset times','nc_findtimes')
       ENDIF
 
       ztime(1)  = time_after_reference(kdate0,hreal_units,ptimes(ktimesid(1)))
 
 
       pweights(1) = (ztime(2) - zmodel_time)/(ztime(2)-ztime(1))
       pweights(2) = (zmodel_time - ztime(1))/(ztime(2)-ztime(1))
 
 
       RETURN
 
    END SUBROUTINE nc_findtimes
!---------------------------------------------------------------------------------


    SUBROUTINE nc_getscalar(nc_ds, hvar, ktimeid, pweights, pvar)
       !!----------------------------------------------------------------------
       !!                 ***  nc_getscalar  ***
       !!
       !! **Purpose  :   read current model date, dataset reference date, times
       !!                and units. Returns time levels that contain model date
       !!                and interpolation weights
       !!
       !!
       !! ** Action  :
       !!
       !! References :   G. De Cillis, first write, 2022-10
       !!----------------------------------------------------------------------
       IMPLICIT NONE
 
       !arguments
       TYPE(nc_data),     INTENT(IN)         :: nc_ds       !netcdf dataset
       CHARACTER(LEN=20), INTENT(IN)         :: hvar        !Variable name
       INTEGER,           INTENT(IN)         :: ktimeid(2)  !Time indices 
       DOUBLE PRECISION,  INTENT(IN)         :: pweights(2) !time interpolation weights
       DOUBLE PRECISION,  INTENT(OUT)        :: pvar        !Interp variable value

       !local
       INTEGER          :: irange(2)
       INTEGER          :: ncid
       DOUBLE PRECISION :: zfill
       DOUBLE PRECISION :: zvar1(1,1), zvar2(1,1)

       ncid = nc_ds%file_id
       irange = (/nc_ds%id_x, nc_ds%id_y/)
       CALL nc_read_var(ncid, hvar, zvar1, irange, zfill, ktime=ktimeid(1))
       CALL nc_read_var(ncid, hvar, zvar2, irange, zfill, ktime=ktimeid(2))

       pvar = pweights(1)*zvar1(1,1) + pweights(2)*zvar2(1,1)

    END SUBROUTINE nc_getscalar

    SUBROUTINE nc_getscalar1d(nc_ds, hvar, ktimeid, pweights, pvar)
       !!----------------------------------------------------------------------
       !!                 ***  nc_getscalar  ***
       !!
       !! **Purpose  :   read current model date, dataset reference date, times
       !!                and units. Returns time levels that contain model date
       !!                and interpolation weights
       !!
       !!
       !! ** Action  :
       !!
       !! References :   G. De Cillis, first write, 2022-10
       !!----------------------------------------------------------------------
       IMPLICIT NONE
 
       !arguments
       TYPE(nc_data1d),     INTENT(IN)         :: nc_ds       !netcdf dataset
       CHARACTER(LEN=20), INTENT(IN)         :: hvar        !Variable name
       INTEGER,           INTENT(IN)         :: ktimeid(2)  !Time indices 
       DOUBLE PRECISION,  INTENT(IN)         :: pweights(2) !time interpolation weights
       DOUBLE PRECISION,  INTENT(OUT)        :: pvar        !Interp variable value

       !local
       INTEGER          :: irange(2)
       INTEGER          :: ncid
       DOUBLE PRECISION :: zfill
       DOUBLE PRECISION :: zvar1, zvar2

       ncid = nc_ds%file_id
       CALL nc_read_var(ncid, hvar, zvar1, ktime=ktimeid(1))
       CALL nc_read_var(ncid, hvar, zvar2, ktime=ktimeid(2))

       pvar = pweights(1)*zvar1 + pweights(2)*zvar2

    END SUBROUTINE nc_getscalar1d



    SUBROUTINE nc_open_read(kncid,hfile)
        !!----------------------------------------------------------------------
        !!                 subroutine nc_open_read
        !!
        !! **Purpose  :   open a netcdf file in read mode
        !!
        !! ** Method  :   opening a netcdf file
        !!
        !! ** Action  : - first action (call nf90_open function)
        !!----------------------------------------------------------------------
  
        IMPLICIT NONE
  
        !arguments
        INTEGER, INTENT( out)             ::      kncid   ! netcdf file id (return)
        CHARACTER*(*), INTENT(in )        ::      hfile   ! netcdf file name
  
        !local
        INTEGER ierr
        CHARACTER*(256)                   ::      ystring ! error message
  
        ierr = nf90_open(hfile, NF90_NOWRITE, kncid)
  
        IF( ierr .eq. nf90_noerr ) RETURN
  
        WRITE(ystring,*) "Error: nf90_open unable to open netcdf file: file not found"
        CALL handle_err(ierr,TRIM(ystring),'nc_open_read')
  
        RETURN

    END SUBROUTINE nc_open_read
    !!-------------------------------------------------------------------------

    SUBROUTINE nc_read_var_1d_i(kncid,hvarname,kvarval,ktime)
        !!----------------------------------------------------------------------
        !!                 ***  nc_read_var_1d_i  ***
        !!
        !! **Purpose  :   retrieve values of a given netcdf variable (integer one-dimensional)
        !!
        !! ** Method  :   call netcdf function
        !!
        !! ** Action  : - first action (gets the id of the netcdf variable)
        !!            : - second action (gets the values of the variable)
        !!----------------------------------------------------------------------
 
        IMPLICIT NONE
 
        !arguments
        INTEGER           , INTENT(in   )                 ::      kncid       ! netcdf file identifier
        CHARACTER(*)      , INTENT(in   )                 ::      hvarname    ! variable name
        INTEGER           , INTENT(inout)                 ::      kvarval(:)  ! variable values
        INTEGER           , INTENT(in   ),OPTIONAL        ::      ktime       ! time level
 
        !local
        INTEGER                   ::      ivarid, ierr
        INTEGER                   ::      itype, indims
        INTEGER                   ::      idimids(2)
        CHARACTER(len=80)         ::      yvarname    ! variable name
 
        ierr = nf90_inq_varid(kncid,hvarname,ivarid)
        IF (ierr /= nf90_noerr) CALL handle_err(ierr,'ERROR: nf90_inq_varid generated an error:'//nf90_strerror(ierr),'nc_read_var_1d_i')
 
        ierr = nf90_inquire_variable(kncid, ivarid,yvarname, itype, indims, idimids)
        IF (ierr /= nf90_noerr) CALL handle_err(ierr,'ERROR: nf90_inquire_variable generated an error','nc_read_var_1d_i')

        IF(itype.LT.4) CALL handle_err(1,'ERROR: the type of data is different, ','nc_read_var_1d_i')
 
        IF (.NOT. PRESENT(ktime)) THEN
           IF(indims.ne.1) CALL handle_err(2,'ERROR: the number of dimension is greater than 1','nc_read_var_1d_i')
           ierr =  nf90_get_var(kncid, ivarid, kvarval(:))
           IF (ierr /= nf90_noerr) CALL handle_err(ierr,'ERROR: nf90_get_var generated an error','nc_read_var_1d_i')
        ELSE
           IF(indims.ne.2) CALL handle_err(2,'ERROR: the number of dimension is greater than 2','nc_read_var_1d_i')
           ierr =  nf90_get_var(kncid, ivarid, kvarval(:), start=(/1,ktime/))
           IF (ierr /= nf90_noerr) CALL handle_err(ierr,'ERROR: nf90_get_var generated an error','nc_read_var_1d_i')
        END IF
 
        RETURN
 
    END SUBROUTINE nc_read_var_1d_i

    SUBROUTINE nc_read_var_0d_d(kncid,hvarname,pvarval,ktime)
        !!----------------------------------------------------------------------
        !!                 ***  nc_read_var_1d_d  ***
        !!
        !! **Purpose  :   retrieve values of a given netcdf variable (double one-dimensional)
        !!
        !! ** Method  :   call netcdf function
        !!
        !! ** Action  : - first action (gets the id of the netcdf variable)
        !!            : - second action (gets the values of the variable)
        !!----------------------------------------------------------------------
 
        IMPLICIT NONE
 
        !arguments
        INTEGER           , INTENT(in   )                 ::      kncid       ! netcdf file identifier
        CHARACTER(*)      , INTENT(in   )                 ::      hvarname    ! variable name
        DOUBLE PRECISION  , INTENT(inout)                 ::      pvarval  ! variable values
        INTEGER           , INTENT(in   ),OPTIONAL        ::      ktime       ! time level
 
        !local
        INTEGER                   ::      ivarid, ierr
        INTEGER                   ::      itype, indims
        INTEGER                   ::      idimids(2)
        CHARACTER(len=80)         ::      yvarname    ! variable name
        DOUBLE PRECISION          ::      zfillval    !: variable FillValue
 
        ierr = nf90_inq_varid(kncid,hvarname,ivarid)
 
        IF (ierr /= nf90_noerr) CALL handle_err(ierr,'ERROR: nf90_inq_varid generated an error: '//nf90_strerror(ierr),'nc_read_var_1d_d')
 
        ierr = nf90_inquire_variable(kncid, ivarid,yvarname, itype, indims, idimids)
        IF (ierr /= nf90_noerr) CALL handle_err(ierr,'ERROR: nf90_inquire_variable generated an error','nc_read_var_1d_d')
 
        IF(itype.LT.5) CALL handle_err(1,'ERROR: the type of data is different, ','nc_read_var_1d_d')
 
        !! retrieve fillvalue
        ierr = nf90_get_att(kncid,ivarid,"_FillValue",zfillval)
        IF(ierr.NE.nf90_noerr) zfillval = -9.D+33
 
        IF (.NOT. PRESENT(ktime)) THEN
           IF(indims.ne.1) CALL handle_err(2,'ERROR: the number of dimension is greater than 1','nc_read_var_0d_d')
           ierr =  nf90_get_var(kncid, ivarid, pvarval)
           IF (ierr /= nf90_noerr) CALL handle_err(ierr,'ERROR: nf90_get_var generated an error','nc_read_var_0d_d')
        ELSE
           IF(indims.ne.1) CALL handle_err(2,'ERROR: the number of dimension is greater than 1','nc_read_var_0d_d')
           ierr =  nf90_get_var(kncid, ivarid, pvarval, start=(/ktime/))
           IF (ierr /= nf90_noerr) CALL handle_err(ierr,'ERROR: nf90_get_var generated an error','nc_read_var_0d_d')
        END IF
 
        !WHERE(ISNAN(pvarval)) pvarval = zfillval
 
        RETURN
 
    END SUBROUTINE nc_read_var_0d_d


    SUBROUTINE nc_read_var_1d_d(kncid,hvarname,pvarval,ktime)
        !!----------------------------------------------------------------------
        !!                 ***  nc_read_var_1d_d  ***
        !!
        !! **Purpose  :   retrieve values of a given netcdf variable (double one-dimensional)
        !!
        !! ** Method  :   call netcdf function
        !!
        !! ** Action  : - first action (gets the id of the netcdf variable)
        !!            : - second action (gets the values of the variable)
        !!----------------------------------------------------------------------
 
        IMPLICIT NONE
 
        !arguments
        INTEGER           , INTENT(in   )                 ::      kncid       ! netcdf file identifier
        CHARACTER(*)      , INTENT(in   )                 ::      hvarname    ! variable name
        DOUBLE PRECISION  , INTENT(inout)                 ::      pvarval(:)  ! variable values
        INTEGER           , INTENT(in   ),OPTIONAL        ::      ktime       ! time level
 
        !local
        INTEGER                   ::      ivarid, ierr
        INTEGER                   ::      itype, indims
        INTEGER                   ::      idimids(2)
        CHARACTER(len=80)         ::      yvarname    ! variable name
        DOUBLE PRECISION          ::      zfillval    !: variable FillValue
 
        ierr = nf90_inq_varid(kncid,hvarname,ivarid)
 
        IF (ierr /= nf90_noerr) CALL handle_err(ierr,'ERROR: nf90_inq_varid generated an error: '//nf90_strerror(ierr),'nc_read_var_1d_d')
 
        ierr = nf90_inquire_variable(kncid, ivarid,yvarname, itype, indims, idimids)
        IF (ierr /= nf90_noerr) CALL handle_err(ierr,'ERROR: nf90_inquire_variable generated an error','nc_read_var_1d_d')
 
        IF(itype.LT.5) CALL handle_err(1,'ERROR: the type of data is different, ','nc_read_var_1d_d')
 
        !! retrieve fillvalue
        ierr = nf90_get_att(kncid,ivarid,"_FillValue",zfillval)
        IF(ierr.NE.nf90_noerr) zfillval = -9.D+33
 
        IF (.NOT. PRESENT(ktime)) THEN
           IF(indims.ne.1) CALL handle_err(2,'ERROR: the number of dimension is greater than 1','nc_read_var_1d_d')
           ierr =  nf90_get_var(kncid, ivarid, pvarval(:))
           IF (ierr /= nf90_noerr) CALL handle_err(ierr,'ERROR: nf90_get_var generated an error','nc_read_var_1d_d')
        ELSE
           IF(indims.ne.2) CALL handle_err(2,'ERROR: the number of dimension is greater than 2','nc_read_var_1d_d')
           ierr =  nf90_get_var(kncid, ivarid, pvarval(:), start=(/1,ktime/))
           IF (ierr /= nf90_noerr) CALL handle_err(ierr,'ERROR: nf90_get_var generated an error','nc_read_var_1d_d')
        END IF
 
        !WHERE(ISNAN(pvarval)) pvarval = zfillval
 
        RETURN
 
    END SUBROUTINE nc_read_var_1d_d

    SUBROUTINE nc_read_var_2d_d(kncid,hvarname,pvarval,irange,pfillval,ktime)

        !!----------------------------------------------------------------------
        !!                 ***  nc_read_var_d ***
        !!
        !! **Purpose  :   retrieve value of a given netcdf variable at one point 
        !!                from regular datasets
        !!
        !! ** Method  :   call netcdf function
        !!
        !! ** Action  : - first action (gets the id of the netcdf variable)
        !!            : - second action (gets the values of the variable)
        !!
        !!----------------------------------------------------------------------
        
        IMPLICIT NONE
        
        !arguments
        INTEGER           , INTENT(in   )                 ::      kncid       !: netcdf file identifier
        CHARACTER(*)      , INTENT(in   )                 ::      hvarname    !: variable name
        DOUBLE PRECISION  , INTENT(inout)                 ::      pvarval(:,:)!: variable values
        INTEGER           , INTENT(in   )                 ::      irange(2)   !: start for x,y
        DOUBLE PRECISION  , INTENT(  out)                 ::      pfillval    !: variable FillValue
        INTEGER           , INTENT(in   ),OPTIONAL        ::      ktime       !: time level
        
        !local
        INTEGER                   ::      ivarid, ierr
        INTEGER                   ::      itype, indims
        INTEGER                   ::      idimids(3)
        CHARACTER(len=80)         ::      yvarname    ! variable name
        DOUBLE PRECISION          ::      zscale,zscale_nc,zoffset,zoffset_nc
        CHARACTER(len=*), PARAMETER :: subname='(nc_read_var_2d_d)'
        
        ierr = nf90_inq_varid(kncid,hvarname,ivarid)
        IF (ierr /= nf90_noerr) CALL handle_err(ierr,'ERROR: nf90_inq_varid generated an error -'//TRIM(hvarname),subname)
        
        ierr = nf90_inquire_variable(kncid, ivarid,yvarname, itype, indims, idimids)
        IF (ierr /= nf90_noerr) CALL handle_err(ierr,'ERROR: nf90_inquire_variable generated an error',subname)
        
        !! retrieve fillvalue
        ierr = nf90_get_att(kncid,ivarid,"_FillValue",pfillval)
        IF(ierr.NE.nf90_noerr) then
            IF(ierr.EQ.nf90_enotatt) then 
                pfillval = -9.D+33
            ELSE
                WRITE(nu_diag,*) nf90_strerror(ierr)
            ENDIF
        ENDIF
        
        !! retrieve offset (if any)
        zoffset   = 0.0D0
        ierr = nf90_get_att(kncid,ivarid,"add_offset",zoffset_nc)
        
        IF ((ierr /= nf90_noerr).AND.(ierr /= nf90_enotatt)) WRITE(nu_diag,*) nf90_strerror(ierr)
        
        IF ( ierr .EQ. 0 ) zoffset = zoffset_nc
        
        !! retrieve scale_factor (if any)
        zscale    = 1.0D0
        ierr = nf90_get_att(kncid,ivarid,"scale_factor",zscale_nc)
        
        IF ((ierr /= nf90_noerr) .AND. (ierr /= nf90_enotatt)) WRITE(nu_diag,*) nf90_strerror(ierr)

        IF ( ierr .EQ. 0 ) zscale = zscale_nc
        
        
        IF(itype.LT.3) CALL handle_err(1,'ERROR: the type of data is different, ', subname) 
        
        IF (.NOT. PRESENT(ktime)) THEN
           IF(indims.ne.2) CALL handle_err(2, 'ERROR: the number of dimension is different than 2', subname)

           ierr =  nf90_get_var(kncid, ivarid, pvarval(:,:), start=(/irange(1),irange(2)/), count=(/1,1/)  )
           IF (ierr /= nf90_noerr) CALL handle_err(ierr,'ERROR: nf90_get_var generated an error: '//nf90_strerror(ierr), subname) 
        ELSE
           IF(indims.ne.3) CALL handle_err(2,'ERROR: the number of dimension is different than 3', subname) 

           ierr =  nf90_get_var(kncid, ivarid, pvarval(:,:), start=(/irange(1),irange(2),ktime/), count=(/1,1,1/))
           IF (ierr /= nf90_noerr) CALL handle_err(ierr,'ERROR: nf90_get_var generated an error: '//nf90_strerror(ierr), subname) 
        END IF
        
        IF ( zscale .NE. 1.0D0 .OR. zoffset .NE. 0.0D0 ) THEN
           WHERE ( pvarval .NE. pfillval ) pvarval = pvarval * zscale + zoffset
        END IF

        
        RETURN

    END SUBROUTINE nc_read_var_2d_d
    !!----------------------------------------------------------------------
    

    SUBROUTINE nc_read_dim(kncid,hdimname,kdimval)
        !!----------------------------------------------------------------------
        !!                 ***  nc_read_dim  ***
        !!
        !! **Purpose  :   retrieve value of a given netcdf dimension
        !!
        !! ** Method  :   call netcdf function
        !!
        !! ** Action  : - first action (gets the id of the netcdf dimension)
        !!            : - second action (gets the value of the dimension)
        !!----------------------------------------------------------------------
 
        IMPLICIT NONE
 
        !arguments
        INTEGER, INTENT(in )              ::      kncid           ! netcdf file identifier
        CHARACTER(*), INTENT(in )         ::      hdimname        ! dimension name
        INTEGER, INTENT( out)             ::      kdimval         ! dimension value (return)
 
        !local
        INTEGER idimid, ierr
 
        ierr = nf90_inq_dimid(kncid,hdimname,idimid)
        IF (ierr /= nf90_noerr) CALL handle_err(ierr,hdimname//TRIM(nf90_strerror(ierr)),'nc_read_dim')
 
        ierr = nf90_inquire_dimension(kncid,idimid,len=kdimval)
        IF (ierr /= nf90_noerr) CALL handle_err(ierr,'nf90_inquire_dimension generated an error','nc_read_dim: 2')
 
        RETURN
 
    END SUBROUTINE nc_read_dim

   !!----------------------------------------------------------------------

   FUNCTION string2date_nc(hstring)
      !!----------------------------------------------------------------------
      !!                 ***  string2date_nc  ***
      !!
      !! **Purpose  :   search for a date (YYYY-m(m)-D(D) HH:MM:SS) in a units string
      !!                from netcdf dataset
      !!
      !! ** Method  :   use IACHAR function
      !!
      !! ** Action  : -
      !!            : -
      !!----------------------------------------------------------------------
      IMPLICIT NONE

      CHARACTER(LEN=*)  , INTENT(in   )         :: hstring

      INTEGER(KIND=8)                           :: idatehour(7)
      INTEGER(KIND=8)                           :: string2date_nc(7)

      !! LOCAL VARIABLES
      INTEGER                                   :: inum, ji,istat,icount
      INTEGER                                   :: js
      CHARACTER(len=4)                          :: ystring_y
      CHARACTER(len=2)                          :: ystring_m
      CHARACTER(len=2)                          :: ystring_d
      CHARACTER(len=2)                          :: ystring_md
      CHARACTER(len=6)                          :: ystring_hms
      CHARACTER(len=8)                          :: ystring_ymd
      INTEGER                                   :: kymd, khms
      INTEGER                                   :: kmd

      icount     = 0

      !! loop on input string
      ! READ YEAR
      DO ji = 1, LEN_TRIM(hstring)
         READ(hstring(ji:ji),*,iostat=istat) inum
         IF (istat .EQ. 0 ) THEN
             icount = icount + 1
             ystring_y(icount:icount) = hstring(ji:ji)
         END IF
         IF(icount .EQ. 4) THEN
             js = ji + 1
             EXIT
         ENDIF
      END DO


      !READ MONTH
      icount = 0
      ystring_md = " "
      DO ji = js, LEN_TRIM(hstring)
         READ(hstring(ji:ji),*,iostat=istat) inum
         IF (istat .EQ. 0 ) THEN
             icount = icount + 1
             ystring_md(icount:icount) = hstring(ji:ji)
         ELSE
             IF (icount .GE. 1) THEN
                 READ(ystring_md,*)        kmd
                 WRITE (ystring_m,FMT='(I0.2)') kmd
                 js = ji + 1
                 EXIT
             ENDIF
         END IF
      END DO

      !READ DAY
      icount = 0
      ystring_md = " "
      DO ji = js, LEN_TRIM(hstring)
         READ(hstring(ji:ji),*,iostat=istat) inum
         IF (istat .EQ. 0 ) THEN
             icount = icount + 1
             ystring_md(icount:icount) = hstring(ji:ji)
         ELSE
             IF (icount .GE. 1) THEN
                 READ(ystring_md,*)       kmd
                 WRITE (ystring_d,FMT='(I0.2)') kmd
                 js = ji + 1
                 EXIT
             ENDIF
         END IF
      END DO

      IF (ji.GT.LEN_TRIM(hstring)) THEN
          !This happens if hours:minutes:seconds are not provided in time units
          READ(ystring_md,*)       kmd
          WRITE (ystring_d,FMT='(I0.2)') kmd
          js = ji
          ystring_hms = "000000"  !default 
      ENDIF


      ! READ HH:MM:SS
      icount = 0
      DO ji = js, LEN_TRIM(hstring)
         READ(hstring(ji:ji),*,iostat=istat) inum
         IF (istat .EQ. 0 ) THEN
            icount = icount + 1
            ystring_hms(icount:icount) = hstring(ji:ji)
         END IF
         IF (icount.EQ.6) EXIT    !Milliseconds not read
      END DO

      ystring_ymd(1:4)  = ystring_y
      ystring_ymd(5:6)  = ystring_m
      ystring_ymd(7:8)  = ystring_d



      READ(ystring_ymd,'(I8)')       kymd
      READ(ystring_hms,'(I6)')       khms

      !! calculate components of date / time array

      idatehour(1)      = kymd / 10**4
      idatehour(2)      = (kymd - idatehour(1)*10**4) / 10**2
      idatehour(3)      = kymd - idatehour(1)*10**4 - idatehour(2)*10**2

      idatehour(4)      = khms / 10**4
      idatehour(5)      = (khms - idatehour(4)*10**4) / 10**2
      idatehour(6)      = khms - idatehour(4)*10**4 - idatehour(5)*10**2
      !! Milliseconds
      idatehour(7)      = 0

      !! DO SOME CHECKS
      IF ( idatehour(2) .GT. 12 ) CALL handle_err(1,'Wrong date month ','string2date_nc')
      IF ( idatehour(3) .GT. 31 ) CALL handle_err(2,'Wrong date day   ','string2date_nc')
      IF ( idatehour(4) .GT. 23 ) CALL handle_err(3,'Wrong date hour  ','string2date_nc')
      IF ( idatehour(5) .GT. 59 ) CALL handle_err(4,'Wrong date minute','string2date_nc')
      IF ( idatehour(6) .GT. 59 ) CALL handle_err(5,'Wrong date second','string2date_nc')

      IF ( .NOT. isleap(idatehour(1)) .AND. idatehour(2) .EQ. 2 .AND. idatehour(3) .GT. 28) THEN
         CALL handle_err(6,'Year is not leap, February cannot have > 28 days','string2date_nc')
      END IF

      string2date_nc = idatehour

      RETURN

   END FUNCTION string2date_nc

   FUNCTION isleap(kyear)
      !!----------------------------------------------------------------------
      !!                 ***  isleap  ***
      !!
      !! **Purpose  :   Return TRUE if kyear is leap
      !!
      !! ** Method  :
      !!
      !! ** Action  : -
      !!            : -
      !!----------------------------------------------------------------------

      LOGICAL           :: isleap
      INTEGER(KIND=8)   :: kyear

      isleap = .FALSE.

      IF(  ( mod(kyear,4) .eq. 0 .and. mod(kyear,100) .ne. 0 ) .or. mod(kyear,400) .eq. 0 ) isleap = .TRUE.

   END FUNCTION isleap


END MODULE icedrv_inputnc   
