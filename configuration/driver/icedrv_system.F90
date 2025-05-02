!=======================================================================

! Diagnostic information output during run
!
! author: Tony Craig

      module icedrv_system

      use icedrv_kinds
      use icedrv_constants, only: nu_diag
      use icedrv_state, only: aice
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted

      implicit none
      private
      public :: icedrv_system_abort
      public :: handle_err

      type, public :: nc_data
          character (char_len_long), public    ::    y_filename = " "
          character (10), public               ::    y_lonname = " "
          character (10), public               ::    y_latname = " "
          integer(kind=int_kind), public       ::    file_id   !netcdf id
          real(kind=dbl_kind), public          ::    zpoint_lon, zpoint_lat
          integer(kind=int_kind), public       ::    id_x      !point lon index
          integer(kind=int_kind), public       ::    id_y      !point lat index
          character (10), public               ::    y_tunits   !time units
          integer(kind=int8_kind), public       ::    idate_0(7)   !dataset reference date
          real(kind=dbl_kind), public, allocatable   ::    ztimes(:)    !: dataset times in time units (y_tunits) from ref date (idate_0)
           
      end type nc_data

      type, public :: nc_data1d
          character (char_len_long), public    ::    y_filename = " "
          integer(kind=int_kind), public       ::    file_id   !netcdf id
          character (10), public               ::    y_tunits   !time units
          integer(kind=int8_kind), public       ::    idate_0(7)   !dataset reference date
          real(kind=dbl_kind), public, allocatable   ::    ztimes(:)    !: dataset times in time units (y_tunits) from ref date (idate_0)
           
      end type nc_data1d


!=======================================================================

      contains

!=======================================================================
! prints error information prior to aborting

      subroutine icedrv_system_abort(icell, istep, string, file, line)

      integer (kind=int_kind), intent(in), optional :: &
         icell       , & ! indices of grid cell where model aborts
         istep       , & ! time step number
         line            ! line number

      character (len=*), intent(in), optional :: string, file

      ! local variables

      character(len=*), parameter :: subname='(icedrv_system_abort)'

      write(nu_diag,*) ' '

      call icepack_warnings_flush(nu_diag)

      write(nu_diag,*) ' '
      write(nu_diag,*) subname,' ABORTED: '
      if (present(file))   write (nu_diag,*) subname,' called from ',trim(file)
      if (present(line))   write (nu_diag,*) subname,' line number ',line
      if (present(istep))  write (nu_diag,*) subname,' istep =', istep
      if (present(icell))  write (nu_diag,*) subname,' i, aice =', icell, aice(icell)
      if (present(string)) write (nu_diag,*) subname,' string = ',trim(string)
      stop

      end subroutine icedrv_system_abort

!=======================================================================

    SUBROUTINE handle_err(kerr,herrstr,hsubname)
       !!----------------------------------------------------------------------
       !!                 ***  handle_err  ***
       !!
       !! **Purpose  :   print a error string that allows identification of
       !!                the error source
       !!
       !! ** Method  :   write a string and stop the program
       !!
       !! ** Action  : - first action
       !!----------------------------------------------------------------------
 
       IMPLICIT NONE
 
       !arguments
       INTEGER, INTENT(in )             ::      kerr            ! error id
       CHARACTER(*), INTENT(in )        ::      herrstr       ! error string
       CHARACTER(*), INTENT(in )        ::      hsubname        ! subroutine or function name
 
       WRITE(nu_diag,*) herrstr,' ', hsubname
       WRITE(nu_diag,*) 'Error code: ', kerr
 
       STOP ' Error stop handle_err'
 
    END SUBROUTINE handle_err

      end module icedrv_system

!=======================================================================
