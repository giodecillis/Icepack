!=======================================================================
!
! Reads and interpolates forcing data for atmosphere and ocean quantities.
!
! authors: Elizabeth C. Hunke and William H. Lipscomb, LANL
!
      module icedrv_forcing

      use icedrv_kinds
      use icedrv_domain_size, only: nx
      use icedrv_calendar, only: time, nyr, dayyr, mday, month, hour, secday
      use icedrv_calendar, only: daymo, daycal, dt, yday, sec
      use icedrv_calendar, only: get_model_date 
      use icedrv_constants, only: nu_diag, nu_forcing, nu_open_clos
      use icedrv_constants, only: c0, c1, c2, c10, c100, p5, c4, c24
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_query_parameters
      use icepack_intfc, only: icepack_sea_freezing_temperature
      use icepack_intfc, only: icepack_init_wave
      use icedrv_system
      use icedrv_flux, only: zlvl, Tair, potT, rhoa, uatm, vatm, wind, &
         strax, stray, fsw, swvdr, swvdf, swidr, swidf, Qa, flw, frain, &
         fsnow, sst, sss, uocn, vocn, qdp, hmix, Tf, opening, closing, sstdat, frzmlt
      use icedrv_inputnc
      use icedrv_phy
      use icedrv_state, only: aice

      implicit none
      private
      public :: init_forcing, get_forcing, interp_coeff, &
                interp_coeff_monthly, get_wave_spec

      integer (kind=int_kind), parameter :: &
         ntime = 8760        ! number of data points in time

      integer (kind=int_kind), public :: &
         ycycle          , & ! number of years in forcing cycle
         fyear_init      , & ! first year of data in forcing cycle
         fyear           , & ! current year in forcing cycle
         fyear_final         ! last year in cycle

      real (kind=dbl_kind), dimension(ntime) :: &
            fsw_data, & ! field values at temporal data points
           cldf_data, &
          fsnow_data, &
           Tair_data, &
           uatm_data, &
           vatm_data, &
           wind_data, &
          strax_data, &
          stray_data, &
           rhum_data, &
             Qa_data, &
           rhoa_data, &
           potT_data, &
            flw_data, &
            qdp_data, &
            sst_data, &
            sss_data, &
           uocn_data, &
           vocn_data, &
          frain_data, &
          swvdr_data, &
          swvdf_data, &
          swidr_data, &
          swidf_data, &
           zlvl_data, &
           hmix_data

      real (kind=dbl_kind), dimension(nx) :: &
          sst_temp

      real (kind=dbl_kind), dimension(ntime) :: &
           open_data, &
           clos_data

      character(char_len), public :: &
         atm_data_format, & ! 'bin'=binary or 'nc'=netcdf
         ocn_data_format, & ! 'bin'=binary or 'nc'=netcdf
         bgc_data_format, & ! 'bin'=binary or 'nc'=netcdf
         atm_data_type,   & ! 'default', 'clim', 'CFS', 'ECMWF'
         ocn_data_type,   & ! 'default', 'SHEBA'
         bgc_data_type,   & ! 'default', 'ISPOL', 'NICE'
         lateral_flux_type,   & ! 'uniform_ice', 'open_water'
         atm_data_file,   & ! atmospheric forcing data file
         atm_pert_file,   & ! atmospheric forcing perturbation file
         ocn_data_file,   & ! ocean forcing data file
         ice_data_file,   & ! ice forcing data file
         bgc_data_file,   & ! biogeochemistry forcing data file
         precip_units       ! 'mm_per_month', 'mm_per_sec', 'mks'

      real (kind=dbl_kind), public :: &
         atm_longitude,   & ! atmospheric longitude for ECMWF
         atm_latitude       ! atmospheric latitude for ECMWF


      character(char_len_long), public :: &
         data_dir           ! top directory for forcing data

      real (kind=dbl_kind), parameter, public :: &
         frcvdr = 0.28_dbl_kind, & ! frac of incoming sw in vis direct band
         frcvdf = 0.24_dbl_kind, & ! frac of incoming sw in vis diffuse band
         frcidr = 0.31_dbl_kind, & ! frac of incoming sw in near IR direct band
         frcidf = 0.17_dbl_kind    ! frac of incoming sw in near IR diffuse band

      logical (kind=log_kind), public :: &
         oceanmixed_ice , & ! if true, use internal ocean mixed layer
         restore_ocn        ! restore sst if true

      real (kind=dbl_kind), public :: &
         trest, &           ! restoring time scale (sec)
         trestore           ! restoring time scale (days)

      character (len=char_len_long), public :: &
         snw_ssp_table      ! snow table type 'test', 'snicar'

      type(nc_data), public :: &
         ecmwf_data, &         ! ecmwf netcdf dataset 
         cglors_data           ! cglors netcdf dataset 

      type(nc_data1d), public :: &
         pert_atm_data         ! pointwise atm pert file

!=======================================================================

      contains

!=======================================================================

      subroutine init_forcing

! Determine the current and final year of the forcing cycle based on
! namelist input; initialize the forcing data.

      integer (kind=int_kind) :: &
         i                ! index

      character(len=*), parameter :: subname='(init_forcing)'

      fyear       = fyear_init + mod(nyr-1,ycycle) ! current year
      fyear_final = fyear_init + ycycle - 1 ! last year in forcing cycle

      write (nu_diag,*) ' Initial forcing data year = ',fyear_init
      write (nu_diag,*) ' Final   forcing data year = ',fyear_final

     !-------------------------------------------------------------------
     ! Initialize forcing data to default values
     !-------------------------------------------------------------------

      ! many default forcing values are set in init_flux_atm
      i = 1 ! use first grid box value

          zlvl_data(:) = zlvl (i)    ! atmospheric data level (m)
          Tair_data(:) = Tair (i)    ! air temperature  (K)
          potT_data(:) = potT (i)    ! air potential temperature  (K)
          rhoa_data(:) = rhoa (i)    ! air density (kg/m^3)
          uatm_data(:) = uatm (i)    ! wind velocity components (m/s)
          vatm_data(:) = vatm (i)
          wind_data(:) = wind (i)    ! wind speed (m/s)
         strax_data(:) = strax(i)    ! wind stress components (N/m^2)
         stray_data(:) = stray(i)
           fsw_data(:) = fsw  (i)    ! incoming shortwave radiation (W/m^2)
         swvdr_data(:) = swvdr(i)    ! sw down, visible, direct  (W/m^2)
         swvdf_data(:) = swvdf(i)    ! sw down, visible, diffuse (W/m^2)
         swidr_data(:) = swidr(i)    ! sw down, near IR, direct  (W/m^2)
         swidf_data(:) = swidf(i)    ! sw down, near IR, diffuse (W/m^2)
            Qa_data(:) = Qa   (i)    ! specific humidity (kg/kg)
           flw_data(:) = flw  (i)    ! incoming longwave radiation (W/m^2)
         frain_data(:) = frain(i)    ! rainfall rate (kg/m^2 s)
         fsnow_data(:) = fsnow(i)    ! snowfall rate (kg/m^2 s)
           qdp_data(:) = qdp  (i)    ! deep ocean heat flux (W/m^2)
           sss_data(:) = sss  (i)    ! sea surface salinity
          uocn_data(:) = uocn (i)    ! ocean current components (m/s)
          vocn_data(:) = vocn (i)
          cldf_data(:) = c0          ! cloud fraction

      if (trim(atm_data_type(1:4)) == 'CFS')   call atm_CFS
      if (trim(atm_data_type(1:5)) == 'ECMWF') call atm_ECMWF(0)
      if (trim(atm_data_type(1:4)) == 'clim')  call atm_climatological
      if (trim(atm_data_type(1:5)) == 'ISPOL') call atm_ISPOL
      if (trim(atm_data_type(1:4)) == 'NICE')  call atm_NICE
      if (trim(ocn_data_type(1:5)) == 'SHEBA') call ice_open_clos

      if (trim(atm_pert_file)/='') call atm_pert(0)

      if (restore_ocn) then
        if (trestore == 0) then
          trest = dt        ! use data instantaneously
        else
          trest = real(trestore,kind=dbl_kind) * secday ! seconds
        end if
        sst_data(:) = sstdat(i)    ! default may be overwritten below
      else
        sst_data(:) = sst   (i)    ! default or restart value if not restoring
      endif

      if (trim(ocn_data_type(1:5)) == 'ISPOL') call ocn_ISPOL
      if (trim(ocn_data_type(1:4)) == 'NICE')  call ocn_NICE
      if (trim(ocn_data_type(1:4)) == 'CGLO')  call ocn_CGLO(0)

      call prepare_forcing (Tair_data,     fsw_data,      &
                            cldf_data,     &
                            frain_data,    fsnow_data,    &
                            Qa_data,       rhoa_data,     &
                            uatm_data,     vatm_data,     &
                            strax_data,    stray_data,    &
                            zlvl_data,     wind_data,     &
                            swvdr_data,    swvdf_data,    &
                            swidr_data,    swidf_data,    &
                            potT_data)

      end subroutine init_forcing

!=======================================================================

      subroutine get_forcing(timestep)

!ECH notes
! We will probably need to send in the time and work out what the data
! time slice is, instead of sending in the timestep.  This currently assumes
! the time step and the data both start Jan 1.

      integer (kind=int_kind), intent(in) :: &
         timestep         ! time step index

      integer (kind=int_kind) :: &
         i            , & ! data index
         recslot      , & ! spline slot for current record
         midmonth         ! middle day of month

      integer (kind=int_kind) :: &
         mlast, mnext     ! indices of bracketing time slices

      real (kind=dbl_kind) :: &
         c1intp, c2intp   ! interpolation coefficients

      integer (kind=int_kind) :: &
          recnum      , & ! record number for current data value
          maxrec      , & ! maximum number of data records
          dataloc     , & ! = 1 for data located in middle of time interval
                          ! = 2 for date located at end of time interval
          offndy          ! Julian day of first data record

      real (kind=dbl_kind) :: &
          sec6hr      , & ! number of seconds in 6 hours
          sec1hr      , & ! number of seconds in 1 hour
          offset          ! time to first data record since 1 Jan (s)

      character(len=*), parameter :: subname='(get_forcing)'

      if (trim(atm_data_type) == 'CFS') then
         ! calculate data index corresponding to current timestep
         i = mod(timestep-1,ntime)+1 ! repeat forcing cycle
         mlast = i
         mnext = mlast
         c1intp = c1
         c2intp = c0

         ! fill all grid boxes with the same forcing data
         Tair (:) = c1intp *  Tair_data(mlast) + c2intp *  Tair_data(mnext)
         Qa   (:) = c1intp *    Qa_data(mlast) + c2intp *    Qa_data(mnext)
         uatm (:) = c1intp *  uatm_data(mlast) + c2intp *  uatm_data(mnext)
         vatm (:) = c1intp *  vatm_data(mlast) + c2intp *  vatm_data(mnext)
         fsnow(:) = c1intp * fsnow_data(mlast) + c2intp * fsnow_data(mnext)
         flw  (:) = c1intp *   flw_data(mlast) + c2intp *   flw_data(mnext)
         fsw  (:) = c1intp *   fsw_data(mlast) + c2intp *   fsw_data(mnext)

         ! derived (or not otherwise set)
         potT (:) = c1intp *  potT_data(mlast) + c2intp *  potT_data(mnext)
         wind (:) = c1intp *  wind_data(mlast) + c2intp *  wind_data(mnext)
         strax(:) = c1intp * strax_data(mlast) + c2intp * strax_data(mnext)
         stray(:) = c1intp * stray_data(mlast) + c2intp * stray_data(mnext)
         rhoa (:) = c1intp *  rhoa_data(mlast) + c2intp *  rhoa_data(mnext)
         frain(:) = c1intp * frain_data(mlast) + c2intp * frain_data(mnext)
         swvdr(:) = c1intp * swvdr_data(mlast) + c2intp * swvdr_data(mnext)
         swvdf(:) = c1intp * swvdf_data(mlast) + c2intp * swvdf_data(mnext)
         swidr(:) = c1intp * swidr_data(mlast) + c2intp * swidr_data(mnext)
         swidf(:) = c1intp * swidf_data(mlast) + c2intp * swidf_data(mnext)

      elseif (trim(atm_data_type) == 'ECMWF') then
          CALL atm_ECMWF(1)

      elseif (trim(atm_data_type) == 'clim') then
         midmonth = 15  ! assume data is given on 15th of every month
         recslot = 1                             ! latter half of month
         if (mday < midmonth) recslot = 2        ! first half of month
         if (recslot == 1) then
            mlast = month
            mnext = mod(month   ,12) + 1
         else ! recslot = 2
            mlast = mod(month+10,12) + 1
            mnext = month
         endif
         call interp_coeff_monthly(recslot, c1intp, c2intp)

         ! fill all grid boxes with the same forcing data
         Tair (:) = c1intp *  Tair_data(mlast) + c2intp *  Tair_data(mnext)
         Qa   (:) = c1intp *    Qa_data(mlast) + c2intp *    Qa_data(mnext)
         uatm (:) = c1intp *  uatm_data(mlast) + c2intp *  uatm_data(mnext)
         vatm (:) = c1intp *  vatm_data(mlast) + c2intp *  vatm_data(mnext)
         fsnow(:) = c1intp * fsnow_data(mlast) + c2intp * fsnow_data(mnext)
         flw  (:) = c1intp *   flw_data(mlast) + c2intp *   flw_data(mnext)
         fsw  (:) = c1intp *   fsw_data(mlast) + c2intp *   fsw_data(mnext)

         ! derived (or not otherwise set)
         potT (:) = c1intp *  potT_data(mlast) + c2intp *  potT_data(mnext)
         wind (:) = c1intp *  wind_data(mlast) + c2intp *  wind_data(mnext)
         strax(:) = c1intp * strax_data(mlast) + c2intp * strax_data(mnext)
         stray(:) = c1intp * stray_data(mlast) + c2intp * stray_data(mnext)
         rhoa (:) = c1intp *  rhoa_data(mlast) + c2intp *  rhoa_data(mnext)
         frain(:) = c1intp * frain_data(mlast) + c2intp * frain_data(mnext)
         swvdr(:) = c1intp * swvdr_data(mlast) + c2intp * swvdr_data(mnext)
         swvdf(:) = c1intp * swvdf_data(mlast) + c2intp * swvdf_data(mnext)
         swidr(:) = c1intp * swidr_data(mlast) + c2intp * swidr_data(mnext)
         swidf(:) = c1intp * swidf_data(mlast) + c2intp * swidf_data(mnext)

         !call atm_ECMWF(1)

      elseif (trim(atm_data_type) == 'ISPOL') then

        offndy = 0                              ! first data record (Julian day)
        offset = real(offndy,dbl_kind)*secday
        dataloc = 1                             ! data located at middle of interval
        maxrec = 365
        recslot = 2
        recnum = mod(int(yday)+maxrec-offndy-1,maxrec)+1
        mlast = mod(recnum+maxrec-2,maxrec) + 1
        mnext = mod(recnum-1,       maxrec) + 1
        call interp_coeff (recnum, recslot, secday, dataloc, &
                           c1intp, c2intp, offset)

        Tair (:) = c1intp *  Tair_data(mlast) + c2intp *  Tair_data(mnext)
        Qa   (:) = c1intp *    Qa_data(mlast) + c2intp *    Qa_data(mnext)
        uatm (:) = c1intp *  uatm_data(mlast) + c2intp *  uatm_data(mnext)
        vatm (:) = c1intp *  vatm_data(mlast) + c2intp *  vatm_data(mnext)
        fsnow(:) = c1intp * fsnow_data(mlast) + c2intp * fsnow_data(mnext)

         ! derived (or not otherwise set)
         potT (:) = c1intp *  potT_data(mlast) + c2intp *  potT_data(mnext)
         wind (:) = c1intp *  wind_data(mlast) + c2intp *  wind_data(mnext)
         strax(:) = c1intp * strax_data(mlast) + c2intp * strax_data(mnext)
         stray(:) = c1intp * stray_data(mlast) + c2intp * stray_data(mnext)
         rhoa (:) = c1intp *  rhoa_data(mlast) + c2intp *  rhoa_data(mnext)
         frain(:) = c1intp * frain_data(mlast) + c2intp * frain_data(mnext)

        sec6hr = secday/c4;                      ! seconds in 6 hours
        offndy = 0
        maxrec = 1460
        recnum = 4*int(yday) - 3 + int(real(sec,kind=dbl_kind)/sec6hr)
        recnum = mod(recnum+maxrec-4*offndy-1,maxrec)+1 ! data begins on 16 June 2004
        recslot = 2
        mlast = mod(recnum+maxrec-2,maxrec) + 1
        mnext = mod(recnum-1,       maxrec) + 1
        call interp_coeff (recnum, recslot, sec6hr, dataloc, &
                           c1intp, c2intp, offset)

        fsw  (:) = c1intp *   fsw_data(mlast) + c2intp *   fsw_data(mnext)
        flw  (:) = c1intp *   flw_data(mlast) + c2intp *   flw_data(mnext)

         ! derived
         swvdr(:) = c1intp * swvdr_data(mlast) + c2intp * swvdr_data(mnext)
         swvdf(:) = c1intp * swvdf_data(mlast) + c2intp * swvdf_data(mnext)
         swidr(:) = c1intp * swidr_data(mlast) + c2intp * swidr_data(mnext)
         swidf(:) = c1intp * swidf_data(mlast) + c2intp * swidf_data(mnext)

      elseif (trim(atm_data_type) == 'NICE') then

        offndy = 0                              ! first data record (Julian day)
        offset = real(offndy,dbl_kind)*secday
        dataloc = 1                          ! data located in middle of interval
        maxrec = 365
        recslot = 2
        recnum = mod(int(yday)+maxrec-offndy-1,maxrec)+1
        mlast = mod(recnum+maxrec-2,maxrec) + 1
        mnext = mod(recnum-1,       maxrec) + 1
        call interp_coeff (recnum, recslot, secday, dataloc, &
                           c1intp, c2intp, offset)

        Tair (:) = c1intp *  Tair_data(mlast) + c2intp *  Tair_data(mnext)
        Qa   (:) = c1intp *    Qa_data(mlast) + c2intp *    Qa_data(mnext)
        uatm (:) = c1intp *  uatm_data(mlast) + c2intp *  uatm_data(mnext)
        vatm (:) = c1intp *  vatm_data(mlast) + c2intp *  vatm_data(mnext)
        fsnow(:) = c1intp * fsnow_data(mlast) + c2intp * fsnow_data(mnext)

         ! derived (or not otherwise set)
         potT (:) = c1intp *  potT_data(mlast) + c2intp *  potT_data(mnext)
         wind (:) = c1intp *  wind_data(mlast) + c2intp *  wind_data(mnext)
         strax(:) = c1intp * strax_data(mlast) + c2intp * strax_data(mnext)
         stray(:) = c1intp * stray_data(mlast) + c2intp * stray_data(mnext)
         rhoa (:) = c1intp *  rhoa_data(mlast) + c2intp *  rhoa_data(mnext)
         frain(:) = c1intp * frain_data(mlast) + c2intp * frain_data(mnext)

        sec6hr = secday/c4;                      ! seconds in 6 hours
        maxrec = 1460
        dataloc = 2                              ! data located at end of interval
        recnum = 4*int(yday) - 3 + int(real(sec,kind=dbl_kind)/sec6hr)
        recnum = mod(recnum+maxrec-4*offndy-1,maxrec)+1
        recslot = 2
        mlast = mod(recnum+maxrec-2,maxrec) + 1
        mnext = mod(recnum-1,       maxrec) + 1
        call interp_coeff (recnum, recslot, sec6hr, dataloc, &
                           c1intp, c2intp, offset)

        fsw  (:) = c1intp *   fsw_data(mlast) + c2intp *   fsw_data(mnext)
        flw  (:) = c1intp *   flw_data(mlast) + c2intp *   flw_data(mnext)

         ! derived
         swvdr(:) = c1intp * swvdr_data(mlast) + c2intp * swvdr_data(mnext)
         swvdf(:) = c1intp * swvdf_data(mlast) + c2intp * swvdf_data(mnext)
         swidr(:) = c1intp * swidr_data(mlast) + c2intp * swidr_data(mnext)
         swidf(:) = c1intp * swidf_data(mlast) + c2intp * swidf_data(mnext)

      endif

      if (trim(atm_pert_file)/='') call atm_pert(1)  !add pert to temperature

! possible bug:  is the ocean data also offset to the beginning of the field campaigns?

      if (trim(ocn_data_type) == 'ISPOL') then

         midmonth = 15  ! assume data is given on 15th of every month
         recslot = 1                             ! latter half of month
         if (mday < midmonth) recslot = 2        ! first half of month
         if (recslot == 1) then
            mlast = month
            mnext = mod(month   ,12) + 1
         else ! recslot = 2
            mlast = mod(month+10,12) + 1
            mnext = month
         endif
         call interp_coeff_monthly(recslot, c1intp, c2intp)

         sst_temp(:) = c1intp *  sst_data(mlast) + c2intp *  sst_data(mnext)
         sss     (:) = c1intp *  sss_data(mlast) + c2intp *  sss_data(mnext)
         uocn    (:) = c1intp * uocn_data(mlast) + c2intp * uocn_data(mnext)
         vocn    (:) = c1intp * vocn_data(mlast) + c2intp * vocn_data(mnext)
         qdp     (:) = c1intp *  qdp_data(mlast) + c2intp *  qdp_data(mnext)

      elseif (trim(ocn_data_type) == 'NICE') then

         dataloc = 2                          ! data located at end of interval
         maxrec = 365
         recslot = 2
         recnum = int(yday)
         mlast = mod(recnum+maxrec-2,maxrec) + 1
         mnext = mod(recnum-1,       maxrec) + 1
         call interp_coeff ( recnum, recslot, secday, dataloc, c1intp, c2intp)

         sst_temp(:) = c1intp *  sst_data(mlast) + c2intp *  sst_data(mnext)
         sss     (:) = c1intp *  sss_data(mlast) + c2intp *  sss_data(mnext)
         uocn    (:) = c1intp * uocn_data(mlast) + c2intp * uocn_data(mnext)
         vocn    (:) = c1intp * vocn_data(mlast) + c2intp * vocn_data(mnext)
         qdp     (:) = c1intp *  qdp_data(mlast) + c2intp *  qdp_data(mnext)
         hmix    (:) = c1intp * hmix_data(mlast) + c2intp * hmix_data(mnext)

      elseif (trim(ocn_data_type) == 'CGLO') then
         CALL ocn_CGLO(1)
      else

         ! use default values for all other data fields
         i = mod(timestep-1,ntime)+1 ! repeat forcing cycle
         mlast = i
         mnext = mlast
         c1intp = c1
         c2intp = c0

         sst_temp(:) = c1intp *  sst_data(mlast) + c2intp *  sst_data(mnext)
         sss     (:) = c1intp *  sss_data(mlast) + c2intp *  sss_data(mnext)
         uocn    (:) = c1intp * uocn_data(mlast) + c2intp * uocn_data(mnext)
         vocn    (:) = c1intp * vocn_data(mlast) + c2intp * vocn_data(mnext)
         qdp     (:) = c1intp *  qdp_data(mlast) + c2intp *  qdp_data(mnext)

      endif

      call finish_ocn_forcing(sst_temp)

      ! Lindsay SHEBA open/close dataset is hourly
      if (trim(ocn_data_type) == 'SHEBA') then

        sec1hr = secday/c24                      ! seconds in 1 hour
        maxrec = ntime
        recnum = 24*int(yday) - 23 + int(real(sec,kind=dbl_kind)/sec1hr)
        recslot = 2
        dataloc = 1                          ! data located at middle of interval
        mlast = mod(recnum+maxrec-2,maxrec) + 1
        mnext = mod(recnum-1,       maxrec) + 1
        call interp_coeff ( recnum, recslot, sec1hr, dataloc, c1intp, c2intp)

        opening(:) =   c1intp * open_data(mlast) + c2intp * open_data(mnext)
        closing(:) = -(c1intp * clos_data(mlast) + c2intp * clos_data(mnext))

      endif

      end subroutine get_forcing

!=======================================================================

      subroutine atm_climatological

      real (kind=dbl_kind), dimension(12) :: &
            fsw_clim, & ! field values at temporal data points
            flw_clim, &
           Tair_clim, &
           wind_clim, &
           rhum_clim, &
          fsnow_clim

      real (kind=dbl_kind) :: &
          Tffresh, qqqice, TTTice, rhos

      character(len=*), parameter :: subname='(atm_climatological)'

      ! Ice station meteorology from Lindsay (1998, J. Climate), Table 1, p. 325
      ! zlvl = c2 ! 2-m temperatures and wind speed

      data  fsw_clim /  0.0d0,   1.2d0,  31.5d0, 146.0d0, 263.3d0, 307.9d0, &
                      230.6d0, 134.7d0,  44.2d0,   2.6d0,   0.0d0,   0.0d0  /
      data  flw_clim /164.0d0, 160.5d0, 164.1d0, 188.1d0, 245.2d0, 291.2d0, &
                      303.9d0, 297.0d0, 263.8d0, 210.9d0, 177.0d0, 166.0d0  /
      data Tair_clim /-31.4d0, -32.8d0, -31.6d0, -24.1d0, -11.0d0,  -1.8d0, &
                       -0.1d0,  -1.4d0,  -8.0d0, -19.5d0, -27.6d0, -31.1d0  /
      data rhum_clim / 78.7d0,  78.4d0,  79.6d0,  82.1d0,  86.5d0,  91.7d0, &
                       95.1d0,  94.3d0,  90.7d0,  83.8d0,  80.1d0,  78.7d0  /
      data wind_clim /  4.4d0,   4.0d0,   4.0d0,   3.9d0,   3.9d0,   4.2d0, &
                        4.1d0,   4.2d0,   4.5d0,   4.2d0,   3.9d0,   4.0d0  /
!      data  shf_clim /  9.9d0,   8.4d0,   6.6d0,   0.1d0,  -5.8d0,  -1.6d0, &
!                        2.2d0,   1.2d0,   0.5d0,   2.0d0,   5.6d0,   7.0d0  /
!      data  lhf_clim /  1.3d0,   1.1d0,   1.1d0,   0.0d0,  -5.9d0, -10.3d0, &
!                       -6.5d0,  -6.7d0,  -3.9d0,  -0.1d0,   1.0d0,   1.1d0  /

      ! Semtner (1976, JPO) snowfall spec., p. 383 in m/s snow volume (.4 m/yr)
      data fsnow_clim/ 3.17d-9, 3.17d-9, 3.17d-9, 3.17d-9, 1.90d-8,    0.0d0, &
                           0.0d0, 1.63d-8, 4.89d-8, 4.89d-8, 3.17d-9, 3.17d-9 /

      !-----------------------------------------------------------------
      ! query icepack values
      !-----------------------------------------------------------------

      call icepack_query_parameters(Tffresh_out=Tffresh, qqqice_out=qqqice, &
           TTTice_out=TTTice, rhos_out=rhos)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      !-----------------------------------------------------------------

       fsw_data (1:12) =  fsw_clim (1:12)
       flw_data (1:12) =  flw_clim (1:12)
      rhum_data (1:12) = rhum_clim (1:12)
      wind_data (1:12) = wind_clim (1:12)

      rhoa_data (1:12) = 1.275_dbl_kind ! air density (kg/m^3)
      Tair_data (1:12) = Tair_clim (1:12) + Tffresh
      uatm_data (1:12) = wind_clim (1:12)
      vatm_data (1:12) = c0

      ! Qa = rhum * saturation humidity (1.275 kg/m^3 = air density)
        Qa_data (1:12) = (rhum_clim(1:12)/c100)*qqqice &
                       * exp(-TTTice/Tair_data(1:12))/rhoa_data(1:12)

      fsnow_data(1:12) = rhos*fsnow_clim(1:12) ! convert vol -> mass flux
      frain_data(1:12) = c0

      ! 6 W/m2 warming of mixed layer from deep ocean
        qdp_data(:) = -6.0 ! 2 W/m2 from deep + 4 W/m2 counteracting larger
                              ! SH+LH with bulk transfer than in MU 71

      end subroutine atm_climatological

!=======================================================================

      subroutine atm_CFS

      integer (kind=int_kind) :: &
         nt             ! loop index

      real (kind=dbl_kind) :: &
         dlwsfc,  &     ! downwelling longwave (W/m2)
         dswsfc,  &     ! downwelling shortwave (W/m2)
         windu10, &     ! wind components (m/s)
         windv10, &     !
         temp2m,  &     ! 2m air temperature (K)
         spechum ,&     ! specific humidity (kg/kg)
         precip         ! precipitation (kg/m2/s)

      character (char_len_long) string1
      character (char_len_long) filename
      character(len=*), parameter :: subname='(atm_CFS)'

!      atm_data_file = 'cfsv2_2015_220_70_01hr.txt'
      filename = trim(data_dir)//'/CFS/'//trim(atm_data_file)

      write (nu_diag,*) 'Reading ',filename

      open (nu_forcing, file=filename, form='formatted')
      read (nu_forcing, *) string1 ! headers
      read (nu_forcing, *) string1 ! units

      do nt = 1, ntime
         read (nu_forcing, '(6(f10.5,1x),2(f10.8,1x))') &
         dswsfc, dlwsfc, windu10, windv10, temp2m, spechum, precip

           flw_data(nt) = dlwsfc
           fsw_data(nt) = dswsfc
          uatm_data(nt) = windu10
          vatm_data(nt) = windv10
          Tair_data(nt) = temp2m
          potT_data(nt) = temp2m
            Qa_data(nt) = spechum
         fsnow_data(nt) = precip
      enddo

      close (nu_forcing)

      end subroutine atm_CFS

!=======================================================================

       subroutine setup_nc_forcing(nc_dataset)
       USE netcdf

       TYPE(nc_data), INTENT(INOUT) :: nc_dataset
 
       integer (kind=int_kind) :: &
           nt, &             ! time coordinate dimension
           nx, &             ! longitude dimension
           ny, &             ! latitude dimension
           id_x, &           ! longitude point index
           id_y, &           ! latitude point index
           ierr, &           ! error code netcdf calls
           ivarid, &         ! netcdf variable id
           idimids(1), &        ! netcdf variable dimensions' ids
           indates, &        ! time variable len
           itype,   &        ! time variable type
           id                ! netcdf file id
       integer (kind=int8_kind) :: &
           idate0(7)               ! reference date in dataset 
       integer (kind=int8_kind), allocatable :: &
           idates(:)      ! integer time variable
 
       real (kind=dbl_kind) :: &
           dlwsfc,  &     ! downwelling longwave (W/m2)
           dswsfc,  &     ! downwelling shortwave (W/m2)
           windu10, &     ! wind components (m/s)
           windv10, &     !
           temp2m,  &     ! 2m air temperature (K)
           spechum ,&     ! specific humidity (kg/kg)
           precip,  &     ! precipitation (kg/m2/s)
           lon_p,   &     ! point longitude 
           lat_p,   &     ! point latitude 
           zx_min,   &     ! min netcdf x (longitude)  
           zy_min,   &     ! min netcdf y (latitude) 
           zx_max,   &     ! max netcdf x (longitude) 
           zy_max,   &     ! max netcdf y (latitude) 
           rdummy
 
       real (kind=dbl_kind), allocatable :: &
           zx_reg(:), &   ! x coordinate netcdf 
           zy_reg(:), &   ! y coordinate netcdf 
           zdates(:)

       logical :: &
           gout_x, &      ! point outside domain flag
           gout_y, &      ! point outside domain flag
           ldummy
 
       character (char_len_long) :: filename
       character(len=10) :: ylonname, ylatname
       character(len=*), parameter :: subname='(setup_nc_forcing)'
       character(len=*), parameter :: ytimevar='time'
       character(len=80) :: ytunits, yreal_units

       filename = nc_dataset%y_filename 
       ylonname = nc_dataset%y_lonname 
       ylatname = nc_dataset%y_latname 

!      open netcdf dataset
       call nc_open_read(id,filename)  
       nc_dataset%file_id = id

!      read lon lat dimesions
       call nc_read_dim(id,ylonname,nx)
       call nc_read_dim(id,ylatname,ny)

       allocate(zx_reg(nx))
       allocate(zy_reg(ny))

!      read netcdf lon lat
       call nc_read_var(id,ylonname,zx_reg)
       call nc_read_var(id,ylatname,zy_reg)
       zx_min = MINVAL(zx_reg)
       zy_min = MINVAL(zy_reg)
       zx_max = MAXVAL(zx_reg)
       zy_max = MAXVAL(zy_reg)

!      test point
       lon_p = nc_dataset%zpoint_lon
       lat_p = nc_dataset%zpoint_lat

!      verify that netcdf domain includes sim coordinates 
       gout_x = (lon_p .LT. zx_min) .OR. (lon_p .GT. zx_max)          
       gout_y = (lat_p .LT. zy_min) .OR. (lat_p .GT. zy_max)          

       IF(gout_x) call handle_err(0,'Error: point outside atm forcing regular domain (LON)', subname)
       IF(gout_y) call handle_err(0,'Error: point outside atm forcing regular domain (LAT)', subname)

!      determine lon lat indices 
       id_x = FINDLOC(zx_reg, lon_p, 1)    
       id_y = FINDLOC(zy_reg, lat_p, 1)    
       nc_dataset%id_x = id_x
       nc_dataset%id_y = id_y

!      Time coordinate management
       ierr = nf90_inq_varid(id,ytimevar,ivarid)
       IF (ierr /= nf90_noerr) CALL handle_err(ierr,'Error: nf90_inq_varid generated an error', subname)

       ierr = nf90_inquire_variable(id, ivarid, xtype=itype, dimids=idimids)
       IF (ierr /= nf90_noerr) CALL handle_err(ierr, 'Error: nf90_inquire_variable generated an error', subname)

       ierr = nf90_get_att(id,ivarid,'units',ytunits)
       IF (ierr /= nf90_noerr) CALL handle_err(ierr, 'Error: nf90_get_att generated an error',subname)

       IF (ytunits(1:7) .EQ. "seconds") yreal_units = "seconds"
       IF (ytunits(1:7) .EQ. "minutes") yreal_units = "minutes"
       IF (ytunits(1:5) .EQ. "hours"  ) yreal_units = "hours"
       IF (ytunits(1:4) .EQ. "days"   ) yreal_units = "days"

       !! get reference date of dataset
       idate0 =  string2date_nc(ytunits)

!       CALL nc_read_dim(id,'time',indates)
       ierr = nf90_inquire_dimension(id,idimids(1),len=indates)
       IF (ierr /= nf90_noerr) CALL handle_err(ierr, 'nf90_inquire_dimension generated an error', subname)

       IF ( itype .EQ. 4 .OR. itype .EQ. 10) THEN                !! type 10 for INT64
          ALLOCATE( idates(indates) )
          ALLOCATE( zdates(indates) )
          ierr = nf90_get_var(id,ivarid,idates)
          zdates = DBLE(idates)
       ELSE IF ( itype .EQ. 5 .OR. itype .EQ. 6 ) THEN           !! float time
          ALLOCATE( zdates(indates))
          ierr = nf90_get_var(id,ivarid,zdates)
       END IF

       IF (ierr /= nf90_noerr) CALL handle_err(ierr, 'Error: nf90_get_var generated an error: '//nf90_strerror(ierr), subname)

       nc_dataset%idate_0 = idate0
       nc_dataset%y_tunits = yreal_units
       nc_dataset%ztimes = zdates


       DEALLOCATE(zdates)
       IF(ALLOCATED(idates)) DEALLOCATE(idates)


       end subroutine setup_nc_forcing

!=======================================================================


       subroutine setup_nc1d_forcing(nc_dataset)
       USE netcdf

       TYPE(nc_data1d), INTENT(INOUT) :: nc_dataset
 
       integer (kind=int_kind) :: &
           nt, &             ! time coordinate dimension
           ierr, &           ! error code netcdf calls
           ivarid, &         ! netcdf variable id
           idimids(1), &        ! netcdf variable dimensions' ids
           indates, &        ! time variable len
           itype,   &        ! time variable type
           id                ! netcdf file id
       integer (kind=int8_kind) :: &
           idate0(7)               ! reference date in dataset 
       integer (kind=int8_kind), allocatable :: &
           idates(:)      ! integer time variable
 
       real (kind=dbl_kind) :: &
           rdummy
 
       real (kind=dbl_kind), allocatable :: &
           zdates(:)

       character (char_len_long) :: filename
       character(len=*), parameter :: subname='(setup_nc1d_forcing)'
       character(len=*), parameter :: ytimevar='time'
       character(len=80) :: ytunits, yreal_units

       filename = nc_dataset%y_filename 

!      open netcdf dataset
       call nc_open_read(id,filename)  
       nc_dataset%file_id = id

!      Time coordinate management
       ierr = nf90_inq_varid(id,ytimevar,ivarid)
       IF (ierr /= nf90_noerr) CALL handle_err(ierr,'Error: nf90_inq_varid generated an error', subname)

       ierr = nf90_inquire_variable(id, ivarid, xtype=itype, dimids=idimids)
       IF (ierr /= nf90_noerr) CALL handle_err(ierr, 'Error: nf90_inquire_variable generated an error', subname)

       ierr = nf90_get_att(id,ivarid,'units',ytunits)
       IF (ierr /= nf90_noerr) CALL handle_err(ierr, 'Error: nf90_get_att generated an error',subname)

       IF (ytunits(1:7) .EQ. "seconds") yreal_units = "seconds"
       IF (ytunits(1:7) .EQ. "minutes") yreal_units = "minutes"
       IF (ytunits(1:5) .EQ. "hours"  ) yreal_units = "hours"
       IF (ytunits(1:4) .EQ. "days"   ) yreal_units = "days"

       !! get reference date of dataset
       idate0 =  string2date_nc(ytunits)

!       CALL nc_read_dim(id,'time',indates)
       ierr = nf90_inquire_dimension(id,idimids(1),len=indates)
       IF (ierr /= nf90_noerr) CALL handle_err(ierr, 'nf90_inquire_dimension generated an error', subname)

       IF ( itype .EQ. 4 .OR. itype .EQ. 10) THEN                !! type 10 for INT64
          ALLOCATE( idates(indates) )
          ALLOCATE( zdates(indates) )
          ierr = nf90_get_var(id,ivarid,idates)
          zdates = DBLE(idates)
       ELSE IF ( itype .EQ. 5 .OR. itype .EQ. 6 ) THEN           !! float time
          ALLOCATE( zdates(indates))
          ierr = nf90_get_var(id,ivarid,zdates)
       END IF

       IF (ierr /= nf90_noerr) CALL handle_err(ierr, 'Error: nf90_get_var generated an error: '//nf90_strerror(ierr), subname)

       nc_dataset%idate_0 = idate0
       nc_dataset%y_tunits = yreal_units
       nc_dataset%ztimes = zdates


       DEALLOCATE(zdates)
       IF(ALLOCATED(idates)) DEALLOCATE(idates)


       end subroutine setup_nc1d_forcing

!=======================================================================



!=======================================================================
!   setup ECMWF nc dataset and get data
       SUBROUTINE atm_ECMWF(mode)

       IMPLICIT NONE

       !arguments
       integer (kind=int_kind), intent(in) :: mode

       ! local
       integer(kind=8)  :: idate(7) !model date 
       integer          :: itimesid(2)
       DOUBLE PRECISION :: zweights(2)
       DOUBLE PRECISION :: zu10     !10m U wind component [m/s]
       DOUBLE PRECISION :: zv10     !10m V wind component [m/s]
       DOUBLE PRECISION :: zt2m     !2m temperature [K]
       DOUBLE PRECISION :: zd2m     !2m dew point temperature [K]
       DOUBLE PRECISION :: zssrd    !short-wave radiation downward [J/m**2]
       DOUBLE PRECISION :: zstrd    !thermal (long-wave) radiation downward [J/m**2]
       DOUBLE PRECISION :: ztp      !total precipitation [m]
       DOUBLE PRECISION :: zmslp    !mean sea level pressure [Pa]
       LOGICAL          :: lice(nx)
       INTEGER          :: i        !do loop index 
       real (kind=dbl_kind) :: &
           Tffresh, &
           rhofresh

          SELECT CASE (mode)
              CASE(0) ! initialization 
                  ecmwf_data%y_filename = trim(data_dir)//'/'//trim(atm_data_file)
                  ecmwf_data%y_lonname = "lon"
                  ecmwf_data%y_latname = "lat"
                  ecmwf_data%zpoint_lon = atm_longitude   !read coordinates from namelist 
                  ecmwf_data%zpoint_lat = atm_latitude
                  CALL setup_nc_forcing(ecmwf_data)

!                  IF (atm_latitude.GE.80) qdp_data(:) = -3.0 ! 2 W/m2 from deep + 4 W/m2 counteracting larger !used by def slab ocean
!                  IF (atm_latitude.GE.85) qdp_data(:) = -6.0 ! 2 W/m2 from deep + 4 W/m2 counteracting larger !used by def slab ocean

              CASE(1) ! get data and interpolate in time

                  call icepack_query_parameters(Tffresh_out=Tffresh, rhofresh_out=rhofresh)
                  idate = get_model_date(ycycle)

                  CALL nc_findtimes(idate, ecmwf_data%idate_0, ecmwf_data%y_tunits, ecmwf_data%ztimes, itimesid, zweights) 
                  CALL nc_getscalar(ecmwf_data, 'u10', itimesid, zweights, zu10) 
                  CALL nc_getscalar(ecmwf_data, 'v10', itimesid, zweights, zv10) 
                  CALL nc_getscalar(ecmwf_data, 't2m', itimesid, zweights, zt2m) 
                  CALL nc_getscalar(ecmwf_data, 'd2m', itimesid, zweights, zd2m) 
                  CALL nc_getscalar(ecmwf_data, 'ssrd', itimesid, zweights, zssrd) 
                  CALL nc_getscalar(ecmwf_data, 'strd', itimesid, zweights, zstrd) 
                  CALL nc_getscalar(ecmwf_data, 'msl', itimesid, zweights, zmslp) 
                  CALL nc_getscalar(ecmwf_data, 'tp', itimesid, zweights, ztp) 

                  
                  uatm(:) = zu10
                  vatm(:) = zv10
                  Tair(:) = zt2m

                  !convert radiations from J/m**2 to W/m**2
                  fsw(:) = zssrd/3600.D0 
                  flw(:) = zstrd/3600.D0 

                  ! divide shortwave into spectral bands
                  swvdr(:) = fsw(:)*frcvdr        ! visible direct
                  swvdf(:) = fsw(:)*frcvdf        ! visible diffuse
                  swidr(:) = fsw(:)*frcidr        ! near IR direct
                  swidf(:) = fsw(:)*frcidf        ! near IR diffuse

                  !convert 2m dew point temperature to specific humidity
                  lice = aice.GT.c0   !what is tmask?
                  DO i = 1, nx
                     Qa(i) = q_sat(zd2m, zmslp, lice(i) )  
                  ENDDO

                  ! Convert total precipitation from m to kg/m**2/s
                  fsnow(:) = ztp*rhofresh/3600.D0

                  ! split precipitation into rain and snow
                  frain(:) = c0
                  if (zt2m >= Tffresh) then
                      frain(:) = fsnow(:)
                      fsnow(:) = c0
                  endif
                  ! 
                  zlvl(:) = c10

                  ! air density
                  DO i = 1, nx
                     rhoa(i) = rho_air(zt2m, Qa(i), zmslp) 
                  ENDDO  
                  !rhoa(:) = 1.3_dbl_kind

                  ! potential temperature
                  potT(:) = Tair(:)*(1.D5/zmslp)**0.286    !
                  !potT(:) = Tair(:)    !

                  ! wind speed
                  wind(:) = SQRT(uatm(:)**2 + vatm(:)**2)

                  ! wind stress
                  strax(:) = c0
                  stray(:) = c0

          END SELECT

      end subroutine atm_ECMWF

!=======================================================================
!   setup atm_pert nc dataset and get data
       SUBROUTINE atm_pert(mode)

       IMPLICIT NONE

       !arguments
       integer (kind=int_kind), intent(in) :: mode

       ! local
       integer(kind=8)  :: idate(7) !model date 
       integer          :: itimesid(2)
       DOUBLE PRECISION :: zweights(2)
       DOUBLE PRECISION :: z_pt2m     !2m temperature perturbation[K]
       DOUBLE PRECISION :: z_pU10     !10m wind speed perturbation[m/s]
       INTEGER          :: i        !do loop index 
          
          IF (trim(atm_pert_file).EQ."") RETURN

          SELECT CASE (mode)
              CASE(0) ! initialization 
                  pert_atm_data%y_filename = trim(data_dir)//'/'//trim(atm_pert_file)
                  CALL setup_nc1d_forcing(pert_atm_data)

              CASE(1) ! get data and interpolate in time

                  idate = get_model_date(ycycle)

                  CALL nc_findtimes(idate, pert_atm_data%idate_0, pert_atm_data%y_tunits, pert_atm_data%ztimes, itimesid, zweights) 
                  CALL nc_getscalar1d(pert_atm_data, 't2m_pert', itimesid, zweights, z_pt2m) 
                  CALL nc_getscalar1d(pert_atm_data, 'U10_pert', itimesid, zweights, z_pU10) 

                  Tair(:) = Tair(:) + z_pt2m  !add perturbation to temperature

                  uatm(:) = z_pU10 * uatm(:)/wind(:)  !add pert to wind speeds
                  vatm(:) = z_pU10 * vatm(:)/wind(:)

                  ! wind speed
                  wind(:) = SQRT(uatm(:)**2 + vatm(:)**2)

                  
          END SELECT

      end subroutine atm_pert


      subroutine prepare_forcing (Tair,     fsw,      &
                                  cldf,     &
                                  frain,    fsnow,    &
                                  Qa,       rhoa,     &
                                  uatm,     vatm,     &
                                  strax,    stray,    &
                                  zlvl,     wind,     &
                                  swvdr,    swvdf,    &
                                  swidr,    swidf,    &
                                  potT)

      ! this routine acts on the data fields prior to interpolation

      real (kind=dbl_kind), dimension(ntime), &
         intent(inout) :: &
         Tair    , & ! air temperature  (K)
         fsw     , & ! incoming shortwave radiation (W/m^2)
         cldf    , & ! cloud fraction
         frain   , & ! rainfall rate (kg/m^2 s)
         fsnow   , & ! snowfall rate (kg/m^2 s)
         Qa      , & ! specific humidity (kg/kg)
         rhoa    , & ! air density (kg/m^3)
         uatm    , & ! wind velocity components (m/s)
         vatm    , &
         strax   , & ! wind stress components (N/m^2)
         stray   , &
         zlvl    , & ! atm level height (m)
         wind    , & ! wind speed (m/s)
!        flw     , & ! incoming longwave radiation (W/m^2)
         swvdr   , & ! sw down, visible, direct  (W/m^2)
         swvdf   , & ! sw down, visible, diffuse (W/m^2)
         swidr   , & ! sw down, near IR, direct  (W/m^2)
         swidf   , & ! sw down, near IR, diffuse (W/m^2)
         potT        ! air potential temperature  (K)

      ! local variables

      integer (kind=int_kind) :: &
         nt

      logical (kind=log_kind) :: &
         calc_strair

      real (kind=dbl_kind), parameter :: &
         lapse_rate = 0.0065_dbl_kind      ! (K/m) lapse rate over sea level

      real (kind=dbl_kind) :: &
         precip_factor, zlvl0, &
         Tffresh

      character(len=*), parameter :: subname='(prepare_forcing)'

      zlvl0 = c10 ! default

      !-----------------------------------------------------------------
      ! query icepack values
      !-----------------------------------------------------------------

      call icepack_query_parameters(calc_strair_out=calc_strair)
      call icepack_query_parameters(Tffresh_out=Tffresh)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      !-----------------------------------------------------------------
      ! convert precipitation units to kg/m^2 s
      !-----------------------------------------------------------------
      if (trim(precip_units) == 'mm_per_month') then
         precip_factor = 12._dbl_kind/(secday*dayyr)
      elseif (trim(precip_units) == 'mm_per_day') then
         precip_factor = c1/secday
      elseif (trim(precip_units) == 'mm_per_sec' .or. &
              trim(precip_units) == 'mks') then
         precip_factor = c1    ! mm/sec = kg/m^2 s
      endif

      do nt = 1, ntime

      !-----------------------------------------------------------------
      ! make sure interpolated values are physically realistic
      !-----------------------------------------------------------------
         cldf (nt) = max(min(cldf(nt),c1),c0)
         fsw  (nt) = max(fsw(nt),c0)
         fsnow(nt) = max(fsnow(nt),c0)
         rhoa (nt) = max(rhoa(nt),c0)
         Qa   (nt) = max(Qa(nt),c0)

      !-----------------------------------------------------------------
      ! calculations specific to datasets
      !-----------------------------------------------------------------

         if (trim(atm_data_type) == 'CFS') then
            ! precip is in kg/m^2/s
            zlvl0 = c10
            ! downward longwave as in Parkinson and Washington (1979)
!            call longwave_parkinson_washington(Tair(nt), cldf(nt), flw(nt))

         elseif (trim(atm_data_type) == 'clim') then
            ! precip is in kg/m^2/s
            zlvl0 = c2

         elseif (trim(atm_data_type) == 'ISPOL' .or. &
                 trim(atm_data_type) == 'NICE') then
            Tair(nt) = Tair(nt) -  lapse_rate*8.0_dbl_kind

         endif                     ! atm_data_type

! this longwave depends on the current ice aice and sst and so can not be
! computed ahead of time
!            ! longwave based on Rosati and Miyakoda, JPO 18, p. 1607 (1988)
!            call longwave_rosati_miyakoda(cldf(i,j), Tsfc(i,j), &
!                                          aice(i,j), sst(i,j),  &
!                                          Qa(i,j),   Tair(i,j), &
!                                          hm(i,j),   flw(i,j))

      !-----------------------------------------------------------------
      ! Compute other fields needed by model
      !-----------------------------------------------------------------

         zlvl(nt) = zlvl0
         potT(nt) = Tair(nt)

         ! divide shortwave into spectral bands
         swvdr(nt) = fsw(nt)*frcvdr        ! visible direct
         swvdf(nt) = fsw(nt)*frcvdf        ! visible diffuse
         swidr(nt) = fsw(nt)*frcidr        ! near IR direct
         swidf(nt) = fsw(nt)*frcidf        ! near IR diffuse

         ! precipitation
         fsnow(nt) = fsnow(nt) * precip_factor

         ! determine whether precip is rain or snow
         ! HadGEM forcing provides separate snowfall and rainfall rather
         ! than total precipitation
!         if (trim(atm_data_type) /= 'hadgem') then
            frain(nt) = c0
            if (Tair(nt) >= Tffresh) then
                frain(nt) = fsnow(nt)
                fsnow(nt) = c0
            endif
!         endif

         if (calc_strair) then
            wind (nt) = sqrt(uatm(nt)**2 + vatm(nt)**2)
            strax(nt) = c0
            stray(nt) = c0
         ! else  ! strax, stray, wind are read from files
         endif                   ! calc_strair

      enddo ! ntime

      end subroutine prepare_forcing

!=======================================================================

      subroutine interp_coeff_monthly (recslot, c1intp, c2intp)

! Compute coefficients for interpolating monthly data to current time step.

      integer (kind=int_kind), intent(in) :: &
          recslot         ! slot (1 or 2) for current record

      real (kind=dbl_kind), intent(inout) :: &
         c1intp, c2intp   ! interpolation coefficients

      ! local variables

      real (kind=dbl_kind) :: &
          tt           , & ! seconds elapsed in current year
          t1, t2           ! seconds elapsed at month midpoint

      real (kind=dbl_kind) :: &
          daymid(0:13)     ! month mid-points

      character(len=*), parameter :: subname='(interp_coeff_monthly)'

      daymid(1:13) = 14._dbl_kind   ! time frame ends 0 sec into day 15
      daymid(0)    = 14._dbl_kind - daymo(12)  ! Dec 15, 0 sec

      ! make time cyclic
      tt = mod(time/secday,dayyr)

      ! Find neighboring times

      if (recslot==2) then      ! first half of month
        t2 = daycal(month) + daymid(month)   ! midpoint, current month
        if (month == 1) then
          t1 = daymid(0)                 ! Dec 15 (0 sec)
        else
          t1 = daycal(month-1) + daymid(month-1) ! midpoint, previous month
        endif
      else                      ! second half of month
        t1 = daycal(month) + daymid(month)    ! midpoint, current month
        t2 = daycal(month+1) + daymid(month+1)! day 15 of next month (0 sec)
      endif

      ! Compute coefficients
      c1intp = (t2 - tt) / (t2 - t1)
      c2intp =  c1 - c1intp

      end subroutine interp_coeff_monthly

!=======================================================================

      subroutine interp_coeff (recnum, recslot, secint, dataloc, &
                               c1intp, c2intp, offset)

! Compute coefficients for interpolating data to current time step.
! Works for any data interval that divides evenly into a
!  year (daily, 6-hourly, etc.)
! Use interp_coef_monthly for monthly data.

      integer (kind=int_kind), intent(in) :: &
          recnum      , & ! record number for current data value
          recslot     , & ! spline slot for current record
          dataloc         ! = 1 for data located in middle of time interval
                          ! = 2 for date located at end of time interval

      real (kind=dbl_kind), intent(in) :: &
          secint          ! seconds in data interval

      real (kind=dbl_kind), intent(inout) :: &
         c1intp, c2intp   ! interpolation coefficients

      real (kind=dbl_kind), intent(in), optional :: &
          offset          ! amount of time data is offset (s)

      ! local variables

      real (kind=dbl_kind) :: &
          secyr            ! seconds in a year

      real (kind=dbl_kind) :: &
          tt           , & ! seconds elapsed in current year
          t1, t2       , & ! seconds elapsed at data points
          rcnum            ! recnum => dbl_kind

      character(len=*), parameter :: subname='(interp_coeff)'

      secyr = dayyr * secday         ! seconds in a year
      if (present(offset)) then
         tt = mod(time-offset,secyr)
      else
         tt = mod(time,secyr)
      endif

      ! Find neighboring times
      rcnum = real(recnum,kind=dbl_kind)
      if (recslot==2) then           ! current record goes in slot 2
         if (dataloc==1) then        ! data located at middle of interval
            t2 = (rcnum-p5)*secint
         else                        !  data located at end of interval
            t2 = rcnum*secint
         endif
         t1 = t2 - secint            !  - 1 interval
      else                           ! recslot = 1
         if (dataloc==1) then        ! data located at middle of interval
            t1 = (rcnum-p5)*secint
         else
            t1 = rcnum*secint        ! data located at end of interval
         endif
         t2 = t1 + secint            !  + 1 interval
      endif

      ! Compute coefficients
      c1intp =  abs((t2 - tt) / (t2 - t1))
      c2intp =  c1 - c1intp

      end subroutine interp_coeff

!=======================================================================

      subroutine atm_ISPOL
      ! there is data for 366 days, but we only use 365

      integer (kind=int_kind) :: &
         i        ! index

      real (kind=dbl_kind), dimension(366) :: &
         tair , & ! air temperature
         qa   , & ! specific humidity
         uatm , & ! wind velocity components
         vatm , &
         fsnow    ! snowfall rate
         ! aday   ! not used

      real (kind=dbl_kind), dimension(1464) :: &
         fsw  , & ! shortwave
         flw      ! longwave
         ! atime  ! not used

      character (char_len_long) filename

      character(len=*), parameter :: subname='(atm_ISPOL)'

!      atm_data_file = 'ISPOL_atm_forcing.txt'
      filename = trim(data_dir)//'/ISPOL_2004/'//trim(atm_data_file)

      write (nu_diag,*) 'Reading ',filename

      open (nu_forcing, file=filename, form='formatted')

      read(nu_forcing,*) tair
      read(nu_forcing,*) qa
      read(nu_forcing,*) fsw
      read(nu_forcing,*) flw
      read(nu_forcing,*) uatm
      read(nu_forcing,*) vatm
      read(nu_forcing,*) fsnow
!      read(nu_forcing,*) aday
!      read(nu_forcing,*) atime

      do i = 1, 366 ! daily
         Tair_data (i) = tair (i)
         Qa_data   (i) = qa   (i)
         uatm_data (i) = uatm (i)
         vatm_data (i) = vatm (i)
         fsnow_data(i) = fsnow(i)
      end do
      do i = 1, 1464 ! 6hr, 1464/4=366 days
         fsw_data  (i) = fsw  (i)
         flw_data  (i) = flw  (i)
      end do

      close(nu_forcing)

      end subroutine atm_ISPOL

!=======================================================================

      subroutine atm_NICE
      ! there is data for 366 days, but we only use 365

      integer (kind=int_kind) :: &
         i        ! index

      real (kind=dbl_kind), dimension(366) :: &
         tair , & ! air temperature
         qa   , & ! specific humidity
         uatm , & ! wind velocity components
         vatm , &
         fsnow    ! snowfall rate
         ! aday   ! not used

      real (kind=dbl_kind), dimension(1464) :: &
         fsw  , & ! shortwave
         flw      ! longwave
         ! atime  ! not used

      character (char_len_long) filename

      character(len=*), parameter :: subname='(atm_NICE)'

!      atm_data_file = 'NICE_atm_forcing.txt'
      filename = trim(data_dir)//'/NICE_2015/'//trim(atm_data_file)

      write (nu_diag,*) 'Reading ',filename

      open (nu_forcing, file=filename, form='formatted')

      read(nu_forcing,*) tair
      read(nu_forcing,*) qa
      read(nu_forcing,*) fsw
      read(nu_forcing,*) flw
      read(nu_forcing,*) uatm
      read(nu_forcing,*) vatm
      read(nu_forcing,*) fsnow
!      read(nu_forcing,*) aday
!      read(nu_forcing,*) atime

      do i = 1, 366
         Tair_data (i) = tair (i)
         Qa_data   (i) = qa   (i)
         uatm_data (i) = uatm (i)
         vatm_data (i) = vatm (i)
         fsnow_data(i) = fsnow(i)
      end do
      do i = 1, 1464
         fsw_data  (i) = fsw  (i)
         flw_data  (i) = flw  (i)
      end do

      close(nu_forcing)

      end subroutine atm_NICE

!=======================================================================

      subroutine ocn_NICE

      integer (kind=int_kind) :: &
         i

      real (kind=dbl_kind), dimension(365) :: &
         t   , &  ! sea surface temperature
         s   , &  ! sea surface salinity
         hblt, &  ! mixed layer depth
         u   , &  ! ocean current, x
         v   , &  ! ocean current, y
         dhdx, &  ! sea surface slope
         dhdy, &  ! sea surface slope
         qdp      ! deep ocean heat flux

      character (char_len_long) filename

      character(len=*), parameter :: subname='(ocn_NICE)'

!      ocn_data_file = 'oceanmixed_daily_3.txt'
      filename = trim(data_dir)//'/NICE_2015/'//trim(ocn_data_file)

      write (nu_diag,*) 'Reading ',filename

      open (nu_forcing, file=filename, form='formatted')

      read(nu_forcing,*) t
      read(nu_forcing,*) s
      read(nu_forcing,*) hblt
      read(nu_forcing,*) u
      read(nu_forcing,*) v
      read(nu_forcing,*) dhdx  ! not used for Icepack
      read(nu_forcing,*) dhdy  ! not used for Icepack
      read(nu_forcing,*) qdp

      close(nu_forcing)

      do i = 1, 365 ! daily
         sst_data (i) = t   (i)
         sss_data (i) = s   (i)
         hmix_data(i) = hblt(i)
         uocn_data(i) = u   (i)
         vocn_data(i) = v   (i)
         qdp_data (i) = qdp (i)
      end do

      end subroutine ocn_NICE

!=======================================================================
!   setup ECMWF nc dataset and get data
       SUBROUTINE ocn_CGLO(mode)

       IMPLICIT NONE

       !arguments
       integer (kind=int_kind), intent(in) :: mode

       ! local
       integer(kind=8)  :: idate(7) !model date 
       integer          :: itimesid(2)
       DOUBLE PRECISION :: zweights(2)
       DOUBLE PRECISION :: zuocn     !ocean surface U current [m/s]
       DOUBLE PRECISION :: zvocn     !ocean surface V current [m/s]
       DOUBLE PRECISION :: zsst      !sea surface potential temperature [C]
       DOUBLE PRECISION :: zsss      !sea surface salinity [psu]
       INTEGER          :: i        !do loop index 
       real (kind=dbl_kind) :: &
           Tffresh, &
           rhofresh

          SELECT CASE (mode)
              CASE(0) ! initialization 
                    
                  !if(oceanmixed_ice) then
                  !    oceanmixed_ice = .false.
                  !    WRITE(nu_diag, *) 'Turning off mixed layer model when forcing with CGLORS'
                  !endif

                  cglors_data%y_filename = trim(data_dir)//'/'//trim(ocn_data_file)
                  cglors_data%y_lonname = "longitude"
                  cglors_data%y_latname = "latitude"
                  cglors_data%zpoint_lon = 150.0   !read coordinates from namelist 
                  cglors_data%zpoint_lat = 80.0
                  CALL setup_nc_forcing(cglors_data)

              CASE(1) ! get data and interpolate in time

                  idate = get_model_date()

                  CALL nc_findtimes(idate, cglors_data%idate_0, cglors_data%y_tunits, cglors_data%ztimes, itimesid, zweights) 
                  CALL nc_getscalar(cglors_data, 'uo_cglo', itimesid, zweights, zuocn) 
                  CALL nc_getscalar(cglors_data, 'vo_cglo', itimesid, zweights, zvocn) 
                  CALL nc_getscalar(cglors_data, 'thetao_cglo', itimesid, zweights, zsst) 
                  CALL nc_getscalar(cglors_data, 'so_cglo', itimesid, zweights, zsss) 

                  sst_temp(:) = zsst
                  !sst     (:) = zsst
                  sss     (:) = zsss
                  uocn    (:) = zuocn
                  vocn    (:) = zvocn

          END SELECT

      end subroutine ocn_CGLO


!=======================================================================

      subroutine ocn_ISPOL

      integer (kind=int_kind) :: &
         i

      real (kind=dbl_kind), dimension(12) :: &
         t   , &  ! sea surface temperature
         s   , &  ! sea surface salinity
         hblt, &  ! mixed layer depth
         u   , &  ! ocean current, x
         v   , &  ! ocean current, y
         dhdx, &  ! sea surface slope
         dhdy, &  ! sea surface slope
         qdp      ! deep ocean heat flux

      character (char_len_long) filename

      character(len=*), parameter :: subname='(ocn_ISPOL)'

!      ocn_data_file = 'pop_frc.gx1v3.051202_but_hblt_from_010815_ispol.txt'
      filename = trim(data_dir)//'/ISPOL_2004/'//trim(ocn_data_file)

      write (nu_diag,*) 'Reading ',filename

      open (nu_forcing, file=filename, form='formatted')

      read(nu_forcing,*) t
      read(nu_forcing,*) s
      read(nu_forcing,*) hblt
      read(nu_forcing,*) u
      read(nu_forcing,*) v
      read(nu_forcing,*) dhdx  ! not used for Icepack
      read(nu_forcing,*) dhdy  ! not used for Icepack
      read(nu_forcing,*) qdp

      close(nu_forcing)

      do i = 1, 12 ! monthly
         sst_data (i) = t   (i)
         sss_data (i) = s   (i)
         hmix_data(i) = hblt(i)
         uocn_data(i) = u   (i)
         vocn_data(i) = v   (i)
         qdp_data (i) = qdp (i)
      end do

     end subroutine ocn_ISPOL

!=======================================================================

      subroutine finish_ocn_forcing(sst_temp)

! Compute ocean freezing temperature Tf based on tfrz_option
! 'minus1p8'         Tf = -1.8 C
! 'constant'         Tf = Tocnfrz
! 'linear_salt'      Tf = -depressT * sss
! 'mushy'            Tf conforms with mushy layer thermo (ktherm=2)

      real (kind=dbl_kind), dimension(nx), intent(in)  :: &
          sst_temp

      ! local variables

      integer (kind=int_kind) :: &
         i           ! horizontal indices

      real (kind=dbl_kind)   :: &
          cprho

      real (kind=dbl_kind), parameter :: &
         frzmlt_max = 1000.0_dbl_kind   ! max magnitude of frzmlt (W/m^2)

      character(len=*), parameter :: subname='(finish_ocn_forcing)'

      call icepack_query_parameters(cprho_out=cprho)

      do i = 1, nx
         sss (i) = max (sss(i), c0)
         hmix(i) = max(hmix(i), c0)
         Tf  (i) = icepack_sea_freezing_temperature(sss(i))
         !if ocn=CGLORS sst = sst_temp ==> sst is read from data in any case
         if (restore_ocn) sst(i) = sst(i) + (sst_temp(i)-sst(i))*dt/trest
         if (.not.oceanmixed_ice) then 
             frzmlt(i) = (Tf(i)-sst(i)) * cprho * hmix(i) / dt 
             frzmlt(i) = min(max(frzmlt(i),-frzmlt_max), frzmlt_max)
         endif
         
      enddo
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      end subroutine finish_ocn_forcing

!=======================================================================

     subroutine ice_open_clos


      integer (kind=int_kind) :: i

      real (kind=dbl_kind) :: xtime

      character (char_len_long) filename

!      ice_data_file = 'open_clos_lindsay.dat'
      filename = trim(data_dir)//'/SHEBA/'//trim(ice_data_file)

      write (nu_diag,*) 'Reading ',filename

      open (nu_open_clos, file=filename, form='formatted')

      ! hourly data
      do i=1,ntime
         read(nu_open_clos,*) xtime, open_data(i), clos_data(i)
      enddo

      close (nu_open_clos)

     end subroutine ice_open_clos

!=======================================================================

      subroutine get_wave_spec

      use icedrv_arrays_column, only: wave_spectrum, wave_sig_ht, &
                                   dwavefreq, wavefreq
      use icedrv_domain_size, only: nfreq

      ! local variables
      integer (kind=int_kind) :: &
         k

      real(kind=dbl_kind), dimension(nfreq) :: &
         wave_spectrum_profile  ! wave spectrum

       wave_spectrum(:,:) = c0

      ! wave spectrum and frequencies
      ! get hardwired frequency bin info and a dummy wave spectrum profile
      call icepack_init_wave(nfreq=nfreq,                 &
                             wave_spectrum_profile=wave_spectrum_profile, &
                             wavefreq=wavefreq, dwavefreq=dwavefreq)

      do k = 1, nfreq
          wave_spectrum(:,k) = wave_spectrum_profile(k)
      enddo

      end subroutine get_wave_spec

!=======================================================================

      end module icedrv_forcing

!=======================================================================
