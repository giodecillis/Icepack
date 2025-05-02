MODULE icedrv_phy

    USE icedrv_constants

    IMPLICIT NONE
    PRIVATE

    PUBLIC :: cp_air
    PUBLIC :: q_sat
    PUBLIC :: rho_air

    CONTAINS

    FUNCTION rho_air( ptak, pqa, ppa )
       !!-------------------------------------------------------------------------------
       !!                           ***  FUNCTION rho_air_sclr  ***
       !!
       !! ** Purpose : compute density of (moist) air using the eq. of state of the atmosphere
       !!
       !! ** Author: L. Brodeau, June 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
       !!-------------------------------------------------------------------------------
       DOUBLE PRECISION, INTENT(in) :: ptak           ! air temperature             [K]
       DOUBLE PRECISION, INTENT(in) :: pqa            ! air specific humidity   [kg/kg]
       DOUBLE PRECISION, INTENT(in) :: ppa           ! pressure in                [Pa]
       DOUBLE PRECISION             :: rho_air   ! density of moist air   [kg/m^3]
       !!-------------------------------------------------------------------------------
       rho_air = MAX( ppa / (ppR_dry*ptak * ( 1.D0 + pp_ctv0*pqa )) , 0.8D0 )
 
    END FUNCTION rho_air


    FUNCTION cp_air( pqa )
       !!-------------------------------------------------------------------------------
       !!                           ***  FUNCTION cp_air_sclr  ***
       !!
       !! ** Purpose : Compute specific heat (Cp) of moist air
       !!
       !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
       !!-------------------------------------------------------------------------------
       DOUBLE PRECISION, INTENT(in) :: pqa           ! air specific humidity         [kg/kg]
       DOUBLE PRECISION             :: cp_air   ! specific heat of moist air   [J/K/kg]
       !!-------------------------------------------------------------------------------
 
       cp_air = ppCp_dry + ppCp_vap * pqa
 
    END FUNCTION cp_air

    FUNCTION q_sat( pta, ppa,  l_ice )
       !!---------------------------------------------------------------------------------
       !!                           ***  FUNCTION q_sat_sclr  ***
       !!
       !! ** Purpose : Conputes specific humidity of air at saturation
       !!
       !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
       !!----------------------------------------------------------------------------------
       DOUBLE PRECISION :: q_sat
       DOUBLE PRECISION, INTENT(in) :: pta  !: absolute temperature of air [K]
       DOUBLE PRECISION, INTENT(in) :: ppa  !: atmospheric pressure        [Pa]
       LOGICAL,  INTENT(in), OPTIONAL :: l_ice  !: we are above ice
       DOUBLE PRECISION :: ze_s
       LOGICAL  :: lice
       !!----------------------------------------------------------------------------------
       lice = .FALSE.
       IF( PRESENT(l_ice) ) lice = l_ice
       IF( lice ) THEN
          ze_s = e_sat_ice( pta )     
          !ze_s = e_sat_ice_tet( pta )   !tetens formula
       ELSE
          ze_s = e_sat( pta ) ! Vapour pressure at saturation (Goff) :
          !ze_s = e_sat_tet( pta ) ! Vapour pressure at saturation (Teten) :
       END IF
       q_sat = pp_eps0*ze_s/(ppa - (1.D0 - pp_eps0)*ze_s)
 
    END FUNCTION q_sat

    FUNCTION e_sat( ptak )
       !!----------------------------------------------------------------------------------
       !!                   ***  FUNCTION e_sat_sclr  ***
       !!                  < SCALAR argument version >
       !! ** Purpose : water vapor at saturation in [Pa]
       !!              Based on accurate estimate by Goff, 1957
       !!
       !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
       !!
       !!----------------------------------------------------------------------------------
       DOUBLE PRECISION             ::   e_sat   ! water vapor at saturation   [kg/kg]
       DOUBLE PRECISION, INTENT(in) ::   ptak    ! air temperature                  [K]
       DOUBLE PRECISION ::   zta, ztmp   ! local scalar
       !!----------------------------------------------------------------------------------
       zta = MAX( ptak , 180.D0 )   ! air temp., prevents fpe0 errors dute to unrealistically low values over masked regions...
       ztmp = pp_tt0 / zta
       !
       ! Vapour pressure at saturation [Pa] : WMO, (Goff, 1957)
       e_sat   = 100.*( 10.**( 10.79574*(1. - ztmp) - 5.028*LOG10(zta/pp_tt0)        &
          &    + 1.50475*10.**(-4)*(1. - 10.**(-8.2969*(zta/pp_tt0 - 1.)) )  &
          &    + 0.42873*10.**(-3)*(10.**(4.76955*(1. - ztmp)) - 1.) + 0.78614) )
       !
    END FUNCTION e_sat

    FUNCTION e_sat_tet( ptak )
       !!----------------------------------------------------------------------------------
       !!                   ***  FUNCTION e_sat_sclr  ***
       !!                  < SCALAR argument version >
       !! ** Purpose : water vapor at saturation in [Pa]
       !!              Based on Teten's formula
       !!
       !!----------------------------------------------------------------------------------
       DOUBLE PRECISION             ::   e_sat_tet   ! water vapor at saturation   [kg/kg]
       DOUBLE PRECISION, INTENT(in) ::   ptak    ! air temperature                  [K]

       DOUBLE PRECISION ::   zta, ztmp   ! local scalar
       DOUBLE PRECISION, PARAMETER :: za1 = 611.21
       DOUBLE PRECISION, PARAMETER :: za3 = 17.502
       DOUBLE PRECISION, PARAMETER :: za4 = 32.19
       !!----------------------------------------------------------------------------------
       zta = MAX( ptak , 180.D0 )   ! air temp., prevents fpe0 errors dute to unrealistically low values over masked regions...
       ztmp = pp_tt0 / zta
       !
       ! Vapour pressure at saturation [Pa]
       e_sat_tet   = za1 * EXP(za3 * (ptak - pp_tt0)/(ptak - za4)) 
       !
    END FUNCTION e_sat_tet


    FUNCTION e_sat_ice(ptak)
       !!---------------------------------------------------------------------------------
       !! Same as "e_sat" but over ice rather than water!
       !!---------------------------------------------------------------------------------
       DOUBLE PRECISION             :: e_sat_ice !: vapour pressure at saturation in presence of ice [Pa]
       DOUBLE PRECISION, INTENT(in) :: ptak
       !!
       DOUBLE PRECISION :: zta, zle, ztmp
       !!---------------------------------------------------------------------------------
       zta = MAX( ptak , 180.D0 )   ! air temp., prevents fpe0 errors dute to unrealistically low values over masked regions...
       ztmp = pp_tt0/zta
       !!
       zle  = rAg_i*(ztmp - 1.D0) + rBg_i*LOG10(ztmp) + rCg_i*(1.D0 - zta/pp_tt0) + rDg_i
       !!
       e_sat_ice = 100.D0 * 10.D0**zle
 
    END FUNCTION e_sat_ice

    FUNCTION e_sat_ice_tet( ptak )
       !!----------------------------------------------------------------------------------
       !!                   ***  FUNCTION e_sat_sclr  ***
       !!                  < SCALAR argument version >
       !! ** Purpose : water vapor at saturation in [Pa]
       !!              Based on Teten's formula
       !!
       !!----------------------------------------------------------------------------------
       DOUBLE PRECISION             ::   e_sat_ice_tet   ! water vapor at saturation   [kg/kg]
       DOUBLE PRECISION, INTENT(in) ::   ptak    ! air temperature                  [K]

       DOUBLE PRECISION ::   zta, ztmp   ! local scalar
       DOUBLE PRECISION, PARAMETER :: za1 = 611.21
       DOUBLE PRECISION, PARAMETER :: za3 = 22.587
       DOUBLE PRECISION, PARAMETER :: za4 = -0.7
       !!----------------------------------------------------------------------------------
       zta = MAX( ptak , 180.D0 )   ! air temp., prevents fpe0 errors dute to unrealistically low values over masked regions...
       ztmp = pp_tt0 / zta
       !
       ! Vapour pressure at saturation [Pa] 
       e_sat_ice_tet   = za1 * EXP(za3 * (ptak - pp_tt0)/(ptak - za4)) 
       !
    END FUNCTION e_sat_ice_tet


END MODULE
