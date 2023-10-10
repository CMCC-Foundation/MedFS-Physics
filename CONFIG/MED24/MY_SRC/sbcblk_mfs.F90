MODULE sbcblk_mfs
   !!======================================================================
   !!                       ***  MODULE  sbcblk_mfs  ***
   !! Ocean forcing:  momentum, heat and freshwater flux formulation
   !!=====================================================================
   !! History :  3.3  !   2010-05 (P. Oddo) Original Code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_blk_mfs  : bulk formulation as ocean surface boundary condition
   !!                   (forced mode, mfs bulk formulae)
   !!   blk_oce_mfs  : ocean: computes momentum, heat and freshwater fluxes
   !!   turb_cd_2z   : Computes iturb surface drag at 2m
   !!   psi_mc           : universal profile stability function for momentum
   !!   psi_hc           : universal profile stability function for temperature and humidity
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE fldread         ! read input fields
   USE sbc_oce         ! Surface boundary condition: ocean fields
   USE iom             ! I/O manager library
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distribued memory computing library
   USE wrk_nemo        ! work arrays
   USE timing          ! Timing
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE prtctl          ! Print control
#ifndef  key_wwiii
   USE sbcwave,ONLY : cdn_wave !wave module
#endif

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_blk_mfs       ! routine called in sbcmod module
   PUBLIC   turb_cd_2z        ! routine called in sbcblk_mfs module
      
   INTEGER , PARAMETER ::   jpfld   = 7         ! maximum number of files to read 
   INTEGER , PARAMETER ::   jp_wndi = 1         ! index of 10m wind velocity (i-component) (m/s) at T-point
   INTEGER , PARAMETER ::   jp_wndj = 2         ! index of 10m wind velocity (j-component) (m/s) at T-point
   INTEGER , PARAMETER ::   jp_clc  = 3         ! index of total cloud cover               ( % )
   INTEGER , PARAMETER ::   jp_msl  = 4         ! index of mean sea level pressure         (Pa)
   INTEGER , PARAMETER ::   jp_tair = 5         ! index of 10m air temperature             (Kelvin)
   INTEGER , PARAMETER ::   jp_rhm  = 6         ! index of dew point temperature           (Kelvin)
   INTEGER , PARAMETER ::   jp_prec = 7         ! index of total precipitation (rain+snow) (Kg/m2/s)
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf ! structure of input fields (file informations, fields read)
         
   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.2 , LOCEAN-IPSL (2009) 
   !! $Id: sbcblk_mfs.F90 5215 2015-04-15 16:11:56Z nicolasmartin $
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS


   SUBROUTINE sbc_blk_mfs( kt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_blk_mfs  ***
      !!                   
      !! ** Purpose :   provide at each time step the surface ocean fluxes
      !!      (momentum, heat, freshwater, runoff is added later in the code) 
      !!
      !! ** Method  : (1) READ Atmospheric data from NetCDF files:
      !!      the 10m wind velocity (i-component) (m/s) at T-point
      !!      the 10m wind velocity (j-component) (m/s) at T-point
      !!      the 2m Dew point Temperature        (k)
      !!      the Cloud COver                     (%)
      !!      the 2m air temperature              (Kelvin)
      !!      the Mean Sea Level Preesure         (hPa)
      !!      the Climatological Precipitation    (kg/m2/s)
      !!              (2) CALL blk_oce_mfs
      !!
      !!      Computes:
      !!      Solar Radiation using Reed formula (1975, 1977)
      !!      Net Long wave radiation using Bignami et al. (1995)
      !!      Latent and Sensible heat using Kondo (1975)
      !!      Drag coeff using Hllerman and Rosenstein (1983)
      !!      C A U T I O N : never mask the surface stress fields
      !!                      the stress is assumed to be in the mesh referential
      !!                      i.e. the (i,j) referential
      !!
      !! ** Action  :   defined at each time-step at the air-sea interface
      !!              - utau, vtau  i- and j-component of the wind stress
      !!              - taum        wind stress module at T-point
      !!              - wndm        10m wind module at T-point over free ocean or leads in presence of sea-ice
      !!              - qns, qsr    non-slor and solar heat flux
      !!              - emp         evaporation minus precipitation
      !!----------------------------------------------------------------------
      REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::  sh_now   ! specific humidity at T-point 
      REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::  catm     ! Cover 
      REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::  alonl    ! Lon 
      REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::  alatl    ! Lat 
      REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::  gsst     ! SST 
      !!---------------------------------------------------------------------
      !! Local fluxes variables
      !!---------------------------------------------------------------------
      REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::  qbw     ! Net Long wave 
      REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::  ha      ! Sesnible 
      REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::  elat    ! Latent 
      REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::  evap    ! evaporation rate 
      REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::  prec    ! precipitation rate

      INTEGER, INTENT( in  ) ::   kt   ! ocean time step
      !!
      INTEGER  :: ierror                          ! return error code
      INTEGER  :: ifpr     ! dummy loop indice
      INTEGER  :: jj,ji    ! dummy loop arguments
      INTEGER  ::   ios    ! Local integer output status for namelist read
      REAL(wp) :: act_hour
      !!--------------------------------------------------------------------
      !! Variables for specific humidity computation
      !!--------------------------------------------------------------------
      REAL(wp) :: onsea,par1,par2
      DATA onsea,par1,par2 / 0.98, 640380., -5107.4 /
      !!                      par1 [Kg/m3], par2 [K]

      CHARACTER(len=100) ::  cn_dir                           ! Root directory for location of Atmospheric forcing files
      TYPE(FLD_N), DIMENSION(jpfld) ::   slf_i                ! array of namelist informations on the fields to read
      TYPE(FLD_N) ::   sn_wndi, sn_wndj, sn_clc, sn_msl       ! informations about the fields to be read
      TYPE(FLD_N) ::   sn_tair , sn_rhm, sn_prec              !   "                                 "
      !!---------------------------------------------------------------------
      NAMELIST/namsbc_mfs/ cn_dir ,                                          &
         &                  sn_wndi , sn_wndj, sn_clc   , sn_msl ,           &
         &                  sn_tair , sn_rhm , sn_prec 
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('sbc_blk_mfs')
      !
      !                                         ! ====================== !
      IF( kt == nit000 ) THEN                   !  First call kt=nit000  !
         !                                      ! ====================== !
         ALLOCATE( sh_now(jpi,jpj), catm(jpi,jpj), alonl(jpi,jpj), alatl(jpi,jpj),     &
         &        gsst(jpi,jpj),  qbw(jpi,jpj),    ha(jpi,jpj),  elat(jpi,jpj),     &
         &        evap(jpi,jpj), prec(jpi,jpj),    STAT=ierror )

         IF( ierror /= 0 )   CALL ctl_warn('sbc_blk_mfs: failed to allocate arrays')

         REWIND( numnam_ref )              ! Namelist namsbc_msf in reference namelist : MFS files
         READ  ( numnam_ref, namsbc_mfs, IOSTAT = ios, ERR = 901)
901      IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsbc_mfs in reference namelist', lwp )

         REWIND( numnam_cfg )              ! Namelist namsbc_msf in configuration namelist : MFS files
         READ  ( numnam_cfg, namsbc_mfs, IOSTAT = ios, ERR = 902 )
902      IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsbc_mfs in configuration namelist', lwp )
         IF(lwm) WRITE ( numond, namsbc_mfs )
         !
         ! store namelist information in an array
         slf_i(jp_wndi) = sn_wndi   ;   slf_i(jp_wndj) = sn_wndj
         slf_i(jp_clc ) = sn_clc    ;   slf_i(jp_msl ) = sn_msl
         slf_i(jp_tair) = sn_tair   ;   slf_i(jp_rhm)  = sn_rhm
         slf_i(jp_prec) = sn_prec   ;  

         IF(lwp) WRITE(numout,*) 'sn_prec frequency is ',sn_prec%nfreqh

         !
         ALLOCATE( sf(jpfld), STAT=ierror )         ! set sf structure
         IF( ierror > 0 ) THEN
            CALL ctl_stop( 'sbc_blk_mfs: unable to allocate sf structure' )   ;   RETURN
         ENDIF
         DO ifpr= 1, jpfld
            ALLOCATE( sf(ifpr)%fnow(jpi,jpj,1) )
            IF( slf_i(ifpr)%ln_tint ) ALLOCATE( sf(ifpr)%fdta(jpi,jpj,1,2) )
         END DO
         ! fill sf with slf_i and control print
         CALL fld_fill( sf, slf_i, cn_dir,'sbc_blk_mfs','bulk formulation for ocean SBC', 'namsbc_mfs' )
            !
      ENDIF

      CALL fld_read( kt, nn_fsbc, sf )                   ! input fields provided at the current time-step

      catm(:,:)   = 0.0    ! initializze cloud cover variable
      sh_now(:,:) = 0.0    ! initializze specifif humidity variable

      DO jj = 1, jpj
         DO ji = 1, jpi   ! vector opt.

         ! Calculate Specific Humidity 
         !-------------------------------------------------
         sh_now(ji,jj) = (1/1.22) * onsea * par1 * EXP(par2/sf(jp_rhm)%fnow(ji,jj,1))

         ! Normalize Clouds
         !-------------------------------------------------
         catm(ji,jj)   = max(0.0,min(1.0,sf(jp_clc)%fnow(ji,jj,1)*0.01))

         ! Convert precipitations from Total prec [m] to
         ! rain fall rate in kg/m2/s
         ! rainfall rate = tot prec [m]* water dens [1000kg/m3] / time period
         ! [s]

         if (sf(jp_prec)%fnow(ji,jj,1).lt.0) sf(jp_prec)%fnow(ji,jj,1)=0.0
         sf(jp_prec)%fnow(ji,jj,1) = sf(jp_prec)%fnow(ji,jj,1) *1000/(sn_prec%nfreqh*3600)

         END DO
      END DO
   


      ! wind module at 10m 
      !--------------------------------------------------
      wndm(:,:) = SQRT(  sf(jp_wndi)%fnow(:,:,1) * sf(jp_wndi)%fnow(:,:,1)   &
           &             + sf(jp_wndj)%fnow(:,:,1) * sf(jp_wndj)%fnow(:,:,1)  )



      ! Some conv for fluxes computation
      !-------------------------------------------------
      alonl(:,:) = glamt(:,:) * rad
      alatl(:,:) = gphit(:,:) * rad
      gsst(:,:)  = tsn(:,:,1,jp_tem)  * tmask(:,:,1)

      IF( MOD( kt - 1, nn_fsbc ) == 0 ) THEN

         ! Force to zero the output of fluxes 
         !-------------------------------------------------
          qsr(:,:)  = 0.0 ; qbw(:,:)  = 0.0 ; 
          ha(:,:)   = 0.0 ; elat(:,:) = 0.0 ; 
          evap(:,:) = 0.0 ; utau(:,:) = 0.0 ; 
          vtau(:,:) = 0.0 ; prec(:,:) = 0.0

          CALL lbc_lnk( sf(jp_wndi)%fnow(:,:,1), 'T', -1. )
          CALL lbc_lnk( sf(jp_wndj)%fnow(:,:,1), 'T', -1. )

          act_hour = (( nsec_year / rday ) - INT (nsec_year / rday)) * rjjhh

          CALL fluxes_mfs(alatl,alonl,act_hour,                                &     ! input static
                            gsst(:,:),sf(jp_tair)%fnow(:,:,1),sh_now(:,:),     &     ! input dynamic
                            sf(jp_wndi)%fnow(:,:,1), sf(jp_wndj)%fnow(:,:,1),  &     ! input dynamic
                            sf(jp_msl)%fnow(:,:,1) , catm(:,:) ,               &     ! input dynamic
                            qsr,qbw,ha,elat,evap,utau,vtau)                          ! output

         ! Shortwave radiation
         !--------------------------------------------------
          qsr(:,:) = qsr(:,:) * tmask(:,:,1)

         ! total non solar heat flux over water
         !--------------------------------------------------
          qns(:,:) = -1. * ( qbw(:,:) + ha(:,:) + elat(:,:) )
          qns(:,:) = qns(:,:)*tmask(:,:,1)

         ! mask the wind module at 10m 
         !--------------------------------------------------
          wndm(:,:) = wndm(:,:) * tmask(:,:,1)

         !   wind stress module (taum) into T-grid
         !--------------------------------------------------
          taum(:,:) = SQRT( utau(:,:) * utau(:,:) + vtau(:,:) * vtau(:,:) ) * tmask(:,:,1)

          CALL lbc_lnk( taum, 'T', 1. )

         ! Interpolate utau, vtau into the grid_V and grid_V
         !-------------------------------------------------
         ! Note the use of 0.5*(2-umask) in order to unmask the stress along coastlines
         ! Note the use of MAX(tmask(i,j),tmask(i+1,j) is to mask tau over ice shelves
          DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1
               utau(ji,jj) = 0.5 * ( 2. - umask(ji,jj,1) ) * ( utau(ji,jj) * tmask(ji,jj,1) &
               &                                + utau(ji+1,jj) * tmask(ji+1,jj,1) )        &
               &                 * MAX(tmask(ji,jj,1),tmask(ji+1,jj  ,1))
               vtau(ji,jj) = 0.5 * ( 2. - vmask(ji,jj,1) ) * ( vtau(ji,jj) * tmask(ji,jj,1) &
               &                                + vtau(ji,jj+1) * tmask(ji,jj+1,1) )        &
               &                 * MAX(tmask(ji,jj,1),tmask(ji  ,jj+1,1))
            END DO
          END DO

          CALL lbc_lnk( utau(:,:), 'U', -1. )
          CALL lbc_lnk( vtau(:,:), 'V', -1. )

         ! for basin budget and cooerence
         !--------------------------------------------------
!CDIR COLLAPSE
           emp (:,:) = evap(:,:) - sf(jp_prec)%fnow(:,:,1) * tmask(:,:,1)
           prec(:,:) = sf(jp_prec)%fnow(:,:,1) * tmask(:,:,1)
!CDIR COLLAPSE

           CALL iom_put( "qlw_oce",  -qbw  )                 ! output downward longwave heat over the ocean
           CALL iom_put( "qsb_oce",   -ha   )                ! output downward sensible heat over the ocean
           CALL iom_put( "qla_oce",   -elat )                ! output downward latent   heat over the ocean
           CALL iom_put( "qns_oce",   qns  )                 ! output downward non solar heat over the ocean
           CALL iom_put( "evap_oce",  evap  )                ! output water evaporation flux
           CALL iom_put( "prec_oce",  prec  )                ! output precipitation flux

      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('sbc_blk_mfs')
      !
   END SUBROUTINE sbc_blk_mfs
  
 
   SUBROUTINE fluxes_mfs(alat,alon,hour,                               &
        sst,tnow,shnow,unow,vnow,mslnow,cldnow,qsw,qbw,ha,elat,        &
        evap,taux,tauy)
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE fluxes_mfs  ***
      !!
      !! --- it provides SURFACE HEAT and MOMENTUM FLUXES in MKS :
      !!
      !!  1)   Water flux (WFLUX)                 [ watt/m*m ]
      !!  2)   Short wave flux (QSW)              [ watt/m*m ] Reed 1977
      !!  3)   Long wave flux backward (QBW)      [ watt/m*m ]
      !!  4)   Latent heat of evaporation (ELAT)  [ watt/m*m ]
      !!  5)   Sensible heat flux   (HA)          [ watt/m*m ]
      !!  6)   Wind stress x-component   (TAUX)   [ newton/m*m ]
      !!  7)   Wind stress y-component   (TAUY)   [ newton/m*m ]
      !!
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in   ) :: hour
      REAL(wp), INTENT(in   ), DIMENSION (:,:) :: sst, unow, alat , alon
      REAL(wp), INTENT(in   ), DIMENSION (:,:) :: vnow, cldnow, mslnow
      REAL(wp), INTENT(out  ), DIMENSION (:,:) :: qsw, qbw, ha, elat
      REAL(wp), INTENT(out  ), DIMENSION (:,:) :: evap,taux,tauy
      REAL(wp), INTENT(inout), DIMENSION (:,:) :: tnow , shnow

      INTEGER :: ji,jj 
      REAL(wp)  :: wair, vtnow, ea, deltemp, s, stp , fh , fe
      REAL(wp)  :: esre, cseep

      REAL(wp), DIMENSION (:,:), POINTER ::   rspeed, cdx, ce, shms
      REAL(wp), DIMENSION (:,:), POINTER ::   rhom, sstk, ch, rel_windu, rel_windv
      !!----------------------------------------------------------------------
      !!     coefficients ( in MKS )  :
      !!----------------------------------------------------------------------

      REAL(wp), PARAMETER ::  ps = 1013.    ! --- surface air pressure
      REAL(wp), PARAMETER ::  expsi=0.622   ! --- expsi
      REAL(wp), PARAMETER ::  rd=287.       ! --- dry air gas constant
      REAL(wp), PARAMETER ::  cp=1005.      ! --- specific heat capacity
      REAL(wp), PARAMETER ::  onsea=0.98    ! --- specific humidity factors
      REAL(wp), PARAMETER ::  par1=640380.  ! [Kg/m3]
      REAL(wp), PARAMETER ::  par2=-5107.4  ! [K]

      !---------------------------------------------------------------------
      !--- define Kondo parameters
      !---------------------------------------------------------------------

      REAL(wp), DIMENSION(5) :: a_h = (/0.0,0.927,1.15,1.17,1.652/)
      REAL(wp), DIMENSION(5) :: a_e = (/0.0,0.969,1.18,1.196,1.68/)
      REAL(wp), DIMENSION(5) :: b_h = (/1.185,0.0546,0.01,0.0075,-0.017/)
      REAL(wp), DIMENSION(5) :: b_e = (/1.23,0.0521,0.01,0.008,-0.016/)
      REAL(wp), DIMENSION(5) :: c_h = (/0.0,0.0,0.0,-0.00045,0.0/)
      REAL(wp), DIMENSION(5) :: c_e = (/0.0,0.0,0.0,-0.0004,0.0/)
      REAL(wp), DIMENSION(5) :: p_h = (/-0.157,1.0,1.0,1.0,1.0/)
      REAL(wp), DIMENSION(5) :: p_e = (/-0.16,1.0,1.0,1.0,1.0/)
      INTEGER :: kku                        !index varing with wind speed
      !
      IF( nn_timing == 1 )  CALL timing_start('fluxes_mfs')
      !
      CALL wrk_alloc( jpi,jpj, rspeed, cdx, ce, shms )
      CALL wrk_alloc( jpi,jpj, rhom, sstk, ch, rel_windu, rel_windv )

      !!----------------------------------------------------------------------
      !! ------------------ (i)      short wave
      !!----------------------------------------------------------------------

      CALL qshort(hour,alat,alon,cldnow,qsw)

      rel_windu(:,:) = 0.0_wp
      rel_windv(:,:) = 0.0_wp

!       DO jj = 2, jpj
!          DO ji = fs_2, jpi
!           rel_windu(ji,jj) = unow(ji,jj) - 0.5_wp * ( ssu_m(ji-1,jj) + ssu_m(ji,jj) )
!           rel_windv(ji,jj) = vnow(ji,jj) - 0.5_wp * ( ssv_m(ji,jj-1) + ssv_m(ji,jj) )
!          END DO
!       END DO

      ! Relative wind speed on T grid

       DO jj = 2, jpj
          DO ji = fs_2, jpi
           rel_windu(ji,jj) = unow(ji,jj) - ( ssu_m(ji-1,jj) + ssu_m(ji,jj) ) &
                            &  /max(1.0, umask(ji-1,jj,1) + umask(ji,jj,1))
           rel_windv(ji,jj) = vnow(ji,jj) - ( ssv_m(ji,jj-1) + ssv_m(ji,jj) ) &
                            &  /max(1.0, vmask(ji,jj-1,1) + vmask(ji,jj,1))
          END DO
       END DO


       CALL lbc_lnk( rel_windu(:,:), 'T', -1. )
       CALL lbc_lnk( rel_windv(:,:), 'T', -1. )


       rspeed(:,:)= SQRT(rel_windu(:,:)*rel_windu(:,:)   &
         &                   + rel_windv(:,:)*rel_windv(:,:)) 

       sstk(:,:) = sst(:,:) + rtt                          !- SST data in Kelvin degrees
       shms(:,:) = (1/1.22)*onsea*par1*EXP(par2/sstk(:,:)) !- Saturation Specific Humidity

      !!----------------------------------------------------------------------
      !! ------------------ (ii)    net long wave
      !!----------------------------------------------------------------------

      DO jj = 1, jpj
         DO ji = 1, jpi
            wair = shnow(ji,jj) / (1 - shnow(ji,jj))    ! mixing ratio of the air (Wallace and Hobbs)
            vtnow = (tnow(ji,jj)*(expsi+wair))/(expsi*(1.+wair))   ! virtual temperature of air
            rhom(ji,jj) = 100.*(ps/rd)/vtnow                       ! density of the moist air
            ea   = (wair  / (wair  + 0.622 )) * (mslnow(ji,jj)/100.0)  !msl from hPa to Pa

            qbw(ji,jj) = emic*stefan*( sstk(ji,jj)**4. )                    &
                 - ( stefan*( tnow(ji,jj)**4. ) * ( 0.653 + 0.00535*ea ) )  &
                   * ( 1. + 0.1762*( cldnow(ji,jj)**2. ) )

         END DO
      END DO

      DO jj = 1, jpj
         DO ji = 1, jpi
      !!----------------------------------------------------------------------
      !! ------------------ (iii)   sensible heat
      !!----------------------------------------------------------------------

      !! --- calculates the term :      ( Ts - Ta )
      !!----------------------------------------------------------------------

            deltemp = sstk(ji,jj) - tnow (ji,jj)

      !!----------------------------------------------------------------------
      !! --- variable turbulent exchange coefficients ( from Kondo 1975 )
      !! --- calculate the Neutral Transfer Coefficent using an empiric formula
      !! --- by Kondo et al. Then it applies the diabatic approximation.
      !!----------------------------------------------------------------------

            s = deltemp/(wndm(ji,jj)**2.)   !! --- calculate S 
            stp = s*abs(s)/(abs(s)+0.01)    !! --- calculate the Stability Parameter 

      !!----------------------------------------------------------------------
      !! --- for stable condition (sst-t_air < 0):
      !!----------------------------------------------------------------------

            IF (s.lt.0. .and. ((stp.gt.-3.3).and.(stp.lt.0.))) THEN
                fh = 0.1_wp+0.03_wp*stp+0.9_wp*exp(4.8_wp*stp)
                fe = fh
            ELSE IF (s.lt.0. .and. stp.le.-3.3) THEN
                fh = 0._wp
                fe = fh
            ELSE                                       ! --- for unstable condition 
                fh = 1.0_wp+0.63_wp*sqrt(stp)
                fe = fh
            ENDIF

      !!----------------------------------------------------------------------
      !! --- calculate the coefficient CH,CE,CD
      !!----------------------------------------------------------------------

            IF (wndm(ji,jj) >= 0. .AND. wndm(ji,jj) <= 2.2)       THEN
                kku=1
            ELSE IF (wndm(ji,jj) > 2.2 .AND. wndm(ji,jj) <= 5.0)  THEN
                kku=2
            ELSE IF (wndm(ji,jj) > 5.0 .AND. wndm(ji,jj) <= 8.0)  THEN
                kku=3
            ELSE IF (wndm(ji,jj) > 8.0 .AND. wndm(ji,jj) <= 25.0) THEN
                kku=4
            ELSE IF (wndm(ji,jj) > 25.0 )                         THEN
                kku=5
            ENDIF

            ch(ji,jj) = ( a_h(kku) + b_h(kku) * wndm(ji,jj) ** p_h(kku)      &
                        + c_h(kku) * (wndm(ji,jj)-8 ) **2) * fh

            ce(ji,jj) = ( a_e(kku) + b_e(kku) * wndm(ji,jj) ** p_e(kku)      &
                        + c_e(kku) * (wndm(ji,jj)-8 ) **2) * fe

            ch(ji,jj) = ch(ji,jj) / 1000.0
            ce(ji,jj) = ce(ji,jj) / 1000.0

            IF (wndm(ji,jj)<0.3) THEN
               ch(ji,jj) = 1.3e-03 * fh
               ce(ji,jj) = 1.5e-03 * fe
            ELSE IF(wndm(ji,jj)>50.0) THEN
               ch(ji,jj) = 1.25e-03 * fh
               ce(ji,jj) = 1.30e-03 * fe
            ENDIF

      !!----------------------------------------------------------------------
      !! --- calculates the SENSIBLE HEAT FLUX in MKS ( watt/m*m )
      !!----------------------------------------------------------------------

            HA(ji,jj) = rhom(ji,jj)*cp*ch(ji,jj)*wndm(ji,jj)*deltemp

      !!----------------------------------------------------------------------
      !! ------------------ (iv)  latent heat
      !! --- calculates the LATENT HEAT FLUX  ( watt/m*m )
      !! --- ELAT = L*rho*Ce*|V|*[qs(Ts)-qa(t2d)]
      !!----------------------------------------------------------------------

            esre  = shms(ji,jj) - shnow(ji,jj)   ! --- calculates the term : qs(Ta)-qa(t2d)

            cseep = ce(ji,jj) * wndm(ji,jj) * esre     ! --- calculates the term : Ce*|V|*[qs(Ts)-qa(t2d)]

            evap(ji,jj) = (cseep * rhom(ji,jj))  ! in [kg/m2/sec] !! --- calculates the EVAPORATION RATE [m/yr]

            elat(ji,jj) = rhom(ji,jj) * cseep * heatlat(sst(ji,jj))

      !!----------------------------------------------------------------------
      !! --- calculates the Drag Coefficient
      !!----------------------------------------------------------------------

      !!----------------------------------------------------------------------
      !! --- deltemp should be (Ts - Ta) in the formula estimating
      !! --- drag coefficient
      !!----------------------------------------------------------------------

              IF ( .NOT. ln_cdgw ) THEN
                 cdx(ji,jj) = cd_HR(wndm(ji,jj),deltemp)
              ENDIF

          END DO
      END DO

      IF (ln_cdgw ) THEN
              CALL turb_cd_2z(2.,10.,sstk,tnow+2*0.0098,shms,shnow,rspeed,cdx)
      ENDIF
   

      !!----------------------------------------------------------------------
      !! --- calculates the wind stresses in MKS ( newton/m*m )
      !! ---            taux= rho*Cd*|V|u      tauy= rho*Cd*|V|v
      !!----------------------------------------------------------------------

      taux(:,:)= rhom(:,:) * cdx(:,:) * rspeed(:,:) * rel_windu(:,:)
      tauy(:,:)= rhom(:,:) * cdx(:,:) * rspeed(:,:) * rel_windv(:,:)

      CALL wrk_dealloc( jpi,jpj, rspeed, cdx, ce, shms )
      CALL wrk_dealloc( jpi,jpj, rhom, sstk, ch, rel_windu, rel_windv )
      !
      IF( nn_timing == 1 )  CALL timing_stop('fluxes_mfs')
      !
   END SUBROUTINE fluxes_mfs

   SUBROUTINE turb_cd_2z( zt, zu, sst, T_zt, q_sat, q_zt, dU, Cd )
      !!----------------------------------------------------------------------
      !!               ***  modified from ROUTINE  turb_core  ***
      !!
      !! ** Purpose :   Computes turbulent transfert coefficients of surface
      !!                fluxes according to Large & Yeager (2004) and Large & Yeager (2008)
      !!                If relevant (zt /= zu), adjust temperature and humidity from height zt to zu
      !!
      !! ** Method : Monin Obukhov Similarity Theory
      !!             + Large & Yeager (2004,2008) closure: CD_n10 = f(U_n10)
      !!
      !! ** References :   Large & Yeager, 2004 / Large & Yeager, 2008
      !!
      !! ** update: Laurent Brodeau, June 2014:
      !!    - handles both cases zt=zu and zt/=zu
      !!    - optimized: less 2D arrays allocated and less operations
      !!    - better first guess of stability by checking air-sea difference of virtual temperature
      !!       rather than temperature difference only...
      !!      Large & Yeager 2008. Drag-coefficient reduction for Cyclone conditions
      !!    - using code-wide physical constants defined into "phycst.mod" rather than redifining them
      !!      => 'vkarmn' and 'grav'
      !! ** Last update: Emanuela Clementi, Nov 2018:
      !!    - removed case not ln_cdgw 
      !!    - removed function cd_neutral_10m (not needed)
      !!    - removed Ce, Ch, ,t_zu, sh_zu from output fields: not needed for mfs
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in   )                     ::   zt       ! height for T_zt and q_zt                  [m]
      REAL(wp), INTENT(in   )                     ::   zu       ! height for dU                             [m]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   sst      ! sea surface temperature              [Kelvin]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   T_zt     ! potential air temperature            [Kelvin]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   q_sat    ! sea surface specific humidity         [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   q_zt     ! specific air humidity                 [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   dU       ! relative wind module at zu              [m/s]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Cd       ! transfer coefficient for momentum       (tau)
      !
      INTEGER ::   j_itt
      INTEGER , PARAMETER ::   nb_itt = 5       ! number of itterations
      LOGICAL ::   l_zt_equal_zu = .FALSE.      ! if q and t are given a different height than U
      !
      REAL(wp), DIMENSION(:,:), POINTER ::   U_zu          ! relative wind at zu [m/s]
      REAL(wp), DIMENSION(:,:), POINTER ::   Ce_n10        ! 10m neutral latent coefficient
      REAL(wp), DIMENSION(:,:), POINTER ::   Ch_n10        ! 10m neutral sensible coefficient
      REAL(wp), DIMENSION(:,:), POINTER ::   Ch            ! transfer coefficient for sensible heat (Q_sens)
      REAL(wp), DIMENSION(:,:), POINTER ::   Ce            ! transfert coefficient for evaporation   (Q_lat)
      REAL(wp), DIMENSION(:,:), POINTER ::   T_zu          ! air temp. shifted at zu                     [K]
      REAL(wp), DIMENSION(:,:), POINTER ::   q_zu          ! spec. hum.shifted at zu               [kg/kg]
      REAL(wp), DIMENSION(:,:), POINTER ::   sqrt_Cd_n10   ! root square of Cd_n10
      REAL(wp), DIMENSION(:,:), POINTER ::   sqrt_Cd       ! root square of Cd
      REAL(wp), DIMENSION(:,:), POINTER ::   zeta_u        ! stability parameter at height zu
      REAL(wp), DIMENSION(:,:), POINTER ::   zeta_t        ! stability parameter at height zt
      REAL(wp), DIMENSION(:,:), POINTER ::   zpsi_h_u, zpsi_m_u
      REAL(wp), DIMENSION(:,:), POINTER ::   ztmp0, ztmp1, ztmp2
      REAL(wp), DIMENSION(:,:), POINTER ::   stab          ! 1st stability test integer
      !!----------------------------------------------------------------------

      IF( nn_timing == 1 )  CALL timing_start('turb_cd_2z')

      CALL wrk_alloc( jpi,jpj, U_zu, Ce_n10, Ch_n10, sqrt_Cd_n10, sqrt_Cd )
      CALL wrk_alloc( jpi,jpj, Ch, Ce, T_zu, q_zu )
      CALL wrk_alloc( jpi,jpj, zeta_u, stab )
      CALL wrk_alloc( jpi,jpj, zpsi_h_u, zpsi_m_u, ztmp0, ztmp1, ztmp2 )

      l_zt_equal_zu = .FALSE.
      IF( ABS(zu - zt) < 0.01 ) l_zt_equal_zu = .TRUE.    ! testing "zu == zt" is risky with double precision

      IF( .NOT. l_zt_equal_zu )   CALL wrk_alloc( jpi,jpj, zeta_t )

      U_zu = MAX( 0.5 , dU )   !  relative wind speed at zu (normally 10m), we don't want to fall under 0.5 m/s

#ifdef  key_wwiii
      ! emanuela needed for online wave coupling
      dT_AS=T_zt-sst
#endif

      !! First guess of stability:
      ztmp0 = T_zt*(1. + 0.608*q_zt) - sst*(1. + 0.608*q_sat) ! air-sea difference of virtual pot. temp. at zt
      stab  = 0.5 + sign(0.5,ztmp0)                           ! stab = 1 if dTv  > 0  => STABLE, 0 if unstable

      !! Neutral coefficients at 10m read from wave output:
      cdn_wave(:,:) = cdn_wave(:,:) + rsmall * ( 1._wp - tmask(:,:,1) )
      ztmp0   (:,:) = cdn_wave(:,:)

      sqrt_Cd_n10 = SQRT( ztmp0 )
      Ce_n10  = 1.e-3*( 34.6 * sqrt_Cd_n10 )
      Ch_n10  = 1.e-3*sqrt_Cd_n10*(18.*stab + 32.7*(1. - stab))

      !! Initializing transf. coeff. with their first guess neutral equivalents:
      Cd = ztmp0   ;   Ce = Ce_n10   ;   Ch = Ch_n10   ;   sqrt_Cd = sqrt_Cd_n10

      !! Initializing values at z_u with z_t values:
      T_zu = T_zt   ;   q_zu = q_zt

      !!  * Now starting iteration loop
      DO j_itt=1, nb_itt
         !
         ztmp1 = T_zu - sst   ! Updating air/sea differences
         ztmp2 = q_zu - q_sat

         ! Updating turbulent scales :   (L&Y 2004 eq. (7))
         ztmp1  = Ch/sqrt_Cd*ztmp1    ! theta*
         ztmp2  = Ce/sqrt_Cd*ztmp2    ! q*

         ztmp0 = T_zu*(1. + 0.608*q_zu) ! virtual potential temperature at zu

         ! Estimate the inverse of Monin-Obukov length (1/L) at height zu:
         ztmp0 =  (vkarmn*grav/ztmp0*(ztmp1*(1.+0.608*q_zu) + 0.608*T_zu*ztmp2)) / (Cd*U_zu*U_zu)
         !                                                           ( Cd*U_zu*U_zu is U*^2 at zu)
         !! Stability parameters :
         zeta_u   = zu*ztmp0   ;  zeta_u = sign( min(abs(zeta_u),10.0), zeta_u )
         zpsi_h_u = psi_hc( zeta_u )
         zpsi_m_u = psi_mc( zeta_u )

         !! Shifting temperature and humidity at zu (L&Y 2004 eq. (9b-9c))
         IF ( .NOT. l_zt_equal_zu ) THEN
            zeta_t = zt*ztmp0 ;  zeta_t = sign( min(abs(zeta_t),10.0), zeta_t )
            stab = LOG(zu/zt) - zpsi_h_u + psi_hc(zeta_t)  ! stab just used as temp array!!!
            T_zu = T_zt + ztmp1/vkarmn*stab    ! ztmp1 is still theta*
            q_zu = q_zt + ztmp2/vkarmn*stab    ! ztmp2 is still q*
            q_zu = max(0., q_zu)
         END IF

         sqrt_Cd = vkarmn / ( vkarmn / sqrt_Cd_n10 - zpsi_m_u )
         Cd      = sqrt_Cd * sqrt_Cd
         !
         ztmp0 = (LOG(zu/10.) - zpsi_h_u) / vkarmn / sqrt_Cd_n10
         ztmp2 = sqrt_Cd / sqrt_Cd_n10 
         ztmp1 = 1. + Ch_n10*ztmp0
         Ch  = Ch_n10*ztmp2 / ztmp1  ! L&Y 2004 eq. (10b)
         !
         ztmp1 = 1. + Ce_n10*ztmp0
         Ce  = Ce_n10*ztmp2 / ztmp1  ! L&Y 2004 eq. (10c)
         !
      END DO

      CALL wrk_dealloc( jpi,jpj, U_zu, Ce_n10, Ch_n10, sqrt_Cd_n10, sqrt_Cd )
      CALL wrk_dealloc( jpi,jpj, zeta_u, stab )
      CALL wrk_dealloc( jpi,jpj, Ch, Ce, T_zu, q_zu )
      CALL wrk_dealloc( jpi,jpj, zpsi_h_u, zpsi_m_u, ztmp0, ztmp1, ztmp2 )

      IF( .NOT. l_zt_equal_zu ) CALL wrk_dealloc( jpi,jpj, zeta_t )

      IF( nn_timing == 1 )  CALL timing_stop('turb_cd_2z')
      !
   END SUBROUTINE turb_cd_2z


   FUNCTION psi_mc(pta)   !! Psis, L&Y 2004 eq. (8c), (8d), (8e)
      !-------------------------------------------------------------------------------
      ! universal profile stability function for momentum
      !-------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pta
      !
      REAL(wp), DIMENSION(jpi,jpj)             :: psi_mc
      REAL(wp), DIMENSION(:,:), POINTER        :: X2, X, stabit
      !-------------------------------------------------------------------------------
      !
      CALL wrk_alloc( jpi,jpj, X2, X, stabit )
      !
      X2 = SQRT( ABS( 1. - 16.*pta ) )  ;  X2 = MAX( X2 , 1. )   ;   X = SQRT( X2 ) 
      stabit = 0.5 + SIGN( 0.5 , pta )
      psi_mc = -5.*pta*stabit  & ! Stable
         &    + (1. - stabit)*(2.*LOG((1. + X)*0.5) + LOG((1. + X2)*0.5) - 2.*ATAN(X) + rpi*0.5)  ! Unstable
      !
      CALL wrk_dealloc( jpi,jpj, X2, X, stabit )
      !
   END FUNCTION psi_mc



   FUNCTION psi_hc( pta )    !! Psis, L&Y 2004 eq. (8c), (8d), (8e)
      !-------------------------------------------------------------------------------
      ! universal profile stability function for temperature and humidity
      !-------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   pta
      !
      REAL(wp), DIMENSION(jpi,jpj)             ::   psi_hc
      REAL(wp), DIMENSION(:,:), POINTER        ::   X2, X, stabit
      !-------------------------------------------------------------------------------
      !
      CALL wrk_alloc( jpi,jpj, X2, X, stabit )
      !
      X2 = SQRT( ABS( 1. - 16.*pta ) )   ;   X2 = MAX( X2 , 1. )   ;   X = SQRT( X2 )
      stabit = 0.5 + SIGN( 0.5 , pta )
      psi_hc = -5.*pta*stabit   &                                       ! Stable
         &    + (1. - stabit)*(2.*LOG( (1. + X2)*0.5 ))                ! Unstable
      !
      CALL wrk_dealloc( jpi,jpj, X2, X, stabit )
      !
   END FUNCTION psi_hc


      REAL(wp) FUNCTION cd_HR(speed,delt)
      !!----------------------------------------------------------------------
      !! --- calculates the Drag Coefficient as a function of the abs. value
      !! --- of the wind velocity ( Hellermann and Rosenstein )
      !!----------------------------------------------------------------------

       REAL(wp), INTENT(in) :: speed,delt
       REAL(wp), PARAMETER  :: a1=0.934e-3 , a2=0.788e-4, a3=0.868e-4     
       REAL(wp), PARAMETER  :: a4=-0.616e-6, a5=-.120e-5, a6=-.214e-5

        cd_HR = a1 + a2*speed + a3*delt + a4*speed*speed        &
           + a5*delt*delt  + a6*speed*delt

      END FUNCTION cd_HR

      REAL(wp) function HEATLAT(t)
      !!----------------------------------------------------------------------
      !! --- calculates the Latent Heat of Vaporization ( J/kg ) as function of
      !! --- the temperature ( Celsius degrees )
      !! --- ( from A. Gill  pag. 607 )
      !!
      !! --- Constant Latent Heat of Vaporization ( Rosati,Miyakoda 1988 )
      !!     L = 2.501e+6  (MKS)
      !!----------------------------------------------------------------------

      REAL(wp) , intent(in) :: t

      heatlat = 2.5008e+6 -2.3e+3*t

      END FUNCTION HEATLAT


   SUBROUTINE qshort(hour,alat,alon,cldnow,qsw)
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE qshort  ***
      !!
      !! ** Purpose :   Compute Solar Radiation
      !!
      !! ** Method  :   Compute Solar Radiation according Astronomical
      !!                formulae
      !!
      !! References :   Reed RK (1975) and Reed RK (1977)
      !!
      !! Note: alat,alon - (lat, lon)  in radians
      !!----------------------------------------------------------------------
        REAL(wp), INTENT (in) :: hour

        REAL(wp), INTENT(in ), DIMENSION(:,:) :: alat,alon
        REAL(wp), INTENT(in ), DIMENSION(:,:) :: cldnow
        REAL(wp), INTENT(out), DIMENSION(:,:) :: qsw
        REAL(wp), DIMENSION(12) :: alpham

        REAL(wp), PARAMETER ::   eclips=23.439* (3.141592653589793_wp / 180._wp)
        REAL(wp), PARAMETER ::   solar = 1350.
        REAL(wp), PARAMETER ::   tau = 0.7
        REAL(wp), PARAMETER ::   aozone = 0.09
        REAL(wp), PARAMETER ::   yrdays = 360.
        REAL(wp) :: days, th0,th02,th03, sundec, thsun, coszen, qatten
        REAL(wp) :: qzer, qdir,qdiff,qtot,tjul,sunbet
        REAL(wp) :: albedo
        INTEGER :: jj, ji

      !!----------------------------------------------------------------------
      !! --- albedo monthly values from Payne (1972) as means of the values
      !! --- at 40N and 30N for the Atlantic Ocean ( hence the same latitudinal
      !! --- band of the Mediterranean Sea ) :
      !!----------------------------------------------------------------------

        data alpham /0.095,0.08,0.065,0.065,0.06,0.06,0.06,0.06,        &
                    0.065,0.075,0.09,0.10/

      !!----------------------------------------------------------------------
      !!   days is the number of days elapsed until the day=nday_year
      !!----------------------------------------------------------------------
        days = nday_year -1.
        th0  = 2.*rpi*days/yrdays
        th02 = 2.*th0
        th03 = 3.*th0

      !! --- sun declination :
      !!----------------------------------------------------------------------
        sundec = 0.006918 - 0.399912*cos(th0) + 0.070257*sin(th0) -   &
                          0.006758*cos(th02) + 0.000907*sin(th02) -   &
                          0.002697*cos(th03) + 0.001480*sin(th03)

      DO jj = 1, jpj
         DO ji = 1, jpi

      !! --- sun hour angle :
      !!----------------------------------------------------------------------
          thsun = (hour -12.)*15.*rad + alon(ji,jj)

      !! --- cosine of the solar zenith angle :
      !!----------------------------------------------------------------------
          coszen =sin(alat(ji,jj))*sin(sundec)                 &
                    +cos(alat(ji,jj))*cos(sundec)*cos(thsun)

          IF(coszen .LE. 5.035D-04) THEN
            coszen = 0.0
            qatten = 0.0
          ELSE
            qatten = tau**(1./coszen)
          END IF

          qzer  = coszen * solar *                                 &
                  (1.0+1.67E-2*cos(rpi*2.*(days-3.0)/365.0))**2
          qdir  = qzer * qatten
          qdiff = ((1.-aozone)*qzer - qdir) * 0.5
          qtot  =  qdir + qdiff
          tjul = (days -81.)*rad

      !! --- sin of the solar noon altitude in radians :
      !!----------------------------------------------------------------------
          sunbet=sin(alat(ji,jj))*sin(eclips*sin(tjul)) +   &
                 cos(alat(ji,jj))*cos(eclips*sin(tjul))

      !! --- solar noon altitude in degrees :
      !!----------------------------------------------------------------------

         sunbet = asin(sunbet)/rad

      !!----------------------------------------------------------------------
      !! --- calculates the albedo according to Payne (1972)
      !!----------------------------------------------------------------------

         albedo = alpham(nmonth)

      !!----------------------------------------------------------------------
      !! --- ( radiation as from Reed(1977), Simpson and Paulson(1979) )
      !! --- calculates SHORT WAVE FLUX ( watt/m*m )
      !! --- ( Rosati,Miyakoda 1988 ; eq. 3.8 )
      !!----------------------------------------------------------------------

          IF(cldnow(ji,jj).LT.0.3) THEN
             qsw(ji,jj) = qtot * (1.-albedo)
          ELSE
             qsw(ji,jj) = qtot*(1.-0.62*cldnow(ji,jj)              &
                                + .0019*sunbet)*(1.-albedo)
          ENDIF

         END DO
      END DO

   END SUBROUTINE qshort


   !!======================================================================

END MODULE sbcblk_mfs
