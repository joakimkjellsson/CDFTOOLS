PROGRAM cdfthermohaline
  !!======================================================================
  !!                     ***  PROGRAM  cdfmocsig  ***
  !!=====================================================================
  !!  ** Purpose : Compute the thermohaline streamfunction 
  !!               using temperature and salinity bins. 
  !!
  !!  ** Method  : The thermohaline streamfunction (THS) is computed from the 
  !!               U,V,W velocity field, collected in temperature and salinity
  !!               bins and integrated throughout the salinity bins. 
  !!               Masking for oceanic basins is available.
  !!               In the present version the masking corresponds to the global
  !!               configuration. Streamfunctions for Global, Atlantic, 
  !!               Indo-Pacific, Indian, and Pacific ocean. 
  !!               Results are saved on psi_thermohaline.nc file with 
  !!               variable names respectively 
  !!               zothsglo, zothsatl, zothsinp, zothsind, zothspac.
  !!               If no new_maskglo.nc file found, then the mask.nc file is used and
  !!               only zothsglo is computed.
  !!
  !!
  !! History : 1.0  : 08/2019  : J. Kjellsson  : Original code from cdfmocsig
  !!                                             and github.com/joakimkjellsson/hydrocode.git
  !!
  !!----------------------------------------------------------------------
  USE cdfio
  USE eos
  USE modcdfnames
  USE modutils
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class transport
  !!----------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER(KIND=2), DIMENSION (:,:,:), ALLOCATABLE ::  ibmask              ! nbasins x npiglo x npjglo
  INTEGER(KIND=2), DIMENSION (:,:,:), ALLOCATABLE ::  itmask              ! tmask from salinity field

  INTEGER(KIND=4)                                 :: jbasin, jj, jk       ! dummy loop index
  INTEGER(KIND=4)                                 :: ji, jt, jn, jbin     ! dummy loop index
  INTEGER(KIND=4)                                 :: im, ip, jm, jp       ! dummy loop index
  INTEGER(KIND=4)                                 :: ikc, ikp, ikn, ikt   ! dummy loop index
  INTEGER(KIND=4)                                 :: m1, mk, ir           ! dummy loop index
  INTEGER(KIND=4)                                 :: it                   ! time index for vvl
  INTEGER(KIND=4)                                 :: item, isal
  INTEGER(KIND=4)                                 :: nbins                ! number of  density  bins
  INTEGER(KIND=4)                                 :: npglo, npatl, npinp  ! basins index (mnemonics)
  INTEGER(KIND=4)                                 :: npind, nppac, npsoc  !  "      "
  INTEGER(KIND=4)                                 :: nbasins              ! number of basins
  INTEGER(KIND=4)                                 :: ierr                 ! working integer
  INTEGER(KIND=4)                                 :: narg, iargc, iarg    ! command line  browsing 
  INTEGER(KIND=4)                                 :: ijarg, ii            !  "             "
  INTEGER(KIND=4)                                 :: ib                   ! current bin number
  INTEGER(KIND=4), DIMENSION(0:6)                 :: mm1, mm2             ! bin indices for each 
                                                                          ! surrounding grid cell
  INTEGER(KIND=4)                                 :: ii1, ii2             ! current I index
  INTEGER(KIND=4)                                 :: ij1, ij2             ! current J index
  INTEGER(KIND=4)                                 :: npiglo,npjglo        ! size of the domain
  INTEGER(KIND=4)                                 :: npk, npt             ! size of the domain
  INTEGER(KIND=4)                                 :: ncout,ncout2         ! ncid of output file
  INTEGER(KIND=4)                                 :: nvaro                ! number of output variables
  INTEGER(KIND=4), DIMENSION(2)                   :: iloc                 ! working array
  INTEGER(KIND=4), DIMENSION(:),      ALLOCATABLE :: ipk, id_varout       ! output variable levels and id
  INTEGER(KIND=4), DIMENSION(:),      ALLOCATABLE :: ipr, id_varout2           ! output variable levels and id
  INTEGER(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: ibin                 ! remaping density in bin number
  INTEGER(KIND=4), DIMENSION(:,:,:),  ALLOCATABLE :: zm1,zm2              ! remaping density in bin number
  INTEGER(KIND=4)                                 :: ithread, nthreads    ! OMP thread identifiers
  
  REAL(KIND=4), PARAMETER                         :: rp_spval=99999.      !
  REAL(KIND=4)                                    :: pref=0.              ! depth reference for pot. density 
  REAL(KIND=4)                                    :: sigmin               ! minimum density for bining
  REAL(KIND=4)                                    :: sigstp               ! density step for bining
  REAL(KIND=4)                                    :: temstp               ! temperature step for bining
  REAL(KIND=4)                                    :: salstp               ! salinity step for bining
  REAL(KIND=4)                                    :: temmin               ! minimum temperature for binning
  REAL(KIND=4)                                    :: salmin               ! minimum salinity for binning
  REAL(KIND=4)                                    :: temmax               ! maximum temperature for binning
  REAL(KIND=4)                                    :: salmax               ! maximum salinity for binning
  REAL(KIND=4)                                    :: zsps                 ! Salinity Missing value
  REAL(KIND=4)                                    :: zspt                 ! Temperature Missing value
  REAL(KIND=4)                                    :: zspu                 ! Zonal Vel.  Missing value
  REAL(KIND=4)                                    :: zspv                 ! Merid. Vel.  Missing value
  REAL(KIND=4)                                    :: zspw                 ! Vert. Vel.  Missing value
  REAL(KIND=4), DIMENSION (:),        ALLOCATABLE :: sigma                ! density coordinate 
  REAL(KIND=4), DIMENSION (:),        ALLOCATABLE :: temperature          ! temperature coordinate 
  REAL(KIND=4), DIMENSION (:),        ALLOCATABLE :: salinity             ! salinity coordinate 
  REAL(KIND=4), DIMENSION (:),        ALLOCATABLE :: e31d                 ! vertical level (full step)
  REAL(KIND=4), DIMENSION (:),        ALLOCATABLE :: gdep                 ! depth of T layers ( full step)
  REAL(KIND=4), DIMENSION (:,:),      ALLOCATABLE :: e1t, e2t             ! horizontal metrics, latitude
  REAL(KIND=4), DIMENSION (:,:),      ALLOCATABLE :: e1v, e2u, gphiv      ! horizontal metrics, latitude
  REAL(KIND=4), DIMENSION (:,:,:),    ALLOCATABLE :: zt, zs               ! temperature, salinity
  REAL(KIND=4), DIMENSION (:,:,:),    ALLOCATABLE :: zu, zv, zw           ! velocities 
  REAL(KIND=4), DIMENSION (:,:),      ALLOCATABLE :: zueiv, zveiv         ! bolus velocities
  REAL(KIND=4), DIMENSION (:,:,:),    ALLOCATABLE :: e3u,e3v              ! vertical metrics
  REAL(KIND=4), DIMENSION (:,:),      ALLOCATABLE :: rdumlon              ! dummy longitude = 0.
  REAL(KIND=4), DIMENSION (:,:),      ALLOCATABLE :: rdumlat              ! latitude for i = north pole
  REAL(KIND=4), DIMENSION (:,:,:),    ALLOCATABLE :: zutmp                ! temporary U array (volume flux)
  REAL(KIND=4), DIMENSION (:,:,:),    ALLOCATABLE :: zvtmp                ! temporary V array (volume flux)
  REAL(KIND=4), DIMENSION (:,:,:),    ALLOCATABLE :: zwtmp                ! temporary W array (volume flux)
  REAL(KIND=4), DIMENSION (:,:,:),    ALLOCATABLE :: zttmp                ! temporary T array
  REAL(KIND=4), DIMENSION (:,:,:),    ALLOCATABLE :: zstmp                ! temporary S array
  REAL(KIND=4), DIMENSION (:,:),      ALLOCATABLE :: zarea                ! product e1v * e3v
  REAL(KIND=4), DIMENSION (:),        ALLOCATABLE :: zflux                ! temporary flux array

  REAL(KIND=8), DIMENSION (:),        ALLOCATABLE :: dtim                 ! time counter
  REAL(KIND=8), DIMENSION(:,:,:),     ALLOCATABLE :: dens                 ! density
  REAL(KIND=8), DIMENSION(:,:,:),     ALLOCATABLE :: psirr_tmp            ! temporary transport array
  REAL(KIND=8), DIMENSION(:,:,:),     ALLOCATABLE :: volrr_tmp            ! temporary cumulated volume array
  REAL(KIND=8), DIMENSION(:,:,:),     ALLOCATABLE :: psiyz                ! nbasins x npjglo x npk
  REAL(KIND=8), DIMENSION(:,:,:),     ALLOCATABLE :: psirr                ! nbasins x npjglo x nbins
  REAL(KIND=8), DIMENSION(:,:,:),     ALLOCATABLE :: volrr                ! cumulated volume 

  CHARACTER(LEN=256)                              :: cf_ufil              ! zonal velocity file
  CHARACTER(LEN=256)                              :: cf_vfil              ! meridional velocity file
  CHARACTER(LEN=256)                              :: cf_wfil              ! vertical velocity file
  CHARACTER(LEN=256)                              :: cf_tfil              ! temperature/salinity file
  CHARACTER(LEN=256)                              :: cf_sfil              ! salinity file (option)
  CHARACTER(LEN=256)                              :: cf_msf='psimsf.nc'   ! output file
  CHARACTER(LEN=256)                              :: cf_ths='psiths.nc'   ! output file
  CHARACTER(LEN=255)                              :: cglobal              ! Global attribute
  CHARACTER(LEN=256)                              :: cldum                ! dummy char variable

  TYPE(variable), DIMENSION(:), ALLOCATABLE       :: styrvar              ! output var properties
  TYPE(variable), DIMENSION(:), ALLOCATABLE       :: strrvar              ! output var properties

  LOGICAL, DIMENSION(3)                           :: lbin                 ! flag for bin specifications
  LOGICAL                                         :: lntr                 ! flag for neutral density
  LOGICAL                                         :: lbas   = .FALSE.     ! flag for basins file
  LOGICAL                                         :: lisodep= .FALSE.     ! flag for isopycnal zonal mean
  LOGICAL                                         :: lprint = .FALSE.     ! flag for extra print
  LOGICAL                                         :: leiv   = .FALSE.     ! flag for Eddy Induced Velocity (GM)
  LOGICAL                                         :: lfull  = .FALSE.     ! flag for full step
  LOGICAL                                         :: lchk   = .FALSE.     ! flag for missing file
  LOGICAL                                         :: lreadW = .FALSE.     ! flag for reading w velocities
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfthermohaline  -u U-file -v V-file -t T-file -s S-file | -ntr [-eiv] [-full] ...'
     PRINT *,'        ... [-sigmin sigmin] [-sigstp sigstp] [-nbins nbins] [-isodep] ...'
     PRINT *,'        ... [-w W-file ] [-o OUT-file] [-vvl] [-verbose]'
     PRINT *,'      '
     PRINT *,'     PURPOSE : '
     PRINT *,'       Compute the MOC in density-latitude coordinates. The global value is '
     PRINT *,'       always computed. Values for oceanic sub-basins are calculated if the '
     PRINT *,'       ', TRIM(cn_fbasins), ' file is provided.'
     PRINT *,'      '
     PRINT *,'       The reference depth for potential density is given with ''-D'' option.'
     PRINT *,'       Density ranges and number of bins to use are pre-defined only for three'
     PRINT *,'       reference depth (0, 1000 and 2000 m). For other reference depth, the '
     PRINT *,'       density binning must be specified using the relevant options for setting'
     PRINT *,'       the minimum density, the density step and the number of bins to use.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'        -v V-file  : Netcdf gridV file.'
     PRINT *,'        -t T-file  : Netcdf gridT file with temperature and salinity.'
     PRINT *,'             If salinity not in T-file use -s option.'
     PRINT *,'        -r ref-depth : reference depth for density. '
     PRINT *,'            For depth values of 0 1000 or 2000 m, pre-defined limits for '
     PRINT *,'            minimum density, number of density bins and width of density '
     PRINT *,'            bins are provided. For other reference depth, you must use the'
     PRINT *,'            options ''-sigmin'', ''-sigstp'' and ''-nbins'' (see below).'
     PRINT *,'        or '
     PRINT *,'        -ntr : uses neutral density (no default bin defined so far), no ''-r'''
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-s S-file   ] : Specify salinity file if not T-file.'
     PRINT *,'       [-eiv ] : takes into account VEIV Meridional eddy induced velocity.'
     PRINT *,'                 -> To be used only if Gent and McWilliams parameterization '
     PRINT *,'                    has been used. '
     PRINT *,'       [-full       ] : Works with full step instead of standard partial steps.'
     PRINT *,'       [-sigmin     ] : Specify minimum of density for bining.'
     PRINT *,'       [-sigstp     ] : Specify density step for bining.'
     PRINT *,'       [-nbins      ] : Specify the number of density bins you want.'
     PRINT *,'       [-isodep     ] : Compute the zonal mean of isopycnal depths used for '
     PRINT *,'                        mocsig.'
     PRINT *,'       [-o OUT-file ] : Specify output file name instead of ', TRIM(cf_ths)
     PRINT *,'       [-vvl        ] : Use time-varying vertical metrics.'
     PRINT *,'       [-verbose    ] : Verbose option for more info during execution.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        Files ', TRIM(cn_fzgr),', ',TRIM(cn_fhgr),', ', TRIM(cn_fmsk)
     PRINT *,'        File ', TRIM(cn_fbasins),' is optional [sub basins masks]'
     PRINT *,'      '
     PRINT *,'     OPENMP SUPPORT : yes '
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_ths) 
     PRINT *,'       variables ',TRIM( cn_zothsglo),' : Global ocean '
     PRINT *,'       variables ',TRIM( cn_zothsatl),' : Atlantic Ocean '
     PRINT *,'       variables ',TRIM( cn_zothsinp),' : Indo Pacific '
     PRINT *,'       variables ',TRIM( cn_zothsind),' : Indian Ocean alone'
     PRINT *,'       variables ',TRIM( cn_zothspac),' : Pacific Ocean alone'
     PRINT *,'       If file ',TRIM(cn_fbasins),' is not present, ',TRIM(cn_fmsk),' file is used and'
     PRINT *,'       only ',TRIM( cn_zothsglo),' is produced.'
     PRINT *,'       If option -isodep is used, each MOC variable is complemented by a iso'
     PRINT *,'       variable, giving the zonal mean of ispycnal depth (e.g.',TRIM(cn_zoisoglo),').'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'       cdfmoc '
     PRINT *,'      '
     STOP 
  ENDIF

  cglobal = 'Partial step computation'
  lbin=(/.TRUE.,.TRUE.,.TRUE./)
  ijarg = 1 ; ii = 0 ! ii is used to count mandatory arguments
  cf_sfil = 'none'
  DO WHILE ( ijarg <= narg )
     CALL getarg (ijarg, cldum) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-u'     ) ;  CALL getarg (ijarg, cf_ufil) ; ijarg=ijarg+1 ; ii=ii+1
     CASE ( '-v'     ) ;  CALL getarg (ijarg, cf_vfil) ; ijarg=ijarg+1 ; ii=ii+1
     CASE ( '-t'     ) ;  CALL getarg (ijarg, cf_tfil) ; ijarg=ijarg+1 ; ii=ii+1
     CASE ( '-s'     ) ;  CALL getarg (ijarg, cf_sfil) ; ijarg=ijarg+1 ; ii=ii+1
     CASE ( '-r'     ) ;  CALL getarg (ijarg, cldum  ) ; ijarg=ijarg+1 ; READ(cldum,*) pref
     CASE ( '-ntr'   ) ;  lntr = .TRUE. ; ii=ii+1
        ! options
     CASE ('-full'   ) ; lfull   = .TRUE. ; cglobal = 'Full step computation'
     CASE ('-eiv'    ) ; leiv    = .TRUE.
     CASE ('-sigmin' ) ; CALL getarg (ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) sigmin ; lbin(1) = .FALSE.
     CASE ('-nbins'  ) ; CALL getarg (ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) nbins  ; lbin(2) = .FALSE.
     CASE ('-sigstp' ) ; CALL getarg (ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) sigstp ; lbin(3) = .FALSE.
     CASE ('-w'      ) ; CALL getarg (ijarg, cf_wfil) ; ijarg=ijarg+1 ; lreadW = .TRUE.
     CASE ('-o'      ) ; CALL getarg (ijarg, cf_ths) ; ijarg=ijarg+1 
     CASE ('-vvl'    ) ; lg_vvl  = .TRUE.
     CASE ('-isodep' ) ; lisodep = .TRUE.
     CASE ('-verbose') ; lprint  = .TRUE.
     CASE DEFAULT     ; PRINT *,' ERROR : ',TRIM(cldum), ' : unknown option.'  ; STOP 99
     END SELECT
  END DO

  IF ( ii /= 4  ) THEN ; PRINT *,' ERROR : mandatory arguments missing, see usage please !'  ; STOP 99
  ENDIF

  IF ( cf_sfil == 'none' ) cf_sfil=cf_tfil

  ! check file existence
  lchk = lchk .OR. chkfile ( cn_fhgr )
  lchk = lchk .OR. chkfile ( cn_fzgr )
  lchk = lchk .OR. chkfile ( cn_fmsk )
  lchk = lchk .OR. chkfile ( cf_vfil )
  lchk = lchk .OR. chkfile ( cf_tfil )
  lchk = lchk .OR. chkfile ( cf_sfil )
  IF ( lchk ) STOP 99  ! missing file(s)

  ! Look for salinity spval
  zsps = getspval(cf_sfil, cn_vosaline)
  zspt = getspval(cf_tfil, cn_votemper)
  zspu = getspval(cf_vfil, cn_vozocrtx)
  zspv = getspval(cf_vfil, cn_vomecrty)
  IF (lreadW) THEN
     zspw = getspval(cf_vfil, cn_vovecrtz)
  ELSE
     zspw = 0.
  END IF
  
  IF ( lg_vvl )  THEN
     cn_fe3v = cf_vfil
     cn_ve3v = cn_ve3vvvl
  ENDIF

  ! re-use lchk for binning control : TRUE if no particular binning specified
  lchk = lbin(1) .OR. lbin(2) .OR. lbin(3) 

  npiglo = getdim (cf_vfil,cn_x)
  npjglo = getdim (cf_vfil,cn_y)
  npk    = getdim (cf_vfil,cn_z)
  npt    = getdim (cf_vfil,cn_t)

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  !setting up the building command in global attribute
  CALL SetGlobalAtt(cglobal, 'A')  ! append command name to global attribute

  !  Detects newmaskglo file
  lbas = .NOT. chkfile (cn_fbasins )

  IF (lbas) THEN ; nbasins = 5
  ELSE           ; nbasins = 1
  ENDIF

  IF ( lisodep ) THEN ; nvaro = 2 * nbasins
  ELSE                ; nvaro =     nbasins
  ENDIF

  ALLOCATE ( styrvar(nvaro), strrvar(nvaro), ipk(nvaro), ipr(nvaro), id_varout(nvaro), id_varout2(nvaro) )

  IF ( lchk )  THEN  ! use default bins definition according to pref 
     ! Define parameters
     nbins  = 151
     
     temmin = -2.
     temmax = 38.
     salmin = 26.
     salmax = 40.
     
     ! define temperature and salinity steps 
     temstp = (temmax-temmin) / float(nbins-1)
     salstp = (salmax-salmin) / float(nbins-1)
     
  END IF   
  
  PRINT '(a,f5.2,a,f5.2,a,i3)', '  You are using -temmin ', temmin,' -temstp ', temstp,' -nbins ', nbins
  PRINT '(a,f5.2,a,f5.2,a,i3)', '  You are using -salmin ', salmin,' -salstp ', salstp,' -nbins ', nbins

  ALLOCATE ( sigma(nbins), salinity(nbins), temperature(nbins) )  

  ! define densities at middle of bins
  DO ji=1,nbins
     !sigma(ji)  = sigmin +(ji-0.5)*sigstp
     temperature(ji) = temmin + (ji-0.5)*temstp
     salinity(ji) = salmin + (ji-0.5)*salstp 
  ENDDO
  
  IF (lprint) PRINT *, ' min temperature: ', temmin, 'max temperature: ',temmax
  IF (lprint) PRINT *, ' min salinity: ', salmin, 'max salinity: ',salmax
  
  ! Allocate arrays
  ALLOCATE ( ibmask(nbasins,npiglo,npjglo) )
  ALLOCATE ( zu (npiglo,npjglo,3), zv (npiglo,npjglo,3), zw (npiglo, npjglo,3) )
  ALLOCATE ( zt (npiglo,npjglo,3), zs (npiglo,npjglo,3), zarea(npiglo, npjglo) )
  ALLOCATE ( e3u(npiglo,npjglo,3), e3v(npiglo,npjglo,3) )
  ALLOCATE ( ibin(npiglo, npjglo) )
  ALLOCATE ( zm1(npiglo, npjglo,0:6), zm2(npiglo, npjglo,0:6) )
  ALLOCATE ( e1t(npiglo, npjglo), e2t(npiglo, npjglo) )
  ALLOCATE ( e1v(npiglo,npjglo), e2u(npiglo,npjglo), gphiv(npiglo,npjglo) )
  !ALLOCATE ( dmoc(nvaro, nbins, npjglo ) )
  ALLOCATE ( rdumlon(1,npjglo) , rdumlat(1,npjglo))
  ALLOCATE ( dens(npiglo,npjglo,3) )
  ALLOCATE ( zutmp(npiglo,npjglo,3), zvtmp(npiglo,npjglo,3), zwtmp(npiglo,npjglo,3) )
  ALLOCATE ( itmask(npiglo,npjglo,3), zttmp(npiglo,npjglo,3), zstmp(npiglo,npjglo,3) )
  ALLOCATE ( dtim(npt), e31d(npk)  )
  ALLOCATE ( zflux(0:6) )
  
  !ALLOCATE ( psiyz(nvaro, npjglo, npk) )
  !ALLOCATE ( psirr(nvaro, nbins, nbins), volrr(nbasins,nbins,nbins) )
  ALLOCATE ( psiyz(npjglo, npk, nvaro) )
  ALLOCATE ( psirr(nbins, nbins, nvaro), volrr(nbasins,nbins,nbins) )
  ALLOCATE ( psirr_tmp(nbasins,nbins,nbins), volrr_tmp(nbasins,nbins,nbins) )
  
  !IF ( lisodep) THEN 
  !   ALLOCATE ( depi(nvaro, nbins, npjglo), gdep(npk))
  !   ALLOCATE ( wdep(nvaro, nbins, npjglo)           )
  !ENDIF
  IF ( leiv   ) ALLOCATE ( zueiv (npiglo,npjglo), zveiv (npiglo,npjglo))

  IF (lprint) PRINT*,' Get horizontal mesh '
  e1v(:,:)   = getvar(cn_fhgr,   cn_ve1v,  1, npiglo, npjglo) 
  e2u(:,:)   = getvar(cn_fhgr,   cn_ve2u,  1, npiglo, npjglo) 
  e1t(:,:)   = getvar(cn_fhgr,   cn_ve1t,  1, npiglo, npjglo) 
  e2t(:,:)   = getvar(cn_fhgr,   cn_ve2t,  1, npiglo, npjglo) 
  
  IF ( lfull  ) e31d(:) = getvare3(cn_fzgr, cn_ve3t1d,  npk )
  IF ( lisodep) gdep(:) = -getvare3(cn_fzgr, cn_gdept, npk )  ! take negative value
  ! to be compliant with zonal mean

  IF ( npjglo > 1 ) THEN 
     gphiv(:,:)   = getvar(cn_fhgr, cn_gphiv, 1, npiglo, npjglo)
     iloc         = MAXLOC(gphiv)
     rdumlat(1,:) = gphiv(iloc(1),:)
  ELSE
     rdumlat(1,:) = 0.
  ENDIF
  rdumlon(:,:) = 0.               ! set the dummy longitude to 0

  ! create output fileset
  !global   ; Atlantic  ; Indo-Pacif ; Indian  ; Pacif    ; Southern Ocean
  npglo= 1  ; npatl=2   ;  npinp=3   ; npind=4 ; nppac=5  ; npsoc=6

  IF (lprint) PRINT*,' Create output file'
  CALL CreateOutputFile

  ! reading the masks
  ibmask(npglo,:,:) = getvar(cn_fmsk, cn_vmask, 1, npiglo, npjglo)

  IF ( lbas ) THEN
     ibmask(npatl,:,:) = getvar(cn_fbasins, cn_tmaskatl, 1, npiglo, npjglo)
     ibmask(npind,:,:) = getvar(cn_fbasins, cn_tmaskind, 1, npiglo, npjglo)
     ibmask(nppac,:,:) = getvar(cn_fbasins, cn_tmaskpac, 1, npiglo, npjglo)
     ibmask(npinp,:,:) = ibmask(nppac,:,:) + ibmask(npind,:,:)
     ! ensure that there are no overlapping on the masks
     WHERE(ibmask(npinp,:,:) > 0 ) ibmask(npinp,:,:) = 1
     ! change global mask for GLOBAL periodic condition
     ibmask(1,1,     :) = 0.
     ibmask(1,npiglo,:) = 0.
  ENDIF

  timeLoop: DO jt=1, npt
     
     IF (lprint) PRINT *,' Reading time step: ',jt
     
     IF ( lg_vvl ) THEN 
        it=jt
     ELSE              
        it=1
     ENDIF
     
     ! initialize streamfunction to 0
     psiyz(:,:,:) = 0.d0 
     psirr(:,:,:) = 0.d0 
     !IF ( lisodep ) THEN ; depi(:,:,:) = 0.d0 ; wdep(:,:,:) = 0.d0
     !ENDIF
     
     ! We always keep three levels in memory
     ! ikc is the center k-index, i.e. the one we are working on 
     ! ikp is the previous (below)
     ! ikn is the next (above)
     !
     ! This is because we need to know the temperature and salinity above and below
     ! when calculating the fluxes.
     ! But for global 1/12 or higher resolution, this becomes really memory-consuming
     ! so its better to keep as little as possible in memory at any time. 
     ! At the end of the levelLoop, we then permute so that
     ! tmp = ikp 
     ! ikc = ikn
     ! ikp = ikc
     ! ikn = ikp 
     ikp = 1
     ikc = 2
     ikn = 3
     
     zutmp(:,:,:) = 0.
     zvtmp(:,:,:) = 0.
     zwtmp(:,:,:) = 0.
     
     ! loop from bottom to top 
     ! so we can assume no vertical flow at the bottom
     levelLoop: DO jk=npk-1,1,-1 ! level loop 
        
        !               for testing purposes only loop from 2 to 400
        
        IF (lprint) PRINT*,' Reading depth level ',jk-1,' into ikn '
        
        IF (jk > 1) THEN
           
           zu(:,:,ikn) = 0.
           zv(:,:,ikn) = 0.
           zw(:,:,ikn) = 0.
           
           ! Read the level above
           ! We read the current level last step in the loop, 
           ! This way we always keep 3 levels in memory
           zu(:,:,ikn) = getvar(cf_ufil, cn_vozocrtx, jk-1, npiglo, npjglo, ktime = jt)
           zv(:,:,ikn) = getvar(cf_vfil, cn_vomecrty, jk-1, npiglo, npjglo, ktime = jt)
           ! We should add an option to get w from u, v
           IF (lreadW) THEN
              zw(:,:,ikn) = getvar(cf_wfil, cn_vovecrtz, jk-1, npiglo, npjglo, ktime = jt)
           END IF
           
           ! Apply missing values
           WHERE( zu(:,:,ikn) == zspu ) zu(:,:,ikn) = 0.
           WHERE( zv(:,:,ikn) == zspv ) zv(:,:,ikn) = 0.
           WHERE( zw(:,:,ikn) == zspw ) zw(:,:,ikn) = 0.
           ! Add eddy-induced velocities
           IF ( leiv ) THEN
              IF (lprint) PRINT*,' Adding EIV '
              zueiv(:,:) = getvar(cf_ufil, cn_vozoeivu, jk-1, npiglo,npjglo, ktime = jt)
              zveiv(:,:) = getvar(cf_vfil, cn_vomeeivv, jk-1, npiglo,npjglo, ktime = jt)
              
              zu(:,:,ikn)    = zu(:,:,ikn) + zueiv(:,:)
              zv(:,:,ikn)    = zv(:,:,ikn) + zveiv(:,:)
           END IF
        
           ! Get T, S and apply missing values
           zt(:,:,ikn) = getvar(cf_tfil, cn_votemper, jk-1, npiglo, npjglo, ktime = jt)
           zs(:,:,ikn) = getvar(cf_sfil, cn_vosaline, jk-1, npiglo, npjglo, ktime = jt)
           WHERE( zt(:,:,ikn) == zspt ) zt(:,:,ikn) = 0.
           WHERE( zs(:,:,ikn) == zsps ) zs(:,:,ikn) = 0.
           
           IF ( lfull ) THEN 
              e3u(:,:,ikn) = e31d(jk-1)
              e3v(:,:,ikn) = e31d(jk-1)
           ELSE              
              e3u(:,:,ikn) = getvar(cn_fe3u, cn_ve3u, jk-1, npiglo, npjglo, ktime=it, ldiom=.NOT.lg_vvl )
              e3v(:,:,ikn) = getvar(cn_fe3v, cn_ve3v, jk-1, npiglo, npjglo, ktime=it, ldiom=.NOT.lg_vvl )
           ENDIF
           
        ELSE
           
           zu(:,:,ikn) = 0.
           zv(:,:,ikn) = 0.
           zw(:,:,ikn) = 0.
        
        END IF
        
        IF (jk == npk-1) THEN
           
           IF (lprint) PRINT*,' Bottom: reading depth level ',jk,' into ikc '
           
           ! For the first step, also read the current level
           zu(:,:,ikc) = getvar(cf_ufil, cn_vozocrtx, jk, npiglo, npjglo, ktime = jt)
           zv(:,:,ikc) = getvar(cf_vfil, cn_vomecrty, jk, npiglo, npjglo, ktime = jt)
           zu(:,:,ikp) = 0.
           zv(:,:,ikp) = 0.
           IF (lreadW) THEN
              zw(:,:,ikc) = getvar(cf_wfil, cn_vovecrtz, jk, npiglo, npjglo, ktime = jt)
           END IF
           ! and assume no flow below
           zw(:,:,ikp) = 0.
           
           WHERE( zu(:,:,ikc) == zspu ) zu(:,:,ikc) = 0.
           WHERE( zv(:,:,ikc) == zspv ) zv(:,:,ikc) = 0.
           WHERE( zw(:,:,ikc) == zspw ) zw(:,:,ikc) = 0.
           
           IF ( leiv ) THEN
              zueiv(:,:) = getvar(cf_ufil, cn_vozoeivu, jk, npiglo,npjglo, ktime = jt)
              zveiv(:,:) = getvar(cf_vfil, cn_vomeeivv, jk, npiglo,npjglo, ktime = jt)
              
              zu(:,:,ikc)    = zu(:,:,ikc) + zueiv(:,:)
              zv(:,:,ikc)    = zv(:,:,ikc) + zveiv(:,:)
           END IF
           
           zt(:,:,ikc) = getvar(cf_tfil, cn_votemper, jk, npiglo, npjglo, ktime = jt)
           zs(:,:,ikc) = getvar(cf_sfil, cn_vosaline, jk, npiglo, npjglo, ktime = jt)
           WHERE( zt(:,:,ikc) == zspt ) zt(:,:,ikc) = 0.
           WHERE( zs(:,:,ikc) == zsps ) zs(:,:,ikc) = 0.
           
           IF ( lfull ) THEN 
              e3u(:,:,ikc) = e31d(jk)
              e3v(:,:,ikc) = e31d(jk)
           ELSE              
              e3u(:,:,ikc) = getvar(cn_fe3u, cn_ve3u, jk, npiglo, npjglo, ktime=it, ldiom=.NOT.lg_vvl )
              e3v(:,:,ikc) = getvar(cn_fe3v, cn_ve3v, jk, npiglo, npjglo, ktime=it, ldiom=.NOT.lg_vvl )
           ENDIF
           
        END IF
        
        ! Apply missing values
        !WHERE( zu == zspu ) zu = 0.
        !WHERE( zv == zspv ) zv = 0.
        !WHERE( zw == zspw ) zw = 0.
        
        ! Get cell face areas 
        ! For u: e2u * e3u
        ! For v: e1v * e3v 
        ! For w: e1t * e2t
        
        
        !zarea(:,:) = e1v(:,:) * e3v(:,:)
        
        !
        ! Calculate volume fluxes
        ! 
        zutmp(:,:,ikn) = zu(:,:,ikn) * e2u(:,:) * e3u(:,:,ikn)
        zvtmp(:,:,ikn) = zv(:,:,ikn) * e1v(:,:) * e3v(:,:,ikn)
        IF (lreadW) THEN
           zwtmp(:,:,ikn) = zw(:,:,ikn) * e1t(:,:) * e2t(:,:) 
        ELSE
           PRINT*,' Calculate vertical flux ' 
           zwtmp(2:npiglo-1,2:npjglo-1,ikc) = zwtmp(2:npiglo-1,2:npjglo-1,ikp) - &
                                            & ( zutmp(2:npiglo-1,2:npjglo-1,ikc) - zutmp(1:npiglo-2,2:npjglo-1,ikc) + &
                                            &   zvtmp(2:npiglo-1,2:npjglo-1,ikc) - zvtmp(2:npiglo-1,1:npjglo-2,ikc) )
           
           IF (lg_vvl) THEN
              PRINT*,' WARNING: Calculating vertical flux with vvl not yet implemented' 
           END IF
           
           ! Cyclic boundaries (is this ever used?)                                 
           ! I dont think we want to use i=1 and npiglo since we dont want to double-count 
           zwtmp(1,     2:npjglo-1,ikc) = zwtmp(npiglo-1,2:npjglo-1,ikc)
           zwtmp(npiglo,2:npjglo-1,ikc) = zwtmp(2,       2:npjglo-1,ikc)
           
           ! Assume Antarctica for j=1 and j=npjglo same as j=npjglo-1
           ! Again, Im not sure we ever use these anyway
           zwtmp(2:npiglo-1,1,ikc) = 0.
           zwtmp(2:npiglo-1,npjglo,ikc) = zwtmp(2:npiglo-1,npjglo-1,ikc)
                                            
        END IF
        
        !
        !  finds density 
        itmask =  1
        WHERE ( zs == zsps ) itmask = 0
        IF ( lntr ) THEN 
           dens(:,:,ikn)  = sigmantr(zt(:,:,ikn), zs(:,:,ikn),       npiglo, npjglo)
        ELSE             
           dens(:,:,ikn)  = sigmai  (zt(:,:,ikn), zs(:,:,ikn), pref, npiglo, npjglo)
        ENDIF
        
        ! Apply mask to T, S, rho
        !zrtmp = dens* itmask 
        zttmp = zt  * itmask
        zstmp = zs  * itmask
        
        ! Find bin numbers for T, S
        !ibin(:,:) = INT( (zttmp-sigmin)/sigstp )
        !ibin(:,:) = MAX( ibin(:,:), 1    )
        !ibin(:,:) = MIN( ibin(:,:), nbins)
        
        IF ( npjglo > 1 ) THEN 
           ij1 = 2 
           ij2 = npjglo-1
        ELSE                  
           ij1 = 1 
           ij2 = 1        ! input file has only one j ( case of extracted broken lines) 
        ENDIF
        
        IF ( npiglo > 1 ) THEN 
           ii1 = 2
           ii2 = npiglo-1
        ELSE
           ii1 = 1
           ii2 = 1
        END IF
        
        zm1(:,:,:) = 0
        zm2(:,:,:) = 0
        
        zm1(ii1:ii2,ij1:ij2,0) = NINT( (zstmp(ii1:ii2    ,ij1:ij2    ,ikc)-salmin)/salstp ) + 1 
        zm1(ii1:ii2,ij1:ij2,1) = NINT( (zstmp(ii1+1:ii2+1,ij1:ij2    ,ikc)-salmin)/salstp ) + 1 
        zm1(ii1:ii2,ij1:ij2,2) = NINT( (zstmp(ii1-1:ii2-1,ij1:ij2    ,ikc)-salmin)/salstp ) + 1 
        zm1(ii1:ii2,ij1:ij2,3) = NINT( (zstmp(ii1:ii2    ,ij1+1:ij2+1,ikc)-salmin)/salstp ) + 1 
        zm1(ii1:ii2,ij1:ij2,4) = NINT( (zstmp(ii1:ii2    ,ij1-1:ij2-1,ikc)-salmin)/salstp ) + 1 
        zm1(ii1:ii2,ij1:ij2,5) = NINT( (zstmp(ii1:ii2    ,ij1:ij2    ,ikn)-salmin)/salstp ) + 1 
        zm1(ii1:ii2,ij1:ij2,6) = NINT( (zstmp(ii1:ii2    ,ij1:ij2    ,ikp)-salmin)/salstp ) + 1
        
        zm2(ii1:ii2,ij1:ij2,0) = NINT( (zttmp(ii1:ii2    ,ij1:ij2    ,ikc)-temmin)/temstp ) + 1 
        zm2(ii1:ii2,ij1:ij2,1) = NINT( (zttmp(ii1+1:ii2+1,ij1:ij2    ,ikc)-temmin)/temstp ) + 1 
        zm2(ii1:ii2,ij1:ij2,2) = NINT( (zttmp(ii1-1:ii2-1,ij1:ij2    ,ikc)-temmin)/temstp ) + 1 
        zm2(ii1:ii2,ij1:ij2,3) = NINT( (zttmp(ii1:ii2    ,ij1+1:ij2+1,ikc)-temmin)/temstp ) + 1 
        zm2(ii1:ii2,ij1:ij2,4) = NINT( (zttmp(ii1:ii2    ,ij1-1:ij2-1,ikc)-temmin)/temstp ) + 1 
        zm2(ii1:ii2,ij1:ij2,5) = NINT( (zttmp(ii1:ii2    ,ij1:ij2    ,ikn)-temmin)/temstp ) + 1 
        zm2(ii1:ii2,ij1:ij2,6) = NINT( (zttmp(ii1:ii2    ,ij1:ij2    ,ikp)-temmin)/temstp ) + 1
        
        zm1(:,:,:) = MAX(zm1(:,:,:), 1)
        zm1(:,:,:) = MIN(zm1(:,:,:), nbins)
        zm2(:,:,:) = MAX(zm2(:,:,:), 1)
        zm2(:,:,:) = MIN(zm2(:,:,:), nbins)
        
        !IF ( lisodep ) ALLOCATE ( depi_tmp(nbins,npiglo) )
        !IF ( lisodep ) ALLOCATE ( wdep_tmp(nbins,npiglo) )
        
        !$OMP PARALLEL  PRIVATE(psirr_tmp, volrr_tmp, ji, im, ip, jj, jp, jm, jn, ir, mm1, mm2, mk, m1, zflux)
        IF (lprint) PRINT*,' Entering main basin loop '
        basinLoop: DO jbasin = 1,nbasins
        
        IF (lprint) PRINT*,' Entering main J loop '
        !$OMP DO REDUCTION(+:psirr,volrr,psiyz) SCHEDULE(RUNTIME)
        latLoop: DO jj= ij1, ij2
           
           jp = jj+1
           jm = jj-1
           ! Not sure how to deal with these boundary conditions
           ! I don't think we need to worry about North fold
           IF (jm == 0) THEN
              jm = 1
           END IF
           IF (jp == ij2+1) THEN
              jp = ij2
           END IF
           
           !PRINT*,' Entering main I loop '
           lonLoop: DO ji= ii1, ii2 
              
              im = ji-1
              ip = ji+1
              ! We always work from 2 to npiglo-1
              ! so I don't think these are needed
              IF (im == 0) THEN
                 im = ii2
              END IF
              IF (ip == ii2+1) THEN
                 ip = ii1
              END IF
               
              !dmoc_tmp = 0.d0 
              psirr_tmp = 0.d0 
              
              !
              ! To calculate the thermohaline streamfunction we need 
              ! the bin indices for the T-point and the neighbouring 
              ! T-points. 
              ! m(0) - i,j
              ! m(1) - i+1,j
              ! m(2) - i-1,j
              ! m(3) - i,j+1
              ! m(4) - i,j-1
              ! m(5) - i,j,k+1 
              ! m(6) - i,j,k-1
              !
              
              ! Find bins for the first tracer coordinate (salinity)
              !mm1(0) = NINT( (zstmp(ji,jj,ikc)-salmin)/salstp ) + 1 
              !mm1(1) = NINT( (zstmp(ip,jj,ikc)-salmin)/salstp ) + 1 
              !mm1(2) = NINT( (zstmp(im,jj,ikc)-salmin)/salstp ) + 1 
              !mm1(3) = NINT( (zstmp(ji,jp,ikc)-salmin)/salstp ) + 1 
              !mm1(4) = NINT( (zstmp(ji,jm,ikc)-salmin)/salstp ) + 1 
              !mm1(5) = NINT( (zstmp(ji,jj,ikn)-salmin)/salstp ) + 1 
              !mm1(6) = NINT( (zstmp(ji,jj,ikp)-salmin)/salstp ) + 1
              
              ! Find bins for the second tracer coordinate (temperature)
              !mm2(0) = NINT( (zttmp(ji,jj,ikc)-temmin)/temstp ) + 1 
              !mm2(1) = NINT( (zttmp(ip,jj,ikc)-temmin)/temstp ) + 1 
              !mm2(2) = NINT( (zttmp(im,jj,ikc)-temmin)/temstp ) + 1 
              !mm2(3) = NINT( (zttmp(ji,jp,ikc)-temmin)/temstp ) + 1 
              !mm2(4) = NINT( (zttmp(ji,jm,ikc)-temmin)/temstp ) + 1 
              !mm2(5) = NINT( (zttmp(ji,jj,ikn)-temmin)/temstp ) + 1 
              !mm2(6) = NINT( (zttmp(ji,jj,ikp)-temmin)/temstp ) + 1
              !
              !WHERE(mm1 < 1) mm1 = 1
              !WHERE(mm1 > nbins) mm1 = nbins
              !WHERE(mm2 < 1) mm2 = 1
              !WHERE(mm2 > nbins) mm2 = nbins
              
              !
              ! Mostly for debug purposes: Calculate meridional overturning in tracer coord
              !
              !psiyz(jbasin,jj,jk)    = psiyz(jbasin,jj,jk)    + zvtmp(ji,jj,ikc) * ibmask(jbasin,ji,jj)
              !psiyz(jbasin,jj,jk)    = psiyz(jbasin,jj,jk)    + zwtmp(ji,jj,ikc) * ibmask(jbasin,ji,jj)
              psiyz(jj,jk,jbasin)    = psiyz(jj,jk,jbasin)    + zwtmp(ji,jj,ikc) * ibmask(jbasin,ji,jj)
              !psiyr(jbasin,jj,mm(0)) = psiyr(jbasin,jj,mm(0)) + vflux(ji,jj,ikc) * ibmask(jbasin,ji,jj)
                  
              !
              ! Stream function using two generalized coordinates (psrr)
              !
              ! The mass flux between two adjacent grid boxes in space or
              ! time is treated as a flux between two grid boxes in 
              ! r-r space (not necessarily adjacent). 
              ! The flux is put along a linear line between the two boxes
              ! 
              
              ! Loop over each direction (east, west, south, north, up, down)
              !
              ! NOTE: Can we vectorize this? Might speed up calculations a lot!
              !
              ! E.g. mk(1:6) = (DBLE(mm2(1:6))-DBLE(mm2(0))) / (DBLE(mm1(1:6))-DBLE(mm1(0)))
              !      
              
              zflux(0) = 0.
              zflux(1) =  zutmp(ji,jj,ikc) * ibmask(jbasin,ji,jj)
              zflux(2) = -zutmp(im,jj,ikc) * ibmask(jbasin,ji,jj)
              zflux(3) =  zvtmp(ji,jj,ikc) * ibmask(jbasin,ji,jj)
              zflux(4) = -zvtmp(ji,jm,ikc) * ibmask(jbasin,ji,jj)
              zflux(5) =  zwtmp(ji,jj,ikc) * ibmask(jbasin,ji,jj)
              zflux(6) = -zwtmp(ji,jj,ikp) * ibmask(jbasin,ji,jj)
              !
              DO jn=1,6
                 IF (zm2(ji,jj,jn) > zm2(ji,jj,0)) THEN
                    !DO m1=mm2(0),mm2(jn)-1 
                    !   psirr(jbasin,ir,m1) = psirr(jbasin,ir,m1) + zflux(jn)
                    !END DO
                    !psirr(jbasin,ir,zm2(ji,jj,jn):zm2(ji,jj,jn)-1) = psirr(jbasin,ir,zm2(ji,jj,jn):zm2(ji,jj,jn)-1) + zflux(jn)
                    psirr(zm1(ji,jj,0),zm2(ji,jj,0):zm2(ji,jj,jn)-1,jbasin) = psirr(zm1(ji,jj,0),zm2(ji,jj,0):zm2(ji,jj,jn)-1,jbasin) + zflux(jn)
                 END IF
              END DO
              !
              ! or this?
              !WHERE (zm2(ji,jj,:) > mm2(0))
              !   psirr(jbasin,zm2(ji,jj,0),:) = psirr(jbasin,zm2(ji,jj,0),:) + zflux(:) 
              !END WHERE
              
              
              
              
           END DO lonLoop ! end of loop on longitude
              
        END DO latLoop ! end of loop on latitude 
        !$OMP END DO
        !DEALLOCATE (dmoc_tmp)
        
        !IF ( lisodep ) DEALLOCATE (depi_tmp)
        !IF ( lisodep ) DEALLOCATE (wdep_tmp)
        
        
        END DO basinLoop ! end loop over basins
        
        !$OMP END PARALLEL
        
        ! Permute k indices
        IF (lprint) PRINT*,' kp,kc,kn were ',ikp,ikc,ikn
        ikt = ikp 
        ikp = ikc
        ikc = ikn
        ikn = ikt 
        IF (lprint) PRINT*,' kp,kc,kn are  ',ikp,ikc,ikn
        
     END DO levelLoop ! end of loop on depths for calculating transports     
     
     !IF ( lisodep ) THEN
     !   WHERE ( wdep(:,:,:) /= 0.d0 ) 
     !      depi(:,:,:) = depi(:,:,:) / wdep (:,:,:)
     !   ELSEWHERE
     !      depi(:,:,:) = rp_spval
     !   END WHERE
     !ENDIF

     ! integrates across bins from highest to lowest density
     !dmoc(:,nbins,:) = dmoc(:,nbins,:)/1.e6
     !DO jbin=nbins-1, 1, -1
     !   dmoc(:,jbin,:) = dmoc(:,jbin+1,:) + dmoc(:,jbin,:)/1.e6
     !END DO  ! loop to next bin
     
     ! integrate in the 'horizontal', i.e. salinity
     
     !IF (lprint) PRINT*,' Entering main bin loop ',jk
     !DO jbasin = 1, nbasins
     !   DO jbin = 2, nbins
     !      psirr(jbasin, jbin, :) = psirr(jbasin, jbin-1, :) + psirr(jbasin,jbin, :)
     !   END DO
     !END DO 
     
     ! Scale to Sverdrup
     psiyz = psiyz / 1e6
     psirr = psirr / 1e6 
     
     ! netcdf output  
     DO jbasin = 1, nbasins
        IF (lprint) PRINT*,' Write basin',jbasin,' to output at step',jt
        DO jk = 1, npk
           !ierr = putvar (ncout, id_varout(jbasin), REAL(psiyz(jbasin,:,jk)), jk, 1, npjglo, ktime = jt)
           ierr = putvar (ncout, id_varout(jbasin), REAL(psiyz(:,jk,jbasin)), jk, 1, npjglo, ktime = jt)
        END DO 
        DO jbin = 1, nbins
           !ierr = putvar (ncout2, id_varout2(jbasin), REAL(psirr(jbasin,:,jbin)), jbin, 1, nbins, ktime = jt)
           ierr = putvar (ncout2, id_varout2(jbasin), REAL(psirr(:,jbin,jbasin)), jbin, 1, nbins, ktime = jt)
        END DO
     END DO

  END DO timeLoop ! time loop
  
  DEALLOCATE( psirr_tmp, volrr_tmp )
  DEALLOCATE( psirr, volrr )
  
  ierr = closeout(ncout)
  ierr = closeout(ncout2)

CONTAINS

  SUBROUTINE CreateOutputFile
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutputFile ***
    !!
    !! ** Purpose :  Initialize and create output files 
    !!
    !! ** Method  :  Check the number of sub_basin, and options 
    !!
    !!----------------------------------------------------------------------

    ! Common to all variables :
    styrvar%cunits            = 'Sverdrup'
    styrvar%rmissing_value    = rp_spval
    styrvar%valid_min         = -1000.
    styrvar%valid_max         =  1000.
    styrvar%scale_factor      = 1.
    styrvar%add_offset        = 0.
    styrvar%savelog10         = 0.
    styrvar%conline_operation = 'N/A'
    styrvar%caxis             = 'TZY'
    
    strrvar%cunits            = 'Sverdrup'
    strrvar%rmissing_value    = rp_spval
    strrvar%valid_min         = -1000.
    strrvar%valid_max         =  1000.
    strrvar%scale_factor      = 1.
    strrvar%add_offset        = 0.
    strrvar%savelog10         = 0.
    strrvar%conline_operation = 'N/A'
    strrvar%caxis             = 'TZZ'

    ipk(:) = npk
    ipr(:) = nbins

    ! Global basin 
    styrvar(npglo)%cname       = cn_zomsfglo
    styrvar(npglo)%clong_name  = 'Meridional.Cell_Global'
    styrvar(npglo)%cshort_name = cn_zomsfglo
    
    strrvar(npglo)%cname       = cn_zothsglo
    strrvar(npglo)%clong_name  = 'Thermohaline.Cell_Global'
    strrvar(npglo)%cshort_name = cn_zothsglo

    IF (lbas) THEN
       strrvar(npatl)%cname       = cn_zothsatl
       strrvar(npatl)%clong_name  = 'Thermohaline.Cell_Atlantic'
       strrvar(npatl)%cshort_name = cn_zothsatl

       strrvar(npinp)%cname       = cn_zothsinp
       strrvar(npinp)%clong_name  = 'Thermohaline.Cell_IndoPacif'
       strrvar(npinp)%cshort_name = cn_zothsinp

       strrvar(npind)%cname       = cn_zothsind
       strrvar(npind)%clong_name  = 'Thermohaline.Cell_Indian'
       strrvar(npind)%cshort_name = cn_zothsind

       strrvar(nppac)%cname       = cn_zothspac
       strrvar(nppac)%clong_name  = 'Thermohaline.Cell_pacif'
       strrvar(nppac)%cshort_name = cn_zothspac
    ENDIF

    IF ( lisodep ) THEN
       ! Global basin
       !strrvar(npglo+nbasins)%cunits      = 'm'
       !strrvar(npglo+nbasins)%cname       = cn_zoisoglo
       !strrvar(npglo+nbasins)%clong_name  = 'Zonal_mean_isopycnal_depth_Global'
       !strrvar(npglo+nbasins)%cshort_name = cn_zoisoglo
       !strrvar(npglo+nbasins)%valid_min   = 0.
       !strrvar(npglo+nbasins)%valid_max   = 8000.
       
       strrvar(npglo+nbasins)%cunits      = 'm3'
       strrvar(npglo+nbasins)%cname       = cn_zovolglo
       strrvar(npglo+nbasins)%clong_name  = 'Volume_Thermohaline_Global'
       strrvar(npglo+nbasins)%cshort_name = cn_zovolglo
       strrvar(npglo+nbasins)%valid_min   = 0.
       strrvar(npglo+nbasins)%valid_max   = 1e10
       
       IF ( lbas ) THEN
          !strrvar(npatl+nbasins)%cunits      = 'm'
          !strrvar(npatl+nbasins)%cname       = cn_zoisoatl
          !strrvar(npatl+nbasins)%clong_name  = 'Zonal_mean_isopycnal_depth_Atlantic'
          !strrvar(npatl+nbasins)%cshort_name = cn_zoisoatl
          !strrvar(npatl+nbasins)%valid_min   = 0.
          !strrvar(npatl+nbasins)%valid_max   = 8000.
          !
          !strrvar(npinp+nbasins)%cunits      = 'm'
          !strrvar(npinp+nbasins)%cname       = cn_zoisoinp
          !strrvar(npinp+nbasins)%clong_name  = 'Zonal_mean_isopycnal_depth_IndoPacif'
          !strrvar(npinp+nbasins)%cshort_name = cn_zoisoinp
          !strrvar(npinp+nbasins)%valid_min   = 0.
          !strrvar(npinp+nbasins)%valid_max   = 8000.
          ! 
          !strrvar(npind+nbasins)%cunits      = 'm'
          !strrvar(npind+nbasins)%cname       = cn_zoisoind
          !strrvar(npind+nbasins)%clong_name  = 'Zonal_mean_isopycnal_depth_Indian'
          !strrvar(npind+nbasins)%cshort_name = cn_zoisoind
          !strrvar(npind+nbasins)%valid_min   = 0.
          !strrvar(npind+nbasins)%valid_max   = 8000.
          !
          !strrvar(nppac+nbasins)%cunits      = 'm'
          !strrvar(nppac+nbasins)%cname       = cn_zoisopac
          !strrvar(nppac+nbasins)%clong_name  = 'Zonal_mean_isopycnal_depth_pacif'
          !strrvar(nppac+nbasins)%cshort_name = cn_zoisopac
          !strrvar(nppac+nbasins)%valid_min   = 0.
          !strrvar(nppac+nbasins)%valid_max   = 8000.
          
          strrvar(npatl+nbasins)%cunits      = 'm3'
          strrvar(npatl+nbasins)%cname       = cn_zovolatl
          strrvar(npatl+nbasins)%clong_name  = 'Volume_Thermohaline_Atlantic'
          strrvar(npatl+nbasins)%cshort_name = cn_zovolglo
          strrvar(npatl+nbasins)%valid_min   = 0.
          strrvar(npatl+nbasins)%valid_max   = 1e10
          
          strrvar(npinp+nbasins)%cunits      = 'm3'
          strrvar(npinp+nbasins)%cname       = cn_zovolinp
          strrvar(npinp+nbasins)%clong_name  = 'Volume_Thermohaline_IndoPacif'
          strrvar(npinp+nbasins)%cshort_name = cn_zovolinp
          strrvar(npinp+nbasins)%valid_min   = 0.
          strrvar(npinp+nbasins)%valid_max   = 1e10
          
          strrvar(npind+nbasins)%cunits      = 'm3'
          strrvar(npind+nbasins)%cname       = cn_zovolind
          strrvar(npind+nbasins)%clong_name  = 'Volume_Thermohaline_Indian'
          strrvar(npind+nbasins)%cshort_name = cn_zovolind
          strrvar(npind+nbasins)%valid_min   = 0.
          strrvar(npind+nbasins)%valid_max   = 1e10
          
          strrvar(nppac+nbasins)%cunits      = 'm3'
          strrvar(nppac+nbasins)%cname       = cn_zovolpac
          strrvar(nppac+nbasins)%clong_name  = 'Volume_Thermohaline_Pacif'
          strrvar(nppac+nbasins)%cshort_name = cn_zovolpac
          strrvar(nppac+nbasins)%valid_min   = 0.
          strrvar(nppac+nbasins)%valid_max   = 1e10
          
       ENDIF
    ENDIF
    
    PRINT*,' Create msf file '
    ncout = create      (cf_msf, 'none', 1,      npjglo, npk,  cdep='sigma')
    ierr  = createvar   (ncout,  styrvar, nvaro, ipk ,id_varout, cdglobal=cglobal)
    ierr  = putheadervar(ncout,  cf_vfil, 1,     npjglo, npk,  pnavlon=rdumlon, pnavlat=rdumlat, pdep=sigma)
    
    PRINT*,' Create ths file '
    ncout2 = create      (cf_ths, 'none', 1,      nbins, nbins,  cdep='sigma')
    ierr  = createvar   (ncout2,  strrvar, nvaro, ipr ,id_varout2, cdglobal=cglobal)
    ierr  = putheadervar(ncout2,  cf_vfil, 1,     nbins, nbins,  pnavlon=rdumlon, pnavlat=rdumlat, pdep=sigma)

    dtim = getvar1d(cf_vfil, cn_vtimec, npt     )
    ierr = putvar1d(ncout,  dtim,       npt, 'T')

  END SUBROUTINE CreateOutputFile

END PROGRAM cdfthermohaline
   
