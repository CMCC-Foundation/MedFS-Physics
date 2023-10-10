!/ ------------------------------------------------------------------- /
      MODULE W3WAVEMD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-May-2009 |
!/                  +-----------------------------------+
!/
!/    04-Feb-2000 : Origination.                        ( version 2.00 )
!/                  For upgrades see subroutines.
!/    14-Feb-2000 : Exact-NL added.                     ( version 2.01 )
!/    05-Jan-2001 : Bug fix to allow model to run       ( version 2.05 )
!/                  without output.
!/    24-Jan-2001 : Flat grid version.                  ( version 2.06 )
!/    09-Feb-2001 : Third propagation scheme added.     ( version 2.08 )
!/    23-Feb-2001 : Check for barrier after source
!/                  terms added ( W3NMIN ).     ( delayed version 2.07 )
!/    16-Mar-2001 : Fourth propagation scheme added.    ( version 2.09 )
!/    30-Mar-2001 : Sub-grid obstacles added.           ( version 2.10 )
!/    23-May-2001 : Clean up and bug fixes.             ( version 2.11 )
!/    10-Dec-2001 : Sub-grid obstacles for UQ schemes.  ( version 2.14 )
!/    11-Jan-2002 : Sub-grid ice.                       ( version 2.15 )
!/    24-Jan-2002 : Zero time step dor data ass.        ( version 2.17 )
!/    18-Feb-2002 : Point output diagnostics added.     ( version 2.18 )
!/    30-Apr-2002 : Add field output types 17-18.       ( version 2.20 )
!/    09-May-2002 : Switch clean up.                    ( version 2.21 )
!/    13-Nov-2002 : Add stress vector.                  ( version 3.00 )
!/    26-Dec-2002 : Moving grid version.                ( version 3.02 )
!/    01-Aug-2003 : Moving grid GSE correction.         ( version 3.03 )
!/    20-Aug-2003 : Output server options added.        ( version 3.04 )
!/    07-Oct-2003 : Output options for NN training.     ( version 3.05 )
!/    29-Dec-2004 : Multiple grid version.              ( version 3.06 )
!/                  W3INIT, W3MPII-O and WWVER moved to w3initmd.ftn
!/    04-Feb-2005 : Add STAMP to par list of W3WAVE.    ( version 3.07 )
!/    04-May-2005 : Change to MPI_COMM_WAVE.            ( version 3.07 )
!/    28-Jun-2005 : Adding map recalc for W3ULEV call.  ( version 3.07 )
!/    07-Sep-2005 : Updated boundary conditions.        ( version 3.08 )
!/                  Fix NRQSG1/2 = 0 array bound issue.
!/    13-Jun-2006 : Split STORE in G/SSTORE             ( version 3.09 )
!/    26-Jun-2006 : Add output type 6.                  ( version 3.09 )
!/    04-Jul-2006 : Consolidate stress arrays.          ( version 3.09 )
!/    18-Oct-2006 : Partitioned spectral data output.   ( version 3.10 )
!/    02-Feb-2007 : Add FLAGST test.                    ( version 3.10 )
!/    02-Apr-2007 : Add partitioned field data.         ( version 3.11 )
!/    07-May-2007 : Bug fix SKIP_O treatment.           ( version 3.11 )
!/    17-May-2007 : Adding NTPROC/NAPROC separation.    ( version 3.11 )
!/    08-Oct-2007 : Adding AS CX-Y to W3SRCE par. list. ( version 3.13 )
!/    22-Feb-2008 : Initialize VGX-Y properly.          ( version 3.13 )
!/    10-Apr-2008 : Bug fix writing log file (MPI).     ( version 3.13 )
!/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
!/
!/    Copyright 2009 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS. 
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!  2. Variables and types :
!
!  3. Subroutines and functions :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      W3WAVE    Subr. Public   Actual wave model.
!      W3GATH    Subr. Public   Data transpose before propagation.
!      W3SCAT    Subr. Public   Data transpose after propagation.
!      W3NMIN    Subr. Public   Calculate minimum number of sea
!                               points per processor.
!     ----------------------------------------------------------------
!
!  4. Subroutines and functions used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3SETx    Subr. W3xDATMD Point to data structure.
!
!      W3UCUR    Subr. W3UPDTMD Interpolate current fields in time.
!      W3UWND    Subr. W3UPDTMD Interpolate wind fields in time.
!      W3UINI    Subr. W3UPDTMD Update initial conditions if init.
!                               with initial wind conditions.
!      W3UBPT    Subr. W3UPDTMD Update boundary points.
!      W3UICE    Subr. W3UPDTMD Update ice coverage.
!      W3ULEV    Subr. W3UPDTMD Transform the wavenumber grid.
!      W3DDXY    Subr. W3UPDTMD Calculate dirivatives of the depth.
!      W3DCXY    Subr. W3UPDTMD Calculate dirivatives of the current.
!
!      W3MAPn    Subr. W3PROnMD Preparation for  ropagation schemes.
!      W3XYPn    Subr. W3PROnMD Longitude-latitude ("XY") propagation.
!      W3KTPn    Subr. W3PROnMD Intra-spectral ("k-theta") propagation.
!
!      W3SRCE    Subr. W3SRCEMD Source term integration and calculation.
!
!      W3IOGR    Subr. W3IOGRMD Reading/writing model definition file.
!      W3OUTG    Subr. W3IOGOMD Generate gridded output fields.
!      W3IOGO    Subr. W3IOGOMD Read/write gridded output.
!      W3IOPE    Subr. W3IOPOMD Extract point output.
!      W3IOPO    Subr. W3IOPOMD Read/write point output.
!      W3IOTR    Subr. W3IOTRMD Process spectral output along tracks.
!      W3IORS    Subr. W3IORSMD Read/write restart files.
!      W3IOBC    Subr. W3IOBCMD Read/write boundary conditions.
!      W3CPRT    Subr. W3IOSFMD Partition spectra.
!      W3IOSF    Subr.   Id.    Write partitioned spectral data.
!
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!      WWTIME    Subr.   Id.    System time in readable format.
!      EXTCDE    Subr.   Id.    Program abort.
!
!      TICK21    Subr. W3TIMEMD Advance the clock.
!      DSEC21    Func.   Id.    Difference between times.
!      STME21    Subr.   Id.    Time in readable format.
!
!      MPI_BARRIER, MPI_STARTALL, MPI_WAITALL
!                Subr.          Basic MPI routines.
!     ----------------------------------------------------------------
!
!  5. Remarks :
!
!  6. Switches :
!
!       !/SHRD  Switch for shared / distributed memory architecture.
!       !/DIST  Id.
!       !/MPI   Id.
!
!       !/LLG   Spherical grid.
!       !/XYG   Cartesian grid.
!
!       !/PR1   First order propagation schemes.
!       !/PR2   ULTIMATE QUICKEST scheme.
!       !/PR3   Averaged ULTIMATE QUICKEST scheme.
!       !/PRX   User-defined scheme.
!
!       !/S     Enable subroutine tracing.
!       !/T     Test output.
!       !/MPIT  Test output for MPI specific code.
!
!  7. Source code :
!
!/ ------------------------------------------------------------------- /
!/MPI      USE W3ADATMD, ONLY: MPIBUF
!
      PUBLIC
!/MPI INTEGER :: FINALTIME(2) !Required to couple with NEMO [DRU]
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3WAVE ( IMOD, TEND, STAMP, NO_OUT )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         10-Apr-2008 |
!/                  +-----------------------------------+
!/
!/    17-Mar-1999 : Distributed FORTRAN 77 version.     ( version 1.18 )
!/    04-Feb-2000 : Upgrade to FORTRAN 90               ( version 2.00 )
!/                  Major changes to logistics.
!/    05-Jan-2001 : Bug fix to allow model to run       ( version 2.05 )
!/                  without output.
!/    24-Jan-2001 : Flat grid version.                  ( version 2.06 )
!/    09-Feb-2001 : Third propagation scheme added.     ( version 2.08 )
!/    23-Feb-2001 : Check for barrier after source
!/                  terms added ( W3NMIN ).     ( delayed version 2.07 )
!/    16-Mar-2001 : Fourth propagation scheme added.    ( version 2.09 )
!/    30-Mar-2001 : Sub-grid obstacles added.           ( version 2.10 )
!/    23-May-2001 : Barrier added for dry run, changed  ( version 2.10 )
!/                  declaration of FLIWND.
!/    10-Dec-2001 : Sub-grid obstacles for UQ schemes.  ( version 2.14 )
!/    11-Jan-2002 : Sub-grid ice.                       ( version 2.15 )
!/    24-Jan-2002 : Zero time step dor data ass.        ( version 2.17 )
!/    09-May-2002 : Switch clean up.                    ( version 2.21 )
!/    13-Nov-2002 : Add stress vector.                  ( version 3.00 )
!/    26-Dec-2002 : Moving grid version.                ( version 3.02 )
!/    01-Aug-2003 : Moving grid GSE correction.         ( version 3.03 )
!/    07-Oct-2003 : Output options for NN training.     ( version 3.05 )
!/    29-Dec-2004 : Multiple grid version.              ( version 3.06 )
!/    04-Feb-2005 : Add STAMP to par list.              ( version 3.07 )
!/    04-May-2005 : Change to MPI_COMM_WAVE.            ( version 3.07 )
!/    28-Jun-2005 : Adding map recalc for W3ULEV call.  ( version 3.07 )
!/    07-Sep-2005 : Updated boundary conditions.        ( version 3.08 )
!/    26-Jun-2006 : Add output type 6.                  ( version 3.09 )
!/    04-Jul-2006 : Consolidate stress arrays.          ( version 3.09 )
!/    18-Oct-2006 : Partitioned spectral data output.   ( version 3.10 )
!/    02-Feb-2007 : Add FLAGST test.                    ( version 3.10 )
!/    02-Apr-2007 : Add partitioned field data.         ( version 3.11 )
!/                  Improve MPI_WAITALL call tests/allocations.
!/    07-May-2007 : Bug fix SKIP_O treatment.           ( version 3.11 )
!/    17-May-2007 : Adding NTPROC/NAPROC separation.    ( version 3.11 )
!/    08-Oct-2007 : Adding AS CX-Y to W3SRCE par. list. ( version 3.13 )
!/    22-Feb-2008 : Initialize VGX-Y properly.          ( version 3.13 )
!/    10-Apr-2008 : Bug fix writing log file (MPI).     ( version 3.13 )
!/
!  1. Purpose :
!
!     Run WAVEWATCH III for a given time interval.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IMOD    Int.   I   Model number.
!       TEND    I.A.   I   Ending time of integration.
!       STAMP   Log.   I   Print time stamp (optional, defaults to T).
!       NO_OUT  Log.   I   Skip output (optional, defaults to F).
!                          Skip at ending time only!
!     ----------------------------------------------------------------
!
!     Local parameters : Flags
!     ----------------------------------------------------------------
!       FLOUTG  Log.  Flag for running W3OUTG.
!       FLPART  Log.  Flag for running W3CPRT.
!       FLZERO  Log.  Flag for zero time interval.
!       FLAG0   Log.  Flag for processors without tasks.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!     See module documentation.
!
!  5. Called by :
!
!     Any program shell or integrated model which uses WAVEWATCH III.
!
!  6. Error messages :
!
!  7. Remarks :
!
!     - Currents are updated before winds as currents are used in wind
!       and USTAR processing.
!     - Ice and water levels can be updated only once per call.
!     - If ice or water level time are undefined, the update
!       takes place asap, otherwise around the "half-way point"
!       betweem the old and new times.
!     - To increase accuracy, the calculation o fthe intra-sprectral
!       propagation is performed around the spatial propagation.
!
!  8. Structure :
!
!     -----------------------------------------------------------
!       0.  Initializations
!         a Point to data structures
!         b Subroutine tracing
!         c Local parameter initialization
!         d Test output
!       1.  Check the consistency of the input.
!         a Ending time versus initial time.
!         b Water level time.
!         c Current time interval.
!         d Wind time interval.
!         e Ice time.
!       2.  Determine next time from ending and output
!           time and get corresponding time step.
!       3.  Loop over time steps (see below).
!       4.  Perform output to file if requested.
!         a Check if time is output time.
!         b Processing and MPP preparations.  ( W3CPRT, W3OUTG )
!         c Reset next output time.
!        -------------- loop over output types ------------------
!         d Perform output.                           ( W3IOxx )
!         e Update next output time.
!        -------------------- end loop --------------------------
!       5.  Update log file.
!       6.  If time is not ending time, branch back to 2.
!     -----------------------------------------------------------
!
!      Section 3.
!     ----------------------------------------------------------
!       3.1  Interpolate winds and currents. ( W3UCUR, W3DCXY )
!                                                    ( W3UWND )
!                                                    ( W3UINI )
!       3.2  Update boundary conditions.     ( W3IOBC, W3UBPT )
!       3.3  Update ice coverage (if new ice map).   ( W3UICE )
!       3.4  Transform grid (if new water level).    ( W3ULEV )
!       3.5  Update maps and dirivatives.    ( W3MAPn, W3DDXY )
!                                            ( W3NMIN, W3UTRN )
!            Update grid advection vector.
!       3.6  Perform propagation
!          a Preparations.
!          b Intra spectral part 1.                  ( W3KTPn )
!          c Longitude-latitude       ( W3GATH, W3XYPn W3SCAT )
!          b Intra spectral part 2.                  ( W3KTPn )
!       3.7  Calculate and integrate source terms.   ( W3SRCE )
!       3.8  Update global time step.
!     ----------------------------------------------------------
!
!  9. Switches :
!
!     See module documentation.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE CONSTANTS
!/
      USE W3GDATMD
      USE W3WDATMD
      USE W3ADATMD
      USE W3IDATMD
      USE W3ODATMD
!/
      USE W3UPDTMD
      USE W3SRCEMD
!/PR1      USE W3PRO1MD
!/PR2      USE W3PRO2MD
!/PR3      USE W3PRO3MD
!/PRX      USE W3PROXMD
!/
      USE W3IOGRMD
      USE W3IOGOMD
      USE W3IOPOMD
      USE W3IOTRMD
      USE W3IORSMD
      USE W3IOBCMD
      USE W3IOSFMD
!/
      USE W3SERVMD
      USE W3TIMEMD
!/MPI      USE NemoWwCoupling    !Required to couple with NEMO [DRU]
!emanuela
      USE W3INITMD, ONLY: FLCOUP 
!emanuela
!
      IMPLICIT NONE
!
!/MPI      INCLUDE "mpif.h"
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)           :: IMOD, TEND(2)
      INTEGER                      :: IERR
      LOGICAL, INTENT(IN), OPTIONAL :: STAMP, NO_OUT
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters :
!/
!/T      INTEGER                 :: ILEN
!/S      INTEGER, SAVE           :: IENT = 0
      INTEGER                 :: TCALC(2), IT, IT0, NT, ITEST,        &
                                 ITLOC, ITLOCH, NTLOC, ISEA, JSEA,    &
                                 IX, IY, ISPEC, J, TOUT(2), TLST(2)
!/MPI      INTEGER                 :: IERR_MPI, NRQMAX
!/MPI      INTEGER, ALLOCATABLE    :: STATCO(:,:), STATIO(:,:)
      REAL                    :: DTTST, DTTST1, DTTST2, DTTST3,       &
                                 DTL0, DTI0, DTGA, DTG, DTRES,        &
                                 FAC, VGX, VGY, FACK, FACTH,          &
                                 FACX, FACY, SXM, SYM, XXX
      REAL, ALLOCATABLE       :: FIELD(:)
      LOGICAL                 :: FLACT, FLZERO, FLFRST, FLMAP, TSTAMP,&
                                 SKIP_O, FLAG_O, FLDDIR, READBC,      &
                                 FLAG0=.FALSE., FLOUTG, FLPFLD,       &
                                 FLPART, LOCAL
!/MPI      LOGICAL                 :: FLGMPI(6)
      CHARACTER(LEN=8)        :: STTIME
      CHARACTER(LEN=11)       :: IDACT, OUTID
      CHARACTER(LEN=23)       :: IDTIME
!
!/MPI REAL(8), DIMENSION(:), ALLOCATABLE :: dT_AS !Required to couple with NEMO [DRU]
!/MPI REAL(8), DIMENSION(:), ALLOCATABLE :: sdu !Required to couple with NEMO [DRU]
!/MPI REAL(8), DIMENSION(:), ALLOCATABLE :: sdv !Required to couple with NEMO [DRU]
!/
!/ ------------------------------------------------------------------- /
! 0.  Initializations
! 0.a Set pointers to data structure
!
!/MPI ALLOCATE(dT_AS(NSEA)) !Required to couple with NEMO [DRU]
!/MPI ALLOCATE(sdu(NSEA)) !Required to couple with NEMO [DRU]
!/MPI ALLOCATE(sdv(NSEA)) !Required to couple with NEMO [DRU]
      IF ( IOUTP  .NE. IMOD ) CALL W3SETO ( IMOD, NDSE, NDST )
      IF ( IGRID  .NE. IMOD ) CALL W3SETG ( IMOD, NDSE, NDST )
      IF ( IWDATA .NE. IMOD ) CALL W3SETW ( IMOD, NDSE, NDST )
      IF ( IADATA .NE. IMOD ) CALL W3SETA ( IMOD, NDSE, NDST )
      IF ( IIDATA .NE. IMOD ) CALL W3SETI ( IMOD, NDSE, NDST )
!
      IF ( PRESENT(STAMP) ) THEN
          TSTAMP = STAMP
        ELSE
          TSTAMP = .TRUE.
        END IF
!
      IF ( PRESENT(NO_OUT) ) THEN
          SKIP_O = NO_OUT
        ELSE
          SKIP_O = .FALSE.
        END IF
!
! 0.b Subroutine tracing
!
!/S      CALL STRACE (IENT, 'W3WAVE')
!
! 0.c Local parameter initialization
!
      IPASS  = IPASS + 1
      IDACT  = '           '
      OUTID  = '           '
      FLACT  = ITIME .EQ. 0
      FLMAP  = ITIME .EQ. 0
      FLDDIR = ITIME .EQ. 0 .AND. ( FLCTH .OR. FLCK )
!
      FLPFLD = .FALSE.
      DO J=15,20
        FLPFLD = FLPFLD .OR. FLOGRD(J)
        END DO
!
      IF ( IAPROC .EQ. NAPLOG ) BACKSPACE ( NDSO )
!
      IF ( FLCOLD ) THEN
          DTDYN = 0.
          FCUT  = SIG(NK) * TPIINV
        END IF
!
      ALLOCATE ( FIELD(1-NY:NY*(NX+2)) )
      FIELD = 0.
!
      LOCAL   = IAPROC .LE. NAPROC
!
! 0.c Test output
!
!/T      ILEN   = LEN_TRIM(FILEXT)
!/T      WRITE (NDST,9000) IMOD, FILEXT(:ILEN), TEND
!
! 1.  Check the consistency of the input ----------------------------- /
! 1.a Ending time versus initial time
!
      DTTST  = DSEC21 ( TIME , TEND )
      FLZERO = DTTST .EQ. 0.
!/T      WRITE (NDST,9010) DTTST, FLZERO
      IF ( DTTST .LT. 0. ) THEN
          IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1000)
          CALL EXTCDE ( 1 )
        END IF
!
! 1.b Water level time
!
      IF ( FLLEV ) THEN
          IF ( TLEV(1) .GE. 0. ) THEN
              DTL0   = DSEC21 ( TLEV , TLN )
            ELSE
              DTL0   = 1.
            END IF
!/T          WRITE (NDST,9011) DTL0
          IF ( DTL0 .LT. 0. ) THEN
              IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1001)
              CALL EXTCDE ( 2 )
            END IF
        ELSE
          DTL0   = 0.
        END IF
!
! 1.c Current interval
!
      IF ( FLCUR ) THEN
          DTTST1 = DSEC21 ( TC0 , TCN )
          DTTST2 = DSEC21 ( TC0 , TIME )
          DTTST3 = DSEC21 ( TEND , TCN )
!/T          WRITE (NDST,9012) DTTST1, DTTST2, DTTST3
          IF ( DTTST1.LT.0. .OR. DTTST2.LT.0. .OR. DTTST3.LT.0. ) THEN
              IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1002)
              CALL EXTCDE ( 3 )
            END IF
          IF ( DTTST2.EQ.0..AND. ITIME.EQ.0 ) THEN
              IDACT(7:7) = 'F'
              TOFRST = TIME
            END IF
        END IF
!
! 1.d Wind interval
!
      IF ( FLWIND ) THEN
          DTTST1 = DSEC21 ( TW0 , TWN )
          DTTST2 = DSEC21 ( TW0 , TIME )
          DTTST3 = DSEC21 ( TEND , TWN )
!/T          WRITE (NDST,9013) DTTST1, DTTST2, DTTST3
          IF ( DTTST1.LT.0. .OR. DTTST2.LT.0. .OR. DTTST3.LT.0. ) THEN
              IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1003)
              CALL EXTCDE ( 4 )
            END IF
          IF ( DTTST2.EQ.0..AND. ITIME.EQ.0 ) THEN
              IDACT(3:3) = 'F'
              TOFRST = TIME
            END IF
        END IF
!
! 1.e Ice time
!
      IF ( FLICE ) THEN
          IF ( TICE(1) .GE. 0 ) THEN
              DTI0   = DSEC21 ( TICE , TIN )
            ELSE
              DTI0   = 1.
            END IF
!/T          WRITE (NDST,9014) DTI0
          IF ( DTI0 .LT. 0. ) THEN
              IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1004)
              CALL EXTCDE ( 5 )
            END IF
        ELSE
          DTI0   = 0.
        END IF
!
! 2.  Determine next time from ending and output --------------------- /
!     time and get corresponding time step.
!
      FLFRST = .TRUE.
      DO
!
        IF ( TOFRST(1) .GT. 0 ) THEN
            DTTST  = DSEC21 ( TEND , TOFRST )
          ELSE
            DTTST  = 0.
          ENDIF
!
        IF ( DTTST.GE.0. ) THEN
            TCALC = TEND
          ELSE
            TCALC = TOFRST
          END IF
!
        DTTST  = DSEC21 ( TIME , TCALC )
        NT     = 1 + INT ( DTTST / DTMAX - 0.001 )
        DTGA   = DTTST / REAL(NT)
        IF ( DTTST .EQ. 0. ) THEN
            IT0    = 0
            IF ( .NOT.FLZERO ) ITIME  = ITIME - 1
            NT     = 0
          ELSE
            IT0    = 1
          END IF
!
!/T        WRITE (NDST,9020) IT0, NT, DTGA
!
! ==================================================================== /
! 3.  Loop over time steps
!
        DTRES  = 0.
!
        DO IT=IT0, NT
!
          ITIME  = ITIME + 1
          DTG    = REAL(NINT(DTGA+DTRES+0.0001))
          DTRES  = DTRES + DTGA - DTG
          IF ( ABS(DTRES) .LT. 0.001 ) DTRES  = 0.
          CALL TICK21 ( TIME , DTG )
!
          IF ( TSTAMP .AND. SCREEN.NE.NDSO .AND. IAPROC.EQ.NAPOUT ) THEN
              CALL WWTIME ( STTIME )
              CALL STME21 ( TIME , IDTIME )
              WRITE (SCREEN,*) "DRU loop on time", IAPROC, TIME, TEND, FINALTIME
              WRITE (SCREEN,950) IDTIME, STTIME
            END IF
!
          VGX = 0.
          VGY = 0. 
          IF(FLAGS(8)) THEN
              DTTST1 = DSEC21 ( TIME, TGN ) 
              DTTST2 = DSEC21 ( TG0, TGN )
              FAC    = DTTST1 / MAX ( 1. , DTTST2 )
              VGX    = (FAC*GA0+(1.-FAC)*GAN) *                       &
                            COS(FAC*GD0+(1.-FAC)*GDN)
              VGY    = (FAC*GA0+(1.-FAC)*GAN) *                       &
                            SIN(FAC*GD0+(1.-FAC)*GDN)
            END IF
!
!/T        WRITE (NDST,9021) ITIME, IT, TIME, FLMAP, FLDDIR,          &
!/T                          VGX, VGY, DTG, DTRES
!
! 3.1 Interpolate winds and currents.
!     (Initialize wave fields with winds)
!
              IF ( IAPROC.EQ.NAPOUT ) THEN
                WRITE (SCREEN,*) "DRU before loop on see", IAPROC, TIME, ITIME, IT, MOD(TIME(2),10000),NSEA,NAPROC,NSEAL
                END IF
!/MPI              IF ( ITIME.EQ.0 ) THEN
!/MPI                    IF ( IAPROC.EQ.NAPOUT ) WRITE (SCREEN,*) "[DRU] NemoWwHandshake to receive t,u,v at itime",itime
!/MPI                    CALL NemoWwHandshake(sdt_w=dT_AS,sdu_w=sdu,sdv_w=sdv,cds_w=CDS,             &
!/MPI                                         gdxd=871,gdyd=253,                                   &
!/MPI                                         nsea=NSEA,iaproc=IAPROC,naproc=NAPROC,mapsf=MAPSF)
!!ini  emanuela
!/MPI            IF (FLCOUP) THEN
!/MPI              CX(1:NSEA)=sdu
!/MPI              CY(1:NSEA)=sdv
!/MPI              AS(1:NSEA)=dT_AS
!/MPI              END IF
!/MPI                   END IF
!!emanuela
          IF ( IAPROC.EQ.NAPOUT ) THEN
                WRITE (SCREEN,*) "ITIME=",ITIME
                WRITE (SCREEN,*) "CDSMIN1=",MINVAL(CDS)
                WRITE (SCREEN,*) "CDSMAX1=",MAXVAL(CDS)
                WRITE (SCREEN,*) "CXMIN1=",MINVAL(CX)
                WRITE (SCREEN,*) "CXMAX1=",MAXVAL(CX)
                WRITE (SCREEN,*) "ASMIN1=",MINVAL(AS)
                WRITE (SCREEN,*) "ASMAX1=",MAXVAL(AS)
                WRITE (SCREEN,*) "NSEA,NSEAL,NSEALM,JSEA,NAPROC="  &
                                 ,NSEA,NSEAL,NSEALM,JSEA,NAPROC
          END IF
          IF ( FLCUR  ) THEN
              CALL W3UCUR ( FLFRST ) 
               CALL W3DCXY
!ini emanuela
!            ELSE IF ( FLFRST) THEN
             ELSE IF (( FLFRST).AND. .NOT.(FLCOUP) ) THEN
              CX = 0.
              CY = 0.
            END IF
!
          IF ( FLWIND ) THEN
               IF ( FLFRST ) ASF = 1.
               CALL W3UWND ( FLFRST, VGX, VGY )
            ELSE IF ( FLFRST ) THEN
              U10    = 0.01
              U10D   = 0.
              UST    = 0.05
              USTDIR = 0.05
            END IF
!
          IF ( FLIWND .AND. LOCAL ) CALL W3UINI ( VA )
!
! 3.2 Update boundary conditions.
!
          IF ( FLBPI .AND. LOCAL ) THEN
!
              DO
                IF ( TBPIN(1) .EQ. -1 ) THEN
                    READBC = .TRUE.
                    IDACT(1:1) = 'F'
                  ELSE
                    READBC = DSEC21(TIME,TBPIN).LT.0.
                    IF (READBC.AND.IDACT(1:1).EQ.' ') IDACT(1:1) = 'X'
                  END IF
                FLACT  = READBC .OR. FLACT
                IF ( READBC ) THEN
                    CALL W3IOBC ( 'READ', NDS(9), TBPI0, TBPIN,       &
                                  ITEST, IMOD )
                    IF ( ITEST .NE. 1 ) CALL W3UBPT 
                  END IF
                IF ( ITEST .LT. 0 ) IDACT(1:1) = 'L'
                IF ( ITEST .GT. 0 ) IDACT(1:1) = ' '
                IF ( .NOT. (READBC.AND.FLBPI) ) EXIT
                END DO
!
            END IF
!
! 3.3 Update ice coverage (if new ice map).
!     Need to be run on output nodes too, to update MAPSTx
!
          IF ( FLICE .AND. DTI0.NE.0. ) THEN
!
              IF ( TICE(1).GE.0 ) THEN
                  IF ( DTI0 .LT. 0. ) THEN
                      IDACT(9:9) = 'B'
                    ELSE
                      DTTST  = DSEC21 ( TIME, TIN )
                      IF ( DTTST .LE. 0.5*DTI0 ) IDACT(9:9) = 'U'
                    END IF
                ELSE
                  IDACT(9:9) = 'I'
                END IF
!
              IF ( IDACT(9:9).NE.' ' ) THEN
                  CALL W3UICE ( VA, VA )
                  DTI0   = 0.
                  FLACT  = .TRUE.
                  FLMAP  = .TRUE.
                END IF
!
            END IF
!
! 3.4 Transform grid (if new water level).
!
          IF ( FLLEV .AND. DTL0.NE.0. ) THEN
!
              IF ( TLEV(1).GE.0 ) THEN
                  IF ( DTL0 .LT. 0. ) THEN
                      IDACT(5:5) = 'B'
                    ELSE
                      DTTST  = DSEC21 ( TIME, TLN )
                      IF ( DTTST .LE. 0.5*DTL0 ) IDACT(5:5) = 'U'
                    END IF
                ELSE
                  IDACT(5:5) = 'I'
                END IF
!
              IF ( IDACT(5:5).NE.' ' ) THEN
                  CALL W3ULEV ( VA, VA )
                  DTL0   = 0.
                  FLACT  = .TRUE.
                  FLMAP  = .TRUE.
                  FLDDIR = FLDDIR .OR. FLCTH .OR. FLCK
                END IF
!
            END IF
!
! 3.5 Update maps and dirivatives.
!
          IF ( FLMAP ) THEN
!/PR1              CALL W3MAP1 ( MAPSTA )
!/PR2              CALL W3MAP2
!/PR3              CALL W3MAP3
!/PRX              CALL W3MAPX
              CALL W3UTRN ( TRNX, TRNY )
!/PR3              CALL W3MAPT
              CALL W3NMIN ( MAPSTA, FLAG0 )
              IF ( FLAG0 .AND. IAPROC.EQ.NAPERR ) WRITE (NDSE,1030) IMOD
              FLMAP  = .FALSE.
            END IF
!
          IF ( FLDDIR ) THEN
              CALL W3DDXY
              FLDDIR = .FALSE.
            END IF
!
          FLIWND = .FALSE.
          FLFRST = .FALSE.
!
          IF ( FLZERO ) THEN
!/T              WRITE (NDST,9022)
              GOTO 400
            END IF
          IF ( FLDRY .OR. IT.EQ.0 .OR. IAPROC.GT.NAPROC ) THEN
!/T              WRITE (NDST,9023)
              GOTO 380
            END IF
!
! 3.6 Perform Propagation = = = = = = = = = = = = = = = = = = = = = = =
! 3.6.1 Preparations
!
          NTLOC  = 1 + INT( DTG/DTCFLI - 0.001 )
          ITLOCH = ( NTLOC + 1 - MOD(ITIME,2) ) / 2
!
          FACTH  = DTG / DTH / REAL(NTLOC)
          FACK   = DTG / REAL(NTLOC)
!
! 3.6.2 Intra-spectral part 1
!
          IF ( FLCTH .OR. FLCK ) THEN
!
              DO ITLOC=1, ITLOCH
!/OMP0/!$OMP PARALLEL DO PRIVATE (JSEA,ISEA,IX,IY) SCHEDULE (DYNAMIC,1)
                DO JSEA=1, NSEAL
!/DIST                  ISEA   = IAPROC + (JSEA-1)*NAPROC
!/SHRD                  ISEA   = JSEA
                  IX     = MAPSF(ISEA,1)
                  IY     = MAPSF(ISEA,2)
                  IF ( MAPSTA(IY,IX) .EQ. 1 ) THEN
!/PR1                      CALL W3KTP1 ( ISEA, FACTH, FACK, CTHG0(IY), &
!/PR1                           CG(:,ISEA), WN(:,ISEA), DW(ISEA),      &
!/PR1                           DDDX(IY,IX), DDDY(IY,IX), CX(ISEA),    &
!/PR1                           CY(ISEA), DCXDX(IY,IX), DCXDY(IY,IX),  &
!/PR1                           DCYDX(IY,IX), DCYDY(IY,IX), VA(:,JSEA) )
!/PR2                      CALL W3KTP2 ( ISEA, FACTH, FACK, CTHG0(IY), &
!/PR2                           CG(:,ISEA), WN(:,ISEA), DW(ISEA),      &
!/PR2                           DDDX(IY,IX), DDDY(IY,IX), CX(ISEA),    &
!/PR2                           CY(ISEA), DCXDX(IY,IX), DCXDY(IY,IX),  &
!/PR2                           DCYDX(IY,IX), DCYDY(IY,IX), VA(:,JSEA) )
!/PR3                      CALL W3KTP3 ( ISEA, FACTH, FACK, CTHG0(IY), &
!/PR3                           CG(:,ISEA), WN(:,ISEA), DW(ISEA),      &
!/PR3                           DDDX(IY,IX), DDDY(IY,IX), CX(ISEA),    &
!/PR3                           CY(ISEA), DCXDX(IY,IX), DCXDY(IY,IX),  &
!/PR3                           DCYDX(IY,IX), DCYDY(IY,IX), VA(:,JSEA) )
!/PRX                      CALL W3KTPX
                    END IF
                  END DO
                END DO
            END IF
!
! 3.6.3 Longitude-latitude
!       (time step correction in routine)
!
          IF ( FLCX .OR. FLCY ) THEN
!
              IF ( FLCX ) THEN
!/LLG                  FACX   = DTG / ( SX * DERA * RADIUS )
!/XYG                  FACX   = DTG / SX
!/LLG                  SXM    = SX * DERA * RADIUS
!/XYG                  SXM    = SX
                ELSE
                  FACX   = 0.
                END IF
              IF ( FLCY ) THEN
!/LLG                  FACY   = DTG / ( SY * DERA * RADIUS )
!/XYG                  FACY   = DTG / SY
!/LLG                  SYM    = SY * DERA * RADIUS
!/XYG                  SYM    = SY
                ELSE
                  FACY   = 0.
                END IF
!
!/MPI              IF ( NRQSG1 .GT. 0 ) THEN
!/MPI                  CALL MPI_STARTALL (NRQSG1, IRQSG1(1,1), IERR_MPI)
!/MPI                  CALL MPI_STARTALL (NRQSG1, IRQSG1(1,2), IERR_MPI)
!/MPI                END IF
!
!/OMP0/!$OMP PARALLEL DO PRIVATE (ISPEC,FIELD) SCHEDULE (DYNAMIC,1)
              DO ISPEC=1, NSPEC
                IF ( IAPPRO(ISPEC) .EQ. IAPROC ) THEN
                    CALL W3GATH ( ISPEC, FIELD )
!/PR1                    CALL W3XYP1 ( ISPEC, FACX, FACY, DTG, MAPSTA, &
!/PR1                                  FIELD, VGX, VGY )
!/PR2                    CALL W3XYP2 ( ISPEC, FACX, FACY, DTG, MAPSTA, &
!/PR2                                  MAPFS, FIELD, VGX, VGY )
!/PR3                    CALL W3XYP3 ( ISPEC, FACX, FACY, DTG, MAPSTA, &
!/PR3                                  MAPFS,  FIELD, VGX, VGY )
!/PRX                    CALL W3XYPX
                     CALL W3SCAT ( ISPEC, MAPSTA, FIELD )
                  END IF
                END DO
!
!/MPI              IF ( NRQSG1 .GT. 0 ) THEN
!/MPI                  ALLOCATE ( STATCO(MPI_STATUS_SIZE,NRQSG1) )
!/MPI                  CALL MPI_WAITALL (NRQSG1, IRQSG1(1,1), STATCO, &
!/MPI                                    IERR_MPI)
!/MPI                  CALL MPI_WAITALL (NRQSG1, IRQSG1(1,2), STATCO, &
!/MPI                                    IERR_MPI)
!/MPI                  DEALLOCATE ( STATCO )
!/MPI                END IF
!
            END IF
!
! 3.6.4 Intra-spectral part 2
!
          IF ( FLCTH .OR. FLCK ) THEN
!
              DO ITLOC=ITLOCH+1, NTLOC
!/OMP0/!$OMP PARALLEL DO PRIVATE (JSEA,ISEA,IX,IY) SCHEDULE (DYNAMIC,1)
                DO JSEA=1, NSEAL
!/DIST                  ISEA   = IAPROC + (JSEA-1)*NAPROC
!/SHRD                  ISEA   = JSEA
                  IX     = MAPSF(ISEA,1)
                  IY     = MAPSF(ISEA,2)
                  IF ( MAPSTA(IY,IX) .EQ. 1 ) THEN
!/PR1                      CALL W3KTP1 ( ISEA, FACTH, FACK, CTHG0(IY), &
!/PR1                           CG(:,ISEA), WN(:,ISEA), DW(ISEA),      &
!/PR1                           DDDX(IY,IX), DDDY(IY,IX), CX(ISEA),    &
!/PR1                           CY(ISEA), DCXDX(IY,IX), DCXDY(IY,IX),  &
!/PR1                           DCYDX(IY,IX), DCYDY(IY,IX), VA(:,JSEA) )
!/PR2                      CALL W3KTP2 ( ISEA, FACTH, FACK, CTHG0(IY), &
!/PR2                           CG(:,ISEA), WN(:,ISEA), DW(ISEA),      &
!/PR2                           DDDX(IY,IX), DDDY(IY,IX), CX(ISEA),    &
!/PR2                           CY(ISEA), DCXDX(IY,IX), DCXDY(IY,IX),  &
!/PR2                           DCYDX(IY,IX), DCYDY(IY,IX), VA(:,JSEA) )
!/PR3                      CALL W3KTP3 ( ISEA, FACTH, FACK, CTHG0(IY), &
!/PR3                           CG(:,ISEA), WN(:,ISEA), DW(ISEA),      &
!/PR3                           DDDX(IY,IX), DDDY(IY,IX), CX(ISEA),    &
!/PR3                           CY(ISEA), DCXDX(IY,IX), DCXDY(IY,IX),  &
!/PR3                           DCYDX(IY,IX), DCYDY(IY,IX), VA(:,JSEA) )
!/PRX                      CALL W3KTPX
                    END IF
                  END DO
                END DO
            END IF
!
! 3.6 End propapgation  = = = = = = = = = = = = = = = = = = = = = = = =
!
! 3.7 Calculate and integrate source terms.
!
!
          IF ( FLSOU ) THEN
!
!/OMP0/!$OMP PARALLEL DO PRIVATE (JSEA,ISEA,IX,IY) SCHEDULE (DYNAMIC,1)
              DO JSEA=1, NSEAL
!/DIST                ISEA   = IAPROC + (JSEA-1)*NAPROC
!/SHRD                ISEA   = JSEA
                IX     = MAPSF(ISEA,1)
                IY     = MAPSF(ISEA,2)
                IF ( MAPSTA(IY,IX).EQ.1 .AND. FLAGST(ISEA) ) THEN
                     CALL W3SRCE ( IX, IY, IMOD, VA(:,JSEA),          &
                          ALPHA(1:NK,JSEA), WN(1:NK,ISEA),            &
                          CG(1:NK,ISEA), DW(ISEA), U10(ISEA),         &
                          U10D(ISEA), AS(ISEA), UST(ISEA),            &
                          USTDIR(ISEA), CX(ISEA), CY(ISEA),           &
                          EMN(ISEA), FMN(ISEA), WNM(ISEA), AMX(ISEA), &
                          FPIS(ISEA), CDS(ISEA), Z0S(ISEA),           &
                          DTDYN(ISEA), FCUT(ISEA), DTG )
                  ELSE
                    UST   (ISEA) = UNDEF
                    USTDIR(ISEA) = UNDEF
                    DTDYN (ISEA) = UNDEF
                    FCUT  (ISEA) = UNDEF
                  END IF
                END DO
!
              IF ( IAPROC.EQ.NAPOUT ) WRITE (SCREEN,*) "DRU after loop on see", IAPROC, TIME, ITIME, IT, &
                                      MOD(TIME(2),10000),NSEA,NAPROC,NSEAL,FINALTIME
!/MPI              IF ( ( MOD(TIME(2),10000).EQ.0 ) .AND. .NOT. ( TIME(1).EQ.FINALTIME(1) .AND. TIME(2).EQ.FINALTIME(2) ) ) THEN
!/MPI                    IF ( IAPROC.EQ.NAPOUT ) WRITE (SCREEN,*) "[DRU] NemoWwHandshake to receive t,u,v and send d at itime",itime
!/MPI                    ! CALL NemoWwHandshake(TIME,TIME,TIME,TIME,2)
!/MPI                    CALL NemoWwHandshake(sdt_w=dT_AS,sdu_w=sdu,sdv_w=sdv,cds_w=CDS,             &
!/MPI                                         gdxd=871,gdyd=253,                                   &
!/MPI                                         nsea=NSEA,iaproc=IAPROC,naproc=NAPROC,mapsf=MAPSF)
!ini emanuela
!/MPI             IF (FLCOUP) THEN
!/MPI               CX(1:NSEA)=sdu
!/MPI               CY(1:NSEA)=sdv
!/MPI               AS(1:NSEA)=dT_AS
!/MPI               END IF
!/MPI                    END IF
!/MPI              IF ( ITIME.EQ.1 ) THEN
!/MPI                 IF ( IAPROC.EQ.NAPOUT ) WRITE (SCREEN,*) "[DRU] NemoWwHandshake to send d at itime",itime
!/MPI                 ! CALL NemoWwHandshake(TIME,TIME,TIME,TIME,2)
!/MPI                    CALL NemoWwHandshake(sdt_w=dT_AS,sdu_w=sdu,sdv_w=sdv,cds_w=CDS,             &
!/MPI                                         gdxd=871,gdyd=253,                                   &
!/MPI                                         nsea=NSEA,iaproc=IAPROC,naproc=NAPROC,mapsf=MAPSF)
!/MPI                 END IF
!
! This barrier is from older code versions. It has been removed in 3.11
! to optimize IO2/3 settings. May be needed on some systems still
!
!!/MPI              IF (FLAG0) CALL MPI_BARRIER (MPI_COMM_WCMP,IERR_MPI)
!!/MPI            ELSE
!!/MPI              CALL MPI_BARRIER (MPI_COMM_WCMP,IERR_MPI)
!
            END IF
 !emanuela
          IF ( IAPROC.EQ.NAPOUT ) THEN
                WRITE (SCREEN,*) "ITIME=",ITIME
                WRITE (SCREEN,*) "CDSMIN3=",MINVAL(CDS)
                WRITE (SCREEN,*) "CDSMAX3=",MAXVAL(CDS)
                WRITE (SCREEN,*) "CXMIN3",MINVAL(CX)
                WRITE (SCREEN,*) "CXMAX3",MAXVAL(CX)
                WRITE (SCREEN,*) "ASMIN3=",MINVAL(AS)
                WRITE (SCREEN,*) "ASMAX3=",MAXVAL(AS)
          END IF
!
! 3.8 Update global time step.
!     (Branch point FLDRY, IT=0)
!
  380     CONTINUE
!
          IF (IT.NE.NT) THEN
              DTTST  = DSEC21 ( TIME , TCALC )
              DTG    = DTTST / REAL(NT-IT)
            END IF
!
          IF ( FLACT .AND. IT.NE.NT .AND. IAPROC.EQ.NAPLOG ) THEN
              CALL STME21 ( TIME , IDTIME )
              IF ( IDLAST .NE. TIME(1) ) THEN
                  WRITE (NDSO,900) ITIME, IPASS, IDTIME(01:19),       &
                                   IDACT, OUTID
                  IDLAST = TIME(1)
                ELSE
                  WRITE (NDSO,901) ITIME, IPASS, IDTIME(12:19),       &
                                   IDACT, OUTID
                END IF
              FLACT  = .FALSE.
              IDACT  = '         '
            END IF
!
          END DO
!
!/T      WRITE (NDST,9030)
!
!     End of loop over time steps
! ==================================================================== /
!
  400 CONTINUE
!
! 4.  Perform output to file if requested ---------------------------- /
! 4.a Check if time is output time
!     Delay if data assimilation time.
!
        IF ( TOFRST(1)  .EQ. -1 ) THEN
            DTTST  = 1.
          ELSE
            DTTST   = DSEC21 ( TIME, TOFRST )
          END IF
!
        IF ( TDN(1)  .EQ. -1 ) THEN
            DTTST1 = 1.
          ELSE
            DTTST1  = DSEC21 ( TIME, TDN )
          END IF
!
        DTTST2 = DSEC21 ( TIME, TEND )
        FLAG_O = .NOT.SKIP_O .OR. ( SKIP_O .AND. DTTST2.NE.0. )
!
!/T        WRITE (NDST,9040) TOFRST, TDN, DTTST, DTTST1, FLAG_O
!
        IF ( DTTST.LE.0. .AND. DTTST1.NE.0. .AND. FLAG_O ) THEN
!
!/T          WRITE (NDST,9041)
!
! 4.b Processing and MPP preparations
!
            IF ( FLOUT(1) ) THEN
                FLOUTG = DSEC21(TIME,TONEXT(:,1)).EQ.0.
              ELSE
                FLOUTG = .FALSE.
              END IF
!
            FLPART = .FALSE.
            IF ( FLOUT(1) .AND. FLPFLD )                               &
                 FLPART = FLPART .OR. DSEC21(TIME,TONEXT(:,1)).EQ.0.
            IF ( FLOUT(6) )                                            &
                 FLPART = FLPART .OR. DSEC21(TIME,TONEXT(:,6)).EQ.0.
!
!/T            WRITE (NDST,9042) LOCAL, FLPART, FLOUTG
!
            IF ( LOCAL .AND. FLPART ) CALL W3CPRT ( IMOD )
            IF ( LOCAL .AND. FLOUTG ) CALL W3OUTG ( VA, FLPFLD )
!
!/MPI            FLGMPI = .FALSE.
!/MPI            NRQMAX = 0
!
!/MPI            IF ( FLOUT(1) .AND. NRQGO.NE.0 ) THEN
!/MPI                IF ( DSEC21(TIME,TONEXT(:,1)).EQ.0. ) THEN
!/MPI                    CALL MPI_STARTALL ( NRQGO, IRQGO , IERR_MPI )
!/MPI                    FLGMPI(1) = .TRUE.
!/MPI                    NRQMAX    = MAX ( NRQMAX , NRQGO )
!/MPIT                    WRITE (NDST,9043) '1 ', NRQGO, NRQMAX, NAPFLD
!/MPI                  END IF
!/MPI              END IF
!
!/MPI            IF ( FLOUT(2) .AND. NRQPO.NE.0 ) THEN
!/MPI                IF ( DSEC21(TIME,TONEXT(:,2)).EQ.0. ) THEN
!/MPI                    CALL MPI_STARTALL ( NRQPO, IRQPO1, IERR_MPI )
!/MPI                    FLGMPI(2) = .TRUE.
!/MPI                    NRQMAX    = MAX ( NRQMAX , NRQPO )
!/MPIT                    WRITE (NDST,9043) '2 ', NRQPO, NRQMAX, NAPPNT
!/MPI                  END IF
!/MPI              END IF
!
!/MPI            IF ( FLOUT(4) .AND. NRQRS.NE.0 ) THEN
!/MPI                IF ( DSEC21(TIME,TONEXT(:,4)).EQ.0. ) THEN
!/MPI                    CALL MPI_STARTALL ( NRQRS, IRQRS , IERR_MPI )
!/MPI                    FLGMPI(4) = .TRUE.
!/MPI                    NRQMAX    = MAX ( NRQMAX , NRQRS )
!/MPIT                    WRITE (NDST,9043) '4 ', NRQRS, NRQMAX, NAPRST
!/MPI                  END IF
!/MPI              END IF
!
!/MPI            IF ( FLOUT(5) .AND. NRQBP.NE.0 ) THEN
!/MPI                IF ( DSEC21(TIME,TONEXT(:,5)).EQ.0. ) THEN
!/MPI                    CALL MPI_STARTALL ( NRQBP , IRQBP1, IERR_MPI )
!/MPI                    FLGMPI(5) = .TRUE.
!/MPI                    NRQMAX    = MAX ( NRQMAX , NRQBP )
!/MPIT                    WRITE (NDST,9043) '5a', NRQBP, NRQMAX, NAPBPT
!/MPI                  END IF
!/MPI              END IF
!
!/MPI            IF ( FLOUT(5) .AND. NRQBP2.NE.0 .AND.                &
!/MPI                 IAPROC.EQ.NAPBPT) THEN
!/MPI                IF ( DSEC21(TIME,TONEXT(:,5)).EQ.0. ) THEN
!/MPI                    CALL MPI_STARTALL (NRQBP2,IRQBP2,IERR_MPI)
!/MPI                    NRQMAX    = MAX ( NRQMAX , NRQBP2 )
!/MPIT                    WRITE (NDST,9043) '5b', NRQBP2, NRQMAX, NAPBPT
!/MPI                  END IF
!/MPI              END IF
!
!/MPI           IF ( NRQMAX .NE. 0 ) ALLOCATE                         &
!/MPI                                 ( STATIO(MPI_STATUS_SIZE,NRQMAX) )
!
! 4.c Reset next output time
!
            TOFRST(1) = -1
            TOFRST(2) =  0
!
            DO J=1, 6
              IF ( FLOUT(J) ) THEN
!
! 4.d Perform output
!
                  TOUT(:) = TONEXT(:,J)
                  DTTST   = DSEC21 ( TIME, TOUT )
!
                  IF ( DTTST .EQ. 0. ) THEN
                      IF ( J .EQ. 1 ) THEN
                          IF ( IAPROC .EQ. NAPFLD ) THEN
!/MPI                              IF ( FLGMPI(1) ) CALL MPI_WAITALL  &
!/MPI                                 ( NRQGO, IRQGO, STATIO, IERR_MPI )
!/MPI                              FLGMPI(1) = .FALSE.
                              CALL W3IOGO                             &
                                 ( 'WRITE', NDS(7), ITEST, IMOD )
                            END IF
                        ELSE IF ( J .EQ. 2 ) THEN
                          IF ( IAPROC .EQ. NAPPNT ) THEN
                              CALL W3IOPE ( VA )
                              CALL W3IOPO                             &
                                 ( 'WRITE', NDS(8), ITEST, IMOD )
                            END IF
                        ELSE IF ( J .EQ. 3 ) THEN
                          CALL W3IOTR ( NDS(11), NDS(12), VA, IMOD )
                        ELSE IF ( J .EQ. 4 ) THEN
                          CALL W3IORS ('HOT', NDS(6), XXX, ITEST, IMOD )
                        ELSE IF ( J .EQ. 5 ) THEN
                          IF ( IAPROC .EQ. NAPBPT ) THEN
!/MPI                              IF (NRQBP2.NE.0) CALL MPI_WAITALL  &
!/MPI                                ( NRQBP2, IRQBP2,STATIO, IERR_MPI )
                              CALL W3IOBC ( 'WRITE', NDS(10),         &
                                            TIME, TIME, ITEST, IMOD )
                            END IF
                        ELSE IF ( J .EQ. 6 ) THEN
                          CALL W3IOSF ( NDS(13), IMOD )
                        END IF
!
                      CALL TICK21 ( TOUT, DTOUT(J) )
                      TONEXT(:,J) = TOUT
                      TLST        = TOLAST(:,J)
                      DTTST       = DSEC21 ( TOUT , TLST )
                      FLOUT(J)    = DTTST.GE.0.
                      IF ( FLOUT(J) ) THEN
                          OUTID(2*J-1:2*J-1) = 'X'
                        ELSE
                          OUTID(2*J-1:2*J-1) = 'L'
                        END IF
                    END IF
!
! 4.e Update next output time
!
                  IF ( FLOUT(J) ) THEN
                      IF ( TOFRST(1).EQ.-1 ) THEN
                          TOFRST = TOUT
                        ELSE
                          DTTST  = DSEC21 ( TOUT , TOFRST )
                          IF ( DTTST.GT.0.) THEN
                              TOFRST = TOUT
                            END IF
                        END IF
                    END IF
!
                END IF
!
              END DO
!
!/MPI            IF ( FLGMPI(1) ) CALL MPI_WAITALL                    &
!/MPI                             ( NRQGO, IRQGO , STATIO, IERR_MPI )
!/MPI            IF ( FLGMPI(2) ) CALL MPI_WAITALL                    &
!/MPI                             ( NRQPO, IRQPO1, STATIO, IERR_MPI )
!/MPI            IF ( FLGMPI(4) ) CALL MPI_WAITALL                    &
!/MPI                             ( NRQRS, IRQRS , STATIO, IERR_MPI )
!/MPI            IF ( FLGMPI(5) ) CALL MPI_WAITALL                    &
!/MPI                             ( NRQBP, IRQBP1, STATIO, IERR_MPI )
!/MPI            IF ( NRQMAX .NE. 0 ) DEALLOCATE ( STATIO )
!
!/T          WRITE (NDST,9044)
!
! This barrier is from older code versions. It has been removed in 3.11
! to optimize IO2/3 settings. May be needed on some systems still
!
!!/MPI            IF (FLDRY) CALL MPI_BARRIER (MPI_COMM_WAVE,IERR_MPI) 
!
          END IF
!
! 5.  Update log file ------------------------------------------------ /
!
        IF ( IAPROC.EQ.NAPLOG ) THEN
!
            CALL STME21 ( TIME , IDTIME )
            IF ( FLCUR ) THEN
                DTTST  = DSEC21 ( TIME , TCN )
                IF ( DTTST .EQ. 0. ) IDACT(7:7) = 'X'
              END IF
            IF ( FLWIND ) THEN
                DTTST  = DSEC21 ( TIME , TWN )
                IF ( DTTST .EQ. 0. ) IDACT(3:3) = 'X'
              END IF
            IF ( TDN(1) .GT. 0  ) THEN
                DTTST  = DSEC21 ( TIME , TDN )
                IF ( DTTST .EQ. 0. ) IDACT(11:11) = 'X'
              END IF
!
            IF ( IDLAST.NE.TIME(1) ) THEN
                WRITE (NDSO,900) ITIME, IPASS, IDTIME(1:19),          &
                                 IDACT, OUTID
                IDLAST = TIME(1)
              ELSE 
                WRITE (NDSO,901) ITIME, IPASS, IDTIME(12:19),         &
                                 IDACT, OUTID
              END IF
!
          END IF
!
        IDACT  = '         '
        OUTID  = '           '
        FLACT  = .FALSE.
!
! 6.  If time is not ending time, branch back to 2 ------------------- /
!
        DTTST  = DSEC21 ( TIME, TEND )
        IF ( DTTST .EQ. 0. ) EXIT
        END DO
!
      IF ( IAPROC .EQ. NAPLOG ) WRITE (NDSO,902)
!
      DEALLOCATE ( FIELD )
!
      RETURN
!
! Formats
!
  900 FORMAT (8X,I6,'  |',I4,'  | ', A19  ,' | ',A,' | ',A,' |')
  901 FORMAT (8X,I6,'  |',I4,'  | ',11X,A8,' | ',A,' | ',A,' |')
  902 FORMAT (8X,'--------+------+---------------------+'             &
                ,'-------------+-------------+')
!
  950 FORMAT ('  WAVEWATCH III calculating for ',A,' at ',A)
!
 1000 FORMAT (/' *** WAVEWATCH III ERROR IN W3WAVE :'/                &
               '     ENDING TIME BEFORE STARTING TIME '/)
 1001 FORMAT (/' *** WAVEWATCH III ERROR IN W3WAVE :'/                &
               '     NEW WATER LEVEL BEFORE OLD WATER LEVEL '/)
 1002 FORMAT (/' *** WAVEWATCH III ERROR IN W3WAVE :'/                &
               '     ILLEGAL CURRENT INTERVAL '/)
 1003 FORMAT (/' *** WAVEWATCH III ERROR IN W3WAVE :'/                &
               '     ILLEGAL WIND INTERVAL '/)
 1004 FORMAT (/' *** WAVEWATCH III ERROR IN W3WAVE :'/                &
               '     NEW ICE FIELD BEFORE OLD ICE FIELD '/)
 1030 FORMAT (/' *** WAVEWATCH III WARING IN W3WAVE :'/               &
               '     AT LEAST ONE PROCESSOR HAS 0 ACTIVE POINTS',     &
                   ' IN GRID',I3)
!
!/T 9000 FORMAT (                                                     &
!/T   '============================================================', &
!/T   '===================='/                                         &
!/T   ' TEST W3WAVE : RUN MODEL',I3,' FILEXT [',A,                    &
!/T                    '] UP TO ',I8.8,I7.6 /                         &
!/T   '====================',                                         &
!/T '============================================================')
!/T 9010 FORMAT (' TEST W3WAVE : DT INT. =',F12.1,'   FLZERO = ',L1)
!/T 9011 FORMAT (' TEST W3WAVE : DT LEV. =',F12.1)
!/T 9012 FORMAT (' TEST W3WAVE : DT CUR. =',F12.1/                    &
!/T              '                        ',F12.1/                    &
!/T              '                        ',F12.1)
!/T 9013 FORMAT (' TEST W3WAVE : DT WIND =',F12.1/                    &
!/T              '                        ',F12.1/                    &
!/T              '                        ',F12.1)
!/T 9014 FORMAT (' TEST W3WAVE : DT ICE  =',F12.1)
!/T 9020 FORMAT (' TEST W3WAVE : IT0, NT, DTG :',2I4,F8.1)
!/T 9021 FORMAT (' TEST W3WAVE : ITIME etc',I6,I4,I10.8,I7.6,1X,2L1,  &
!/T                                         2F6.2,F7.1,F6.2)
!/T 9022 FORMAT (' TEST W3WAVE : SKIP TO 400 IN 3.5')
!/T 9023 FORMAT (' TEST W3WAVE : SKIP TO 380 IN 3.5')
!/T 9030 FORMAT (' TEST W3WAVE : END OF COMPUTATION LOOP')
!/T 9040 FORMAT (' TEST W3WAVE : CHECKING FOR OUTPUT'/                &
!/T              '               TOFRST           :',I9.8,I7.6/       &
!/T              '               TND              :',I9.8,I7.6/       &
!/T              '               DTTST[1], FLAG_O :',2F8.1,L4)
!/T 9041 FORMAT (' TEST W3WAVE : PERFORMING OUTPUT')
!/T 9042 FORMAT (' TEST W3WAVE : OUTPUT COMPUTATION FLAGS: '3L2)
!/MPIT 9043 FORMAT (' TEST W3WAVE : TYPE, NRQ, NRQMAX, NA : ',A2,3I6)
!/T 9044 FORMAT (' TEST W3WAVE : END OF OUTPUT')
!/
!/ End of W3WAVE ----------------------------------------------------- /
!/
      END SUBROUTINE W3WAVE
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3GATH ( ISPEC, FIELD )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         13-Jun-2006 |
!/                  +-----------------------------------+
!/
!/    04-Jan-1999 : Distributed FORTRAN 77 version.     ( version 1.18 )
!/    13-Jan-2000 : Upgrade to FORTRAN 90               ( version 2.00 )
!/                  Major changes to logistics.
!/    29-Dec-2004 : Multiple grid version.              ( version 3.06 )
!/    13-Jun-2006 : Split STORE in G/SSTORE             ( version 3.09 )
!/
!  1. Purpose :
!
!     Gather spectral bin information into a propagation field array.
!
!  2. Method :
!
!     Direct copy or communication calls (MPP version).
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       ISPEC   Int.   I   Spectral bin considered.
!       FIELD   R.A.   O   Full field to be propagated.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!
!      MPI_STARTALL, MPI_WAITALL
!                Subr. mpif.h   MPI persistent comm. routines (!/MPI).
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3WAVE    Subr. W3WAVEMD Actual wave model routine.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!       None.
!
!  7. Remarks :
!
!     - The field is extracted but not converted.
!     - Array FIELD is not initialized.
!     - MPI version requires posing of send and receive calls in
!       W3WAVE to match local calls.
!     - MPI version does not require an MPI_TESTALL call for the
!       posted gather operation as MPI_WAITALL is mandatory to
!       reset persistent communication for next time step.
!     - MPI version allows only two new pre-fetch postings per
!       call to minimize chances to be slowed down by gathers that
!       are not yet needed, while maximizing the pre-loading
!       during the early (low-frequency) calls to the routine
!       where the amount of calculation needed for proagation is
!       the largest.
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/SHRD  Switch for message passing method.
!     !/MPI   Id.
!
!     !/S     Enable subroutine tracing.
!     !/MPIT  MPI test output.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/S      USE W3SERVMD, ONLY: STRACE
!/
      USE W3GDATMD, ONLY: NSPEC, NX, NY, NSEA, NSEAL, MAPSF
      USE W3WDATMD, ONLY: A => VA
!/MPI      USE W3ADATMD, ONLY: MPIBUF, BSTAT, IBFLOC, ISPLOC, BISPL, &
!/MPI                          NSPLOC, NRQSG2, IRQSG2, GSTORE
!/MPI      USE W3ODATMD, ONLY: NDST, IAPROC, NAPROC
!/
      IMPLICIT NONE
!
!/MPI      INCLUDE "mpif.h"
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: ISPEC
      REAL, INTENT(OUT)       :: FIELD(1-NY:NY*(NX+2))
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
!/SHRD      INTEGER                 :: ISEA, IXY
!/MPI      INTEGER                 :: STATUS(MPI_STATUS_SIZE,NSPEC),  &
!/MPI                                 IOFF, IERR_MPI, JSEA, ISEA,     &
!/MPI                                 IXY, IS0, IB0, NPST, J
!/S      INTEGER, SAVE           :: IENT
!/MPIT      CHARACTER(LEN=15)       :: STR(MPIBUF), STRT
!/
!/ ------------------------------------------------------------------- /
!/
!/S      CALL STRACE (IENT, 'W3GATH')
!
!      FIELD  = 0.
!
! 1.  Shared memory version ------------------------------------------ /
!
!/SHRD      DO ISEA=1, NSEA
!/SHRD        IXY        = MAPSF(ISEA,3)
!/SHRD        FIELD(IXY) = A(ISPEC,ISEA)
!/SHRD        END DO
!
!/SHRD      RETURN
!
! 2.  Distributed memory version ( MPI ) ----------------------------- /
! 2.a Update counters
!
!/MPI      ISPLOC = ISPLOC + 1
!/MPI      IBFLOC = IBFLOC + 1
!/MPI      IF ( IBFLOC .GT. MPIBUF ) IBFLOC = 1
!
!/MPIT      IF ( ISPLOC .EQ. 1 ) THEN
!/MPIT          STR = '--------------+'
!/MPIT          WRITE (NDST,9000) STR
!/MPIT        END IF
!/MPIT      STR    = '              |'
!/MPIT      STRT   = STR(IBFLOC)
!/MPIT      STRT(9:9) = 'A'
!
! 2.b Check status of present buffer
! 2.b.1 Scatter (send) still in progress, wait to end
!
!/MPI      IF ( BSTAT(IBFLOC) .EQ. 2 ) THEN
!/MPI          IOFF =  1 + (BISPL(IBFLOC)-1)*NRQSG2
!/MPI          IF ( NRQSG2 .GT. 0 ) CALL                              &
!/MPI               MPI_WAITALL ( NRQSG2, IRQSG2(IOFF,2),             &
!/MPI                             STATUS, IERR_MPI )
!/MPI          BSTAT(IBFLOC) = 0
!/MPIT          STRT(13:13) = 'S'
!/MPI        END IF
!
! 2.b.2 Gather (recv) not yet posted, post now
!
!/MPI      IF ( BSTAT(IBFLOC) .EQ. 0 ) THEN
!/MPI          BSTAT(IBFLOC) = 1
!/MPI          BISPL(IBFLOC) = ISPLOC
!/MPI          IOFF =  1 + (ISPLOC-1)*NRQSG2
!/MPI          IF ( NRQSG2 .GT. 0 ) CALL                              &
!/MPI               MPI_STARTALL ( NRQSG2, IRQSG2(IOFF,1), IERR_MPI )
!/MPIT          STRT(10:10) = 'g'
!/MPI        END IF
!
! 2.c Put local spectral densities in store
!
!/MPI      DO JSEA=1, NSEAL
!/MPI        GSTORE(IAPROC+(JSEA-1)*NAPROC,IBFLOC) = A(ISPEC,JSEA)
!/MPI        END DO
!
! 2.d Wait for remote spectral densities
!
!/MPI      IOFF =  1 + (BISPL(IBFLOC)-1)*NRQSG2
!/MPI      IF ( NRQSG2 .GT. 0 ) CALL                                  &
!/MPI           MPI_WAITALL ( NRQSG2, IRQSG2(IOFF,1), STATUS, IERR_MPI )
!
!/MPIT      STRT(11:11) = 'G'
!/MPIT      WRITE (STRT(1:7),'(I2,I5)') BSTAT(IBFLOC), ISPLOC
!/MPIT      STR(IBFLOC) = STRT
!
! 2.e Convert storage array to field.
!
!/MPI      DO ISEA=1, NSEA
!/MPI        IXY        = MAPSF(ISEA,3)
!/MPI        FIELD(IXY) = GSTORE(ISEA,IBFLOC)
!/MPI        END DO
!
! 2.f Pre-fetch data in available buffers
!
!/MPI      IS0    = ISPLOC
!/MPI      IB0    = IBFLOC
!/MPI      NPST   = 0
!
!/MPI      DO J=1, MPIBUF-1
!/MPI        IS0    = IS0 + 1
!/MPI        IF ( IS0 .GT. NSPLOC ) EXIT
!/MPI        IB0    = 1 + MOD(IB0,MPIBUF)
!/MPI        IF ( BSTAT(IB0) .EQ. 0 ) THEN
!/MPI            BSTAT(IB0) = 1
!/MPI            BISPL(IB0) = IS0
!/MPI            IOFF       = 1 + (IS0-1)*NRQSG2
!/MPI            IF ( NRQSG2 .GT. 0 ) CALL                            &
!/MPI                 MPI_STARTALL ( NRQSG2, IRQSG2(IOFF,1), IERR_MPI )
!/MPI            NPST       = NPST + 1
!/MPIT            STRT        = STR(IB0)
!/MPIT            STRT(10:10) = 'g'
!/MPIT            WRITE (STRT(1:7),'(I2,I5)') BSTAT(IB0), BISPL(IB0)
!/MPIT            STR(IB0)    = STRT
!/MPI          END IF
!/MPI        IF ( NPST .GE. 2 ) EXIT
!/MPI        END DO
!
! 2.g Test output
!
!/MPIT      DO IB0=1, MPIBUF
!/MPIT        STRT   = STR(IB0)
!/MPIT        IF ( STRT(2:2) .EQ. ' ' ) THEN
!/MPIT            IF ( BSTAT(IB0) .EQ. 0 ) THEN
!/MPIT                WRITE (STRT(1:2),'(I2)') BSTAT(IB0)
!/MPIT              ELSE
!/MPIT                WRITE (STRT(1:7),'(I2,I5)') BSTAT(IB0), BISPL(IB0)
!/MPIT              END IF
!/MPIT            STR(IB0) = STRT
!/MPIT          END IF
!/MPIT        END DO
!/MPIT      WRITE (NDST,9010) ISPLOC, STR
!
!/MPI      RETURN
!
! Formats
!
!/MPIT 9000 FORMAT ( ' TEST OF BUFFER MANAGEMENT MPI :'/              &
!/MPIT               ' -------------------------------'/              &
!/MPIT      '      RECORDS ALTERNATELY WRITTEN BY W3GATH AND W3SCAT'/ &
!/MPIT      '      FRIST COLLUMN  : LOCAL ISPEC'/                     &
!/MPIT      '      OTHER COLLUMNS : BUFFER STATUS INDICATOR '/        &
!/MPIT      '                        0 : INACTIVE'/                   &
!/MPIT      '                        1 : RECEIVING'/                  &
!/MPIT      '                        2 : SENDING'/                    &
!/MPIT      '                       LOCAL ISPEC FOR BUFFER'/          &
!/MPIT      '                       A  : ACTIVE BUFFER'/              &
!/MPIT      '                       g/G: START/FINISH RECIEVE'/       &
!/MPIT      '                       s/S: START/FINISH SEND'/          &
!/MPIT      ' +-----+',8A15)
!/MPIT 9010 FORMAT ( ' |',I4,' |',8A15)
!/
!/ End of W3GATH ----------------------------------------------------- /
!/
      END SUBROUTINE W3GATH
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3SCAT ( ISPEC, MAPSTA, FIELD )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         13-Jun-2006 |
!/                  +-----------------------------------+
!/
!/    04-Jan-1999 : Distributed FORTRAN 77 version.     ( version 1.18 )
!/    13-Jan-2000 : Upgrade to FORTRAN 90               ( version 2.00 )
!/                  Major changes to logistics.
!/    28-Dec-2004 : Multiple grid version.              ( version 3.06 )
!/    07-Sep-2005 : Updated boundary conditions.        ( version 3.08 )
!/    13-Jun-2006 : Split STORE in G/SSTORE             ( version 3.09 )
!/
!  1. Purpose :
!
!     'Scatter' data back to spectral storage after propagation.
!
!  2. Method :
!
!     Direct copy or communication calls (MPP version).
!     See also W3GATH.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       ISPEC   Int.   I   Spectral bin considered.
!       MAPSTA  I.A.   I   Status map for spatial grid.
!       FIELD   R.A.   I   Full field to be propagated.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!
!      MPI_STARTALL, MPI_WAITALL, MPI_TESTALL
!                Subr. mpif.h   MPI persistent comm. routines (!/MPI).
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3WAVE    Subr. W3WAVEMD Actual wave model routine.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!     None.
!
!  7. Remarks :
!
!     - The field is put back but not converted !
!     - MPI persistent communication calls initialize in W3MPII.
!     - See W3GATH and W3MPII for additional comments on data
!       buffering.
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/SHRD  Switch for message passing method.
!     !/MPI   Id.
!
!     !/S     Enable subroutine tracing.
!     !/MPIT  MPI test output.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/S      USE W3SERVMD, ONLY: STRACE
!/
      USE W3GDATMD, ONLY: NSPEC, NX, NY, NSEA, NSEAL, MAPSF
      USE W3WDATMD, ONLY: A => VA
!/MPI      USE W3ADATMD, ONLY: MPIBUF, BSTAT, IBFLOC, ISPLOC, BISPL, &
!/MPI                          NSPLOC, NRQSG2, IRQSG2, SSTORE
      USE W3ODATMD, ONLY: NDST
!/MPI      USE W3ODATMD, ONLY: IAPROC, NAPROC
!/
      IMPLICIT NONE
!
!/MPI      INCLUDE "mpif.h"
!/ 
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: ISPEC, MAPSTA(NY*NX)
      REAL, INTENT(IN)        :: FIELD(1-NY:NY*(NX+2))
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
!/SHRD      INTEGER                 :: ISEA, IXY
!/MPI      INTEGER                 :: ISEA, IXY, IOFF, IERR_MPI, J,   &
!/MPI                                 STATUS(MPI_STATUS_SIZE,NSPEC),  &
!/MPI                                 JSEA, IB0
!/S      INTEGER, SAVE           :: IENT
!/MPIT      CHARACTER(LEN=15)       :: STR(MPIBUF), STRT
!/MPI      LOGICAL                 :: DONE
!/
!/ ------------------------------------------------------------------- /
!/
!/S      CALL STRACE (IENT, 'W3SCAT')
!
! 1.  Shared memory version ------------------------------------------ *
!
!/SHRD      DO ISEA=1, NSEA
!/SHRD        IXY           = MAPSF(ISEA,3)
!/SHRD        IF ( MAPSTA(IXY) .GE. 1 ) A(ISPEC,ISEA) = FIELD(IXY)
!/SHRD        END DO
!
!/SHRD      RETURN
!
! 2.  Distributed memory version ( MPI ) ----------------------------- *
! 2.a Initializations
!
!/MPIT      DO IB0=1, MPIBUF
!/MPIT        STR(IB0) = '              |'
!/MPIT        END DO
!
!/MPIT      STRT   = STR(IBFLOC)
!/MPIT      STRT(9:9) = 'A'
!
! 2.b Convert full grid to sea grid, active points only
!
!/MPI      DO ISEA=1, NSEA
!/MPI        IXY    = MAPSF(ISEA,3)
!/MPI        IF ( MAPSTA(IXY) .GE. 1 ) SSTORE(ISEA,IBFLOC) = FIELD(IXY)
!/MPI        END DO
!
! 2.c Send spectral densities to appropriate remote
!
!/MPI      IOFF   = 1 + (ISPLOC-1)*NRQSG2
!/MPI      IF ( NRQSG2 .GT. 0 ) CALL                                  &
!/MPI           MPI_STARTALL ( NRQSG2, IRQSG2(IOFF,2), IERR_MPI )
!/MPI      BSTAT(IBFLOC) = 2
!/MPIT      STRT(12:12) = 's'
!/MPIT      WRITE (STRT(1:7),'(I2,I5)') BSTAT(IBFLOC), ISPLOC
!/MPIT      STR(IBFLOC) = STRT
!
! 2.d Save locally stored results
!
!/MPI      DO JSEA=1, NSEAL
!/MPI        ISEA   = IAPROC+(JSEA-1)*NAPROC
!/MPI        IXY    = MAPSF(ISEA,3)
!/MPI        IF (MAPSTA(IXY) .GE. 1) A(ISPEC,JSEA) = SSTORE(ISEA,IBFLOC)
!/MPI        END DO
!
! 2.e Check if any sends have finished
!
!/MPI      IB0    = IBFLOC
!
!/MPI      DO J=1, MPIBUF
!/MPI        IB0    = 1 + MOD(IB0,MPIBUF)
!/MPI        IF ( BSTAT(IB0) .EQ. 2 ) THEN
!/MPI            IOFF   = 1 + (BISPL(IB0)-1)*NRQSG2
!/MPI            IF ( NRQSG2 .GT. 0 ) THEN
!/MPI               CALL MPI_TESTALL ( NRQSG2, IRQSG2(IOFF,2), DONE,  &
!/MPI                                 STATUS, IERR_MPI )
!/MPI              ELSE
!/MPI                DONE   = .TRUE.
!/MPI              END IF
!/MPI            IF ( DONE .AND. NRQSG2.GT.0 ) CALL                   &
!/MPI                     MPI_WAITALL ( NRQSG2, IRQSG2(IOFF,2),       &
!/MPI                                   STATUS, IERR_MPI )
!/MPI            IF ( DONE ) THEN
!/MPI                BSTAT(IB0) = 0
!/MPIT                STRT        = STR(IB0)
!/MPIT                WRITE (STRT(1:7),'(I2,I5)') BSTAT(IB0), BISPL(IB0)
!/MPIT                STRT(13:13) = 'S'
!/MPIT                STR(IB0)    = STRT
!/MPI              END IF
!/MPI          END IF
!/MPI        END DO
!
! 2.f Last component, finish message passing, reset buffer control
!
!/MPI      IF ( ISPLOC .EQ. NSPLOC ) THEN
!
!/MPI          DO IB0=1, MPIBUF
!/MPI            IF ( BSTAT(IB0) .EQ. 2 ) THEN
!/MPI                IOFF   = 1 + (BISPL(IB0)-1)*NRQSG2
!/MPI                IF ( NRQSG2 .GT. 0 ) CALL                        &
!/MPI                     MPI_WAITALL ( NRQSG2, IRQSG2(IOFF,2),       &
!/MPI                                   STATUS, IERR_MPI )
!/MPI                BSTAT(IB0) = 0
!/MPIT                STRT        = STR(IB0)
!/MPIT                WRITE (STRT(1:7),'(I2,I5)') BSTAT(IB0), BISPL(IB0)
!/MPIT                STRT(13:13) = 'S'
!/MPIT                STR(IB0)    = STRT
!/MPI              END IF
!/MPI            END DO
!
!/MPI          ISPLOC = 0
!/MPI          IBFLOC = 0
!
!/MPI        END IF
!
! 2.g Test output
!
!/MPIT      DO IB0=1, MPIBUF
!/MPIT        STRT   = STR(IB0)
!/MPIT        IF ( STRT(2:2) .EQ. ' ' ) THEN
!/MPIT            IF ( BSTAT(IB0) .EQ. 0 ) THEN
!/MPIT                WRITE (STRT(1:2),'(I2)') BSTAT(IB0)
!/MPIT              ELSE
!/MPIT                WRITE (STRT(1:7),'(I2,I5)') BSTAT(IB0), BISPL(IB0)
!/MPIT              END IF
!/MPIT            STR(IB0) = STRT
!/MPIT          END IF
!/MPIT        END DO
!
!/MPIT      WRITE (NDST,9000) STR
!
!/MPIT      IF ( ISPLOC .EQ. 0 ) THEN
!/MPIT          DO IB0=1, MPIBUF
!/MPIT            STR(IB0) = '--------------+'
!/MPIT            END DO
!/MPIT          WRITE (NDST,9010) STR
!/MPIT          WRITE (NDST,*)
!/MPIT        END IF
!
!/MPI      RETURN
!
! Formats
!
!/MPIT 9000 FORMAT ( ' |     |',8A15)
!/MPIT 9010 FORMAT ( ' +-----+',8A15)
!/
!/ End of W3SCAT ----------------------------------------------------- /
!/
      END SUBROUTINE W3SCAT
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3NMIN ( MAPSTA, FLAG0 )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         28-Dec-2004 |
!/                  +-----------------------------------+
!/
!/    23-Feb-2001 : Origination.                        ( version 2.07 )
!/    28-Dec-2004 : Multiple grid version.              ( version 3.06 )
!/
!  1. Purpose :
!
!     Check minimum number of active sea points at given processor to
!     evaluate the need for a MPI_BARRIER call.
!
!  2. Method :
!
!     Evaluate mapsta.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       MAPSTA  I.A.   I   Status map for spatial grid.
!       FLAG0   log.   O   Flag to identify 0 as minimum.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3WAVE    Subr. W3WAVEMD Actual wave model routine.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!     None.
!
!  7. Remarks :
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/S     Enable subroutine tracing.
!     !/T     Test output.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/S      USE W3SERVMD, ONLY: STRACE
!/
      USE W3GDATMD, ONLY: NX, NY, NSEA, MAPSF
      USE W3ODATMD, ONLY: NDST, NAPROC
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: MAPSTA(NY*NX)
      LOGICAL, INTENT(OUT)    :: FLAG0
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: NMIN, IPROC, NLOC, ISEA, IXY
!/S      INTEGER, SAVE           :: IENT
!/
!/ ------------------------------------------------------------------- /
!/
!/S      CALL STRACE (IENT, 'W3NMIN')
!
      NMIN   = NSEA
!
      DO IPROC=1, NAPROC
        NLOC   = 0
        DO ISEA=IPROC, NSEA, NAPROC
          IXY    = MAPSF(ISEA,3)
          IF ( MAPSTA(IXY) .EQ. 1 ) NLOC = NLOC + 1
          END DO
!/T        WRITE (NDST,9000) IPROC, NLOC
        NMIN   = MIN ( NMIN , NLOC )
        END DO
!
      FLAG0  = NMIN .EQ. 0
!/T      WRITE (NDST,9001) NMIN, FLAG0
!
      RETURN
!
! Formats
!
!/T 9000 FORMAT ( ' TEST W3NMIN : IPROC =',I3,'  NLOC =',I5)
!/T 9001 FORMAT ( ' TEST W3NMIN : NMIN =',I5,'  FLAG0 =',L2)
!/
!/ End of W3NMIN ----------------------------------------------------- /
!/
      END SUBROUTINE W3NMIN
!/
!/ End of module W3WAVEMD -------------------------------------------- /
!/
      END MODULE W3WAVEMD
