!/ ------------------------------------------------------------------- /
      PROGRAM W3OUTF
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-May-2009 |
!/                  +-----------------------------------+
!/
!/    19-Oct-1998 : Final FORTRAN 77                    ( version 1.18 )
!/    19-Jan-2000 : Upgrade to FORTRAN 90               ( version 2.00 )
!/    24-Jan-2001 : Flat grid version                   ( version 2.06 )
!/    23-Apr-2002 : Clean-up                            ( version 2.19 )
!/    29-Apr-2002 : Adding output fields 17-18.         ( version 2.20 )
!/    13-Nov-2002 : Add stress vector.                  ( version 3.00 )
!/    24-Dec-2004 : Multiple grid version.              ( version 3.06 )
!/    21-Jul-2005 : Adding output fields 19-21.         ( version 3.07 )
!/    28-Jun-2006 : Adding file name preamble.          ( version 3.09 )
!/    05-Jul-2006 : Consolidate stress arrays.          ( version 3.09 )
!/    28-Mar-2007 : Adding partitioned output.          ( version 3.11 )
!/                  Adding user slots for outputs.
!/    31-Jul-2007 : Fix file extension errors.          ( version 3.12 )
!/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
!/
!/    Copyright 2009 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS. 
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     Post-processing of grid output.
!
!  2. Method :
!
!     Data is read from the grid output file out_grd.ww3 (raw data)
!     and from the file ww3_outf.inp ( NDSI, output requests ).
!     Model definition and raw data files are read using WAVEWATCH III
!     subroutines.
!
!     Output types :
!      1 : print plots
!      2 : field statistics
!      3 : transfer file
!
!  3. Parameters :
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3NMOD    Subr. W3GDATMD Set number of model.
!      W3SETG    Subr.   Id.    Point to selected model.
!      W3NDAT    Subr. W3WDATMD Set number of model for wave data.
!      W3SETW    Subr.   Id.    Point to selected model for wave data.
!      W2NAUX    Subr. W3ADATMD Set number of model for aux data.
!      W3SETA    Subr.   Id.    Point to selected model for aux data.
!      ITRACE    Subr. W3SERVMD Subroutine tracing initialization.
!      STRACE    Subr.   Id.    Subroutine tracing.
!      NEXTLN    Subr.   Id.    Get next line from input filw
!      EXTCDE    Subr.   Id.    Abort program as graceful as possible.
!      STME21    Subr. W3TIMEMD Convert time to string.
!      TICK21    Subr.   Id.    Advance time.
!      DSEC21    Func.   Id.    Difference between times.
!      W3IOGR    Subr. W3IOGRMD Reading/writing model definition file.
!      W3IOGO    Subr. W3IOGOMD Reading/writing raw gridded data file.
!      W3EXGO    Subr. Internal Execute grid output.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!     None, stand-alone program.
!
!  6. Error messages :
!
!     Checks on input, checks in W3IOxx.
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
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE CONSTANTS
!/
!     USE W3GDATMD, ONLY: W3NMOD, W3SETG
      USE W3WDATMD, ONLY: W3NDAT, W3SETW
      USE W3ADATMD, ONLY: W3NAUX, W3SETA
      USE W3ODATMD, ONLY: W3NOUT, W3SETO
      USE W3SERVMD, ONLY : ITRACE, NEXTLN, EXTCDE
!/S      USE W3SERVMD, ONLY : STRACE
      USE W3TIMEMD
      USE W3IOGRMD, ONLY: W3IOGR
      USE W3IOGOMD, ONLY: W3IOGO
!/
      USE W3GDATMD
      USE W3WDATMD, ONLY: TIME, WLV, ICE, UST, USTDIR
      USE W3ADATMD, ONLY: DW, UA, UD, AS, CX, CY, HS, WLM, TMN, THM,  &
                          THS, FP0, THP0, FP1, THP1, DTDYN, FCUT,     &
                          ABA, ABD, UBA, UBD, SXX, SYY, SXY, USERO,   &
                          PHS, PTP, PLP, PTH, PSI, PWS, PWST, PNR
      USE W3ODATMD, ONLY: NDSO, NDSE, NDST, NOGRD, IDOUT, UNDEF,      &
                          FLOGRD, FNMPRE, NOSWLL
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: NDSI, NDSM, NDSOG, NDSDAT, NDSDT,    &
                                 NDSTRC, NTRACE, IERR, I, J,          &
                                 TOUT(2), TDUM(2), IOTEST, NOUT,      &
                                 ITYPE, IX1, IXN, IXS, IY1, IYN, IYS, &
                                 IDLA, IDFM, IOUT, IPART
!/S      INTEGER, SAVE           :: IENT = 0
      REAL                    :: DTREQ, DTEST
      CHARACTER               :: COMSTR*1, IDTIME*23, IDDDAY*11,      &
                                 TABNME*9
      LOGICAL                 :: FLREQ(NOGRD), SCALE, VECTOR
!/
!/ ------------------------------------------------------------------- /
!/
! 1.  IO set-up.
!
      CALL W3NMOD ( 1, 6, 6 )
      CALL W3SETG ( 1, 6, 6 )
      CALL W3NDAT (    6, 6 )
      CALL W3SETW ( 1, 6, 6 )
      CALL W3NAUX (    6, 6 )
      CALL W3SETA ( 1, 6, 6 )
      CALL W3NOUT (    6, 6 )
      CALL W3SETO ( 1, 6, 6 )
!
      NDSI   = 10
      NDSM   = 20
      NDSOG  = 20
      NDSDAT = 50
!
      NDSTRC =  6
      NTRACE = 10
      CALL ITRACE ( NDSTRC, NTRACE )
!
!/S      CALL STRACE (IENT, 'W3OUTF')
!
      WRITE (NDSO,900)
!
      J      = LEN_TRIM(FNMPRE)
      OPEN (NDSI,FILE=FNMPRE(:J)//'ww3_outf.inp',STATUS='OLD',       &
            ERR=800,IOSTAT=IERR)
      READ (NDSI,'(A)',END=801,ERR=802) COMSTR
      IF (COMSTR.EQ.' ') COMSTR = '$'
      WRITE (NDSO,901) COMSTR
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 2.  Read model definition file.
!
      CALL W3IOGR ( 'READ', NDSM )
      WRITE (NDSO,920) GNAME
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 3.  Read general data and first fields from file
!
      CALL W3IOGO ( 'READ', NDSOG, IOTEST )
!
      WRITE (NDSO,930)
      DO I=1, NOGRD
        IF ( FLOGRD(I) ) WRITE (NDSO,931) IDOUT(I)
        END DO
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 4.  Read requests from input file.
!     Output times
!
      CALL NEXTLN ( COMSTR , NDSI , NDSE )
      READ (NDSI,*,END=801,ERR=802) TOUT, DTREQ, NOUT
      DTREQ  = MAX ( 0. , DTREQ )
      IF ( DTREQ.EQ.0. ) NOUT = 1
      NOUT   = MAX ( 1 , NOUT )
!
      CALL STME21 ( TOUT , IDTIME )
      WRITE (NDSO,940) IDTIME
!
      TDUM = 0
      CALL TICK21 ( TDUM , DTREQ )
      CALL STME21 ( TDUM , IDTIME )
      IF ( DTREQ .GE. 86400. ) THEN
          WRITE (IDDDAY,'(I10,1X)') INT(DTREQ/86400.)
        ELSE
          IDDDAY = '           '
        END IF
      IDTIME(1:11) = IDDDAY
      IDTIME(21:23) = '   '
      WRITE (NDSO,941) IDTIME, NOUT
!
! ... Output fields
!
      CALL NEXTLN ( COMSTR , NDSI , NDSE )
      READ (NDSI,*,END=801,ERR=802) FLREQ
!
! ... Output type
!
      CALL NEXTLN ( COMSTR , NDSI , NDSE )
      READ (NDSI,*,END=801,ERR=802) ITYPE, IPART
      IF ( ITYPE.LT.0 .OR. ITYPE.GT.3 ) THEN
          WRITE (NDSE,1010) ITYPE
          CALL EXTCDE ( 1 )
        END IF
      IPART  = MAX ( 0 , MIN ( NOSWLL , IPART ) )
!
! ... ITYPE = 0
!
      IF ( ITYPE .EQ. 0 ) THEN
          WRITE (NDSO,942) ITYPE, 'Checking contents of file'
          DO
            CALL STME21 ( TIME , IDTIME )
            WRITE (NDSO,943) IDTIME
            CALL W3IOGO ( 'READ', NDSOG, IOTEST )
            IF ( IOTEST .EQ. -1 ) THEN
                WRITE (NDSO,944)
                GOTO 888
              END IF
            END DO
!
! ... ITYPE = 1
!
        ELSE IF (ITYPE .EQ. 1) THEN
          WRITE (NDSO,942) ITYPE, 'Print plots'
          CALL NEXTLN ( COMSTR , NDSI , NDSE )
          READ (NDSI,*,END=801,ERR=802)                               &
                IX1, IXN, IXS, IY1, IYN, IYS, SCALE, VECTOR
          IX1    = MAX ( IX1 , 1 )
          IXN    = MIN ( IXN , NX )
          IXS    = MAX ( IXS , 1 )
          IY1    = MAX ( IY1 , 1 )
          IYN    = MIN ( IYN , NY )
          IYS    = MAX ( IYS , 1 )
          WRITE (NDSO,1940) IX1, IXN, IXS, IY1, IYN, IYS
          IF ( SCALE ) WRITE (NDSO,1941)
!
! ... ITYPE = 2
!
        ELSE IF (ITYPE .EQ. 2) THEN
          WRITE (NDSO,942) ITYPE, 'Field statistics'
          NDSDT  = NDSDAT - 1
          CALL NEXTLN ( COMSTR , NDSI , NDSE )
          READ (NDSI,*,END=801,ERR=802) IX1, IXN, IY1, IYN
          IX1    = MAX ( IX1 , 1 )
          IXN    = MIN ( IXN , NX )
          IY1    = MAX ( IY1 , 1 )
          IYN    = MIN ( IYN , NY )
          WRITE (NDSO,2940) IX1, IXN, IY1, IYN
!
! ... ITYPE = 3
!
        ELSE IF (ITYPE .EQ. 3) THEN
          WRITE (NDSO,942) ITYPE, 'Transfer files'
          CALL NEXTLN ( COMSTR , NDSI , NDSE )
          READ (NDSI,*,END=801,ERR=802)                               &
                IX1, IXN, IY1, IYN, IDLA, IDFM
          IX1    = MAX ( IX1 , 1 )
          IXN    = MIN ( IXN , NX )
          IY1    = MAX ( IY1 , 1 )
          IYN    = MIN ( IYN , NY )
          IF (IDLA.LT.1 .OR. IDLA.GT.5) IDLA   = 1
          IF (IDFM.LT.1 .OR. IDFM.GT.3) IDFM   = 1
          VECTOR = .TRUE.
          WRITE (NDSO,3940) IX1, IXN, IY1, IYN, IDLA, IDFM
!
        END IF
!
! ... Output of output fields
!
      IF ( ITYPE.NE.2 ) THEN
          WRITE (NDSO,945)
        ELSE
          WRITE (NDSO,2945)
        END IF
!
      DO I=1, NOGRD
        IF ( FLREQ(I) ) THEN
            IF ( FLOGRD(I) ) THEN
                IF ( ITYPE.NE.2 ) THEN
                    WRITE (NDSO,946) IDOUT(I), ' '
                  ELSE
                    J      = LEN_TRIM(FNMPRE)
                    NDSDT  = NDSDT + 1
                    WRITE (TABNME,'(A3,I2.2,A4)') 'tab', NDSDT, '.ww3'
                    WRITE (NDSO,2946) TABNME, IDOUT(I)
                    OPEN (NDSDT,FILE=FNMPRE(:J)//TABNME)
                    WRITE (NDSDT,2947) IDOUT(I)
                  END IF
              ELSE
                WRITE (NDSO,946) IDOUT(I), '*** NOT AVAILABLE ***'
                FLREQ(I) = .FALSE.
              END IF
          END IF
        END DO
!
      IF ( IPART .EQ. 0 ) THEN
          WRITE (NDSO,948)
        ELSE
          WRITE (NDSO,949) IPART
        END IF
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 5.  Time management.
!
      IOUT   = 0
      IF (ITYPE.EQ.3) WRITE (NDSO,970)
!
      DO
        DTEST  = DSEC21 ( TIME , TOUT )
        IF ( DTEST .GT. 0. ) THEN
            CALL W3IOGO ( 'READ', NDSOG, IOTEST )
              IF ( IOTEST .EQ. -1 ) THEN
                WRITE (NDSO,944)
                GOTO 888
              END IF
            CYCLE
          END IF
        IF ( DTEST .LT. 0. ) THEN
            CALL TICK21 ( TOUT , DTREQ )
            CYCLE
          END IF
!
        IOUT   = IOUT + 1
        CALL STME21 ( TOUT , IDTIME )
        IF (ITYPE.EQ.1) THEN
            WRITE (NDSO,950) IDTIME
          ELSE IF (ITYPE.EQ.3) THEN
            WRITE (NDSO,971) IDTIME
          END IF
!
        CALL W3EXGO ( NX, NY, NSEA )
        CALL TICK21 ( TOUT , DTREQ )
        IF ( IOUT .GE. NOUT ) EXIT
        END DO
!
      IF (ITYPE.EQ.3) WRITE (NDSO,972)
!
      GOTO 888
!
! Escape locations read errors :
!
  800 CONTINUE
      WRITE (NDSE,1000) IERR
      CALL EXTCDE ( 10 )
!
  801 CONTINUE
      WRITE (NDSE,1001)
      CALL EXTCDE ( 11 )
!
  802 CONTINUE
      WRITE (NDSE,1002) IERR
      CALL EXTCDE ( 12 )
!
  888 CONTINUE
      WRITE (NDSO,999)
!
! Formats
!
  900 FORMAT (/15X,'   *** WAVEWATCH III Field output postp. ***   '/ &
               15X,'==============================================='/)
  901 FORMAT ( '  Comment character is ''',A,''''/)
!
  920 FORMAT ( '  Grid name : ',A/)
!
  930 FORMAT ( '  Fields in file : '/                                 &
               ' --------------------------')
  931 FORMAT ( '      ',A)
!
  940 FORMAT (/'  Output time data : '/                               &
               ' --------------------------------------------------'/ &
               '      First time         : ',A)
  941 FORMAT ( '      Interval           : ',A/                       &
               '      Number of requests : ',I4)
  942 FORMAT (/'  Output type ',I2,' :'/                              &
               ' --------------------------------------------------'/ &
               '      ',A/)
  943 FORMAT ( '      Data for ',A)
  944 FORMAT (/'      End of file reached '/)
!
  945 FORMAT (/'  Requested output fields : '/                        &
               ' --------------------------------------------------')
 2945 FORMAT (/'  Output files and fields : '/                        &
               ' --------------------------------------------------')
  946 FORMAT ( '      ',A,2X,A)
 2946 FORMAT ( '      ',A,' : ',A)
 2947 FORMAT ( ' Statitics of ',A/                                    &
               '   (time, min, max, avg, std)'/)
  948 FORMAT (/'         Partitioned field data for wind seas')
  949 FORMAT (/'         Partitioned field data for swell field',I2)
!
 1940 FORMAT ( '      X range and interval : ',3I5/                   &
               '      Y range and interval : ',3I5)
 1941 FORMAT ( '      Data is normalized ')
!
 2940 FORMAT ( '      X range : ',2I5/                                &
               '      Y range : ',2I5)
!
 3940 FORMAT ( '      X range          : ',2I5/                       &
               '      Y range          : ',2I5/                       &
               '      Layout indicator : ',I5/                        &
               '      Format indicator : ',I5)
!
  950 FORMAT (//'  Output for ',A/                                    &
               ' --------------------------------------------------')
!
  970 FORMAT (//'  Generating files '/                                &
               ' --------------------------------------------------')
  971 FORMAT ( '      Files for ',A)
  972 FORMAT ( ' ')
!
  999 FORMAT (/'  End of program '/                                   &
               ' ========================================='/          &
               '         WAVEWATCH III Field output '/)
!
 1000 FORMAT (/' *** WAVEWATCH III ERROR IN W3OUTF : '/               &
               '     ERROR IN OPENING INPUT FILE'/                    &
               '     IOSTAT =',I5/)
!
 1001 FORMAT (/' *** WAVEWATCH III ERROR IN W3OUTF : '/               &
               '     PREMATURE END OF INPUT FILE'/)
!
 1002 FORMAT (/' *** WAVEWATCH III ERROR IN W3OUTF : '/               &
               '     ERROR IN READING FROM INPUT FILE'/               &
               '     IOSTAT =',I5/)
!
 1010 FORMAT (/' *** WAVEWATCH III ERROR IN W3OUTF : '/               &
               '     ILLEGAL TYPE, ITYPE =',I4/)
!/
!/ Internal subroutine W3EXGO ---------------------------------------- /
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3EXGO ( NX, NY, NSEA )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         31-Jul-2007 |
!/                  +-----------------------------------+
!/
!/    26-Sep-1997 : Final FORTRAN 77                    ( version 1.18 )
!/    19-Jan-2000 : Upgrade to FORTRAN 90               ( version 2.00 )
!/                  Massive changes to logistics
!/    24-Jan-2001 : Flat grid version                   ( version 2.06 )
!/    23-Apr-2002 : Clean-up                            ( version 2.19 )
!/    29-Apr-2002 : Adding output fields 17-18.         ( version 2.20 )
!/    16-Oct-2002 : Fix bound. error for stress output. ( version 3.00 )
!/    16-Oct-2002 : Fix statistical output for UNDEF.   ( version 3.00 )
!/    13-Nov-2002 : Add stress vector.                  ( version 3.00 )
!/    24-Dec-2004 : Multiple grid version.              ( version 3.06 )
!/    21-Jul-2005 : Adding output fields 19-21.         ( version 3.07 )
!/    28-Jun-2006 : Adding file name preamble.          ( version 3.09 )
!/    05-Jul-2006 : Consolidate stress arrays.          ( version 3.09 )
!/    28-Mar-2007 : Adding partitioned output.          ( version 3.11 )
!/                  Adding user slots for outputs.
!/    31-Jul-2007 : Fix file extension errors.          ( version 3.12 )
!/
!  1. Purpose :
!
!     Perform actual grid output.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       NX/Y    Int.  I  Grid dimensions.
!       NSEA    Int.  I  Number of sea points.
!     ----------------------------------------------------------------
!
!     Internal parameters
!     ----------------------------------------------------------------
!       FLONE   Log.  Flags for one or two-dimensional field.
!       X1, X2, XX, XY
!               R.A.  Output fields
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!      EXTCDE    Subr.   Id.    Abort program as graceful as possible.
!      W3S2XY    Subr.   Id.    Convert from storage to spatial grid.
!      PRTBLK    Subr. W3ARRYMD Print plot of array.
!      OUTA2I    Subr.   Id.    Print array of INTEGERS.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!     Main program in which it is contained.
!
!  6. Error messages :
!
!       None.
!
!  7. Remarks :
!
!     - Note that arrays CX and CY of the main program now contain
!       the absolute current speed and direction respectively.
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!       !/LLG  Spherical grid.
!       !/XYG  Spherical grid.
!
!       !/S  Enable subroutine tracing.
!       !/T  Enable test output.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3SERVMD, ONLY : W3S2XY
      USE W3ARRYMD, ONLY : OUTA2I, PRTBLK
      USE OutNetCDF                      ![GIR] write outputs in netCDF
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER                 :: NX, NY, NSEA
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: NXMAX, NXTOT, NBLOK, IH, IM, IS,     &
                                 MFILL, J, ISEA, IX, IY, IXB, IB,     &
                                 IXA, NINGRD, JJ
      INTEGER                 :: MAP(NX+1,NY), MP2(NX+1,NY),          &
                                 MX1(NX,NY), MXX(NX,NY), MYY(NX,NY)
      INTEGER, SAVE           :: IPASS
!     INTEGER, SAVE           :: NCOL   = 80
      INTEGER, SAVE           :: NCOL   = 132
!/S      INTEGER, SAVE           :: IENT   =   0
      REAL                    :: FSC, CABS, UABS, FSCA, XMIN, XMAX,   &
                                 XAVG, XSTD, YGBX, XGBX, AABS
      REAL                    :: X1(NX+1,NY), X2(NX+1,NY),            &
                                 XX(NX+1,NY), XY(NX+1,NY) 
      DOUBLE PRECISION        :: XDS, XDSQ
      LOGICAL                 :: FLONE
      CHARACTER               :: OLDTID*8, FNAME*16, ENAME*4,         &
                                 FORMF*11, UNITS*10
      CHARACTER, SAVE         :: TIMEID*8 = '00000000'
      CHARACTER, SAVE         :: FILEID*13 = 'WAVEWATCH III'
!/
!/ ------------------------------------------------------------------- /
!/
!
!/S      CALL STRACE (IENT, 'W3EXGO')
!
!/T      WRITE (NDST,9000) (FLREQ(J),J=1,NOGRD)
!/T      WRITE (NDST,9001) ITYPE, IX1, IXN, IXS, IY1, IYN, IYS,          &
!/T                        SCALE, VECTOR, NDSDAT
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 1.  Preparations
!
      X1     = UNDEF
      X2     = UNDEF
      XX     = UNDEF
      XY     = UNDEF
!
!     Number of print-plots
!
      IF ( ITYPE .EQ. 1 ) THEN
          IF ( SCALE ) THEN
              NXMAX  = ( NCOL - 10 ) / 2
            ELSE
              NXMAX  = ( NCOL - 10 ) / 5
            END IF
          NXTOT  = 1 + (IXN-IX1)/IXS
          NBLOK  = 1 + (NXTOT-1)/NXMAX
!/T          WRITE (NDST,9012) NXMAX, NXTOT, NBLOK
        END IF
!
!     Output file unit number
!
      IF ( ITYPE .EQ. 2 ) THEN
          NDSDT = NDSDAT - 1
          IH    = TIME(2) / 10000
          IM    = MOD ( TIME(2) , 10000 ) / 100
          IS    = MOD ( TIME(2) , 100 )
        END IF
!
!     Set-up transfer files
!
      IF ( ITYPE .EQ. 3 ) THEN
          MFILL  = -999
          OLDTID = TIMEID
          WRITE (TIMEID,'(I6.6,I2.2)') MOD( TIME(1) , 1000000 ),      &
                                       TIME(2)/10000
          FNAME(05:12) = TIMEID
          FNAME(13:13) = '.'
          IF ( TIMEID .NE. OLDTID ) THEN
              FNAME(1:4) = 'ww3.'
              IPASS    = 1
            ELSE
              WRITE (ENAME,'(A1,I2.2,A1)') 'e', IPASS, '.'
              FNAME(1:4) = ENAME
              IPASS   = IPASS + 1
            END IF
!/T          WRITE (NDST,9014) FNAME(1:13)
        END IF
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 2.  Loop over output fields.
!

      DO J=1, NOGRD
        IF ( FLREQ(J) ) THEN
!
            FORMF  = '(1X,32I4)'
!/T            WRITE (NDST,9020) IDOUT(J)
!
! 2.a Set output arrays and parameters
!
            IF ( J .EQ.  1 ) THEN
                FLONE  = .TRUE.
                FSC    = 1.
                IF ( ITYPE .EQ. 3 ) FSC = 0.1
                UNITS  = 'm'
                ENAME  = '.dpt'
                FORMF  = '(1X,20I6)'
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, DW(1:NSEA)        &
                                                          , MAPSF, X1 )
!
              ELSE IF ( J .EQ.  2 ) THEN
                FLONE  = .FALSE.
                FSC    = 0.01
                ENAME  = '.cur'
                UNITS  = 'm/s'
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, CX(1:NSEA)        &
                                                          , MAPSF, XX )
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, CY(1:NSEA)        &
                                                          , MAPSF, XY )
                DO ISEA=1, NSEA
                  CABS   = SQRT(CX(ISEA)**2+CY(ISEA)**2)
                  IF ( CABS .GT. 0.05 ) THEN
                      CY(ISEA) = MOD ( 630. -                         &
                            RADE*ATAN2(CY(ISEA),CX(ISEA)) , 360. )
                    ELSE
                      CY(ISEA) = UNDEF
                    END IF
                  CX(ISEA) = CABS
                  END DO
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, CX(1:NSEA)        &
                                                          , MAPSF, X1 )
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, CY(1:NSEA)        &
                                                          , MAPSF, X2 )
!
              ELSE IF ( J .EQ.  3 ) THEN
                FLONE  = .FALSE.
                FSC    = 0.1
                ENAME  = '.wnd'
                UNITS  = 'm/s'
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, UA(1:NSEA)        &
                                                          , MAPSF, XX )
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, UD(1:NSEA)        &
                                                          , MAPSF, XY )
                DO ISEA=1, NSEA
                  UABS   = SQRT(UA(ISEA)**2+UD(ISEA)**2)
                  IF ( UABS .GT. 1.0 ) THEN
                      UD(ISEA) = MOD ( 630. -                         &
                            RADE*ATAN2(UD(ISEA),UA(ISEA)) , 360. )
                    ELSE
                      UD(ISEA) = UNDEF
                    END IF
                  UA(ISEA) = UABS
                  END DO
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, UA(1:NSEA)        &
                                                          , MAPSF, X1 )
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, UD(1:NSEA)        &
                                                          , MAPSF, X2 )
!
              ELSE IF ( J .EQ.  4 ) THEN
                FLONE  = .TRUE.
                FSC    = 0.1
                ENAME  = '.dt'
                UNITS  = 'deg'
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, AS(1:NSEA)        &
                                                          , MAPSF, X1 )
!
              ELSE IF ( J .EQ.  5 ) THEN
                FLONE  = .FALSE.
                FSC    = 0.001
                ENAME  = '.ust'
                FORMF  = '(1X,20I6)'
                UNITS  = 'm/s'
                CALL W3S2XY (NSEA,NSEA,NX+1,NY, UST   (1:NSEA)        &
                                                      , MAPSF, XX )
                CALL W3S2XY (NSEA,NSEA,NX+1,NY, USTDIR(1:NSEA)        &
                                                      , MAPSF, XY )
                DO ISEA=1, NSEA
                  UABS   = SQRT(UST(ISEA)**2+USTDIR(ISEA)**2)
                  IF ( UST(ISEA) .EQ. UNDEF ) THEN
                      USTDIR(ISEA) = UNDEF
                      UABS         = UNDEF
                    ELSE IF ( UABS .GT. 0.05 ) THEN
                      USTDIR(ISEA) = MOD ( 630. -                     &
                        RADE*ATAN2(USTDIR(ISEA),UST(ISEA)) , 360. )
                    ELSE
                      USTDIR(ISEA) = UNDEF
                    END IF
                  UST(ISEA) = UABS
                  END DO
                CALL W3S2XY (NSEA,NSEA,NX+1,NY, UST   (1:NSEA)        &
                                                      , MAPSF, X1 )
                CALL W3S2XY (NSEA,NSEA,NX+1,NY, USTDIR(1:NSEA)        &
                                                      , MAPSF, X2 )
!
              ELSE IF ( J .EQ.  6 ) THEN
                FLONE  = .TRUE.
                FSC    = 0.01
                UNITS  = 'm'
                ENAME  = '.hs'
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, HS    , MAPSF, X1 )
!
              ELSE IF ( J .EQ.  7 ) THEN
                FLONE  = .TRUE.
                FSC    = 1.
                UNITS  = 'm'
                ENAME  = '.l'
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, WLM   , MAPSF, X1 )
!
              ELSE IF ( J .EQ.  8 ) THEN
                FLONE  = .TRUE.
                FSC    = 0.01
                UNITS  = 's'
                ENAME  = '.t'
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, TMN   , MAPSF, X1 )
!
              ELSE IF ( J .EQ.  9 ) THEN
                FLONE  = .TRUE.
                FSC    = 1.
                UNITS  = 'degr.'
                ENAME  = '.dir'
                DO ISEA=1, NSEA
                  IF ( THM(ISEA) .NE. UNDEF )                         &
                       THM(ISEA) = MOD ( 630. - RADE*THM(ISEA) , 360. )
                  IF ( THM(ISEA) .NE. UNDEF )                         &
                        THM(ISEA) = MOD (THM(ISEA)+180, 360.)
                  END DO
                
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, THM   , MAPSF, X1 )
!
              ELSE IF ( J .EQ. 10 ) THEN
                FLONE  = .TRUE.
                FSC    = 0.1
                UNITS  = 'degr.'
                ENAME  = '.spr'
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, THS   , MAPSF, X1 )
!
              ELSE IF ( J .EQ. 11 ) THEN
                FLONE  = .TRUE.
                FSC    = 0.001
                UNITS  = 'Hz'
                ENAME  = '.fp'
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, FP0   , MAPSF, X1 )
!
              ELSE IF ( J .EQ. 12 ) THEN
                FLONE  = .TRUE.
                FSC    = 1.
                UNITS  = 'degr.'
                ENAME  = '.dp'
                DO ISEA=1, NSEA
                  IF ( THP0(ISEA) .NE. UNDEF ) THEN
                      THP0(ISEA) = MOD ( 630-RADE*THP0(ISEA) , 360. )
                    END IF
                  END DO
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, THP0  , MAPSF, X1 )
!
              ELSE IF ( J .EQ. 13 ) THEN
                FLONE  = .TRUE.
                FSC    = 0.001
                UNITS  = 'Hz'
                ENAME  = '.fpl'
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, FP1   , MAPSF, X1 )
!
              ELSE IF ( J .EQ. 14 ) THEN
                FLONE  = .TRUE.
                FSC    = 1.
                UNITS  = 'degr.'
                ENAME  = '.dpl'
                DO ISEA=1, NSEA
                  IF ( THP1(ISEA) .NE. UNDEF ) THEN
                      THP1(ISEA) = MOD ( 630-RADE*THP1(ISEA) , 360. )
                    END IF
                  END DO
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, THP1  , MAPSF, X1 )
!
              ELSE IF ( J .EQ. 15 ) THEN
                FLONE  = .TRUE.
                FSC    = 0.01
                UNITS  = 'm'
                ENAME  = '.phs'
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, PHS(:,IPART)      &
                                                      , MAPSF, X1 )
!
              ELSE IF ( J .EQ. 16 ) THEN
                FLONE  = .TRUE.
                FSC    = 0.01
                UNITS  = 's'
                ENAME  = '.ptp'
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, PTP(:,IPART)      &
                                                      , MAPSF, X1 )
!
              ELSE IF ( J .EQ. 17 ) THEN
                FLONE  = .TRUE.
                FSC    = 1.
                UNITS  = 'm'
                ENAME  = '.plp'
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, PLP(:,IPART)      &
                                                      , MAPSF, X1 )
!
              ELSE IF ( J .EQ. 18 ) THEN
                FLONE  = .TRUE.
                FSC    = 1.
                UNITS  = 'degr.'
                ENAME  = '.pth'
                DO ISEA=1, NSEA
                  IF ( PTH(ISEA,IPART) .NE. UNDEF ) THEN
                      PTH(ISEA,IPART) =                               &
                            MOD ( 630-RADE*PTH(ISEA,IPART) , 360. )
                    END IF
                  END DO
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, PTH(:,IPART)      &
                                                      , MAPSF, X1 )
!
              ELSE IF ( J .EQ. 19 ) THEN
                FLONE  = .TRUE.
                FSC    = 0.1
                UNITS  = 'degr.'
                ENAME  = '.psi'
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, PSI(:,IPART)      &
                                                      , MAPSF, X1 )
!
              ELSE IF ( J .EQ. 20 ) THEN
                FLONE  = .TRUE.
                FSC    = 0.001
                UNITS  = ' '
                ENAME  = '.pws'
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, PWS(:,IPART)      &
                                                      , MAPSF, X1 )
!
              ELSE IF ( J .EQ. 21 ) THEN
                FLONE  = .TRUE.
                FSC    = 0.001
                UNITS  = ' '
                ENAME  = '.wsf'
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, PWST(:), MAPSF, X1 )
!
              ELSE IF ( J .EQ. 22 ) THEN
                FLONE  = .TRUE.
                FSC    = 1.
                UNITS  = ' '
                ENAME  = '.pnr'
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, PNR(:), MAPSF, X1 )
!
              ELSE IF ( J .EQ. 23 ) THEN
                FLONE  = .TRUE.
                FSC    = 0.1
                UNITS  = 'min.'
                ENAME  = '.dtd'
                DO ISEA=1, NSEA
                  IF ( DTDYN(ISEA) .NE. UNDEF )                       &
                       DTDYN(ISEA) = DTDYN(ISEA) / 60.
                  END DO
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, DTDYN , MAPSF, X1 )
!
              ELSE IF ( J .EQ. 24 ) THEN
                FLONE  = .TRUE.
                FSC    = 0.001
                UNITS  = 'Hz'
                ENAME  = '.fc'
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, FCUT  , MAPSF, X1 )
!
              ELSE IF ( J .EQ. 25 ) THEN
                FLONE  = .TRUE.
                FSC    = 0.001
                UNITS  = ' '
                ENAME  = '.ice'
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, ICE   , MAPSF, X1 )
!
              ELSE IF ( J .EQ. 26 ) THEN
                FLONE  = .TRUE.
                FSC    = 0.01
                UNITS  = 'm'
                ENAME  = '.wlv'
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, WLV   , MAPSF, X1 )
!
              ELSE IF ( J .EQ. 27 ) THEN
                FLONE  = .FALSE.
                FSC    = 0.01
                ENAME  = '.abr'
                UNITS  = 'm'
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, ABA(1:NSEA)       &
                                                          , MAPSF, XX )
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, ABD(1:NSEA)       &
                                                          , MAPSF, XY )
                DO ISEA=1, NSEA
                  IF ( ABA(ISEA) .NE. UNDEF ) THEN
                      AABS   = SQRT(ABA(ISEA)**2+ABD(ISEA)**2)
                      IF ( AABS .GT. 0.005 ) THEN
                          ABD(ISEA) = MOD ( 630. -                    &
                                RADE*ATAN2(ABD(ISEA),ABA(ISEA)) , 360. )
                        ELSE
                          ABD(ISEA) = UNDEF
                        END IF
                      ABA(ISEA) = AABS
                    END IF
                  END DO
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, ABA(1:NSEA)       &
                                                          , MAPSF, X1 )
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, ABD(1:NSEA)       &
                                                          , MAPSF, X2 )

!
              ELSE IF ( J .EQ. 28 ) THEN
                FLONE  = .FALSE.
                FSC    = 0.01
                ENAME  = '.ubr'
                UNITS  = 'm/s'
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, UBA(1:NSEA)       &
                                                          , MAPSF, XX )
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, UBD(1:NSEA)       &
                                                          , MAPSF, XY )
                DO ISEA=1, NSEA
                  IF ( UBA(ISEA) .NE. UNDEF ) THEN
                      UABS   = SQRT(UBA(ISEA)**2+UBD(ISEA)**2)
                      IF ( UABS .GT. 0.005 ) THEN
                          UBD(ISEA) = MOD ( 630. -                    &
                                RADE*ATAN2(UBD(ISEA),UBA(ISEA)) , 360. )
                        ELSE
                          UBD(ISEA) = UNDEF
                        END IF
                      UBA(ISEA) = UABS
                    END IF
                  END DO
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, UBA(1:NSEA)       &
                                                          , MAPSF, X1 )
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, UBD(1:NSEA)       &
                                                          , MAPSF, X2 )
!
              ELSE IF ( J .EQ. 29 ) THEN
                FLONE  = .TRUE.
                FSC    = 10.
                UNITS  = 'N/m'
                ENAME  = '.Sxy'
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, SXX(1:NSEA)       &
                                                          , MAPSF, X1 )
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, SYY(1:NSEA)       &
                                                          , MAPSF, X2 )
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, SXY(1:NSEA)       &
                                                          , MAPSF, XY )
!
              ELSE IF ( J .EQ. 30 ) THEN
                FLONE  = .TRUE.
                FSC    = 1
                UNITS  = '-'
                ENAME  = '.Cdg'
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, USERO(:,1)        &
                                                          , MAPSF, X1 )
!
              ELSE IF ( J .EQ. 31 ) THEN
                FLONE  = .TRUE.
                FSC    = 1
                UNITS  = 'm/s'
                ENAME  = '.SDx'
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, USERO(:,2)        &
                                                          , MAPSF, X1 )
!
              ELSE IF ( J .EQ. 32 ) THEN
                FLONE  = .TRUE.
                FSC    = 1
                UNITS  = 'm/s'
                ENAME  = '.SDy'
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, USERO(:,3)        &
                                                          , MAPSF, X1 )
!
             ELSE IF ( J .EQ. 33 ) THEN
                FLONE  = .TRUE.
                FSC    = 1
                UNITS  = '%'
                ENAME  = '.nws'
                CALL W3S2XY ( NSEA, NSEA, NX+1, NY, USERO(:,4)        &
                                                          , MAPSF, X1 )             
!
              ELSE
                WRITE (NDSE,999)
                CALL EXTCDE ( 1 )
!
              END IF
!
! 2.b Make map
!
            DO IX=1, NX
              DO IY=1, NY
                IF ( MAPSTA(IY,IX) .EQ. 0 ) THEN
                    X1(IX,IY) = UNDEF
                    X2(IX,IY) = UNDEF
                    XX(IX,IY) = UNDEF
                    XY(IX,IY) = UNDEF
                  END IF
                IF ( X1(IX,IY) .EQ. UNDEF ) THEN
                    MAP(IX,IY) = 0
                  ELSE
                    MAP(IX,IY) = 1
                  END IF
                IF ( X2(IX,IY) .EQ. UNDEF ) THEN
                    MP2(IX,IY) = 0
                  ELSE
                    MP2(IX,IY) = 1
                  END IF
                END DO
              END DO
  

![GIR]:start  Perform output file NetCDF 
                          
 CALL WriteNc(NX,NY,X1(1:NX,1:NY),UNDEF,ENAME,FLREQ,NOGRD, TIME, & 
              X0+SX*REAL(IX1-1),X0+SX*REAL(IXN-1), Y0+SY*REAL(IY1-1),Y0+SY*REAL(IYN-1))
![GIR]:end               


!
! 2.c Perform output type 1 ( print plots )
!
            IF ( ITYPE .EQ. 1 ) THEN
!
                IF ( SCALE ) THEN
                    FSC    = 0.
                    FSCA   = 0.
                  ELSE
                    FSCA   = 1.
                  END IF
                IXB    = IX1 - IXS
!
                DO IB=1, NBLOK
                  IXA    = IXB + IXS
                  IXB    = IXA + (NXMAX-1)*IXS
                  IXB    = MIN ( IXB , IXN )
                  IF ( J .EQ. 29 ) THEN
                      CALL PRTBLK (NDSO, NX, NY, NX+1, X1, MAP, 0, FSC, &
                        IXA, IXB, IXS, IY1, IYN, IYS, IDOUT(J), UNITS)
                      CALL PRTBLK (NDSO, NX, NY, NX+1, X2, MAP, 0, FSC, &
                        IXA, IXB, IXS, IY1, IYN, IYS, IDOUT(J), UNITS)
                      CALL PRTBLK (NDSO, NX, NY, NX+1, XY, MAP, 0, FSC, &
                        IXA, IXB, IXS, IY1, IYN, IYS, IDOUT(J), UNITS)
                    ELSE IF ( FLONE ) THEN
                      CALL PRTBLK (NDSO, NX, NY, NX+1, X1, MAP, 0, FSC, &
                        IXA, IXB, IXS, IY1, IYN, IYS, IDOUT(J), UNITS)
                    ELSE IF ( VECTOR ) THEN
                      CALL PRTBLK (NDSO, NX, NY, NX+1, XX, MAP, 0, FSC, &
                        IXA, IXB, IXS, IY1, IYN, IYS, IDOUT(J), UNITS)
                      CALL PRTBLK (NDSO, NX, NY, NX+1, XY, MAP, 0, FSC, &
                        IXA, IXB, IXS, IY1, IYN, IYS, IDOUT(J), UNITS)
                    ELSE
                      CALL PRTBLK (NDSO, NX, NY, NX+1, X1, MAP, 0, FSC, &
                        IXA, IXB, IXS, IY1, IYN, IYS, IDOUT(J), UNITS)
                      CALL PRTBLK (NDSO, NX, NY, NX+1, X2, MP2, 0,FSCA, &
                        IXA, IXB, IXS, IY1, IYN, IYS, IDOUT(J), 'Deg.')
                    END IF
                  END DO
!
! 2.d Perform output type 2 ( statistics )
!
              ELSE IF ( ITYPE .EQ. 2 ) THEN
                XMIN   =  1.E20
                XMAX   = -1.E20
                XDS    =  0.D0
                XDSQ   =  0.D0
                NINGRD =  0
!
                DO IX=IX1, IXN
                  DO IY=IY1, IYN
                    IF ( MAPSTA(IY,IX) .GT. 0 .AND.                   &
                         X1(IX,IY) .NE. UNDEF ) THEN
                        NINGRD = NINGRD + 1
                        XMIN   = MIN ( XMIN , X1(IX,IY) )
                        XMAX   = MAX ( XMAX , X1(IX,IY) )
                        XDS    = XDS  + DBLE(X1(IX,IY))
                        XDSQ   = XDSQ + DBLE(X1(IX,IY))**2
                      END IF
                    END DO
                  END DO
!
                NDSDT  = NDSDT + 1
!
                IF ( NINGRD .EQ. 0 ) THEN
                    WRITE (NDSDT,940) TIME(1), IH, IM, IS
                  ELSE IF ( NINGRD .LE. 2 ) THEN
                    XAVG   = REAL ( XDS / DBLE(NINGRD) )
                    WRITE (NDSDT,940) TIME(1), IH, IM, IS,            &
                                   XMIN, XMAX
                  ELSE
                    XAVG   = REAL ( XDS / DBLE(NINGRD) )
                    XSTD   = REAL ( ( XDSQ - XDS**2/DBLE(NINGRD) )    &
                                         / DBLE(NINGRD-1) )
                    XSTD   = SQRT ( MAX ( XSTD , 0. ) )
                    WRITE (NDSDT,940) TIME(1), IH, IM, IS,            &
                                   XMIN, XMAX, XAVG, XSTD
                  END IF
!
! 2.e Perform output type 3 ( file )
!
              ELSE IF ( ITYPE .EQ. 3 ) THEN
!
                FNAME(13:16) = ENAME
                IF ( IDFM .EQ. 3 ) THEN
                    JJ     = LEN_TRIM(FNMPRE)
![GIR]                    OPEN (NDSDAT,FILE=FNMPRE(:JJ)//FNAME,             &
![GIR]                          FORM='UNFORMATTED',ERR=800,IOSTAT=IERR)
![GIR]                    WRITE (NDSDAT) FILEID, TIME,                      &
![GIR]                       X0+SX*REAL(IX1-1),X0+SX*REAL(IXN-1),IXN-IX1+1, &
![GIR]                       Y0+SY*REAL(IY1-1),Y0+SY*REAL(IYN-1),IYN-IY1+1, &
![GIR]                       ENAME, FSC, UNITS, IDLA, IDFM, FORMF, MFILL
                  ELSE
                    JJ     = LEN_TRIM(FNMPRE)
![GIR]                    OPEN (NDSDAT,FILE=FNMPRE(:JJ)//FNAME,ERR=800,     &
![GIR]                          IOSTAT=IERR)
![GIR]                    WRITE (NDSDAT,950) FILEID, TIME,                  &
![GIR]                       X0+SX*REAL(IX1-1),X0+SX*REAL(IXN-1),IXN-IX1+1, &
![GIR]                       Y0+SY*REAL(IY1-1),Y0+SY*REAL(IYN-1),IYN-IY1+1, &
![GIR]                       ENAME, FSC, UNITS, IDLA, IDFM, FORMF, MFILL
                  END IF
!
                IF ( J.EQ.2 .OR. J.EQ.3 .OR. J.EQ.5 .OR.              &
                     J.EQ.27 .OR. J.EQ.28 ) THEN
                    DO IX=IX1, IXN
                     DO IY=IY1, IYN
                        IF ( MAPSTA(IY,IX) .LE. 0 .OR.                &
                             XX(IX,IY) .EQ. UNDEF ) THEN
                            MXX(IX,IY) = MFILL
                            MYY(IX,IY) = MFILL
                          ELSE
                            MXX(IX,IY) = NINT(XX(IX,IY)/FSC)
                            MYY(IX,IY) = NINT(XY(IX,IY)/FSC)
                          END IF
                        END DO
                      END DO
                    IF ( IDLA .NE. 5 ) THEN
![GIR]                        CALL OUTA2I ( MXX, NX, NY, IX1, IXN, IY1, IYN, &
![GIR]                             NDSDAT, NDST, NDSE, IDFM, FORMF, IDLA, 1 )
![GIR]                        CALL OUTA2I ( MYY, NX, NY, IX1, IXN, IY1, IYN, &
![GIR]                             NDSDAT, NDST, NDSE, IDFM, FORMF, IDLA, 1 )
                      ELSE
                        DO IY=IY1,IYN
                          YGBX   = Y0 + REAL(IY-1)*SY
                          DO IX=IX1, IXN
                            XGBX   = X0 + REAL(IX-1)*SX
                            IF ( MXX(IX,IY) .NE. MFILL ) THEN
                                IF ( IDFM .EQ. 3 ) THEN
  ![GIR]                                    WRITE (NDSDAT)                    &
  ![GIR]                                      XGBX, YGBX, MXX(IX,IY), MYY(IX,IY)
                                  ELSE
  ![GIR]                                    WRITE (NDSDAT,951)                &
  ![GIR]                                      XGBX, YGBX, MXX(IX,IY), MYY(IX,IY)
                                  END IF
                              END IF
                            END DO
                          END DO
                      END IF
                  ELSE
                    DO IX=IX1, IXN
                      DO IY=IY1, IYN
                        IF ( MAPSTA(IY,IX) .LE. 0 .OR.                &
                             X1(IX,IY) .EQ. UNDEF ) THEN
                            MX1(IX,IY) = MFILL
                          ELSE
                            MX1(IX,IY) = NINT(X1(IX,IY)/FSC)
                          END IF
                        END DO
                      END DO
                    IF ( IDLA .NE. 5 ) THEN
![GIR]                        CALL OUTA2I ( MX1, NX, NY, IX1, IXN, IY1, IYN, &
![GIR]                             NDSDAT, NDST, NDSE, IDFM, FORMF, IDLA, 1 )
                      ELSE
                        DO IY=IY1,IYN
                          YGBX   = Y0 + REAL(IY-1)*SY
                          DO IX=IX1, IXN
                            XGBX   = X0 + REAL(IX-1)*SX
                            IF ( MX1(IX,IY) .NE. MFILL ) THEN
                                IF ( IDFM .EQ. 3 ) THEN
![GIR]                                    WRITE (NDSDAT)                    &
![GIR]                                      XGBX, YGBX, MX1(IX,IY)
                                  ELSE
![GIR]                                    WRITE (NDSDAT,951)                &
![GIR]                                      XGBX, YGBX, MX1(IX,IY)
                                  END IF
                              END IF
                            END DO
                          END DO
                      END IF
                  END IF
!
![GIR]                CLOSE (NDSDAT)
!
              END IF
!
! ... End of fields loop
!
          END IF
        END DO
!
      RETURN
!
! Error escape locations
!
  800 CONTINUE
      WRITE (NDSE,1000) IERR
      CALL EXTCDE (2)
!
! Formats
!
  940 FORMAT (1X,I8,3I3.2,2X,4E12.4)
!/LLG  950 FORMAT (1X,A13,I9.8,I7.6,2(2F8.2,I4),                      &
!/LLG              1X,A4,F8.4,1X,A10,2I2,1X,A11,I4)
!/XYG  950 FORMAT (1X,A13,I9.8,I7.6,2(2E11.3,I4),                     &
!/XYG              1X,A4,F8.4,1X,A10,2I2,1X,A11,I4)
!/LLG  951 FORMAT (1X,2F10.5,2I8)
!/XYG  951 FORMAT (1X,2E12.4,2I8)
!
  999 FORMAT (/' *** WAVEWATCH III ERROR IN W3EXGO :'/                &
               '     PLEASE UPDATE FIELDS !!! '/)
!
 1000 FORMAT (/' *** WAVEWATCH III ERROR IN W3EXGO : '/               &
               '     ERROR IN OPENING OUTPUT FILE'/                   &
               '     IOSTAT =',I5/)
!
!/T 9000 FORMAT (' TEST W3EXGO : FLAGS :',40L2)
!/T 9001 FORMAT (' TEST W3EXGO : ITPYE :',I4/                         &
!/T              '             IX1/N/S :',3I4/                        &
!/T              '             IY1/N/S :',3I4/                        &
!/T              '       SCALE, VECTOR :',2L2/                        &
!/T              '              NDSDAT :',I4)
!
!/T 9012 FORMAT (' TEST W3EXGO : BLOK PARS    : ',3I4)
!/T 9014 FORMAT ('           BASE NAME : ',A)
!
!/T 9020 FORMAT (' TEST W3EXGO : OUTPUT FIELD : ',A)
!/
!/ End of W3EXGO ----------------------------------------------------- /
!/
      END SUBROUTINE W3EXGO
!/
!/ End of W3OUTF ----------------------------------------------------- /
!/
      END PROGRAM W3OUTF