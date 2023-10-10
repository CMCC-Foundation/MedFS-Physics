!/ ------------------------------------------------------------------- /
      PROGRAM GXOUTF
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-May-2009 |
!/                  +-----------------------------------+
!/
!/    30-Jun-1999 : Final FORTRAN 77                    ( version 1.18 )
!/    24-Jan-2000 : Upgrade to FORTRAN 90               ( version 2.00 )
!/    25-Jan-2001 : Cartesian grid version              ( version 2.06 )
!/    29-Jan-2001 : Add output fields 17-18             ( version 2.20 )
!/    13-Nov-2002 : Add stress vector.                  ( version 3.00 )
!/    24-Dec-2004 : Multiple grid version.              ( version 3.06 )
!/    27-Jun-2005 : Adding MAPST2.                      ( version 3.07 )
!/    21-Jul-2005 : Add output fields 19-21.            ( version 3.07 )
!/    15-Dec-2005 : Updating MAPST2 for 2-way nest.     ( version 3.08 )
!/    13-Mar-2006 : MSOUT and MBOUT added.              ( version 3.09 )
!/    29-Jun-2006 : Adding file name preamble.          ( version 3.09 )
!/    05-Jul-2006 : Consolidate stress arrays.          ( version 3.09 )
!/    18-Jan-2007 : Update MSOUT/MBOUT treatment.       ( version 3.10 )
!/    28-Mar-2007 : Adding partitioned output.          ( version 3.11 )
!/                  Adding user slots for outputs.
!/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
!/
!/    Copyright 2009 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS. 
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     Generate GrADS input files from raw WAVEWATCH data file.
!
!  2. Method :
!
!     Data is read from the grid output file out_grd.ww3 (raw data)
!     and from the file gx_outf.inp ( NDSI, output requests ).
!     Model definition and raw data files are read using WAVEWATCH III
!     subroutines.
!
!     Output files are ww3.ctl and ww3.grads. the output files
!     contains a land-sea map, followed by requested fields. See the
!     control file for the names of the fields.
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
!      W3NAUX    Subr. W3ADATMD Set number of model for aux data.
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
!     - For the Cartesian grid version the X and Y increment are
!       artificially converted to long-lat by assuming the 1 degree
!       equals 100 km.
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!        !/LLG   Spherical grid.
!        !/XYG   Cartesian grid.
!
!        !/S     Enable subroutine tracing.
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
      USE W3IOGRMD, ONLY: W3IOGR
      USE W3IOGOMD, ONLY: W3IOGO
      USE W3SERVMD, ONLY : ITRACE, NEXTLN, EXTCDE
!/S      USE W3SERVMD, ONLY : STRACE
      USE W3TIMEMD, ONLY: STME21, TICK21, DSEC21
!/
      USE W3GDATMD
      USE W3WDATMD, ONLY: TIME, WLV, ICE, UST, USTDIR
      USE W3ADATMD, ONLY: DW, UA, UD, AS, CX, CY, HS, WLM, TMN, THM,  &
                          THS, FP0, THP0, FP1, THP1, DTDYN, FCUT,     &
                          ABA, ABD, UBA, UBD, SXX, SYY, SXY, USERO,   &
                          PHS, PTP, PLP, PTH, PSI, PWS, PWST, PNR
      USE W3ODATMD, ONLY: NDSE, NDST, NDSO, NOGRD, IDOUT, UNDEF,      &
                          FLOGRD, FNMPRE, NOSWLL
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: NDSI, NDSM, NDSOG, NDSDAT, NDSCTL,   &
                                 NDSTRC, NTRACE, IERR, I, J,          &
                                 TOUT(2), NOUT, TDUM(2), NVAR, IOUT,  &
                                 IX0, IXN, IY0, IYN, TIME0(2), IH0,   &
                                 IM0, ID0, IID, IJ0, IOTEST, IINC, IU,&
                                 TIMEN(2)
!/S      INTEGER, SAVE           :: IENT = 0
      REAL                    :: DTREQ, DTEST
!/LLG      REAL                    :: FAC = 1.
!/XYG      REAL                    :: FAC, XYMAX
      CHARACTER               :: COMSTR*1, IDTIME*23, IDDDAY*11,      &
                                 CINC*2
      CHARACTER*3             :: MNTH(12)
      CHARACTER*5             :: PARID
      LOGICAL                 :: FLREQ(NOGRD), MSOUT, MBOUT
!/
!/ ------------------------------------------------------------------- /
!/
      DATA TIME0  / -1, 0 /
      DATA MNTH   / 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',         &
                    'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC' /
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
      NDSCTL = 51
!
      NDSTRC =  6
      NTRACE = 10
!
      WRITE (NDSO,900)
!
      CALL ITRACE ( NDSTRC, NTRACE )
!/S      CALL STRACE (IENT, 'GXOUTF')
!
      J      = LEN_TRIM(FNMPRE)
      OPEN (NDSI,FILE=FNMPRE(:J)//'gx_outf.inp',STATUS='OLD',         &
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
      IF ( DTREQ.EQ.0 ) NOUT = 1
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
      IF ( MOD(NINT(DTREQ),60) .NE. 0 ) GOTO 810
!
! ... Output fields
!
      CALL NEXTLN ( COMSTR , NDSI , NDSE )
      READ (NDSI,*,END=801,ERR=802) FLREQ
!
      WRITE (NDSO,945)
!
      NVAR  = 1
      DO I=1, NOGRD
        IF ( FLREQ(I) ) THEN
            IF ( .NOT. FLOGRD(I) ) THEN
                WRITE (NDSO,946) IDOUT(I), '*** DATA NOT AVAILABLE ***'
              ELSE
                WRITE (NDSO,946) IDOUT(I), ' '
              END IF
          END IF
        FLREQ(I) = FLREQ(I) .AND. FLOGRD(I)
        IF ( I.EQ.29 ) THEN
            IF ( FLREQ(I) ) NVAR = NVAR + 3
          ELSE IF ( I.GE.15 .AND. I.LE.20 ) THEN
            IF ( FLREQ(I) ) NVAR = NVAR + NOSWLL + 1
          ELSE IF ( I.EQ.2 .OR. I.EQ.3 .OR. I.EQ.5 .OR. I.EQ.27 .OR.  &
            I.EQ.28 ) THEN
            IF ( FLREQ(I) ) NVAR = NVAR + 2
          ELSE
            IF ( FLREQ(I) ) NVAR = NVAR + 1
          END IF
        END DO
!
! ... Grid range
!
      CALL NEXTLN ( COMSTR , NDSI , NDSE )
      READ (NDSI,*,END=801,ERR=802) IX0, IXN, IY0, IYN, MSOUT, MBOUT
!
      WRITE (NDSO,947)
!
      IX0    = MAX ( 1, IX0 )
      IY0    = MAX ( 1, IY0 )
      IXN    = MIN ( NX, IXN )
      IYN    = MIN ( NY, IYN )
!
      WRITE (NDSO,948) IX0, IXN, IY0, IYN
!
      IF ( MSOUT ) THEN
          WRITE (NDSO,950) 'YES/--'
        ELSE
          WRITE (NDSO,950) '---/NO'
        END IF
!
      IF ( .NOT. MSOUT ) MBOUT = .FALSE.
      IF ( MBOUT ) THEN
          WRITE (NDSO,951) 'YES/--'
        ELSE
          WRITE (NDSO,951) '---/NO'
        END IF
!
      MSOUT  = .NOT. MSOUT
      MBOUT  = .NOT. MBOUT
!
      OPEN (NDSDAT,FILE=FNMPRE(:J)//'ww3.grads',FORM='UNFORMATTED',   &
            ERR=811,IOSTAT=IERR)
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 5.  Time management.
!
      IOUT   = 0
      WRITE (NDSO,970)
!
      DO
        DTEST  = DSEC21 ( TIME , TOUT )
        IF ( DTEST .GT. 0. ) THEN
            CALL W3IOGO ( 'READ', NDSOG, IOTEST )
            IF ( IOTEST .EQ. -1 ) THEN
                WRITE (NDSO,942)
                GOTO 600
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
        WRITE (NDSO,971) IDTIME
!
        CALL GXEXGO ( NX, NY, NSEA )
        TIMEN  = TOUT
!
        IF ( TIME0(1) .EQ. -1 ) TIME0 = TIME
!
        CALL TICK21 ( TOUT , DTREQ )
        IF ( IOUT .GE. NOUT ) EXIT
        END DO
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 6.  Close data file and write control file
!
  600 CONTINUE
      WRITE (NDSO,980)
!
      WRITE (NDSO,981)
      CLOSE (NDSDAT)
      WRITE (NDSO,982)
      OPEN (NDSCTL,FILE=FNMPRE(:J)//'ww3.ctl',ERR=812,IOSTAT=IERR)
!
      IH0    = TIME0(2)/10000
      IM0    = MOD(TIME0(2)/100,100)
      ID0    = MOD(TIME0(1),100)
      IID    = MOD(TIME0(1)/100,100)
      IJ0    = TIME0(1)/10000
!
      IF ( IOUT .GT. 1 ) DTREQ  = DSEC21 ( TIME0, TIMEN ) / REAL(IOUT-1)
      IF ( IOUT .EQ. 1 ) DTREQ  = 3600.
      IF ( DTREQ .GT. 3599. ) THEN
          CINC   = 'HR'
          IINC   = NINT(DTREQ/3600.)
          IF ( MOD(NINT(DTREQ),3600) .NE. 0 ) GOTO 820
        ELSE
          CINC   = 'MN'
          IINC   = NINT(DTREQ/60.)
        END IF
!
      WRITE (NDSO,983) IOUT, IH0, IM0, ID0, MNTH(IID), IJ0, IINC, CINC
!
!/XYG      XYMAX  = MAX ( ABS(X0+REAL(IX0-1)*SX), ABS(X0+REAL(IXN-1)*SX),  &
!/XYG                     ABS(Y0+REAL(IY0-1)*SY), ABS(Y0+REAL(IYN-1)*SY) )
!/XYG      IF ( XYMAX .LT. 1.E3 ) THEN
!/XYG          FAC    = 1.E-1
!/XYG        ELSE IF ( XYMAX .LT. 1.E4 ) THEN
!/XYG          FAC    = 1.E-2
!/XYG        ELSE IF ( XYMAX .LT. 1.E5 ) THEN
!/XYG          FAC    = 1.E-3
!/XYG        ELSE IF ( XYMAX .LT. 1.E6 ) THEN
!/XYG          FAC    = 1.E-4
!/XYG        ELSE
!/XYG          FAC    = 1.E-5
!/XYG        END IF
!
      WRITE (NDSCTL,990) UNDEF,                                       &
                         (1+IXN-IX0), FAC*(X0+REAL(IX0-1)*SX), FAC*SX,&
                         (1+IYN-IY0), FAC*(Y0+REAL(IY0-1)*SY), FAC*SY,&
                         1, 1000., 1.,                                &
                         IOUT, IH0, IM0, ID0, MNTH(IID), IJ0,         &
                         IINC, CINC, NVAR
!
      IU     = 99
      WRITE (NDSCTL,991) 'MAP  ', 0, IU, 'grid use map   '
!
      IF ( FLREQ(01) )                                                &
           WRITE (NDSCTL,991) 'DEPTH', 0, IU, 'Water depth    '
      IF ( FLREQ(02) )                                                &
           WRITE (NDSCTL,991) 'CX   ', 0, IU, 'Current U (m/s)'
      IF ( FLREQ(02) )                                                &
           WRITE (NDSCTL,991) 'CY   ', 0, IU, 'Current V (m/s)'
      IF ( FLREQ(03) )                                                &
           WRITE (NDSCTL,991) 'WU   ', 0, IU, 'Wind U (m/s)   '
      IF ( FLREQ(03) )                                                &
           WRITE (NDSCTL,991) 'WV   ', 0, IU, 'Wind V (m/s)   '
      IF ( FLREQ(04) )                                                &
           WRITE (NDSCTL,991) 'DELT ', 0, IU, 'AT-SST (degr)  '
      IF ( FLREQ(05) )                                                &
           WRITE (NDSCTL,991) 'USTU ', 0, IU, 'Fr.Vel. U(m/s) '
      IF ( FLREQ(05) )                                                &
           WRITE (NDSCTL,991) 'USTV ', 0, IU, 'Fr.Vel. V(m/s) '
      IF ( FLREQ(06) )                                                &
           WRITE (NDSCTL,991) 'HS   ', 0, IU, 'Wave height (m)'
      IF ( FLREQ(07) )                                                &
           WRITE (NDSCTL,991) 'LMN  ', 0, IU, 'Mean L (m)     '
      IF ( FLREQ(08) )                                                &
           WRITE (NDSCTL,991) 'TMN  ', 0, IU, 'Mean Per. (s)  '
      IF ( FLREQ(09) )                                                &
           WRITE (NDSCTL,991) 'DIRMN', 0, IU, 'Mean Dir. (rad)'
      IF ( FLREQ(10) )                                                &
           WRITE (NDSCTL,991) 'DIRSP', 0, IU, 'Dir. spread    '
      IF ( FLREQ(11) )                                                &
           WRITE (NDSCTL,991) 'PEAKP', 0, IU, 'Peak Per. (s)  '
      IF ( FLREQ(12) )                                                &
           WRITE (NDSCTL,991) 'PEAKD', 0, IU, 'Peak Dir. (rad)'
      IF ( FLREQ(13) )                                                &
           WRITE (NDSCTL,991) 'WNDSP', 0, IU, 'Wind Sea Per.  '
      IF ( FLREQ(14) )                                                &
           WRITE (NDSCTL,991) 'WNDSD', 0, IU, 'Wind Sea Dir.  '
      IF ( FLREQ(15) ) THEN
          PARID  = 'PHS  '
          DO I=0, NOSWLL
            WRITE (PARID(4:5),'(I2.2)') I
            WRITE (NDSCTL,991) PARID , 0, IU, 'Part. Hs (m)   '
            END DO
        END IF
      IF ( FLREQ(16) ) THEN
          PARID  = 'PTP  '
          DO I=0, NOSWLL
            WRITE (PARID(4:5),'(I2.2)') I
            WRITE (NDSCTL,991) PARID , 0, IU, 'Part. Tp (s)   '
            END DO
        END IF
      IF ( FLREQ(17) ) THEN
          PARID  = 'PLP  '
          DO I=0, NOSWLL
            WRITE (PARID(4:5),'(I2.2)') I
            WRITE (NDSCTL,991) PARID , 0, IU, 'Part. L  (m)   '
            END DO
        END IF
      IF ( FLREQ(18) ) THEN
          PARID  = 'PTH  '
          DO I=0, NOSWLL
            WRITE (PARID(4:5),'(I2.2)') I
            WRITE (NDSCTL,991) PARID , 0, IU, 'Part. Th (deg.)'
            END DO
        END IF
      IF ( FLREQ(19) ) THEN
          PARID  = 'PSI  '
          DO I=0, NOSWLL
            WRITE (PARID(4:5),'(I2.2)') I
            WRITE (NDSCTL,991) PARID , 0, IU, 'Part. si (deg.)'
            END DO
        END IF
      IF ( FLREQ(20) ) THEN
          PARID  = 'PWS  '
          DO I=0, NOSWLL
            WRITE (PARID(4:5),'(I2.2)') I
            WRITE (NDSCTL,991) PARID , 0, IU, 'Part. ws frac. '
            END DO
        END IF
      IF ( FLREQ(21) )                                                &
           WRITE (NDSCTL,991) 'PWST ', 0, IU, 'Total ws frac. '
      IF ( FLREQ(22) )                                                &
           WRITE (NDSCTL,991) 'PNR  ', 0, IU, 'Number of part.'
      IF ( FLREQ(23) )                                                &
           WRITE (NDSCTL,991) 'DTDYN', 0, IU, 'DTAVG ST (min) '
      IF ( FLREQ(24) )                                                &
           WRITE (NDSCTL,991) 'FCUT ', 0, IU, 'fcut (Hz)      '
      IF ( FLREQ(25) )                                                &
           WRITE (NDSCTL,991) 'ICE  ', 0, IU, 'Ice Conc. (-)  '
      IF ( FLREQ(26) )                                                &
           WRITE (NDSCTL,991) 'WLEV ', 0, IU, 'Water Level (m)'
      IF ( FLREQ(27) )                                                &
           WRITE (NDSCTL,991) 'ABX  ', 0, IU, 'a_b,rms (m)    '
      IF ( FLREQ(27) )                                                &
           WRITE (NDSCTL,991) 'ABY  ', 0, IU, 'a_b,rms (m)    '
      IF ( FLREQ(28) )                                                &
           WRITE (NDSCTL,991) 'UBX  ', 0, IU, 'u_b,rms (m/s)  '
      IF ( FLREQ(28) )                                                &
           WRITE (NDSCTL,991) 'UBY  ', 0, IU, 'u_b,rms (m/s)  '
      IF ( FLREQ(29) )                                                &
           WRITE (NDSCTL,991) 'SXX  ', 0, IU, 'Sxx (N/m)      '
      IF ( FLREQ(29) )                                                &
           WRITE (NDSCTL,991) 'SYY  ', 0, IU, 'Syy (N/m)      '
      IF ( FLREQ(29) )                                                &
           WRITE (NDSCTL,991) 'SXY  ', 0, IU, 'Sxy (N/m)      '
      IF ( FLREQ(30) )                                                &
           WRITE (NDSCTL,991) 'USR1 ', 0, IU, 'Drag Coeff. (-) '
!! ini emanuela
      IF ( FLREQ(31) )                                                &
           WRITE (NDSCTL,991) 'USR2 ', 0, IU, 'X StokesD(m/s) '
      IF ( FLREQ(32) )                                                &
           WRITE (NDSCTL,991) 'USR3 ', 0, IU, 'Y StokesD(m/s) '
!!end emanuela
      IF ( FLREQ(33) )                                                &
           WRITE (NDSCTL,991) 'USR4 ', 0, IU, 'User defined 2 '
!
      WRITE (NDSCTL,992)
!
      GOTO 888
!
! Escape locations read errors :
!
  800 CONTINUE
      WRITE (NDSE,1000) IERR
      CALL EXTCDE ( 1 )
!
  801 CONTINUE
      WRITE (NDSE,1001)
      CALL EXTCDE ( 2 )
!
  802 CONTINUE
      WRITE (NDSE,1002) IERR
      CALL EXTCDE ( 3 )
!
  810 CONTINUE
      WRITE (NDSE,1010)
      CALL EXTCDE ( 10 )
!
  811 CONTINUE
      WRITE (NDSE,1011)
      CALL EXTCDE ( 11 )
!
  812 CONTINUE
      WRITE (NDSE,1012)
      CALL EXTCDE ( 12 )
!
  820 CONTINUE
      WRITE (NDSE,1020) DTREQ
      CALL EXTCDE ( 20 )
!
  821 CONTINUE
      WRITE (NDSE,1021)
      CALL EXTCDE ( 21 )
!
  888 CONTINUE
      WRITE (NDSO,999)
!
! Formats
!
  900 FORMAT (/12X,'   *** WAVEWATCH III GrADS field output postp. ***   '/ &
               12X,'====================================================='/)
  901 FORMAT ( '  Comment character is ''',A,''''/)
!
  920 FORMAT ( '  Grid name : ',A/)
!
  930 FORMAT ( '  Fields in file : '/                                 &
               ' --------------------------')
  931 FORMAT ( '      ',A)
!
  940 FORMAT (/'  Output time data : '/                               &
               ' -----------------------------------------------------'/ &
               '      First time         : ',A)
  941 FORMAT ( '      Interval           : ',A/                       &
               '      Number of requests : ',I4)
  942 FORMAT (/'      End of file reached '/)
!
  945 FORMAT (/'  Requested output fields : '/                        &
               ' -----------------------------------------------------')
  946 FORMAT ( '      ',A,1X,A)
!
  947 FORMAT (/'  Requested discrete grid ranges : '/                 &
               ' -----------------------------------------------------')
  948 FORMAT ( '      Longitudes         : ',2I6/                     &
               '      lattidutes         : ',2I6/                     &
               '      Opening file ww3.grads')
  949 FORMAT ( '      Alternative definition is used ')
  950 FORMAT ( '      Sea points in mask :      ',A)
  951 FORMAT ( '      Bound. pts. in mask:      ',A)
!
  970 FORMAT (//'  Generating file '/                                 &
               ' -----------------------------------------------------')
  971 FORMAT ( '      Data for ',A)
!
  980 FORMAT (//'  Final file management '/                           &
               ' -----------------------------------------------------')
  981 FORMAT ( '      Closing file ww3.grads')
  982 FORMAT ( '      Opening file ww3.ctl')
  983 FORMAT ( '         Number of times : ',I6/                      &
           '         Initial time ID : ',I2.2,':',I2.2,'Z',I2.2,A3,I4/ &
           '         Time step ID    : ',I2,A2)
!
  990 FORMAT ('DSET      ww3.grads'/                                  &
              'TITLE     WAVEWATCH III gridded data'/                 &
              'OPTIONS   sequential'/                                 &
              'UNDEF    ',F8.1/                                       &
              'XDEF     ',I4,'  LINEAR ',2F12.5/                      &
              'YDEF     ',I4,'  LINEAR ',2F12.5/                      &
              'ZDEF     ',I4,'  LINEAR ',2F12.5/                      &
              'TDEF     ',I4,'  LINEAR ',I6.2,':',I2.2,'Z',I2.2,A3,I4, &
               2x,I2,A2/                                              &
              'VARS     ',I4)
  991 FORMAT (A5,2I4,2X,A15)
  992 FORMAT ('ENDVARS')
!
  999 FORMAT (/'  End of program '/                                   &
               ' ========================================='/          &
               '         WAVEWATCH III GrADS field output '/)
!
!/T 9050 FORMAT ( ' TEST GXOUTF : KPDS : ',13I4/                      &
!/T               '                      ',12I4)
!/T 9051 FORMAT ( ' TEST GXOUTF : KGDS : ',8I6/                       &
!/T               '                      ',8I6/                       &
!/T               '                      ',6I6)
!
 1000 FORMAT (/' *** WAVEWATCH III ERROR IN GXOUTF : '/               &
               '     ERROR IN OPENING INPUT FILE'/                    &
               '     IOSTAT =',I5/)
!
 1001 FORMAT (/' *** WAVEWATCH III ERROR IN GXOUTF : '/               &
               '     PREMATURE END OF INPUT FILE'/)
!
 1002 FORMAT (/' *** WAVEWATCH III ERROR IN GXOUTF : '/               &
               '     ERROR IN READING FROM INPUT FILE'/               &
               '     IOSTAT =',I5/)
!
 1010 FORMAT (/' *** WAVEWATCH III ERROR IN GXOUTF : '/               &
               '     SMALLEST OUTPUT INCREMENT IS 60 SEC.'/)
!
 1011 FORMAT (/' *** WAVEWATCH III ERROR IN GXOUTF : '/               &
               '     ERROR IN OPENING OUTPUT FILE ww3.grads'/         &
               '     IOSTAT =',I5/)
!
 1012 FORMAT (/' *** WAVEWATCH III ERROR IN GXOUTF : '/               &
               '     ERROR IN OPENING OUTPUT FILE ww3.ctl'/           &
               '     IOSTAT =',I5/)
!
 1020 FORMAT (/' *** WAVEWATCH III ERROR IN GXOUTF : '/               &
               '     FIELD INCREMENT > 1HR BUT NOT MULTIPLE',F10.0/)
!
 1021 FORMAT (/' *** WAVEWATCH III ERROR IN GXOUTF : '/               &
               '     UPDATE PARS IN LOOP 610 !!!'/)
!/
!/ Internal subroutine GXEXGO ---------------------------------------- /
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE GXEXGO ( NX, NY, NSEA )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         28-Mar-2007 |
!/                  +-----------------------------------+
!/
!/    30-Jun-1999 : Final FORTRAN 77                    ( version 1.18 )
!/    24-Jan-2000 : Upgrade to FORTRAN 90               ( version 2.00 )
!/                  Massive changes to logistics.
!/    29-Jan-2001 : Add output fields 17-18             ( version 2.20 )
!/    24-Dec-2004 : Multiple grid version.              ( version 3.06 )
!/    27-Jun-2005 : Adding MAPST2.                      ( version 3.07 )
!/    21-Jul-2005 : Add output fields 19-21.            ( version 3.07 )
!/    05-Jul-2006 : Consolidate stress arrays.          ( version 3.09 )
!/    18-Jan-2007 : Update MSOUT/MBOUT treatment.       ( version 3.10 )
!/    28-Mar-2007 : Adding partitioned output.          ( version 3.11 )
!/                  Adding user slots for outputs.
!/
!  1. Purpose :
!
!     Perform actual output for GrADS postprocessing.
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
!       X1, X2, XX
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
!     ---------------------------------------------------------------
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
!     - VALLND added to assure that map with only land plots in 
!       GrADS.
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/S  Enable subroutine tracing.
!     !/T  Enable test output.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3SERVMD, ONLY: W3S2XY
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: NX, NY, NSEA
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IX, IY, J, ISEA, IXL, IXR
      INTEGER                 :: MAPXCL(NY,NX), MAPDRY(NY,NX),        &
                                 MAPICE(NY,NX), MAPLND(NY,NX),        &
                                 MAPMSK(NY,NX)
!/S      INTEGER, SAVE           :: IENT   = 0
      REAL                    :: X1(NX,NY), XX(NX,NY), XY(NX,NY),     &
                                 XA(NX,NY,0:NOSWLL)
      REAL                    :: VALLND = 0.001
!/
!/ ------------------------------------------------------------------- /
!/
!/S      CALL STRACE (IENT, 'GXEXGO')
!
!/T      WRITE (NDST,9000) (FLREQ(J),J=1,NOGRD)
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 1.  Preparations
! 1.a Write map to file
!
      MAPXCL = MOD(MAPST2,2)
      MAPICE = MOD(MAPST2,2)
      MAPDRY = MOD(MAPST2/2,2)
      MAPLND = MOD(MAPST2/4,2)
      MAPMSK = MOD(MAPST2/8,2)
!
      DO IY=1, NY
        DO IX=1, NX
          IF ( MAPSTA(IY,IX).EQ.0 ) THEN
              IF ( MAPXCL(IY,IX).EQ.1 ) THEN
                  X1(IX,IY) = UNDEF
                ELSE
                  X1(IX,IY) = VALLND
                END IF
            ELSE IF ( MAPSTA(IY,IX).LT.0 ) THEN
              IF ( MAPMSK(IY,IX).EQ.1 ) THEN
                  X1(IX,IY) = -4.
                ELSE IF ( MAPLND(IY,IX).EQ.1 ) THEN
                  X1(IX,IY) = VALLND
                ELSE IF ( MAPICE(IY,IX).EQ.1 .AND.                    &
                          MAPDRY(IY,IX).EQ.1 ) THEN
                  X1(IX,IY) = -3.
                ELSE IF ( MAPDRY(IY,IX).EQ.1 ) THEN
                  X1(IX,IY) = -2.
                ELSE IF ( MAPICE(IY,IX).EQ.1 ) THEN
                  X1(IX,IY) = -1.
                ELSE
                  X1(IX,IY) = -5.
                END IF
            ELSE
              X1(IX,IY) = REAL(MAPSTA(IY,IX))
              IF ( MSOUT ) THEN
                  IF ( MAPSTA(IY,IX) .GT. 0 ) X1(IX,IY) = UNDEF
                ELSE IF ( MBOUT ) THEN
                  IF ( MAPSTA(IY,IX).EQ.2  .OR.                      &
                       IY.EQ.1 .OR. IY.EQ.NY .OR.                    &
                       ( .NOT.GLOBAL .AND.                           &
                             (IX.EQ.1 .OR. IX.EQ.NX) ) ) THEN
                      X1(IX,IY) = UNDEF
                    ELSE
                      IXl    = 1 + MOD(IX+NX-2,NX)
                      IXR    = 1 + MOD(IX,NX)
                      IF ( MAPSTA(IY+1,IXL).EQ.0 .AND.               &
                           MAPXCL(IY+1,IXL).EQ.1 ) X1(IX,IY) = UNDEF
                      IF ( MAPSTA(IY+1,IX ).EQ.0 .AND.               &
                           MAPXCL(IY+1,IX ).EQ.1 ) X1(IX,IY) = UNDEF
                      IF ( MAPSTA(IY+1,IXR).EQ.0 .AND.               &
                           MAPXCL(IY+1,IXR).EQ.1 ) X1(IX,IY) = UNDEF
                      IF ( MAPSTA( IY ,IXR).EQ.0 .AND.               &
                           MAPXCL( IY ,IXR).EQ.1 ) X1(IX,IY) = UNDEF
                      IF ( MAPSTA(IY-1,IXR).EQ.0 .AND.               &
                           MAPXCL(IY-1,IXR).EQ.1 ) X1(IX,IY) = UNDEF
                      IF ( MAPSTA(IY-1,IX ).EQ.0 .AND.               &
                           MAPXCL(IY-1,IX ).EQ.1 ) X1(IX,IY) = UNDEF
                      IF ( MAPSTA(IY-1,IXL).EQ.0 .AND.               &
                           MAPXCL(IY-1,IXL).EQ.1 ) X1(IX,IY) = UNDEF
                      IF ( MAPSTA( IY ,IXL).EQ.0 .AND.               &
                           MAPXCL( IY ,IXL).EQ.1 ) X1(IX,IY) = UNDEF
                    END IF
                END IF
              IF ( MSOUT .AND. MAPSTA(IY,IX).EQ.1 ) X1(IX,IY) = UNDEF
              IF ( MBOUT .AND. MAPSTA(IY,IX).EQ.2 ) X1(IX,IY) = UNDEF
            END IF
          VALLND = - VALLND
          END DO
        END DO
!
      WRITE (NDSDAT) ((X1(IX,IY),IX=IX0,IXN),IY=IY0,IYN)
!
! 1.b Initialize arrays
!
      X1 = UNDEF
      XX = UNDEF
      XY = UNDEF
      XA = UNDEF
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 2.  Loop over output fields.
!
      DO J=1, NOGRD
        IF ( FLREQ(J) ) THEN
!
!/T            WRITE (NDST,9020) IDOUT(J)
!
! 2.a Set output arrays and parameters
!
!     Water depth
!
            IF ( J .EQ.  1 ) THEN
                CALL W3S2XY ( NSEA, NSEA, NX, NY, DW(1:NSEA)          &
                                                        , MAPSF, X1 )
!
!     Current
!
              ELSE IF ( J .EQ.  2 ) THEN
                CALL W3S2XY ( NSEA, NSEA, NX, NY, CX(1:NSEA)          &
                                                        , MAPSF, XX )
                CALL W3S2XY ( NSEA, NSEA, NX, NY, CY(1:NSEA)          &
                                                        , MAPSF, XY )
!
!     Wind speed
!
              ELSE IF ( J .EQ.  3 ) THEN
                CALL W3S2XY ( NSEA, NSEA, NX, NY, UA(1:NSEA)          &
                                                        , MAPSF, XX )
                CALL W3S2XY ( NSEA, NSEA, NX, NY, UD(1:NSEA)          &
                                                        , MAPSF, XY )
!
!     Air-sea temp. dif.
!
              ELSE IF ( J .EQ.  4 ) THEN
                CALL W3S2XY ( NSEA, NSEA, NX, NY, AS(1:NSEA)          &
                                                        , MAPSF, X1 )
!
!     Friction velocity
!
              ELSE IF ( J .EQ.  5 ) THEN
                CALL W3S2XY ( NSEA, NSEA, NX, NY, UST   (1:NSEA)      &
                                                        , MAPSF, XX )
                CALL W3S2XY ( NSEA, NSEA, NX, NY, USTDIR(1:NSEA)      &
                                                        , MAPSF, XY )
!
!     Significant wave height
!
              ELSE IF ( J .EQ.  6 ) THEN
                CALL W3S2XY ( NSEA, NSEA, NX, NY, HS    , MAPSF, X1 )
!
!     Mean wave length
!
              ELSE IF ( J .EQ.  7 ) THEN
                CALL W3S2XY ( NSEA, NSEA, NX, NY, WLM   , MAPSF, X1 )
!
!     Mean wave period
!
              ELSE IF ( J .EQ.  8 ) THEN
                CALL W3S2XY ( NSEA, NSEA, NX, NY, TMN   , MAPSF, X1 )
!
!     Mean wave direction
!
              ELSE IF ( J .EQ.  9 ) THEN
                CALL W3S2XY ( NSEA, NSEA, NX, NY, THM   , MAPSF, X1 )
!
!     Directional spread
!
              ELSE IF ( J .EQ. 10 ) THEN
                CALL W3S2XY ( NSEA, NSEA, NX, NY, THS   , MAPSF, X1 )
!
!     Peak frequency
!
              ELSE IF ( J .EQ. 11 ) THEN
                DO ISEA=1, NSEA
                  IF ( FP0(ISEA) .NE. UNDEF ) THEN
                      FP0(ISEA) = 1. / FP0(ISEA)
                    END IF
                  END DO
                CALL W3S2XY ( NSEA, NSEA, NX, NY, FP0   , MAPSF, X1 )
!
!     Peak direction
!
              ELSE IF ( J .EQ. 12 ) THEN
                CALL W3S2XY ( NSEA, NSEA, NX, NY, THP0  , MAPSF, X1 )
!
!     Wind sea peak frequency
!
              ELSE IF ( J .EQ. 13 ) THEN
                DO ISEA=1, NSEA
                  IF ( FP1(ISEA) .NE. UNDEF ) THEN
                      FP1(ISEA) = 1. / FP1(ISEA)
                    END IF
                  END DO
                CALL W3S2XY ( NSEA, NSEA, NX, NY, FP1   , MAPSF, X1 )
!
!     Wind sea peak direction
!
              ELSE IF ( J .EQ. 14 ) THEN
                CALL W3S2XY ( NSEA, NSEA, NX, NY, THP1  , MAPSF, X1 )
!
!     Partitioned wave heights
!
              ELSE IF ( J .EQ. 15 ) THEN
                DO I=0, NOSWLL
                  CALL W3S2XY ( NSEA, NSEA, NX, NY, PHS(:,I),         &
                                                   MAPSF, XA(:,:,I) )
                  END DO
!
!     Partitioned peak period
!
              ELSE IF ( J .EQ. 16 ) THEN
                DO I=0, NOSWLL
                  CALL W3S2XY ( NSEA, NSEA, NX, NY, PTP(:,I),         &
                                                   MAPSF, XA(:,:,I) )
                  END DO
!
!     Partitioned wave leangths (peak)
!
              ELSE IF ( J .EQ. 17 ) THEN
                DO I=0, NOSWLL
                  CALL W3S2XY ( NSEA, NSEA, NX, NY, PLP(:,I),         &
                                                   MAPSF, XA(:,:,I) )
                  END DO
!
!     Partitioned directions
!
              ELSE IF ( J .EQ. 18 ) THEN
                DO I=0, NOSWLL
                  CALL W3S2XY ( NSEA, NSEA, NX, NY, PTH(:,I),         &
                                                   MAPSF, XA(:,:,I) )
                  END DO
!
!     Partitioned direstional spread
!
              ELSE IF ( J .EQ. 19 ) THEN
                DO I=0, NOSWLL
                  CALL W3S2XY ( NSEA, NSEA, NX, NY, PSI(:,I),         &
                                                   MAPSF, XA(:,:,I) )
                  END DO
!
!     Partitioned wind sea fraction
!
              ELSE IF ( J .EQ. 20 ) THEN
                DO I=0, NOSWLL
                  CALL W3S2XY ( NSEA, NSEA, NX, NY, PWS(:,I),         &
                                                   MAPSF, XA(:,:,I) )
                  END DO
!
!     Total wind sea fraction
!
              ELSE IF ( J .EQ. 21 ) THEN
                DO I=0, NOSWLL
                  CALL W3S2XY ( NSEA, NSEA, NX, NY, PWST(:),MAPSF, X1 )
                  END DO
!
!     Number of artitions
!
              ELSE IF ( J .EQ. 22 ) THEN
                DO I=0, NOSWLL
                  CALL W3S2XY ( NSEA, NSEA, NX, NY, PNR(:), MAPSF, X1 )
                  END DO
!
!     Average source term time step
!
              ELSE IF ( J .EQ. 23 ) THEN
                DO ISEA=1, NSEA
                  IF ( DTDYN(ISEA) .NE. UNDEF )                       &
                       DTDYN(ISEA) = DTDYN(ISEA) / 60.
                  END DO
                CALL W3S2XY ( NSEA, NSEA, NX, NY, DTDYN , MAPSF, X1 )
!
!     Cut-off frequency
!
              ELSE IF ( J .EQ. 24 ) THEN
                CALL W3S2XY ( NSEA, NSEA, NX, NY, FCUT  , MAPSF, X1 )
!
!     Ice concentration
!
              ELSE IF ( J .EQ. 25 ) THEN
                CALL W3S2XY ( NSEA, NSEA, NX, NY, ICE   , MAPSF, X1 )
!
!     Water level
!
              ELSE IF ( J .EQ. 26 ) THEN
                CALL W3S2XY ( NSEA, NSEA, NX, NY, WLV   , MAPSF, X1 )
!
!     Near-bottom amplitude
!
              ELSE IF ( J .EQ. 27 ) THEN
                CALL W3S2XY ( NSEA, NSEA, NX, NY, ABA   , MAPSF, XX )
                CALL W3S2XY ( NSEA, NSEA, NX, NY, ABD   , MAPSF, XY )
!
!     Near-bottom velocity
!
              ELSE IF ( J .EQ. 28 ) THEN
                CALL W3S2XY ( NSEA, NSEA, NX, NY, UBA   , MAPSF, XX )
                CALL W3S2XY ( NSEA, NSEA, NX, NY, UBD   , MAPSF, XY )
!
!     Radiation stresses
!
              ELSE IF ( J .EQ. 29 ) THEN
                CALL W3S2XY ( NSEA, NSEA, NX, NY, SXX   , MAPSF, X1 )
                CALL W3S2XY ( NSEA, NSEA, NX, NY, SYY   , MAPSF, XX )
                CALL W3S2XY ( NSEA, NSEA, NX, NY, SXY   , MAPSF, XY )
!
!     User defined #1 : Drag Coefficient
!
              ELSE IF ( J .EQ. 30 ) THEN
                CALL W3S2XY ( NSEA, NSEA, NX, NY, USERO(:,1)          &
                                                        , MAPSF, X1 )
!
!     User defined #2 : Stokes Drift
!
              ELSE IF ( J .EQ. 31 ) THEN
                CALL W3S2XY ( NSEA, NSEA, NX, NY, USERO(:,2)          &
                                                        , MAPSF, X1 )
!! ini emanuela
!
!     User defined #3 : Stokes Drift 
!
              ELSE IF ( J .EQ. 32 ) THEN
                CALL W3S2XY ( NSEA, NSEA, NX, NY, USERO(:,3)          &
                                                        , MAPSF, X1 )
!
!     User defined #4
!
              ELSE IF ( J .EQ. 33 ) THEN
                CALL W3S2XY ( NSEA, NSEA, NX, NY, USERO(:,4)          &
                                                        , MAPSF, X1 )
!! end emanuela
!
              ELSE
                WRITE (NDSE,999)
                CALL EXTCDE ( 1 )
!
              END IF
!
! 3   Perform output
!
            IF ( J.EQ.29 ) THEN
                WRITE (NDSDAT)                                        &
                      ((X1(IX,IY),IX=IX0,IXN),IY=IY0,IYN)
                WRITE (NDSDAT)                                        &
                      ((XX(IX,IY),IX=IX0,IXN),IY=IY0,IYN)
                WRITE (NDSDAT)                                        &
                      ((XY(IX,IY),IX=IX0,IXN),IY=IY0,IYN)
              ELSE IF ( J.GE.15 .AND. J.LE.20 ) THEN
                DO I=0, NOSWLL
                  WRITE (NDSDAT)                                      &
                    ((XA(IX,IY,I),IX=IX0,IXN),IY=IY0,IYN)
                  END DO
              ELSE IF ( J.EQ.2 .OR. J.EQ.3 .OR. J.EQ.5.               &
                .OR. J.EQ.27 .OR. J.EQ.28 ) THEN
                WRITE (NDSDAT)                                        &
                      ((XX(IX,IY),IX=IX0,IXN),IY=IY0,IYN)
                WRITE (NDSDAT)                                        &
                      ((XY(IX,IY),IX=IX0,IXN),IY=IY0,IYN)
              ELSE
                WRITE (NDSDAT)                                        &
                      ((X1(IX,IY),IX=IX0,IXN),IY=IY0,IYN)
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
! Formats
!
  940 FORMAT (1X,I8,3I3.2,2X,4E12.4)
  950 FORMAT (1X,A13,I9.8,I7.6,2(2F8.2,I4),                           &
              1X,A4,F8.4,1X,A10,2I2,1X,A11,I4)
  951 FORMAT (1X,2F10.5,2I8)
!
  999 FORMAT (/' *** WAVEWATCH III ERROR IN GXEXGO :'/                &
               '     PLEASE UPDATE FIELDS !!! '/)
!
!/T 9000 FORMAT (' TEST GXEXGO : FLAGS  :',40L2)
!
!/T 9020 FORMAT (' TEST GXEXGO : OUTPUT FIELD : ',A)
!/
!/ End of GXEXGO ----------------------------------------------------- /
!/
      END SUBROUTINE GXEXGO
!/
!/ End of GXOUTF ----------------------------------------------------- /
!/
      END PROGRAM GXOUTF