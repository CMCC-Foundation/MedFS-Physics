!/ ------------------------------------------------------------------- /
      MODULE W3IOGOMD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-May-2009 |
!/                  +-----------------------------------+
!/
!/    04-Jan-2001 : Origination.                        ( version 2.00 )
!/    23-Apr-2002 : Clean up.                           ( version 2.19 )
!/    29-Apr-2002 : Add output parameters 17-18.        ( version 2.20 )
!/    30-May-2002 : Switch clean up.                    ( version 2.21 )
!/    13-Nov-2002 : Add stress vector.                  ( version 3.00 )
!/    25-Oct-2004 : Multiple grid version.              ( version 3.06 )
!/    27-Jun-2005 : Adding MAPST2.                      ( version 3.07 )
!/    21-Jul-2005 : Adding output fields 19-21.         ( version 3.07 )
!/    23-Apr-2006 : Filter for directional spread.      ( version 3.09 )
!/    27-Jun-2006 : Adding file name preamble.          ( version 3.09 )
!/    05-Jul-2006 : Consolidate stress arrays.          ( version 3.09 )
!/    02-Apr-2007 : Adding partitioned output.          ( version 3.11 )
!/                  Adding user slots for outputs.
!/    08-Oct-2007 : Adding ST3 source term option.      ( version 3.13 )
!/                  ( F. Ardhuin )
!/    05-Mar-2008 : Added NEC sxf90 compiler directives
!/                  (Chris Bunney, UK Met Office)       ( version 3.13 )
!/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
!/
!/    Copyright 2009 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS. 
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     Gridded output of mean wave parameters.
!
!  2. Variables and types :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      VEROGR    C*10  Private  Gridded output file version number.
!      IDSTR     C*30  Private  Gridded output file ID string.
!     ----------------------------------------------------------------
!
!  3. Subroutines and functions :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      W3OUTG    Subr. Public   Calculate mean parameters.
!      W3IOGO    Subr. Public   IO to raw gridded fields file.
!     ----------------------------------------------------------------
!
!  4. Subroutines and functions used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3SETO    Subr. W3ODATMD Point to data structure.
!      W3SETG    Subr. W3GDATMD Point to data structure.
!      W3SETW    Subr. W3WDATMD Point to data structure.
!      W3SETA    Subr. W3ADATMD Point to data structure.
!      W3DIMW    Subr. W3WDATMD Allocate data structure.
!      W3DIMA    Subr. W3ADATMD Allocate data structure.
!      STRACE    Subr. W3SERVMD Subroutine tracing.           ( !/S )
!      EXTCDE    Subr. W3SERVMD Program abort with exit code.
!     ----------------------------------------------------------------
!
!  5. Remarks :
!
!     - The different output fields are not folded in with this module
!       due to the different requirements for a element '0' in some of
!       the fields.
!
!  6. Switches :
!
!       !/SHRD  Switch for shared / distributed memory architecture.
!       !/DIST  Id.
!
!       !/OMP1  OpenMP compiler directive for loop splitting.
!       !/C90   Cray FORTRAN 90 compiler directive.
!       !/NEC   NEC SXF90 compiler directives.
!
!       !/O8    Filter for low wave heights ( HSMIN )
!       !/O9    Negative wave height alowed, other mean parameters will
!             not be correct.
!
!       !/ST0   No source terms.
!       !/ST1   Source term set 1 (WAM equiv.)
!       !/ST2   Source term set 2 (Tolman and Chalikov)
!       !/ST3   Source term set 3 (WAM 4+)
!       !/STX   Open source term slot (implemented as ST0).
!
!       !/S     Enable subroutine tracing.
!       !/T     Test output.
!
!  7. Source code :
!
!/ ------------------------------------------------------------------- /
!/S      USE W3SERVMD, ONLY : STRACE
!/
      PUBLIC
!/
!/ Private parameter statements (ID strings)
!/
      CHARACTER(LEN=10), PARAMETER, PRIVATE :: VEROGR = 'III  2.03 '
      CHARACTER(LEN=30), PARAMETER, PRIVATE ::                        &
                            IDSTR = 'WAVEWATCH III GRID OUTPUT FILE'
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3OUTG ( A, FLPART )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         05-Mar-2008 |
!/                  +-----------------------------------+
!/
!/    10-Dec-1998 : Distributed FORTRAN 77 version.     ( version 1.18 )
!/    04-Jan-2000 : Upgrade to FORTRAN 90               ( version 2.00 )
!/                  Major changes to logistics.
!/    09-May-2002 : Switch clean up.                    ( version 2.21 )
!/    19-Oct-2004 : Multiple grid version.              ( version 3.06 )
!/    21-Jul-2005 : Adding output fields 19-21.         ( version 3.07 )
!/    23-Apr-2006 : Filter for directional spread.      ( version 3.09 )
!/    02-Apr-2007 : Adding partitioned output.          ( version 3.11 )
!/                  Adding user slots for outputs.
!/    08-Oct-2007 : Adding ST3 source term option.      ( version 3.13 )
!/                  ( F. Ardhuin )
!/    05-Mar-2008 : Added NEC sxf90 compiler directives
!/                  (Chris Bunney, UK Met Office)       ( version 3.13 )
!/
!  1. Purpose :
!
!     Fill necessary arrays with gridded data for output.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       A       R.A.   I   Input spectra. Left in par list to change
!                          shape.
!       FLPART  Log.   I   Flag for filling fields with part. data.
!     ----------------------------------------------------------------
!
!     Locally saved parameters
!     ----------------------------------------------------------------
!       HSMIN   Real  Filter level in Hs for calculation of mean
!                     wave parameters.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!     See module documentation.
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
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/SHRD  Switch for shared / distributed memory architecture.
!     !/DIST  Id.
!
!     !/C90   Cray FORTRAN 90 compiler directives.
!     !/NEC   NEC SXF90 compiler directives.
!     !/OMP1  OpenMP compiler directive for loop splitting.
!
!     !/O8    Filter for low wave heights ( HSMIN )
!     !/O9    Negative wave height alowed, other mean parameters will
!             not be correct.
!
!     !/ST0   No source terms.
!     !/ST1   Source term set 1 (WAM equiv.)
!     !/ST2   Source term set 2 (Tolman and Chalikov)
!     !/ST3   Source term set 3 (WAM 4+)
!     !/STX   Open source term slot (implemented as ST0).
!
!     !/S     Enable subroutine tracing.
!     !/T     Test output.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE CONSTANTS
      USE W3GDATMD
      USE W3WDATMD, ONLY: UST, FPIS
      USE W3ADATMD, ONLY: CG, WN, DW, HS, WLM, TMN, THM, THS, FP0,    &
                          THP0, FP1, THP1, ABA, ABD, UBA, UBD,        &
                          SXX, SYY, SXY, PHS, PTP, PLP, PTH, PSI, PWS,&
                          PWST, PNR, CDS, USERO
!! ini emanuela
      USE W3ADATMD, ONLY: USTKS, VSTKS,WNM
!! end emanuela
      USE W3ODATMD, ONLY: NDST, UNDEF, IAPROC, NAPROC, ICPRT, DTPRT,  &
                          WSCUT, NOSWLL
!/S      USE W3SERVMD, ONLY: STRACE
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL, INTENT(IN)        :: A(NTH,NK,0:NSEAL)
      LOGICAL, INTENT(IN)     :: FLPART
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IK, ITH, JSEA, ISEA, IX, IY,         &
                                 IKP0(NSEAL), IKP1(NSEAL),            &
                                 ILOW, ICEN, IHGH, I, J
!/S      INTEGER, SAVE           :: IENT = 0
      REAL                    :: FXPMC, FACTOR, EBAND, FKD,           &
                                 FP1STR, FP1TST, FPISTR, AABS, UABS,  &
                                 XL, XH, XL2, XH2, EL, EH, DENOM
      REAL                    :: ET(NSEAL), EWN(NSEAL), ETR(NSEAL),   &
                                 ETX(NSEAL), ETY(NSEAL), AB(NSEAL),   &
                                 ABX(NSEAL), ABY(NSEAL),              &
                                 EBD(NK,NSEAL), EC(NSEAL),            &
                                 ABR(NSEAL), UBR(NSEAL),              &
                                 ABXX(NSEAL), ABYY(NSEAL), ABXY(NSEAL)
!! ini emanuela 
      REAL                    :: USSCO, KD, FACTOR2
      REAL                    :: ETUSCX(NSEAL), ETUSCY(NSEAL)
!! end emanuela
      REAL, SAVE              :: HSMIN = 0.05
!/
!/ ------------------------------------------------------------------- /
!/
!/S      CALL STRACE (IENT, 'W3OUTG')
!
      FXPMC  = 0.66 * GRAV / 28.
      HSMIN  = HSMIN
!
! 1.  Initialize storage arrays -------------------------------------- *
!
      ET     = 0.
      EWN    = 0.
      ETR    = 0.
      ETX    = 0.
      ETY    = 0.
      ABR    = 0.
      ABA    = 0.
      ABD    = 0.
      UBR    = 0.
      UBA    = 0.
      UBD    = 0.
      SXX    = 0.
      SYY    = 0.
      SXY    = 0.
!
!! ini emanuela 
      USTKS   = 0.
      VSTKS   = 0.
      ETUSCX  = 0.
      ETUSCY  = 0.
      USSCO   = 0.
!! end emanuela
!
      HS     = UNDEF
      WLM    = UNDEF
      TMN    = UNDEF
      THM    = UNDEF
      THS    = UNDEF
      FP0    = UNDEF
      FP1    = UNDEF
      THP0   = UNDEF
      THP1   = UNDEF
!
! 2.  Integral over discrete part of spectrum ------------------------ *
!
      DO IK=1, NK
!
! 2.a Initialize energy in band
!
        AB     = 0.
        ABX    = 0.
        ABY    = 0.
        ABXX   = 0.
        ABYY   = 0.
        ABXY   = 0.
!
! 2.b Integrate energy in band
!
        DO ITH=1, NTH
!! 
!/OMP1/!$OMP PARALLEL DO PRIVATE(JSEA)
          DO JSEA=1, NSEAL
            AB (JSEA)  = AB (JSEA) + A(ITH,IK,JSEA)
            ABX(JSEA)  = ABX(JSEA) + A(ITH,IK,JSEA)*ECOS(ITH)
            ABY(JSEA)  = ABY(JSEA) + A(ITH,IK,JSEA)*ESIN(ITH)
!/DIST          ISEA         = IAPROC + (JSEA-1)*NAPROC
!/SHRD          ISEA         = JSEA
            FACTOR     = MAX ( 0.5 , CG(IK,ISEA)/SIG(IK)*WN(IK,ISEA) )
            ABXX(JSEA) = ABXX(JSEA) + ((1.+EC2(ITH))*FACTOR-0.5) *    &
                                     A(ITH,IK,JSEA)
            ABYY(JSEA) = ABYY(JSEA) + ((1.+ES2(ITH))*FACTOR-0.5) *    &
                                     A(ITH,IK,JSEA)
            ABXY(JSEA) = ABXY(JSEA) + ESC(ITH)*FACTOR * A(ITH,IK,JSEA)
          END DO
!
        END DO
!
! 2.c Finalize integration over band and update mean arrays
!
!/OMP1/!$OMP PARALLEL DO PRIVATE(JSEA,ISEA,FACTOR)
        DO JSEA=1, NSEAL
!/DIST          ISEA         = IAPROC + (JSEA-1)*NAPROC
!/SHRD          ISEA         = JSEA
          FACTOR       = DDEN(IK) / CG(IK,ISEA)
          EBD(IK,JSEA) = AB(JSEA) * FACTOR
          ET (JSEA)    = ET (JSEA) + EBD(IK,JSEA)
          EWN(JSEA)    = EWN(JSEA) + EBD(IK,JSEA) / WN(IK,ISEA)
          ETR(JSEA)    = ETR(JSEA) + EBD(IK,JSEA) / SIG(IK)
          ETX(JSEA)    = ETX(JSEA) + ABX(JSEA) * FACTOR
          ETY(JSEA)    = ETY(JSEA) + ABY(JSEA) * FACTOR
!! ini emanuela
! Directional moments in the last freq. band
!
          IF (IK.EQ.NK) THEN
            FACTOR2       = SIG(IK)**5/(GRAV**2)/DSII(IK)
            ETUSCX(JSEA)  = ABX(JSEA)*FACTOR*FACTOR2
            ETUSCY(JSEA)  = ABY(JSEA)*FACTOR*FACTOR2
            END IF
!
          KD    = MAX ( 0.001 , WN(IK,ISEA) * DW(ISEA) )
          IF ( KD .LT. 6. ) THEN
              FKD       = FACTOR / SINH(KD)**2 
!! end emanuela
              ABR(JSEA) = ABR(JSEA) + AB(JSEA) * FKD
              ABA(ISEA) = ABA(ISEA) + ABX(JSEA) * FKD
              ABD(ISEA) = ABD(ISEA) + ABY(JSEA) * FKD
              UBR(JSEA) = UBR(JSEA) + AB(JSEA) * SIG(IK)**2 * FKD
              UBA(ISEA) = UBA(ISEA) + ABX(JSEA) * SIG(IK)**2 * FKD
              UBD(ISEA) = UBD(ISEA) + ABY(JSEA) * SIG(IK)**2 * FKD
!! ini emanuela
              USSCO=FKD*SIG(IK)*WN(IK,ISEA)*COSH(2.*KD )
          ELSE
            USSCO=FACTOR*SIG(IK)*2.*WN(IK,ISEA)
            END IF
          USTKS(JSEA)  = USTKS(JSEA) + ABX(JSEA)*USSCO
          VSTKS(JSEA)  = VSTKS(JSEA) + ABY(JSEA)*USSCO
!! end emanuela
          ABXX(JSEA)   = MAX ( 0. , ABXX(JSEA) ) * FACTOR
          ABYY(JSEA)   = MAX ( 0. , ABYY(JSEA) ) * FACTOR
          ABXY(JSEA)   = ABXY(JSEA) * FACTOR
          SXX(ISEA)    = SXX(ISEA)  + ABXX(JSEA)
          SYY(ISEA)    = SYY(ISEA)  + ABYY(JSEA)
          SXY(ISEA)    = SXY(ISEA)  + ABXY(JSEA)
          EBD(IK,JSEA) = EBD(IK,JSEA) / DSII(IK)
        END DO
!
      END DO
! write(669,*) AB
! close (669)
! stop
!
! 3.  Finalize computation of mean parameters ------------------------ *
! 3.a Add tail
!     ( DTH * SIG absorbed in FTxx )
!
!/OMP1/!$OMP PARALLEL DO PRIVATE(JSEA,ISEA,EBAND)
      DO JSEA=1, NSEAL
!/DIST        ISEA      = IAPROC + (JSEA-1)*NAPROC
!/SHRD        ISEA      = JSEA
        EBAND     = AB(JSEA) / CG(NK,ISEA)
        ET (JSEA) = ET (JSEA) + FTE  * EBAND
        EWN(JSEA) = EWN(JSEA) + FTWL * EBAND
        ETR(JSEA) = ETR(JSEA) + FTTR * EBAND
        ETX(JSEA) = ETX(JSEA) + FTE * ABX(JSEA) / CG(NK,ISEA)
        ETY(JSEA) = ETY(JSEA) + FTE * ABY(JSEA) / CG(NK,ISEA)
        SXX(ISEA) = SXX(ISEA) + FTE * ABXX(JSEA) / CG(NK,ISEA)
        SYY(ISEA) = SYY(ISEA) + FTE * ABYY(JSEA) / CG(NK,ISEA)
        SXY(ISEA) = SXY(ISEA) + FTE * ABXY(JSEA) / CG(NK,ISEA)
!! ini emanuela 
! Tail for surface stokes drift is commented out: very sensitive to tail power
!	
        USTKS(JSEA)  = USTKS(JSEA) + 2*GRAV*ETUSCX(JSEA)/SIG(NK)
        VSTKS(JSEA)  = VSTKS(JSEA) + 2*GRAV*ETUSCY(JSEA)/SIG(NK)
!! end emanuela
        END DO
!
      SXX    = SXX * DWAT * GRAV
      SYY    = SYY * DWAT * GRAV
      SXY    = SXY * DWAT * GRAV
!
!/OMP1/!$OMP PARALLEL DO PRIVATE(JSEA,ISEA,IX,IY)
      DO JSEA=1, NSEAL
!/DIST        ISEA   = IAPROC + (JSEA-1)*NAPROC
!/SHRD        ISEA   = JSEA
        IX     = MAPSF(ISEA,1)
        IY     = MAPSF(ISEA,2)
        IF ( MAPSTA(IY,IX) .GT. 0 ) THEN
!/O9            IF ( ET(JSEA) .GE. 0. ) THEN
            HS (ISEA) = 4. * SQRT ( ET(JSEA) )
!/O9              ELSE
!/O9                HS (ISEA) = - 4. * SQRT ( -ET(JSEA) )
!/O9              END IF
            IF ( ET(JSEA) .GT. 1.E-7 ) THEN
                WLM(ISEA) = EWN(JSEA) / ET(JSEA) * TPI
                TMN(ISEA) = ETR(JSEA) / ET(JSEA) * TPI
                THS(ISEA) = RADE * SQRT ( MAX ( 0. , 2. * ( 1. - SQRT ( &
                MAX(0.,(ETX(JSEA)**2+ETY(JSEA)**2)/ET(JSEA)**2) ) ) ) )
                IF ( THS(ISEA) .LT. 0.01*RADE*DTH ) THS(ISEA) = 0.
!!!!!! ini emanuela
              ELSE 
                TMN(ISEA) = TPI / SIG(NK) 
                WLM(ISEA) = (SIG(NK)**2) / GRAV
!                TMN(ISEA) = 0.
!                WLM(ISEA) = 0.
                THS(ISEA) = 0.
!!!!!! end emanuela
              END IF
            IF ( ABS(ETX(JSEA))+ABS(ETY(JSEA)) .GT. 1.E-7 ) THEN
                THM(ISEA) = ATAN2(ETY(JSEA),ETX(JSEA))
              END IF
            ABR(JSEA) = SQRT ( 2. * MAX ( 0. , ABR(JSEA) ) )
            IF ( ABR(JSEA) .GE. 1.E-7 ) THEN
                ABD(ISEA) = ATAN2(ABD(ISEA),ABA(ISEA))
              ELSE
                ABD(ISEA) = 0.
              ENDIF
            ABA(ISEA) = ABR(JSEA)
            UBR(JSEA) = SQRT ( 2. * MAX ( 0. , UBR(JSEA) ) )
            IF ( UBR(JSEA) .GE. 1.E-7 ) THEN
                UBD(ISEA) = ATAN2(UBD(ISEA),UBA(ISEA))
              ELSE
                UBD(ISEA) = 0.
              ENDIF
            UBA(ISEA) = UBR(JSEA)
!            USERO(ISEA,1) = HS(ISEA) / MAX ( 0.001 , DW(ISEA) )
             USERO(ISEA,1) = CDS(ISEA)
!! ini emanuela
             USERO(ISEA,2) = USTKS(JSEA)
             USERO(ISEA,3) = VSTKS(JSEA)
             USERO(ISEA,4) = WNM(ISEA)
!! end emanuela
          END IF
        END DO
!
! 3.b Clean-up small values if !/O8 switch selected
!
!/O8      DO ISEA=IAPROC, NSEA, NAPROC
!/O8        IF ( HS(ISEA).LE.HSMIN .AND. HS(ISEA).NE.UNDEF) THEN
!/O8            WLM(ISEA) = UNDEF
!/O8            TMN(ISEA) = UNDEF
!/O8            THM(ISEA) = UNDEF
!/O8            THS(ISEA) = UNDEF
!/O8          END IF
!/O8        END DO
!
! 4.  Peak frequencies and directions -------------------------------- *
! 4.a Initialize
!
!/OMP1/!$OMP PARALLEL DO PRIVATE(JSEA,ISEA,FPISTR,FP1STR,FP1TST)
      DO JSEA=1, NSEAL
!/DIST        ISEA       = IAPROC + (JSEA-1)*NAPROC
!/SHRD        ISEA       = JSEA
        EC  (JSEA) = EBD(NK,JSEA)
        FP0 (ISEA) = UNDEF
        IKP0(JSEA) = 0
        THP0(ISEA) = UNDEF
!/ST0        FP1 (ISEA) = UNDEF
!/ST0        IKP1(JSEA) = NK
!/ST1        FP1 (ISEA) = UNDEF
!/ST1        IKP1(JSEA) = 0
!/ST2        FP1 (ISEA) = UNDEF
!/ST2        IKP1(JSEA) = NK
!/ST2        FPISTR     = MAX ( 0.003 , FPIS(ISEA) * UST(ISEA) / GRAV )
!/ST2        FP1STR     = 3.6E-4 + 0.92*FPISTR - 6.3E-10/FPISTR**3
!/ST2        FP1TST     = FP1STR / UST(ISEA) * GRAV
!/ST2        IF ( FP1TST.LE.SIG(NK) .AND. FP1TST.GT.SIG(1) ) THEN
!/ST2            FP1 (ISEA) = TPIINV * FP1TST
!/ST2            IKP1(JSEA) = MAX ( 1 , NINT(FACTI2+FACTI1*LOG(FP1TST)) )
!/ST2          END IF
!/ST3        FP1 (ISEA) = UNDEF
!/ST3        IKP1(JSEA) = 0
!/STX        FP1 (ISEA) = UNDEF
!/STX        IKP1(JSEA) = NK
!/STX        FPISTR     = MAX ( 0.003 , FPIS(ISEA) * UST(ISEA) / GRAV )
!/STX        FP1STR     = 3.6E-4 + 0.92*FPISTR - 6.3E-10/FPISTR**3
!/STX        FP1TST     = FP1STR / UST(ISEA) * GRAV
!/STX        IF ( FP1TST.LE.SIG(NK) .AND. FP1TST.GT.SIG(1) ) THEN
!/STX            FP1 (ISEA) = TPIINV * FP1TST
!/STX            IKP1(JSEA) = MAX ( 1 , NINT(FACTI2+FACTI1*LOG(FP1TST)) )
!/STX          END IF
        THP1(ISEA) = UNDEF
        END DO
!
! 4.b Discrete peak frequencies
!
      DO IK=NK-1, 2, -1
!/OMP1/!$OMP PARALLEL DO PRIVATE(JSEA,ISEA)
        DO JSEA=1, NSEAL
!/DIST          ISEA   = IAPROC + (JSEA-1)*NAPROC
!/SHRD          ISEA   = JSEA
          IF ( EC(JSEA) .LT. EBD(IK,JSEA) ) THEN
              EC  (JSEA) = EBD(IK,JSEA)
              IKP0(JSEA) = IK
            END IF
!/ST1          IF ( IKP1(JSEA).EQ.0                             &
!/ST1                 .AND. EBD(IK-1,JSEA).LT.EBD(IK,JSEA)      &
!/ST1                 .AND. EBD(IK-1,JSEA).LT.EBD(IK+1,JSEA)    &
!/ST1                 .AND. SIG(IK).GT.FXPMC/UST(ISEA)          &
!/ST1                 .AND. SIG(IK).LT.0.75*SIG(NK) )           &
!/ST1              IKP1(JSEA) = IK
!/ST3          IF ( IKP1(JSEA).EQ.0                             &
!/ST3                 .AND. EBD(IK-1,JSEA).LT.EBD(IK,JSEA)      &
!/ST3                 .AND. EBD(IK-1,JSEA).LT.EBD(IK+1,JSEA)    &
!/ST3                 .AND. SIG(IK).GT.FXPMC/UST(ISEA)          &
!/ST3                 .AND. SIG(IK).LT.0.75*SIG(NK) )           &
!/ST3              IKP1(JSEA) = IK
          END DO
        END DO
!
!/OMP1/!$OMP PARALLEL DO PRIVATE(JSEA,ISEA)
      DO JSEA=1, NSEAL
!/DIST        ISEA   = IAPROC + (JSEA-1)*NAPROC
!/SHRD        ISEA   = JSEA
        IF ( IKP0(JSEA) .NE. 0 ) FP0(ISEA) = SIG(IKP0(JSEA)) * TPIINV
!/ST1        IF ( IKP1(JSEA) .NE. 0 ) FP1(ISEA) = SIG(IKP1(JSEA)) * TPIINV
!/ST3        IF ( IKP1(JSEA) .NE. 0 ) FP1(ISEA) = SIG(IKP1(JSEA)) * TPIINV
        END DO
!
! 4.c Continuous peak frequencies
!
      XL     = 1./XFR - 1.
      XH     =  XFR - 1.
      XL2    = XL**2
      XH2    = XH**2
!
!/OMP1/!$OMP PARALLEL DO PRIVATE(JSEA,ISEA,ILOW,ICEN,IHGH,EL,EH,DENOM)
      DO JSEA=1, NSEAL
!/DIST        ISEA   = IAPROC + (JSEA-1)*NAPROC
!/SHRD        ISEA   = JSEA
        ILOW   = MAX (  1 , IKP0(JSEA)-1 )
        ICEN   = MAX (  1 , IKP0(JSEA)   )
        IHGH   = MIN ( NK , IKP0(JSEA)+1 )
        EL     = EBD(ILOW,JSEA) - EBD(ICEN,JSEA)
        EH     = EBD(IHGH,JSEA) - EBD(ICEN,JSEA)
        DENOM  = XL*EH - XH*EL
        FP0(ISEA) = FP0 (ISEA) * ( 1. + 0.5 * ( XL2*EH - XH2*EL )     &
                       / SIGN ( MAX(ABS(DENOM),1.E-15) , DENOM ) )
!/ST1        ILOW   = MAX (  1 , IKP1(JSEA)-1 )
!/ST1        ICEN   = MAX (  1 , IKP1(JSEA)   )
!/ST1        IHGH   = MIN ( NK , IKP1(JSEA)+1 )
!/ST1        EL     = EBD(ILOW,JSEA) - EBD(ICEN,JSEA)
!/ST1        EH     = EBD(IHGH,JSEA) - EBD(ICEN,JSEA)
!/ST1        DENOM  = XL*EH - XH*EL
!/ST1        FP1(ISEA) = FP1(ISEA) * ( 1. + 0.5 * (XL2*EH - XH2*EL )  &
!/ST1                       / SIGN ( MAX(ABS(DENOM),1.E-15) , DENOM ) )
!/ST3        ILOW   = MAX (  1 , IKP1(JSEA)-1 )
!/ST3        ICEN   = MAX (  1 , IKP1(JSEA)   )
!/ST3        IHGH   = MIN ( NK , IKP1(JSEA)+1 )
!/ST3        EL     = EBD(ILOW,JSEA) - EBD(ICEN,JSEA)
!/ST3        EH     = EBD(IHGH,JSEA) - EBD(ICEN,JSEA)
!/ST3        DENOM  = XL*EH - XH*EL
!/ST3        FP1(ISEA) = FP1(ISEA) * ( 1. + 0.5 * (XL2*EH - XH2*EL )  &
!/ST3                       / SIGN ( MAX(ABS(DENOM),1.E-15) , DENOM ) )
        END DO
!
! 4.d Peak directions
!
!/OMP1/!$OMP PARALLEL DO PRIVATE(JSEA)
      DO JSEA=1, NSEAL
        ETX(JSEA) = 0.
        ETY(JSEA) = 0.
        END DO
!
      DO ITH=1, NTH
!/C90/!DIR$ IVDEP
!/NEC/!CDIR NODEP
!/OMP1/!$OMP PARALLEL DO PRIVATE(JSEA,ISEA)
        DO JSEA=1, NSEAL
!/DIST          ISEA   = IAPROC + (JSEA-1)*NAPROC
!/SHRD          ISEA   = JSEA
          IF ( FP0(ISEA).NE.UNDEF) THEN
              ETX(JSEA) = ETX(JSEA) + A(ITH,IKP0(JSEA),JSEA)*ECOS(ITH)
              ETY(JSEA) = ETY(JSEA) + A(ITH,IKP0(JSEA),JSEA)*ESIN(ITH)
            END IF
          END DO
        END DO
!
!/OMP1/!$OMP PARALLEL DO PRIVATE(JSEA,ISEA)
      DO JSEA=1, NSEAL
!/DIST        ISEA   = IAPROC + (JSEA-1)*NAPROC
!/SHRD        ISEA   = JSEA
        IF ( ABS(ETX(JSEA))+ABS(ETY(JSEA)) .GT. 1.E-7 .AND.           &
             FP0(ISEA).NE.UNDEF )                                     &
            THP0(ISEA) = ATAN2(ETY(JSEA),ETX(JSEA))
        ETX(JSEA) = 0.
        ETY(JSEA) = 0.
        IKP1(JSEA) = MAX ( 1 , IKP1(JSEA) )
        END DO
!
      DO ITH=1, NTH
!/C90/!DIR$ IVDEP
!/NEC/!CDIR NODEP
!/OMP1/!$OMP PARALLEL DO PRIVATE(JSEA,ISEA)
        DO JSEA=1, NSEAL
!/DIST          ISEA   = IAPROC + (JSEA-1)*NAPROC
!/SHRD          ISEA   = JSEA
          IF ( FP1(ISEA).NE.UNDEF) THEN
              ETX(JSEA) = ETX(JSEA) + A(ITH,IKP1(JSEA),JSEA)*ECOS(ITH)
              ETY(JSEA) = ETY(JSEA) + A(ITH,IKP1(JSEA),JSEA)*ESIN(ITH)
            END IF
          END DO
        END DO
!
!/OMP1/!$OMP PARALLEL DO PRIVATE(ISEA,IX,IY)
      DO ISEA=IAPROC, NSEA, NAPROC
        IX          = MAPSF(ISEA,1)
        IY          = MAPSF(ISEA,2)
        IF ( MAPSTA(IY,IX) .LE. 0 ) THEN
            FP0 (ISEA) = UNDEF
            THP0(ISEA) = UNDEF
            FP1 (ISEA) = UNDEF
          END IF
        END DO
!
!/OMP1/!$OMP PARALLEL DO PRIVATE(ISEA,JSEA)
      DO JSEA=1, NSEAL
!/DIST        ISEA   = IAPROC + (JSEA-1)*NAPROC
!/SHRD        ISEA   = JSEA
        IF ( ABS(ETX(JSEA))+ABS(ETY(JSEA)) .GT. 1.E-7 .AND.           &
             FP1(ISEA) .NE. UNDEF )                                   &
            THP1(ISEA) = ATAN2(ETY(JSEA),ETX(JSEA))
        END DO
!
! 5.  Test output (local to MPP only)
!
!/T      WRITE (NDST,9050)
!/T      DO ISEA=IAPROC, NSEA, NAPROC
!/T        IX     = MAPSF(ISEA,1)
!/T        IY     = MAPSF(ISEA,2)
!/T        IF ( HS(ISEA) .EQ. UNDEF ) THEN
!/T            WRITE (NDST,9051) ISEA, IX, IY
!/T          ELSE IF ( WLM(ISEA) .EQ. UNDEF ) THEN
!/T            WRITE (NDST,9051) ISEA, IX, IY, HS(ISEA)
!/T          ELSE IF ( FP0(ISEA) .EQ. UNDEF ) THEN
!/T            WRITE (NDST,9051) ISEA, IX, IY, HS(ISEA), WLM(ISEA),   &
!/T                   TMN(ISEA), RADE*THM(ISEA), THS(ISEA)
!/T          ELSE IF ( FP1(ISEA) .EQ. UNDEF ) THEN
!/T            WRITE (NDST,9051) ISEA, IX, IY, HS(ISEA), WLM(ISEA),   &
!/T                   TMN(ISEA), RADE*THM(ISEA), THS(ISEA), FP0(ISEA),&
!/T                   THP0(ISEA)
!/T          ELSE
!/T            WRITE (NDST,9051) ISEA, IX, IY, HS(ISEA), WLM(ISEA),   &
!/T                   TMN(ISEA), RADE*THM(ISEA), THS(ISEA), FP0(ISEA),&
!/T                   THP0(ISEA), FP1(ISEA), THP1(ISEA)
!/T          END IF
!/T        END DO
!
! 6.  Fill arrays wth partitioned data
!
      IF ( FLPART ) THEN
!
! 6.a Initializations
!
          PHS    = UNDEF
          PTP    = UNDEF
          PLP    = UNDEF
          PTH    = UNDEF
          PSI    = UNDEF
          PWS    = UNDEF
          PWST   = UNDEF
          PNR    = UNDEF
!
! 6.b Loop over local sea points
!
          DO JSEA=1, NSEAL
!/DIST            ISEA   = IAPROC + (JSEA-1)*NAPROC
!/SHRD            ISEA   = JSEA
            IX          = MAPSF(ISEA,1)
            IY          = MAPSF(ISEA,2)
!
            IF ( MAPSTA(IY,IX).GT.0 ) THEN
                I         = ICPRT(JSEA,2)
                PNR(ISEA) = MAX ( 0. , REAL(ICPRT(JSEA,1)-1) )
                IF ( ICPRT(JSEA,1).GE.1 ) PWST(ISEA) = DTPRT(6,I)
              END IF
!
            IF ( MAPSTA(IY,IX).GT.0 .AND. ICPRT(JSEA,1).GT.1 ) THEN
                I      = ICPRT(JSEA,2) + 1
                IF ( DTPRT(6,I) .GE. WSCUT ) THEN
                    PHS(ISEA,0) = DTPRT(1,I)
                    PTP(ISEA,0) = DTPRT(2,I)
                    PLP(ISEA,0) = DTPRT(3,I)
                    PTH(ISEA,0) = ( 270. - DTPRT(4,I) ) * DERA
                    PSI(ISEA,0) = DTPRT(5,I)
                    PWS(ISEA,0) = DTPRT(6,I)
                    I      = I + 1
                  END IF
                DO J=1, NOSWLL
                  IF ( I .GT.  ICPRT(JSEA,2)+ICPRT(JSEA,1)-1 ) EXIT
                  PHS(ISEA,J) = DTPRT(1,I)
                  PTP(ISEA,J) = DTPRT(2,I)
                  PLP(ISEA,J) = DTPRT(3,I)
                  PTH(ISEA,J) = ( 270. - DTPRT(4,I) ) * DERA
                  PSI(ISEA,J) = DTPRT(5,I)
                  PWS(ISEA,J) = DTPRT(6,I)
                  I      = I + 1
                  END DO
              END IF
!
            END DO
!
        END IF
!
      RETURN
!
! Formats
!
!/T 9050 FORMAT (' TEST W3OUTG : ISEA, IX, IY, HS, L, Tm, THm, THs',     &
!/T              ', FP0, THP0, FP1, THP1')
!/T 9051 FORMAT (2X,I8,2I4,F6.2,F7.1,F6.2,2F6.1,2(F6.3,F6.0))
!/
!/ End of W3OUTG ----------------------------------------------------- /
!/
      END SUBROUTINE W3OUTG
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3IOGO ( INXOUT, NDSOG, IOTST, IMOD )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         02-Apr-2007 |
!/                  +-----------------------------------+
!/
!/    17-Mar-1999 : Distributed FORTRAN 77 version.     ( version 1.18 )
!/    04-Jan-2000 : Upgrade to FORTRAN 90               ( version 2.00 )
!/                  Major changes to logistics.
!/    24-Jan-2001 : Flat grid version (formats only)    ( version 2.06 )
!/    23-Apr-2002 : Clean up                            ( version 2.19 )
!/    29-Apr-2002 : Add output types 17-18.             ( version 2.20 )
!/    13-Nov-2002 : Add stress vector.                  ( version 3.00 )
!/    25-Oct-2004 : Multiple grid version.              ( version 3.06 )
!/    27-Jun-2005 : Adding MAPST2.                      ( version 3.07 )
!/    21-Jul-2005 : Adding output fields 19-21.         ( version 3.07 )
!/    27-Jun-2006 : Adding file name preamble.          ( version 3.09 )
!/    05-Jul-2006 : Consolidate stress arrays.          ( version 3.09 )
!/    02-Apr-2007 : Adding partitioned output.          ( version 3.11 )
!/                  Adding user slots for outputs.
!/
!  1. Purpose :
!
!     Read/write gridded output.
!
!  2. Method :
!
!     Fields in file are determined by flags in FLOGRD in W3ODATMD.
!
!         Nr  Identifies
!       --------------------------------------------
!          1  Water depth.
!          2  Current velocity.
!          3  Wind speed.
!          4  Air-sea temperature difference.
!          5  Friction velocity.
!          6  Wave height.
!          7  Mean wave length.
!          8  Mean wave period.
!          9  Mean wave direction.
!         10  Mean dirextional spread.
!         11  Peak frequency.
!         12  Peak direction.
!         13  Lowest peak frequeny / wind sea peak frequency.
!         14  Direction for 13.
!         15  Partitioned wave heights.
!         16  Partitioned peak period.
!         17  Partitioned peak wave length.
!         18  Partitioned mean direction.
!         19  Partitioned mean directional spread.
!         20  Total wind sea fraction.
!         21  Number of partitions.
!         22  Partitioned wind sea fraction.
!         23  Average time step in integration.
!         24  Cut-off frequency.
!         25  Ice concentration.
!         26  Water levels.
!         27  Near bottom rms amplitides.
!         28  Near bottom rms velocities.
!         29  Radiation stresses.
!         30  User defined #1. (requires coding ...)
!         31  User defined #1. (requires coding ...)
!
!     15-20 consist of a set of fields, index 0 = wind sea, index
!     1:NOSWLL are first NOSWLL swell fields.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       INXOUT  C*(*)  I   Test string for read/write, valid are:
!                          'READ' and 'WRITE'.
!       NDSOG   Int.   I   File unit number.
!       IOTST   Int.   O   Test indictor for reading.
!                           0 : Fields read.
!                          -1 : Past end of file.
!       IMOD    Int.   I   Model number for W3GDAT etc.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!       See module documentation above.
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3WAVE    Subr. W3WAVEMD Actual wave model routine.
!      WW3_OUTF  Prog.   N/A    Ouput postprocessor.
!      WW3_GRIB  Prog.   N/A    Ouput postprocessor.
!      GX_OUTF   Prog.   N/A    Ouput postprocessor.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!       Tests on INXOUT, file status and on array dimensions.
!
!  7. Remarks :
!
!     - MAPSTA is dumped as it contains information on the ice edge.
!       Dynamic ice edges require MAPSTA to be dumped every time step.
!     - The output file has the pre-defined name 'out_grd.FILEXT'.
!     - The current components CX and CY are converted to the absolute
!       value and direction in this routine.
!     - All written direction are in degrees, nautical convention,
!       but in reading, all is convered back to radians and cartesian
!       conventions.
!     - Before writing, wind and current directions are converted,
!       wave directions are already in correct convention (see W3OUTG).
!     - In MPP version of model data is supposed to be gatherd at the
!       correct processor before the routine is called.
!     - In MPP version routine is called by only one process, therefore
!       no test on process for error messages is needed.
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/ST1   First source term package (WAM3).
!     !/ST2   Second source term package (TC96).
!     !/S     Enable subroutine tracing.
!     !/T     Test output.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3GDATMD
!/
      USE W3WDATMD, ONLY: W3SETW, W3DIMW
      USE W3ADATMD, ONLY: W3SETA, W3DIMA
      USE W3ODATMD, ONLY: W3SETO
!/
      USE W3WDATMD, ONLY: TIME, DINIT, WLV, ICE, UST, USTDIR, ASF
      USE W3ADATMD, ONLY: AINIT, DW, UA, UD, AS, CX, CY, HS, WLM,     &
                          TMN, THM, THS, FP0, THP0, FP1, THP1, DTDYN, &
                          FCUT, ABA, ABD, UBA, UBD, SXX, SYY, SXY,    &
                          PHS, PTP, PLP, PTH, PSI, PWS, PWST, PNR,    &
                          USERO
      USE W3ODATMD, ONLY: NOGRD, IDOUT, UNDEF, NDST, NDSE, FLOGRD,    &
                          IPASS => IPASS1, WRITE => WRITE1, FNMPRE,   &
                          NOSWLL, NOEXTR
!/
      USE W3SERVMD, ONLY: EXTCDE
!/S      USE W3SERVMD, ONLY: STRACE
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(INOUT)        :: IOTST
      INTEGER, INTENT(IN)           :: NDSOG
      INTEGER, INTENT(IN), OPTIONAL :: IMOD
      CHARACTER, INTENT(IN)         :: INXOUT*(*)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IGRD, IERR, I, J, IX, IY, IO, MOGRD, &
                                 ISEA, MOSWLL
      INTEGER, ALLOCATABLE    :: MAPTMP(:,:)
!/S      INTEGER, SAVE           :: IENT = 0
      REAL                    :: AUX1(NSEA), AUX2(NSEA)
      CHARACTER(LEN=30)       :: IDTST, TNAME
      CHARACTER(LEN=10)       :: VERTST
!/
!/ ------------------------------------------------------------------- /
!/
!/S      CALL STRACE (IENT, 'W3IOGO')
!
! test input parameters ---------------------------------------------- *
!
      IF ( PRESENT(IMOD) ) THEN
          IGRD   = IMOD
        ELSE
          IGRD   = 1
        END IF
!
      CALL W3SETO ( IGRD, NDSE, NDST )
      CALL W3SETG ( IGRD, NDSE, NDST )
      CALL W3SETA ( IGRD, NDSE, NDST )
      CALL W3SETW ( IGRD, NDSE, NDST )
!
      IPASS  = IPASS + 1
      IOTST  = 0
!
      IF (INXOUT.NE.'READ' .AND. INXOUT.NE.'WRITE' ) THEN
          WRITE (NDSE,900) INXOUT
          CALL EXTCDE ( 1 )
        END IF
!
      IF ( IPASS.EQ.1 ) THEN
          WRITE  = INXOUT.EQ.'WRITE'
        ELSE
          IF ( WRITE .AND. INXOUT.EQ.'READ' ) THEN
              WRITE (NDSE,901) INXOUT
              CALL EXTCDE ( 2 )
            END IF
        END IF
!
!/T      WRITE (NDST,9000) IPASS, INXOUT, WRITE, NDSOG, IGRD, FILEXT
!
! open file ---------------------------------------------------------- *
! ( IPASS = 1 )
!
      IF ( IPASS.EQ.1 ) THEN
!
          I      = LEN_TRIM(FILEXT)
          J      = LEN_TRIM(FNMPRE)
!
!/T          WRITE (NDST,9001) FNMPRE(:J)//'out_grd.'//FILEXT(:I)
          IF ( WRITE ) THEN
              OPEN (NDSOG,FILE=FNMPRE(:J)//'out_grd.'//FILEXT(:I),    &
                    FORM='UNFORMATTED',ERR=800,IOSTAT=IERR)
            ELSE
              OPEN (NDSOG,FILE=FNMPRE(:J)//'out_grd.'//FILEXT(:I),    &
                    FORM='UNFORMATTED',ERR=800,IOSTAT=IERR,STATUS='OLD')
            END IF
!
          REWIND ( NDSOG )
!
! test info --------------------------------------------------------- *
! ( IPASS = 1 )
!
          IF ( WRITE ) THEN
              WRITE (NDSOG)                                           &
                IDSTR, VEROGR, GNAME, NOGRD, NSEA, NX, NY,            &
                X0, Y0, SX, SY, UNDEF, NOSWLL
            ELSE
              READ (NDSOG,END=801,ERR=802,IOSTAT=IERR)                &
                IDTST, VERTST, TNAME, MOGRD, NSEA, NX, NY,            &
                X0, Y0, SX, SY, UNDEF, MOSWLL
!
              IF ( IDTST .NE. IDSTR ) THEN
                  WRITE (NDSE,902) IDTST, IDSTR
                  CALL EXTCDE ( 20 )
                END IF
              IF ( VERTST .NE. VEROGR ) THEN
                  WRITE (NDSE,903) VERTST, VEROGR
                  CALL EXTCDE ( 21 )
                END IF
              IF ( NOGRD .NE. MOGRD ) THEN
                  WRITE (NDSE,904) MOGRD, NOGRD
                  CALL EXTCDE ( 22 )
                END IF
              IF ( TNAME .NE. GNAME ) THEN
                  WRITE (NDSE,905) TNAME, GNAME
                END IF
              IF ( NOSWLL .NE. MOSWLL ) THEN
                  WRITE (NDSE,906) MOSWLL, NOSWLL
                  CALL EXTCDE ( 24 )
                END IF
!
            END IF
!
!/T          WRITE (NDST,9002) IDSTR, VEROGR, GNAME, NSEA, NX, NY,    &
!/T                            X0, Y0, SX, SY, UNDEF
!
        END IF
!
! TIME and flags ----------------------------------------------------- *
!
      IF ( WRITE ) THEN
          WRITE (NDSOG)                            TIME, FLOGRD
        ELSE
          READ (NDSOG,END=803,ERR=802,IOSTAT=IERR) TIME, FLOGRD
        END IF
!
!/T      WRITE (NDST,9003) TIME, FLOGRD
!
! MAPSTA ------------------------------------------------------------- *
!
      ALLOCATE ( MAPTMP(NY,NX) )
      IF ( WRITE ) THEN
          MAPTMP = MAPSTA + 8*MAPST2
          WRITE (NDSOG)                                               &
               ((MAPTMP(IY,IX),IX=1,NX),IY=1,NY)
        ELSE
          READ (NDSOG,END=801,ERR=802,IOSTAT=IERR)                    &
               ((MAPTMP(IY,IX),IX=1,NX),IY=1,NY)
          MAPSTA = MOD(MAPTMP+2,8) - 2
          MAPST2 = (MAPTMP-MAPSTA) / 8
        END IF
      DEALLOCATE ( MAPTMP )
!
! Fields ------------------------------------------------------------- *
!
      IF ( WRITE ) THEN
          DO ISEA=1, NSEA
            IF ( MAPSTA(MAPSF(ISEA,2),MAPSF(ISEA,1)) .LT. 0 ) THEN
                UST   (ISEA) = UNDEF
                USTDIR(ISEA) = UNDEF
                DTDYN (ISEA) = UNDEF
                FCUT  (ISEA) = UNDEF
                WLM   (ISEA) = UNDEF
                TMN   (ISEA) = UNDEF
                THM   (ISEA) = UNDEF
                THS   (ISEA) = UNDEF
                ABA   (ISEA) = UNDEF
                ABD   (ISEA) = UNDEF
                UBA   (ISEA) = UNDEF
                UBD   (ISEA) = UNDEF
                SXX   (ISEA) = UNDEF
                SYY   (ISEA) = UNDEF
                SXY   (ISEA) = UNDEF
              END IF
            END DO
        ELSE
          IF (.NOT.DINIT) CALL W3DIMW ( IGRD, NDSE, NDST, .TRUE. )
          IF (.NOT.AINIT) CALL W3DIMA ( IGRD, NDSE, NDST, .TRUE. )
        END IF
!
      DO IO=1, NOGRD
        IF ( FLOGRD(IO) ) THEN
!
!/T            WRITE (NDST,9010) FLOGRD(IO), IDOUT(IO)
!
            IF ( WRITE ) THEN
!
                IF ( IO .EQ.  1 ) THEN
                    WRITE ( NDSOG ) DW(1:NSEA)
                  ELSE IF ( IO .EQ.  2 ) THEN
                    WRITE ( NDSOG ) CX(1:NSEA)
                    WRITE ( NDSOG ) CY(1:NSEA)
                  ELSE IF ( IO .EQ.  3 ) THEN
                    DO ISEA=1, NSEA
                      AUX1(ISEA) = UA(ISEA)*COS(UD(ISEA))
                      AUX2(ISEA) = UA(ISEA)*SIN(UD(ISEA))
                      END DO
                    WRITE ( NDSOG ) AUX1
                    WRITE ( NDSOG ) AUX2
                  ELSE IF ( IO .EQ.  4 ) THEN
                    WRITE ( NDSOG ) AS(1:NSEA)
                  ELSE IF ( IO .EQ.  5 ) THEN
                    DO ISEA=1, NSEA
                      IX     = MAPSF(ISEA,1)
                      IY     = MAPSF(ISEA,2)
                      IF ( MAPSTA(IY,IX) .EQ. 1 ) THEN
                          AUX1(ISEA) = UST(ISEA) * ASF(ISEA) *        &
                                                      COS(USTDIR(ISEA))
                          AUX2(ISEA) = UST(ISEA) * ASF(ISEA) *        &
                                                      SIN(USTDIR(ISEA))
                        ELSE
                          AUX1(ISEA) = UNDEF
                          AUX2(ISEA) = UNDEF
                        END IF
                      END DO
                    WRITE ( NDSOG ) AUX1
                    WRITE ( NDSOG ) AUX2
                  ELSE IF ( IO .EQ.  6 ) THEN
                    WRITE ( NDSOG ) HS(1:NSEA)
                  ELSE IF ( IO .EQ.  7 ) THEN
                    WRITE ( NDSOG ) WLM(1:NSEA)
                  ELSE IF ( IO .EQ.  8 ) THEN
                    WRITE ( NDSOG ) TMN(1:NSEA)
                  ELSE IF ( IO .EQ.  9 ) THEN
                    WRITE ( NDSOG ) THM(1:NSEA)
                  ELSE IF ( IO .EQ. 10 ) THEN
                    WRITE ( NDSOG ) THS(1:NSEA)
                  ELSE IF ( IO .EQ. 11 ) THEN
                    WRITE ( NDSOG ) FP0(1:NSEA)
                  ELSE IF ( IO .EQ. 12 ) THEN
                    WRITE ( NDSOG ) THP0(1:NSEA)
                  ELSE IF ( IO .EQ. 13 ) THEN
                    WRITE ( NDSOG ) FP1(1:NSEA)
                  ELSE IF ( IO .EQ. 14 ) THEN
                    WRITE ( NDSOG ) THP1(1:NSEA)
                  ELSE IF ( IO .EQ. 15 ) THEN
                    WRITE ( NDSOG ) PHS(1:NSEA,0:NOSWLL)
                  ELSE IF ( IO .EQ. 16 ) THEN
                    WRITE ( NDSOG ) PTP(1:NSEA,0:NOSWLL)
                  ELSE IF ( IO .EQ. 17 ) THEN
                    WRITE ( NDSOG ) PLP(1:NSEA,0:NOSWLL)
                  ELSE IF ( IO .EQ. 18 ) THEN
                    WRITE ( NDSOG ) PTH(1:NSEA,0:NOSWLL)
                  ELSE IF ( IO .EQ. 19 ) THEN
                    WRITE ( NDSOG ) PSI(1:NSEA,0:NOSWLL)
                  ELSE IF ( IO .EQ. 20 ) THEN
                    WRITE ( NDSOG ) PWS(1:NSEA,0:NOSWLL)
                  ELSE IF ( IO .EQ. 21 ) THEN
                    WRITE ( NDSOG ) PWST(1:NSEA)
                  ELSE IF ( IO .EQ. 22 ) THEN
                    WRITE ( NDSOG ) PNR(1:NSEA)
                  ELSE IF ( IO .EQ. 23 ) THEN
                    WRITE ( NDSOG ) DTDYN(1:NSEA)
                  ELSE IF ( IO .EQ. 24 ) THEN
                    WRITE ( NDSOG ) FCUT(1:NSEA)
                  ELSE IF ( IO .EQ. 25 ) THEN
                    WRITE ( NDSOG ) ICE(1:NSEA)
                  ELSE IF ( IO .EQ. 26 ) THEN
                    WRITE ( NDSOG ) WLV(1:NSEA)
                  ELSE IF ( IO .EQ. 27 ) THEN
                    DO ISEA=1, NSEA
                      IF ( ABA(ISEA) .NE. UNDEF ) THEN
                          AUX1(ISEA) = ABA(ISEA)*COS(ABD(ISEA))
                          AUX2(ISEA) = ABA(ISEA)*SIN(ABD(ISEA))
                        ELSE
                          AUX1(ISEA) = UNDEF
                          AUX2(ISEA) = UNDEF
                        END IF
                      END DO
                    WRITE ( NDSOG ) AUX1
                    WRITE ( NDSOG ) AUX2
                  ELSE IF ( IO .EQ. 28 ) THEN
                    DO ISEA=1, NSEA
                      IF ( ABA(ISEA) .NE. UNDEF ) THEN
                          AUX1(ISEA) = UBA(ISEA)*COS(UBD(ISEA))
                          AUX2(ISEA) = UBA(ISEA)*SIN(UBD(ISEA))
                        ELSE
                          AUX1(ISEA) = UNDEF
                          AUX2(ISEA) = UNDEF
                        END IF
                      END DO
                    WRITE ( NDSOG ) AUX1
                    WRITE ( NDSOG ) AUX2
                  ELSE IF ( IO .EQ. 29 ) THEN
                    WRITE ( NDSOG ) SXX(1:NSEA)
                    WRITE ( NDSOG ) SYY(1:NSEA)
                    WRITE ( NDSOG ) SXY(1:NSEA)
                  ELSE IF ( IO .EQ. 30 ) THEN
                    WRITE ( NDSOG ) USERO(1:NSEA,1)
                  ELSE IF ( IO .EQ. 31 ) THEN
                    WRITE ( NDSOG ) USERO(1:NSEA,2)
!! ini emanuela
                  ELSE IF ( IO .EQ. 32 ) THEN
                    WRITE ( NDSOG ) USERO(1:NSEA,3)
                  ELSE IF ( IO .EQ. 33 ) THEN
                    WRITE ( NDSOG ) USERO(1:NSEA,4)
!! end emanuela
                  ELSE
                    WRITE (NDSE,999)
                    CALL EXTCDE ( 30 )
                  END IF
!
              ELSE
!
                IF ( IO .EQ.  1 ) THEN
                    READ (NDSOG,END=801,ERR=802,IOSTAT=IERR) DW(1:NSEA)
                  ELSE IF ( IO .EQ.  2 ) THEN
                    READ (NDSOG,END=801,ERR=802,IOSTAT=IERR) CX(1:NSEA)
                    READ (NDSOG,END=801,ERR=802,IOSTAT=IERR) CY(1:NSEA)
                  ELSE IF ( IO .EQ.  3 ) THEN
                    READ (NDSOG,END=801,ERR=802,IOSTAT=IERR) UA(1:NSEA)
                    READ (NDSOG,END=801,ERR=802,IOSTAT=IERR) UD(1:NSEA)
                  ELSE IF ( IO .EQ.  4 ) THEN
                    READ (NDSOG,END=801,ERR=802,IOSTAT=IERR) AS(1:NSEA)
                  ELSE IF ( IO .EQ.  5 ) THEN
                    READ (NDSOG,END=801,ERR=802,IOSTAT=IERR)          &
                                                         UST(1:NSEA)
                    READ (NDSOG,END=801,ERR=802,IOSTAT=IERR)          &
                                                         USTDIR(1:NSEA)
                  ELSE IF ( IO .EQ.  6 ) THEN
                    READ (NDSOG,END=801,ERR=802,IOSTAT=IERR) HS(1:NSEA)
                  ELSE IF ( IO .EQ.  7 ) THEN
                    READ (NDSOG,END=801,ERR=802,IOSTAT=IERR) WLM(1:NSEA)
                  ELSE IF ( IO .EQ.  8 ) THEN
                    READ (NDSOG,END=801,ERR=802,IOSTAT=IERR) TMN(1:NSEA)
                  ELSE IF ( IO .EQ.  9 ) THEN
                    READ (NDSOG,END=801,ERR=802,IOSTAT=IERR) THM(1:NSEA)
                  ELSE IF ( IO .EQ. 10 ) THEN
                    READ (NDSOG,END=801,ERR=802,IOSTAT=IERR) THS(1:NSEA)
                  ELSE IF ( IO .EQ. 11 ) THEN
                    READ (NDSOG,END=801,ERR=802,IOSTAT=IERR) FP0(1:NSEA)
                  ELSE IF ( IO .EQ. 12 ) THEN
                    READ (NDSOG,END=801,ERR=802,IOSTAT=IERR)         &
                                                         THP0(1:NSEA)
                  ELSE IF ( IO .EQ. 13 ) THEN
                    READ (NDSOG,END=801,ERR=802,IOSTAT=IERR) FP1(1:NSEA)
                  ELSE IF ( IO .EQ. 14 ) THEN
                    READ (NDSOG,END=801,ERR=802,IOSTAT=IERR)         &
                                                         THP1(1:NSEA)
                  ELSE IF ( IO .EQ. 15 ) THEN
                    READ (NDSOG,END=801,ERR=802,IOSTAT=IERR)         &
                                                  PHS(1:NSEA,0:NOSWLL)
                  ELSE IF ( IO .EQ. 16 ) THEN
                    READ (NDSOG,END=801,ERR=802,IOSTAT=IERR)         &
                                                  PTP(1:NSEA,0:NOSWLL)
                  ELSE IF ( IO .EQ. 17 ) THEN
                    READ (NDSOG,END=801,ERR=802,IOSTAT=IERR)         &
                                                  PLP(1:NSEA,0:NOSWLL)
                  ELSE IF ( IO .EQ. 18 ) THEN
                    READ (NDSOG,END=801,ERR=802,IOSTAT=IERR)         &
                                                  PTH(1:NSEA,0:NOSWLL)
                  ELSE IF ( IO .EQ. 19 ) THEN
                    READ (NDSOG,END=801,ERR=802,IOSTAT=IERR)         &
                                                  PSI(1:NSEA,0:NOSWLL)
                  ELSE IF ( IO .EQ. 20 ) THEN
                    READ (NDSOG,END=801,ERR=802,IOSTAT=IERR)         &
                                                  PWS(1:NSEA,0:NOSWLL)
                  ELSE IF ( IO .EQ. 21 ) THEN
                    READ (NDSOG,END=801,ERR=802,IOSTAT=IERR)         &
                                                          PWST(1:NSEA)
                  ELSE IF ( IO .EQ. 22 ) THEN
                    READ (NDSOG,END=801,ERR=802,IOSTAT=IERR) PNR(1:NSEA)
                  ELSE IF ( IO .EQ. 23 ) THEN
                    READ (NDSOG,END=801,ERR=802,IOSTAT=IERR)         &
                                                         DTDYN(1:NSEA)
                  ELSE IF ( IO .EQ. 24 ) THEN
                    READ (NDSOG,END=801,ERR=802,IOSTAT=IERR)         &
                                                         FCUT(1:NSEA)
                  ELSE IF ( IO .EQ. 25 ) THEN
                    READ (NDSOG,END=801,ERR=802,IOSTAT=IERR) ICE(1:NSEA)
                  ELSE IF ( IO .EQ. 26 ) THEN
                    READ (NDSOG,END=801,ERR=802,IOSTAT=IERR) WLV(1:NSEA)
                  ELSE IF ( IO .EQ. 27 ) THEN
                    READ (NDSOG,END=801,ERR=802,IOSTAT=IERR) ABA(1:NSEA)
                    READ (NDSOG,END=801,ERR=802,IOSTAT=IERR) ABD(1:NSEA)
                  ELSE IF ( IO .EQ. 28 ) THEN
                    READ (NDSOG,END=801,ERR=802,IOSTAT=IERR) UBA(1:NSEA)
                    READ (NDSOG,END=801,ERR=802,IOSTAT=IERR) UBD(1:NSEA)
                  ELSE IF ( IO .EQ. 29 ) THEN
                    READ (NDSOG,END=801,ERR=802,IOSTAT=IERR) SXX(1:NSEA)
                    READ (NDSOG,END=801,ERR=802,IOSTAT=IERR) SYY(1:NSEA)
                    READ (NDSOG,END=801,ERR=802,IOSTAT=IERR) SXY(1:NSEA)
                  ELSE IF ( IO .EQ. 30 ) THEN
                    READ (NDSOG,END=801,ERR=802,IOSTAT=IERR)         &
                                                       USERO(1:NSEA,1)
                  ELSE IF ( IO .EQ. 31 ) THEN
                    READ (NDSOG,END=801,ERR=802,IOSTAT=IERR)         &
                                                       USERO(1:NSEA,2)
!! ini emanuela
                  ELSE IF ( IO .EQ. 32 ) THEN
                    READ (NDSOG,END=801,ERR=802,IOSTAT=IERR)         &
                                                       USERO(1:NSEA,3)
                  ELSE IF ( IO .EQ. 33 ) THEN
                    READ (NDSOG,END=801,ERR=802,IOSTAT=IERR)         &
                                                       USERO(1:NSEA,4)
!! end emanuela
                  ELSE
                    WRITE (NDSE,999)
                    CALL EXTCDE ( 30 )
                  END IF
!
              END IF
!
          END IF
        END DO
!
      RETURN
!
! Escape locations read errors
!
  800 CONTINUE
      WRITE (NDSE,1000) IERR
      CALL EXTCDE ( 41 )
!
  801 CONTINUE
      WRITE (NDSE,1001)
      CALL EXTCDE ( 42 )
!
  802 CONTINUE
      WRITE (NDSE,1002) IERR
      CALL EXTCDE ( 43 )
!
  803 CONTINUE
      IOTST  = -1
!/T      WRITE (NDST,9020)
      RETURN
!
! Formats
!
  900 FORMAT (/' *** WAVEWATCH III ERROR IN W3IOGO :'/                &
               '     ILEGAL INXOUT VALUE: ',A/)
  901 FORMAT (/' *** WAVEWATCH III ERROR IN W3IOGO :'/                &
               '     MIXED READ/WRITE, LAST REQUEST: ',A/)
  902 FORMAT (/' *** WAVEWATCH III ERROR IN W3IOGO :'/                &
               '     ILEGAL IDSTR, READ : ',A/                        &
               '                  CHECK : ',A/)
  903 FORMAT (/' *** WAVEWATCH III ERROR IN W3IOGO :'/                &
               '     ILEGAL VEROGR, READ : ',A/                       &
               '                   CHECK : ',A/)
  904 FORMAT (/' *** WAVEWATCH III ERROR IN W3IOGO :'/                &
               '     DIFFERENT NUMBER OF FIELDS, FILE :',I8/          &
               '                              PROGRAM :',I8/)
  905 FORMAT (/' *** WAVEWATCH III WARNING IN W3IOGO :'/              &
               '     ILEGAL GNAME, READ : ',A/                        &
               '                  CHECK : ',A/)
  906 FORMAT (/' *** WAVEWATCH III ERROR IN W3IOGO :'/                &
               '     ILEGAL NOSWLL, READ : ',I4/                      &
               '                   CHECK : ',I4/)
!
  999 FORMAT (/' *** WAVEWATCH III ERROR IN W3IOGO :'/                &
               '     PLEASE UPDATE FIELDS !!! '/)
!
 1000 FORMAT (/' *** WAVEWATCH III ERROR IN W3IOGO : '/               &
               '     ERROR IN OPENING FILE'/                          &
               '     IOSTAT =',I5/)
 1001 FORMAT (/' *** WAVEWATCH III ERROR IN W3IOGO : '/               &
               '     PREMATURE END OF FILE'/)
 1002 FORMAT (/' *** WAVEWATCH III ERROR IN W3IOGO : '/               &
               '     ERROR IN READING FROM FILE'/                     &
               '     IOSTAT =',I5/)
!
!/T 9000 FORMAT (' TEST W3IOGO : IPASS =',I4,' INXOUT = ',A,          &
!/T              ' WRITE = ',L1,' UNIT =',I3/                         &
!/T              '               IGRD =',I3,' FEXT = ',A)
!/T 9001 FORMAT (' TEST W3IOGO : OPENING NEW FILE [',A,']')
!/T 9002 FORMAT (' TEST W3IOGO : TEST PARAMETERS:'/                   &
!/T              '       IDSTR :  ',A/                                &
!/T              '      VEROGR :  ',A/                                &
!/T              '       GNAME :  ',A/                                &
!/T              '        NSEA :',I6/                                 &
!/T              '       NX,NY : ',I9,I12/                            &
!/T              '       X0,Y0 : ',2F12.2/                            &
!/T              '       SX,SY : ',2F12.2/                            &
!/T              '       UNDEF : ',F8.2)
!/T 9003 FORMAT (' TEST W3IOGO : TIME  :',I9.8,I7.6/                  &
!/T              '               FLAGS :',31L2)
!/T 9010 FORMAT (' TEST W3IOGO : PROC = ',L1,' FOR ',A)
!/T 9020 FORMAT (' TEST W3IOGO : END OF FILE REACHED')
!/
!/ End of W3IOGO ----------------------------------------------------- /
!/
      END SUBROUTINE W3IOGO
!/
!/ End of module W3IOGOMD -------------------------------------------- /
!/
      END MODULE W3IOGOMD