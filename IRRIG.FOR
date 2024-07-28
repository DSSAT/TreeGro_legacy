C=======================================================================
C  IRRIG, Subroutine
C  Determines when irrigation occurs
C-----------------------------------------------------------------------
C  REVISION HISTORY
C  09/01/1988 BB  Restructured from WATBAL
C  02/07/1993 PWW Header revision and minor changes
C  02/25/1994 BDB Added Fixed amount auto irrigation
C  04/19/1994 BDB Added NAP calc. to auto irrigations
C  09/19/1995 GH  Modified to match GRO code
C  07/11/1996 GH  Separated irrigation and precipitation
C  03/23/1997 JWJ Added P and W options
C  10/17/1997 CHP Updated for modular format.
C  09/01/1999 GH  Incorporated into CROPGRO
C  10/15/2001 CHP IRRIG now called from Operations Management module
C                   (MGMTOPS). Output variables are no longer exported.
C  04/16/2002 GH  Adjusted for crop rotations
C  06/12/2002 CHP/US Added flooded irrigation option from rice model
C  08/01/2002 CHP Merged RUNINIT and SEASINIT into INIT section
C  08/19/2002 GH  Modified for Y2K
C  08/12/2003 CHP Added I/O error checking
C  08/12/2003 CHP Fixed problem with automatic irrigation
C  10/28/2004 CHP Fixed problem with multiple applications on same day.
!  06/06/2006 CHP Export TIL_IRR, the irrigation amount which affects 
!                 soil dynamics (excludes drip irrigation).
!  01/11/2007 CHP Changed GETPUT calls to GET and PUT
!  03/30/2011 CHP Add Automatic Management for IIRRI
!             Use (DYNAMIC .EQ. INIT) for both RUNINIT and SEASINIT
!  08/15/2011 error handling for very small drip duration
!             remove using flooded conditions with Century method 
!             Now PUT('MGMT', Irr info) are in daily rate which were in daily integration
C-----------------------------------------------------------------------
C  Called by: WATBAL
C  Calls  : None
C=======================================================================
      SUBROUTINE IRRIG(CONTROL, ISWITCH,
     &    RAIN, SOILPROP, SW, MDATE, YRPLT,               !Input
     &    FLOODWAT, IIRRI, IRRAMT, NAP, TIL_IRR, TOTIR)   !Output

!-----------------------------------------------------------------------
      USE ModuleDefs 
      USE ModuleData
      USE FloodModule
      IMPLICIT NONE
      SAVE
!-----------------------------------------------------------------------
      CHARACTER*6 ERRKEY
!      CHARACTER*70 IrrText
      PARAMETER (ERRKEY = 'IRRIG')

      CHARACTER*1  IIRRI, ISWWAT, PLME, RNMODE, MEHYD   !, MESOM
      CHARACTER*6  SECTION
      CHARACTER*30 FILEIO
      CHARACTER*98 CHAR
      CHARACTER*78 MSG(10)

      INTEGER AIRRCOD, DAP, DYNAMIC, ERRNUM, FOUND, I, IDATE
      INTEGER LUNIO, LINC, LNUM
      INTEGER MULTI, NAP, NDRIP, NAPW, JIRR, NLAYR, NTBL
      INTEGER YR, MDATE, RUN, YRDOY, YRPLT, YRSIM, TIMDIF
      INTEGER YRDIF
      INTEGER, DIMENSION(NAPPL) :: IDLAPL, IRRCOD
      INTEGER, DIMENSION(NAPPL) :: JULAPL, JULWTB, JWTBRD

      REAL AIRAMT, AIRAMX, ATHETA, DEPIR, DSOIL, DSOILX
      REAL EFFIRR, EFFIRX, IRRAMT
      REAL SWDEF, THETAC, THETCX, TOTAPW, TOTEFFIRR, TOTIR
      REAL DLAYR(NL), DS(NL), DUL(NL), LL(NL), SW(NL)
      REAL, DIMENSION(NAPPL) :: AMIR, AMT, WTABL

!  Added for flooded field management
      LOGICAL PUDDLED
      INTEGER NBUND, NCOND, NPERC
      INTEGER IBDAT(NAPPL), IIRRCV(NAPPL), IPDAT(NAPPL) !, IIRRP(100)
      INTEGER CONDAT(NAPPL)   !, IIRRC(NAPPL)
      REAL BUND(NAPPL), IPERC(NAPPL), COND(NAPPL)
      REAL RAIN, IRRAPL, TIL_IRR

!  Added for drip irrigation
      INTEGER, DIMENSION(NAPPL) :: DripDat, DripNum
      INTEGER, DIMENSION(NAPPL) :: IRSTRH, IRSTRM, IRNUM
      REAL DripSpc, DripOfset
      REAL, DIMENSION(NAPPL) :: DripDur, DripInt, DripRate, DripStart 
      REAL, DIMENSION(NAPPL) :: IRDUR, IRINT
      LOGICAL AUTO   !indicator that automatic irrigation will be used
      LOGICAL DRIP   !indicator that drip irrigation will be used

!  Added for water table management
      REAL MgmtWTD, ICWD

!-----------------------------------------------------------------------
      TYPE (ControlType)  CONTROL
      TYPE (SwitchType)   ISWITCH
      TYPE (SoilType)     SOILPROP
      TYPE (FloodWatType) FLOODWAT

!     Transfer values from constructed data types into local variables.
      DYNAMIC = CONTROL % DYNAMIC
      YRDOY   = CONTROL % YRDOY

      DLAYR  = SOILPROP % DLAYR  
      DS     = SOILPROP % DS
      DUL    = SOILPROP % DUL    
      LL     = SOILPROP % LL     
      NLAYR  = SOILPROP % NLAYR  

      IIRRI  = ISWITCH % IIRRI

      PUDDLED= FLOODWAT % PUDDLED

C***********************************************************************
C***********************************************************************
C    Input and Initialization 
C***********************************************************************
      IF (DYNAMIC .EQ. INIT) THEN
C-----------------------------------------------------------------------
      FILEIO  = CONTROL % FILEIO
      LUNIO   = CONTROL % LUNIO
      MULTI   = CONTROL % MULTI
      YRDIF   = CONTROL % YRDIF
      RNMODE  = CONTROL % RNMODE
      RUN     = CONTROL % RUN
      YRSIM   = CONTROL % YRSIM

      ISWWAT  = ISWITCH % ISWWAT
      MEHYD   = ISWITCH % MEHYD

      TOTAPW = 0
      AMT    = 0. !irrigation amounts
      NAP    = 0
      NAPW   = 0  !irrigation application
      NBUND  = 0  !# bunds
      NTBL   = 0  !# water tables
      NCOND  = 0  !# irrigation applications (same as NAPW??)
      NPERC  = 0  !# percs
      NDRIP  = 0  !# drip irrigation days

      IRRAMT = 0.0
      TOTIR  = 0.
      TOTEFFIRR = 0.
      TIL_IRR = 0.0
      IRRAPL  = 0.0

      AUTO = .FALSE.

!     Water table depth (-99 indicates no water table present)
      MgmtWTD = -99.  

      IF (ISWWAT .EQ. 'Y') THEN
      !Data is read if not sequenced or seasonal run or for first
      !  season of sequence or seasonal runs.
!        IF ((INDEX('PQF',RNMODE) .EQ. 0 .OR. RUN .EQ. 1) 
!     &        .AND. MULTI .LE. 1) THEN

          JIRR = 0.0
          OPEN (LUNIO, FILE = FILEIO, STATUS = 'OLD', IOSTAT=ERRNUM)
          IF (ERRNUM .NE. 0) CALL ERROR(ERRKEY,ERRNUM,FILEIO,0)
          LNUM = 0

C-----------------------------------------------------------------------
!     May be needed for MgmtOps output later (OPOPS)
!      SELECT CASE (IIRRI)
!      CASE ('R')
!        IrrText = 'As reported in field schedule in YYDDD format'
!      CASE ('D')
!        IrrText = 'As reported in field schedule in days after ' //
!     &       'planting)'
!      CASE ('A')
!        IrrText = 'Automatic irrigation based on soil water deficit'
!      CASE ('F')
!        IrrText = 'Automatic irrigation with fixed amount (AIRAMT)'
!      CASE ('P')
!        IrrText = 'As reported through last reported day, then ' //
!     &       'automatic to re-fill profile (as in option A)'
!      CASE ('W')
!        IrrText = 'As reported through last reported day, then ' //
!     &       'automatic with fixed amount.'
!      CASE ('T')
!         IrrText = 'Repeat irrigation amount every day until a new record is input.'
!       CASE ('N')    !No irrigation
!      END SELECT

C-----------------------------------------------------------------------
C      Read Automatic Management
C-----------------------------------------------------------------------
          IF (INDEX('AFPW', ISWITCH % IIRRI) > 0) THEN
            SECTION = '!AUTOM'
            CALL FIND(LUNIO, SECTION, LINC, FOUND) ; LNUM = LNUM + LINC
            IF (FOUND .EQ. 0) CALL ERROR(SECTION, 42, FILEIO, LNUM)
            READ(LUNIO,'(/,14X,2(1X,F5.0),16X,I2,2(1X,F5.0))',
     &        IOSTAT=ERRNUM) DSOIL, THETAC, AIRRCOD, AIRAMT, EFFIRR
            LNUM = LNUM + 2
            IF (ERRNUM .NE. 0) CALL ERROR(ERRKEY,ERRNUM,FILEIO,LNUM)
          ENDIF

C-----------------------------------------------------------------------
C         Find and Read Initial Conditions Section
C-----------------------------------------------------------------------
          IF (INDEX('FQ',RNMODE) .LE. 0 .OR. RUN == 1) THEN
            REWIND(LUNIO)
            SECTION = '*INITI'
            CALL FIND(LUNIO, SECTION, LINC, FOUND) ; LNUM = LINC
            IF (FOUND .EQ. 0) THEN
              CALL ERROR(SECTION, 42, FILEIO, LNUM)
            ELSE
              READ(LUNIO,'(40X,F6.0)',IOSTAT=ERRNUM) ICWD ; LNUM =LNUM+1
              IF (ERRNUM .NE. 0) CALL ERROR(ERRKEY,ERRNUM,FILEIO,LNUM)
              MgmtWTD = ICWD
              CALL PUT('MGMT','WATTAB',MgmtWTD)
            ENDIF
          ENDIF

C-----------------------------------------------------------------------
!     Find and Read Planting Details Section
C-----------------------------------------------------------------------
          SECTION = '*PLANT'
          CALL FIND(LUNIO, SECTION, LINC, FOUND) ; LNUM = LNUM + LINC
          IF (FOUND .EQ. 0) CALL ERROR (SECTION, 42, FILEIO,LNUM)
          READ(LUNIO,'(35X,A1)', IOSTAT=ERRNUM) PLME
          LNUM = LNUM + 1
          IF (ERRNUM .NE. 0) CALL ERROR(ERRKEY,ERRNUM,FILEIO,LNUM)

C-----------------------------------------------------------------------
C    Find and Read Irrigation Section
C-----------------------------------------------------------------------
          SECTION = '*IRRIG'
          CALL FIND(LUNIO, SECTION, LINC, FOUND) ; LNUM = LNUM + LINC
          IF (FOUND .EQ. 0) CALL ERROR(SECTION, 42, FILEIO, LNUM)
          READ(LUNIO,'(2X,F6.3,2(1X,F5.0),18X,F6.1,2F6.0)'
     &      ,IOSTAT=ERRNUM)
     &      EFFIRX, DSOILX, THETCX, AIRAMX, DripSpc, DripOfset
          LNUM = LNUM + 1
          IF (ERRNUM .NE. 0) CALL ERROR(ERRKEY,ERRNUM,FILEIO,LNUM)
          JIRR = 0
          DO I = 1,NAPPL
!           READ(LUNIO,'(3X,I7,3X,I3,1X,F5.0,1X,I5)',IOSTAT=ERRNUM,
!     &        ERR=50)  IDLAPL(I), IRRCOD(I), AMT(I)   !, IIRRC(I)
            READ(LUNIO,'(3X,I7,3X,A90)',ERR=50, END=50) IDLAPL(I),CHAR
            LNUM = LNUM + 1

            READ(CHAR,'(I3,1X,F5.0)',IOSTAT=ERRNUM) IRRCOD(I), AMT(I)
            !READ(CHAR,'(I3,1X,F5.0,1X,I5)',IOSTAT=ERRNUM) 
            IF (ERRNUM .NE. 0) CALL ERROR(ERRKEY,ERRNUM,FILEIO,LNUM)

            IF (INDEX('GC',MEHYD) > 0 .AND. IRRCOD(I) == 5) THEN
              READ(CHAR,'(9X,2(1X,I2),2F6.0,I6)',IOSTAT=ERRNUM) 
     &          IRSTRH(I),IRSTRM(I),IRDUR(I),IRINT(I),IRNUM(I)
              IF (ERRNUM .NE. 0) CALL ERROR(ERRKEY,59,FILEIO,LNUM)

!             Check for small drip duration
              IF (IRDUR(I) < 1. .AND. IRNUM(I) > 0) THEN
                MSG(1) = "Drip irrigation duration less than 1 minute"
                WRITE(MSG(2),'(A,I4)') "Irrigation number ",I
                WRITE(MSG(3),'(I7,1X,A)') IDLAPL(I), CHAR(1:70)
                MSG(4) = "Program will stop."
                CALL WARNING(4,ERRKEY,MSG)
                CALL ERROR(ERRKEY,20,FILEIO,LNUM)
              ENDIF
            ENDIF

            JIRR = JIRR + 1
          ENDDO

   50     CONTINUE
          CLOSE (LUNIO)
!        ENDIF

C-----------------------------------------------------------------------
C     Set Irrigation Management
C-----------------------------------------------------------------------
        IF (IIRRI .EQ. 'R' .OR. IIRRI .EQ. 'D') THEN
          IF (EFFIRX .GT. 0) EFFIRR = EFFIRX
          IF (THETCX .GT. 0) THETAC = THETCX
          IF (DSOILX .GT. 0)  DSOIL = DSOILX
          IF (AIRAMX .GT. 0) AIRAMT = AIRAMX
        ENDIF
        
        IF (EFFIRR < 1.E-3) EFFIRR = 1.0

        IF (INDEX('AF',IIRRI) > 0) THEN
          AUTO = .TRUE.
        ENDIF

!     AMTMIN was not being used -- should it be used in place
!             of AIRAMT?  CHP
        !IF (AIRAMT .GT. 0.0) THEN   
        !  AMTMIN = AIRAMT
        !ELSE
        !  AMTMIN = 5.0
        !ENDIF

C       DSOIL in the CROPGRO model has units of cm
C       THETA has units of %
C
        IF (DSOIL .LE. 0.0) THEN
          DSOIL = (DLAYR(1) + DLAYR(2) + DLAYR(3))
        ENDIF
      ENDIF

!-----------------------------------------------------------------------
!     Sort time series records into water table and irrigation
!-----------------------------------------------------------------------
      CONDAT = 0
      IBDAT  = 0
      IIRRCV = 0
      IPDAT  = 0
      JULAPL = 0
      JWTBRD = 0
      JULWTB = 0
      DripDat= 0

      AMIR  = 0.0
      BUND  = 0.0
      COND  = 0.0
      IPERC = 0.0
      WTABL = -99.

      PUDDLED = .FALSE.
      DRIP = .FALSE.

!-----------------------------------------------------------------------
!     Irrigation Codes: IRRCOD
!         1:  Furrow irrigation of specified amount (mm)
!         2:  Alternating furrows; irrigation of specified amount (mm)
!         3:  Flood irrigation of specified amount (mm)
!         4:  Sprinkler irrigation of specified amount (mm)
!         5:  Drip or trickle irrigation of specified amount (mm)
!         6:  Single irrigation to specified total flood depth (mm)
!         7:  Water table depth (cm)
!         8:  Percolation rate (mm/d)
!         9:  Bund height (mm)
!        10:  Puddling (Puddled if IRRCOD = 10 record is present)
!        11:  Maintain constant specified flood depth (mm)

!     Transfer irrigation dates to appropriate arrays.
      IF (JIRR .GT. 0) THEN
        DO I = 1,JIRR
          SELECT CASE (IRRCOD(I))

          !------------------------------
           CASE (1:6)    
          !Regular irrigation (bunded or upland)

             IF (INDEX('GC',MEHYD) > 0 .AND. IRRCOD(I) == 5) THEN
                !------------------------------
                !Drip irrigation 
                NDRIP = NDRIP + 1
                DripDat(NDRIP) = IDLAPL(I)
                DripRate(NDRIP)= AMT(I)          !drip rate ml/s
                DripStart(NDRIP)= IRSTRH(I) + IRSTRM(I) / 60. !hours
                DripDur(NDRIP)  = IRDUR(I) / 60. !duration (hr)
                DripInt(NDRIP)  = IRINT(I) / 60. !interval (hr)
                DripNum(NDRIP)  = IRNUM(I)       !# per day

             ELSE
                NCOND         = NCOND + 1        
                CONDAT(NCOND) = IDLAPL(I)
                IIRRCV(NCOND) = IRRCOD(I)
                COND(NCOND)   = AMT(I)
             
                ! Regular irrigation upland fields
                NAPW = NAPW + 1
                JULAPL(NAPW) = IDLAPL(I)
                AMIR(NAPW)   = AMT(I)
             ENDIF

          !------------------------------
           CASE (7)
             !
             ! Water table
             !
             NTBL = NTBL + 1
!            JWTBRD stores original values
             JWTBRD(NTBL) = IDLAPL(I)  
!            JULWTB can be modified for sequenced or multi-year runs  
             JULWTB(NTBL) = IDLAPL(I)
             WTABL(NTBL)   = AMT(I)         !cm
!            PWAT(NTBL)    = AMT(I)/10.0

          !------------------------------
           CASE (8)
           ! Percolation

             NPERC         = NPERC + 1
             IPDAT(NPERC)  = IDLAPL(I)
             IPERC(NPERC)  = AMT(I)/10.0

          !------------------------------
           CASE (9)
             !
             ! Bunding
             !
             NBUND         = NBUND + 1
             IBDAT(NBUND)  = IDLAPL(I)
             BUND(NBUND)   = AMT(I)

          !------------------------------
           CASE (10)
             !
             ! Puddling for rice added
             !
             PUDDLED = .TRUE.

          !------------------------------
           CASE (11)
             !
             ! Maintain constant specified flood depth 
             !   until next irrigation record 
             NCOND         = NCOND + 1        
             CONDAT(NCOND) = IDLAPL(I)
             IIRRCV(NCOND) = IRRCOD(I)
             COND(NCOND)   = AMT(I)

         END SELECT
        ENDDO
      ENDIF

!!     Check for using flooded conditions with Century method of 
!!     soil organic matter.  
!      MESOM = ISWITCH % MESOM
!      IF (NBUND .GT. 0 .AND. MESOM .EQ. 'P') THEN
!!        MSG(1) = 'Flooded field not available with Century SOM model.'
!!        MSG(2) = 'Change either management or SOM method in FILEX.'
!!        MSG(3) = 'Program will stop.'
!
!        MSG(1) = 
!     &'Century model was not intended to be used with flooded fields.'
!        MSG(2) = 'Change either management or SOM method in FILEX.'
!        CALL WARNING(2, ERRKEY, MSG)
!!        WRITE(*,'(/,A78,/,A78,/,A78,/)') MSG(1), MSG(2), MSG(3)
!!        CALL ERROR(ERRKEY, 10, "", 0)
!      ENDIF
!

      IF (INDEX('GC',MEHYD) > 0 .AND. NDRIP .GT. 0) THEN
        DRIP = .TRUE.
      ENDIF

!     Transfer data to ModuleData
      CALL PUT('MGMT','TOTIR',TOTIR)
      CALL PUT('MGMT','EFFIRR',EFFIRR)
      CALL PUT('MGMT','IRRAMT',IRRAMT)
      CALL PUT('MGMT','DEPIR', DEPIR)
!-----------------------------------------------------------------------
!     Adjust irrigation dates for multi-year simulations
!     This section was taken from the MRUN subroutine.
!-----------------------------------------------------------------------
      IF (MULTI .GT. 1 .AND. JIRR .GT. 0 .AND. IIRRI .NE. 'D') THEN
        IF (NAPW .GT. 0 .AND. JULAPL(1) .LT. YRSIM) THEN
          DO I = 1, NAPW
            CALL YR_DOY(JULAPL(I),YR,IDATE)
            JULAPL(I) = (YR + MULTI - 1) * 1000 + IDATE
          ENDDO
        ENDIF

        IF (NDRIP .GT. 0 .AND. DripDat(1) .LT. YRSIM) THEN
          DO I = 1, NDRIP
            CALL YR_DOY(DripDat(I),YR,IDATE)
            DripDat(I) = (YR + MULTI - 1) * 1000 + IDATE
          ENDDO
        ENDIF

        IF (NTBL .GT. 0 .AND. JULWTB(1) .LT. YRSIM) THEN
          DO I = 1, NTBL
            CALL YR_DOY(JULWTB(I),YR,IDATE)
            JULWTB(I) = (YR + MULTI - 1) * 1000 + IDATE
          ENDDO
        ENDIF

        IF (NPERC .GT. 0 .AND. IPDAT(1) .LT. YRSIM) THEN
          DO I = 1, NPERC
            CALL YR_DOY(IPDAT(I),YR,IDATE)
            IPDAT(I) = (YR + MULTI - 1) * 1000 + IDATE
          ENDDO
        ENDIF

        IF (NBUND .GT. 0 .AND. IBDAT(1) .LT. YRSIM) THEN
          DO I = 1, NBUND
            CALL YR_DOY(IBDAT(I),YR,IDATE)
            IBDAT(I) = (YR + MULTI - 1) * 1000 + IDATE
          ENDDO
        ENDIF

        IF (NCOND .GT. 0 .AND. CONDAT(1) .LT. YRSIM) THEN
          DO I = 1, NCOND
            CALL YR_DOY(CONDAT(I),YR,IDATE)
            CONDAT(I) = (YR + MULTI - 1) * 1000 + IDATE
          ENDDO
        ENDIF
      ENDIF

C-----------------------------------------------------------------------
C     Adjust for crop rotations
C-----------------------------------------------------------------------
      !IF (RNMODE .EQ. 'Q') THEN
      IF (INDEX('PQ',RNMODE) .GT. 0) then 
        IF (NAPW .GT. 0 .AND. JULAPL(1) .LT. YRSIM .AND. IIRRI.NE.'D')
     &      THEN
          DO I = 1, NAPW
            CALL YR_DOY(JULAPL(I),YR,IDATE)
            JULAPL(I) = (YR + YRDIF) * 1000 + IDATE
          END DO
        ENDIF

        IF (NDRIP .GT. 0 .AND. DripDat(1) .LT. YRSIM .AND. IIRRI.NE.'D')
     &      THEN
          DO I = 1, NDRIP
            CALL YR_DOY(DripDat(I),YR,IDATE)
            DripDat(I) = (YR + YRDIF) * 1000 + IDATE
          END DO
        ENDIF

        IF (NTBL .GT. 0 .AND. JULWTB(1) .LT. YRSIM .AND. IIRRI.NE.'D')
     &      THEN
          DO I = 1, NTBL
            CALL YR_DOY(JULWTB(I),YR,IDATE)
            JULWTB(I) = (YR + YRDIF) * 1000 + IDATE
          END DO
        ENDIF

        IF (NPERC.GT. 0 .AND. IPDAT(1) .LT. YRSIM .AND. IIRRI.NE.'D')
     &      THEN
          DO I = 1, NPERC
            CALL YR_DOY(IPDAT(I),YR,IDATE)
            IPDAT(I) = (YR + YRDIF) * 1000 + IDATE
          END DO
        ENDIF

        IF (NBUND .GT. 0 .AND. IBDAT(1) .LT. YRSIM .AND. IIRRI.NE.'D')
     &      THEN
          DO I = 1, NBUND
            CALL YR_DOY(IBDAT(I),YR,IDATE)
            IBDAT(I) = (YR + YRDIF) * 1000 + IDATE
          END DO
        ENDIF

        IF (NCOND .GT. 0 .AND. CONDAT(1) .LT. YRSIM .AND. IIRRI.NE.'D')
     &      THEN
          DO I = 1, NCOND
            CALL YR_DOY(CONDAT(I),YR,IDATE)
            CONDAT(I) = (YR + YRDIF) * 1000 + IDATE
          END DO
        ENDIF
      ENDIF

!     Adjust IPERC values if necessary.
      DO I = 1, NPERC
        IF (DS(NLAYR) .LE. 20.0 .AND. IPERC(I) .GT. 0.) THEN
          IPERC(I) = IPERC(I) * 20.
        ENDIF
        IF (IPERC(I) .LE. 0.0) THEN
          IF (DS(NLAYR) .LE. 20.) THEN
            IPERC(I) = IPERC(1) * 20.0
          ELSE
            IPERC(I) = IPERC(1)
          ENDIF
        ENDIF
      END DO

!     Check for water table depths compatible with irrigation code
      IF (NTBL > 0) THEN
        SELECT CASE (IIRRI)
        CASE ('R','D'); CONTINUE
        CASE DEFAULT
          MSG(1)=
     &      "Water table depths require irrigation switch 'D' or 'R'."
          MSG(2)="Water table input records will be ignored."
          CALL WARNING(2,ERRKEY,MSG)
        END SELECT
      ENDIF

!      IF (NBUND .GT. 0) THEN
        CALL FLOOD_IRRIG (SEASINIT, 
     &    BUND, COND, CONDAT, IBDAT, IIRRCV, IIRRI,       !Input
     &    IPDAT, IPERC, NBUND, NCOND, NPERC,              !Input
     &    PUDDLED, RAIN, SOILPROP, SW, YRDOY, YRPLT,      !Input
     &    FLOODWAT,                                       !I/O
     &    DEPIR)                                          !Output
!      ENDIF

!     Store NBUND in composite variable. Used as a trigger for 
!       potential flooding.
      FLOODWAT % NBUND   = NBUND

      CALL DRIP_IRRIG (CONTROL,  
     &    AUTO, IRRAPL,  
     &    DripDat, DripDur, DripInt, DripNum, DripOfset,  !Input
     &    DripRate, DripSpc, DripStart, EFFIRR, IIRRI,    !Input
     &    MEHYD, NDRIP, YRDOY, YRPLT,                     !Input
     &    DEPIR)                                          !Output

      IF (INDEX('GC',MEHYD) > 0 .AND. NDRIP .GT. 0) THEN
        DRIP = .TRUE.
      ENDIF

!     Transfer data to ModuleData
      CALL PUT('MGMT','TOTIR',TOTIR)
      CALL PUT('MGMT','EFFIRR',EFFIRR)
      CALL PUT('MGMT','IRRAMT',IRRAMT)
      CALL PUT('MGMT','DEPIR', DEPIR)

!***********************************************************************
!***********************************************************************
!     DAILY RATE CALCULATIONS
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. RATE) THEN
C-----------------------------------------------------------------------
      DEPIR  = 0.
      IRRAMT = 0.
      IRRAPL = 0.0

!     Irrigation amount that affects soil properties after a tillage
!     event, expressed as equivalent rainfall depth (mm).
      TIL_IRR = 0.0

!-----------------------------------------------------------------------
!     Check to see if flood irrigation is done today
!-----------------------------------------------------------------------
      IF (NBUND .GT. 0) THEN
        CALL FLOOD_IRRIG (RATE, 
     &    BUND, COND, CONDAT, IBDAT, IIRRCV, IIRRI,       !Input
     &    IPDAT, IPERC, NBUND, NCOND, NPERC,              !Input
     &    PUDDLED, RAIN, SOILPROP, SW, YRDOY, YRPLT,      !Input
     &    FLOODWAT,                                       !I/O
     &    DEPIR)                                          !Output
        IF (DEPIR > 1.E-3) NAP = NAP + 1

!-----------------------------------------------------------------------
!     Now handle regular irrigation events
!-----------------------------------------------------------------------
      ELSE

      SELECT CASE (IIRRI)

C-----------------------------------------------------------------------
C** IIRRI = R - As Reported in FIELD SCHEDULE. If day of irrigation then.
C-----------------------------------------------------------------------
      CASE ('R')
        IF (NAPW .GT. 0) THEN
          LOOP1: DO I = 1, NAPW
            IF (JULAPL(I) .EQ. YRDOY) THEN
              DEPIR = DEPIR + AMIR(I)
              SELECT CASE(IRRCOD(I))
                CASE(1:4,6); TIL_IRR = TIL_IRR + AMIR(I)
              END SELECT
              NAP = NAP + 1
            ELSEIF (JULAPL(I) .GT. YRDOY) THEN
              EXIT LOOP1
            ENDIF
          END DO LOOP1
        ENDIF

C-----------------------------------------------------------------------
C** IIRRI = D - As Reported in FIELD SCHEDULE, days after planting
C-----------------------------------------------------------------------
      CASE ('D')
        IF (NAPW .GT. 0) THEN
          DAP = MAX(0,TIMDIF(YRPLT,YRDOY))
          LOOP2: DO I = 1, NAPW
            IF (JULAPL(I) .EQ. DAP) THEN
              DEPIR = DEPIR + AMIR(I)
              SELECT CASE(IRRCOD(I))
                CASE(1:4,6); TIL_IRR = TIL_IRR + AMIR(I)
              END SELECT
              NAP = NAP + 1
            ELSEIF (JULAPL(I) .GT. DAP) THEN
              EXIT LOOP2
            ENDIF
          END DO LOOP2
        ENDIF

C-----------------------------------------------------------------------
C** IIRRI = P - As Reported through last reported day (YYYYDDD), then automatic
C          to re-fill profile (as in option A)
C   IIRRI = W - As Reported through last reported day (YYYYDDD), then automatic
C**         adding AIRAMT each time
C-----------------------------------------------------------------------
      CASE ('P', 'W')

!       Not drip irrig
        IF (NAPW > 0) THEN
!       If Today's date is after the last record in the IRRIG section
!         check to see if automatic irrigation is needed (P or W option)
          IF (YRDOY .GT. JULAPL(NAPW))THEN
            AUTO = .TRUE.
          ELSE
            LOOP3: DO I = 1, NAPW
              IF (JULAPL(I) .EQ. YRDOY) THEN
                DEPIR = DEPIR + AMIR(I)
                NAP = NAP + 1
              ENDIF
            END DO LOOP3
          ENDIF

!       Drip irrigation
        ELSEIF (NDRIP > 0 .AND. YRDOY .GT. DripDat(NDRIP))THEN
          AUTO = .TRUE.
        ENDIF

      END SELECT
      ENDIF

!-----------------------------------------------------------------------
!   Automatic irrigation
!   AUTO = .TRUE.
!   IIRRI = A - Automatic irrigation or 
!           F - Fixed Amount Automatic Irrigation
!           P - after last reported date, refill profile
!           W - after last reported date, fixed amount
!-----------------------------------------------------------------------
      IF (AUTO) THEN
!       Only apply irrigation between planting and harvest
        IF ((YRDOY .GE. YRPLT .AND. YRDOY .LE. MDATE ).OR. 
     &      (YRDOY .GE. YRPLT .AND. MDATE .LE.  -99)) THEN

          IF (INDEX('GC',MEHYD) > 0) THEN
!           SWDEFICIT_bed looks for bedded system, which has a different
!           method for calculating volume of water deficit
            CALL SWDEFICIT_bed(
     &        DSOIL, DLAYR, DUL, LL, NLAYR, SW,           !Input
     &        ATHETA, SWDEF)                              !Output
          ELSE
            CALL SWDEFICIT(
     &        DSOIL, DLAYR, DUL, LL, NLAYR, SW,           !Input
     &        ATHETA, SWDEF)                              !Output
          ENDIF

          IF (ATHETA .LE. THETAC*0.01) THEN
!         A soil water deficit exists - automatic irrigation today.

!           IF (IIRRI .EQ. 'A') THEN
            IF (INDEX('AP',IIRRI) > 0) THEN
C             Determine supplemental irrigation amount.
C             Compensate for expected water loss due to soil evaporation
C             and transpiration today.
C             Estimate that an average of 5 mm of water will be lost.
              IRRAPL = SWDEF*10 + 5.0
              IRRAPL = MAX(0.,IRRAPL)

!           ELSE IF (IIRRI .EQ. 'F') THEN
            ELSE IF (INDEX('FW',IIRRI) > 0) THEN
C             Apply fixed irrigation amount
              IRRAPL = AIRAMT
            ENDIF

            SELECT CASE(AIRRCOD)
              CASE(1:4,6); TIL_IRR = TIL_IRR + IRRAPL
            END SELECT

            IF (.NOT. DRIP) THEN
              DEPIR = DEPIR + IRRAPL
              NAP = NAP + 1
              JULAPL(NAP) = YRDOY
              AMIR(NAP)   = IRRAPL
            ENDIF
          ENDIF
        ENDIF
      ENDIF

!-----------------------------------------------------------------------
!     Check to see if drip irrigation is done today
!-----------------------------------------------------------------------
      IF (DRIP) THEN
        CALL DRIP_IRRIG (CONTROL, 
     &    AUTO, IRRAPL,  
     &    DripDat, DripDur, DripInt, DripNum, DripOfset,  !Input
     &    DripRate, DripSpc, DripStart, EFFIRR, IIRRI,    !Input
     &    MEHYD, NDRIP, YRDOY, YRPLT,                     !Input
     &    DEPIR)                                          !Output

          IF (DEPIR > 1.E-4) THEN
            NAP = NAP + 1
          ENDIF
      ENDIF
!-----------------------------------------------------------------------
!     Water table management
      IF (NTBL .GT. 0) THEN
        DAP = MAX(0,TIMDIF(YRPLT,YRDOY))
        DO I = 1, NTBL
   !       IF ((IIRRI == 'R' .AND. YRDOY == JULWTB(I)) .OR.
   !  &        (IIRRI == 'D' .AND. DAP == JULWTB(I))) THEN
          IF ((IIRRI == 'R' .AND. YRDOY == JULWTB(I)) .OR.
     &        (IIRRI == 'D' .AND. DAP == JULWTB(I))
     &          .OR. ((IIRRI == 'A').AND. YRDOY == JULWTB(I))) THEN
            MgmtWTD =  WTABL(I) 
            CALL PUT('MGMT','WATTAB',MgmtWTD)
            EXIT
          ENDIF
        END DO
      ENDIF

C-----------------------------------------------------------------------
C    *********    IRRIGATE     **********
C-----------------------------------------------------------------------
      IF (EFFIRR .GT. 0.0) THEN
        IRRAMT = DEPIR*EFFIRR
      ELSE
        IRRAMT = DEPIR
      ENDIF
      IF (DEPIR .GT. 0.0) THEN
        !NAP    = NAP + 1
        TOTIR  = TOTIR + DEPIR
        TOTEFFIRR = TOTEFFIRR + IRRAMT
      ENDIF
!     Transfer data to ModuleData
      CALL PUT('MGMT','DEPIR', DEPIR)   !Total irrig amt today (mm)
      CALL PUT('MGMT','TOTIR', TOTIR)   !Total applied irrigAtion (mm)
      CALL PUT('MGMT','EFFIRR',EFFIRR)  !Effective irrigation %
      CALL PUT('MGMT','IRRAMT',IRRAMT)  !Effective irrig amt today (mm)
      
!***********************************************************************
!***********************************************************************
!     DAILY INTEGRATION CALCULATIONS
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. INTEGR) THEN
C-----------------------------------------------------------------------
    ! IF (DEPIR .GT. 0.0) THEN  ! Jin move to Daily Rate
        !!NAP    = NAP + 1
    !    TOTIR  = TOTIR + DEPIR
    !    TOTEFFIRR = TOTEFFIRR + IRRAMT
    !  ENDIF

!     Transfer data to ModuleData
 !     CALL PUT('MGMT','DEPIR', DEPIR)   !Total irrig amt today (mm)
 !     CALL PUT('MGMT','TOTIR', TOTIR)   !Total applied irrigAtion (mm)
 !     CALL PUT('MGMT','EFFIRR',EFFIRR)  !Effective irrigation %
 !     CALL PUT('MGMT','IRRAMT',IRRAMT)  !Effective irrig amt today (mm)

!***********************************************************************
      ENDIF
!***********************************************************************
!***********************************************************************
!     END OF DYNAMIC IF CONSTRUCT
!***********************************************************************
      FLOODWAT % PUDDLED = PUDDLED
      
      RETURN
      END SUBROUTINE IRRIG
C=======================================================================
!***********************************************************************
!***********************************************************************
!     IRRIG VARIABLE DEFINITIONS:     
!-----------------------------------------------------------------------
! AIRAMT    Amount of irrigation applied, if fixed, for automatic 
!             irrigation (mm)
! AIRAMX    Fixed irrigation amount (mm)
! AMIR(I)   Irrigation depth of the Ith application (mm)
! AMT(I)    Irrigation amount, depth of water/water table, bund height, or 
!             percolation rate  (mm or mm/day)
!           For Dripper: it is drip rate ml/s. Read from IRVAL of X file
! ATHETA    Average available water in top irrigation management depth of 
!             soil (%)
! DAP       Number of days after planting (d)
! DEPIR     Total irrigation depth for today (mm)
! DLAYR(L)  Soil thickness in layer L (cm)
! DripDur(NDRIP)  Dripper irrigation duration (hr)
! DripInt(NDRIP)  Dripper irrigation interval (hr)
! DripNum(NDRIP)  =        !# per day
! DripStart(NDRIP)= IRSTRH(I) + IRSTRM(I) / 60. !hours
! DSOIL     Irrigation management depth (cm)
! DSOILX    Management depth for automatic irrigation (cm)
! EFFIRR    Irrigation application efficiency (fraction) or  Effective irrigation fraction
! EFFIRX    Irrigation application efficiency reading from X file
! ERRKEY    Subroutine name for error file 
! ERRNUM    Error number for input 
! FILEIO    Filename for input file (e.g., IBSNAT35.INP) 
! FOUND     Indicator that good data was read from file by subroutine FIND 
!             (0 - End-of-file encountered, 1 - NAME was found) 
! IDATE     Day of irrigation or fertilizer application (d)
! IDLAPL(I) Irrigation or water table dates read from input file. 
! IIRRI     Irrigation switch R=on reported dates, D=as reported, days 
!             after planting, A=automatic, when rqd., F=automatic w/ fixed 
!             amt, P=as reported thru last reported day then automatic, 
!             W=as reported thru last reported day then fixed amount, N=not 
!             irrigated 
! IRDUR(I)   Irrigation duration
! IRINT(I)   Irrigation interval
! IRNUM(I)  
! IRRAMT    Irrigation amount (mm) or Effective irrig amt today (mm)
! IRRAPL    Supplemental irrigation amount. Compensate for expected water loss due to soil evaporation
!             and transpiration today.
! IRRCOD(I) Irrigation operation code: 1=Furrow, 2=Alternating furrows, 
!             3=flood, 4=Sprinkler, 5=Drip or trickle, 6=Flood depth, 
!             7=Water table depth, 8=Percolation rate, 9=Bund height
!             (mm or mm/day)
! IRSTRH(I) 
! IRSTRM(I) 
! JULAPL    Julian date for scheduled irrigation application (DOY)
! JULWTB    Julian date for scheduled irrigation application (YYDDD)
! JWTBRD    Recorded water table dates, saved for date modification on 
!             seasonal or sequenced runs 
! LNUM      Current line number of input file 
! LUNIO     Logical unit number for FILEIO 
! MULTI     Current simulation year (=1 for first or single simulation, 
!             =NYRS for last seasonal simulation) 
! NAP       Number of irrigation applications 
! NAPPL     Maximum number of applications (Here is refer to irrigation)
! NAPW      Number of irrigation values read from input file 
! JIRR      Number of irrigation records 
! NAPW       irrigation application
! NBUND      # bunds
! NCOND      # irrigation applications (same as NAPW??)
! NDRIP      # drip irrigation days    
! NL        Maximum number of soil layers = 20 
! NPERC      # percs
! NTBL      Number of water table values read 
! SECTION   Section name in input file 
! SWDEF     Soil water deficit (cm)
! THETAC    Threshold, % of maximum available water triggering irrigation
!             (%)
! THETCX    Threshold for automatic aplication, % of maximum avail water.
!             (%)
! TIL_IRR, the irrigation amount which affects 
!                 soil dynamics (excludes drip irrigation).   
! TOTAPW    Cumulative irrigation applied (mm)
! TOTIR     Total seasonal irrigation (mm)
!            Total applied irrigAtion (mm)
! WTABL     Water table record associated with JULWTB (cm)
! YR        Year portion of date 
! YRDIF     Function subroutine which calculates number of days between two 
!             dates (da)
! YRDOY     Current day of simulation (YYDDD)
! MDATE     Harvest maturity date (YYDDD)
! YRPLT     Planting date (YYDDD)
! YRSIM     Start of simulation date (YYDDD)
!-----------------------------------------------------------------------
!     END SUB-MODULE IRRIG
!-----------------------------------------------------------------------   
!=======================================================================

C=======================================================================
C  SWDEFICIT, Subroutine
C  Determines soil water deficit for automatic irrigation requirments
C-----------------------------------------------------------------------

      SUBROUTINE SWDEFICIT(
     &    DSOIL, DLAYR, DUL, LL, NLAYR, SW,               !Input
     &    ATHETA, SWDEF)                                  !Output

      USE ModuleDefs
      IMPLICIT NONE

      INTENT(IN) DSOIL, DLAYR, DUL, LL, NLAYR, SW
      INTENT(OUT) ATHETA, SWDEF

      INTEGER L, NLAYR
      REAL, DIMENSION(NL) :: DLAYR, DUL, LL, SW
      REAL ATHETA, DEPMAX, DSOIL, SWDEF, TSWTOP, WET1, XDEP, XDEPL

      WET1 = 0.0
      DEPMAX = 0.0
      TSWTOP = 0.0
      DEPMAX = 0.0

      DO L = 1,NLAYR
        IF (DEPMAX .LT. DSOIL) THEN
          XDEPL  = DEPMAX
          DEPMAX = DEPMAX + DLAYR(L)
          IF (DEPMAX .GT. DSOIL) THEN
            XDEP = (DSOIL - XDEPL)
          ELSE
            XDEP = DLAYR(L)
          ENDIF
          WET1 = WET1   + (DUL(L) - LL(L)) * XDEP
          TSWTOP = TSWTOP + (SW(L) - LL(L)) * XDEP
        ENDIF
      ENDDO

      ATHETA = TSWTOP / WET1
      SWDEF  = MAX(0.0,(WET1 - TSWTOP))

      RETURN
      END SUBROUTINE SWDEFICIT

C=======================================================================

