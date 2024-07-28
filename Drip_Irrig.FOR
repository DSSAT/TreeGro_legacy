!=======================================================================
!  Drip_IRRIG, Subroutine
!
!  Determines drip irrigation - works with MGAR water balance routine
!-----------------------------------------------------------------------
!  Revision history
!
!  09/04/2008 CHP Written
!  07/27/2010 CHP Drip irrigation emitter can be offset from centerline.
!  03/30/2011 CHP Add Automatic Management for IIRRI
!=======================================================================

      SUBROUTINE Drip_IRRIG (CONTROL,  
     &    AUTO, IRRAPL,  
     &    DripDat, DripDur, DripInt, DripNum, DripOfset,  !Input
     &    DripRate, DripSpc, DripStart, EFFIRR, IIRRI,    !Input
     &    MEHYD, NDRIP, YRDOY, YRPLT,                     !Input
     &    DEPIR)                                          !Output

      USE ModuleDefs
      USE ModuleData
      IMPLICIT NONE
      SAVE

      CHARACTER*1 IIRRI, MEHYD
      CHARACTER*7, PARAMETER :: ERRKEY = 'DRIPIRR'
      INTEGER DAP, DYNAMIC, I
      INTEGER NDRIP, TIMDIF, YRDOY, YRPLT
      INTEGER, DIMENSION(NAPPL) :: DripDat, DripNum
      REAL DEPIR, DripSpc, DripOfset, EFFIRR, IrrRate, ROWSPC_CM, SPD
      REAL, DIMENSION(NAPPL) :: DripInt, DripRate, DripDur, DripStart
      REAL IRRAPL

      CHARACTER*6 SECTION
      INTEGER ERR, FOUND, LINC, LNUM, LUNIO
      LOGICAL AUTO

      TYPE (DripIrrType) DripIrrig
      TYPE (ControlType) CONTROL

      DYNAMIC = CONTROL % DYNAMIC

!***********************************************************************
!***********************************************************************
!     Seasonal initialization - run once per season
!***********************************************************************
      IF (DYNAMIC .EQ. RUNINIT .OR. DYNAMIC .EQ. SEASINIT) THEN
C-----------------------------------------------------------------------
      DripIrrig % DripRate  = 0.0
      DripIrrig % DripNum   = 0
      DripIrrig % DripStart = 0.0  
      DripIrrig % DripDur   = 0.0    
      DripIrrig % DripInt   = 0.0    
      DripIrrig % DripSpc   = DripSpc
      DripIrrig % NDrip     = NDrip
      DripIrrig % DripOfset = DripOfset
      CALL PUT(DripIrrig)

      IrrRate = 0.0
      DEPIR = 0.0

      IF (INDEX('GC',MEHYD) < 1 .OR. NDRIP < 1) RETURN

!-----------------------------------------------------------------------
!    Get row spacing from FILEIO
      LUNIO = CONTROL % LUNIO
      OPEN (LUNIO, FILE = CONTROL%FILEIO,STATUS = 'OLD',IOSTAT=ERR)
      IF (ERR .NE. 0) CALL ERROR(ERRKEY,ERR,CONTROL%FILEIO,0)
      LNUM = 0
      SECTION = '*PLANT'
      CALL FIND(LUNIO, SECTION, LINC, FOUND) ; LNUM = LNUM + LINC
      IF (FOUND == 0) THEN
        CALL ERROR(SECTION, 42, CONTROL%FILEIO, LNUM)
      ELSE
        READ(LUNIO,'(42X,F6.0,12X,F6.0)',IOSTAT=ERR)ROWSPC_CM
        LNUM = LNUM + 1
        IF (ERR .NE. 0) CALL ERROR(ERRKEY,ERR,CONTROL%FILEIO,LNUM)
      ENDIF

!***********************************************************************
!***********************************************************************
!     DAILY RATE CALCULATIONS
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. RATE) THEN
      DAP = MAX(0,TIMDIF(YRPLT,YRDOY))

      IF (AUTO) THEN
        DripIrrig % DripRate  = DripRate(NDRIP) * EFFIRR
        DripIrrig % DripStart = DripStart(NDRIP)  
        DripIrrig % DripInt   = DripInt(NDRIP) 
   
        IF (IRRAPL > 1.E-6) THEN
!         Calculate drip duration based on known quantity and drip rate
          DripIrrig % DripDur = 
     &      IRRAPL / DripRate(NDRIP) * DripSpc * ROWSPC_CM / 10. / 3600.
!         hr          mm/d        cm/emitter * cm
!         -- =    ------------- * ---------------
!          d      cm3/s-emitter     mm/cm * s/hr

          IrrRate = IRRAPL
          DripIrrig % IrrRate = IRRAPL    
          DripIrrig % DripNum   = 1

        ELSE
          DripIrrig % IrrRate = 0.0    
          DripIrrig % DripNum   = 0
          DripIrrig % DripDur = 0.0
          IrrRate = 0.0
        ENDIF

      ELSE
!-----------------------------------------------------------------------
!** IIRRI = R - As Reported in FIELD SCHEDULE. If day of irrigation then.
!** IIRRI = D - As Reported in FIELD SCHEDULE, days after planting
!** IIRRI = P - As Reported through last reported day (YYYYDDD), then automatic
!               to re-fill profile (as in option A)
!** IIRRI = W - As Reported through last reported day (YYYYDDD), then automatic
!               adding AIRAMT each time
!-----------------------------------------------------------------------
        DO I = 1, NDRIP
          IF ((INDEX('RPW',IIRRI) > 0 .AND. DripDat(I) .EQ. YRDOY) .OR.
     &        (INDEX('D',IIRRI)   > 0 .AND. DripDat(I) .EQ. DAP)) THEN
      
!           New drip irrigation schedule today
            DripIrrig % DripRate  = DripRate(I) * EFFIRR
            DripIrrig % DripNum   = DripNum(I)
            DripIrrig % DripStart = DripStart(I)  
            DripIrrig % DripDur   = DripDur(I)    
            DripIrrig % DripInt   = DripInt(I)    
      
            SPD = DripNum(I) * DripDur(I) * 3600. 
!           s/d =   #/day    *   hr       * s/hr 
      
            IrrRate = DripRate(I) / (DripSpc * ROWSPC_CM) * 10. * SPD
!                        cm3/s              emitter         mm     s
!             mm/d  =   -------   /        -------        * --  *  -
!                       emitter             cm * cm         cm     d
      
            DripIrrig % IrrRate   = IrrRate    
            EXIT
          ENDIF
        ENDDO
      ENDIF
     
!     Drip irrigation schedule maintained until new schedule specified
      DEPIR = MAX(0.0, IrrRate)

!***********************************************************************
      ENDIF
!***********************************************************************
!***********************************************************************
!     END OF DYNAMIC IF CONSTRUCT
!***********************************************************************

      CALL PUT(DripIrrig)

      RETURN
      END SUBROUTINE Drip_IRRIG

!=======================================================================

C=======================================================================
C  SWDEFICIT_bed, Subroutine
C  Determines soil water deficit for automatic irrigation requirments
!    for bedded system
C-----------------------------------------------------------------------

      SUBROUTINE SWDEFICIT_bed(
     &    DSOIL, DLAYR, DUL, LL, NLAYR, SW,               !Input
     &    ATHETA, SWDEF)                                  !Output

      USE ModuleDefs
      USE Cells_2D
      IMPLICIT NONE

      INTENT(IN) DSOIL, DLAYR, DUL, LL, NLAYR, SW
      INTENT(OUT) ATHETA, SWDEF

      INTEGER L, NLAYR
      REAL, DIMENSION(NL) :: DLAYR, DUL, LL, SW
      REAL ATHETA, DEPMAX, DSOIL, SWDEF, TSWTOP, WET1, XDEP, XDEPL
      REAL MULT, MULT_bed

      IF (BedDimension % RaisedBed) THEN
!       If raised bed, then volume of water deficit is smaller due to
!       presence of furrow.
        MULT_bed = BedDimension % BEDWD / BedDimension % ROWSPC_cm

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
            IF (L < BedDimension % FurRow1) THEN
              MULT = MULT_bed
            ELSE
              MULT = 1.0
            ENDIF
            WET1   = WET1   + (DUL(L) - LL(L)) * XDEP * MULT
            TSWTOP = TSWTOP + (SW(L)  - LL(L)) * XDEP * MULT
          ENDIF
        ENDDO
      
        ATHETA = TSWTOP / WET1
        SWDEF  = MAX(0.0,(WET1 - TSWTOP))

      ELSE
!       If no raised bed, just call the original routine.
        CALL SWDEFICIT(
     &    DSOIL, DLAYR, DUL, LL, NLAYR, SW,               !Input
     &    ATHETA, SWDEF)                                  !Output
      ENDIF

      RETURN
      END SUBROUTINE SWDEFICIT_bed

C=======================================================================


