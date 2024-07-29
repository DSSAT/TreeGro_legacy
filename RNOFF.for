C=======================================================================
C  RNOFF, Subroutine, J.T. Ritchie
C  Calculate runoff by Williams-SCS curve number (CN) technique.
C-----------------------------------------------------------------------
C  REVISION HISTORY
C  01/01/1989 JR  Written
C  10/01/1987 GH  Correction of RUNOFF calc. according to SCS curve no.
C  12/05/1993 NBP Made into subroutine
C  07/12/1996 GH  Changed Precipitation to RAIN
C  07/01/1997 BDB Changed CN and SMX initialization (removed CN1, CN2,
C                 CN3, WF, WX, XX) from old INSOIL code
C  07/01/1997 BDB Simplified CN calculations (Removed C1, C2, C3, WF, DUL
C                  and made SMX input variables) (old WBSUBS code)
C  10/11/1997 CHP Updated for modular format.
C  09/01/1999  GH Incorporated into CROPGRO
!  06/12/2007 CHP Increase initial abstraction if mulch layer present.
!  08/19/2008 CHP/JTR Include Salus runoff option
!-----------------------------------------------------------------------
!  Called by: WATBAL
!  Calls:     None
C=======================================================================
      SUBROUTINE RNOFF(DYNAMIC,
     &    CN, LL, MEINF, MULCH, WEATHER,                  !Input
     &    SAT, SW, SWCN, TMAX, WATAVL,                    !Input
     &    RUNOFF)                                         !Output

C-----------------------------------------------------------------------
      USE ModuleDefs   ; USE MODULEDATA  !TEMP CHP
      IMPLICIT NONE
      SAVE

      INTEGER DYNAMIC

      CHARACTER*1 MEINF
      CHARACTER*5 Method
      CHARACTER*6 ERRKEY
      PARAMETER (ERRKEY = 'RNOFF')

      REAL CN, IABS
      REAL Ksat, Tmax
      REAL PB, WATAVL, SMX
      REAL RUNOFF, SWABI

      REAL RRC, D3, MeanPrecip, Offset, D1, D2, Slope

!     Maximum initial abstraction ratio
      REAL, PARAMETER :: MAXIABS = 0.6

!     Max number of soil layers defined in ModuleDefs.for
      REAL LL(NL), SAT(NL), SW(NL), SWCN(NL)

!     Mulch layer
      Type (MulchType) MULCH
      Type (WeatherType) WEATHER

!***********************************************************************
!***********************************************************************
!     SEASONAL INITIALIZATION
!***********************************************************************
      IF (DYNAMIC == SEASINIT) THEN
!-----------------------------------------------------------------------
      SELECT CASE (MEINF)
      CASE ('L')    !Salus runoff method
!       Default value for rainfall intenstiy factor
!       RainInt value can be read from Weather file (Header "RINT")
        Method = 'Salus'
!       Must have Ksat for this method to work
        IF (SWCN(1) < 0.0) THEN
          Method = 'SCSCN'
        ENDIF
      CASE Default !SCS Runoff curve number method 
        Method = 'SCSCN'
      END SELECT 

      IF (Method == 'Salus') THEN
!       Regional runoff coefficient, RRC
!       RRC is spatially variable, map available for USA prepared by Bruno and Joe
        RRC = WEATHER % RRC
        IF (RRC < 0.0) THEN
          RRC = 0.063
        ENDIF

!       Ranges of values for RRC based on JTR analysis
!       0.048 North Dakota
!       0.052 New Orleans
!       0.054 Lake Charles, LA
!       0.065-0.066 most FL locations
!       0.066 Jacksonville
!       0.072 Boston MA, Wilmington DE, ME
     
!       6/2/2010 JTR - D3 not dependant on Mean precip
!       No need for mean precip anymore
        D3 = 0.12

!!       MeanPrecip = WEATHER % MeanPrecip
!        SELECT CASE (NINT(MeanPrecip))
!        CASE (0:520); D3 = 0.24     !JR 6/2/2010
!        CASE DEFAULT; D3 = 0.12
!        END SELECT
      ENDIF
      
!***********************************************************************
!***********************************************************************
!     DAILY RATE CALCULATIONS
!***********************************************************************
      ELSEIF (DYNAMIC == RATE) THEN
!-----------------------------------------------------------------------
      IF (WATAVL < 1.E-3) THEN
        RUNOFF = 0.0
        RETURN
      ENDIF

      SELECT CASE (METHOD)

!     ----------------------------------------------
      CASE ('Salus')
        KSAT = SWCN(1)* 24.    !cm/d
      
!       Offset for runoff:rainfall line
        IF (TMAX < 10.) THEN
          Offset = -RRC * Ksat**0.76 * (-42.)         
        ELSEIF (TMAX < 52.) THEN
          Offset = -RRC * Ksat**0.76 * (TMAX - 52.)   
        ELSE
          Offset = 0.0  					        
        ENDIF
      
!       Slope of runoff:rainfall line
        D1 = exp(-0.015 * Ksat)                !Maximum slope 
        D2 = 0.04 * (1. - exp(-D3 * Ksat))     !slope of the slope 
        IF (TMAX < 30.) THEN
          Slope = MAX(D1 - D2 * (30. - TMAX), 0.0)
        ELSE
          Slope = D1
        ENDIF 
      
        IF (WATAVL > Offset) THEN
          RUNOFF = Slope * (WATAVL - Offset)   !mm/d
        ELSE
          RUNOFF= 0.0
        ENDIF

!     ----------------------------------------------
      CASE DEFAULT  !SCS CN method
        SMX = 254.0 * (100.0/CN - 1.0)
!       Initial abstraction ratio
!       Runoff is related to the average soil water content of the top
!       two layers of soil
        SWABI = 0.15 * ((SAT(1) - SW(1)) / (SAT(1) - LL(1) * 0.5) +
     &                (SAT(2) - SW(2)) / (SAT(2) - LL(2) * 0.5))
        SWABI = MAX(0.0, SWABI)

!       05/08/2006 CHP increase initial abstraction if surface mulch is
!       present. Initial abstraction ratio increases from SWABI
!       (no-mulch conditions) to 0.6 (MAXIABS) at full mulch coverage. 
!       IF (INDEX('RSN',MEINF) .LE. 0) THEN
        IF (INDEX('RSML',MEINF) > 0) THEN   
!         Model effects of mulch on runoff
          IABS = SWABI + (MAXIABS - SWABI) * MULCH % MULCHCOVER
          IABS = MAX(SWABI, IABS)
        ELSE
!         No mulch effects on runoff
          IABS = SWABI
        ENDIF
      
        PB = WATAVL - IABS * SMX
      
        IF (WATAVL .GT. 0.001) THEN
          IF (PB .GT. 0) THEN
            RUNOFF = PB**2/(WATAVL + (1.0-IABS) * SMX)
          ELSE
            RUNOFF = 0.0
          END IF
        ELSE
          RUNOFF = 0.0
        ENDIF
      END SELECT 

!***********************************************************************
      ENDIF
!***********************************************************************
      RETURN
      END SUBROUTINE RNOFF

!-----------------------------------------------------------------------
!     Variable definitions for RNOFF Module (updated 2/12/2004)
!-----------------------------------------------------------------------
! CN       Runoff Curve Number - measure of runoff potential based on soil 
!            type and current soil water content. 
! IABS     Initial abstraction ratio, modified for surface mulch layer effects.
! LL(L)    Volumetric soil water content in soil layer L at lower limit
!           (cm3 [water] / cm3 [soil])
! PB       Determines threshold amount of rainfall that will occur before 
!            runoff starts (mm/d)
! RUNOFF   Calculated runoff (mm/d)
! SAT(L)   Volumetric soil water content in layer L at saturation
!           (cm3 [water] / cm3 [soil])
! SMX      Soil storage available for surface water based on CN formula
!           (mm)
! SW(L)    Volumetric soil water content in layer L
!           (cm3 [water] / cm3 [soil])
! SWABI(L) A soil water abstraction index, a unitless indicator of the soil 
!            water condition at the time of a rainfall event.  This affects 
!            the intercept of the runoff axis when runoff starts to 
!            occur--later when drier and sooner when wetter.
! WATAVL   Water available for infiltration or runoff (rainfall plus 
!            irrigation) (mm/d)
!-----------------------------------------------------------------------
!     END SUBROUTINE RNOFF
!-----------------------------------------------------------------------


