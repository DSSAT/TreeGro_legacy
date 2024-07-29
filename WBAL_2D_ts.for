!=====================================================================
!  Wbal_2D_ts, Subroutine, Gerrit Hoogenboom
!  Modified for 2-D drip irrigation model
!  Seasonally: Provides output Water balance.  Prints file SoilWat_ts.OUT
!  Data is obtained from WATBAL, SPAM and IRRIG modules daily.  
!  Data from SPAM and IRRIG are sent via GETPUT routines.
!-----------------------------------------------------------------------
!  REVISION HISTORY
!  09/05/2008 CHP adapted WBAL for 2-D model 
!  08/25/2009 CHP modified for sub-daily time step
!  08/15/2011 Make detail available for cell(1,DripCol) and for cell(FurRow1,j), add INF_vol for cell detail
!             Add handling of LIMIT_2D for WBALAN
!             For 1st timestep, LatFlow include the portion which is calculated in WaterTable_2D
!-----------------------------------------------------------------------
!  Called by: WATBAL
!=====================================================================
      SUBROUTINE Wbal_2D_ts(CONTROL, ISWITCH, Time, TimeIncr,   !Input
     &    DRAIN, RUNOFF, IRRAMT, RAIN,                          !Input
     &    TES, TEP, TSW, CritCell, Diffus, Kunsat, LatFlow_ts,  !Input
     &    Count, LatFlow,                                       !Input
!         Temp chp
     &    CellArea, SWV_D, EP_vf, ES_vf_ts, IrrVol, INF_vol_dtal) !Input

!     ------------------------------------------------------------------
!      USE ModuleDefs 
      USE Cells_2D  !temp chp

      USE ModuleData
      IMPLICIT NONE
      SAVE

      CHARACTER*14, PARAMETER :: SWBAL = 'SoilWat_ts.OUT'
      INTEGER DAS, DYNAMIC, LUNWBL, I, Count
      INTEGER YRDOY
      INTEGER YR2, DY2, CritCell(2)

      REAL WBALAN, Time, TimeIncr, LatFlow_ts, LatFlow
      REAL CUMWBAL, Diffus1, Kunsat1

      REAL, DIMENSION(MaxRows,MaxCols) :: Kunsat, Diffus

      LOGICAL FEXIST, DOPRINT

      integer clun, detailRow, detailCol, DripCol
      real swij, swimj, swijm, swijp, swipj, swijcm2y, swijcm2
      real cbalij, rwuij, wbalij, esij
      REAL, DIMENSION(MaxRows,MaxCols) :: CellArea

      Double Precision DRAIN, IRRAMT, RAIN, RUNOFF
      Double Precision TEP, TES, TSW, TSWY, es_vf_ts
      Double Precision, DIMENSION(MaxRows,MaxCols) :: SWV_D, ep_vf
      Double Precision IrrVol, INF_vol_dtal
      TYPE (ControlType)  CONTROL
      TYPE (SwitchType)   ISWITCH

!     ------------------------------------------------------------------
      DYNAMIC = CONTROL % DYNAMIC
      YRDOY   = CONTROL % YRDOY
      DAS     = CONTROL % DAS
      DripCol = BedDimension % DripCol
!***********************************************************************
!***********************************************************************
!     Seasonal initialization - run once per season
!***********************************************************************
      IF (DYNAMIC .EQ. SEASINIT) THEN
!-----------------------------------------------------------------------
      DOPRINT=.TRUE.
      IF (ISWITCH % IDETW .EQ. 'N') THEN
        DOPRINT=.FALSE.
      ENDIF
      IF (ISWITCH % ISWWAT .EQ. 'N') THEN
        DOPRINT=.FALSE.
      ENDIF
      IF (ISWITCH % IDETL /= 'D') THEN
        DOPRINT=.FALSE.
      ENDIF
      IF (.NOT. DOPRINT) RETURN

!     Open output file SoilWat_ts.OUT
      CALL GETLUN(SWBAL, LUNWBL)
      INQUIRE (FILE = SWBAL, EXIST = FEXIST)
      IF (FEXIST) THEN
        OPEN (UNIT = LUNWBL, FILE = SWBAL, STATUS = 'OLD',
     &    POSITION = 'APPEND')
      ELSE
        OPEN (UNIT = LUNWBL, FILE = SWBAL, STATUS = 'NEW')
        WRITE(LUNWBL,'("*WATER BALANCE OUTPUT FILE")')
      ENDIF

      CALL HEADER(SEASINIT, LUNWBL, CONTROL % RUN)

!     Write header for daily output
      WRITE (LUNWBL,1120,ADVANCE='NO')
 1120 FORMAT('@YEAR DOY   DAS   TIME   INCR  Diffus  Kunsat  Row  Col',
     & '      SWTD',                               !State vars
     & '     IRRD     PRED     LAFD',              !Inflows
     & '     DRND     ROFD     ESAD     EPAD',     !Outflows
     & '     WBAL    CUMWBAL')                     !Balance

!     Soil water content for 1D simulations
      IF (NColsTot == 1) THEN
        WRITE(LUNWBL,1121) ("SW",I,"T",I=1,NRowsTot)
 1121   FORMAT(50(5X,A2,I2.2,A1))
      ELSE
        WRITE(LUNWBL,'(" ")')
      ENDIF

      TSWY   = TSW
      CUMWBAL = 0.0

      CALL YR_DOY(YRDOY, YR2, DY2)
      WRITE (LUNWBL,1300,ADVANCE='NO') YR2, DY2, DAS, Time, TimeIncr, 
     &    0.0, 0.0, 0, 0,
     &    TSW,                                     !State variables
     &    0.0, 0.0, 0.0,                           !Inflows
     &    0.0, 0.0, 0.0, 0.0,                      !Outflows
     &    0.0, CUMWBAL                             !Balance
      
!     Soil water content for 1D simulations
      IF (NColsTot == 1) THEN
          WRITE(LUNWBL,'(50F10.4)') (SWV_D(I,1),I=1,NRowsTot)
      ELSE
        WRITE(LUNWBL,'(" ")')
      ENDIF

!     ------------------------------------------------------------------
!     temp chp
!     water balance for single cell as defined in Cell_detail
      detailRow = Cell_detail%Row
      detailCol = Cell_detail%Col

      if (detailRow == 1 .and. detailCol == DripCol) then
      !   print *, "no cell detail available for cell (1,DripCol)"
      !   return
      endif

      CALL GETLUN("CellDetail.OUT", CLun)
      OPEN (UNIT = CLun, FILE = "CellDetail.OUT", STATUS = 'REPLACE')
      WRITE(CLun,'("*WATER BALANCE FOR CELL(",I2,",",I2,")")') 
     &      Cell_detail%row, Cell_detail%col
      CALL HEADER(SEASINIT, CLun, CONTROL % RUN)
      WRITE (CLun,1130)
 1130 FORMAT('@YEAR DOY   DAS   TIME   INCR',
     & '   SW(i,j)',                               !State vars
     & '      H_in      V_in',                     !Inflows
     & '     H_out     V_out     RWUij      ESij', !Outflows
     & '      WBAL   CUMWBAL',                     !Balance
     & ' SW(i-1,j) SW(i,j-1)   SW(i,j) SW(i,j+1) SW(i+1,j)', !Extra info
     & '  D(i-1,j)  D(i,j-1)    D(i,j)',           !Extra info
     & '  K(i-1,j)    K(i,j)     VoutD     VoutG') !Extra info

      SWimj = 0.0
      SWijm = 0.0
      SWijp = 0.0
      SWipj = 0.0

      SWijcm2  =SWV_D(detailRow,detailCol)*CellArea(detailRow,detailCol)   !cm2
      SWij     = SWV_D(detailRow,detailCol)
      IF (detailRow > 1)        SWimj = SWV_D(detailRow-1,detailCol)
      IF (detailCol > 1)        SWijm = SWV_D(detailRow,detailCol-1)
      IF (detailCol < NColsTot) SWijp = SWV_D(detailRow,detailCol+1)
      IF (detailRow < NRowsTot) SWipj = SWV_D(detailRow+1,detailCol)

      WRITE (CLUN,1320) YR2, DY2, DAS, 0.0, 0.0, 
     &      SWijcm2,                       !State variables
     &      0.0, 0.0,                   !Inflows
     &      0.0, 0.0, 0.0, 0.0,         !Outflows
     &      0.0, 0.0,                   !Balance
     &      SWimj, SWijm, SWij, SWijp, SWipj, !Extras
     &      0.0, 0.0, 0.0,              !Extras
     &      0.0, 0.0, 0.0, 0.0

      SWijcm2y = SWijcm2
      CBALij = 0.0

!***********************************************************************
!***********************************************************************
!     DAILY OUTPUT 
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. RATE) THEN
!-----------------------------------------------------------------------
      IF (.NOT. DOPRINT) RETURN
!       Change in storage = Inflows - Outflows
!       Balance = Inflows - Outflows - Change in storage
        if (BedDimension % LIMIT_2D .GE. NRowsTot) then 
          WBALAN = 
     &         + IRRAMT + RAIN + LatFlow_ts      !Inflows
     &         - DRAIN - RUNOFF - TES - TEP   !Outflows
     &         - (TSW - TSWY)                 !Change in soil water 
        else ! Drain is part of LatFlow_ts
          WBALAN = 
     &         + IRRAMT + RAIN + LatFlow_ts      !Inflows
     &         - RUNOFF - TES - TEP   !Outflows
     &         - (TSW - TSWY)                 !Change in soil water 
        Endif
        ! for 1st timestep, LatFlow include the portion which is calculated in WaterTable_2D
        If (Count .eq. 1)  WBALAN =  WBALAN - LatFlow_ts + LatFlow
        CUMWBAL = CUMWBAL + WBALAN

        IF (CritCell(1) > 0 .and. CritCell(1) <= NRowsTot .and. 
     &      CritCell(2) > 0 .and. CritCell(2) <= NColsTot) THEN
          Diffus1 = Diffus(CritCell(1),CritCell(2))
          Kunsat1 = Kunsat(CritCell(1),CritCell(2))
        ELSE
          Diffus1 = Diffus(1,1)
          Kunsat1 = Kunsat(1,1)
        ENDIF

        CALL YR_DOY(YRDOY, YR2, DY2)
        if (Count .eq. 1) then 
          WRITE (LUNWBL,1300,ADVANCE='NO') YR2, DY2, DAS, Time,TimeIncr,
     &      Diffus1, Kunsat1, 
     &      CritCell(1), CritCell(2),
     &      TSW,                                        !State variables
     &      IRRAMT, RAIN, LatFlow,                      !Inflows
     &      DRAIN, RUNOFF, TES, TEP,                    !Outflows
     &      WBALAN, CUMWBAL                             !Balance
        else
          WRITE (LUNWBL,1300,ADVANCE='NO') YR2, DY2, DAS, Time,TimeIncr,
     &      Diffus1, Kunsat1, 
     &      CritCell(1), CritCell(2),
     &      TSW,                                        !State variables
     &      IRRAMT, RAIN, LatFlow_ts,                      !Inflows
     &      DRAIN, RUNOFF, TES, TEP,                    !Outflows
     &      WBALAN, CUMWBAL                             !Balance 
        endif
 1300  FORMAT
     &    (1X,I4,1X,I3.3,1X,I5,2F7.3,   !Time
     &    2F8.1,                        !D, K
     &    2I5,                          !CritCells
     &    F10.4,                        !TSW
     &    3F9.4,                        !Inflows
     &    4F9.4,                        !Outflows
     &    F9.4,F11.4)                   !Balances, TSRadFrac

!       Soil water content for 1D simulations
        IF (NColsTot == 1) THEN
          WRITE(LUNWBL,'(50F10.4)') (SWV_D(I,1),I=1,NRowsTot)
        ELSE
          WRITE(LUNWBL,'(" ")')
        ENDIF

        !Save values for comparison tomorrow
        TSWY   = TSW

!     --------------------------------------------------------------------
!     temp chp
!     Water balance detail for one cell in cm2
      if (detailRow == 1 .and. detailCol==DripCol) then 
         Cell_detail%v_in = IrrVol
      !SWV_avail(1,DripCol) = SWV_avail(1,DripCol) + IrrVol
    ! &                                          / CellArea(1,DripCol)
      endif
      if (detailRow ==BedDimension % FurRow1) then 
!           Within furrow, add infiltration to top cells
        Cell_detail%v_in = INF_vol_dtal * CellArea(detailRow, detailCol)
            !INF_vol = WINF_col(j) * 0.1 * DayIncr / Thick(FurRow1,j)
!           cm3[water]   mm[water]   cm         1           
!           ---------- = --------- * -- * d * --------
!            cm3[soil]       d       mm       cm[soil]   
        
           !! SWV_avail(FurRow1,j) = SWV_avail(FurRow1,j) + INF_vol
      Endif
      if (BedDimension % LIMIT_2D . LT. detailRow) then
        Cell_detail%v_in = -999999
      Endif
      if (BedDimension % LIMIT_2D . EQ. detailRow) then
       ! Cell_detail%v_out = -999999
      Endif
      SWimj = 0.0
      SWijm = 0.0
      SWijp = 0.0
      SWipj = 0.0

      SWijcm2  =SWV_D(detailRow,detailCol)*CellArea(detailRow,detailCol)   !cm2
      SWij     = SWV_D(detailRow,detailCol) !cm3/cm3
      IF (detailRow > 1)        SWimj = SWV_D(detailRow-1,detailCol) 
      IF (detailCol > 1)        SWijm = SWV_D(detailRow,detailCol-1) 
      IF (detailCol < NColsTot) SWijp = SWV_D(detailRow,detailCol+1) 
      IF (detailRow < NRowsTot) SWipj = SWV_D(detailRow+1,detailCol) 

      RWUij = ep_vf(detailRow,detailCol) * CellArea(detailRow,detailCol)
      ESij = es_vf_ts * CellArea(detailRow,detailCol)
      if (abs(Cell_detail%v_in - -999999) .LT. 0.00001) then 
        WBALij = -999999
      elseif (abs(Cell_detail%v_out - -999999) .LT. 0.00001) then 
        WBALij = -999999
      else 
        WBALij = 
     &     + Cell_detail%h_in + Cell_detail%v_in     !Inflows
     &     - Cell_detail%h_out - Cell_detail%v_out   !Outflows
     &     - RWUij - ESij                            !Outflows
     &     - (SWijcm2 - SWijcm2y)                    !Change in soil water
        CBALij = CBALij + WBALij 
      endif

      if (detailCol >1 .and. detailRow > 1) then
        WRITE (CLUN,1320) YR2, DY2, DAS, Time, TimeIncr, 
     &      SWijcm2,                                      !State vars
     &      Cell_detail%h_in, Cell_detail%v_in,           !Inflows
     &      Cell_detail%h_out, Cell_detail%v_out,         !Outflows
     &      RWUij, ESij,                                  !Outflows
     &      WBALij, CBALij,                               !Balance
     &      SWimj, SWijm, SWij, SWijp, SWipj, 
     &      Diffus(detailRow-1,detailCol), Diffus(detailRow,detailCol-1)
     &      , Diffus(detailRow,detailCol), 
     &      Kunsat(detailRow-1,detailCol), Kunsat(detailRow,detailCol), 
     &      Cell_detail%vdiff, Cell_detail%vgrav
 1320   FORMAT(1X,I4,1X,I3.3,1X,I5,2F7.3,14F10.6,3F10.2,4F10.6)
      elseif (detailRow > 1) then
        WRITE (CLUN,1320) YR2, DY2, DAS, Time, TimeIncr, 
     &      SWijcm2,                                      !State vars
     &      Cell_detail%h_in, Cell_detail%v_in,           !Inflows
     &      Cell_detail%h_out, Cell_detail%v_out,         !Outflows
     &      RWUij, ESij,                                  !Outflows
     &      WBALij, CBALij,                               !Balance
     &      SWimj, SWijm, SWij, SWijp, SWipj, 
     &      Diffus(detailRow-1,detailCol), 0.0, 
     &      Diffus(detailRow,detailCol), 
     &      Kunsat(detailRow-1,detailCol), Kunsat(detailRow,detailCol), 
     &      Cell_detail%vdiff, Cell_detail%vgrav
      elseif (detailCol >1) then !DetailRow=1
        WRITE (CLUN,1320) YR2, DY2, DAS, Time, TimeIncr, 
     &      SWijcm2,                                      !State vars
     &      Cell_detail%h_in, Cell_detail%v_in,           !Inflows
     &      Cell_detail%h_out, Cell_detail%v_out,         !Outflows
     &      RWUij, ESij,                                  !Outflows
     &      WBALij, CBALij,                               !Balance
     &      SWimj, SWijm, SWij, SWijp, SWipj, 
     &      0.0, Diffus(detailRow,detailCol-1), 
     &      Diffus(detailRow,detailCol), 
     &      0.0, Kunsat(detailRow,detailCol), 
     &      Cell_detail%vdiff, Cell_detail%vgrav
      else !DetailRow=1 & DetailCol=1
        WRITE (CLUN,1320) YR2, DY2, DAS, Time, TimeIncr, 
     &      SWijcm2,                                      !State vars
     &      Cell_detail%h_in, Cell_detail%v_in,           !Inflows
     &      Cell_detail%h_out, Cell_detail%v_out,         !Outflows
     &      RWUij, ESij,                                  !Outflows
     &      WBALij, CBALij,                               !Balance
     &      SWimj, SWijm, SWij, SWijp, SWipj, 
     &      0.0, 0.0, 
     &      Diffus(detailRow,detailCol), 
     &      0.0, Kunsat(detailRow,detailCol), 
     &      Cell_detail%vdiff, Cell_detail%vgrav
      endif

      SWijcm2y = SWijcm2

!***********************************************************************
!***********************************************************************
!     SEASEND - Seasonal output
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. SEASEND) THEN
C-----------------------------------------------------------------------
      IF (.NOT. DOPRINT) RETURN

      CLOSE(LUNWBL)    
      close(clun)   

!***********************************************************************
!***********************************************************************
!     END OF DYNAMIC IF CONSTRUCT
!***********************************************************************
      ENDIF
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE Wbal_2D_ts
C=======================================================================
C=====================================================================
!     Wbal_2D_ts VARIABLE DEFINITIONS:
!-----------------------------------------------------------------------
! CritCell(2) The # of row and # of Col of cell which require smallest time step
! Diffus  Diffusivity
! DRAIN   Drain water from the bottom of simulation depth
! EP_vf(Row,Col) In a time step, cell actual plant transpiration rate volumetric fraction in mm/mm
! ES_vf_ts(Row,Col) In a time step, cell soil evaporation volumetric fraction in mm/mm
! H_in(Row,Col)    2D Horizontal water amount into cell in current time step in cm2
! H_out(Row,Col)   2D Horizontal water amount out of cell in current time step in cm2
! INF_vol_dtal INF_vol for cell detail
! IRRAMT  Irrigation amount (mm) 
! RAIN
! RUNOFF   Run off water in mm
! SWijcm2  Total cell water in cm2
! TEP     Total potential root uptake water in current time step
! TES     Total soil evaporation in current time step(mm)
! Time    Time of the day in hour
! TimeIncr length of time step
! TSW     Total soil water in profile (cm?)
!-----------------------------------------------------------------------
!     END SUBROUTINE Wbal_2D_ts
!=======================================================================
