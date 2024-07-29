!=======================================================================
C  SoilNiBal_2D, Subroutine
C 
C  Purpose: Provides seasonal inorganic soil N balance.  
C     Based on SoilNBal.FOR.
C
C  REVISION   HISTORY
C  03/04/2005 CHP wrote based on SoilNBal
!  08/15/2011 floodN is comment out
!=======================================================================

      SUBROUTINE SoilNiBal_2D (CONTROL, ISWITCH, Cells, DENITRIF, IMM,  
     &    MNR, ALGFIX, CIMMOBN, CMINERN, CUMFNRO, FERTDATA, NBUND, TLCH,
     &    TNH4, TNO3, TNOX, TOTAML, TUREA, WTNUP) 
!     ------------------------------------------------------------------
      USE Cells_2D
      USE ModuleDefs     !Definitions of constructed variable types, 
                         !which contain control information, soil
                         !parameters, hourly weather data.
      USE Interface_SoilNBalSum
      IMPLICIT NONE
      SAVE
!     ------------------------------------------------------------------
      LOGICAL FEXIST

      CHARACTER*1  IDETN, IDETL, ISWNIT
      CHARACTER*13, PARAMETER :: SNiBAL = 'SoilNiBal.OUT'

      INTEGER DAS, DYNAMIC, INCDAT, YRDOY
      INTEGER YRSIM, RUN, LUNSNC, NBUND, Clunn
      INTEGER YR, DOY, YRI, DOYI, detailRow, detailCol, DripCol

      REAL ALGFIX, ALGFIXI, AMTFER, TALLN, TALLNI, TLCH, TNH4, TNH4I,
     &  TNO3, TNO3I,  TNOX, TUREA, TUREAI, WTNUP
      REAL TOTAML, CUMFNRO !, TOTFLOODN, TOTFLOODNI
      REAL STATEN, BALANCE
      REAL DENITRIF(MaxRows,MaxCols)
      REAL LCHTODAY, NOXTODAY, IMMOBTODAY, MINERTODAY
      REAL WTNUPTODAY, AMLTODAY, FNROTODAY, AMTFERTODAY
      REAL TLCHY, TNOXY, WTNUPY, CIMMOBY, CMINERY
      REAL TOTAMLY, CUMFNROY, AMTFERY
      REAL TOTSTATE, TOTADD, TOTSUB, DAYBAL, TOTSTATY, CUMBAL 
      REAL CIMMOBN, CMINERN, CelNtot, CelWD, CelHT, CelNBal, CCelNBal
      REAL YCelNBal, CelNtotI, CumNH4U, CumNO3U, FertAdd
      REAL IMM(0:NL,NELEM), MNR(0:NL,NELEM), CelMNR, CelIMM
!     ------------------------------------------------------------------
      TYPE (ControlType) CONTROL
      TYPE (SwitchType)  ISWITCH
      TYPE (FertType)    FertData
      TYPE (CellType) CELLS(MaxRows,MaxCols), CellDetail
!     ------------------------------------------------------------------
!     Return if detail not requested.
      IDETL   = ISWITCH % IDETL
      IDETN   = ISWITCH % IDETN
      ISWNIT  = ISWITCH % ISWNIT
      IF (IDETL  == 'N' .OR. 
     &    IDETL  == '0' .OR.    !zero
     &    IDETN  == 'N' .OR. 
     &    ISWNIT == 'N') RETURN

!     Transfer values from constructed data types into local variables.
      DAS     = CONTROL % DAS
      DYNAMIC = CONTROL % DYNAMIC
      RUN     = CONTROL % RUN
      YRDOY   = CONTROL % YRDOY
      YRSIM   = CONTROL % YRSIM

      AMTFER = FERTDATA % AMTFER(N)
      detailRow = Cell_detail%Row
      detailCol = Cell_detail%Col
      CellDetail = Cells(detailRow, detailCol)
      CelWD =  CellDetail % STRUC % WIDTH
      CelHT  = CellDetail % STRUC % THICK
!***********************************************************************
!***********************************************************************
!     Seasonal Initialization phase
!***********************************************************************
      IF (DYNAMIC == SEASINIT) THEN
!     ------------------------------------------------------------------
!     Initialize output file
      CALL GETLUN(SNiBAL, LUNSNC)
      INQUIRE (FILE = SNiBAL, EXIST = FEXIST)
      IF (FEXIST) THEN
        OPEN (UNIT = LUNSNC, FILE = SNiBAL, STATUS = 'OLD',
     &    POSITION = 'APPEND')
        WRITE(LUNSNC,'(/,"!",79("="))') 
      ELSE
        OPEN (UNIT = LUNSNC, FILE = SNiBAL, STATUS = 'NEW')
        WRITE(LUNSNC,'("*SOIL INORGANIC N BALANCE")')
      ENDIF

      CALL HEADER(SEASINIT, LUNSNC, RUN)

!     Initial value for extractable N summed across all soil layers.
      TNO3I  = TNO3
      TNH4I  = TNH4
      TUREAI = TUREA
      CelNtotI = CelNtot
      CumNH4U = 0.
      CumNO3U = 0.

!     Sum the initial value of all abiotic N pools (soil, air)
      TALLNI = TNO3I + TNH4I

   !   TOTFLOODNI = TOTFLOODN
      ALGFIXI = ALGFIX

!     If detailed printout requested, print daily soil N balance
      IF (IDETL == 'D') THEN
!       Cumulative values (yesterday)
!       Save today's cumulative values for use tomorrow
        TLCHY    = 0.0
        TNOXY    = 0.0
        WTNUPY   = 0.0
        TOTAMLY  = 0.0
        CUMFNROY = 0.0
        AMTFERY  = 0.0
        CMINERY = 0.0
        CIMMOBY = 0.0
        TOTSTATY = TNO3I + TNH4I + TUREAI + ALGFIXI ! + TOTFLOODNI
        !FLOODNY  = TOTFLOODNI

        CUMBAL   = 0.0
        DAYBAL = 0.0
        CALL YR_DOY(INCDAT(YRDOY,-1), YR, DOY)
        
        WRITE(LUNSNC,10)
   10     FORMAT('!',15X,'------- STATE VARIABLES -------',
     &    '  ---- ADDED ----   --------------- REMOVED TODAY -----',
     &    '----------     DAILY     CUMUL',/,
     &    '@YEAR DOY  DAS      NO3     NH4   TUREA  ALGFIX',
    !  &    '@YEAR DOY  DAS      NO3     NH4   TUREA  FLOODN  ALGFIX',
     &    '   AFERT  AMINER    RLCH    RNOX    RNUP    RAML    RNRO',
     &    '  RIMMOB      DBAL      CBAL')
        WRITE (LUNSNC,50) YR, DOY, 0, 
     &    TNO3, TNH4, TUREA, ALGFIX, 
    ! &    TNO3, TNH4, TUREA, TOTFLOODN, ALGFIX,
     &    0.0, 0.0, 
     &    0.0, 0.0, 0.0, 
     &    0.0, 0.0, 0.0, 0.0, 0.0
       ! For Cell Detail N
        CCelNBal = 0.
        CALL GETLUN("CellDetailN.OUT", CLunn)
   !   OPEN (UNIT = CLunn, FILE = "CellDetailN.OUT", STATUS = 'OLD',
   ! &    POSITION = 'APPEND')
        OPEN (UNIT = CLunn, FILE ="CellDetailN.OUT", STATUS = 'REPLACE')
        WRITE(CLunn,'("*Nitrate BALANCE FOR CELL(",I2,",",I2,")")') 
     &      Cell_detail%row, Cell_detail%col
        CALL HEADER(SEASINIT, CLunn, CONTROL % RUN)
        WRITE (CLunn,1130)
        CelNtot =10.* (1/CelHT) * CelWD * CelHT * (
     &  CellDetail % state % SNO3 +
     &  CellDetail % state % SNH4 +
     &  CellDetail % state % UREA )
        WRITE (CLunn,1325) YR, DOY, 0, 
     &    TNO3, TNH4, TUREA, 0.0, 0.0, 0.0, 
    ! &    TNO3, TNH4, TUREA, TOTFLOODN, ALGFIX,
     &    0.0, 0.0,  0.0, 0.0, 0.0
      ENDIF
     
      
      ! JZW NHFlux, NVFlux represen Nitrogen; WHFlux, WVFlux, represent water
 1130 FORMAT('!',15X,'   STATE   --------------- ADDED ',
     & '---------------   --------------------- REMOVED TODAY ',
     & '--------------------     DAILY     CUMUL',/,
     & '@YEAR DOY   DAS   TotalN   AFERT     CelNom', ! NH4(i,j)    NO3(i,j)   
     & '    NFluxL    NFluxD    NFluxR    NFluxU  ', !State vars
     & ' NH4UpTak  NO3UpTak   NOX       CelIMM   CellNBal  CumNBal')  !Inflows

!      WRITE (CLUNN,1320) YR2, DY2, DAS, 0.0, 0.0, 
  !   &      SWijcm2,                       !State variables
   !  &      0.0, 0.0,                   !Inflows
   !  &      0.0, 0.0, 0.0, 0.0,         !Outflows
   !  &      0.0, 0.0,                   !Balance
  !   &      SWimj, SWijm, SWij, SWijp, SWipj, !Extras
  !   &      0.0, 0.0, 0.0,              !Extras
   !  &      0.0, 0.0, 0.0, 0.0
      CALL SoilNBalSum (CONTROL, TNH4=TNH4, TNO3=TNO3)

!***********************************************************************
!***********************************************************************
!     Daily output
!***********************************************************************
      ELSEIF (DYNAMIC == OUTPUT) THEN
!     ------------------------------------------------------------------
      IF (IDETL == 'D') THEN
        !Compute daily rates from cumulative values

!       Additions:
        AMTFERTODAY = AMTFER  - AMTFERY
        MINERTODAY  = CMINERN - CMINERY

!       Subtractions:
        LCHTODAY   = TLCH - TLCHY
        NOXTODAY   = TNOX - TNOXY
        WTNUPTODAY = (WTNUP - WTNUPY) * 10.
        AMLTODAY   = TOTAML - TOTAMLY
        FNROTODAY  = CUMFNRO - CUMFNROY
        IMMOBTODAY = CIMMOBN - CIMMOBY

        !FLOODNTODAY = TOTFLOODN - FLOODNY

        TOTSTATE = TNO3 + TNH4 + TUREA + ALGFIX ! + TOTFLOODN
        TOTADD   = AMTFERTODAY + MINERTODAY
        TOTSUB   = LCHTODAY  + NOXTODAY  + WTNUPTODAY + AMLTODAY + 
     &                FNROTODAY + IMMOBTODAY
        DAYBAL = TOTSTATE - TOTSTATY - TOTADD + TOTSUB
        CUMBAL   = CUMBAL + DAYBAL

        CALL YR_DOY(YRDOY, YR, DOY)
!       Write daily output to SoilNiBal.OUT.
        WRITE (LUNSNC,50) YR, DOY, DAS, 
     &    TNO3, TNH4, TUREA, ALGFIX, 
   ! &    TNO3, TNH4, TUREA, TOTFLOODN, ALGFIX,   
     &    AMTFERTODAY, MINERTODAY, 
     &    LCHTODAY, NOXTODAY, WTNUPTODAY, 
     &    AMLTODAY, FNROTODAY, IMMOBTODAY, DAYBAL, CUMBAL
   50     FORMAT(I5, I4.3, I5, 1X, 12F8.3, 2F10.3)
  !50     FORMAT(I5, I4.3, I5, 1X, 13F8.3, 2F10.3)
 !       WRITE (CLUN,1320) YR2, DY2, DAS, Time, TimeIncr, 
  !   &      SWijcm2,                                      !State vars
  !   &      Cell_detail%h_in, Cell_detail%v_in,           !Inflows
  !   &      Cell_detail%h_out, Cell_detail%v_out,         !Outflows
  !   &      RWUij, ESij,                                  !Outflows
  !   &      WBALij, CBALij,                               !Balance
  !   &      SWimj, SWijm, SWij, SWijp, SWipj, 
  !   &      0.0, Diffus(Row,Col-1), Diffus(Row,Col), 
  !   &      0.0, Kunsat(Row,Col), 
  !   &      Cell_detail%vdiff, Cell_detail%vgrav
 1320   FORMAT(1X,I4,1X,I3.3,1X,I5,2F7.3,14F10.6,3F10.2,4F10.6) 
      !Cells(L, J) % State % SNO3; Cells(L, J) % State % SNH4;
      !CELLS % RATE % SWFlux_L
      
      CelNtot =10.* (1/CelHT) * CelWD * CelHT * (
     & CellDetail % state % SNO3 +
     & CellDetail % state % SNH4 +
     & CellDetail % state % UREA )
      !10 * ug/cm = (1/cm) * cm * cm * kg/ha
      ! TNH4  = TNH4  + SNH4_2D(L, J) * ColFrac(j)
      ! ! SNH4_2D(L,J)   Total extractable ammonium N in soil layer L (kg [N] / ha)
      FertAdd = 10.* (1/CelHT) * CelWD * CelHT * (
     & FERTDATA % ADDSNO3(detailRow) + 
     &    FERTDATA % ADDSNH4(detailRow) + FERTDATA % ADDUREA(detailRow))
      If (detailRow .LE. NELEM) then 
         CelMNR = 10.* (1/CelHT) * CelWD * CelHT * 
     &     MNR(detailRow,detailCol)
         CelIMM = 10.* (1/CelHT) * CelWD * CelHT * 
     &     IMM(detailRow,detailCol) 
      else
        CelMNR = 0.  
      Endif
      CelNBal = CelNBal - YCelNBal +  10.* (1/CelHT) * CelWD * CelHT *(
     & Cell_Ndetail % NFlux_L +
     &  FertAdd+CelMNR+Cell_Ndetail % NFlux_D -Cell_Ndetail % NFlux_R - 
     &  Cell_Ndetail % NFlux_U -
     & CellDetail % Rate % NH4Uptake -
     & CellDetail % Rate % NO3Uptake -DENITRIF(detailRow, detailCol) 
     & - CelIMM)
      CCelNBal = CCelNBal + CelNBal
        WRITE (CLunn,1325) YR, DOY, DAS, CelNtot, FertAdd, CelMNR,
     &    Cell_Ndetail % NFlux_L,
     &    Cell_Ndetail% NFlux_D, 
     &    Cell_Ndetail % NFlux_R, 
     &    Cell_Ndetail%NFlux_U,
     &    CellDetail % Rate % NH4Uptake,
     &    CellDetail % Rate % NO3Uptake,
     &    DENITRIF(detailRow, detailCol), CelIMM,
     &    CelNBal, CCelNBal
 1325   FORMAT(1X,I4,1X,I3.3,1X,I5,F9.2,14F10.5,3F10.2,F10.6, 4F10.6)  
!       Save today's cumulative values for use tomorrow
        AMTFERY  = AMTFER
        CUMFNROY = CUMFNRO
        !FLOODNY  = TOTFLOODN
        CIMMOBY  = CIMMOBN
        CMINERY  = CMINERN
        TLCHY    = TLCH
        TNOXY    = TNOX
        TOTAMLY  = TOTAML
        WTNUPY   = WTNUP

        TOTSTATY = TOTSTATE
        YCelNBal = CelNBal
        CumNH4U = CumNH4U + CellDetail % Rate % NH4Uptake
        CumNO3U = CumNO3U + CellDetail % Rate % NO3Uptake
      ENDIF

!***********************************************************************
!***********************************************************************
!     Seasonal Output
!***********************************************************************
      ELSEIF (DYNAMIC == SEASEND) THEN
!     ------------------------------------------------------------------
!     N balance will be off by the amount of N uptake on the last day 
!       of season because this amount has been added to the N uptake by
!       plant, but has not been subtracted from soil.  This is because
!       soil processes are computed before plant processes.
!     May want to subtract this amount from balance?

        CALL YR_DOY(YRSIM, YRI, DOYI)
        CALL YR_DOY(YRDOY, YR, DOY)
!       Add the fertilizer N to the initial N pool. Also add the N from
!       organic residues applied during the growth period and sum this
!       with the initial TALLNI to make the balance fit with the final
!       TALLN. SEEDNI is not needed, because for the plant the NBAL
!       only deals with N uptake from the soil.
        TALLNI  = TALLNI + AMTFER + CMINERN !Initial + add

!       Sum state N at end of season
        STATEN = TNO3 + TNH4 + TUREA          !State end of day

!       Sum the initial value of all abiotic N pools (soil, air,
!       fertilizer), SOM and N uptake (WTNUP multiplied by 10 to
!       convert from g/m2 to kg/ha). Deduct the N in senesced material
!       from the N removed by the plant, because senesced material has
!       been returned to the soil.
        TALLN = STATEN +                          !State end of day
     &          TLCH + TNOX + WTNUP * 10. +       !Losses
     &          TOTAML + CIMMOBN           !Losses

!       Write output to NBAL.OUT.
        WRITE (LUNSNC,100) YRI, DOYI, YR, DOY

        WRITE (LUNSNC, 200) TNO3I, TNO3, TNH4I, TNH4, TUREAI, TUREA

        IF (NBUND .GT. 0) THEN
          WRITE(LUNSNC,300) ALGFIXI, ALGFIX
          ! WRITE(LUNSNC,300) TOTFLOODNI, TOTFLOODN, ALGFIXI, ALGFIX
          TALLN = TALLN +  ALGFIX !+TOTFLOODN +
          !JZW, JW, do we need TALLN?  TALLN = TALLN + TOTFLOODN + ALGFIX
        ENDIF

        WRITE (LUNSNC,600) AMTFER, CMINERN, TLCH, TNOX, 
     &        WTNUP * 10., TOTAML, CIMMOBN

        IF (NBUND .GT. 0) THEN
          TALLN = TALLN + CUMFNRO     !Factor in runoff over bund
          WRITE(LUNSNC,700) CUMFNRO
        ENDIF

        WRITE(LUNSNC,800) TALLNI, TALLN
        BALANCE = TALLN - TALLNI
        WRITE(LUNSNC,900) BALANCE

        CLOSE (UNIT = LUNSNC)
       
100     FORMAT (//,'!',T42,'INITIAL Year/Doy', T64, 'FINAL Year/Doy', 
     &   /,'!',T49, I5,'/',I3.3, T69,I5,'/',I3.3,
     &    /,'!','SOIL INORGANIC N BALANCE',T49,
     &        '-----------kg N/ha-----------')

200     FORMAT (
     &     /,'!', 3X, 'Soil NO3',             T48, F10.2, T68, F10.2,
     &     /,'!', 3X, 'Soil NH4',             T48, F10.2, T68, F10.2,
     &     /,'!', 3X, 'Soil Urea',            T48, F10.2, T68, F10.2)

300     FORMAT (
     &       '!', 3X, 'Flood water N',        T48, F10.2, T68, F10.2, ! JW, JZW how to change?
    !    &       '!', 3X, 'Flood water N',        T48, F10.2, T68, F10.2,
     &     /,'!', 3X, 'Algae N',              T48, F10.2, T68, F10.2)

600     FORMAT (
     &     /,'!', 3X, 'Added to / Removed from Soil:'
     &     /,'!', 3X, 'Fertilizer N',         T48, F10.2,
     &     /,'!', 3X, 'Mineralized N',        T48, F10.2,

     &     /,'!', 3X, 'Leached NO3',                      T68, F10.2,
     &     /,'!', 3X, 'N Denitrified',                    T68, F10.2,
     &     /,'!', 3X, 'N Uptake From Soil',               T68, F10.2,
     &     /,'!', 3X, 'Ammonia volatilization',           T68, F10.2,
     &     /,'!', 3X, 'N immobilized',                    T68, F10.2)

700     FORMAT ('!',3X,'N in runoff over bund',           T68, F10.2)

800     FORMAT (/,'!',3X, 'Total N balance',     T48, F10.2, T68, F10.2)
900     FORMAT ('!',3X,'Balance',                            T68, F10.3)
       !       Write output to CellDetailN.OUT.
        WRITE (CLunn,100) YRI, DOYI, YR, DOY
        WRITE (CLunn, 201) CelNtotI, CelNtot
        WRITE (CLunn,601) CumNH4U + CumNO3U
        !WRITE(CLunn,800) TALLNI, TALLN
        BALANCE = CelNtot - CelNtotI
        WRITE(CLunn,900) BALANCE
201     FORMAT (
     &     /,'!', 3X, 'Soil TotN',            T48, F10.2, T68, F10.2)
601     FORMAT (
     &     /,'!', 3X, 'Total N Uptake From Soil',         T68, F10.2)
        CLOSE (UNIT = CLunn)
        CALL SoilNBalSum (CONTROL, 
     &    AMTFER, Balance, 
     &    TLCH=TLCH, TNH4=TNH4, TNO3=TNO3, 
     &    TNOX=TNOX, TOTAML=TOTAML, WTNUP=WTNUP*10.)

!***********************************************************************
!***********************************************************************
!     End of DYNAMIC IF construct
!***********************************************************************
      END IF

      RETURN
      END SUBROUTINE SoilNiBal_2D

!=======================================================================
! SoilNiBal_2D VARIABLE DEFINITIONS:
!-----------------------------------------------------------------------
! ALGFIX        N in algae (kg [N] / ha)
! AMTFER     Cumulative amount of N in fertilizer applications
! BD(L,j)   Bulk density, soil layer L (g [soil] / cm3 [soil])
! cellNtot  Nitrogen in soil cell (µg[N] / g[soil])
! CMINERN   Cumulative seasonal mineralization of N in soil profile (kg[N]/ha)
! CUMFNRO   Cumulative N lost in runoff over bund (kg [N] / ha)
! KG2PPM(L) Conversion factor to switch from kg [N] / ha to ug [N] / g 
!           KG2PPM(L) = 1.0/(BD*1.E-01*DLAYR(L))
! NFlux_D   Downward movement of nitrogen with the water flow. (kg [N] / ha / d)
! NFlux_U   upward movement of nitrogen with the water flow.
! NH4Uptake      in kg[N]/ha
! NO3Uptake      in kg[N]/ha
! SNO3(L,j) Total extractable nitrate N in soil layer L (kg [N] / ha) 
! TNOX      Season cumulative denitrification across the total soil profile adding to 
!                 the nitrous oxide (NOx) pool of the air (kg [N] / ha)  
! TNOXY     Yesterday's TNOX
! TOTAML         Cumulative ammonia volatilization (kg [N] / ha)
!-----------------------------------------------------------------------
! END SUBROUTINE SoilNiBal_2D
!=======================================================================

