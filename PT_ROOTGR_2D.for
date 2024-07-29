C=======================================================================
C  PT_ROOTGR_2D, Subroutine
C
C  Determines root growth
C-----------------------------------------------------------------------
C  Revision history
C
C  Written
C  09/  /1988 EA & BB Modified by E. Alocilja & B. Baer 
C  04/  /1989 TJ Modified by T. Jou
C  02/08/1989 PWW Header revision and minor changes 
C  02/08/1993 PWW Added switch block, etc. 
C  08/23/2001 CHP Modified for modular format
C-----------------------------------------------------------------------
C  INPUT  : None
C
C  LOCAL  : RLDF,RNFAC,RLNEW,RLVF,SWDF,TRLDF,RNLF,L,L1
C
C  OUTPUT : None
C-----------------------------------------------------------------------
C  Called : WATBAL
C
C  Calls  : None
C-----------------------------------------------------------------------
C                         DEFINITIONS
C
C  RLDF(i, j) : A root length density factor for soil layer L used to calculate
C           new root growth distribution - unitless
C  RNFAC  : Zero to unity factor describing mineral N availability effect on
C           root growth in Layer L
C  RLNEW  : New root length to be added to the total root system length -
C           cm.  root per sq. cm. ground
C  RLVF   :
C  SWDF   : Soil water deficit factor for Layer L used to calculate root
C           growth and water uptake - unitless value between 0 and 1
C  TRLDF  : An intermediate calculation used to calculate distribution of
C           new root growth in soil
C  RNLF   : Intermediate factor used to calculate distribution of new root
C           growth in the soil - unitless value between 0 and 1
C  L,L1   : Loop counter
C=======================================================================

      SUBROUTINE PT_ROOTGR_2D (DYNAMIC, CELLS, YRDOY,
     &    DLAYR, DS, DTT, FILEIO, GRORT, ISWNIT,     !Input
     &    NH4, NLAYR, NO3, PLTPOP, SHF, SWFAC,    !Input
     &    CUMDEP, RLV, RTDEP)                             !Output

!-----------------------------------------------------------------------
      USE Cells_2D
      USE ModuleDefs     !Definitions of constructed variable types, 
                         ! which contain control information, soil
                         ! parameters, hourly weather data.
      IMPLICIT  NONE
      SAVE

      LOGICAL FIRST
      CHARACTER*1   ISWNIT
      CHARACTER*30 FILEIO

      INTEGER DYNAMIC, L, J, NLAYR, YRDOY0, YRDOY, DAS
      INTEGER ROW, Col, LastCol, LastRow
      INTEGER CSTIMDIF      ! Time difference function
      REAL HalfRow, ROWSPC

      REAL CUMDEP, DEP, DEPMAX, DTT, GRORT, PLTPOP
      REAL RLINIT, RLNEW, RLWR, RNFAC, RNLF, RTDEP, RTDEPI, RTWIDI
      REAL SDEPTH, SWDF, SWFAC, TRLDF, TRLV
      REAL CumWid, LastCumDep, LastCumWid
      REAL RTDEPnew, RTWID
      REAL RTWIDr(MaxRows), RTWIDnew(MaxRows), WidMax(MaxRows)
      REAL WidFrac(MaxRows,MaxCols), DepFrac(MaxRows,MaxCols) 
      REAL TotRootMass, RFAC3 
      ! above two variables are required in Aggregate_Roots. Here do not need

      REAL, DIMENSION(NL) :: DLAYR, DS
      REAL, DIMENSION(NL) :: NO3, NH4, RLV, SHF ! SW, RLDF
      !REAL, DIMENSION(MaxRows,MaxCols) :: NO3_2D, NH4_2D
      REAL, DIMENSION(MaxRows,MaxCols) :: RLV_2D, RLDF, DUL, LL, SWV
      REAL, DIMENSION(MaxRows,MaxCols) :: Thick, Width, CellArea, ESW 
      
      INTEGER, DIMENSION(MaxRows,MaxCols) :: TypeCell
      TYPE (CellType) CELLS(MaxRows,MaxCols)
      TYPE (CellStrucType) Struc(MaxRows,MaxCols)


!***********************************************************************
!***********************************************************************
!     Seasonal Initialization - Called once per season
!***********************************************************************
      IF (DYNAMIC .EQ. SEASINIT) THEN
!-----------------------------------------------------------------------
      ! We only need RTWIDI to calculate LastCol
      YRDOY0 = YRDOY
      STRUC = CELLS%STRUC
      Thick = STRUC%THICK
      Width = STRUC%WIDTH
      CellArea = STRUC%CellArea
      TypeCell = STRUC%CellType
      DUL = CELLS%STATE%DUL
      LL  = CELLS%STATE%LL
      !SAT = CELLS%STATE%SAT
      RLV_2D   = 0.0
      
      CALL PT_IPROOT_2D(FILEIO,                    !Input
     &               RLWR, SDEPTH, ROWSPC)              !Output

!********* TEMPORARY CHP *********************************
!     RLWR Sensitivity
!      SELECT CASE(RUN)
!        CASE(1); RLWR = 0.50
!        CASE(2); RLWR = 0.75
!        CASE(3); RLWR = 2.5
!        CASE(4); RLWR = 5.0
!        CASE(5); RLWR = 7.5
!        CASE(6); RLWR = 10.0
!      END SELECT
!*********************************************************

      FIRST = .TRUE.

      ! DO L = 1, NL
      !   RLV_2D(L) = 0.0
      ! END DO
      RLV_2D = 0.0
      DEPMAX = DS(NRowsTot) !DEPMAX = DS(NLAYR)
      CUMDEP = 0.0
      RTDEP  = 0.0 
      
      WidMax = 0.0
      
      HalfRow = ROWSPC * 100. / 2. 

      CALL PT_Aggregate_Roots(
     &    DLAYR, HalfRow,                               !Input
     &    NLAYR, RLV_2D, Struc,                         !Input
     &    RLV, TRLV)                                    !Output
      LastRow = 1
      LastCol = 1
      CALL PT_OPRoots_2D(TotRootMass, RFAC3, RLV_2D, Thick, Width)
       WRITE (92,1110)
 1110 FORMAT('Row, Col,',
     & ' RTDEP, CUMDEP, LastCumDep, RTDEPnew, Depfrac, LastRow, ',    
     & ' RTWIDr,CumWid,LastCumWid,RTWIDnew,Widfrac,LastCol,RLDF,',
     &   'RLV_2D in cm root / cm3 soil,RNFAC,NH4(Row), NO3(Row)')   

!***********************************************************************
!***********************************************************************
!     Daily rate calculations 
!***********************************************************************
      ! ELSEIF (DYNAMIC .EQ. RATE) THEN
      ELSEIF (DYNAMIC .EQ. INTEGR) THEN
!-----------------------------------------------------------------------
      SWV = CELLS%STATE%SWV  
!     Initial root distribution:  
      IF (FIRST) THEN

!********* TEMPORARY CHP *********************************
!     RLWR Sensitivity
!     Write to Overview.out file - can't do it when value is 
!       set because file is not open yet.
!      CALL GETLUN('OUTO',L)   !Get unit # for Overview.out
!      WRITE(L,*) ' Sensitivity analysis. RLWR = ',RLWR  
!*********************************************************

        !RTDEPI = SDEPTH  
        RTDEPI = MIN(20.0,DS(NLAYR))     !CHP per JWJ  
        ! Initial root width (specify half because we are modeling half a row)             
        RTWIDI = MIN(BedDimension%BEDWD / 2.0, RTDEPI / 4.0)     
        ! Tomato 2D use *.spe to give RTWIDI, YRTFACH, XRTFACH. Potato 2D does not need to change *.spe
        ! ROOOTS_2D use YRTFACH, XRTFACH to calculate RFAC2H which is used to calculate RTWIDnew 

        FIRST  = .FALSE.
        
C-------------------------------------------------------------------------
!       CHP 5/29/03 - Added this section based on CROPGRO initialization
!           at emergence. 
C       INITIALIZE ROOT DEPTH AT EMERGENCE
C       DISTRIBUTE ROOT LENGTH EVENLY IN ALL LAYERS TO A DEPTH OF
C       RTDEPTI (ROOT DEPTH AT EMERGENCE)
C-------------------------------------------------------------------------
        CUMDEP = 0.
        RLINIT = GRORT * RLWR * PLTPOP 
        ! JZW GROUT is previous days data
        !JZW in the begining days & lest days of plant life time, the RLNNEW might be zero 
        CALL PT_INROOT_2D(
     &  DepMax, HalfRow, RLINIT,                           !Input
     &  RTDEPI, RTWIDI, Thick, WidMax, Width,             !Input
     &  RLV_2D, RTDEP, RTWID, RTWIDr)                      !Output 
          
!        DO L = 1,NLAYR
!          DEP = MIN(RTDEPI - CUMDEP, DLAYR(L))
!   !       RLINIT = WTNEW * FRRT * PLTPOP * RFAC1 * DEP / ( RTDEP *
!   !    &     10000 )     
!          CUMDEP = CUMDEP + DEP
!          DO J = 1,NColsTot
!            IF (TypeCell(L, J) < 3 .OR. TypeCell(L,J) > 5) CYCLE
!            RLV_2D(L, J) = RLINIT / DLAYR(L)
!          Enddo
!            IF (CUMDEP .GE. RTDEPI) EXIT
!        ENDDO
!
!        RTDEP = RTDEPI 

!***********************************************************************
      ELSE !if not first, i.e not initial
!     Daily root growth and distribution

        RLNEW  = GRORT * RLWR * PLTPOP  !CHP   
        TRLDF  = 0.0
        CUMDEP = 0.0
        RNFAC  = 1.0
        RTDEPnew = RTDEP
        RTWIDnew = RTWIDr !it is array
        
!     First, root expansion.
!     Root depth is calculated in column 1 only.
!     Root width is calculated for each row. 
        RowLoop: DO Row = 1, NRowsTot
          LastCumdep = CUMDEP      
          CUMDEP = CUMDEP + Thick(Row,1)
          CumWid = 0.0
          ColLoop: Do Col = 1, NColsTot
            IF (TypeCell(Row,Col) < 3 .OR. TypeCell(Row,Col) > 5) CYCLE
            IF (TypeCell(Row,Col) .EQ. 3) 
     &          WIDMAX(Row) = BedDimension % BEDWD / 2
            IF (TypeCell(Row,Col) .EQ. 4 .OR. TypeCell(Row,Col) .EQ. 5) 
     &          WIDMAX(Row) = HalfRow
            LastCumWid = CumWid
            CumWid = CumWid + Width(Row,Col)

            ESW(ROW, Col) = DUL(ROW, Col) - LL(ROW, Col)
            SWDF   = 1.0
            IF (SWV(ROW, Col)-LL(ROW, Col) .LT. 0.25*ESW(ROW, Col)) THEN
              SWDF = 4.0*(SWV(ROW, Col)-LL(ROW, Col))/ESW(ROW, Col)
            ENDIF
            SWDF = AMAX1 (SWDF,0.0)      
            IF (ISWNIT .NE. 'N') THEN
!             RNFAC = 1.0 - (1.17 * EXP(-0.15 * TOTIN)
!             RNFAC = 1.0 - (1.17 * EXP(-0.15 * (SNH4(L) + SNO3(L))))
              RNFAC = 1.0 - (1.17 * EXP(-0.15 * (NH4(ROW) + NO3(ROW))))
              ! JZW NH4 and NO3 need to be 2D !!!
              RNFAC = AMAX1 (RNFAC,0.01)
            ENDIF
            ! Weighting factor for each cell
            RLDF(Row,Col) =AMIN1(SWDF,RNFAC)*SHF(Row)*CellArea(Row,Col)
           
!           Calculate new vertical growth in column 1 only
            IF (COL == 1) THEN
              WidFrac = 1.0
              IF (RTDEP >= CUMDEP) THEN
                DepFrac(Row,Col) = 1.0
              ELSEIF (RTDEP >= LastCumDep)THEN
!               Roots have partially filled the depth of this cell
                IF (CELLS(Row,Col)%STATE%WR > 0.0 .AND. RLNEW >0.0) THEN
                    !CHP: RTWID = RTWID + DTT*0.65*AMIN1((SWFAC*2.0),SWDF) ! Where 0.65 is an assumed data
                  RTDEPnew = RTDEP + DTT * 1.3 *
     &                       AMIN1((SWFAC * 2.0 ), SWDF)
                  RTDEPnew = MIN(RTDEPnew, DEPMAX)
                ENDIF
                  DepFrac(Row,Col) = MIN(1.0, 1. - (CUMDEP - RTDEPnew)/ 
     &                         Thick(Row,Col))
                  ! if the new root is more than one row, take minum
                  IF (Row > LastRow) LastRow = Row
                  ! JZW: we'd better to add exit statement
               ELSE
!                 No roots in this cell 
                 !JZW this is equivalent exit the do loop of row
                  DepFrac(Row,Col) = 0.0
               ENDIF
          
!             Check for new roots in this cell
              IF (RTDEPnew > LastCumDep .AND. 
     &          RTDEP <= LastCumDep) THEN
!               New roots have just grown into this cell
                RTWIDnew(Row) = Width(Row,Col)
                DepFrac(Row,Col) = MIN(1.0, 1. - (CUMDEP - RTDEPnew) / 
     &                         Thick(Row,Col))
              ENDIF
          
            Else ! if col!=1
  !           Calculate new horizontal growth in this cell (RTWIDnew) 
!             horizontal portion of cell occupied by roots (WidFrac)
!           Horizontal root growth only occurs when DepFrac of adjacent 
!             cell is > 0.99.  No need to calculate for Column 1, since
!             width fraction is initialized to 1.0 there.
              IF (RTWIDr(Row) >= CumWid) THEN
                WidFrac(Row,Col) = 1.0
                DepFrac(Row,Col) = 1.0
              ELSEIF (RTWIDr(Row) >= LastCumWid) THEN
!             Roots have partially filled the width of this cell
                IF (CELLS(Row,Col)%STATE%WR > 0.0 .AND. RLNEW >0.0) THEN
                    !CHP: RTWID = RTWID + DTT*0.65*AMIN1((SWFAC*2.0),SWDF) ! Where 0.65 is an assumed data
                    RTWIDnew(Row) = RTWIDr(Row) + DTT * 1.3 *
     &                          AMIN1((SWFAC*2.0),SWDF) 
                  RTWIDnew(Row) = MIN(RTWIDnew(Row), WIDMAX(Row))
                ENDIF
                WidFrac(Row,Col) = MIN(1.0, 1. - (CumWid -RTWIDnew(Row))
     &                        / Width(Row,Col))    
                
                IF (Col > LastCol) LastCol = Col
                  DepFrac(Row,Col) = 1.0
                ELSE
!               No roots in this cell
                  WidFrac(Row,Col) = 0.0
                  DepFrac(Row,Col) = 0.0
                ENDIF
             
    !           Check for new roots in this cell
                IF (RTWIDnew(Row) > LastCumWid .AND. 
     &              RTWIDr(Row) <= LastCumWid) THEN
!             New roots have just grown into this cell
                  WidFrac(Row,Col) = MIN(1.0, 1. -(CumWid-RTWIDnew(Row))
     &                        / Width(Row,Col))
      
                  DepFrac(Row,Col) = 1.0 
                  IF (Col > LastCol) LastCol = Col 
                ENDIF
            ENDIF !! end of  col!=1  
            
!-----------------------------------------------------------------------
!         Apply factor for this cell
            RLDF(Row,Col) =
     &               RLDF(Row,Col)*DepFrac(Row,Col)*WidFrac(Row,Col)
!         Sum of all factors
            TRLDF = TRLDF + RLDF(Row,Col)
            WRITE (92,1120)Row, Col,
     &        RTDEP, CUMDEP, LastCumDep, 
     &        RTDEPnew, Depfrac(Row,Col), LastRow,  
     &        RTWIDr(Row),CumWid, LastCumWid,RTWIDnew(Row),
     &        Widfrac(Row,Col),LastCol,RLDF(Row,Col), RLV_2D(Row,Col),
     &        RNFAC, NH4(Row), NO3(ROW)
 1120 FORMAT(2(I4,","),5(F6.2,","),I2,",",
     &    5(F6.2,","),I2,",",F6.2,4(",",F6.2))   
       !JZW day 1980159 RLDF is wrong for row=5,col=1,due to RNFAC is wrong
            IF (RTWIDnew(Row) < CumWid) EXIT ColLoop    
          ENDDO ColLoop 
        ENDDO RowLoop
        RTDEP  = RTDEPnew 
        RTWIDr = RTWIDnew ! it is array

!-------------------------------------------------------------------------
        IF (TRLDF .GE. RLNEW*0.00001) THEN
           RNLF = RLNEW/TRLDF
           !DO L = 1, L1
           DO Row = 1, LastRow !JZW LastRow is the last row of root
             DO Col = 1, LastCol 
               IF (TypeCell(Row,Col)<3 .OR. TypeCell(Row,Col) > 5) CYCLE
             ! To calculate LastCol need RTWIDr(Row), LastCumWid, RTWIDI
               RLV_2D(Row,Col) = RLV_2D(Row,Col)
     &             +RLDF(Row,Col)*RNLF/DLAYR(ROW)-0.005*RLV_2D(Row,Col) 
             ! JZW minus is due to root senescence
               RLV_2D(Row,Col) = AMAX1 (RLV_2D(Row,Col),0.0)
               RLV_2D(Row,Col) = AMIN1 (RLV_2D(Row,Col),5.0)
             END DO
           ENDDO
        END IF
      ENDIF ! end of IF not (FIRST)
      TRLV = 0.0
      DO Row = 1, NRowsTot
        Do Col = 1, NColsTot
          IF (TypeCell(Row,Col) < 3 .OR. TypeCell(Row,Col) > 5) CYCLE
            TRLV = TRLV + RLV_2D(Row,Col) * CellArea(Row,Col) 
            ! JZW, TRLV is calculated in PT_Aggregate_Roots, we do not need to calculate here
        End do
      ENDDO
      IF (RTWIDr(Row) > RTWID) RTWID = RTWIDr(Row) 
      ! RTWID is not used, it can be as output of this subroutine for watch variable
      Write(92,1125) TRLV
 1125 Format("TRLV=", F8.1)
      CALL PT_Aggregate_Roots(
     &    DLAYR, HalfRow,                              !Input
     &    NLAYR, RLV_2D, Struc,                        !Input
     &    RLV, TRLV)                                   !Output

      CELLS%STATE%RLV = RLV_2D
      
       DAS = MAX(0,CSTIMDIF(YRDOY0,YRDOY))
       write(92,1130) YRDOY, DAS, RTDEP, TRLV
 1130  Format("YRDAY=", I8, ",DAS=", I4,", RTDEP=",F6.2, ",TRLV=",F8.1) 
       DO Row = 1, LastRow 
          write(92,1140)Row, RTWIDr(ROW), RLV(ROW)
       Enddo
 1140  Format("Row=", I2, ", RIWIDr=", F6.2, ",  RLV=", F6.2, 
     &  "cm3[root]/cm3[ground)")

       
!***********************************************************************
      ELSEIF (DYNAMIC == OUTPUT .OR. DYNAMIC == SEASEND) THEN
!-----------------------------------------------------------------------
      CALL PT_OPRoots_2D(TotRootMass, RFAC3, RLV_2D, Thick, Width)
     
!***********************************************************************
!     END OF DYNAMIC IF CONSTRUCT
!***********************************************************************
      ENDIF
!***********************************************************************
      RETURN
      END SUBROUTINE PT_ROOTGR_2D
C=======================================================================


C=======================================================================
C  PT_IPROOT_2D, Subroutine
C
C  Input data for potato root module
C-----------------------------------------------------------------------
C  Revision history
C
C  08/23/2001 CHP Written
C  10/25/2002 CHP Modified read format for Y2K
C  08/12/2003 CHP Added I/O error checking
C-----------------------------------------------------------------------

      SUBROUTINE PT_IPROOT_2D(FILEIO,                    !Input
     &                     RLWR, SDEPTH, ROWSPC)              !Output

!     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER LUNIO, LUNCRP
      CHARACTER*1, PARAMETER :: BLANK = ' '
      CHARACTER*6, PARAMETER :: ERRKEY = 'ROOTGR'

      CHARACTER*6   SECTION
      CHARACTER*12  FILEC
      CHARACTER*30  FILEIO
      CHARACTER*80  PATHCR
      CHARACTER*92  FILECC
      CHARACTER*180 CHAR

      INTEGER ERR, FOUND, ISECT, LINC, LNUM, PATHL

      REAL RLWR, SDEPTH, ROWSPC
!     LOGICAL EOF
!-----------------------------------------------------------------------
!     Read data from FILEIO for use in ROOTGR module
      CALL GETLUN('FILEIO', LUNIO)
      OPEN (LUNIO, FILE = FILEIO, STATUS = 'OLD', IOSTAT=ERR)
      IF (ERR .NE. 0) CALL ERROR(ERRKEY,ERR,FILEIO,0)

      READ(LUNIO,'(6(/),15X,A12,1X,A80)', IOSTAT=ERR) FILEC, PATHCR
      LNUM = 7
      IF (ERR .NE. 0) CALL ERROR(ERRKEY,ERR,FILEIO,LNUM)

C-----------------------------------------------------------------------
C    Read Planting Details Section
C-----------------------------------------------------------------------
      SECTION = '*PLANT'
      CALL FIND(LUNIO, SECTION, LINC, FOUND) ; LNUM = LNUM + LINC
      IF (FOUND .EQ. 0) THEN
        CALL ERROR(SECTION, 42, FILEIO, LNUM)
      ELSE
        !READ (LUNIO,'(55X,F5.1)', IOSTAT=ERR) SDEPTH ; LNUM = LNUM + 1
         READ (LUNIO,'(43X, F5.1, 6X, F5.1)', IOSTAT=ERR)  
     &    ROWSPC, SDEPTH ; LNUM = LNUM + 1  
        IF (ERR .NE. 0) CALL ERROR(ERRKEY,ERR,FILEIO,LNUM)
      ENDIF

      CLOSE (LUNIO)

C-----------------------------------------------------------------------
C     Read Crop Parameters from FILEC
C-----------------------------------------------------------------------
      LNUM   = 0
      PATHL  = INDEX (PATHCR,BLANK)
      IF (PATHL .LE. 1) THEN
         FILECC = FILEC
       ELSE
         FILECC = PATHCR(1:(PATHL-1)) // FILEC
      ENDIF
      CALL GETLUN('FILEC', LUNCRP)
      OPEN (LUNCRP,FILE = FILECC, STATUS = 'OLD',IOSTAT=ERR)
      IF (ERR .NE. 0) CALL ERROR(ERRKEY,ERR,FILEC,0)

!     EOF not portable. CHP 7/24/2007
!     DO WHILE (.NOT. EOF (LUNCRP))
      DO WHILE (ERR == 0)
        CALL IGNORE(LUNCRP,LNUM,ISECT,CHAR)
!       IF (ISECT .EQ. 0) CALL ERROR(ERRKEY,33,FILECC,LNUM)
        IF (ISECT .EQ. 0) EXIT
        IF (ISECT .EQ. 2) CYCLE
        IF (CHAR(10:13) .EQ. 'RLWR') THEN
          READ (CHAR,'(14X,F6.0)',IOSTAT=ERR) RLWR
          IF (ERR .NE. 0) CALL ERROR(ERRKEY,ERR,FILEC,LNUM)
          EXIT
        ENDIF
      ENDDO

      CLOSE (LUNCRP)

C-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE PT_IPROOT_2D
!=======================================================================
!  PT_INROOT Subroutine
!  Initializes root variables at emergence.
!----------------------------------------------------------------------
!  REVISION HISTORY
!  04/01/1991 GH  Adapted for CROPGRO
!  06/17/1998 CHP Modified for modular format
!  05/11/1999 GH  Incorporated in CROPGRO
!  02/21/2009 CHP Adapted for 2D roots
!-----------------------------------------------------------------------
!  Called : CROPGRO
!  Calls  : None
!=======================================================================
      SUBROUTINE PT_INROOT_2D(
     &  DepMax, HalfRow, RLINIT,     !Input
     &  RTDEPI, RTWIDI, Thick, WidMax, Width,              !Input
     &  RLV_2D, RTDEP, RTWID, RTWIDr)                      !Output

!     ------------------------------------------------------------------
      USE Cells_2D
      IMPLICIT NONE

      INTEGER Row, Col
      REAL DepMax, RLINIT, WidMax(MaxRows)
      REAL HalfRow, X, Z
      REAL RTDEPI, RTDEP, LastCumDep, CumDep
      REAL RTWIDI, RTWID, LastCumWid, CumWid, RTWIDr(MaxRows)
      REAL TotRootArea
      REAL, DIMENSION(MaxRows,MaxCols) :: Thick, Width, CellArea
      REAL, DIMENSION(MaxRows,MaxCols) :: RLV_2D, RootArea
!-----------------------------------------------------------------------
      RTDEPI = MAX(MIN(RTDEPI, DepMax), Thick(1,1))
      RTWIDI = MAX(MIN(RTWIDI, WidMax(1)), Width(1,1))
      RLV_2D = 0.0
      RootArea = 0.  !cell area containing roots
      TotRootArea = 0.0
      RTDEP = RTDEPI
      RTWID = RTWIDI
      RTWIDr = 0.0
      X = 0.0
      Z = 0.0

!     Distribute root length and width evenly thru cells
      CUMDEP = 0.
      RowLoop: DO Row = 1, NRowsTot
        LastCumDep = CUMDEP
        CUMDEP = CUMDEP + Thick(Row,1)
        IF (RTDEPI >= CUMDEP) THEN
          Z = Thick(Row,1)
        ELSEIF (RTDEPI > LastCumDep) THEN
          Z = RTDEPI - LastCumDep 
        ELSE
          Z = 0.0
          EXIT RowLoop
        ENDIF
        
        IF (Row == 1 .OR. Z > 0.98 * Thick(Row,1)) THEN
          RTWIDr(Row) = RTWIDI
        ELSEIF (Z > 0.0) THEN
          RTWIDr(Row) = WIDTH(Row,1)
        ENDIF
        
        CumWid = 0.
        ColLoop: DO Col = 1,NColsTot
          LastCumWid = CumWid
          CumWid = CumWid + Width(Row,Col)
          CellArea(Row,Col) = Width(Row,Col) * Thick(Row,Col)
          IF (RTWIDI >= CumWid) THEN
            X = Width(Row,Col)
          ELSEIF (RTWIDI > LastCumWid) THEN
            X = RTWIDI - LastCumWid
          ELSE
            X = 0.0
            EXIT ColLoop
          ENDIF
          
          IF (ROW == 1 .OR. COL == 1 .OR. Z > 0.98 * Thick(Row,Col))THEN
            RootArea(Row,Col) = X * Z
          ENDIF
          TotRootArea = TotRootArea + RootArea(Row,Col)
        ENDDO ColLoop
      ENDDO RowLoop
          
      DO Row = 1, NRowsTot
        DO Col = 1, NColsTot
          IF (RootArea(Row,Col) > 1.E-6) THEN
            RLV_2D(Row,Col) = RLINIT  * RootArea(Row,Col) / TotRootArea
            RLV_2D(Row,Col) = RLV_2D(Row,Col) / CellArea(Row,Col)
!            cm[root]         cm[root]      1  
!           ----------- = -------------- * ----
!            cm3[soil]    cm[row length]   cm2 
          ENDIF
        ENDDO
      ENDDO
      
!***********************************************************************
      RETURN
      END SUBROUTINE PT_INROOT_2D
!=======================================================================

!=======================================================================
!  OPRoots_2D, Subroutine, C.H.Porter from Soil Water portions of OPDAY
!  Generates output for daily soil water data
!-----------------------------------------------------------------------
!  REVISION       HISTORY
!  07/02/2009 CHP Written
!-----------------------------------------------------------------------
!  Called from:   WatBal2D
!  Calls:         None
!=======================================================================
      SUBROUTINE PT_OPRoots_2D(TotRootMass, RFAC3, RLV_2D, Thick, Width)
!                                 kg/ha,  cm/g, cm/cm3,   cm , cm
!-----------------------------------------------------------------------
      USE Cells_2D
      USE ModuleData
      IMPLICIT NONE
      SAVE

      REAL, DIMENSION(MaxRows,MaxCols), INTENT(IN) :: RLV_2D,Thick,Width
      REAL, INTENT(IN) :: TotRootMass, RFAC3

      CHARACTER*1 IDETG, IDETL, RNMODE
      CHARACTER*13 OUTRoot
      PARAMETER (OUTRoot = 'PT_Root2D.OUT')
      CHARACTER*17 FMT

      INTEGER COL, DAS, DOY, DYNAMIC, ERRNUM, FROP
      INTEGER NOUTDW, ROW, RUN
      INTEGER YEAR, YRDOY, REPNO, YRSTART, INCDAT

      LOGICAL FEXIST, DOPRINT

!-----------------------------------------------------------------------
!     Define constructed variable types based on definitions in
!     ModuleDefs.for.
      TYPE (ControlType) CONTROL
      TYPE (SwitchType)  ISWITCH
      
      CALL GET(CONTROL)

      DAS     = CONTROL % DAS
      DYNAMIC = CONTROL % DYNAMIC
      FROP    = CONTROL % FROP
      RUN     = CONTROL % RUN
      RNMODE  = CONTROL % RNMODE
      REPNO   = CONTROL % REPNO
      YRDOY   = CONTROL % YRDOY

      CALL YR_DOY(YRDOY, YEAR, DOY) 

!***********************************************************************
!***********************************************************************
!     Seasonal initialization - run once per season
!***********************************************************************
      IF (DYNAMIC == SEASINIT) THEN
!-----------------------------------------------------------------------
!   Set initial values to calculate average values
!-----------------------------------------------------------------------
      CALL GET(ISWITCH)
      IDETL   = ISWITCH % IDETL
      IDETG   = ISWITCH % IDETG

      IF (IDETG == 'N' .OR. IDETL == '0') THEN
        DOPRINT = .FALSE.
      ELSE
        DOPRINT = .TRUE.
      ENDIF
      IF (.NOT. DOPRINT) RETURN

!-----------------------------------------------------------------------
!   Generate headings for output file
!-----------------------------------------------------------------------
      CALL GETLUN('OUTRoot', NOUTDW)
      INQUIRE (FILE = OUTRoot, EXIST = FEXIST)
      IF (FEXIST) THEN
        OPEN (UNIT = NOUTDW, FILE = OUTRoot, STATUS = 'OLD',
     &    IOSTAT = ERRNUM, POSITION = 'APPEND')
      ELSE
        OPEN (UNIT = NOUTDW, FILE = OUTRoot, STATUS = 'NEW',
     &    IOSTAT = ERRNUM)
        WRITE(NOUTDW,'("*2D ROOTS DAILY OUTPUT FILE")')
      ENDIF

!-----------------------------------------------------------------------
!     Variable heading for WATER.OUT
!-----------------------------------------------------------------------
      IF (RNMODE .NE. 'Q' .OR. RUN .EQ. 1) THEN
        IF (RNMODE .EQ. 'Q') THEN
          CALL HEADER(SEASINIT, NOUTDW, REPNO)
        ELSE
          CALL HEADER(SEASINIT, NOUTDW, RUN)
        ENDIF

        YRSTART = YRDOY
        CALL YR_DOY(INCDAT(YRSTART,-1),YEAR,DOY)
      ENDIF

!***********************************************************************
!***********************************************************************
      ENDIF !DYNAMIC CONTROL
!***********************************************************************
!***********************************************************************
!     Daily Output
!***********************************************************************
      IF (DYNAMIC == SEASINIT .OR. DYNAMIC == OUTPUT .OR. 
     &      DYNAMIC == SEASEND) THEN
!-----------------------------------------------------------------------
      IF (DOPRINT) THEN
!           Print initial conditions, 
        IF (DYNAMIC == SEASINIT .OR.
!           Print every FROP days, and
     &     (DYNAMIC .EQ. OUTPUT .AND. MOD(DAS, FROP) .EQ. 0) .OR. 
!           Print on last day if not already done.
     &     (DYNAMIC .EQ. SEASEND  .AND. MOD(DAS, FROP) .NE. 0)) THEN

          Write(NOUTDW,'(/,"Year DOY:",I5,I4.3)') YEAR, DOY
          Write(NOUTDW,'("Root Mass =     ",F10.2," kg/ha")')TotRootMass
          Write(NOUTDW,'("Root L:M ratio =",F10.2," cm/g")') RFAC3

          Write(NOUTDW,'("  Column ->",20I10)') (Col, Col=1, NColsTOT)
          Write(NOUTDW,'("Width(cm)->",20F10.3)') 
     &                  (width(1,Col),Col = 1, NColsTOT)
          Write(NOUTDW,'("      Thick")') 
          Write(NOUTDW,'("Lyr    (cm)   ------- ",
     &  "RLV (cm[root]/cm3[soil] -------")')
          WRITE(FMT,'("(I3,F8.1,",I2,"F10.4)")') NColsTot 
          DO Row = 1, NRowsTot  
            Write(NOUTDW,FMT)     
     &      Row, Thick(Row,1), (RLV_2D(Row,Col),Col = 1, NColsTOT) 
          Enddo 

        ENDIF
      ENDIF

!***********************************************************************
!***********************************************************************
!     SEASEND - Sesaonal Output
!***********************************************************************
        IF (DYNAMIC .EQ. SEASEND) THEN
!-----------------------------------------------------------------------
            !Close daily output files.
            CLOSE (NOUTDW)
        ENDIF
!***********************************************************************
!***********************************************************************
!     END OF DYNAMIC IF CONSTRUCT
!***********************************************************************
      ENDIF
!***********************************************************************
      RETURN
      END SUBROUTINE PT_OPRoots_2D
!=======================================================================
      SUBROUTINE PT_Aggregate_Roots(
     &    DLAYR, HalfRow,                                 !Input
     &    NLAYR, RLV_2D, Struc,                           !Input
     &    RLV, TRLV)                                      !Output

      Use Cells_2D
      IMPLICIT NONE
      SAVE

      INTEGER Row, Col, L, NLAYR
      REAL HalfRow, RFAC3, TRLV
      REAL, DIMENSION(NL) :: DLAYR, RLV
      REAL, DIMENSION(MaxRows,MaxCols) :: RLV_2D, Width, Thick, RtLen
      TYPE (CellStrucType) Struc(MaxRows,MaxCols)

      Width = Struc%Width
      Thick = Struc%Thick

      TRLV = 0.0
      DO Row = 1, NRowsTot
        DO Col = 1, NColsTot
          RtLen(Row,Col) =RLV_2D(Row,Col)*THICK(Row,Col)*Width(Row,Col)
!             cm[root]         cm[root]
!          -------------- =   ----------- * cm[cell depth] * cm[cell width]
!          cm[row length]     cm3[ground]

          TRLV = TRLV + RtLen(Row,Col)
        ENDDO
      ENDDO

!     Aggregate cells across a row to get layer total.  Units for layers
!     are in cm[root]/cm[row length]
      CALL Cell2Layer_2D(
     &   RtLen, Struc, NLAYR,                 !Input
     &   RLV)                                  !Output

      DO L = 1, NLAYR
        RLV(L) = RLV(L) / HalfRow / DLAYR(L)
!    cm[root]       cm[root]            1               1
!  ----------- = -------------- * ------------- * -----------------
!  cm3[ground]   cm[row length]   cm[row width]   cm[row thickness]
      ENDDO

      RETURN
      END Subroutine PT_Aggregate_Roots
!==============================================================================
!-----------------------------------------------------------------------
! Variable definitions
!-----------------------------------------------------------------------
! CUMDEP       The buttom of current row
! CumWid       The width of the right side of the current column
! DepFrac      Fracton of the root in a row thickness
! DTT          Growing degree days today, degrees C 
! ESW(Row,Col) Plant extractable soil water by layer (= DUL - LL) (cm3/cm3)
! GRORT        Root growth rate, g/plant/day
! LastCumdep   The top of current layer
! LastCumWif   The width of the left side of the current column
! LastRow      the deepest of the row which is occupied by the root
! NH4(L)       Ammonium N in soil layer L (µg[N] / g[soil])
! PLTPOP       Plant population (# plants / m2)
! RLDF(Row,Col)Combined weighting factor to determine root growth distribution
! RLNEW        New root growth added (cm[root]/cm2[ground]/d)
! RLV_2D(Row, Col) Root length density for soil cell in cm root / cm3 soil
! RLWR         Root length to weight ration, cm/g   
! RNFAC        Zero to unity factor describing mineral N availability effect on
!              root growth in Layer L
! RNLF         Intermediate factor used to calculate distribution of new root(1/cm2[ground]/d)
! RTDEP        Root length in col=1 at the begining of the day
! RTDEPnew     Root length in col=1 at the end of the day
! RTWID        Maximum width used for watch variable
! RTWIDr(Row)  Root width for each row
! RTWIDnew(row)Root width for specific row at the end of the day (calculated up to last col)
! SHF          Soil hospitality factor 0-1  
! SWDF         Soil water deficit factor for layer with deepest roots (0-1)    
! SWFAC        Effect of soil-water stress on photosynthesis, 1.0=no stress,0.0=max stress 
! TRLDF        Total root length density factor for root depth (cm)
!***********************************************************************
! END SUBROUTINES PT_ROOTGR_2D, PT_IPROOT_2D
!=======================================================================
! JZW report to Cheryl, We need variable root lenth density cm/cm3
