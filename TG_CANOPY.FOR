C=======================================================================
C  TG_CANOPY, Subroutine, G. Hoogenboom, K.J. Boote, J.W. Jones
C  Calculates canopy height and canopy width as a function of V-Stage,
C  air temperature, drought stress, daylength, and radiation.
C-----------------------------------------------------------------------
C  REVISION HISTORY
C  05/01/1989 Written.
C  04/24/1994 NBP Changed TAIRHR to TGRO.  Edited indentation.
C  01/19/1996 KJB Include PAR effect on expansion.
C  07/15/1998 CHP Modified for modular format
C  05/15/1999 GH  Incorporated into CROPGRO
C  01/22/2003 KJB Add checks for minimum canopy height and width.
C  08/12/2003 CHP Revised I/O error checking
C  06/30/2004 CHP/CDM Added KC_SLOPE to SPE file and KC_ECO to ECO file.
C                 Added optional KCAN to ECO file.
!  03/14/2012 JZW Add the input argument PLME for read I File when PLME = "T" 
C-----------------------------------------------------------------------
C  Called : VEGGR
C  Calls  : ERROR, FIND, IGNORE
C========================================================================

      SUBROUTINE TG_CANOPY(DYNAMIC, CONTROL,PLME,
     &    ECONO, FILECC, FILEGC, KCAN, PAR, ROWSPC,       !Input
     &     RVSTGE, TGRO, TURFAC, VSTAGE, XLAI,             !Input
     &    CANHT, CANWH)                                   !Output
       !If read IFile VStage will be output JZW debug
C-----------------------------------------------------------------------
      USE ModuleDefs     !Definitions of constructed variable types, 
                         ! which contain control information, soil
                         ! parameters, hourly weather data.
      USE Interface_Read_IFile
      IMPLICIT NONE
      SAVE

      CHARACTER*1 PLME
      CHARACTER*6 ERRKEY
      PARAMETER (ERRKEY = 'TG_CANOPY')

      CHARACTER*6   SECTION
      CHARACTER*6   ECOTYP, ECONO
      CHARACTER*92  FILECC, FILEGC
      CHARACTER*255 C255

      INTEGER I, II, LUNCRP, LUNECO, ERR, LINC, LNUM, ISECT
      INTEGER DYNAMIC
      INTEGER FOUND

      REAL PAR, ROWSPC, RVSTGE, TURFAC, VSTAGE
      REAL CANHT, CANWH, XLAI
      REAL KCAN, RHGHT, RWIDTH
      REAL HWTEM, RCANHT, RCANWH, PARNOD, HPAR, WPAR
      REAL TABEX
      REAL XHWPAR(10), XHWTEM(10), YHWPAR(10), YHWTEM(10)
      REAL XVSHT(15), YVSHT(15), YVSWH(15)
      REAL TGRO(TS)
      TYPE (ControlType) CONTROL

!***********************************************************************
!***********************************************************************
!     Run Initialization - Called once per simulation
!***********************************************************************
      IF (DYNAMIC .EQ. RUNINIT) THEN
!-----------------------------------------------------------------------
!     Read in values from input file, which were previously input
!       in Subroutine IPCROP.
!-----------------------------------------------------------------------
      CALL GETLUN('FILEC', LUNCRP)      
      OPEN (LUNCRP,FILE = FILECC, STATUS = 'OLD',IOSTAT=ERR)
      IF (ERR .NE. 0) CALL ERROR(ERRKEY,ERR,FILECC,0)
      LNUM = 0
!-----------------------------------------------------------------------
!    Find and Read Photosynthesis Section
!-----------------------------------------------------------------------
!     Subroutine FIND finds appropriate SECTION in a file by
!     searching for the specified 6-character string at beginning
!     of each line.
!-----------------------------------------------------------------------
!CHP 7/30/2004 - Get KCAN from main routine.
!                May be overriden by value in ECOTYPE file.
!      SECTION = '!*PHOT'
!      CALL FIND(LUNCRP, SECTION, LINC, FOUND); LNUM = LNUM + LINC
!      IF (FOUND .EQ. 0) THEN
!        CALL ERROR(SECTION, 42, FILECC, LNUM)
!      ELSE
!        ISECT = 2
!        DO WHILE (ISECT .NE. 1)
!          CALL IGNORE(LUNCRP,LNUM,ISECT,C255)
!          IF (ISECT .EQ. 0) CALL ERROR(ERRKEY,ERR,FILECC,LNUM)
!        ENDDO
!        READ(C255,'(12X,F6.0)',IOSTAT=ERR) KCAN
!        IF (ERR .NE. 0) CALL ERROR(ERRKEY,ERR,FILECC,LNUM)
!      ENDIF
!-----------------------------------------------------------------------
C     ***** READ CANOPY HEIGHT & WIDTH PARAMETERS ******************
C-----------------------------------------------------------------------
      SECTION = '!*CANO'
      CALL FIND(LUNCRP, SECTION, LINC, FOUND); LNUM = LNUM + LINC
      IF (FOUND .EQ. 0) THEN
        CALL ERROR(SECTION, 42, FILECC, LNUM)
      ELSE
        CALL IGNORE(LUNCRP,LNUM,ISECT,C255)
        READ(C255,'(10F6.0)',IOSTAT=ERR)(XVSHT(II),II = 1,10)
        IF (ERR .NE. 0) CALL ERROR(ERRKEY,ERR,FILECC,LNUM)

        CALL IGNORE(LUNCRP,LNUM,ISECT,C255)
        READ(C255,'(10F6.0)',IOSTAT=ERR)(YVSHT(II),II = 1,10)
        IF (ERR .NE. 0) CALL ERROR(ERRKEY,ERR,FILECC,LNUM)

        CALL IGNORE(LUNCRP,LNUM,ISECT,C255)
        READ(C255,'(10F6.0)',IOSTAT=ERR)(YVSWH(II),II = 1,10)
        IF (ERR .NE. 0) CALL ERROR(ERRKEY,ERR,FILECC,LNUM)

        CALL IGNORE(LUNCRP,LNUM,ISECT,C255)
        READ(C255,'(5F6.0)',IOSTAT=ERR)(XHWTEM(II),II = 1,5)
        IF (ERR .NE. 0) CALL ERROR(ERRKEY,ERR,FILECC,LNUM)

        CALL IGNORE(LUNCRP,LNUM,ISECT,C255)
        READ(C255,'(5F6.0)',IOSTAT=ERR)(YHWTEM(II),II = 1,5)
        IF (ERR .NE. 0) CALL ERROR(ERRKEY,ERR,FILECC,LNUM)

        CALL IGNORE(LUNCRP,LNUM,ISECT,C255)
        READ(C255,'(8F6.0)',IOSTAT=ERR)(XHWPAR(II),II = 1,8)
        IF (ERR .NE. 0) CALL ERROR(ERRKEY,ERR,FILECC,LNUM)

        CALL IGNORE(LUNCRP,LNUM,ISECT,C255)
        READ(C255,'(8F6.0)',IOSTAT=ERR)(YHWPAR(II),II = 1,8)
        IF (ERR .NE. 0) CALL ERROR(ERRKEY,ERR,FILECC,LNUM)
      ENDIF

      CLOSE (LUNCRP)

C-----------------------------------------------------------------------
C    Read Ecotype Parameter File
C-----------------------------------------------------------------------
      CALL GETLUN('FILEE', LUNECO)
      OPEN (LUNECO,FILE = FILEGC,STATUS = 'OLD',IOSTAT=ERR)
      IF (ERR .NE. 0) CALL ERROR(ERRKEY,ERR,FILEGC,0)
      ECOTYP = '      '
      LNUM = 0
      DO WHILE (ECOTYP .NE. ECONO)
        CALL IGNORE(LUNECO, LNUM, ISECT, C255)
          IF ((ISECT .EQ. 1) .AND. (C255(1:1) .NE. ' ') .AND.
     &        (C255(1:1) .NE. '*')) THEN
          READ (C255,'(A6,90X,2(1X,F5.0))',IOSTAT=ERR)
     &        ECOTYP, RWIDTH, RHGHT
          IF (ERR .NE. 0) CALL ERROR(ERRKEY,ERR,FILEGC,LNUM)
          IF (ECOTYP .EQ. ECONO) THEN
              EXIT
          ENDIF

        ELSE IF (ISECT .EQ. 0) THEN
          IF (ECONO .EQ. 'DFAULT') CALL ERROR(ERRKEY,35,FILEGC,LNUM)
          ECONO = 'DFAULT'
          REWIND(LUNECO)
          LNUM = 0
        ENDIF
      ENDDO

      CLOSE (LUNECO)

      CANHT = 0.0
      CANWH = 0.0
      ! In RUNINIT of TG_VEGGR.for, only call TG_Canopy when RUN=1 
!***********************************************************************
!***********************************************************************
!     SEASONAL INITIALIZATION 
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. SEASINIT) THEN
!-----------------------------------------------------------------------
      IF (CONTROL%RUN .EQ. 1 .OR. INDEX('QPF',CONTROL%RNMODE).LE.0) THEN
        CANHT = 0.0
        CANWH = 0.0
      Endif

!***********************************************************************
!***********************************************************************
!     EMERGENCE CALCULATIONS - Performed once per season upon emergence
!         or transplanting of plants
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. EMERG) THEN
!-----------------------------------------------------------------------
        IF (CONTROL %RUN == 1 .OR. INDEX('FPQ',CONTROL%RNMODE) < 0) THEN
          CANHT  = TABEX(YVSHT,XVSHT,VSTAGE,10)       
          CANWH  = TABEX(YVSWH,XVSHT,VSTAGE,10)
          If ( PLME .EQ. 'T')
           ! The data of I File will over write the data of *.spe 
     &    CALL Read_IFile (CONTROL, CHTI=CANHT, CWIDI=CANWH)
          ! CHTI = Initial canopy height m
          ! CWIDI= Initial canopy width m
        endif       
!***********************************************************************
!***********************************************************************
!     DAILY RATE/INTEGRATION
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. INTEGR) THEN
!-----------------------------------------------------------------------
C     Calculate effect of temperature on canopy expansion, HWTEM
C-----------------------------------------------------------------------
      HWTEM = 0.0
      DO I = 1, 24
        HWTEM = HWTEM + TABEX(YHWTEM,XHWTEM,TGRO(I),5)
      ENDDO
      HWTEM = HWTEM / 24.
      ! Jin Debug: RUN=2, day 1, HWTEM = 0.763, Day2, HWTEM=0.608
C-----------------------------------------------------------------------
C     Calculate effect of day's PAR on canopy expansion, HPAR.
C     ASSUME THAT UPPER 30% OF CANOPY SHADES THE GROWING POINT
C     WPAR IS EFFECT ON WIDTH.  SHADE DOES NOT MAKE MUCH WIDER. LOWER K?
C-----------------------------------------------------------------------
C     IF (XLAI .GT. 0.1) THEN
C        PARNOD = PAR * EXP(-KCAN*0.3*(XLAI-0.1))
C     ELSE
C        PARNOD = PAR
C     ENDIF
C-----------------------------------------------------------------------
      PARNOD = PAR * EXP(-KCAN*(0.3*XLAI))
      ! JZW debug, RUN=2, day 1, PAR = 20, KCAN = 0.6, XLAI=6.12, PARNOD = 6.65
      ! HPAR = 1.67, WPAR =1, RCANHT=0, RVSTAGE=0, CANHT=2.4
      ! Day 2. PAR = 5.5,KCAN = 0.6, XLAI=6.12, PARNOD =1.828
      !  HPAR = 3.26, WPAR = 1.9, RVSTGE = 13.97, RCANHT = 3.8, CANHT =6.29
      HPAR = TABEX(YHWPAR,XHWPAR,PARNOD,8)
      WPAR = TABEX(YHWPAR,XHWPAR,PAR,8)
C-----------------------------------------------------------------------
C     Calculate rate of increase in canopy height and update height, CANHT
C-----------------------------------------------------------------------
      !JZW, VStage should be continue from this year to next year, except when prone
   !   if ((CONTROL %RUN .GE. 2) .and. (CONTROL % DAS .eq. 2)) 
   !  &  RVSTGE = 0. ! JZW temp to fix bug in TG_VSTAGES, The canopy hight up or down with shap change on 2nd day of each year
      RCANHT= RVSTGE * TABEX(YVSHT,XVSHT,VSTAGE,10) * HWTEM *
     &  TURFAC * HPAR * RHGHT
      CANHT = CANHT + RCANHT

!     Set minimum Canopy height based on lookup function
      CANHT = MAX(CANHT, TABEX(YVSHT,XVSHT, 0.0, 10))

C-----------------------------------------------------------------------
C     Calculate rate of increase in canopy width and update width, CANWH
C     RWIDTH,RHGHT are used to normalize other crops to the values in tables
C     Values of RHGHT and RWIDTH = 1.00 are for Florunner peanut variety
C     1/22/03 KJB - Don't allow reduction in vstage to reduce canopy
C       width.
!-----------------------------------------------------------------------
      RCANWH= MAX(0.0,RVSTGE) * TABEX(YVSWH,XVSHT,VSTAGE,10) * HWTEM *
     &  TURFAC * WPAR * RWIDTH
      CANWH = CANWH + RCANWH

!     Set minimum Canopy width based on lookup function
      CANWH = MAX(CANWH, TABEX(YVSWH, XVSHT, 0.0, 10))  
      CANWH = MIN(CANWH,ROWSPC)

!***********************************************************************
!***********************************************************************
!     END OF DYNAMIC IF CONSTRUCT
!***********************************************************************
      ENDIF
!***********************************************************************
      RETURN
      END ! SUBROUTINE TG_CANOPY
!=======================================================================
! TG_CANOPY Definitions:  updated 25 Feb 2004
!-----------------------------------------------------------------
! C255      255-character record read from file 
! CANHT     Canopy height (m)
! CANWH     Canopy width normal to row (m)
! ECONO     Ecotype code - used to match ECOTYP in .ECO file 
! ECOTYP    Ecotype code for this simulation 
! ERR       Error code for file operation 
! FILECC    Path plus filename for species file (*.spe) 
! FILEGC    Pathname plus filename for ECO file 
! FOUND     Indicator that good data was read from file by subroutine FIND 
!             (0 - End-of-file encountered, 1 - NAME was found) 
! HPAR      Effect of day's PAR on canopy expansion 
! HWTEM     Effect of temperature on canopy expansion 
! ISECT     Indicator of completion of IGNORE routine: 0 - End of file 
!             encountered, 1 - Found a good line to read, 2 - End of 
!             Section in file encountered denoted by * in column 1. 
! KCAN      Canopy light extinction coefficient for daily PAR, for 
!             equidistant plant spacing, modified when in-row and between 
!             row spacing are not equal 
! LINC      Line number of input file 
! LNUM      Current line number of input file 
! LUNCRP    Logical unit number for FILEC (*.spe file) 
! LUNECO    Logical unit number for FILEE (*.eco file) 
! PAR       Daily photosynthetically active radiation or photon flux 
!             density (moles[quanta]/m2-d)
! PARNOD    Effective PAR at growth point (moles[quanta]/m2-d)
! RCANHT    Rate of increase in canopy height (m/d)
! RCANWH    Rate of increase in canopy width (m/d)
! RHGHT     Relative height of this ecotype in comparison to the standard 
!             height per node (YVSHT) defined in the species file (*.SPE) 
! ROWSPC    Row spacing (m)
! RVSTGE    Rate of VSTAGE change (nodes/day)
! RWIDTH    Relative width of this ecotype in comparison to the standard 
!             width per node (YVSWH) defined in the species file (*.SPE) (m)
! SECTION   Section name in input file 
! TGRO(I)   Hourly canopy temperature (°C)
! TURFAC    Water stress factor for expansion (0 - 1) 
! VSTAGE    Number of nodes on main stem of plant (nodes)
! WPAR      Effect of PAR on canopy width 
! XHWPAR(I) PAR values for table look-up for modifying height and width 
!             growth rate, particularly to allow etiliolation at low PAR 
!             values (mol/day)
! XHWTEM    Temperatures in a table look-up function for modifying height 
!             and width growth rates (°C)
! XLAI      Leaf area (one side) per unit of ground area
!            (m2[leaf] / m2[ground])
! XVSHT     Node number on main stem for use in computing height and width 
!             growth rates 
! YHWPAR(I) Relative increase in height and width growth rates with low PAR 
!             as given in XHWPAR 
! YHWTEM(I) Relative (0-1) expansion in height and width with temperatures 
!             given in XHWTEM 
! YVSHT     Length of internode (m) Vs position on the main stem defined by 
!             XVSHT (m/node)
! YVSWH     Increase in canopy width per node developed on the main stem
!            (m/node)
!***********************************************************************
!      END SUBROUTINE CANOPY
!=======================================================================

