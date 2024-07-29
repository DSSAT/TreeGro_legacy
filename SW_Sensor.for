!=======================================================================
!  SW_SensorD, Subroutine
!  Work with PraxSoft sensor depth
!  2D water balance for drip irrigation bed with
!     plastic mulch cover.   
!-----------------------------------------------------------------------
!=======================================================================

      SUBROUTINE SW_SensorD(SOILPROP, CONTROL, Cells, SWV)       !Input
!-----------------------------------------------------------------------
      USE Cells_2D
      USE ModuleData
      IMPLICIT NONE
      SAVE

      TYPE (ControlType), INTENT(IN) :: CONTROL
      Type (CellType) Cells(MaxRows,MaxCols)
      TYPE (SoilType) SOILPROP
      
      REAL, DIMENSION(MaxRows,MaxCols) :: SWV
      CHARACTER*8, PARAMETER :: ERRKEY = 'SensorD'
      CHARACTER*78 MSG(30)
      INTEGER i, j, L, DYNAMIC, NLAYR, LIMIT_2D, iPDAT, iHDAT
      INTEGER, DIMENSION(MaxRows,MaxCols) :: Cell_Type

      REAL, DIMENSION(NL) :: BD, DLAYR, DS
      REAL, DIMENSION(MaxRows,MaxCols) :: CellArea
      REAL, DIMENSION(MaxRows,MaxCols) :: Thick, Width
      LOGICAL FEXIST
      
      REAL, DIMENSION(4,MaxCols) :: SS_SW
      REAL, DIMENSION(4):: SS_SWmax, SS_SWmin, SS_SWavg 
      Integer ERRNUM,LNUM, LUNIO, ISENS, LUNWBLxmlD, YR2, DY2, FurCol1
      CHARACTER*6  SSLayerText(4)
      CHARACTER*12 FILEA
      CHARACTER*30 FILEIO
	CHARACTER*80 PATHEX
	CHARACTER*6, DIMENSION(EvaluateNum) :: OLAB
      CHARACTER*6 X(EvaluateNum)
      INTEGER TRTNUM
      
	CHARACTER*17, PARAMETER :: SWBALxmlD = 'SWSensorD.OUT'
	Real ssd1, ssd2, ssd3, ssd4, SSDS(4)
!     Define headings for observed data file (FILEA)
      DATA OLAB / !Pred.          Obs.   Definition 
  !   & '       ', !1                
  !   & '       ', !2                 
  !   & '       ', !3               
  !   & '       ', !4
  !   & '       ', !5                
  !   & '       ', !6                  
     & 'PDAT   ', !7                 
     & 'HDAT   ', !8 
     & 'SWSD1  ', !9 Soil dep 1 (cm) Depth of first soil water sensor (cm)                    .
     & 'SWSD2  ', !10  Soil dep 2 (cm) Depth of second soil water sensor (cm)                   .
     & 'SWSD3  ', !11  Soil dep 3 (cm) Depth of third soil water sensor (cm)                    .
     & 'SWSD4  ', !12  Soil dep 4 (cm) Depth of fourth soil water sensor (cm) 
     &  34 * '      '/  
      Common / SensorDepth/SSDS
      DYNAMIC = CONTROL % DYNAMIC
      TRTNUM =  CONTROL % TRTNUM
      FILEIO  = CONTROL % FILEIO
!***********************************************************************      
!     Seasonal initialization - run once per season
!***********************************************************************
      IF (DYNAMIC .EQ. SEASINIT) THEN
C-----------------------------------------------------------------------
!       Read FILEIO
        !     Read FILEIO
 !     OPEN (LUNIO, FILE = FILEIO, STATUS = 'OLD', IOSTAT=ERRNUM)
 !     IF (ERRNUM .NE. 0) CALL ERROR(ERRKEY,ERRNUM,FILEIO,0)

 !     READ (LUNIO,'(55X,I5)',IOSTAT=ERRNUM) ISENS; LNUM = 1   
 !     IF (ERRNUM .NE. 0) CALL ERROR(ERRKEY,ERRNUM,FILEIO,LNUM)
 !     READ (LUNIO,'(3(/),15X,A12,1X,A80)',IOSTAT=ERRNUM) FILEA,
 !    &     PATHEX
!      CLOSE (LUNIO)
      !CALL GETDESC(ACOUNT, OLAB, DESCRIP)
!-----------------------------------------------------------------------
!     Read Measured (measured) data from FILEA
!-----------------------------------------------------------------------
      ! CALL READA (FILEA, PATHEX, OLAB, TRT_ROT, YRSIM, X)
      CALL READASensor (FILEA, PATHEX, OLAB, TRTNUM, X)
       READ(X(1),'(I7)')iPDAT
	 READ(X(2),'(I7)')iHDAT
	 READ(X(3),'(F10.0)') ssd1
	 READ(X(4),'(F10.0)') ssd2
	 READ(X(5),'(F10.0)') ssd3
       READ(X(6),'(F10.0)') ssd4
      SSDS(1) = ssd1
      SSDS(2) = ssd2
      SSDS(3) = ssd3
      SSDS(4) = ssd4
      FurCol1 = BedDimension % FurCol1
!***********************************************************************
!     Open output file
!-----------------------------------------------------------------------------
      CALL GETLUN('SWBALxmlD', LUNWBLxmlD)
      INQUIRE (FILE = SWBALxmlD, EXIST = FEXIST)
      IF (FEXIST) THEN
        OPEN (UNIT = LUNWBLxmlD, FILE = SWBALxmlD, STATUS = 'OLD',
     &    POSITION = 'APPEND')
      ELSE
        OPEN (UNIT = LUNWBLxmlD, FILE = SWBALxmlD, STATUS = 'NEW')
  !      OPEN (UNIT = LUNWBLxmlD, FILE = SWBALxmlD, STATUS = 'REPLACE')
        WRITE(LUNWBLxmlD,'("*Daily Soil Water BALANCE OUTPUT FILE in xml
     & project")')
      ENDIF
      CALL HEADER(SEASINIT, LUNWBLxmlD, CONTROL % RUN)
      if (CONTROL % RUN .eq. 1) then
        WRITE (LUNWBLxmlD, '(4(A,I1,A4,F7.2,/))')
     &       (('! Depth of Sensor',L, ' is ', SSDS(L)), L=1,4) 
      Endif
      WRITE (LUNWBLxmlD, '("!",T21,"Maximum Soil Water Content",
     &      T54,"Minimum Soil Water Content",T87, 
     &      "Average Soil Water Content")') ! (" ",L=1,N_LYR-3)
      WRITE (LUNWBLxmlD,1123)
 1123   FORMAT('@YEAR DOY   DAS',
     &    '   LYR1max LYR2max LYR3max LYR4max  LYR1min LYR2min '
     &    'LYR3min LYR4min  LYR1avg LYR2avg LYR3avg LYR4avg')


!***********************************************************************
!     OUTPUT - Daily output
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. OUTPUT) THEN
C-----------------------------------------------------------------------
!-----------------------------------------------------------------
      !  SWV    = CELLS % State % SWV
        Call CalSS_SW(SOILPROP, SWV, SSDS, SS_SW )
        DO i = 1, 4
          SS_SWavg(i) = SS_SW(i,1)
          DO j = 2, FurCol1-1 !NColsTOT = N_Bed_Cols + N_Fur_Cols
            SS_SWmax(i)=AMAX1(SS_SW(i,j-1), SS_SW(i,j))
            SS_SWmin(i) =AMIN1(SS_SW(i,j-1), SS_SW(i,j))
            SS_SWavg(i) = SS_SWavg(i)+ SS_SW(i,j)
          Enddo
          SS_SWavg(i) = SS_SWavg(i)/ (FurCol1 -1)
        Enddo
        !     IF (IDETL .EQ. 'D') THEN
        !Write header for daily output
   !     WRITE (LUNWBLxmlD, '("!",T7,"Soil Layer depths (cm):",7A8,
   !  &      " Soil Layer depths (cm):")') (" ",L=1,N_LYR-3)
        

 
!          WRITE (LUNWBLxmlD,1122,ADVANCE='NO') ("SW",L,"D",L=1,9),"    SW10"
! 1122     FORMAT(9("    ",A2,I1,A1),A8)
!          WRITE(LUNWBLxmlD,'(9("    ",A2,I1,A1),A8)')("TCSW",L,"D",L=1,9)
!     &              ,"    TCSW10"
  !    ENDIF
1300   FORMAT(1X,I4,1X,I3.3,1X,I5, 2X, 3(4F8.3,1x))       !Inflows
        CALL YR_DOY(CONTROL % YRDOY, YR2, DY2)
        WRITE (LUNWBLxmlD,1300) YR2, DY2, CONTROL % DAS, 
     &    SS_SWmax(1), SS_SWmax(2), SS_SWmax(3), SS_SWmax(4),
     &    SS_SWmin(1), SS_SWmin(2), SS_SWmin(3), SS_SWmin(4),    
     &    SS_SWavg(1), SS_SWavg(2), SS_SWavg(3), SS_SWavg(4)
!------------------------------------------------------------
!***********************************************************************
!     END OF DYNAMIC IF CONSTRUCT
!***********************************************************************
      ENDIF
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE SW_SensorD

!=======================================================================
C=====================================================================
!     SW_SensorD VARIABLE DEFINITIONS:
!-----------------------------------------------------------------------
! EvaluateNum Number of evaluation variables. It equals 40 in ModuleDefs
! SWV         Single precision cell soil water content in mm3/mm3  

! Width(Row,Col)Cell width in cm[soil]
!-----------------------------------------------------------------------
!     END SUBROUTINE SW_SensorD
!=======================================================================
!=======================================================================
!  SW_SensorH, Subroutine
!  Work with PraxSoft sensor depth
!  2D water balance for drip irrigation bed with
!     plastic mulch cover.  
!-----------------------------------------------------------------------
!=======================================================================

      SUBROUTINE SW_SensorH(SOILPROP, CONTROL, Cells,SWV,iHr)       !Input
!-----------------------------------------------------------------------
      USE Cells_2D
      USE ModuleData
      IMPLICIT NONE
      SAVE

      TYPE (ControlType), INTENT(IN) :: CONTROL
      Type (CellType) Cells(MaxRows,MaxCols)
      TYPE (SoilType) SOILPROP
      
      REAL, DIMENSION(MaxRows,MaxCols) :: SWV
      CHARACTER*8, PARAMETER :: ERRKEY = 'SensorH'
      CHARACTER*78 MSG(30)
      INTEGER i, j, L, DYNAMIC, NLAYR, LIMIT_2D, iHr, iPDAT, iHDAT
      INTEGER, DIMENSION(MaxRows,MaxCols) :: Cell_Type

      REAL, DIMENSION(NL) :: BD, DLAYR, DS
      REAL, DIMENSION(MaxRows,MaxCols) :: CellArea
      REAL, DIMENSION(MaxRows,MaxCols) :: Thick, Width
      LOGICAL FEXIST
      
      REAL, DIMENSION(4,MaxCols) :: SS_SW
      REAL, DIMENSION(4):: SS_SWmax, SS_SWmin, SS_SWavg 
      Integer ERRNUM,LNUM, LUNIO, ISENS, LUNWBLxmlH, YR2, DY2, FurCol1
      CHARACTER*6  SSLayerText(4)
      CHARACTER*12 FILEA
      CHARACTER*30 FILEIO
	CHARACTER*80 PATHEX
	CHARACTER*6, DIMENSION(EvaluateNum) :: OLAB
      CHARACTER*6 X(EvaluateNum)
      INTEGER TRTNUM
      
	CHARACTER*17, PARAMETER :: SWBALxmlH = 'SWSensorH.OUT'
	Real ssd1, ssd2, ssd3, ssd4, SSDS(4)
!     Define headings for observed data file (FILEA)
      DATA OLAB / !Pred.          Obs.   Definition 
  !   & '       ', !1                
  !   & '       ', !2                 
  !   & '       ', !3               
  !   & '       ', !4
  !   & '       ', !5                
  !   & '       ', !6                  
     & 'PDAT   ', !7                 
     & 'HDAT   ', !8   
     & 'SWSD1  ', !9 Soil dep 1 (cm) Depth of first soil water sensor (cm)                    .
     & 'SWSD2  ', !10  Soil dep 2 (cm) Depth of second soil water sensor (cm)                   .
     & 'SWSD3  ', !11  Soil dep 3 (cm) Depth of third soil water sensor (cm)                    .
     & 'SWSD4  ', !12  Soil dep 4 (cm) Depth of fourth soil water sensor (cm) 
     &  34 * '      '/  
      Common / SensorDepth/SSDS
      DYNAMIC = CONTROL % DYNAMIC
      TRTNUM =  CONTROL % TRTNUM
      FILEIO  = CONTROL % FILEIO
     
!***********************************************************************      
!     Seasonal initialization - run once per season
!***********************************************************************
      IF (DYNAMIC .EQ. SEASINIT) THEN
C-----------------------------------------------------------------------
!       Read FILEIO
        !     Read FILEIO
      OPEN (LUNIO, FILE = FILEIO, STATUS = 'OLD', IOSTAT=ERRNUM)
      IF (ERRNUM .NE. 0) CALL ERROR(ERRKEY,ERRNUM,FILEIO,0)

      READ (LUNIO,'(55X,I5)',IOSTAT=ERRNUM) ISENS; LNUM = 1   
      IF (ERRNUM .NE. 0) CALL ERROR(ERRKEY,ERRNUM,FILEIO,LNUM)
      READ (LUNIO,'(3(/),15X,A12,1X,A80)',IOSTAT=ERRNUM) FILEA,
     &     PATHEX
      CLOSE (LUNIO)
      !CALL GETDESC(ACOUNT, OLAB, DESCRIP)
!-----------------------------------------------------------------------
!     Read Measured (measured) data from FILEA
!-----------------------------------------------------------------------
      ! CALL READA (FILEA, PATHEX, OLAB, TRT_ROT, YRSIM, X)
      CALL READASensor (FILEA, PATHEX, OLAB, TRTNUM, X)
   !    READ(X,'(4F10.0)') SSDS(1), SSDS(2), SSDS(3), SSDS(4)
       READ(X(1),'(I7)')iPDAT
       READ(X(2),'(I7)')iHDAT
       READ(X(3),'(F10.0)') ssd1
       READ(X(4),'(F10.0)') ssd2
       READ(X(5),'(F10.0)') ssd3
       READ(X(6),'(F10.0)') ssd4
      SSDS(1) = ssd1
      SSDS(2) = ssd2
      SSDS(3) = ssd3
      SSDS(4) = ssd4
   !   SSDS(1) = (ssd1 + ssd2 ) /2. !Depth to bottom of sensor layer 1	cm
   !   SSDS(2) = (ssd2 + ssd3 ) /2. !Depth to bottom of sensor layer 2	cm
   !   SSDS(3) = (ssd3 + ssd4 ) /2. !Depth to bottom of sensor layer 3	cm
   !   SSDS(4) = ssd4 + (ssd4 - SSDS(3) ) !Depth to bottom of sensor Layer 4	cm
!      CALL SoilLayerText(SSDS, 4, SSLayerText)
      FurCol1 = BedDimension % FurCol1
!***********************************************************************
!     Open output file
!-----------------------------------------------------------------------------
      CALL GETLUN('SWBALxmlH', LUNWBLxmlH)
      INQUIRE (FILE = SWBALxmlH, EXIST = FEXIST)
      IF (FEXIST) THEN
       OPEN (UNIT = LUNWBLxmlH, FILE = SWBALxmlH, STATUS = 'OLD',
     &    POSITION = 'APPEND')
      ELSE
  !      OPEN (UNIT = LUNWBLxmlH, FILE = SWBALxmlH, STATUS = 'REPLACE')
        OPEN (UNIT = LUNWBLxmlH, FILE = SWBALxmlH, STATUS = 'NEW')
        WRITE(LUNWBLxmlH,'("*Hourly and daily Soil WATER BALANCE OUTPUT 
     &FILE in xml project")')
      ENDIF
      CALL HEADER(SEASINIT, LUNWBLxmlH, CONTROL % RUN)
      if (CONTROL % RUN .EQ.1) then
        WRITE (LUNWBLxmlH, '(A, I7,/))')'!Farmers predicted harvest date
     & is ',iHDAT
        WRITE (LUNWBLxmlH, '(4(A,I1,A4,F7.2,/))')
     &       (('! Depth of Sensor',L, ' is ', SSDS(L)), L=1,4) 
      Endif
      WRITE (LUNWBLxmlH, '("!",T24,"Maximum Soil Water Content",
     &      T57,"Minimum Soil Water Content",T90, 
     &      "Average Soil Water Content")') ! (" ",L=1,N_LYR-3)
      WRITE (LUNWBLxmlH,1123)
 1123   FORMAT('@YEAR DOY   DAS iHr',
     &    '  LYR1max LYR2max LYR3max LYR4max  LYR1min LYR2min '
     &    'LYR3min LYR4min  LYR1avg LYR2avg LYR3avg LYR4avg')


!***********************************************************************
!     DAILY RATE CALCULATIONS
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. RATE) THEN
C-----------------------------------------------------------------------
!-----------------------------------------------------------------
        Call CalSS_SW(SOILPROP, SWV, SSDS, SS_SW )
        DO i = 1, 4
          SS_SWavg(i) = SS_SW(i,1)
          DO j = 2, FurCol1-1 !NColsTOT = N_Bed_Cols + N_Fur_Cols
            SS_SWmax(i)=AMAX1(SS_SW(i,j-1), SS_SW(i,j))
            SS_SWmin(i) =AMIN1(SS_SW(i,j-1), SS_SW(i,j))
            SS_SWavg(i) = SS_SWavg(i)+ SS_SW(i,j)
          Enddo
          SS_SWavg(i) = SS_SWavg(i)/(FurCol1 -1)
        Enddo
        !     IF (IDETL .EQ. 'D') THEN
        !Write header for daily output
   !     WRITE (LUNWBLxmlH, '("!",T7,"Soil Layer depths (cm):",7A8,
   !  &      " Soil Layer depths (cm):")') (" ",L=1,N_LYR-3)
        

 
!          WRITE (LUNWBLxmlH,1122,ADVANCE='NO') ("SW",L,"D",L=1,9),"    SW10"
! 1122     FORMAT(9("    ",A2,I1,A1),A8)
!          WRITE(LUNWBLxmlH,'(9("    ",A2,I1,A1),A8)')("TCSW",L,"D",L=1,9)
!     &              ,"    TCSW10"
  !    ENDIF
1300   FORMAT(1X,I4,1X,I3.3,1X,I5, 2X,I2,1X, 3(4F8.3,1x))       !Inflows
        CALL YR_DOY(CONTROL % YRDOY, YR2, DY2)
        !if (iHr. eq. 0) iHr=24
        WRITE (LUNWBLxmlH,1300) YR2, DY2, CONTROL % DAS, iHr,
     &    SS_SWmax(1), SS_SWmax(2), SS_SWmax(3), SS_SWmax(4),
     &    SS_SWmin(1), SS_SWmin(2), SS_SWmin(3), SS_SWmin(4),    
     &    SS_SWavg(1), SS_SWavg(2), SS_SWavg(3), SS_SWavg(4)
!------------------------------------------------------------
!***********************************************************************
!     END OF DYNAMIC IF CONSTRUCT
!***********************************************************************
      ENDIF
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE SW_SensorH

!=======================================================================
C=====================================================================
!     SW_SensorH VARIABLE DEFINITIONS:
!-----------------------------------------------------------------------
! EvaluateNum Number of evaluation variables. It equals 40 in ModuleDefs
! SWV         Single precision cell soil water content in mm3/mm3  

! Width(Row,Col)Cell width in cm[soil]
!-----------------------------------------------------------------------
!     END SUBROUTINE SW_SensorH
!=======================================================================


!=======================================================================
!  CalSS_SW, Subroutine
!  Work with PraxSoft sensor depth
!  Calculate the soil water content for sensor layers
!=======================================================================
      Subroutine CalSS_SW(SOILPROP, SWV, SSDS, SS_SW )
      Use Cells_2D
      USE ModuleData
      Implicit None
      
      TYPE (SoilType) SOILPROP
      Real SSDS(4), SS_SW(4, NColsTot), Top, Bottom !, Thick
      INTEGER i, j, L, NLAYR, FurCol1 
      REAL HalfFurrow, HalfRow
      Real, Dimension(MaxRows,MaxCols) :: CellArea, SWV
      REAL, DIMENSION(NL) :: DLAYR, DS
      
   !   CellArea = CELLS % Struc % CellArea !CELLS % Struc % Width , CELLS % Struc % Thick
      DS    = SOILPROP % DS 
      NLAYR = SOILPROP % NLAYR 
      FurCol1 = BedDimension % FurCol1
   !   DLAYR = SOILPROP % DLAYR
     
      DO i = 1, 4
         DO j = 1, Furcol1-1
           DO L = 1, NLAYR 
              IF (L == 1) THEN 
                Top = 0.
              ELSE
                Top = DS(L-1)
              ENDIF
              Bottom = DS(L)
              IF ((SSDS(i) > TOP) .and. (SSDS(i) < Bottom) )THEN
!               This sensor depth is within a simulation layer layer; done.
                SS_SW(i, j) = SWV(L, j)
              ELSEIF (SSDS(i) == Bottom) THEN
!               This sensor layer is in between of two simulation layers;  
                SS_SW(i,j) = (SWV(L, j) + SWV(L+1, j))/2
              ELSE
                Cycle
!               This sensor layer is out of current simulation layer;
              Endif
           Enddo
        ENDDO
      ENDDO

      RETURN
      END Subroutine CalSS_SW
!=======================================================================
C=====================================================================
!     SS_SW VARIABLE DEFINITIONS:
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     END SUBROUTINE SS_SW
!=======================================================================
!=======================================================================
!  Subroutine READASeasor
!   Reads measured development and final harvest data from FILEA 
!   and maps measured data to appropriate headers for output to 
!   OVERVEIW.OUT
!-----------------------------------------------------------------------
!  Revision history:
!  08/12/2005 CHP Modified to read "alias" headers for some variables
C  02/09/2007 GH  Add path for fileA
!=======================================================================
      SUBROUTINE READASensor(FILEA, PATHEX, OLAB, TRTNUM, X)

!-----------------------------------------------------------------------
!     READ DEVELOPMENT AND FINAL HARVEST DATA FROM  FILEA
!-----------------------------------------------------------------------
      USE ModuleDefs
      IMPLICIT NONE

      INTEGER TRTNUM,ERRNUM,LUNA,LINEXP,ISECT,NTR,I, J
      INTEGER YR,ISIM
      INTEGER COUNT
!     Headers with aliases -- save column
      INTEGER HWAM, HWAH, BWAM, BWAH, PDFT, R5AT  

      REAL TESTVAL

      CHARACTER*6   OLAB(EvaluateNum), HD
      CHARACTER*6   HEAD(EvaluateNum)
      CHARACTER*6   DAT(EvaluateNum), X(EvaluateNum)  !, ERRKEY
      CHARACTER*12  FILEA
      CHARACTER*78  MSG(10)
	CHARACTER*80  PATHEX
	CHARACTER*92  FILEA_P
      CHARACTER*255 C255

      LOGICAL FEXIST

      FILEA_P = TRIM(PATHEX)//FILEA

C-----------------------------------------------------------------------
C     Initialize measured values to -99 before reading values
C
      X = '   -99'

      CALL GETLUN('FILEA', LUNA)
      LINEXP = 0

      INQUIRE (FILE = FILEA_P, EXIST = FEXIST)

      IF (FEXIST) THEN
        OPEN (LUNA,FILE = FILEA_P,STATUS = 'OLD',IOSTAT=ERRNUM)
        IF (ERRNUM .NE. 0) GOTO 5000
      !  CALL YR_DOY(YRSIM,YR,ISIM)

C       FIND THE HEADER LINE, DESIGNATED BY @TRNO
        DO WHILE (.TRUE.)
          READ(LUNA,'(A)',END=5000) C255
          LINEXP = LINEXP + 1
          IF (C255(1:1) .EQ. '@') EXIT    
        ENDDO

C       FOUND HEADER LINE, SAVE IT IN HEAD AND SEARCH FOR TREATMENT
        DO I = 1,EvaluateNum
          READ(C255,'(1X,A5)') HEAD(I)
          IF (HEAD(I) .EQ. '     ') THEN
            COUNT = I - 1
            EXIT
          ENDIF
          C255 = C255(7:255)
        ENDDO

C       FIND THE RIGHT TREATMENT LINE OF DATA
        DO I = 1,1000
          CALL IGNORE(LUNA,LINEXP,ISECT,C255)
C
C    Return if no matching treatment is found in file FILA
C    No field measured data are necessary to be able to run the
C    model

          IF (ISECT .EQ. 0) GO TO 100
          READ(C255(1:6),'(2X,I4)',IOSTAT=ERRNUM) NTR
          IF (ERRNUM .NE. 0) GOTO 5000
          IF(NTR .EQ. TRTNUM) GO TO 60
        ENDDO
  
  60    CONTINUE
  
C       READ DATA LINE
        DO I = 1,COUNT
          READ(C255,'(A6)',IOSTAT=ERRNUM) DAT(I)
          IF (ERRNUM .NE. 0) GOTO 5000

          !Test for numeric value -- set non-numeric values to -99
          READ(C255,'(F6.0)',IOSTAT=ERRNUM) TESTVAL
          IF (ERRNUM .NE. 0 .AND. 
     &        TRIM(ADJUSTL(HEAD(I))) .NE. 'TNAM') THEN
            DAT(I) = '   -99'
          ENDIF 

          C255 = C255(7:255)
        ENDDO


C       MATCH HEADER WITH DATA
        DO I = 2, COUNT   !Loop thru FILEA headers
          HD = ADJUSTL(HEAD(I))

!         
            DO J = 1, EvaluateNum    !Loop thru crop-specific headers
              IF (OLAB(J) == HD) THEN
                X(J) = DAT(I)
                EXIT
              ENDIF
            ENDDO

        ENDDO
  
 100    CONTINUE
        CLOSE(LUNA)
        RETURN

!       Error handling
 5000   CONTINUE
        X = '   -99'
        WRITE (MSG(1),'(" Error in FILEA - Measured data not used")')
        CALL INFO(1, "READA ", MSG)
      ENDIF

      CLOSE (LUNA)
      RETURN
      END SUBROUTINE READASensor
C=====================================================================
!     READASensor VARIABLE DEFINITIONS:
!-----------------------------------------------------------------------
! ISECT   Data record code (0 - End of file encountered, 1 - Found a good 
!           line to read, 2 - End of Section in file encountered, denoted 
!           by * in column 1
! LNUM    Current line number of input file 
!-----------------------------------------------------------------------
!     END SUBROUTINE READASensor
!=======================================================================