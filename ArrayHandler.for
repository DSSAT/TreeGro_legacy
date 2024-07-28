!=======================================================================

      SUBROUTINE ArrayHandler(CELLS, CONTROL, SOILPROP, CellArray, 
     &      ArrayName, MinVal, MaxVal)

!---------------------------------------------------------------------
      USE Cells_2D
      IMPLICIT NONE
      SAVE

      TYPE (ControlType), INTENT(IN) :: CONTROL
      TYPE (SoilType)   , INTENT(IN) :: SOILPROP
      CHARACTER*3,        INTENT(IN) :: ArrayName
      REAL,               INTENT(IN) :: MinVal, MaxVal
      Type (CellType) CELLS(MaxRows,MaxCols)

      CHARACTER*8, PARAMETER :: ERRKEY = 'ArrayOut'
      INTEGER DYNAMIC

      REAL, DIMENSION(MaxRows,MaxCols) :: CellArray   !, RLV_2D

!     For scalable array
      REAL, ALLOCATABLE :: ArrayPlot(:,:)
      INTEGER NR, NC

!The following are needed to control timing for 3D graphics.  Change
!random number seed: higher number, slower graph
      CHARACTER*1 RNMODE
      CHARACTER*6 SECTION
      CHARACTER*30 FILEIO
      INTEGER RSEED1
      INTEGER LUNIO, LINC, LNUM, FOUND
      INTEGER ERRNUM, RUN
      LOGICAL ShowGraph

      DYNAMIC = CONTROL % DYNAMIC

!-----------------------------------------------------------------------
      IF (DYNAMIC .EQ. SEASINIT) THEN
!       Transfer values from constructed data types into local variables.
        FILEIO  = CONTROL % FILEIO
        RUN     = CONTROL % RUN
        RNMODE  = CONTROL % RNMODE
      
!       Read RSEED1, used to control array viewing
        CALL GETLUN('FILEIO', LUNIO)
        OPEN (LUNIO, FILE = FILEIO, STATUS = 'OLD', IOSTAT=ERRNUM)
        IF (ERRNUM .NE. 0) CALL ERROR(ERRKEY,ERRNUM,FILEIO,0)
        SECTION = '*SIMUL'
        CALL FIND(LUNIO, SECTION, LINC, FOUND) ; LNUM = LNUM + LINC
        IF (FOUND .EQ. 0) THEN
          CALL ERROR(SECTION, 42, FILEIO, LNUM)
        ELSE
          READ (LUNIO,'(40X,I6)', IOSTAT=ERRNUM) RSEED1 ; LNUM = LNUM +1
          IF (ERRNUM .NE. 0) CALL ERROR(ERRKEY,ERRNUM,FILEIO,LNUM)
        ENDIF
        CLOSE (LUNIO)
        
        IF (RSEED1 < 10) THEN
          ShowGraph = .FALSE.
        ELSE
          ShowGraph = .TRUE.
        ENDIF
        
        IF (.NOT. ShowGraph) RETURN
      
!       Set up array structure
        NR = INT(soilprop%DS(soilprop%NLAYR) / 5.0)
        NC = INT(BedDimension%ROWSPC_cm / 2.0 / 5.0)
        ALLOCATE (ArrayPlot(0:NR,0:NC))
      ENDIF
!-----------------------------------------------------------------------

      IF (.NOT. ShowGraph) RETURN

!     RLV_2D = CELLS % STATE % RLV

      CALL ScaledArray(DYNAMIC, SOILPROP, CellArray, Cells, NR, NC, 
     &                      ArrayPlot)
      CALL ShowArray(DYNAMIC, ArrayPlot, MinVal, MaxVal, ArrayName, 
     &          NR, NC, RSEED1)

      IF (DYNAMIC == SEASEND) DEALLOCATE (ArrayPlot)

      RETURN
      END SUBROUTINE ArrayHandler

!=======================================================================
C=====================================================================
!    ArrayHandler VARIABLE DEFINITIONS:
!-----------------------------------------------------------------------
! Allocatable arrays are those explicitly declared ALLOCATABLE. An allocatable array may be local to a procedure or may be placed in a module and effectively be global to all procedures of the application. An allocatable array is explicitly allocated with the ALLOCATE statement, and deallocated either explicitly with the DEALLOCATE statement or, if it is a local array for which SAVE has not been specified, automatically upon exit from the procedure. (If SAVE has been specified, local allocatable arrays can persist from one execution of the procedure to the next - they must be explicitly deallocated with a DEALLOCATE statement.) 
! The arrays declared B(10:19,50:100) and C(10,51) have the same shape. The arrays D(100) and E(10,10) do not have the same shape, even though they both contain 100 elements. 
! ERRNUM    Error number for input 
! ERRKEY    Subroutine name for error file 
! FILEIO    Filename for input file (e.g., IBSNAT35.INP) 
! FOUND   Indicator that good data was read from file by subroutine FIND 
! LINC    Line number of input file, Indicates if a line is a good line
! LNUM    Current line number of input file 
! LUNIO     Logical unit number for FILEIO 
! RSEED1  Random number generator seed- user input 
! SECTION   Section name in input file
! NC   # of rows
! NR,  # of columns
!-----------------------------------------------------------------------
!     END SUBROUTINE ArrayHandler
!=======================================================================


!=======================================================================
! Subroutine ScaledArray 
!     Takes 2D soil water array and produces an array which is 
!     spaced in spatially even increments of 5cm width and depth.  This 
!     is so the array viewer will show the animation to scale.
 
      Subroutine ScaledArray(DYNAMIC, SOILPROP, CellArray, Cells, NR,NC,
     &                       ArrayPlot)

      USE Cells_2D
      IMPLICIT NONE
      SAVE

      Type (CellType) Cells(MaxRows,MaxCols)
      TYPE (SoilType) SOILPROP
      REAL, DIMENSION(MaxRows,MaxCols) :: CellArray
      INTEGER DYNAMIC, i, ii, j, jj, NR, NC, NLAYR
      REAL, DIMENSION(NL) :: DS
      REAL ArrayPlot(0:NR,0:NC)
      REAL WidTot
      REAL CumWid(NColsTot)
      INTEGER RowAssoc(0:50)
      INTEGER ColAssoc(0:50)

      DS    = SOILPROP%DS
      NLAYR = SOILPROP%NLAYR

!----------------------------------------------------------------------
      IF (DYNAMIC == SEASINIT) THEN
        WidTot = 0.0
        DO j = 1, NColsTot
          WidTot = WidTot + Cells(BedDimension%FurRow1,j)%Struc%Width
          CumWid(j) = WidTot
        ENDDO

        RowAssoc = 0
        ColAssoc = 0

        RowAssoc(0) = 1
        i = 1
        DO ii = 1, NR
          IF (float(ii) * 5.0 > DS(i)) THEN
            i = i + 1
            i = min(i,NRowsTot)
          ENDIF
          RowAssoc(ii) = i
        ENDDO

        ColAssoc(0) = 1
        j = 1
        DO jj = 1, NC
          IF (float(jj) * 5.0 > CumWid(j)) THEN
            j = j + 1
            j = min(j, NColsTot)
          ENDIF
          ColAssoc(jj) = j
        ENDDO
      ENDIF

!----------------------------------------------------------------------
      ArrayPlot = 0

      DO ii = 0, NR
        DO jj = 0, NC
          ArrayPlot(ii,jj) = CellArray(RowAssoc(ii),ColAssoc(jj))
        ENDDO
      ENDDO

      Return
      End Subroutine ScaledArray

!=======================================================================
C=====================================================================
!     ScaledArray VARIABLE DEFINITIONS:
!-----------------------------------------------------------------------
! NC   # of rows
! NR,  # of columns
!-----------------------------------------------------------------------
!     END SUBROUTINE ScaledArray
!=======================================================================
