!=================================================================================
module ArrayData
implicit none
real(4), allocatable, public :: Var2D(:,:)
!DEC$ATTRIBUTES array_visualizer :: Var2D
end module 

!=================================================================================
Subroutine ShowArray(Dynamic, ArrayPlot, MinVal, MaxVal, ArrayName, NR, NC, RSEED1) 
use IfCore
use AvFRT    ! AvFRT is the AV module file
USE Cells_2D
Use ArrayData
use ModuleData

implicit none
SAVE

real(4), dimension(0:NR,0:NC) :: ArrayPlot 
real(4) MinVal, MaxVal
integer :: dim2d(2) = 0
integer i,j, dynamic, nr, nc, RSEED1
integer :: viewerId, status
character*3 ArrayName

!=================================================================================
 SELECT CASE(DYNAMIC)
!=================================================================================
 CASE (SEASINIT)
!=================================================================================
!   allocate arrays
    allocate (Var2D(0:NR, 0:NC))

!   Call StartWatch to let the Av lib know we're interested in viewing M.
!   avStartWatch needs to know the location of M, the rank, the shape, and the type.
!   The string "M" is use to identify the array in the viewer.
    dim2d = shape(Var2D)
    call avStartWatch(LOC(Var2D), 2, dim2d, AV_REAL4, "Var2D", status)

!   Set min and max values before creating the graphs.
    Var2D(0,0) = MinVal
    Var2D(NR,NC) = MaxVal

!   Create a heightplot graph
    call avCreateGraph3D("graph1", "plot:HEIGHTPLOT, xysource:Var2D, colorsource:Var2D", status)

!   Create a viewer instance
    print *, "Starting Array Viewer for variable: ", ArrayName

    call avNewViewer(viewerId)
    call avViewerUpdate(viewerId, status)

!   point viewer to the first graph
    call avSetViewerPath(viewerId, "graph:/graph1", status)

    do i = 0, NR
      do j = 0, NC
        Var2D(i,j) = ArrayPlot(i,j)
      enddo
    enddo

    call avUpdate(LOC(Var2D), status)

!   Make it visible
    call avVisible(viewerId, AV_TRUE, status)

!=================================================================================
 CASE (RATE, INTEGR, OUTPUT) 
!=================================================================================

    do i = 0, NR
      do j = 0, NC
        Var2D(i,j) = ArrayPlot(i,j)
      enddo
    enddo

    call avUpdate(LOC(Var2D), status)

    do while (i < RSEED1*1500)
      i = i + 1
    enddo

    
!    pause "Press ENTER to continue"

!=================================================================================
 CASE (SEASEND) 
!=================================================================================
   
   call avCloseViewer(viewerId, status)

   call avEndWatch(Loc(Var2D), status)

   DEALLOCATE (Var2D)

!=================================================================================
 END SELECT
!=================================================================================
Return
End Subroutine ShowArray
!=================================================================================
!=================================================================================
