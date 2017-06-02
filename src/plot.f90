MODULE plot
  USE decimal
  USE util
  IMPLICIT NONE
CONTAINS
  SUBROUTINE save_matrix(data, names, name, dir)
    ! Save data matrix in file
    CHARACTER(LEN=32)             :: name
    CHARACTER(LEN=8)              :: names(:)
    REAL(kind=dp)                 :: data(:,:)
    CHARACTER(LEN=32)             :: name_dat
    INTEGER                       :: iunit1, N, i
    ! dir = 0: series in columns, dir = 1: series in rows
    INTEGER                       :: dir 
    name_dat = TRIM(ADJUSTL(name))//".txt"
    CALL util_get_unit(iunit1)
    N = size(data,1)
    OPEN(iunit1,file=name_dat,   status='replace', action='write')
      if (dir == 0) then
        WRITE(iunit1,*) names
      end if
      DO i = 1,N
        if (dir == 0) then
          WRITE(iunit1,*) data(i,:)
        else
          WRITE(iunit1,*) names(i), data(i,:)
        end if
      END DO
    CLOSE(iunit1)
  END SUBROUTINE

  SUBROUTINE read_matrix(name, data, rows, cols)
    ! Read data matrix in file
    CHARACTER(LEN=32)             :: name
    CHARACTER(LEN=5), ALLOCATABLE :: names(:)
    REAL(kind=dp), ALLOCATABLE    :: data(:,:)
    CHARACTER(LEN=32)             :: name_dat
    INTEGER                       :: iunit1, rows, cols, i
    name_dat = TRIM(ADJUSTL(name))//".txt"
    CALL util_get_unit(iunit1)
    ALLOCATE(data(rows, cols), names(cols))
    data = 0.0_dp
    OPEN(iunit1,file=name_dat,   status='old', action='read')
      READ(iunit1,*) names
      DO i = 1,rows
        READ(iunit1,*) data(i,:)
      END DO
    CLOSE(iunit1)
    DEALLOCATE(names)
  END SUBROUTINE 

  SUBROUTINE plot_results(uu, uinit, xx, name)
    !
    ! Save data in paraview format

    REAL(KIND=dp),INTENT(IN)     :: xx(:)
    REAL(KIND=dp),INTENT(IN)     :: uu(:)
    REAL(KIND=dp),INTENT(IN)     :: uinit(:)
    CHARACTER(LEN=32)             :: name
    !

     CALL plot_paraview1D(uu, uinit, xx, name)

  END SUBROUTINE plot_results

  SUBROUTINE plot_paraview1D(uu, uinit, x, name)
    REAL(KIND=dp),INTENT(IN)     :: x(:)
    REAL(KIND=dp),INTENT(IN)     :: uu(:)
    REAL(KIND=dp),INTENT(IN)     :: uinit(:)
    CHARACTER(LEN=32)             :: name
    CHARACTER(LEN=32)            :: name_dat
    INTEGER                      :: i,j,nod,iunit1
    !
    nod = SIZE(x, 1)
    name_dat = TRIM(ADJUSTL(name))//".vtk"  
    !
    CALL util_get_unit(iunit1)
    OPEN(iunit1,file=name_dat,   status='replace', action='write' )
    !
    ! .vtk file header
    !
    !
    !
    WRITE(iunit1,'(A)') '# vtk DataFile Version 2.0'
    WRITE(iunit1,'(A)') 'Solution'
    WRITE(iunit1,'(A)') 'ASCII'
    !
    WRITE(iunit1,'(A)')'DATASET STRUCTURED_GRID'
    WRITE(iunit1,'(A,1x,i12,1x,i12,1x,i12)')'DIMENSIONS',size(x,1), 1, 1
    !
    ! save mesh
    !
    nod = size(x,1)
    WRITE(iunit1,'(A,1x,i12,1x,A)')'POINTS',nod,'float'
    !
    DO i = 1,size(x,1)
      WRITE(iunit1,'(3(F30.15,2x))') x(i),0.0,0.0
    END DO
    !
    ! save output
    !
    WRITE(iunit1,'(A)') 
    WRITE(iunit1,'(A,1x,i12)')'POINT_DATA', nod
    !
    WRITE(iunit1,'(A)')'SCALARS Uo float  1'
    WRITE(iunit1,'(A)')'LOOKUP_TABLE default'
    DO j=1,size(x,1)
        WRITE(iunit1,'(F30.15)') uinit(j)
    END DO
    !
    WRITE(iunit1,'(A)')'SCALARS ref float  1'
    WRITE(iunit1,'(A)')'LOOKUP_TABLE default'
    DO j=1,size(x,1)
        WRITE(iunit1,'(F30.15)') uu(j)
    END DO
    CLOSE(iunit1)
  END SUBROUTINE plot_paraview1D

END MODULE plot
