PROGRAM DiffCorrect
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  Finite Volume Methods                               !!
  !!                                                      !!
  !!                                                      !!
  !! Author: Paul Mendez                                  !!
  !! Date: 2/June/2017                                    !!
  !!                                                      !!
  !! Version: 0.4                                         !!  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  USE test1
  USE test2
  USE test3
  IMPLICIT NONE

  ! Init threads
  integer t
  call omp_set_num_threads( 4 )
  t = omp_get_max_threads()
  write(*,*)'Number of threads:',t

  ! Run Test(run_reference, run_error)
  CALL test1_run()
  !CALL test2_run()
  !CALL test3_run()
END PROGRAM DiffCorrect