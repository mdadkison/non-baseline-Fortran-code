Subroutine RANDIR(nr,avec,yvec)
  !given Dirchlet parameters avec = a1,a2,...a(nr) returns a Dirichlet distributed random vector
  !calls RNGAM from IMSL, RNSET needs to be called first outside of the subroutine to set the random seed
  ! See Gelman et al text 3rd ed. p. 583 for justification  
  !***** 5/28/2015 Milo Adkison **********  
  use RNGAM_INT  
  implicit none
  integer :: nr,i
  real :: a,y(1)
  real :: avec(nr),yvec(nr)
  do i = 1,nr
    a=avec(i) 
    call RNGAM(a,y)
    yvec(i)=y(1)
  enddo
  yvec=yvec/sum(yvec)
  return
end !subroutine RANDIR
    