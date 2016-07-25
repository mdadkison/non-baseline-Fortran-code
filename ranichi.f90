Subroutine RANICHI(df,s2,x)
  !given scaled inverse chi-square parameters df and s2, returns an inverse chi-square distributed random y
  !to create an inverse chi-square, set s2 = 1.0
  !calls RNGCHI from IMSL, RNSET needs to be called first outside of the subroutine to set the random seed
  ! See Gelman et al text 3rd ed. p. 583 for justification  
  !***** 5/28/2015 Milo Adkison **********  
  use RNCHI_INT  
  implicit none
  real :: x,y(1)
  real :: df,s2
  call RNCHI(df,y)
  x=df*s2/y(1) !from Gelman et al.
  write(*,*) df,s2,y(1),x
return
end !subroutine RANICHI
