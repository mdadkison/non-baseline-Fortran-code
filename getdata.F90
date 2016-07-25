  subroutine getdata()  
    implicit none
    include 'include.txt'
    character*8 fname
!********* stuff to change ****** (change nchars in include.txt) ************    
!    data basechar/1,1,0/ !if non-baseline characters included
!    data nalleles/4,2,4/ !if non-baseline characters included
    data basechar/1,1/ !if not
    data nalleles/4,2/ !if not
!    4-70 has a mixture of 70% stock 1, 10% each of other stocks. 
!       requires change in setting mixnum in line 24 of this subroutine    
    fname='4-70/HHN' !1st letter =L/M/H for nbasesamples = 20/100/500, 2nd L/H base contrast, 3rd N/L/H/V/P for non-base contrast
    basefile='4-70/HL.txt' ! first low/high is baseline contrast, second low/high/very is non-baseline character contrast
    nmix=500 !low,medium,high base samples =20,100,500
!*******************************************************************    
    
    do ic=1,nchars
      if(basechar(ic).eq.1) then  
        nbasesamples(:,ic)=nmix 
      else
        nbasesamples(:,ic)=1 !non-baseline characters have negligible prior weight
      endif
    enddo !ic  
 
    !mixnum = nmix/nstocks !mixture has equal proportions of each stock
    mixnum(2:nstocks) = nint(float(nmix)/10.)
    mixnum(1) = nmix - mixnum(2)*(nstocks-1)
    mixpriors = 1./float(nstocks) !equal priors
!output file names
    indfile=fname//'.ind' !stock assignments for individuals
    propfile=fname//'.pro' !proportion -> output
    freqfile=fname//'.frq' !character frequencies
    histfile=fname//'.hst' !histogram of proportions
    b1file=fname//'.b1' !freq state 1 character 1 stock 1
    nbfile=fname//'.nb' !freq state 1 nb character stock 1
   
    
!read baseline frequencies
    open(16,file=basefile,ACTION='READ')
    do ic =1,nchars
      do is =1,nstocks
        read(16,*,ERR=100) i,j,(truefreqs(is,ic,k),k=1,nalleles(ic))  
      enddo !is  
    enddo !ic
100 close(16)
 

!    return
    end !subroutine getdata