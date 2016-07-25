!*******************************************
Module my_mod

contains
    
    
subroutine gendat(idtrue,mixvals)
    use RNMTN_INT 
    implicit none
    integer, intent(out) :: idtrue(:)
    integer, intent(out) :: mixvals(:,:)
    include 'include.txt'
    integer :: loci,maxl
    integer, allocatable :: locimat(:,:)
    real,allocatable :: locjunk(:),junka(:),junkb(:)
    
    !create baseline sample counts of character states
    do ic =1,nchars
      do is =1,nstocks          
        k=nalleles(ic)
        allocate(junka(k))
        allocate(junkb(k))
        basefreqs(is,ic,:) = 0.

         if(basechar(ic).eq.1) then
          junka=truefreqs(is,ic,:)*nbasesamples(is,ic) !Dirichlet parameters = true frequencies times sample size
          call RANDIR(k,junka,junkb)
          basefreqs(is,ic,1:k)=junkb
        else
          basefreqs(is,ic,1:k) =1./float(k) !Dirichlet parameters are equal and sum to 1.0
        endif  
          
        deallocate(junka)
        deallocate(junkb)
        
      enddo !is  
    enddo !ic
100 close(16)
    

!************** This section replaces reading the mixture states in from an external file (commented out below) ********    
    !reandomly draw genotypes for the individuals specified in mixnum()
    do i = 1,maxalleles
      ijunkm(i)=i !array of integer values
    enddo !i  
    
    im=0
    do is = 1,nstocks
      do ii = 1,mixnum(is)
        im=im+1
        idtrue(im)=is
!       !select allele for each character
        do ic = 1,nchars
          k=nalleles(ic)
          allocate(locimat(1,k))
          allocate(locjunk(k))
            
          locjunk=0.  
          locjunk = truefreqs(is,ic,:)
          maxl = maxval(locjunk)
          if(maxl.gt.0.999) then
            locimat(1,:)=loci*locjunk  
          else  
          loci=1
            call RNMTN(loci,locjunk,locimat)
          endif  
          !save result
          mixvals(im,ic) = dot_product(ijunkm(1:k),locimat(1,1:k))
          !write(20,'(4i6,f20.2)') is,im,ii,ic,mixvals(im,ic)
          deallocate(locimat)
          deallocate(locjunk)
        enddo !ic = character
        !output individual
      enddo !im = individual
    enddo !is = stock  
return
end subroutine gendat

end module my_mod