!****************************************************************************
!
!  PROGRAM: Pella_plus
!
!  PURPOSE:  Modified Pella-Masuda stock mixture estimates
!
!****************************************************************************

program Pella_plus
    !include 'link_fnl_static.h'
    include 'link_fnl_shared.h'
    use RNSET_INT
    use RNGDA_INT
    use UVSTA_INT
    use my_mod
    implicit none
    
! Variables
    !baseline & mixture files
    !integer :: nbasesamples(nstocks,nchars),mixstates(nmix,nchars,maxalleles)
    !real :: basefreqs(nstocks,nchars,maxalleles) !only allows 10 character nstates
    include 'include.txt'
    integer :: ichain,ijunks(nstocks),ijunk1(1),isim,nsims
    integer :: iscount(nstocks,nstocks),i5,i95,imed
    real :: total,stat(15,nstocks),xmed(nstocks),avgqfreqs(nstocks,nchars,maxalleles)
    real(kind=8) :: pmix(nstocks),junks(nstocks)
    real,parameter :: bsum=100.
    integer, allocatable :: mixvals(:,:)
    integer, allocatable :: idtrue(:)
    real :: df,s2,x,xx5(nstocks),xx95(nstocks)
    
! Body of Pella_plus
    call RNSET(12726) !Random seed for IMSL routines 1st files
!    call RNSET(88746) !Random seed for IMSL routines 1st files
  
    call getdata
    write(*,*) nstocks,nchars,nmix
    open(15,FILE=indfile,STATUS='UNKNOWN')
    open(14,FILE=propfile,STATUS='UNKNOWN')
    open(13,FILE=freqfile,STATUS='UNKNOWN')
    open(12,FILE=histfile,STATUS='UNKNOWN')
    open(11,FILE=b1file,STATUS='UNKNOWN')
    open(10,FILE=nbfile,STATUS='UNKNOWN')
    

    allocate(mixstates(nmix,nchars,maxalleles))
    allocate(mixvals(nmix,nchars))
    allocate(idstock(nmix))
    allocate(idtrue(nmix))

    
nsims=1000    
do isim=1,nsims !main simulation loop    
    call gendat(idtrue,mixvals) !generate random baseline and mixture data
    
!Computation
    !initial p's and q's
    qfreqs=basefreqs !initial values for allele frequencies
    pvec=1/float(nstocks) !initial mixture frequencies
    !parameters for priors on baseline allele frequencies
    betapriors=0.
    do ic=1,nchars
      do ia=1,nalleles(ic)  
        junks=basefreqs(:,ic,ia)
        betapriors(ic,ia)=sum(junks)/float(nstocks)
      enddo !ia
    enddo !ic  
         
  !*******************MCMC loop********************   
    allocate(xxx(nsamples,nstocks))
    allocate(xx(nsamples))
    iscount=0
    avgqfreqs=0.
    prophist=0.

    do ichain = 1,2*nsamples !MCMC chain loop
     
      !assign identities to mixture individuals
       avec=1/float(nstocks) !Dirichlet parameters for p's - prior = 1/#stocks  
      do im = 1,nmix !loop through mixture individuals 
        !calculate p(i)*f(X(im)|stock i)
        total=0.  
        do is = 1,nstocks !stock loop
         pmix(is) = pvec(is)
         do ic = 1,nchars !character loop
            i=mixvals(im,ic)  
            pmix(is)=pmix(is)*qfreqs(is,ic,i)  
          enddo !ic
          total=total+pmix(is)
        enddo !is
        !assign a random stock identity to the mixture individual
        pmix=pmix/total !rescale probabilities to sum to 1.0
        i=0
        j=1
        call RNGDA(i,j,pmix,ijunks,junks,ijunk1)
        idstock(im)=ijunk1(1) !stock identity of individual im
        if(ichain.EQ.nsamples+1) write(21,'(6i6,20f8.4)') isim,im,idstock(im),mixvals(im,:),pvec 
        if(ichain.gt.nsamples) iscount(idtrue(im),ijunk1(1))=iscount(idtrue(im),ijunk1(1))+1
        avec(ijunk1(1))=avec(ijunk1(1))+1 !cumulate # fish assigned to each stock for Dirichlet draw of p's
      enddo !im  
      !write(15,'(i10,100i3)') ichain,idstock 
      
      !draw p's for stock frequency in mixture
      call RANDIR(nstocks,avec,pvec)      
      !write(14,'(i10,20f5.2)') ichain,pvec
      !write(*,'(i10,20f5.2)') ichain,pvec
      if(ichain.gt.nsamples) then
        xxx(ichain-nsamples,:)=pvec
        x = 1./(float(nhist-1))
        do is=1,nstocks
          i=int(pvec(is)/x)+1
          if(i.gt.nhist) i=nhist
          prophist(i)=prophist(i)+1.
        enddo !is
      endif !ichain 
      
      !draw q's for allele frequencies in stocks
      do is=1,nstocks
        do ic=1,nchars
          !betas = prior parameters = ybar*bsum
          bvec=betapriors(ic,:)
          k=nalleles(ic)
          do ia=1,k
            bvec(ia)=bvec(ia)+nbasesamples(is,ic)*basefreqs(is,ic,ia) !count of allele in baseline
          enddo !ia  
          !count alleles from this stock in the mixture
          do im=1,nmix
            if(idstock(im).EQ.is) then
              i=mixvals(im,ic)  
              bvec(i)=bvec(i)+1
            endif !idstock=is  
          enddo ! im   
          qvec=0.
          call RANDIR(k,bvec,qvec) !Dirichlet draw
          qfreqs(is,ic,:)=qvec !update stock allele frequencies
          
          !save all values of state 1 for character 1 and non-baseline character for stock 1
          if((is.eq.1).and.(ichain.gt.nsamples)) then
            if(ic.eq.1) then  
              b1freq(ichain-nsamples)=qvec(1)
            elseif(ic.eq.3) then
              nbfreq(ichain-nsamples)=qvec(1)
            endif
          endif  
          !*****************************************************************************
          
          if(ichain.gt.nsamples) avgqfreqs(is,ic,:)=avgqfreqs(is,ic,:)+qvec/float(nsamples)
          !write(13,'(i10,2i4,10f7.4)') ichain,is,ic,qvec
        enddo !ic
      enddo !is  
        
    enddo !ichain    
  !***************** end MCMC loop ************************  
    call UVSTA(xxx,stat)
    write(*,*) 'iter = ',isim
    write(*,'(20f7.2)') stat(1,:)
    write(*,'(20f7.4)') stat(3,:)
    
    imed=nsamples/2
    i5=nsamples/20
    i95=nsamples-nsamples/20
    do j=1,nstocks
      xx=xxx(:,j)  
      call sort(nsamples,xx)
      xmed(j)=xx(imed)
      xx5(j)=xx(i5)
      xx95(j)=xx(i95)
    enddo !j
     write(*,*) 'Median =',xmed
     write(14,'(i10,20f5.2)') isim,xmed,xx5,xx95
     deallocate(xxx)
     deallocate(xx)
     
     do j=1,nstocks !j = true stock
       write(15,'(2i10,4i15)') isim,j,iscount(j,:)
     enddo  
     
     do i=1,nstocks
       do j=1,nchars  
         write(13,'(i10,2i4,10f7.4)') isim,i,j,avgqfreqs(i,j,:)
       enddo !j
     enddo !i  
    
    call sort(nsamples,b1freq)
    junk=b1freq(i95)-b1freq(i5)
    write(11,'(i6,2f10.7)') isim,b1freq(imed),junk
    
    call sort(nsamples,nbfreq)
    junk=nbfreq(i95)-nbfreq(i5)
    write(10,'(i6,2f10.7)') isim,nbfreq(imed),junk
     
  enddo !isim    
  
    prophist=prophist/float(nsims*nstocks)
     do i=1,nhist
       write(12,'(f10.7)') prophist(i)
     enddo !i    
     
    
    deallocate(mixstates)
    deallocate(mixvals)
    deallocate(idstock)
    deallocate(idtrue)

 !pause   
end program Pella_plus

