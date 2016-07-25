      subroutine Metro(width,opte,optc,fkept,ff,pp)
	include 'defs.inc'
c Metropolis algorithm to calculate effort producing largest expected catch
	real*8 totlike !cumulates total likelihood of chain
	real*8 zlike,olike,r
	nburn=200
	nkeep=1000
	totlike=0.
	dprob=0.

c	width=.03 !width of hypercube used in proposal distribution

	do i=1,nd
	  deffort(i)=float(i)*1.
      enddo
	dcatch=0.

	  pkeep(1)=alpha
	  pkeep(2)=beta*10000.
	  pkeep(3)=q*100.
	  pkeep(4)=zM
	  pkeep(5)=rho
	  pkeep(6)=z1+1.
	  pkeep(7)=z2+1.

	  pvals=pkeep !old values for likelihood and catch
	  call ff(ssq)
	  olike=-ssq
c	  olike=exp(-0.5*ssq/sigb**2)
	  do i=1,nd
	    pe=deffort(i)
          oldcatch(i)=-1.*pp(pe)
	  enddo	    


	fkept=0. !fraction of steps accepted
	do iter =1,nburn+nkeep !chain values that are kept
c        write(*,*) iter
	  do i=1,np !generate proposal from uniform hypercube
	    if(inout(i).eq.1) then
	      call rnun(1,x)
	      pvals(i)=pkeep(i)-width/2.+x*width
          endif
        enddo
	  
	  call ff(ssq)
	  zlike=-ssq
	  itake=0
	  if(zlike.gt.olike) then
	    itake=1
	  else
	    r=zlike-olike
		call rnun(1,u)
	    u=log(u)
		if(u.lt.r) itake=1
	  endif
	  
	  if(itake.eq.1) then !take proposed point
	    if(iter.gt.nburn) fkept=fkept+1./float(nkeep)
		olike=zlike	
		pkeep=pvals	 
	  	        
		if(iter.gt.nburn) totlike=totlike+zlike

	    do i=1,nd
	      pe=deffort(i)
            c=-1.*pp(pe)
            if(c.lt.0.001) exit !if catch has gone to zero, don't bother checking larger efforts	      
            oldcatch(i)=c
		  if(iter.gt.nburn) then
		    dcatch(i)=dcatch(i)+c*zlike
            endif
	    enddo
c        if(ny.eq.nym-5) write(70,'(200f12.2)') olike,pvals(1:7),oldcatch

	  else !keep previous point in chain
		if(iter.gt.nburn) then
	      totlike=totlike+olike
	      do i=1,nd
		    dcatch(i)=dcatch(i)+oldcatch(i)*olike
		  enddo
          endif

	  endif !itake	  		    
c        write(*,'(i12,9f12.2,i3)') iter,pvals(1:7),zlike,olike,itake
c        write(50,'(i12,9f12.4,i3)') iter,pvals(1:7),zlike,olike,itake
	enddo !iter

      optc=0. !calculate optimal effort as that producing largest expected catch
	do i=1,nd
	  dcatch(i)=dcatch(i)/totlike
	  if(dcatch(i).gt.optc) then
	    optc=dcatch(i)
		opte=deffort(i)
	  endif	
	enddo  
	
c	write(40,'(200f12.2)') fkept,dcatch 

	return
	end
