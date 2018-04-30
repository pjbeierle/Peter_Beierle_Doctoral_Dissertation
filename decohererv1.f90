	PROGRAM Path_Integral_Decoherence
	!Author:Peter Beierle
	!Last Revised 10/13/2015
	implicit none

	integer Nf,Nl, Ng,Ns,Nscr,Ngg,Nd,endrun
	parameter(Nf=24,Nl=1500,Ng=12000,Ns=1500,Nscr=1500,Nd=5000)
	parameter(Ngg=5000,endrun=100)
   	integer j,i,inc,k,m,numberofslits,nn,i2,incend,jend
	integer(2) ihr,imin,isec,i100th,starttime,endtime
	real programtime
	real*8 hbar,h,c,position,d,norm,me,Pi,gam,ee ,D0
	real*8 lambdae,sourcesize,dsource,ve,D2,D3,phaseslope,phasestep
	real*8 Dsg,Dgg, gndglsssize,dgndglss
	real*8 Energy,kernel,screen3size,dscreen3,lenssize,dlens
	real*8 dist1,dist2,lngh,image(Ngg)
	complex*8 kern3(Nl,Ng),kern2(Ng,Ngg),kern(Ngg,Ns),kern0(Ns,1)
	complex*8 kern4(Nscr,Nl)
	complex*8 phiin(1),phi2(Nl),amp(Ng),phi0(Ns)
	complex*8 phigndglss(Ngg) 
	complex*8 cmI,phioutx,imagephase(Ngg)
	complex*8 phiscreen(Nscr),phi1(Ng)
	real*8 probscreen(Nscr),probscreenmem(Nscr,Nf)
	real*8 probscreenaverage(Nscr)
	real*8 probgrating(Ng),probgratingmem(Ng,Nf)
	real*8 probgratingaverage(Ng),finalprob2(Nl,endrun)
	real*8 x,prob2(Nl)
	real*8 dsecondslit,secondslit
	real*8 phasecurveglobal
	real*8 dfirstslit,firstslit
	complex*8 matrix(Nscr,Nscr),finalmatrix(Nscr,Nscr)
	complex*8 matrix2(Nl,Nl),finalmatrix2(Nl,Nl)
	complex*8 matrix3(Ns,Ns),finalmatrix3(Ns,Ns)
	real*8 sigma2(1000),acumphase(1000),positionphase,returnvalue
	integer runner,runnumber,idum,nnmodnumofpnts

	!for random phase
	integer Q,p,numofpoints
	parameter(Q=123456)
	real*8 dephasingsize,value(1000),sigma(1000)
	!call rnset(Q)

	dephasingsize=0.7
	runnumber=1


	!do i=1,Nscr
	! do j=1,Nscr
	!    finalmatrix(i,j)=0d0
	! enddo
	!enddo

	!do i=1,Nl
	! do j=1,Nl
	!    finalmatrix2(i,j)=0d0
	! enddo
	!enddo

	!do i=1,Ns
	! do j=1,Ns
	!    finalmatrix3(i,j)=0d0
    ! enddo
	!enddo

	open(unit=17,file='phase.dat') !phase written



	ee=1.6d-19
	Pi=3.14159265d0
	h=6.626d-34
	hbar=h/2d0/Pi
	c=3.d8
	d=100.d-9
	me=9.11d-31
	Energy=1670d0*ee+me*c*c
	!ve=dsqrt(2.*Energy/(.511E6/c**2.))	 !non relativistic leftover
	gam=Energy/(me*c*c)
	ve=sqrt(1d0-1d0/(gam*gam))*c
	lambdae=h/(gam*me*ve)
	write(6,*) "The gamma factor is  ",gam
	write(6,*) "The electron velocity is  ",ve

	firstslit=12.7d-6
	write(6,*) "The first slit is ",firstslit, " wide."
	secondslit=2.5d-6 !10.22E-6 !14.59E-6
	write(6,*) "The second slit is  ",secondslit, "   wide."
	lenssize=20000d-9
	gndglsssize=20d-6
	write(6,*) "the ground glass size is ",gndglsssize," wide"
	write(6,*) "The lenssize is  ",lenssize, "   wide."
	sourcesize=15000d-9			 !YOU NEED TO CHANGE THE NUMBER OF SLITS!
	write(6,*) "The grating source is  ",sourcesize, "   wide."
	screen3size=500d-6
	write(6,*) "The detecttion screen is  ",screen3size, "   wide."

	dfirstslit=firstslit/Nf
	dsecondslit=secondslit/Ns
	dlens=lenssize/Nl
	dsource=sourcesize/Ng
	dscreen3=screen3size/Nscr
	dgndglss=gndglsssize/Ngg
	D0=.24d0		!distance from the first slit to the second slit
	Dsg=.06d0-.005d0		!0.3
	write(6,*) "The dis from the 2nd slit to ground glass is ",Dsg
	Dgg= .005d0

	D2=.005d0	!distance at UNL
	write(6,*) "The distance from the grating to the lens is  ",D2

	D3=0.24d0	!distance at UNL   0.231
	write(6,*) "The distance from the lens to the screen is  ",D3


	!call gettim(ihr,imin,isec,i100th)
	!starttime=0.
	!starttime=ihr*60*60+imin*60+isec

	cmI=(0,1)	!the complex number I




!				The next subroutines sets up the path integral kernels
   write(6,*) "01"
   call firstkernalsubroutine(kern,Ngg,Ns,gndglsssize,dgndglss,secondslit,dsecondslit,Dsg,lambdae)
   call secondkernalsubroutine(kern2,Ng,Ngg,sourcesize,dsource,gndglsssize,dgndglss,Dgg,lambdae)
   call thirdkernalsubroutine(kern3,Nl,Ng,lenssize,dlens,sourcesize,dsource,D2,lambdae)
   call fourthkernalsubroutine(kern4,Nscr,Nl,screen3size,dscreen3,lenssize,dlens,D3,lambdae)
   

!				The next subroutine sets up the amplitude at the grating
   call gratingamplitude (amp,sourcesize,dsource,Ng)




	do runnumber=1,endrun
	do i=1,Nl
	finalprob2(i,runnumber)=0d0
	enddo
	enddo

	do runnumber=1,endrun

	do i=1,200 !value is an array of 100 random numbers from 0 to 1
	call ran0(idum,returnvalue)
	 value(i)=returnvalue
	enddo


	numofpoints=25 !should be a factor of 5000
	p=1

	do k=1,Ngg
	 nn=k
	 image(nn)=dephasingsize*2.0*Pi*value(p)
	 imagephase(nn)=cmI*dsin(image(nn))+dcos(image(nn))
	 nnmodnumofpnts=MOD(nn,numofpoints)
	  if (nnmodnumofpnts.eq.0) then
	   p=p+1
	  endif
	enddo

		


	do k=1,Ngg
       nn=k
	 x=(gndglsssize/2d0-dgndglss*nn)
	 write(17,*) x, image(k)
	enddo

	
   
     
     	! ---------start of large incoherent loop-----------------------------
	incend=Nf

	do inc=1,incend

	call zerothkernel(kern0,Ns,firstslit,dfirstslit,inc,secondslit,dsecondslit,D0,lambdae)	

    !This routine propogates from the first slit to the second slit
    call wavefunction1(Ns,kern0,phi0)

    !The next routine propagates from the source (second slit) to the gndglss
	call wavefunction2(phi0,Ngg,Ns,kern,phigndglss,imagephase)											  
	
	!The next routine propagates from the gndglss to the grating
	call wavefunction3(phigndglss,Ng,Ngg,kern2,phi1,amp)

    !The next routine propagates from the grating to the "lens" (just the nearfield)
    call wavefunction4(phi1,Nl,Ng,kern3,phi2)										

	!The next routine propagates from the lens to the screen
	call wavefunction5(phi2,Nscr,Nl,kern4,phiscreen)
    

	!the probability is calculated	next
    call probabilitycalculation(phiscreen,probscreen,probscreenmem,prob2,finalprob2,Nscr,Nl,runnumber,inc,incend,endrun,phi2)	


!	 calculation of matrices
!	call matrixcalculation(phiscreen,phi2,phi0,matrix,finalmatrix,matrix2,finalmatrix2,matrix3,finalmatrix3,Nscr,Nl,Ns)
	
      enddo   !close incoherent loop
write(6,*) 'end of runnumber',runnumber
 enddo !end of runnumber loop

	
	 
!	  the incoherent sum at the screen is calculated next
call screenincoherentsumming(Nscr,incend,inc,probscreenaverage,probscreenmem)

     
     
     	
!	  next write data to file
call writeprobabilities(Nscr,Nl,screen3size,dscreen3,probscreenaverage,runnumber,endrun,finalprob2)

!call writematrices(Nscr,Nl,Ns,finalmatrix,finalmatrix2,finalmatrix3)


!       next find runtime

!	call gettim(ihr,imin,isec,i100th)
!	endtime=0.
!	endtime=ihr*60*60+imin*60+isec

!	programtime=(endtime-starttime)/60.
!	write(6,*) 'Program run time =  ',programtime

666	format(E12.6,X,E12.6) !there was a 666 there
667   format(E12.6,X,E12.6)

	END


!---------------------------------------------------------------------	  
subroutine firstkernalsubroutine(kern,Ngg,Ns,gndglsssize,dgndglss,secondslit,dsecondslit,Dsg,lambdae)
implicit none 
  integer i,j,Ngg,Ns
  real*8 gndglsssize,dgndglss,secondslit,dsecondslit,Dsg,lambdae
  real*8 dist1,dist2,lngh,kernel,Pi
  complex*8 cmI,kern(Ngg,Ns)
  Pi=3.14159265d0
  cmI=(0,1)	!the complex number I

	do i=1,Ngg
		do j=1,Ns
		dist1=(gndglsssize/2d0-dgndglss*i)	 !source stand for grating
		dist2=(secondslit/2d0-dsecondslit*j)
		lngh=sqrt((dist1-dist2)**2d0+(Dsg)**2d0)
		kernel=0d0
		kernel=2d0*Pi*lngh/(lambdae)
		kern(i,j)=cmI*dsin(kernel)+dcos(kernel)
 		enddo
	enddo 
	write(6,*) "02"
return
end

!----------------------------------------------------------------------
 subroutine secondkernalsubroutine(kern2,Ng,Ngg,sourcesize,dsource,gndglsssize,dgndglss,Dgg,lambdae)
 implicit none
 integer i,j,Ng,Ngg
 real*8 sourcesize,dsource,gndglsssize,dgndglss,Dgg,lambdae
 real*8 dist1,dist2,lngh,kernel,Pi
 complex*8 cmI,kern2(Ng,Ngg)
 Pi=3.14159265d0
 cmI=(0,1)	!the complex number I

	do i=1,Ng
	    do j=1,Ngg
	    dist1=(sourcesize/2d0-dsource*i)
	    dist2=(gndglsssize/2d0-dgndglss*j)
		lngh=sqrt((dist1-dist2)**2d0+(Dgg)**2d0)
		kernel=0d0
		kernel=2d0*Pi*lngh/(lambdae)
		kern2(i,j)=cmI*dsin(kernel)+dcos(kernel)
		enddo
      enddo
	write(6,*) "03"
return
end
!---------------------------------------------------------------------
subroutine thirdkernalsubroutine(kern3,Nl,Ng,lenssize,dlens,sourcesize,dsource,D2,lambdae)
implicit none
integer i,j,Nl,Ng
real*8 lenssize,dlens,sourcesize,dsource,D2,lambdae
real*8 dist1,dist2,lngh,kernel,Pi
complex*8 cmI,kern3(Nl,Ng)
Pi=3.14159265d0
cmI=(0,1)	!the complex number I

	do i=1,Nl
		do j=1,Ng
		dist1=(lenssize/2d0-dlens*i)
		dist2=(sourcesize/2d0-dsource*j)
		lngh=sqrt((dist1-dist2)**2d0+(D2)**2d0)
		kernel=0
		kernel=2d0*Pi*lngh/(lambdae)
		kern3(i,j)=cmI*dsin(kernel)+dcos(kernel)
 		enddo
	enddo
	write(6,*) "04"
return
end
!-------------------------------------------------------------------
subroutine fourthkernalsubroutine(kern4,Nscr,Nl,screen3size,dscreen3,lenssize,dlens,D3,lambdae)
implicit none
integer i,j,Nscr,Nl
real*8 screen3size,dscreen3,lenssize,dlens,D3,lambdae
real*8 dist1,dist2,lngh,kernel,Pi
complex*8 cmI,kern4(Nscr,Nl)
Pi=3.14159265d0
cmI=(0,1)	!the complex number I

	do i=1,Nscr
		do j=1,Nl
		dist1=(screen3size/2d0-dscreen3*i)
		dist2=(lenssize/2d0-dlens*j)
		lngh=sqrt((dist1-dist2)**2d0+(D3)**2d0)
		kernel=0d0
		kernel=2d0*Pi*lngh/(lambdae)
		kern4(i,j)=cmI*dsin(kernel)+dcos(kernel)
 		enddo
	enddo	
 	write(6,*) "05"
return
end

!-------------------------------------------------------------------


subroutine gratingamplitude (amp,sourcesize,dsource,Ng)
implicit none
integer numberofslits,jend,k,m,j,nn,Ng
real*8 x,sourcesize,dsource
complex*8 amp(Ng)
	numberofslits=150  !sourcesize/d
	jend=int(Ng/numberofslits/2)
	do k=0,numberofslits-1
	do m=0,1
	do j=1,jend
		nn=j+k*jend*2+m*jend
		x=(sourcesize/2d0-dsource*nn)
		amp(nn)=m  !1. for pure phasegrating
	enddo
	enddo
	enddo
return
end

!-------------------------------------------------------------------

subroutine zerothkernel(kern0,Ns,firstslit,dfirstslit,inc,secondslit,dsecondslit,D0,lambdae)	
implicit none
integer Ns,i,j,inc
real*8 firstslit,dfirstslit,secondslit,dsecondslit,D0,lambdae
real*8 kernel,dist1,dist2,lngh,Pi
complex*8 cmI,kern0(Ns,1)
Pi=3.14159265d0
cmI=(0,1)	!the complex number I
	do i=1,Ns
		do j=1,1
		dist1=(firstslit/2d0-dfirstslit*inc)			  
		dist2=(secondslit/2d0-dsecondslit*i)
		lngh=sqrt((dist1-dist2)**2d0+(D0)**2d0)
		kernel=0d0
		kernel=2d0*Pi*lngh/(lambdae)
		kern0(i,j)=cmI*dsin(kernel)+dcos(kernel)
 		enddo
	enddo
return
end
!-------------------------------------------------------------------

subroutine wavefunction1(Ns,kern0,phi0)
integer Ns,i,i2
complex*8 phiin(1),phi0(Ns),phioutx
complex*8 kern0(Ns,1)
	phiin(:)=1d0
	phi0(:)=0d0
	write(6,*) "1"
	do i=1,Ns
		phioutx=0d0
		do i2=1,1
 		phioutx=phioutx+kern0(i,i2)*phiin(i2)	
		enddo
	phi0(i)=phioutx
	enddo
return
end

!-------------------------------------------------------------------

subroutine wavefunction2(phi0,Ngg,Ns,kern,phigndglss,imagephase)											  
implicit none
integer Ngg,Ns,i,i2
complex*8 phi0(Ns),phigndglss(Ngg),phioutx
complex*8 kern(Ngg,Ns),imagephase(Ngg)
	phigndglss(:)=0d0
	write(6,*) "2"
	do i=1,Ngg
		phioutx=0d0
		do i2=1,Ns
 		phioutx=phioutx+kern(i,i2)*phi0(i2)	
		enddo
	phigndglss(i)=phioutx*imagephase(i)
!	write(6,*) phi1(i)
      enddo
return
end

!-------------------------------------------------------------------
subroutine wavefunction3(phigndglss,Ng,Ngg,kern2,phi1,amp)
implicit none
integer Ngg,Ng,i,i2
complex*8 phigndglss(Ngg),phioutx,phi1(Ng)
complex*8 kern2(Ng,Ngg),amp(Ng)

	phi1(:)=0d0
	do i=1,Ng
	    phioutx=0d0
	    do i2=1,Ngg
	    phioutx=phioutx+kern2(i,i2)*phigndglss(i2)
	    enddo
      phi1(i)=phioutx*amp(i)

    enddo
return
end


!-------------------------------------------------------------------

subroutine wavefunction4(phi1,Nl,Ng,kern3,phi2)
implicit none
integer Nl,Ng,i,j
complex*8 phi2(Nl),phioutx,kern3(Nl,Ng),phi1(Ng)
	write(6,*) "2"
 	phi2(:)=0d0
	do i=1,Nl
	phioutx=0d0
		do j=1,Ng
		phioutx=phioutx+kern3(i,j)*phi1(j)
		enddo
	phi2(i)=phioutx
!	write(6,*) phi2(i)
	enddo
   write(6,*) "3"
return
end

!-------------------------------------------------------------------

subroutine wavefunction5(phi2,Nscr,Nl,kern4,phiscreen)
implicit none
integer Nscr,Nl,i,j
complex*8 phiscreen(Nscr),phi2(Nl),kern4(Nscr,Nl),phioutx

   
 	phiscreen(:)=0d0
	do i=1,Nscr
	phioutx=0d0
		do j=1,Nl
 		phioutx=phioutx+kern4(i,j)*phi2(j)
		enddo
	phiscreen(i)=phioutx
	enddo
	
      write(6,*) "4"
return
end


!-------------------------------------------------------------------

subroutine probabilitycalculation(phiscreen,probscreen,probscreenmem,prob2,finalprob2,Nscr,Nl,runnumber,inc,incend,endrun,phi2)
implicit none 
integer i,Nscr,inc,incend,runnumber,Nl,endrun
real*8 probscreen(Nscr),norm,probscreenmem(Nscr,incend)
complex*8 phiscreen(Nscr),phi2(Nl)
real*8 prob2(Nl),finalprob2(Nl,endrun)


	do i=1,Nscr
		probscreen(i)=phiscreen(i)*conjg(phiscreen(i))
	enddo

	norm=0
	do i=1,Nscr
		norm=norm+probscreen(i)
	enddo

	write(6,*) "before"

	do i=1,Nscr
		probscreen(i)=1/norm*phiscreen(i)*conjg(phiscreen(i))
	enddo

	write(6,*) "after"

	do i=1,Nscr
	probscreenmem(i,inc)=probscreen(i)
	enddo     
	
	!not normalized correctly
    prob2(:)=0d0
	do i=1,Nl
	  prob2(i)=phi2(i)*conjg(phi2(i))
	enddo


	do i=1,Nl
	 finalprob2(i,runnumber)=finalprob2(i,runnumber)+prob2(i)
	enddo

 return
 end

 !----------------------------------------------------------------

subroutine matrixcalculation(phiscreen,phi2,phi0,matrix,finalmatrix,matrix2,finalmatrix2,matrix3,finalmatrix3,Nscr,Nl,Ns)
implicit none	 

integer i,j,Nscr,Nl,Ns
complex*8 phiscreen(Nscr),matrix(Nscr,Nscr),finalmatrix(Nscr,Nscr)
complex*8 phi2(Nl),matrix2(Nl,Nl),finalmatrix2(Nl,Nl)
complex*8 phi0(Ns),matrix3(Ns,Ns),finalmatrix3(Ns,Ns)

    do i=1,Nscr
     do j=1,Nscr
	   matrix(i,j)=phiscreen(i)*conjg(phiscreen(j))
	 enddo
	enddo
	
	do i=1,Nscr
	  do j=1,Nscr
	   finalmatrix(i,j)=finalmatrix(i,j)+matrix(i,j)
	  enddo
	enddo

	do i=1,Nl
	 do j=1,Nl
	   matrix2(i,j)=phi2(i)*conjg(phi2(j))
	 enddo
	enddo

	do i=1,Nl
	 do j=1,Nl
	   finalmatrix2(i,j)=finalmatrix2(i,j)+matrix2(i,j)
	 enddo
	enddo

	do i=1,Ns
	 do j=1,Ns
	    matrix3(i,j)=phi0(i)*conjg(phi0(j))
	 enddo
	enddo

	do i=1,Ns
	 do j=1,Ns
	    finalmatrix3(i,j)=finalmatrix3(i,j)+matrix3(i,j)
     enddo
	enddo
return
end

 !----------------------------------------------------------------

subroutine screenincoherentsumming(Nscr,incend,inc,probscreenaverage,probscreenmem)
implicit none

integer i,inc,incend,Nscr
real*8 probscreenaverage(Nscr),probscreenmem(Nscr,incend)

	  probscreenaverage(:)=0d0
	  do i=1,Nscr
	  do inc=1,incend
		probscreenaverage(i)=probscreenaverage(i)+probscreenmem(i,inc)
	  enddo
	  enddo
return
end

 !----------------------------------------------------------------

subroutine writeprobabilities(Nscr,Nl,screen3size,dscreen3,probscreenaverage,runnumber,endrun,finalprob2)
implicit none
integer i,j,Nscr,Nl,runnumber,endrun
real*8 position,screen3size,dscreen3
real*8 probscreenaverage(Nscr),finalprob2(Nl,endrun)

	open(unit=12,file='PX.dat')	   !probability distribution
	open(unit=13,file='probabilitydephasor.dat') !density matrix at dephasor


	do j=1,Nscr
		position=(screen3size/2d0-dscreen3*j)
		write(12,666) position,probscreenaverage(j)
    enddo

	do runnumber=1,endrun
	 do i=1,Nl
	  write(13,*) runnumber, finalprob2(i,runnumber)
	 enddo
	enddo

666	format(E12.6,X,E12.6) !there was a 666 there

return
end

 !----------------------------------------------------------------

subroutine writematrices(Nscr,Nl,Ns,finalmatrix,finalmatrix2,finalmatrix3)
implicit none

integer i,j,Nscr,Nl,Ns
complex*8 finalmatrix(Nscr,Nscr),finalmatrix2(Nl,Nl),finalmatrix3(Ns,Ns)

	open(unit=14,file='matrixscreen.dat') !density matrix at screen
	open(unit=15,file='matrixlens.dat') !density matrix at lens
	open(unit=16,file='matrixsecondslit.dat')!density matrix at 2nd slit
	

	do i=1,Nscr
	 do j=1,Nscr
	    write(14,667) finalmatrix(i,j)
	 enddo
	enddo


	do i=1,Nl
	 do j=1,Nl
	    write(15,667) finalmatrix2(i,j)
	 enddo
	enddo

	do i=1,Ns
	 do j=1,Ns
	    write(16,667) finalmatrix3(i,j)
	 enddo
	enddo
 
 667   format(E12.6,X,E12.6)

 return
 end


 !----------------------------------------------------------------


subroutine bubblepot(acumphase,sigma2,idum,runnumber)
!use msimsl
implicit none
!this program has its own random number generator to avoid
!using msimsl
integer i,endset,j,runnumber
parameter(endset=40000)
real*8 devmean,devdeviation,x1(endset),x2(endset),sigma(endset),z1(1000),z2(1000)
real*8 ampmean,ampdeviation,y1(endset),y2(endset),amplitude(endset)
real*8 Pi,phase(1000),velocity,hbar,totallength(1000),sigma2(1000)
integer idum
real*8 value
real*8 acumphase(1000)
!integer Q !,p,numofpoints(1000)
!parameter(Q=987654)
!real*8 dephasingsize,value(1000)
!call rnset(Q)

open(unit=11,file='phasetest.dat') !testingthephase




do i=1,1000
acumphase(i)=0.0
enddo

!open(unit=12,file='random.dat')
hbar=6.582e-16
Pi=3.14159265
velocity=2.422E7


devmean=250.E-9
devdeviation=250.E-9
ampdeviation=.35
ampmean=0.

do i=1,1000
  phase(i)=0.0
  totallength(i)=0.0
enddo

!using the box-mueller method
do j=1,1000
do i=1,endset
   call ran0(idum,value)
   x1(i)=value
   call ran0(idum,value)   
   x2(i)=value
   call ran0(idum,value)
   y1(i)=value
   call ran0(idum,value)
   y2(i)=value
   sigma(i)=(devdeviation/(2*sqrt(2.0)))*(sqrt(-2.0*log(x1(i)))*cos(2*Pi*x2(i)))+devmean  
   amplitude(i)=(ampdeviation/(2*sqrt(2.0)))*(sqrt(-2.0*log(y1(i)))*cos(2*Pi*y2(i)))+ampmean 
   phase(j)=phase(j)+(sqrt(Pi)*abs(sigma(i))*amplitude(i)/sqrt(log(16.0)))
   totallength(j)=totallength(j)+sigma(i)
enddo
!write(*,*) j
enddo

do i=1,1000
   call ran0(idum,value)
   z1(i)=value
   call ran0(idum,value)
   z2(i)=value
   sigma2(i)=(devdeviation/(2*sqrt(2.0)))*(sqrt(-2.0*log(z1(i)))*cos(2*Pi*z2(i)))+devmean
enddo

do i=1,1000
 acumphase(i)=phase(i)/(velocity*hbar)
enddo
!do i=1,1000
! write(12,*)  acumphase(i),sigma2(i)
!enddo
if (runnumber.eq.1) then
do i=1,1000
	write(11,*) acumphase(i),sigma2(i)
enddo
endif

return
end

!----------------------------------------------------------------

subroutine randomnumberdistros(acumphase,sigma2,idum,runnumber)
implicit none
! this is a patch subroutine that uses the output of bubblepot and sets the
!fwhm of the accumulated phase as ampdeviation
integer i,runnumber
real*8 devmean,devdeviation,ampmean,ampdeviation,Pi
real*8 acumphase(1000),sigma2(1000),value
integer idum
real*8 x1(1000),x2(1000),y1(1000),y2(1000)

open(unit=11,file='phasetest2.dat') !testingthephase
 devmean=250.0d-9
 devdeviation=250.0d-9
 ampmean=0.0d0
 ampdeviation=1.0335d3*1.1935d0
 Pi=3.14159265d0

do i=1,1000
 call ran0(idum,value)
 x1(i)=value
 call ran0(idum,value)
 x2(i)=value
 call ran0(idum,value)
 y1(i)=value
 call ran0(idum,value)
 y2(i)=value
 sigma2(i)=(devdeviation/(2.0d0*sqrt(2.0d0)))*(sqrt(-2.0d0*log(x1(i)))*cos(2.0d0*Pi*x2(i)))+devmean
 acumphase(i)=(ampdeviation/(2.0d0*sqrt(2.0d0)))*(sqrt(-2.0d0*log(x1(i)))*cos(2.0d0*Pi*x2(i)))+ampmean
enddo

if (runnumber.eq.1) then
do i=1,1000
	write(11,*) acumphase(i),sigma2(i)
enddo
endif

!return
end


!----------------------------------------------------------------
subroutine ran0(idum,returnvalue)
implicit none
integer idum,ia,im,iq,ir,mask
real*8 returnvalue,am
parameter(ia=16807,im=2147483647,am=1.0d0/im)
parameter(iq=127773,ir=2836,mask=123459876)

integer k
idum=ieor(idum,mask)
k=idum/iq
idum=ia*(idum-k*iq)-ir*k
if(idum.lt.0)idum=idum+im
returnvalue=am*idum
idum=ieor(idum,mask)
return
end
