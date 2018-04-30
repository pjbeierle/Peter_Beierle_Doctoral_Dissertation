	PROGRAM Dephaser
	!includes incoherent summing at first slit
	! phaseslope, phasestep can be set at a value, the amp can be turned on and off
	!use msimsl
	!turned off dephas

	implicit none

	integer Nsource,Nin,Nout,Ndet,endrun,p,q,r,Nx,idum
	integer blksum,blknum1,blknum,blknumm !if block exceed range,block number will be changed
	parameter(Nsource=1500,Nin=1000,Nout=1000,Ndet=1500,endrun=500,blknum=2)	

	integer blkwid(blknum),blkwid1 !the first one is for grating, the second one is for source
	real*8 ee,Pi,h,hbar,c,d,a,me,Energy,gam,ve,lambdae
	real*8 source,dsource,phiin,dphiin,phiout,dphiout
	real*8 Dsd,Dds,dephase,image(Nin,endrun),delx,w
	real*8 screen,dscreen
	complex*8 kern2(Ndet,Nout),kern(Nin,Nsource),decfactor(Nin)
	complex*8 amp(Nin) 
	complex*8 cmI,imagephase(Nin)
	complex*8 wavephidet(Ndet),wavephiin(Nin),wavephiout(Nout)
	real*8 probscreenmem(Ndet,endrun)
	integer runnumber


	real*8 dephasingsize
	!call rnset(Q)
	dephasingsize=0.7d0
	
	
	idum=6
	
	
	ee=1.6d-19
	Pi=3.14159265d0
	h=6.626d-34
	hbar=h/2d0/Pi
	c=3.d8
	
	
	d=150d-9 !periodicity of grating
	a=50d-9 !slit width of grating
	
	me=9.11d-31
	Energy=1670d0*ee+me*c*c
	gam=Energy/(me*c*c)
	ve=sqrt(1d0-1d0/(gam*gam))*c
	lambdae=h/(gam*me*ve)
	write(6,*) "The gamma factor is  ",gam
	write(6,*) "The electron velocity is  ",ve

	source=15d-6
	write(6,*) "The source is ",source, " wide."
	screen=8d-4
	write(6,*) "The screen is  ",screen, "   wide."
	dephase=0.5d-6
	write(6,*) "the dephasor ",dephase," wide"

	phiin=dephase
	phiout=dephase

	dsource=source/Nsource
	dscreen=screen/Ndet
	dphiin=phiin/Nin
	dphiout=phiout/Nout
		
	delx=0.25d0*a !shift of gaussian
	w=2d0*a !width of gaussian
																		
	blknum1=1 !seperate source into 25 parts
	blkwid1=Nsource/blknum1 !make sure it is an integer
	Nx=idnint(dephase/(delx))!this is how many gaussians are seperated on the grating,if want to close decoherer, make it 1.
	Nx=1

	Dsd=.24d0		!distanace from the source to the dephasor	
	write(6,*) "The distance from the source to dephasor is ",Dsd

	Dds=.24d0	!distance from the dephasor to the detectorscreen
	write(6,*) "The distance from the dephasor to the detectorscreen is  ",Dds

	cmI=(0,1)	!the complex number I

!				The next subroutines sets up the path integral kernels
   write(6,*) "01"

   call firstkernalsubroutine(kern,Nin,Nsource,phiin,dphiin,source,dsource,Dsd,lambdae)
   call secondkernalsubroutine(kern2,Ndet,Nout,screen,dscreen,phiout,dphiout,Dds,lambdae)  
!				The next subroutine sets up the amplitude at the grating
   call gratingamplitude(amp,phiin,Nin,d,a)	
   call gratingwrite(amp,Nin,phiin,dphiin)
   

   
	do runnumber=1,endrun
		wavephiin(:)=0d0
		probscreenmem(:,runnumber)=0d0
		
!------------random widthblock on grating ------------------------	 
		!this routine create random with block
		call randomblock(idum,blknum,blknumm,Nin,blkwid)!if close blocked source, make the blkwid very big
		

		    
!------------smooth phase-------------------------------!choose one of these potential
		call smoothpotential(dephasingsize,idum,Nin,image,imagephase,runnumber,endrun)
		
!-------------------------------------------------------		

		
		write(6,*) "grating block number is",blknumm
		
			do q=1,blknum1
				
				!This routine propogates from the source to the place right before grating
				write(6,*) "source block no.",q
				
				call wavefunction1(Nin,Nsource,kern,wavephiin,q,blkwid1)
				!plot wavephiin
				call plotwavephiin(Nin,wavephiin)
!------------tiltphase potential ------------------------ !choose one of these potential
		!create tilt potential in each block
				!call tiltphase(idum,imagephase,image,Nin,runnumber,endrun,blknum,blknumm,blkwid)  
				!write phase to check
				call phasewrite(endrun,Nin,image)				
				!dephasor
				call dephasor(wavephiin,Nin,imagephase,Nsource)
				!grating			

				call grating(wavephiin,Nin,amp)
				
				
!There starts the decoherent loop
!------------start of decoherence loop---------------------
!wavephiout will change in every loop			
					
					
				do r=1,Nx
					write(6,*) r,Nx,q,blknum1
					blksum=0 !this number is accumulating blockwidth
					do p=1,blknumm	
						call decoheregenerator(Nin,dephase,decfactor,r,delx,w)

						call decohere(decfactor,wavephiout,wavephiin,Nin,Nout)!to cancel decoherer, you need go into this route and comment out defactor
						!The next routine propagates from after the grating to the detector
						call plotwavephiout(Nout,wavephiout)
						call wavefunction2(Ndet,Nout,kern2,wavephiout,wavephidet,p,blkwid,blksum,blknum)
	
						!the probability is calculated	next
						call probabilitycalculation(wavephidet,probscreenmem,Ndet,runnumber,endrun)
						write(6,*) 'block ',p
					enddo
					write(6,*) blksum 
				enddo !close decoherence loop
			enddo

		write(6,*) 'end of runnumber',runnumber

      		
	enddo
        !next write data to file
		call writeprobabilities(probscreenmem,Ndet,endrun)
		
		
666	format(E20.6,X,E20.6) !there was a 666 there
667 format(E20.6,X,E20.6)

	END







!---------------------------------------------------------------------
	  
subroutine firstkernalsubroutine(kern,Nin,Nsource,phiin,dphiin,source,dsource,Dsd,lambdae)
implicit none 
  integer i,j,Nin,Nsource
  real*8 phiin,dphiin,Dsd,lambdae
  real*8 dist1,dist2,lngh,kernel,Pi,source,dsource
  complex*8 cmI,kern(Nin,Nsource)
  Pi=3.14159265d0
  cmI=(0d0,1d0)	!the complex number I


	do i=1,Nin
		do j=1,Nsource
			dist1=phiin/2d0-dphiin*i+dphiin/2d0						!checked
			dist2=source/2d0-dsource*(j)+dsource/2d0	 !source	!checked
			!write(6,*) 'the coordinates of',j,' source point is ',dist2	
			lngh=sqrt((dist1-dist2)**2d0+(Dsd)**2d0)
			!write(6,*) 'the length of source',j,' to phiin',i,' is',lngh 
			kernel=2d0*Pi*lngh/(lambdae)
			!write(6,*) i,'to',j,kernel
			kern(i,j)=cmI*dsin(kernel)+dcos(kernel)
			!write(6,*) j,' to ',i, real(kern(i,j)),imag(kern(i,j))
 		enddo
		    !write(6,*) 'the coordinates of',i,' phiin point is ',dist1
	enddo 
	write(6,*) "02"
return
end

!----------------------------------------------------------------------

 subroutine secondkernalsubroutine(kern2,Ndet,Nout,screen,dscreen,phiout,dphiout,Dds,lambdae)
 implicit none
 integer i,j,Nout,Ndet
 real*8 screen,dscreen,phiout,dphiout,Dds,lambdae
 real*8 dist1,dist2,lngh,kernel,Pi
 complex*8 cmI,kern2(Ndet,Nout)
 Pi=3.14159265d0
 cmI=(0d0,1d0)	!the complex number I


	do i=1,Ndet
	    do j=1,Nout
			dist1=(screen/2d0-dscreen*i+dscreen/2d0)				
			dist2=(phiout/2d0-dphiout*j+dphiout/2d0)
			lngh=sqrt((dist1-dist2)**2d0+(Dds)**2d0)
			kernel=2d0*Pi*lngh/(lambdae)
			kern2(i,j)=cmI*dsin(kernel)+dcos(kernel)
		enddo
      enddo
	write(6,*) "03"
return
end

!---------------------------------------------------------------------

subroutine gratingamplitude(amp,phiin,Nin,d,a)
implicit none
integer hfdist,hfsltwid,Nin,j
real*8 a,d,phiin
complex*8 amp(Nin)
	amp(:)=0
	
	hfdist=idnint(d/phiin/2*Nin)
	hfsltwid=idnint(a/phiin/2*Nin)
	
	
	
	do j=(Nin/2-hfdist)-hfsltwid,(Nin/2-hfdist)+hfsltwid
		 amp(j)=1
	 enddo
	 do j=(Nin/2+hfdist)-hfsltwid,(Nin/2+hfdist)+hfsltwid
		 amp(j)=1
	 enddo	
	
	
return
end


!---------------------------------------------------------------------
!modified grating to remove unsymmetric phenomenon

! subroutine gratingamplitude (amp,Nin)
! implicit none
! integer	jend,k,m,j,nn,Nin
! complex*8 amp(Nin)

	! jend=25
	! do k=0,199
		! do m=0,1
			! do j=1,jend
				! nn=j+k*jend*2+m*jend
				
				! amp(nn)=m
			! enddo
		! enddo
	! enddo
	! do j=1,jend
		! nn=j+10000
	    ! amp(nn)=0d0
	! enddo
! return
! end
!----------------------------------------------------------------
subroutine gratingwrite(amp,Nin,phiin,dphiin)
implicit none
integer Nin,i
real*8 phiin,dphiin
complex*8 amp(Nin)
			
open(unit=15,file='grating.dat') 
	do i=1,Nin
		write(15,666) (phiin/2d0-dphiin*i), Real(amp(i))
	enddo
666	format(E20.6,X,E20.6)
close(15)

return
end


!----------------------------------------------------------------

subroutine phasewrite(endrun,Nin,image)
integer i,runnumber,Nin,endrun
real*8 image(Nin,endrun)

open(unit=14,file='phase.dat') 
do runnumber=1,endrun
	do i=1,Nin
			write(14,*) runnumber, image(i,runnumber)
	enddo
enddo
666	format(E20.6,X,E20.6)
close(14)

return
end

!-------------------------------------------------------------------

subroutine wavefunction1(Nin,Nsource,kern,wavephiin,q,blkwid1)
implicit none
integer Nin,Nsource,i,i2,q,blkwid1
complex*8 wavephiin(Nin),phisource(Nsource)
complex*8 kern(Nin,Nsource),phiinsum
	phisource(:)=1d0
	wavephiin(:)=0d0
	write(6,*) "1"
		
	do i=1,Nin
		phiinsum=0d0
		do i2=(q-1)*blkwid1+1,q*blkwid1
 			phiinsum=phiinsum+kern(i,i2)*phisource(i2)
			!write(6,*) i2,'to',i,kern(i,i2)	
		enddo
		wavephiin(i)=phiinsum
		!write(6,*) wavephiin(i)		
	enddo
return
end

!-------------------------------------------------------------------
subroutine plotwavephiin(Nin,wavephiin)
implicit none
integer Nin,i 
real*8 phiin_probability(Nin)
complex*8 wavephiin(Nin)
	do i=1,Nin
		phiin_probability(i)=wavephiin(i)*conjg(wavephiin(i))
	enddo
open(unit=18,file='wavephiin.dat') 
		do i=1,Nin
			write(18,666) real(wavephiin(i)), phiin_probability(i)
		enddo
	666	format(E20.6,X,E20.6)
	close(18)

return
end	



!------------------------------------------------------------------- 
subroutine dephasor(wavephiin,Nin,imagephase,Nsource)
implicit none
integer i,Nin,Nsource
complex*8 wavephiin(Nin),imagephase(Nsource)
	do i=1,Nin
		wavephiin(i)=wavephiin(i)*imagephase(i)
	enddo
return 
end


!-------------------------------------------------------------------

subroutine grating(wavephiin,Nin,amp)
implicit none
integer Nin,i
complex*8 wavephiin(Nin),amp(Nin)
    do i=1,Nin
		wavephiin(i)=wavephiin(i)*amp(i)		
		!write(6,*) wavephiout(i)
	enddo
	

	
return
end
!----------------------------------------------------------------
subroutine decoheregenerator(Nin,dephase,decfactor,r,delx,w)
implicit none
integer Nin,i,r
real*8 A,w,Pi,l,dephase,delx
complex*8 decfactor(Nin)

Pi=3.14159265d0

A=delx/sqrt(w*Pi)

l=dephase/Nin


do i=1,Nin
	!if (dexp(-(((i-1)*l-r*w/10)/w)**2).le.1d-6) then
 
	decfactor(i)=A*dexp(-(((i-1)*l-r*delx)/w)**2)

enddo



		 

return
end
!------------------------------------------------------------------- 
subroutine decohere(decfactor,wavephiout,wavephiin,Nin,Nout)
implicit none
integer i,Nin,Nout
complex*8 wavephiin(Nin),wavephiout(Nout),decfactor(Nin)

	


do i=1,Nout
	wavephiout(i)=wavephiin(i)!*decfactor(i)!comment out decfactor to cancel decoherer
	!write(6,*) wavephiout(i)
enddo


	


return
end
!-------------------------------------------------------------------
!plot the wave right after the grating
subroutine plotwavephiout(Nout,wavephiout)
implicit none
integer Nout,i
real*8 phiout_probability(Nout)
complex*8 wavephiout(Nout)
	do i=1,Nout
		phiout_probability(i)=wavephiout(i)*conjg(wavephiout(i))
	enddo

	open(unit=17,file='wavephiout.dat') 
		do i=1,Nout
			write(17,666) Imag(wavephiout(i)), phiout_probability(i)
		enddo
	666	format(E20.6,X,E20.6)
	close(17)

return
end



!-------------------------------------------------------------------
subroutine randomblock(idum,blknum,blknumm,Nin,blkwid)
implicit none 
integer idum,Nin,blknum,i,blknumm
integer blkwid(blknum),blkwidsum
real*8 returnvalue,Pi,y1,x1,x2
blkwidsum=0
blkwid(:)=0
Pi=3.14159265d0

do i=1,blknum-1
  !this is box muller method
  call ran0(idum,returnvalue)
  x1=returnvalue
  call ran0(idum,returnvalue)
  x2=returnvalue
  y1=dsqrt(-2.0d0*dlog(x1))*dcos(2.0d0*Pi*x2)
 !y2=dsqrt(-2.0d0*dlog(x1))*dsin(2.0d0*Pi*x2)
  blkwid(i)=idnint((300.0d0*y1)+20000.0d0)
  blkwidsum=blkwidsum+blkwid(i)
 !in case of blocks exceed the width of grating
	if (blkwidsum.ge.Nin) then	 
		blkwidsum=blkwidsum-blkwid(i)
		blknumm=i
		exit
	endif
  !write(6,*) blkwid(i)
!  write(28,*) amplitude(i), x0(i), sigma(i)  
enddo
  blkwid(blknumm)=Nin-blkwidsum   !the rest points will form the last block
  write(6,*) blkwid


return
end
!----------------------------------------------------------------

subroutine tiltphase(idum,imagephase,image,Nin,runnumber,endrun,blknum,blknumm,blkwid)
implicit none
integer idum,p,i,Nin,runnumber,endrun,blknum,blknumm,blkwid(blknum),x,blkwidsum
real*8 k(blknum),image(Nin,endrun),returnvalue,center
complex*8 imagephase(Nin),cmI

cmI=(0d0,1d0)       !the complex number I
blkwidsum=0
do i=1,blknumm
	call ran0(idum,returnvalue)
	k(i)=(-8d-1)*returnvalue+4d-1
	write(6,*) returnvalue   !slope
enddo

do p=1,blknumm
	do i=blkwidsum+1,blkwidsum+blkwid(p)
	center=(blkwid(p)+1)/2
	x=i-blkwidsum !so it starts from 1 again
	image(i,runnumber)=(x-center)*k(p)
	imagephase(i)=cmI*dsin(image(i,runnumber))+dcos(image(i,runnumber))
	enddo
	blkwidsum=blkwidsum+blkwid(p)
enddo

write(6,*) "tiltphase is created" 
return
end

!----------------------------------------------------------------

subroutine smoothpotential(dephasingsize,idum,Nin,image,imagephase,runnumber,endrun)
implicit none
integer Ngauss,idum,runnumber,endrun
parameter(Ngauss=500)
integer Nin,i,j,x0(Ngauss)
real*8 dephasingsize,image(Nin,endrun),Pi,returnvalue,amplitude(Ngauss)
real*8 x1,x2,y1,sigma(Nin),sumimage(Nin),maxsum
complex*8 imagephase(Nin),cmI
Pi=3.14159265d0
cmI=(0,1)       !the complex number I

do i=1,Nin
 sumimage(i)=0d0
enddo

!open(unit=28,file='debug.dat') !phase written


do i=1,Ngauss
  call ran0(idum,returnvalue)
  amplitude(i)=dephasingsize*2.0d0*Pi*returnvalue
  call ran0(idum,returnvalue)
  x0(i)=idnint(returnvalue*Nin)

  call ran0(idum,returnvalue)
  x1=returnvalue
  call ran0(idum,returnvalue)
  x2=returnvalue
  y1=dsqrt(-2.0d0*dlog(x1))*dcos(2.0d0*Pi*x2)
 !y2=dsqrt(-2.0d0*dlog(x1))*dsin(2.0d0*Pi*x2)
  sigma(i)=(1d0*y1)+8d0
!  write(28,*) amplitude(i), x0(i), sigma(i)  
enddo

do i=1,Ngauss
 do j=1,Nin
   sumimage(j)=sumimage(j)+amplitude(i)*dexp(-((j-x0(i))/(dsqrt(2.0d0)*sigma(i)))**2)
 enddo
enddo

!find the maximum value of sumimage
maxsum=sumimage(1)
do i=1,Nin-1
 if (maxsum.le.sumimage(i+1)) then
 maxsum=sumimage(i+1)
 endif
enddo

do j=1,Nin
   image(j,runnumber)=sumimage(j)*dephasingsize*2*Pi
   imagephase(j)=cmI*dsin(image(j,runnumber))+dcos(image(j,runnumber))
enddo

return
end
!-------------------------------------------------------------------

subroutine wavefunction2(Ndet,Nout,kern2,wavephiout,wavephidet,p,blkwid,blksum,blknum)											  
implicit none
integer Ndet,Nout,i,i2,p,blknum,blkwid(blknum),blksum
complex*8 phidetsum				 
complex*8 kern2(Ndet,Nout),wavephidet(Ndet),wavephiout(Nout)

	
	
	wavephidet(:)=0d0
	write(6,*) "2"
	do i=1,Ndet
		phidetsum=0d0
		do i2=blksum+1,blksum+blkwid(p)
 			phidetsum=phidetsum+kern2(i,i2)*wavephiout(i2)
		enddo
		wavephidet(i)=phidetsum
		!write(6,*) wavephiout(i)
    enddo
	blksum=blksum+blkwid(p)
return
end

!-------------------------------------------------------------------

subroutine probabilitycalculation(wavephidet,probscreenmem,Ndet,runnumber,endrun)
implicit none 
integer i,Ndet,runnumber,endrun
real*8 probscreen(Ndet),norm,probscreenmem(Ndet,endrun)
complex*8 wavephidet(Ndet)
	    probscreen(:)=0d0
	do i=1,Ndet
		probscreen(i)=wavephidet(i)*conjg(wavephidet(i))
	enddo

	norm=0d0
	do i=1,Ndet
		norm=norm+probscreen(i)
	enddo
		write(6,*) norm
	write(6,*) "before"

	do i=1,Ndet
		probscreen(i)=(1/norm)*probscreen(i)		
	enddo
		
	write(6,*) "after"

	do i=1,Ndet
		probscreenmem(i,runnumber)=probscreen(i)+probscreenmem(i,runnumber)
		
	enddo     
		 
 return
 end

 !----------------------------------------------------------------

subroutine writeprobabilities(probscreenmem,Ndet,endrun)
implicit none
integer i,Ndet,runnumber,endrun
real*8 probscreenmem(Ndet,endrun)

open(unit=13,file='probabilitydetector.dat') !probability on detector
	do runnumber=1,endrun
		do i=1,Ndet
			write(13,*) runnumber, probscreenmem(i,runnumber)
		enddo
	enddo
666	format(E20.6,X,E20.6) !there was a 666 there

close(13)

return
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



!----------------------------------------------------------------