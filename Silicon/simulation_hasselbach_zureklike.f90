PROGRAM Euler
!====== MPI =====
    use mpi
!================

        IMPLICIT NONE
        common/par/ me,q,k,v0,theta,length
        real*8 me,q,k,v0,theta,finalheight,trajectory(500)
        real*8 length
        real*8 rho_1p5(2001,2001),near_field_x(2001),reduced_diagonal(2001)
        real ratio
        integer s1,s2,flag
        integer counter,sucesses
        integer Nbin
        integer, parameter :: N = 1281
        parameter(Nbin=1100)
        real*8 binvalue(Nbin),gratingarray(32)
        real*8 accum_reduced_diagonal(Nbin,2001)
        real*8 binlower(Nbin),binupper(Nbin),binposition(Nbin)

!============================== MPI =================================
    integer ind
    real*8, dimension(:), allocatable :: y_local
    integer numnodes,myid,rc,ierr,start_local,end_local,N_local
    real*8 allsum

    character(len=8) :: fmt !format descriptor
    character(5) x1
    fmt = '(I5.5)' ! an integer of width 5 with zeros at the left

!====================================================================
!============================== MPI =================================
    call mpi_init( ierr )
    call mpi_comm_rank ( mpi_comm_world, myid, ierr )
    call mpi_comm_size ( mpi_comm_world, numnodes, ierr )

    write(x1,fmt) myid !converting integer to string using an 'internal file'

!    OPEN(UNIT=15,FILE='off_diagonal_dist_scheel'//trim(x1)//'.dat')



    N_local = N/numnodes
!    allocate ( y_local(N_local) )
    start_local = N_local*myid + 1
    end_local =  N_local*myid + N_local
!====================================================================



		OPEN(UNIT=11,FILE='propogation_curvedsurface.dat')
		OPEN(UNIT=12,FILE='trajectory_curvedsurface.dat')
		OPEN(UNIT=13,FILE='debug.dat')
                OPEN(UNIT=14,FILE='electrondistribution_zureklike_lb'//trim(x1)//'.dat')
                OPEN(UNIT=15,FILE='off_diagonal_dist_zureklike_lb'//trim(x1)//'.dat')
		OPEN(UNIT=16,FILE='paths_determination.dat')
		OPEN(UNIT=17,FILE='initial_coherence_distribution.dat')

        call parameters
	call gratinginitialization(gratingarray)
	call initialdensitymatrix(rho_1p5,near_field_x)
	call bininitialization(Nbin,binlower,binupper,binvalue,binposition,accum_reduced_diagonal)
	call write_rho(rho_1p5,near_field_x)

	counter=0
	sucesses=0
!        do s1=1,1901
!        do s2=1,1281
	 do s1=1,1901
         do s2=start_local,end_local
	  call propagate(trajectory,finalheight,s1,s2,counter,sucesses,flag,gratingarray)
	   if (flag.eq.0) then
	    call hasslebachdecoherence(trajectory,reduced_diagonal,rho_1p5)
	    call binning(finalheight,reduced_diagonal,Nbin,binlower,binupper,binvalue,accum_reduced_diagonal,s1,s2)
           endif
	enddo !end of slit2 do loop
        enddo !end of slit1 do loop
       write(6,*) 'total runs= ', counter
       write(6,*) 'total landing hits= ', sucesses
!       ratio=sucesses*100.0/counter
	ratio=sucesses*100.0/1623398.0
       write(6,*) 'percent of beam= ', ratio
	call writeelectrondistribution(binvalue,binposition,Nbin)
	call writeoffdiagonal(binposition,accum_reduced_diagonal,Nbin)
		close(11)
		close(13)
		close(14)
		close(15)
		close(16)
		close(17)

!============================================== MPI ==================================
    call mpi_finalize(rc)
!=====================================================================================

Stop
End Program


!====================================================================================
!====================================================================================
!====================================================================================
!====================================================================================



	
	subroutine parameters
	implicit none
	common/par/ me,q,k,v0,theta,length
	real*8 me,q,k,v0,theta,length
	
	me=9.11d-31
	q=(1.6d-19)*.8425d0
	k=8.988d9
	!v0=1.676d7	!corresponds to a .8keV beam
	v0=2.422d7  !corresponds to a 1.67keV beam
        length=1.00d-2


	return
	end
	 
	
	subroutine propagate(trajectory,finalheight,s1,s2,counter,sucesses,flag,gratingarray)
	implicit none

	real*8 endoftim,delt
	integer ido,tel,i,s1,s2,flag,counter,j,sucesses
	real*8 fcn,t,tend,tft,y(4)
	common/par/ me,q,k,v0,theta,length
	real*8 me,q,k,v0,yprime(4),tfinal
	real*8 length,surfheight,theta,finalheight
	real*8 slit1position,slit2position,cube
	external fcn,divprk,sset
	real*8 trajectory(500)!,sucesstraj(40000,500)
	real*8 gratingarray(32)
	real*8 chargeforce,deflection

	do i=1,500
	   trajectory(i)=0d0
	enddo


	t=0.0
	do i=1,4
	   y(i)=0d0
	enddo

	y(1)=0d0 		!initial x coordinate
	
	counter=counter+1 
	flag=0
	
	slit1position=(-9.5d-6)+(s1-1)*(1d-8)
	slit2position=(-6.4d-6)+(s2-1)*(1d-8)

                         
	theta=datan((slit2position-slit1position)/(25d-2))

!	deflection=
!     	surfheight=-.59d-6 !-1d0+dsqrt(1d0+((y(1)-5.0d-4)*(y(1)-5.0d-4)))+3d-6
	surfheight=-.81d-6


      y(1)=0d0
      y(2)=v0*dcos(theta)		!initial x velocity of beam	
      y(3)=slit1position+(31d-2)*((slit2position-slit1position)/(25d-2)) 
	!initial y coordinate of beam,curved (1.67keV)
      y(4)=v0*dsin(theta) 		!initial y velocity of beam
	
!	write(11,111) y(3)

!	do i=1,16
!	 if (y(3).ge.gratingarray(2*i-1)) then
!	  if (y(3).le.gratingarray(2*i)) then
!		flag=1
!	  endif
!	 endif
!	enddo

!	chargeforce=deflection*me*v0*v0/(length*23.0d-2)

	yprime(1)=y(2)   !y(1)=x, y(2)=vx
	yprime(2)=0d0
	yprime(3)=y(4)   !y(3)=y, y(4)=vy	 
	yprime(4)=-k*q*q/(4*(y(3)+surfheight)*(y(3)+surfheight)*me)
	
	

	     
	ido=1
	tel=0   
	tft=0.0d0
	endoftim=length/v0 

if (s1.eq.1) then
  if (s2.eq.1) then
        write(6,*) 'initial=', surfheight
  endif
endif
	
      delt=endoftim/500d0 
	do i=1,500

          tft=tft+endoftim/500d0 
	  tend=tft						  
	 
	    y(1)=y(1)+(y(2)*delt)
	   !y(2)=y(2)				
	    y(3)=y(3)+(y(4)*delt)
	    y(4)=y(4)+(yprime(4)*delt)
          yprime(4)=-k*q*q/(4*(y(3)+surfheight)*(y(3)+surfheight)*me)
	    surfheight=-.81d-6 !+(4*3.33333d-6)*(tft/endoftim)
		trajectory(i)=y(3)+surfheight
	  if (y(3)+surfheight.le.5d-8) then  
	    flag=1
	  endif
	enddo
if (s1.eq.1) then
  if (s2.eq.1) then
	write(6,*) 'final=', surfheight 
  endif
endif
	!travel another 23cm, free propogation
	
       	tfinal=(23d-2)/v0
	  
	   y(1)=y(1)+(y(2)*tfinal)
	   y(3)=y(3)+(y(4)*tfinal)

	if (flag.eq.0) then	         

!	 write(11,111) y(3)
	 finalheight=y(3)
	 sucesses=sucesses+1
!          do i=1,500
!           write(12,112) trajectory(i)
!          enddo
	endif

!	write (*,'(A,2x,I5)') '+ Number= ',  counter 
	
      !write(6,*) counter

!	write(6,*) 'finished loop'
!	write(6,*) 'total runs= ', counter
!	write(6,*) 'total landing hits= ', sucesses
!	ratio=sucesses*100.0/counter
!	write(6,*) 'percent of beam= ', ratio
	
      


111   format(1E15.8)
112	format(1E15.8)
113   format(3E15.8,1I12.6)   	
      return
	end

!----------------------------------------------------------------------

subroutine initialdensitymatrix(rho_1p5,near_field_x)
implicit none


real*8 rho_1(2001,2001),b,x,y,near_field_x(2001)
real*8 rho_1p5(2001,2001),fwhm1p67,beta,w,sigma
integer i,j

do i=1,2001
 do j=1,2001
   rho_1(i,j)=0d0
   rho_1p5(i,j)=0d0
 enddo
enddo

w=2.2745d0*dsqrt(2d0)
sigma=w/(2d0*dsqrt(2d0*dlog(2d0)))

b=1/(2d0*(sigma**2)) 

do i=1,2001
 do j=1,2001
  x=i*10d0/2000d0
  y=j*10d0/2000d0
  rho_1(i,j)=dexp(-b*((x-5.005)**2+(y-5.005)**2))
 enddo
enddo

do i=1,2001
 near_field_x(i)=-5d0+(i-1)*.005d0
enddo


fwhm1p67=215.44d0
beta=fwhm1p67/(2d0*dsqrt(dlog(2d0)))

do i=1,2001
  do j=1,2001
   rho_1p5(i,j)=rho_1(i,j)*exp(-((i-j)/beta)**2)
  enddo
enddo


return
end

!-------------------------------------------------------------------------------
subroutine hasslebachdecoherence(trajectory,reduced_diagonal,rho_1p5)
implicit none

real*8 trajectory(500),reduced_diagonal(2001),rho_1p5(2001,2001)
common/par/ me,q,k,v0,theta,length
real*8 me,q,k,v0,theta,length
real*8 delt,b,a
integer i,j,m
real*8 Pi,hbar,h,kbolt,Tkelvin
real*8 resistivity,resistivitygold,resistivitysilicon
real*8 inversehasstime,hassexp,numerator

delt=length/v0/500d0 	!assuming 500 steps (see propogation)

a=(10.9d0-5.8d0)*(1d-6)/(500d0)


do i=1,2001
 reduced_diagonal(i)=0d0
enddo

do m=1,1
 do i=890,1112 !it is unecessary to calculate the change at the tails, suff small
   j=2002-i
	b=a*((i*(5d-6)/1000d0)-(j*(5d-6)/1000d0))**2
	hassexp=dexp(-b/(trajectory(m)**3))
	reduced_diagonal(i)=rho_1p5(i,j)*hassexp
 enddo
enddo

do m=2,500
 do i=890,1112 !it is unecessary to calculate the change at the tails, suff small
   j=2002-i
        b=a*((i*(5d-6)/1000d0)-(j*(5d-6)/1000d0))**2
        hassexp=dexp(-b/(trajectory(m)**3))
	reduced_diagonal(i)=reduced_diagonal(i)*hassexp
 enddo
enddo


return
end
!----------------------------------------------------------------------------------

subroutine bininitialization(Nbin,binlower,binupper,binvalue,binposition,accum_reduced_diagonal)
implicit none

integer Nbin,i,j
real*8 accum_reduced_diagonal(Nbin,2001),binvalue(Nbin)
real*8 binlower(Nbin),binupper(Nbin),binposition(Nbin)
real*8 lowerbound,binsize

lowerbound=-6d-4
binsize=(6d-4)/1000d0

do i=1,Nbin
 binlower(i)=lowerbound+binsize*(i-1)
 binupper(i)=lowerbound+binsize*i
 binposition(i)=(binlower(i)+binupper(i))/2.0d0
 binvalue(i)=0d0
enddo

do i=1,Nbin
 do j=1,2001
  accum_reduced_diagonal(i,j)=0d0
 enddo
enddo

return
end

!--------------------------------------------------------------------------------

subroutine binning(finalheight,reduced_diagonal,Nbin,binlower,binupper,binvalue,accum_reduced_diagonal,s1,s2)
implicit none

integer Nbin,i,j,s1,s2
real*8 accum_reduced_diagonal(Nbin,2001)
real*8 binlower(Nbin),binupper(Nbin)
real*8 lowerbound,binsize,binvalue(Nbin)
real*8 reduced_diagonal(2001),finalheight,norm
norm=0d0

do i=1,2001
  norm=norm+reduced_diagonal(i)
enddo

do i=1,2001
  reduced_diagonal(i)=reduced_diagonal(i)/norm
enddo

do i=1,Nbin
 if (binlower(i).le.finalheight) then
  if (finalheight.le.binupper(i)) then
   binvalue(i)=binvalue(i)+1d0
!	if (binlower(i).ge.-1.82d-5) then
!	 if (binupper(i).le.1.82d-5) then
!	   write(16,*) s1,' ',s2
!	 endif
!	endif
   do j=1,2001
    accum_reduced_diagonal(i,j)=accum_reduced_diagonal(i,j)+reduced_diagonal(j)
   enddo
  endif
 endif
enddo

return
end

!--------------------------------------------------------------------------------

subroutine writeelectrondistribution(binvalue,binposition,Nbin)
implicit none

integer Nbin,i
real*8 binposition(Nbin),binvalue(Nbin)

do i=1,Nbin
 write(14,114) binposition(i), binvalue(i)
enddo

114     format(100E16.6E4)

return
end

!-----------------------------------------------------------------

subroutine writeoffdiagonal(binposition,accum_reduced_diagonal,Nbin)
implicit none

integer Nbin,i,j
real*8 binposition(Nbin),accum_reduced_diagonal(Nbin,2001)
 
do i=1,Nbin
 do j=890,1112
  write(15,115) binposition(i), accum_reduced_diagonal(i,j)
 enddo
enddo

115	format(100E16.6E4)

return
end

!-------------------------------------------------------------

subroutine gratinginitialization(gratingarray)
implicit none

real*8 gratingarray(32),offset
integer i
offset=(-11d-6)+.5d-6


do i=1,31,2
gratingarray(i)=(.75d-6)*(i-1)+offset
enddo

do i=2,32,2
gratingarray(i)=(.75d-6)*(i-2)+.5d-6+offset
enddo

do i=1,32
write(6,*) gratingarray(i)
enddo

return
end

!-------------------------------------------------------------


subroutine write_rho(rho_1p5,near_field_x)
implicit none

real*8 rho_1p5(2001,2001),near_field_x(2001)
integer i,j

do i=1,2001
 j=2002-i
 write(17,117) near_field_x(i),rho_1p5(i,j)
enddo


117      format(100E16.6E4)

return
end
