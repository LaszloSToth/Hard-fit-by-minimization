c     
c
c  fit of "Toth" strain hardening function 
c                               Laszlo S. Toth June 18, Metz, France
c input data: a given number (N) of stress - strain values from the experimental curve
c solution by: minimum value of the sum of the absolute values of the N differences between the experiment and teh fitted function
c output data: the simulated curve together with the input data, and the Kocks-Mecking strain hardening rate curve 
   																											 
                
      implicit real*8 (a-h,o-z)    
	dimension u(3),P(10000,2)

	open(2,file='hard.ctl',status='old')
	read(2,*)
	read(2,*) sig0
	read(2,*) u1min,u1max    ! initial interval for h      
	read(2,*) u2min,u2max    ! initial interval for tsat      
	read(2,*) u3min,u3max    ! initial interval for a
	read(2,*) ndiv           ! number of intervals between max and min      
	read(2,*) tol            ! tolerance in precision of minimization, 
	read(2,*) nmax1          ! maximum number of iterations in minimization 
	read(2,*) delta,epsmax   ! increment in strain for output data file, and maximum strain
	 
	close (2)
	
      open(4,file='result.dat',status='unknown')
      open(3,file='running.dat',status='unknown')
	open(8,file='exp-curve.dat',status='old')

	read(8,*) ndata
	do i=1,ndata
	read(8,*) P(i,1),P(i,2)		! x and y coordinates of input data
	end do

	write(3,*) ' iteration    h           sig-sat         a          fmin    average deviation'
	write(*,*) ' iteration    h           sig-sat         a          fmin    average deviation'
	write(4,*) '      strain stress-exp. stress-fit rate of hardening'
      u2min=1.00001*sig0
	fmin=1.d+20
      nfit=0
10    continue
      nfit=nfit+1
   	if(nfit.gt.nmax1) goto 20
	write(*,*)
 
      step1=(u1max-u1min)/ndiv
      step2=(u2max-u2min)/ndiv
      step3=(u3max-u3min)/ndiv

	call fit(sig0,P,h,taus,aa,u1min,u1max,u2min,u2max,u3min,u3max,fmin,step1,step2,step3,ndata)

	u1min=h-step1
	u1max=h+step1
	u2min=taus-step2
	u2max=taus+step2
	u3min=aa-step3
	u3max=aa+step3

	u(1)=h
	u(2)=taus
	u(3)=aa
	write(3,'(i4,7x,5e12.4)') nfit,u(1),u(2),u(3),fmin,(fmin/ndata)

	write(*,'(i4,7x,5e12.4)') nfit,u(1),u(2),u(3),fmin,(fmin/ndata)

	if(fmin.gt.tol) goto 10

20	continue
      write(*,*)
      if(nfit.gt.nmax1) write(*,*) ' maximum number of iterations reached in minimization'
	write(3,*)
      if(nfit.gt.nmax1) write(3,*) ' maximum number of iterations reached in minimization'

              
c write result into output file 'result.dat':
	do i=1,ndata
	eps=P(i,1)
      PP=-eps*u(1)*(1.-u(3))/(u(2)**u(3))+(u(2)-sig0)**(1.-u(3))
	stress=u(2)-PP**(1./(1.-u(3)))
	rate=u(1)*(1.-stress/u(2))**u(3)
	stress0=P(i,2)
	write(4,'(4f12.4)') eps,stress0,stress,rate 	 ! calculated stress
	end do	
      close (4)
	close (3)

c write result into output file for plotting simulation result only:
      open(4,file='hardplot.dat',status='unknown')
	write(4,*) '  strain      stress-fit rate of hardening'
	do eps=0,epsmax,delta
	PP=-eps*u(1)*(1.-u(3))/(u(2)**u(3))+(u(2)-sig0)**(1.-u(3))
	stress=u(2)-PP**(1./(1.-u(3)))			 ! simulated curve
	rate=u(1)*(1.-stress/u(2))**u(3)
	write(4,'(3f12.4)') eps,stress,rate 	 ! Kocks-Mecking plot
	end do
	close (4)


      stop
	end

      subroutine fit(sig0,P,h,taus,a,u1min,u1max,u2min,u2max,u3min,u3max,fmin,step1,step2,step3,ndata)
	implicit real*8 (a-h,o-z)
	dimension P(10000,2),f(10000)

	fmin=1.d+10
	do 300 u1=u1min,u1max,step1
	do 300 u2=u2min,u2max,step2
	do 300 u3=u3min,u3max,step3

      f0=0.
      do 200 i=1,ndata
	if(u2-sig0.lt.0.d0) goto 300
      PP=P(i,1)*u1*(u3-1.d0)/(u2**u3)+(u2-sig0)**(1.d0-u3)
	if(PP.lt.0.) goto 300
	f(i)=P(i,2)-u2+PP**(1.d0/(1.d0-u3))
	f0=f0+dabs(f(i))
200   continue

	if(f0.lt.fmin) then
	fmin=f0
	h=u1
	taus=u2
	a=u3
	endif

300   continue
                 
	return
	end

