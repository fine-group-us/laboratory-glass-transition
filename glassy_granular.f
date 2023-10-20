      program glassy
c     This program simulates a cooling process for a granular gas with restitution coefficient alpha. It integrates the Boltzmann-Foikker-Planck equation
c     that governs the gas dynamics by using the Direct Simulation Monte Carlo method.
c     T0=initial temperature
c     tempf= final temperature
c     N=number of particles
c     rteo= adimensional cooling rate, while rc=dT/dt
c     stepkick=number of collisions after which the gas particles are kicked at random (thermostat) 
c     v(N,d)= velocity in d dimensions
c     vcm(d)=velocity of the mass center
c     fdv=velocity distribution funcion of c=v/v_th, being v_th=2 k_b T/m
c     Deltat=unit of time
c     z0=collision rate for hard spheres
c     a2s= stationary kurtosis
c     tf=final time
c     tau=adimensional time

      implicit none  
    
      integer i,j,id,N,d,iseed,ncol,stepkick
      integer nkick,jmax
            
      real ran3
    
      real*8 vmod,varg,pi,rN,x,y,alpha,q
      real*8 phi,scprod,scprodmax,prob,tf
      real*8 time,a2s,z0,chi,rc,rteo,cj,c,vth,deltac
      real*8 timeold,timeinterval,val
      real*8 v2kk,Ts,v2,v4,a2,temp,tempf,logT,ctrlogT
      real*8 Deltat,tau,taueq,T0,deltalogT
      
      parameter (alpha=0.1d0,T0=5.d0,tempf=1.d-6)
      parameter (N=10**6,d=3,rteo=0.001d0)
      parameter (stepkick=N/2000,taueq=6.d0)
      parameter (jmax=100,deltac=0.04d0)	

      double precision v(N,d),vcm(d),sigma(d),fdv(jmax)
      
      open(30,file='cooling_0001rc_01a_1e6N.dat',status="new")
      open(1,file='0001rc01a.dat',status="new")
c      open(20,file='data_001rc_09a_1e6N.dat',status="new")
      open(2,file='fdv_0001rc_01a_1e6N.dat',status="new")


!     Initialize some parameters
      pi=4.d0*datan(1.d0)
      rN=dble(N)          
      iseed=-923456
      scprodmax=23.d0
      Deltat=(4.d0-d)/pi/rN/scprodmax
      a2s=16.d0*(1.d0-alpha)*(1.d0-2.d0*alpha**2)
      q=73.d0+56.d0*d-24.d0*d*alpha-105.d0*alpha
      q=q+30.d0*(1.d0-alpha)*alpha**2
      a2s=a2s/q
      if (d.eq.2) then 
        z0=(1.d0-alpha**2)*dsqrt(pi)
      else 
        z0=4.d0*(1.d0-alpha**2)*dsqrt(pi)/3.d0 
      endif
      rc=rteo*T0**1.5*z0
      tf=(T0-tempf)/rc
       
      ncol=0
      nkick=0
      time=0.d0
      timeold=0.d0  
      deltalogT=(dlog(T0)-dlog(tempf))/100.d0
      ctrlogT=dlog(T0)-deltalogT

!---------------------------------------------------                               
!     Gaussian velocity distribution with temp=T0
!------------------------------------------------------------
      do i=1,N
13      continue        
        x=ran3(iseed)
        if (x.eq.0) goto 13         
        y=ran3(iseed)
        vmod=dsqrt(-2.d0*T0*dlog(x))
        varg=2.d0*pi*y
        v(i,1)=vmod*dcos(varg)
        v(i,2)=vmod*dsin(varg)
14      x=ran3(iseed)
        if (x.eq.0) goto 14         
        y=ran3(iseed)
        vmod=dsqrt(-2.d0*T0*dlog(x))
        varg=2.d0*pi*y
        v(i,3)=vmod*dcos(varg)
      enddo   
      
!----------------------------------------------------      
!     Prepare the initial stationary state 
!----------------------------------------------------  
      do 99 while (tau.lt.taueq)
        chi=z0*T0**1.5d0*(1.d0+3.d0/16.d0*a2s)
!     We choose a random direction sigma in d-dimensions     
        x=ran3(iseed)
        phi=2.d0*pi*x
        sigma(1)=dcos(phi)
        sigma(2)=dsin(phi)
        if (d.eq.3) then
        x=ran3(iseed)
        phi=2.d0*pi*x
        sigma(1)=sigma(1)*dsin(phi)
        sigma(2)=sigma(2)*dsin(phi)
        sigma(3)=dcos(phi)        
       endif     
         
!     We choose the pair to collide: there is always a collision, because either
!     (i,j) or (j,i) verifies that the scalar product of the relative velocity with
!     sigma is positive 
       i=1+int(N*ran3(iseed))
10     continue
       j=1+int(N*ran3(iseed))      
       if (j.eq.i) goto 10
       if (i.gt.N) i=N
       if (j.gt.N) j=N
       scprod=0.d0
       do id=1,d
         scprod=scprod+(v(i,id)-v(j,id))*sigma(id)
       enddo
                  
!     Time is increased either a collision is accepted or not
       time=time+deltat
     
 !    Accept the collision with a probability proportional to the scalar product
 !    of its relative velocity and change the velocities
       x=ran3(iseed)
       prob=dabs(scprod)/scprodmax
       if (dabs(scprod).gt.scprodmax) then 
           write(1,199) T0,scprod
           scprodmax=dabs(scprod)
           Deltat=(4.d0-d)/pi/rN/scprodmax
       endif
       if (x.lt.prob) then
            do id=1,d
                v(i,id)=v(i,id)-(1.d0+alpha)*scprod*sigma(id)/2.d0
                v(j,id)=v(j,id)+(1.d0+alpha)*scprod*sigma(id)/2.d0
            enddo 

!     updates the collision counters
         ncol=ncol+1 
	 nkick=1       
       endif
                  
!     Every given number of collisions, we kick all the particles at random
        if (mod(ncol,stepkick).eq.0 .and.(nkick.eq.1)) then
          timeinterval=time-timeold
          val=sqrt(chi*timeinterval)

          do i=1,N
15          continue
            x=ran3(iseed)
            if (x.eq.0.d0) goto 15
            y=ran3(iseed)
            v(i,1)=v(i,1)+val*dsqrt(-2.d0*dlog(x))*dcos(2.d0*pi*y)
            v(i,2)=v(i,2)+val*dsqrt(-2.d0*dlog(x))*dsin(2.d0*pi*y)
            if (d.eq.3) then
16            continue
              x=ran3(iseed)
              if (x.eq.0.d0) goto 16
              y=ran3(iseed)
              v(i,3)=v(i,3)+val*dsqrt(-2.d0*dlog(x))*dcos(2.d0*pi*y)
            endif
          enddo

	  nkick=0
          timeold=time
        endif
       tau=0.5d0*z0*dsqrt(T0)*time
99     continue 
     
c     Initial Temperature after eliminating a nonvanishing CM velocity
      do id=1,d
         vcm(id)=0.d0
         do i=1,N
            vcm(id)=vcm(id)+v(i,id)
         enddo
         vcm(id)=vcm(id)/rN
      enddo
      
      v2=0.d0
      v4=0.d0
      do i=1,N
            v2kk=0.d0
            do id=1,d
              v(i,id)=v(i,id)-vcm(id)
              v2kk=v2kk+v(i,id)**2.d0            
            enddo
            v2=v2+v2kk
            v4=v4+v2kk**2
      enddo        
      v2=v2/rN
      temp=v2/d
      v4=v4/rN
      a2=v4/temp**2.d0/dble(d*(d+2))-1.d0         
      time=0.d0
      write(30,199) time,T0,temp,a2

!----------------------------------------------------      
!     Cooling process
!----------------------------------------------------
      ncol=0
      nkick=0
      timeold=0.d0 
      do 100 while (time.lt.tf) 
       Ts=T0-rc*time
       logT=dlog(Ts)
       chi=z0*Ts**1.5d0*(1.d0+3.d0/16.d0*a2s)

!     We choose a random direction sigma in d-dimensions     
        x=ran3(iseed)
        phi=2.d0*pi*x
        sigma(1)=dcos(phi)
        sigma(2)=dsin(phi)
        if (d.eq.3) then
         x=ran3(iseed)
         phi=2.d0*pi*x
         sigma(1)=sigma(1)*dsin(phi)
         sigma(2)=sigma(2)*dsin(phi)
         sigma(3)=dcos(phi)        
        endif     
         
!     We choose the pair to collide: there is always a collision, because either
!     (i,j) or (j,i) verifies that the scalar product of the relative velocity with
!     sigma is positive 
       i=1+int(N*ran3(iseed))
1     continue
       j=1+int(N*ran3(iseed))      
       if (j.eq.i) goto 1
       if (i.gt.N) i=N
       if (j.gt.N) j=N
       scprod=0.d0
       do id=1,d
         scprod=scprod+(v(i,id)-v(j,id))*sigma(id)
       enddo
                  
!     Time is increased either a collision is accepted or not
       time=time+deltat
     
 !    Accept the collision with a probability proportional to the scalar product
 !    of its relative velocity and change the velocities
       x=ran3(iseed)
       prob=dabs(scprod)/scprodmax
       if (dabs(scprod).gt.scprodmax) then 
           write(1,*) Ts,scprod
           scprodmax=dabs(scprod)
           Deltat=(4.d0-d)/pi/rN/scprodmax
       endif
       if (x.lt.prob) then
            do id=1,d
                v(i,id)=v(i,id)-(1.d0+alpha)*scprod*sigma(id)/2.d0
                v(j,id)=v(j,id)+(1.d0+alpha)*scprod*sigma(id)/2.d0
            enddo 
!     updates the collision counters
         ncol=ncol+1 
	 nkick=1       
       endif
                  
!     Every given number of collisions, we kick all the particles at random
        if (mod(ncol,stepkick).eq.0 .and.(nkick.eq.1)) then
          timeinterval=time-timeold
          val=sqrt(chi*timeinterval)
          do i=1,N
5           continue
            x=ran3(iseed)
            if (x.eq.0.d0) goto 5
            y=ran3(iseed)
            v(i,1)=v(i,1)+val*dsqrt(-2.d0*dlog(x))*dcos(2.d0*pi*y)
            v(i,2)=v(i,2)+val*dsqrt(-2.d0*dlog(x))*dsin(2.d0*pi*y)
            if (d.eq.3) then
6             continue
              x=ran3(iseed)
              if (x.eq.0.d0) goto 6
              y=ran3(iseed)
              v(i,3)=v(i,3)+val*dsqrt(-2.d0*dlog(x))*dcos(2.d0*pi*y)
            endif
          enddo
	  nkick=0
          timeold=time
        endif
     
c     Compute granular temperature every DeltalogT
        if (logT.le.ctrlogT) then 
          do id=1,d
            vcm(id)=0.d0
            do i=1,N
            vcm(id)=vcm(id)+v(i,id)
            enddo
            vcm(id)=vcm(id)/rN
          enddo
      
          v2=0.d0
          v4=0.d0
          do i=1,N
            v2kk=0.d0
            do id=1,d
              v(i,id)=v(i,id)-vcm(id)
              v2kk=v2kk+v(i,id)**2.d0            
            enddo
            v2=v2+v2kk
            v4=v4+v2kk**2
          enddo        
          v2=v2/rN
          temp=v2/d
          v4=v4/rN
          a2=v4/temp**2.d0/dble(d*(d+2))-1.d0   
          write(30,199) time,Ts,temp,a2
          ctrlogT=ctrlogT-deltalogT
        endif
100   continue

c     Write the frozen state
c      do i=1,N
c        write(20,199) (v(i,j), j=1,d)
c      enddo

c     Compute fdv
      fdv=0.d0
      vth=2.d0*v2/d
      do i=1,N
            v2kk=0.d0
            do id=1,d
              v2kk=v2kk+v(i,id)**2.d0            
            enddo
            c=dsqrt(v2kk/vth)
            j=int(c/deltac)+1
            if (j.le.jmax) fdv(j)=fdv(j)+1.d0
      enddo        

      fdv=fdv/rN/deltac
      do j=1,jmax
        cj=deltac*(j-0.5d0)
      write(2,199)  cj, fdv(j)
      enddo
      
      close(30)
199   format(10(e14.7,2x))
      
      stop
      end
    
   
c====================================================================
!     Subrutina Ran3 numeros aleatorios	
C====================================================================
C       PROGRAM: ran3.f
C       TYPE   : function
C       PURPOSE: generate random numbers
C       VERSION: 17 June 94
C       COMMENT: Initialize idum with negative integer
C======================================================================
        real function ran3(idum)
        Parameter (mbig=1000000000,Mseed=161803398,Mz=0,fac=1./Mbig)
        Dimension MA(55)
        save
        if (idum.lt.0.or.iff.eq.0) then
	   iff=1
	   mj=mseed-iabs(idum)
	   mj=mod(mj,mbig)
	   ma(55)=mj
	   mk=1
	   do 11 i=1,54
	      ii=mod(21*i,55)
	      ma(ii)=mk
	      mk=mj-mk
	      if (mk.lt.mz) mk=mk+mbig
	      mj=ma(ii)
 11	   continue
	   do 13 k=1,4
	      do 12 i=1,55
		 ma(i)=ma(i)-ma(1+mod(i+30,55))
		 if (ma(i).lt.mz) ma(i)=ma(i)+mbig
 12	      continue
 13	   continue
	   inext=0
	   inextp=31
	   idum=1
	end if
	inext=inext+1
	if (inext.eq.56) inext=1
	inextp=inextp+1
	if (inextp.eq.56) inextp=1
	mj=ma(inext)-ma(inextp)
	if (mj.lt.mz) mj=mj+mbig
	ma(inext)=mj
	ran3=mj*fac
	return
	end











    
    
