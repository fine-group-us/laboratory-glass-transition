      program dsmc
c     Integrate the langevin equation equivalent to the FP equation of a molecular gas with nonlineal drag
c     The gas relaxes because the background temperature decresases as T(t)=T_0-r t 

      implicit none  
    
      integer id,d,iseed
      integer N,i,j,l,jmax
     
      real ran3
    
      real(8) vmod,varg,pi,rN,x,y,chi2,z0,zeff
      real*8 phi,DT,r,ctrlT,tmf,y1,y2,y3,sumv
      real*8 time,a2,a3,cj,c,vth,deltac,tf
      real*8 tboltz2d,val,deltalogT,ctrlogT,logT
      real*8 v2kk,v2,v4,temp,temp0,g,v6,tempf
      real*8 Deltat,tau,tmed
     
      parameter (g=0.1d0,r=0.005d0,DT=0.001d0,deltat=1.d-5)
      parameter (N=1*10**5,d=3,temp0=1.d0,tempf=1.d-6)
      parameter (jmax=200,deltac=0.01d0,tf=Temp0/r-2.d0*deltat)

      double precision v(N,d),vcm(d) 
      double precision fdv(jmax),cmoments(10)
       
c      open(1,file='kk2.dat',status="unknown")
      open(1,file='cool01g_0005r_1e5N.dat',status="unknown")
      open(2,file='fdv01g_0005r_1e5N.dat',status="unknown")
      open(3,file='Moments01g_0005r_1e5N.dat',status="unknown")
      open(4,file='Data01g_0005r_1e5N.dat',status="unknown")

      pi=4.d0*datan(1.d0)
      rN=dble(N)        
      iseed=-123456

!     Initialize some parameters 
      time=0.d0
      ctrlT=temp0-DT
      deltalogT=(dlog(Temp0)-dlog(tempf))/100.d0
      ctrlogT=dlog(temp0)-deltalogT
                        
!     initially: gaussian velocity distribution with T=temp0 
      do i=1,N
13      continue        
        x=ran3(iseed)
        if (x.eq.0) goto 13         
        y=ran3(iseed)
        vmod=dsqrt(-2.d0*temp0*dlog(x))
        varg=2.d0*pi*y
        v(i,1)=vmod*dcos(varg)
        v(i,2)=vmod*dsin(varg)
14      x=ran3(iseed)
        if (x.eq.0) goto 14         
        y=ran3(iseed)
        vmod=dsqrt(-2.d0*temp0*dlog(x))
        varg=2.d0*pi*y
        v(i,3)=vmod*dcos(varg)
      enddo   
      
!     Velocity of the Center of mass and change to the CM frame     
      do id=1,d
        vcm(id)=0.d0
      enddo    
      
      do id=1,d
        do i=1,N
            vcm(id)=vcm(id)+v(i,id)
        enddo
        vcm(id)=vcm(id)/rN
      enddo
           
      do i=1,N
        do id=1,d
            v(i,id)=v(i,id)-vcm(id)
        enddo
      enddo

c     Compute Temperature of the molecular fluid
      v2=0.d0
      v4=0.d0
      v6=0.d0
      do i=1,N
            v2kk=0.d0
            do id=1,d
              v(i,id)=v(i,id)-vcm(id)
              v2kk=v2kk+v(i,id)**2.d0            
            enddo
            v2=v2+v2kk
            v4=v4+v2kk**2
            v6=v6+v2kk**3
      enddo        
      v2=v2/rN
      tmf=v2/d
      v4=v4/rN
      v6=v6/rN/(2.d0*tmf)**3
      a2=v4/tmf**2.d0/dble(d*(d+2))-1.d0
      a3=1.d0+3.d0*a2-8.d0/d*v6/(d+2.d0)/(d+4.d0)
      logT=dlog(temp0)
      write(1,199) time,temp0,logT,tmf,a2,a3

!----------------------------------------------------      
!     Cooling process
!----------------------------------------------------
      do while (time.lt.tf) 
       time=time+deltat
       temp=temp0-r*time
       logT=dlog(temp)
                     
!     Apply the nonlinear drag and the thermostat
         do i=1,N
           v2=v(i,1)**2+v(i,2)**2+v(i,3)**2
           sumv=v(i,1)+v(i,2)+v(i,3)
           z0=2.d0*1.772453851d0*dsqrt(temp)
           zeff=z0*(1.d0-2.d0*g+g*v2/temp)
11         continue
           v(i,1)=v(i,1)*(1.d0-zeff*deltat) 
           v(i,2)=v(i,2)*(1.d0-zeff*deltat)
           v(i,3)=v(i,3)*(1.d0-zeff*deltat)
           chi2=2*z0*temp*(1.d0+g*v2/temp)
           val=dsqrt(chi2*deltat)
           x=ran3(iseed)
           if (x.eq.0.d0) goto 11
           y=ran3(iseed)
           y1=dsqrt(-2.d0*dlog(x))*dcos(2.d0*pi*y)
           y2=dsqrt(-2.d0*dlog(x))*dsin(2.d0*pi*y)
           v(i,1)=v(i,1)+val*y1+z0*g*sumv*(y1**2-1.d0)*deltat
           v(i,2)=v(i,2)+val*y2+z0*g*sumv*(y2**2-1.d0)*deltat
12         x=ran3(iseed)
           if (x.eq.0) goto 12         
           y=ran3(iseed)
           y3=dsqrt(-2.d0*dlog(x))*dcos(2.d0*pi*y)
           v(i,3)=v(i,3)+val*y3+z0*g*sumv*(y3**2-1.d0)*deltat
         enddo

c     Compute fluid temperature every Deltalog
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
         v6=0.d0
         do i=1,N
            v2kk=0.d0
            do id=1,d
              v(i,id)=v(i,id)-vcm(id)
              v2kk=v2kk+v(i,id)**2.d0            
            enddo
            v2=v2+v2kk
            v4=v4+v2kk**2
            v6=v6+v2kk**3
         enddo        
         v2=v2/rN
         tmf=v2/d
         v4=v4/rN
         v6=v6/rN/(2.d0*tmf)**3
         a2=v4/tmf**2.d0/dble(d*(d+2))-1.d0
         a3=1.d0+3.d0*a2-8.d0/d*v6/(d+2.d0)/(d+4.d0)
         write(1,199) time,temp,logT,tmf,a2,a3
         ctrlT=ctrlT-DT
         ctrlogT=ctrlogT-deltalogT
        endif 

      enddo
C     End of  relaxation


c     Compute fdv and moments
      fdv=0.d0
      cmoments=0.d0
      vth=2.d0*v2/d
      do i=1,N
            v2kk=0.d0
            do id=1,d
              v2kk=v2kk+v(i,id)**2.d0            
            enddo
            c=dsqrt(v2kk/vth)
            j=int(c/deltac)+1
            if (j.le.jmax) fdv(j)=fdv(j)+1.d0

c           Moments
            do 22 l=1,10
            cmoments(l)=c**l+cmoments(l)
22          continue
        write(4,199) (v(i,j),j=1,d) 
      enddo        

      fdv=fdv/rN/deltac
      do j=1,jmax
        cj=deltac*(j-0.5d0)
      write(2,199)  cj, fdv(j)
      enddo

      cmoments=cmoments/rN
      do l=1,10
       write(3,*)  l,cmoments(l)
      enddo


      close(30)
199   format(10(e14.7,2x))
20    format(5(e14.7,2x),2x,i8)
      
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
        integer mbig,mseed,mz,ma,mj,mk,i,ii,k,inext,inextp,iff,idum
        double precision fac
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











    
    
