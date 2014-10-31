c
c     freeze out
c
      implicit none
      integer imax,tmax,i,j,k,l,ndim,npart,ii,jmax,m
      parameter (ndim=2,npart=4492,imax=100,tmax=200,jmax=10)
      real z,z0,z1,x1
      double precision xm,Tfo,Tc,pmax,dp,dtheta,ans,
     &                 ans1,ans2,ratio,thet,px,py,x2,i1,i2,i1d,
     &                 u0,u1,u2,dp0,unp,press,entalp,dwds,c2,p,e,dndp,
     &                 nor,np,nu,f
      double precision un0(npart)
      double precision un(ndim,npart),kotae(imax,tmax)
      double precision temp(npart)
      double precision t(npart),entrop(npart),theta(npart),
     &                 sphbulk(npart),sigma(npart),s(npart),ref(npart)
      double precision r(ndim,npart),v(ndim,npart),
     &                 gradP(ndim,npart),f1(ndim,npart)
      double precision du(ndim,ndim,npart),sphshear(ndim,ndim,npart)
      double precision t1(npart),entrop1(npart),theta1(npart),
     &                 sphbulk1(npart),sigma1(npart),s1(npart),
     &                 ref1(npart)
      double precision r1(ndim,npart),v1(ndim,npart),
     &                 gradP1(ndim,npart),f11(ndim,npart)
      double precision du1(ndim,ndim,npart),sphshear1(ndim,ndim,npart)
      double precision gm(3,3),gmd(3,3)
      double precision gammaloren(npart),gammaloren1(npart)
      double precision weight1(npart),weight2(npart)
      double precision pi,hc
      real bessk

      common /unit/ hc
      common /num/ pi

      open(unit=11,file='hyper_130.dat')
c      open(unit=11,file='hyper_150.dat')
c      open(unit=11,file='hyper_170.dat')


      open(unit=13,file='cut_100.dat')
      

      open(unit=12,file='FONO_130_Tc=100_ver2.dat')
c      open(unit=12,file='FO_150.dat')
c      open(unit=12,file='FO_170_withoutstep.dat')


      open(unit=14,file='test.dat')



      call params

      xm = 139.0d0/hc
      Tfo = 130.0d0/hc
      Tc = 100.0d0/hc

      pmax = 3000.0d0/hc
      ymax = 10.0d0
      ymin = -ymax


      dp = pmax/dfloat(imax) 
      dy = (ymax-ymin)/dfloat(jmax)
      dtheta = 2.0d0*pi/dfloat(tmax)

      call condini

      do i=1,npart
      weight1(i) = 0.0d0
      weight2(i) = 0.0d0
      end do

c
c     v:index up, others : index down 
c


10    format (25(1x,E12.4))


c
c     calculation of weight
c



      do i=1,npart
      


      read(unit=11,fmt=10)  t(i),r(1,i),r(2,i),r(3,i),v(1,i),v(2,i),
     &                v(3,i),entrop(i),gradP(1,i),gradP(2,i),gradP(3,i),
     &                du(1,1,i),du(1,2,i),du(2,1,i),du(2,2,i),
     &                f1(1,i),f1(2,i),theta(i),
     &                sphshear(1,1,i),sphshear(1,2,i),sphshear(2,1,i),
     &                sphshear(2,2,i),sphbulk(i),sigma(i),s(i),ref(i)

      read(unit=13,fmt=10)  t1(i),r1(1,i),r1(2,i),r1(3,i),v1(1,i),v1(2,i),
     &                v1(3,i),entrop1(i),gradP1(1,i),gradP1(2,i),gradP1(3,i),
     &                du1(1,1,i),du1(1,2,i),du1(2,1,i),du1(2,2,i),
     &                f11(1,i),f11(2,i),theta1(i),
     &               sphshear1(1,1,i),sphshear1(1,2,i),sphshear1(2,1,i),
     &               sphshear1(2,2,i),sphbulk1(i),sigma1(i),
     &               s1(i),ref1(i)


      do k=1,3
	   do l=1,3
	     gm(k,l) = 0.0d0
	   enddo
      enddo

	gm(1,1) = -1.0d0
	gm(2,2) = -1.0d0
	gm(3,3) = -1.0d0/t(i)**2

      do k=1,3
	   do l=1,3
	     gmd(k,l) = 0.0d0
	   enddo
      enddo

      gmd(1,1) = -1.0d0
      gmd(2,2) = -1.0d0
      gmd(3,3) = -t(i)**2


      gammaloren(i) = 1.0d0
      gammaloren1(i) = 1.0d0
      do j=1,ndim
        gammaloren(i) = gammaloren(i) + v(j,i)*v(j,i)*gmd(j,j)
        gammaloren1(i) = gammaloren1(i) + v1(j,i)*v1(j,i)*gmd(j,j)
      end do
      gammaloren(i) = 1.0d0/dsqrt(gammaloren(i))
      gammaloren1(i) = 1.0d0/dsqrt(gammaloren1(i))

c
c     calculation of weight
c

 
      weight1(i) = weight1(i)*sigma(i)/ref(i)*t(i)
      

      z1 = xm/Tc

      weight2(i) = ref1(i)/sigma1(i)*weight1(i)
     &    /(bessk(2,z1)*xm**2*Tc/(2.0d0*pi**2))/t1(i)


        dp0 = 0.0d0
      do k=1,ndim
        dp0 =  dp0 - v(k,i)*gradP(k,i)
      end do

      call eqstate(entrop(i),press,temp(i),entalp,dwds)
      c2 = dwds/Tc - 1.0d0


      dp0 = dp0 - (entalp + sphbulk(i))*(theta(i)+gammaloren(i)/t(i))
     &      *c2/gammaloren(i)

      un0(i) = - dp0

      do j=1,ndim
        un(j,i) = - gradP(j,i)
      end do


      enddo
      
      
c
c     after weight
c
      
      do k=1,imax
          p=dp*dfloat(k)
          write(6,*) k
        do j=1,tmax
          thet = dtheta*dfloat(j-1)
          
          do m=1,jmax
            y = -ymax + dy*dfloat(m)


          px=p*dcos(thet)
          py=p*dsin(thet)
          E = dsqrt(xm**2 + px**2 + py**2)*dcosh(y)
          pl = dsqrt(xm**2 + px**2 + py**2)*dsinh(y)

      
      kotae(k,j,m) = 0.0d0


          do i=1, npart


      u0 = gammaloren(i)
      u1 = v(1,i)*gammaloren(i)
      u2 = v(2,i)*gammaloren(i)
      u3 = v(3,i)*gammaloren(i)


c      nor = un0(i)**2 - un(1,i)**2 - un(2,i)**2 +gmd(3,3)*un(3,i)**2

c      if(nor.le.0.0d0)then
      
c      kotae(k,j) = kotae(k,j) + sigma1(i)/ref1(i)*t1(i)
c     &       *weight2(i)*(E*i1d)/(2.0d0*pi)**3
     
c      else

      np gammaloren(i)*E - u1*px - u2*py - gmd(3,3)*u3*pl
      nu = un0(i)*u0 + un(1,i)*u1 + un(2,i)*u2 + un(3,i)*u3

      f = ( un0(i)*E*dcosh(y-x(3,i)) + (px*un(1,i)+py*un(2,i)) 
     &   + E/t(i)*un(3,i)*dsinh(y-x(3,i)))
     &      /(2.0d0*pi)**3
     &   *dexp(-ne/temp(i))/(1.0d0 - dexp(-ne/temp(i)))

      kotae(k,j,m) = kotae(k,j,m)
     & + sigma(i)*t(i)*gammaloren(i)/(ref(i)*nu)*f
      
c      end if

          end do

          end do
          

       end do

      
      end do


      do m=1,jmax
          y = -ymax + dy*dfloat(m)

      do k=1,imax
        p = dp*dfloat(k)
        dndp = 0.0d0

      do j=1,tmax

         dndp = dndp + kotae(k,j,m)*dtheta

      end do

      write(12,10) y, p*hc/1000.0d0,
     &       dndp/hc/hc*1000000.0d0/2.0d0/pi

      end do
      end do



      stop
      end 


c==========================================
      subroutine params
c==========================================
      implicit none
      integer ndim,npart
      doubleprecision pi,hc
      common /unit/ hc
      common /num/ pi

	pi = 4.0d0*atan(1.0d0)
        hc = 197.3d0

	return

	end



c==========================================
      subroutine condini
c========================================== 
      implicit none
      integer ndim,npart,imax
      integer i,j,k
      double precision table_s(50000),table_e(50000),table_P(50000),
     &                 table_cs2(50000), table_T(50000)
      common /table/ table_s,table_e,table_P,
     &               table_cs2,table_T


c
c     eos
c

      open(unit=35,file='EoS_pasi_uni.dat')

	do i=33001,1,-1
	read(35,*) table_T(i),table_cs2(i),table_e(i),
     &           table_P(i),table_s(i)
	end do

	do i=1,33001
	table_T(i) = table_T(i)*1000.0d0
	table_e(i) = table_e(i)*1000.0d0
        table_P(i) = table_P(i)*1000.0d0
	end do

	return

	end

!     ================================================================
      SUBROUTINE eqstate(s,p,t,entalp,dwds)
!     Equation of state
!     Input:  entrop = entropy density in local rest frame.
!     Output: pressure, temperature and enthalpy density in LRF.
!     ================================================================ 
c     Calculate P,e and T as a function of s 
c     Input:  s
c     output:  P, e, T
      implicit none
      integer i,j
      double precision ds, dT, cs2 ,logs, entalp, dwds, e, p, t, s
      double precision table_s(50000),table_e(50000),table_P(50000),
     &                 table_cs2(50000),table_T(50000)
      common /table/ table_s,table_e,table_P,
     &               table_cs2,table_T
 
      ds=0.0005d0

      logs = dlog(s)
   
      i = int(( logs - table_s(1) )/ds + 1)    !must be integer


      T = ( table_T(i+1) - table_T(i))/ds*(logs - table_s(i)) 
     &	+ table_T(i)
      P = ( table_P(i+1) - table_P(i))/ds*(logs - table_s(i)) 
     &	+ table_P(i)
      e = ( table_e(i+1) - table_e(i))/ds*(logs - table_s(i)) 
     &	+ table_e(i)
      cs2 = ( table_cs2(i+1) - table_cs2(i))/ds*(logs - table_s(i)) 
     &	+ table_cs2(i)



      dwds = (cs2+1.0d0)*T

      entalp = e + P

      T=T/197.3d0
      P=P/197.3d0
      e=e/197.3d0
      dwds=dwds/197.3d0
      entalp=entalp/197.3d0


	return

	 
      end

c
      subroutine sumi1(x1,x2,i1)
c
c
      implicit none
      integer ii
      double precision i1,ans1,ans2,ratio,x2,x2d
      real x1,x1d
      real bessk1
      external bessk1


      ans1 = 2.0d0*dexp(x2)*bessk1(x1)


      i1 = ans1

      ii=1

 102  continue


      x1d= x1*dfloat(ii+1)
      x2d = x2*dfloat(ii+1)

      ans2 = 2.0d0*dexp(x2d)*bessk1(x1d)

      ratio = dabs(ans2/ans1)

      if(ratio.ge.0.0001d0)then
      i1 = i1 + ans2
      ii = ii +1
      ans1 = ans2
      go to 102
      else
      return
      end if


      end

c
      subroutine sumi2(x1,x2,i2)
c
c
      implicit none
      integer ii
      double precision ans1,ans2,ratio,i2,x2,x2d
      real x1,x1d
      real bessk0
      external bessk0


      ans1 = 2.0d0*dexp(x2)*bessk0(x1)

      i2 = ans1

      ii=1

 103  continue


      x1d= x1*dfloat(ii+1)
      x2d = x2*dfloat(ii+1)

      ans2 = 2.0d0*dexp(x2d)*bessk0(x1d)

      ratio = dabs(ans2/ans1)

      if(ratio.ge.0.0001d0)then
      i2 = i2 + ans2
      ii = ii +1
      ans1 = ans2
      go to 103
      else
      return
      end if


      end

c
c
c
      FUNCTION bessk(n,x)  
      INTEGER n
      REAL bessk,x
CU    USES bessk0,bessk1
      INTEGER j
      REAL bk,bkm,bkp,tox,bessk0,bessk1

	if (n.lt.2) pause 'bad argument n in bessk'
	      tox=2.0/x  
	      bkm=bessk0(x)  
	      bk=bessk1(x)  
	do 11 j=1,n-1
	        bkp=bkm+j*tox*bk  
	        bkm=bk  
	        bk=bkp  
11    continue
	      bessk=bk  
	return
	END
c
c
c
      FUNCTION bessk0(x)  
      REAL bessk0,x
CU    USES bessi0
      REAL bessi0
      DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,y
      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7
      DATA p1,p2,p3,p4,p5,p6,p7/-0.57721566d0,0.42278420d0,0.23069756d0,
     *0.3488590d-1,0.262698d-2,0.10750d-3,0.74d-5/
      DATA q1,q2,q3,q4,q5,q6,q7/1.25331414d0,-0.7832358d-1,0.2189568d-1,
     *-0.1062446d-1,0.587872d-2,-0.251540d-2,0.53208d-3/
	if (x.le.2.0) then
	y=x*x/4.0
      bessk0=(-log(x/2.0)*bessi0(x))+(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*
     *(p6+y*p7))))))
	else
	y=(2.0/x)
      bessk0=(exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*
     *q7))))))
	endif
	return
	END
c
c
c
      FUNCTION bessk1(x)  
      REAL bessk1,x
CU    USES bessi1
      REAL bessi1
      DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,y
      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7
      DATA p1,p2,p3,p4,p5,p6,p7/1.0d0,0.15443144d0,-0.67278579d0,
     *-0.18156897d0,-0.1919402d-1,-0.110404d-2,-0.4686d-4/
      DATA q1,q2,q3,q4,q5,q6,q7/1.25331414d0,0.23498619d0,-0.3655620d-1,
     *0.1504268d-1,-0.780353d-2,0.325614d-2,-0.68245d-3/
	if (x.le.2.0) then
	y=x*x/4.0
	bessk1=(log(x/2.0)*bessi1(x))+(1.0/x)*(p1+y*(p2+y*(p3+y*(p4+y*
     *(p5+y*(p6+y*p7))))))
	else
        y=2.0/x
	bessk1=(exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*
     *q7))))))
	endif
	return
	END
c
c
c
      FUNCTION bessi0(x)
      REAL bessi0,x
      REAL ax
      DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y
      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
      DATA p1,p2,p3,p4,p5,p6,p7/1.0d0,3.5156229d0,3.0899424d0,
     *1.2067492d0,0.2659732d0,0.360768d-1,0.45813d-2/
      DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,0.1328592d-1,
     *0.225319d-2,-0.157565d-2,0.916281d-2,-0.2057706d-1,0.2635537d-1,
     *-0.1647633d-1,0.392377d-2/
	      if (abs(x).lt.3.75) then
	        y=(x/3.75)**2
	        bessi0=p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))
	      else
	        ax=abs(x)
	        y=3.75/ax
      bessi0=(exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*
     *(q7+y*(q8+y*q9))))))))
	      endif
	      return
	      END
c
c
c
      FUNCTION bessi1(x)
      REAL bessi1,x
      REAL ax
      DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y
      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
      DATA p1,p2,p3,p4,p5,p6,p7/0.5d0,0.87890594d0,0.51498869d0,
     *0.15084934d0,0.2658733d-1,0.301532d-2,0.32411d-3/
      DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,-0.3988024d-1,
     *-0.362018d-2,0.163801d-2,-0.1031555d-1,0.2282967d-1,-0.2895312d-1,
     *0.1787654d-1,-0.420059d-2/
	      if (abs(x).lt.3.75) then
	        y=(x/3.75)**2
	        bessi1=x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
	      else
	        ax=abs(x)
	        y=3.75/ax
      bessi1=(exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*
     *(q7+y*(q8+y*q9))))))))
	        if(x.lt.0.)bessi1=-bessi1
	      endif
	      return
	      END

