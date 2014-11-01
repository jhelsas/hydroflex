c
c     freeze out
c
      implicit none
      integer i,j,npart,ndim,k,l
      integer kx,ky,kz,ix,iy,iz,ixk,iyk,izk
      integer imax,tmax
      parameter (ndim=2,npart=6844)
      double precision press,entalp,dwds,c2,dp0,nor,dndp
      double precision hw,w0,hc
      double precision entrop(npart),y(ndim),shear00(npart),un0(npart),
                     & sphbulk(npart),theta(npart),t(npart)
      double precision r(ndim,npart),v(ndim,npart),gradP(ndim,npart),
                     & f1(ndim,npart),gm(3,3),gmd(3,3),
                     & shear0(ndim,npart),ans(800,800)
      double precision sphshear(3,3,npart),du(ndim,ndim,npart)
      double precision un(ndim,npart)
      double precision pmax,dp,px,py,u0,u1,u2,np,nu,ne,f,e,xm,pi,
                     & dtheta,thet,p
      double precision sigma(npart),s(npart),ref(npart),
                     & gammaloren(npart),temp(npart)
      common /kernel/ hw,w0
      common /unit/ hc
      common /num/ pi

      open(unit=11,file='hyper.dat')
      open(unit=12,file='feezout_exclpult0.dat')
      open(unit=13,file='test_3.dat')
      open(unit=14,file='nu_distribution.dat')
    


      call params

      xm = 139.0d0/hc
      pmax = 3000.0d0/hc
      imax = 100
      dp = pmax/dfloat(imax) 
      tmax = 400
      dtheta = 2.0d0*pi/dfloat(tmax)

      call condini
c
c     v:index up, others : index down 
c

10    format (25(1x,E12.4))
      
      do i=1,npart
        read(unit=11,fmt=10)  t(i),r(1,i),r(2,i),v(1,i),v(2,i),
                            & entrop(i),gradP(1,i),gradP(2,i),
                            & du(1,1,i),du(1,2,i),du(2,1,i),du(2,2,i),
                            & f1(1,i),f1(2,i),theta(i),
                      & sphshear(1,1,i),sphshear(1,2,i),sphshear(2,1,i),
                       & sphshear(2,2,i),sphbulk(i),sigma(i),s(i),ref(i)
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
        do j=1,ndim
          gammaloren(i) = gammaloren(i) + v(j,i)*v(j,i)*gmd(j,j)
        end do
        gammaloren(i) = 1.0d0/dsqrt(gammaloren(i))

        write(14,*)  gammaloren(i), sigma(i)/t(i)

c
c     calculation of other components of sphshear
c

        shear00(i) = 0.0d0 
        do j=1,ndim
	        shear0(j,i) = 0.0d0
        end do

        do j=1,ndim
          do k=1,ndim
             shear00(i) = shear00(i)
           &	+ v(j,i)*v(k,i)*sphshear(j,k,i)
             shear0(j,i) = shear0(j,i) 
           &    - v(k,i)*sphshear(j,k,i)
          end do
        end do
c
c    calculation of time derivative of pressure
c

        dp0 = 0.0d0
        
        call eqstate(entrop(i),press,temp(i),entalp,dwds)
        c2 = dwds/temp(i) - 1.0d0
        
        do k=1,ndim
          dp0 =  dp0 - v(k,i)*gradP(k,i)
        end do

        dp0 = dp0 - (entalp + sphbulk(i))*theta(i)*c2/gammaloren(i) 
c
c     definition of unit vector
c
        nor = dp0**2 
        do k=1,ndim
          nor = nor + gradP(k,i)*gradP(k,i)*gm(k,k)
        end do
        un0(i) = - dp0
        
        do j=1,ndim
          un(j,i) = - gradP(j,i)
        end do
      
        write(13,*) i, un0(i), un(1,i), un(2,i), c2

      end do


c
c     calcualtion of Ed^3 N/dp^3
c


      do k=1,imax
        p=dp*dfloat(k)
        write(6,*) k 
        do j=1,tmax+1
          thet = dtheta*dfloat(j-1)
          px=p*dcos(thet)
          py=p*dsin(thet)
          ans(k,j) = 0.0d0
          do i = 1,npart
            u0 = gammaloren(i)
            u1 = v(1,i)*gammaloren(i)
            u2 = v(2,i)*gammaloren(i)
            E = dsqrt(xm**2 + px**2 + py**2)
            np = un0(i)*E + un(1,i)*px + un(2,i)*py
            nu = un0(i)*u0 + un(1,i)*u1 + un(2,i)*u2 
            ne = u0*E - u1*px -u2*py
            
            f = dexp(- ne/temp(i) )/(1.0d0 - dexp(- ne/temp(i) ) )
                 &      /(2.0d0*pi)**3
            
            if(nor.ge.0.0d0) then
              if(np.le.0.0d0)then
                ans(k,j) = ans(k,j)
              else
                ans(k,j) = ans(k,j)
                & + sigma(i)*np/( t(i)*ref(i)/u0*nu )*f
              end if
            else
              if(np.le.0.0d0) then
                ans(k,j) = ans(k,j)
              else
                ans(k,j) = ans(k,j)
                & + sigma(i)*np/( t(i)*ref(i)/u0)/dsqrt(-nor)*f
              end if
            endif
          end do
        end do
      end do

c
c     calculation of dN/pdp
c
      do k=1,imax
        p = dp*dfloat(k)
        e = p**2 + xm**2
        e = dsqrt(e)
        dndp = 0.0d0
        do j=1,tmax+1
         dndp = dndp + ans(k,j)*dtheta
        end do
  
        write(12,10) p*hc/1000.0d0,dndp/hc/hc*1000000.0d0
      end do

      stop
      end 


c========================================== Set common constants 
      subroutine params
c========================================== \hbar c,\pi,hw?(kernel), w0 
      implicit none
      integer ndim,npart
	
      doubleprecision hw,w0,pi,hc
	    
      common /kernel/ hw,w0
	    common /unit/ hc
      common /num/ pi

      pi = 4.0d0*atan(1.0d0)
      hc = 197.3d0
      hw = 0.33941d0

      w0 = 10.d0/7.d0/pi/hw**2

      return
      end



c========================================== Initializes EOS tables
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
                 & table_P(i),table_s(i)
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
        & +table_T(i)
      P = ( table_P(i+1) - table_P(i))/ds*(logs - table_s(i)) 
        & +table_P(i)
      e = ( table_e(i+1) - table_e(i))/ds*(logs - table_s(i)) 
        & + table_e(i)
      cs2 = ( table_cs2(i+1) - table_cs2(i))/ds*(logs - table_s(i)) 
        & + table_cs2(i)

      dwds = (cs2+1.0d0)*T
      entalp = e + P
      
      T=T/197.3d0
      P=P/197.3d0
      e=e/197.3d0
      dwds=dwds/197.3d0
      entalp=entalp/197.3d0
      
      return
      end


