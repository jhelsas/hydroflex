c     Relativistic Hydrodynamics SPH Program from Tomoi Koide
c     tau/eta coordinates, freezeout and simplified shear+bulk viscosity
c     freezeout uses quadratric interpolation
c     Equation of state is a lattice based, input file EoS
      implicit none
      integer ndim,npart,imax
      parameter (ndim=3,npart=10648,imax=100)
      integer ntstep,ndiv,i,j,ioutp,ii
      integer frag(npart),frag2(npart),frag3(npart)
      double precision t,ez,esum,esum0
      double precision x(ndim,npart),u(ndim,npart),v(ndim,npart)
      double precision gradP2(ndim,npart)
      double precision s(npart),sigma(npart),entrop(npart),
     &                 gammaloren(npart),ref(npart),rho(npart)
      double precision bulk(npart),theta(npart)
      double precision shear(3,3,npart)
      double precision du(3,3,npart)
      double precision f1(3,npart)
      double precision link(npart),rmin(ndim),rmax(ndim),dr(ndim)
      double precision lead(-3:imax,-3:imax,-3:imax)
c     double precision lead(-3:imax,-3:imax)

      character*7 outfile(60)
      data outfile/'t01.dat','t02.dat','t03.dat','t04.dat','t05.dat',
     &             't06.dat','t07.dat','t08.dat','t09.dat','t10.dat',
     &             't11.dat','t12.dat','t13.dat','t14.dat','t15.dat',
     &             't16.dat','t17.dat','t18.dat','t19.dat','t20.dat',
     &             't21.dat','t22.dat','t23.dat','t24.dat','t25.dat',
     &             't26.dat','t27.dat','t28.dat','t29.dat','t30.dat',
     &             't31.dat','t32.dat','t33.dat','t34.dat','t35.dat',
     &             't36.dat','t37.dat','t38.dat','t39.dat','t40.dat',
     &             't41.dat','t42.dat','t43.dat','t44.dat','t45.dat',
     &             't46.dat','t47.dat','t48.dat','t49.dat','t50.dat',
     &             't51.dat','t52.dat','t53.dat','t54.dat','t55.dat',
     &             't56.dat','t57.dat','t58.dat','t59.dat','t60.dat'/
      common /steps/ ntstep,ndiv 

      call params
      call condini(x,u,s,sigma,bulk,shear,t)

      ez = 0.0d0
      ioutp = 0
      do i=1,npart
        frag(i) = 0
        frag2(i) = 0
        frag3(i) = 0
      end do
      
      esum0=0.0d0
      ii = 0
      call hyper0(x,u,s,sigma,entrop,t,rmin,rmax,dr,link,lead,
     &            frag,frag2,frag3,ii)
      call out(x,u,entrop,sigma,ref,t,esum,esum0,
     &         ez,bulk,shear,rmin,rmax,dr,link,lead,'t00.dat')
             
      do i=1,ntstep
        write(6,*) i
        call step(x,u,s,sigma,ez,bulk,shear,t,gradP2,f1,du,theta)
        call hyper(x,u,s,sigma,entrop,t,rmin,rmax,dr,link,lead
     &             ,gradP2,f1,du,shear,bulk,theta,frag,frag2,frag3,ii)
        if(mod(i,ndiv).eq.0) then
          ioutp = ioutp + 1
          call out(x,u,entrop,sigma,ref,t,esum,esum0,
     &             ez,bulk,shear,rmin,rmax,dr,link,lead,outfile(ioutp))
          write(6,*) 't=',t
        end if
      end do
       
      stop
      end

c==========================================
      subroutine params
c==========================================
      implicit none
      integer ndim,npart
      parameter (ndim=3,npart=10648)
      integer ntstep,nout,ndiv,i,j
      double precision hwt,hwl,dt,w0,pi,hc,w1,w2
      double precision covisb,corelb
      common /kernel/ hwt,hwl,w0
      common /time/ dt
      common /steps/ ntstep,ndiv 
      common /unit/ hc
      common /bulk/ covisb,corelb

      pi = 4.0d0*atan(1.0d0)
      hc = 197.3d0
      ntstep = 1600
      dt = 0.2d0
      hwt = 0.4d0
      hwl = 0.4d0
      nout = 20
      ndiv = ntstep/nout
      w2 = 10.d0/7.d0/pi/(hwt**2) !w0 for 2D
      w1 = 2.0d0/3.0d0/hwl        !w0 for 1D
      w0 = w1*w2
      covisb = 0.1d0
      corelb = 6.0d0
      
      return
      end

c==========================================
      subroutine condini(x,u,s,sigma,bulk,shear,t)
c========================================== 
      implicit none
      integer ndim,npart,imax
      parameter (ndim=3,npart=10648,imax=100)
      integer i,j,k
      double precision x(ndim,npart),u(ndim,npart)
      double precision lead(-3:imax,-3:imax)
      double precision s(npart),sigma(npart),entrop(npart)
      double precision bulk(npart)
      double precision shear(3,3,npart)
      double precision t
      double precision table_s(50000),table_e(50000),table_P(50000),
     &                 table_cs2(50000) 
      double precision table_T(50000),table_T2(50000)
      common /table/ table_s,table_e,table_P,
     &               table_cs2,table_T

      open(unit=31,file='ini_koide_2.dat') ! Initial condition

      do i=1,npart
        read(31,*) x(1,i), x(2,i), x(3,i), s(i) ! sigma(i)
      end do
      
      do i=1,npart
        sigma(i) = 1.0d0 ! s(i)
        do j=1,3
          u(j,i) = 0.0d0
        end do
        
        bulk(i) = 0.0d0
        do j=1,3
          do k=1,3
            shear(j,k,i) = 0.0d0
          end do
        end do
      end do

      t = 0.6d0

      open(unit=35,file='EoS_pasi_uni.dat') ! eos
      
      do i=33001,1,-1
        read(35,*) table_T(i),table_cs2(i),table_e(i),
     &             table_P(i),table_s(i)
      end do
      
      do i=1,33001
        table_T(i) = table_T(i)*1000.0d0
        table_e(i) = table_e(i)*1000.0d0
        table_P(i) = table_P(i)*1000.0d0
      end do
      
c     open(unit=51,file='test1.dat')
c     do i= 1,33001
c	    write(51,*) table_T(i),table_cs2(i),table_e(i),
c               & table_P(i),table_s(i)
c	    end do
      
      return
      end

c     ================================================
      subroutine wkernel(xx,wk,dwk)
c
c     kernel and its derivative
c     ================================================      
      implicit none
      integer ndim,k
      parameter (ndim=3)
      double precision xx(ndim),dwk(ndim),dwkt(ndim),dwkl(ndim)
      double precision hwt,hwl,w0,wk,wkt,wkl
      common /kernel/ hwt,hwl,w0

      call wkernel_t(xx,wkt,dwkt)
      call wkernel_l(xx,wkl,dwkl)

      wk = w0*wkt*wkl
      do k=1,2
        dwk(k) = dwkt(k)*wkl*w0 
      end do 
      dwk(3) = dwkl(3)*wkt*w0

      return
      end

c     ===============================================================
      subroutine wkernel_t(xx,wk,dwk)
c
c	    Kernal in transverse direction
c     ===============================================================
      implicit none
      integer ndim,k
      parameter (ndim=3)
      double precision xx(ndim),dwk(ndim)
      double precision r,q,hwt,hwl,w0,qq,qq2,q2,dw,wk,hwl2,hwt2
      common /kernel/ hwt,hwl,w0

      hwl2 = hwl**2
      hwt2 = hwt**2
      do k=1,ndim
       dwk(k) = 0.0d0
      end do

c     r = 0.0d0
c	    do k=1,2
c       r = r + xx(k)*xx(k)
c	    enddo

      r = dsqrt(xx(1)*xx(1)+xx(2)*xx(2))
      q=r/hwt
      if (q.ge.2.d0) then
        wk = 0.d0
        do k=1,2
          dwk(k) = 0.d0
        end do
        return
      end if
      
      if (q.ge.1.d0) then
        qq = 2.d0 - q
        qq2 = qq*qq
        wk = qq2*qq * 0.25d0
        dw = - qq2 * 0.75d0 / hwt
        do k=1,2
          dwk(k) = xx(k)/r * dw
        end do
        return
      end if

      wk = (1.d0 - 1.5d0*q**2 + 0.75d0*q**3)
      dw = (-3.d0 + 2.25d0*q) / hwt2

      do k=1,2        
        dwk(k) = xx(k) * dw
      end do

      return 
      end


c     ================================================================
      subroutine wkernel_l(xx,wk,dwk)
c	kernal in langitudinal direction
c     ================================================================
      implicit none
      integer ndim,k
      parameter (ndim=3)
      double precision xx(ndim),dwk(ndim)
      double precision r,q,hwt,hwl,hwt2,hwl2,qq,qq2,q2,dw,w0,wk
      common /kernel/ hwt,hwl,w0

 
      hwl2 = hwl**2
      hwt2 = hwt**2       
      do k=1,ndim
        dwk(k) = 0.0d0
      end do

      r = dsqrt( xx(3)*xx(3) )

      q=r/hwl
       
      if (q.ge.2.d0) then
        wk = 0.d0
        dwk(3) = 0.d0
        return
      end if
      
      if (q.ge.1.d0) then
        qq = 2.d0 - q
        qq2 = qq*qq
        wk = qq2*qq * 0.25d0
        dw = - qq2 * 0.75d0 / hwl
        dwk(3) = xx(3)/r * dw
        return
      end if

      wk = (1.d0 - 1.5d0*q**2 + 0.75d0*q**3)
      dw = (-3.d0 + 2.25d0*q) / hwl2
      dwk(3) = xx(3)*dw


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
     &                 table_cs2,table_T
 
      ds=0.0005d0

      logs = dlog(s)
      
      if(logs.ge.table_s(1)) then
        i = int(( logs - table_s(1) )/ds + 1)    !must be integer
        T = ( table_T(i+1) - table_T(i))/ds*(logs - table_s(i)) 
     &      + table_T(i)
        P = ( table_P(i+1) - table_P(i))/ds*(logs - table_s(i)) 
     &      + table_P(i)
        e = ( table_e(i+1) - table_e(i))/ds*(logs - table_s(i)) 
     &      + table_e(i)
        cs2 = ( table_cs2(i+1) - table_cs2(i))/ds*(logs - table_s(i)) 
     &      + table_cs2(i)
      else
        T = table_T(1)
        P = table_P(1)
        e = table_e(1)
        cs2 = table_cs2(1)
      endif

      dwds = (cs2+1.0d0)*T
      entalp = e + P
      T=T/197.3d0
      P=P/197.3d0
      e=e/197.3d0
      dwds=dwds/197.3d0
      entalp=entalp/197.3d0

      return
      end

c=====================================================================
      subroutine tcoeff(temp,zeta)
c=====================================================================
      implicit none
      integer i
      double precision temp,zeta
      double precision x,hc,tc,a1,a2,a3,l1,l2,l3,l4,s1,s2,s3,s4
      common /unit/ hc
 
      tc = 200.0d0/hc
      a1 = -13.77d0
      a2 = 27.55d0
      a3 = 13.45d0

      l1 = 0.9d0
      l2 = 0.25d0
      l3 = 0.9d0
      l4 = 0.22d0

      s1 = 0.025d0
      s2 = 0.13d0
      s3 = 0.0025d0
      s4 = 0.022d0

      x=temp/tc

      if(temp.gt.1.05d0*tc)then
        zeta = l1*dexp(-(x-1.0d0)/s1) + l2*dexp(-(x-1.0d0)/s2) 
     &         + 0.001d0
        return
      end if

      if(temp.lt.0.995d0*tc)then
        zeta = l3*dexp((x-1.0d0)/s3) + l4*dexp((x-1.0d0)/s4)
     &         + 0.03d0
        return
      end if

      zeta = a1*x**2 + a2*x - a3
c        zeta = 0.0d0
      return
      end

c     ================================================================
      SUBROUTINE derives(r,u,s,sigma,ez,bulk,shear,v,f,ds,dez,dbulk,
     &                   dshear,t,gradP2,du,theta)
c     Computes velocity v, force f, and ds=ds/dt for each particle.
c     s(i) : entropy density per reference density
c     sigma = reference density per SPH particle
c     ds: time derivative of s(i)/ref(i)
c     dbulk: time derivative of bulk(i)/ref(i)
c     cur: current
c     dsig: gradient
c     dcur: divergence
c     dsigmadt: 
c     ================================================================
      implicit none
      integer ndim,npart,i,j,k,l,imax,ix,iy,iz,kx,ky,kz,ixk,iyk,izk,kk
      parameter(ndim=3,npart=10648,imax=100)
      double precision t,u2,hwt,hwl,ssig,wk,omega,AA,f12,
     & covisb,corelb,dcurd,thetad,bb,w0,w12,
     & ez,dez
      double precision s(npart),theta(npart),ds(npart),bulk(npart),
     &                 dbulk(npart),bulk2(npart),dbulk2(npart),
     &                 rmax(ndim),rmin(ndim),dr(ndim),link(npart),
     &                 gammaloren(npart),wgt(npart),rho(npart),
     &                 ref(npart),cur(ndim,npart),sphbulk(npart),
     &                 sphbulk2(npart),dwk(ndim),entrop(npart),
     &                 press(npart),temp(npart),entalp(npart),
     &                 dwds(npart),zeta(npart),
     &                 relaxtb(npart),prho2(npart),
     &                 atrixd1(npart),atrixd2(npart),force2(npart),
     &                 y(ndim),dsigdt(npart),dcur(npart),dw12(ndim),
     &                 sigma(npart)
      double precision r(ndim,npart),u(ndim,npart),v(ndim,npart),
     &                 f(ndim,npart),wgtv(ndim,npart),dsig(ndim,npart),
     &                 mass(ndim,ndim),fd(ndim,npart),gm(3,3),gmd(3,3),
     &                 gradP(ndim,npart),gradP2(ndim,npart)
      double precision shear(3,3,npart),dshear(3,3,npart),
     &                 sphshear(3,3,npart),btrixd(ndim,ndim,npart)
      double precision eta(npart),relaxts(npart),
     &                 shear00(npart)
      double precision shear0(3,npart),force3(ndim,npart),
     &                 dpi(ndim,npart),vdpi(ndim,npart)
      double precision atrixinv(ndim,ndim,npart),dv(ndim,ndim,npart),
     &                 du(ndim,ndim,npart)
      double precision lead(-3:imax,-3:imax,-3:imax)!
c	  double precision lead(-3:imax,-3:imax)
c     double precision lead(-3:imax)
      common /kernel/ hwt,hwl,w0
      common /bulk/ covisb,corelb
      
      do i=1,3 ! definition of metric
        do j=1,3
          gm(i,j) = 0.0d0
        enddo
      enddo

      gm(1,1) = -1.0d0
      gm(2,2) = -1.0d0
      gm(3,3) = -1.0d0/t**2

      do i=1,3
        do j=1,3
          gmd(i,j) = 0.0d0
        enddo
      enddo

      gmd(1,1) = -1.0d0
      gmd(2,2) = -1.0d0
      gmd(3,3) = -t**2

      do i=1,npart ! calculating gamma and velocity
        u2=0.0d0
        do k=1,3
          u2 = u2 + u(k,i)*u(k,i)*gm(k,k)
        enddo
        gammaloren(i) = dsqrt(1.0d0 - u2)
        
        do k=1,3
          v(k,i) = u(k,i)/gammaloren(i)*gm(k,k)
        enddo
      enddo
      
      do k=1,3 ! Linklist
        rmin(k)=1.e10
        rmax(k)=-1.e10
      end do

      do i=1,npart
        do k=1,3
          if(rmin(k).gt.r(k,i)) rmin(k)=r(k,i)
          if(rmax(k).lt.r(k,i)) rmax(k)=r(k,i)
        enddo
      enddo

      do i=1,2
        dr(i) = 2.0d0*hwt
      enddo
      dr(3) = 2.0d0*hwl
      
      call link_3d(r,rmax,rmin,dr,link,lead)
c------------------------------------------------
c     begginning of first loop
c------------------------------------------------
      do i=1,npart ! SPH 1 interpolation sigma*,s*,grad sigma*, J, div J
          
        rho(i) = 0.0d0    ! setting the initial values
        ref(i) = 0.0d0
        sphbulk(i) = 0.0d0
        shear00(i) = 0.0d0 
        do j=1,3
          shear0(j,i) = 0.0d0
          do k=1,3
            dv(j,k,i) = 0.0d0
            sphshear(j,k,i) = 0.0d0
            du(j,k,i) = 0.0d0
          end do
        end do
            
        ix=(r(1,i)-rmin(1))/dr(1) ! box coordinates
        iy=(r(2,i)-rmin(2))/dr(2)
        iz=(r(3,i)-rmin(3))/dr(3)

        do kx=1,5     ! Loops for the boxes next to the central box of i
          do ky=1,5   !   9 units cubics centered at i
            do kz=1,5
              ixk=ix-3+kx
              iyk=iy-3+ky
              izk=iz-3+kz
              
              j=lead(ixk,iyk,izk) ! j=lead(ixk,iyk)	j = lead(ixk)
              
              if(j.ne.0) then
                
 1              continue
                
                do k=1,3
                  y(k) = r(k,i) - r(k,j)
                enddo
                
                call wkernel(y,w12,dw12)
                
                rho(i) = rho(i) + sigma(j)*s(j)*w12
                ref(i) = ref(i) + sigma(j)*w12
                sphbulk(i) = sphbulk(i) 
     &                       + sigma(j)*bulk(j)*w12/gammaloren(j)/t
                do k=1,3
                  do l=1,3
                    dv(k,l,i) = dv(k,l,i) 
     &                          + sigma(j)*(v(l,j)-v(l,i))*dw12(k) 
                    sphshear(k,l,i) = sphshear(k,l,i)
     &                       + sigma(j)*shear(k,l,j)*w12/gammaloren(j)/t
                  end do
                end do
                
                j=link(j)
                
                if(j.ne.0) then ! MODIFIED
                  go to 1
                end if
              end if
            end do
          end do
        end do
c-----------------------------------------------
c     SPH 1 INTERPOLATION FINISHED
c-----------------------------------------------
        do j=1,3
          do k=1,3
            dv(j,k,i) = dv(j,k,i)/ref(i)
          end do
        end do
      
        dsigdt(i) = 0.0d0
        do j=1,3
          dsigdt(i) = dsigdt(i) - ref(i)*dv(j,j,i)
        end do
        
        do j=1,3
          do k=1,3
            do l=1,3
              du(j,k,i) = du(j,k,i) + gammaloren(i)
     &                    *( gmd(k,l) - u(k,i)*u(l,i) )*dv(j,l,i)
            end do
          end do
        end do
      
        shear00(i) = -(v(3,i)*t)**2*( sphshear(1,1,i)+sphshear(2,2,i) )
        do j=1,2
          shear00(i) = shear00(i)
     &                 + 2.0d0*v(3,i)*v(j,i)*sphshear(j,3,i)
        end do
          
        do j=1,2
          do k=1,2
            shear00(i) = shear00(i)
     &                   + v(j,i)*v(k,i)*sphshear(j,k,i)
          end do
        end do
        shear00(i) = shear00(i)/(1.0d0 - (t*v(3,i))**2)
      
        do j=1,3
          do k=1,3
            shear0(j,i) = shear0(j,i)
     &                    - v(k,i)*sphshear(j,k,i)
          end do
        end do
        
        ! Calculating the Thermodynamic Variables
        entrop(i) = rho(i)/gammaloren(i)/t 
        call eqstate (entrop(i),press(i),temp(i),entalp(i),dwds(i))
        call tcoeff(temp(i),zeta(i))
        zeta(i) =zeta(i)*entrop(i)

        ! relaxtb(i) = corelb*zeta(i)/entalp(i) relaxtb(i) = 1.0d0
        relaxtb(i) = 9.0d0*zeta(i)/(entalp(i)-4.0d0*press(i))   
        eta(i) = 2.0d0*0.08d0*entrop(i)
        relaxts(i) = eta(i)/press(i)/2.0d0
        zeta(i) = 0.0d0
        eta(i) = 0.0d0
      enddo
      
c------------------------------------------------
c     finished first particle loop
c     begginning of second particle loop
c     SPH 2 interpolation gradP + gradBulk
c-----------------------------------------------
      dez = 0.0d0
      do i=1,npart  ! setting the initial values
        sphshear(3,3,i) = shear00(i) - sphshear(1,1,i) - sphshear(2,2,i)
        sphshear(3,3,i) = sphshear(3,3,i)*t**2
        do k=1,ndim
          gradP(k,i) = 0.0d0
          gradP2(k,i) = 0.0d0
          dpi(k,i) = 0.0d0
          vdpi(k,i) = 0.0d0
        enddo

        ix=(r(1,i)-rmin(1))/dr(1)  ! box coordinates
        iy=(r(2,i)-rmin(2))/dr(2)
        iz=(r(3,i)-rmin(3))/dr(3)

        do kx=1,5 ! Loops for the boxes adjacent to the central box of i
          do ky=1,5 ! 9 units cubics centered at i
            do kz=1,5
              ixk=ix-3+kx
              iyk=iy-3+ky
              izk=iz-3+kz
              
              j=lead(ixk,iyk,izk) ! j=lead(ixk,iyk) j = lead(ixk)
              
              if(j.ne.0) then
 2              continue
                
                do k=1,ndim
                  y(k) = r(k,i)-r(k,j)
                enddo
                
                call wkernel(y,w12,dw12)
                do k=1,3
                  gradP(k,i) = gradP(k,i) 
     &                         + sigma(j)*ref(i)*(
     &	                       ( press(i) + sphbulk(i) )/ref(i)**2
     &                         + ( press(j) + sphbulk(j) )/ref(j)**2 
     &                         )*dw12(k)
                  gradP2(k,i) = gradP2(k,i) 
     &                          + sigma(j)*ref(i)*(
     &                          ( press(i) + sphbulk(i) )/ref(i)**2
     &                          + ( press(j) + sphbulk(j) )/ref(j)**2 
     &                          )*dw12(k)
                
                  do l=1,3
                    dpi(k,i) = dpi(k,i) + ref(i)*sigma(j)
     &                          *( sphshear(k,l,j)/ref(j)**2 
     &                     + sphshear(k,l,i)/ref(i)**2 )*dw12(l)*gm(l,l)
                    vdpi(k,i) = vdpi(k,i) + ref(i)*sigma(j)
     &                          *( shear0(k,j)/ref(j)**2
     &                      + shear0(k,i)/ref(i)**2 )*v(l,i)*dw12(l)
                  end do
                end do
                j=link(j)
                if(j.ne.0) go to 2
              end if
            end do
          end do
        end do

c------------------------------------------------
c     SPH 2 INTERPOLATION FINISHED
c-----------------------------------------------

        omega = entalp(i) + sphbulk(i) ! mass matrix coefficients
     &          - 0.5d0*eta(i)/relaxts(i)*(gammaloren(i)**2 - 1.0d0)
     &          /gammaloren(i)**2
        
        atrixd1(i) = entalp(i) 
     &               - dwds(i)*( entrop(i) + sphbulk(i)/temp(i) ) 
     &               - zeta(i)/relaxtb(i)

        atrixd2(i) = eta(i)/relaxts(i)*(0.5d0 - 1.0d0/3.0d0)
     &               - atrixd1(i) 
     &               + shear00(i)*( 1.0d0 - dwds(i)/temp(i) ) 

        atrixd2(i) = atrixd2(i)/gammaloren(i)

        do j=1,3
          do k=1,3
            btrixd(j,k,i) = ( -gammaloren(i)*sphshear(j,k,i)
     &                        - u(j,i)*shear0(k,i)
     &             + (1.0d0 + 1.0d0/gammaloren(i)**2)*shear0(j,i)*u(k,i)
     &                    + dwds(i)/temp(i)*shear0(k,i)*u(j,i) )*gm(k,k)
          end do 
        end do

        force2(i)=( atrixd1(i) + eta(i)/3.0d0/relaxts(i) )
     &               *( gammaloren(i)*dsigdt(i)/ref(i) - gammaloren(i)/t 
     &               + u(3,i)**2/gammaloren(i)/t**3 )
     &               + sphbulk(i)/relaxtb(i)
     &               + dwds(i)/temp(i)/t**3
     &       *( gammaloren(i)*sphshear(3,3,i) - 2.0d0*u(3,i)*shear0(3,i) 
     &               + u(3,i)**2/gammaloren(i)*shear00(i) )
        do j=1,3 ! Calculating the force coefficient /\ and \/
          do k=1,3
            force2(i) = force2(i) 
     &               - dwds(i)/temp(i)*( sphshear(j,k,i)*gm(j,j)*gm(k,k) 
     &               + shear00(i)*v(j,i)*v(k,i)
     &               - shear0(j,i)*gm(j,j)*v(k,i)
     &               - shear0(k,i)*gm(k,k)*v(j,i) )
     &               *du(j,k,i)
          end do
        end do

        do j=1,3
          force3(j,i)= -shear0(j,i)/gammaloren(i)
     &              *( gammaloren(i)*dsigdt(i)/ref(i) - 1.0d0/relaxts(i)  
     &              + u(3,i)**2/gammaloren(i)/t**3 )
     &              + ( shear0(3,i) - u(3,i)/gammaloren(i)*shear00(i) 
     &          - eta(i)*u(3,i)/gammaloren(i)/relaxts(i) )*gmd(3,j)/t**3
     &              + u(3,i)/gammaloren(i)/t**3*sphshear(j,3,i) 
          do k=1,3
	        force3(j,i) = force3(j,i) 
     &                  + 0.5d0*eta(i)/gammaloren(i)/relaxts(i)*v(k,i)
     &                  *( du(j,k,i) + du(k,j,i) )
          end do
        end do

        do k=1,3 ! Calculating the force
          f(k,i) = force2(i)*u(k,i) + gradP(k,i) + force3(k,i) 
     &             - dpi(k,i) + vdpi(k,i)
          fd(k,i) = f(k,i)
        enddo 
        
        mass(1,1) = gammaloren(i)*omega ! Calculating the mass matrix
     &              + atrixd2(i)*u(1,i)*u(1,i)*gm(1,1) 
     &              + btrixd(1,1,i)
        mass(2,1) = atrixd2(i)*u(2,i)*u(1,i)*gm(1,1) 
     &              + btrixd(2,1,i)
        mass(3,1) = atrixd2(i)*u(3,i)*u(1,i)*gm(1,1) 
     &              + btrixd(3,1,i)
        mass(1,2) = atrixd2(i)*u(1,i)*u(2,i)*gm(2,2) 
     &              + btrixd(1,2,i)
        mass(2,2) = gammaloren(i)*omega
     &              + atrixd2(i)*u(2,i)*u(2,i)*gm(2,2) 
     &              + btrixd(2,2,i)
        mass(3,2) = atrixd2(i)*u(3,i)*u(2,i)*gm(2,2) 
     &              + btrixd(3,2,i)
        mass(1,3) = atrixd2(i)*u(1,i)*u(3,i)*gm(3,3) 
     &              + btrixd(1,3,i)
	    mass(2,3) = atrixd2(i)*u(2,i)*u(3,i)*gm(3,3) 
     &              + btrixd(2,3,i)
	    mass(3,3) = gammaloren(i)*omega
     &              + atrixd2(i)*u(3,i)*u(3,i)*gm(3,3) 
     &              + btrixd(3,3,i)
      
        call gaussj(mass, 3, 3, fd, 1, 1) !Inverting the mass matrix 
	  
        do j=1,3 ! Calculating the acceleration
          fd(j,i) = 0.0d0
          do k=1,3
            fd(j,i) = fd(j,i) + mass(j,k)*f(k,i)
          end do
        end do 
        
        do j=1,3
          f(j,i) = fd(j,i)
        end do

        thetad = 0.0d0 ! Calculating the divergence of the velocity
        
        do k=1,3
          thetad = thetad - u(k,i)*f(k,i)/gammaloren(i)*gm(k,k)
        enddo

        theta(i) = thetad - gammaloren(i)/ref(i)*dsigdt(i)
     & -u(3,i)**2/gammaloren(i)/t**3


        ds(i) =	- sphbulk(i) ! Calculating ds_dt and dBulk_dt
     &            *t/temp(i)/ref(i)*( theta(i) + gammaloren(i)/t )
     &           +t/temp(i)/ref(i)*( -gammaloren(i)/t**3*sphshear(3,3,i)
     &           + 2.0d0*u(3,i)/t**3*shear0(3,i) 
     &           - u(3,i)**2/gammaloren(i)/t**3*shear00(i) )
        do j=1,3
          ds(i) = ds(i) + t/temp(i)/ref(i)
     &            *( shear0(j,i)*gm(j,j) - shear00(i)*v(j,i) )*f(j,i)
          do k=1,3
            ds(i) = ds(i) + t/temp(i)/ref(i)
     &              *( sphshear(j,k,i)*gm(j,j)*gm(k,k)
     &              + shear00(i)*v(j,i)*v(k,i)
     &              - shear0(j,i)*gm(j,j)*v(k,i)
     &              -shear0(k,i)*gm(k,k)*v(j,i) )*du(j,k,i)
          end do
        end do
        
        dbulk(i) = - t/relaxtb(i)/ref(i)*(
     &               sphbulk(i)
     &               + zeta(i)
     &               *(theta(i) + gammaloren(i)/t)
     &               )
                
        do j=1,3
          do k=1,3
            dshear(j,k,i) = - sphshear(j,k,i)*t/ref(i)/relaxts(i) 
     &           +0.5d0*eta(i)*t/ref(i)/relaxts(i)*(du(j,k,i)+du(k,j,i))
     &           -0.5d0*eta(i)*t*gammaloren(i)/ref(i)/relaxts(i)
     &           *(u(j,i)*f(k,i)+u(k,i)*f(j,i)) 
     &           -eta(i)/3.0d0*t/ref(i)/relaxts(i)
     &          *( gmd(j,k)-u(j,i)*u(k,i) )*(theta(i) + gammaloren(i)/t)
     &           -eta(i)*t**2*gammaloren(i)/ref(i)/relaxts(i)
     &           *(-t**2*gm(j,3))*(-t**2*gm(k,3))
     &           +gammaloren(i)/ref(i)*( (-t**2*gm(j,3))*sphshear(k,3,i) 
     &           +(-t**2*gm(k,3))*sphshear(j,3,i) )
     &          +t**2/ref(i)*gm(3,3)*u(3,i)*((-t**2*gm(j,3))*shear0(k,i) 
     &           +(-t**2*gm(k,3))*shear0(j,i) )
            do l=1,3
              dshear(j,k,i) = dshear(j,k,i) 
     &                        -(u(j,i)*sphshear(k,l,i)*gm(l,l)
     &    +u(k,i)*sphshear(j,l,i)*gm(l,l))*f(l,i)/ref(i)*t*gammaloren(i)
     &               +( u(j,i)*shear0(k,i)+u(k,i)*shear0(j,i) )*t/ref(i)
     &               *u(l,i)*gm(l,l)*f(l,i)
            end do
          end do
        end do
        
        dez = dez + sigma(i)*t**2/ref(i)
     &                 *( 
     &            (entalp(i) + sphbulk(i))*u(3,i)*u(3,i)*gm(3,3)*gm(3,3) 
     &            +( press(i) + sphbulk(i) )/t**2 
     &             +sphshear(3,3,i)*gm(3,3)*gm(3,3) )
      end do
      return
      END
      
c==========================================
      subroutine step(x,u,s,sigma,ez,bulk,shear,t,gradP2,f1,du,theta)
c========================================== 
      implicit none
      integer ndim,npart,imax
      parameter (ndim=3,npart=10648,imax=100)
      integer i,j,k,n
      double precision x(ndim,npart),u(ndim,npart),v(ndim,npart),
     &                 f(ndim,npart)
      double precision x1(ndim,npart),u1(ndim,npart),v1(ndim,npart),
     &                 f1(ndim,npart)
      double precision gradP2(ndim,npart)
      double precision gammaloren(npart),entrop(npart),sigma(npart),
     &                 s(npart),ds(npart),rho(npart)
      double precision s1(npart),ds1(npart)
      double precision theta(npart)
      double precision ez,t,dez,dt
      double precision ez1,dez1,hwt,hwl,w0
      double precision link(npart),rmin(ndim),rmax(ndim),wgt(npart),
     &                 dr(ndim)
      double precision lead(-3:imax,-3:imax,-3:imax)
                     ! lead(-3:imax,-3:imax)
      double precision bulk(npart),bulk1(npart),
     &                 dbulk(npart),dbulk1(npart)
      double precision shear(3,3,npart),shear1(3,3,npart),
     &                 dshear(3,3,npart),dshear1(3,3,npart)
      double precision du(3,3,npart)
      common /time/ dt
      common /kernel/ hwt,hwl,w0

      call derives(x,u,s,sigma,ez,bulk,shear,v,f,ds,dez,dbulk,dshear,t
     &             ,gradP2,du,theta)

      do i=1,npart
        do j=1,3
          x1(j,i) = x(j,i) + 0.5d0*dt*v(j,i)
          u1(j,i) = u(j,i) + 0.5d0*dt*f(j,i)
        end do
        
        do j=1,3
          do k=1,3
            shear1(k,j,i) = shear(k,j,i) + 0.5d0*dt*dshear(k,j,i)
          end do
        end do
        
c       do j=1,3
c	      shear1(1,j,i) = shear(1,j,i) + 0.5d0*dt*dshear(1,j,i)
c	    end do
c	    do j=2,3
c	      shear1(2,j,i) = shear(2,j,i) + 0.5d0*dt*dshear(2,j,i)	  
c	    end do
c       shear1(3,3,i) = shear(3,3,i) + 0.5d0*dt*dshear(3,3,i)
      
c	    shear1(2,1,i) = shear1(1,2,i)
c	    do j=1,2
c         shear1(3,j,i) = shear1(j,3,1)
c	    end do
        s1(i) = s(i) + 0.5d0*dt*ds(i)
        bulk1(i) = bulk(i) + 0.5d0*dt*dbulk(i)
      end do
      
      ez1 = ez + 0.5d0*dt*dez
      t = t + 0.5d0*dt
      
      call derives(x1,u1,s1,sigma,ez1,bulk1,shear1,v1,f1,ds1,dez1,
     &             dbulk1,dshear1,t,gradP2,du,theta)
      
      do i=1,npart
        do j=1,ndim
          x(j,i) = x(j,i) + dt*v1(j,i)
	        u(j,i) = u(j,i) + dt*f1(j,i)
        end do
        
        do j=1,3
	      do k=1,3
	        shear(k,j,i) = shear(k,j,i) + dt*dshear1(k,j,i)
          end do
        end do
c       do j=1,3
c         shear(1,j,i) = shear(1,j,i) + dt*dshear1(1,j,i)
c	    end do
c       do j=2,3
c         shear(2,j,i) = shear(2,j,i) + dt*dshear1(2,j,i)	  
c       end do
c       shear(3,3,i) = shear(3,3,i) + dt*dshear1(3,3,i)
c	    shear(2,1,i) = shear(1,2,i)
c	    do j=1,2
c         shear(3,j,i) = shear(j,3,1)
c       end do
        s(i) = s(i) + dt*ds1(i)
        bulk(i) = bulk(i) + dt*dbulk1(i)
      end do
      
      ez = ez + dt*dez1
      t = t + 0.5d0*dt
      return
      end

c     ================================================
 
      subroutine link_3d(r,rmax,rmin,dr,link,lead)

c     This routine generates just generates the table for linklist,
c     that is,  lead(ix,iy,iz) and link(i)
c
c     rmax, rmin are vectors which define the limit of the coordinates.
c                   rmin(k)< r(i,k) < rmax(k)  for all i  and k=1,2,3.
c
c     dr is a vector which defines the range of the binary force.
c     For a given particle i, those particles j which satisfy 
c                          |r(i,k)-r(j,k)| < dr(k) 
c     are considered in the sum of the forces.  
c     
c     Caution:
c     
c                |rmax(k)-rmin(k)|/dx(k) < imax
c     and always should update the value of rmax and rmin whenever the
c     coordinates are modified.
c     
c     common line
c 
c   	common lead(-3:imax,-3:imax,-3:imax),link(n)  
c
c     should be copied in the program which
c     uses these lists, together with the parameter, imax and n values. 
c 

      implicit none 
      integer npart,ndim,imax,ix,iy,iz,i,kx,ky,kz,ixk,iyk,izk
      parameter(ndim=3,npart=10648,imax=100)
      double precision rmin(ndim),rmax(ndim),dr(ndim),link(npart)
      double precision r(ndim,npart)
      double precision lead(-3:imax,-3:imax,-3:imax)
c     double precision lead(-3:imax,-3:imax)
c     double precision lead(-3:imax)
      
      do i=1,npart
        ix=(r(1,i)-rmin(1))/dr(1)
        iy=(r(2,i)-rmin(2))/dr(2)
        iz=(r(3,i)-rmin(3))/dr(3)
        
        do kx=1,5
          do ky=1,5
            do kz=1,5
              ixk=ix-3+kx
              iyk=iy-3+ky
              izk=iz-3+kz
              
              lead(ixk,iyk,izk)=0
c             lead(ixk,iyk)=0
c             lead(ixk)=0
            end do
          end do
        end do
      end do

      do i=1,npart
        ix=(r(1,i)-rmin(1))/dr(1)
        iy=(r(2,i)-rmin(2))/dr(2)
        iz=(r(3,i)-rmin(3))/dr(3)
	    
        link(i)=lead(ix,iy,iz)
        lead(ix,iy,iz)=i

c	    link(i)=lead(ix,iy)
c	    lead(ix,iy)=i
c       link(i)=lead(ix)
c	    lead(ix)=i
      end do
      
      return
      end

C     =====================================================================
C 
       SUBROUTINE gaussj(a,n,np,b,m,mp)  
C
C     =====================================================================
 
      INTEGER m,mp,n,np,NMAX
      double precision a(np,np),b(np,mp)
      PARAMETER (NMAX=50)
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX) 
      double precision big,dum,pivinv  
 
      do 11 j=1,n
         ipiv(j)=0   
11    continue  
       
      do 22 i=1,n
        big=0.  
        do 13 j=1,n 
          if(ipiv(j).ne.1)then 
            do 12 k=1,n 
              if (ipiv(k).eq.0) then 
                if (abs(a(j,k)).ge.big)then 
                  big=abs(a(j,k)) 
                  irow=j 
                  icol=k 
                endif  
                   
                else if (ipiv(k).gt.1) then 
                  pause 'singular matrix in gaussj'
                endif  
12          continue
          endif  
13      continue  
         
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)  
            a(irow,l)=a(icol,l)
            a(icol,l)=dum 
14        continue  
          do 15 l=1,m  
            dum=b(irow,l)
            b(irow,l)=b(icol,l) 
            b(icol,l)=dum
15        continue  
        endif  
        indxr(i)=irow   
        indxc(i)=icol  
        
        if (a(icol,icol).eq.0.) pause 'singular matrix in gaussj'  
        
        pivinv=1./a(icol,icol)
        a(icol,icol)=1.
        
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
16      continue 
        
        do 17 l=1,m 
          b(icol,l)=b(icol,l)*pivinv 
17      continue

        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0.
            do 18 l=1,n 
              a(ll,l)=a(ll,l)-a(icol,l)*dum 
18          continue 
            do 19 l=1,m 
              b(ll,l)=b(ll,l)-b(icol,l)*dum  
19          continue 
          endif
21      continue
22    continue 

      do 24 l=n,1,-1 
        if(indxr(l).ne.indxc(l)) then 
          do 23 k=1,n 
            dum=a(k,indxr(l)) 
            a(k,indxr(l))=a(k,indxc(l)) 
            a(k,indxc(l))=dum 
23        continue 
        endif
24    continue 
      
      return 
      END  

c     ===============================================
      subroutine sum(f,wgt,r,rmin,rmax,dr,link,lead)
c     summation of the sph particle using the linklist
c     ===============================================
      implicit none
      integer npart,ndim,i,j,imax,k,ix,iy,iz,kx,ky,kz,ixk,iyk,izk
      parameter (ndim=3,npart=10648,imax=100)
      double precision w12,w0,hwt,hwl
      double precision f(npart),wgt(npart),y(ndim),rmin(ndim),
     &                 rmax(ndim),dr(ndim),link(npart),dw12(ndim)
      double precision r(ndim,npart)
      double precision lead(-3:imax,-3:imax,-3:imax)
c     double precision lead(-3:imax,-3:imax)
c	    double precision lead(-3:imax)
      common /kernel/ hwt,hwl,w0

      do i=1,npart
        f(i) = wgt(i)*w0
      enddo

      do i=1,npart
        ix=(r(1,i)-rmin(1))/dr(1)
        iy=(r(2,i)-rmin(2))/dr(2)
        iz=(r(3,i)-rmin(3))/dr(3)
          
        do kx=1,5 ! Loops for the boxes adjacent to the central box of i
          do ky=1,5 ! 9 units cubics centered at i
            do kz=1,5
              ixk=ix-3+kx
              iyk=iy-3+ky
              izk=iz-3+kz
              
              j=lead(ixk,iyk,izk) ! j=lead(ixk,iyk) j = lead(ixk)
              if(j.ne.0) then
2               continue
                if(i.ne.j) then  ! found one j within distance condition
                  do k=1,3
                    y(k) = r(k,i)-r(k,j)
                  enddo
                  call wkernel(y,w12,dw12)
                  
                  if(w12.gt.0.d0)then
                    f(i) = f(i) + wgt(j)*w12
                  endif
c                 Sum = Sum + f(i,j)  
                end if
	              j=link(j)
                if(j.ne.0) go to 2
              end if
            end do
          end do
        end do
      end do
      
      return
      end

c==========================================
      subroutine hyper0(x,u,s,sigma,entrop,t,rmin,rmax,dr,link,lead,
     &                  frag,frag2,frag3,ii)
c==========================================
      implicit none
      integer ndim,npart,imax
      parameter (ndim=3,npart=10648,imax=100)
      integer i,j,ii
      integer frag(npart),frag2(npart),frag3(npart)
      integer ix,iy,iz,k,l,kx,ky,kz,ixk,iyk,izk
      double precision x(ndim,npart),u(ndim,npart),entrop(npart),
     &                 gm(3,3),gmd(3,3),v(ndim,npart)
      double precision press(npart),temp(npart),entalp(npart),
     &                 dwds(npart),gammaloren(npart),sigma(npart),
     &                 s(npart),rho(npart)
      double precision t,hc,Tout,Tcut
      double precision link(npart),rmin(ndim),rmax(ndim),wgt(npart),
     &                 dr(ndim)
      double precision w12
      double precision ref(npart),y(ndim),dw12(ndim)
      double precision gradP2(ndim,npart)
      double precision lead(-3:imax,-3:imax,-3:imax) !
c     double precision lead(-3:imax,-3:imax)
      double precision dv(ndim,ndim,npart)
      double precision hwt,hwl,w0
      common /kernel/ hwt,hwl,w0
      common /unit/ hc
c     save frag,frag2,frag3
c     save ii

      do i=1,3
        do j=1,3
          gm(i,j) = 0.0d0
        enddo
      enddo

      gm(1,1) = -1.0d0
      gm(2,2) = -1.0d0
      gm(3,3) = -1.0d0/t**2

      do i=1,3
        do j=1,3
          gmd(i,j) = 0.0d0
        enddo
      enddo

      gmd(1,1) = -1.0d0
      gmd(2,2) = -1.0d0
      gmd(3,3) = -t**2

      do i=1,npart
	    gammaloren(i) = 1.0d0
        do j=1,ndim
          gammaloren(i)=gammaloren(i)-u(j,i)*u(j,i)*gm(j,j)
        end do
        gammaloren(i) = dsqrt(gammaloren(i))
        do j= 1,ndim
          v(j,i) = u(j,i)/gammaloren(i)*gm(j,j)
        end do
      end do

      do i=1,ndim
        rmin(i)=1.e10
        rmax(i)=-1.e10
      end do

      do i=1,npart
        do j=1,ndim
          if(rmin(j).gt.x(j,i)) rmin(j)=x(j,i)
          if(rmax(j).lt.x(j,i)) rmax(j)=x(j,i)
        enddo
      enddo
      
      do i=1,2
        dr(i) = 2.0d0*hwt
      enddo
      dr(3) = 2.0d0*hwl
      
      call link_3d(x,rmax,rmin,dr,link,lead)
      
      do i= 1,npart
        wgt(i) = sigma(i)*s(i)
      end do
      call sum(rho,wgt,x,rmin,rmax,dr,link,lead)
      
c     *******************************
c      open(unit=88,file='test.dat')
c      do i=1,npart
c      write(88,*) i, sigma(i), s(i)
c      end do
c      pause
c     *******************************

      do i=1,npart
        entrop(i) = rho(i)/gammaloren(i)/t
        call eqstate(entrop(i),press(i),temp(i),entalp(i),dwds(i))
        Tout = 130.0d0/hc
        if(temp(i).le.Tout)then
          if(frag(i).eq.0)then
            frag(i) = 2
            frag2(i) = 2
            ii = ii + 1
            write(*,*) ii
          else
          end if
        else
        end if
      end do
      
      return
      end

c==========================================
      subroutine hyper(x,u,s,sigma,entrop,t,rmin,rmax,dr,link,lead
     &           ,gradP2,f1,du,shear,bulk,theta,frag,frag2,frag3,ii)
c==========================================
      implicit none
      integer ndim,npart,imax
      parameter (ndim=3,npart=10648,imax=100)
      integer i,j,k,l,ii
      integer frag(npart),frag2(npart),frag3(npart)
      integer ix,iy,iz,kx,ky,kz,ixk,iyk,izk,ifirst
      double precision x(ndim,npart),u(ndim,npart),entrop(npart),
     &                 gm(3,3),gmd(3,3),v(ndim,npart)
      double precision press(npart),temp(npart),entalp(npart),
     &                 dwds(npart),gammaloren(npart),sigma(npart),
     &                 s(npart),rho(npart)
      double precision hc,Tout,t,Tcut
      double precision link(npart),rmin(ndim),rmax(ndim),wgt(npart),
     &                 dr(ndim)
      double precision w12
      double precision ref(npart),y(ndim),dw12(ndim)
      double precision f1(ndim,npart),gradP2(ndim,npart)
      double precision shaa(npart),sphbulk(npart),bulk(npart),
     &                 theta(npart)
      double precision shear(3,3,npart),sphshear(3,3,npart)
      double precision lead(-3:imax,-3:imax,-3:imax)
c     double precision lead(-3:imax,-3:imax)
      double precision du(3,3,npart)
      double precision hwt,hwl,w0
      double precision tempo(npart),entropo(npart),so(npart),refo(npart)
     &                 ,to(npart),thetao(npart)
      double precision tempo2(npart),entropo2(npart),so2(npart),
     &                 refo2(npart),to2(npart),thetao2(npart)
      double precision xo(ndim,npart),vo(ndim,npart),
     &                 gradP2o(ndim,npart),f1o(ndim,npart)
      double precision xo2(ndim,npart),vo2(ndim,npart),
     &                 gradP2o2(ndim,npart),f1o2(ndim,npart)
      double precision duo(ndim,ndim,npart)
      double precision duo2(ndim,ndim,npart)
      double precision entropn(npart),sn(npart),refn(npart),tn(npart),
     &                 thetan(npart)
      double precision entropn2(npart),sn2(npart),refn2(npart),
     &                 tn2(npart),thetan2(npart)
      double precision xn(ndim,npart),vn(ndim,npart),
     &                 gradP2n(ndim,npart),f1n(ndim,npart)
      double precision xn2(ndim,npart),vn2(ndim,npart),
     &                 gradP2n2(ndim,npart),f1n2(ndim,npart)
      double precision dun(ndim,ndim,npart)
      double precision dun2(ndim,ndim,npart)
      common /kernel/ hwt,hwl,w0
      common /unit/ hc
      data ifirst/0/
      save ifirst

      if(ifirst.eq.0) then
        open(unit=92,file='hyper_130.dat')
        open(unit=95,file='cut_100.dat')
        open(unit=96,file='deiri.dat')
        ifirst=1
      else
      end if

      do i=1,3
        do j=1,3
          gm(i,j) = 0.0d0
        enddo
      enddo

      gm(1,1) = -1.0d0
      gm(2,2) = -1.0d0
      gm(3,3) = -1.0d0/t**2

      do i=1,3
        do j=1,3
          gmd(i,j) = 0.0d0
        enddo
      enddo

      gmd(1,1) = -1.0d0
      gmd(2,2) = -1.0d0
      gmd(3,3) = -t**2

      do i=1,npart
        gammaloren(i) = 1.0d0
        do j=1,ndim
          gammaloren(i)=gammaloren(i)-u(j,i)*u(j,i)*gm(j,j)
        end do
        gammaloren(i) = dsqrt(gammaloren(i))
        do j= 1,ndim
          v(j,i) = u(j,i)/gammaloren(i)*gm(j,j)
        end do
      end do
      
      do i=1,ndim
        rmin(i)=1.e10
        rmax(i)=-1.e10
      end do

      do i=1,npart
        do j=1,ndim
          if(rmin(j).gt.x(j,i)) rmin(j)=x(j,i)
          if(rmax(j).lt.x(j,i)) rmax(j)=x(j,i)
        enddo
      enddo

      do i=1,2
        dr(i) = 2.0d0*hwt
      enddo
      dr(3) = 2.0d0*hwl

      call link_3d(x,rmax,rmin,dr,link,lead)

      do i=1,npart ! preparation for the interpolation of hypersurface
        if(frag(i).eq.0) then
          tempo(i) = 0.0d0
          to(i) = t
          entropo(i) = entrop(i)
          so(i) = s(i)
          refo(i) = ref(i)
          thetao(i) = theta(i)
          
          do k=1,ndim
            xo(k,i) = x(k,i)
            vo(k,i) = v(k,i)
            gradP2o(k,i) = gradP2(k,i)
            f1o(k,i) = f1(k,i)
            do l=1,ndim
              duo(k,l,i) = du(k,l,i)
            end do
          end do
                
          frag(i)=1
        end if
          
        if(frag2(i).eq.0) then
          tempo2(i) = 0.0d0
          to2(i) = t
          entropo2(i) = entrop(i)
          so2(i) = s(i)
          refo2(i) = ref(i)
          thetao2(i) = theta(i)
          do k=1,ndim
            xo2(k,i) = x(k,i)
            vo2(k,i) = v(k,i)
            gradP2o2(k,i) = gradP2(k,i)
            f1o2(k,i) = f1(k,i)
            do l=1,ndim
              duo2(k,l,i) = du(k,l,i)
            end do
          end do
          
          frag2(i)=1
        end if
      end do

      do i= 1,npart
        wgt(i) = sigma(i)*s(i)
      end do
      call sum(rho,wgt,x,rmin,rmax,dr,link,lead)
      
      do l=1,npart
        wgt(l) = sigma(l)
      enddo
      call sum(ref,wgt,x,rmin,rmax,dr,link,lead)
      
      do i=1,npart ! freezeout & cut
        
        entrop(i) = rho(i)/gammaloren(i)/t
        call eqstate(entrop(i),press(i),temp(i),entalp(i),dwds(i))
        Tout = 130.0d0/hc ! freezeout
        if(frag(i).eq.3.and.temp(i).gt.Tout) then
          frag(i) = 1
          write(96,*) i, 'pass1'
        else
        endif
        
        if(frag(i).eq.1)then
          if(temp(i).le.Tout) then
            tn(i) = to(i) + (Tout-tempo(i))/(temp(i)-tempo(i))*(t-to(i))
            entropn(i) = entropo(i) + (Tout-tempo(i))/(temp(i)-tempo(i))
     &                   *(entrop(i)-entropo(i))
            sn(i)=so(i) +(Tout-tempo(i))/(temp(i)-tempo(i))*(s(i)-so(i))
            refn(i) = refo(i) + (Tout-tempo(i))/(temp(i)-tempo(i))
     &                   *(ref(i)-refo(i))
            thetan(i) = thetao(i) + (Tout-tempo(i))/(temp(i)-tempo(i))
     &                   *(theta(i)-thetao(i))
            
            do k=1,ndim
              xn(k,i) = xo(k,i) 
     &              +(Tout-tempo(i))/(temp(i)-tempo(i))*(x(k,i)-xo(k,i))
              vn(k,i) = vo(k,i) 
     &              +(Tout-tempo(i))/(temp(i)-tempo(i))*(v(k,i)-vo(k,i))
              gradP2n(k,i) = gradP2o(k,i) + 
     &     (Tout-tempo(i))/(temp(i)-tempo(i))*(gradP2(k,i)-gradP2o(k,i))
              f1n(k,i) = f1o(k,i) 
     &            +(Tout-tempo(i))/(temp(i)-tempo(i))*(f1(k,i)-f1o(k,i))
           
              do l=1,ndim
                dun(k,l,i) = duo(k,l,i) 
     &       + (Tout-tempo(i))/(temp(i)-tempo(i))*(du(k,l,i)-duo(k,l,i))
              end do
            end do
            
            frag(i) = 3
          else
            tempo(i) = temp(i)
            to(i) = t
            entropo(i) = entrop(i)
            so(i) = s(i)
            refo(i) = ref(i)
            thetao(i) = theta(i)
            do k=1,ndim
              xo(k,i) = x(k,i)
              vo(k,i) = v(k,i)
              gradP2o(k,i) = gradP2(k,i)
              f1o(k,i) = f1(k,i)
              do l=1,ndim
                duo(k,l,i) = du(k,l,i)
              end do
            end do
          endif
        else
        end if
            
        Tcut = 100.0d0/hc ! cut
c        
c     search the sph part which pass the hypersurface at the last time
c       
        if(frag2(i).eq.3.and.temp(i).gt.Tout) then
          frag2(i) = 1
          ii = ii -1
          write(96,*) i, 'pass2'
        else
        endif
      
        if(frag2(i).eq.1)then
          if(temp(i).le.Tcut) then
            tn2(i) = to2(i) + (Tcut-tempo2(i))/(temp(i)-tempo2(i))
     &               *(t-to2(i))
          entropn2(i) = entropo2(i)+(Tcut-tempo2(i))/(temp(i)-tempo2(i))
     &               *(entrop(i)-entropo2(i))
            sn2(i) = so2(i) + (Tcut-tempo2(i))/(temp(i)-tempo2(i))
     &               *(s(i)-so2(i))
            refn2(i) = refo2(i) + (Tcut-tempo2(i))/(temp(i)-tempo2(i))
     &               *(ref(i)-refo2(i))
            thetan2(i) = thetao2(i)+(Tout-tempo2(i))/(temp(i)-tempo2(i))
     &               *(theta(i)-thetao2(i))
            do k=1,ndim
              xn2(k,i) = xo2(k,i)
     &           +(Tcut-tempo2(i))/(temp(i)-tempo2(i))*(x(k,i)-xo2(k,i))
              vn2(k,i) = vo2(k,i)
     &           +(Tcut-tempo2(i))/(temp(i)-tempo2(i))*(v(k,i)-vo2(k,i))
              gradP2n2(k,i) = gradP2o2(k,i) + (Tout-tempo2(i))
     &                  /(temp(i)-tempo2(i))*(gradP2(k,i)-gradP2o2(k,i))
              f1n2(k,i) = f1o2(k,i)
     &         +(Tcut-tempo2(i))/(temp(i)-tempo2(i))*(f1(k,i)-f1o2(k,i))
              do l=1,ndim
                dun2(k,l,i) = duo2(k,l,i)
     &     +(Tout-tempo2(i))/(temp(i)-tempo2(i))*(du(k,l,i)-duo2(k,l,i))
              end do
            end do
            
            frag2(i) = 3
            ii=ii+1
            write(*,*)ii
          else
            tempo2(i) = temp(i)
            to2(i) = t
            entropo2(i) = entrop(i)
            so2(i) = s(i)
            refo2(i) = ref(i)
            thetao2(i) = theta(i)
            do k=1,ndim
              xo2(k,i) = x(k,i)
              vo2(k,i) = v(k,i)
              gradP2o2(k,i) = gradP2(k,i)
              f1o2(k,i) = f1(k,i)
              do l=1,ndim
                duo2(k,l,i) = du(k,l,i)
              end do
            end do
          endif
        else
        end if
      
        if(ii.eq.npart) goto 1111
      end do
      
1111  continue
       
      if(ii.eq.npart) then
        do i=1,npart
          if(frag(i).eq.3) then
            write(unit=92,fmt=10) tn(i),xn(1,i),xn(2,i),xn(3,i),vn(1,i),
     &                vn(2,i),vn(3,i),
     &                entropn(i),gradP2n(1,i),gradP2n(2,i),gradP2n(3,i),
     &                dun(1,1,i),dun(1,2,i),dun(2,1,i),dun(2,2,i),
     &                f1n(1,i),f1n(2,i),thetan(i),0.0d0,0.0d0,0.0d0,
     &                0.0d0,0.0d0,sigma(i),sn(i),refn(i)
          else
          endif
          
          if(frag2(i).eq.3) then
            write(unit=95,fmt=10) tn2(i),xn2(1,i),xn2(2,i),
     &            xn2(3,i),vn2(1,i),vn2(2,i),vn2(3,i),
     &            entropn2(i),gradP2n2(1,i),gradP2n2(2,i),gradP2n2(3,i),
     &            dun2(1,1,i),dun2(1,2,i),dun2(2,1,i),dun2(2,2,i),
     &            f1n2(1,i),f1n2(2,i),thetan2(i),0.0d0,0.0d0,0.0d0,
     &            0.0d0,0.0d0,sigma(i),sn2(i),refn2(i)
          else
          endif
        enddo
        
        stop
      else
      endif

10    format (30(1x,E12.4))
      return
      end
       
c==========================================
      subroutine out(x,u,entrop,sigma,ref,t,esum,esum0,ez,bulk,shear,
     &               rmin,rmax,dr,link,lead,outfile)
c==========================================
      implicit none
      character*7 outfile
      integer ndim,npart,imax
      parameter (ndim=3,npart=10648,imax=100)
      integer i,j,k
      double precision x(ndim,npart),v(ndim,npart),u(ndim,npart),
     &                 gm(3,3),gmd(3,3)
      double precision entrop(npart),entalp(npart),press(npart),
     &                 gammaloren(npart),temp(npart),dwds(npart),
     &                 sigma(npart),ref(npart)
      double precision t,esum,esum0,ez,hc,gam
      double precision link(npart),rmin(ndim),rmax(ndim),wgt(npart),
     &                 dr(ndim),shaa(npart)
      double precision lead(-3:imax,-3:imax,-3:imax)
c     double precision lead(-3:imax,-3:imax)
      double precision bulk(npart),sphbulk(npart)
      double precision shear(3,3,npart),sphshear(3,3,npart)
      double precision shear00(npart)
      double precision hwt,hwl,w0
      common /unit/ hc
      common /kernel/ hwt,hwl,w0

      write(6,*) 'Writing output...'

      open(unit=33,file=outfile)
 
      do i=1,3
        do j=1,3
          gm(i,j) = 0.0d0
        enddo
      enddo
      
      gm(1,1) = -1.0d0
      gm(2,2) = -1.0d0
      gm(3,3) = -1.0d0/t**2
      
      do i=1,3
        do j=1,3
	        gmd(i,j) = 0.0d0
        enddo
      enddo
      
      gmd(1,1) = -1.0d0
      gmd(2,2) = -1.0d0
      gmd(3,3) = -t**2
      
      do i=1,npart
        gammaloren(i) = 1.0d0
        do j=1,3
          gammaloren(i) = gammaloren(i)-u(j,i)*u(j,i)*gm(j,j)
        end do
        gammaloren(i) = dsqrt(gammaloren(i))
	      do j= 1,3
	        v(j,i) = u(j,i)/gammaloren(i)*gm(j,j)
	      end do
      end do

      do i=1,npart
        wgt(i) = sigma(i)
      end do
      call sum(ref,wgt,x,rmin,rmax,dr,link,lead)
      
      do i=1,npart
        wgt(i) = sigma(i)*bulk(i)/gammaloren(i)/t
      end do
      call sum(sphbulk,wgt,x,rmin,rmax,dr,link,lead)

      do j=1,3
        do k=1,3
          do i=1,npart
            wgt(i) = sigma(i)*shear(j,k,i)/gammaloren(i)/t
          end do
          call sum(shaa,wgt,x,rmin,rmax,dr,link,lead)
          do i=1,npart
            sphshear(j,k,i) = shaa(i)
          end do
        end do
      end do

      esum=0.0d0
      do i=1,npart
        shear00(i) = -(v(3,i)*t)**2*( sphshear(1,1,i)+sphshear(2,2,i) )
        do j=1,2
          shear00(i) = shear00(i)
     &               + 2.0d0*v(3,i)*v(j,i)*sphshear(j,3,i)
          do k=1,2
            shear00(i) = shear00(i)
     &               + v(j,i)*v(k,i)*sphshear(j,k,i)
          end do
        end do
        shear00(i) = shear00(i)/(1.0d0 - (t*v(3,i))**2)
        call eqstate(entrop(i),press(i),temp(i),entalp(i),dwds(i))
        
        esum = esum + sigma(i)*t/ref(i)*
     &         ( (entalp(i) + sphbulk(i) )*gammaloren(i)*gammaloren(i) 
     &           - press(i) - sphbulk(i) + shear00(i) )

        write(unit=33,fmt=10) x(1,i),x(2,i),x(3,i),v(1,i),v(2,i),v(3,i),
     &                   entrop(i),temp(i)*hc,press(i)*hc,sphbulk(i)*hc,
     &                   shear00(i)*hc,shear(1,2,i)*hc,shear(1,1,i)*hc,
     &                   shear(2,2,i)*hc
      end do
      
      if(esum0.eq.0.0d0) esum0=esum

      write(6,*) esum+ez,(esum0-esum-ez)/esum0
      open(unit=34,file='ener_con.dat')
      write(34,10) t,esum+ez,(esum0-esum-ez)/esum0

10    format (20(1x,E12.4))
      return
      end 
