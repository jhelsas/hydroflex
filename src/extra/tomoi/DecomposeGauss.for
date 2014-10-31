c     Test for 2D Gauss fit by N small Gauss
      implicit real*8 (a-h,o-z)
	parameter (ndim=500)
	parameter (ndim2=ndim/2)
	dimension r(2,Ndim2), gzi(2)
      
      open(1, file='Firstfit_X_Y_W.dat')
	open(2, file='SPHparticles.dat')
	open(3, file='finalfit.dat')
	open(4, file='ini.dat')
      open(5, file='ini_v.dat')

      read(1,*) NSPH
	read(1,*) H0, Wmax,Wmin

      q=3.
	h=h0/q

      Npart=0

      do NP=1,NSPH

      read(1,*) X0,Y0,W
    
      N=W/Wmin+0.1
      Wi=W/N

      if(NP.eq.1) write(4,*)Wi*100.0d0

	write(*,*) '***', NP, ' of ', NSPH, ' ***'
	write(*,*) 'N =', N, 'Nacc=', Npart

      call MinGaussfit(n,r,H0,h)


      do i=1,N
	Npart=Npart+1
	write(2,100) Npart,x0+r(1,i), y0+r(2,i), Wi, h0*1.6
      write(4,*) x0+r(1,i), y0+r(2,i)
	write(5,*) 0.0d0, 0.0d0
	end do

      end do

100    format(I5, 5E15.5)

	stop
	end

      subroutine MinGaussfit(n,r,HH0,hh)
	implicit real*8 (a-h,o-z)
	parameter (ndim=500)
	parameter (ndim2=ndim/2)
	dimension p(ndim), r(2,ndim2), pcom(ndim),xicom(ndim)
	common /fparam/ h0,h
	COMMON /f1com/ pcom,xicom,ncom
      npam=2*n
	h0=HH0
	h=hh
      ncom=n*2
      ftol=1.d-08
	iter=150
	l=0
      do i=1,n
	cost=2.*rand()-1.
	rad=0.5*h0*rand()
	l=l+1
	r(1,i)=rad*cost
	p(l)=r(1,i)
	l=l+1
	r(2,i)=rad*dsqrt(1.d0-cost**2)
	p(l)=r(2,i)

	end do
	fret0=1000.
1     call frprmn(p,npam,ftol,iter,fret)
      write(*,*) fret
	dfret=fret-fret0
	fret0=fret
      if(dabs(dfret).gt.0.00000*dabs(fret)) go to 1
	l=0
      do i=1,n
	do k=1,2
	l=l+1
	r(k,i)=p(l)
      end do
	end do
      return
	end 	

      double precision function func(p)
	implicit real*8 (A-H, O-Z)
	parameter (ndim=500)
	parameter (ndim2=ndim/2)
	dimension p(ndim), r(2,ndim2)
	dimension pcom(ndim), xicom(ndim)
	common /fparam/ h0,h
	COMMON /f1com/ pcom,xicom,ncom
	n=ncom/2
      l=0
      do i=1,n
      l=l+1
	r(1,i)=p(l)
	l=l+1
	r(2,i)=p(l)
	end do
      call SumGauss2D(N,r,h0,h, fdist2)
      func=fdist2+10.
	return
	end

      subroutine dfunc(p,xi)
	implicit real*8 (A-H, O-Z)
	parameter (ndim=500)
	parameter (ndim2=ndim/2)
	dimension p(ndim), xi(ndim),r(2,ndim2), grad(2,ndim2)
	dimension pcom(ndim), xicom(ndim)
	common /fparam/ h0,h
	COMMON /f1com/ pcom,xicom,ncom
	n=ncom/2

      l=0
      do i=1,n
      l=l+1
	r(1,i)=p(l)
	l=l+1
	r(2,i)=p(l)
	end do
      call Grad_SumGauss2D(N,r,h0,h, grad)
      l=0
      do i=1,n
      l=l+1
	xi(l)=grad(1,i)
	l=l+1
	xi(l)=grad(2,i)
	end do
   	return
	end

      subroutine SumGauss2D(N,r,h0,h, fdist2)
      implicit real*8 (a-h,o-z)
      parameter (ndim=500)
	parameter (ndim2=ndim/2)
	dimension r(2,ndim2), gzi(2)

      data pi/3.14159265358979d0/

      hhat=dsqrt(h0**2+h**2)
	hd=dsqrt(2.d0)*h

	A=0.d0
	B=0.d0

      do i=1,N

	do j=i,N

	do k=1,2
	gzi(k)=r(k,i)-r(k,j)
	end do

	A=A+Gauss2D(gzi,hd)
	end do                    ! end do j

      hanare=0.d0
	do k=1,2
	gzi(k)=r(k,i)
	hanare=hanare+r(k,i)**2
	end do               

      B=B+Gauss2D(gzi,hhat)

      adition=0.d0
	if(hanare.gt.8.d0*h0**2) then
      adition=N*((hanare-8.d0*h0**2))**2
      end if

      end do                  !   end do i
 
      fdist2=A/N**2-2.d0*B/N+adition

      return
	end

      subroutine Grad_SumGauss2D(N,r,h0,h, grad)

      implicit real*8 (a-h,o-z)
      parameter (ndim=500)
	parameter (ndim2=ndim/2)
	dimension r(2,ndim2), gzi(2), grad(2,ndim2)
	dimension  Agrad(2), Bgrad(2), Adgrad(2)
      dimension adition(2)
      data pi/3.14159265358979d0/

      hhat=dsqrt(h0**2+h**2)
	hd=dsqrt(2.d0)*h

      do i=1,N

      do k=1,2
	Agrad(k)=0.d0
	Bgrad(k)=0.d0
	end do

	do j=1,N

		do k=1,2
		gzi(k)=r(k,i)-r(k,j)
		end do

	      do k=1,2
	Agrad(k)=Agrad(k)-2.d0*Gauss2D(gzi,hd)*gzi(k)/hd**2
	      end do

      end do           ! end do of j


      hanare=0.d0
	do k=1,2
	gzi(k)=r(k,i)
	hanare=hanare+r(k,i)**2
	end do  
	             
      do k=1,2   
      Bgrad(k)=Bgrad(k)-2.d0*Gauss2D(gzi,hhat)*gzi(k)/hhat**2

      adition(k)=0.d0
	if(hanare.gt.8.d0*h0**2) then
      adition(k)=4.d0*10.d0*((hanare-8.d0*h0**2))*gzi(k)
      end if
      end do

      do k=1,2
      grad(k,i)=2.d0*(Agrad(k)/N**2-Bgrad(k)/N)+adition(k)
	end do

      end do                  !   end do i
 
      return
	end


	double precision function Gauss2D(gzi,hg)
	implicit real*8(a-h,o-z)
      dimension gzi(2)
      data pi/3.14159265358979d0/

      Gauss2D=0.d0
	q2=(gzi(1)**2+gzi(2)**2)/hg**2
      if(q2.gt.50.d0) return
	Gauss2D=dexp(-q2)/(pi*hg**2) 
	return
	end



c     Conjugate Gradient Method with Derivatives

      SUBROUTINE frprmn(p,n,ftol,iter,fret)
      INTEGER iter,n,NMAX,ITMAX
      REAL*8 fret,ftol,p(n),EPS,func
      EXTERNAL func
      PARAMETER (NMAX=500,ITMAX=4*Nmax,EPS=1.e-10)
CU    USES dfunc,func,linmin
      INTEGER its,j
      REAL*8 dgg,fp,gam,gg,g(NMAX),h(NMAX),xi(NMAX)
      fp=func(p)
      call dfunc(p,xi)
      do 11 j=1,n
        g(j)=-xi(j)
        h(j)=g(j)
        xi(j)=h(j)
11    continue
      do 14 its=1,ITMAX
        iter=its
        call linmin(p,xi,n,fret)
        if(2.d0*dabs(fret-fp).le.ftol*(dabs(fret)+dabs(fp)+EPS))return
        fp=func(p)
c	write(*,*) fp
        call dfunc(p,xi)
        gg=0.d0
        dgg=0.d0
        do 12 j=1,n
          gg=gg+g(j)**2
C         dgg=dgg+xi(j)**2
          dgg=dgg+(xi(j)+g(j))*xi(j)
12      continue
        if(gg.eq.0.d0)return
        gam=dgg/gg
        do 13 j=1,n
          g(j)=-xi(j)
          h(j)=g(j)+gam*h(j)
          xi(j)=h(j)
13      continue
14    continue
c      pause 'frprmn maximum iterations exceeded'
      return
      END

      SUBROUTINE linmin(p,xi,n,fret)
      INTEGER n,NMAX
      REAL*8 fret,p(n),xi(n),TOL
      PARAMETER (NMAX=500,TOL=1.e-5)
CU    USES dbrent,f1dim,df1dim,mnbrak
      INTEGER j,ncom
      REAL*8 ax,bx,fa,fb,fx,xmin,xx,pcom(NMAX),xicom(NMAX),brent
      COMMON /f1com/ pcom,xicom,ncom
      EXTERNAL f1dim, df1dim
      ncom=n
      do 11 j=1,n
        pcom(j)=p(j)
        xicom(j)=xi(j)
11    continue
      ax=0.d0
      xx=1.d0
      call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
      fret=dbrent(ax,xx,bx,f1dim,df1dim,TOL,xmin)
      do 12 j=1,n
        xi(j)=xmin*xi(j)
        p(j)=p(j)+xi(j)
12    continue
      return
      END

      SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,func)
      REAL*8 ax,bx,cx,fa,fb,fc,func,GOLD,GLIMIT,TINY
      EXTERNAL func
      PARAMETER (GOLD=1.618034, GLIMIT=100., TINY=1.e-20)
      REAL*8 dum,fu,q,r,u,ulim
      fa=func(ax)
      fb=func(bx)
      if(fb.gt.fa)then
        dum=ax
        ax=bx
        bx=dum
        dum=fb
        fb=fa
        fa=dum
      endif
      cx=bx+GOLD*(bx-ax)
      fc=func(cx)
1     if(fb.ge.fc)then
        r=(bx-ax)*(fb-fc)
        q=(bx-cx)*(fb-fa)
        u=bx-((bx-cx)*q-(bx-ax)*r)/(2.d0*dsign(max(dabs(q-r),TINY),q-r))
        ulim=bx+GLIMIT*(cx-bx)
        if((bx-u)*(u-cx).gt.0.d0)then
          fu=func(u)
          if(fu.lt.fc)then
            ax=bx
            fa=fb
            bx=u
            fb=fu
            return
          else if(fu.gt.fb)then
            cx=u
            fc=fu
            return
          endif
          u=cx+GOLD*(cx-bx)
          fu=func(u)
        else if((cx-u)*(u-ulim).gt.0.)then
          fu=func(u)
          if(fu.lt.fc)then
            bx=cx
            cx=u
            u=cx+GOLD*(cx-bx)
            fb=fc
            fc=fu
            fu=func(u)
          endif
        else if((u-ulim)*(ulim-cx).ge.0.)then
          u=ulim
          fu=func(u)
        else
          u=cx+GOLD*(cx-bx)
          fu=func(u)
        endif
        ax=bx
        bx=cx
        cx=u
        fa=fb
        fb=fc
        fc=fu
        goto 1
      endif
      return
      END

      FUNCTION dbrent(ax,bx,cx,f,df,tol,xmin)
      INTEGER ITMAX
      REAL*8 dbrent,ax,bx,cx,tol,xmin,df,f,ZEPS
      EXTERNAL df,f
      PARAMETER (ITMAX=500,ZEPS=1.0e-10)
      INTEGER iter
      REAL*8 a,b,d,d1,d2,du,dv,dw,dx,e,fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,
     *v,w,x,xm
      LOGICAL ok1,ok2
      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=0.d0
      fx=f(x)
      fv=fx
      fw=fx
      dx=df(x)
      dv=dx
      dw=dx
      do 11 iter=1,ITMAX
        xm=0.5d0*(a+b)
        tol1=tol*abs(x)+ZEPS
        tol2=2.d0*tol1
        if(dabs(x-xm).le.(tol2-.5*(b-a))) goto 3
        if(dabs(e).gt.tol1) then
          d1=2.d0*(b-a)
          d2=d1
          if(dw.ne.dx) d1=(w-x)*dx/(dx-dw)
          if(dv.ne.dx) d2=(v-x)*dx/(dx-dv)
          u1=x+d1
          u2=x+d2
          ok1=((a-u1)*(u1-b).gt.0.).and.(dx*d1.le.0.)
          ok2=((a-u2)*(u2-b).gt.0.).and.(dx*d2.le.0.)
          olde=e
          e=d
          if(.not.(ok1.or.ok2))then
            goto 1
          else if (ok1.and.ok2)then
            if(abs(d1).lt.abs(d2))then
              d=d1
            else
              d=d2
            endif
          else if (ok1)then
            d=d1
          else
            d=d2
          endif
          if(abs(d).gt.abs(0.5*olde))goto 1
          u=x+d
          if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
          goto 2
        endif
1       if(dx.ge.0.) then
          e=a-x
        else
          e=b-x
        endif
        d=0.5d0*e
2       if(abs(d).ge.tol1) then
          u=x+d
          fu=f(u)
        else
          u=x+sign(tol1,d)
          fu=f(u)
          if(fu.gt.fx)goto 3
        endif
        du=df(u)
        if(fu.le.fx) then
          if(u.ge.x) then
            a=x
          else
            b=x
          endif
          v=w
          fv=fw
          dv=dw
          w=x
          fw=fx
          dw=dx
          x=u
          fx=fu
          dx=du
        else
          if(u.lt.x) then
            a=u
          else
            b=u
          endif
          if(fu.le.fw .or. w.eq.x) then
            v=w
            fv=fw
            dv=dw
            w=u
            fw=fu
            dw=du
          else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
            v=u
            fv=fu
            dv=du
          endif
        endif
11    continue
      pause 'dbrent exceeded maximum iterations'
3     xmin=x
      dbrent=fx
      return
      END

      FUNCTION f1dim(x)
      INTEGER NMAX
      REAL*8 f1dim,func,x
      PARAMETER (NMAX=500)
CU    USES func
      INTEGER j,ncom
      REAL*8 pcom(NMAX),xicom(NMAX),xt(NMAX)
      COMMON /f1com/ pcom,xicom,ncom
      do 11 j=1,ncom
        xt(j)=pcom(j)+x*xicom(j)
11    continue
      f1dim=func(xt)
      return
      END

      FUNCTION df1dim(x)
      INTEGER NMAX
      REAL*8 df1dim,x
      PARAMETER (NMAX=500)
CU    USES dfunc
      INTEGER j,ncom
      REAL*8 df(NMAX),pcom(NMAX),xicom(NMAX),xt(NMAX)
      COMMON /f1com/ pcom,xicom,ncom
      do 11 j=1,ncom
        xt(j)=pcom(j)+x*xicom(j)
11    continue
      call dfunc(xt,df)
      df1dim=0.
      do 12 j=1,ncom
        df1dim=df1dim+df(j)*xicom(j)
12    continue
      return
      END


