      implicit real*8 (a-h,o-z)

c       2Dimensional case

c       Fit of a general function  by N Gaussian superposition on 
c       grid, first by ajusting the weight, then decompose them into
C       smaller, equal (approximately) weighted gaussian distribution.
C       Results are written into files, 2Dfit.dat
c       first line,                  Ntot 
c       then follows lines with      x(i), y(i), h(i), w(i)
c        
C       Before the use, fix the number of grid N1 in the parameter statement
C       In calling subroutine, r, H and W should be appropriately declared dimension.
C       On return, NSPH is the number of effective gaussians, and r, H and W are their parameters.

	parameter(N1=20, ndim=3)
	parameter (N=N1*N1*N1)
	dimension r(ndim,N),gzi(ndim),Xmax(ndim), Xmin(ndim),dx(ndim)
      dimension W(N)
      dimension gzi0(ndim)

      open(1, file='Firstfit_X_Y_W.dat')
	open(3, file='first_plots.dat')


      call FitGaussF_2DSimple(r,H,W,NSPH)

      wmax=0.d0
      wmin=1.d+50
	hmin=1.d+50
	do i=1,NSPH
	if(W(i).lt.wmin) wmin=W(i)
	if(W(i).gt.wmax) wmax=W(i)
	end do
      write(1,*) NSPH
	write(1, 100) H, Wmax,Wmin
      write(*,*) NSPH, H
	write(*, 100) H, Wmax/Wmin
      do i=1,NSPH
	write(1,100) r(1,i),r(2,i),r(3,i),W(i)
      end do
    
100    format(5E15.5)

      do k=1,ndim
	gzi(k)=0.d0
	end do
      fdum=funct(gzi,Xmax,Xmin,ndim) 
  
      dx(1)=(Xmax(1)-Xmin(1))/N1
      dx(2)=(Xmax(2)-Xmin(2))/N1
	dx(3)=(Xmax(3)-Xmin(3))/N1
      write(*,*) dx(1), dx(2), dx(3),h

      do i=1,3*N1
	do l=1,3*N1
	do m=1,3*N1
      gzi0(1)=Xmin(1)+(i-0.5d0)*dx(1)/3.d0
      gzi0(2)=Xmin(2)+(l-0.5d0)*dx(2)/3.d0
	gzi0(3)=Xmin(3)+(m-0.5d0)*dx(3)/3.d0
      ft=funct(gzi0,Xmax,Xmin,ndim)
      fcal=0.d0

	do j=1,NSPH
	do k=1,ndim
      gzi(k)=gzi0(k)-r(k,j)
	end do
	fcal=fcal+Gauss(gzi,H,ndim)*W(j)
	end do

	write(3,100) gzi0(1),gzi0(2), gzi0(3), ft, fcal
      end do
	end do
	end do

200   format(I10, 4E15.6)
	stop
	end

      subroutine FitGaussF_2DSimple(rSPH,H,WSPH,NSPH)
      implicit real*8 (a-h,o-z)
c       Fit of a general function  by N Gaussian superposition on 
c       grid, by ajusting the weight.
C       Before the use, fix the number of grid N1 in the parameter statement
C       In calling subroutine, r, H and W should be appropriately declared dimension.
C       On return, NSPH is the number of effective gaussians, and r, H and W are their parameters.

	parameter(N1=20, ndim=3)
	parameter (N=N1*N1*N1)

     	dimension rSPH(ndim,N),WSPH(N)
     
	dimension r(ndim,N),gzi(ndim),Xmax(ndim), Xmin(ndim),dx(ndim)
      dimension F(n), W(N) 
      dimension gzi0(ndim)

      Call Integral(Xmax, Xmin, S, fmax,ndim)
  
      dx(1)=(Xmax(1)-Xmin(1))/N1
      dx(2)=(Xmax(2)-Xmin(2))/N1
	dx(3)=(Xmax(3)-Xmin(3))/N1
      H=dx(1)
	if(H.lt.dx(2)) H=dx(2)
	if(H.lt.dx(3)) H=dx(3)
	H=H*dsqrt(2.d0)
      i=0
      do ix=1,N1
	do iy=1,N1
	i=i+1
      r(1,i)=Xmin(1)+(ix-0.5d0)*dx(1)
      r(2,i)=Xmin(2)+(iy-0.5d0)*dx(2)
	r(3,i)=Xmin(3)+(iz-0.5d0)*dx(3)
      
	do k=1,ndim
	gzi(k)=r(k,i)
	end do

	f(i)=funct(gzi,Xmax,Xmin,ndim)

      end do
      end do

      sum=0.d0
      do i=1,N
	sum=sum+F(i)
      end do
      do i=1,N
	W(i)=S*f(i)/sum
	end do
      write(*,*) sum, S

c     reduction process
c     discard those which contributes less than 0.5%

      NSPH=0
      Wmax=0.
      do i=1,N
	WSPH(i)=0.d0
      if(W(i).gt.Wmax) Wmax=W(i)
	end do

	do i=1,N
	if(W(i).gt.0.05d0*Wmax) then

	NSPH=NSPH+1
      WSPH(NSPH)=W(i)
	do k=1,ndim
	RSPH(k,NSPH)=r(k,i)
	end do
	end if

	end do
      sum=0.d0
	do i=1,NSPH
      sum=sum+WSPH(i)
	end do

	do i=1,NSPH
	WSPH(i)=WSPH(i)*S/sum
	end do

      return

	end


      subroutine Integral(Xmax,Xmin,S, fmax, ndim)
	implicit real*8(a-h, o-z)
      dimension gzi(2), dx(2), Xmax(ndim), Xmin(ndim)

      N=100
	fmax=0.d0
      do k=1,ndim
	gzi(k)=0.d0
	end do
      fdum= funct(gzi,Xmax,Xmin,ndim) 
      if(fdum.gt.fmax) fmax=fdum
      do k=1,ndim
	dx(k)=(Xmax(k)-Xmin(k))/(2*N)
	end do

	S=0.d0
	gzi(1)=Xmin(1)
      do i=1,N

	gzi(1)=gzi(1)+dx(1)
      sy=0.d0
	gzi(2)=Xmin(2)

      do j=1,N
	gzi(2)=gzi(2)+dx(2)

c      sz=0.0d0
c	gzi(3)=Xmin(3)
c      do k=1,N
c	gzi(3)=gzi(3)+dx(3)

	fdum=funct(gzi, Xmax, Xmin, ndim)
      if(fdum.gt.fmax) fmax=fdum
	sy=sy+4.d0*fdum
	gzi(2)=gzi(2)+dx(2)
	fdum=funct(gzi, Xmax, Xmin, ndim)
      if(fdum.gt.fmax) fmax=fdum
	sy=sy+2.d0*fdum

      end do
	S=S+4.d0*sy

	gzi(1)=gzi(1)+dx(1)
      sy=0.d0
	gzi(2)=Xmin(2)
      do j=1,N
	gzi(2)=gzi(2)+dx(2)	
	fdum=funct(gzi, Xmax, Xmin, ndim)
      if(fdum.gt.fmax) fmax=fdum
	sy=sy+4.d0*fdum
	gzi(2)=gzi(2)+dx(2)
	fdum=funct(gzi, Xmax, Xmin, ndim)
      if(fdum.gt.fmax) fmax=fdum
	sy=sy+2.d0*fdum
      end do

	S=S+2.d0*sy
      end do


	S=S*dx(1)*dx(2)/9.d0
	return
	end

	double precision function Gauss(gzi,hg,ndim)
	implicit real*8(a-h,o-z)
      dimension gzi(ndim)
      data pi/3.14159265358979d0/

      Gauss=0.d0
	q2=0.d0
	do k=1,ndim
	q2=q2+gzi(k)**2
	end do
	q2=q2/hg**2
      if(q2.gt.50.d0) return
	Gauss=dexp(-q2)/dsqrt(pi*hg**2)**ndim
	return
	end

      double precision function funct(gzi,Xmax,Xmin,ndim)
	implicit real*8 (a-h,o-z)
	dimension gzi(ndim),Xmax(ndim),Xmin(ndim)
	data a/0.35d0/,r0/4.d0/

	Xmin(1)=-3.d0
	Xmax(1)=3.d0
	Xmin(2)=-2.d0
	Xmax(2)=2.d0
	Xmin(3)=-1.0d0
	Xmax(3)=1.0d0

	funct=0.d0

	do k=1,ndim
	if(gzi(k).gt.Xmax(k) .or. gzi(k).lt.Xmin(k)) return
	end do

c	r2=gzi(1)**2+(gzi(2)*Xmax(1)/Xmax(2))**2

c	r=dsqrt(r2)
c	funct=1.d0/(1+dexp((r-r0)/a))
      funct=exp(-15.d0*((gzi(1)/Xmax(1))**2+(gzi(2)/Xmax(2))**2
     &    	+(gzi(3)/Xmax(3))**2))
	return
	end
