      implicit real*8 (a-h,o-z)

c       2Dimensional case

	parameter(Ntot=100000, ndim=2)
	dimension r(ndim,Ntot),rgrid(ndim,Ntot), ft(Ntot)
      dimension	gzi(ndim)
      dimension W(Ntot)
      dimension gzi0(ndim)

      open(1, file='first_plots.dat')
	open(2, file='SPHparticles.dat')
	open(3, file='finalfit.dat')
	Nt=0
1     read(1,*,end=2) x0,y0,ff
      Nt=Nt+1
      rgrid(1,Nt)=x0
	rgrid(2,Nt)=y0
	ft(Nt)=ff
	go to 1

2     Np=0
3     read(2,*,end=4) NSPH, x0,y0,ww,h
      Np=Np+1
      r(1,Np)=x0
	r(2,Np)=y0
	w(Np)=ww
	go to 3
4     continue
    
100    format(5E15.5)

	do it=1,Nt
      write(*,*) '--- ', it, '  of  ',Nt
	fcal=0.d0

	do j=1,Np
	do k=1,ndim
      gzi(k)=rgrid(k,it)-r(k,j)
	end do
	fcal=fcal+Gauss(gzi,H,ndim)*W(j)
	end do

	write(3,100) rgrid(1,it),rgrid(2,it), ft(it), fcal
      end do

200   format(I10, 4E15.6)
	stop
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
