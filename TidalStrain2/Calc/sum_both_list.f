c
c	**************************************************
c	*             ----- sum_both_list -----          *
c       *                                                *
c       *                         by fuyuki 2015/02/17   *
c       *                         by fhirose 2015/07/16: time   *
c       *                                                *
c	**************************************************
c

c--
	double precision sExx,sEyy,sEzz,sExy,sExz,sEyz
!	double precision time,oExx,oEyy,oEzz,oExy,oExz,oEyz
	double precision oExx,oEyy,oEzz,oExy,oExz,oEyz
	double precision Exx,Eyy,Ezz,Exy,Exz,Eyz
	integer in
	character*18     UTC
	double precision sumtime,tinterval

c-- read parameter data
	open(10,file='./solid_earthtide.output9',
     &       status='old')
	open(11,file='./ocean_loading.out.2',
     &       status='old')
	open(20,file='./solid-ocean.out',
     &       status='unknown')
c--
	sumtime=0.0
!	tinterval=15.0  ! [minutes]
	tinterval=60.0
!	do i=1,97  ! 1day
!	do i=1,193  ! 2days
!	do i=1,289  ! 3days
!	do i=1,961  ! 2days, dt=3min
!	do i=1,131761
!	do i=1,35137  ! 366days (~1 year)
	do i=1,1405441  ! 14640days (~40 year)
	read(10,100) in,sExx,sEyy,sEzz,sExy,sExz,sEyz
!	read(11,111) time,oExx,oEyy,oEzz,oExy,oExz,oEyz
	read(11,111) UTC,oExx,oEyy,oEzz,oExy,oExz,oEyz

	Exx=sExx+oExx
	Eyy=sEyy+oEyy
	Ezz=sEzz+oEzz
	Exy=sExy+oExy
	Exz=sExz+oExz
	Eyz=sEyz+oEyz

	write(20,112) UTC,Exx,Eyy,Ezz,Exy,Exz,Eyz,sumtime
	sumtime=sumtime+tinterval

	end do

 100	format(i18,e12.4,e12.4,e12.4,e12.4,e12.4,e12.4)
 111	format(a18,e12.4,e12.4,e12.4,e12.4,e12.4,e12.4)
 112	format(a18,e12.4,e12.4,e12.4,e12.4,e12.4,e12.4,1x,f10.1)
c--
	close(10)
	close(11)
	close(20)

	end
