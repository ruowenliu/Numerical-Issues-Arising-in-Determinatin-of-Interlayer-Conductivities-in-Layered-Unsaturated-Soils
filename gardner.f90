! Richards' equation for 1D with two distinct soil layers
! Gardner Model
! Boundary Conditions: Dirichlet

program gardner
use dvode_f90_m
implicit none
character :: prefix*25, fileout*30
real(8), dimension(2) :: a, b, thr, ths ! parameters
integer :: neq, maxiter, initfile, ix, i, iint, k, iters, pl
integer :: ISTATE, IERROR, ITASK
integer :: numtime=0, num=0, multinum=0, pos_up_num=0, pos_low_num=0, dif_flux_num=0
integer :: p1, p2
real(8) :: cdint, wvlint(1), hvlint(1)
real(8), dimension(1) :: ht, hb, tht, thb
real(8) :: t0, tf, dt, rtol, atol, tol, tspan(8), hmx
real(8) :: dz, hu, hl, ku, kl, r=1.0d0, rsucc, lambda, mu, kminus, kplus, yd(3), qd(3)
real(8), allocatable :: y(:), h(:)
logical :: upperx_sat, lowerx_sat, multi_root, incons_fl
type(VODE_OPTS) :: OPTIONS
real(8) :: start=-1.0d0, ending=5.0d0
real(8) :: tstart1, tend2
!**********************************
read*, prefix
read*, a
read*, b
read*, thr
read*, ths
read*, neq
read*, ht
read*, hb
read*, tspan
read*, rtol
read*, atol
read*, maxiter
read*, tol
read*, initfile
ix=index(prefix//' ',' ')-1
dz = 1.0d0/dble(neq+1)
!h = (/(hb(1)-dble(neq-i+1)/dble(neq+1),i=1,neq)/)
cdint=(50d0/2d0+1d0)/51d0
p1 = floor(cdint/dz)
p2 = ceiling(cdint/dz)
pl = floor((0.05d0-0.5d0*dz)/dz)
allocate(h(neq))
h(neq/2+pl+1:neq) = hb(1)
!h(1:neq/2)=hb(1)-dz
do k=1,neq/2+pl
  h(neq/2+pl+1-k)=hb(1)-dble(k)*dz
end do
allocate(y(neq))
call THETA_H(1,1,ht,tht)
call THETA_H(2,1,hb,thb)
call THETA_H(1,neq/2,h(1:neq/2),y(1:neq/2))
call THETA_H(2,neq/2,h(neq/2+1:neq),y(neq/2+1:neq))
! write initial
fileout(1:ix+2)=prefix(1:ix)//'h.'
write(fileout(ix+3:ix+6),'(i4.4)') initfile
open(unit=10,file=fileout(1:ix+6),form='unformatted',access='stream')
write(10) ht(1)
do k=1,neq
  write(10) h(k)
end do
write(10) hb(1)
close(10)
fileout(1:ix+2)=prefix(1:ix)//'w.'
write(fileout(ix+3:ix+6),'(i4.4)') initfile
open(unit=10,file=fileout(1:ix+6),form='unformatted',access='stream')
write(10) tht(1)
do k=1,neq
  write(10) y(k)
end do
write(10) thb(1)
close(10)
! end write initial
open(unit=11,file=prefix(1:ix)//'sat.txt')
open(unit=21,file=prefix(1:ix)//'phu.bin',form='unformatted',access='stream')
open(unit=22,file=prefix(1:ix)//'phl.bin',form='unformatted',access='stream')
open(unit=23,file=prefix(1:ix)//'time.bin',form='unformatted',access='stream')
open(unit=24,file=prefix(1:ix)//'lam.bin',form='unformatted',access='stream')
open(unit=25,file=prefix(1:ix)//'mu.bin',form='unformatted',access='stream')
open(unit=26,file=prefix(1:ix)//'r.bin',form='unformatted',access='stream')
open(unit=28,file=prefix(1:ix)//'pwu.bin',form='unformatted',access='stream')
open(unit=29,file=prefix(1:ix)//'pwl.bin',form='unformatted',access='stream')
open(unit=30,file=prefix(1:ix)//'hbelow_neq50.bin',form='unformatted',access='stream')
open(unit=31,file=prefix(1:ix)//'wbelow_neq50.bin',form='unformatted',access='stream')
open(unit=32,file=prefix(1:ix)//'multitime.bin',form='unformatted',access='stream')
open(unit=41,file=prefix(1:ix)//'kminus.bin',form='unformatted',access='stream')
open(unit=42,file=prefix(1:ix)//'kplus.bin',form='unformatted',access='stream')
IERROR = 0
ITASK = 5
numtime = size(tspan)
call cpu_time(tstart1)
do iint = 1,numtime-1
  dt = tspan(iint+1)-tspan(iint)
  t0 = tspan(iint)
  tf = tspan(iint+1)
  write(*,*) 'Time from ', t0,' to ', tf
  ISTATE = 1  ! re-initialzation for the DVODE
  if (iint==1) then
     hmx=1d-8
  else
     hmx = (1d-3)*dt
  end if
  OPTIONS = SET_OPTS(METHOD_FLAG=22,RELERR=rtol,ABSERR=atol,TCRIT=tf,MXSTEP=100000,HMAX=hmx)
  do while (t0<tf)
    call DVODE_F90(dthetadt,neq,y,t0,tf,ITASK,ISTATE,OPTIONS)
!    if (t0<1e-7) then
!      hmx=1e-8
!      OPTIONS = SET_OPTS(METHOD_FLAG=22,RELERR=rtol,ABSERR=atol,TCRIT=tf,MXSTEP=100000,HMAX=hmx)
!      ISTATE = 3
!    else
!      hmx = 1d-3*dt
!    end if
    rsucc = r
    write(11,*) 'At time = ', t0
    write(11,*) 'upper > 0 ', upperx_sat
    write(11,*) 'lower > 0 ', lowerx_sat
    write(11,*) 'multiple roots ', multi_root
    write(11,*) 'inconsistent flux ', incons_fl
    write(11,*) rsucc, kminus, kplus, hu, hl
    write(11,*) yd
    write(11,*) qd
    write(21) hu
    write(22) hl
    write(23) t0
    write(24) lambda
    write(25) mu
    write(26) rsucc
    write(28) y(neq/2)
    write(29) y(neq/2+1)
    write(41) kminus
    write(42) kplus
    if (p1 == p2) then
      wvlint(1) = y(p1)
    else
      wvlint(1) = y(p1)+(y(p2)-y(p1))/dz*(cdint-p1*dz)
    end if
    call H_THETA(2,1,wvlint,hvlint) ! point below interface
    write(30) hvlint(1)
    write(31) wvlint(1)
    if (multi_root) then
      write(32) t0
      multinum = multinum+1
    end if
    if (upperx_sat) then
      pos_up_num = pos_up_num+1
    end if
    if (lowerx_sat) then
      pos_low_num = pos_low_num+1
    end if
    if (incons_fl) then
      dif_flux_num = dif_flux_num+1
    end if
    num = num+1
  end do
  write(*,*) 'After step ',iint,' number of integrations is ',num
  call cpu_time(tend2)
  write(*,*) 'Elapsed CPU time (in seconds) = ', tend2-tstart1
  write(*,*) '# of time steps = ', num
  write(*,*) '# of multiple = ', multinum
  write(*,*) '# of positive upper h = ', pos_up_num
  write(*,*) '# of positive lower h = ', pos_low_num
  write(*,*) '# of inconsistent flux = ', dif_flux_num

  ! **** record every step pressure heads and water contents ***
  initfile = initfile+1
  fileout(1:ix+2)=prefix(1:ix)//'w.'
  write(fileout(ix+3:ix+6),'(i4.4)') initfile
  open(unit=10,file=fileout(1:ix+6),form='unformatted',access='stream')
  write(10) tht(1)
  do k=1,neq
    write(10) y(k)
  end do
  write(10) thb(1)
  close(10)
  call H_THETA(1,neq/2,y(1:neq/2),h(1:neq/2))
  call H_THETA(2,neq/2,y(neq/2+1:neq),h(neq/2+1:neq))
  fileout(1:ix+2)=prefix(1:ix)//'h.'
  write(fileout(ix+3:ix+6),'(i4.4)') initfile
  open(unit=10,file=fileout(1:ix+6),form='unformatted',access='stream')
  write(10) ht(1)
  do k=1,neq
    write(10) h(k)
  end do
  write(10) hb(1)
  close(10)
  ! **** end record every step pressure heads and water contents ***
end do
!write(*,*) '# of time steps = ', num
!write(*,*) '# of multiple = ', multinum
!write(*,*) '# of positive upper h = ', pos_up_num
!write(*,*) '# of positive lower h = ', pos_low_num
!write(*,*) '# of inconsistent flux = ', dif_flux_num
close(11)
close(21)
close(22)
close(23)
close(24)
close(25)
close(26)
close(28)
close(29)
close(30)
close(31)
close(32)
close(41)
close(42)

contains  !******** pack all the subroutines and functions

subroutine DTHETADT(neq,t,y,ydot)
integer,intent(in) :: neq
real(8),intent(in) :: t, y(neq)
real(8),intent(out) :: ydot(neq)
real(8) :: theta(0:neq+1),hj(0:neq+1),kj(0:neq+1),kjm(0:neq),qjm(0:neq)
real(8),dimension(1) :: hux, hlx, kux, klx
real(8) :: dh, qminus, qplus
upperx_sat = .false.
lowerx_sat = .false.
multi_root = .false.
incons_fl = .false.
hj(0) = ht(1)
hj(neq+1) = hb(1)
theta(1:neq) = y
call H_THETA(1,neq/2,y(1:neq/2),hj(1:neq/2))
call H_THETA(2,neq/2,y(neq/2+1:neq),hj(neq/2+1:neq))
call K_H(1,neq/2+1,hj(0:neq/2),kj(0:neq/2))
call K_H(2,neq/2+1,hj(neq/2+1:neq+1),kj(neq/2+1:neq+1))
kjm = dsqrt(kj(0:neq)*kj(1:neq+1))
qjm = kjm*(1.0d0-(hj(1:neq+1)-hj(0:neq))/dz)

hu = hj(neq/2)
hl = hj(neq/2+1)
ku = kj(neq/2)
kl = kj(neq/2+1)

call LAMBDA_G(lambda)
call MU_G(mu)
if ((mu<=-2.0d0).AND.&
(dabs(dlog(lambda))<=dsqrt(mu*(2.0d0+mu))+dlog(dabs(1.0d0+mu+dsqrt(mu*(2.0d0+mu)))))) then
  multi_root = .true.
  write(*,*) 'Warning: multiple roots at time = ', T
end if

r=rsucc

call fixedpoint(G,r,iters)
if (iters >= maxiter) then
    r = rsucc
    call anewton(F,r,iters)
  if (iters >= maxiter) then
   write(*,*) 'Neither Fixed-point nor Newton method converges. Need to reduce the tolerance.'
  end if
end if

dh = (dz-(hl-hu))*(1.0d0-r)/(1.0d0+r)

hux(1) = hu-dh
hlx(1) = hl-dh
if (hux(1)>=0.0d0) then
  upperx_sat = .true.
end if
if (hlx(1)>=0.0d0) then
  lowerx_sat = .true.
end if
call K_H(1,1,hlx,klx)
call K_H(2,1,hux,kux)
kminus = dsqrt(ku*klx(1))
kplus = dsqrt(kl*kux(1))
qminus = kminus*(1.0d0-((hl-dh)-hu)/dz)
qplus = kplus*(1.0d0-(hl-(hu-dh))/dz)
if (dabs(qminus-qplus)>10.0d0*atol) then
  incons_fl = .true.
end if
qjm(neq/2) = (qminus+qplus)/2.0d0
ydot = -(qjm(1:neq)-qjm(0:neq-1))/dz
yd = ydot(neq/2-1:neq/2+1)
qd = qjm(neq/2-1:neq/2+1)
end subroutine DTHETADT

subroutine S_THETA(layer,dim,theta,s)
integer,intent(in) :: layer,dim
real(8),intent(in) :: theta(dim)
real(8),intent(out) :: s(dim)
s = dmax1(0.0d0,dmin1(1.0d0, (theta-thr(layer))/(ths(layer)-thr(layer))))
return
end subroutine S_THETA

subroutine THETA_S(layer,dim,s,theta)
integer,intent(in) :: layer,dim
real(8),intent(in) :: s(dim)
real(8),intent(out) :: theta(dim)
theta = dmax1(thr(layer),dmin1(ths(layer),s*(ths(layer)-thr(layer))+thr(layer)))
return
end subroutine THETA_S

subroutine S_H(layer,dim,h,s)
implicit none
integer,intent(in) :: layer,dim
real(8),intent(in) :: h(dim)
real(8),intent(out) :: s(dim)
s = dexp(a(layer)*dmin1(h,0.0d0))
return
end subroutine S_H

subroutine H_S(layer,dim,s,h)
integer,intent(in) :: layer,dim
real(8),intent(in) :: s(dim)
real(8),intent(out) :: h(dim)
h = dmin1(0.0d0,(-dlog(1.0d0/dmax1(s,0.0d0)))/a(layer))
return
end subroutine H_S

subroutine THETA_H(layer,dim,h,theta)
integer,intent(in) :: layer,dim
real(8),intent(in) :: h(dim)
real(8),intent(out) :: theta(dim)
real(8) :: s(dim)
call S_H(layer,dim,h,s)
call THETA_S(layer,dim,s,theta)
end subroutine THETA_H

subroutine H_THETA(layer,dim,theta,h)
integer,intent(in) :: layer,dim
real(8),intent(in) :: theta(dim)
real(8),intent(out) :: h(dim)
real(8) :: s(dim)
call S_THETA(layer,dim,theta,s)
call H_S(layer,dim,s,h)
end subroutine H_THETA

subroutine K_S(layer,dim,s,k)
integer, intent(in) :: layer,dim
real(8), intent(in) :: s(dim)
real(8), intent(out) :: k(dim)
k = dmax1(0.0d0,b(layer)*dmin1(s,1.0d0))
return
end subroutine K_S

subroutine K_H(layer,dim,h,k)
integer, intent(in) :: layer,dim
real(8), intent(in) :: h(dim)
real(8), intent(out) :: k(dim)
real(8) :: s(dim)
call S_H(layer,dim,h,s)
call K_S(layer,dim,s,k)
return
end subroutine K_H

subroutine fixedpoint(f,x,iters)
real(8),intent(inout) :: x
real(8),external :: f
integer,intent(out) :: iters
real(8) :: xold
xold = x
x = f(x)
iters = 1
do while ( ((dabs(x-xold)/dabs(xold) > tol).and.(dabs(x-xold) > tol)) .and. iters < maxiter)
  xold = x
  x = f(x)
  iters = iters+1
enddo
if (iters >= maxiter) then
  xold = x
  x = f(x)
!  if (DABS(x-xold) > tol) then
!    print*, '***Warning: Fixed-point iteration has not yet converged.'
!    print*, 'root, second last = '
!    print*, xold
!    print*, 'root, last = '
!    print*, x
!    print*
!  endif
endif
return
end subroutine fixedpoint

subroutine anewton(f,x,iters)
real(8), intent(inout) :: x
real(8), external :: f
integer, intent(out) :: iters
real(8) :: macheps, hstep, fx, xold
iters = 0
macheps = EPSILON(0d0)
fx = f(x)
if (DABS(fx) > 0.0D0) then
hstep = (macheps*DABS(fx))**(1.0D0/3.0D0)
xold = x
x = x-fx/((f(x+hstep)-f(x-hstep))/(2.0D0*hstep))
iters = iters+1
do while (((dabs(x-xold)/dabs(xold) > tol).and.(dabs(x-xold) > tol)).and.iters < maxiter)
fx = f(x)
hstep = (macheps*DABS(fx))**(1.0D0/3.0D0)
xold = x
x = x-fx/((f(x+hstep)-f(x-hstep))/(2.0D0*hstep))
iters = iters+1
end do
end if
if (iters >= maxiter) then
  fx = f(x)
!  if (DABS(fx) > tol) then
!    write (*,*) '***Warning: Newton iteration has not yet converged'
!    print*, 'root, second last = '
!    print*, xold
!    print*, 'root, last = '
!    print*, x
!    print*
!  endif
endif
return
end subroutine anewton

real(8) function G(x)
real(8),intent(in) :: x
real(8),dimension(1) :: hlex, huex, klex, kuex
hlex(1) = hl-(dz-(hl-hu))*(1.0d0-x)/(1.0d0+x)
huex(1) = hu-(dz-(hl-hu))*(1.0d0-x)/(1.0d0+x)
call K_H(1,1,hlex,klex)
call K_H(2,1,huex,kuex)
G = dsqrt((ku*klex(1))/(kl*kuex(1)))
return
end function G

real(8) function F(x)
real(8),intent(in) :: x
F = G(x)-x
return
end function F

subroutine LAMBDA_G(lam)
real(8),intent(out) :: lam
lam = b(1)/b(2)*dexp((a(1)-a(2))*(hl+hu)/2.0d0)
return
end subroutine LAMBDA_G

subroutine MU_G(muv)
real(8),intent(out) :: muv
muv = 0.5d0*(a(2)-a(1))*(dz-(hl-hu))
return
end subroutine MU_G

end program gardner