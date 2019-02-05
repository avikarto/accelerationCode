module ODEModule
    use mathModule
    use fieldsModule
    use IOModule
    implicit none

    contains

      SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,hmin,errmax,ntraj) ! Integrate the starting values ystart from x1 to x2 with accuracy eps
      implicit none

      INTEGER nbad,nok,nvar,NMAX,ntraj
      real(kind=8) eps,h1,hmin,x1,x2,ystart(nvar),TINY,rhoConvergence,rho
!      EXTERNAL derivs,bsstep
      PARAMETER (NMAX=6,TINY=1.d-30)
      INTEGER i,nstp
      real(kind=8) h,hdid,hnext,odeX,dydx(NMAX),odeY(NMAX),yscal(NMAX),errmax

      odeX=x1
      h=sign(h1,x2-x1)
      nok=0
      nbad=0
      nstp=0

      do 11 i=1,nvar
        odeY(i)=ystart(i)
11    enddo

      do while ( (odeX-x2)*(x2-x1)<0.d0 )

        rhoConvergence=( (1.d0+odeY(3)**2/zr**2)**(1.5d0) * w0**6*k**2 )**0.25d0
        rho=sqrt( odeY(1)**2 + odeY(2)**2 )
        if (rho.ge.rhoConvergence) then
            ystart(3)=-1.d5 ! sets z<0 such that the particle is filtered out of the main results
            !write(*,*) -1 ! for counting removed particles
            goto 99 ! ends trajectory calcualtion and returns results
        endif

        nstp = nstp + 1

       call derivs(odeX,odeY,dydx)
        do 12 i=1,nvar
          yscal(i)=abs(odeY(i))+abs(h*dydx(i))+TINY
12      enddo

        if((odeX+h-x2)*(odeX+h-x1)>0.d0) h=x2-odeX
        call bsstep(odeY,dydx,nvar,odeX,h,eps,yscal,hdid,hnext,errmax)

          if(nctt==3) then
            do i=1,nvar
              ystart(i)=odeY(i)
            enddo
            return
          endif

        if(hdid==h)then
          nok=nok+1
        else
          nbad=nbad+1
        endif

        if(abs(hnext)<hmin) then
        write(*,*) 'stepsize smaller than minimum in odeint'
          nctt = 4
          return
        endif

        h=hnext

        if (ntraj==1) then
            call makeFields(odeY(1),odeY(2),odeY(3),odeX)
            write(91,*) Real(fields(1))
            write(92,*) Real(fields(2))
            write(93,*) Real(fields(3))
            write(95,*) Real(fields(4))
            write(96,*) Real(fields(5))
            write(97,*) Real(fields(6))
            write(61,*) odeY(1)
            write(62,*) odeY(2)
            write(63,*) odeY(3)
            write(65,*) odeY(4)
            write(66,*) odeY(5)
            write(67,*) odeY(6)
            write(64,*) odeX
            write(71,*) sqrt( 1.d0 + (odeY(4)**2+odeY(5)**2+odeY(6)**2)/c2 )*c2*eVperHartree
        !elseif (ntraj.gt.1) then
            ! p*dr is infinitesimal phase accumulation
        !    distancePhase = distancePhase + sqrt( (odeY(1)-oldPosition(1))**2 + (odeY(2)-oldPosition(2))**2 + (odeY(3)-oldPosition(3))**2 ) * sqrt(odeY(4)**2+odeY(5)**2+odeY(6)**2)
        endif
      enddo ! while ( (odeX-x2)*(x2-x1)<0.d0 )

      do 14 i=1,nvar
        ystart(i)=odeY(i)
14    enddo

99    continue
      return
    END

! adapted from Numerical Recipes, chapter 16
      SUBROUTINE bsstep(bssY,dydx,nv,bssX,htry,eps,yscal,hdid,hnext,errmax)
      implicit none
      INTEGER nv,NMAX,KMAXX,IMAX
      REAL(kind=8) eps,hdid,hnext,htry,bssX,dydx(nv),bssY(nv),yscal(nv),SAFE1,&
        SAFE2,REDMAX,REDMIN,TINY,SCALMX
      PARAMETER (NMAX=6,KMAXX=8,IMAX=KMAXX+1,SAFE1=.25d0,SAFE2=.7d0,&
        REDMAX=1.d-5,REDMIN=.7d0,TINY=1.d-30,SCALMX=.1d0) !------------------SCALMX=.1d0 (my .9d0)------------------------
!CU    USES mmid,pzextr
!Bulirsch-Stoer step with monitoring of local truncation error to ensure accuracy and adjust
!stepsize. Input are the dependent variable vector y(1:nv) and its derivative dydx(1:nv)
!at the starting value of the independent variable x. Also input are the stepsize to be attempted
!htry, the required accuracy eps, and the vector yscal(1:nv) against which the
!error is scaled. On output, y and x are replaced by their new values, hdid is the stepsize
!that was actually accomplished, and hnext is the estimated next stepsize. derivs is the
!user-supplied subroutine that computes the right-hand side derivatives. Be sure to set htry
!on successive steps to the value of hnext returned from the previous step, as is the case
!if the routine is called by odeint.

!Parameters: NMAX is the maximum value of nv; KMAXX is the maximum row number used
!in the extrapolation; IMAX is the next row number; SAFE1 and SAFE2 are safety factors;
!REDMAX is the maximum factor used when a stepsize is reduced, REDMIN the minimum;
!TINY prevents division by zero; 1/SCALMX is the maximum factor by which a stepsize can
!be increased.

      INTEGER i,iq,k,kk,km,kmax,kopt,nseq(IMAX)
      REAL(kind=8) eps1,epsold,errmax,fact,h,red,scale,work,wrkmin,xest,&
        xnew,a(IMAX),alf(KMAXX,KMAXX),err(KMAXX),yerr(NMAX),ysav(NMAX),&
        yseq(NMAX)
      LOGICAL first,reduct
      SAVE a,alf,epsold,first,kmax,kopt,nseq,xnew
      DATA first/.true./,epsold/-1.d0/
      DATA nseq /2,4,6,8,10,12,14,16,18/
      if(eps.ne.epsold)then
        hnext=-1.d29
        xnew=-1.d29
        eps1=SAFE1*eps
        a(1)=nseq(1)+1
        do 11 k=1,KMAXX
          a(k+1)=a(k)+nseq(k+1)
11      continue
        do 13 iq=2,KMAXX
          do 12 k=1,iq-1
            alf(k,iq)=eps1**((a(k+1)-a(iq+1))/((a(iq+1)-a(1)+1.)*(2*k+1)))
12        continue
13      continue
        epsold=eps
        do 14 kopt=2,KMAXX-1
          if(a(kopt+1).gt.a(kopt)*alf(kopt-1,kopt))goto 1
14      continue
1       kmax=kopt
      endif
      h=htry
      do 15 i=1,nv
        ysav(i)=bssY(i)
15    continue
      if(h.ne.hnext.or.bssX.ne.xnew)then
        first=.true.
        kopt=kmax
      endif
      reduct=.false.
2     do 17 k=1,kmax
        xnew=bssX+h

!        if(xnew.eq.bssX)pause 'step size underflow in bsstep'
        if(xnew==bssX) then
           write(*,*)      'step size underflow in bsstep'
           write(*,*) "failed step size: h=",h
          nctt = 3
          return
        endif

        call mmid(ysav,dydx,nv,bssX,h,nseq(k),yseq)
        xest=(h/nseq(k))**2
        call pzextr(k,xest,yseq,bssY,yerr,nv)
        if(k.ne.1)then
          errmax=TINY
          do 16 i=1,nv
            errmax=max(errmax,abs(yerr(i)/yscal(i)))
16        continue
          errmax=errmax/eps
          km=k-1
          err(km)=(errmax/SAFE1)**(1.d0/(2*km+1))
        endif
        if(k.ne.1.and.(k.ge.kopt-1.or.first))then
          if(errmax.lt.1.)goto 4
          if(k.eq.kmax.or.k.eq.kopt+1)then
            red=SAFE2/err(km)
            goto 3
          else if(k.eq.kopt)then
            if(alf(kopt-1,kopt).lt.err(km))then
              red=1.d0/err(km)
              goto 3
            endif
          else if(kopt.eq.kmax)then
            if(alf(km,kmax-1).lt.err(km))then
              red=alf(km,kmax-1)*SAFE2/err(km)
              goto 3
            endif
          else if(alf(km,kopt).lt.err(km))then
            red=alf(km,kopt-1)/err(km)
            goto 3
          endif
        endif
17    continue
3     red=min(red,REDMIN)
      red=max(red,REDMAX)
      h=h*red
      reduct=.true.
      goto 2
4     bssX=xnew
      hdid=h
      first=.false.
      wrkmin=1.d35
      do 18 kk=1,km
        fact=max(err(kk),SCALMX)
        work=fact*a(kk+1)
        if(work.lt.wrkmin)then
          scale=fact
          wrkmin=work
          kopt=kk+1
        endif
18    continue
      hnext=h/scale
      if(kopt.ge.k.and.kopt.ne.kmax.and..not.reduct)then
        fact=max(scale/alf(kopt-1,kopt),SCALMX)
        if(a(kopt+1)*fact.le.wrkmin)then
          hnext=h/fact
          kopt=kopt+1
        endif
      endif
      return
      END

      SUBROUTINE mmid(mmY,dydx,nvar,xs,htot,nstep,yout)
      implicit none
      INTEGER nstep,nvar,NMAX
      REAL(kind=8) htot,xs,dydx(nvar),mmY(nvar),yout(nvar)
!      EXTERNAL derivs
      PARAMETER (NMAX=6)
      INTEGER i,n
      REAL(kind=8) h,h2,swap,mmX,ym(NMAX),yn(NMAX)
      h=htot/nstep
      do 11 i=1,nvar
        ym(i)=mmY(i)
        yn(i)=mmY(i)+h*dydx(i)
11    continue
      mmX=xs+h
      call derivs(mmX,yn,yout)
      h2=2.d0*h
      do 13 n=2,nstep
        do 12 i=1,nvar
          swap=ym(i)+h2*yout(i)
          ym(i)=yn(i)
          yn(i)=swap
12      continue
        mmX=mmX+h
        call derivs(mmX,yn,yout)
13    continue
      do 14 i=1,nvar
        yout(i)=0.5d0*(ym(i)+yn(i)+h*yout(i))
14    continue
      return
      END

      SUBROUTINE pzextr(iest,xest,yest,yz,dy,nv)
      implicit none
      INTEGER iest,nv,IMAX,NMAX
      REAL(kind=8) xest,dy(nv),yest(nv),yz(nv)
      PARAMETER (IMAX=13,NMAX=6)
      INTEGER j,k1
      REAL(kind=8) delta,f1,f2,q,d(NMAX),qcol(NMAX,IMAX),pzX(IMAX)
      SAVE qcol,pzX
      pzX(iest)=xest
      do 11 j=1,nv
        dy(j)=yest(j)
        yz(j)=yest(j)
11    continue
      if(iest.eq.1) then
        do 12 j=1,nv
          qcol(j,1)=yest(j)
12      continue
      else
        do 13 j=1,nv
          d(j)=yest(j)
13      continue
        do 15 k1=1,iest-1
          delta=1.d0/(pzX(iest-k1)-xest)
          f1=xest*delta
          f2=pzX(iest-k1)*delta
          do 14 j=1,nv
            q=qcol(j,k1)
            qcol(j,k1)=dy(j)
            delta=d(j)-q
            dy(j)=f1*delta
            d(j)=f2*delta
            yz(j)=yz(j)+dy(j)
14        continue
15      continue
        do 16 j=1,nv
          qcol(j,iest)=dy(j)
16      continue
      endif
      return
      END

    Subroutine derivs(myX,myY,dydx)
        IMPLICIT none

        REAL(kind=8) myX,myY(6),DYDX(6)
        real(kind=8) ex,ey,ez,bx,by,bz,gama,regama

        call makeFields(myY(1),myY(2),myY(3),myX)

        if(polar<1 .or. polar>4) then
            write(*,*) "Invalid Polarization.  See ODEModule."
            stop
        endif

        if (polar==2) then ! radial polarization.  fields={Er,Ep,Ez,Br,Bp,Bz}
          ex=real(fields(1))*cp-sp*real(fields(2))
          ey=real(fields(1))*sp+cp*real(fields(2))
          ez=real(fields(3))
          bx=real(fields(4))*cp-sp*real(fields(5))
          by=real(fields(4))*sp+cp*real(fields(5))
          bz=real(fields(6))
        else if (polar==1 .or. polar==3 .or. polar==4) then ! lin/circ/azim polarization.  fields={Ex,Ey,Ez,Bx,By,Bz}
          ex=real(fields(1))
          ey=real(fields(2))
          ez=real(fields(3))
          bx=real(fields(4))
          by=real(fields(5))
          bz=real(fields(6))
        endif

        gama =sqrt(1.d0+(myY(4)**2+myY(5)**2+myY(6)**2)/c2) ! sqrt(1+P^2/c^2)
        regama=1.d0/gama

        ! r-dot = v = p/(gamma*m) = p/gamma  (in a.u.)
        dydx(1) = myY(4)*regama ! x-dot = Px/gamma
        dydx(2) = myY(5)*regama
        dydx(3) = myY(6)*regama

        ! F = m*a = m*v-dot = p-dot
        dydx(4) = -ex -(myY(5)*bz-myY(6)*by)*regama ! Px-dot = F_x = -E_x-(VxB)_x = -E_x-((P/gamma)xB)_x = -E_x-(PyBz-PzBy)/gamma
        dydx(5) = -ey -(myY(6)*bx-myY(4)*bz)*regama
        dydx(6) = -ez -(myY(4)*by-myY(5)*bx)*regama

        return
    END subroutine

end module ODEModule
