module fieldsModule
    use constants
    use IOModule
    use mathModule
    implicit none

contains

    subroutine checkIonization(ppp,t,ionized)
        implicit none
        real*8 ppp(1:3),t,rateT,rateMP,prob,rand,IPeV,keldysh,IPHartree,Emin
        integer nPhotonsK ! multiphoton variable
        real*8 sigma_nK ! multiphoton variable
        real*8 tmu,txi,txi2,txi3,teps! semiclassical rel. tunneling variables
        real*8 rateADK,rateBSI,nstar,lstar,Ec,potvar ! BSI variables
        logical ionized

        call makeFields(ppp(1),ppp(2),ppp(3),t)
        ETotal=sqrt( real(fields(1))**2 + real(fields(2))**2 + real(fields(3))**2 )

        if (electronNum==1) then ! from "Calculated Ionization Potentials for Multiply Charged Ions"
            if (znum==1) then
                IPeV=13.6 ! H
                Emin=0.d0
            elseif (znum==2) then
                IPeV=5.433d01 ! He+
                Emin=0.d0
            elseif (znum==10) then
                IPeV=1.296d03 ! Ne 9+
                Emin=0.d0
            elseif (znum==18) then
                IPeV=4.264d03 ! Ar 17+
                Emin=0.d0
            elseif (znum==19) then
                IPeV=4.759d03 ! K 18+
                Emin=0.d0
            elseif (znum==20) then
                IPeV=5.280d03 ! Ca 19+
                Emin=0.d0
            elseif (znum==21) then
                IPeV=5.827d03 ! Sc 20+
                Emin=0.d0
            elseif (znum==36) then
                IPeV=1.740d04 ! Kr 35+
                Emin=0.d0
            else
                write(*,*) "No given IPeV for chosen znum.  Manually add from resources."
                stop
            endif
        else
            write(*,*) "Error: given value of electronNum is not supported."
            stop
        endif

        if (ETotal.le.Emin) then ! skip calculations for small ionization probability
            ionized=.false.
            return
        endif

        rateT=0
        rateMP=0
        rateADK=0
        rateBSI=0

        IPHartree=IPeV/eVperHartree

        if (ionizationMode==0 .and. keldysh .le. 2.d0) then ! tunnel ionization regime (<0.5), and intermediate regime (0.5<keldysh<2)
            ! Joachain "Atoms in Intense Laser Fields", EQ (6.22).  OK for lower atomic number ions (znum <= 20 with maximum ~20% error)
            rateT=(4.d0*znum**5/ETotal)*exp(-2.d0*znum**3/(3.d0*ETotal))
        elseif (ionizationMode==1) then

            nstar=znum/sqrt(2.d0*IPHartree)
            lstar=nstar-1.d0

            if (omega0>IPHartree) then ! E = hbar*omega > IP, BSI regime
                ! Alistair Lawrence-Douglas "Ionization Effects for Laser-Plasma Interactions by Particle-in-Cell Code", Eq. 2.10
                potvar=(2.d0*IPHartree)**1.5d0
                Ec=znum**3/(16.d0*nstar**4)
                rateADK=IPhartree*2.d0**(2.d0*nstar)/(nstar*gamma(nstar+lstar+1)*gamma(nstar-lstar))&
                    *sqrt( 3.d0*Ec/(pi*potvar) )&
                    *(potvar/Ec)**(2.d0*nstar-1.d0)&
                    *exp(-2.d0*potvar/(3.d0*Ec))
                rateBSI=(2.d0*IPHartree)**(3.d0/2.d0)*(1.d0-IPhartree**2/(4.d0*znum*ETotal))/(4.d0*pi*znum)
            else
                keldysh=(c*omega0/ETotal)*sqrt( 1.d0-(IPHartree/c2)**2 ) ! relativistic parameter

                if (keldysh .le. 2.d0) then ! tunnel ionization regime (<0.5), and intermediate regime (0.5<keldysh<2)
                    ! Joachain "Atoms in Intense Laser Fields", EQ (6.22).  OK for lower atomic number ions (znum <= 20 with maximum ~20% error)
                    !                rateT=(4.d0*znum**5/ETotal)*exp(-2.d0*znum**3/(3.d0*ETotal))

                    ! relativistic semiclassical tunneling rate (Milosevic - "Semiclassical Dirac Theory of Tunnel Ionization", Eq. 8)
                    tmu=znum
                    teps=sqrt(1.d0-tmu**2)
                    txi=sqrt( 1.d0-0.5d0*teps*(sqrt(teps**2+8.d0)-teps) )
                    txi2=txi**2
                    txi3=txi**3
                    rateT=ETotal**(1.d0-2.d0*teps)/(2.d0*sqrt(3.d0)*txi*gamma(2.d0*teps+1.d0))&
                        *sqrt( (3.d0-txi2)/(3.d0+txi2) )&
                        *( (4.d0*txi3*(3.d0-txi2)**2)/(sqrt(3.d0)*(1+txi2)) )**(2.d0*teps)&
                        *exp( 6.d0*tmu*asin(txi/sqrt(3.d0)) - (2.d0*sqrt(3.d0)*txi3)/(ETotal*(1.d0+txi2)) )
                elseif (keldysh .ge. 0.5d0) then ! multiphoton ionization regime (>2), and intermediate regime (0.5<keldysh<2)
                    ! Lawrence-Douglas "Ionization Effects for Laser-Plasma Interactions by Particle-in-Cell Code", Eq.
                    nPhotonsK=floor(IPHartree/omega0+1.d0)
                    sigma_nK=4.8d0*(1.3d0)**(2.d0*nPhotonsK)*ETotal**(2.d0*nPhotonsK-2.d0)/(c*fact(nPhotonsK)**2*omega0**((10.d0*nPhotonsK-1.d0)/3.d0)*sqrt(real(nPhotonsK))*((2.d0*nPhotonsK-1.d0)))
                    rateMP=sigma_nK*(c*ETotal**2/(8.d0*pi*omega0))**nPhotonsK
                endif
            endif !ionizationMode==1

            prob=(rateT+rateMP+rateADK+rateBSI)*dt
            call random_number(rand)
            if (rand.le.prob) then
                ionized=.true.
                write(*,*) "Rel. Keldysh value is",keldysh
                if (omega0>IPHartree) then
                    write(*,*) "Rate dominates in BSI regime."
                elseif (keldysh .ge. 0.5d0) then
                    write(*,*) "Rate dominates in MPI or mixed regime."
                elseif (keldysh .le. 2.d0) then
                    write(*,*) "Rate dominates in T or mixed regime."
                endif
            else
                ionized=.false.
            endif
        endif ! ionized?

        return
    end subroutine

    subroutine makeFields(x,y,z,t)
        implicit none
        real*8 x,y,z,t,smp,cmp
        complex*16 dx2U,dy2U,dz2U,dxdyU,dxdzU,dydzU,dxdtU,dydtU,dzdtU ! cartesian deriv (local variables)
        complex*16 drU,dpU,dr2U,dp2U,drdpU,drdzU,dpdzU,drdtU,dpdtU ! cylindrical deriv (local variables)
        complex*16 p00,p10,p01 ! repeated polynomials (local variables)
        complex*16 tempCirc(1:6)

        ! for BGV field types only...
        real*8 bEps !c
        complex*16 bSum0,bSum1,bf0,bf1,dvbf0,dvbf1,d2vbf0,d2vbf1,m2v1 !c
        complex*16 dhdz,dvdh,dvdr,d2vdrdh,d2vdr2,dzU,dUdh,drdhU !c
        complex*16 bU,bU0,bUv,bUvr,bh,bh2,bv,bPhase !c
            ! specific things for BGV w/ Poisson-like spectrum (whosFields=3)
        integer j
        complex*16 bT,bxi,bbeta,bSum2,&
            bSum42,bSum43,bSum44,bSum63,bSum64,bSum65,bSum66,&
            bbInv,bbInv2,bbInv3,bbInv4,bTemp,bTemp2,bTemp3,bTemp4,bTemp6
        complex*16 dUdbbeta,dUdbT,dUdbxi,&
            dUbetadbxi,dUbetadbT,dUxidbT,dUxidbxi,d2UdbT2
            !dUbetadr,dbxidr,dUdt,d2Udrdbxi,d2UdrdbT

        ! plasma coordinates
        real*8 CGSz,CGSr

        fields=0.d0

        if(x==0.d0 .and. y==0.d0) then ! to catch potential divide by r=0 errors
            x=1.0d-5
            y=1.0d-5
        endif

        ! Calcluate ONLY the pump/wake fields?  (injection pulse lagging)
        plasmaOnly=0
        if(plasmaWake==1) then
            ! 1.d2*z*a0=CGSz, and likewise for rho.  t [CGS] = t [SI] = 2.41888E-17 [au].
            CGSz=1.d2*z*a0
            CGSr=1.d2*r*a0

            bigPsi = (c*t - z)*1.d2*a0 ! in CGS
            if( bigPsi .lt. bigPsiInjection) plasmaOnly=1 !if ahead of injection pulse
        endif


        ! "r" is the cylindrical radius rho, not to be confused with the spherical radius r.
        r2=x**2+y**2
        r=sqrt(r2)
        phi=atan2(y,x) ! atan2 finds the angle in the proper quadrant, since atan is periodic and can't do this with a single argument

        ! these are used in April's fields, and in ODEModule (and therefore trajectories for all fields)
        sp=sin(phi)
        cp=cos(phi)

        if (whosFields==1.and.plasmaOnly==0) then

            !!!!!!!!!!!!!!!!!!!!!!! ~~~~~~~~~~~~~~~~~~~~ !!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!! ~~~~~~~~~~~~~~~~~~~~ !!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!! April's exact fields !!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!! ~~~~~~~~~~~~~~~~~~~~ !!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!! ~~~~~~~~~~~~~~~~~~~~ !!!!!!!!!!!!!!!!!!!!!!!

            smp=sin(mReal*phi)
            cmp=cos(mReal*phi)

            zia=z+ii*a
            zia2=zia**2

            if(rtType==1) then
                Rt=sqrt(r2+zia2) ! choice 1, expect discontinuities at z=0 for odd m with complex legendre
            elseif(rtType==2) then
                Rt=ii*sqrt(-r2-zia2) ! choice 2, should not be discontinuous for complex legendre
            else
                write(*,*) "Invalid rtType!"
                stop
            endif

            if( fixSingularity==0 .or. ( (a-r)**2 + z**2 .gt. 0.1d0*w0 ) ) then ! we don't worry about the numerical ring singularity
                crt=c/Rt
                rrt=r/Rt
                Tp=1.d0+(ws)*(-ii*t+a/c+ii*Rt/c)
                Tm=1.d0+(ws)*(-ii*t+a/c-ii*Rt/c)

                call makeFABs() ! calculate the T+- pieces(A/Bs) for the fields
                call makeUs() ! calculate the phasor sums and their derivatives

                ! some polynomial calculations that are used multiple times
                p00=Pmm(z,0,0)
                p10=Pmm(z,1,0)
                p01=Pmm(z,0,1)

                ! cylindrical derivatives
                    !U=Gpr/2 cmp p00 Upr
                drU   = Gpr/2.d0*cmp*( p10*Upr + p00*Ur )
                dr2U  = Gpr/2.d0*cmp*( Pmm(z,2,0)*Upr + 2.d0*p10*Ur + p00*Urr )
                dpU   = Gpr/2.d0*(-mReal)*smp*p00*Upr
                dp2U  = -Gpr/2.d0*mReal**2*cmp*p00*Upr
                dz2U  = Gpr/2.d0*cmp*( Pmm(z,0,2)*Upr + 2.d0*p01*Uz + p00*Uzz )
                drdpU = -mReal*Gpr/2.d0*smp*( p10*Upr + p00*Ur )
                drdzU = Gpr/2.d0*cmp*( p10*Uz + Pmm(z,1,1)*Upr + p00*Urz + p01*Ur )
                dpdzU = -mReal*Gpr/2.d0*smp*( p01*Upr + p00*Uz )
                drdtU = Gpr/2.d0*cmp*( p10*Uprt + p00*Urt )
                dpdtU = Gpr/2.d0*(-mReal)*smp*p00*Uprt
                dzdtU = Gpr/2.d0*cmp*( p01*Uprt + p00*Uzt )

            else ! take care of the numerical singularity
                ATilde=ws**(3.d0*mReal/2.d0+1.d0)*(1.d0+ws*(a/c-ii*t))**(-s-2.d0-3.d0*mReal/2.d0)
                ATildeT=ii*ws**(3.d0*mReal/2.d0+2.d0)*(s+2.d0+3.d0*mReal/2.d0)*(1.d0+ws*(a/c-ii*t))**(-s-3.d0-3.d0*mReal/2.d0)

                drU   = GprTilde*cmp*ATilde*mReal*r**(m-1)
                dr2U  = GprTilde*cmp*ATilde*mReal*(mReal-1.d0)*r**(m-2)
                dpU   = -mReal*GprTilde*smp*ATilde*r**m
                dp2U  = -(mReal**2)*GprTilde*cmp*ATilde*r**m
                dz2U  = 0.d0
                drdpU = -mReal*GprTilde*smp*ATilde*mReal*r**(m-1)
                drdzU = 0.d0
                dpdzU = 0.d0
                drdtU = GprTilde*cmp*ATildeT*mReal*r**(m-1)
                dpdtU = -mReal*GprTilde*smp*ATildeT*r**m
                dzdtU = 0.d0
            endif

            !calculating and writing the complex fields
            if(polar==1 .or. polar==3 .or. polar==4) then
                ! writing cartesian derivatives in cylindrical CS
                dx2U=cp**2*dr2U + sp**2/r2*dp2U + sp**2/r*drU + 2.d0*cp*sp/r2*dpU - 2.d0*sp*cp/r*drdpU ! c,d
                dy2U=sp**2*dr2U + cp**2/r2*dp2U - 2.d0*sp*cp/r2*dpU + cp**2/r*drU + 2.d0*sp*cp/r*drdpU ! c,d
                dxdyU=cp*sp*dr2U-cp*sp/r2*dp2U-sp*cp/r*drU+(sp**2-cp**2)/r2*dpU+(cp**2-sp**2)/r*drdpU
                dxdzU=cp*drdzU-sp/r*dpdzU
                dydzU=sp*drdzU+cp/r*dpdzU
                dxdtU=cp*drdtU-sp/r*dpdtU
                dydtU=sp*drdtU+cp/r*dpdtU
            endif
            if(polar==1 .or. polar==4) then ! c,d
                fields(1)=( -dz2U-dy2U+dzdtU*mu0/eta0 )             !Ex
                fields(2)=( dxdyU )                                 !Ey
                fields(3)=( dxdzU-mu0/eta0*dxdtU )                  !Ez
                fields(4)=( (dxdyU/eta0)*mu0 )                      !Bx
                fields(5)=( (-dz2U/eta0-dx2U/eta0+eps0*dzdtU)*mu0 ) !By
                fields(6)=( (dydzU/eta0-eps0*dydtU)*mu0 )           !Bz
            elseif(polar==2) then ! c,d
                fields(1)=( drdzU )                                 !Er
                fields(2)=( dpdzU/r )                               !Ep
                fields(3)=( (-drU-r*dr2U-dp2U/r)/r )                !Ez
                fields(4)=( (dpdtU/r)/c2 )                          !Br
                fields(5)=( (-drdtU)/c2 )                           !Bp
                fields(6)=0.d0                                      !Bz
            elseif(polar==3) then ! c,d
                fields(1)=( -dydtU*mu0/eta0 )                       !Ex
                fields(2)=( dxdtU*mu0/eta0 )                        !Ey
                fields(3)=( 0.d0 )                                  !Ez
                fields(4)=( dxdzU*mu0/eta0 )                        !Bx
                fields(5)=( dydzU*mu0/eta0 )                        !By
                fields(6)=( (-dx2U-dy2U)*mu0/eta0 )                 !Bz
            endif
            if (polar==4) then
                ! assuming RHC polarization.  Change "ii" phases to "-ii" for LHC polarization.
                tempCirc(1)=fields(1)-ii*fields(2)
                tempCirc(2)=fields(2)+ii*fields(1)
                tempCirc(3)=fields(3)+ii*fields(3)
                tempCirc(4)=fields(4)-ii*fields(5)
                tempCirc(5)=fields(5)+ii*fields(4)
                tempCirc(6)=fields(6)+ii*fields(6)
                fields=tempCirc
            endif

            fields=fields*norm
            ETotal=sqrt(real(fields(1))**2+real(fields(2))**2+real(fields(3))**2)

            ! these won't be good for circular.
            phasor(1)=abs(Gpr/2.d0*cmp*Pmm(z,0,0)*Upr)**2
            phasor(2)=real(Gpr/2.d0*cmp*Pmm(z,0,0)*Upr)
            phasor(3)=aimag(Gpr/2.d0*cmp*Pmm(z,0,0)*Upr)

        elseif (whosFields==2 .and. polar==2 .and. plasmaOnly==0) then

            !!!!!!!!!!!!!!!!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!! BGV radial fields with temporal Gaussian !!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!!!!!!!!!!!!!!!!!

            !bConst=(-1.d0)**(nn+m)*2.d0**(2*nn+m)*exp(ii*phi0) ! this is in fields.f95
            bh=(1.d0+ii*z/zr)**(-0.5d0) !c,d
            bh2=1.d0/(1.d0+ii*z/zr) !c,d
            bv=bh**2*r2/w0**2 !c,d
            bEps=1.d0/(k*w0) !c,d
            bPhase=exp(ii*(k*z+mReal*phi)) !c,d
            m2v1=(mReal/(2.d0*bv)-1.d0) !c,d

            ! computing terms in the sum, including f(2j) = bfj.  ex: bf1 = f(j=1) = f(2)
            bf0= fact(nn)  * bAssLagPoly(nn,m,bv) !c,d
            dvbf0= fact(nn)*( -bAssLagPoly(nn-1,m+1,bv) ) !c,d
            d2vbf0=fact(nn)*( bAssLagPoly(nn-2,m+2,bv) ) !c,d
            bf1= 2.d0*fact(nn+1)  * bAssLagPoly(nn+1,m,bv)     - fact(nn+2)* bAssLagPoly(nn+2,m,bv) !c,d
            dvbf1= 2.d0*fact(nn+1)*( -bAssLagPoly(nn,m+1,bv) ) - fact(nn+2)*( -bAssLagPoly(nn+1,m+1,bv) ) !c,d
            d2vbf1=2.d0*fact(nn+1)*( bAssLagPoly(nn-1,m+2,bv) ) - fact(nn+2)*( bAssLagPoly(nn,m+2,bv) ) !c,d
            bSum0=bf0 !c,d
            bSum1=bh2*bEps**2*bf1 !c,d

            ! simple derivatives
            dhdz=-ii*bh**3/(2.d0*zr) !c,d
            dvdh=2.d0*bv/bh !c,d
            dvdr=2.d0*bv/r !c,d
            d2vdrdh=4.d0*bv/(r*bh) !c,d
            d2vdr2=2.d0*bv/r2 !c,d

            ! the phasor and derivative pieces
            bU0=bConst*bPhase*bh**(2*nn+m+2)*bv**(mReal/2.d0)*exp(-bv)*&
                exp(-ii*omega0*t)*exp(-omega0**2*(t-z/c)**2/(2.d0*s)) !c,d
            bU=bU0*(bSum0+bSum1) !c,d

            bUv=m2v1*bU + bU0*(dvbf0+bh2*bEps**2*dvbf1) !c,d, notes Eq 5
            !singular            bUv=bU*( mReal/(2.d0*bv) - 1.d0 + (dvbf0+(bh*bEps)**2*dvbf1)/(bSum0+bSum1) ) ! notes Eq 5
            drU=dvdr*bUv !c,d, notes Eq 4

            bUvr=drU*m2v1 + dvdr*&
                ( bU0*(m2v1*(dvbf0+bh2*bEps**2*dvbf1)+d2vbf0+bh2*bEps**2*d2vbf1)-mReal*bU/(2.d0*bv**2) ) !c,d notes Eq 7
            !singular            bUvr=bUv/bU*drU + bU*dvdr*( -mReal/(2.d0*bv**2) + (d2vbf0+(bh*bEps)**2*d2vbf1)/(bSum0+bSum1)-((dvbf0+(bh*bEps)**2*dvbf1)/(bSum0+bSum1))**2 ) ! notes Eq 7
            dr2U=d2vdr2*bUv + dvdr*bUvr !c,d, notes Eq 6

            dUdh=bU0*( ((2.d0*nnReal+mReal+2.d0)/bh-dvdh)*&
                (bSum0+bSum1) + dvdh*(dvbf0+bh2*bEps**2*dvbf1) +&
                2.d0*bh*bEps**2*bf1 ) !c,d, notes Eq 9
            !singular            dUdh=bU*( (2.d0*nnReal+mReal+2.d0)/bh + 2.d0*bh*bEps**2*bf1/(bSum0+bSum1) ) + dvdh*bUv ! notes Eq 9
            dzU=dUdh*dhdz + bUv*dvdh*dhdz + ii*k*bU !c,d, notes Eq 8

            drdhU=dvdr*dUdh*m2v1  + bU0*( d2vdrdh*&
                (dvbf0+bh2*bEps**2*dvbf1-bf0-bh2*bEps**2*bf1)+&
                (((2.d0*nnReal+mReal+2.d0)/bh-dvdh)*(dvbf0+bh2*bEps**2*dvbf1)+&
                dvdh*(d2vbf0+bh2*bEps**2*d2vbf1)+2.d0*bh*bEps**2*dvbf1 )*dvdr) !c,d, notes Eq 12
            !singular            drdhU=bUvr*dvdh + bUv*d2vdrdh + drU*((2.d0*nnReal+mReal+2.d0)/bh+2.d0*bh*bEps**2*bf1/(bSum0+bSum1)) + &
            !                bU*dvdr*( 2.d0*bh*bEps**2*dvbf1/(bSum0+bSum1) - &
            !                2.d0*bh*bEps**2*bf1/(dvbf0+(bh*bEps)**2*dvbf1) ) ! notes Eq 12
            drdzU=ii*k*bUv*dvdr + dhdz*( drdhU + bUv*d2vdrdh + dvdh*bUvr ) !c,d, notes Eq 11

            ! radial polarization required for these fields
            fields(1)=( drdzU )                                 !Er !c
            fields(2)=( ii*mReal/r * dzU )                      !Ep !c
            fields(3)=( (-drU - r*dr2U + mReal**2*bU/r)/r )     !Ez !c
            fields(4)=( (mReal*omega0*bU/r)/c2 )                !Br !c
            fields(5)=( (ii*omega0*drU)/c2 )                    !Bp !c
            fields(6)=0.d0                                      !Bz !c

            fields=fields*norm
            ETotal=sqrt(real(fields(1))**2+real(fields(2))**2+real(fields(3))**2)

        elseif (whosFields==3 .and. polar==2 .and. plasmaOnly==0) then

            !!!!!!!!!!!!!!!!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!! BGV radial fields with Poisson spectrum !!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!!!!!!!!!!!!!!!!!

            bbeta=ii*z+zr ! phasor p.1, c
            bbInv=1/bbeta
            bbInv2=bbInv**2
            bbInv3=bbInv**3
            bbInv4=bbInv**4

            bxi=r2/(2.d0*c*bbeta) ! phasor p.1, c
            bT=1.d0+ws*(-ii*z/c + bxi + ii*t) ! phasor p.4, c

            bU0=exp(ii*mReal*phi)*(bbeta/zr)**(-nnReal-mReal/2.d0-1.d0)*bxi**(mReal/2.d0) ! phasor p.1, c


            if(mod(pertOrder,2)==1) then
                write(*,*) "Invalid perturbative order.  Even values required (odd orders are zero)."
                write(*,*) "Stopping."
                stop
            endif

            !

            if(pertOrder==4) then



            bSum0=0.d0
            bSum1=0.d0
            bSum2=0.d0
            bSum42=0.d0
            bSum43=0.d0
            bSum44=0.d0
            do j=0,nn+4 ! nothing down
                bTemp=bxi**j
                bTemp2=bTemp*bT**(-bgj(j)+1.d0)
                if(j.le.nn) bSum0=bSum0 + bc0(j)*bTemp2/bT**2!bTemp*bT**(-bgj(j)-1.d0)
                if(j.le.nn+1) bSum1=bSum1 + bc1(j)*bTemp2/bT!bTemp*bT**(-bgj(j))
                if(j.le.nn+2) then
                    bSum2=bSum2 +   bc2(j) *bTemp2/bT!bTemp*bT**(-bgj(j))
                    bSum42=bSum42 + bc42(j)*bTemp2
                endif
                if(j.le.nn+3) bSum43=bSum43 + bc43(j)*bTemp2
                bSum44=bSum44 + bc44(j)*bTemp2
            enddo
            bU=Cmn*bU0*( bSum0 + bbInv*(bSum1-bSum2) + bbInv2*(bSum42-bSum43+bSum44) ) !c
            dUdbbeta=Cmn*bU0*( -bbInv2*(bSum1-bSum2) - 2.d0*bbInv3*(bSum42-bSum43+bSum44) ) !c

            bSum0=0.d0
            bSum1=0.d0
            bSum2=0.d0
            bSum42=0.d0
            bSum43=0.d0
            bSum44=0.d0
            do j=0,nn+4 ! gamma down, bT-deriv
                bTemp=bxi**j
                bTemp2=bTemp*bT**(-bgj(j))
                bTemp3=(1.d0-bgj(j))*bTemp2
                if(j.le.nn) bSum0=bSum0 + bc0(j)*(-bgj(j)-1.d0)*bTemp2/bT**2!bTemp*bT**(-bgj(j)-2.d0)
                if(j.le.nn+1) bSum1=bSum1 + bc1(j)*(-bgj(j))*bTemp2/bT!bTemp*bT**(-bgj(j)-1.d0)
                if(j.le.nn+2) then
                    bSum2=bSum2 +   bc2(j) *(-bgj(j))*bTemp2/bT!bTemp*bT**(-bgj(j)-1.d0)
                    bSum42=bSum42 + bc42(j)*bTemp3
                endif
                if(j.le.nn+3) bSum43=bSum43 + bc43(j)*bTemp3
                bSum44=bSum44 + bc44(j)*bTemp3
            enddo
            dUdbT=Cmn*bU0*( bSum0 + bbInv*(bSum1-bSum2) + bbInv2*(bSum42-bSum43+bSum44) ) !c
            dUbetadbT=Cmn*bU0*( -bbInv2*(bSum1-bSum2) - 2.d0*bbInv3*(bSum42-bSum43+bSum44) ) !c

            bSum0=0.d0
            bSum1=0.d0
            bSum2=0.d0
            bSum42=0.d0
            bSum43=0.d0
            bSum44=0.d0
            do j=0,nn+4 ! j down, bXi-deriv
                bTemp=real(j,kind=8)*bxi**(j-1)
                bTemp2=bTemp*bT**(-bgj(j)+1.d0)
                if(j.le.nn) bSum0=bSum0 + bc0(j)*bTemp2/bT**2!bTemp*bT**(-bgj(j)-1.d0)
                if(j.le.nn+1) bSum1=bSum1 + bc1(j)*bTemp2/bT!bTemp*bT**(-bgj(j))
                if(j.le.nn+2) then
                    bSum2=bSum2 +   bc2(j) *bTemp2/bT!bTemp*bT**(-bgj(j))
                    bSum42=bSum42 + bc42(j)*bTemp2
                endif
                if(j.le.nn+3) bSum43=bSum43 + bc43(j)*bTemp2
                bSum44=bSum44 + bc44(j)*bTemp2
            enddo
            dUdbxi=Cmn*bU0*( bSum0 + bbInv*(bSum1-bSum2) + bbInv2*(bSum42-bSum43+bSum44) ) !c
            dUbetadbxi=Cmn*bU0*( -bbInv2*(bSum1-bSum2) - 2.d0*bbInv3*(bSum42-bSum43+bSum44) ) !c

            bSum0=0.d0
            bSum1=0.d0
            bSum2=0.d0
            bSum42=0.d0
            bSum43=0.d0
            bSum44=0.d0
            do j=0,nn+4 ! gamma down twice
                bTemp=bxi**j
                bTemp2=bTemp*bT**(-bgj(j)-1.d0)
                bTemp3=(1.d0-bgj(j))*(-bgj(j))*bTemp2
                if(j.le.nn) bSum0=bSum0 + bc0(j)*(bgj(j)+1.d0)*(bgj(j)+2.d0)*bTemp2/bT**2!bTemp*bT**(-bgj(j)-3.d0)
                if(j.le.nn+1) bSum1=bSum1 + bc1(j)*bgj(j)*(bgj(j)+1.d0)*bTemp2/bT!bTemp*bT**(-bgj(j)-2.d0)
                if(j.le.nn+2) then
                    bSum2=bSum2 +   bc2(j) *bgj(j)*(bgj(j)+1.d0)   *bTemp2/bT!bTemp*bT**(-bgj(j)-2.d0)
                    bSum42=bSum42 + bc42(j)*bTemp3
                endif
                if(j.le.nn+3) bSum43=bSum43 + bc43(j)*bTemp3
                bSum44=bSum44 + bc44(j)*bTemp3
            enddo
            d2UdbT2=Cmn*bU0*( bSum0 + bbInv*(bSum1-bSum2) + bbInv2*(bSum42-bSum43+bSum44) ) !c

            bSum0=0.d0
            bSum1=0.d0
            bSum2=0.d0
            bSum42=0.d0
            bSum43=0.d0
            bSum44=0.d0
            do j=0,nn+4 ! j and gamma down
                bTemp=real(j,kind=8)*bxi**(j-1)
                bTemp2=bTemp*bT**(-bgj(j))
                bTemp3=(1.d0-bgj(j))*bTemp2
                if(j.le.nn) bSum0=bSum0 + bc0(j)*(-bgj(j)-1.d0)*bTemp2/bT**2!bTemp*bT**(-bgj(j)-2.d0)
                if(j.le.nn+1) bSum1=bSum1 + bc1(j)*(-bgj(j))*bTemp2/bT!bTemp*bT**(-bgj(j)-1.d0)
                if(j.le.nn+2) then
                    bSum2=bSum2 +   bc2(j) *(-bgj(j))*bTemp2/bT!bTemp*bT**(-bgj(j)-1.d0)
                    bSum42=bSum42 + bc42(j)*bTemp3
                endif
                if(j.le.nn+3) bSum43=bSum43 + bc43(j)*bTemp3
                bSum44=bSum44 + bc44(j)*bTemp3
            enddo
            dUxidbT=Cmn*bU0*( bSum0 + bbInv*(bSum1-bSum2) + bbInv2*(bSum42-bSum43+bSum44) ) !c

            bSum0=0.d0
            bSum1=0.d0
            bSum2=0.d0
            bSum42=0.d0
            bSum43=0.d0
            bSum44=0.d0
            do j=0,nn+4 ! j down twice
                bTemp=real(j,kind=8)*(real(j,kind=8)-1.d0)*bxi**(j-2)
                bTemp2=bTemp*bT**(-bgj(j)+1.d0)
                if(j.le.nn) bSum0=bSum0 + bc0(j)*bTemp2/bT**2!bTemp*bT**(-bgj(j)-1.d0)
                if(j.le.nn+1) bSum1=bSum1 + bc1(j)*bTemp2/bT!bTemp*bT**(-bgj(j))
                if(j.le.nn+2) then
                    bSum2=bSum2 +   bc2(j) *bTemp2/bT!bTemp*bT**(-bgj(j))
                    bSum42=bSum42 + bc42(j)*bTemp2
                endif
                if(j.le.nn+3) bSum43=bSum43 + bc43(j)*bTemp2
                bSum44=bSum44 + bc44(j)*bTemp2
            enddo
            dUxidbxi=Cmn*bU0*( bSum0 + bbInv*(bSum1-bSum2) + bbInv2*(bSum42-bSum43+bSum44) ) !c



            elseif(pertOrder==6) then



            bSum0=0.d0
            bSum1=0.d0
            bSum2=0.d0
            bSum42=0.d0
            bSum43=0.d0
            bSum44=0.d0
            bSum63=0.d0
            bSum64=0.d0
            bSum65=0.d0
            bSum66=0.d0
            do j=0,nnMax ! nothing down
                bTemp=bxi**j
                bTemp2=bTemp*bT**(-bgj(j))
                bTemp4=bTemp*bT**(-bgj(j)+1.d0)
                bTemp6=bTemp*bT**(-bgj(j)+2.d0)
                if(j.le.nn) bSum0=bSum0 + bc0(j)*bTemp*bT**(-bgj(j)-1.d0)
                if(j.le.nn+1) bSum1=bSum1 + bc1(j)*bTemp2
                if(j.le.nn+2) then
                    bSum2=bSum2 +   bc2(j) *bTemp2
                    bSum42=bSum42 + bc42(j)*bTemp4
                endif
                if(j.le.nn+3) then
                    bSum43=bSum43 + bc43(j)*bTemp4
                    bSum63=bSum63 + bc63(j)*bTemp6
                endif
                if(j.le.nn+4) then
                    bSum44=bSum44 + bc44(j)*bTemp4
                    bSum64=bSum64 + bc64(j)*bTemp6
                endif
                if(j.le.nn+5) then
                    bSum65=bSum65 + bc65(j)*bTemp6
                endif
                bSum66=bSum66 + bc66(j)*bTemp6
            enddo
            bU=Cmn*bU0*( bSum0 + bbInv*(bSum1-bSum2) + bbInv2*(bSum42-bSum43+bSum44)  +&
                 bbInv3*(bSum63-bSum64+bSum65-bSum66) )
            dUdbbeta=Cmn*bU0*( -bbInv2*(bSum1-bSum2) - 2.d0*bbInv3*(bSum42-bSum43+bSum44) -&
                3.d0*bbInv4*(bSum63-bSum64+bSum65-bSum66) )

            bSum0=0.d0
            bSum1=0.d0
            bSum2=0.d0
            bSum42=0.d0
            bSum43=0.d0
            bSum44=0.d0
            bSum63=0.d0
            bSum64=0.d0
            bSum65=0.d0
            bSum66=0.d0
            do j=0,nnMax ! gamma down, bT-deriv
                bTemp=bxi**j
                bTemp2=bTemp*(-bgj(j))*bT**(-bgj(j)-1.d0)
                bTemp4=(1.d0-bgj(j))*bT**(-bgj(j))
                bTemp6=(2.d0-bgj(j))*bT**(-bgj(j)+1.d0)
                if(j.le.nn) bSum0=bSum0 + bc0(j)*(-bgj(j)-1.d0)*bTemp*bT**(-bgj(j)-2.d0)
                if(j.le.nn+1) bSum1=bSum1 + bc1(j)*bTemp2
                if(j.le.nn+2) then
                    bSum2=bSum2 +   bc2(j) *bTemp2
                    bSum42=bSum42 + bc42(j)*bTemp4
                endif
                if(j.le.nn+3) then
                    bSum43=bSum43 + bc43(j)*bTemp4
                    bSum63=bSum63 + bc63(j)*bTemp6
                endif
                if(j.le.nn+4) then
                    bSum44=bSum44 + bc44(j)*bTemp4
                    bSum64=bSum64 + bc64(j)*bTemp6
                endif
                if(j.le.nn+5) then
                    bSum65=bSum65 + bc65(j)*bTemp6
                endif
                bSum66=bSum66 + bc66(j)*bTemp6
            enddo
            dUdbT=Cmn*bU0*( bSum0 + bbInv*(bSum1-bSum2) + bbInv2*(bSum42-bSum43+bSum44) +&
                bbInv3*(bSum63-bSum64+bSum65-bSum66) )
            dUbetadbT=Cmn*bU0*( -bbInv2*(bSum1-bSum2) - 2.d0*bbInv3*(bSum42-bSum43+bSum44) -&
                3.d0*bbInv4*(bSum63-bSum64+bSum65-bSum66) )

            bSum0=0.d0
            bSum1=0.d0
            bSum2=0.d0
            bSum42=0.d0
            bSum43=0.d0
            bSum44=0.d0
            bSum63=0.d0
            bSum64=0.d0
            bSum65=0.d0
            bSum66=0.d0
            do j=0,nnMax ! j down, bXi-deriv
                bTemp=real(j,kind=8)*bxi**(j-1)
                bTemp2=bTemp*bT**(-bgj(j))
                bTemp4=bTemp*bT**(-bgj(j)+1.d0)
                bTemp6=bTemp*bT**(-bgj(j)+2.d0)
                if(j.le.nn) bSum0=bSum0 + bc0(j)*bTemp*bT**(-bgj(j)-1.d0)
                if(j.le.nn+1) bSum1=bSum1 + bc1(j)*bTemp2
                if(j.le.nn+2) then
                    bSum2=bSum2 +   bc2(j) *bTemp2
                    bSum42=bSum42 + bc42(j)*bTemp4
                endif
                if(j.le.nn+3) then
                    bSum43=bSum43 + bc43(j)*bTemp4
                    bSum63=bSum63 + bc63(j)*bTemp6
                endif
                if(j.le.nn+4) then
                    bSum44=bSum44 + bc44(j)*bTemp4
                    bSum64=bSum64 + bc64(j)*bTemp6
                endif
                if(j.le.nn+5) then
                    bSum65=bSum65 + bc65(j)*bTemp6
                endif
                bSum66=bSum66 + bc66(j)*bTemp6
            enddo
            dUdbxi=Cmn*bU0*( bSum0 + bbInv*(bSum1-bSum2) + bbInv2*(bSum42-bSum43+bSum44) +&
                bbInv3*(bSum63-bSum64+bSum65-bSum66) )
            dUbetadbxi=Cmn*bU0*( -bbInv2*(bSum1-bSum2) - 2.d0*bbInv3*(bSum42-bSum43+bSum44) -&
                3.d0*bbInv4*(bSum63-bSum64+bSum65-bSum66) )

            bSum0=0.d0
            bSum1=0.d0
            bSum2=0.d0
            bSum42=0.d0
            bSum43=0.d0
            bSum44=0.d0
            bSum63=0.d0
            bSum64=0.d0
            bSum65=0.d0
            bSum66=0.d0
            do j=0,nnMax ! gamma down twice
                bTemp=bxi**j
                bTemp2=bTemp*(-bgj(j))*(-bgj(j)-1.d0)*bT**(-bgj(j)-2.d0)
                bTemp4=(1.d0-bgj(j))*(-bgj(j))*bT**(-bgj(j)-1.d0)
                bTemp6=(2.d0-bgj(j))*(-bgj(j)+1.d0)*bT**(-bgj(j))
                if(j.le.nn) bSum0=bSum0 + bc0(j)*(-bgj(j)-1.d0)*(-bgj(j)-2.d0)*bTemp*bT**(-bgj(j)-3.d0)
                if(j.le.nn+1) bSum1=bSum1 + bc1(j)*bTemp2
                if(j.le.nn+2) then
                    bSum2=bSum2 +   bc2(j) *bTemp2
                    bSum42=bSum42 + bc42(j)*bTemp4
                endif
                if(j.le.nn+3) then
                    bSum43=bSum43 + bc43(j)*bTemp4
                    bSum63=bSum63 + bc63(j)*bTemp6
                endif
                if(j.le.nn+4) then
                    bSum44=bSum44 + bc44(j)*bTemp4
                    bSum64=bSum64 + bc64(j)*bTemp6
                endif
                if(j.le.nn+5) then
                    bSum65=bSum65 + bc65(j)*bTemp6
                endif
                bSum66=bSum66 + bc66(j)*bTemp6
            enddo
            d2UdbT2=Cmn*bU0*( bSum0 + bbInv*(bSum1-bSum2) + bbInv2*(bSum42-bSum43+bSum44) +&
                bbInv3*(bSum63-bSum64+bSum65-bSum66) )

            bSum0=0.d0
            bSum1=0.d0
            bSum2=0.d0
            bSum42=0.d0
            bSum43=0.d0
            bSum44=0.d0
            bSum63=0.d0
            bSum64=0.d0
            bSum65=0.d0
            bSum66=0.d0
            do j=0,nnMax ! j and gamma down
                bTemp=real(j,kind=8)*bxi**(j-1)
                bTemp2=bTemp*(-bgj(j))*bT**(-bgj(j)-1.d0)
                bTemp4=(1.d0-bgj(j))*bT**(-bgj(j))
                bTemp6=(2.d0-bgj(j))*bT**(-bgj(j)+1.d0)
                if(j.le.nn) bSum0=bSum0 + bc0(j)*(-bgj(j)-1.d0)*bTemp*bT**(-bgj(j)-2.d0)
                if(j.le.nn+1) bSum1=bSum1 + bc1(j)*bTemp2
                if(j.le.nn+2) then
                    bSum2=bSum2 +   bc2(j) *bTemp2
                    bSum42=bSum42 + bc42(j)*bTemp4
                endif
                if(j.le.nn+3) then
                    bSum43=bSum43 + bc43(j)*bTemp4
                    bSum63=bSum63 + bc63(j)*bTemp6
                endif
                if(j.le.nn+4) then
                    bSum44=bSum44 + bc44(j)*bTemp4
                    bSum64=bSum64 + bc64(j)*bTemp6
                endif
                if(j.le.nn+5) then
                    bSum65=bSum65 + bc65(j)*bTemp6
                endif
                bSum66=bSum66 + bc66(j)*bTemp6
            enddo
            dUxidbT=Cmn*bU0*( bSum0 + bbInv*(bSum1-bSum2) + bbInv2*(bSum42-bSum43+bSum44) +&
                bbInv3*(bSum63-bSum64+bSum65-bSum66) )

            bSum0=0.d0
            bSum1=0.d0
            bSum2=0.d0
            bSum42=0.d0
            bSum43=0.d0
            bSum44=0.d0
            bSum63=0.d0
            bSum64=0.d0
            bSum65=0.d0
            bSum66=0.d0
            do j=0,nnMax ! j down twice
                bTemp=real(j,kind=8)*real(j-1,kind=8)*bxi**(j-2)
                bTemp2=bTemp*bT**(-bgj(j))
                bTemp4=bTemp*bT**(-bgj(j)+1.d0)
                bTemp6=bTemp*bT**(-bgj(j)+2.d0)
                if(j.le.nn) bSum0=bSum0 + bc0(j)*bTemp*bT**(-bgj(j)-1.d0)
                if(j.le.nn+1) bSum1=bSum1 + bc1(j)*bTemp2
                if(j.le.nn+2) then
                    bSum2=bSum2 +   bc2(j) *bTemp2
                    bSum42=bSum42 + bc42(j)*bTemp4
                endif
                if(j.le.nn+3) then
                    bSum43=bSum43 + bc43(j)*bTemp4
                    bSum63=bSum63 + bc63(j)*bTemp6
                endif
                if(j.le.nn+4) then
                    bSum44=bSum44 + bc44(j)*bTemp4
                    bSum64=bSum64 + bc64(j)*bTemp6
                endif
                if(j.le.nn+5) then
                    bSum65=bSum65 + bc65(j)*bTemp6
                endif
                bSum66=bSum66 + bc66(j)*bTemp6
            enddo
            dUxidbxi=Cmn*bU0*( bSum0 + bbInv*(bSum1-bSum2) + bbInv2*(bSum42-bSum43+bSum44) +&
                bbInv3*(bSum63-bSum64+bSum65-bSum66) )

            else!pertOrder!=valid
                write(*,*) "Undefined perturbative order.  use pertOrder=(2j)_max >=4.  Stopping."
                stop
            endif

            fields(1)=-ii/r*( mReal*(nnReal+mReal+1.d0)/bbeta*bU - 2.d0*ws*bxi*dUbetadbT + &
                        ws*(2.d0*bxi*(nnReal+mReal+2.d0)/bbeta+mReal*(bxi/bbeta+1.d0/c))*dUdbT + &
                        bxi*(2.d0*nnReal+3.d0*mReal+4.d0)/bbeta*dUdbxi + 2.d0*bxi*ws**2*(bxi/bbeta+1.d0/c)*d2UdbT2 + &
                        2.d0*bxi**2/bbeta*dUxidbxi + 2.d0*ws*bxi*(2.d0*bxi/bbeta+1.d0/c)*dUxidbT &
                        - 2.d0*bxi*dUbetadbxi - mReal*dUdbbeta ) ! c,d

            fields(2)=mReal/r*( (nnReal+mReal+1.d0)/bbeta*bU  + ws*(bxi/bbeta+1.d0/c)*dUdbT + &
                        bxi/bbeta*dUdbxi - dUdbbeta ) ! c,d

            fields(3)=bxi/r2*( - 4.d0*ws*(mReal+1.d0)*dUdbT - 4.d0*(mReal+1.d0)*dUdbxi &
                        - 4.d0*ws**2*bxi*d2UdbT2 - 4.d0*bxi*dUxidbxi - 8.d0*ws*bxi*dUxidbT ) ! c,d

            fields(4)=-mReal*ws/(c2*r)*dUdbT ! c,d

            fields(5)=-ii*ws/(c2*r)*( mReal*dUdbT + 2.d0*bxi*(ws*d2UdbT2 + dUxidbT) ) ! c,d

            fields(6)=0.d0 ! c,d

            fields=fields*norm
            ETotal=sqrt(real(fields(1))**2+real(fields(2))**2+real(fields(3))**2)

            phasor(1)=abs(bU)**2
            !phasor(2)=real(bU)
            !phasor(3)=aimag(bU)

        else !whosFields =/= 1,2,3, or polarization other than radial chosen for BGV
            write(*,*) "Invalid field choice (April vs B&G-V).  Stopping."
            stop
        endif !whosFields

        if(plasmaWake==1.and.normed==1) then
            if(polar.ne.2) then
                write(*,*) "Stopping.  Invalid polarization for including the plasma wake."
                stop
            endif
            plasmaFields=0.d0

            ! t [CGS] = t [SI] = 2.41888E-17 [au].
            psip=omegap*(CGSz/(CGSc*betaphi)-2.41888d-17*t)

            ! fields from Eq.(15-16) of "Generation of ultrashort electron bunches..."
            ! PRE 1999: Schroeder, Esarey, Leemans, Lee, Wurtele

            if(bigPsi .gt. lambdap ) then !behind the pump pulse, use osc. fields
                plasmaFields(1)=(pi*a0hat**2)*(r/w0pump**2)*cos(psip)
                plasmaFields(3)=(0.25d0*pi*a0hat**2)*(omegap/CGSc)*sin(psip)
            elseif( (bigPsi .le. lambdap) .and. (bigPsi.ge.0.d0) ) then !interaction region
                plasmaFields(1)=(a0hat**2)*(r/w0pump**2)*(1.d0+sin(psip)+cos(psip)*(3.d0*pi*0.25d0-psip*0.5d0))
                plasmaFields(3)=-(0.25d0*a0hat**2)*(omegap/CGSc)*&
                                ( cos(psip) - sin(psip)*(3.d0*pi*0.25d0-psip*0.5d0) - 0.5d0*cos(psip) )
            !else...ahead of pump, plasma wake fields already set to zero.
            endif!bigPsi
            plasmaFields=plasmaFields*(CGSme*CGSc**2/CGSqe)*exp(-2.d0*CGSr**2/w0pump**2)*EAUperECGS

            fields=fields+plasmaFields
        endif ! plasmaWake

        return
    end subroutine

    subroutine makeBGVPvars() ! c,d
        implicit none
        integer jindex

        do jindex=0,nnMax
            bgj(jindex)=bgamma(jindex)
            if(jindex.le.nn) bc0(jindex)=fbc0(jindex)
            if(jindex.le.nn+1)bc1(jindex)=fbc1(jindex)
            if(jindex.le.2) then
                bc2(jindex)=fbc2(jindex)
                bc42(jindex)=fbc42(jindex)
            endif
            if(jindex.le.3) then
                bc43(jindex)=fbc43(jindex)
                bc63(jindex)=fbc63(jindex)
            endif
            if(jindex.le.4) then
                bc44(jindex)=fbc44(jindex)
                bc64(jindex)=fbc64(jindex)
            endif
            if(jindex.le.5) then
                bc65(jindex)=fbc65(jindex)
            endif
            bc66(jindex)=fbc66(jindex)
        enddo
        return
    end subroutine

    subroutine makeF01() ! c,d
        implicit none
        integer i
        real*8 m0,p1
        complex*16 mim,mim1

        mim=(-ii)**m
        mim1=(-ii)**(m+1)

        do i=0,kmax1
            m0=2.d0*i-spr
            p1=m0+1.d0

            F0(i)=(-1.d0)**i*aterm(2*i,m)*mim1*ws**((-1.d0)*(s+p1)) ! c,d
            F0(i)=F0(i)*exp( log_gamma(-m0)-Log_gamma(s+1.d0) )

            if(i.le.kmax2) then
                F1(i)=(-1.d0)**i*aterm(2*i+1,m)*mim*(ws)**((-1.d0)*(s+1.d0+p1)) ! c,d
                F1(i)=F1(i)*exp( log_gamma(-p1)-Log_gamma(s+1.d0) )
            endif
        enddo
        return
    end subroutine

    subroutine makeFABs() ! c
        implicit none
        integer i
        real*8 m0,m1,m2,p1,posneg
        complex*16 tm0,tm1,tm2,tp0,tp1,tp2

        posneg=(-1.d0)**m

        do i=0,kmax1
            m0=2.d0*i-spr
            m1=m0-1.d0
            m2=m0-2.d0
            p1=m0+1.d0
            tm0=tm**m0
            tm1=tm**m1
            tm2=tm**m2
            tp0=tp**m0
            tp1=tp**m1
            tp2=tp**m2

            FA(i)=tm0-posneg*tp0 ! c,d

            FAr(i)=m0*iwsc*( -tm1 - posneg*tp1 )*rrt ! c,d
            FAz(i)=m0*iwsc*( -tm1 - posneg*tp1 )*zia/Rt ! c,d
            FAt(i)=-ii*ws*m0*( tm1 - posneg*tp1 ) ! c

            FArr(i)=m0*( m1*tm2*(-iwsc)**2*rrt**2      + tm1*(-iwsc)/Rt + tm1*iwsc*r2/(Rt**3) )   - posneg*m0*iwsc*( m1*tp2*iwsc*rrt**2      + tp1*( 1.d0/Rt - r2/(Rt**3) ) ) ! c
            FAzz(i)=m0*( m1*tm2*(-iwsc)**2*(zia/Rt)**2 + tm1*(-iwsc)/Rt + tm1*iwsc*zia2/(Rt**3) ) - posneg*m0*iwsc*( m1*tp2*iwsc*(zia/Rt)**2 + tp1*( 1.d0/Rt - zia2/(Rt**3) ) ) ! c
            FArz(i)=m0*( m1*tm2*(-iwsc)**2*r*zia/Rt**2 +                  tm1*iwsc*r*zia/(Rt**3) ) - posneg*m0*iwsc*( m1*tp2*iwsc*r*zia/Rt**2 + tp1*(        -r*zia/(Rt**3) ) ) ! c
            FArt(i)=(-ii*ws*m1)*(m0*rrt*iwsc)*( -tm2 - posneg*tp2 ) ! c
            FAzt(i)=(-ii*ws*m1)*(m0*zia/Rt*iwsc)*( -tm2 - posneg*tp2 ) ! c

            if(i.le.kmax2) then
                FB(i)=Tm**p1+posneg*Tp**p1 ! c,d

                FBr(i)=p1*iwsc*( -tm0 + posneg*tp0 )*rrt ! c,d
                FBz(i)=p1*iwsc*( -tm0 + posneg*tp0 )*zia/Rt ! c,d
                FBt(i)=-ii*ws*p1*( tm0 + posneg*tp0 ) ! c

                FBrr(i)=p1*( m0*tm1*(-iwsc)**2*rrt**2      + tm0*(-iwsc)/Rt + tm0*iwsc*r2/(Rt**3) )   + posneg*p1*iwsc*( m0*tp1*iwsc*rrt**2      + tp0*( 1.d0/Rt - r2/(Rt**3) ) ) ! c
                FBzz(i)=p1*( m0*tm1*(-iwsc)**2*(zia/Rt)**2 + tm0*(-iwsc)/Rt + tm0*iwsc*zia2/(Rt**3) ) + posneg*p1*iwsc*( m0*tp1*iwsc*(zia/Rt)**2 + tp0*( 1.d0/Rt - zia2/(Rt**3) ) ) ! c
                FBrz(i)=p1*( m0*tm1*(-iwsc)**2*r*zia/Rt**2 +                  tm0*iwsc*r*zia/(Rt**3) ) + posneg*p1*iwsc*( m0*tp1*iwsc*r*zia/Rt**2 + tp0*( -r*zia/(Rt**3) ) ) ! c
                FBrt(i)=(rrt*iwsc*p1)*(-ii*ws*m0)*( -tm1 + posneg*tp1 ) ! c
                FBzt(i)=(zia/Rt*iwsc*p1)*(-ii*ws*m0)*( -tm1 + posneg*tp1 ) ! c
            endif
        enddo
        return
    end subroutine

    subroutine makeUs()
        implicit none
        complex*16 sum1(1:9),sum2(1:9),crtPow,rtPow,rtPowP,rtPowM
        real*8 cPow,cPowMult,plus
        integer i

        sum1=0.d0
        sum2=0.d0

        do i=0,kmax1
            plus=2.d0*i+3.d0
            cPow=c**(plus-2.d0)
            cPowMult=cPow*(plus-2.d0)
            crtPow=crt**(plus-2.d0)
            rtPow=Rt**(plus)
            rtPowP=Rt**(plus+2.d0)
            rtPowM=Rt**(plus-2.d0)

            sum1(1)=sum1(1)+F0(i)*FA(i)*crtPow
            sum1(2)=sum1(2)+F0(i)*( -FA(i)*cPowMult*r/rtPow + crtPow*FAr(i) ) !dU'/dr sum 1 ... c,d
            sum1(3)=sum1(3)+F0(i)*( -FA(i)*cPowMult*zia/rtPow + crtPow*FAz(i) )
            sum1(4)=sum1(4)+F0(i)*FAt(i)*crtPow
            sum1(5)=sum1(5)+F0(i)*( -2.d0*FAr(i)*cPowMult*r/rtPow + FA(i)*cPowMult*(-1.d0/rtPow + plus*r2/rtPowP) + cPow*FArr(i)/rtPowM )
            sum1(6)=sum1(6)+F0(i)*( -2.d0*FAz(i)*cPowMult*zia/rtPow + FA(i)*cPowMult*(-1.d0/rtPow + plus*zia2/rtPowP) + cPow*FAzz(i)/rtPowM )
            sum1(7)=sum1(7)+F0(i)*( cPowMult*r*(-FAz(i)/rtPow+FA(i)*zia*plus/rtPowP) + (cPow*FArz(i)/rtPowM-FAr(i)*cPowMult*zia/rtPow) )
            sum1(8)=sum1(8)+F0(i)*( -FAt(i)*cPowMult*r/rtPow + crtPow*FArt(i) )
            sum1(9)=sum1(9)+F0(i)*( -FAt(i)*cPowMult*zia/rtPow + crtPow*FAzt(i) )
        enddo

        do i=0,kmax2
            plus=2.d0*i+4.d0
            cPow=c**(plus-2.d0)
            cPowMult=cPow*(plus-2.d0)
            crtPow=crt**(plus-2.d0)
            rtPow=Rt**(plus)
            rtPowP=Rt**(plus+2.d0)
            rtPowM=Rt**(plus-2.d0)

            sum2(1)=sum2(1)+F1(i)*FB(i)*crtPow
            sum2(2)=sum2(2)+F1(i)*( -FB(i)*cPowMult*r/rtPow + crtPow*FBr(i) ) !dU'/dr sum 2 ... c,d
            sum2(3)=sum2(3)+F1(i)*( -FB(i)*cPowMult*zia/rtPow + crtPow*FBz(i) )
            sum2(4)=sum2(4)+F1(i)*FBt(i)*crtPow
            sum2(5)=sum2(5)+F1(i)*( -2.d0*FBr(i)*cPowMult*r/rtPow + FB(i)*cPowMult*(-1.d0/rtPow + plus*r2/rtPowP) + cPow*FBrr(i)/rtPowM )
            sum2(6)=sum2(6)+F1(i)*( -2.d0*FBz(i)*cPowMult*zia/rtPow + FB(i)*cPowMult*(-1.d0/rtPow + plus*zia2/rtPowP) + cPow*FBzz(i)/rtPowM )
            sum2(7)=sum2(7)+F1(i)*( cPowMult*r*(-FBz(i)/rtPow+FB(i)*zia*plus/rtPowP) + (cPow*FBrz(i)/rtPowM-FBr(i)*cPowMult*zia/rtPow) )
            sum2(8)=sum2(8)+F1(i)*( -FBt(i)*cPowMult*r/rtPow + crtPow*FBrt(i) )
            sum2(9)=sum2(9)+F1(i)*( -FBt(i)*cPowMult*zia/rtPow + crtPow*FBzt(i) )
        enddo

        Upr=sum1(1)+sum2(1) ! U'        c
        Ur=sum1(2)+sum2(2) ! dU'/dr     c,d
        Uz=sum1(3)+sum2(3) ! dU'/dz     c,d
        Uprt=sum1(4)+sum2(4) ! dU'/dt   c,d
        Urr=sum1(5)+sum2(5) ! d2U'/dr2  c,d
        Uzz=sum1(6)+sum2(6) ! d2U'/dz2  c,d
        Urz=sum1(7)+sum2(7) ! d2U'/drdz c,d
        Urt=sum1(8)+sum2(8) ! d2U'/drdt c,d
        Uzt=sum1(9)+sum2(9) ! d2U'/dzdt c,d

        return
    end subroutine

end module fieldsModule
