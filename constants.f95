module constants
    implicit none

    !----------physical constants----------
    real*8,parameter::&
        pi=4.d0*atan(1.d0),&
        a0=5.2917721092d-11,& ! (SI) bohr radius in meters
        GEVperHartree=27.21138505d-9,& ! GeV per Hartree
        eVperHartree=27.21138505d0,& ! eV per Hartree
        secPerAU=2.41888d-17,& ! seconds per au time
        c=137.035999074d0,c2=c**2,c3=c**3,& ! speed of light
        mu0=4.d0*pi/c2,&
        eps0=1.d0/(4.d0*pi),&
        eta0=sqrt(mu0/eps0) ! impedence of free space
    complex*16,parameter::&
        ii=(0.d0,1.d0) ! imaginary i

!----------global parameters (in a.u. when applicable)----------
    integer,parameter::&
        m=0,& ! angular L-G index
        nn=2,& ! radial L-G index.

        polar=2,& ! 1=linear polarization (x) (cartesian cs),
                  ! 2=radial polarization (cylindrical cs),
                  ! 3=azimuthal polarization (cartesian cs),
                  ! 4= Circular polarization (cartesian cs).  Assuming RHC.  See fieldsModule.f95
                  !!!! also see ODEModule if changing coordinate systems for a polarization !!!!

        polyType=2,& ! 1=real argument, 2=complex argument (for associated Legendre polynomial definitions)
        rtType=2,& ! choice of R-tilde.  1 or 2.
        kmax1=floor(m/2.d0),&
        kmax2=floor((m-1.d0)/2.d0),&
        grid=200,& ! number of x,y grid cells.  for plotting and normalizing.  2D mesh is (grid x grid).
        ntraj=500000,&
        znum=1,& ! atomic number Z, for ionization from a gas
        electronNum=1,& ! number of electrons on each atom, for ionization from a gas
        writeSelection=0,& ! pick which files to create.  0 writes all fields, 1/2/3 writes only Ex/Ey/Ez, -1 writes no fields
        ionizationMode=0,& ! 0 for Joacain tunneling, 1 for SC tunneling + multiphoton
        injection=0,& ! 0 for no injection, 1 for injection
        fixSingularity=1,& ! correct for the numerical singularity at Rt=0? (in April's fields)
        whosFields=3,& ! 1=April, 2=BGV w/ Gaussian, 3=BGV w/ Poisson
            pertOrder=4,& ! max value of (2j).  pert order to include in BGVP field calc.
            ! check required order with myBGV.nb
            nnMax=nn+pertOrder,&
        plasmaWake=0,& ! include the electric fields of the plasma wake?
        normType=2 ! normalize by fixed peak power (1), or fixed pulse energy (2).  Set values in fields.95
    real*8,parameter::&
        phi0=0.d0*pi/4.d0,& ! initial phase
        lambda=800.d-9/a0,& ! laser wavelength (800 nm)
        s=70,&! cycles = s-value:          3=64, 5=178, 6=256, 7=349, 8=456, 9=577, 10=712
        w0=785.d-9/a0,&

        omega0=2.d0*pi*c/lambda,& ! laser frequency (angular)
        Tlas=lambda/c,& ! laser period
        TFWHM=2.d0*sqrt(2.d0*s*log(2.d0))/omega0*2.41888d-17,& ! in seconds
        RMSt=s/(omega0*sqrt(2.d0*s-1.d0)),& ! RMS of temporal distribution (see April, EQ. 23a - "will be used as the expression of the pulse duration")
        RMSc=s/(2.d0*pi*sqrt(2.d0*s-1.d0)),& ! RMS in cycles
        hpulse=((4.d0*sqrt(log(2.d0))*RMSc)*1.3d0+0.5d0)*Tlas,& ! Liangwen's half pulse duration estimate
        k=2.d0*pi/lambda,& ! wave number
        zr=(1.d0/2.d0)*k*w0**2,& ! Rayleigh range
        a=zr*sqrt(1.d0+2.d0/(k*zr)),& ! confocal parameter

        mReal=real(m,kind=8),&
        nnReal=real(nn,kind=8),&

        ! plasma variables, in Gaussian units
        EAUperECGS=5.83407d-8,&
        CGSme=9.1094d-28,&
        CGSqe=4.8032d-10,&
        CGSc=2.9979d10,&
        CGSlambda0pump=800.d-7,&! 800E-7cm=800nm
        CGSomega0=2.d0*pi*CGSc/CGSlambda0pump,&
        w0pump=15.d-4,& !15E-4cm...15E-6m
        omegap=CGSomega0/10.d0,&
        lambdap=2.d0*pi*CGSc/omegap,&
        betaphi=1.d0-0.5d0*omegap**2/CGSomega0**2,&
        a0hat=0.94d0,&

        bigPsiInjection=1.5d0*lambdap,&

    !----------for printing fields, and normalizing----------
        xmin=-3.d0*w0, ymin=xmin, xmax=-xmin, ymax=xmax,zmin=-2.d0*zr,zmax=-zmin,& ! grid dimensions
        dx=(xmax-xmin)/grid, dy=dx, dz=(zmax-zmin)/grid ! grid mesh size

    !----------MPI variables----------
    integer ierr,PID,ncores

    !----------other global variables----------
    integer nionized,randSeed,normed,test,plasmaOnly
    integer eval_num ! for the norm integrator
    real*8 pw,pulseEnergyJ,intensityWcm2
    real*8 sphericalPhi
    real*8 ws,norm,spr ! common factors in calculations
    real*8 bgj(0:nnMax),bc0(0:nn),bc1(0:nn+1),bc2(0:nn+2),&
        bc42(0:nn+2),bc43(0:nn+3),bc44(0:nn+4),&
        bc63(0:nn+3),bc64(0:nn+4),bc65(0:nn+5),bc66(0:nn+6) ! common factors in calculations (BGV_P)
    real*8 x,y,z,r,phi,t,r2,cp,sp!,& ! coordinates
        !distancePhase, oldPosition(1:3) ! for tracking phase information
    real*8 phasor(1:3),& ! for plotting the phasor (mod2, realpart, imagpart)
        ETotal ! (real) E-field magnitude
    real*8 psip,bigPsi ! plasma varaibles
    complex*16 Gpr, GprTilde, iwsc ! global constant phasor pieces
    complex*16 zia,zia2,Rt,Tm,Tp,crt,rrt
    complex*16 F0(0:kmax1),F1(0:kmax2),FA(0:kmax1),FB(0:kmax2),&
        FAr(0:kmax1),FBr(0:kmax2),FAz(0:kmax1),FBz(0:kmax2),FAt(0:kmax1),FBt(0:kmax2),&
        FArr(0:kmax1),FBrr(0:kmax2),FAzz(0:kmax1),FBzz(0:kmax2),&
        FArz(0:kmax1),FBrz(0:kmax2),FArt(0:kmax1),FBrt(0:kmax2),FAzt(0:kmax1),FBzt(0:kmax2)
    complex*16 Upr,Ur,Uz,Uprt,Urr,Uzz,Urz,Urt,Uzt
    complex*16 bConst, Cmn
    complex*16 ATilde,ATildeT
    complex*16 fields(1:6), plasmaFields(1:6) ! (Ex,Ey,Ez,Bx,By,Bz) for linear/azimuthal, (Er,Ephi,Ez,Br,Bphi,Bz) for radial

    !----------for DE solver----------
    integer nctt
    real*8,parameter::h1=1.d-2,hmin=0.d0,& ! initial and minimum step sizes
        eps=1.d-9,& ! This is the local error, given by the user.
        dt=Tlas/1.d3,& ! temporal grid for integrator
        !tmax=100.d0*Tlas ! to some number of laser cycles
        tmax=1.2d4/2.41888d0 ! 120 fs, as in Varin's 2013 results
        !tmax=2.75d2*Tlas
    integer,parameter::ntsteps=int(2.d0*hpulse/dt) ! number of time steps per pulse
    real*8 errmax

end module constants
