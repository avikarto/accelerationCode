! Code for calculating fields from April's [1] phasor via Hertz potentials for a (p=0,m) L-G beam
! [1] A. April, in Coherence Ultrashort Pulse Laser Emiss.
!
! Author: Andrew Vikartofsky, May 2016 (avikarto@gmail.com)

program main
    use mathModule
    use constants
    use IOModule
    use fieldsModule
    use normModule
    use ODEModule
    use mpi
    implicit none
    integer ix,iy,iz
    integer i,status(MPI_STATUS_SIZE),nsent,nrec,temp1,temp2
    real*8 ppp(1:6),& ! particle position and momentum (X0x,X0y,X0z,P0x,...)
        energy,p2,pPerp,gama,angle,r0Temp(1:3)
    logical ionized
    real*8 injAngle,injEeVrnd,Emax,injPosZ
    real*8 injEeV, injEHart
!----------begin MPI workflow----------
    call MPI_INIT(ierr) ! start MPI
    call MPI_COMM_RANK(MPI_COMM_WORLD,PID,ierr) ! get process ID, PID
    call MPI_COMM_SIZE(MPI_COMM_WORLD,ncores,ierr) ! get total number of processes

    call init_random_seed()
    test=0

    if(normType==1) then
        pw=300.d9
        pulseEnergyJ=TFWHM*pw
    elseif(normType==2) then
        pulseEnergyJ=2.5d-3!8.5d0
        pw=pulseEnergyJ/TFWHM
    else
        write(*,*) "Invalid norm type.  Stopping."
        stop
    endif
    intensityWcm2=pw/(pi*(w0*a0*100.d0)**2) ! W/cm^2 (100 cm/m)

    ws=omega0/s
    !----------global constant phasor pieces----------
    if(whosFields==1) then ! April
        Gpr=exp(ii*phi0)*4.d0/fac2(2*m-1)*(a/(2.d0*c))**(1.d0+mReal/2.d0)
        iwsc=ii*ws/c
        spr=s+1.d0+mReal/2.d0
        GprTilde=exp(ii*phi0)*2.d0*a/c*(cmplx(-a/(2.d0*c),kind=16))**(mReal/2.d0)/(c**m*fac2(2*m+1))*exp( log_gamma(spr+mReal+1.d0)-log_gamma(s+1.d0) )
        call makeF01() ! calculate the sum coefficients(Fs) for the fields
    elseif(whosFields==2) then ! BGV_TG
        bConst=(-1.d0)**(nn+m)*2.d0**(2*nn+m)*exp(ii*phi0)
    elseif(whosFields==3) then ! BGV_P
        Cmn=cmplx(sqrt(2.d0*pi)*(-1.d0)**(nn+m)*2.d0**(2*nn+m),kind=16)*exp(ii*phi0)*cmplx(fact(nn),kind=16) ! phasor p. 1&3 ,c
        call makeBGVPvars()
    else
        write(*,*) "Invalid whosFields choice"
        stop
    endif !whosFields
    z=1.d-5
    t=z/c !real time
    norm=1.d0
    normed=0
    norm=normalize() ! calculate normalization factor for fields.
    normed=1

    if(PID==0) then ! the master process
        call printParameters()
!        if (injection==1) write(*,*) "Injection Energy (eV):",injEeV

        if (polar.lt.1 .or. polar.gt.4) then ! checking polarization settings
            write(*,*) "ERROR!! INVALID POLARIZATION!!"
            stop
        endif

        call openFiles()

        nsent=0
        nionized=0
        do while (nsent<ntraj)
            temp1=0 ! for knowing when to cycle back to slave 1
            do i=1,ncores-1 ! sends one job to each slave
                if(nsent<ntraj) then
                    if (injection==0) then
                        ! initial postions x,y,z of ions in range ...
                        if(ntraj==1) then
                            ppp(1)=0.57d0*lambda
                            ppp(2)=0.d0!1.d-5
                            ppp(3)=2.91d0*lambda
                        else ! ntraj>1
                            ! random positions
                            call random_number(ppp(4))
                            call random_number(ppp(5))
                            call random_number(ppp(6))
                            ppp(1)=( ppp(4)*2.d0 - 1.d0 )*lambda
                            ppp(2)=0.d0!( ppp(5)*2.d0 - 1.d0 )*lambda
                            ppp(3)=( ppp(6)*6.d0 - 3.d0 )*lambda
                        endif
                    else ! injection=1
                        ! give initial postions x,y,z for injected electrons
                        if (ntraj==1) then
                            ppp(1)=0.1d0*w0
                            ppp(2)=0.1d0*w0
                            ppp(3)=-1.d0*w0
                        else ! ntraj > 1 and injection=1
                            call random_number(ppp(4))
                            call random_number(ppp(5))
                            call random_number(ppp(6))
                            ppp(1)=2.d0*w0 + (2.d0*ppp(4)-1.d0)*0.05d0*w0
                            ppp(2)=1.d-5 + (2.d0*ppp(5)-1.d0)*0.05d0*w0
                            ppp(3)=1.d-5 + (2.d0*ppp(6)-1.d0)*0.05d0*w0
                            ! write initial ensemble positions
                            write(61,*) ppp(1)
                            write(62,*) ppp(2)
                            write(63,*) ppp(3)
                        endif ! ntraj>1 and injection == 1
                    endif ! injection == 0/1
                    call mpi_send(ppp(1:3),3,mpi_real8,i,1,mpi_comm_world,ierr) ! sending work to slaves: tag 1
                    nsent=nsent+1
                    temp1=temp1+1
                else
                    call mpi_send(ppp(1:3),3,mpi_real8,i,99,mpi_comm_world,ierr) ! sending kill code to slaves: tag 99
                endif
            enddo ! i,1,ncores-1
            temp2=0 ! for knowing when to stop collecting and start sending next set of jobs
            do i=1,ncores-1 ! collect results from all slaves and write
                if (temp2<temp1) then
                    call mpi_recv(ppp,6,mpi_real8,i,mpi_any_tag,mpi_comm_world,status,ierr) ! getting results from slaves: tag 2 or 3
                    if (ntraj.gt.1) then
                        call mpi_recv(r0Temp,3,mpi_real8,i,mpi_any_tag,mpi_comm_world,status,ierr) ! getting results from slaves: tag 4
                    endif

                    if(status(mpi_tag)==2) then ! if ionized
                        nionized=nionized+1

                        if( ppp(3).ge.0.d0 ) then ! zf filtering
                        !if( abs(ppp(1))/lambda.le.2.d0 ) then ! xf filtering
                        !if( abs(ppp(2))/lambda.le.0.5d0 ) then ! yf filtering

                            p2=ppp(4)**2+ppp(5)**2+ppp(6)**2
                            pPerp=sqrt(ppp(4)**2+ppp(5)**2) ! P_perp = sqrt( P_x^2 + P_y^2 )
                            angle=atan(pPerp/ppp(6)) ! ejection angle = atan( P_perp / P_z )
                            gama=sqrt(1.d0+p2/c2)
                            ! E = sqrt( m^2 c^4 + p^2 c^2 )
                            energy=gama*c2*eVperHartree

                            if(ntraj>1) then
                                write(61,*) r0Temp(1) ! initial x
                                write(62,*) r0Temp(2) ! initial y
                                write(63,*) r0Temp(3) ! initial z
                                write(41,*) ppp(1) ! final x
                                write(42,*) ppp(2) ! final y
                                write(43,*) ppp(3) ! final z
                                write(44,*) energy ! total energy, includes rest energy
                                write(45,*) ppp(4) ! final Px
                                write(46,*) ppp(5) ! final Py
                                write(47,*) ppp(6) ! final Pz
                                !write(48,*) distancePhase
                            endif
                            write(51,*) angle

                        !endif ! yf filtering
                        !endif ! xf filtering
                        endif ! zf filtering
                    endif ! ionized

                    temp2=temp2+1
                endif
            enddo
        enddo ! while loop

        ! write the fields to their respective files
!        Emax=0.d0
!        do iz=0,2*grid
!            z=-3.d0*lambda + real(iz,kind=8)*(6.d0*lambda/(2.d0*real(grid,kind=8)))
!            z=1.d-5
!            y=0.d0
!            do ix=0,grid
!                do iy=0,grid
!                    y=ymin+iy*dy
!                    x=0.d0+real(ix,kind=8)*(1.5d0*lambda/real(grid,kind=8))
!                    t=z/c !real time
!                    test=1
!                    call makeFields(x,y,z,t)
!                    call printFields(writeSelection)
!                    if (Etotal>Emax) Emax=Etotal
!                    test=0
!                enddo !y-loop
!            enddo !x-loop
!        enddo !z loop
!        write(*,*) "Emax=",Emax,"au"
!        if (injection==0) then
!            write(*,*) real(nionized)/real(ntraj)*1.d2,"% of electrons were ionized for Z=",znum
!            write(*,*) nionized,"were ionized out of",ntraj
!            write(*,*) "----------------------------------------------------"
!        endif
        call closeFiles()
        write(*,*) "Tmax (cycles)=",tmax/Tlas
        write(*,*) "Tmax (fs)=",tmax*2.41888d-2

!-----------------------------------------------------------
!----------end master process, start slave process----------
!-----------------------------------------------------------

    else ! a slave process
        if(ntraj==1) then
            open(unit=61,file='results/trajPosX.txt')
            open(unit=62,file='results/trajPosY.txt')
            open(unit=63,file='results/trajPosZ.txt')
            open(unit=64,file='results/trajPosT.txt')
            open(unit=65,file='results/pxt.txt')
            open(unit=66,file='results/pyt.txt')
            open(unit=67,file='results/pzt.txt')
            open(unit=71,file='results/EtotAlongTraj_eV.txt')
            open(unit=91,file='results/ext.txt')
            open(unit=92,file='results/eyt.txt')
            open(unit=93,file='results/ezt.txt')
            open(unit=94,file='results/phasort.txt')
            open(unit=95,file='results/bxt.txt')
            open(unit=96,file='results/byt.txt')
            open(unit=97,file='results/bzt.txt')

        endif
        nrec=0
        do while(nrec<ntraj)
            call mpi_recv(ppp,3,mpi_real8,0,mpi_any_tag,mpi_comm_world,status,ierr)
            if(status(mpi_tag)==99) then ! kill code recieved from master: tag 99
!                write(*,*) PID,"reached finalize by getting a kill command"
                call MPI_FINALIZE(ierr)
                stop
            elseif(status(mpi_tag)==1) then ! get inital position from master: tag 1

                !distancePhase=0.d0
                !oldPosition(1:3)=ppp(1:3)
                r0Temp(1)=ppp(1)
                r0Temp(2)=ppp(2)
                r0Temp(3)=ppp(3)
                nrec=nrec+ncores-1
                if(injection==0) then ! electron source is ionization
                    ppp(4:6)=1.d-5 ! zero initial momentum, next line also assumes this
                    ppp(6)=sqrt( ( (0.511d6/eVperHartree)**2-c**4 )/c2 ) ! temp storage of total
                else ! electron source is external injection (injection==1)

                    injPosZ=zr
                    injEeV=1.d6

                    ! spherical angle off of Z
                    injAngle=atan(ppp(1),abs(ppp(3)-injPosZ)) ! distance from z0 to injPosZ

                    ! spherical angle in X-Y plane
                    sphericalPhi=atan(ppp(2),ppp(1))

                    injEeVrnd=injEeV
                    if (ntraj>1) then ! and injection = 1
                        call random_number(injEeVrnd)
                        injEeVrnd=(2.d0*injEeVrnd-1.d0)*(0.215d0*InjEeV)+injEeV
                    endif
                    ppp(6)=sqrt( (( (injEeVrnd+0.511d6)/eVperHartree )**2-c**4)/c2 ) ! temporary storage of total momentum
                    if (ntraj==1) write(71,*) sqrt(1.d0+ppp(6)**2/c2)*c2*eVperHartree ! write total energy (eV)

                    ppp(4)=-ppp(6)*sin(injAngle)*cos(sphericalPhi)
                    ppp(5)=-ppp(6)*sin(injAngle)*sin(sphericalPhi)
                    ppp(6)=( -ppp(5) / (sin(injAngle)*sin(sphericalPhi)) ) * cos(injAngle) !=ppp(6)*cos(injAngle)
                endif ! injection == ?

                ! do ionization/trajectory calculations

                !t=ppp(3)/c-hpulse ! initial location (in time) of pulse's peak amplitude
                t=ppp(3)/c

                if(plasmaWake==1) t=t-bigPsiInjection/CGSc ! delay the injection pulse.  Pump@z=0@t=0.

                ionized=.false.
                iy=3 ! not ionized tag
                do i=1,ntsteps
                    !t=t+dt
                    !if (injection==0) call checkIonization(ppp(1:3),t,ionized)
                    ionized=.true. ! to force immediate ionization
                    if (ionized .or. injection==1) then
                        iy=2 ! ionized tag
                        if (ntraj==1) then
                            write(61,*) ppp(1)
                            write(62,*) ppp(2)
                            write(63,*) ppp(3)
                            write(65,*) ppp(4)
                            write(66,*) ppp(5)
                            write(67,*) ppp(6)
                            write(64,*) t
                            call makeFields(ppp(1),ppp(2),ppp(3),t)
                            write(71,*) sqrt( 1 + ppp(6)**2/c2 )*c2*eVperHartree
                            write(91,*) Real(fields(1))
                            write(92,*) Real(fields(2))
                            write(93,*) Real(fields(3))
                            write(95,*) Real(fields(4))
                            write(96,*) Real(fields(5))
                            write(97,*) Real(fields(6))
                        endif
                        ppp(6)=1.d-5 ! reset ppp(6) to P_z instead of P_tot
                        call odeint(ppp,6,t,tmax,eps,h1,hmin,errmax,ntraj)
                        exit
                    endif
                enddo ! time loop
                call mpi_send(ppp,6,mpi_real8,0,iy,mpi_comm_world,ierr) ! returning results to master: tag 2 or 3
                if (ntraj.gt.1) then
                    call mpi_send(r0Temp,3,mpi_real8,0,iy,mpi_comm_world,ierr) ! returning to master: tag 2
                endif
            endif ! mpi tag check
        enddo ! while true
        if(ntraj==1) then ! injection trajectories
            close(61)
            close(62)
            close(63)
            close(64)
            close(71)
        endif
    endif ! slave

    call MPI_FINALIZE(ierr)
end program
