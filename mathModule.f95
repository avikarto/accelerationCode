module mathModule
    use constants
    use mpi
    implicit none

    contains

    subroutine init_random_seed()

        INTEGER :: n, clock!,i
        INTEGER, DIMENSION(:), ALLOCATABLE :: seed

        call MPI_COMM_RANK(MPI_COMM_WORLD,PID,ierr) ! get process ID, PID

        CALL RANDOM_SEED(size = n)
        ALLOCATE(seed(n))

        CALL SYSTEM_CLOCK(COUNT=clock)

!        seed = clock + 37 * (/ (i - 1, i = 1, n) /)
        seed = abs( mod((clock*181)*((PID-83)*359), 104729) ) ! http://stackoverflow.com/questions/8920411/possible-sources-for-random-number-seeds
        CALL RANDOM_SEED(PUT = seed)

        DEALLOCATE(seed)
    end

    function Pmm(z,dr,dz) ! calculating Pmm(z,dr-count,dz-count)
        implicit none
        integer dr,dz
        complex*16 Pmm,part,amiz,amiz2,root
        real*8 z

        if(dr<0 .or. dz<0) then
            write(*,*) "Error: Pmm derivative order can't be negative."
            stop
        endif

        if(polyType.ne.1 .and. polyType.ne.2) then
            write(*,*) "Error: invalid choice between real/complex argument for Legendre.  See polyType variable."
            stop
        endif

        amiz=(a-ii*z)
        amiz2=amiz**2
        part=(r2-amiz2)
        root=sqrt(-1.d0/part)

        pmm=0.d0

        if (m==0) then ! real and complex valued
            if(dz==0 .and. dr==0) then
                Pmm=1.d0
            else
                Pmm=0.d0
            endif
        elseif (m==1 .and. polyType==1) then ! real valued
            if(dz==0 .and. dr==0) then
                Pmm=-sqrt(r2/part)
            elseif(dz==0 .and. dr==1) then
                Pmm=r*amiz2/( sqrt(r2/part)*part**2 )
            elseif(dz==1 .and. dr==0) then
                Pmm=sqrt(r2/part)**3*zia/r2
            elseif(dz==0 .and. dr==2) then
                Pmm=-3.d0*sqrt(r2/part)*amiz2/part**2
            elseif(dz==2 .and. dr==0) then
                Pmm=sqrt(r2/part)*(r2+2.d0*amiz2)/part**2
            elseif(dz==1 .and. dr==1) then
                Pmm=-ii*r*(2.d0*r2+amiz2)*amiz/( sqrt(r2/part)*part**3 )
            else
                write(*,*) "Invalid deriv counter in Pmm.  See mathModule."
                stop
            endif
        elseif (m==1 .and. polyType==2) then ! complex valued
            if(dz==0 .and. dr==0) then
                Pmm=r*root
            elseif(dz==0 .and. dr==1) then
                Pmm=root**3*amiz2
            elseif(dz==1 .and. dr==0) then
                Pmm=r*root**3*zia
            elseif(dz==0 .and. dr==2) then
                Pmm=3.d0*r*root*amiz2/part**2
            elseif(dz==2 .and. dr==0) then
                Pmm=-r*root**5*(r2+2.d0*amiz2)
            elseif(dz==1 .and. dr==1) then
                Pmm=root**5*(2.d0*r2+amiz2)*zia
            else
                write(*,*) "Invalid deriv counter in Pmm.  See mathModule."
                stop
            endif
        elseif (m==2 .and. polyType==1) then ! real valued
            if(dz==0 .and. dr==0) then
                Pmm=3.d0*r2/part
            elseif(dz==0 .and. dr==1) then
                Pmm=-6.d0*r*amiz2/part**2
            elseif(dz==1 .and. dr==0) then
                Pmm=-6.d0*ii*r2*amiz/part**2
            elseif(dz==0 .and. dr==2) then
                Pmm=6.d0*(3.d0*r2+amiz2)*amiz2/part**3
            elseif(dz==2 .and. dr==0) then
                Pmm=-6.d0*r2*(r2+3.d0*amiz2)/part**3
            elseif(dz==1 .and. dr==1) then
                Pmm=12.d0*r*(r2+amiz2)*zia/part**3
            else
                write(*,*) "Invalid deriv counter in Pmm.  See mathModule."
                stop
            endif
        elseif (m==2 .and. polyType==2) then ! complex valued
            if(dz==0 .and. dr==0) then
                Pmm=-3.d0*r2/part
            elseif(dz==0 .and. dr==1) then
                Pmm=6.d0*r*amiz2/part**2
            elseif(dz==1 .and. dr==0) then
                Pmm=6.d0*r2*zia/part**2
            elseif(dz==0 .and. dr==2) then
                Pmm=-6.d0*(3.d0*r2+amiz2)*amiz2/part**3
            elseif(dz==2 .and. dr==0) then
                Pmm=6.d0*r2*(r2+3.d0*amiz2)/part**3
            elseif(dz==1 .and. dr==1) then
                Pmm=-12.d0*r*(r2+amiz2)*zia/part**3
            else
                write(*,*) "Invalid deriv counter in Pmm.  See mathModule."
                stop
            endif
        elseif (m==3 .and. polyType==1) then ! real valued
            if(dz==0 .and. dr==0) then
                Pmm=-15.d0*sqrt(r2/part)**3
            elseif(dz==0 .and. dr==1) then
                Pmm=45.d0*r*amiz2*sqrt(r2/part)/part**2
            elseif(dz==1 .and. dr==0) then
                Pmm=45.d0*zia*sqrt(r2/part)**5/r2
            elseif(dz==0 .and. dr==2) then
                Pmm=-45.d0*sqrt(r2/part)*(3.d0*r2+2.d0*amiz2)*amiz2/part**3
            elseif(dz==2 .and. dr==0) then
                Pmm=45.d0*sqrt(r2/part)**3*(r2+4.d0*amiz2)/part**2
            elseif(dz==1 .and. dr==1) then
                Pmm=-45.d0*ii*r*sqrt(r2/part)*(2.d0*r2+3.d0*amiz2)*amiz/part**3
            else
                write(*,*) "Invalid deriv counter in Pmm.  See mathModule."
                stop
            endif
        elseif (m==3 .and. polyType==2) then ! complex valued
            if(dz==0 .and. dr==0) then
                Pmm=15.d0*sqrt(-r2/part)**3
            elseif(dz==0 .and. dr==1) then
                Pmm=45.d0*r2*root*amiz2/part**2
            elseif(dz==1 .and. dr==0) then
                Pmm=45.d0*r**3*root*zia/part**2
            elseif(dz==0 .and. dr==2) then
                Pmm=45.d0*r*root**7*(3.d0*r2+2.d0*amiz2)*amiz2
            elseif(dz==2 .and. dr==0) then
                Pmm=-45.d0*r**3*root**7*(r2+4.d0*amiz2)
            elseif(dz==1 .and. dr==1) then
                Pmm=45.d0*r2*root**7*(2.d0*r2+3.d0*amiz2)*zia
            else
                write(*,*) "Invalid deriv counter in Pmm.  See mathModule."
                stop
            endif
        elseif (m==4 .and. polyType==2) then ! complex valued
            if(dz==0 .and. dr==0) then
                Pmm=105.d0*r**4/part**2
            elseif(dz==0 .and. dr==1) then
                Pmm=-420.d0*r**3*amiz2/part**3
            elseif(dz==1 .and. dr==0) then
                Pmm=-420.d0*r**4*zia/part**3
            elseif(dz==0 .and. dr==2) then
                Pmm=1260.d0*r2*(r2+amiz2)*amiz2/part**4
            elseif(dz==2 .and. dr==0) then
                Pmm=-420.d0*r**4*(r2+5.d0*amiz2)/part**4
            elseif(dz==1 .and. dr==1) then
                Pmm=840.d0*r**3*(r2+2.d0*amiz2)*zia/part**4
            else
                write(*,*) "Invalid deriv counter in Pmm.  See mathModule."
                stop
            endif
        elseif (m==5 .and. polyType==2) then ! complex valued
            if(dz==0 .and. dr==0) then
                Pmm=945.d0*r**5*root
            elseif(dz==0 .and. dr==1) then
                Pmm=-4725.d0*root*r**4*amiz2/part**3
            elseif(dz==1 .and. dr==0) then
                Pmm=4725.d0*ii*root**7*r**5*amiz
            elseif(dz==0 .and. dr==2) then
                Pmm=4725.d0*root**9*r**3*(3.d0*r2+4.d0*amiz2)*amiz2
            elseif(dz==2 .and. dr==0) then
                Pmm=-4725.d0*root*r**5*(r2+6.d0*amiz2)/part**4
            elseif(dz==1 .and. dr==1) then
                Pmm=4725.d0*ii*root**9*r**4*(2.d0*r2+5.d0*amiz2)*amiz
            else
                write(*,*) "Invalid deriv counter in Pmm.  See mathModule."
                stop
            endif
        elseif (m==6 .and. polyType==2) then ! complex valued
            if(dz==0 .and. dr==0) then
                Pmm=-10395.d0*r2/part**3
            elseif(dz==0 .and. dr==1) then
                Pmm=62370.d0*r**5*amiz2/part**4
            elseif(dz==1 .and. dr==0) then
                Pmm=62370.d0*r**6*zia/part**4
            elseif(dz==0 .and. dr==2) then
                Pmm=-62370*r**4*(3.d0*r2+5.d0*amiz2)*amiz2/part**5
            elseif(dz==2 .and. dr==0) then
                Pmm=62370*r**6*(r2+7.d0*amiz2)/part**5
            elseif(dz==1 .and. dr==1) then
                Pmm=-124740.d0*r**5*(r2+3.d0*amiz2)*zia/part**5
            else
                write(*,*) "Invalid deriv counter in Pmm.  See mathModule."
                stop
            endif
        elseif (m==15 .and. polyType==2) then ! complex valued
            if(dz==0 .and. dr==0) then
                Pmm=6190283353629375.d0*sqrt(-r2/part)**15
            elseif(dz==0 .and. dr==1) then
                Pmm=92854250304440625.d0*r**14*root*amiz2/part**8
            elseif(dz==1 .and. dr==0) then
                Pmm=92854250304440625.d0*r**15*root*zia/part**8
            elseif(dz==0 .and. dr==2) then
                Pmm=92854250304440625.d0*r**13*root**19*(3.d0*r2+14.d0*amiz2)*amiz2
            elseif(dz==2 .and. dr==0) then
                Pmm=-92854250304440625.d0*r**15*root**19*(r2+16.d0*amiz2)
            elseif(dz==1 .and. dr==1) then
                Pmm=92854250304440625.d0*r**14*root**19*(2.d0*r2+15.d0*amiz2)*zia
            else
                write(*,*) "Invalid deriv counter in Pmm.  See mathModule."
                stop
            endif
        else
            write(*,*) "Invalid m value in Pmm.  See mathModule."
            stop
        endif
        if(polyType==1) pmm=pmm/ii**m ! pmm definitions from Mathematica were for real arguments. This translates to complex argument pmm's. (not really?)
        return
    end function

    function aterm(theta,m) ! a_2k(m+1/2) ,&c
        implicit none
        real*8 aterm
        integer theta,m

        aterm=fact(theta+m)/(2.d0**theta*fact(theta)*fact(m-theta))

        return
    end function

    function bAssLagPoly(bn,bm,bx) ! associated Laguerre polynomial L_bn^bm (bx)
        implicit none
        integer bn,bm
        real*8 rex,imx
        complex*16 bx,bAssLagPoly

        if (bn<0) then
            bAssLagPoly=(0.d0,0.d0)
            return
        elseif (bn==0) then
            bAssLagPoly=(1.d0,0.d0)
            return
        else

            rex=realpart(bx)
            imx=imagpart(bx)

            if(bm==0) then
                if(bn==1) then
                    bAssLagPoly=1.d0 - ii*imx - rex
                elseif(bn==2) then
                    bAssLagPoly=(2.d0 - 4.d0*ii*imx - imx**2 - 4.d0*rex + 2.d0*ii*imx*rex + rex**2)/2.d0
                elseif(bn==3) then
                    bAssLagPoly=(6.d0 - 18.d0*ii*imx - 9.d0*imx**2 + ii*imx**3 -18.d0*rex + &
                        18.d0*ii*imx*rex + 3.d0*imx**2*rex + 9.d0*rex**2 - &
                        3.d0*ii*imx*rex**2 - rex**3)/6.d0
                elseif(bn==4) then
                    bAssLagPoly=(24.d0 - 96.d0*ii*imx - 72.d0*imx**2 + 16.d0*ii*imx**3 + imx**4 - &
                        96.d0*rex + 144.d0*ii*imx*rex + 48.d0*imx**2*rex - 4.d0*ii*imx**3*rex + &
                        72.d0*rex**2 - 48.d0*ii*imx*rex**2 - 6.d0*imx**2*rex**2 - 16.d0*rex**3 + &
                        4.d0*ii*imx*rex**3 + rex**4)/24.d0
                elseif(bn==5) then
                    bAssLagPoly=(120.d0 - 600.d0*ii*imx - 600.d0*imx**2 + 200.d0*ii*imx**3 + 25.d0*imx**4 - &
                        ii*imx**5 - 600.d0*rex + 1200.d0*ii*imx*rex + 600.d0*imx**2*rex - 100.d0*ii*imx**3*rex - &
                        5.d0*imx**4*rex + 600.d0*rex**2 - 600.d0*ii*imx*rex**2 - 150.d0*imx**2*rex**2 + &
                        10.d0*ii*imx**3*rex**2 - 200.d0*rex**3 + 100.d0*ii*imx*rex**3 + 10.d0*imx**2*rex**3 + &
                        25.d0*rex**4 - 5.d0*ii*imx*rex**4 - rex**5)/120.d0
                else
                    write(*,*) "Undefined Laguerre poly index.  Stopping.",bn,bm
                    stop
                endif
            elseif (bm==1) then
                if(bn==1) then
                    bAssLagPoly=2.d0 - ii*imx - rex
                elseif(bn==2) then
                    bAssLagPoly=(6.d0 - 6.d0*ii*imx - imx**2 - 6.d0*rex + 2.d0*ii*imx*rex + rex**2)/2.d0
                elseif(bn==3) then
                    bAssLagPoly=(24.d0 - 36.d0*ii*imx - 12.d0*imx**2 + ii*imx**3 -36.d0*rex +&
                         24.d0*ii*imx*rex + 3.d0*imx**2*rex + 12.d0*rex**2 - 3.d0*ii*imx*rex**2 - rex**3)/6.d0
                elseif(bn==4) then
                    bAssLagPoly=(120.d0 - 240.d0*ii*imx - 120.d0*imx**2 + 20.d0*ii*imx**3 + imx**4 - 240.d0*rex + &
                        240.d0*ii*imx*rex + 60.d0*imx**2*rex - 4.d0*ii*imx**3*rex + 120.d0*rex**2 - &
                        60.d0*ii*imx*rex**2 - 6.d0*imx**2*rex**2 - 20.d0*rex**3 + 4.d0*ii*imx*rex**3 + rex**4)/24.d0
                elseif(bn==5) then
                    bAssLagPoly=(720.d0 - 1800.d0*ii*imx - 1200.d0*imx**2 + 300.d0*ii*imx**3 + 30.d0*imx**4 - &
                        ii*imx**5 - 1800.d0*rex + 2400.d0*ii*imx*rex + 900.d0*imx**2*rex - 120.d0*ii*imx**3*rex - &
                        5.d0*imx**4*rex + 1200.d0*rex**2 - 900.d0*ii*imx*rex**2 - 180.d0*imx**2*rex**2 + &
                        10.d0*ii*imx**3*rex**2 - 300.d0*rex**3 + 120.d0*ii*imx*rex**3 + 10.d0*imx**2*rex**3 + &
                        30.d0*rex**4 - 5.d0*ii*imx*rex**4 - rex**5)/120.d0
                else
                    write(*,*) "Undefined Laguerre poly index.  Stopping.",bn,bm
                    stop
                endif
            elseif (bm==2) then
                if(bn==1) then
                    bAssLagPoly=3.d0 - ii*imx - rex
                elseif(bn==2) then
                    bAssLagPoly=(12.d0 - 8.d0*ii*imx - imx**2 - 8.d0*rex + 2.d0*ii*imx*rex + rex**2)/2.d0
                elseif(bn==3) then
                    bAssLagPoly=(60.d0 - 60.d0*ii*imx - 15.d0*imx**2 + ii*imx**3 - 60.d0*rex + &
                        30.d0*ii*imx*rex + 3.d0*imx**2*rex + 15.d0*rex**2 - 3.d0*ii*imx*rex**2 - rex**3)/6.d0
                elseif(bn==4) then
                    bAssLagPoly=(360.d0 - 480.d0*ii*imx - 180.d0*imx**2 + 24.d0*ii*imx**3 + imx**4 - 480.d0*rex + &
                        360.d0*ii*imx*rex + 72.d0*imx**2*rex - 4.d0*ii*imx**3*rex + 180.d0*rex**2 - &
                        72.d0*ii*imx*rex**2 - 6.d0*imx**2*rex**2 - 24.d0*rex**3 + 4.d0*ii*imx*rex**3 + rex**4)/24.d0
                elseif(bn==5) then
                    bAssLagPoly=(2520.d0 - 4200.d0*ii*imx - 2100.d0*imx**2 + 420.d0*ii*imx**3 + 35.d0*imx**4 - &
                        ii*imx**5 - 4200.d0*rex + 4200.d0*ii*imx*rex + 1260.d0*imx**2*rex - 140.d0*ii*imx**3*rex - &
                        5.d0*imx**4*rex + 2100.d0*rex**2 - 1260.d0*ii*imx*rex**2 - 210.d0*imx**2*rex**2 + &
                        10.d0*ii*imx**3*rex**2 - 420.d0*rex**3 + 140.d0*ii*imx*rex**3 + 10.d0*imx**2*rex**3 + &
                        35.d0*rex**4 - 5.d0*ii*imx*rex**4 - rex**5)/120.d0
                else
                    write(*,*) "Undefined Laguerre poly index.  Stopping.",bn,bm
                    stop
                endif
            elseif (bm==3) then
                if(bn==1) then
                    bAssLagPoly=4.d0 - ii*imx - rex
                elseif(bn==2) then
                    bAssLagPoly=(20.d0 - 10.d0*ii*imx - imx**2 - 10.d0*rex + 2.d0*ii*imx*rex + rex**2)/2.d0
                elseif(bn==3) then
                    bAssLagPoly=(120.d0 - 90.d0*ii*imx - 18.d0*imx**2 + ii*imx**3 - 90.d0*rex + &
                        36.d0*ii*imx*rex + 3.d0*imx**2*rex + 18.d0*rex**2 - 3.d0*ii*imx*rex**2 - rex**3)/6.d0
                elseif(bn==4) then
                    bAssLagPoly=(840.d0 - 840.d0*ii*imx - 252.d0*imx**2 + 28.d0*ii*imx**3 + imx**4 - 840.d0*rex + &
                        504.d0*ii*imx*rex + 84.d0*imx**2*rex - 4.d0*ii*imx**3*rex + 252.d0*rex**2 - &
                        84.d0*ii*imx*rex**2 - 6.d0*imx**2*rex**2 - 28.d0*rex**3 + 4.d0*ii*imx*rex**3 + rex**4)/24.d0
                elseif(bn==5) then
                    bAssLagPoly=(6720.d0 - 8400.d0*ii*imx - 3360.d0*imx**2 + 560.d0*ii*imx**3 + 40.d0*imx**4 - &
                        ii*imx**5 - 8400.d0*rex + 6720.d0*ii*imx*rex + 1680.d0*imx**2*rex - 160.d0*ii*imx**3*rex - &
                        5.d0*imx**4*rex + 3360.d0*rex**2 - 1680.d0*ii*imx*rex**2 - 240.d0*imx**2*rex**2 + &
                        10.d0*ii*imx**3*rex**2 - 560.d0*rex**3 + 160.d0*ii*imx*rex**3 + 10.d0*imx**2*rex**3 + &
                        40.d0*rex**4 - 5.d0*ii*imx*rex**4 - rex**5)/120.d0
                else
                    write(*,*) "Undefined Laguerre poly index.  Stopping.",bn,bm
                    stop
                endif
            elseif (bm==4) then
                if(bn==1) then
                    bAssLagPoly=5.d0 - ii*imx - rex
                elseif(bn==2) then
                    bAssLagPoly=(30.d0 - 12.d0*ii*imx - imx**2 - 12.d0*rex + 2.d0*ii*imx*rex + rex**2)/2.d0
                elseif(bn==3) then
                    bAssLagPoly=(210.d0 - 126.d0*ii*imx - 21.d0*imx**2 + ii*imx**3 - 126.d0*rex + 42.d0*ii*imx*rex + &
                        3.d0*imx**2*rex + 21.d0*rex**2 - 3.d0*ii*imx*rex**2 - rex**3)/6.d0
                elseif(bn==4) then
                    bAssLagPoly=(1680.d0 - 1344.d0*ii*imx - 336.d0*imx**2 + 32.d0*ii*imx**3 + imx**4 - 1344.d0*rex + &
                        672.d0*ii*imx*rex + 96.d0*imx**2*rex - 4.d0*ii*imx**3*rex + 336.d0*rex**2 - &
                        96.d0*ii*imx*rex**2 - 6.d0*imx**2*rex**2 - 32.d0*rex**3 + 4.d0*ii*imx*rex**3 + rex**4)/24.d0
                elseif(bn==5) then
                    bAssLagPoly=(15120.d0 - 15120.d0*ii*imx - 5040.d0*imx**2 + 720.d0*ii*imx**3 + 45.d0*imx**4 - &
                        ii*imx**5 - 15120.d0*rex + 10080.d0*ii*imx*rex + 2160.d0*imx**2*rex - 180.d0*ii*imx**3*rex - &
                        5.d0*imx**4*rex + 5040.d0*rex**2 - 2160.d0*ii*imx*rex**2 - 270.d0*imx**2*rex**2 + &
                        10.d0*ii*imx**3*rex**2 - 720.d0*rex**3 + 180.d0*ii*imx*rex**3 + 10.d0*imx**2*rex**3 + &
                        45.d0*rex**4 - 5.d0*ii*imx*rex**4 - rex**5)/120.d0
                else
                    write(*,*) "Undefined Laguerre poly index.  Stopping.",bn,bm
                    stop
                endif
            elseif (bm==5) then
                if(bn==1) then
                    bAssLagPoly=6.d0 - ii*imx - rex
                elseif(bn==2) then
                    bAssLagPoly=(42.d0 - 14.d0*ii*imx - imx**2 - 14.d0*rex + 2.d0*ii*imx*rex + rex**2)/2.d0
                elseif(bn==3) then
                    bAssLagPoly=(336.d0 - 168.d0*ii*imx - 24.d0*imx**2 + ii*imx**3 - 168.d0*rex + 48.d0*ii*imx*rex + &
                        3.d0*imx**2*rex + 24.d0*rex**2 - 3.d0*ii*imx*rex**2 - rex**3)/6.d0
                elseif(bn==4) then
                    bAssLagPoly=(3024.d0 - 2016.d0*ii*imx - 432.d0*imx**2 + 36.d0*ii*imx**3 + imx**4 - 2016.d0*rex + &
                        864.d0*ii*imx*rex + 108.d0*imx**2*rex - 4.d0*ii*imx**3*rex + 432.d0*rex**2 - &
                        108.d0*ii*imx*rex**2 - 6.d0*imx**2*rex**2 - 36.d0*rex**3 + 4.d0*ii*imx*rex**3 + rex**4)/24.d0
                elseif(bn==5) then
                    bAssLagPoly=(30240.d0 - 25200.d0*ii*imx - 7200.d0*imx**2 + 900.d0*ii*imx**3 + 50.d0*imx**4 - &
                        ii*imx**5 - 25200.d0*rex + 14400.d0*ii*imx*rex + 2700.d0*imx**2*rex - 200.d0*ii*imx**3*rex - &
                        5.d0*imx**4*rex + 7200.d0*rex**2 - 2700.d0*ii*imx*rex**2 - 300.d0*imx**2*rex**2 + &
                        10.d0*ii*imx**3*rex**2 - 900.d0*rex**3 + 200.d0*ii*imx*rex**3 + 10.d0*imx**2*rex**3 + &
                        50.d0*rex**4 - 5.d0*ii*imx*rex**4 - rex**5)/120.d0
                else
                    write(*,*) "Undefined Laguerre poly index.  Stopping.",bn,bm
                    stop
                endif
            else
                write(*,*) "Stopping.  Invalid or undefined Laguerre index M=",bm
                stop
            endif !bm==?
        endif ! bn<=0 or else?
        return
    end function

    function fact(n) ! factorial
        implicit none
        integer fact,n,i

        fact = 1

        do i=n,1,-1
            fact = fact*i
        enddo

        return
    end function

    function fac2 (n) ! double factorial
          implicit none

          integer n
          double precision fac2
          double precision r8_n

          if ( n .lt. 1 ) then
            fac2 = 1.0D+00
            return
          end if

          r8_n = dble ( n )
          fac2 = 1.0D+00

990       continue

          if ( 1.0D+00 .lt. r8_n ) then
            fac2 = fac2 * r8_n
            r8_n = r8_n - 2.0D+00
            go to 990
          end if

          return
    end function

    function bgamma(jjj) ! phasor p.3, c
        implicit none
        real*8 bgamma
        integer jjj

        bgamma=mReal/2.d0+s+real(jjj,kind=8)

        return
    end function

    function fbc0(jjj) ! phasor p.3, c
        implicit none
        real*8 fbc0
        integer jjj

        fbc0=bgabj(nn,m,jjj)*ws**(bgj(jjj)-s)*&
            exp( log_gamma(bgj(jjj)+1.d0)-Log_gamma(s+1.d0) ) ! phasor p.3&4, c

        return
    end function

    function fbc1(jjj) ! phasor p.3, c
        implicit none
        real*8 fbc1
        integer jjj

        fbc1=c*(nnReal+1.d0)*bgabj(nn+1,m,jjj)*ws**(bgj(jjj)-s-1.d0)*&
            exp( log_gamma(bgj(jjj))-Log_gamma(s+1.d0) ) ! phasor p.2&3&4, c

        return
    end function

    function fbc2(jjj) ! phasor p.3, c
        implicit none
        real*8 fbc2
        integer jjj

        fbc2=(c/2.d0)*(nnReal+1.d0)*(nnReal+2.d0)*bgabj(nn+2,m,jjj)*ws**(bgj(jjj)-s-1.d0)*&
            exp( log_gamma(bgj(jjj))-Log_gamma(s+1.d0) ) ! phasor p.2&3&4, c

        return
    end function

    function fbc42(jjj) ! 4th order constant, in the sum to (n+2)
        implicit none
        real*8 fbc42
        integer jjj

        fbc42=(3.d0*c2/2.d0)*(nnReal+1.d0)*(nnReal+2.d0)*bgabj(nn+2,m,jjj)*ws**(bgj(jjj)-s-2.d0)*&
            exp( log_gamma(bgj(jjj)-1.d0)-Log_gamma(s+1.d0) )

        return
    end function

    function fbc43(jjj) ! 4th order constant, in the sum to (n+3)
        implicit none
        real*8 fbc43
        integer jjj

        fbc43=c2*(nnReal+1.d0)*(nnReal+2.d0)*(nnReal+3.d0)*bgabj(nn+3,m,jjj)*ws**(bgj(jjj)-s-2.d0)*&
            exp( log_gamma(bgj(jjj)-1.d0)-Log_gamma(s+1.d0) )

        return
    end function

    function fbc44(jjj) ! 4th order constant, in the sum to (n+4)
        implicit none
        real*8 fbc44
        integer jjj

        fbc44=(c2/8.d0)*(nnReal+1.d0)*(nnReal+2.d0)*(nnReal+3.d0)*(nnReal+4.d0)*bgabj(nn+4,m,jjj)*ws**(bgj(jjj)-s-2.d0)*&
            exp( log_gamma(bgj(jjj)-1.d0)-Log_gamma(s+1.d0) )

        return
    end function

    function fbc63(jjj)
        implicit none
        real*8 fbc63
        integer jjj

        fbc63=(5.d0*c3/2.d0)*(nnReal+1.d0)*(nnReal+2.d0)*(nnReal+3.d0)*bgabj(nn+3,m,jjj)*ws**(bgj(jjj)-s-3.d0)*&
            exp( log_gamma(bgj(jjj)-2.d0)-Log_gamma(s+1.d0) )

        return
    end function

    function fbc64(jjj)
        implicit none
        real*8 fbc64
        integer jjj

        fbc64=(15.d0*c3/8.d0)*(nnReal+1.d0)*(nnReal+2.d0)*(nnReal+3.d0)*(nnReal+4.d0)*bgabj(nn+4,m,jjj)*ws**(bgj(jjj)-s-3.d0)*&
            exp( log_gamma(bgj(jjj)-2.d0)-Log_gamma(s+1.d0) )

        return
    end function

    function fbc65(jjj)
        implicit none
        real*8 fbc65
        integer jjj

        fbc65=(3.d0*c3/8.d0)*(nnReal+1.d0)*(nnReal+2.d0)*(nnReal+3.d0)*(nnReal+4.d0)*(nnReal+5.d0)*bgabj(nn+5,m,jjj)*ws**(bgj(jjj)-s-3.d0)*&
            exp( log_gamma(bgj(jjj)-2.d0)-Log_gamma(s+1.d0) )

        return
    end function

    function fbc66(jjj)
        implicit none
        real*8 fbc66
        integer jjj

        fbc66=(c3/48.d0)*(nnReal+1.d0)*(nnReal+2.d0)*(nnReal+3.d0)*(nnReal+4.d0)*(nnReal+5.d0)*(nnReal+6.d0)*bgabj(nn+6,m,jjj)*ws**(bgj(jjj)-s-3.d0)*&
            exp( log_gamma(bgj(jjj)-2.d0)-Log_gamma(s+1.d0) )

        return
    end function

    function fbcad(bcalpha,bcdelta,jjj) ! general BGVP coefficients c_{\alpha,\delta}
        implicit none
        real*8 fbcad
        integer jjj,bcalpha,bcdelta,chi

        if(bcalpha==0) then ! the kronecker delta term
            fbcad=1.d0
        else
            fbcad=2.d0
        endif
        fbcad=fbcad*fact(bcalpha)/(fact(bcdelta-bcalpha)*fact(2*bcalpha-bcdelta)) ! the binomial coefficient
        fbcad=fbcad/fact(bcdelta)*(-1.d0)**(bcalpha+bcdelta)*omega0**bcalpha*bgabj(nn+bcdelta,m,jjj)*&
            ws**(bcalpha-s+bgj(jjj))*exp( log_gamma(bgj(jjj)+1.d0-bcalpha)-Log_gamma(s+1.d0) )
        do chi=2,bcalpha
            fbcad=fbcad*(4.d0*real(chi,kind=8)-2.d0)
        enddo
        do chi=1,bcdelta
            fbcad=fbcad*(nnReal+real(chi,kind=8))
        enddo

        return
    end function

    function bgabj(aaa,bbb,jjj) ! phasor p.2, c
        implicit none
        real*8 bgabj
        integer aaa,bbb,jjj

        bgabj=((-1.d0)**jjj) * real(fact(aaa+bbb),kind=8)/real(fact(aaa-jjj)*fact(bbb+jjj)*fact(jjj),kind=8)

        return
    end function

end module mathModule
