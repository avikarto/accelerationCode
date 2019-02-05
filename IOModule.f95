module IOModule
    use constants
    implicit none

    contains

    subroutine openFiles()
        implicit none

        if(writeSelection .ne. -1) then
            open(unit=11,file='fields/xs.txt')
            open(unit=12,file='fields/ys.txt')
            open(unit=13,file='fields/zs.txt')
        endif
        if (polar==1 .or. polar==3) then
            if(writeSelection==1.or.writeSelection==0) then
                open(unit=21,file='fields/ExReal.txt')
                open(unit=22,file='fields/EyReal.txt')
                open(unit=23,file='fields/EzReal.txt')
                open(unit=24,file='fields/ExImag.txt')
                open(unit=25,file='fields/EyImag.txt')
                open(unit=26,file='fields/EzImag.txt')
                !open(unit=27,file='fields/transvInten.txt')
            endif
            if(writeSelection==0) then
                !open(unit=31,file='fields/Bx.txt')
                !open(unit=32,file='fields/By.txt')
                !open(unit=33,file='fields/Bz.txt')
            endif
        elseif (polar==2) then
            if(writeSelection==1.or.writeSelection==0) then
                open(unit=21,file='fields/ErReal.txt')
                open(unit=22,file='fields/EpReal.txt')
                open(unit=23,file='fields/EzReal.txt')
                open(unit=24,file='fields/ErImag.txt')
                open(unit=25,file='fields/EpImag.txt')
                open(unit=26,file='fields/EzImag.txt')
                !open(unit=27,file='fields/transvInten.txt')
            endif
            if(writeSelection==0) then
                open(unit=31,file='fields/Br.txt')
                open(unit=32,file='fields/Bp.txt')
            endif
        endif
        if(ntraj>1) then
            open(unit=41,file='results/xpos.txt')
            open(unit=42,file='results/ypos.txt')
            open(unit=43,file='results/zpos.txt')
            open(unit=44,file='results/energyEV.txt')
            open(unit=45,file='results/xmom.txt')
            open(unit=46,file='results/ymom.txt')
            open(unit=47,file='results/zmom.txt')
            !open(unit=48,file='results/distancePhases.txt')
        endif
        open(unit=51,file='results/angleRad.txt')
        if(writeSelection==0) then
            open(unit=52,file='fields/phasorMod2.txt')
            !open(unit=53,file='fields/phasorReal.txt')
            !open(unit=54,file='fields/phasorImag.txt')
            !open(unit=55,file='fields/poynting.txt')
        endif
        if(ntraj>1) then
            open(unit=61,file='results/initialX.txt')
            open(unit=62,file='results/initialY.txt')
            open(unit=63,file='results/initialZ.txt')
        endif

        return
    end subroutine

    subroutine closeFiles()
        if (writeSelection .ne. -1) then
            close(11) ! xs
            close(12) ! ys
            close(13) ! zs
        endif
        if(writeSelection==1.or.writeSelection==0) then
            close(21) ! Ex/Er Real
            close(22) ! Ey/Ep Real
            close(23) ! Ez Real
            close(24) ! Ex/Er Imag
            close(25) ! Ey/Ep Imag
            close(26) ! Ez Imag
            !close(27) ! transverse intensity
        endif
        if(writeSelection==0) close(31) ! Bx/Br
        if(writeSelection==0) close(32) ! By/Bp
        !if((polar==1.and.writeSelection==0) .or. (polar==3.and.writeSelection==0)) close(33) ! Bz
        if(ntraj>1) then
            close(41) ! xpos
            close(42) ! ypos
            close(43) ! zpos
            close(44) ! energy
            close(45) ! xmom
            close(46) ! ymom
            close(47) ! zmom
            !close(48) ! phases
        endif
        close(51) ! angle
        if(writeSelection==0) then
            close(52) ! phasorMod2
            !close(53) ! phasorReal
            !close(54) ! phasorImag
            !close(55) ! poynting
        endif
        if(ntraj>1) then
            close(61)
            close(62)
            close(63)
        endif

        return
    end subroutine

    subroutine printParameters()
        if(polar==1) then
            write(*,*) "Laser parameters (linearly polarized):"
        elseif(polar==2) then
            write(*,*) "Laser parameters (radially polarized):"
        elseif(polar==3) then
            write(*,*) "Laser parameters (azimuthally polarized):"
        elseif(polar==4) then
            write(*,*) "Laser parameters (circularly polarized):"
        endif
        if(whosFields==1) write(*,*) "---April's Fields---"
        if(whosFields==2) write(*,*) "---BGV Fields w/ temporal Gaussian---"
        if(whosFields==3) write(*,*) "---BGV Fields w/ Poisson-like spectrum---"
        if(whosFields.ne.1) write(*,*) "perturbative terms retained (2j):",pertOrder
        write(*,*) "----------------------------------------------------"
        write(*,*) "Angular mode:",m
        if(whosFields==2.or.whosFields==3) write(*,*) "Radial mode:",nn
        write(*,*) "s=",s
        write(*,*) "phi0=",phi0/pi,"pi"
        write(*,*) "beam waist=",w0*a0*1.d6,"um=",w0/lambda,"lambda"
        write(*,*) "wavelength=",lambda*a0*1.d9,"nm"
        write(*,*) "frequency=",omega0,"a.u."
        write(*,*) "wave number k=",k,"au"
        if(whosFields==1) write(*,*) "confocal parameter=",a*a0*1.d6,"um=",a/lambda,"lambda"
        if(whosFields==1) write(*,*) "k*a=",k*a
        write(*,*) "zR=",zr*a0*1.d6,"um=",zr,"au=",zr/lambda,"lambda"
        write(*,*) "Intensity (W/cm^2) =",intensityWcm2
        write(*,*) "Field normalization calculated at",norm,"with",eval_num,"func calls."
        write(*,*)
        write(*,*) "----------------------------------------------------"
!        write(*,*) "Pulse duration calculations:"
        write(*,*) "T_fwhm =",TFWHM*1.d15,"fs"
        write(*,*) "T_fwhm =",2.d0*sqrt(2.d0*s*log(2.d0))/omega0/Tlas,"cycles"
        write(*,*) "Peak power (GW)=",pw*1.d-9
        write(*,*) "Pulse energy (J)=",pulseEnergyJ
!        write(*,*) "sigma_t=",RMSt*secPerAU*1.d15,"fs"
 !       write(*,*) "LWs=",2.d0*hpulse*secPerAU*1.d15,"fs"
!        if(rtType==1) then
!            if(whosFields==1) write(*,*) "R-tilde choice: 1"
!        else
!            if(whosFields==1) write(*,*) "R-tilde choice: 2"
!        endif
!        if(polyType==1) then
!            if(whosFields==1) write(*,*) "Associated Legendre polynomial choice: real argument"
!        else
!            if(whosFields==1) write(*,*) "Associated Legendre polynomial choice: complex argument"
!        endif
        write(*,*)
        write(*,*) "----------------------------------------------------"
        write(*,*)
        if(plasmaWake==1) then
            write(*,*) "w0pump=",w0pump,"cm"
            write(*,*) "CGS_omega_0/omega_plasma=",CGSomega0/omegap
            write(*,*) "beta_phi=",betaphi
            write(*,*) "a0hat=",a0hat
        endif
        write(*,*)
        return
    end subroutine

    subroutine printFields(selection)
        implicit none
        integer selection

        if (selection .ne. -1) then
            write(11,*) x
            write(12,*) y
            write(13,*) z
        endif
        if(selection==1.or.selection==0) write(21,"(ES12.5E2)") real(fields(1)) ! Ex or Er
        if(selection==2.or.selection==0) write(22,"(ES12.5E2)") real(fields(2)) ! Ey or Ep
        if(selection==3.or.selection==0) write(23,"(ES12.5E2)") real(fields(3)) ! Ez
        if(selection==1.or.selection==0) write(24,"(ES12.5E2)") aimag(fields(1)) ! Ex or Er
        if(selection==2.or.selection==0) write(25,"(ES12.5E2)") aimag(fields(2)) ! Ey or Ep
        if(selection==3.or.selection==0) write(26,"(ES12.5E2)") aimag(fields(3)) ! Ez
        !if(selection==3.or.selection==0) write(27,"(ES12.5E2)") sqrt( abs(fields(1))**2 + abs(fields(2))**2 ) ! Transverse Electric intensity (instant.)
        if(selection==0) write(31,"(ES12.5E2)") real(fields(4)) ! Bx or Br
        if(selection==0) write(32,"(ES12.5E2)") real(fields(5)) ! By or Bp
        !if((polar==1.or.polar==3).and.selection==0) write(33,"(ES12.5E2)") real(fields(6)) ! Bz
        if(selection==0) write(52,*) phasor(1) ! mod squared
        !if(selection==0) write(53,*) phasor(2) ! real
        !if(selection==0) write(54,*) phasor(3) ! imag
        !if(selection==0) write(55,"(ES12.5E2)") real( fields(1)*conjg(fields(5))-fields(2)*conjg(fields(4)) )*c2/(8.d0*pi) ! z-comp of complex poynting vec
        return
    end subroutine

end module IOModule
