module normModule
    use constants
    use fieldsModule
    implicit none

    contains

    function normalize()
        implicit none
        integer xygrid(1:2),iterationDepth,err,id
        real*8 mins(1:2),maxs(1:2),errTol,normalize

        errTol=2.d-3 ! error tolerance

        xygrid=grid

        mins(1)=xmin
        mins(2)=ymin
        maxs(1)=xmax
        maxs(2)=ymax

        do id=10,10
            iterationDepth=id
            call romberg_nd (mins,maxs,2,xygrid,iterationDepth,errTol,normalize,err,eval_num)
            !write(*,*) "raw norm:",normalize
    !           (1 au energy unit in J)/(1 au time unit in s)=(1 au power unit in J/s)=0.180238 [W per au]
            normalize=normalize*0.180238d0 ! (SI power before normalizing fields)=(power in au before normalizing fields)*(W per au)

            ! Power~Intesnsity~|E|^2, so multiply fields by normalize=sqrt(goalPower/calculatedPower)
            normalize=sqrt(pw/normalize)

            if(err==-1) then
                write(*,*) "Normalization failed: the romberg integral's error tolerance could not be achieved for depth",iterationDepth
    !                stop
            endif
            !write(*,*) "scaled norm=:",normalize
        enddo

        return
    end function

subroutine romberg_nd ( aaa, bbb, dim_num, sub_num, it_max, tol, result, &
  ind, eval_num )

!*****************************************************************************80
!
!! ROMBERG_ND estimates a multidimensional integral using Romberg integration.
!
!  Discussion:
!
!    The routine uses a Romberg method based on the midpoint rule.
!
!    In the reference, this routine is called "NDIMRI".
!
!    Thanks to Barak Bringoltz for pointing out problems in a previous
!    FORTRAN90 implementation of this routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 September 2006
!
!  Author:
!
!    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, real ( kind = 8 ), external FUNC, a routine which evaluates
!    the function to be integrated, of the form:
!      function func ( dim_num, x )
!      integer ( kind = 4 ) dim_num
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x(dim_num)
!      func = ...
!      return
!      end
!
!    Input, real ( kind = 8 ) AAA(DIM_NUM), BBB(DIM_NUM), the integration limits.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) SUB_NUM(DIM_NUM), the number of subintervals
!    into which the I-th integration interval (AAA(I), BBB(I)) is
!    initially subdivided.  SUB_NUM(I) must be greater than 0.
!
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations to
!    be performed.  The number of function evaluations on
!    iteration J is at least J**DIM_NUM, which grows very rapidly.
!    IT_MAX should be small!
!
!    Input, real ( kind = 8 ) TOL, an error tolerance for the approximation
!    of the integral.
!
!    Output, real ( kind = 8 ) RESULT, the approximate value of the integral.
!
!    Output, integer ( kind = 4 ) IND, error return flag.
!    IND = -1 if the error tolerance could not be achieved.
!    IND = 1 if the error tolerance was achieved.
!
!    Output, integer ( kind = 4 ) EVAL_NUM, the number of function evaluations.
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) IWORK(DIM_NUM), a pointer used to generate
!    all the points XX in the product region.
!
!    Local, integer ( kind = 4 ) IWORK2(IT_MAX), a counter of the number of
!    points used at each step of the Romberg iteration.
!
!    Local, integer ( kind = 4 ) SUB_NUM2(DIM_NUM), the number of subintervals
!    used in each direction, a refinement of the user's input SUB_NUM.
!
!    Local, real ( kind = 8 ) TABLE(IT_MAX), the difference table.
!
!    Local, real ( kind = 8 ) XX(DIM_NUM), an evaluation point.
!
  implicit none

  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) aaa(dim_num)
  real ( kind = 8 ) bbb(dim_num)
!  real ( kind = 8 ) en
  integer ( kind = 4 ) eval_num
  real ( kind = 8 ) factor
!  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) it
  integer ( kind = 4 ) iwork(dim_num)
  integer ( kind = 4 ) iwork2(it_max)
  integer ( kind = 4 ) kdim
  integer ( kind = 4 ) ll
  integer ( kind = 4 ) sub_num(dim_num)
  integer ( kind = 4 ) sub_num2(dim_num)
  real ( kind = 8 ) result
  real ( kind = 8 ) result_old
  real ( kind = 8 ) rnderr
!  real ( kind = 8 ) submid
  real ( kind = 8 ) sum1
  real ( kind = 8 ) weight
  real ( kind = 8 ) table(it_max)
  real ( kind = 8 ) tol
  real ( kind = 8 ) xx(dim_num)

  eval_num = 0

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ROMBERG_ND - Fatal error!'
    write ( *, '(a,i8)' ) '  DIM_NUM is less than 1.  DIM_NUM = ', dim_num
    stop
  end if

  if ( it_max < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ROMBERG_ND - Fatal error!'
    write ( *, '(a,i8)' ) '  IT_MAX is less than 1.  IT_MAX = ', it_max
    stop
  end if

  do i = 1, dim_num
    if ( sub_num(i) <= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ROMBERG_ND - Fatal error!'
      write ( *, '(a)' ) '  SUB_NUM(I) is less than 1.'
      write ( *, '(a,i8)' ) '  for I = ', i
      write ( *, '(a,i8)' ) '  SUB_NUM(I) = ', sub_num(i)
      stop
    end if
  end do

  ind = 0
  rnderr = epsilon ( 1.0D+00 )
  iwork2(1) = 1
  sub_num2(1:dim_num) = sub_num(1:dim_num)

  if ( 1 < it_max ) then
    iwork2(2) = 2
  end if

  it = 1

  !if (PID==0) then
  !  open(unit=94,file='results/rawerho.txt')
  !  open(unit=95,file='results/rawephi.txt')
  !  open(unit=96,file='results/rawez.txt')
  !  open(unit=97,file='results/funcx.txt')
  !  open(unit=98,file='results/funcy.txt')
  !  open(unit=99,file='results/funcval.txt')
  !endif
  do! loop for augmenting 'it'
  !  if (PID==0) then
  !      close(94,status="delete")
  !      close(95,status="delete")
  !      close(96,status="delete")
  !      close(97,status="delete")
  !      close(98,status="delete")
  !      close(99,status="delete")
  !      open(unit=94,file='results/rawerho.txt')
  !      open(unit=95,file='results/rawephi.txt')
  !      open(unit=96,file='results/rawez.txt')
  !      open(unit=97,file='results/funcx.txt')
  !      open(unit=98,file='results/funcy.txt')
  !      open(unit=99,file='results/funcval.txt')
  !  endif

    sum1 = 0.0D+00
!
!  Generate every point XX in the product region, and evaluate F(XX).
!
    iwork(1:dim_num) = 1

    do! loop over points xx, summing

      xx(1:dim_num) = &
       ( real ( 2 * sub_num2(1:dim_num) - 2 * iwork(1:dim_num) + 1, kind = 8 ) &
       * aaa(1:dim_num)   &
       + real (                         + 2 * iwork(1:dim_num) - 1, kind = 8 ) &
       * bbb(1:dim_num) ) &
       / real ( 2 * sub_num2(1:dim_num),                            kind = 8 )

      sum1 = sum1 + func ( xx )
      eval_num = eval_num + 1

      kdim = dim_num

      do while ( 0 < kdim )

        if ( iwork(kdim) < sub_num2(kdim) ) then
          iwork(kdim) = iwork(kdim) + 1
          exit
        end if

        iwork(kdim) = 1

        kdim = kdim - 1

      end do

      if ( kdim == 0 ) then
        exit
      end if

    end do! end loop over points xx, summing

    weight = product ( ( bbb(1:dim_num) - aaa(1:dim_num) ) &
    / ( real(sub_num2(1:dim_num),kind=8) )         )

!
!  Done with summing.
!
    table(it) = weight * sum1

    if ( it <= 1 ) then

      result = table(1)
      result_old = result

      if ( it_max <= it ) then
        ind = 1
        exit
      end if

      it = it + 1

      sub_num2(1:dim_num) = iwork2(it) * sub_num2(1:dim_num)

      cycle

    end if
!
!  Compute the difference table for Richardson extrapolation.
!
    do ll = 2, it
      i = it + 1 - ll
      factor = real ( iwork2(i)**2, kind = 8 ) &
             / real ( iwork2(it)**2 - iwork2(i)**2, kind = 8 )
      table(i) = table(i+1) + ( table(i+1) - table(i) ) * factor
    end do

    result = table(1)

!
!  Terminate successfully if the estimated error is acceptable.
!
    if ( abs ( result - result_old ) <= abs ( result * ( tol + rnderr ) ) ) then
      ind = 1
      exit
    end if
!
!  Terminate unsuccessfully if the iteration limit has been reached.
!
    if ( it_max <= it ) then
      ind = -1
      exit
    end if
!
!  Prepare for another step.
!
    result_old = result

    it = it + 1

    iwork2(it) = int ( 1.5D+00 * real ( iwork2(it-1), kind = 8 ) )

    sub_num2(1:dim_num) = &
      int ( 1.5D+00 * real ( sub_num2(1:dim_num), kind = 8 ) )

  end do ! loop for augmenting 'it'

  return
end

function func ( xx )
    ! calculates z-comp of complex poynting vec.
    ! See Jackson Eq 6.132, replace H with mu0 and B, replace mu0 with atomic unit rep (mu0=4pi/c2)

    implicit none
    real ( kind = 8 ) func, xx(2)

    call makeFields(xx(1),xx(2),z,t)
    func=real( fields(1)*conjg(fields(5))-fields(2)*conjg(fields(4)) )*c2/(8.d0*pi)

!    if (PID==0) then
!        write(94,*) Real(fields(1))
!        write(95,*) Real(fields(2))
!        write(96,*) Real(fields(3))
!        write(97,*) xx(1)
!        write(98,*) xx(2)
!        write(99,*) func
!    endif

    return
end function

end module normModule
