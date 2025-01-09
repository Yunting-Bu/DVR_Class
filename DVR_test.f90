program DVR_Test
  use DVR1d
  use, intrinsic :: iso_fortran_env, only : f8 => real64
  implicit none

  real(kind=f8), parameter :: au2ev = 27.2114_f8
  real(kind=f8), parameter :: dalton2au = 1822.888486209_f8
  real(kind=f8), parameter :: m = 0.50391_f8

    interface
        function harmonic(x) result(res)
            import :: f8
            real(kind=f8), intent(in) :: x
            real(kind=f8) :: res
        end function 

        function morse(x) result(res)
            import :: f8
            real(kind=f8), intent(in) :: x
            real(kind=f8) :: res
            real(kind=f8) :: De, re, a
        end function 

        function heiles_x(x) result(res)
            import :: f8
            real(kind=f8), intent(in) :: x
            real(kind=f8) :: res
        end function

        function heiles_y(y) result(res)
            import :: f8
            real(kind=f8), intent(in) :: y
            real(kind=f8) :: res
        end function
    end interface

  type(DVR_Class) :: dvr_harmonic
  type(DVR_Class) :: dvr_morse
  integer :: ngrid
  integer :: i, j

  !--for 2d Henon-Heiles--!
  type(DVR_Class) :: dvr_heiles_x
  type(DVR_Class) :: dvr_heiles_y

!------Harmonic------!
  ngrid = 50
  lefta = -10.0_f8
  rightb = 10.0_f8

  call dvr_harmonic%DVR_init(ngrid, lefta, rightb, 1.0_f8, harmonic)
  call dvr_harmonic%DVR_matrix_calc()
  call dvr_harmonic%DVR_diagH()
  
  write(*,'(/,a)') 'Harmonic Oscillator Eigenenergies in Hartree:'
  do i = 1, 6
    write(*,'(a,i4,2x,a,2x,f12.5)') 'n = ',i-1,'E_DVR = ',dvr_harmonic%E(i)
  end do 

  open(100,file='harmonic.log',status='replace')
  do i = 1, ngrid
      write(100,'(*(f15.6))') (dvr_harmonic%C(i,j),j=1,ngrid)
  end do 
  close(100)

!------Morse------!
  ngrid = 100
  lefta = -1.0_f8
  rightb = 6.0_f8
  call dvr_morse%DVR_init(ngrid, lefta, rightb, m*dalton2au, morse)
  call dvr_morse%DVR_matrix_calc()
  call dvr_morse%DVR_diagH()

  write(*,'(/,a)') 'Morse Potential Eigenenergies for H2 in eV:'
  do i = 1, 6
    write(*,'(a,i4,2x,a,2x,f12.5)') 'n = ',i-1,'E_DVR = ',dvr_morse%E(i)*au2ev
  end do 

  open(101,file='morse.log',status='replace')
  do i = 1, ngrid
      write(101,'(*(f15.6))') (dvr_morse%C(i,j),j=1,ngrid)
  end do 
  close(101)

!------Henon-Heiles------!

  call dvr_heiles_x%DVR_init(50, -5.0_f8, 5.0_f8, 1.0_f8, heiles_x)
  call dvr_heiles_x%DVR_matrix_calc()
  call dvr_heiles_x%DVR_diagH()

  call dvr_heiles_y%DVR_init(50, -5.0_f8, 5.0_f8, 1.0_f8, heiles_y)
  call dvr_heiles_y%DVR_matrix_calc()
  call dvr_heiles_y%DVR_diagH()

  write(*,'(/,a)') 'Henon-Heiles Potential Eigenenergies in Hartree:'
  do i = 1, 2
    do j = 1, 2
      write(*,'(a,i4,2x,a,i4,2x,a,2x,f12.5)') 'n1 = ',i-1,'n2 = ',j-1,'E_DVR = ',&
                                              dvr_heiles_x%E(i)+dvr_heiles_y%E(j)
    end do
  end do 

  open(102,file='heiles_x.log',status='replace')
  do i = 1, 50
      write(102,'(*(f15.6))') (dvr_heiles_x%C(i,j),j=1,50)
  end do 
  close(102)

  open(103,file='heiles_y.log',status='replace')
  do i = 1, 50
      write(103,'(*(f15.6))') (dvr_heiles_y%C(i,j),j=1,50)
  end do 
  close(103)

end program DVR_Test

function harmonic(x) result(res)
  use, intrinsic :: iso_fortran_env, only : f8 => real64
    real(kind=f8), intent(in) :: x
    real(kind=f8) :: res
    res = 0.5_f8 * x**2
end function 

function morse(x) result(res)
  use, intrinsic :: iso_fortran_env, only : f8 => real64
    real(kind=f8), intent(in) :: x
    real(kind=f8) :: De, re, a
    real(kind=f8) :: res
    De = 4.7446_f8 / 27.2114_f8
    re = 0.7416_f8 * 1.88973_f8
    a = 1.9426_f8 / 1.88973_f8
    res = De * (1.0_f8 - exp(-a * (x - re)))**2
end function

function heiles_x(x) result(res)
  use, intrinsic :: iso_fortran_env, only : f8 => real64
    real(kind=f8), intent(in) :: x
    real(kind=f8) :: res
    res = 0.5_f8 * x**2
end function

function heiles_y(y) result(res)
  use, intrinsic :: iso_fortran_env, only : f8 => real64
    real(kind=f8), intent(in) :: y
    real(kind=f8) :: res
    res = 0.5_f8 * y**2 + sqrt(0.0125_f8)*(-0.333_f8*y**3)
end function

