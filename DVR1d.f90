module DVR1d
  implicit none
  
  real(8), parameter :: pi = 4.0D0*atan(1.0D0)
  real(8) :: lefta, rightb   ! [a,b]

  abstract interface
        function PotentialFunction(x) result(V)
            real(8), intent(in) :: x
            real(8) :: V
        end function PotentialFunction
    end interface

  type :: DVR_Class
    integer :: ngrid        ! The number of grids
    real(8) :: weight       ! The DVR weight
    real(8) :: mass
    real(8), allocatable, dimension(:) :: grid          ! Grids
    real(8), allocatable, dimension(:,:) :: TransMat    ! The transMat for FBR->DVR
    real(8), allocatable, dimension(:,:) :: T           ! Kinetic Matrix
    real(8), allocatable, dimension(:,:) :: V           ! Potential Matrix
    real(8), allocatable, dimension(:,:) :: H           ! Hamiltonian Matrix
    real(8), allocatable, dimension(:,:) :: C           ! Coefficient Matrix
    real(8), allocatable, dimension(:) :: E             ! Eigenenergies
    procedure(PotentialFunction), pointer, nopass :: potential => null()    ! Potential Function

    contains
    procedure, public :: DVR_init
    procedure, public :: DVR_matrix_calc
    procedure, public :: DVR_diagH
  end type DVR_Class

contains

  subroutine DVR_init(this, GridsNumber, lefta, rightb, MassInAU, potential)
    class(DVR_Class) :: this
    integer, intent(in) :: GridsNumber
    real(8), intent(in) :: lefta, rightb, MassInAU
    procedure(PotentialFunction) :: potential
    integer :: i

    this%ngrid = GridsNumber
    this%mass = MassInAU
    this%potential => potential

    allocate(this%grid(this%ngrid))
    allocate(this%TransMat(this%ngrid, this%ngrid))
    allocate(this%T(this%ngrid, this%ngrid))
    allocate(this%V(this%ngrid, this%ngrid))
    allocate(this%H(this%ngrid, this%ngrid))
    allocate(this%C(this%ngrid, this%ngrid))
    allocate(this%E(this%ngrid))

    do i = 1, this%ngrid
        this%grid(i) = lefta + i * (rightb - lefta) / (this%ngrid + 1)
    end do
    this%weight = (rightb - lefta) / (this%ngrid + 1)
  end subroutine 

  subroutine DVR_matrix_calc(this)
    class(DVR_Class) :: this
    integer :: i, j

    associate(n => this%ngrid, T => this%T, V => this%V, H => this%H, B => this%TransMat, m => this%mass)
      do i = 1, n
        do j = 1, n
          B(i,j) = sqrt(2.0D0 / real(n + 1, kind=8)) * sin(i * j * pi / real(n + 1, kind=8))
          if (i == j) then
            T(i,j) = (pi**2) / (4.0D0 * m * (rightb - lefta)**2) & 
                     * ((2.0D0 * (real(n + 1, kind=8)**2) + 1.0D0) / 3.0D0 &
                     - 1.0D0 / max(sin(pi * i / real(n + 1, kind=8))**2, 1.0d-10))
            V(i,j) = this%potential(this%grid(i))
          else
            T(i,j) = ((-1.0)**(i - j)) * pi**2 / (4.0D0 * m * (rightb - lefta)**2) & 
                     * (1.0D0 / max(sin(pi * (i - j) / (2.0D0 * real(n + 1, kind=8)))**2, 1.0d-10) &
                     - 1.0D0 / max(sin(pi * (i + j) / (2.0D0 * real(n + 1, kind=8)))**2, 1.0d-10))
            V(i,j) = 0.0D0
          end if
        end do
      end do
      H = T + V
    end associate
  end subroutine

  subroutine DVR_diagH(this)
    class(DVR_Class) :: this
    real(8), dimension(3 * this%ngrid - 1) :: w
    real(8), allocatable, dimension(:,:) :: tempH
    integer :: info

    allocate(tempH(this%ngrid, this%ngrid))
    associate(n => this%ngrid, H => this%H, C => this%C, E => this%E)
      tempH = H
      call dsyev('V', 'L', n, tempH, n, E, w, 3 * n - 1, info)
      C = tempH
      deallocate(tempH)
      C = C/sqrt(this%weight)
    end associate
  end subroutine

end module DVR1d
