! 2D Solver for Heat Diffusion through Different Geological Materials
!==================================================================== 

! Discretization method: Finite Differences (FD)

! module for different data types
module DataTypes
    implicit none

    integer, parameter :: &

    ! REALS
    !------ 
    ! single precision
    R_SP = selected_real_kind(6, 36), &
    ! double precision
    R_DP = selected_real_kind(15, 307), &
    ! quadruple-precision
    R_QP = selected_real_kind(33, 4931)

    ! default; whats used
    integer, parameter :: RK = R_DP

    ! INTEGERS
    !--------- 
    I1 = selected_int_kind(2), &
    I2 = selected_int_kind(4), &
    I3 = selected_int_kind(9), &
    I4 = selected_int_kind(18)

    ! default; whats used
    integer, parameter :: IK = I3

    ! For formatting
    ! --------------
    ! REALS
    character(len=*), parameter :: &
    
    ! single precision
    R_SP_FMT = "f0.6", &
    ! double precision
    R_DP_FMT = "f0.15", &
    ! quadruple precision
    R_QP_FMT = "f0.33"

    ! default; whats used 
    character(len=*), parameter :: RK_FMT = R_DP_FMT

end module DataTypes

!------------------------------------------------------

! Software Configuration; useful for OOP purposes
! Note: data members &/or methods for a class are prefixed with m
module Class_Configuration
    ! need access to different data types
    use DataTypes
    implicit none 
    ! No access allowed! 
    private

    type, public :: Configuration

        ! Note: The rock (or rather, domain) is a square
        ! Thermal diffusivity coefficients: units: (m^2/s); each one to be run individually with differing numbers of threads; e.g. 2, 4, 6, 8
        real(RK) :: m_ThermDiffCoeff = 1.3e-4 ! sandstone
        ! real(RK) :: m_ThermDiffCoeff = 9e-5! basalt
        ! real(RK) :: m_ThermDiffCoeff = 2.6e-4 ! quartzite

        ! temperature of rock type; Boundary Conditions (BC); temperature at boundaries is uniform 
        real(RK) :: m_TempTopLeft = 100._RK, &
        real(RK) :: m_TempTopRight = 85._RK, &
        real(RK) :: m_TempBottomLeft = 40.0_RK, &
        real(RK) :: m_TempBottomRight = 35.0_RK

        ! domain properties and geometry
        real(RK) :: m_DomainLength = 50._RK, &
        integer(IK) :: m_xNumNodes = 300 ! number of nodes (or points) in the x-direction of the domain

    end type Configuration

    ! interface for a constructor; user-defined
    interface Configuration
        module procedure createConfiguration
    end interface Configuration

!
contains
    type(Configuration) function createConfiguration(configFilePath)

        character(len=*), intent(in) :: configFilePath
        integer :: configFileID

        ! # related to the mantissa of a single-/double-precision IEEE real
        ! => int(REAL-TYPE-NUMBER, IK) == MISS, True, providing:
        ! a) REAL-TYPE-NUMBER was initialised with MISS
        ! b) other instructions (e.g. NAMELIST-I/O here) did not modify the value of REAL-TYPE-NUMBER

        ! Basically, reals can't be true/false because they hold floating-point numbers

        integer(IK), parameter :: trueValue = - 9999 ! true value; used for initialisation 

        ! local variables to mirror the ones in the NAMELIST
        ! length of domain
        real :: DomainLength = trueValue
        ! number of nodes/points in the domain
        integer :: xNode = trueValue

        ! thermal diffusivity coefficient
        real :: ThermDiffCoeff = trueValue

        ! temperatures across the domain
        real :: TempTopLeft = trueValue, &
                TempTopRight = trueValue, &
                TempBottomLeft = trueValue, &
                TempBottomRight = trueValue

        ! Name list; so we can combine related variables together which are then read with a single statement
        NameList/heat_diffusion_parameters/ DomainLength, ThermDiffCoeff, xNode, &
                                            TempTopLeft, TempTopRight, TempBottomLeft, TempBottomRight

        ! I/O operations; open the file on the unit number specified by configFileID
        open(newunit = configFileID, file = trim(configFilePath), status = 'old', action = 'read')
        ! read the name list file 
        read(configFileID, nml = heat_diffusion_parameters)
        ! close name list file
        close(configFileID)

        ! Note: this isn't essential; good to let the user know about (if any) issues r.e. the name list file
        write(*,'(">> START: Name List Read <<")')
        write(*, nml = heat_diffusion_parameters)
        write(*,'(">> END: Name List Read <<")')

        ! assimilate data read from name list into new object's internal state.
        ! Note: the safe-guard-constant is used so that default values from the type-definition are overwritten only if 
        ! the user provides replacement values (in the name list)

        ! length of material (or domain)
        if(int(DomainLength, IK) /= trueValue) createConfiguration % m_DomainLength = DomainLength
        ! thermal diffusivity 
        if(int(ThermDiffCoeff, IK) /= trueValue) createConfiguration % m_ThermDiffCoeff = ThermDiffCoeff
        ! number nodes/points in domain
        if(xNode /= trueValue) createConfiguration % m_xNumNodes
        ! boundary temperatures of the domain
        if(int(TempTopLeft, IK) /= trueValue) createConfiguration % m_TempTopLeft = TempTopLeft
        if(int(TempTopRight, IK) /= trueValue) createConfiguration % m_TempTopRight = TempTopRight
        if(int(TempBottomLeft, IK) /= trueValue) createConfiguration % m_TempBottomLeft = TempBottomLeft
        if(int(TempBottomRight, IK) /= trueValue) createConfiguration % m_TempBottomRight = TempBottomRight

    end function createConfiguration
end module Class_Configuration

! solver for the heat diffusion equation
module Class_Solver
    ! bring in other required classes/modules
    use DataTypes
    use Class_Configuration
    implicit none
    private

    type, public :: Solver

        private
        type(Configuration) :: m_Configuration

        ! number of time iterations
        real(RK) :: m_tSteps 
        ! spacing between nodes of domain
        real(RK) :: m_dx
        ! time step/increment
        real(RK) :: m_dt
        ! these coefficients hold some algebra; Method of Least Squares
        real(RK) :: m_a, &
                    m_b
        ! horizontal (U) component and vertical (V) component of solution arrays
        real(RK), allocatable, dimension(:,:) :: m_U, m_V
        ! iterations: max. number of iterations and current iteration
        integer(IK) :: m_MaxNumIterations, m_CurrentIteration = 0

    contains
        ! hide methods
        private
        procedure, public :: initialisation 
        procedure, public :: runSolver
        procedure, public :: writeToASCII ! method for writing data to ASCII file; default
        ! procedure, public :: writeToDAT ! method for writing to a .dat file
        procedure, public :: getTemperature ! get method
        ! these are kept private; user doesn't need to know about these!
        procedure :: advanceU ! for horizontal component of temperature field advanced from 1 point to the next via time marching method
        procedure :: advanceV ! for vertical component of temperature field advanced from 1 point to the next via time marching method

        ! clean-up; memory management; fortran destructor; gfortran compiler doesn't destroy objects inside of arrays; might have to comment-out the following line of code
        ! final :: cleanup 
    end type Solver

contains
    ! initialisation sub-routine 
    subroutine initalisation(this, configFilePath, simulationTime)

        ! initialise objects
        class(Solver), intent(inout) :: this
        character(len=*), intent(in) :: configFilePath
        real(RK), intent(in) :: simulationTime
        integer(IK) :: xNode
        integer(IK) :: i, j ! iteration counters
        ! amplification factor; relates to the stability of the equations
        real(RK) :: lambda

        ! call constructor 
        this % m_Configuration = Configuration(configFilePath)
        ! point in the domain
        xNode = this % m_Configuration % m_xNumNodes
        ! geometry specific; relationship of characteristic time to characteristic length for the discretized domain
        this % m_tSteps = xNode**2
        ! evaluate derived parameter _mlambda (stability of solution)
        this % m_MaxNumIterations = int(simulationTime * this % m_tSteps)
        ! value for space increment/differential
        this % m_dx = 1. / xNode
        ! value for time increment/differential 
        this % m_dt = 1. / this % m_tSteps
        ! Note: see paper noted in README.md under references; link provided to this open access paper
        ! amplitude factor
        lambda = (2.*this % m_dt) / (this % m_dx**2)
        ! coefficients a & b from the Method of Least Squares
        this % m_a = (1. - lambda) / (1. + lambda)
        this % m_b = lambda / (2.*(1. + lambda))

        ! shape of domain

                !           TOP
                !      --------------
                !      |            | 
                ! LEFT |            |  RIGHT
                !      |            | 
                !      --------------
                !          BOTTOM
       
        ! allocate memory for solution arrays
        allocate(this % m_U(0:xNode, 0:xNode), this % m_V(0:xNode, 0:xNode))
        ! initialise temperature across the domain
        ! initial horizontal component of temperature field
        this % m_U = 1.
        ! boundary conditions:
        ! TOP
        this % m_U(:, xNode) = [(1./3.*(i/real(xNode, RK)) + 2./3., i = 0, xNode)]
        ! LEFT
        this % m_U(0, :) = [(1./3.*(j/real(xNode, RK)) + 1./3., j = 0, xNode)]
        ! BOTTOM
        this % m_U(:, 0) = [(-1./3.*(i/real(xNode, RK)) + 1./3., i = 0, xNode)]
        ! RIGHT
        this % m_U(:, 0) = [(j/real(xNode, RK), j = 0, xNode)]

        ! use the horizontal component setup of the heat field to initialise the vertical component array
        this % m_V = this % m_U

    end subroutine initalisation

    ! getter method for temperature
    real(RK) function getTemperature(this, i, j)

        class(Solver), intent(in) :: this
        integer(IK), intent(in) :: i, j
        ! get the temperature of the horizontal and vertical components at point i,j in the domain
        getTemperature = 0.5*(this % m_U(i, j) + this % m_V(i, j))

    end function getTemperature

    ! method for time marching
    subroutine runSolver(this)

        class(Solver), intent(inout) ::this
        integer(IK) :: K ! iterative counter for time marching

        ! time marching: computationally resource-heavy so OpenMP will be used; i.e. share work across available threads
        ! pragmas indicate when we use the OpenMP API for multithreading the time marching technique 

        do k = 1, this % m_MaxNumIterations
            if(mod(k - 1, (this % m_MaxNumIterations - 1)/ 10) == 0) then
                write(*, '(i5, a)') nint((k*100.0) / this % m_MaxNumIterations), "%"
            end if

            ! now using OpenMP: 2 threads to be used: 1 for each component of the temperature field

            !$omp parallel num_threads(2)
                !$omp sections
                !$omp section
                    ! horizontal component: 1st thread
                    call this % advanceU()

                !$omp section
                    ! vertical component: 2nd thread
                    call this % advanceV()

                !$omp end sections
            !$omp end parallel
            this % m_CurrentNumIterations = this % m_CurrentNumIterations + 1 ! tracking each time step

        end do
    end subroutine runSolver

    ! method finding horizontal component of temperature field; advance from 1 point to the next
    subroutine advanceU(this)

        class(Solver), intent(inout) :: this
        integer(IK) :: i, j ! iteration counters

        ! update m_U in the TOP-RIGHT direction
        do j = 1, this % m_Configuration % m_xNumNodes - 1
            do i = 1, this % m_Configuration % m_xNumNodes - 1 ! boundaries

                this % m_U(i, j) = this % m_a*this % m_U(i, j) + this % m_b*(&
                                   this % m_U(i - 1, j) + this % m_U(i + 1, j) + this % m_U(i, j - 1) + this % m_U(i, j + 1))

            end do
        end do
    end subroutine advanceU

    ! now perform the same for the vertical component
    ! method finding horizontal component of temperature field; advance from 1 point to the next
    subroutine advanceV(this)

        class(Solver), intent(inout) :: this 
        integer(IK) :: i, j ! iteration counters

        ! update m_V in the BOTTOM-LEFT direction
        do j = 1, this % m_Configuration % m_xNumNodes - 1
            do i = 1, this % m_Configuration % m_xNumNodes - 1 ! boundaries

                this % m_V(i, j) = this % m_a*this % m_V(i, j) + this % m_b*(&
                                   this % m_U(i - 1, j) + this % m_U(i + 1, j) + this % m_U(i, j - 1) + this % m_U(i, j + 1))

            end do
        end do
    end subroutine advanceV

    ! Output file
    ! writing data to ASCII file
    subroutine writeToASCII(this, outputFilePath)

        class(Solver), intent(in) :: this
        character(len=*), intent(in) :: outputFilePath

        ! coordinates and file I.D.
        integer(IK) :: x, y, outputFileID
        ! open file
        open(newunit = outputFileID, file = trim(outputFilePath), status = 'replace', action = 'write')

        ! write data to file
        write(outputFileID, '(a)') "# Output Data for the 2D Heat Diffusion Equation for Various Geological Materials"
        write(outputFileID, '(a, 2x, a)') '"s"', "# Time Unit"
        write(outputFileID, '(f0.8, 2x, a)') &
             (this % m_CurrentIteration * this % m_Configuration % m_DomainLength ** 2) / &
             (this % m_Configuration % m_ThermDiffCoeff * this % m_tSteps), &
             "# Current Time"

        write(outputFileID, '(a, 2x, a)') '"m"', "# X unit"
        write(outputFileID, '(i0, 2x, a)') this % m_Configuration % m_xNumNodes, "# Nx"
        write(outputFileID, '(a,2x,a)') '"m"', "# Y unit"
        write(outputFileID, '(i0, 2x, a)') this % m_Configuration % m_xNumNodes, "# Ny"
        write(outputFileID, '(a,2x,a)') '"Degree C"', "# Temperature Unit"

        ! x-axis
        do x = 0, this % m_Configuration % m_xNumNodes

            write(outputFileID, '(f0.8, 2x)', advance = 'no') this % m_dx * this % m_Configuration % m_DomainLength * x 

        end do

        write(outputFileID, '(a)') "X Values"

        ! y-axis
        do x = 0, this % m_Configuration % m_xNumNodes

            write(outputFileID, '(f0.8, 2x)', advance = 'no') this % m_dx * this % m_Configuration % m_DomainLength * y 

        end do

        write(outputFileID, '(a)') "Y Values"

        ! simulation results
        write(outputFileID, '(a)') "From Next Line to End: Simulated Heat Diffusion Data"

        do y = 0, this % m_Configuration % m_xNumNodes
            do x = 0, this % m_Configuration % m_xNumNodes

                write(outputFileID, '(f0.8, 2x)', advance = 'no') &

                    this % m_Configuration % m_TempBottomRight + this % getTemperature(x, y) * (this % m_Configuration % m_TempBottomLeft - this % m_Configuration % m_TempBottomRight)
            end do

            write(outputFileID, *) ! new line for separating rows in visualisation script
        end do

        ! close the file
        close(outputFileID)

    end subroutine writeToASCII

    ! destructor method; need to destroy objects and free memory
    subroutine freeMemory(this)

        type(Solver), intent(inout) :: this 

        ! free memory; destroy arrays which hold the horizontal and vertical components of the heat field
        deallocate(this % m_U, this % m_V)

    end subroutine freeMemory

end module Class_Solver

! main: might be worth separating this from the modules/classes; these can be in their own (module) files
program 2DHeatDiffusion

    ! bring in modules
    use DataTypes
    use Class_Solver
    implicit none

    type(Solver) :: square ! domain shape we apply the solver to
    real(RK) :: simulationTime = 0.1 ! number of time intervals to simulate

    ! bring in the name list file to configure the simulation and write resulting data to a .dat file
    character(len = 200) :: ConfigurationFile = "heat_diff_namelist.nml", & ! the name list file containing temperatures, thermal diffusivity values etc
                            outputFile = "Thermal_Diffusion_Field.dat"

    ! call initialiser
    call square % initalisation(ConfigurationFile, simulationTime)
    ! call and run the solver on the square configuration
    call square % runSolver()

    ! call the ASCII method and write the data to the .dat file
    call square % writeToASCII(outputFile)

    ! program finished
end program 2DHeatDiffusion
