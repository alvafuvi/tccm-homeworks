program MD_code
    implicit none    
    integer :: num_atoms, i, j, k
    double precision :: epsilon, sigma, total_v, total_T, res, dt
    character(len=100) :: input_file  ! Name of the input file
    double precision, allocatable :: coord(:,:), mass(:), distance(:,:), acceleration(:,:), velocity(:,:)

    ! Constants
    epsilon = 0.0661d0
    sigma = 0.3345d0
    dt = 0.2d0
    write(*,*) "Put the name of the input file: "
    read (*,*) input_file

    ! Read the number of atoms
    num_atoms = read_Natoms(input_file)
    print *, 'Number of atoms: ', num_atoms

    ! Dynamically allocate memory for the matrices
    allocate(coord(num_atoms, 3))
    allocate(mass(num_atoms))
    allocate(distance(num_atoms, num_atoms))
    allocate(acceleration(num_atoms, 3))
    allocate(velocity(num_atoms, 3))

    ! Initialize velocity and acceleration
    velocity = 0.0d0
    acceleration = 0.0d0

    ! Read and display the molecule data
    call read_molecule(input_file, num_atoms, coord, mass)

    ! Compute and display interatomic distances
    call compute_distances(num_atoms, coord, distance)

    ! Compute the total Lennard-Jones potential
    total_v = V(epsilon, sigma, num_atoms, distance)
    print *, 'Total potential V: ', total_v

    ! Compute accelerations
    call compute_acc(num_atoms, coord, mass, distance, acceleration)

    ! Print initial accelerations
    print *, 'Initial accelerations:'
    do i = 1, num_atoms
        print *, acceleration(i, 1), acceleration(i, 2), acceleration(i, 3)
    end do

    ! Perform Verlet integration
    call verlet_algorithm(num_atoms, coord, velocity, acceleration, mass, distance, epsilon, sigma, dt)

    

contains

    ! -------------------------------------------------------------------
    ! Function to read the number of atoms
    ! -------------------------------------------------------------------
    integer function read_Natoms(input_file) result(result_num_atoms)
        implicit none
        character(len=*), intent(in) :: input_file
        integer :: unit_number

        unit_number = 10  ! Logical unit for the file

        ! Open the file
        open(unit=unit_number, file=input_file, status='old', action='read')

        ! Read the number of atoms from the first line
        read(unit_number, *) result_num_atoms

        ! Close the file
        close(unit_number)
    end function read_Natoms      

    ! -------------------------------------------------------------------
    ! Subroutine to read coordinates and masses
    ! -------------------------------------------------------------------
    subroutine read_molecule(input_file, num_atoms, coord, mass)
        implicit none
        character(len=*), intent(in) :: input_file
        integer, intent(in) :: num_atoms
        double precision, intent(out) :: coord(num_atoms, 3)
        double precision, intent(out) :: mass(num_atoms)
        integer :: i, unit_number

        unit_number = 10  ! Logical unit for the file

        ! Open the file
        open(unit=unit_number, file=input_file, status='old', action='read')

        ! Skip the first line (number of atoms)
        read(unit_number, *)

        ! Read coordinates and masses
        do i = 1, num_atoms
            read(unit_number, *) coord(i, 1), coord(i, 2), coord(i, 3), mass(i)
        end do

        ! Close the file
        close(unit_number)
    end subroutine read_molecule

    ! -------------------------------------------------------------------
    ! Subroutine to compute distances between atoms
    ! -------------------------------------------------------------------
     ! Subroutine to compute distances between atoms
    subroutine compute_distances(num_atoms, coord, distance)
    implicit none
    integer, intent(in) :: num_atoms
    double precision, intent(in) :: coord(num_atoms, 3)
    double precision, intent(out) :: distance(num_atoms, num_atoms)  ! Correct intent
    integer :: i, j

    ! Compute interatomic distances
        do i = 1, num_atoms
            do j = i + 1, num_atoms
            distance(i, j) = sqrt((coord(i, 1) - coord(j, 1))**2 + &
                                  (coord(i, 2) - coord(j, 2))**2 + &
                                  (coord(i, 3) - coord(j, 3))**2)
                        distance(j, i) = distance(i, j)  ! Symmetry
            end do
        end do

        
    end subroutine compute_distances

    ! -------------------------------------------------------------------
    ! Function to compute Lennard-Jones potential
    ! -------------------------------------------------------------------
    double precision function V(epsilon, sigma, num_atoms, distance) result(total_potential)
        implicit none
        double precision, intent(in) :: epsilon, sigma
        integer, intent(in) :: num_atoms
        double precision, intent(in) :: distance(num_atoms, num_atoms)
        integer :: i, j

        total_potential = 0.0d0

        ! Compute Lennard-Jones potential
        do i = 1, num_atoms
            do j = i + 1, num_atoms
                if (distance(i, j) > 0.0d0) then
                    total_potential = total_potential + 4.0d0 * epsilon * &
                        ((sigma / distance(i, j))**12 - (sigma / distance(i, j))**6)
                end if
            end do
        end do
    end function V

    double precision function T(num_atoms, velocity, mass) result(total_T)
        implicit none
        integer, intent(in) :: num_atoms
        double precision, intent(in) :: velocity(num_atoms, 3), mass(num_atoms)
        integer :: i

        total_T = 0.0d0
        do i = 1, num_atoms
            total_T = total_T + 0.5d0 * mass(i) * (velocity(i, 1)**2 + velocity(i, 2)**2 + velocity(i, 3)**2)
        end do
    end function T

    double precision function total_energy(V, T) result(total)
        implicit none
        double precision, intent(in) :: V, T  ! Inputs
                

        total = V + T  ! Compute the total energy

    end function total_energy


    ! -------------------------------------------------------------------
    ! Subroutine to compute accelerations
    ! -------------------------------------------------------------------
    subroutine compute_acc(num_atoms, coord, mass, distance, acceleration)
        implicit none
        integer, intent(in) :: num_atoms
        double precision, intent(in) :: coord(num_atoms, 3)
        double precision, intent(in) :: mass(num_atoms)
        double precision, intent(in) :: distance(num_atoms, num_atoms)
        double precision, intent(out) :: acceleration(num_atoms, 3)
        double precision :: epsilon, sigma, r, force
        integer :: i, j, k

        epsilon = 0.0661d0
        sigma = 0.3345d0

        ! Initialize acceleration to zero
        acceleration = 0.0d0

        ! Compute forces and update accelerations
        do i = 1, num_atoms
            do j = 1, num_atoms
                if (i /= j) then
                    r = distance(i, j)
                    if (r > 0.0d0) then
                        force = (24.0d0 * epsilon / r) * ((2.0d0 * (sigma / r)**12) - (sigma / r)**6)
                        do k = 1, 3
                            acceleration(i, k) = acceleration(i, k) - (force / mass(i)) * &
                                (coord(i, k) - coord(j, k)) / r
                        end do
                    end if
                end if
            end do
        end do
    end subroutine compute_acc

    ! -------------------------------------------------------------------
    ! Subroutine for Verlet integration
    ! -------------------------------------------------------------------
    subroutine verlet_algorithm(num_atoms, coord, velocity, acceleration, mass, distance, epsilon, sigma, dt)
        implicit none
        integer, intent(in) :: num_atoms
        double precision, intent(inout) :: coord(num_atoms, 3), velocity(num_atoms, 3), acceleration(num_atoms, 3)
        double precision, intent(in) :: mass(num_atoms)
        double precision, intent(inout) :: distance(num_atoms, num_atoms)
        double precision, intent(in) :: epsilon, sigma, dt
        integer :: j, k
        double precision :: total_V, total_T, final_energy
        integer :: step
        character(len=100) :: traj_file
    
        traj_file = 'trajectories.xyz'
    
        ! Open the trajectory file
        open(unit=20, file=traj_file, status='replace', action='write')
        write(20, '(A,I5)')'Number of atoms: ', num_atoms
        do step = 1, 1000
            ! Update coordinates
            do j = 1, num_atoms
                do k = 1, 3
                    coord(j, k) = coord(j, k) + velocity(j, k) * dt + 0.5d0 * acceleration(j, k) * dt**2
                end do
            end do
    
            ! Half-step velocity update
            do j = 1, num_atoms
                do k = 1, 3
                    velocity(j, k) = velocity(j, k) + 0.5d0 * acceleration(j, k) * dt
                end do
            end do
    
            ! Recompute distances and accelerations
            call compute_distances(num_atoms, coord, distance)
            call compute_acc(num_atoms, coord, mass, distance, acceleration)
    
            ! Full-step velocity update
            do j = 1, num_atoms
                do k = 1, 3
                    velocity(j, k) = velocity(j, k) + 0.5d0 * acceleration(j, k) * dt
                end do
            end do
    
            ! Compute total energy
            total_V = V(epsilon, sigma, num_atoms, distance)
            total_T = T(num_atoms, velocity, mass)
            final_energy = total_energy(total_V, total_T)
    
            ! Print energy at each step
            print *, 'Step: ', step, 'V = ', total_V, 'T = ', total_T, 'total energy = ', final_energy
    
            ! Write coordinates to the trajectory file in XYZ format
            
            write(20, '(A, E10.4, A, E10.4, A, E10.4)') '# T = ', total_T, ' V = ', total_V, ' T + V = ', final_energy
            write(20, '(A)') 'Step: ' // trim(adjustl(itoa(step)))
            do j = 1, num_atoms
                write(20, '(A, 3F10.4)') 'Ar', coord(j, 1), coord(j, 2), coord(j, 3)
            end do
            write(20, '(A)') ' '
        end do
    
        ! Close the trajectory file
        close(20)
    end subroutine verlet_algorithm
    

    ! -------------------------------------------------------------------
    ! Function to convert integer to string (itoa)
    ! -------------------------------------------------------------------
    function itoa(value) result(str)
        implicit none
        integer, intent(in) :: value
        character(len=20) :: str
        write(str, '(I0)') value
    end function itoa
  

end program MD_code
