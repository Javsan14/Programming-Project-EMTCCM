program geom_opt
    
    ! Main part of the code
    implicit none
    integer*4 :: i, j, natoms, nbonds, ncatoms, ncbonds, nangles, ndihedrals, npairs, nq
    integer*4 :: bending_terms(300,3), dihedral_terms(500,4), step, cartesian_count, ext_pos
    character(len=1), allocatable :: symbol(:)
    character(len=100) :: input_file, output_file, out_xyz_file
    logical :: file_exists
    integer*4, allocatable :: list_bonds(:,:), atom_pairs(:,:)
    real*8, parameter :: threshold = 0.001d0, cartesian_thresh = 0.00001d0, pi = acos(-1.0d0)
    real*8 :: total_stretch_E, total_bend_E, total_torsion_E, total_LJ_E, total_E, new_total_E, rms_grad, s_rms
    real*8, allocatable :: coord_mat(:,:), coord(:), new_coord(:), bond_lengths(:), angle_values(:), dihedral_values(:)
    real*8, allocatable :: int_coord(:), new_int_coord(:), update_int_coord(:), diff_coord(:)
    real*8, allocatable :: stretch_G(:), bend_G(:), torsion_G(:), LJ_G(:), total_G(:), total_G_int(:), new_total_G_int(:)
    real*8, allocatable :: inv_hess(:,:), new_inv_hess(:,:)
    real*8, allocatable :: big_vec(:), p_vec(:), s_vec(:), v_vec(:), y_vec(:), outer_big_s(:,:), outer_v_s(:,:), outer_s_v(:,:)
    real*8, allocatable :: B_mat(:,:), G_inv_mat(:,:)
    
    ! Read the value of the filename (.mol2 file) variable from the user
    write(*,*) 'Welcome, please enter the name of the mol2 file of the molecule to optimize:'
    write(*,*)
    read(*,*) input_file
    write(*,*)

    ! Find the position of the ".mol2" extension
    ext_pos = index(input_file, '.mol2')
    if (ext_pos > 0) then
        ! Extract the base filename without ".mol2"
        output_file = trim(input_file(1:ext_pos-1)) // '_opt.out'
    else
        ! Default output if ".mol2" not found
        output_file = trim(input_file) // '_opt.out'
    endif

    ! Read the mol2 file
    call read_mol2_file(input_file, symbol, coord_mat, list_bonds, natoms, nbonds, ncatoms, ncbonds)
    
    ! Store the nq value for the Wilson B matrix
    nq = nbonds + 6*ncatoms + 9*ncbonds
    
    ! We call the three subroutines to obtain the internal coordinates
    call get_bond_lenghts_and_stretching(coord_mat, symbol, list_bonds, nbonds, bond_lengths, total_stretch_E)
    call get_bending_terms_and_energy(coord_mat, symbol, list_bonds, nbonds, &
    bending_terms, angle_values, nangles, total_bend_E)
    call get_dihedral_terms_and_energy(coord_mat, symbol, bending_terms, nangles, &
    dihedral_terms, dihedral_values, ndihedrals, total_torsion_E)
    allocate(int_coord(nq))
    
    ! Initialize the Hessian and the internal coordinates vector by looping over the internal coordinates
    allocate(inv_hess(nq, nq))
    inv_hess = 0.0d0
    ! First streching terms
    do i = 1, nbonds
        inv_hess(i,i) = 1.0d0/600.0d0
        int_coord(i) = bond_lengths(i)
    enddo
    ! Now the bending ones
    do i = nbonds+1,nbonds+nangles
        inv_hess(i,i) = 1.0d0/150.0d0
        int_coord(i) = angle_values(i-nbonds)
    enddo
    ! And finally the dihedral terms
    do i = nbonds+nangles+1, nbonds+nangles+ndihedrals
        inv_hess(i,i) = 1.0d0/80.0d0
        int_coord(i) = dihedral_values(i-nbonds-nangles)
    enddo

    ! Initialize the RMS gradient (convergence criteria) and the step number
    rms_grad = 1.0d0
    step = 0
    
    ! Allocate also the necessary matrices
    allocate(coord(natoms*3))
    allocate(new_coord(natoms*3))
    allocate(stretch_G(natoms*3))
    allocate(bend_G(natoms*3))
    allocate(torsion_G(natoms*3))
    allocate(LJ_G(natoms*3))
    allocate(total_G(natoms*3))
    allocate(new_inv_hess(nq, nq))
    allocate(big_vec(nq))
    allocate(s_vec(nq))
    allocate(v_vec(nq))
    allocate(y_vec(nq))
    allocate(outer_big_s(nq, nq))
    allocate(outer_v_s(nq, nq))
    allocate(outer_s_v(nq, nq))
    allocate(B_mat(nq,3*natoms), G_inv_mat(nq,nq))
    allocate(total_G_int(nq))
    allocate(new_total_G_int(nq))
    allocate(new_int_coord(nq))
    allocate(update_int_coord(nq))
    allocate(diff_coord(3*natoms))

    write(*,*) '*drum roll*'
    write(*,*)

    ! Open the output file for writing
    open(unit=10, file=output_file, status='replace', action='write')
    
    ! Main loop of the optimization algorithm: BFGS in internal coordinates
    do while (rms_grad .ge. threshold)
        step = step + 1
        write(10, '(a,i3,a)') '*****  Optimization step: ', step, ' *****'

        ! 1. Calculate all internal coordinates, energies, and gradients.
        ! Internal coordinates and associated energies:
        ! - Bond lengths and stretching energy
        call get_bond_lenghts_and_stretching(coord_mat, symbol, list_bonds, nbonds, bond_lengths, total_stretch_E)
        ! - Bending angles and bending energy
        call get_bending_terms_and_energy(coord_mat, symbol, list_bonds, nbonds, &
        bending_terms, angle_values, nangles, total_bend_E)
        ! - Dihedral angles and torsional energy
        call get_dihedral_terms_and_energy(coord_mat, symbol, bending_terms, nangles, &
        dihedral_terms, dihedral_values, ndihedrals, total_torsion_E)
        ! - Atom pair interactions and Lennard-Jones energy
        call get_atom_pairs_and_energy(coord_mat, symbol, nbonds, list_bonds, nangles, bending_terms, natoms, &
        atom_pairs, npairs, total_LJ_E)

        ! Gradients with respect to Cartesian coordinates:
        call get_stretching_gradient(coord_mat, symbol, natoms, list_bonds, nbonds, bond_lengths, stretch_G)
        call get_bending_gradient(coord_mat, symbol, natoms, bending_terms, angle_values, nangles, bend_G)
        call get_dihedral_gradient(coord_mat, dihedral_terms, dihedral_values, ndihedrals, torsion_G)
        call get_pairs_gradient(coord_mat, symbol, atom_pairs, npairs, LJ_G)
        ! Total energy and gradient in Cartesian coordinates
        total_E = total_stretch_E + total_bend_E + total_torsion_E + total_LJ_E
        total_G = stretch_G + bend_G + torsion_G + LJ_G
        ! 2. Compute Wilson's B matrix and its inverse G matrix for the transformation to internal coordinates
        call get_wilson_b_mat(coord_mat, natoms, nq, nbonds, list_bonds, bond_lengths, &
        nangles, bending_terms, ndihedrals, dihedral_terms, B_mat, G_inv_mat)
        
        ! Transform gradient to internal coordinates
        total_G_int = matmul(matmul(G_inv_mat, B_mat), total_G)

        ! 3. Compute search direction using the inverse Hessian
        p_vec = -1.0d0 * matmul(inv_hess,total_G_int)

        ! 4. Step size control: Limit the RMS of the step vector to 0.02
        s_vec = p_vec
        call get_rms(nq, s_vec, s_rms)
        if (s_rms .gt. 0.02d0) then
            s_vec = s_vec * 0.02d0/s_rms
        endif
        ! Update internal coordinates with the scaled step vector
        new_int_coord = int_coord + s_vec

        ! 5. Convert updated internal coordinates back to Cartesian coordinates iteratively
        coord = reshape(transpose(coord_mat), [3*natoms])  ! Flatten Cartesian coordinate matrix for processing
        ! Initialize the convergence criteria and step count
        diff_coord = 1.0d0
        cartesian_count = 0

        do while (any(abs(diff_coord) .gt. cartesian_thresh))
            diff_coord = matmul(matmul(transpose(B_mat), G_inv_mat), s_vec) ! Difference vector
            new_coord = coord + diff_coord ! Update Cartesian coordinates
            coord_mat = transpose(reshape(new_coord, [3,natoms])) ! Reshape back to matrix form
            ! Update internal coordinates based on the new Cartesian coordinates
            call get_bond_lenghts_and_stretching(coord_mat, symbol, list_bonds, nbonds, bond_lengths, total_stretch_E)
            call get_bending_terms_and_energy(coord_mat, symbol, list_bonds, nbonds, &
            bending_terms, angle_values, nangles, total_bend_E)
            call get_dihedral_terms_and_energy(coord_mat, symbol, bending_terms, nangles, &
            dihedral_terms, dihedral_values, ndihedrals, total_torsion_E)

            ! Assign values to the updated internal coordinate vector
            do i = 1, nbonds
                update_int_coord(i) = bond_lengths(i)
            enddo
            do i = nbonds+1,nbonds+nangles
                update_int_coord(i) = angle_values(i-nbonds)
            enddo
            do i = nbonds+nangles+1, nbonds+nangles+ndihedrals
                update_int_coord(i) = dihedral_values(i-nbonds-nangles)
            enddo
            ! Compute the difference in internal coordinates for the next step
            s_vec = new_int_coord - update_int_coord
            ! Ensure dihedral angles remain within bounds (-π to π)
            do i = nbonds+nangles+1, nbonds+nangles+ndihedrals
                if (s_vec(i) .gt. pi) then
                    s_vec(i) = s_vec(i) - 2*pi
                elseif (s_vec(i) .lt. -pi) then
                    s_vec(i) = s_vec(i) + 2*pi
                endif
            enddo
            coord = new_coord ! Update Cartesian coordinates for next step
            cartesian_count = cartesian_count + 1
            !write(10, '(a,i3)') "Cartesian iteration number: ", cartesian_count
            !write(10, '(a,3f12.6)') "dx: ", diff_coord(1), diff_coord(2), diff_coord(3)
            !write(10, '(a,3f12.6)') "New s vector: ", update_int_coord(1) - int_coord(1), update_int_coord(2) - &
            !                                         int_coord(2), update_int_coord(3) - int_coord(3)
            !write(10, '(a,3f12.6)') "New cartesians: ", new_coord(1), new_coord(2), new_coord(3)
            !write(10, '(a,f12.6)') "diff_coord max: ", maxval(abs(diff_coord))
            !write(10, *)
        enddo

        ! 6. Recalculate step vector based on the final optimal internal coordinates
        ! and ensure dihedrals remain within bounds
        s_vec = update_int_coord - int_coord
        do i = nbonds+nangles+1, nbonds+nangles+ndihedrals
            if (s_vec(i) .gt. pi) then
                s_vec(i) = s_vec(i) - 2*pi
            elseif (s_vec(i) .lt. -pi) then
                s_vec(i) = s_vec(i) + 2*pi
            endif
        enddo

        ! 7. Calculate updated energies and gradients for the new coordinates
        call get_stretching_gradient(coord_mat, symbol, natoms, list_bonds, nbonds, bond_lengths, stretch_G)
        call get_bending_gradient(coord_mat, symbol, natoms, bending_terms, angle_values, nangles, bend_G)
        call get_dihedral_gradient(coord_mat, dihedral_terms, dihedral_values, ndihedrals, torsion_G)
        call get_pairs_gradient(coord_mat, symbol, atom_pairs, npairs, LJ_G)
        new_total_E = total_stretch_E + total_bend_E + total_torsion_E + total_LJ_E
        total_G = stretch_G + bend_G + torsion_G + LJ_G

        call get_wilson_b_mat(coord_mat, natoms, nq, nbonds, list_bonds, bond_lengths, &
        nangles, bending_terms, ndihedrals, dihedral_terms, B_mat, G_inv_mat)
        new_total_G_int = matmul(matmul(G_inv_mat, B_mat), total_G)

        ! 8. Compute the y vector (change in gradients in internal coordinates)
        y_vec = new_total_G_int - total_G_int

        ! 9. Update the inverse Hessian using the BFGS formula
        v_vec = matmul(inv_hess,y_vec)
        big_vec = (dot_product(s_vec, y_vec) + dot_product(y_vec, v_vec))*s_vec
        
        call get_outer_product(big_vec, s_vec, outer_big_s)
        call get_outer_product(v_vec, s_vec, outer_v_s)
        call get_outer_product(s_vec, v_vec, outer_s_v)

        new_inv_hess = inv_hess + (outer_big_s)/((dot_product(s_vec,y_vec))**2) - &
        (outer_v_s + outer_s_v)/(dot_product(s_vec,y_vec))

        ! 10. Update variables and compute convergence criteria
        coord = new_coord
        int_coord = update_int_coord
        inv_hess = new_inv_hess
        call get_rms(3*natoms, total_G, rms_grad) ! Recalculate RMS gradient for convergence check
        
        ! Print information about the step in the output file
        !write(10, '(a)') "p vector:"
        !write(10, '(10f12.6)') (p_vec(i), i = 1, size(p_vec))

        write(10, '(a)') "Search direction vector after cartesian iteration (s_k):"
        write(10, '(10f12.6)') (s_vec(i), i = 1, size(s_vec))

        write(10, '(a)') "New cartesian coordinates:"
        do i = 1, natoms
            write(10, '(3f12.6)') (new_coord(3*i - 2:3*i))
        end do

        write(10, '(a)') "New internal coordinates:"
        write(10, '(10f10.4)') (update_int_coord(i), i = 1, size(update_int_coord))

        write(10, '(a)') "New gradient (in internal coordinates):"
        write(10, '(10f12.6)') (new_total_G_int(i), i = 1, size(new_total_G_int))

        write(10, '(a,f12.6,a)') "Old energy: ", total_E, ' kcal/mol'
        write(10, '(a,f12.6,a)') "New energy: ", new_total_E, ' kcal/mol'
        write(10, '(a,f10.4)') "GRMS: ", rms_grad

        !write(10, '(a)') "Updated inverse Hessian:"
        !do i = 1, nq
        !    write(10, '(10f10.4)') (new_inv_hess(i, j), j = 1, nq)
        !end do
        write(10, *)  ! Blank line for readability
    enddo

    ! Close the output file
    close(10)

    ! Name the xyz output file
    if (ext_pos > 0) then
        ! Extract the base filename without ".mol2"
        out_xyz_file = trim(input_file(1:ext_pos-1)) // '_opt.xyz'
    else
        ! Default output if ".mol2" not found
        out_xyz_file = trim(input_file) // '_opt.xyz'
    endif

    ! Open the output xyz file
    open(unit=20, file=out_xyz_file, status="replace", action="write")

    ! Write in XYZ format
    write(20, '(i0)') natoms  ! First line: number of atoms
    write(20, '(a)') "Optimized structure by programming_project.f90"  ! Second line: comment
    do i = 1, natoms
        write(20, '(a2, 3f12.6)') symbol(i), new_coord(3*i-2:3*i)
    end do

    ! Close the xyz file
    close(20)

    write(*,*) 'Tadaaaaa'
    write(*,*)
    write(*, '(a,i3,a)') 'Optimization succesfully converged after ', step, ' steps!'
    write(*, '(a)') "Final cartesian coordinates:"
    do i = 1, natoms
        write(*, '(3f12.6)') (new_coord(3*i - 2:3*i))
    end do
    write(*, '(a,f12.6,a)') 'Energy at minimum: ', new_total_E, ' kcal/mol'
    write(*,*)
    write(*, '(A, A, A)') 'Info about the optimization saved in ', trim(adjustl(output_file)), ','
    write(*, '(A, A, A)') 'and final coordinates saved in ', trim(adjustl(out_xyz_file)), '.'

    ! The last part of the program is deallocating the arrays
    deallocate(symbol, coord, list_bonds, bond_lengths, angle_values, dihedral_values, atom_pairs)
    deallocate(stretch_G, bend_G, torsion_G, LJ_G, total_G, new_total_G_int)
    deallocate(int_coord, new_int_coord, update_int_coord, diff_coord)
    deallocate(inv_hess, new_inv_hess, new_coord)
    deallocate(big_vec, s_vec, v_vec, y_vec, outer_big_s, outer_v_s, outer_s_v)
    deallocate(B_mat, G_inv_mat)

contains

    subroutine read_mol2_file(mol2_file, symbol_array, coord_mat, bonds_list, natoms, nbonds, ncatoms, ncbonds)
        ! This subroutine reads a .mol2 file and extracts atomic symbols, Cartesian coordinates, 
        ! bond information, and counts of atoms and bonds.
        implicit none
        character(len=100), intent(in) :: mol2_file
        integer*4, intent(out) :: natoms, nbonds, ncatoms, ncbonds
        character(len=1), allocatable, intent(out) :: symbol_array(:)
        integer*4, allocatable, intent(out) :: bonds_list(:,:)
        real*8, allocatable, intent(out) :: coord_mat(:,:)

        ! Check if the file exists
        inquire(file=mol2_file, exist=file_exists)
        if (.not. file_exists) then
            write(*,*) 'Error: File does not exist: ', trim(mol2_file)
            stop
        endif

        ! Open the file
        open(2, file=mol2_file, status='old', action='read')

        ! Read data from the file and allocate the different arrays
        read(2,*) natoms, nbonds, ncatoms, ncbonds
        allocate(symbol_array(natoms))
        allocate(coord_mat(natoms,3))
        allocate(bonds_list(nbonds,2))
        do i=1,natoms
            read(2,*) (coord_mat(i,j), j=1,3), symbol_array(i)
        enddo

        do i=1,nbonds
            read(2,*) (bonds_list(i,j), j=1,2)
        enddo

        ! Close the file
        close(2)

        ! Remember to deallocate the corresponding arrays in the main program

    end subroutine read_mol2_file

    subroutine get_bond_lenghts_and_stretching(coord_mat, symbol_array, bonds_list, nbonds, bond_lengths, total_stretch_E)
        ! This subroutine calculates the bond lengths and bond stretching forces for a given molecule.
        ! The bond length is determined by the Euclidean distance between two atoms involved in a bond.
        ! The bond stretching force is computed based on a typical force field of MM.
        ! The inputs are the cartesian coordinates, the atomic symbols (for parameter assignment),
        ! the bonds list and number of bonds.
        ! The output is the total stretching energy
        implicit none
        character(len=1), intent(in) :: symbol_array(:)
        integer*4, intent(in) :: nbonds, bonds_list(:,:)
        real*8, intent(in) :: coord_mat(:,:)
        integer*4 :: atom1, atom2
        real*8, parameter :: cc_r = 1.5300d0, ch_r = 1.1100d0, cc_k = 300d0, ch_k = 350d0
        real*8, allocatable :: stretch_E(:)
        real*8, allocatable, intent(out) ::  bond_lengths(:)
        real*8, intent(out) :: total_stretch_E

        ! Allocate arrays
        allocate(bond_lengths(nbonds))
        allocate(stretch_E(nbonds))
    
        ! Initialize the total stretching energy
        total_stretch_E = 0.0d0

        ! Loop through bonds to compute bond lengths and the energy
        do i = 1, nbonds
            atom1 = bonds_list(i, 1)  ! First atom of the bond
            atom2 = bonds_list(i, 2)  ! Second atom of the bond

            ! Calculate Euclidean distance between atom1 and atom2
            bond_lengths(i) = norm2(coord_mat(atom2,:) - coord_mat(atom1,:))
            ! Calculate the stretching energy
            ! As there are only C-C or C-H bonds, assign the respective parameter
            ! depending on the symbol (if they share it, it's a C-C; if not, a C-H)
            if (symbol_array(atom1) == symbol_array(atom2)) then
                stretch_E(i) = cc_k*(bond_lengths(i) - cc_r)**2
            else 
                stretch_E(i) = ch_k*(bond_lengths(i) - ch_r)**2
            endif

            ! Accumulate the stretching energy for all bonds
            total_stretch_E = total_stretch_E + stretch_E(i)
        end do

        ! Deallocate arrays no longer useful
        deallocate(stretch_E)

    end subroutine get_bond_lenghts_and_stretching

    subroutine get_stretching_gradient(coord_mat, symbol_array, natoms, bonds_list, nbonds, bond_lengths, vec_stretch_G)
        ! This subroutine calculates the gradient of the bond stretching energy with respect to 
        ! the Cartesian coordinates of the atoms. The gradient is computed using the bond lengths 
        ! and atomic symbols for all bonds in the molecule. The resulting gradient is returned as a vector.
        implicit none
        character(len=1), intent(in) :: symbol_array(:)
        integer*4, intent(in) :: natoms, nbonds, bonds_list(:,:)
        real*8, intent(in) :: coord_mat(:,:), bond_lengths(:)
        integer*4 :: atom1, atom2
        real*8, parameter :: cc_r = 1.5300d0, ch_r = 1.1100d0, cc_k = 300d0, ch_k = 350d0
        real*8 :: gradient_mat(natoms,3), stretch_G(natoms,3)
        real*8, intent(out) ::  vec_stretch_G(3*natoms)
    
        ! Initialize the total stretching gradient
        stretch_G = 0.0d0

        ! Loop through bonds and calculate gradients
        do i = 1, nbonds
            atom1 = bonds_list(i, 1)  ! First atom of the bond
            atom2 = bonds_list(i, 2)  ! Second atom of the bond

            ! Obtain the gradient matrix for this bond
            gradient_mat = 0.0d0
            gradient_mat(atom1,:) = (coord_mat(atom1,:) - coord_mat(atom2,:))/bond_lengths(i)
            gradient_mat(atom2,:) = -1.d0 * gradient_mat(atom1,:)
            
            ! Calculate the stretching gradient (same parameter criteria as in energy)
            if (symbol_array(atom1) == symbol_array(atom2)) then
                stretch_G(:,:) = stretch_G(:,:) + 2*cc_k*(bond_lengths(i) - cc_r) * gradient_mat(:,:)
            else 
                stretch_G(:,:) = stretch_G(:,:) + 2*ch_k*(bond_lengths(i) - ch_r) * gradient_mat(:,:)
            endif

        end do

        ! Reshape the stretching gradient for the optimization (convert it to a vector)
        vec_stretch_G = reshape(transpose(stretch_G), [3*natoms])

    end subroutine get_stretching_gradient

    subroutine get_bending_terms_and_energy(coord_mat, symbol_array, list_bonds, nbonds, &
        bending_terms, bending_angles, nangles, total_bend_E)
        ! This subroutine identifies all bending angles in a molecular structure and calculates 
        ! the bending angles and their associated energies based on the atomic coordinates and 
        ! bond connectivity. It returns the total bending energy, individual angles, and the 
        ! bending terms (three-atom groupings forming angles).
        implicit none
        character(len=1), intent(in) :: symbol_array(:)
        integer*4, intent(in) :: nbonds, list_bonds(:,:)
        real*8, intent(in) :: coord_mat(:,:)
        integer*4 :: i, j, atom1, atom2, atom3
        real*8, allocatable :: bend_E(:)
        real*8, parameter :: pi = acos(-1.0d0), c_k = 60.0d0, h_k = 35.0d0, ang_eq = 109.5d0*pi/180.0d0
        real*8 :: norm_1, norm_2, vec1(3), vec2(3)
        integer*4, intent(out) :: bending_terms(300,3)
        integer*4, intent(out) :: nangles
        real*8, intent(out) :: total_bend_E
        real*8, allocatable, intent(out) :: bending_angles(:)
    
        ! Find the number of bending terms (angles)
        bending_terms = 0
        nangles = 0
        do i = 1, nbonds
            atom1 = list_bonds(i, 1)
            atom2 = list_bonds(i, 2)

            ! Loop through bonds again to find pairs that share one atom
            ! The shared atom is the central or 2, and the others are assigned as 1 and 3
            do j = i+1, nbonds
                if (list_bonds(j, 1) == atom1) then
                    nangles = nangles + 1
                    atom3 = list_bonds(j, 2)
                    bending_terms(nangles, 1) = atom2
                    bending_terms(nangles, 2) = atom1
                    bending_terms(nangles, 3) = atom3
                elseif (list_bonds(j, 2) == atom1) then
                    nangles = nangles + 1
                    atom3 = list_bonds(j, 1)
                    bending_terms(nangles, 1) = atom2
                    bending_terms(nangles, 2) = atom1
                    bending_terms(nangles, 3) = atom3
                elseif (list_bonds(j, 1) == atom2) then
                    nangles = nangles + 1
                    atom3 = list_bonds(j, 2)
                    bending_terms(nangles, 1) = atom1
                    bending_terms(nangles, 2) = atom2
                    bending_terms(nangles, 3) = atom3
                elseif (list_bonds(j, 2) == atom2) then
                    nangles = nangles + 1
                    atom3 = list_bonds(j, 1)
                    bending_terms(nangles, 1) = atom1
                    bending_terms(nangles, 2) = atom2
                    bending_terms(nangles, 3) = atom3
                end if
            end do
        end do

        ! With the bending terms, we can calculate the bending energy
        ! Initialize the total energy
        total_bend_E = 0.0d0
        allocate(bending_angles(nangles))
        allocate(bend_E(nangles))

        ! Loop through terms and calculate bending angles
        do i = 1, nangles
            atom1 = bending_terms(i, 1)  ! First atom of the angle
            atom2 = bending_terms(i, 2)  ! Second atom of the angle (central)
            atom3 = bending_terms(i, 3)  ! Third atom of the angle

            vec1(:) = coord_mat(atom1,:) - coord_mat(atom2,:)
            vec2(:) = coord_mat(atom3,:) - coord_mat(atom2,:)
            norm_1 = norm2(vec1)
            norm_2 = norm2(vec2)
            bending_angles(i) = (acos((dot_product(vec1,vec2))/(norm_1*norm_2)))
            
            ! Calculate the bending energy
            ! H-C-C and H-C-H parameters are equal, so if at least one atom is H that parameter is used
            if ((symbol_array(atom1) == 'H') .or. (symbol_array(atom2) == 'H') .or. (symbol_array(atom3) == 'H')) then
                bend_E(i) = h_k*(bending_angles(i) - ang_eq)**2
            else 
                bend_E(i) = c_k*(bending_angles(i) - ang_eq)**2
            endif
            ! Accumualte the total energy for every angle
            total_bend_E = total_bend_E + bend_E(i)
        end do

        ! Deallocate arrays no longer useful
        deallocate(bend_E)

    end subroutine get_bending_terms_and_energy

    subroutine get_bending_gradient(coord_mat, symbol_array, natoms, bending_terms, bending_angles, nangles, vec_bend_G)
        ! This subroutine calculates the bending gradient for each angle in the molecular structure. 
        ! The bending gradient is computed by considering the cross products of vectors formed by
        ! the atoms involved in the bending angle. The result is stored in a vector format for 
        ! optimization purposes.
        implicit none
        character(len=1), intent(in) :: symbol_array(:)
        integer*4, intent(in) :: nangles, natoms, bending_terms(:,:)
        real*8, intent(in) :: coord_mat(:,:), bending_angles(:)
        integer*4 :: i, atom1, atom2, atom3
        real*8, parameter :: pi = acos(-1.0d0), c_k = 60.0d0, h_k = 35.0d0, ang_eq = 109.5d0*pi/180.0d0
        real*8 :: norm_1, norm_2, norm_p, vec1(3), vec2(3), vec_p(3), cross_ba_p(3), cross_bc_p(3)
        real*8 :: gradient_mat(natoms,3), bend_G(natoms,3)
        real*8, intent(out) :: vec_bend_G(natoms*3)
    
        ! Initialize the total bending gradient
        bend_G = 0.0d0

        ! Loop through bending terms to calculate the gradient
        do i = 1, nangles
            atom1 = bending_terms(i, 1)  ! First atom of the angle
            atom2 = bending_terms(i, 2)  ! Second atom of the angle
            atom3 = bending_terms(i, 3)  ! Third atom of the angle

            ! Obtain the vectors that define an angle, and calculate their norm
            vec1(:) = coord_mat(atom1,:) - coord_mat(atom2,:)
            vec2(:) = coord_mat(atom3,:) - coord_mat(atom2,:)
            call get_cross_product(vec1, vec2, vec_p)
            norm_1 = norm2(vec1)
            norm_2 = norm2(vec2)
            norm_p = norm2(vec_p)
            call get_cross_product(vec1, vec_p, cross_ba_p)
            call get_cross_product(vec2, vec_p, cross_bc_p)
            ! Obtain the gradient matrix with these vectors
            gradient_mat = 0.0d0
            gradient_mat(atom1,:) = cross_ba_p(:)/(norm_1**2 * norm_p)
            gradient_mat(atom2,:) = -1.0d0*cross_ba_p(:)/(norm_1**2 * norm_p) + cross_bc_p(:)/(norm_2**2 * norm_p)
            gradient_mat(atom3,:) = -1.0d0*cross_bc_p(:)/(norm_2**2 * norm_p)
            
            ! Calculate the bending gradient (same parameter criteria as in energy)
            if ((symbol_array(atom1) == 'H') .or. (symbol_array(atom2) == 'H') .or. (symbol_array(atom3) == 'H')) then
                bend_G(:,:) = bend_G(:,:) + 2*h_k*(bending_angles(i) - ang_eq)*gradient_mat(:,:)
            else 
                bend_G(:,:) = bend_G(:,:) + 2*c_k*(bending_angles(i) - ang_eq)*gradient_mat(:,:)
            endif
        end do

        ! Reshape the bending gradient for the optimization (convert it to a vector)
        vec_bend_G = reshape(transpose(bend_G), [3*natoms])

    end subroutine get_bending_gradient

    subroutine get_dihedral_terms_and_energy(coord_mat, symbol_array, bending_terms, & 
        nangles, dihedral_terms, dihedral_angles, ndihedrals, total_torsion_E)
        ! This subroutine calculates the dihedral (torsional) terms and energy for a molecular system. 
        ! The torsional energy is computed based on the dihedral angles between planes formed by four atoms
        ! and the standard torsion potential formula is applied.
        implicit none
        character(len=1), intent(in) :: symbol_array(:)
        integer*4, intent(in) :: nangles, bending_terms(:,:)
        real*8, intent(in) :: coord_mat(:,:)
        integer*4, intent(out) :: dihedral_terms(500,4), ndihedrals
        real*8, intent(out) :: total_torsion_E
        real*8, allocatable, intent(out) :: dihedral_angles(:)
        integer*4 :: i, j, atom1, atom2, atom3, atom4, count_shared, shared_atoms(2), remaining_atoms(2)
        real*8, allocatable :: torsion_E(:)
        real*8, parameter :: pi = acos(-1.0d0), torsion_A = 0.3d0, torsion_n = 3.0d0
        real*8 :: vec1(3), vec2(3), vec3(3), normt, normu, norm_2, cross_t(3), cross_u(3), cross_v(3), cos_phi, sin_phi
    
        ! Initialize dihedral terms and counter for the number of dihedrals
        dihedral_terms = 0
        ndihedrals = 0
        ! Loop through all angles and check for shared atoms between angle pairs to form dihedrals
        do i = 1, nangles
            do j = i + 1, nangles
                ! Count the number of shared atoms between angle i and angle j
                count_shared = 0

                ! Check for the shared atoms, between atoms in diferent columns.
                ! If you compare all of them, you will get dihedrals like 3-2-1-4, that does not exist in ethane f.i.
                ! Also, we are only interested in C atoms, because only them can be central atoms in carbohydrates
                ! First check for shared atom 1
                if ((bending_terms(i, 1) == bending_terms(j, 2) .or. & 
                    bending_terms(i, 1) == bending_terms(j, 3)) .and. & 
                    symbol_array(bending_terms(i, 1)) == 'C') then
                    count_shared = count_shared + 1
                    shared_atoms(count_shared) = bending_terms(i, 1)
                endif

                ! Check for shared atom 2
                if ((bending_terms(i, 2) == bending_terms(j, 1) .or. & 
                    bending_terms(i, 2) == bending_terms(j, 3)) .and. & 
                    symbol_array(bending_terms(i, 2)) == 'C') then
                    count_shared = count_shared + 1
                    shared_atoms(count_shared) = bending_terms(i, 2)
                endif

                ! Check for shared atom 3
                if ((bending_terms(i, 3) == bending_terms(j, 1) .or. & 
                    bending_terms(i, 3) == bending_terms(j, 2)) .and. & 
                    symbol_array(bending_terms(i, 3)) == 'C') then
                    count_shared = count_shared + 1
                    shared_atoms(count_shared) = bending_terms(i, 3)
                endif

                ! If exactly 2 atoms are shared, we have a valid dihedral
                if (count_shared == 2) then
                    if (shared_atoms(1) /= shared_atoms(2)) then
                        ndihedrals = ndihedrals + 1

                        ! Identify the remaining atoms from both bending terms
                        ! We check that they are different to any central atom in order to be accepted
                        if (bending_terms(i, 1) /= shared_atoms(1) .and. bending_terms(i, 1) /= shared_atoms(2)) then
                            remaining_atoms(1) = bending_terms(i, 1)
                        elseif (bending_terms(i, 2) /= shared_atoms(1) .and. bending_terms(i, 2) /= shared_atoms(2)) then
                            remaining_atoms(1) = bending_terms(i, 2)
                        else
                            remaining_atoms(1) = bending_terms(i, 3)
                        endif

                        if (bending_terms(j, 1) /= shared_atoms(1) .and. bending_terms(j, 1) /= shared_atoms(2)) then
                            remaining_atoms(2) = bending_terms(j, 1)
                        elseif (bending_terms(j, 2) /= shared_atoms(1) .and. bending_terms(j, 2) /= shared_atoms(2)) then
                            remaining_atoms(2) = bending_terms(j, 2)
                        else
                            remaining_atoms(2) = bending_terms(j, 3)
                        endif

                        ! Form the dihedral term: (remaining_atoms(1) - shared_atoms - remaining_atoms(2))
                        dihedral_terms(ndihedrals, 1) = remaining_atoms(1)
                        dihedral_terms(ndihedrals, 2) = shared_atoms(2)
                        dihedral_terms(ndihedrals, 3) = shared_atoms(1)
                        dihedral_terms(ndihedrals, 4) = remaining_atoms(2)
                    endif
                endif
            enddo
        enddo

        ! With the torsional terms, we can calculate the torsional energy
        ! Initialize the total energy
        total_torsion_E = 0.0d0
        allocate(dihedral_angles(ndihedrals))
        allocate(torsion_E(ndihedrals))

        do i=1,ndihedrals
            atom1 = dihedral_terms(i, 1)
            atom2 = dihedral_terms(i, 2)
            atom3 = dihedral_terms(i, 3)
            atom4 = dihedral_terms(i, 4)
            ! Compute the torsional angle between the planes (AB-C) and (BC-D)
            ! First get the vectors of the dihedral
            vec1(:) = coord_mat(atom2,:) - coord_mat(atom1,:)
            vec2(:) = coord_mat(atom3,:) - coord_mat(atom2,:)
            vec3(:) = coord_mat(atom4,:) - coord_mat(atom3,:)

            ! Cross products of the planes
            call get_cross_product(vec1, vec2, cross_t)
            call get_cross_product(vec2, vec3, cross_u)
            call get_cross_product(cross_t, cross_u, cross_v)

            ! Norms of the necesary vectors
            normt = norm2(cross_t)
            normu = norm2(cross_u)
            norm_2 = sqrt(dot_product(vec2, vec2))

            cos_phi = dot_product(cross_t, cross_u) / (normt * normu)
            sin_phi = dot_product(vec2, cross_v) / (norm_2 * normt * normu)

            ! Calculate the torsional angle (dihedral angle) using the atan2 function
            dihedral_angles(i) = atan2(sin_phi, cos_phi)

            ! Compute the torsional energy: E_torsion = k(1 + cos(n * phi))
            ! No distinction regarding parameters
            torsion_E(i) = torsion_A * (1.0d0 + cos(torsion_n * dihedral_angles(i)))

            ! Accumulate the total torsional energy
            total_torsion_E = total_torsion_E + torsion_E(i)
        end do
    
        ! Deallocate arrays no longer useful
        deallocate(torsion_E)

    end subroutine get_dihedral_terms_and_energy

    subroutine get_dihedral_gradient(coord_mat, dihedral_terms, dihedral_angles, ndihedrals, vec_torsion_G)
        ! This subroutine computes the dihedral gradient, which is the derivative of the torsional energy with respect to cartesian coordinates.
        ! The torsional energy is calculated based on the dihedral angles of the system, and the gradient is used for optimization purposes.
        !
        ! For each dihedral term, the subroutine computes the necessary bond vectors (e.g., vec12, vec23, vec34) that define the geometry
        ! of the system. Cross products of these vectors are calculated to determine the directions of the torsional forces. 
        ! The gradient for each atom is then computed based on these vectors and the cross products.
        implicit none
        integer*4, intent(in) :: ndihedrals, dihedral_terms(:,:)
        real*8, intent(in) :: coord_mat(:,:), dihedral_angles(:)
        integer*4 :: i, atom1, atom2, atom3, atom4
        real*8, parameter :: pi = acos(-1.0d0), torsion_A = 0.3d0, torsion_n = 3.0d0
        real*8 :: vec12(3), vec23(3), vec34(3), vec13(3), vec24(3), vec_t(3), vec_u(3)
        real*8 :: tmp1(3), tmp21(3), tmp22(3), tmp31(3), tmp32(3), tmp4(3), cross_t_23(3), cross_u_23(3)
        real*8 :: normt, normu, norm23, gradient_mat(natoms,3), torsion_G(natoms,3)
        real*8, intent(out) :: vec_torsion_G(natoms*3)
    
        ! Initialize the total dihedral gradient
        torsion_G = 0.0d0

        ! Loop through each dihedral and compute the gradient
        do i=1,ndihedrals
            atom1 = dihedral_terms(i, 1)   ! First atom in the dihedral
            atom2 = dihedral_terms(i, 2)   ! Second (central) atom in the dihedral
            atom3 = dihedral_terms(i, 3)   ! Third (central) atom in the dihedral
            atom4 = dihedral_terms(i, 4)   ! Fourth atom in the dihedral
            ! Compute bond vectors between atoms (atom1-atom2), (atom2-atom3), and (atom3-atom4)
            vec12(:) = coord_mat(atom2,:) - coord_mat(atom1,:)
            vec23(:) = coord_mat(atom3,:) - coord_mat(atom2,:)
            vec34(:) = coord_mat(atom4,:) - coord_mat(atom3,:)
            vec13(:) = coord_mat(atom3,:) - coord_mat(atom1,:)
            vec24(:) = coord_mat(atom4,:) - coord_mat(atom2,:)

            ! Compute the cross products of the vectors to find the direction of the torsional force
            call get_cross_product(vec12, vec23, vec_t)
            call get_cross_product(vec23, vec34, vec_u)
            call get_cross_product(vec_t, vec23, cross_t_23)
            call get_cross_product(vec_u, vec23, cross_u_23)
            ! Compute the norms (magnitudes) of the cross products
            normt = norm2(vec_t)
            normu = norm2(vec_u)
            norm23 = norm2(vec23)

            ! To build the gradient matrix properly, let's decompose the process for every atom (row)
            gradient_mat = 0.0d0
            ! Atom 1
            call get_cross_product(cross_t_23, vec23, tmp1)
            gradient_mat(atom1,:) = tmp1/(normt**2 * norm23)
            ! Atom 2
            call get_cross_product(vec13, cross_t_23, tmp21)
            call get_cross_product(cross_u_23, vec34, tmp22)
            gradient_mat(atom2,:) = tmp21/(normt**2 * norm23) - tmp22/(normu**2 * norm23)
            ! Atom 2
            call get_cross_product(cross_t_23, vec12, tmp31)
            call get_cross_product(vec24, cross_u_23, tmp32)
            gradient_mat(atom3,:) = tmp31/(normt**2 * norm23) - tmp32/(normu**2 * norm23)
            ! Atom 4
            call get_cross_product(cross_u_23, vec23, tmp4)
            gradient_mat(atom4,:) = -1.0d0*tmp4/(normu**2 * norm23)

            ! Calculate the torsional energy gradient and update the total torsional gradient for the system
            torsion_G = torsion_G - torsion_n*torsion_A*sin(torsion_n*dihedral_angles(i)) * gradient_mat(:,:)
        end do

        ! Reshape the torsional gradient for the optimization (convert it to a vector)
        vec_torsion_G = reshape(transpose(torsion_G), [3*natoms])

    end subroutine get_dihedral_gradient

    subroutine get_atom_pairs_and_energy(coord_mat, symbol_array, nbonds, bond_list, &
        nangles, bending_terms, natoms, atom_pairs, npairs, total_LJ_E)
        ! This subroutine calculates the nonbonded atom pairs in the system and computes their Lennard-Jones (LJ) interaction energies.
        ! The goal is to identify pairs of atoms that are not directly bonded or part of an angle, then compute their interaction energies based on the 
        ! Lennard-Jones potential.
        !
        ! The Lennard-Jones potential for a pair of atoms i and j is given by the formula:
        !   V_LJ(r) = A_ij / r^12 - B_ij / r^6,
        ! where A_ij and B_ij are constants specific to the atom types, and r is the distance between the atoms.
        implicit none
        character(len=1), intent(in) :: symbol_array(:)
        integer*4, intent(in) :: nbonds, nangles, natoms, bond_list(:,:), bending_terms(:,:)
        real*8, intent(in) :: coord_mat(:,:)
        integer*4 :: i, j, k, is_bonded, is_angle, atom1, atom2, nvdws
        real*8, parameter :: A_HH = 4382.44d0, B_HH = 22.932d0
        real*8, parameter :: A_HC = 64393.99d0, B_HC = 108.644d0
        real*8, parameter :: A_CC = 946181.74d0, B_CC = 514.714d0
        real*8 :: rij(3), dist, Aij, Bij, r2_inv, r6_inv, r12_inv
        real*8, allocatable :: LJ_E(:)
        integer*4, intent(out) :: npairs
        integer*4, allocatable, intent(out) :: atom_pairs(:,:)
        real*8, intent(out) :: total_LJ_E

        ! Calculate the total number of possible atom pairs for the nonbonded interactions
        nvdws = natoms*(natoms-1)/2
        allocate(atom_pairs(nvdws,2))
    
        ! Initialize the number of nonbonded pairs
        npairs = 0
        atom_pairs = 0
    
        ! Loop over all pairs of atoms (i, j) with i < j
        do i = 1, natoms
            do j = i + 1, natoms
                ! Assume the atoms are not bonded or part of an angle by default
                is_bonded = 0
                is_angle = 0
    
                ! Check if the atoms form a bond (bond list check)
                do k = 1, nbonds
                    if ((bond_list(k, 1) == i .and. bond_list(k, 2) == j) .or. &
                        (bond_list(k, 1) == j .and. bond_list(k, 2) == i)) then
                        is_bonded = 1
                        exit ! Exit the loop if a bond is found
                    endif
                enddo
    
                ! Check if the atoms are bonded to the same atom (bending term check)
                do k = 1, nangles
                    if ((bending_terms(k, 1) == i .and. bending_terms(k, 3) == j) .or. &
                        (bending_terms(k, 1) == j .and. bending_terms(k, 3) == i)) then
                        is_angle = 1
                        exit ! Exit the loop if an angle term is found
                    endif
                enddo
    
                ! If atoms i and j are neither bonded nor part of an angle, add them to the nonbonded pair list
                if (is_bonded == 0 .and. is_angle == 0) then
                    npairs = npairs + 1
                    atom_pairs(npairs, 1) = i
                    atom_pairs(npairs, 2) = j
                endif
            enddo
        enddo

        ! Initialize the total Lennard-Jones energy
        total_LJ_E = 0.0d0
        allocate(LJ_E(npairs))

        ! Loop over all nonbonded pairs to compute the Lennard-Jones potential
        do i = 1, npairs
            atom1 = atom_pairs(i, 1)
            atom2 = atom_pairs(i, 2)

            ! Compute the distance between atoms atom1 and atom2
            rij(:) = coord_mat(atom1, :) - coord_mat(atom2, :)
            dist = norm2(rij)

            ! Select the correct Lennard-Jones parameters (Aij, Bij) based on atom types
            if (symbol_array(atom1) == 'H' .and. symbol_array(atom2) == 'H') then
                Aij = A_HH
                Bij = B_HH
            elseif (symbol_array(atom1) == 'C' .and. symbol_array(atom2) == 'C') then
                Aij = A_CC
                Bij = B_CC
            else
                Aij = A_HC
                Bij = B_HC
            endif

            ! Compute the Lennard-Jones energy for the pair (atom1, atom2)
            r2_inv = 1.0d0 / dist**2  ! Inverse of the squared distance
            r6_inv = r2_inv**3        ! Inverse of the distance raised to the 6th power
            r12_inv = r6_inv**2       ! Inverse of the distance raised to the 12th power

            LJ_E(i) = Aij * r12_inv - Bij * r6_inv

            ! Accumulate the total Lennard-Jones energy
            total_LJ_E = total_LJ_E + LJ_E(i)
        enddo
    
        ! Deallocate arrays no longer useful
        deallocate(LJ_E)

    end subroutine get_atom_pairs_and_energy

    subroutine get_pairs_gradient(coord_mat, symbol_array, atom_pairs, npairs, vec_LJ_G)
        ! This subroutine calculates the gradient of the Lennard-Jones (LJ) potential for nonbonded atom pairs in a molecular system.
        ! The gradient represents how the LJ potential changes with respect to cartesian coordinates and is essential for optimization and energy minimization.
        implicit none
        character(len=1), intent(in) :: symbol_array(:)
        integer*4, intent(in) :: atom_pairs(:,:), npairs
        real*8, intent(in) :: coord_mat(:,:)
        integer*4 :: i, atom1, atom2
        real*8, parameter :: A_HH = 4382.44d0, B_HH = 22.932d0
        real*8, parameter :: A_HC = 64393.99d0, B_HC = 108.644d0
        real*8, parameter :: A_CC = 946181.74d0, B_CC = 514.714d0
        real*8 :: rij(3), dist, Aij, Bij, r2_inv, r8_inv, r14_inv, gradient_mat(natoms,3), LJ_G(natoms,3)
        real*8, intent(out) :: vec_LJ_G(natoms*3)

        ! Initialize the LJ gradient
        LJ_G = 0.0d0
        ! Loop over all nonbonded pairs to compute the LJ gradient
        do i = 1, npairs
            atom1 = atom_pairs(i, 1)
            atom2 = atom_pairs(i, 2)

            ! Compute the distance between atoms atom1 and atom2
            rij(:) = coord_mat(atom1, :) - coord_mat(atom2, :)
            dist = norm2(rij)

            ! Here, there is no gradient matrix, but we can use the vector between both atoms as it
            gradient_mat = 0.0d0
            gradient_mat(atom1,:) = coord_mat(atom1,:) - coord_mat(atom2,:)
            gradient_mat(atom2,:) = -1.0d0*gradient_mat(atom1,:)

            ! Select the correct Lennard-Jones parameters (Aij, Bij) based on atom types
            if (symbol_array(atom1) == 'H' .and. symbol_array(atom2) == 'H') then
                Aij = A_HH
                Bij = B_HH
            elseif (symbol_array(atom1) == 'C' .and. symbol_array(atom2) == 'C') then
                Aij = A_CC
                Bij = B_CC
            else
                Aij = A_HC
                Bij = B_HC
            endif

            ! Compute the Lennard-Jones gradient for the pair (atom1, atom2)
            r2_inv = 1.0d0 / dist**2
            r8_inv = r2_inv**4           ! (1 / r^8)
            r14_inv = r8_inv*r2_inv**3   ! (1 / r^14)
            ! Calculate the gradient of the Lennard-Jones potential for the current pair
            LJ_G(:,:) = LJ_G(:,:) + (-12.0d0 * Aij * r14_inv + 6 * Bij * r8_inv)*gradient_mat(:,:)
        enddo

        ! Reshape the pairs gradient for the optimization (convert it to a vector)
        vec_LJ_G = reshape(transpose(LJ_G), [3*natoms])
    
    end subroutine get_pairs_gradient

    subroutine get_wilson_b_mat(coord_mat, natoms, nq, nbonds, bonds_list, bond_lengths, &
        nangles, bending_terms, ndihedrals, dihedral_terms, B_mat, G_inv_mat)
        ! This subroutine calculates the Wilson B matrix, which is used in the construction 
        ! of the generalized coordinate matrix for a molecular system. The B matrix relates the changes in 
        ! Cartesian coordinates to the internal coordinates (bond lengths, bond angles, and dihedrals) in the 
        ! system. This matrix is critical in normal mode analysis and force constant calculations. The subroutine 
        ! performs the following steps:
        !
        ! 1. **Bond Stretching Terms**: For each bond in the system, the B matrix is populated based on the 
        !    displacement between the bonded atoms, normalized by the bond length.
        !
        ! 2. **Bending Terms (Angles)**: For each bond angle, the subroutine computes the cross products of 
        !    vectors formed by the angle atoms, updating the B matrix with terms that describe how atomic 
        !    displacements affect the angle.
        !
        ! 3. **Dihedral Terms (Torsions)**: For each dihedral angle, the subroutine calculates cross products 
        !    involving the vectors between atoms forming the dihedral. These terms are then incorporated into 
        !    the B matrix, capturing how atomic displacements influence torsional motion.
        !
        ! 4. **Matrix Construction**: The B matrix is assembled from the above terms. The resulting matrix is 
        !    symmetric and reflects the relationships between the internal coordinates and the atomic Cartesian 
        !    coordinates.
        !
        ! 5. **Generalized Matrix (G) Calculation**: The matrix \( G \) is computed as the product of the B 
        !    matrix and its transpose. This matrix is used in subsequent calculations of the system's properties.
        !
        ! 6. **Eigenvalue Decomposition**: The matrix G is diagonalized using LAPACK’s "dsyev" routine to 
        !    obtain its eigenvalues and eigenvectors. These eigenvalues are used to compute the inverse of G (G-), 
        !    which is essential for calculating force constants and Hessians.
        !
        ! 7. **Numerical Stability**: A threshold is applied to the eigenvalues to ensure numerical stability 
        !    during the inversion process, avoiding divisions by very small values.
        !
        ! The subroutine outputs the Wilson B matrix and the inverse of the generalized matrix G-.
        implicit none
        external LAPACK
        integer*4, intent(in) :: natoms, nq, nbonds, bonds_list(:,:)
        integer*4, intent(in) :: nangles, bending_terms(:,:)
        integer*4, intent(in) :: ndihedrals, dihedral_terms(:,:)
        real*8, intent(in) :: coord_mat(:,:), bond_lengths(:)
        integer*4 :: i, atom1, atom2, atom3, atom4, lwork, info
        real*8 :: norm_1, norm_2, normp, vec1(3), vec2(3), vec_p(3), cross_ba_p(3), cross_bc_p(3)
        real*8 :: vec12(3), vec23(3), vec34(3), vec13(3), vec24(3), vec_t(3), vec_u(3)
        real*8 :: tmp1(3), tmp21(3), tmp22(3), tmp31(3), tmp32(3), tmp4(3), cross_t_23(3), cross_u_23(3)
        real*8 :: normt, normu, norm23
        real*8, allocatable :: eigenvalues(:), eigenvectors(:,:), work(:), G_mat(:,:), lambda(:,:)
        real*8, allocatable, intent(out) ::  B_mat(:,:), G_inv_mat(:,:)
    
        ! Allocate the B matrix to hold the Wilson B matrix (nq x 3*natoms)
        allocate(B_mat(nq,3*natoms))
        B_mat = 0.0d0
        ! First, do the first nbonds terms (streching internal terms)
        do i = 1, nbonds
            atom1 = bonds_list(i, 1)  ! First atom of the bond
            atom2 = bonds_list(i, 2)  ! Second atom of the bond
            B_mat(i,atom1*3-2:atom1*3) = (coord_mat(atom1,:) - coord_mat(atom2,:))/bond_lengths(i)
            B_mat(i,atom2*3-2:atom2*3) = -1.d0 * B_mat(i,atom1*3-2:atom1*3)
        enddo

        ! Now, the bending terms
        do i = 1, nangles
            atom1 = bending_terms(i, 1)  ! First atom of the angle
            atom2 = bending_terms(i, 2)  ! Second atom of the angle
            atom3 = bending_terms(i, 3)  ! Third atom of the angle

            vec1(:) = coord_mat(atom1,:) - coord_mat(atom2,:)
            vec2(:) = coord_mat(atom3,:) - coord_mat(atom2,:)
            norm_1 = norm2(vec1)
            norm_2 = norm2(vec2)
            call get_cross_product(vec1, vec2, vec_p)
            normp = norm2(vec_p)
            call get_cross_product(vec1, vec_p, cross_ba_p)
            call get_cross_product(vec2, vec_p, cross_bc_p)
            ! Here, the indices of the rows are the actual bending term plus all bonds
            B_mat(nbonds+i,atom1*3-2:atom1*3) = cross_ba_p(:)/(norm_1**2 * normp)
            B_mat(nbonds+i,atom2*3-2:atom2*3) = -1.0d0*cross_ba_p(:)/(norm_1**2 * normp) + &
            cross_bc_p(:)/(norm_2**2 * normp)
            B_mat(nbonds+i,atom3*3-2:atom3*3) = -1.0d0*cross_bc_p(:)/(norm_2**2 * normp)
        end do

        ! Finally, the dihedral
        do i=1,ndihedrals
            atom1 = dihedral_terms(i, 1)
            atom2 = dihedral_terms(i, 2)
            atom3 = dihedral_terms(i, 3)
            atom4 = dihedral_terms(i, 4)
            ! Obtain the vectors for the B matrix evaluation
            vec12(:) = coord_mat(atom2,:) - coord_mat(atom1,:)
            vec23(:) = coord_mat(atom3,:) - coord_mat(atom2,:)
            vec34(:) = coord_mat(atom4,:) - coord_mat(atom3,:)
            vec13(:) = coord_mat(atom3,:) - coord_mat(atom1,:)
            vec24(:) = coord_mat(atom4,:) - coord_mat(atom2,:)

            ! Cross products
            call get_cross_product(vec12, vec23, vec_t)
            call get_cross_product(vec23, vec34, vec_u)
            call get_cross_product(vec_t, vec23, cross_t_23)
            call get_cross_product(vec_u, vec23, cross_u_23)

            normt = norm2(vec_t)
            normu = norm2(vec_u)
            norm23 = norm2(vec23)

            ! Atom 1's contribution to the B matrix
            call get_cross_product(cross_t_23, vec23, tmp1)
            B_mat(nbonds+nangles+i,atom1*3-2:atom1*3) = tmp1/(normt**2 * norm23)
            ! Atom 2
            call get_cross_product(vec13, cross_t_23, tmp21)
            call get_cross_product(cross_u_23, vec34, tmp22)
            B_mat(nbonds+nangles+i,atom2*3-2:atom2*3) = tmp21/(normt**2 * norm23) - tmp22/(normu**2 * norm23)
            ! Atom 3
            call get_cross_product(cross_t_23, vec12, tmp31)
            call get_cross_product(vec24, cross_u_23, tmp32)
            B_mat(nbonds+nangles+i,atom3*3-2:atom3*3) = tmp31/(normt**2 * norm23) - tmp32/(normu**2 * norm23)
            ! Atom 4
            call get_cross_product(cross_u_23, vec23, tmp4)
            B_mat(nbonds+nangles+i,atom4*3-2:atom4*3) = -1.0d0*tmp4/(normu**2 * norm23)
        end do

        ! From the Wilson B matrix, we can compute the G and G- matrices
        ! Allocate arrays
        allocate(G_mat(nq,nq), G_inv_mat(nq,nq), eigenvalues(nq), eigenvectors(nq,nq))
        allocate(lambda(nq,nq))
        ! Set the workspace size for LAPACK's diagonalization routine
        lwork = max(1, 3*nq - 1)
        allocate(work(lwork))

        ! Calculate G
        G_mat = matmul(B_mat, transpose(B_mat))
        eigenvectors = G_mat

        ! Call DSYEV to diagonalize G
        call dsyev('V', 'U', nq, eigenvectors, nq, eigenvalues, work, lwork, info)
        ! Check if diagonalization was successful
        if (info /= 0) then
            print *, "DSYEV failed with error code: ", info
            stop
        endif

        ! Construct lambda matrix (inverse of eigenvalues) with a small cutoff for numerical stability
        lambda = 0.0d0
        do i = 1, nq
            if (eigenvalues(i) .gt. 1.0E-12) then
                lambda(i,i) = 1.0d0 / eigenvalues(i)
            endif
        enddo

        ! Finally, calculate the inverse general matrix
        G_inv_mat = matmul(matmul(eigenvectors, lambda), transpose(eigenvectors))

        ! Deallocate arrays no longer useful
        deallocate(G_mat, eigenvalues, eigenvectors, work, lambda)
        
    end subroutine get_wilson_b_mat

    subroutine get_rms(n, vector, rms)
        ! This subroutine computes the Root Mean Square (RMS) of a given vector.
        ! The RMS is calculated by summing the squares of all elements in the vector, 
        ! dividing the sum by the number of elements (n), and then taking the square root of that result.
        integer*4, intent(in) :: n
        real*8, intent(in) :: vector(n)
        integer*4 :: i
        real*8 :: sum
        real*8, intent(out) ::  rms
    
        ! Initialize the total sum
        sum = 0.0d0
        ! Loop over all elements in the vector, squaring them and adding to the sum
        do i = 1,n
            sum = sum + vector(i)**2
        enddo

        ! Calculate the RMS
        rms = sqrt(sum/n)
        
    end subroutine get_rms

    subroutine get_cross_product(vec1, vec2, cross_prod)
        ! This subroutine calculates the cross product of two vectors.
        real*8, intent(in) :: vec1(:), vec2(:)
        real*8, intent(out) :: cross_prod(3)
    
        cross_prod(1) = vec1(2) * vec2(3) - vec1(3) * vec2(2)
        cross_prod(2) = vec1(3) * vec2(1) - vec1(1) * vec2(3)
        cross_prod(3) = vec1(1) * vec2(2) - vec1(2) * vec2(1)
        
    end subroutine get_cross_product

    subroutine get_outer_product(vec1, vec2, outer)
        ! This subroutine computes the outer product of two vectors.
        ! The outer product of two vectors vec1 (of size m) and vec2 (of size n) results in an m x n matrix,
        ! where each element "outer(i, j)" is the product of "vec1(i)" and "vec2(j)". 
        integer*4 :: i, j, m, n
        real*8, intent(in) :: vec1(:), vec2(:)
        real*8, intent(out) ::  outer(:,:)

        ! Get the size of the input vectors
        m = size(vec1)
        n = size(vec2)
        ! Loop over all combinations of elements from vec1 and vec2 to fill the outer product matrix
        do i = 1, m
            do j = 1, n
                outer(i, j) = vec1(i) * vec2(j)
            end do
        end do    
        
    end subroutine get_outer_product

end program geom_opt