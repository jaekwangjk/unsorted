program read_data
    implicit none
    
    integer :: n_atoms, i, atoms, id, ierror
    real :: lattice(3, 3), atom_positions(3)
    character(100) :: line, temp
    
    ! Open the file
    open(unit=10, file='pre_minimization.dump', status='old')

    ! Skip "ITEM: TIMESTEP"
    read(10, *)
	read(10, *)
	! Skip line "ITEM: NUMBER OF ATOMS"
	read(10, *)
    ! Read the number of atoms
    read(10, *) n_atoms
	write(*,*) n_atoms
    
    ! Skip the line "ITEM: BOX BOUNDS pp pp pp"
    read(10, *) 
	read(10, *) ! box x vector
	read(10, *) ! box y vector
	read(10, *) ! box z vector
	
    ! Skip the line "ITEM: ATOMS id type xs ys zs"
    read(10, *) 
	
    ! Read and output atom positions
    print *, "Atom positions:"
    do i = 1, n_atoms
        read(10, *, iostat=ierror) atoms, id, atom_positions(:)
		if (ierror .ne. 0) then
			write(*,*) "Premature end of file"
		endif
		
        write(*,*) atom_positions
    end do
    
	write(*,*) "saved at 8:33 AM"
	
    ! Close the file
    close(unit=10)
    
end program read_data