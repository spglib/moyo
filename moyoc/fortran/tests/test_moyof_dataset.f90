program test_moyof_dataset
    use moyo
    use, intrinsic :: iso_c_binding
    implicit none

    type(c_ptr) :: version_ptr
    type(c_ptr) :: dataset_ptr
    type(c_ptr) :: sgt_ptr
    type(MoyoDataset), pointer :: dataset
    type(MoyoSpaceGroupType), pointer :: sgt
    integer(c_int32_t), pointer :: orbits(:)
    integer(c_int32_t), pointer :: mapping_std_prim(:)
    integer(c_int32_t), pointer :: rot(:, :, :)
    integer(c_int32_t), parameter :: identity(3, 3) = reshape( &
        [1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])

    real(c_double) :: basis(3, 3)
    real(c_double) :: positions(3, 2)
    integer(c_int32_t) :: numbers(2)
    integer(c_int32_t) :: num_atoms
    real(c_double) :: symprec, angle_tolerance
    integer(c_int32_t) :: setting, hall_number
    logical(c_bool) :: rotate_basis

    real(c_double) :: a, c
    integer :: i, j

    ! Version
    print *, "moyo_version: ", moyo_to_string(moyo_version())
    if (len(moyo_to_string(moyo_version())) == 0) then
        error stop "moyo_version returned an empty string"
    end if

    ! hcp structure
    a = 3.17d0
    c = 5.14d0
    basis(:, 1) = [a, 0.0d0, 0.0d0]
    basis(:, 2) = [-a / 2.0d0, a * sqrt(3.0d0) / 2.0d0, 0.0d0]
    basis(:, 3) = [0.0d0, 0.0d0, c]
    positions(:, 1) = [1.0d0 / 3.0d0, 2.0d0 / 3.0d0, 1.0d0 / 4.0d0]
    positions(:, 2) = [2.0d0 / 3.0d0, 1.0d0 / 3.0d0, 3.0d0 / 4.0d0]
    numbers = [0_c_int32_t, 0_c_int32_t]
    num_atoms = 2

    symprec = 1d-4
    angle_tolerance = -1d0
    setting = MOYO_SETTING_SPGLIB
    hall_number = -1
    rotate_basis = .true._c_bool

    dataset_ptr = moyo_dataset_new(basis, positions, numbers, num_atoms, symprec, &
                                   angle_tolerance, setting, hall_number, rotate_basis)
    if (.not. c_associated(dataset_ptr)) then
        error stop "moyo_dataset_new returned NULL"
    end if

    call c_f_pointer(dataset_ptr, dataset)

    ! Identification
    print *, "dataset%number: ", dataset%number
    print *, "dataset%hall_number: ", dataset%hall_number
    print *, "dataset%hm_symbol: ", moyo_to_string(dataset%hm_symbol)
    if (dataset%number /= 194) then
        error stop "dataset%number /= 194"
    end if
    if (dataset%hall_number /= 488) then
        error stop "dataset%hall_number /= 488"
    end if
    if (moyo_to_string(dataset%hm_symbol) /= "P 6_3/m m c") then
        error stop "dataset%hm_symbol /= P 6_3/m m c"
    end if

    ! Input cell
    print *, "dataset%num_atoms: ", dataset%num_atoms
    if (dataset%num_atoms /= 2) then
        error stop "dataset%num_atoms /= 2"
    end if

    ! Symmetry operations in the input cell
    print *, "dataset%operations%num_operations: ", dataset%operations%num_operations
    if (dataset%operations%num_operations /= 24) then
        error stop "dataset%operations%num_operations /= 24"
    end if

    ! First operation must be the identity
    call c_f_pointer(dataset%operations%rotations, rot, &
                     [3, 3, int(dataset%operations%num_operations)])
    print *, "dataset%operations%rotations(:, :, 1):"
    do i = 1, 3
        print *, (rot(i, j, 1), j = 1, 3)
    end do
    if (any(rot(:, :, 1) /= identity)) then
        error stop "first operation is not the identity matrix"
    end if

    ! Site symmetry
    call c_f_pointer(dataset%orbits, orbits, [2])
    print *, "dataset%orbits: ", orbits
    if (any(orbits /= [0_c_int32_t, 0_c_int32_t])) then
        error stop "dataset%orbits /= [0, 0]"
    end if

    print *, "dataset%wyckoffs: ", moyo_to_string(dataset%wyckoffs)
    if (moyo_to_string(dataset%wyckoffs) /= "cc") then
        error stop "dataset%wyckoffs /= cc"
    end if

    ! Standardized cell
    print *, "dataset%pearson_symbol: ", moyo_to_string(dataset%pearson_symbol)
    if (moyo_to_string(dataset%pearson_symbol) /= "hP2") then
        error stop "dataset%pearson_symbol /= hP2"
    end if

    print *, "dataset%std_cell%num_atoms: ", dataset%std_cell%num_atoms
    if (dataset%std_cell%num_atoms /= 2) then
        error stop "dataset%std_cell%num_atoms /= 2"
    end if

    ! Primitive standardized cell
    call c_f_pointer(dataset%mapping_std_prim, mapping_std_prim, [2])
    print *, "dataset%mapping_std_prim: ", mapping_std_prim
    if (any(mapping_std_prim /= [0_c_int32_t, 1_c_int32_t])) then
        error stop "dataset%mapping_std_prim /= [0, 1]"
    end if

    call moyo_dataset_free(dataset_ptr)

    ! Space-group type
    sgt_ptr = moyo_space_group_type_new(194_c_int32_t)
    if (.not. c_associated(sgt_ptr)) then
        error stop "moyo_space_group_type_new(194) returned NULL"
    end if
    call c_f_pointer(sgt_ptr, sgt)
    print *, "sgt%number: ", sgt%number
    print *, "sgt%crystal_system: ", moyo_to_string(sgt%crystal_system)
    if (sgt%number /= 194) then
        error stop "sgt%number /= 194"
    end if
    if (moyo_to_string(sgt%crystal_system) /= "Hexagonal") then
        error stop "sgt%crystal_system /= Hexagonal"
    end if
    call moyo_space_group_type_free(sgt_ptr)

    if (c_associated(moyo_space_group_type_new(231_c_int32_t))) then
        error stop "moyo_space_group_type_new(231) did not return NULL"
    end if
    if (c_associated(moyo_space_group_type_new(0_c_int32_t))) then
        error stop "moyo_space_group_type_new(0) did not return NULL"
    end if
    if (c_associated(moyo_space_group_type_new(-huge(1_c_int32_t) - 1_c_int32_t))) then
        error stop "moyo_space_group_type_new(INT32_MIN) did not return NULL"
    end if
    if (c_associated(moyo_hall_symbol_entry_new(-huge(1_c_int32_t) - 1_c_int32_t))) then
        error stop "moyo_hall_symbol_entry_new(INT32_MIN) did not return NULL"
    end if

    print *, "test_moyof_dataset: all checks passed"
end program test_moyof_dataset
