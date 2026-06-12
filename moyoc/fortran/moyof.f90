! Fortran interface to moyoc, the C bindings for moyo.
!
! Conventions for crossing the C boundary:
! - `basis(:, i)` is the `i`th basis vector (the Fortran array is the
!   transpose view of the row-major C `basis[3][3]`, so columns are vectors).
! - `positions(:, i)` are the fractional coordinates of the `i`th site.
! - 3x3 matrices stored in datasets (`std_linear`, `rotations(:, :, k)`, ...)
!   are likewise transposed with respect to their C/Rust definition.
! - Functions returning `type(c_ptr)` give NULL (`c_associated` is false) on
!   failure; map the result with `c_f_pointer` and release it with the
!   matching `*_free` routine.
! - C strings (`type(c_ptr)` components) are converted with `moyo_to_string`.
module moyo
    use, intrinsic :: iso_c_binding, only: c_ptr, c_char, c_double, c_bool, &
                                           c_int32_t, c_size_t, c_null_char, &
                                           c_associated, c_f_pointer
    implicit none
    private

    ! ------------------------------------------------------------------------
    ! Setting preferences (values of the C `MoyoSetting` enum)
    ! ------------------------------------------------------------------------
    !> Hall number specified by the `hall_number` argument
    integer(c_int32_t), parameter, public :: MOYO_SETTING_HALL_NUMBER = 0
    !> The setting of the smallest Hall number
    integer(c_int32_t), parameter, public :: MOYO_SETTING_SPGLIB = 1
    !> Unique axis b, cell choice 1 for monoclinic, hexagonal axes for
    !> rhombohedral, and origin choice 2 for centrosymmetric space groups
    integer(c_int32_t), parameter, public :: MOYO_SETTING_STANDARD = 2

    ! Values of the C `MoyoLayerSetting` enum
    !> Layer Hall number specified by the `hall_number` argument (1 - 116)
    integer(c_int32_t), parameter, public :: MOYO_LAYER_SETTING_HALL_NUMBER = 0
    !> The setting of the smallest layer Hall number for each layer group
    integer(c_int32_t), parameter, public :: MOYO_LAYER_SETTING_SPGLIB = 1
    !> BCS standard choice (origin choice 2 for layer groups 52, 62, 64)
    integer(c_int32_t), parameter, public :: MOYO_LAYER_SETTING_STANDARD = 2

    ! ------------------------------------------------------------------------
    ! Interoperable derived types mirroring moyoc.h
    ! ------------------------------------------------------------------------
    !> Symmetry operations. `rotations` points to `num_operations` int32 3x3
    !> matrices and `translations` to `num_operations` real vectors; map them
    !> with `c_f_pointer(ops%rotations, rot, [3, 3, int(ops%num_operations)])`.
    type, bind(c), public :: MoyoOperations
        type(c_ptr) :: rotations
        type(c_ptr) :: translations
        integer(c_int32_t) :: num_operations
    end type MoyoOperations

    !> Magnetic symmetry operations; `time_reversals` points to
    !> `num_operations` `logical(c_bool)` values.
    type, bind(c), public :: MoyoMagneticOperations
        type(c_ptr) :: rotations
        type(c_ptr) :: translations
        type(c_ptr) :: time_reversals
        integer(c_int32_t) :: num_operations
    end type MoyoMagneticOperations

    !> Crystal structure. `basis(:, i)` is the `i`th basis vector;
    !> `positions` points to `num_atoms` fractional coordinate triplets and
    !> `numbers` to `num_atoms` atomic numbers.
    type, bind(c), public :: MoyoCell
        real(c_double) :: basis(3, 3)
        type(c_ptr) :: positions
        type(c_ptr) :: numbers
        integer(c_int32_t) :: num_atoms
    end type MoyoCell

    !> Collinear magnetic structure; `magnetic_moments` points to
    !> `cell%num_atoms` scalar moments.
    type, bind(c), public :: MoyoCollinearMagneticCell
        type(MoyoCell) :: cell
        type(c_ptr) :: magnetic_moments
    end type MoyoCollinearMagneticCell

    !> Non-collinear magnetic structure; `magnetic_moments` points to
    !> `cell%num_atoms` Cartesian moment vectors.
    type, bind(c), public :: MoyoNonCollinearMagneticCell
        type(MoyoCell) :: cell
        type(c_ptr) :: magnetic_moments
    end type MoyoNonCollinearMagneticCell

    !> A dataset of the symmetry analysis, created by `moyo_dataset_new` and
    !> freed by `moyo_dataset_free`. See moyoc.h for field documentation.
    type, bind(c), public :: MoyoDataset
        integer(c_int32_t) :: number
        integer(c_int32_t) :: hall_number
        type(c_ptr) :: hm_symbol
        integer(c_int32_t) :: num_atoms
        type(MoyoOperations) :: operations
        type(c_ptr) :: orbits
        type(c_ptr) :: wyckoffs
        type(c_ptr) :: site_symmetry_symbols
        type(MoyoCell) :: std_cell
        real(c_double) :: std_linear(3, 3)
        real(c_double) :: std_origin_shift(3)
        real(c_double) :: std_rotation_matrix(3, 3)
        type(c_ptr) :: pearson_symbol
        type(MoyoCell) :: prim_std_cell
        real(c_double) :: prim_std_linear(3, 3)
        real(c_double) :: prim_std_origin_shift(3)
        type(c_ptr) :: mapping_std_prim
        real(c_double) :: symprec
        real(c_double) :: angle_tolerance
    end type MoyoDataset

    !> A dataset of the layer-group symmetry analysis, created by
    !> `moyo_layer_dataset_new` and freed by `moyo_layer_dataset_free`.
    type, bind(c), public :: MoyoLayerDataset
        integer(c_int32_t) :: number
        integer(c_int32_t) :: hall_number
        type(c_ptr) :: hm_symbol
        integer(c_int32_t) :: num_atoms
        type(MoyoOperations) :: operations
        type(c_ptr) :: orbits
        type(c_ptr) :: wyckoffs
        type(c_ptr) :: site_symmetry_symbols
        type(MoyoCell) :: std_cell
        real(c_double) :: std_linear(3, 3)
        real(c_double) :: std_origin_shift(3)
        real(c_double) :: std_rotation_matrix(3, 3)
        type(c_ptr) :: pearson_symbol
        type(MoyoCell) :: prim_std_cell
        real(c_double) :: prim_std_linear(3, 3)
        real(c_double) :: prim_std_origin_shift(3)
        type(c_ptr) :: mapping_std_prim
        real(c_double) :: symprec
        real(c_double) :: angle_tolerance
    end type MoyoLayerDataset

    !> A dataset of the magnetic symmetry analysis of a collinear magnetic
    !> structure, created by `moyo_collinear_magnetic_dataset_new` and freed
    !> by `moyo_collinear_magnetic_dataset_free`.
    type, bind(c), public :: MoyoCollinearMagneticDataset
        integer(c_int32_t) :: uni_number
        integer(c_int32_t) :: num_atoms
        type(MoyoMagneticOperations) :: magnetic_operations
        type(c_ptr) :: orbits
        type(MoyoCollinearMagneticCell) :: std_mag_cell
        real(c_double) :: std_linear(3, 3)
        real(c_double) :: std_origin_shift(3)
        real(c_double) :: std_rotation_matrix(3, 3)
        type(MoyoCollinearMagneticCell) :: prim_std_mag_cell
        real(c_double) :: prim_std_linear(3, 3)
        real(c_double) :: prim_std_origin_shift(3)
        type(c_ptr) :: mapping_std_prim
        real(c_double) :: symprec
        real(c_double) :: angle_tolerance
        real(c_double) :: mag_symprec
    end type MoyoCollinearMagneticDataset

    !> A dataset of the magnetic symmetry analysis of a non-collinear magnetic
    !> structure, created by `moyo_noncollinear_magnetic_dataset_new` and
    !> freed by `moyo_noncollinear_magnetic_dataset_free`.
    type, bind(c), public :: MoyoNonCollinearMagneticDataset
        integer(c_int32_t) :: uni_number
        integer(c_int32_t) :: num_atoms
        type(MoyoMagneticOperations) :: magnetic_operations
        type(c_ptr) :: orbits
        type(MoyoNonCollinearMagneticCell) :: std_mag_cell
        real(c_double) :: std_linear(3, 3)
        real(c_double) :: std_origin_shift(3)
        real(c_double) :: std_rotation_matrix(3, 3)
        type(MoyoNonCollinearMagneticCell) :: prim_std_mag_cell
        real(c_double) :: prim_std_linear(3, 3)
        real(c_double) :: prim_std_origin_shift(3)
        type(c_ptr) :: mapping_std_prim
        real(c_double) :: symprec
        real(c_double) :: angle_tolerance
        real(c_double) :: mag_symprec
    end type MoyoNonCollinearMagneticDataset

    !> Hall symbol database entry, created by `moyo_hall_symbol_entry_new`
    !> and freed by `moyo_hall_symbol_entry_free`.
    type, bind(c), public :: MoyoHallSymbolEntry
        integer(c_int32_t) :: hall_number
        integer(c_int32_t) :: number
        integer(c_int32_t) :: arithmetic_number
        type(c_ptr) :: setting
        type(c_ptr) :: hall_symbol
        type(c_ptr) :: hm_short
        type(c_ptr) :: hm_full
        type(c_ptr) :: centering
    end type MoyoHallSymbolEntry

    !> Layer Hall symbol database entry, created by
    !> `moyo_layer_hall_symbol_entry_new` and freed by
    !> `moyo_layer_hall_symbol_entry_free`.
    type, bind(c), public :: MoyoLayerHallSymbolEntry
        integer(c_int32_t) :: hall_number
        integer(c_int32_t) :: number
        integer(c_int32_t) :: arithmetic_number
        type(c_ptr) :: setting
        type(c_ptr) :: hall_symbol
        type(c_ptr) :: hm_short
        type(c_ptr) :: hm_full
        type(c_ptr) :: centering
    end type MoyoLayerHallSymbolEntry

    !> Space-group type information, created by `moyo_space_group_type_new`
    !> and freed by `moyo_space_group_type_free`.
    type, bind(c), public :: MoyoSpaceGroupType
        integer(c_int32_t) :: number
        type(c_ptr) :: hm_short
        type(c_ptr) :: hm_full
        integer(c_int32_t) :: arithmetic_number
        type(c_ptr) :: arithmetic_symbol
        type(c_ptr) :: geometric_crystal_class
        type(c_ptr) :: crystal_system
        type(c_ptr) :: bravais_class
        type(c_ptr) :: lattice_system
        type(c_ptr) :: crystal_family
    end type MoyoSpaceGroupType

    !> Layer-group type information, created by `moyo_layer_group_type_new`
    !> and freed by `moyo_layer_group_type_free`.
    type, bind(c), public :: MoyoLayerGroupType
        integer(c_int32_t) :: number
        type(c_ptr) :: hm_short
        type(c_ptr) :: hm_full
        integer(c_int32_t) :: arithmetic_number
        type(c_ptr) :: arithmetic_symbol
        type(c_ptr) :: geometric_crystal_class
        type(c_ptr) :: bravais_class
        type(c_ptr) :: lattice_system
    end type MoyoLayerGroupType

    !> Magnetic space-group type information, created by
    !> `moyo_magnetic_space_group_type_new` and freed by
    !> `moyo_magnetic_space_group_type_free`.
    type, bind(c), public :: MoyoMagneticSpaceGroupType
        integer(c_int32_t) :: uni_number
        integer(c_int32_t) :: litvin_number
        type(c_ptr) :: bns_number
        type(c_ptr) :: og_number
        integer(c_int32_t) :: number
        integer(c_int32_t) :: construct_type
    end type MoyoMagneticSpaceGroupType

    public :: moyo_version
    public :: moyo_dataset_new, moyo_dataset_free
    public :: moyo_layer_dataset_new, moyo_layer_dataset_free
    public :: moyo_collinear_magnetic_dataset_new, moyo_collinear_magnetic_dataset_free
    public :: moyo_noncollinear_magnetic_dataset_new, moyo_noncollinear_magnetic_dataset_free
    public :: moyo_operations_from_number, moyo_operations_from_layer_number
    public :: moyo_operations_free
    public :: moyo_magnetic_operations_from_uni_number, moyo_magnetic_operations_free
    public :: moyo_hall_symbol_entry_new, moyo_hall_symbol_entry_free
    public :: moyo_layer_hall_symbol_entry_new, moyo_layer_hall_symbol_entry_free
    public :: moyo_space_group_type_new, moyo_space_group_type_free
    public :: moyo_layer_group_type_new, moyo_layer_group_type_free
    public :: moyo_magnetic_space_group_type_new, moyo_magnetic_space_group_type_free
    public :: moyo_to_string

    interface
        !> Version string of moyoc; convert with `moyo_to_string`. The
        !> returned string is statically allocated and must not be freed.
        function moyo_version() bind(c, name="moyo_version") result(version)
            import :: c_ptr
            type(c_ptr) :: version
        end function moyo_version

        !> Analyze the symmetry of the given cell; returns a pointer to a
        !> `MoyoDataset` or NULL on failure. Pass a negative
        !> `angle_tolerance` to use the default tolerance; `hall_number` is
        !> only used with `MOYO_SETTING_HALL_NUMBER`.
        function moyo_dataset_new(basis, positions, numbers, num_atoms, symprec, &
                                  angle_tolerance, setting, hall_number, rotate_basis) &
            bind(c, name="moyo_dataset_new") result(dataset)
            import :: c_ptr, c_double, c_int32_t, c_bool
            real(c_double), intent(in) :: basis(3, 3)
            real(c_double), intent(in) :: positions(3, *)
            integer(c_int32_t), intent(in) :: numbers(*)
            integer(c_int32_t), value :: num_atoms
            real(c_double), value :: symprec
            real(c_double), value :: angle_tolerance
            integer(c_int32_t), value :: setting
            integer(c_int32_t), value :: hall_number
            logical(c_bool), value :: rotate_basis
            type(c_ptr) :: dataset
        end function moyo_dataset_new

        !> Free a dataset created by `moyo_dataset_new`.
        subroutine moyo_dataset_free(dataset) bind(c, name="moyo_dataset_free")
            import :: c_ptr
            type(c_ptr), value :: dataset
        end subroutine moyo_dataset_free

        !> Analyze the layer-group symmetry of the given cell; returns a
        !> pointer to a `MoyoLayerDataset` or NULL on failure. The third
        !> basis vector must be the aperiodic stacking direction.
        function moyo_layer_dataset_new(basis, positions, numbers, num_atoms, symprec, &
                                        angle_tolerance, setting, hall_number, rotate_basis) &
            bind(c, name="moyo_layer_dataset_new") result(dataset)
            import :: c_ptr, c_double, c_int32_t, c_bool
            real(c_double), intent(in) :: basis(3, 3)
            real(c_double), intent(in) :: positions(3, *)
            integer(c_int32_t), intent(in) :: numbers(*)
            integer(c_int32_t), value :: num_atoms
            real(c_double), value :: symprec
            real(c_double), value :: angle_tolerance
            integer(c_int32_t), value :: setting
            integer(c_int32_t), value :: hall_number
            logical(c_bool), value :: rotate_basis
            type(c_ptr) :: dataset
        end function moyo_layer_dataset_new

        !> Free a dataset created by `moyo_layer_dataset_new`.
        subroutine moyo_layer_dataset_free(dataset) bind(c, name="moyo_layer_dataset_free")
            import :: c_ptr
            type(c_ptr), value :: dataset
        end subroutine moyo_layer_dataset_free

        !> Analyze the magnetic symmetry of a collinear magnetic cell;
        !> returns a pointer to a `MoyoCollinearMagneticDataset` or NULL on
        !> failure. Pass a negative `mag_symprec` to reuse `symprec`.
        function moyo_collinear_magnetic_dataset_new(basis, positions, numbers, &
                                                     magnetic_moments, num_atoms, symprec, &
                                                     angle_tolerance, mag_symprec, is_axial, &
                                                     rotate_basis) &
            bind(c, name="moyo_collinear_magnetic_dataset_new") result(dataset)
            import :: c_ptr, c_double, c_int32_t, c_bool
            real(c_double), intent(in) :: basis(3, 3)
            real(c_double), intent(in) :: positions(3, *)
            integer(c_int32_t), intent(in) :: numbers(*)
            real(c_double), intent(in) :: magnetic_moments(*)
            integer(c_int32_t), value :: num_atoms
            real(c_double), value :: symprec
            real(c_double), value :: angle_tolerance
            real(c_double), value :: mag_symprec
            logical(c_bool), value :: is_axial
            logical(c_bool), value :: rotate_basis
            type(c_ptr) :: dataset
        end function moyo_collinear_magnetic_dataset_new

        !> Free a dataset created by `moyo_collinear_magnetic_dataset_new`.
        subroutine moyo_collinear_magnetic_dataset_free(dataset) &
            bind(c, name="moyo_collinear_magnetic_dataset_free")
            import :: c_ptr
            type(c_ptr), value :: dataset
        end subroutine moyo_collinear_magnetic_dataset_free

        !> Analyze the magnetic symmetry of a non-collinear magnetic cell
        !> (`magnetic_moments(:, i)` is the Cartesian moment of site `i`);
        !> returns a pointer to a `MoyoNonCollinearMagneticDataset` or NULL
        !> on failure.
        function moyo_noncollinear_magnetic_dataset_new(basis, positions, numbers, &
                                                        magnetic_moments, num_atoms, symprec, &
                                                        angle_tolerance, mag_symprec, is_axial, &
                                                        rotate_basis) &
            bind(c, name="moyo_noncollinear_magnetic_dataset_new") result(dataset)
            import :: c_ptr, c_double, c_int32_t, c_bool
            real(c_double), intent(in) :: basis(3, 3)
            real(c_double), intent(in) :: positions(3, *)
            integer(c_int32_t), intent(in) :: numbers(*)
            real(c_double), intent(in) :: magnetic_moments(3, *)
            integer(c_int32_t), value :: num_atoms
            real(c_double), value :: symprec
            real(c_double), value :: angle_tolerance
            real(c_double), value :: mag_symprec
            logical(c_bool), value :: is_axial
            logical(c_bool), value :: rotate_basis
            type(c_ptr) :: dataset
        end function moyo_noncollinear_magnetic_dataset_new

        !> Free a dataset created by `moyo_noncollinear_magnetic_dataset_new`.
        subroutine moyo_noncollinear_magnetic_dataset_free(dataset) &
            bind(c, name="moyo_noncollinear_magnetic_dataset_free")
            import :: c_ptr
            type(c_ptr), value :: dataset
        end subroutine moyo_noncollinear_magnetic_dataset_free

        !> Symmetry operations for a space-group ITA number (1 - 230);
        !> returns a pointer to `MoyoOperations` or NULL on invalid input.
        function moyo_operations_from_number(number, setting, hall_number, primitive) &
            bind(c, name="moyo_operations_from_number") result(operations)
            import :: c_ptr, c_int32_t, c_bool
            integer(c_int32_t), value :: number
            integer(c_int32_t), value :: setting
            integer(c_int32_t), value :: hall_number
            logical(c_bool), value :: primitive
            type(c_ptr) :: operations
        end function moyo_operations_from_number

        !> Symmetry operations for a layer-group number (1 - 80);
        !> returns a pointer to `MoyoOperations` or NULL on invalid input.
        function moyo_operations_from_layer_number(number, setting, hall_number, primitive) &
            bind(c, name="moyo_operations_from_layer_number") result(operations)
            import :: c_ptr, c_int32_t, c_bool
            integer(c_int32_t), value :: number
            integer(c_int32_t), value :: setting
            integer(c_int32_t), value :: hall_number
            logical(c_bool), value :: primitive
            type(c_ptr) :: operations
        end function moyo_operations_from_layer_number

        !> Free operations returned by `moyo_operations_from_number` or
        !> `moyo_operations_from_layer_number` (not operations embedded in a
        !> dataset).
        subroutine moyo_operations_free(operations) bind(c, name="moyo_operations_free")
            import :: c_ptr
            type(c_ptr), value :: operations
        end subroutine moyo_operations_free

        !> Magnetic symmetry operations for a UNI number (1 - 1651); returns
        !> a pointer to `MoyoMagneticOperations` or NULL on invalid input.
        function moyo_magnetic_operations_from_uni_number(uni_number, primitive) &
            bind(c, name="moyo_magnetic_operations_from_uni_number") result(magnetic_operations)
            import :: c_ptr, c_int32_t, c_bool
            integer(c_int32_t), value :: uni_number
            logical(c_bool), value :: primitive
            type(c_ptr) :: magnetic_operations
        end function moyo_magnetic_operations_from_uni_number

        !> Free magnetic operations returned by
        !> `moyo_magnetic_operations_from_uni_number` (not operations
        !> embedded in a dataset).
        subroutine moyo_magnetic_operations_free(magnetic_operations) &
            bind(c, name="moyo_magnetic_operations_free")
            import :: c_ptr
            type(c_ptr), value :: magnetic_operations
        end subroutine moyo_magnetic_operations_free

        !> Hall symbol entry for `hall_number` (1 - 530); returns a pointer
        !> to `MoyoHallSymbolEntry` or NULL if out of range.
        function moyo_hall_symbol_entry_new(hall_number) &
            bind(c, name="moyo_hall_symbol_entry_new") result(entry)
            import :: c_ptr, c_int32_t
            integer(c_int32_t), value :: hall_number
            type(c_ptr) :: entry
        end function moyo_hall_symbol_entry_new

        !> Free an entry created by `moyo_hall_symbol_entry_new`.
        subroutine moyo_hall_symbol_entry_free(entry) &
            bind(c, name="moyo_hall_symbol_entry_free")
            import :: c_ptr
            type(c_ptr), value :: entry
        end subroutine moyo_hall_symbol_entry_free

        !> Layer Hall symbol entry for `hall_number` (1 - 116); returns a
        !> pointer to `MoyoLayerHallSymbolEntry` or NULL if out of range.
        function moyo_layer_hall_symbol_entry_new(hall_number) &
            bind(c, name="moyo_layer_hall_symbol_entry_new") result(entry)
            import :: c_ptr, c_int32_t
            integer(c_int32_t), value :: hall_number
            type(c_ptr) :: entry
        end function moyo_layer_hall_symbol_entry_new

        !> Free an entry created by `moyo_layer_hall_symbol_entry_new`.
        subroutine moyo_layer_hall_symbol_entry_free(entry) &
            bind(c, name="moyo_layer_hall_symbol_entry_free")
            import :: c_ptr
            type(c_ptr), value :: entry
        end subroutine moyo_layer_hall_symbol_entry_free

        !> Space-group type information for ITA `number` (1 - 230); returns
        !> a pointer to `MoyoSpaceGroupType` or NULL if out of range.
        function moyo_space_group_type_new(number) &
            bind(c, name="moyo_space_group_type_new") result(space_group_type)
            import :: c_ptr, c_int32_t
            integer(c_int32_t), value :: number
            type(c_ptr) :: space_group_type
        end function moyo_space_group_type_new

        !> Free a space-group type created by `moyo_space_group_type_new`.
        subroutine moyo_space_group_type_free(space_group_type) &
            bind(c, name="moyo_space_group_type_free")
            import :: c_ptr
            type(c_ptr), value :: space_group_type
        end subroutine moyo_space_group_type_free

        !> Layer-group type information for layer-group `number` (1 - 80);
        !> returns a pointer to `MoyoLayerGroupType` or NULL if out of range.
        function moyo_layer_group_type_new(number) &
            bind(c, name="moyo_layer_group_type_new") result(layer_group_type)
            import :: c_ptr, c_int32_t
            integer(c_int32_t), value :: number
            type(c_ptr) :: layer_group_type
        end function moyo_layer_group_type_new

        !> Free a layer-group type created by `moyo_layer_group_type_new`.
        subroutine moyo_layer_group_type_free(layer_group_type) &
            bind(c, name="moyo_layer_group_type_free")
            import :: c_ptr
            type(c_ptr), value :: layer_group_type
        end subroutine moyo_layer_group_type_free

        !> Magnetic space-group type information for `uni_number`
        !> (1 - 1651); returns a pointer to `MoyoMagneticSpaceGroupType` or
        !> NULL if out of range.
        function moyo_magnetic_space_group_type_new(uni_number) &
            bind(c, name="moyo_magnetic_space_group_type_new") result(magnetic_space_group_type)
            import :: c_ptr, c_int32_t
            integer(c_int32_t), value :: uni_number
            type(c_ptr) :: magnetic_space_group_type
        end function moyo_magnetic_space_group_type_new

        !> Free a magnetic space-group type created by
        !> `moyo_magnetic_space_group_type_new`.
        subroutine moyo_magnetic_space_group_type_free(magnetic_space_group_type) &
            bind(c, name="moyo_magnetic_space_group_type_free")
            import :: c_ptr
            type(c_ptr), value :: magnetic_space_group_type
        end subroutine moyo_magnetic_space_group_type_free
    end interface

    interface
        function moyo_c_strlen(s) bind(c, name="strlen") result(length)
            import :: c_ptr, c_size_t
            type(c_ptr), value :: s
            integer(c_size_t) :: length
        end function moyo_c_strlen
    end interface

contains

    !> Convert a NUL-terminated C string (`type(c_ptr)`) into a Fortran
    !> character value. Returns an empty string for a NULL pointer.
    function moyo_to_string(c_str) result(f_str)
        type(c_ptr), intent(in) :: c_str
        character(len=:), allocatable :: f_str
        character(kind=c_char), pointer :: chars(:)
        integer :: i, length

        if (.not. c_associated(c_str)) then
            f_str = ""
            return
        end if

        length = int(moyo_c_strlen(c_str))
        call c_f_pointer(c_str, chars, [length])
        allocate (character(len=length) :: f_str)
        do i = 1, length
            f_str(i:i) = chars(i)
        end do
    end function moyo_to_string

end module moyo
