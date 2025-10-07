module hdf5_utils_m
   !> A set of high level wrapper subroutines for HDF5, taken directly from
   !> [https://github.com/jterwin/HDF5_utils](https://github.com/jterwin/HDF5_utils)
   use, intrinsic :: iso_fortran_env, only: dp => real64, sp => real32
   ! HDF5 types
   use hdf5, only: SIZE_T, HID_T, HSIZE_T
   ! HDF5 datatypes
   use hdf5, only: H5T_NATIVE_INTEGER, H5T_NATIVE_DOUBLE, H5T_NATIVE_CHARACTER
   use hdf5, only: h5tclose_f, h5tcopy_f, h5tset_size_f
   ! HDF5 general
   use hdf5, only: h5open_f, h5close_f
   ! HDF5 constants
   use hdf5, only: H5F_ACC_RDONLY_F, H5F_ACC_RDWR_F, H5F_ACC_TRUNC_F, H5F_ACC_EXCL_F
   ! HDF5 objects
   use hdf5, only: h5oopen_f, h5oclose_f
   ! HDF5 files
   use hdf5, only: h5fopen_f, h5fcreate_f, h5fclose_f
   ! HDF5 groups
   use hdf5, only: h5gopen_f, h5gcreate_f, h5gclose_f
   ! HDF5 attributes
   use hdf5, only: h5aopen_f, h5acreate_f, h5awrite_f, h5aread_f, h5aclose_f
   ! HDF5 links
   use hdf5, only: h5lexists_f
   ! HDF5 datasets
   use hdf5, only: h5dopen_f, h5dcreate_f, h5dget_space_f, h5dread_f, h5dwrite_f, h5dclose_f
   ! HDF5 dataspaces
   use hdf5, only: H5S_SCALAR_F, H5S_SELECT_SET_F
   use hdf5, only: h5screate_f, h5screate_simple_f, h5sget_simple_extent_dims_f, h5sget_simple_extent_ndims_f, &
      h5sget_simple_extent_npoints_f, h5sselect_hyperslab_f, h5sclose_f
   ! HDF5 dimension scales
   use h5ds, only: h5dsset_scale_f, h5dsattach_scale_f, h5dsset_label_f
   implicit none

   private

   public :: hdf_open_file, hdf_close_file
   public :: hdf_create_group, hdf_open_group, hdf_close_group
   public :: hdf_exists, hdf_get_rank, hdf_get_dims
   public :: hdf_set_data_scale, hdf_attach_data_scale, hdf_label_dim
   public :: hdf_write_dataset, hdf_read_dataset
   public :: hdf_write_attribute, hdf_read_attribute
   public :: hdf_create_dataset
   public :: hdf_write_vector_to_dataset, hdf_read_vector_from_dataset
   public :: HID_T, hdf_set_print_messages

   interface hdf_write_dataset
      !> Generic interface to write a dataset
      !>  Supported types
      !>   - integers (scalar and 1d-6d arrays)
      !>   - doubles (scalar and 1d-6d arrays)
      !>   - reals (scalar and 1d-6d arrays)
      !>  - string (scalar and 1d-6d arrays)
      module procedure hdf_write_dataset_integer_0
      module procedure hdf_write_dataset_integer_1
      module procedure hdf_write_dataset_integer_2
      module procedure hdf_write_dataset_integer_3
      module procedure hdf_write_dataset_integer_4
      module procedure hdf_write_dataset_integer_5
      module procedure hdf_write_dataset_integer_6
      module procedure hdf_write_dataset_double_0
      module procedure hdf_write_dataset_double_1
      module procedure hdf_write_dataset_double_2
      module procedure hdf_write_dataset_double_3
      module procedure hdf_write_dataset_double_4
      module procedure hdf_write_dataset_double_5
      module procedure hdf_write_dataset_double_6
   end interface hdf_write_dataset

   interface hdf_write_vector_to_dataset
      !> Generic interface to write a vector to dataset
      !>  The vector is is written out in along the fast dimension
      !>  (column oriented in FORTRAN, row oriented in the HDF5 file).
      !>  So the vector should have the same length as the first dimension
      !>  of the dataset, and the offset should agree with dims(2:rank) of
      !>  of the dataset.
      !>  Supported types
      !>   - integers (scalar and 1d-6d arrays)
      !>   - doubles (scalar and 1d-6d arrays)
      !>   - reals (scalar and 1d-6d arrays)
      !>   - string (scalar and 1d-6d arrays)
      module procedure hdf_write_vector_to_dataset_double
      module procedure hdf_write_vector_to_dataset_integer
   end interface hdf_write_vector_to_dataset


   interface hdf_read_dataset
      !> Generic interface to read a dataset of doubles
      !>  Supported types:
      !>   - integers (scalar and 1d-6d arrays)
      !>   - doubles (scalar and 1d-6d arrays)
      !>   - reals (scalar and 1d-6d arrays)
      !>   - string (scalar and 1d-6d arrays)
      module procedure hdf_read_dataset_integer_0
      module procedure hdf_read_dataset_integer_1
      module procedure hdf_read_dataset_integer_2
      module procedure hdf_read_dataset_integer_3
      module procedure hdf_read_dataset_integer_4
      module procedure hdf_read_dataset_integer_5
      module procedure hdf_read_dataset_integer_6
      module procedure hdf_read_dataset_double_0
      module procedure hdf_read_dataset_double_1
      module procedure hdf_read_dataset_double_2
      module procedure hdf_read_dataset_double_3
      module procedure hdf_read_dataset_double_4
      module procedure hdf_read_dataset_double_5
      module procedure hdf_read_dataset_double_6
   end interface hdf_read_dataset


   interface hdf_read_vector_from_dataset
      !> Generic interface to read a vector from a dataset
      !>  The vector is is read in along the fast dimension
      !>  (column oriented in FORTRAN, row oriented in the HDF5 file).
      !>  So the vector should have the same length as the first dimension
      !>  of the dataset, and the offset should agree with dims(2:rank) of
      !>  of the dataset.
      !>  Supported types
      !>   - integers (scalar and 1d-6d arrays)
      !>   - doubles (scalar and 1d-6d arrays)
      !>   - reals (scalar and 1d-6d arrays)
      !>   - string (scalar and 1d-6d arrays)
      module procedure hdf_read_vector_from_dataset_double
      module procedure hdf_read_vector_from_dataset_integer
   end interface hdf_read_vector_from_dataset

   interface hdf_write_attribute
      !>   Generic interface to write attribute
      !>  Supported types
      !>   - integers (scalar and 1d-6d arrays)
      !>   - doubles (scalar and 1d-6d arrays)
      !>   - string (scalar)
      !>  - reals (scalar and 1d-6d arrays)
      module procedure hdf_write_attr_integer_0
      module procedure hdf_write_attr_integer_1
      module procedure hdf_write_attr_double_0
      module procedure hdf_write_attr_double_1
      module procedure hdf_write_attr_string
   end interface hdf_write_attribute

   interface hdf_read_attribute
      !>   Generic interface to read attribute
      !>  Supported types
      !>   - integers (scalar and 1d-6d arrays)
      !>   - doubles (scalar and 1d-6d arrays)
      !>   - string (scalar and 1d-6d arrays)
      !>  - reals (scalar)
      module procedure hdf_read_attr_integer_0
      module procedure hdf_read_attr_integer_1
      module procedure hdf_read_attr_double_0
      module procedure hdf_read_attr_double_1
      module procedure hdf_read_attr_string
   end interface hdf_read_attribute

   !
   logical :: hdf_print_messages = .false.

contains


subroutine hdf_set_print_messages(val_print_messages)
   !>   Sets the value of hdf_print_messages
   !>  By default, hdf_print_messages = .false. By setting it
   !>  to .true., some messages are printed detailing what hdf_utils
   !>  is doing.
   logical, intent(in) :: val_print_messages  !>  new value for hdf_print_messages

   hdf_print_messages = val_print_messages

end subroutine hdf_set_print_messages


function hdf_exists(loc_id, obj_name) result(exists)
   !> Check if location exists.
   !> Also checks is intemediate paths exists in a safe way.

   integer(HID_T), intent(in) :: loc_id       !> local id
   character(len=*), intent(in) :: obj_name   !> relative path to object

   logical :: exists  !> .TRUE. if everything exists, .FALSE. otherwise

   integer :: hdferror, pos, cpos, str_len

   if (hdf_print_messages) then
      write(*,'(A,A)') "->hdf_exists: " // obj_name
   end if

   ! check intermediate paths (subgroups)
   str_len = len_trim(obj_name)
   cpos = 0
   do
      !start = cpos + 1
      !write(*,*) start, str_len, obj_name(start:str_len)

      pos = index(obj_name(cpos+1:str_len), "/")

      ! no subgroup found
      if (pos == 0) exit

      ! check subgroup
      cpos = cpos + pos
      call h5lexists_f(loc_id, obj_name(1:cpos-1), exists, hdferror)
      !write(*,*) obj_name(1:cpos-1), exists

      ! return if intermediate path fails
      if (exists .eqv. .false.) then
         if (hdf_print_messages) then
            write(*,'(A,A,A)') "--->hdf_exists: subpath '", obj_name(1:cpos-1), "' does not exist, return false"
         end if
         exists = .false.
         return
      end if

   end do

   ! check object (unless obj_name ended with "/"
   if (cpos /= str_len) then
      call h5lexists_f(loc_id, obj_name, exists, hdferror)
      !write(*,*) obj_name, exists
      if (exists .eqv. .false.) then
         if (hdf_print_messages) then
            write(*,'(A,A,A)') "--->hdf_exists: object '", obj_name, "' does not exist, return false"
         end if
         exists = .false.
         return
      end if
   end if

   exists = .true.
   return

end function hdf_exists


subroutine hdf_open_file(file_id, filename, STATUS, ACTION)
   !> Opens file and return identifier
   !>  | STATUS  | ACTION    | Description                          |
   !>  | :-----: | :-------: | :----------------------------------- |
   !>  | NEW     | na        | calls h5fcreate with H5F_ACC_TRUNC_F |
   !>  | REPLACE | na        | calls h5fcreate with H5F_ACC_EXCL_F  |
   !>  | OLD     | READ      | calls h5fopen with H5F_ACC_RDONLY_F  |
   !>  | OLD     | WRITE     | calls h5fopen with H5F_ACC_RDWR_F    |
   !>  | OLD     | READWRITE | calls h5fopen with H5F_ACC_RDWR_F    |
   integer(HID_T), intent(out) :: file_id            !> HDF5 id of the file
   character(len=*), intent(in) :: filename          !> the HDF5 filename
   character(len=*), optional, intent(in) :: STATUS  !> file status (OLD, NEW, REPLACE)
   character(len=*), optional, intent(in) :: ACTION  !> file action (READ, WRITE, READWRITE)

   integer :: hdferror
   character(len=16) :: status2, action2

   if (hdf_print_messages) then
      write(*,'(A)') "->hdf_open_file: " // trim(filename)
   end if

   ! open hdf5 interface
   call h5open_f(hdferror)
   !write(*,'(A20,I0)') "h5open: ", hdferror

   ! set defaults
   status2 = 'NEW'
   if (present(STATUS)) status2 = STATUS
   action2 = 'READWRITE'
   if (present(STATUS)) action2 = ACTION

   ! open/create hdf5 file
   if (status2 == 'OLD') then
      if (action2 == 'READ') then
         call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, hdferror)
      elseif ( (action2 == 'WRITE') .or. (action2 == 'READWRITE') ) then
         call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, hdferror)
      else
         write(*,*) "hdf_open: action = ", action2, " not supported."
         stop
      end if
   elseif (status2 == 'NEW') then
      call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, hdferror)
   elseif (status2 == 'REPLACE') then
      call EXECUTE_COMMAND_LINE("rm -f " // filename)
      call h5fcreate_f(filename, H5F_ACC_EXCL_F, file_id, hdferror)
   else
      write(*,*) "hdf_open: status = ", status2, " not supported."
      stop
   end if

   !write(*,'(A20,I0)') "h5fcreate: ", hdferror

end subroutine hdf_open_file



subroutine hdf_close_file(file_id)
   !>   Closes a hdf5 file
   integer(HID_T), intent(in) :: file_id  !> file id to be closed

   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "->hdf_close_file"
   end if

   call h5fclose_f(file_id, hdferror)
   !write(*,'(A20,I0)') "h5fclose: ", hdferror

   call h5close_f(hdferror)

end subroutine hdf_close_file


subroutine hdf_create_group(loc_id, group_name)
   !>   Create a new group
   integer(HID_T), intent(in) :: loc_id         !> location id where to put the group
   character(len=*), intent(in) :: group_name   !> name of the group

   integer(HID_T) :: grp_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "->hdf_create_group: " // trim(group_name)
   end if

   call h5gcreate_f(loc_id, group_name, grp_id, hdferror)
   !write(*,'(A20,I0)') "h5gcreate: ", hdferror

   call h5gclose_f(grp_id, hdferror)
   !write(*,'(A20,I0)') "h5gclose: ", hdferror

end subroutine hdf_create_group


subroutine hdf_open_group(loc_id, group_name, group_id)
   !>   Opens a group and returns the identifier
   integer(HID_T), intent(in) :: loc_id         !> location id where to put the group
   character(len=*), intent(in) :: group_name   !> name of the group
   integer(HID_T), intent(out) :: group_id      !> id for the group

   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A,A,A)') "->hdf_open_group: '" // trim(group_name) // "'"
   end if

   if (hdf_exists(loc_id, group_name)) then
      if (hdf_print_messages) then
         write(*,'(A,A,A)') "->hdf_open_group: opening group '" // trim(group_name) // "'"
      end if
      call h5gopen_f(loc_id, group_name, group_id, hdferror)
   else
      if (hdf_print_messages) then
         write(*,'(A,A,A)') "->hdf_open_group: group '" // trim(group_name) // "' does not exist, return with error"
      end if
      hdferror = -1
   end if

end subroutine hdf_open_group


subroutine hdf_close_group(group_id)
   !>   Close a group by identifier
   integer(HID_T), intent(in) :: group_id   !> id for the group

   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "->hdf_close_group"
   end if

   call h5gclose_f(group_id, hdferror)
   !write(*,'(A20,I0)') "h5gclose: ", hdferror

end subroutine hdf_close_group


subroutine hdf_get_rank(loc_id, dset_name, rank)
   !>   Get the rank of a dataset
   integer(HID_T), intent(in) :: loc_id        !> location id
   character(len=*), intent(in) :: dset_name   !> dataset name
   integer, intent(out) :: rank                !> rank of the dataset

   integer(HID_T) :: dset_id, dspace_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "->hdf_get_rank"
   end if

   ! open dataset
   call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

   ! get dataspace
   call h5dget_space_f(dset_id, dspace_id, hdferror)

   ! get rank (ndims)
   call h5sget_simple_extent_ndims_f(dspace_id, rank, hdferror)

   ! close id's
   call h5sclose_f(dspace_id, hdferror)
   call h5dclose_f(dset_id, hdferror)

end subroutine hdf_get_rank


subroutine hdf_get_dims(loc_id, dset_name, dims)
   !>   get the dimensions of a dataset
   integer(HID_T), intent(in) :: loc_id        !> location id
   character(len=*), intent(in) :: dset_name   !> name of dataset
   integer, intent(out) :: dims(:)             !> dimensions of the dataset

   integer(HID_T) :: dset_id, dspace_id
   integer :: rank
   integer(HSIZE_T) :: dset_dims(6), max_dims(6)
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "->hdf_get_dims"
   end if

   ! open dataset
   call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

   ! get dataspace
   call h5dget_space_f(dset_id, dspace_id, hdferror)

   ! get rank (ndims)
   call h5sget_simple_extent_ndims_f(dspace_id, rank, hdferror)

   ! get dims
   call h5sget_simple_extent_dims_f(dspace_id, dset_dims(1:rank), max_dims(1:rank), hdferror)
   dims(1:rank) = int(dset_dims(1:rank))

   ! close id's
   call h5sclose_f(dspace_id, hdferror)
   call h5dclose_f(dset_id, hdferror)

end subroutine hdf_get_dims

subroutine hdf_label_dim(loc_id, dset_name, dim_idx, dim_label)
   !>   Label a dimension of a dataset
   integer(HID_T), intent(in) :: loc_id        !> local id in file
   character(len=*), intent(in) :: dset_name   !> name of dataset
   integer, intent(in) :: dim_idx               !> index of dimension to label (0-based)
   character(len=*), intent(in) :: dim_label    !> label for the dimension

   integer(HID_T) :: dset_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A,A,A,I0,A)') "--->hdf_label_dim: ", trim(dset_name), " dim ", dim_idx, " to ", trim(dim_label)
   end if

   ! open dataset
   call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

   ! set label for dimension
   call h5dsset_label_f(dset_id, dim_idx, dim_label, hdferror)

   ! close id's
   call h5dclose_f(dset_id, hdferror)
end subroutine hdf_label_dim

subroutine hdf_set_data_scale(loc_id, dset_name, dscale_name)
   !>   Convert a dataset to a dimension scale, with optional name `dscale_name`
   integer(HID_T), intent(in) :: loc_id        !> local id in file
   character(len=*), intent(in) :: dset_name   !> name of dataset
   character(len=*), intent(in), optional :: dscale_name  !> name of dimension scale

   integer(HID_T) :: dset_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_set_data_scale: " // trim(dset_name)
   end if

   ! open dataset
   call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

   ! set as dimension scale
   if (present(dscale_name)) then
      call h5dsset_scale_f(dset_id, hdferror, dscale_name)
   else
      call h5dsset_scale_f(dset_id, hdferror, '')
   end if

   ! close id's
   call h5dclose_f(dset_id, hdferror)

end subroutine hdf_set_data_scale

subroutine hdf_attach_data_scale(dscale_loc_id, dscale_name, dset_loc_id, dset_name, dim_idx)
   !>   Attach a dimension scale to a dataset
   integer(HID_T), intent(in) :: dscale_loc_id        !> local id in file
   character(len=*), intent(in) :: dscale_name  !> name of dimension scale
   integer(HID_T), intent(in) :: dset_loc_id        !> local id in file
   character(len=*), intent(in) :: dset_name   !> name of dataset
   integer, intent(in), optional :: dim_idx  !> index of dimension to attach scale to (default: 1)

   integer(HID_T) :: dset_id, dscale_id
   integer :: idx
   integer :: hdferror

   if (.not. present(dim_idx)) then
      idx = 0
   else
      idx = dim_idx
   end if

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf5_attach_data_scale: " // trim(dset_name) // " to " // trim(dscale_name)
   end if

   ! open dataset
   call h5dopen_f(dset_loc_id, dset_name, dset_id, hdferror)

   ! open dimension scale
   call h5dopen_f(dscale_loc_id, dscale_name, dscale_id, hdferror)

   ! attach scale to dataset
   call h5dsattach_scale_f(dset_id, dscale_id, idx, hdferror)

   ! close id's
   call h5dclose_f(dscale_id, hdferror)
   call h5dclose_f(dset_id, hdferror)

end subroutine hdf_attach_data_scale



 !----------------------------------------------------------------------------------------
 !--------------------------------hdf_write_vector_to_dataset-----------------------------
 !----------------------------------------------------------------------------------------



subroutine hdf_create_dataset(loc_id, dset_name, dset_dims, dset_type)
   !>   create a dataset 'dset_name' with shape 'dset_dims' of type 'dset_type'
   integer(HID_T), intent(in) :: loc_id        !> local id in file
   character(len=*), intent(in) :: dset_name   !> name of dataset
   integer, intent(in) :: dset_dims(:)         !> dimensions of the dataset
   character(len=*), intent(in) :: dset_type   !> type of dataset (integer or double)

   integer :: rank
   integer(SIZE_T) :: dims(8)
   integer(HID_T) :: dset_id, dspace_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_create_dataset: " // trim(dset_name)
   end if

   ! set rank and dims
   rank = size(dset_dims, 1)
   dims(1:rank) = int(dset_dims, SIZE_T)

   ! create dataspace
   call h5screate_simple_f(rank, dims, dspace_id, hdferror)
   !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

   ! create dataset
   select case (dset_type)
    case('integer')
      call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_INTEGER, dspace_id, dset_id, hdferror)
      call h5dclose_f(dset_id, hdferror)
    case('double')
      call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferror)
      call h5dclose_f(dset_id, hdferror)
    case default
      write(*,'(A,A,A)') "---> ERROR: dset_type ", dset_type," not supported"
   end select

   ! close all id's
   call h5sclose_f(dspace_id, hdferror)

end subroutine hdf_create_dataset

 !
subroutine hdf_write_vector_to_dataset_double(loc_id, dset_name, offset, vector)

   integer(HID_T), intent(in) :: loc_id        !> local id in file
   character(len=*), intent(in) :: dset_name   !> name of dataset
   integer, intent(in) :: offset(:)            !> position within dataset
   real(dp), intent(in) :: vector(:)           !> data to be written

   integer(HID_T) :: dset_id, dspace_id, mspace_id
   integer :: rank
   integer(HSIZE_T) :: dset_dims(6), max_dims(6), mdims(1)
   integer(HSIZE_T) :: hs_count(6), hs_offset(6)
   integer :: hdferror

   integer :: i
   character(len=32) :: format_string

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_write_vector_to_dataset_double: " // trim(dset_name)
   end if

   ! open dataset
   call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

   ! get dataspace
   call h5dget_space_f(dset_id, dspace_id, hdferror)

   ! get rank (ndims), and dims
   call h5sget_simple_extent_ndims_f(dspace_id, rank, hdferror)
   call h5sget_simple_extent_dims_f(dspace_id, dset_dims(1:rank), max_dims(1:rank), hdferror)

   ! check size and offset
   if (size(vector,1) == dset_dims(1)) then

      if (all(offset(1:rank-1) <= dset_dims(2:rank)) .and. all(offset(1:rank-1) > 0)) then

         ! select hyperslab in dataset (note: convert FORTRAN offset to C offset by subtracting 1)
         hs_count(1) = dset_dims(1)
         hs_count(2:rank) = 1
         hs_offset(1) = 0
         hs_offset(2:rank) = offset(1:rank-1) - 1
         call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, hs_offset(1:rank), hs_count(1:rank), hdferror)

         ! set mspace to a vector
         mdims(1) = size(vector,1)
         call h5screate_simple_f(1, mdims, mspace_id, hdferror)

         ! write out vector
         call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, vector, mdims, hdferror, mspace_id, dspace_id)

         ! close mspace_id
         call h5sclose_f(mspace_id, hdferror)

      else
         write(format_string, '(A,I0,A,I0,A)') '(A,', rank-1, '(I0,A),A,', rank-1, '(I0,A),A)'
         write(*,format_string) "--->ERROR: offset=(", (offset(i), ',', i=1,rank-1) , &
            "), is not constent with dset_dims(2:rank)=(", (dset_dims(i), ',', i=2,rank),")"
      end if

   else
      write(*,'(A,I0,A,I0)') "--->ERROR: size(vector)=", size(vector), &
         ", is not constent with dset_dims(1)=", dset_dims(1)
   endif

   ! close id's
   call h5sclose_f(dspace_id, hdferror)
   call h5dclose_f(dset_id, hdferror)

end subroutine hdf_write_vector_to_dataset_double

 !
subroutine hdf_write_vector_to_dataset_integer(loc_id, dset_name, offset, vector)

   integer(HID_T), intent(in) :: loc_id        !> local id in file
   character(len=*), intent(in) :: dset_name   !> name of dataset
   integer, intent(in) :: offset(:)            !> position within dataset
   integer, intent(out) :: vector(:)           !> data to be written

   integer(HID_T) :: dset_id, dspace_id, mspace_id
   integer :: rank
   integer(HSIZE_T) :: dset_dims(6), max_dims(6), mdims(1)
   integer(HSIZE_T) :: hs_count(6), hs_offset(6)  ! Hyperslab offset
   integer :: hdferror

   integer :: i
   character(len=32) :: format_string

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_write_vector_to_dataset_integer: " // trim(dset_name)
   end if

   ! open dataset
   call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

   ! get dataspace
   call h5dget_space_f(dset_id, dspace_id, hdferror)

   ! get rank (ndims), and dims
   call h5sget_simple_extent_ndims_f(dspace_id, rank, hdferror)
   call h5sget_simple_extent_dims_f(dspace_id, dset_dims(1:rank), max_dims(1:rank), hdferror)

   ! check size and offset
   if (size(vector,1) == dset_dims(1)) then

      if (all(offset(1:rank-1) <= dset_dims(2:rank)) .and. all(offset(1:rank-1) > 0)) then

         ! select hyperslab in dataset (note: convert FORTRAN offset to C offset by subtracting 1)
         hs_count(1) = dset_dims(1)
         hs_count(2:rank) = 1
         hs_offset(1) = 0
         hs_offset(2:rank) = offset(1:rank-1) - 1
         call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, hs_offset(1:rank), hs_count(1:rank), hdferror)

         ! set mspace to a vector
         mdims(1) = size(vector,1)
         call h5screate_simple_f(1, mdims, mspace_id, hdferror)

         ! write out vector
         call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, vector, mdims, hdferror, mspace_id, dspace_id)

         ! close mspace_id
         call h5sclose_f(mspace_id, hdferror)

      else
         write(format_string, '(A,I0,A,I0,A)') '(A,', rank-1, '(I0,A),A,', rank-1, '(I0,A),A)'
         write(*,format_string) "--->ERROR: offset=(", (offset(i), ',', i=1,rank-1) , &
            "), is not constent with dset_dims(2:rank)=(", (dset_dims(i), ',', i=2,rank),")"
      end if

   else
      write(*,'(A,I0,A,I0)') "--->ERROR: size(vector)=", size(vector), &
         ", is not constent with dset_dims(1)=", dset_dims(1)
   endif

   ! close id's
   call h5sclose_f(dspace_id, hdferror)
   call h5dclose_f(dset_id, hdferror)

end subroutine hdf_write_vector_to_dataset_integer

 !
subroutine hdf_read_vector_from_dataset_double(loc_id, dset_name, offset, vector)

   integer(HID_T), intent(in) :: loc_id        !> local id in file
   character(len=*), intent(in) :: dset_name   !> name of dataset
   integer, intent(in) :: offset(:)            !> position within dataset
   real(dp), intent(out) :: vector(:)          !> data to be written

   integer(HID_T) :: dset_id, dspace_id, mspace_id
   integer :: rank
   integer(HSIZE_T) :: dset_dims(6), max_dims(6), mdims(1)
   integer(HSIZE_T) :: hs_count(6), hs_offset(6)
   integer :: hdferror

   integer :: i
   character(len=32) :: format_string

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_read_vector_from_dataset_double: " // trim(dset_name)
   end if

   ! open dataset
   call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

   ! get dataspace
   call h5dget_space_f(dset_id, dspace_id, hdferror)

   ! get rank (ndims), and dims
   call h5sget_simple_extent_ndims_f(dspace_id, rank, hdferror)
   call h5sget_simple_extent_dims_f(dspace_id, dset_dims(1:rank), max_dims(1:rank), hdferror)

   ! check size and offset
   if (size(vector,1) == dset_dims(1)) then

      if (all(offset(1:rank-1) <= dset_dims(2:rank)) .and. all(offset(1:rank-1) > 0)) then

         ! select hyperslab in dataset (note: convert FORTRAN offset to C offset by subtracting 1)
         hs_count(1) = dset_dims(1)
         hs_count(2:rank) = 1
         hs_offset(1) = 0
         hs_offset(2:rank) = offset(1:rank-1) - 1
         call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, hs_offset(1:rank), hs_count(1:rank), hdferror)

         ! set mspace to a vector
         mdims(1) = size(vector,1)
         call h5screate_simple_f(1, mdims, mspace_id, hdferror)

         ! write out vector
         call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, vector, mdims, hdferror, mspace_id, dspace_id)

         ! close mspace_id
         call h5sclose_f(mspace_id, hdferror)

      else
         write(format_string, '(A,I0,A,I0,A)') '(A,', rank-1, '(I0,A),A,', rank-1, '(I0,A),A)'
         write(*,format_string) "--->ERROR: offset=(", (offset(i), ',', i=1,rank-1) , &
            "), is not constent with dset_dims(2:rank)=(", (dset_dims(i), ',', i=2,rank),")"
      end if

   else
      write(*,'(A,I0,A,I0)') "--->ERROR: size(vector)=", size(vector), &
         ", is not constent with dset_dims(1)=", dset_dims(1)
   endif

   ! close id's
   call h5sclose_f(dspace_id, hdferror)
   call h5dclose_f(dset_id, hdferror)

end subroutine hdf_read_vector_from_dataset_double

 !
subroutine hdf_read_vector_from_dataset_integer(loc_id, dset_name, offset, vector)

   integer(HID_T), intent(in) :: loc_id        !> local id in file
   character(len=*), intent(in) :: dset_name   !> name of dataset
   integer, intent(in) :: offset(:)            !> position within dataset
   integer, intent(out) :: vector(:)           !> data to be written

   integer(HID_T) :: dset_id, dspace_id, mspace_id
   integer :: rank
   integer(HSIZE_T) :: dset_dims(6), max_dims(6), mdims(1)
   integer(HSIZE_T) :: hs_count(6), hs_offset(6)  ! Hyperslab offset
   integer :: hdferror

   integer :: i
   character(len=32) :: format_string

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_read_vector_from_dataset_integer: " // trim(dset_name)
   end if

   ! open dataset
   call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

   ! get dataspace
   call h5dget_space_f(dset_id, dspace_id, hdferror)

   ! get rank (ndims), and dims
   call h5sget_simple_extent_ndims_f(dspace_id, rank, hdferror)
   call h5sget_simple_extent_dims_f(dspace_id, dset_dims(1:rank), max_dims(1:rank), hdferror)

   ! check size and offset
   if (size(vector,1) == dset_dims(1)) then

      if (all(offset(1:rank-1) <= dset_dims(2:rank)) .and. all(offset(1:rank-1) > 0)) then

         ! select hyperslab in dataset (note: convert FORTRAN offset to C offset by subtracting 1)
         hs_count(1) = dset_dims(1)
         hs_count(2:rank) = 1
         hs_offset(1) = 0
         hs_offset(2:rank) = offset(1:rank-1) - 1
         call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, hs_offset(1:rank), hs_count(1:rank), hdferror)

         ! set mspace to a vector
         mdims(1) = size(vector,1)
         call h5screate_simple_f(1, mdims, mspace_id, hdferror)

         ! write out vector
         call h5dread_f(dset_id, H5T_NATIVE_INTEGER, vector, mdims, hdferror, mspace_id, dspace_id)

         ! close mspace_id
         call h5sclose_f(mspace_id, hdferror)

      else
         write(format_string, '(A,I0,A,I0,A)') '(A,', rank-1, '(I0,A),A,', rank-1, '(I0,A),A)'
         write(*,format_string) "--->ERROR: offset=(", (offset(i), ',', i=1,rank-1) , &
            "), is not constent with dset_dims(2:rank)=(", (dset_dims(i), ',', i=2,rank),")"
      end if

   else
      write(*,'(A,I0,A,I0)') "--->ERROR: size(vector)=", size(vector), &
         ", is not constent with dset_dims(1)=", dset_dims(1)
   endif

   ! close id's
   call h5sclose_f(dspace_id, hdferror)
   call h5dclose_f(dset_id, hdferror)

end subroutine hdf_read_vector_from_dataset_integer

 !----------------------------------------------------------------------------------------
 !--------------------------------hdf_write_dataset_double--------------------------------
 !----------------------------------------------------------------------------------------

subroutine hdf_write_dataset_double_0(loc_id, dset_name, data)
   !> writes a scalar to an hdf5 file
   integer(HID_T), intent(in) :: loc_id        !> local id in file
   character(len=*), intent(in) :: dset_name   !> name of dataset
   real(dp), intent(in) :: data                !> data to be written

   integer(SIZE_T) :: dims(1)
   integer(HID_T) :: dset_id, dspace_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_write_dataset_double_0: " // trim(dset_name)
   end if

   ! set rank and dims
   dims = (/ 0 /)

   ! create dataspace
   call h5screate_f(H5S_SCALAR_F, dspace_id, hdferror)
   !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

   ! create dataset
   call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dcreate: ", hdferror

   ! write dataset
   call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferror)
   !write(*,'(A20,I0)') "h5dwrite: ", hdferror

   ! close all id's
   call h5dclose_f(dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dclose: ", hdferror
   call h5sclose_f(dspace_id, hdferror)
   !write(*,'(A20,I0)') "h5sclose: ", hdferror

end subroutine hdf_write_dataset_double_0


subroutine hdf_write_dataset_double_1(loc_id, dset_name, data)
   !> writes a 1d array to an hdf5 file
   integer(HID_T), intent(in) :: loc_id        !> local id in file
   character(len=*), intent(in) :: dset_name   !> name of dataset
   real(dp), intent(in) :: data(:)             !> data to be written

   integer :: rank
   integer(SIZE_T) :: dims(1)
   integer(HID_T) :: dset_id, dspace_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_write_dataset_double_1: " // trim(dset_name)
   end if

   ! set rank and dims
   rank = 1
   dims = shape(data, KIND=HID_T)

   ! create dataspace
   call h5screate_simple_f(rank, dims, dspace_id, hdferror)
   !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

   ! create dataset
   call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dcreate: ", hdferror

   ! write dataset
   call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferror)
   !write(*,'(A20,I0)') "h5dwrite: ", hdferror

   ! close all id's
   call h5dclose_f(dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dclose: ", hdferror
   call h5sclose_f(dspace_id, hdferror)
   !write(*,'(A20,I0)') "h5sclose: ", hdferror

end subroutine hdf_write_dataset_double_1


subroutine hdf_write_dataset_double_2(loc_id, dset_name, data)
   !> writes a 2d array to an hdf5 file
   integer(HID_T), intent(in) :: loc_id        !> local id in file
   character(len=*), intent(in) :: dset_name   !> name of dataset
   real(dp), intent(in) :: data(:,:)           !> data to be written

   integer :: rank
   integer(SIZE_T) :: dims(2)
   integer(HID_T) :: dset_id, dspace_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_write_dataset_double_2: " // trim(dset_name)
   end if

   ! set rank and dims
   rank = 2
   dims = shape(data, KIND=HID_T)

   ! create dataspace
   call h5screate_simple_f(rank, dims, dspace_id, hdferror)
   !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

   ! create dataset
   call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dcreate: ", hdferror

   ! write dataset
   call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferror)
   !write(*,'(A20,I0)') "h5dwrite: ", hdferror

   ! close all id's
   call h5dclose_f(dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dclose: ", hdferror
   call h5sclose_f(dspace_id, hdferror)
   !write(*,'(A20,I0)') "h5sclose: ", hdferror

end subroutine hdf_write_dataset_double_2

subroutine hdf_write_dataset_double_3(loc_id, dset_name, data)
   !> writes a 3d array to an hdf5 file
   integer(HID_T), intent(in) :: loc_id        !> local id in file
   character(len=*), intent(in) :: dset_name   !> name of dataset
   real(dp), intent(in) :: data(:,:,:)         !> data to be written

   integer :: rank
   integer(SIZE_T) :: dims(3)
   integer(HID_T) :: dset_id, dspace_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_write_dataset_double_3: " // trim(dset_name)
   end if

   ! set rank and dims
   rank = 3
   dims = shape(data, KIND=HID_T)

   ! create dataspace
   call h5screate_simple_f(rank, dims, dspace_id, hdferror)
   !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

   ! create dataset
   call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dcreate: ", hdferror

   ! write dataset
   call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferror)
   !write(*,'(A20,I0)') "h5dwrite: ", hdferror

   ! close all id's
   call h5dclose_f(dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dclose: ", hdferror
   call h5sclose_f(dspace_id, hdferror)
   !write(*,'(A20,I0)') "h5sclose: ", hdferror

end subroutine hdf_write_dataset_double_3


subroutine hdf_write_dataset_double_4(loc_id, dset_name, data)
   !>   writes a 4d array to an hdf5 file
   integer(HID_T), intent(in) :: loc_id        !> local id in file
   character(len=*), intent(in) :: dset_name   !> name of dataset
   real(dp), intent(in) :: data(:,:,:,:)       !> data to be written

   integer :: rank
   integer(SIZE_T) :: dims(4)
   integer(HID_T) :: dset_id, dspace_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_write_dataset_double_4: " // trim(dset_name)
   end if

   ! set rank and dims
   rank = 4
   dims = shape(data, KIND=HID_T)

   ! create dataspace
   call h5screate_simple_f(rank, dims, dspace_id, hdferror)
   !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

   ! create dataset
   call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dcreate: ", hdferror

   ! write dataset
   call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferror)
   !write(*,'(A20,I0)') "h5dwrite: ", hdferror

   ! close all id's
   call h5dclose_f(dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dclose: ", hdferror
   call h5sclose_f(dspace_id, hdferror)
   !write(*,'(A20,I0)') "h5sclose: ", hdferror

end subroutine hdf_write_dataset_double_4


subroutine hdf_write_dataset_double_5(loc_id, dset_name, data)
   !>  writes a 5d array to an hdf5 file
   integer(HID_T), intent(in) :: loc_id        !> local id in file
   character(len=*), intent(in) :: dset_name   !> name of dataset
   real(dp), intent(in) :: data(:,:,:,:,:)     !> data to be written

   integer :: rank
   integer(SIZE_T) :: dims(5)
   integer(HID_T) :: dset_id, dspace_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_write_dataset_double_5: " // trim(dset_name)
   end if

   ! set rank and dims
   rank = 5
   dims = shape(data, KIND=HID_T)

   ! create dataspace
   call h5screate_simple_f(rank, dims, dspace_id, hdferror)
   !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

   ! create dataset
   call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dcreate: ", hdferror

   ! write dataset
   call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferror)
   !write(*,'(A20,I0)') "h5dwrite: ", hdferror

   ! close all id's
   call h5dclose_f(dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dclose: ", hdferror
   call h5sclose_f(dspace_id, hdferror)
   !write(*,'(A20,I0)') "h5sclose: ", hdferror

end subroutine hdf_write_dataset_double_5


subroutine hdf_write_dataset_double_6(loc_id, dset_name, data)
   !>   writes a 6d array to an hdf5 file
   integer(HID_T), intent(in) :: loc_id        !> local id in file
   character(len=*), intent(in) :: dset_name   !> name of dataset
   real(dp), intent(in) :: data(:,:,:,:,:,:)   !> data to be written

   integer :: rank
   integer(SIZE_T) :: dims(6)
   integer(HID_T) :: dset_id, dspace_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_write_dataset_double_6: " // trim(dset_name)
   end if

   ! set rank and dims
   rank = 6
   dims = shape(data, KIND=HID_T)

   ! create dataspace
   call h5screate_simple_f(rank, dims, dspace_id, hdferror)
   !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

   ! create dataset
   call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dcreate: ", hdferror

   ! write dataset
   call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferror)
   !write(*,'(A20,I0)') "h5dwrite: ", hdferror

   ! close all id's
   call h5dclose_f(dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dclose: ", hdferror
   call h5sclose_f(dspace_id, hdferror)
   !write(*,'(A20,I0)') "h5sclose: ", hdferror

end subroutine hdf_write_dataset_double_6


 !----------------------------------------------------------------------------------------
 !--------------------------------hdf_write_dataset_double--------------------------------
 !----------------------------------------------------------------------------------------

subroutine hdf_write_dataset_integer_0(loc_id, dset_name, data)
   !>   writes a scalar to an hdf5 file
   integer(HID_T), intent(in) :: loc_id        !> local id in file
   character(len=*), intent(in) :: dset_name   !> name of dataset
   integer, intent(in) :: data                 !> data to be written

   integer(SIZE_T) :: dims(1)
   integer(HID_T) :: dset_id, dspace_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_write_dataset_integer_0: " // trim(dset_name)
   end if

   ! set rank and dims
   dims = (/ 0 /)

   ! create dataspace
   call h5screate_f(H5S_SCALAR_F, dspace_id, hdferror)
   !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

   ! create dataset
   call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_INTEGER, dspace_id, dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dcreate: ", hdferror

   ! write dataset
   call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data, dims, hdferror)
   !write(*,'(A20,I0)') "h5dwrite: ", hdferror

   ! close all id's
   call h5dclose_f(dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dclose: ", hdferror
   call h5sclose_f(dspace_id, hdferror)
   !write(*,'(A20,I0)') "h5sclose: ", hdferror

end subroutine hdf_write_dataset_integer_0

subroutine hdf_write_dataset_integer_1(loc_id, dset_name, data)
   !>   writes a 1d array to an hdf5 file
   integer(HID_T), intent(in) :: loc_id        !> local id in file
   character(len=*), intent(in) :: dset_name   !> name of dataset
   integer, intent(in) :: data(:)              !> data to be written

   integer :: rank
   integer(SIZE_T) :: dims(1)
   integer(HID_T) :: dset_id, dspace_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_write_dataset_integer_1: " // trim(dset_name)
   end if

   ! set rank and dims
   rank = 1
   dims = shape(data, KIND=HID_T)

   ! create dataspace
   call h5screate_simple_f(rank, dims, dspace_id, hdferror)
   !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

   ! create dataset
   call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_INTEGER, dspace_id, dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dcreate: ", hdferror

   ! write dataset
   call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data, dims, hdferror)
   !write(*,'(A20,I0)') "h5dwrite: ", hdferror

   ! close all id's
   call h5dclose_f(dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dclose: ", hdferror
   call h5sclose_f(dspace_id, hdferror)
   !write(*,'(A20,I0)') "h5sclose: ", hdferror

end subroutine hdf_write_dataset_integer_1

subroutine hdf_write_dataset_integer_2(loc_id, dset_name, data)
   !>   writes a 2d array to an hdf5 file
   integer(HID_T), intent(in) :: loc_id        !> local id in file
   character(len=*), intent(in) :: dset_name   !> name of dataset
   integer, intent(in) :: data(:,:)            !> data to be written

   integer :: rank
   integer(SIZE_T) :: dims(2)
   integer(HID_T) :: dset_id, dspace_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_write_dataset_integer_2: " // trim(dset_name)
   end if

   ! set rank and dims
   rank = 2
   dims = shape(data, KIND=HID_T)

   ! create dataspace
   call h5screate_simple_f(rank, dims, dspace_id, hdferror)
   !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

   ! create dataset
   call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_INTEGER, dspace_id, dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dcreate: ", hdferror

   ! write dataset
   call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data, dims, hdferror)
   !write(*,'(A20,I0)') "h5dwrite: ", hdferror

   ! close all id's
   call h5dclose_f(dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dclose: ", hdferror
   call h5sclose_f(dspace_id, hdferror)
   !write(*,'(A20,I0)') "h5sclose: ", hdferror

end subroutine hdf_write_dataset_integer_2

subroutine hdf_write_dataset_integer_3(loc_id, dset_name, data)
   !>   writes a 3d array to an hdf5 file
   integer(HID_T), intent(in) :: loc_id        !> local id in file
   character(len=*), intent(in) :: dset_name   !> name of dataset
   integer, intent(in) :: data(:,:,:)          !> data to be written

   integer :: rank
   integer(SIZE_T) :: dims(3)
   integer(HID_T) :: dset_id, dspace_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_write_dataset_integer_3: " // trim(dset_name)
   end if

   ! set rank and dims
   rank = 3
   dims = shape(data, KIND=HID_T)

   ! create dataspace
   call h5screate_simple_f(rank, dims, dspace_id, hdferror)
   !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

   ! create dataset
   call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_INTEGER, dspace_id, dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dcreate: ", hdferror

   ! write dataset
   call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data, dims, hdferror)
   !write(*,'(A20,I0)') "h5dwrite: ", hdferror

   ! close all id's
   call h5dclose_f(dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dclose: ", hdferror
   call h5sclose_f(dspace_id, hdferror)
   !write(*,'(A20,I0)') "h5sclose: ", hdferror

end subroutine hdf_write_dataset_integer_3

subroutine hdf_write_dataset_integer_4(loc_id, dset_name, data)
   !>   writes a 4d array to an hdf5 file
   integer(HID_T), intent(in) :: loc_id        !> local id in file
   character(len=*), intent(in) :: dset_name   !> name of dataset
   integer, intent(in) :: data(:,:,:,:)        !> data to be written

   integer :: rank
   integer(SIZE_T) :: dims(4)
   integer(HID_T) :: dset_id, dspace_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_write_dataset_integer_4: " // trim(dset_name)
   end if

   ! set rank and dims
   rank = 4
   dims = shape(data, KIND=HID_T)

   ! create dataspace
   call h5screate_simple_f(rank, dims, dspace_id, hdferror)
   !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

   ! create dataset
   call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_INTEGER, dspace_id, dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dcreate: ", hdferror

   ! write dataset
   call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data, dims, hdferror)
   !write(*,'(A20,I0)') "h5dwrite: ", hdferror

   ! close all id's
   call h5dclose_f(dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dclose: ", hdferror
   call h5sclose_f(dspace_id, hdferror)
   !write(*,'(A20,I0)') "h5sclose: ", hdferror

end subroutine hdf_write_dataset_integer_4

subroutine hdf_write_dataset_integer_5(loc_id, dset_name, data)
   !>   writes a 5d array to an hdf5 file
   integer(HID_T), intent(in) :: loc_id        !> local id in file
   character(len=*), intent(in) :: dset_name   !> name of dataset
   integer, intent(in) :: data(:,:,:,:,:)      !> data to be written

   integer :: rank
   integer(SIZE_T) :: dims(5)
   integer(HID_T) :: dset_id, dspace_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_write_dataset_integer_5: " // trim(dset_name)
   end if

   ! set rank and dims
   rank = 5
   dims = shape(data, KIND=HID_T)

   ! create dataspace
   call h5screate_simple_f(rank, dims, dspace_id, hdferror)
   !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

   ! create dataset
   call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_INTEGER, dspace_id, dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dcreate: ", hdferror

   ! write dataset
   call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data, dims, hdferror)
   !write(*,'(A20,I0)') "h5dwrite: ", hdferror

   ! close all id's
   call h5dclose_f(dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dclose: ", hdferror
   call h5sclose_f(dspace_id, hdferror)
   !write(*,'(A20,I0)') "h5sclose: ", hdferror

end subroutine hdf_write_dataset_integer_5

subroutine hdf_write_dataset_integer_6(loc_id, dset_name, data)
   !>   writes a 6d array to an hdf5 file
   integer(HID_T), intent(in) :: loc_id        !> local id in file
   character(len=*), intent(in) :: dset_name   !> name of dataset
   integer, intent(in) :: data(:,:,:,:,:,:)    !> data to be written

   integer :: rank
   integer(SIZE_T) :: dims(6)
   integer(HID_T) :: dset_id, dspace_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_write_dataset_integer_6: " // trim(dset_name)
   end if

   ! set rank and dims
   rank = 6
   dims = shape(data, KIND=HID_T)

   ! create dataspace
   call h5screate_simple_f(rank, dims, dspace_id, hdferror)
   !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

   ! create dataset
   call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_INTEGER, dspace_id, dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dcreate: ", hdferror

   ! write dataset
   call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data, dims, hdferror)
   !write(*,'(A20,I0)') "h5dwrite: ", hdferror

   ! close all id's
   call h5dclose_f(dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dclose: ", hdferror
   call h5sclose_f(dspace_id, hdferror)
   !write(*,'(A20,I0)') "h5sclose: ", hdferror

end subroutine hdf_write_dataset_integer_6


 !---------------------------------------------------------------------------------------
 !--------------------------------hdf_read_dataset_integer--------------------------------
 !---------------------------------------------------------------------------------------



subroutine hdf_read_dataset_integer_0(loc_id, dset_name, data)
   !>   reads a scalar from an hdf5 file
   integer(HID_T), intent(in) :: loc_id        !> local id in file
   character(len=*), intent(in) :: dset_name   !> name of dataset
   integer, intent(out) :: data                !> data to be written


   integer(SIZE_T) :: dims(1)
   integer(HID_T) :: dset_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_read_dataset_integer_0: " // trim(dset_name)
   end if

   ! set rank and dims
   dims = (/ 0 /)

   ! open dataset
   call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dcreate: ", hdferror

   ! read dataset
   call h5dread_f(dset_id, H5T_NATIVE_INTEGER, data, dims, hdferror)
   !write(*,'(A20,I0)') "h5dread: ", hdferror

   ! close all id's
   call h5dclose_f(dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dclose: ", hdferror

end subroutine hdf_read_dataset_integer_0

subroutine hdf_read_dataset_integer_1(loc_id, dset_name, data)
   !>   reads a 1d array from an hdf5 file
   integer(HID_T), intent(in) :: loc_id       !> local id in file
   character(len=*), intent(in) :: dset_name  !> name of dataset
   integer, intent(out) :: data(:)            !> data to be written

   integer :: rank
   integer(SIZE_T) :: dims(1)
   integer(HID_T) :: dset_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_read_dataset_integer_1: " // trim(dset_name)
   end if

   ! set rank and dims
   rank = 1
   dims = shape(data, KIND=HID_T)

   ! open dataset
   call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dcreate: ", hdferror

   ! read dataset
   call h5dread_f(dset_id, H5T_NATIVE_INTEGER, data, dims, hdferror)
   !write(*,'(A20,I0)') "h5dread: ", hdferror

   ! close all id's
   call h5dclose_f(dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dclose: ", hdferror

end subroutine hdf_read_dataset_integer_1


subroutine hdf_read_dataset_integer_2(loc_id, dset_name, data)
   !>   reads a 2d array from an hdf5 file
   integer(HID_T), intent(in) :: loc_id       !> local id in file
   character(len=*), intent(in) :: dset_name  !> name of dataset
   integer, intent(out) :: data(:,:)          !> data to be written

   integer :: rank
   integer(SIZE_T) :: dims(2)
   integer(HID_T) :: dset_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_read_dataset_integer_2: " // trim(dset_name)
   end if

   ! set rank and dims
   rank = 2
   dims = shape(data, KIND=HID_T)

   ! open dataset
   call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dcreate: ", hdferror

   ! read dataset
   call h5dread_f(dset_id, H5T_NATIVE_INTEGER, data, dims, hdferror)
   !write(*,'(A20,I0)') "h5dread: ", hdferror

   ! close all id's
   call h5dclose_f(dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dclose: ", hdferror

end subroutine hdf_read_dataset_integer_2


subroutine hdf_read_dataset_integer_3(loc_id, dset_name, data)
   !>   reads a 3d array from an hdf5 file
   integer(HID_T), intent(in) :: loc_id       !> local id in file
   character(len=*), intent(in) :: dset_name  !> name of dataset
   integer, intent(out) :: data(:,:,:)        !> data to be written

   integer :: rank
   integer(SIZE_T) :: dims(3)
   integer(HID_T) :: dset_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_read_dataset_integer_3: " // trim(dset_name)
   end if

   ! set rank and dims
   rank = 3
   dims = shape(data, KIND=HID_T)

   ! open dataset
   call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dcreate: ", hdferror

   ! read dataset
   call h5dread_f(dset_id, H5T_NATIVE_INTEGER, data, dims, hdferror)
   !write(*,'(A20,I0)') "h5dread: ", hdferror

   ! close all id's
   call h5dclose_f(dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dclose: ", hdferror

end subroutine hdf_read_dataset_integer_3


subroutine hdf_read_dataset_integer_4(loc_id, dset_name, data)
   !>   reads a 4d array from an hdf5 file
   integer(HID_T), intent(in) :: loc_id       !> local id in file
   character(len=*), intent(in) :: dset_name  !> name of dataset
   integer, intent(out) :: data(:,:,:,:)      !> data to be written

   integer :: rank
   integer(SIZE_T) :: dims(4)
   integer(HID_T) :: dset_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_read_dataset_integer_4: " // trim(dset_name)
   end if

   ! set rank and dims
   rank = 4
   dims = shape(data, KIND=HID_T)

   ! open dataset
   call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dcreate: ", hdferror

   ! read dataset
   call h5dread_f(dset_id, H5T_NATIVE_INTEGER, data, dims, hdferror)
   !write(*,'(A20,I0)') "h5dread: ", hdferror

   ! close all id's
   call h5dclose_f(dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dclose: ", hdferror

end subroutine hdf_read_dataset_integer_4


subroutine hdf_read_dataset_integer_5(loc_id, dset_name, data)
   !>   reads a 5d array from an hdf5 file
   integer(HID_T), intent(in) :: loc_id       !> local id in file
   character(len=*), intent(in) :: dset_name  !> name of dataset
   integer, intent(out) :: data(:,:,:,:,:)    !> data to be written

   integer :: rank
   integer(SIZE_T) :: dims(5)
   integer(HID_T) :: dset_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_read_dataset_integer_5: " // trim(dset_name)
   end if

   ! set rank and dims
   rank = 5
   dims = shape(data, KIND=HID_T)

   ! open dataset
   call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dcreate: ", hdferror

   ! read dataset
   call h5dread_f(dset_id, H5T_NATIVE_INTEGER, data, dims, hdferror)
   !write(*,'(A20,I0)') "h5dread: ", hdferror

   ! close all id's
   call h5dclose_f(dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dclose: ", hdferror

end subroutine hdf_read_dataset_integer_5


subroutine hdf_read_dataset_integer_6(loc_id, dset_name, data)
   !>   reads a 6d array from an hdf5 file
   integer(HID_T), intent(in) :: loc_id       !> local id in file
   character(len=*), intent(in) :: dset_name  !> name of dataset
   integer, intent(out) :: data(:,:,:,:,:,:)  !> data to be written

   integer :: rank
   integer(SIZE_T) :: dims(6)
   integer(HID_T) :: dset_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_read_dataset_integer_6: " // trim(dset_name)
   end if

   ! set rank and dims
   rank = 6
   dims = shape(data, KIND=HID_T)

   ! open dataset
   call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dcreate: ", hdferror

   ! read dataset
   call h5dread_f(dset_id, H5T_NATIVE_INTEGER, data, dims, hdferror)
   !write(*,'(A20,I0)') "h5dread: ", hdferror

   ! close all id's
   call h5dclose_f(dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dclose: ", hdferror

end subroutine hdf_read_dataset_integer_6


 !---------------------------------------------------------------------------------------
 !--------------------------------hdf_read_dataset_double--------------------------------
 !---------------------------------------------------------------------------------------



subroutine hdf_read_dataset_double_0(loc_id, dset_name, data)
   !>   reads a scalar from an hdf5 file
   integer(HID_T), intent(in) :: loc_id       !> local id in file
   character(len=*), intent(in) :: dset_name  !> name of dataset
   real(dp), intent(out) :: data              !> data to be written


   integer(SIZE_T) :: dims(1)
   integer(HID_T) :: dset_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_read_dataset_double_0: " // trim(dset_name)
   end if

   ! set rank and dims
   dims = (/ 0 /)

   ! open dataset
   call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dcreate: ", hdferror

   ! read dataset
   call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferror)
   !write(*,'(A20,I0)') "h5dread: ", hdferror

   ! close all id's
   call h5dclose_f(dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dclose: ", hdferror

end subroutine hdf_read_dataset_double_0


subroutine hdf_read_dataset_double_1(loc_id, dset_name, data)
   !>   reads a 1d array from an hdf5 file
   integer(HID_T), intent(in) :: loc_id       !> local id in file
   character(len=*), intent(in) :: dset_name  !> name of dataset
   real(dp), intent(out) :: data(:)           !> data to be written

   integer :: rank
   integer(SIZE_T) :: dims(1)
   integer(HID_T) :: dset_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_read_dataset_double_1: " // trim(dset_name)
   end if

   ! set rank and dims
   rank = 1
   dims = shape(data, KIND=HID_T)

   ! open dataset
   call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dcreate: ", hdferror

   ! read dataset
   call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferror)
   !write(*,'(A20,I0)') "h5dread: ", hdferror

   ! close all id's
   call h5dclose_f(dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dclose: ", hdferror

end subroutine hdf_read_dataset_double_1


subroutine hdf_read_dataset_double_2(loc_id, dset_name, data)
   !>   reads a 2d array from an hdf5 file
   integer(HID_T), intent(in) :: loc_id       !> local id in file
   character(len=*), intent(in) :: dset_name  !> name of dataset
   real(dp), intent(out) :: data(:,:)         !> data to be written

   integer :: rank
   integer(SIZE_T) :: dims(2)
   integer(HID_T) :: dset_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_read_dataset_double_1: " // trim(dset_name)
   end if

   ! set rank and dims
   rank = 2
   dims = shape(data, KIND=HID_T)

   ! open dataset
   call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dcreate: ", hdferror

   ! read dataset
   call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferror)
   !write(*,'(A20,I0)') "h5dread: ", hdferror

   ! close all id's
   call h5dclose_f(dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dclose: ", hdferror

end subroutine hdf_read_dataset_double_2


subroutine hdf_read_dataset_double_3(loc_id, dset_name, data)
   !>   reads a 3d array from an hdf5 file
   integer(HID_T), intent(in) :: loc_id       !> local id in file
   character(len=*), intent(in) :: dset_name  !> name of dataset
   real(dp), intent(out) :: data(:,:,:)       !> data to be written

   integer :: rank
   integer(SIZE_T) :: dims(3)
   integer(HID_T) :: dset_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_read_dataset_double_3: " // trim(dset_name)
   end if

   ! set rank and dims
   rank = 3
   dims = shape(data, KIND=HID_T)

   ! open dataset
   call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dcreate: ", hdferror

   ! read dataset
   call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferror)
   !write(*,'(A20,I0)') "h5dread: ", hdferror

   ! close all id's
   call h5dclose_f(dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dclose: ", hdferror

end subroutine hdf_read_dataset_double_3


subroutine hdf_read_dataset_double_4(loc_id, dset_name, data)
   !>   reads a 4d array from an hdf5 file
   integer(HID_T), intent(in) :: loc_id       !> local id in file
   character(len=*), intent(in) :: dset_name  !> name of dataset
   real(dp), intent(out) :: data(:,:,:,:)     !> data to be written

   integer :: rank
   integer(SIZE_T) :: dims(4)
   integer(HID_T) :: dset_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_read_dataset_double_4: " // trim(dset_name)
   end if

   ! set rank and dims
   rank = 4
   dims = shape(data, KIND=HID_T)

   ! open dataset
   call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dcreate: ", hdferror

   ! read dataset
   call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferror)
   !write(*,'(A20,I0)') "h5dread: ", hdferror

   ! close all id's
   call h5dclose_f(dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dclose: ", hdferror

end subroutine hdf_read_dataset_double_4


subroutine hdf_read_dataset_double_5(loc_id, dset_name, data)
   !>   reads a 5d array from an hdf5 file
   integer(HID_T), intent(in) :: loc_id       !> local id in file
   character(len=*), intent(in) :: dset_name  !> name of dataset
   real(dp), intent(out) :: data(:,:,:,:,:)   !> data to be written

   integer :: rank
   integer(SIZE_T) :: dims(5)
   integer(HID_T) :: dset_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_read_dataset_double_5: " // trim(dset_name)
   end if

   ! set rank and dims
   rank = 5
   dims = shape(data, KIND=HID_T)

   ! open dataset
   call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dcreate: ", hdferror

   ! read dataset
   call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferror)
   !write(*,'(A20,I0)') "h5dread: ", hdferror

   ! close all id's
   call h5dclose_f(dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dclose: ", hdferror

end subroutine hdf_read_dataset_double_5


subroutine hdf_read_dataset_double_6(loc_id, dset_name, data)
   !>   reads a 6d array from an hdf5 file
   integer(HID_T), intent(in) :: loc_id       !> local id in file
   character(len=*), intent(in) :: dset_name  !> name of dataset
   real(dp), intent(out) :: data(:,:,:,:,:,:)  !> data to be written

   integer :: rank
   integer(SIZE_T) :: dims(6)
   integer(HID_T) :: dset_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_read_dataset_double_6: " // trim(dset_name)
   end if

   ! set rank and dims
   rank = 6
   dims = shape(data, KIND=HID_T)

   ! open dataset
   call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dcreate: ", hdferror

   ! read dataset
   call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferror)
   !write(*,'(A20,I0)') "h5dread: ", hdferror

   ! close all id's
   call h5dclose_f(dset_id, hdferror)
   !write(*,'(A20,I0)') "h5dclose: ", hdferror

end subroutine hdf_read_dataset_double_6


 !---------------------------------------------------------------------------------------
 !---------------------------------------hdf_write_attr*---------------------------------
 !---------------------------------------------------------------------------------------



subroutine hdf_write_attr_double_0(loc_id, obj_name, attr_name, data)
   !>   writes a scalar attribute
   integer(HID_T), intent(in) :: loc_id       !> local id in file
   character(len=*), intent(in) :: obj_name   !> object name attribute will be attached to (if "" use loc_id)
   character(len=*), intent(in) :: attr_name  !> name of attribute
   real(dp), intent(in) :: data               !> data to write to attribute

   !integer :: rank
   integer(SIZE_T) :: dims(1)
   integer(HID_T) :: obj_id, aspace_id, attr_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_write_attr_double_0: " // trim(obj_name) // "/" //trim(attr_name)
   end if

   ! open object
   if (obj_name == "") then
      obj_id = loc_id
   else
      call h5oopen_f(loc_id, obj_name, obj_id, hdferror)
   end if

   ! create dataspace
   dims = (/ 0 /)
   call h5screate_f(H5S_SCALAR_F, aspace_id, hdferror)
   !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

   ! create attribute
   call h5acreate_f(obj_id, attr_name, H5T_NATIVE_DOUBLE, aspace_id, attr_id, hdferror)
   !write(*,'(A20,I0)') "h5dcreate: ", hdferror

   ! write dataset
   call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, data, dims, hdferror)
   !write(*,'(A20,I0)') "h5dwrite: ", hdferror

   ! close all id's
   call h5aclose_f(attr_id, hdferror)
   !write(*,'(A20,I0)') "h5dclose: ", hdferror
   call h5sclose_f(aspace_id, hdferror)
   !write(*,'(A20,I0)') "h5sclose: ", hdferror
   if (obj_name /= "") then
      call h5oclose_f(obj_id, hdferror)
   end if

end subroutine hdf_write_attr_double_0


subroutine hdf_write_attr_double_1(loc_id, obj_name, attr_name, data)
   !>   writes 1d array attribute
   integer(HID_T), intent(in) :: loc_id       !> local id in file
   character(len=*), intent(in) :: obj_name   !> object name attribute will be attached to (if "" use loc_id)
   character(len=*), intent(in) :: attr_name  !> name of attribute
   real(dp), intent(in) :: data(:)            !> data to write to attribute

   integer :: rank
   integer(SIZE_T) :: dims(1)
   integer(HID_T) :: obj_id, aspace_id, attr_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_write_attr_double_1: " // trim(obj_name) // "/" //trim(attr_name)
   end if

   ! open object
   if (obj_name == "") then
      obj_id = loc_id
   else
      call h5oopen_f(loc_id, obj_name, obj_id, hdferror)
   end if

   ! create dataspace
   rank = 1
   dims = shape(data, KIND=HID_T)
   call h5screate_simple_f(rank, dims, aspace_id, hdferror)
   !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

   ! create attribute
   call h5acreate_f(obj_id, attr_name, H5T_NATIVE_DOUBLE, aspace_id, attr_id, hdferror)
   !write(*,'(A20,I0)') "h5dcreate: ", hdferror

   ! write dataset
   call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, data, dims, hdferror)
   !write(*,'(A20,I0)') "h5dwrite: ", hdferror

   ! close all id's
   call h5aclose_f(attr_id, hdferror)
   !write(*,'(A20,I0)') "h5dclose: ", hdferror
   call h5sclose_f(aspace_id, hdferror)
   !write(*,'(A20,I0)') "h5sclose: ", hdferror
   if (obj_name /= "") then
      call h5oclose_f(obj_id, hdferror)
   end if

end subroutine hdf_write_attr_double_1


subroutine hdf_write_attr_integer_0(loc_id, obj_name, attr_name, data)
   !>   writes a scalar attribute
   integer(HID_T), intent(in) :: loc_id       !> local id in file
   character(len=*), intent(in) :: obj_name   !> object name attribute will be attached to (if "" use loc_id)
   character(len=*), intent(in) :: attr_name  !> name of attribute
   integer, intent(in) :: data                !> data to write to attribute

   integer(SIZE_T) :: dims(1)
   integer(HID_T) :: obj_id, aspace_id, attr_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_write_attr_integer_0: " // trim(obj_name) // "/" // trim(attr_name)
   end if

   ! open object
   if (obj_name == "") then
      obj_id = loc_id
   else
      call h5oopen_f(loc_id, obj_name, obj_id, hdferror)
   end if

   ! create dataspace
   dims = (/ 0 /)
   call h5screate_f(H5S_SCALAR_F, aspace_id, hdferror)
   !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

   ! create attribute
   call h5acreate_f(obj_id, attr_name, H5T_NATIVE_INTEGER, aspace_id, attr_id, hdferror)
   !write(*,'(A20,I0)') "h5dcreate: ", hdferror

   ! write dataset
   call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, data, dims, hdferror)
   !write(*,'(A20,I0)') "h5dwrite: ", hdferror

   ! close all id's
   call h5aclose_f(attr_id, hdferror)
   !write(*,'(A20,I0)') "h5dclose: ", hdferror
   call h5sclose_f(aspace_id, hdferror)
   !write(*,'(A20,I0)') "h5sclose: ", hdferror
   if (obj_name /= "") then
      call h5oclose_f(obj_id, hdferror)
   end if

end subroutine hdf_write_attr_integer_0


subroutine hdf_write_attr_integer_1(loc_id, obj_name, attr_name, data)
   !>   writes 1d array attribute
   integer(HID_T), intent(in) :: loc_id       !> local id in file
   character(len=*), intent(in) :: obj_name   !> object name attribute will be attached to (if "" use loc_id)
   character(len=*), intent(in) :: attr_name  !> name of attribute
   integer, intent(in) :: data(:)             !> data to write to attribute

   integer :: rank
   integer(SIZE_T) :: dims(1)
   integer(HID_T) :: obj_id, aspace_id, attr_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_write_attr_integer_1: " // trim(obj_name) // "/" //trim(attr_name)
   end if

   ! open object
   if (obj_name == "") then
      obj_id = loc_id
   else
      call h5oopen_f(loc_id, obj_name, obj_id, hdferror)
   end if

   ! create dataspace
   rank = 1
   dims = shape(data, KIND=HID_T)
   call h5screate_simple_f(rank, dims, aspace_id, hdferror)
   !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

   ! create attribute
   call h5acreate_f(obj_id, attr_name, H5T_NATIVE_INTEGER, aspace_id, attr_id, hdferror)
   !write(*,'(A20,I0)') "h5dcreate: ", hdferror

   ! write dataset
   call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, data, dims, hdferror)
   !write(*,'(A20,I0)') "h5dwrite: ", hdferror

   ! close all id's
   call h5aclose_f(attr_id, hdferror)
   !write(*,'(A20,I0)') "h5dclose: ", hdferror
   call h5sclose_f(aspace_id, hdferror)
   !write(*,'(A20,I0)') "h5sclose: ", hdferror
   if (obj_name /= "") then
      call h5oclose_f(obj_id, hdferror)
   end if

end subroutine hdf_write_attr_integer_1


subroutine hdf_write_attr_string(loc_id, obj_name, attr_name, data)
   !>   writes a string attribute
   integer(HID_T), intent(in) :: loc_id       !> local id in file
   character(len=*), intent(in) :: obj_name   !> object name attribute will be attached to (if "" use loc_id)
   character(len=*), intent(in) :: attr_name  !> name of attribute
   character(len=*), intent(in) :: data       !> data to write to attribute

   integer(SIZE_T) :: dims(1)
   integer(HID_T) :: obj_id, type_id, aspace_id, attr_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_write_attr_string: " // trim(obj_name) // "/" //trim(attr_name)
   end if

   ! open object
   if (obj_name == "") then
      obj_id = loc_id
   else
      call h5oopen_f(loc_id, obj_name, obj_id, hdferror)
   end if

   ! create type_id and aspace_id
   dims(1) = len(data, KIND=HID_T)
   call h5tcopy_f (H5T_NATIVE_CHARACTER, type_id, hdferror)
   !write(*,*) 'h5tcopy_f returns', type_id
   call h5tset_size_f (type_id, dims(1), hdferror)
   !write(*,*) 'h5tset_size_f returns', hdferror
   call h5screate_f (H5S_SCALAR_F, aspace_id, hdferror)

   ! create attribute
   call h5acreate_f(obj_id, attr_name, type_id, aspace_id, attr_id, hdferror)
   !write(*,'(A20,I0)') "h5dcreate: ", hdferror

   ! write dataset
   call h5awrite_f(attr_id, type_id, data, dims, hdferror)
   !write(*,'(A20,I0)') "h5dwrite: ", hdferror

   ! close all id's
   call h5tclose_f(type_id, hdferror)
   call h5aclose_f(attr_id, hdferror)
   !write(*,'(A20,I0)') "h5dclose: ", hdferror
   call h5sclose_f(aspace_id, hdferror)
   !write(*,'(A20,I0)') "h5sclose: ", hdferror
   if (obj_name /= "") then
      call h5oclose_f(obj_id, hdferror)
   end if

end subroutine hdf_write_attr_string


 !---------------------------------------------------------------------------------------
 !-------------------------------------hdf_read_attr-------------------------------------
 !---------------------------------------------------------------------------------------



subroutine hdf_read_attr_double_0(loc_id, obj_name, attr_name, data)
   !>   writes a scalar attribute
   integer(HID_T), intent(in) :: loc_id       !> local id in file
   character(len=*), intent(in) :: obj_name   !> object name attribute will be attached to (if "" use loc_id)
   character(len=*), intent(in) :: attr_name  !> name of attribute
   real(dp), intent(out) :: data              !> data to write to attribute

   integer(SIZE_T) :: dims(1)
   integer(HID_T) :: obj_id, attr_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_read_attr_double_0: " // trim(obj_name) // "/" //trim(attr_name)
   end if

   ! open object
   if (obj_name == "") then
      obj_id = loc_id
   else
      call h5oopen_f(loc_id, obj_name, obj_id, hdferror)
   end if

   ! create attribute
   call h5aopen_f(obj_id, attr_name, attr_id, hdferror)
   !write(*,'(A20,I0)') "h5dcreate: ", hdferror

   ! read attribute
   call h5aread_f(attr_id, H5T_NATIVE_DOUBLE, data, dims, hdferror)
   !write(*,'(A20,I0)') "h5dwrite: ", hdferror

   ! close all id's
   call h5aclose_f(attr_id, hdferror)
   !write(*,'(A20,I0)') "h5dclose: ", hdferror
   if (obj_name /= "") then
      call h5oclose_f(obj_id, hdferror)
   end if

end subroutine hdf_read_attr_double_0


subroutine hdf_read_attr_double_1(loc_id, obj_name, attr_name, data)
   !>   reads 1d array attribute
   integer(HID_T), intent(in) :: loc_id       !> local id in file
   character(len=*), intent(in) :: obj_name   !> object name attribute will be attached to (if "" use loc_id)
   character(len=*), intent(in) :: attr_name  !> name of attribute
   real(dp), intent(out) :: data(:)           !> data to write to attribute

   integer :: rank
   integer(SIZE_T) :: dims(1)
   integer(HID_T) :: obj_id, attr_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_read_attr_double_1: " // trim(obj_name) // "/" //trim(attr_name)
   end if

   ! open object
   if (obj_name == "") then
      obj_id = loc_id
   else
      call h5oopen_f(loc_id, obj_name, obj_id, hdferror)
   end if

   ! create dataspace
   rank = 1
   dims = shape(data, KIND=HID_T)

   ! create attribute
   call h5aopen_f(obj_id, attr_name, attr_id, hdferror)
   !write(*,'(A20,I0)') "h5dcreate: ", hdferror

   ! read attribute
   call h5aread_f(attr_id, H5T_NATIVE_DOUBLE, data, dims, hdferror)
   !write(*,'(A20,I0)') "h5dwrite: ", hdferror

   ! close all id's
   call h5aclose_f(attr_id, hdferror)
   !write(*,'(A20,I0)') "h5dclose: ", hdferror
   if (obj_name /= "") then
      call h5oclose_f(obj_id, hdferror)
   end if

end subroutine hdf_read_attr_double_1


subroutine hdf_read_attr_integer_0(loc_id, obj_name, attr_name, data)
   !>   writes a scalar attribute
   integer(HID_T), intent(in) :: loc_id       !> local id in file
   character(len=*), intent(in) :: obj_name   !> object name attribute will be attached to (if "" use loc_id)
   character(len=*), intent(in) :: attr_name  !> name of attribute
   integer, intent(out) :: data               !> data to write to attribute

   integer(SIZE_T) :: dims(1)
   integer(HID_T) :: obj_id, attr_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_read_attr_integer_0: " // trim(obj_name) // "/" //trim(attr_name)
   end if

   ! open object
   if (obj_name == "") then
      obj_id = loc_id
   else
      call h5oopen_f(loc_id, obj_name, obj_id, hdferror)
   end if

   ! create attribute
   call h5aopen_f(obj_id, attr_name, attr_id, hdferror)
   !write(*,'(A20,I0)') "h5dcreate: ", hdferror

   ! read attribute
   call h5aread_f(attr_id, H5T_NATIVE_INTEGER, data, dims, hdferror)
   !write(*,'(A20,I0)') "h5dwrite: ", hdferror

   ! close all id's
   call h5aclose_f(attr_id, hdferror)
   !write(*,'(A20,I0)') "h5dclose: ", hdferror
   if (obj_name /= "") then
      call h5oclose_f(obj_id, hdferror)
   end if

end subroutine hdf_read_attr_integer_0


subroutine hdf_read_attr_integer_1(loc_id, obj_name, attr_name, data)
   !>   writes a scalar attribute
   integer(HID_T), intent(in) :: loc_id       !> local id in file
   character(len=*), intent(in) :: obj_name   !> object name attribute will be attached to (if "" use loc_id)
   character(len=*), intent(in) :: attr_name  !> name of attribute
   integer, intent(out) :: data(:)            !> data to write to attribute

   integer :: rank
   integer(SIZE_T) :: dims(1)
   integer(HID_T) :: obj_id, attr_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_read_attr_integer_1: " // trim(obj_name) // "/" //trim(attr_name)
   end if

   ! open object
   if (obj_name == "") then
      obj_id = loc_id
   else
      call h5oopen_f(loc_id, obj_name, obj_id, hdferror)
   end if

   ! create dataspace
   rank = 1
   dims = shape(data, KIND=HID_T)

   ! create attribute
   call h5aopen_f(obj_id, attr_name, attr_id, hdferror)
   !write(*,'(A20,I0)') "h5dcreate: ", hdferror

   ! read attribute
   call h5aread_f(attr_id, H5T_NATIVE_INTEGER, data, dims, hdferror)
   !write(*,'(A20,I0)') "h5dwrite: ", hdferror

   ! close all id's
   call h5aclose_f(attr_id, hdferror)
   !write(*,'(A20,I0)') "h5dclose: ", hdferror
   if (obj_name /= "") then
      call h5oclose_f(obj_id, hdferror)
   end if

end subroutine hdf_read_attr_integer_1


subroutine hdf_read_attr_string(loc_id, obj_name, attr_name, data)
   !>   writes a string attribute
   integer(HID_T), intent(in) :: loc_id        !> local id in file
   character(len=*), intent(in) :: obj_name    !> object name attribute will be attached to (if "" use loc_id)
   character(len=*), intent(in) :: attr_name   !> name of attribute
   character(len=*), intent(out) :: data       !> data to write to attribute

   integer(SIZE_T) :: dims(1)
   integer(HID_T) :: obj_id, type_id, attr_id
   integer :: hdferror

   if (hdf_print_messages) then
      write(*,'(A)') "--->hdf_read_attr_string: " // trim(obj_name) // "/" //trim(attr_name)
   end if

   ! open object
   if (obj_name == "") then
      obj_id = loc_id
   else
      call h5oopen_f(loc_id, obj_name, obj_id, hdferror)
   end if

   ! create type_id
   dims(1) = len(data, KIND=HID_T)
   call h5tcopy_f (H5T_NATIVE_CHARACTER, type_id, hdferror)
   !write(*,*) 'h5tcopy_f returns', type_id
   call h5tset_size_f (type_id, dims(1), hdferror)

   ! create attribute
   call h5aopen_f(obj_id, attr_name, attr_id, hdferror)
   !write(*,'(A20,I0)') "h5dcreate: ", hdferror

   ! read attribute
   call h5aread_f(attr_id, type_id, data, dims, hdferror)
   !write(*,'(A20,I0)') "h5dwrite: ", hdferror

   ! close all id's
   call h5aclose_f(attr_id, hdferror)
   !write(*,'(A20,I0)') "h5dclose: ", hdferror
   if (obj_name /= "") then
      call h5oclose_f(obj_id, hdferror)
   end if

end subroutine hdf_read_attr_string

end module hdf5_utils_m
