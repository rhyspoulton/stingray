! Fortran HDF5 Module
! Standard functions for reading and writing HDF5 data
! Author: Danail Obreschkow, 2019
   
module module_hdf5

   use hdf5
   
   private 
   public   :: hdf5_create, hdf5_open, hdf5_close
   public   :: hdf5_dataset_size
   public   :: hdf5_add_group
   public   :: hdf5_read_data, hdf5_write_data
   public   :: hdf5_read_attribute
   public   :: hdf5_read_dataset_compound

   interface hdf5_read_attribute
      module procedure read_attribute_int4
      module procedure read_attribute_real4
   end interface hdf5_read_attribute

   interface hdf5_read_data
      module procedure read_dataset_0d_string
      module procedure read_dataset_0d_int4
      module procedure read_dataset_0d_int8
      module procedure read_dataset_0d_real4
      module procedure read_dataset_0d_real8
      module procedure read_dataset_1d_int4
      module procedure read_dataset_1d_int8
      module procedure read_dataset_1d_real4
      module procedure read_dataset_1d_real8
      module procedure read_dataset_2d_real4
   end interface hdf5_read_data
   
   interface hdf5_write_data
      module procedure write_dataset_0d_int4
      module procedure write_dataset_0d_int8
      module procedure write_dataset_0d_real4
      module procedure write_dataset_0d_real8
      module procedure write_dataset_0d_string
      module procedure write_dataset_1d_int4
      module procedure write_dataset_1d_int8
      module procedure write_dataset_1d_real4
      module procedure write_dataset_1d_real8
      module procedure write_dataset_2d_real4
   end interface hdf5_write_data
   
   integer(hid_t) :: file_id
   
   integer*4,allocatable   :: i4(:)
   integer*8,allocatable   :: i8(:)
   real*4,allocatable      :: r4(:)
   real*8,allocatable      :: r8(:)
   
contains

   function exists(filename,do_not_stop) result(res)
      implicit none
      character(len=*),intent(in)   :: filename
      logical,intent(in),optional   :: do_not_stop
      logical                       :: res
      inquire(file=trim(filename), exist=res)
      if ((.not.res).and.(.not.present(do_not_stop))) then
         write(*,'(A)') 'HDF5 ERROR: File does not exist: '//trim(filename)
         stop
      end if
   end function exists
   
   subroutine error(txt)
      implicit none
      character(*),intent(in) :: txt
      write(*,'(A)') 'HDF5 ERROR: '//txt
      stop
   end subroutine error

   subroutine hdf5_open(filename,write_access) ! opens an existing HDF5 file for read and write
   
      implicit none
      character(*),intent(in)       :: filename
      logical,optional,intent(in)   :: write_access
      integer*4                     :: status
      logical                       :: wa
              
      wa = .false.
      if (present(write_access)) then
         wa = write_access
      end if

      if (exists(filename)) then
         call h5open_f(status) ! Initialize the Fortran interface
         if (wa) then
            call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,status) ! Open an hdf5 file
         else
            call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,status) ! Open an hdf5 file
         end if
      end if
      
   end subroutine hdf5_open
   
   subroutine hdf5_create(filename) ! creates a new HDF5 file and closes it
   
      implicit none
      character(*),intent(in)    :: filename
      integer*4                  :: status

      call h5open_f(status) ! Initialize the Fortran interface
      call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,status) ! Create an hdf5 file
      call hdf5_close()
      
   end subroutine hdf5_create
   
   subroutine hdf5_add_group(groupname) ! adds a new group to the open file (=> hdf5_open must be called first)
   
      implicit none
      character(*),intent(in)    :: groupname
      integer(hid_t)             :: group_id
      integer*4                  :: error
   
      call h5gcreate_f(file_id, groupname, group_id, error)
      call h5gclose_f(group_id, error)
      
   end subroutine hdf5_add_group
   
   subroutine hdf5_close()
   
      implicit none
      integer*4                  :: status
      integer*4                  :: err
      call h5fclose_f(file_id,err) ! Terminate access to the file
      call h5close_f(status) ! Close the Fortran interface
      
   end subroutine hdf5_close
   
   integer*8 function hdf5_dataset_size(dataset)
   
      implicit none
      character(*),intent(in)          :: dataset
      integer*4                        :: err,rankr
      integer(hid_t)                   :: dataset_id
      integer(hid_t)                   :: dataspace
      integer(hsize_t), dimension(2)   :: dimsr, maxdimsr
   
      call h5dopen_f(file_id, dataset, dataset_id, err) ! Open dataset
      call h5dget_space_f(dataset_id, dataspace, err)
      call h5sget_simple_extent_ndims_f(dataspace, rankr, err)
      call h5sget_simple_extent_dims_f(dataspace, dimsr, maxdimsr, err)
      call h5dclose_f(dataset_id, err)
      hdf5_dataset_size = dimsr(1)
      
   end function hdf5_dataset_size
   
   ! ===========================================================================================================
   ! subroutines used by interface hdf5_write_dataset
   ! ===========================================================================================================
   
   subroutine write_dataset_0d_int4(datasetname,dat,attribute)
   
      implicit none
      character(*),intent(in)          :: datasetname
      integer*4,intent(in)             :: dat
      character(*),intent(in),optional :: attribute ! optional explanation
      
      call write_dataset_1d_int4(datasetname,(/dat/),attribute)
      
   end subroutine write_dataset_0d_int4
   
   subroutine write_dataset_0d_int8(datasetname,dat,attribute)
   
      implicit none
      character(*),intent(in)          :: datasetname
      integer*8,intent(in)             :: dat
      character(*),intent(in),optional :: attribute ! optional explanation
      
      call write_dataset_1d_int8(datasetname,(/dat/),attribute)
      
   end subroutine write_dataset_0d_int8
   
   subroutine write_dataset_0d_real4(datasetname,dat,attribute)
   
      implicit none
      character(*),intent(in)          :: datasetname
      real*4,intent(in)                :: dat
      character(*),intent(in),optional :: attribute ! optional explanation
      
      call write_dataset_1d_real4(datasetname,(/dat/),attribute)
      
   end subroutine write_dataset_0d_real4
   
   subroutine write_dataset_0d_real8(datasetname,dat,attribute)
   
      implicit none
      character(*),intent(in)          :: datasetname
      real*8,intent(in)                :: dat
      character(*),intent(in),optional :: attribute ! optional explanation
      
      call write_dataset_1d_real8(datasetname,(/dat/),attribute)
      
   end subroutine write_dataset_0d_real8
   
   subroutine write_dataset_1d_int4(datasetname,dat,attribute)
   
      implicit none
      character(*),intent(in)          :: datasetname
      integer*4,intent(in)             :: dat(:)
      character(*),intent(in),optional :: attribute ! optional explanation
      integer(hid_t)                   :: space_id
      integer(hid_t)                   :: dataset_id
      integer*4                        :: error
      integer*4,parameter              :: rank = 1
      integer*8                        :: n(1)
      
      n(1) = size(dat,1)
      call h5screate_simple_f(rank,n,space_id,error)  ! Create the space for the dataset
      call h5dcreate_f(file_id,datasetname,H5T_STD_I32LE,space_id,dataset_id,error) ! create dataset
      call h5dwrite_f(dataset_id,H5T_STD_I32LE,dat,n,error) ! write data set
      call h5sclose_f(space_id,error) ! close data space
      if (present(attribute)) call write_attribute(attribute,dataset_id)
      call h5dclose_f(dataset_id,error) ! close data set
      
   end subroutine write_dataset_1d_int4
   
   subroutine write_dataset_1d_int8(datasetname,dat,attribute)
   
      implicit none
      character(*),intent(in)          :: datasetname
      integer*8,intent(in)             :: dat(:)
      character(*),intent(in),optional :: attribute ! optional explanation
      integer(hid_t)                   :: space_id
      integer(hid_t)                   :: dataset_id
      integer*4                        :: error
      integer*4,parameter              :: rank = 1
      integer*8                        :: n(1)
      
      n(1) = size(dat,1)
      call h5screate_simple_f(rank,n,space_id,error)  ! Create the space for the dataset
      call h5dcreate_f(file_id,datasetname,H5T_STD_I64LE,space_id,dataset_id,error) ! create dataset
      call h5dwrite_f(dataset_id,H5T_STD_I64LE,dat,n,error) ! write data set
      call h5sclose_f(space_id,error) ! close data space
      if (present(attribute)) call write_attribute(attribute,dataset_id)
      call h5dclose_f(dataset_id,error) ! close data set
      
   end subroutine write_dataset_1d_int8
           
   subroutine write_dataset_1d_real4(datasetname,dat,attribute)
   
      implicit none
      character(*),intent(in)          :: datasetname
      real*4,intent(in)                :: dat(:)
      character(*),intent(in),optional :: attribute ! optional explanation
      integer(hid_t)                   :: space_id
      integer(hid_t)                   :: dataset_id
      integer*4                        :: error
      integer*4,parameter              :: rank = 1
      integer*8                        :: n(1)
      
      n(1) = size(dat,1)
      call h5screate_simple_f(rank,n,space_id,error)  ! Create the space for the dataset
      call h5dcreate_f(file_id,datasetname,H5T_IEEE_F32LE,space_id,dataset_id,error) ! create dataset
      call h5dwrite_f(dataset_id,H5T_IEEE_F32LE,dat,n,error) ! write data set
      call h5sclose_f(space_id,error) ! close data space
      if (present(attribute)) call write_attribute(attribute,dataset_id)
      call h5dclose_f(dataset_id,error) ! close data set
      
   end subroutine write_dataset_1d_real4
   
   subroutine write_dataset_1d_real8(datasetname,dat,attribute)
   
      implicit none
      character(*),intent(in)          :: datasetname
      real*8,intent(in)                :: dat(:)
      character(*),intent(in),optional :: attribute ! optional explanation
      integer(hid_t)                   :: space_id
      integer(hid_t)                   :: dataset_id
      integer*4                        :: error
      integer*4,parameter              :: rank = 1
      integer*8                        :: n(1)
      
      n(1) = size(dat,1)
      call h5screate_simple_f(rank,n,space_id,error)  ! Create the space for the dataset
      call h5dcreate_f(file_id,datasetname,H5T_IEEE_F64LE,space_id,dataset_id,error) ! create dataset
      call h5dwrite_f(dataset_id,H5T_IEEE_F64LE,dat,n,error) ! write data set
      call h5sclose_f(space_id,error) ! close data space
      if (present(attribute)) call write_attribute(attribute,dataset_id)
      call h5dclose_f(dataset_id,error) ! close data set
      
   end subroutine write_dataset_1d_real8
   
   subroutine write_dataset_0d_string(datasetname,dat,attribute)
   
      implicit none
      character(*),intent(in)          :: datasetname
      character(*),intent(in)          :: dat
      character(*),intent(in),optional :: attribute ! optional explanation
      integer(hid_t)                   :: space_id
      integer(hid_t)                   :: dataset_id,filetype
      integer*4                        :: error,hdferr
      integer*4,parameter              :: rank = 1
      integer*8                        :: n(1)
      
      n(1) = len(dat)
      call h5screate_simple_f(rank,(/1_8/),space_id,error)  ! Create the space for the dataset
      call H5Tcopy_f(H5T_FORTRAN_S1, filetype, hdferr)
      call H5Tset_size_f(filetype, n(1), hdferr)
      call h5dcreate_f(file_id,datasetname,filetype,space_id,dataset_id,error) ! create dataset
      call h5dwrite_f(dataset_id,filetype,dat,n,error) ! write data set
      call h5sclose_f(space_id,error) ! close data space
      if (present(attribute)) call write_attribute(attribute,dataset_id)
      call h5dclose_f(dataset_id,error) ! close data set
      
   end subroutine write_dataset_0d_string
   
   subroutine write_attribute(attribute,dataset_id)
   
      implicit none
      character(*),intent(in)       :: attribute
      integer(hid_t),intent(in)     :: dataset_id
      integer(hid_t)                :: attr_id ! Attribute identifier
      integer(hid_t)                :: aspace_id ! Attribute Dataspace identifier
      integer(hid_t)                :: atype_id ! Attribute Dataspace identifier
      integer*8,dimension(1)        :: adims = (/1/) ! Attribute dimension
      integer*4,parameter           :: arank = 1 ! Attribure rank
      integer(size_t)               :: attrlen ! Length of the attribute string
      character(*),parameter        :: aname = "Comment" ! Attribute name
      integer(hsize_t),dimension(1) :: data_dims = (/1/)
      integer*4                     :: error
      
      attrlen = len(attribute)
      call h5screate_simple_f(arank,adims,aspace_id,error) ! Create scalar data space for the attribute.
      call h5tcopy_f(H5T_NATIVE_CHARACTER,atype_id,error) ! Create datatype for the attribute.
      call h5tset_size_f(atype_id,attrlen,error)
      call h5acreate_f(dataset_id,aname,atype_id,aspace_id,attr_id,error)
      call h5awrite_f(attr_id,atype_id,attribute,data_dims,error) ! Write the attribute data.
      call h5aclose_f(attr_id,error) ! Close the attribute.
      
   end subroutine write_attribute
   
   subroutine write_dataset_2d_real4(datasetname,dat,attribute)
   
      implicit none
      character(*),intent(in)          :: datasetname
      real*4,intent(in)                :: dat(:,:)
      character(*),intent(in),optional :: attribute ! optional explanation
      integer(hid_t)                   :: space_id
      integer(hid_t)                   :: dataset_id
      integer*4                        :: error
      integer*4,parameter              :: rank = 2
      integer*8                        :: n(2)
      
      n(1) = size(dat,1)
      n(2) = size(dat,2)
      call h5screate_simple_f(rank,n,space_id,error)  ! Create the space for the dataset
      call h5dcreate_f(file_id,datasetname,H5T_IEEE_F32LE,space_id,dataset_id,error) ! create dataset
      call h5dwrite_f(dataset_id,H5T_IEEE_F32LE,dat,n,error) ! write data set
      call h5sclose_f(space_id,error) ! close data space
      if (present(attribute)) call write_attribute(attribute,dataset_id)
      call h5dclose_f(dataset_id,error) ! close data set
      
   end subroutine write_dataset_2d_real4

   ! ===========================================================================================================
   ! subroutines used by interface hdf5_read_attribute
   ! ===========================================================================================================

   subroutine read_attribute_int4(attribute,buff)

      implicit none
      character(*),intent(in)             :: attribute
      integer(hid_t)                      :: attribute_id
      integer*4                           :: err
      integer*4                           :: buff
      integer(hsize_t), dimension(1)      :: dims = (/1/)

      CALL h5aopen_f(file_id, attribute, attribute_id, err)
      CALL h5aread_f(attribute_id, H5T_NATIVE_INTEGER, buff,dims, err)
      CALL h5aclose_f(attribute_id , err)

   end subroutine

   subroutine read_attribute_real4(group,attribute,buff)

      implicit none
      character(*),intent(in)             :: attribute,group
      integer(hid_t)                      :: attribute_id, group_id
      integer*4                           :: err
      real*4                              :: buff
      real*8                              :: buff2
      integer(hsize_t), dimension(1)      :: dims = (/1/)

      CALL h5gopen_f(file_id, group, group_id, err)
      CALL h5aopen_f(group_id, attribute, attribute_id, err)
      CALL h5aread_f(attribute_id, H5T_NATIVE_DOUBLE, buff2,dims, err)
      CALL h5aclose_f(attribute_id , err)
      CALL h5gclose_f(group_id, err)

      ! Cast this to a real4 for the moment
      buff=buff2

   end subroutine

   
   
   ! ===========================================================================================================
   ! subroutines used by interface hdf5_read_dataset
   ! ===========================================================================================================
   
   subroutine read_dataset_0d_string(dataset, dat, convert)
   
      implicit none
      character(*),intent(in)       :: dataset
      character(255),intent(inout)  :: dat
      logical,intent(in),optional   :: convert
      integer*4                     :: err,hdferr
      integer(hsize_t)              :: size(2)
      integer(hid_t)                :: dataset_id,filetype
      integer*8                     :: n(1)
      
      if (present(convert)) call error('String variables cannot be converted in HDF5.')
      
      n(1) = 255
      
      call h5dopen_f(file_id, dataset, dataset_id, err)
      call H5Tcopy_f(H5T_FORTRAN_S1, filetype, hdferr)
      call H5Tset_size_f(filetype, n(1), hdferr)
      call h5dread_f(dataset_id, filetype, dat, size, err)
      call h5dclose_f(dataset_id, err)
   
   end subroutine read_dataset_0d_string
   
   subroutine read_dataset_0d_int4(dataset, dat, field_name, convert)
   
      implicit none
      character(*),intent(in)          :: dataset
      integer*4,intent(inout)          :: dat
      logical,intent(in),optional      :: convert
      character(*),intent(in),optional :: field_name ! Required if reading a compound datatype
      character(2),parameter           :: expected_type = 'i4'
      
      call read_numeric(dataset, 1, expected_type, field_name, convert)
      dat = i4(1)
      deallocate(i4)
      
   end subroutine read_dataset_0d_int4
   
   subroutine read_dataset_0d_int8(dataset, dat, field_name, convert)
   
      implicit none
      character(*),intent(in)          :: dataset
      integer*8,intent(inout)          :: dat
      logical,intent(in),optional      :: convert
      character(*),intent(in),optional :: field_name ! Required if reading a compound datatype
      character(2),parameter           :: expected_type = 'i8'
      
      call read_numeric(dataset, 1, expected_type, field_name, convert)
      dat = i8(1)
      deallocate(i8)
      
   end subroutine read_dataset_0d_int8
   
   subroutine read_dataset_0d_real4(dataset, dat, field_name, convert)
   
      implicit none
      character(*),intent(in)          :: dataset
      real*4,intent(inout)             :: dat
      logical,intent(in),optional      :: convert
      character(*),intent(in),optional :: field_name ! Required if reading a compound datatype
      character(2),parameter           :: expected_type = 'r4'
      
      call read_numeric(dataset, 1, expected_type, field_name, convert)
      dat = r4(1)
      deallocate(r4)
      
   end subroutine read_dataset_0d_real4
   
   subroutine read_dataset_0d_real8(dataset, dat, field_name, convert)
   
      implicit none
      character(*),intent(in)          :: dataset
      real*8,intent(inout)             :: dat
      logical,intent(in),optional      :: convert
      character(*),intent(in),optional :: field_name ! Required if reading a compound datatype
      character(2),parameter           :: expected_type = 'r8'
      
      call read_numeric(dataset, 1, expected_type, field_name, convert)
      dat = r8(1)
      deallocate(r8)
      
   end subroutine read_dataset_0d_real8
   
   subroutine read_dataset_1d_int4(dataset, dat, field_name, convert)
   
      implicit none
      character(*),intent(in)          :: dataset
      integer*4,intent(inout)          :: dat(:)
      logical,intent(in),optional      :: convert
      character(*),intent(in),optional :: field_name ! Required if reading a compound datatype
      character(2),parameter           :: expected_type = 'i4'
      
      call read_numeric(dataset, int(size(dat),4), expected_type, field_name, convert)
      dat = i4
      deallocate(i4)
      
   end subroutine read_dataset_1d_int4
   
   subroutine read_dataset_1d_int8(dataset, dat, field_name, convert)
   
      implicit none
      character(*),intent(in)          :: dataset
      integer*8,intent(inout)          :: dat(:)
      logical,intent(in),optional      :: convert
      character(*),intent(in),optional :: field_name ! Required if reading a compound datatype
      character(2),parameter           :: expected_type = 'i8'
      
      call read_numeric(dataset, int(size(dat),4), expected_type, field_name, convert)
      dat = i8
      deallocate(i8)
      
   end subroutine read_dataset_1d_int8
   
   subroutine read_dataset_1d_real4(dataset, dat, field_name, convert)
   
      implicit none
      character(*),intent(in)          :: dataset
      real*4,intent(inout)             :: dat(:)
      logical,intent(in),optional      :: convert
      character(*),intent(in),optional :: field_name ! Required if reading a compound datatype
      character(2),parameter           :: expected_type = 'r4'
      
      call read_numeric(dataset, int(size(dat),4), expected_type, field_name, convert)
      dat = r4
      deallocate(r4)
      
   end subroutine read_dataset_1d_real4
   
   subroutine read_dataset_1d_real8(dataset, dat, field_name, convert)
   
      implicit none
      character(*),intent(in)          :: dataset
      real*8,intent(inout)             :: dat(:)
      logical,intent(in),optional      :: convert
      character(*),intent(in),optional :: field_name ! Required if reading a compound datatype
      character(2),parameter           :: expected_type = 'r8'

      call read_numeric(dataset, int(size(dat),4), expected_type, field_name, convert)
      dat = r8
      deallocate(r8)
      
   end subroutine read_dataset_1d_real8

   subroutine read_dataset_2d_real4(dataset, dat, field_name, convert)
   
      implicit none
      character(*),intent(in)          :: dataset
      real*4,intent(inout)             :: dat(:,:)
      logical,intent(in),optional      :: convert
      character(*),intent(in),optional :: field_name ! Required if reading a compound datatype
      character(2),parameter           :: expected_type = 'r4'
      character(2)                     :: detected_type
      integer*4                        :: err, class, sign, index
      integer(hsize_t)                 :: array_dims(2)
      integer(hid_t)                   :: dataset_id, datatype_id,compound_datatype_id
      integer(hid_t)                   :: hdf_type
      character(255)                   :: msg
      integer(size_t)                  :: size
      integer(size_t)                  :: offset=0
      TYPE(C_PTR) :: f_ptr

      call h5dopen_f(file_id, dataset, dataset_id, err)

      ! Open up the dataset and extract its datatype
      call h5dget_type_f(dataset_id,compound_datatype_id,err)

      ! Find the index of the member we want to read and find its type
      call h5tget_member_index_f(compound_datatype_id,field_name,index,err)
      call h5tget_member_type_f(compound_datatype_id,index,datatype_id,err)

      call h5tget_array_dims_f(datatype_id, array_dims, err)

      ! Get the array datatype
      call h5tget_size_f(datatype_id, size, err)

      call h5tcreate_f(H5T_COMPOUND_F, size, hdf_type, err)
      call h5tinsert_f(hdf_type, field_name, offset, datatype_id, err)

      call h5dread_f(dataset_id, hdf_type, dat, array_dims, err)

      call h5dclose_f(dataset_id, err)

   end subroutine read_dataset_2d_real4

   subroutine read_numeric(dataset, n, expected_type, field_name, convert)
   
      character(*),intent(in)             :: dataset
      integer*4,intent(in)                :: n
      character(2),intent(in)             :: expected_type
      logical,intent(in),optional         :: convert
      character(*),intent(in),optional    :: field_name ! Required if reading a compound datatype
      character(2)                        :: detected_type
      integer*4                           :: err, class, sign, index
      integer(hsize_t)                    :: size(2)
      integer(hid_t)                      :: dataset_id
      integer(hid_t)                      :: hdf_type
      character(255)                      :: msg

      call h5dopen_f(file_id, dataset, dataset_id, err)
      call get_mem_type_id(dataset_id, hdf_type, detected_type, field_name)
      
      if (detected_type.ne.expected_type) then
         write(msg,'(7A)') 'Expected type ',expected_type,' but detected type ', &
         & detected_type,' in dataset ',dataset, field_name
         if (present(convert)) then
            if (.not.convert) then
               call error(trim(msg))
            end if
         else 
            call error(trim(msg))
         end if
      end if
      
      select case(detected_type)
      case('i4')
         allocate(i4(n))
         call h5dread_f(dataset_id, hdf_type, i4, size, err)
         if (expected_type.ne.detected_type) then
            if (expected_type=='i8') i8 = int(i4,8)
            if (expected_type=='r4') r4 = real(i4,4)
            if (expected_type=='r8') r8 = real(i4,8)
            deallocate(i4)
         end if
      case('i8')
         allocate(i8(n))
         call h5dread_f(dataset_id, hdf_type, i8, size, err)
         if (expected_type.ne.detected_type) then
            if (expected_type=='i4') i4 = int(i8,4)
            if (expected_type=='r4') r4 = real(i8,4)
            if (expected_type=='r8') r8 = real(i8,8)
            deallocate(i8)
         end if
      case('r4')
         allocate(r4(n))
         call h5dread_f(dataset_id, hdf_type, r4, size, err)
         if (expected_type.ne.detected_type) then
            if (expected_type=='i4') i4 = int(r4,4)
            if (expected_type=='i8') i8 = int(r4,8)
            if (expected_type=='r8') r8 = real(r4,8)
            deallocate(r4)
         end if
      case('r8')
         allocate(r8(n))
         call h5dread_f(dataset_id, hdf_type, r8, size, err)
         if (expected_type.ne.detected_type) then
            if (expected_type=='i4') i4 = int(r8,4)
            if (expected_type=='i8') i8 = int(r8,8)
            if (expected_type=='r4') r4 = real(r8,4)
            deallocate(r8)
         end if
      case default
         write(msg,'(4A)') 'Unknown type ',detected_type,' in dataset ',dataset
         call error(trim(msg))
      end select
      
      call h5dclose_f(dataset_id, err)
      
   end subroutine read_numeric

   
   subroutine get_mem_type_id(dataset_id, hdf_type, fortran_type, field_name)

      implicit none
   
      integer(hid_t),intent(in)        :: dataset_id
      integer(hid_t),intent(out)       :: hdf_type
      character(2),intent(out)         :: fortran_type
      character(*),intent(in),optional :: field_name ! Required if reading a compound datatype
      integer*4                        :: err, class, sign, index
      integer(size_t)                  :: size
      integer(hid_t)                   :: datatype_id
      integer(hid_t)                   :: compound_datatype_id
      integer(hid_t)                   :: base_datatype_id
      integer(size_t)                  :: array_dims(1)
      integer(size_t)                  :: offset=0

      if(field_name=='') then
         call h5dget_type_f(dataset_id, datatype_id, err)
         call h5tget_class_f(datatype_id, class, err)
         call h5tget_size_f(datatype_id, size, err)
         call h5tget_sign_f(datatype_id, sign, err)
      else
         ! Open up the dataset and extract its datatype
         call h5dget_type_f(dataset_id,compound_datatype_id,err)

         ! Find the index of the member we want to read and find its type
         call h5tget_member_index_f(compound_datatype_id,field_name,index,err)
         call h5tget_member_type_f(compound_datatype_id,index,datatype_id,err)

         ! Get the info for the member
         call h5tget_class_f(datatype_id, class, err)
         call h5tget_size_f(datatype_id, size, err)
         call h5tget_sign_f(datatype_id, sign, err)
      end if

      ! First check if the class is an array, if so we need get the datatype in the array
      if (class == h5t_array_f) then
         call h5tget_super_f(datatype_id,base_datatype_id,err)
         call h5tget_class_f(base_datatype_id, class, err)
         call h5tget_size_f(base_datatype_id, size, err)
         call h5tget_sign_f(base_datatype_id, sign, err)
      end if

      if (class == h5t_integer_f) then ! integer
         if (size == 4) then
            if (sign == 0) then
               hdf_type = H5T_STD_U32LE
            else if (sign == 1) then
               hdf_type = H5T_STD_I32LE
            else
               call error('Unknown type')
            end if
         else if (size == 8) then
            if (sign == 0) then
               hdf_type = H5T_STD_U64LE
            else if (sign == 1) then
               hdf_type = H5T_STD_I64LE
            else
               call error('Unknown type')
            end if
         else
            call error('Unknown type')
         end if
      else if (class == h5t_float_f) then ! floating point
         if (size == 4) then
            hdf_type = H5T_IEEE_F32LE
         else if (size == 8) then
            hdf_type = H5T_IEEE_F64LE
         else
            call error('Unknown type')
         end if
      else
         call error('Unknown type')
      end if
      
      if (class == 0) then
         write(fortran_type,'(A,I0)') 'i',size
      else
         write(fortran_type,'(A,I0)') 'r',size
      end if

      ! If there is a field name we need to create a compound datatype to know which field to read in
      if(field_name.ne.'') then
         ! Create a new compound datatype to contain the field we want to read in
         call h5tcreate_f(H5T_COMPOUND_F, size, hdf_type, err)
         call h5tinsert_f(hdf_type, field_name, offset, datatype_id, err)
      end if

   end subroutine get_mem_type_id


   subroutine hdf5_read_dataset_compound(dataset,field_name,dat)

      implicit none
      character(*),intent(in)             :: dataset
      character(*),intent(in)             :: field_name
      integer*4                           :: index
      integer*4                           :: err
      integer*4                           :: nmemb
      integer(hsize_t) ,DIMENSION(1)      :: isize
      integer(hid_t)                      :: dataset_id
      integer(hid_t)                      :: datatype_id
      integer(hid_t)                      :: feild_datatype_id
      integer(hid_t)                      :: compound_datatype_id
      integer(hid_t)                      :: hdf_type
      integer(size_t)                     :: type_size
      integer(size_t)                     :: typesize
      integer(size_t)                     :: offset=0
      character(255)                      :: member_name
      integer*8,intent(inout)             :: dat(:)
      integer*4                           :: n
      n=int(size(dat),4)

      isize(1)=n

      ! Open up the dataset and extract its datatype
      call h5dopen_f(file_id, dataset, dataset_id, err)
      call h5dget_type_f(dataset_id,datatype_id,err)

      ! Find the index of the member we want to read and find its type
      call h5tget_member_index_f(datatype_id,field_name,index,err)
      call h5tget_member_type_f(datatype_id,index,feild_datatype_id,err)

      ! Get the size of the datatype
      call h5tget_size_f(feild_datatype_id,type_size,err)

      ! Create a new compound datatype to contain the field we want to read in
      call h5tcreate_f(H5T_COMPOUND_F, type_size, compound_datatype_id, err)
      call h5tinsert_f(compound_datatype_id, field_name, offset, feild_datatype_id, err)

      allocate(i8(n))
      call h5dread_f(dataset_id, compound_datatype_id, i8 , isize, err)

      call h5tclose_f(datatype_id,err)
      call h5tclose_f(compound_datatype_id,err)
      call h5dclose_f(dataset_id,err)

      dat = i8
      deallocate(i8)

   end subroutine hdf5_read_dataset_compound
   
end module module_hdf5