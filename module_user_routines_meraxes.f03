module module_user_routines

! ==============================================================================================================
! IMPORT MODULES
! ==============================================================================================================

! default modules, do not edit
use module_constants
use module_system
use module_types
use module_io
use module_linalg
use module_cosmology
use module_conversion
use module_emission_lines

! custom modules
use module_hdf5

private

public :: type_sam   ! SAM class, requires procedures get_position, get_groupid, is_selected
public :: type_sky_galaxy
public :: type_sky_group
public :: parameter_filename_default
public :: make_automatic_parameters
public :: make_redshifts
public :: load_sam_snapshot
public :: make_hdf5

! ==============================================================================================================
! SET DEFAULT PARAMETER FILENAME
! ==============================================================================================================

! Here, specify the default parameter filename, used when stingray is called without the optional argument
! "parameterfile"

character(len=255),parameter  :: parameter_filename_default = &
   & '/Users/do/Dropbox/Code/Fortran/stingray/stingray/parameters.txt'


! ==============================================================================================================
! TYPE DECLARATIONS
! ==============================================================================================================

! Here, specify the class of SAM-properties to be loaded for making the mock sky. How these properties are
! read from files is specified in the subroutine load_sam_snapshot below
! The class requires four class functions:
! get_position: returns xyz-position of the galaxy
! get_groupid: returns the group id of the galaxy
! is_group_center: logical function specifying if the galaxy is the group_center

type type_sam

   integer*8   :: ID             ! unique galaxy ID
   integer*8   :: HaloID         ! unique ID of parent halo
   integer*4   :: snapshot       ! snapshot ID
   integer*4   :: core           ! core  index
   integer*4   :: type           ! galaxy type ?
   real*4      :: pos(3)         ! [Mpc/h] position of galaxy centre in simulation box
   real*4      :: vel(3)         ! [proper km/s] peculiar velocity
   real*4      :: J(3)           ! [proper Msun/h pMpc/h km/s] angular momentum
   real*4      :: stellarmass    ! [Msun/h] stellar mass 
   real*4      :: hotgas         ! [Msun/h] hot gass mass
   real*4      :: coldgas        ! [Msun/h] cold gass mass
   real*4      :: HIMass         ! [Msun/h] atomic gas mass disk
   real*4      :: H2Mass         ! [Msun/h] molecular gas mass bulge
   real*4      :: mvir           ! [Msun/h]
   real*4      :: vvir           ! [km/s] virial velocity of subhalo
   real*4      :: vmax           ! [km/s] maximum circular velocity of subhalo
   
contains

   procedure   :: get_position      => sam_get_position     ! required function
   procedure   :: get_groupid       => sam_get_groupid      ! required function               
   procedure   :: is_group_center   => sam_is_group_center  ! required function

end type type_sam

! Here, specify the class or classes of mock object properties in the sky, such as apparent galaxy properties.
! See the existing routines below for clarifications. The class must contain the following class functions:
! make_from_sam
! write_to_file
! is_selected
   
type type_sky

   integer*4   :: snapshot       ! snapshot ID
   integer*4   :: core      ! core index
   integer*4   :: tile           ! tile ID
   real*4      :: zobs           ! redshift in observer-frame
   real*4      :: zcmb           ! redshift in CMB frame
   real*4      :: zcos           ! cosmological redshift without peculiar motions
   real*4      :: dc             ! [simulation units = Mpc/h] comoving distance
   real*4      :: ra             ! [rad] right ascension
   real*4      :: dec            ! [rad] declination
   real*4      :: vpec(3)        ! [proper km/s] 3D peculiar velocity
   real*4      :: vpecrad        ! [proper km/s] radial peculiar velocity
   integer*8   :: HaloID_sam    ! galaxy parent halo ID in the SAM
   integer*8   :: id_group_sky   ! unique group ID in the mock sky
   
   contains

   procedure   :: write_to_file  => sky_write_to_file    ! required subroutine
   
end type type_sky

type,extends(type_sky) :: type_sky_galaxy ! must exist

   integer*8   :: id_galaxy_sky           ! unique ID in the mock sky
   integer*8   :: id_galaxy_sam           ! galaxy ID in the SAM
   integer*4   :: type                    ! galaxy type (0=central, 1=satellite, 2=orphan)
   
   real*4      :: inclination             ! [rad] inclination
   real*4      :: pa                      ! [rad] position angle from North to East
   real*4      :: mag                     ! apparent magnitude (generic: M/L ratio of 1, no k-correction)

   real*4      :: stellarmass             ! [Msun/h] hot gass mass   
   real*4      :: hotgas                  ! [Msun/h] hot gass mass
   real*4      :: coldgas                 ! [Msun/h] cold gass mass
   real*4      :: HIMass                  ! [Msun/h] atomic gas mass disk
   real*4      :: H2Mass                  ! [Msun/h] molecular gas mass bulge
   
   real*4      :: J(3)                    ! [Msun pMpc km/s] total angular momentum
   
   real*4      :: mvir                    ! [Msun/h] mass 
   
   real*4      :: vvir           ! [km/s] virial velocity of subhalo
   real*4      :: vmax           ! [km/s] maximum circular velocity of subhalo
   
   ! HI line
   real*4      :: hiline_flux_int         ! [W/m^2] integrated HI line flux
   real*4      :: hiline_flux_int_vel     ! [Jy km/s] velocity-integrated HI line flux
   real*4      :: hiline_flux_peak        ! [s/km] normalised peak HI line flux density of inclined galaxy (multiply by hiline_flux_int_vel to get Jy)
   real*4      :: hiline_flux_central     ! [s/km] normalised central HI line flux density of inclined galaxy
   real*4      :: hiline_width_peak       ! [km/s] line width of inclined galaxy in rest-frame velocity at peak flux
   real*4      :: hiline_width_50         ! [km/s] line width of inclined galaxy in rest-frame velocity at 50% or peak flux
   real*4      :: hiline_width_20         ! [km/s] line width of inclined galaxy in rest-frame velocity at 20% of peak flux
   real*4      :: hiline_flux_peak_eo     ! [s/km] peak HI line flux density of edge-on galaxy
   real*4      :: hiline_flux_central_eo  ! [s/km] central HI line flux density of edge-on galaxy
   real*4      :: hiline_width_peak_eo    ! [km/s] line width of edge-on galaxy in rest-frame velocity at peak flux
   real*4      :: hiline_width_50_eo      ! [km/s] line width of edge-on galaxy in rest-frame velocity at 50% or peak flux
   real*4      :: hiline_width_20_eo      ! [km/s] line width of edge-on galaxy in rest-frame velocity at 20% of peak flux
   
   ! CO (1-0) line
   real*4      :: coline_flux_int         ! [W/m^2] integrated CO(1-0) line flux
   real*4      :: coline_flux_int_vel     ! [Jy km/s] velocity-integrated CO(1-0) line flux
   real*4      :: coline_flux_peak        ! [s/km] peak CO line flux density of inclined galaxy (multiply by coline_flux_int_vel to get Jy)
   real*4      :: coline_flux_central     ! [s/km] central CO line flux density of inclined galaxy
   real*4      :: coline_width_peak       ! [km/s] line width of inclined galaxy in rest-frame velocity at peak flux
   real*4      :: coline_width_50         ! [km/s] line width of inclined galaxy in rest-frame velocity at 50% or peak flux
   real*4      :: coline_width_20         ! [km/s] line width of inclined galaxy in rest-frame velocity at 20% of peak flux
   real*4      :: coline_flux_peak_eo     ! [s/km] peak CO line flux density of edge-on galaxy
   real*4      :: coline_flux_central_eo  ! [s/km] central CO line flux density of edge-on galaxy
   real*4      :: coline_width_peak_eo    ! [km/s] line width of edge-on galaxy in rest-frame velocity at peak flux
   real*4      :: coline_width_50_eo      ! [km/s] line width of edge-on galaxy in rest-frame velocity at 50% or peak flux
   real*4      :: coline_width_20_eo      ! [km/s] line width of edge-on galaxy in rest-frame velocity at 20% of peak flux
   
   contains
   
   procedure   :: make_from_sam  => make_sky_galaxy   ! required subroutine
   
end type type_sky_galaxy


type,extends(type_sky) :: type_sky_group ! must exist
   
   real*4      :: mvir                 ! [Msun/h] virial mass of group
   real*4      :: vvir                 ! [km/s]	virial velocity of group halo
   real*4      :: vmax                 ! [km/s]	maximum circular velocity of group halo
   real*4      :: cnfw                 ! [-] concentration of NFW fit to group halo
   integer*4   :: group_ntot           ! total number of galaxies in group
   integer*4   :: group_nsel           ! number of selected galaxies in group
   integer*4   :: group_flag           ! 0=group complete, >0 group truncated
   real*4      :: sigma_los_detected   ! [km/s] line-of-sight velocity dispersion of selected galaxies
   real*4      :: sigma_3D_all         ! [km/s] 3D velocity dispersion of all galaxies in group
   
   contains
   
   procedure   :: make_from_sam  => make_sky_group   ! required subroutine
   
end type type_sky_group


! In order to place the objects in the mock sky, the class type_sam must have the following functions
! enabling stingray to extract the position and group id of each object.

contains

function sam_get_position(self) result(pos)
   class(type_sam) :: self
   real*4 :: pos(3)
   pos = self%pos ! [length_unit of simulation] position if the galaxy in the box
end function sam_get_position

integer*8 function sam_get_groupid(self)
   class(type_sam) :: self
   sam_get_groupid = self%HaloID ! unique group identifier
end function sam_get_groupid

logical function sam_is_group_center(self)
   class(type_sam) :: self
   sam_is_group_center = self%type==0
end function sam_is_group_center

! For instances of the user-defined sky-types to be saved, add all sky types in the following subroutine

subroutine sky_write_to_file(self,fileID)
   class(type_sky) :: self
   integer*4,intent(in) :: fileID
   select type (self)
   type is (type_sky_galaxy); write(fileID) self
   type is (type_sky_group); write(fileID) self
   class default
   call error('Unknown class for I/O.')
   end select
end subroutine sky_write_to_file

! ==============================================================================================================
! CONVERSION BETWEEN INTRINSIC AND APPARENT GALAXY PROPERTIES
! ==============================================================================================================

subroutine make_sky_object(sky_object,sam,base,groupid)

   class(type_sky),intent(out)            :: sky_object
   type(type_sam),intent(in)              :: sam
   type(type_base),intent(in)             :: base                 ! basic properties of the position of this galaxy in the sky
   integer*8,intent(in)                   :: groupid              ! unique group in sky index
   real*4                                 :: vector_rotation(3,3) ! rotation matrix for vectors
   real*4                                 :: pos(3)               ! [simulation length units] position vector of galaxy
   real*4                                 :: elos(3)              ! unit vector pointing from the observer to the object in comoving space
   
   call nil(sky_object,sam,base,groupid) ! dummy statement to avoid compiler warnings

   ! sky coordinates
   sky_object%dc  = base%dc*para%box_side  ! [Mpc/h]
   sky_object%ra  = base%ra         ! [rad]
   sky_object%dec = base%dec        ! [rad]
   
   ! make redshift, provided the galaxy position [simulation units] galaxy-velocity [km/s]
   call sph2car(sky_object%dc,sky_object%ra,sky_object%dec,pos)
   elos = pos/norm(pos) ! line-of-sight vector
   call make_redshift(pos*(para%length_unit/Mpc),sam%vel,&
   &zobs=sky_object%zobs,zcmb=sky_object%zcmb,zcos=sky_object%zcos)
   
   ! make IDs
   sky_object%tile          = base%tile
   sky_object%snapshot      = sam%snapshot
   sky_object%core          = sam%core
   sky_object%HaloID_sam    = sam%HaloID
   sky_object%id_group_sky  = groupid         ! unique group id
   
   ! peculiar velocity
   vector_rotation      = tile(base%tile)%Rvector
   sky_object%vpec      = rotate(vector_rotation,sam%vel)
   sky_object%vpecrad   = sum(sky_object%vpec*elos)
   
end subroutine make_sky_object

subroutine make_sky_galaxy(sky_galaxy,sam,base,groupid,galaxyid)

   implicit none
   class(type_sky_galaxy),intent(out)  :: sky_galaxy
   type(type_sam),intent(in)           :: sam
   type(type_base),intent(in)          :: base                 ! basic properties of the position of this galaxy in the sky
   integer*8,intent(in)                :: groupid              ! unique group in sky index
   integer*8,intent(in)                :: galaxyid             ! (only for galaxy types) unique galaxy in sky index
   real*4                              :: pseudo_rotation(3,3) ! rotation matrix for pseudo-vectors (axial vectors)
   real*4                              :: pos(3)               ! [simulation length units] position vector of galaxy
   real*4                              :: dl                   ! [simulation length units] luminosity distance to observer
   real*4                              :: elos(3)              ! unit vector pointing from the observer to the object in comoving space
   real*4                              :: mstars,mHI,mH2,LCO
   
   call nil(sky_galaxy,sam,base,groupid,galaxyid) ! dummy statement to avoid compiler warnings
   
   ! basics
   call make_sky_object(sky_galaxy,sam,base,groupid)
   
   
   ! INTRINSIC PROPERTIES
   
   ! make IDs
   sky_galaxy%id_galaxy_sam   = sam%ID          ! copy other properties
   sky_galaxy%id_galaxy_sky   = galaxyid        ! unique galaxy id
      
   ! basic properties
   sky_galaxy%type             = sam%type
   
   ! intrinsic halo properties
   sky_galaxy%mvir   = sam%mvir
   sky_galaxy%vvir   = sam%vvir
   sky_galaxy%vmax   = sam%vmax
   
   ! intrinsic masses
   sky_galaxy%stellarmass     = sam%stellarmass
   sky_galaxy%HotGas          = sam%HotGas
   sky_galaxy%ColdGas         = sam%ColdGas 
   sky_galaxy%HIMass          = sam%HIMass
   sky_galaxy%H2Mass          = sam%H2Mass
   
   ! intrinsic angular momentum
   pseudo_rotation   = tile(base%tile)%Rpseudo
   sky_galaxy%J      = rotate(pseudo_rotation,sam%J)
   
      
   ! APPARENT PROPERTIES
   
   ! inclination and position angle
   call sph2car(sky_galaxy%dc,sky_galaxy%ra,sky_galaxy%dec,pos)
   elos = pos/norm(pos) ! position vector
   call make_inclination_and_pa(pos,sky_galaxy%J,inclination=sky_galaxy%inclination,pa=sky_galaxy%pa)
   
   ! rough generic optical magnitude
   dl = sky_galaxy%dc*(1+sky_galaxy%zobs) ! [Mpc/h]
   mstars = (sam%stellarmass+sam%stellarmass)/para%h ! [Msun]
   sky_galaxy%mag = convert_absmag2appmag(convert_stellarmass2absmag(mstars,1.0),dl/para%h)
      
   ! cold gas emission lines
   mHI = (sam%HIMass)/1.35 ! [Msun/h] HI mass  
   mH2 = (sam%H2Mass)/1.35 ! [Msun/h] H2 mass
   sky_galaxy%hiline_flux_int = real(convert_luminosity2flux(real(mHI/para%h,8)*real(L2MHI,8)*Lsun,dl/para%h),4) ! [W/m^2]
   sky_galaxy%hiline_flux_int_vel = convert_intflux2velintflux(sky_galaxy%hiline_flux_int,0.21106114,sky_galaxy%zobs)
   LCO = mH2/para%h/(313*X_CO) ! [Jy km/s Mpc^2] CO(1-0) luminosity (note this equation has a J^2 dependence)
   sky_galaxy%coline_flux_int = LCO/(4.0*pi*(dl/para%h)**2)/0.00260075761e23 ! [W/m^2] integrated flux
   sky_galaxy%coline_flux_int_vel = convert_intflux2velintflux(sky_galaxy%hiline_flux_int,0.00260075761,sky_galaxy%zobs)
   if (para%line_parameters==1) call make_line_profiles
   
contains
   
   subroutine make_line_profiles
   
      implicit none
      type(type_lineinput)                :: line_input
      type(type_line)                     :: line(4)
      real*4                              :: a,c_halo,cMpch2pkpc,rvir
      real,parameter                      :: G = 4.3022682e-6     ! [(km/s)^2 kpc/Msun] gravitational constant
   
         
   end subroutine make_line_profiles
   
end subroutine make_sky_galaxy

subroutine make_sky_group(sky_group,sam,sky_galaxy,selected,base,groupid,group_nselected)

   ! the first object in sam, sky_galaxy, selected is the central galaxy of the group;
   ! sky_galaxy only exists for selected objects

   implicit none
   class(type_sky_group),intent(out)   :: sky_group
   type(type_sam),intent(in)           :: sam(:)
   type(type_sky_galaxy),intent(in)    :: sky_galaxy(:)
   logical,intent(in)                  :: selected(:)
   type(type_base),intent(in)          :: base                 ! basic properties of the position of this galaxy in the sky
   integer*8,intent(in)                :: groupid              ! unique group in sky index
   integer*4,intent(in)                :: group_nselected
   
   integer*4                           :: i,n,count
   real*4                              :: dv,v0(3),vr
   
   call nil(sky_group,sam(1),sky_galaxy(1),selected(1),base,groupid,group_nselected) ! dummy statement to avoid compiler warnings
   
end subroutine make_sky_group


! ==============================================================================================================
! IO routines
! ==============================================================================================================

! Set parameters automatically (e.g. from information provided with the snapshot files).
! These automatic parameters are only adopted if the values in the parameter files are set to 'auto'.
! Otherwise the parameter file overwrites these values.

subroutine make_automatic_parameters
   
   implicit none
   
   character(len=255)           :: filename
   integer*4                    :: ncores
   character(*),parameter       :: g = '/InputParams/'  ! Dataset name
   real*4                       :: baryon_frac
   
   filename = ''
   para%snapshot_min = 0
   para%snapshot_max = 163

   write(filename,'(A,A)') trim(para%path_input),'/meraxes.hdf5'


   call out('File of automatic parameters: '//trim(filename))
   call hdf5_open(filename)
   call hdf5_read_attribute('NCores',ncores)
   para%subvolume_min = 0
   para%subvolume_max = 2 !ncores-1
   call hdf5_read_attribute(g,'Hubble_h',para%h)
   call hdf5_read_attribute(g,'OmegaLambda',para%omega_l)
   call hdf5_read_attribute(g,'OmegaM',para%omega_m)
   call hdf5_read_attribute(g,'BaryonFrac',baryon_frac) 
   para%omega_b = para%omega_m * baryon_frac
   call hdf5_read_attribute(g,'BoxSize',para%box_side)
   para%length_unit = Mpc/para%h
   call hdf5_close()
   
end subroutine make_automatic_parameters

! load redshifts
! this routine must write the redshift of each snapshot into the real*4-array
! snapshot(isnapshot)%redshift

subroutine make_redshifts

   implicit none
   character(len=255)                              :: filename,groupname
   integer*4                                       :: isnapshot
   real*4                                          :: z

   write(filename,'(A,A)') trim(para%path_input),'meraxes.hdf5'
   call hdf5_open(filename) ! NB: this routine also checks if the file exists
   
   do isnapshot = para%snapshot_min,para%snapshot_max
      write(groupname,'(A,I0.3)') "Snap",isnapshot
      call hdf5_read_attribute(groupname,'Redshift',z)
      snapshot(isnapshot)%redshift = z
   end do
   call hdf5_close()
   
end subroutine make_redshifts

! load SAM snapshot file

subroutine load_sam_snapshot(index,subindex,sam)

   ! variable declaration
   implicit none
   integer*4,intent(in)                            :: index             ! snapshot index
   integer*4,intent(in)                            :: subindex          ! subindex, if the snapshot is split into several files
   type(type_sam),allocatable,intent(out)          :: sam(:)            ! derived type of all the relevant SAM properties
   character(len=255)                              :: filename
   character(len=17)                               :: g
   integer*8                                       :: n
   real*4,allocatable                              :: buff(:,:)
   
   ! open file
   write(filename,'(A,A,I0,A)') trim(para%path_input),'/meraxes_',subindex,'.hdf5'
   call hdf5_open(filename) ! NB: this routine also checks if the file exists

   write(g,'(A,I0.3,A)') '/Snap',index,"/Galaxies"
   
   ! determine number of galaxies in this (sub)snapshot
   n =  hdf5_dataset_size(g)

   ! allocate array
   if (allocated(sam)) deallocate(sam)
   allocate(sam(n))

   if (allocated(buff)) deallocate(sam)
   allocate(buff(n,3))

   
   ! ! read file
   call hdf5_read_data(g,sam%ID,field_name='ID',convert=.true.)
   call hdf5_read_data(g,sam%HaloID,field_name='HaloID')
   call hdf5_read_data(g,sam%type,field_name='Type')

   !! Extract the Position
   call hdf5_read_data(g,buff,field_name='Pos') 
   sam%pos(1) = buff(:,1)
   sam%pos(2) = buff(:,2)
   sam%pos(3) = buff(:,3)

   !! Extract the Velocity
   call hdf5_read_data(g,buff,field_name='Vel') 
   sam%vel(1) = buff(:,1)
   sam%vel(2) = buff(:,2)
   sam%vel(3) = buff(:,3)

   call hdf5_read_data(g,sam%stellarmass,field_name='StellarMass')
   call hdf5_read_data(g,sam%hotgas,field_name='HotGas')
   call hdf5_read_data(g,sam%coldgas,field_name='ColdGas')
   call hdf5_read_data(g,sam%HIMass,field_name='HIMass')
   call hdf5_read_data(g,sam%H2Mass,field_name='H2Mass')
   call hdf5_read_data(g,sam%mvir,field_name='Mvir')
   call hdf5_read_data(g,sam%vvir,field_name="Vvir")
   call hdf5_read_data(g,sam%vmax,field_name="Vmax")
   
   ! assign other properties
   sam%snapshot = index
   sam%core = subindex

   deallocate(buff)
   
   ! close file
   call hdf5_close()
   
end subroutine load_sam_snapshot

! convert binary files into HDF5

subroutine make_hdf5
   
   implicit none
   logical,parameter                   :: show_test_sum = .false.
   character(len=255)                  :: filename_bin
   character(len=255)                  :: filename_hdf5
   type(type_sky_galaxy),allocatable   :: sky_galaxy(:)
   integer*8                           :: n,i,n_galaxies,n_groups
   character(len=255)                  :: name,str,seedstr
   character(len=255)                  :: filename
   character(len=255)                  :: shark_version,shark_git_revision,shark_timestamp
   real*8                              :: test(0:10),expected_sum
   real*4                              :: n_replica_mean
   integer*4                           :: n_replica_max
   
   call tic
   call out('CONVERT BINARY FILES TO HDF5')
   
   ! load auxilary data
   call load_parameters
   call load_tile_list
   call load_snapshot_list
   
   ! load auxilary data from shark output
   write(filename,'(A,I0,A)') trim(para%path_input),para%snapshot_min,'/0/galaxies.hdf5'
   ! call hdf5_open(filename)
   ! call hdf5_read_data('/run_info/shark_version',shark_version)
   ! call hdf5_read_data('/run_info/shark_git_revision',shark_git_revision)
   ! call hdf5_read_data('/run_info/timestamp',shark_timestamp)
   ! call hdf5_close()

   ! create HDF5 file
   write(seedstr,'(I0)') para%seed
   filename_hdf5 = trim(para%path_output)//'mocksky_'//trim(para%survey)//'_'//trim(seedstr)//'.hdf5'
   call hdf5_create(filename_hdf5)
   
   ! open HDF5 file
   call hdf5_open(filename_hdf5,.true.)
   
   ! Group "Parameters"
   call hdf5_add_group('parameters')
   call hdf5_write_data('parameters/survey',para%survey,'name of simulated survey')
   call hdf5_write_data('parameters/path_output',para%path_output)
   call hdf5_write_data('parameters/path_input',para%path_input)
   call hdf5_write_data('parameters/box_side',para%box_side, &
   & '[length_unit] comoving side-length of simulation box in multiples of length_unit')
   call hdf5_write_data('parameters/length_unit',para%length_unit,'[m] SI-value of comoving length unit')
   call hdf5_write_data('parameters/snapshot_min',para%snapshot_min, &
   & 'index of earliest snapshot used for the mock sky')
   call hdf5_write_data('parameters/snapshot_max',para%snapshot_max, &
   & 'index of latest snapshot used for the mock sky')
   call hdf5_write_data('parameters/core_min',para%subvolume_min, &
   & 'index of first core used for the mock sky')
   call hdf5_write_data('parameters/core_max',para%subvolume_min, &
   & 'index of last core used for the mock sky')
   call hdf5_write_data('parameters/h',para%h,'[-] Hubble parameter H0=h*100 km/s/Mpc')
   call hdf5_write_data('parameters/omega_l',para%omega_l, &
   & 'energy density of dark energy relative to closure density')
   call hdf5_write_data('parameters/omega_m',para%omega_m, &
   & 'energy density of all matter relative to closure density')
   call hdf5_write_data('parameters/omega_b',para%omega_b, &
   & 'energy density of baryonic matter relative to closure density')
   call hdf5_write_data('parameters/dc_min',para%dc_min, &
   & '[length_unit] minimal comoving distance of mock sky')
   call hdf5_write_data('parameters/dc_max',para%dc_max, &
   & '[length_unit] maximal comoving distance of mock sky')
   call hdf5_write_data('parameters/ra_min',para%ra_min/degree, &
   & '[deg] min right ascension of mock sky')
   call hdf5_write_data('parameters/ra_max',para%ra_max/degree, &
   & '[deg] max right ascension of mock sky')
   call hdf5_write_data('parameters/dec_min',para%dec_min/degree, &
   & '[deg] min declination of mock sky')
   call hdf5_write_data('parameters/dec_max',para%dec_max/degree, &
   & '[deg] max declination of mock sky')
   call hdf5_write_data('parameters/zaxis_ra',para%zaxis_ra/degree, &
   & '[deg] RA coordinate of the SAM z-axis in spherical survey coordinates')
   call hdf5_write_data('parameters/zaxis_dec',para%zaxis_dec/degree, &
   & '[deg] Dec coordinate the SAM z-axis in spherical survey coordinates')
   call hdf5_write_data('parameters/xy_angle',para%xy_angle/degree, &
   & '[deg] Rotation of the SAM (x,y)-plane on the sky')
   call hdf5_write_data('parameters/seed',para%seed, &
   & 'seed for the random number generator of symmetry operations')
   call hdf5_write_data('parameters/translate',para%translate, &
   & 'logical flag specifying if random translations are applied (0=false, 1=true)')
   call hdf5_write_data('parameters/rotate',para%rotate, &
   & 'logical flag specifying if random rotations are applied (0=false, 1=true)')
   call hdf5_write_data('parameters/invert',para%invert, &
   & 'logical flag specifying if random inversions are applied (0=false, 1=true)')
   call hdf5_write_data('parameters/velocity_norm',para%velocity_norm, &
   & '[km/s] observer velocity relative to CMB rest-frame')
   call hdf5_write_data('parameters/velocity_ra',para%velocity_ra, &
   & '[deg] RA coordinate to which the observer is moving relative to CMB rest-frame')
   call hdf5_write_data('parameters/velocity_dec',para%velocity_dec, &
   & '[deg] Dec coordinate to which the observer is moving relative to CMB rest-frame')
   call hdf5_write_data('parameters/search_angle',para%search_angle, &
   & '[deg] typical angle in which overlaps between survey volume and tiling grid are searched')
   call hdf5_write_data('parameters/volume_search_level',para%volume_search_level, &
   & 'specifies the number of search points (2^#)^3 inside each tile')
   call hdf5_write_data('parameters/line_parameters',para%line_parameters, &
   & 'logical flag specifying if global emission line parameters are saved (0=false, 1=true)')
   call hdf5_write_data('parameters/keep_binaries',para%keep_binaries, &
   & 'logical flag specifying if binary output files are kept in additino to this HDF5 (0=false, 1=true)')
   call hdf5_write_data('parameters/keep_log',para%keep_binaries, &
   & 'logical flag specifying if the logfile is kept after successful runs (0=false, 1=true)')
   call hdf5_write_data('parameters/skyrotation',para%sky_rotation, &
   & 'Rotation matrix to map xyz-coordinates of the tiling structure onto sky coordinates.')
   
   ! Group "galaxies"
   name = 'galaxies'
   filename_bin = trim(para%path_output)//fn_galaxies
   open(1,file=trim(filename_bin),action='read',form='unformatted',access='stream')
   read(1) n,n_replica_mean,n_replica_max
   n_galaxies = n
   allocate(sky_galaxy(n))
   read(1) sky_galaxy
   close(1)
   call hdf5_add_group(trim(name))
   call hdf5_write_data(trim(name)//'/snapshot',sky_galaxy%snapshot,'snapshot index')
   call hdf5_write_data(trim(name)//'/core',sky_galaxy%core,'core index')
   call hdf5_write_data(trim(name)//'/tile',sky_galaxy%tile,'tile index in tiling array')
   call hdf5_write_data(trim(name)//'/zobs',sky_galaxy%zobs,'redshift in observer-frame')
   call hdf5_write_data(trim(name)//'/zcmb',sky_galaxy%zcmb,'redshift in CMB frame')
   call hdf5_write_data(trim(name)//'/zcos',sky_galaxy%zcos,'cosmological redshift without peculiar motions')
   call hdf5_write_data(trim(name)//'/dc',sky_galaxy%dc,'[Mpc/h] comoving distance')
   call hdf5_write_data(trim(name)//'/ra',sky_galaxy%ra/degree,'[deg] right ascension')
   call hdf5_write_data(trim(name)//'/dec',sky_galaxy%dec/degree,'[deg] declination')
   call hdf5_write_data(trim(name)//'/id_galaxy_sky',sky_galaxy%id_galaxy_sky,'unique galaxy ID in mock sky')
   call hdf5_write_data(trim(name)//'/id_galaxy_sam',sky_galaxy%id_galaxy_sam,'galaxy ID in SAM')
   call hdf5_write_data(trim(name)//'/HaloID_sam',sky_galaxy%HaloID_sam,'host halo ID in SAM')
   call hdf5_write_data(trim(name)//'/type',sky_galaxy%type,'galaxy type (0=central, 1=satellite in halo, 2=orphan)')
   call hdf5_write_data(trim(name)//'/inclination',sky_galaxy%inclination/degree, &
   & '[deg] inclination = angle between line-of-sight and spin axis')
   call hdf5_write_data(trim(name)//'/pa',sky_galaxy%pa/degree,'[deg] position angle from north to east')
   call hdf5_write_data(trim(name)//'/mag',sky_galaxy%mag, &
   & 'apparent magnitude (generic: M/L ratio of 1, no k-correction)')
   call hdf5_write_data(trim(name)//'/vpec_x',sky_galaxy%vpec(1),'[proper km/s] x-component of peculiar velocity')
   call hdf5_write_data(trim(name)//'/vpec_y',sky_galaxy%vpec(2),'[proper km/s] y-component of peculiar velocity')
   call hdf5_write_data(trim(name)//'/vpec_z',sky_galaxy%vpec(3),'[proper km/s] z-component of peculiar velocity')
   call hdf5_write_data(trim(name)//'/vpec_r',sky_galaxy%vpecrad,'[proper km/s] line-of-sight peculiar velocity')
   call hdf5_write_data(trim(name)//'/stellarmass',sky_galaxy%stellarmass,'[Msun/h] stellar mass')
   call hdf5_write_data(trim(name)//'/HotGas',sky_galaxy%HotGas,'[Msun/h] hot gas mass')
   call hdf5_write_data(trim(name)//'/ColdGas',sky_galaxy%ColdGas,'[Msun/h] cold gas mass')
   call hdf5_write_data(trim(name)//'/HIMass',sky_galaxy%HIMass,'[Msun/h] HI mass')
   call hdf5_write_data(trim(name)//'/H2Mass',sky_galaxy%H2Mass,'[Msun/h] H2 mass')
   call hdf5_write_data(trim(name)//'/l_x',sky_galaxy%J(1),'[Msun pMpc km/s] x-component of total angular momentum')
   call hdf5_write_data(trim(name)//'/l_y',sky_galaxy%J(2),'[Msun pMpc km/s] y-component of total angular momentum')
   call hdf5_write_data(trim(name)//'/l_z',sky_galaxy%J(3),'[Msun pMpc km/s] z-component of total angular momentum')
   call hdf5_write_data(trim(name)//'/Mvir',sky_galaxy%Mvir,'[Msun/h] subhalo mass')
   
   if (para%line_parameters==1) then
      
      test(0) = 0
   
   else
   
      test(0) = 0
      
   end if
   
   test(1) = n
   test(3) = sum(sky_galaxy%tile)
   test(4) = sum(sky_galaxy%inclination)
   test(5) = sum(sky_galaxy%zobs)
   deallocate(sky_galaxy)
   
   ! close HDF5 file
   call hdf5_close()
   
   ! check test sum
   if (show_test_sum) then
      expected_sum = 831651.793945312
      if (abs(sum(test)/expected_sum-1)>1e-5) then
         write(str,'(F20.9)') sum(test)
         call out('Test sum FAILED, found = '//trim(str))
      else
         call out('Test sum OK.')
      end if
   end if
   
   call toc

end subroutine make_hdf5

end module module_user_routines