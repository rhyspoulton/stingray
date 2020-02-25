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

   integer*8   :: ID                               ! unique galaxy ID
   integer*8   :: HaloID                           ! unique ID of parent halo
   integer*4   :: snapshot                         ! snapshot ID
   integer*4   :: core                             ! core  index
   integer*4   :: type                             ! galaxy type ?
   real*4      :: pos(3)                           ! [Mpc/h] position of galaxy centre in simulation box
   real*4      :: vel(3)                           ! [proper km/s] peculiar velocity
   real*4      :: StellarMass                      ! [Msun/h yr] Stellar Mass of galaxy
   real*4      :: Sfr                              ! [Msun/h yr] Star formation rate
   real*4      :: BlackHoleMass                    ! [Msun/h] Black hole mass
   real*4      :: BlackHoleAccretedHotMass         ! [Msun/h] Black hole accreted hot mass
   real*4      :: dt                               ! [yr] Time between outputs
   real*4      :: Spin                             ! [Msun/h] Spin of black hole
   real*4      :: mvir                             ! [Msun/h] Virial mass of subhalo
   
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
   integer*4   :: core           ! core index
   integer*4   :: tile           ! tile ID
   real*4      :: zobs           ! redshift in observer-frame
   real*4      :: zcmb           ! redshift in CMB frame
   real*4      :: zcos           ! cosmological redshift without peculiar motions
   real*4      :: dc             ! [simulation units = Mpc/h] comoving distance
   real*4      :: ra             ! [rad] right ascension
   real*4      :: dec            ! [rad] declination
   real*4      :: vpec(3)        ! [proper km/s] 3D peculiar velocity
   real*4      :: vpecrad        ! [proper km/s] radial peculiar velocity
   integer*8   :: HaloID_sam     ! galaxy parent halo ID in the SAM
   integer*8   :: id_group_sky   ! unique group ID in the mock sky
   
   contains

   procedure   :: write_to_file  => sky_write_to_file    ! required subroutine
   
end type type_sky

type,extends(type_sky) :: type_sky_galaxy ! must exist

   integer*8   :: id_galaxy_sky           ! unique ID in the mock sky
   integer*8   :: id_galaxy_sam           ! galaxy ID in the SAM
   integer*4   :: type                    ! galaxy type (0=central, 1=satellite, 2=orphan)

   real*8      :: StellarMass                      ! [Msun/h yr] Stellar Mass of galaxy
   real*8      :: Sfr                              ! [Msun/h yr] Star formation rate
   real*8      :: BlackHoleMass                    ! [Msun/h] Black hole mass
   real*8      :: BlackHoleAccretedHotMass         ! [Msun/h] Black hole accreted hot mass
   real*8      :: dt                               ! [yr] Time between outputs
   real*8      :: Spin                             ! [Msun/h] Spin of black

   real*8      :: fluxBH         ! Black Hole Flux
   real*8      :: fluxStars      ! Stars Flux
   
   real*4      :: mvir           ! [Msun/h] mass
   
   contains
   
   procedure   :: make_from_sam  => make_sky_galaxy   ! required subroutine
   
end type type_sky_galaxy


type,extends(type_sky) :: type_sky_group ! must exist
   
   real*4      :: mvir                 ! [Msun/h] virial mass of group
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
   real*4                              :: LBH,LStars,fluxBH,fluxStars
   real*4                              :: frequency0 = 1420e6
   real*4                              :: frequency = 1500e6
   real*4                              :: mu = 0.8
   real*4                              :: sigma = 0.09

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
   
   ! intrinsic masses
   sky_galaxy%Sfr                      = sam%Sfr
   sky_galaxy%BlackHoleMass            = sam%BlackHoleMass
   sky_galaxy%BlackHoleAccretedHotMass = sam%BlackHoleAccretedHotMass
   sky_galaxy%dt                       = sam%dt
   sky_galaxy%Spin                     = sam%Spin
   
   ! intrinsic angular momentum
   pseudo_rotation   = tile(base%tile)%Rpseudo

   ! APPARENT PROPERTIES

   ! Compute the luminosity of the BH and Stars

   LBH = Lnu_BH(sky_galaxy,1.4e9)

   fluxBH = fluxDensity(LBH,snapshot(sam%snapshot)%redshift)

   sky_galaxy%fluxBH = scaledFluxDensity(fluxBH,frequency,frequency0,mu,sigma,snapshot(sam%snapshot)%redshift,.false.)

   LStars = Lnu_SFR(sky_galaxy%SFR)

   fluxStars = fluxDensity(LBH,snapshot(sam%snapshot)%redshift)

   sky_galaxy%fluxStars = scaledFluxDensity(fluxStars,frequency,frequency0,mu,sigma,snapshot(sam%snapshot)%redshift,.true.)
   
contains

   function M_Eddington(BHMass) result(rate_Eddington)

      implicit none
      real*8,intent(in)                       :: BHMass
      real*8                                  :: rate_Eddington
      real*8,parameter                        :: eScatteringOpacity=0.3
      real*8,parameter                        :: accretionEfficiency=0.1
      real*8                                  :: L_Eddington

      L_Eddington = 4 * pi * G * BHMass * c /  (eScatteringOpacity* 0.01**2 / 0.001) ! [g] [m]^4 / ([cm]^2 [s]3)

      rate_Eddington = L_Eddington / c**2 * (100 * 0.01)**2 * (0.001 / 0.001 ) / accretionEfficiency ! [kg]/[s]

   end function M_Eddington

   function Lnu_BH(sky_galaxy,nu) result(Lnu)
   
      implicit none
      type(type_sky_galaxy),intent(inout)    :: sky_galaxy
      real*4,intent(in)                      :: nu
      real*8                                 :: accretion_rate
      real*8                                 :: M_Eddinton_rate
      real*8                                 :: BHMass
      real*8                                 :: dm
      real*8                                 :: dt
      real*8                                 :: normalization
      real*8                                 :: Ljet_radiomode
      real*8                                 :: Ljet_quasar
      real*8                                 :: Lnu

      BHMass = sky_galaxy%BlackHoleMass * 1e10 * Msun ! kg

      dm = (sky_galaxy%BlackHoleAccretedHotMass+sky_galaxy%BlackHoleAccretedHotMass)*1e10 * Msun
      dt = sky_galaxy%dt*1e6*365*24*60*60
      accretion_rate = dm/dt
      accretion_rate = accretion_rate/M_Eddington(BHMass)

      if (accretion_rate>0) then

         if(accretion_rate<=0.01) then ! radio mode
            normalization = 8e-5

            Ljet_radiomode = 2.0_8 * 1E45_8 * BHMass / 10**9_8 / Msun * accretion_rate / 0.01_8 * sky_galaxy%Spin**2

            Lnu = normalization * (BHMass / 10**9 / Msun * accretion_rate / 0.01)**0.42 * Ljet_radiomode / nu

         else ! quasar mode

            normalization = 5e-2

            Ljet_quasar = 2.5 * 1E43_8 * (BHMass / 10**9 / Msun)**1.1 * (accretion_rate / 0.01)**1.2 * sky_galaxy%Spin**2

            Lnu = normalization * (BHMass / 10**9 / Msun)**0.32 * (accretion_rate / 0.01)**(-1.2) * Ljet_quasar / nu

         end if

      end if
   
   end function Lnu_BH


   function Lnu_SFR(SFR) result(Lnu)

      implicit none
      real*8,intent(in)                      :: SFR
      real*8                                 :: Lnu

      Lnu =  SFR / (0.75 * 10**(-21)) ! [erg] /[s] /[Hz]

   end function Lnu_SFR

   function fluxDensity(Luminosity,redshift) result(flux)

      implicit none
      real*4,intent(in)                      :: Luminosity
      real*4,intent(in)                      :: redshift
      real*8                                 :: distance
      real*8                                 :: prefactor
      real*8                                 :: flux

      distance = redshift_to_dc(redshift) * Mpc * 100 ! [m]

      prefactor = 10e26_8 / (4 * pi *  distance**2 )

      flux =  prefactor * Luminosity ! [erg] /[s] /[cm]^2 /[Hz] or [Jy]

   end function fluxDensity

   function frequencyRedshift(redshift,basefrequency) result(frequency)
      implicit none
      real*4,intent(in)                      :: redshift
      real*4,intent(in)                      :: basefrequency
      real*4                                 :: frequency

      frequency = basefrequency / (redshift + 1)

   end function frequencyRedshift

   function scaledFluxDensity(fluxDensity,frequency,frequency0,mu,sigma,redshift,redshiftFreq) result(scaledFlux)

      implicit none
      real*4,intent(in)                      :: fluxDensity
      real*4,intent(in)                      :: frequency
      real*4,intent(inout)                   :: frequency0
      real*4,intent(in)                      :: mu
      real*4,intent(in)                      :: sigma
      real*4,intent(in)                      :: redshift
      logical,intent(in)                     :: redshiftFreq
      real*4                                 :: gamma
      real*4                                 :: scaledFlux

      gamma = get_normal_random_number(mu,sigma)

      if (redshiftFreq) then
         frequency0 = frequencyRedshift(redshift,frequency0)
      end if

      scaledFlux = fluxDensity * (frequency/frequency0)**(-gamma)

   end function scaledFluxDensity
   
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
   para%subvolume_max = ncores-1
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

   call hdf5_read_data(g,sam%StellarMass,field_name='StellarMass')
   call hdf5_read_data(g,sam%Sfr,field_name='Sfr')
   call hdf5_read_data(g,sam%BlackHoleMass,field_name='BlackHoleMass')
   call hdf5_read_data(g,sam%BlackHoleAccretedHotMass,field_name='BlackHoleAccretedHotMass')
   call hdf5_read_data(g,sam%dt,field_name='dt')
   call hdf5_read_data(g,sam%Spin,field_name='Spin')
   call hdf5_read_data(g,sam%mvir,field_name='Mvir')
   
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
   call hdf5_write_data(trim(name)//'/icore',sky_galaxy%core,'core index')
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
   call hdf5_write_data(trim(name)//'/vpec_x',sky_galaxy%vpec(1),'[proper km/s] x-component of peculiar velocity')
   call hdf5_write_data(trim(name)//'/vpec_y',sky_galaxy%vpec(2),'[proper km/s] y-component of peculiar velocity')
   call hdf5_write_data(trim(name)//'/vpec_z',sky_galaxy%vpec(3),'[proper km/s] z-component of peculiar velocity')
   call hdf5_write_data(trim(name)//'/vpec_r',sky_galaxy%vpecrad,'[proper km/s] line-of-sight peculiar velocity')
   ! call hdf5_write_data(trim(name)//'/fluxBH',sky_galaxy%fluxBH,'[ergs/s /Hz] BlackHole flux')
   ! call hdf5_write_data(trim(name)//'/fluxStars',sky_galaxy%fluxStars,'[ergs/s /Hz] Stars flux')
   call hdf5_write_data(trim(name)//'/stellarmass',sky_galaxy%stellarmass,'[Msun/h] stellar mass')
   call hdf5_write_data(trim(name)//'/Sfr',sky_galaxy%Sfr,'[Msun/h yr] Star formation rate')
   call hdf5_write_data(trim(name)//'/BlackHoleMass',sky_galaxy%BlackHoleMass,'[Msun/h] Black Hole Mass')
   call hdf5_write_data(trim(name)//'/BlackHoleAccretedHotMass',sky_galaxy%BlackHoleAccretedHotMass, &
   & '[Msun/h] Black hole accreted hot mass')
   call hdf5_write_data(trim(name)//'/dt',sky_galaxy%dt,'[yr] Time between snapshot')
   call hdf5_write_data(trim(name)//'/Spin',sky_galaxy%Spin,'Spin of the black home')
   call hdf5_write_data(trim(name)//'/Mvir',sky_galaxy%Mvir,'[Msun/h] subhalo mass')
   
   if (para%line_parameters==1) then
      
      test(0) = 0
   
   else
   
      test(0) = 0
      
   end if
   
   test(1) = n
   test(3) = sum(sky_galaxy%tile)
   test(4) = sum(sky_galaxy%dc)
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