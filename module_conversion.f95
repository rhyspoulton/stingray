module module_conversion

! standard functions for the conversion from intrinsic to apparent galaxy properties

   use module_constants
   use module_linalg
   use module_system
   use module_cosmology
   
   public
   
   real*4   :: Rvector(3,3), Rpseudo(3,3) ! argument used by function rotate, defined as global variable to
                                          ! avoid passing them through the use module

contains

function rotate(x,pseudovector) result(y)

   implicit none
   real*4,intent(in)    :: x(3)
   logical,intent(in)   :: pseudovector
   real*4               :: y(3)
   
   if (pseudovector) then
      y = matmul(Rpseudo,x)
   else
      y = matmul(Rvector,x)
   end if
   
end function rotate

subroutine make_redshift(x,v,z,zcos)

   implicit none
   real*4,intent(in)             :: x(3)     ! [box side-length] position vector
   real*4,intent(in)             :: v(3)     ! [km/s] velocity of galaxy
   real*4,intent(out),optional   :: z        ! total redshift
   real*4,intent(out),optional   :: zcos     ! cosmological redshift due to Hubble flow without accounting for relative velocity
   real*4                        :: zpec     ! redshift due to the peculiar motion of object relative to Hubble flow
   real*4                        :: zobs     ! redshift due to the peculiar motion of the observer relative to the Hubble flow
   real*4                        :: elos(3)  ! unit vector pointing along the line of sight
   real*4                        :: dc       ! [box side-length] comoving distance
   real*4                        :: z_,zcos_
   
   dc = norm(x)
   if (dc<=epsilon(dc)) call error('make_distances: norm of x is zero.')
   
   elos = x/dc ! unit-vector along the line of slight
   zcos_ = dc_to_redshift(dc*para%L*(para%length_unit/Mpc))
   zpec = min(0.1,max(-0.1,sum(v*elos)/c*1e3)) ! limited to 0.1 to avoid relativistic regime, ok for all practical purposes
   zobs = min(0.1,max(-0.1,-sum(para%velocity*elos)/c*1e3)) ! limited to 0.1 to avoid relativistic regime, ok for all practical purposes
   z_ = (1+zcos_)*(1+zpec)*(1+zobs)-1 ! following Davis & Scrimgeour 2014
   
   if (present(z))      z     = z_
   if (present(zcos))   zcos  = zcos_
   
end subroutine make_redshift

subroutine make_inclination_and_pa(x,J,inclination,pa)

   implicit none
   real*4,intent(in)                :: x(3)           ! [arbitrary units] position in cone
   real*4,intent(in)                :: J(3)           ! [arbitrary units] angular momentum in cone
   real*4,intent(out),optional      :: inclination    ! [rad] inclination (angle between LOS and galaxy axis)
   real*4,intent(out),optional      :: pa             ! [rad] position angle from North to East
   real*4                           :: eLOS(3)        ! unit vector along x = line-of-sight (LOS)
   real*4                           :: eJ(3)          ! unit vector pointing along galaxy axis
   real*4                           :: eMajor(3)      ! unit vector pointing along the major axis (orthogonal to LOS)
   real*4                           :: eNorth(3)      ! unit vector pointing north (or south) (orthoconal to LOS)
   real*4                           :: eEast(3)       ! unit vector pointing east (or west) (orthoconal to LOS)
   real*4                           :: normx,normJ
   
   normx = norm(x)
   if (normx<=epsilon(normx)) call error('make_inclination_and_pa: norm of x is zero.')
   normJ = norm(J)
   if (normJ<=epsilon(normJ)) call error('make_inclination_and_pa: norm of J is zero.')
   
   eLOS = x/normx
   eJ = J/normJ
   
   inclination = acos(sum(eLOS*eJ))
   if (inclination>pi) inclination = 2*pi-inclination
   
   eMajor = cross_product(eLOS,eJ)
   eMajor = eMajor/norm(eMajor)
   
   eNorth = (/0,0,1/)-eLOS(3)*eLOS
   eNorth = eNorth/norm(eNorth)
   eEast = cross_product(eNorth,eLOS)
   
   pa = atan2(sum(eMajor*eEast),sum(eMajor*eNorth))
   
end subroutine make_inclination_and_pa

subroutine make_sky_coordinates(x,dc,ra,dec)

   implicit none
   real*4,intent(in)                :: x(3)  ! [box side-length] position vector
   real*4,intent(out),optional      :: dc    ! [length units of simulation] comoving distance
   real*4,intent(out),optional      :: ra    ! [rad] right ascension
   real*4,intent(out),optional      :: dec   ! [rad] declination
   real*4                           :: normx
   
   normx = norm(x)
   if (normx<=epsilon(normx)) call error('make_sky_coordinates: norm of x is zero.')
   
   if (present(dc))  dc  = normx*para%L
   if (present(ra))  ra  = atan2(x(1),x(3))
   if (present(dec)) dec = asin(x(2)/normx)
   
end subroutine make_sky_coordinates

function convert_stellarmass2absmag(M,M2L) result(mag)
      
   implicit none
   real*4,intent(in)    :: M     ! [Msun] stellar mass
   real*4,intent(in)    :: M2L   ! [Msun/Lsun] stellar mass-to-light ratio
   real*4               :: mag   ! absolute magnitude
   
   mag = convert_luminosity2absmag(M/M2L)

end function convert_stellarmass2absmag

function convert_luminosity2absmag(L) result(mag)
      
   implicit none
   real*4,intent(in)    :: L              ! [Lsun] luminosity
   real*4               :: mag            ! absolute magnitude
   real*4,parameter     :: magSun = 4.83  ! absolute magnitude of sun
   
   mag = magSun-2.5*log10(L)

end function convert_luminosity2absmag

function convert_luminosity2flux(L,dl) result(S)

   implicit none
   real*8,intent(in) :: L     ! [W] Luminosity
   real*4,intent(in) :: dl    ! [Mpc] luminosity distance
   real*8            :: S     ! [W/m^2] Flux

   S = L/real(dl,8)**2/ASphereMpc

end function convert_luminosity2flux

function convert_absmag2appmag(absmag,dl) result(appmag)

   implicit none
   real*4,intent(in) :: absmag   ! absolute magnitude
   real*4,intent(in) :: dl       ! [Mpc] luminosity distance
   real*4            :: appmag   ! apparent magnitude

   appmag = absmag+5*log10(dl)+25

end function convert_absmag2appmag

function convert_vector(x,rotation) result(y)

   implicit none
   real*4,intent(in)    :: x(3)
   integer*4,intent(in) :: rotation
   real*4               :: y(3)

   select case (abs(rotation))
      case (1) ! identity
         y = x
      case(2) ! invert x-axis, while permuting y and z
         y = (/-x(1),x(3),x(2)/)
      case(3) ! invert y-axis, while permuting z and x
         y = (/x(3),-x(2),x(1)/)
      case(4) ! invert z-axis, while permuting x and y
         y = (/x(2),x(1),-x(3)/)
      case(5) ! permute (x,y,z) -> (y,z,x)
         y = (/x(2),x(3),x(1)/)
      case(6) ! permute (x,y,z) -> (z,x,y)
         y = (/x(3),x(1),x(2)/)
      case default
         call error('Unknown rotation')
   end select

   if (rotation<0) y = -y ! inversion

end function convert_vector

function convert_pseudovector(x,rotation) result(y)

   implicit none
   real*4,intent(in)    :: x(3)
   integer*4,intent(in) :: rotation
   real*4               :: y(3)

   y = convert_vector(x,abs(rotation))

end function convert_pseudovector

end module module_conversion