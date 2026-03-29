!
! This module contains general purpose utilities and WRF wrappers because some
! applications require this module to run standalone. No physics here.
! Some are dependent on WRF indexing scheme. Some violate WRF conventions but these
! are not called from the WRF dependent code. Some are not called at all.
! 

module module_fr_fire_util

! various method selection parameters
! 1. add the parameter and its static default here
! optional:
! 2. add copy from config_flags in module_fr_fire_driver%%set_flags
! 3. add a line in Registry.EM to define the variable and set default value
! 4. add the variable and value in namelist.input


contains

end module module_fr_fire_util
