
! SFIRE - Spread fire model in WRF-Fire
!
! This module is the entry point from WRF ARW to the wildland 
! fire module. The call to fire_driver advances the fire module by 
! one timestep. The fire module inputs the wind and outputs 
! temperature and humidity tendencies. The fire module also inputs a 
! number of constant arrays (fuel data, topography). Additional 
! arguments are model state (for data assimilation) and constant arrays 
! the model gives to WRF for safekeeping because it is not allowed 
! to save anything.
!
! Contributions to this wildland fire module have come from individuals at
! NCAR, the U.S.D.A. Forest Service, the Australian Bureau of Meteorology, 
! and the University of Colorado at Denver. 
!

module module_fr_fire_driver
! use this module for standalone call, you only need to provide some mock-up wrf modules  


USE module_domain, only: domain
USE module_configure, only: grid_config_rec_type
use module_model_constants, only: reradius,    & ! 1/earth radiusw
                                  g,           & ! gravitation acceleration
                                  pi2            ! 2*pi

implicit none



contains

 

end module module_fr_fire_driver
