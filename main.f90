Program main

  use precision, only: dp
  use input, only: pi
  use input, only: source_information, partition_setup, layer_setup, setup_input
  use array, only: pyramid_creation, cartesian_creation, domain_creation, &
                   pyramid_destruction, cartiesian_destruction, domain_destruction
  use radiation, only: photo_table_setup  
  use pyramid_initialization, only: pyramid_volume, pyramid_cellsize, &
                                    pyramid_solid_angle, root_setup 
  use case_division, only: pyramid_case, cartesian_case
  use field, only: pyramid_field_generation, short_field_generation, &
                   long_field_generation, pyramid_analytical_field_generation, &
                   cartesian_analytical_field_generation
  use geometry_transformation, only: volume_partition
  use pyramid_global_source_transformation, only: pyramid_global_to_source_species_density_transformation, &
                                                  pyramid_source_to_global_photo_rate_transformation
  use pyramid_column_density, only: pyramid_domain_column_density_calculation
  use pyramid_photo_rate, only: pyramid_domain_photo_rate_calculation 
  use pyramid_source_domain_transformation, only: pyramid_source_to_domain_species_density_transformation, &
                                                  pyramid_domain_to_source_photo_rate_transformation
  use short_coordinate_transformation, only: short_global_to_source_transformation, &
                                             short_source_to_global_transformation
  use short, only: short_characteristic
  use short_column_density, only: short_column_density_calculation
  use short_photo_rate, only: short_photo_rate_calculation
  use long_coordinate_transformation, only: long_global_to_source_transformation, &
                                            long_source_to_global_transformation
  use long_column_density, only: long_column_density_calculation
  use long_photo_rate, only: long_photo_rate_calculation
  use pyramid_analytical_column_density, only: pyramid_analytical_column_density_calculation
  use pyramid_analytical_photo_rate, only: pyramid_analytical_photo_rate_calculation
  use cartesian_analytical_column_density, only: cartesian_analytical_column_density_calculation
  use cartesian_analytical_photo_rate, only: cartesian_analytical_photo_rate_calculation
  use output, only: output_global
  use trilinear, only: trilinear_interpolation

  implicit none

!!! Basic initialization !!!

  call setup_input ()  

  call source_information()
    
  call partition_setup()

  call layer_setup()
  
  call pyramid_creation()

  call cartesian_creation()
    
  call domain_creation()

  call photo_table_setup()

  call pyramid_volume()

  call pyramid_cellsize()
  
  call pyramid_solid_angle()

  call root_setup()
  
  call pyramid_case()
  
  call cartesian_case()
  
  call volume_partition()
  
  !!! pyramidal ray tracing !!!
  call pyramid_field_generation()
    
  call pyramid_global_to_source_species_density_transformation() 
   
  call pyramid_source_to_domain_species_density_transformation()

  call pyramid_domain_column_density_calculation()

  call pyramid_domain_photo_rate_calculation()
  
  call pyramid_domain_to_source_photo_rate_transformation() 
  
  call pyramid_source_to_global_photo_rate_transformation()
    
  !!! short ray tracing !!!  
  call short_field_generation()

  call short_global_to_source_transformation()
      
  call short_column_density_calculation()
  
  call short_photo_rate_calculation() 
   
  call short_source_to_global_transformation()
  
  !!! long ray tracing !!!
  call long_field_generation()

  call long_global_to_source_transformation()
  
  call long_column_density_calculation()

  call long_photo_rate_calculation()
     
  call long_source_to_global_transformation()
   
  !!! trilinear interpolation from (2N+1)^3 to (2N)^3
  call trilinear_interpolation ()

  !!! pyramid analytical solution
  call pyramid_analytical_field_generation()

  call pyramid_analytical_column_density_calculation()

  call pyramid_analytical_photo_rate_calculation()

  !!! cartesian analytical solution
  call cartesian_analytical_field_generation()

  call cartesian_analytical_column_density_calculation()

  call cartesian_analytical_photo_rate_calculation()

  !!!
  call output_global()
  
  call pyramid_destruction()
  
  call cartiesian_destruction()

  call domain_destruction()

end Program main
