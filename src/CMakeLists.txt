# Generate placeholder "ersem_version" module (versions all set to "unknown" - updated at build time if FABM_EMBED_VERSION is set)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/version.F90.in ${CMAKE_CURRENT_BINARY_DIR}/version.F90)

add_library(fabm_models_gpm OBJECT
            gpm_model_library.F90
            ${CMAKE_CURRENT_BINARY_DIR}/version.F90
            common.F90
            phytoplankton.F90
            zooplankton.F90     
           )

add_dependencies(fabm_models_gpm fabm_base)

if(FABM_EMBED_VERSION)
  # Dynamically generate "ersem_version" module with git commit/branch information at every build.
  add_custom_target(fabm_models_gpm_version
     ${CMAKE_COMMAND} -DINFILE=${CMAKE_CURRENT_SOURCE_DIR}/version.F90.in -DOUTFILE=${CMAKE_CURRENT_BINARY_DIR}/version.F90 -P ${GET_GIT_INFO_SCRIPT}
     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
     COMMENT "Retrieving description of last GPM commit..." 
     SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/version.F90.in ${GET_GIT_INFO_SCRIPT}
     VERBATIM
  )
  add_dependencies(fabm_models_gpm fabm_models_gpm_version)
endif()

