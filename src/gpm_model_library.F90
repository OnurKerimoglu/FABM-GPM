module gpm_model_library

   use fabm_types, only: type_base_model_factory, type_base_model

   use gpm_plankton
   use gpm_abio_pel
   use gpm_abio_sed
   use gpm_version
   ! Add new components modules here

   implicit none

   private

   type,extends(type_base_model_factory) :: type_factory
      contains
      procedure :: initialize
      procedure :: create
   end type

   type (type_factory),save,target,public :: gpm_model_factory

contains
   
   subroutine initialize(self)
      class (type_factory), intent(inout) :: self
      call self%register_version('GPM',git_commit_id//' ('//git_branch_name//' branch)')
   end subroutine initialize

   subroutine create(self,name,model)
      class (type_factory),intent(in) :: self
      character(*),        intent(in) :: name
      class (type_base_model),pointer :: model

      select case (name)
         case ('plankton'); allocate(type_gpm_plankton::model)
         case ('abio_pel'); allocate(type_gpm_abio_pel::model)
         case ('abio_sed'); allocate(type_gpm_abio_sed::model)
         ! Add new components models here
      end select
   end subroutine create

end module gpm_model_library
