#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: gpm_abio_sed --- GPM basic sediment abiotic model: DI(N,P,Si),det(C,N,P,Si)
!
! !INTERFACE:
   module gpm_abio_sed
!
! !DESCRIPTION:
! This is a very simple model for remineralization of pom to dim in sediment
!
! !USES:
   use fabm_types
!   use fabm_driver

   implicit none

!  default: all is private.
   private
!!
! !PUBLIC TYPES:
   type,extends(type_base_model),public    ::   type_gpm_abio_sed
!     Variable identifiers
      type (type_bottom_state_variable_id) :: id_detC_sed,id_detP_sed,id_DIP_sed,id_detN_sed,id_DIN_sed,id_detSi_sed,id_DISi_sed
      type (type_state_variable_id)        :: id_detC_pel,id_detP_pel,id_DIP_pel,id_detN_pel,id_DIN_pel,id_detSi_pel,id_DISi_pel
      type (type_horizontal_diagnostic_variable_id)   ::id_dsedrC,id_ddecC,id_drhsdetC
      type (type_horizontal_diagnostic_variable_id)   ::id_detQP,id_drhsDIP,id_drhsdetP,id_dDIP_diff,id_dsedrP,id_ddecP
      type (type_horizontal_diagnostic_variable_id)   ::id_detQN,id_drhsDIN,id_drhsdetN,id_dDIN_diff,id_dsedrN,id_ddecN
      type (type_horizontal_diagnostic_variable_id)   ::id_detQSi,id_drhsDISi,id_drhsdetSi,id_dDISi_diff,id_dsedrSi,id_ddecSi
      type (type_dependency_id)         :: id_temp
      real(rk) :: det_wsed,thick_sed,DIM_diffcoef,det_decr,Q10,Tref
      logical  :: coupled2water
      logical  :: resolve_P,resolve_N,resolve_Si
      real(rk) :: const_detC,const_detP,const_DIP,const_detN,const_DIN,const_detSi,const_DISi
      
      !     Model procedures
      contains
      procedure :: initialize
      procedure :: do_bottom      
   end type type_gpm_abio_sed
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the sediment model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the parameters and variables exported by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_gpm_abio_sed), intent(inout), target :: self
   integer,                          intent(in)            :: configunit
!
! !LOCAL VARIABLES:
   real(rk), parameter :: secs_pr_day = 86400.
                                    
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   call self%get_parameter(self%resolve_P,     'resolve_P',    '-',       'whether to resolve P cycle',          default=.true.) 
   call self%get_parameter(self%resolve_N,     'resolve_N',    '-',       'whether to resolve N cycle',          default=.false.) 
   call self%get_parameter(self%resolve_Si,    'resolve_Si',    '-',      'whether to resolve Si cycle',         default=.false.)
   call self%get_parameter(self%det_wsed,      'det_wsed',     'm/d',     'sedimentation rate',                  default= 0.5_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%thick_sed,     'thick_sed',    'm',       'thickness of the sediment layer',     default=0.2_rk)
   call self%get_parameter(self%DIM_diffcoef,  'DIM_diffcoef', 'm^2/d',    'DIM dif. coeff. at w-s interface',    default=5e-4_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%det_decr,      'det_decr',     'd-1',     'decay rate constant',                 default=0.1_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%Q10,           'Q10',          '-',       'Q10 for sediment processes',          default=2.0_rk)    
   call self%get_parameter(self%Tref,          'Tref',         'celcius', 'reference temperature for bac. proc.',default=20.0_rk)
   call self%get_parameter(self%coupled2water, 'coupled2water','-',       'whether coupled with pelagic model',  default=.true.)    
   !C
   call self%get_parameter(self%const_detC,    'const_detC',   'mmol/m^3','constant pel. detC, if not coupled', default=116.0_rk)
   !P
   if (self%resolve_P) then
     call self%get_parameter(self%const_detP,  'const_detP',   'mmol/m^3','constant pel. detP, if not coupled', default=1.0_rk)
     call self%get_parameter(self%const_DIP,   'const_DIP',    'mmol/m^3','constant pel. DIP, if not coupled',  default=2.0_rk)
   end if
   !N
   if (self%resolve_N) then
     call self%get_parameter(self%const_detN,  'const_detN',   'mmol/m^3','constant pel. detN, if not coupled', default=14.0_rk)
     call self%get_parameter(self%const_DIN,   'const_DIN',    'mmol/m^3','constant pel. DIN, if not coupled',  default=28.0_rk)
   end if
   !Si
   if (self%resolve_Si) then
     call self%get_parameter(self%const_detSi,  'const_detSi',   'mmol/m^3','constant pel. detSi, if not coupled', default=14.0_rk)
     call self%get_parameter(self%const_DISi,   'const_DISi',    'mmol/m^3','constant pel. DISi, if not coupled',  default=28.0_rk)
   end if
   
   ! Register state variables
   !C
   call self%register_state_variable(self%id_detC_sed,'detC','mmol/m^2','detrital carbon in sediment', minimum=_ZERO_)
   !P
   if (self%resolve_P) then
     call self%register_state_variable(self%id_detP_sed,'detP','mmol/m^2','detrital phosphorus in sediment', minimum=_ZERO_)
     call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_detP_sed)
     call self%register_state_variable(self%id_DIP_sed,'DIP','mmol/m^2','dissolved inorganic phosphorus in sediment', minimum=_ZERO_)
     call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_DIP_sed)
   end if
   !N
   if (self%resolve_N) then
     call self%register_state_variable(self%id_detN_sed,'detN','mmol/m^2','detrital nitrogen in sediment', minimum=_ZERO_)
     call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_detN_sed)
     call self%register_state_variable(self%id_DIN_sed,'DIN','mmol/m^2','dissolved inorganic nitrogen in sediment', minimum=_ZERO_)
     call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_DIN_sed)
   end if
   !Si
   if (self%resolve_Si) then
     call self%register_state_variable(self%id_detSi_sed,'detSi','mmol/m^2','detrital silica in sediment', minimum=_ZERO_)
     call self%add_to_aggregate_variable(standard_variables%total_silicate,self%id_detSi_sed)
     call self%register_state_variable(self%id_DISi_sed,'DISi','mmol/m^2','dissolved inorganic silica in sediment', minimum=_ZERO_)
     call self%add_to_aggregate_variable(standard_variables%total_silicate,self%id_DISi_sed)
   end if
   
   !Regester state dependencies
   
   
   if (self%coupled2water) then
    !C
    call self%register_state_dependency(self%id_detC_pel, 'pelagic_detC') !,   'mmol m-3','detrital C variable in pelagic')
    !P
    if (self%resolve_P) then
      call self%register_state_dependency(self%id_detP_pel, 'pelagic_detP') !,   'mmol m-3','detrital P variable in pelagic')
      call self%register_state_dependency(self%id_DIP_pel, 'pelagic_DIP') !,   'mmol m-3','DIP variable in pelagic')
    end if
    !N
    if (self%resolve_N) then
      call self%register_state_dependency(self%id_detN_pel, 'pelagic_detN') !,   'mmol m-3','detrital N variable in pelagic')
      call self%register_state_dependency(self%id_DIN_pel, 'pelagic_DIN') !,   'mmol m-3','DIN variable in pelagic')
    end if
    !Si
    if (self%resolve_Si) then
      call self%register_state_dependency(self%id_detSi_pel, 'pelagic_detSi') !,   'mmol m-3','detrital Si variable in pelagic')
      call self%register_state_dependency(self%id_DISi_pel, 'pelagic_DISi') !,   'mmol m-3','DISi variable in pelagic')
    end if
   end if
   
   ! Register diagnostic variables
   !C
   call self%register_diagnostic_variable(self%id_dsedrC,'detCsed','mmol/m^2/d',  'detC sedimentation flux',             &
                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_ddecC,'detCrem','mmol/m^2/d',  'detC decay rate',                 &
                     output=output_instantaneous) !output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_drhsdetC,'ddetC','mmol/m^2/d',  'detC_RHS',                          &
                     output=output_time_step_averaged)
   !P
   if (self%resolve_P) then
     call self%register_diagnostic_variable(self%id_detQP,'detQP','molP/molC',  'molar P:C ratio of detritus',             &
                     output=output_time_step_averaged)
     call self%register_diagnostic_variable(self%id_dsedrP,'detPsed','mmol/m^2/d',  'detP sedimentation flux',             &
                     output=output_time_step_averaged)
     call self%register_diagnostic_variable(self%id_dDIP_diff,'DIPdif','mmol/m^2/d',  'DIP diffusive flux',                  &
                     output=output_time_step_averaged)
     call self%register_diagnostic_variable(self%id_ddecP,'detPrem','mmol/m^2/d',  'detP decay rate',                 &
                     output=output_instantaneous) !output_time_step_averaged)
     call self%register_diagnostic_variable(self%id_drhsdetP,'ddetP','mmol/m^2/d',  'detP_RHS',                          &
                     output=output_time_step_averaged)
     call self%register_diagnostic_variable(self%id_drhsDIP,'dDIP','mmol/m^2/d',  'DIP_RHS',                          &
                     output=output_time_step_averaged)
   end if
   !N
   if (self%resolve_N) then
     call self%register_diagnostic_variable(self%id_detQN,'detQN','molN/molC',  'molar N:C ratio of detritus',             &
                     output=output_time_step_averaged)
     call self%register_diagnostic_variable(self%id_dsedrN,'detNsed','mmol/m^2/d',  'detN sedimentation flux',             &
                     output=output_time_step_averaged)
     call self%register_diagnostic_variable(self%id_dDIN_diff,'DINdif','mmol/m^2/d',  'DIN diffusive flux',                  &
                     output=output_time_step_averaged)
     call self%register_diagnostic_variable(self%id_ddecN,'detNrem','mmol/m^2/d',  'detN decay rate',                 &
                     output=output_instantaneous) !output_time_step_averaged)
     call self%register_diagnostic_variable(self%id_drhsdetN,'ddetN','mmol/m^2/d',  'detN_RHS',                          &
                     output=output_time_step_averaged)
     call self%register_diagnostic_variable(self%id_drhsDIN,'dDIN','mmol/m^2/d',  'DIN_RHS',                          &
                     output=output_time_step_averaged)
   end if
   !Si
   if (self%resolve_Si) then
     call self%register_diagnostic_variable(self%id_detQSi,'detQSi','molSi/molC',  'molar Si:C ratio of detritus',             &
                     output=output_time_step_averaged)
     call self%register_diagnostic_variable(self%id_dsedrSi,'detSised','mmol/m^2/d',  'detSi sedimentation flux',             &
                     output=output_time_step_averaged)
     call self%register_diagnostic_variable(self%id_dDISi_diff,'DISidif','mmol/m^2/d',  'DISi diffusive flux',                  &
                     output=output_time_step_averaged)
     call self%register_diagnostic_variable(self%id_ddecSi,'detSirem','mmol/m^2/d',  'detSi decay rate',                 &
                     output=output_instantaneous) !output_time_step_averaged)
     call self%register_diagnostic_variable(self%id_drhsdetSi,'ddetSi','mmol/m^2/d',  'detSi_RHS',                          &
                     output=output_time_step_averaged)
     call self%register_diagnostic_variable(self%id_drhsDISi,'dDISi','mmol/m^2/d',  'DISi_RHS',                          &
                     output=output_time_step_averaged)
   end if

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of the sediment model
!
! !INTERFACE:
   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
!
! !DESCRIPTION:
! This routine calculates the benthic remineralization and
! (matching) bottom fluxes for pelagic variables. 
! All flux and change terms have units mmol/m**2/s.
!
! !INPUT PARAMETERS:
   class (type_gpm_abio_sed),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
! !LOCAL VARIABLES:
   real(rk)                   :: fT,temp,spdecr
   real(rk)                   :: detC_sed,detC_pel
   real(rk)                   :: detC_sedr,detC_decr,ddetC
   real(rk)                   :: DIP_sed,detP_sed,DIP_pel,detP_pel
   real(rk)                   :: detP_sedr,detP_decr,ddetP,dDIP,DIP_loss2pel
   real(rk)                   :: DIN_sed,detN_sed,DIN_pel,detN_pel
   real(rk)                   :: detN_sedr,detN_decr,ddetN,dDIN,DIN_loss2pel
   real(rk)                   :: DISi_sed,detSi_sed,DISi_pel,detSi_pel
   real(rk)                   :: detSi_sedr,detSi_decr,ddetSi,dDISi,DISi_loss2pel
   real(rk), parameter        :: secs_pr_day = 86400.
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops over the horizontal domain (if any).
   _HORIZONTAL_LOOP_BEGIN_
    
   ! Retrieve current environmental conditions.
   _GET_(self%id_temp,temp)  ! temperature

   ! Retrieve current (local) state variable values.
   if (self%coupled2water) then
     !C
     _GET_(self%id_detC_pel,detC_pel)    ! det-C concentration - pelagic
     !P
     if (self%resolve_P) then
       _GET_(self%id_detP_pel,detP_pel)    ! det-P concentration - pelagic
       _GET_(self%id_DIP_pel,DIP_pel)      ! DIP concentration - pelagic
     end if
     !N
     if (self%resolve_N) then
       _GET_(self%id_detN_pel,detN_pel)    ! det-N concentration - pelagic
       _GET_(self%id_DIN_pel,DIN_pel)      ! DIN concentration - pelagic
     end if
     !Si
     if (self%resolve_Si) then
       _GET_(self%id_detSi_pel,detSi_pel)    ! det-Si concentration - pelagic
       _GET_(self%id_DISi_pel,DISi_pel)      ! DISi concentration - pelagic
     end if
   else
     !C
     detC_pel = self%const_detC          ! no coupling - constant pelagic detC
     !P
     if (self%resolve_P) then
       detP_pel = self%const_detP          ! no coupling - constant pelagic detP
       DIP_pel = self%const_DIP            ! no coupling - constant pelagic DIP
     end if
     !N
     if (self%resolve_N) then
       detN_pel = self%const_detN          ! no coupling - constant pelagic detN
       DIN_pel = self%const_DIN            ! no coupling - constant pelagic DIN
     end if
     !Si
     if (self%resolve_Si) then
       detSi_pel = self%const_detSi          ! no coupling - constant pelagic detSi
       DISi_pel = self%const_DISi            ! no coupling - constant pelagic DISi
     end if
   end if
   
   !C
   _GET_HORIZONTAL_(self%id_detC_sed,detC_sed) ! det-C density - sediment
   !P
   if (self%resolve_P) then
     _GET_HORIZONTAL_(self%id_detP_sed,detP_sed) ! det-P density - sediment
     _GET_HORIZONTAL_(self%id_DIP_sed,DIP_sed)  ! DIP density - sediment
   end if 
   !N
   if (self%resolve_N) then
     _GET_HORIZONTAL_(self%id_detN_sed,detN_sed) ! det-N density - sediment
     _GET_HORIZONTAL_(self%id_DIN_sed,DIN_sed)  ! DIN density - sediment
   end if
   !Si
   if (self%resolve_Si) then
     _GET_HORIZONTAL_(self%id_detSi_sed,detSi_sed) ! det-Si density - sediment
     _GET_HORIZONTAL_(self%id_DISi_sed,DISi_sed)  ! DISi density - sediment
   end if
   
   ! Decay rate of sediment detritus
   fT=self%Q10**((temp-self%Tref)/10)
   spdecr = fT*self%det_decr  
   !C
   detC_decr = spdecr*detC_sed
   !P
   if (self%resolve_P) then
     detP_decr = spdecr*detP_sed
   end if
   !N
   if (self%resolve_N) then
     detN_decr = spdecr*detN_sed
   end if
   !Si
   if (self%resolve_Si) then
     detSi_decr = spdecr*detSi_sed
   end if
   
   ! sedimentation rate 
   !C
   detC_sedr = self%det_wsed*detC_pel
   !P
   if (self%resolve_P) then
     detP_sedr = self%det_wsed*detP_pel
   end if
   !N
   if (self%resolve_N) then
     detN_sedr = self%det_wsed*detN_pel
   end if
   !Si
   if (self%resolve_Si) then
     detSi_sedr = self%det_wsed*detSi_pel
   end if
   
   ! loss rate of DIM to pelagic. First the sediment variable (DIM_sed. areal units) has to be converted to concentration using the thickness of the sediment (thick_sed) for being able to calculate the gradient. Then this gradient is assumed to be taking place within a boundary layer, thickness of which is equal to the thickness of the sediment layer itself:
   !C: not calculated
   !P
   if (self%resolve_P) then
     DIP_loss2pel = self%DIM_diffcoef*(DIP_sed/self%thick_sed-DIP_pel)/self%thick_sed
     ! m2/d * mmol/m3 *1/m = mmol/m2/d
   end if
   !N
   if (self%resolve_N) then
     DIN_loss2pel = self%DIM_diffcoef*(DIN_sed/self%thick_sed-DIN_pel)/self%thick_sed
     ! m2/d * mmol/m3 *1/m = mmol/m2/d
   end if
   !Si
   if (self%resolve_Si) then
     DISi_loss2pel = self%DIM_diffcoef*(DISi_sed/self%thick_sed-DISi_pel)/self%thick_sed
     ! m2/d * mmol/m3 *1/m = mmol/m2/d
   end if
   
   ! rhs
   !C
   ddetC = detC_sedr - detC_decr
   !P
   if (self%resolve_P) then
     ddetP = detP_sedr - detP_decr
     dDIP = detP_decr - DIP_loss2pel 
     !write(*,*),'r,d:',remin*detP_sed,-1.0*DIM_diffusion 
   end if
   !N
   if (self%resolve_N) then
     ddetN = detN_sedr - detN_decr
     dDIN = detN_decr - DIN_loss2pel 
     !write(*,*),'r,d:',remin*detN_sed,-1.0*DIM_diffusion 
   end if
   !Si
   if (self%resolve_Si) then
     ddetSi = detSi_sedr - detSi_decr
     dDISi = detSi_decr - DISI_loss2pel 
     !write(*,*),'r,d:',remin*detSi_sed,-1.0*DIM_diffusion 
   end if
   
   ! Set local temporal derivatives of sediment variables
   !C
   _SET_ODE_BEN_(self%id_detC_sed,ddetC)
   !P
   if (self%resolve_P) then
     _SET_ODE_BEN_(self%id_detP_sed,ddetP)
     _SET_ODE_BEN_(self%id_DIP_sed,dDIP)
   end if
   !N
   if (self%resolve_N) then
     _SET_ODE_BEN_(self%id_detN_sed,ddetN)
     _SET_ODE_BEN_(self%id_DIN_sed,dDIN)
   end if
   !Si
   if (self%resolve_Si) then
     _SET_ODE_BEN_(self%id_detSi_sed,ddetSi)
     _SET_ODE_BEN_(self%id_DISi_sed,dDISi)
   end if
   
   ! Set bottom fluxes of pelagic variables (these mirror local sediment derivatives), if coupled
   if (self%coupled2water) then
     !C
     _SET_BOTTOM_EXCHANGE_(self%id_detC_pel,-detC_sedr)
     !P
     if (self%resolve_P) then
       _SET_BOTTOM_EXCHANGE_(self%id_detP_pel,-detP_sedr)
       _SET_BOTTOM_EXCHANGE_(self%id_DIP_pel,DIP_loss2pel)
     end if
     !N
     if (self%resolve_N) then
       _SET_BOTTOM_EXCHANGE_(self%id_detN_pel,-detN_sedr)
       _SET_BOTTOM_EXCHANGE_(self%id_DIN_pel,DIN_loss2pel)
     end if
     !Si
     if (self%resolve_Si) then
       _SET_BOTTOM_EXCHANGE_(self%id_detSi_pel,-detSi_sedr)
       _SET_BOTTOM_EXCHANGE_(self%id_DISI_pel,DISi_loss2pel)
     end if
   end if 

   ! Export diagnostic variables
   !C
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_dsedrC,detC_sedr*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ddecC,detC_decr*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_drhsdetC, ddetC*secs_pr_day)
   !P
   if (self%resolve_P) then
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_detQP,detP_sed/detC_sed)
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_dsedrP,detP_sedr*secs_pr_day)
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_dDIP_diff,DIP_loss2pel*secs_pr_day)
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ddecP,detP_decr*secs_pr_day)
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_drhsdetP, ddetP*secs_pr_day)
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_drhsDIP, dDIP*secs_pr_day)
   end if
   !N
   if (self%resolve_N) then
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_detQN,detN_sed/detC_sed)
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_dsedrN,detN_sedr*secs_pr_day)
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_dDIN_diff,DIN_loss2pel*secs_pr_day)
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ddecN,detN_decr*secs_pr_day)
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_drhsdetN, ddetN*secs_pr_day)
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_drhsDIN, dDIN*secs_pr_day)
   end if
   !Si
   if (self%resolve_Si) then
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_detQSi,detSi_sed/detSi_sed)
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_dsedrSi,detSi_sedr*secs_pr_day)
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_dDISi_diff,DISi_loss2pel*secs_pr_day)
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ddecSi,detSi_decr*secs_pr_day)
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_drhsdetSi, ddetSi*secs_pr_day)
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_drhsDISi, dDISi*secs_pr_day)
   end if
   
   ! Leave spatial loops over the horizontal domain (if any).
   _HORIZONTAL_LOOP_END_

   end subroutine do_bottom
!EOC

!-----------------------------------------------------------------------

   end module gpm_abio_sed
