#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: gpm_abio_sed_EH --- GPM sediment abiotic model based on ECOHAMv5: DI(N,P,Si),det(C,N,P,Si)
!
! !INTERFACE:
   module gpm_abio_sed_EH
!
! !DESCRIPTION:
! This is a very simple model for remineralization of pom to dim in sediment
!
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
   
   real(rk), parameter :: pi=acos(-1._rk)
   real(rk), parameter :: s2d = 86400._rk
   real(rk), parameter :: Dsw= 5e-4_rk/s2d ![m2/s] diffusion coefficient at the sediment/water interface
   real(rk), parameter :: Dasl=0.1_rk ![m] depth of the active sediment layer
   
   ! !PUBLIC DERIVED TYPES:
   
   !elements resolved
   type,public :: type_elms 
     real(rk) :: C
     real(rk) :: P     
     real(rk) :: N
     real(rk) :: Si
     !real(rk) :: Ccal
   end type
   
   ! dissolved inorganic material species resolved
   type, public :: type_dim 
     real(rk) :: C
     real(rk) :: P
     real(rk) :: NO3
     real(rk) :: NH4
     real(rk) :: Si
     !real(rk) :: Ccal
   end type

   ! gpm_abio_sed_EH 
   type,extends(type_base_model),public    ::   type_gpm_abio_sed_EH
!     Variable identifiers
      type (type_bottom_state_variable_id) :: id_POMsed_C,id_POMsed_P,id_POMsed_N,id_POMsed_Si
      type (type_bottom_state_variable_id) :: id_DIMsed_C,id_DIMsed_P,id_DIMsed_N,id_DIMsed_Si
      type (type_state_variable_id)        :: id_detpel1_C,id_detpel1_P,id_detpel1_N !slow sinking det in pel
      type (type_state_variable_id)        :: id_detpel2_C,id_detpel2_P,id_detpel2_N,id_detpel2_Si !fast sinking det in pel
      type (type_state_variable_id)        :: id_O2_pel,id_DIMpel_C,id_DIMpel_P,id_DIMpel_NO3,id_DIMpel_NH4,id_DIMpel_Si
      type (type_horizontal_diagnostic_variable_id)   :: id_det1_sed_C,id_det2_sed_C,id_sed_wdim_C,id_sed_rem_C
      type (type_horizontal_diagnostic_variable_id)   :: id_det1_sed_P,id_det2_sed_P,id_sed_wdim_P,id_sed_rem_P,id_POM_QP
      type (type_horizontal_diagnostic_variable_id)   :: id_det1_sed_N,id_det2_sed_N,id_sed_wdim_N,id_sed_rem_N,id_POM_QN
      type (type_horizontal_diagnostic_variable_id)   :: id_det2_sed_Si,id_sed_wdim_Si,id_sed_rem_Si,id_POM_QSi
      type (type_horizontal_diagnostic_variable_id)   :: id_o2o_brm,id_sed_o3c,id_sed_nn2
      type (type_horizontal_diagnostic_variable_id)   :: id_fsorbed,id_sedO2
      type (type_dependency_id)         :: id_temp
      real(rk) :: det1_wsed,det2_wsed,Q10,Tref
      real(rk) :: brc,brp,brn,brsi
      real(rk) :: do_sorpeq,K_T2do
      integer :: sorpmeth,sedO2meth
      logical  :: coupled2water
      logical  :: resolve_Si,resolve_DIC,resolve_sedDIM
      real(rk) :: const_det1C,const_det1P,const_det1N
      real(rk) :: const_det2C,const_det2P,const_det2N,const_det2Si
      real(rk) :: const_pelC,const_pelP,const_pelNH4,const_pelNO3,const_pelSi
      real(rk) :: const_O2,const_NO3!const_DIC,const_DIP,const_NH4,const_DISi

      !     Model procedures
      contains
      procedure :: initialize
      procedure :: do_bottom      
   end type type_gpm_abio_sed_EH
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
   class (type_gpm_abio_sed_EH), intent(inout), target :: self
   integer,                          intent(in)            :: configunit
!
! !LOCAL VARIABLES:
   real(rk), parameter :: s2d = 86400.
                                    
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   call self%get_parameter(self%resolve_Si,    'resolve_Si',    '-',      'whether to resolve Si cycle',         default=.false.)
   call self%get_parameter(self%resolve_DIC,  'resolve_DIC',  '-',         'whether to resolve DIC',          default=.false.)
   call self%get_parameter(self%resolve_sedDIM,  'resolve_sedDIM',  '-',   'whether to resolve sediment-DIM',  default=.false.)
   call self%get_parameter(self%sedO2meth, 'sedO2meth',   '-',          'method to estimate the O2 in sediment',  default=0)
   call self%get_parameter(self%sorpmeth, 'sorpmeth',   '-',  'method to calculate sorbed-fraction of sediment-P. 1: instantaneous, based on sed.O2 [only if resolve_sedDIM]',    default=0)
   call self%get_parameter(self%det1_wsed,      'det1_wsed',     'm/d',     'sedimentation rate of the slow sinking det',                  default= 0.5_rk, scale_factor=0.4_rk/s2d)
   call self%get_parameter(self%det2_wsed,      'det2_wsed',     'm/d',     'sedimentation rate of the fast sinking det',                  default= 0.5_rk, scale_factor=10.0_rk/s2d)
   !call self%get_parameter(self%DIM_diffcoef,  'DIM_diffcoef', 'm^2/d',    'DIM dif. coeff. at w-s interface',    default=5e-4_rk, scale_factor=1.0_rk/s2d)
   call self%get_parameter(self%brc,      'brc',     'd-1',     'benthic rem rate for C',                 default=0.028_rk, scale_factor=1.0_rk/s2d)
   call self%get_parameter(self%brp,      'brp',     'd-1',     'benthic rem rate for P',                 default=0.0333_rk, scale_factor=1.0_rk/s2d)
   call self%get_parameter(self%brn,      'brn',     'd-1',     'benthic rem rate for N',                 default=0.0333_rk, scale_factor=1.0_rk/s2d)
   call self%get_parameter(self%brsi,     'brsi',     'd-1',     'benthic rem rate for Si',                 default=0.0130_rk, scale_factor=1.0_rk/s2d)
   call self%get_parameter(self%K_T2do,    'K_T2do',   'mmolO2/C,', '[if sedO2meth=1] constant for estimating DO from T',          default=20.0_rk)
   call self%get_parameter(self%do_sorpeq, 'do_sorpeq','mmolO2/m^3','[if sorpmeth=1] DO threshold at which sorption exceeds desorption', default=200.0_rk)
   call self%get_parameter(self%Q10,       'Q10',       '-',       'Q10 for sediment processes',          default=2.0_rk)    
   call self%get_parameter(self%Tref,      'Tref',     'celcius', 'reference temperature for bac. proc.',default=20.0_rk)
   
   call self%get_parameter(self%coupled2water, 'coupled2water','-',       'whether coupled with pelagic model',  default=.true.)    
   call self%get_parameter(self%const_O2,    'const_O2',   'mmol/m^3','constant sed. O2, if not coupled', default=14.0_rk)
   call self%get_parameter(self%const_det1C,  'const_det1C',   'mmol/m^3','constant pel. det1C, if not coupled', default=116.0_rk)
   call self%get_parameter(self%const_det1P,  'const_det1P',   'mmol/m^3','constant pel. det1P, if not coupled', default=1.0_rk)
   call self%get_parameter(self%const_det1N,  'const_det1N',   'mmol/m^3','constant pel. det1N, if not coupled', default=14.0_rk)
   call self%get_parameter(self%const_det2C,  'const_det2C',   'mmol/m^3','constant pel. det2C, if not coupled', default=116.0_rk)
   call self%get_parameter(self%const_det2P,  'const_det2P',   'mmol/m^3','constant pel. det2P, if not coupled', default=1.0_rk)
   call self%get_parameter(self%const_det2N,  'const_det2N',   'mmol/m^3','constant pel. det2N, if not coupled', default=14.0_rk)
   call self%get_parameter(self%const_pelC,   'const_pelDIC',   'mmol/m^3','constant pel. DIC, if not coupled', default=116.0_rk)
   call self%get_parameter(self%const_pelP,   'const_pelDIP',   'mmol/m^3','constant pel. DIP, if not coupled', default=1.0_rk)
   call self%get_parameter(self%const_pelNH4,   'const_pelNH4',   'mmol/m^3','constant pel. NH4, if not coupled', default=1.6_rk)
   call self%get_parameter(self%const_pelNO3,   'const_pelNO3',   'mmol/m^3','constant pel. NO3, if not coupled', default=16.0_rk)
   call self%get_parameter(self%const_pelSi,   'const_pelSi',   'mmol/m^3','constant pel. Si, if not coupled', default=15.0_rk)
   if (self%resolve_Si) then
     call self%get_parameter(self%const_det2Si,  'const_det2Si',   'mmol/m^3','constant pel. det2Si, if not coupled', default=15.0_rk)
     !call self%get_parameter(self%const_DISi,   'const_DISi',    'mmol/m^3','constant pel. DISi, if not coupled', default=15.0_rk) 
   end if
   
   ! Register state variables
   !POM pool in sediment:
   !C
   call self%register_state_variable(self%id_POMsed_C,'POC','mmol C/m^2','particulate organic carbon in sediment', minimum=0.0_rk)
   call self%add_to_aggregate_variable(standard_variables%total_carbon,self%id_POMsed_C)
   !P
   call self%register_state_variable(self%id_POMsed_P,'POP','mmol P/m^2','particulate organic detrital phosphorus in sediment', minimum=0.0_rk)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_POMsed_P)
   !N
   call self%register_state_variable(self%id_POMsed_N,'PON','mmol N/m^2','particulate organic detrital nitrogen in sediment', minimum=0.0_rk)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_POMsed_N)
   !Si
   if (self%resolve_Si) then
     call self%register_state_variable(self%id_POMsed_Si,'POSi','mmol Si/m^2','particulate organic detrital silica in sediment', minimum=0.0_rk)
     call self%add_to_aggregate_variable(standard_variables%total_silicate,self%id_POMsed_Si)
   end if
   !(optional) DIM pool in sediment:
   if (self%resolve_sedDIM) then
     !C
     if (self%resolve_DIC) then
       call self%register_state_variable(self%id_DIMsed_C,'DIC','mmol C/m^2','dissolved inorganic carbon in sediment', minimum=0.0_rk)
       call self%add_to_aggregate_variable(standard_variables%total_carbon,self%id_DIMsed_C)
     end if
     !P
     call self%register_state_variable(self%id_DIMsed_P,'DIP','mmol P/m^2','dissolved inorganic phosphorus in sediment', minimum=0.0_rk)
     call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_DIMsed_P)
     !N
     call self%register_state_variable(self%id_DIMsed_N,'DIN','mmol N/m^2','dissolved inorganic nitrogen in sediment', minimum=0.0_rk)
     call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_DIMsed_N)
     !Si
     if (self%resolve_Si) then
       call self%register_state_variable(self%id_DIMsed_Si,'DISi','mmol Si/m^2','dissolved inorganic silica in sediment', minimum=0.0_rk)
       call self%add_to_aggregate_variable(standard_variables%total_silicate,self%id_DIMsed_Si)
     end if 
   end if
   
   !Regester state dependencies
   
   
   if (self%coupled2water) then
    !O2
    call self%register_state_dependency(self%id_O2_pel, 'pelagic_O2') !'mmol m-3','O2 in pelagic')
    !C
    call self%register_state_dependency(self%id_detpel1_C, 'pelagic_det1C') ! 'mmol m-3','slow sinking detrital C variable in pelagic')
    call self%register_state_dependency(self%id_detpel2_C, 'pelagic_det2C') ! 'mmol m-3','fast sinking detrital C variable in pelagic')
    if (self%resolve_DIC) then
        call self%register_state_dependency(self%id_DIMpel_C, 'pelagic_DIC') !'mmol m-3','DIC in pelagic')
      end if
    !P
    call self%register_state_dependency(self%id_detpel1_P, 'pelagic_det1P') !'mmol m-3','slow sinking detrital P variable in pelagic')
    call self%register_state_dependency(self%id_detpel2_P, 'pelagic_det2P') !'mmol m-3','fast sinking detrital P variable in pelagic')
    call self%register_state_dependency(self%id_DIMpel_P, 'pelagic_DIP') !'mmol m-3','DIP in pelagic')
    !N
    call self%register_state_dependency(self%id_detpel1_N, 'pelagic_det1N') !'mmol m-3','slow sinking detrital N variable in pelagic')
    call self%register_state_dependency(self%id_detpel2_N, 'pelagic_det2N') !'mmol m-3','fast sinking detrital N variable in pelagic')
    call self%register_state_dependency(self%id_DIMpel_NO3, 'pelagic_NO3') !'mmol m-3','NO3 in pelagic')
    call self%register_state_dependency(self%id_DIMpel_NH4, 'pelagic_NH4') !'mmol m-3','NH4 in pelagic')
    !Si
    if (self%resolve_Si) then
      call self%register_state_dependency(self%id_detpel2_Si, 'pelagic_det2Si') !,   'mmol m-3','fast sinking detrital Si variable in pelagic')
      call self%register_state_dependency(self%id_DIMpel_Si, 'pelagic_DISi') !'mmol m-3','DISi in pelagic')
    end if
   end if
   
   ! Register diagnostic variables
   !O2
   call self%register_diagnostic_variable(self%id_o2o_brm,'o2o_brm','mmol/m^2/d',  'O2 removal in benthos',          &
                     output=output_time_step_averaged)
   !C
   call self%register_diagnostic_variable(self%id_det1_sed_C,'det1_sed_C','mmol/m^2/d',  'det1C sedimentation flux',    &
                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_det2_sed_C,'det2_sed_C','mmol/m^2/d',  'det2C sedimentation flux',    &
                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_sed_rem_C,'sed_dim_C','mmol/m^2/d',  'benthic C remineralization', &
                     output=output_time_step_averaged)
   !P
   call self%register_diagnostic_variable(self%id_det1_sed_P,'det1_sed_P','mmol/m^2/d',  'det1P sedimentation flux',    &
                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_det2_sed_P,'det2_sed_P','mmol/m^2/d',  'det2P sedimentation flux',    &
                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_sed_rem_P,'sed_dim_P','mmol/m^2/d',  'benthic P remineralization', &
                     output=output_time_step_averaged)                  
   call self%register_diagnostic_variable(self%id_POM_QP,'POM_QP','molP/molC',  'molar P:C ratio of detritus',       &
                     output=output_time_step_averaged)
   if (self%sedO2meth .ne. 0) then 
     call self%register_diagnostic_variable(self%id_sedO2,'sedO2','-',  'estimated O2 in sediment',       &
                     output=output_time_step_averaged)
   end if
   if ((self%sorpmeth .eq. 1) .and. self%resolve_sedDIM) then
      call self%register_diagnostic_variable(self%id_fsorbed,'fsorbed','-',  'fraction of P sorbed in iron complexes',       &
                     output=output_time_step_averaged)
   end if

   !N
   call self%register_diagnostic_variable(self%id_det1_sed_N,'det1_sed_N','mmol/m^2/d',  'det1N sedimentation flux',   &
                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_det2_sed_N,'det2_sed_N','mmol/m^2/d',  'det2N sedimentation flux',   &
                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_sed_rem_N,'sed_dim_N','mmol/m^2/d', 'benthic N remineralization', &
                     output=output_time_step_averaged)  
   call self%register_diagnostic_variable(self%id_POM_QN,'POM_QN','molN/molC',  'molar N:C ratio of detritus',      &
                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_sed_nn2,'sed_nn2','molN/m^2/d',  'benthic denitrification rate',  &
                     output=output_time_step_averaged)
   !Si
   if (self%resolve_Si) then
     call self%register_diagnostic_variable(self%id_det2_sed_Si,'det2_sed_Si','mmol/m^2/d',  'det2Si sedimentation flux',  &
                     output=output_time_step_averaged)
     call self%register_diagnostic_variable(self%id_sed_rem_Si,'sed_dim_Si','mmol/m^2/d', 'benthic Si remineralization',&
                     output=output_time_step_averaged) 
     call self%register_diagnostic_variable(self%id_POM_QSi,'POM_QSi','molSi/molC',  'molar Si:C ratio of detritus',    &
                     output=output_time_step_averaged)
   end if
   
   !if DIM is resolved, save exchange rates of DIM between sed <-> pel
   if (self%resolve_sedDIM) then
     if (self%resolve_DIC) then
       call self%register_diagnostic_variable(self%id_sed_wdim_C,'sed_wdim_C','mmol/m^2/d',  'DIM-C flux towards water', &
                     output=output_time_step_averaged)
     end if
     call self%register_diagnostic_variable(self%id_sed_wdim_P,'sed_wdim_P','mmol/m^2/d',  'DIM-P flux towards water', &
                     output=output_time_step_averaged)
     call self%register_diagnostic_variable(self%id_sed_wdim_N,'sed_wdim_N','mmol/m^2/d', 'DIM-N flux towards water', &
                     output=output_time_step_averaged)
     if (self%resolve_Si) then
       call self%register_diagnostic_variable(self%id_sed_wdim_Si,'sed_wdim_Si','mmol/m^2/d', 'DIM-Si flux towards water', &
                     output=output_time_step_averaged)
     end if
   end if

   ! Register environmental dependencies
   call self%register_dependency(self%id_temp,standard_variables%temperature)

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
   class (type_gpm_abio_sed_EH),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
! !LOCAL VARIABLES:
   real(rk)                   :: fT,temp,spdecr
   real(rk)                   :: fsorbed,trapP,freeP
   !svars
   type(type_elms)            :: POMsed,detpel1,detpel2
   real(rk)                   :: O2_sed,NO3_pel
   !rates
   type(type_elms)            :: det1_sed,det2_sed,dPOMs
   type(type_dim)             :: sed_rem,sed_wdim,DIMsed,DIMpel
   real(rk)                   :: sed_nn2,sed_n4n,o2o_brm,no3_brm,bdnf_basic,osw,nsw
! !PARAMETERS
   real(rk), parameter        :: p_seitz=0.116 !conversion factor f_o2o_brm to bdnf_basic
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops over the horizontal domain (if any).
   _HORIZONTAL_LOOP_BEGIN_
   
   !-------------------------------------------------------------------------
   !PREPARE: retrieve variables
   !
   ! Retrieve current environmental conditions.
   _GET_(self%id_temp,temp)  ! temperature
   
   ! Retrieve local state variables:
   _GET_HORIZONTAL_(self%id_POMsed_C,POMsed%C) ! POC density - sediment
   _GET_HORIZONTAL_(self%id_POMsed_P,POMsed%P) ! POP density - sediment
   _GET_HORIZONTAL_(self%id_POMsed_N,POMsed%N) ! PON density - sediment
   !Si
   if (self%resolve_Si) then
     _GET_HORIZONTAL_(self%id_POMsed_Si,POMsed%Si) ! det-Si density - sediment
   end if
   
   ! Retrieve local state variables:
   if (self%resolve_sedDIM) then
     if (self%resolve_DIC) then
       _GET_HORIZONTAL_(self%id_DIMsed_C,DIMsed%C) ! POC density - sediment
     end if
     _GET_HORIZONTAL_(self%id_DIMsed_P,DIMsed%P) ! POP density - sediment
     _GET_HORIZONTAL_(self%id_DIMsed_N,DIMsed%NH4) ! PON density - sediment
     !Si
     if (self%resolve_Si) then
       _GET_HORIZONTAL_(self%id_DIMsed_Si,DIMsed%Si) ! det-Si density - sediment
     end if
   end if
   
   ! to do: move the sedimentation calculation to the pelagic modules
   ! Pelagic variables:
   if (self%coupled2water) then
     !O2
     if (self%sedO2meth .eq. 0) then !O2 in the sediment is equal to the O2 in the pelagic bottom layer 
       _GET_(self%id_O2_pel,O2_sed)
     else if (self%sedO2meth .eq. 1) then !estimate first the sed-o2 from temperature (as in Kerimoglu et al 2017)
       O2_sed=max(0.0, 300.0_rk-self%K_T2do*temp)
     !else if !there may be a relationship between pel-o2 vs benthic-POM or salinity?
     end if
     
     !C
     _GET_(self%id_detpel1_C,detpel1%C)    ! slow sinking det-C concentration - pelagic
     _GET_(self%id_detpel2_C,detpel2%C)    ! fast sinking det-C concentration - pelagic
     !P
     _GET_(self%id_detpel1_P,detpel1%P)    ! slow sinking det-P concentration - pelagic
     _GET_(self%id_detpel2_P,detpel2%P)    ! fast sinking det-P concentration - pelagic
     !N
     _GET_(self%id_detpel1_N,detpel1%N)    ! slow sinking det-N concentration - pelagic
     _GET_(self%id_detpel2_N,detpel2%N)    ! fast sinking det-N concentration - pelagic
     _GET_(self%id_DIMpel_NO3,DIMpel%NO3)  ! NO3 concentration - pelagic
     !Si
     if (self%resolve_Si) then
       _GET_(self%id_detpel2_Si,detpel2%Si)    ! fast sinking det-Si concentration - pelagic
     end if
     ! Retrieve local state variables:
     if (self%resolve_sedDIM) then
       if (self%resolve_DIC) then
         _GET_(self%id_DIMpel_C,DIMpel%C) ! POC density - sediment
       end if
       _GET_(self%id_DIMpel_P,DIMpel%P) ! POP density - sediment
       _GET_(self%id_DIMpel_NH4,DIMpel%NH4) ! PON density - sediment
       !Si
       if (self%resolve_Si) then
         _GET_(self%id_DIMpel_Si,DIMpel%Si) ! det-Si density - sediment
       end if
     end if
   
   else 
     !O2
     O2_sed=self%const_O2                ! no coupling -constant sediment O2
     !C
     detpel1%C = self%const_det1C          ! no coupling - constant pelagic det1C - slow sinking
     detpel2%C = self%const_det2C          ! no coupling - constant pelagic det2C - fast sinking
     !P
     detpel1%P = self%const_det1P          ! no coupling - constant pelagic det1P - slow sinking
     detpel2%P = self%const_det2P          ! no coupling - constant pelagic det2P - fast sinking
     !N
     detpel1%N = self%const_det1N          ! no coupling - constant pelagic det1N - slow sinking
     detpel2%N = self%const_det2N          ! no coupling - constant pelagic det2N - fast sinking
     DIMpel%NO3  = self%const_NO3           ! no coupling -constant pelagic NO3
     !Si
     if (self%resolve_Si) then
       detpel2%Si = self%const_det2Si      ! no coupling - constant pelagic detSi - fast sinking
     end if
     
     if (self%resolve_sedDIM) then
      if (self%resolve_DIC) then
        DIMpel%C=self%const_pelC
      end if
      DIMpel%P=self%const_pelP
      DIMpel%NH4=self%const_pelNH4
      DIMpel%Si=self%const_pelSi
     end if
   end if
   !
   !END OF PREPARE
   !-------------------------------------------------------------------------
   
   
   !-------------------------------------------------------------------------
   !START OF CALCULATIONS: no other calculation outside this box
   !
   ! sedimentation rate !todo: move this to pelagic modules 
   !C
   det1_sed%C = self%det1_wsed*detpel1%C 
   det1_sed%P = self%det1_wsed*detpel1%P
   det1_sed%N = self%det1_wsed*detpel1%N 
   det2_sed%C = self%det2_wsed*detpel2%C
   det2_sed%P = self%det2_wsed*detpel2%P
   det2_sed%N = self%det2_wsed*detpel2%N
   !Si
   if (self%resolve_Si) then
     det2_sed%Si = self%det2_wsed*detpel2%Si
   end if
   
   ! Temp function
   fT=self%Q10**((temp-self%Tref)/10)
   
   !oxic-anoxic switch (EH L.589-592)
   osw=1.0
   if (O2_sed <= 0.) osw=0.0
   nsw=1.0
   if (DIMpel%NO3 <= 0.1) nsw=0.0
   
   ! remineralization
   !C
   sed_rem%C= self%brc*POMsed%C !*fT?
   !O2
   o2o_brm = osw*sed_rem%C &            !oxic
           + (1.-osw)*(1-nsw)*sed_rem%C !anoxic
   !N
   sed_n4n = self%brn*POMsed%N  !in EH this is h_bon_n4n: uncorrected (ok:total?) benthic N remineralization
   bdnf_basic = p_seitz*o2o_brm !potential denitrification based on o2 consumption (Seitzinger&Giblin 96)
   sed_nn2 = bdnf_basic - max(0.0,bdnf_basic-sed_n4n) !sed_nn2 in EH bon_n2n
   sed_rem%NH4 = max(0.0,sed_n4n-bdnf_basic)   !bon_n4n benthic N remineralization
   !ok: it should be sed_n4n - sed_nn2? 
   no3_brm = 0.5*(1.-osw)*nsw*sed_rem%NH4
   !P
   sed_rem%P = self%brp*POMsed%P
   !Si
   if (self%resolve_Si) then
     sed_rem%Si = self%brsi*POMsed%Si
   end if
   
   if (self%resolve_sedDIM) then
     if (self%resolve_DIC) then
       sed_wdim%C=Dsw*(DIMsed%C/Dasl - DIMpel%C)/Dasl !(+: conc higher in sediment -> flux towards water; -: flux towards sediment)
       ! Flux = K*(delC/delZ)
       ! mmol/m2/s = m2/s * mmol/m3 *1/s
     end if
     if (self%sorpmeth .eq. 0) then 
       sed_wdim%P=Dsw*(DIMsed%P/Dasl - DIMpel%P)/Dasl
     else if (self%sorpmeth .eq. 1) then !instantaneous sorption under oxic conditions
       ! based on O2 in sediment, estimate the fraction of benthic-P trapped in iron-complexation phase (fsorbed), 
       ! which is unavailable for pelagic exchange (as in Kerimoglu et al. 2017)
       fsorbed=1.0_rk/(1.0_rk+exp(0.05*(self%do_sorpeq-O2_sed)))
       trapP=fsorbed*DIMsed%P
       freeP=(1.0-fsorbed)*DIMsed%P
       sed_wdim%P = Dsw*(freeP/Dasl - DIMpel%P)/Dasl
     end if
     sed_wdim%NH4=Dsw*(DIMsed%NH4/Dasl - DIMpel%NH4)/Dasl
     if (self%resolve_Si) then
       sed_wdim%Si=Dsw*(DIMsed%Si/Dasl - DIMpel%Si)/Dasl
     end if
   else !if DIM is not explicitly represented, rem is assigned as the transport towards water (sed_wdim)
     if (self%resolve_DIC) then
       sed_wdim%C=sed_rem%C
     end if
     sed_wdim%P=sed_rem%P
     sed_wdim%NH4=sed_rem%NH4
     if (self%resolve_Si) then
       sed_wdim%Si=sed_rem%Si
     end if
   end if
   !
   !END OF CALCULATIONS: no other calculation outside this box
   !-------------------------------------------------------------------------
   
   !-------------------------------------------------------------------------
   !WRITE
   !
   ! Set local temporal derivatives of sediment variables
   !C
   _SET_BOTTOM_ODE_(self%id_POMsed_C,det1_sed%C+det2_sed%C - sed_rem%C)
   _SET_BOTTOM_ODE_(self%id_POMsed_P,det1_sed%P+det2_sed%P - sed_rem%P)
   _SET_BOTTOM_ODE_(self%id_POMsed_N,det1_sed%N+det2_sed%N - sed_rem%NH4 - sed_nn2)
   !Si
   if (self%resolve_Si) then
     _SET_BOTTOM_ODE_(self%id_POMsed_Si,det2_sed%Si - sed_rem%Si)
   end if
   
   if (self%resolve_sedDIM) then
     if (self%resolve_DIC) then
       _SET_BOTTOM_ODE_(self%id_DIMsed_C, sed_rem%C - sed_wdim%C)
     end if
     _SET_BOTTOM_ODE_(self%id_DIMsed_P,sed_rem%P - sed_wdim%P)
     _SET_BOTTOM_ODE_(self%id_DIMsed_N,sed_rem%NH4 - sed_wdim%NH4)
     if (self%resolve_Si) then
       _SET_BOTTOM_ODE_(self%id_DIMsed_Si,sed_rem%Si - sed_wdim%Si)
     end if
   end if
   
   ! Set bottom fluxes of pelagic variables (these mirror local sediment derivatives), if coupled
   if (self%coupled2water) then
     !O2
     _SET_BOTTOM_EXCHANGE_(self%id_O2_pel,-o2o_brm)
     !water -> sediment fluxes
     _SET_BOTTOM_EXCHANGE_(self%id_detpel1_C,-det1_sed%C)
     _SET_BOTTOM_EXCHANGE_(self%id_detpel2_C,-det2_sed%C)
     _SET_BOTTOM_EXCHANGE_(self%id_detpel1_P,-det1_sed%P)
     _SET_BOTTOM_EXCHANGE_(self%id_detpel2_P,-det2_sed%P)
     _SET_BOTTOM_EXCHANGE_(self%id_detpel1_N,-det1_sed%N)
     _SET_BOTTOM_EXCHANGE_(self%id_detpel2_N,-det2_sed%N)
     if (self%resolve_Si) then
       _SET_BOTTOM_EXCHANGE_(self%id_detpel2_Si,-det2_sed%Si)
     end if
     
     !sediment -> water fluxes
     if (self%resolve_DIC) then
       _SET_BOTTOM_EXCHANGE_(self%id_DIMpel_C, +sed_wdim%C)
     end if
     _SET_BOTTOM_EXCHANGE_(self%id_DIMpel_P, +sed_wdim%P)
     _SET_BOTTOM_EXCHANGE_(self%id_DIMpel_NH4, +sed_wdim%NH4)
     _SET_BOTTOM_EXCHANGE_(self%id_DIMpel_NO3, -no3_brm)
     if (self%resolve_Si) then
       _SET_BOTTOM_EXCHANGE_(self%id_DIMpel_Si, +sed_wdim%Si)
     end if
     
   end if 
   
   ! Export diagnostic variables
   !O2
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_o2o_brm,o2o_brm*s2d)
   !C
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_det1_sed_C,det1_sed%C*s2d)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_det2_sed_C,det2_sed%C*s2d)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_sed_rem_C,sed_rem%C*s2d)
   !P
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_det1_sed_P,det1_sed%P*s2d)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_det2_sed_P,det2_sed%P*s2d)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_sed_rem_P,sed_rem%P*s2d)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_POM_QP,POMsed%P/POMsed%C)
   if ((self%sedO2meth .ne. 0) .and. (self%coupled2water)) then
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_sedO2, o2_sed)
   end if
   if ((self%sorpmeth .eq. 1) .and. self%resolve_sedDIM) then
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fsorbed, fsorbed)
   end if
   !N
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_det1_sed_N,det1_sed%N*s2d)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_det2_sed_N,det2_sed%N*s2d)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_sed_rem_N,sed_rem%NH4*s2d)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_POM_QN,POMsed%N/POMsed%C)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_sed_nn2,sed_nn2*s2d)
   
   !Si
   if (self%resolve_Si) then
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_det2_sed_Si,det2_sed%Si*s2d)
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_sed_rem_Si,sed_rem%Si*s2d)
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_POM_QSi,POMsed%Si/POMsed%Si)
   end if
   
   if (self%resolve_sedDIM) then
     if (self%resolve_DIC) then
       _SET_HORIZONTAL_DIAGNOSTIC_(self%id_sed_wdim_C,sed_wdim%C*s2d)
     end if
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_sed_wdim_P,sed_wdim%P*s2d)
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_sed_wdim_N,sed_wdim%NH4*s2d)
     if (self%resolve_Si) then
       _SET_HORIZONTAL_DIAGNOSTIC_(self%id_sed_wdim_Si,sed_wdim%Si*s2d)
     end if
   end if
   !
   !END OF WRITE
   !-------------------------------------------------------------------------
   
   ! Leave spatial loops over the horizontal domain (if any).
   _HORIZONTAL_LOOP_END_

   end subroutine do_bottom
!EOC

!-----------------------------------------------------------------------

   end module gpm_abio_sed_EH
