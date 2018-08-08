#include "fabm_driver.h"
#include "fabm.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: gpm_abio_pel --- GPM pelagic abiotic model: DI(N,P,Si),det(C,N,P,Si)
!
! !INTERFACE:
   module gpm_abio_pel
!
! !DESCRIPTION:
!
! !USES:
   use fabm_types

!  default: all is private.
   private

! !REVISION HISTORY:!
!
! !PUBLIC DERIVED TYPES:
   type, extends (type_base_model),public :: type_gpm_abio_pel
!     Variable identifiers
      type (type_state_variable_id)     :: id_detC,id_detP,id_detN,id_detSi,id_DIP,id_DIN,id_DISi
      type (type_dependency_id)         :: id_temp
      type (type_global_dependency_id)   :: id_doy
      type (type_diagnostic_variable_id)  :: id_detQP,id_detQN,id_detQSi
      type (type_horizontal_diagnostic_variable_id)  :: id_detC_Bflux
      type (type_horizontal_diagnostic_variable_id)  :: id_DIP_Sflux,id_DIP_Bflux,id_detP_Bflux
      type (type_horizontal_diagnostic_variable_id)  :: id_DIN_Sflux,id_DIN_Bflux,id_detN_Bflux
      type (type_horizontal_diagnostic_variable_id)  :: id_DISi_Sflux,id_DISi_Bflux,id_detSi_Bflux

!     Model parameters
      !common
      real(rk) :: Q10,Tref
      logical  :: resolve_P,resolve_N,resolve_Si
      !Det
      real(rk) :: det_w,det_kc,det_decr !,det_c2p
      integer  :: det_velmet,det_BCtype
      !DIM
      real(rk) :: relrate_dir,DIP_bot,DIP_Sflux,DIN_bot,DIN_Sflux,DISi_bot,DISi_Sflux
      integer  :: DIM_BCtype
      
      contains
      
      procedure  :: initialize
      procedure  :: do
      procedure  :: do_bottom
      procedure  :: do_surface
      procedure  :: get_light_extinction
      procedure  :: get_vertical_movement
      
   end type 
   
   !Local variables
   LOGICAL :: debug = .false.
   
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the mat-model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the parameters and variables are registered det_with FABM.
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   class (type_gpm_abio_pel), intent(inout), target  :: self
   integer,                  intent(in)              :: configunit
!
! !LOCAL VARIABLES:
real(rk), parameter :: d_per_s = 1._rk/86400.
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Store parameter values in our odet_wn derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   !general
   call self%get_parameter(self%resolve_P,  'resolve_P',   '-',           'whether to resolve P cycle',          default=.true.) 
   call self%get_parameter(self%resolve_N,  'resolve_N',   '-',           'whether to resolve N cycle',          default=.false.)  
   call self%get_parameter(self%resolve_Si,  'resolve_Si',  '-',           'whether to resolve Si cycle',          default=.false.)  
   call self%get_parameter(self%Q10,        'Q10',         '-',           'Q10 for bacterial processes',         default=2.0_rk)    
   call self%get_parameter(self%Tref,       'Tref',        'celcius',     'reference temperature for bac. proc.',default=20.0_rk)
   !det
   call self%get_parameter(self%det_kc,     'det_kc',      'm^2/mmolC',  'specific light extinction',           default=0.03_rk)
   call self%get_parameter(self%det_velmet, 'det_velmet',  '-',           'velocity method of detritus',         default=1)
   call self%get_parameter(self%det_w,      'det_w',       'm/d',       'vertical velocity (<0 for sinking)',  default=-1.0_rk, scale_factor=d_per_s)
   call self%get_parameter(self%det_BCtype, 'det_BCtype',  '-',           'boundary condition type for detritus',default=0)
   call self%get_parameter(self%det_decr,   'det_decr',    'd-1',         'detritus decay rate',                 default=0.1_rk,  scale_factor=d_per_s)
   !call self%get_parameter(self%det_c2p,   'det_c2p',      'molC/molP',   'det_c2p ratio',                              default=116.0_rk)
   !nut
   call self%get_parameter(self%DIM_BCtype, 'DIM_BCtype',  '-',           'boundary condition type at bottom',   default=0)
   call self%get_parameter(self%relrate_dir,'relrate_dir', 'd-1',         'relaxation rate for Dirichlet bc',    default=1.0_rk, scale_factor=d_per_s)
   !P
   call self%get_parameter(self%DIP_bot,    'DIP_bot',     'mmol/m^3',    'concentration of DIP at bottom to relax',    default=2.0_rk)
   call self%get_parameter(self%DIP_Sflux,  'DIP_Sflux',   'mmol/m^2/d','DIP flux at surface',                 default=0.0_rk, scale_factor=d_per_s)
   !N
   call self%get_parameter(self%DIN_bot,    'DIN_bot',     'mmol/m^3',    'concentration of DIN at bottom to relax',    default=2.0_rk)
   call self%get_parameter(self%DIN_Sflux,  'DIN_Sflux',   'mmol/m^2/d',  'DISi flux at surface',                 default=0.0_rk, scale_factor=d_per_s)
   !Si
   call self%get_parameter(self%DISi_bot,   'DISi_bot',     'mmol/m^3',   'concentration of DISi at bottom to relax',    default=2.0_rk)
   call self%get_parameter(self%DISi_Sflux, 'DISi_Sflux',   'mmol/m^2/d', 'DISi flux at surface',                 default=0.0_rk, scale_factor=d_per_s)
   
   ! Register state variables
   !C
   call self%register_state_variable(self%id_detC,'detC','mmol C/m^3','detrital carbon in water', & 
                                    minimum=_ZERO_, specific_light_extinction=self%det_kc,vertical_movement=self%det_w*d_per_s)
   call self%add_to_aggregate_variable(standard_variables%total_carbon,self%id_detC)
   !call self%register_dependency(self%id_carbon_vmean,vertical_mean(self%id_totcarb))
   !P
   if (self%resolve_P) then
     call self%register_state_variable(self%id_detP,'detP','mmol P/m^3','detrital phosphorus in water', & 
                                    minimum=_ZERO_,vertical_movement=self%det_w*d_per_s)
     call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_detP)
     call self%register_state_variable(self%id_DIP,'DIP','mmol P/m^3','dissolved inorganic phosphorus in water', & 
                                    minimum=_ZERO_)
     call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_DIP)
   end if
   !N
   if (self%resolve_N) then
     call self%register_state_variable(self%id_detN,'detN','mmol N/m^3','detrital nitrogen in water', & 
                                    minimum=_ZERO_,vertical_movement=self%det_w*d_per_s)
     call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_detN)
     call self%register_state_variable(self%id_DIN,'DIN','mmol P/m^3','dissolved inorganic nitrogen in water', & 
                                    minimum=_ZERO_)
     call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_DIN)
   end if
   !Si
   if (self%resolve_Si) then
     call self%register_state_variable(self%id_detSi,'detSi','mmol Si/m^3','detrital silica in water', & 
                                    minimum=_ZERO_,vertical_movement=self%det_w*d_per_s)
     call self%add_to_aggregate_variable(standard_variables%total_silicate,self%id_detSi)
     call self%register_state_variable(self%id_DISi,'DISi','mmol Si/m^3','dissolved inorganic silica in water', & 
                                    minimum=_ZERO_)
     call self%add_to_aggregate_variable(standard_variables%total_silicate,self%id_DISi)
   end if 
   
   ! Register link to external pools

   ! Register diagnostic variables
   !C
   if (self%det_BCtype .eq. 20) then
     call self%register_horizontal_diagnostic_variable(self%id_detC_Bflux,'detC_f_bot','mmolC/m^2/d', ' detC bottom flux',   &
                                          output=output_time_step_averaged)
   end if
   !P
   if (self%resolve_P) then
     call self%register_diagnostic_variable(self%id_detQP,'det_QP','molP/molC', 'molar P:C ratio of detritus',   &
                                          output=output_time_step_averaged)
     call self%register_horizontal_diagnostic_variable(self%id_DIP_Sflux,'dip_f_surf','mmolP/m^2/d', 'DIP surface flux',   &
                                          output=output_time_step_averaged)
     if (self%det_BCtype .eq. 20) then
       call self%register_horizontal_diagnostic_variable(self%id_detP_Bflux,'detP_f_bot','mmolP/m^2/d', ' detP bottom flux',   &
                                          output=output_time_step_averaged)
     end if
     if (self%DIM_BCtype .eq. 10) then
       call self%register_horizontal_diagnostic_variable(self%id_DIP_Bflux,'dip_f_bot','mmolP/m^2/d', ' DIP bottom flux',   &
                                          output=output_time_step_averaged)
     end if                                            
   end if
   !N
   if (self%resolve_N) then
     call self%register_diagnostic_variable(self%id_detQN,'det_QN','molN/molC', 'molar N:C ratio of detritus',   &
                                          output=output_time_step_averaged)
     call self%register_horizontal_diagnostic_variable(self%id_DIN_Sflux,'din_f_surf','mmolN/m^2/d', 'DIN surface flux',   &
                                          output=output_time_step_averaged)
     if (self%det_BCtype .eq. 20) then
       call self%register_horizontal_diagnostic_variable(self%id_detN_Bflux,'detN_f_bot','mmolN/m^2/d', ' detN bottom flux',   &
                                          output=output_time_step_averaged)
     end if
     if (self%DIM_BCtype .eq. 10) then
       call self%register_horizontal_diagnostic_variable(self%id_DIN_Bflux,'din_f_bot','mmolN/m^2/d', ' DIN bottom flux',   &
                                          output=output_time_step_averaged)
     end if
   end if
   !Si
   if (self%resolve_Si) then
     call self%register_diagnostic_variable(self%id_detQSi,'det_QSi','molSi/mol C', 'molar Si:C ratio of detritus',   &
                                          output=output_time_step_averaged)
     call self%register_horizontal_diagnostic_variable(self%id_DISi_Sflux,'DISi_f_surf','mmolSi/m^2/d', 'DISi surface flux',   &
                                          output=output_time_step_averaged)
     if (self%det_BCtype .eq. 20) then
       call self%register_horizontal_diagnostic_variable(self%id_detSi_Bflux,'detSi_f_bot','mmolSi/m^2/d', ' detSi bottom flux',   &
                                          output=output_time_step_averaged)
     end if
     if (self%DIM_BCtype .eq. 10) then
       call self%register_horizontal_diagnostic_variable(self%id_DISi_Bflux,'detSi_f_bot','mmolSi/m^2/d', ' DISi bottom flux',   &
                                          output=output_time_step_averaged)
     end if
   end if
   
   ! Register conserved quantities

   ! Register environmental dependencies
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_global_dependency(self%id_doy,standard_variables%number_of_days_since_start_of_the_year)

   return

!99 call self%fatal_error('gpm_abio_pel_init','Error reading namelist npzd')

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of Detritus model
!
! !INTERFACE:
   subroutine do(self,_FABM_ARGS_DO_RHS_)
!
! !DESCRIPTION:
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   class (type_gpm_abio_pel), intent(in)     :: self
   _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
   real(rk)                   :: detP,detC,detN,detSi
   real(rk)                   :: det_decr,temp,doy
!EOP
!-----------------------------------------------------------------------
!BOC
   if (debug) then
     _GET_GLOBAL_(self%id_doy,doy) !day of year
     write(*,'(A,2F9.5)')'doy:',doy
   end if
   
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_
   
   ! Retrieve current environmental conditions.
   _GET_(self%id_temp,temp)  ! temperature
   !_GET_GLOBAL_(self%id_doy,doy) !day of year
   !write(*,'(A,2F7.3)')'doy,t:',doy,temp
   
   ! Retrieve current (local) state variable values.
   !C
   _GET_(self%id_detC,detC) ! detritus C
   !P
   if (self%resolve_P) then
     _GET_(self%id_detP,detP) ! detritus P
   end if
   !N
   if (self%resolve_N) then
     _GET_(self%id_detN,detN) ! detritus N
   end if 
   !Si
   if (self%resolve_Si) then
     _GET_(self%id_detSi,detSi) ! detritus Si
   end if 
   
   det_decr=self%det_decr*self%Q10**((temp-self%Tref)/10)

   ! Set temporal derivatives
   !C
   _SET_ODE_(self%id_detC,-det_decr*detC)
   !P
   if (self%resolve_P) then
     _SET_ODE_(self%id_detP,-det_decr*detP)
     _SET_ODE_(self%id_DIP, det_decr*detP)
   end if
   !N
   if (self%resolve_N) then
     _SET_ODE_(self%id_detN,-det_decr*detN)
     _SET_ODE_(self%id_DIN, det_decr*detN)
   end if
   !Si
   if (self%resolve_Si) then
     _SET_ODE_(self%id_detSi,-det_decr*detSi)
     _SET_ODE_(self%id_DISi, det_decr*detSi)
   end if
   
   ! Export diagnostic variables
   !P
   if (self%resolve_P) then
     _SET_DIAGNOSTIC_(self%id_detQP,detP/detC)
   end if
   !N
   if (self%resolve_N) then
     _SET_DIAGNOSTIC_(self%id_detQN,detN/detC)
   end if
   !Si
   if (self%resolve_Si) then
     _SET_DIAGNOSTIC_(self%id_detQSi,detSi/detC)
   end if
   
   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

   end subroutine do
!EOC

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: Right hand sides of bottom exchange model
!
! !INTERFACE:
   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
!
! !DESCRIPTION:
! !IROUTINE: Process interaction betdet_ween benthos and bottom layer of the
! pelagic. This calculates the fluxes into all bottom pelagic and benthic variables,
! in variable quantity per surface area per time. This typically imples variable units * m/s
! [bottom fluxes] for the pelagic, and variable units/s [temporal derivatives] for the benthos.
! Positive values denote state variable increases, negative values state variable decreases.
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   class(type_gpm_abio_pel), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
! !LOCAL VARIABLES:
   real(rk)                   ::vert_vel,temp,Bflux
   real(rk)                   ::detC_bot,detP_bot,DIP_bot,detN_bot,DIN_bot,detSi_bot,DISi_bot
!EOP
!-----------------------------------------------------------------------
!BOC

   !Enter spatial loops over the horizontal domain (if any).
   _FABM_HZ_LOOP_BEGIN_
   
   !det
   select case (self%det_BCtype)
     case (0) !do nothing
     case (10) !1x:family of dirichlet bc's. 10= constant value. Not implemented for the detritus
        !Retrieve current (local) state variable values.
        !_GET_STATE_(self%id_DIP,DIP_bot) 
        !_SET_BOTTOM_EXCHANGE_(self%id_DIP,(self%nbot-DIP_bot)*d_per_s) !nutrient restoring : dirichlet b.c
        call self%fatal_error('gpm_abio_pel_do_benthos', 'doBC=10 not implemented')  
     case (20) !2x: family of neumann bc's. 20= constant flux. concentration at the bottom layer * sinking rate
        !calculate the sinking velocity
        select case (self%det_velmet)
          case (1) 
            vert_vel=self%det_w
          case (2) !velocity is corrected by the viscosity, calculated as a function of temperature
            _GET_DEPENDENCY_  (self%id_temp,temp)  ! temperature
            !in the expression belodet_w, denominator (10**..) gives the viscosity ratio at the ambient and reference (T=20 oC) temperatures (Kestin et al. 1978). So the overall expression is: vel*mu(20)/mu(T)
            vert_vel=self%det_w*1/(10**(((20-temp)/(temp+96))*(1.2378 - 1.303e-3*(20-temp) + 3.06e-6*(20-temp)**2 + 2.55e-8*(20-temp)**3))) 
          end select
          !Retrieve current (local) state variable values and set the rates
          !C
          _GET_STATE_(self%id_detC,detC_bot)
          Bflux=vert_vel*detC_bot
          _SET_BOTTOM_EXCHANGE_(self%id_detC,Bflux)             !detritus-C sedimentation
          _SET_HORIZONTAL_DIAGNOSTIC_(self%id_detC_Bflux, Bflux)
          !P
          if (self%resolve_P) then
            _GET_STATE_(self%id_detP,detP_bot)
            Bflux=vert_vel*detP_bot
            _SET_BOTTOM_EXCHANGE_(self%id_detP,Bflux)             !detritus-P sedimentation
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_detP_Bflux, Bflux)
          end if
          !N
          if (self%resolve_N) then
            _GET_STATE_(self%id_detN,detN_bot)
            Bflux=vert_vel*detN_bot
            _SET_BOTTOM_EXCHANGE_(self%id_detN,Bflux)             !detritus-N sedimentation
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_detN_Bflux, Bflux)
          end if
          !N
          if (self%resolve_Si) then
            _GET_STATE_(self%id_detSi,detSi_bot)
            Bflux=vert_vel*detSi_bot
            _SET_BOTTOM_EXCHANGE_(self%id_detSi,Bflux)             !detritus-Si sedimentation
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_detSi_Bflux, Bflux)
          end if
   end select
   
   !nut
   select case (self%DIM_BCtype)
     case (0) !do nothing
     case (10) !1x:family of dirichlet bc's. 10= constant value. Not implemented for the detritus
        !Retrieve current (local) state variable values and set the rates
        !P
        if (self%resolve_P) then
          _GET_STATE_(self%id_DIP,DIP_bot) 
          Bflux=(self%DIP_bot-DIP_bot)*self%relrate_dir
          _SET_BOTTOM_EXCHANGE_(self%id_DIP,Bflux) !nutrient restoring : (fake) dirichlet b.c
          _SET_HORIZONTAL_DIAGNOSTIC_(self%id_DIP_Bflux, Bflux)
        end if
        !N
        if (self%resolve_N) then
          _GET_STATE_(self%id_DIN,DIN_bot) 
          Bflux=(self%DIN_bot-DIN_bot)*self%relrate_dir
          _SET_BOTTOM_EXCHANGE_(self%id_DIN,Bflux) !nutrient restoring : (fake) dirichlet b.c
          _SET_HORIZONTAL_DIAGNOSTIC_(self%id_DIN_Bflux, Bflux)
        end if
        !Si
        if (self%resolve_Si) then
          _GET_STATE_(self%id_DISi,DISi_bot) 
          Bflux=(self%DISi_bot-DISi_bot)*self%relrate_dir
          _SET_BOTTOM_EXCHANGE_(self%id_DISi,Bflux) !nutrient restoring : (fake) dirichlet b.c
          _SET_HORIZONTAL_DIAGNOSTIC_(self%id_DISi_Bflux, Bflux)
        end if
     case (20) !2x: family of neumann bc's. 20= constant flux. concentration at the bottom layer * sinking rate
        !Retrieve current (local) state variable values.
        !_GET_STATE_(self%id_detP,detb)
        !_SET_BOTTOM_EXCHANGE_(self%id_detP,vert_vel*detb)             !detritus burial
        call self%fatal_error('gpm_components_nut_do_benthos', 'doBC=20 not implemented')  
   end select
   
   !Leave spatial loops over the horizontal domain (if any).
   _FABM_HZ_LOOP_END_

   end subroutine do_bottom
!EOC

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: Right hand sides of surface exchange model
!
! !INTERFACE:
   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
!
! !DESCRIPTION:
! !IROUTINE: Process interaction betdet_ween benthos and bottom layer of the
! pelagic. This calculates the fluxes into all bottom pelagic and benthic variables,
! in variable quantity per surface area per time. This typically imples variable units * m/s
! [bottom fluxes] for the pelagic, and variable units/s [temporal derivatives] for the benthos.
! Positive values denote state variable increases, negative values state variable decreases.
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   class (type_gpm_abio_pel), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_
!
! !LOCAL VARIABLES:

   !real(rk) :: n_Sflux
!EOP
!-----------------------------------------------------------------------
!BOC
   
   ! Enter spatial loops (if any)
   _HORIZONTAL_LOOP_BEGIN_

   ! nut flux due to rivers, surface runoff, deposition
   !P
   if (self%resolve_P) then
     _SET_SURFACE_EXCHANGE_(self%id_DIP, self%DIP_Sflux)
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_DIP_Sflux, self%DIP_Sflux)
   end if
   !N
   if (self%resolve_N) then
     _SET_SURFACE_EXCHANGE_(self%id_DIN, self%DIN_Sflux)
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_DIN_Sflux, self%DIN_Sflux)
   end if
   !Si
   if (self%resolve_Si) then
     _SET_SURFACE_EXCHANGE_(self%id_DISi, self%DISi_Sflux)
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_DISi_Sflux, self%DISi_Sflux)
   end if
     
   !Leave spatial loops over the horizontal domain (if any).
   _HORIZONTAL_LOOP_END_

   end subroutine do_surface
!EOC


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the vertical velocity of pelagic biogeochemical variables
!
! !INTERFACE:
   subroutine get_vertical_movement(self,_FABM_ARGS_GET_VERTICAL_MOVEMENT_)
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   class (type_gpm_abio_pel), intent(in) :: self
   _DECLARE_FABM_ARGS_GET_VERTICAL_MOVEMENT_
!
! !LOCAL VARIABLES:
   real(rk)                   :: vert_vel,temp !,yearday,par
!
!EOP
!-----------------------------------------------------------------------
!BOC

! det_way to get the day of year
!_GET_DEPENDENCY_SCALAR_(self%id_yearday,yearday)  ! day of year (diff(d/m/Y - 1/1/Y))

  ! Enter spatial loops (if any)
  _FABM_LOOP_BEGIN_

  select case (self%det_velmet)
    case (1) 
     vert_vel=self%det_w
    case (2) !velocity is corrected by the viscosity, calculated as a function of temperature
     _GET_DEPENDENCY_  (self%id_temp,temp)  ! temperature
     !in the expression belodet_w, denominator (10**..) gives the viscosity ratio at the ambient and reference (T=20 oC) temperatures (Kestin et al. 1978). So the overall expression is: vel*mu(20)/mu(T)
     vert_vel=self%det_w*1/(10**(((20-temp)/(temp+96))*(1.2378 - 1.303e-3*(20-temp) + 3.06e-6*(20-temp)**2 + 2.55e-8*(20-temp)**3))) 
  end select

  !Set these calculated vertical_movement values for the current p 
  _SET_VERTICAL_MOVEMENT_(self%id_detC,vert_vel)
  if (self%resolve_P) then
    _SET_VERTICAL_MOVEMENT_(self%id_detP,vert_vel) 
  end if
  if (self%resolve_N) then
    _SET_VERTICAL_MOVEMENT_(self%id_detN,vert_vel)
  end if
  if (self%resolve_Si) then
    _SET_VERTICAL_MOVEMENT_(self%id_detSi,vert_vel)
  end if
! Leave spatial loops (if any)
_FABM_LOOP_END_

end subroutine get_vertical_movement
!EOC
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the light extinction coefficient due to biogeochemical
! variables
!
! !INTERFACE:
   subroutine get_light_extinction(self,_ARGUMENTS_GET_EXTINCTION_)
!
! !USES:
   implicit none
!   
! !INPUT PARAMETERS:
   class(type_gpm_abio_pel), intent(in)     :: self
   _DECLARE_ARGUMENTS_GET_EXTINCTION_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   real(rk)                     :: detC
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_detC,detC) ! detritus-carbon
   !detC=detP*self%det_c2p !detritus-carbon

   ! Self-shading det_with explicit contribution from detritus-C concentration.
   _SET_EXTINCTION_(self%det_kc*detC)
   !_SET_EXTINCTION_(0.01)

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

   end subroutine get_light_extinction
!EOC
!-----------------------------------------------------------------------

   end module gpm_abio_pel
