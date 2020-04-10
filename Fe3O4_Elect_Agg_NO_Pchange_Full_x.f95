!  1-D, transient, rectangular, electrode problem solved by the method of Newman.
! Nicholas Brady 5/11/2017

!   The user must specify
!     1. parameters for the reaction kinetics, operating conditions, and electrode structure
!     2. an initial guess
!     3. the coefficients aik, bik, dik, gi, pik, eik, fi
!       (see J.S. Newman, Electrochemical Systems, p. 540-)

!________________________________________________
! program description
!       Finite Volume Formulation
!       Diffusion through a rectangular medium
!       First Order Reaction
!       Fixed Concentrations at the boundaries
!
! ********** Governing Equations **********
!       electrolyte concentration (1): D*d2cdx2  + a*irxn/Fconst = dcdt
!       solid-state concentration (2):           - a*irxn/Fconst = dcs/dt
!       solid-state potential     (3): sigma*d2Phi_1dx2 - a*irxn = 0
!       electrolyte potential     (4): kappa*d2Phi_2dx2 + a*irxn = 0
!
! ********** Boundary Conditions **********
!       at x=0,     dcdx=0,   dVdr=0,               Phi_2 =        
!       at x=xmax,  c=1,      -sigma*dVdr=-i_app,   dPhi_2dx = 0
!       at t=0,     c=1,      V=U0,                 Phi_2 = 0,     ci = 0
!________________________________________________
!******************************INPUT MODULE******************************************

module user_input
SAVE    
     parameter(N=4,NJ=22,N_c=5,NJ_c=52,Numbertimesteps=3.6d3*500,tmax=3.6d3*500)
     ! 3.6d3 = 1 hour
     
          ! *********************** ELECTRODE PARAMETERS ***********************
     double precision :: xmax                                     ! Electrode size

          ! *********************** GENERAL  PARAMETERS ***********************
     double precision :: Rigc = 8.314,Temp = 298, Fconst = 96485  ! Ideal gas constant [J/(mol*K)], Temperature [K], Faraday's Constant [C/mol]
     real             :: PI=3.141592654                           ! Geometric constant

          ! *********************** ELECTROLYTE  ***********************
     double precision :: eps_tot = 0.63, eps_elect = 0.47, diff = 1.00d-6 * dfdf              ! electrode porosity, diffusion coefficient
     double precision :: eps_agg = 0.26, diff_agg = 2.25d-13
     double precision :: c0, cbulk=0.001, c0min=1.0d-22           ! electrolyte concentration, bulk electrolyte conc, and minimum conc
     double precision :: kappa
          
          ! *********************** ELECTRODE MATERIAL PHYSICAL PARAMETERS ***********************
     double precision :: sigma=4.269d-2        ! solid-state electronic conductivity
     double precision :: vmin=0.5                                 ! minimum voltage allowed (numerically)
     double precision :: molar_mass_Fe3O4 = 231.533, density_Fe3O4 = 5.15 ! [g/mol],  [g/cm3]
     double precision :: molar_volume_Fe3O4 ! [mol/cm3] - density / molar_mass

          ! *********************** INITIAL CONDITIONS ***********************
     double precision :: Phi_1_init=3.0, Phi_2_init = 0.0, c0_init=0.001, cs_init = 1.0d-5   ! initial conditions: voltage, electrolyte conc, solid-state conc.

          !****************** CHANGES BY EXPERIMENTAL CONDITIONS ***********************
     double precision :: Discharge_Time, Discharge_Relax_Time, Charge_Time, Charge_Relax_Time           !
     double precision :: c_specific, c_density, c_specific_agg, c_density_agg                            ! [A/g], [A/cm2]
     double precision :: spec_a, spec_a_agg    

     ! ********************* ELECTROLYTE CONDUCTIVITY ***************************
     double precision :: z_cat = +1.0, z_an = -1.0            ! cation and anion species charge
     double precision :: diff_cat, diff_an                    ! cation and anion diffusion coefficients
     double precision :: xmax_c, diff_c = 1.00d-13

          !****************** COULOMB COUNTING ******************
     double precision :: mAhg = 0.0, LixFe3O4
     double precision :: c_test = 0.0
     double precision :: to_electrons ! multiply by this to convert mol/cm3 to electrons (mol/mol)

          !****************** VARIABLE: STATE OF EXPERIMENT ***********************
     CHARACTER(LEN=1) :: state
     CHARACTER(LEN=4) :: electron_equiv
     integer(kind=8)  :: Discharge

        ! N: number of dependent variables
        ! NJ: number of node points (for 1-D, 101 is typical)
        ! tmax: duration of the simulation [=] seconds
        ! c_on: time when current is applied [=] seconds

        ! Rigc: ideal gas constant [=] J/mol/K
        ! Temp: temeperature [=] K
        ! c_applied: applied current density [=] A/cm^2
        ! c0min_fail: lower bound on unreacted concentration in pores [=] mol/cm^3
        ! cimax_fail: upper bound on reacted concentration [=] mol/cm^3
        
        ! eps: porosity
        ! diff_agg: diffusivity coefficent [=] cm^2/s
        ! sigma: electronic conductivity [=] S/cm
        ! spec_a: specificy surface area [=] cm^2/cm^3
        ! xmax: radius of the particle [=] cm
        ! cbulk: concentration of lithium ions in the electroltyte [=] mol/cm^3
        
        ! rxn_k: reaction rate [=] 1/(mol*cm)**(0.5)/s
        ! alpha_a/c: charge transfer coefficient
        ! cimax: maximum concentration of reacted lithium in Fe3O4 crystal [=] mol/cm^3
!_______________________

               
end module user_input

!************************ ELECTROCHEMICAL REACTION MODULE **********************************

module MOD_echem_rxn
  use user_input
  SAVE

  double precision :: cimax=1.788337d-1                        ! maximum solid-state concentration [=] mol/cm3
  double precision :: alpha_a = 0.5, alpha_c = 0.5             ! cathodic and anodic charge transfer coefficients

contains 
 FUNCTION OSP(c00, cii)                                        ! Function to calculate the Open-Circuit Potential
  double precision :: OSP, cii, xish, c00                      ! xish is the normalized solid-state concentration
  double precision :: U_ref=1.561677                           ! U_ref, RK_Ax are constants used to calculate the solid-state potential
  double precision :: RK_A0=-6.581135d-1, RK_A1=6.586353d-3, RK_A2=1.224915d-1, RK_A3=2.765079d-1
  double precision :: RK_A4=-5.146964d-1, RK_A5=-1.204923d-4, RK_A6=-4.364893d-8, RK_A7=1.109934d-1
  double precision, DIMENSION(0:7) :: AK

      xish = cii/cimax

      !!! NEW NEW !!!

!       U_ref = 1.518
!       RK_A0 = -6.4722d-1
!       RK_A1 = 1.2827d-1
!       RK_A2 = 3.2177d-1
!       RK_A3 = 1.5766d-1
!       RK_A4 = -6.6700d-1
!       RK_A5 = 3.4446d-1
!       RK_A6 = 1.2277d-1
!       RK_A7 = -1.4046d-1      

      !!!!!!!!!!!!!!!

      AK = (/ RK_A0, RK_A1, RK_A2, RK_A3, RK_A4, RK_A5, RK_A6, RK_A7 /)

      Vint=0.0
  
      do kk=0,7
        Vint=Vint+AK(kk)*((2*xish-1)**(kk+1)-(2*xish*kk*(1-xish))/(2*xish-1)**(1-kk))
      end do

      OSP = U_ref + Rigc*Temp/Fconst*DLOG(c00/cbulk*(cimax-cii)/cii) + Vint - 0.200

RETURN
END FUNCTION OSP

FUNCTION EX_I(c00,cii) ! exchange current
  double precision :: EX_I, c00, cii
  double precision :: rxn_k=5.62d-9
! c00 is the electrolyte concentration
! cii is the solid-state concentration
! rxn_k is the reaction rate [=] 1/(mol*cm)**(0.5)/s

  EX_I = Fconst*rxn_k*(c00**alpha_a)*((cimax-cii)**alpha_a)*(cii**alpha_c)

RETURN
END FUNCTION EX_I


! Electrochemical Reaction
! ex_i * [exp(alpha_a*F*eta/RT) - exp(-alpha_c*F*eta/RT)]

FUNCTION echem_rxn(c00, cii, Phi_1, Phi_2) 
  double precision :: echem_rxn, eta, U0, Phi_1, Phi_2
  double precision :: c00, cii

U0 = OSP(c00, cii)
eta = Phi_1 - Phi_2 - U0

ex_curr = EX_I(c00,cii)

echem_rxn = ex_curr * ( DEXP(alpha_a*Fconst*eta/(Rigc*Temp)) - DEXP(-alpha_c*Fconst*eta/(Rigc*Temp)) )


RETURN
END FUNCTION echem_rxn

end module MOD_echem_rxn

! ______________________________________________________________________________________

! *****************************************************************************
! **************************** PHASE CHANGE MODULE ****************************
! **************************** PHASE CHANGE MODULE ****************************
! *****************************************************************************

module MOD_Phase_Change

  use user_input
  SAVE

  contains
    FUNCTION Rxn_Beta(c_alpha, thet_beta)
      double precision :: c_alpha, thet_beta ! input values: alpha phase conc., beta_phase volume fraction
      double precision :: k_beta, gamma, c_alpha_sat, m_exp

!     m_exp  = 2.0/3.0
    m_exp  = 0.0
    k_beta = 1.0d-4 * 70.0
    gamma  = 4.0
    c_alpha_sat  = 5.2 * molar_volume_Fe3O4

    ! Forming and Unforming Beta Phase are not identically reversible
    ! See Paper: Discharge, Relaxation, and Charge Model for the Lithium Trivanadate Electrode: Reactions, Phase CHange, and Transport

!     if (c_alpha.GE.c_alpha_sat) then ! forming Beta Phase
!       Rxn_Beta = k_beta*(c_alpha - c_alpha_sat)*thet_beta**m_exp*(1.0 - thet_beta)

    
!     elseif (c_alpha.LT.c_alpha_sat) then ! unforming Beta Phase
      
!       if (theta_beta.LE.0.0) then
!         Rxn_Beta = 0.0 ! we can't unform beta phase if there is no beta phase

!       else
!         Rxn_Beta = k_beta*(c_alpha - c_alpha_sat)*(1.0 - thet_beta)**m_exp*(thet_beta)

!       end if

!     end if


  if (c_alpha.GE.c_alpha_sat) then ! do not calculate pchange formation if c_alpha < c_alpha_sat
    Rxn_Beta = k_beta*(c_alpha - c_alpha_sat)*thet_beta**m_exp*(1.0 - thet_beta)
!     Rxn_Beta = k_beta * (c_alpha - c_alpha_sat)

  else
    Rxn_Beta = 0.0

  end if
      
      

    RETURN
  END FUNCTION Rxn_Beta

end module MOD_Phase_Change


!******************************VARIABLE MODULE******************************************
! do not edit this module

module variables
use user_input
SAVE
     integer :: NP1
! Declare array/matrix sizes depending on the user inputted values.
     double precision, dimension(N,N+1,NJ)  :: E
     double precision, dimension(N,2*N+1)   :: D
     double precision, dimension(N,N)       :: A,B,X,Y,smA,smB,smD,smE,smP
     double precision, dimension(N,NJ)      :: cprev, delC
     double precision, dimension(N)         :: dcdx,d2cdx2,G,smG,smF,ID
     double precision, dimension(NJ)        :: xx, delx, i_rxn_0, total_lithium_electrode
     double precision :: delT,time             
     integer :: t1,t2,clock_rate,clock_max
     !****************** FVM VARIABLES ***********************
     double precision, dimension(N,N) :: fE, fW, dE, dW, rj
     double precision, dimension(N)   :: cE, cW, dcdxE, dcdxW
     double precision :: alphaE, alphaW, betaE, betaW

!*****************************************
! ********* Crystal VARIABLES ************
!*****************************************
     integer :: NP1_c
! Declare array/matrix sizes depending on the user inputted values.
     double precision, dimension(N_c,N_c+1,NJ_c)  :: E_c
     double precision, dimension(N_c,2*N_c+1)     :: D_c
     double precision, dimension(N_c,N_c)         :: A_c,B_c,X_c,Y_c,smA_c,smB_c,smD_c,smE_c,smP_c
     double precision, dimension(N_c)             :: dcdx_c,d2cdx2_c,G_c,smG_c,smF_c,ID_c
     double precision, dimension(NJ_c)            :: xx_c,delx_c 
     double precision, dimension(N_c,NJ_c)        :: delC_c
     double precision, dimension(NJ,N_c,NJ_c)     :: cprev_c

     !****************** FVM VARIABLES ***********************
     double precision, dimension(N_c,N_c) :: fE_c, fW_c, dE_c, dW_c, rj_c
     double precision, dimension(N_c)     :: cE_c, cW_c, dcdxE_c, dcdxW_c
     double precision :: alphaE_c, alphaW_c, betaE_c, betaW_c

end module variables

!**********************SUBROUTINES FOR WRITING DATA******************************
!____________________________________________________________________________
!   This section of code stores all of the subroutines used to write data files.
!   It is up to the user to edit the "call" in the main program to get the right data.
!____________________________________________________________________________

subroutine write_all_voltage(it)
    use user_input
    use variables
    use MOD_echem_rxn
    implicit double precision(a-h,o-z)
    double precision :: Phi_1
    CHARACTER(50) :: position

    ! ***** NOTE *********
    ! Headers CANNOT HAVE SPACES - this messes up the code to plot the data
    ! underscores are preferable, but dashes can be used as well
    ! i.e. Solution Pot --> Solution_Pot

    ! Write data at all Electrode positions
    open(56, file = 'Time_Voltage_Position.txt', status = 'unknown')

    ! Write data at all Agglomerate positions at the separator
    open(551, file = 'Time_Voltage_Position_Agg_Separator.txt', status = 'unknown')

    ! Write data at all Agglomerate positions at the current collector
    open(5511, file = 'Time_Voltage_Position_Agg_Middle.txt', status = 'unknown')

    ! Write data at all Agglomerate positions at the current collector
    open(5522, file = 'Time_Voltage_Position_Agg_Collector.txt', status = 'unknown')

    ! Write just edge information
    open(57, file = 'Time_Voltage.txt', status = 'unknown')

    ! Theta Beta information
    open(60, file = 'Time_Theta_Beta.txt', status = 'unknown')

    t_write=time/float(3600)

! Just the  information from node NJ
      c0 = cprev(1,1)
      
! solid-state concentration
      cs = cprev(4,1)
      theta_beta = cprev_c(1,5,NJ_c)
      
      Phi_1 = cprev(2,NJ)
      Phi_2 = cprev(3,NJ)

      U0 = OSP(c0, cs)


      
!         i_rxn_0(j) = echem_rxn(c0,cs,Phi_1, Phi_2)

      LixFe3O4 = mAhg*molar_mass_Fe3O4*3.6/Fconst ! 1 mAh = 3.6 C

      ! [mol Li/ cm3 Fe3O4] --> [mol Li/ mol Fe3O4]
      ! [mol Li / cm3 Fe3O4] * [cm3/g Fe3O4] * [g/mol Fe3O4] = mol Li / mol Fe3O4
      to_electrons = 1.0/density_Fe3O4*molar_mass_Fe3O4

! **** TERMINAL SCREEN INFORMATION ****
      if (it.EQ.1) then
!           write the headers on the first entrance into write all voltage
        write(*,17) 'State', 'Time', 'Voltage', 'Solution_Pot', 'mAh/g', 'LixFe3O4', 'Edge c0', 'Edge cs', &
        &           'edge reaction', 'edge OCV', 'c_density', 'kappa', 'c(SS) NJ_c = 1', 'Theta_beta'
        write(*,17) 'CDR',  'hours', 'Volts',  'Volts',         'mAh/g', 'LixFe3O4', 'mol/L', 'LixFe3O4', &
        &           ','            , 'Volts',    'A/cm2'    , 'S/cm',  'LixFe3O4',  'Fraction'
        ! need to write units 

      end if

!         c0 = cprev_c(1,1,NJ_c)
      
!         ! solid-state concentration
!         cs = cprev_c(1,4,NJ_c)
!         !         theta_beta = cprev_c(1,2,NJ_c)
      
!         Phi_1 = cprev_c(l,2,NJ_c)
!         Phi_2 = cprev_c(l,3,NJ_c)

!         U0 = OSP(c0, cs)


!         write(*,15) state, t_write, Phi_1, Phi_2, mAhg, LixFe3O4, c0*1000, cs*to_electrons, i_rxn_0(1), U0, c_density, kappa, &
!         &           cprev_c(1,1,NJ_c) * to_electrons, theta_beta
      
      difference_lithium = total_lithium_electrode(1) - total_lithium_electrode(NJ)

      write(*,15) state, t_write, Phi_1, Phi_2, mAhg, LixFe3O4, c0*1000, cs*to_electrons, i_rxn_0(1), U0, c_density, kappa, &
      &           cprev_c(1,4,NJ_c) * to_electrons, theta_beta, total_lithium_electrode(1), total_lithium_electrode(NJ), &
      &           difference_lithium, difference_lithium/total_lithium_electrode(1) * 100.0

! **** TEXT FILE INFORMATION ****
! * -- 56 - Positional Information
! * -- 57 - Just Edge Information
      if (it.EQ.1) then
!           write the headers on the first entrance into write all voltage
        write(56,17) 'State', 'Time', 'Voltage', 'Solution_Pot', 'Position', 'mAhg', 'Equivalence', 'Elect_Conc', 'Solid_Conc', &
        &            'Rxn_Rate', 'Theta_beta_Edge', 'Theta_beta_Mid', 'Theta_beta_Cen', 'Vol_avg_Theta_beta', &
        &            'Tot_Lithium_Elect'
        write(56,17) 'CDR',   'hours', 'Volts',  'Volts',        'cm',      'mAh/g', 'LixFe3O4', 'mol/L'  , 'LixFe3O4'   , &
        &            ','  ,      'Fraction'       , 'Fraction'      , 'Fraction'      , 'Fraction'             , &
        &            'mol/cm3'

      end if

      ! Write the information at each electrode position
      do j = 1,NJ
          c0 = cprev(1,j)
          cs = cprev(4,j)
          theta_beta_edge   = cprev_c(j,5,NJ_c)
          theta_beta_mid    = cprev_c(j,5,NJ_c/2)
          theta_beta_center = cprev_c(j,5,1)

          volume_avg_theta_beta = 0.0
          do jj = 2,NJ_c
            volume_avg_theta_beta = volume_avg_theta_beta + 4.0/3.0*PI*( xx_c(jj)**3 - xx_c(jj-1)**3 ) * cprev_c(j,5,jj)
          end do

          volume_avg_theta_beta = volume_avg_theta_beta / (4.0/3.0*PI * xx_c(NJ_c)**3.0)

          
          Phi_1 = cprev(2,j)
          Phi_2 = cprev(3,j)
          U0 = OSP(c0, cs)

          i_rxn_0(j) = echem_rxn(c0,cs,Phi_1, Phi_2)
        
        write(56,15) state, t_write, Phi_1, Phi_2, xx(j), mAhg, LixFe3O4, c0*1000, cs*to_electrons, i_rxn_0(j), theta_beta_edge, &
        &            theta_beta_mid, theta_beta_center, volume_avg_theta_beta, total_lithium_electrode(j)

      
      end do

      !*** Write Electrode Separator Information ***!
      if (it.EQ.1) then
        write(57,17) 'State', 'Time', 'Voltage', 'Solution_Pot', 'mAhg', 'Equivalence', 'Elect_Conc', 'Solid_Conc', 'Rxn_Rate',&
                     'Theta_beta'
        write(57,17) 'CDR',   'hours', 'Volts' ,  'Volts',      'mAh/g', 'LixFe3O4'   , 'mol/L'     , 'LixFe3O4'  , ',' , &
                     'Fraction'
      end if  

      theta_beta = cprev_c(1,5,NJ_c)
      cs  = cprev(4,1)
        write(57,15) state,  t_write,    Phi_1,     Phi_2,         mAhg,   LixFe3O4,     c0*1000,   cs*to_electrons, i_rxn_0(1), &
                     theta_beta



      ! write the agglomerate positional data at the separator and current collector
      ! Separator
      if (it.EQ.1) then
        write(551,17)  'State', 'Time', 'Voltage', 'Solution_Pot', 'Position', 'mAhg', 'Equivalence', 'Elect_Conc', 'Solid_Conc',&
                       'Theta_beta'
        write(551,17) 'CDR',   'hours', 'Volts'  ,  'Volts',      'cm'      , 'mAh/g', 'LixFe3O4'   , 'mol/L'     , 'LixFe3O4'  ,&
                       'Fraction'

        write(5511,17)  'State', 'Time', 'Voltage', 'Solution_Pot', 'Position', 'mAhg', 'Equivalence', 'Elect_Conc', 'Solid_Conc',&
                       'Theta_beta'
        write(5511,17) 'CDR',   'hours', 'Volts'  ,  'Volts',       'cm'      , 'mAh/g', 'LixFe3O4'   , 'mol/L'     , 'LixFe3O4' ,&
                       'Fraction'  


        write(5522,17) 'State', 'Time', 'Voltage', 'Solution_Pot', 'Position', 'mAhg', 'Equivalence', 'Elect_Conc', 'Solid_Conc',&
                       'Theta_beta'
        write(5522,17) 'CDR',   'hours', 'Volts'  ,  'Volts',      'cm'      , 'mAh/g',  'LixFe3O4'   , 'mol/L'     , 'LixFe3O4' ,&
                       'Fraction'   

        write(60,21)   'State', 'Time', 'Voltage', 'Elect_Pos',  'Equivalence', ( j, j=1,NJ_c)
        write(60,17)   'CDR',   'hours', 'Volts'  ,   'cm'    ,  'LixFe3O4',    ('-', j=1,NJ_c)
        
21 format(1(A5,1X), 2(A12,1X),   2(A15,1X), 100(I15, 1X))

      end if 


      do l = 1, NJ
          write(60,15) state, t_write, Phi_1, xx(l), LixFe3O4, (cprev_c(l,5,j), j=1,NJ_c)
      end do
          

      do j = 1,NJ_c
        l = 1 ! electrode separator
        Phi_1 = cprev_c(l,2,j)
        Phi_2 = cprev_c(l,3,j)
        c0    = cprev_c(l,1,j)
        cs    = cprev_c(l,4,j)
        theta_beta = cprev_c(l,5,j)

        write(551,15) state,  t_write,     Phi_1,     Phi_2,    xx_c(j),     mAhg,   LixFe3O4,     c0*1000,   cs*to_electrons,  &
                      theta_beta

        l = NJ/2 ! electrode current collector
        Phi_1 = cprev_c(l,2,j)
        Phi_2 = cprev_c(l,3,j)
        c0    = cprev_c(l,1,j)
        cs    = cprev_c(l,4,j)
        theta_beta = cprev_c(l,5,j)

        write(5511,15) state,  t_write,     Phi_1,     Phi_2,    xx_c(j),     mAhg,   LixFe3O4,     c0*1000,   cs*to_electrons,  &
                      theta_beta                            

        l = NJ ! electrode current collector
        Phi_1 = cprev_c(l,2,j)
        Phi_2 = cprev_c(l,3,j)
        c0    = cprev_c(l,1,j)
        cs    = cprev_c(l,4,j)
        theta_beta = cprev_c(l,5,j)

        write(5522,15) state,  t_write,     Phi_1,     Phi_2,    xx_c(j),     mAhg,   LixFe3O4,     c0*1000,   cs*to_electrons,  &
                      theta_beta     

! 19    format(1(A5,1X), 1(ES15.5,1X), 2(F12.5,1X), 20(ES15.5,1X))                                         
      end do



       

! Format for terminal
! numbers
 15    format(1(A5,1X), 2(F12.5,1X), 60(ES15.5,1X))
! headers
 17    format(1(A5,1X), 2(A12,1X),   60(A15,1X))
 

! Format for text file
 16    format(1(A7,1X), 10(ES15.5,1X))
 18    format()

end subroutine write_all_voltage
      
!****************************** MAIN PROGRAM ******************************************  

program unsteady  
use user_input  
use variables
use MOD_echem_rxn
implicit double precision(a-h,o-z)

! *************************************************************************************************************
! ****************************************** EXPERIMENTAL CONDITIONS ******************************************
! *************************************************************************************************************
! Under this line, conditions that change for each experimental run are defined.
! Using conditional statements the simulated parameters are adjusted for different experiments
! Usually these conditional statements are triggered through the correspond shell script

! Experimental specific parameters include
! xmax_c0         - the crystal size
! xmax            - agglomerate size
! c_applied       - the applied current in A/g (per gram of active material)

! Discharge_ON    - time that current is being passed during discharge
! Discharge_OFF   - time that current is OFF after discharge
! Charge_ON       - time that current is passed during charge; usually this time is irrelevant because charge is done until a cutoff voltage is reached
! Charge_OFF      - time that current is OFF after charge
! Const_Volt_Time - Amount of time for the constant voltage hold (at the end of charge)

! Sometimes, the initial conditions will change as well for different experiments or simulations
! Phi_1_init, c0_init, cs_init, mAhg
! These conditions produce Li1.5Fe3O4
molar_volume_Fe3O4    = density_Fe3O4/molar_mass_Fe3O4
Discharge_Time        = 36.288 * 3600               ! hour * sec/hour
Discharge_Relax_Time  = 5.0                         ! 5.0 seconds
Charge_Time           = 00.0 * 3600.0               ! hour * sec/hour
Charge_Relax_Time     = 30.0 * 24.0 * 3600.0        ! 30 days
c_specific            = 4.77951388889/1000.0 * 2.0  ! A/g
c_specific            = 9.523759/1000.0 * 4.0       ! A/g
c_specific            = 36.7710725244574448 / 1000.0! A/g (taken from the experimental data)
ci_init               = 1.0d-5                      !
cv_init               = 3.3                         !
xmax                  = 1.0d-4 * 450.0              ! Electrode size 400 um
xmax_c                = 1.0d-4 * 1.0 /float(xaxa)               ! Crystal size     1 um


state     = 'D'

spec_a    = 3.0*(1.0 - eps_elect)/(xmax_c) ! check this one
! spec_a    = 3.0 / (xmax_c * density_Fe3O4 * (1.0 - eps_agg))

! convert c_specific to c_density [A/g] -> [A/cm2]
! c_density   = c_specific * 4.0/3.0*PI*(xmax**3)*(1.0-eps)*density_Fe3O4/(4.0*PI*xmax**2) ! g Fe3O4 = agg_vol*(1-eps)*density_Fe3O4;  area of agg = 4*PI*(xmax**2)
c_density   = c_specific * (xmax**3) * (1.0 - eps_tot) * density_Fe3O4 / (xmax**2)

print*, xmax_c
print*, spec_a
print*, c_density
      
call system_clock(t1,clock_rate,clock_max)
call initial_condition(h,h_c)

do j = 1,NJ
  total_lithium_electrode(j) = 0.0
end do

do 1 it=1,Numbertimesteps

    U0 = OSP(cprev(1,1),cprev(4,1))


! !       write one of the early data points
      if (it.EQ.1) then 
          call write_all_voltage(it)  
      
      ! Write 1000 evenly spaced data points
      elseif (MOD(INT(it),INT(Numbertimesteps/1000)).EQ.0) then 
          call write_all_voltage(it)
      
      ! Write a data point near the very end of the program
      elseif (it.GE.(Numbertimesteps)) then 
          call write_all_voltage(it)

      elseif ((cprev(2,NJ).GE.3.0).AND.(state.EQ.'C')) then ! write before exitting;
          call write_all_voltage(it)
          exit

      elseif (ISNAN(cprev(2,NJ))) then ! if the voltage is NaN then write and exit
          call write_all_voltage(it)
          exit

! !       if the the program has to exit for some reason, then write the data before exiting
      elseif (U0.LE.vmin) then
        call write_all_voltage(it)
        exit

      elseif(time.GE.10.0*3600.0) then ! if the total simulation time is 1000 hours or greater, exit.  
        call write_all_voltage(it)
         exit

      end if

!       write(*,15) delC(4,1)*Fconst/density_Fe3O4/delT * 1000.0, delC_c(4,1)*Fconst/density_Fe3O4/delT * 1000.0 


15    format(9(ES20.5,1X))
16    format(9(A20,1X))

! call write_all_voltage(it) 

! dynamic time-stepping so we can quickly move through the relaxation period.  
! Significantly reduces simulation time during long relaxation periods
      
      if (time.EQ.0.0) then
        delT = 1.0/1000.0

      else if (delT.LT.1.0) then
        delT = delT * 1.0002

      else if (state.EQ.'R') then 
        delT = delT * 1.0001 ! * 1.00005

      else
        delT=tmax/float(Numbertimesteps)

      end if

!       delT = 1.0
!       if (cprev(2,1).LE.1.0) then
!       	c_specific = 9.523759/1000.0 * 2.0
!       	c_density  = c_specific * (xmax**3) * (1.0 - eps) * density_Fe3O4 / (xmax**2)
!       end if
  
  call bound_val_e(h)
  call bound_val_c(h_c)

  do j = 1,NJ
    total_lithium_electrode(j) = total_lithium_electrode(j) + delC(4,j)
  end do


      time=time+delT

 1    continue
call system_clock(t2,clock_rate,clock_max)
write ( *, * ) 'Elapsed real time =', real(t2-t1)/real(clock_rate )
end program unsteady  

!*********************************BOUND VAL**************************************      
!____________________________________________________________________________
!   This subroutine assumes solution of a linear boundary-value problem at 
!   a given time step.  There is thus no need for iteration.
!     c=cprev+delC, we solve for delC using forward time central difference method.
!     It is assumed that the time step is small enough that the linearization is exact.
!     Thus, iterations are not necessary.  This assumption is of course tested by 
!     computer experiments on the effect of time-step size.
!____________________________________________________________________________

subroutine bound_val_e(h)
      use user_input
      use variables
      implicit double precision(a-h,o-z)
      
! initialize all of the finite volume (matrix) coefficients to zero

      do 150 j=1,nj
      do 123 ic=1,n
      do 123 kc=1,n
        dE(ic,kc) = 0.0
        dW(ic,kc) = 0.0
        fE(ic,kc) = 0.0
        fW(ic,kc) = 0.0
        rj(ic,kc) = 0.0
123      continue

call fillmat(h,j)
! Fillmat is determining small coefficents based on the user input. Must go to
! fillmat and change coefficent definitions based on governing differential equations.

call ABDGXY(h,j)
! ABDGXY equates the large coefficents based on the small coefficents. These equations
! can be found in Newman appendix C.
 
call BAND(J)
! BAND(J) computes delC and calls MATINV to solve the problem using gaussian elimination.

! for all the dependent variables, the values are update as c = cprev + delC

150   continue
        do 12 k=1,n
         do 12 j=1,nj
             cprev(k,j) = cprev(k,j) + delC(k,j)
12       continue    
      return
end subroutine bound_val_e

!************************************ Crystal ******************************************
!************************************ Crystal ******************************************
!************************************ Crystal ******************************************

subroutine bound_val_c(h_c)
      use user_input
      use variables
      implicit double precision(a-h,o-z)
      
! initialize all of the finite volume (matrix) coefficients to zero
  do 321 l=1,NJ
      do 150 j=1,nj_c
      do 123 ic=1,n_c
      do 123 kc=1,n_c
        dE_c(ic,kc) = 0.0
        dW_c(ic,kc) = 0.0
        fE_c(ic,kc) = 0.0
        fW_c(ic,kc) = 0.0
        rj_c(ic,kc) = 0.0
123      continue

call fillmat_c(h_c,j,l)
! Fillmat_c is determining small coefficents based on the user input. Must go to
! fillmat_c and change coefficent definitions based on governing differential equations.

call ABDGXY_c(h_c,j)
! ABDGXY_c equates the large coefficents based on the small coefficents. These equations
! can be found in Newman appendix C.
 
call BAND_c(J)
! BAND_c(J) computes delC and calls MATINV_c to solve the problem using gaussian elimination.

! for all the dependent variables, the values are update as c = cprev + delC

150   continue
        do 12 ic=1,n_c
         do 12 j=1,nj_c
             cprev_c(l,ic,j) = cprev_c(l,ic,j) + delC_c(ic,j)
12       continue  

321 continue

      return
end subroutine bound_val_c
    


!******************************INITIAL GUESS****************************************

subroutine initial_condition(h, h_c)
      use user_input
      use variables
      implicit double precision(a-h,o-z)

! Define delta t       
!       delT=tmax/float(Numbertimesteps)

!************************************ Electrode ******************************************
!************************************ Electrode ******************************************
!************************************ Electrode ******************************************

      h=xmax/float(nj-2)
      do j=1,NJ

        if (j.EQ.1) then
          xx(j) = 0.0
        else if (j.EQ.NJ) then
          xx(NJ) = xmax
        else
          xx(j)=h*float(j-1) - h/2.0
        end if

      end do

      do j=2,NJ-1
         delx(j)=h
      end do
      ! it is common in the finite volume code to set the control volumes at the boundaries to zero
         delx(1)=0.0d0
         delx(NJ)=0.0d0
          
      do j=1,NJ
        cprev(1,j) = c0_init
        cprev(4,j) = cs_init
        cprev(2,j) = Phi_1_init
        cprev(3,j) = Phi_2_init
      end do

!************************************ Crystal ******************************************
!************************************ Crystal ******************************************
!************************************ Crystal ******************************************

      h_c = xmax_c/float(nj_c-2)
      
      do j=1,NJ_c

        if (j.EQ.1) then
          xx_c(j) = 0.0
        else if (j.EQ.NJ_c) then
          xx_c(NJ_c) = xmax_c
        else
          xx_c(j)=h_c*float(j-1) - h_c/2.0
        end if

      end do

      do j=2,NJ_c-1
         delx_c(j) = h_c
      end do
      ! it is common in the finite volume code to set the control volumes at the boundaries to zero
         delx_c(1) = 0.0d0
         delx_c(NJ_c) = 0.0d0
       
      do l=1,NJ   
        do j=1,NJ_c
          cprev_c(l,1,j) = c0_init
          cprev_c(l,4,j) = cs_init
          cprev_c(l,5,j) = 1.0d-6
          cprev_c(l,2,j) = 3.3 ! Phi_1_init
          cprev_c(l,3,j) = Phi_2_init
        end do
      end do




      return
end subroutine initial_condition


!***********************************FILLMAT****************************************

subroutine fillmat(h,j)
      use user_input
      use variables
      use MOD_echem_rxn
      implicit double precision(a-h,o-z)
      double precision :: cs, Phi_1, Phi_2, U
      double precision :: conc_cat, conc_an, u_cat, u_an
      double precision :: step_c0, step_cs, step_phi_1, step_phi_2
      
!       to_electrons = 1.0/density_Fe3O4*molar_mass_Fe3O4

      c0 = cprev(1,j)     ! electrolyte concentration
      cs = cprev_c(j,4,NJ_c)     ! solid-state concentration
!       cs    = cprev(4,j)
      cprev(4,j) = cs
      Phi_1 = cprev(2,j)  ! solid-state potential
      Phi_2 = cprev(3,j)  ! electrolyte potential

      U = OSP(c0,cs)

      conc_cat = c0
      conc_an  = c0
      diff_cat = diff
      diff_an  = diff
      u_cat    = diff_cat/(Rigc*Temp)
      u_an     = diff_an/(Rigc*Temp)

!       kappa = Fconst**2*(z_cat**2*(diff_cat/(Rigc*Temp))*conc_cat + z_an**2 * (diff_an/(Rigc*Temp))*conc_an)
      kappa = Fconst**2*(z_cat**2*u_cat*conc_cat + z_an**2*u_an*conc_an)

! Butler-Volmer Equation and its partial derivatives

! ________ ELECTROCHEMICAL REACTION RATE [A/cm2]
      
      i_rxn_0(j)    = echem_rxn(c0, cs, Phi_1, Phi_2)

! ________ NUMERICAL DERIVATIVE STEP SIZES
      step_c0    = 1.0d-6
      step_cs    = 1.0d-6
      step_phi_1 = 1.0d-6
      step_phi_2 = 1.0d-6

! _______ NUMERICAL DERIVATIVES WRT c0, ci, phi_1, phi_2
! *** NOTES *** For c0, ci, it is unphysical for those values to be below 0
! For cases, where those values are close to zero we use a simple forward stepping numerical derivative
      if (c0.LE.step_c0) then
        di_rxn_dc0    = ( echem_rxn(c0+step_c0,cs,Phi_1, Phi_2) - echem_rxn(c0,cs,Phi_1, Phi_2) )/(step_c0)
      else
        di_rxn_dc0    = ( echem_rxn(c0+step_c0,cs,Phi_1, Phi_2) - echem_rxn(c0-step_c0,cs,Phi_1, Phi_2) )/(2.0*step_c0)
      end if

      if (cs.LE.step_cs) then
        di_rxn_dcs    = ( echem_rxn(c0,cs+step_cs,Phi_1, Phi_2) - echem_rxn(c0,cs,Phi_1, Phi_2) )/(step_cs)
      else
        di_rxn_dcs    = ( echem_rxn(c0,cs+step_cs,Phi_1, Phi_2) - echem_rxn(c0,cs-step_cs,Phi_1, Phi_2) )/(2.0*step_cs)
      end if
      
      di_rxn_dphi_1 = ( echem_rxn(c0,cs,Phi_1+step_phi_1, Phi_2) - echem_rxn(c0,cs,Phi_1-step_phi_1, Phi_2) )/(2.0*step_phi_1)

      di_rxn_dphi_2 = ( echem_rxn(c0,cs,Phi_1, Phi_2+step_phi_2) - echem_rxn(c0,cs,Phi_1, Phi_2-step_phi_2) )/(2.0*step_phi_2)


!**** NOTES ******
      ! It is important to benchmark the numerical derivatives and have a good way to determine
      ! the appropriate step size that will produce the most accurate results 
      ! In addition this method of multiplying a certain fraction does not work well when the variable
      ! is close to zero
      ! Use a conditional statement when close to zero to use additive steps
      ! if ( ABS(variable).LE.1.0d-10 ) then
      !     step  = 1.0d-10
      !     df_dv = (f(v+step) - f(v-step))/(2*step)
      ! end if

      ! for the concentrations, they are in log() terms and log(-#) is undefined
      ! for these variables it might be more useful to instead do something like:
      ! if ( variable.LE.1.0d-10 ) then
      !     step  = 1.0d-10
      !     df_dv = (f(2*step) - f(step))/(step)
      ! end if

! it is often necessary to have "write" statements inside of fillmat when debugging the code      

!       if (j.EQ.1) then
        
!         if (time.EQ.0.0) then
! !           write(*,19) 'time',   'i_rxn_0', 'di_rxn_dc0', 'di_rxn_dcs', 'di_rxn_dphi_1', 'di_rxn_dphi_2'
!         end if
! !           write(*,20)  time,     i_rxn_0,   di_rxn_dc0,   di_rxn_dcs,   di_rxn_dphi_1,   di_rxn_dphi_2, c0, cs, Phi_1, Phi_2
! !         write(*,21) kappa, Fconst**2, diff**2, (Rigc*Temp)**2
    
!       end if

!  19    format(1(A5,1X),  2(A14,1X),  12(A14,1X))
!  20    format(1(F5.2,1X),2(ES14.5,1X),12(ES14.5,1X))
 21    format(20(ES14.5,1X)) ! all scientific notation

!************************************************************************************
!*********** x = 0 ********************* x = 0 ********************* x = 0 **********
!************************************************************************************

!______________________Boundary-condition at x = 0
! At the boundaries it is usually useful to set the volume to zero (delX(j) = 0)
!     dE(i,k)*dcdxE + rj(i,k)*c = g(i)

      if (j.eq.1) then 

        alphaE=delx(j)/(delx(j+1)+delx(j))
        betaE=2.0/(delx(j)+delx(j+1))

        do ic=1,N
          cE(ic)    = alphaE*cprev(ic,j+1) + (1.d0 - alphaE)*cprev(ic,j)
          dcdxE(ic) = betaE * (cprev(ic,j+1) - cprev(ic,j))
        end do

! *************** Electrolyte Concentration ****************
! Both Fixed Concentration and No Flux work

!  Fixed Concentration Condition
!     c0 = cbulk at x = 0
        rj(1,1) = 1.d0
        smG(1) = cbulk - rj(1,1)*cprev(1,j)

! *************** Solid-State Concentration ****************
! for cs both dcs/dt = rxn and dcs/dx = 0 produce the same results
! because the solid-state concentration does not have any spatial gradients, 
! spatial boundary conditions are also not necessary

! dc/dt = rxn
! dc/dt = i_rxn/Fconst
! dc/dt - (di_rxn_dc * delC + di_rxn_dcs * delCs + di_rxn_dphi_1 * delPhi_1)/Fconst = i_rxn/Fconst

        rj(4,1) = -spec_a*di_rxn_dc0/Fconst
        rj(4,4) = -spec_a*di_rxn_dcs/Fconst - 1.0*(1.0-eps_tot)/delT
!         rj(4,4) = - 1.0*(1.0-eps)/delT
        rj(4,2) = -spec_a*di_rxn_dphi_1/Fconst
        rj(4,3) = -spec_a*di_rxn_dphi_2/Fconst

        smG(4)  = +spec_a*i_rxn_0(j)/Fconst

!         c_test = c_test - spec_a*i_rxn_0(j)/Fconst/(1.0 - eps)

! *************** Solid-State Potential ****************
! At x = 0 this is where the electrode meets the separator
! At this point there cannot be any flow of current through the solid-state

! -(1-eps) * simga * dV/dx = 0 ! no electron flux at the separator
        
        dE(2,2) = -(1.0 - eps_tot)*sigma
        smG(2)  =  0.0 - dE(2,2)*dcdxE(2)    

! *************** Solution Potential ****************    
! At x = 0 this is where the electrode meets the separator
! At this point there cannot be any flow of current through the solid-state
! Because there is no electronic current, all the current must be ionic
! **  -eps * kappa * dPhi_2/dx = c_applied **
! However for some numerical reason this boundary condition is not working so instead
! we are assuming that Phi_2 = 0.0 at x = 0
! Note that this condition still leads to the condition -eps * kappa * dPhi_2/dx = c_applied 
! being met

        rj(3,3) = 1.0
        smG(3)  = 0.0 - rj(3,3)*cprev(3,j)

        return
      end if

!*********************************************************************************
!*********** XMAX ********************* XMAX ********************* XMAX **********
!*********************************************************************************

!______________________Boundary-condition at x=xmax
! At the boundaries it is usually useful to set the volume to zero (delX(j) = 0)
! The general boundary condition equation looks like: dW(i,k)*dcdxW + rj(i,k)*delC = g(i)

      if (j.eq.NJ) then

        alphaW=delx(j-1)/(delx(j-1)+delx(j))
        betaW=2.0/(delx(j-1)+delx(j))

        do ic=1,N
          cW(ic)    = alphaW*cprev(ic,j) + (1.d0 - alphaW)*cprev(ic,j-1)
          dcdxW(ic) = betaW * (cprev(ic,j) - cprev(ic,j-1))
        end do

! *************** Electrolyte Concentration ****************
! at x = xmax this is where the electrode meets the current collector
! at this point there cannot be any flow of electrolyte through the solid electrode
! Therefore there is a no flux condition at xmax

!  No flux Condition
!     dc0/dx = 0 at x = xmax

!         dW(1,1) = 1.0
!         smG(1)  = 0.0 - dW(1,1)*dcdxW(1)

        dW(1,1) = -eps_elect*diff_cat
        fW(1,1) = -eps_elect*z_cat*u_cat*Fconst*dcdxW(3)

        dW(1,3) = -eps_elect*z_cat*u_cat*Fconst*cW(1)
        fW(1,3) = 0.0

        smG(1)  = 0.0 - dW(1,1)*dcdxW(1) - fW(1,1)*cW(1)

! *************** Solid-State Concentration ****************
! No Boundary conditions on cs
! dc/dt = rxn
! dc/dt = i_rxn/Fconst
! dc/dt - (di_rxn_dc * delC + di_rxn_dcs * delCs + di_rxn_dphi_1 * delPhi_1)/Fconst = i_rxn/Fconst

!         rj(2,1) = -spec_a*di_rxn_dc0/Fconst
!         rj(2,2) = -spec_a*di_rxn_dcs/Fconst - (1.0-eps)/delT
!         rj(2,3) = -spec_a*di_rxn_dphi_1/Fconst
!         rj(2,4) = -spec_a*di_rxn_dphi_2/Fconst

!         smG(2)  = +spec_a*i_rxn_0(j)/Fconst
        rj(4,1) = -spec_a*di_rxn_dc0/Fconst
        rj(4,4) = -spec_a*di_rxn_dcs/Fconst - 1.0*(1.0-eps_tot)/delT
!         rj(4,4) = - 1.0*(1.0-eps)/delT
        rj(4,2) = -spec_a*di_rxn_dphi_1/Fconst
        rj(4,3) = -spec_a*di_rxn_dphi_2/Fconst

        smG(4)  = +spec_a*i_rxn_0(j)/Fconst

! *************** Solid-State Potential ****************
! at x = xmax this is where the electrode meets the current collector
! at this point there cannot be any flow of current in the electrolyte
! ie all the current must be in the solid-state (or electronic)

! -(1 - eps) * sigma * dV/dx = c_density at x = xmax

! This boundary condition changes during the experiment:
! Discharge, recovery, charge

! *** NOTE ***
! c_density is defined as a positive (+) number is this code
! in the agglomerate code it was defined as a negative (-) number

      if (time.LE.Discharge_Time) then                                     ! Lithiate
        dW(2,2) = -(1.0-eps_tot)*sigma
        smG(2)  = c_density - dW(2,2)*dcdxW(2)
        state = 'D'
        mAhg  = mAhg + 1000.0*c_specific*delT/3600.0

      else if (time.LE.(Discharge_Time + Discharge_Relax_Time)) then        ! Recover
        dW(2,2) = -(1.0-eps_tot)*sigma
        smG(2)  = 0.0 - dW(2,2)*dcdxW(2)
        state = 'R'
        mAhg  = mAhg 
                                                                            ! Delithiate
      else if ( (time.LE.(Discharge_Time+Discharge_Relax_Time+Charge_Time)) .AND. (cprev(2,NJ).LE.3.0) ) then ! Delithiate
        dW(2,2) = -(1.0-eps_tot)*sigma
        smG(2)  = -c_density - dW(2,2)*dcdxW(2)
        state = 'C'
        mAhg  = mAhg - 1000.0*c_specific*delT/3600.0

      else if (Relax.EQ.1) then                                             ! Recover
        dW(2,2) = -(1.0-eps_tot)*sigma
        smG(2)  = 0.0 - dW(2,2)*dcdxW(2)
        state = 'R'
        mAhg  = mAhg

      else                                                                   ! Recover
        dW(2,2) = -(1.0-eps_tot)*sigma
        smG(2)  = 0.0 - dW(2,2)*dcdxW(2)
        state = 'R'
        mAhg  = mAhg  

      end if

! *************** Solution Potential ****************
! at x = xmax this is where the electrode meets the current collector
! at this point there cannot be any flow of current in the electrolyte
! Therefore a no flux condition is met for ionic current

! @ x = xmax the current carried in the solution equals 0 (all the current is carried in the solid-state)

!         dW(4,4) = -eps*kappa
!         smG(4) = 0.0 - dW(4,4)*dcdxW(4)

!         dW(4,1) = -eps*diff_cat
!         fW(4,1) = -eps*z_cat*u_cat*Fconst*dcdxW(4)

!         dW(4,4) = -eps*z_cat*u_cat*Fconst*cW(1)
!         fW(4,4) = 0.0

!         smG(4)  =
           

        dW(3,1) = -eps_elect*Fconst*(z_cat*diff_cat + z_an*diff_an)
        fW(3,1) = -eps_elect*Fconst**2*(z_cat**2*u_cat + z_an**2*u_an)*dcdxW(3)
        
        dW(3,3) = -eps_elect*Fconst**2*(z_cat**2*u_cat + z_an**2*u_an)*cW(1)
        fW(3,3) = 0.0

        smG(3)  = 0.0 - dW(3,1)*dcdxW(1) - fW(3,1)*cW(1)

        return
      end if


! !______________________Governing Equations
! Finite Difference Method (Control Volume Formulation)
! A material-balance for species i at node k

!       dc(i,k)/dt * delX(k) = Flux_west(i) - Flux_east(i) + R(i)*delX(k)

! where R(i) is the production rate per unit volume, Flux_west(i) is the flux of species i
! at the "west" face of node k, and Flux_east(i) is the flux of species i at the "east" face 
! of node k.

! In general the flux of species i can be affected by all of the variable and their gradients
      
!       Flux_west(i) = d_west(j,k) * dcdx_west(j) + f_west*c_west(j)
!       Flux_east(i) = d_east(j,k) * dcdx_east(j) + f_east*c_east(j)

! The general formula to apply to each node point can be written as:

!      g(i) = rj(i,k)*c + fW(i,k)*cW + dW(i,k)*dcdxW - fE(i,k)*cE - dE(i,k)*dcdxE

! The optimal interpolation formula to use to evaluate the variable c and their derivatives
! depend on the local Peclet number. The general form can be seen below
! For a low Peclet number, a linear interpolation is most appropriate
! Linear interpolation, alpha = 0.5, beta = 1/delX

! *** NOTE *** For large Peclet numbers this central-difference approximation causes an oscillatory
! behavior in the concentration. This can be eliminated by an upwind scheme, in which betaW and betaE
! remain the same, but 
! alphaE = 0, alphaW = 0 when flow is from west to east
! alphaE = 1, alphaW = 1 when flow is from east to west

        alphaW=delx(j-1)/(delx(j-1)+delx(j))
        alphaE=delx(j)/(delx(j+1)+delx(j))
        betaW=2.0/(delx(j-1)+delx(j))
        betaE=2.0/(delx(j)+delx(j+1))


        do ic=1,N
          cW(ic)    = alphaW*cprev(ic,j) + (1.d0 - alphaW)*cprev(ic,j-1)
          cE(ic)    = alphaE*cprev(ic,j+1) + (1.d0 - alphaE)*cprev(ic,j)
          dcdxW(ic) = betaW * (cprev(ic,j) - cprev(ic,j-1))
          dcdxE(ic) = betaE * (cprev(ic,j+1) - cprev(ic,j))
        end do

! *************** Electrolyte Concentration ****************
! dc/dt = D * d2c/dx2 + rxn
! delX * dc/dt = -D * dcdx_W - (-D * dcdx_E) - i_rxn/Fconst * delX
! delX * dc/dt + (di_rxn_dc * delC + di_rxn_dcs * delCs + di_rxn_dphi_1 * delPhi_1 + di_rxn_dphi_2 * delPhi_2)/Fconst * delX
!            = -D * dcdx_W - (-D * dcdx_E) - i_rxn/Fconst*delX

        dW(1,1) = -eps_elect*diff_cat
        dE(1,1) = -eps_elect*diff_cat
        fW(1,1) = -eps_elect*z_cat*u_cat*Fconst*dcdxW(3)
        fE(1,1) = -eps_elect*z_cat*u_cat*Fconst*dcdxE(3)

        dW(1,3) = -eps_elect*z_cat*u_cat*Fconst*cW(1)
        dE(1,3) = -eps_elect*z_cat*u_cat*Fconst*cE(1)
        fW(1,3) = 0.0
        fE(1,3) = 0.0

        rj(1,1) = spec_a*di_rxn_dc0/Fconst*delX(j) - eps_elect/delT*delX(j)
        rj(1,4) = spec_a*di_rxn_dcs/Fconst*delX(j)
        rj(1,2) = spec_a*di_rxn_dphi_1/Fconst*delX(j)
        rj(1,3) = spec_a*di_rxn_dphi_2/Fconst*delX(j)

        smG(1)  = -spec_a*i_rxn_0(j)/Fconst*delX(j) &
        &         - (fW(1,1)*cW(1) + dW(1,1)*dcdxW(1)) &
        &         + (fE(1,1)*cE(1) + dE(1,1)*dcdxE(1))


! *************** Solid-State Concentration ****************
! dc/dt = -rxn
! dc/dt = i_rxn/Fconst
! dc/dt - (di_rxn_dc * delC + di_rxn_dcs * delCs + di_rxn_dphi_1 * delPhi_1 + di_rxn_dphi_2 * delPhi_2)/Fconst = i_rxn/Fconst
        
!         dW(2,2) = 0.0
!         dE(2,2) = 0.0
!         fW(2,2) = 0.0
!         fE(2,2) = 0.0
        
!         rj(2,1) = -spec_a*di_rxn_dc0/Fconst
!         rj(2,2) = -spec_a*di_rxn_dcs/Fconst - (1.0-eps)/delT
!         rj(2,3) = -spec_a*di_rxn_dphi_1/Fconst
!         rj(2,4) = -spec_a*di_rxn_dphi_2/Fconst

!         smG(2)  = +spec_a*i_rxn_0(j)/Fconst
        rj(4,1) = -spec_a*di_rxn_dc0/Fconst
        rj(4,4) = -spec_a*di_rxn_dcs/Fconst - 1.0*(1.0-eps_tot)/delT
!         rj(4,4) = - 1.0*(1.0-eps)/delT
        rj(4,2) = -spec_a*di_rxn_dphi_1/Fconst
        rj(4,3) = -spec_a*di_rxn_dphi_2/Fconst

        smG(4)  = +spec_a*i_rxn_0(j)/Fconst

! *************** Solid-State Potential ****************
! It is assumed that electrons travel very quickly so there is no accumulation of electrons
! (Kirchhoff's laws also assume no accumulation of current)
! dce/dt = sigma * d2V/dx2 + i_rxn 
!      0 = sigma * d2V/dx2 + i_rxn 
!  - delX*(di_rxn_dc * delC + di_rxn_dcs * delCs + di_rxn_dphi_1 * delPhi_1 + di_rxn_dphi_2 * delPhi_2) = -sigma * dcdx_W - (-sigma * dcdx_E) + delX*i_rxn

        dW(2,2) = -(1.0-eps_tot)*sigma
        dE(2,2) = -(1.0-eps_tot)*sigma
        fW(2,2) = 0.0
        fE(2,2) = 0.0

        rj(2,1) = -spec_a*di_rxn_dc0*delX(j)
        rj(2,4) = -spec_a*di_rxn_dcs*delX(j)
        rj(2,2) = -spec_a*di_rxn_dphi_1*delX(j)
        rj(2,3) = -spec_a*di_rxn_dphi_2*delX(j)

        smG(2)  = +spec_a*i_rxn_0(j)*delX(j) & 
        &         - (fW(2,2)*cW(2) + dW(2,2)*dcdxW(2)) &
        &         + (fE(2,2)*cE(2) + dE(2,2)*dcdxE(2))

! *************** Solution Potential **************** 
! Similar to what is assumed in the solid-state, in the electrolyte it is also assumed that ions do not accumulate 
! (also called conservation of charge) therefore dPhi_2/dt = 0 everywhere
! 0 = -eps * kappa * d2Phi_2/dx2 - i_rxn
! + delX*(di_rxn_dc * delC + di_rxn_dcs * delCs + di_rxn_dphi_1 * delPhi_1 + di_rxn_dphi_2 * delPhi_2) = -esp * kappa * dcdx_W - (-eps * kappa * dcdx_E) - delX*i_rxn       

!       kappa = Fconst**2*(z_cat**2*u_cat*conc_cat + z_an**2*u_an*conc_an)
!       kappa = Fconst**2*(z_cat**2*u_cat + z_an**2*u_an)*conc
        
        dW(3,1) = -eps_elect*Fconst*(z_cat*diff_cat + z_an*diff_an)
        dE(3,1) = -eps_elect*Fconst*(z_cat*diff_cat + z_an*diff_an)
        fW(3,1) = -eps_elect*Fconst**2*(z_cat**2*u_cat + z_an**2*u_an)*dcdxW(3)
        fE(3,1) = -eps_elect*Fconst**2*(z_cat**2*u_cat + z_an**2*u_an)*dcdxE(3)
        
        dW(3,3) = -eps_elect*Fconst**2*(z_cat**2*u_cat + z_an**2*u_an)*cW(1)
        dE(3,3) = -eps_elect*Fconst**2*(z_cat**2*u_cat + z_an**2*u_an)*cE(1)
        fW(3,3) = 0.0
        fE(3,3) = 0.0

        rj(3,1) = +spec_a*di_rxn_dc0*delX(j)
        rj(3,4) = +spec_a*di_rxn_dcs*delX(j)
        rj(3,2) = +spec_a*di_rxn_dphi_1*delX(j)
        rj(3,3) = +spec_a*di_rxn_dphi_2*delX(j)

        smG(3)  = -spec_a*i_rxn_0(j)*delX(j) &
        &         - (fW(3,3)*cW(3) + dW(3,3)*dcdxW(3)) &
        &         + (fE(3,3)*cE(3) + dE(3,3)*dcdxE(3))

!         rj(4,4) = 1.0
!         smG(4)  = 0.0 - rj(4,4)*cprev(4,j)

      return

end subroutine fillmat




subroutine fillmat_c(h_c,j,l)
      use user_input
      use variables
      use MOD_echem_rxn
      use MOD_Phase_Change
      implicit double precision(a-h,o-z)
      double precision :: k_a, N_W, N_E, i_rxn_agg
      double precision :: cs, Phi_1, U
      double precision :: conc_cat, conc_an, u_cat, u_an
      double precision :: step_c0, step_cs, step_phi_1, step_phi_2
      double precision :: mA_per_crystal, mA_per_gram

      c0 = cprev_c(l,1,j)           ! electrolyte concentration
      c0_edge = cprev(1,l)
      cs = cprev_c(l,4,j)           ! solid-state concentration
      theta_beta = cprev_c(l,5,j)   ! solid-state theta beta volume fraction
      gamma_conc = 8.0*molar_volume_Fe3O4

      Phi_1 = cprev_c(l,2,j)        ! solid-state potential
      Phi_2 = cprev_c(l,3,j)        ! electrolyte potential
      Phi_2_elect = cprev(3,l)

      xmax_xtal = 4.0d-7
      spec_a_agg    = 3.0*(1.0d0-eps_agg)/(xmax_xtal)


      c_specific_agg  = delC(4,l)/delT * Fconst * (1.0-eps_agg) / density_Fe3O4           ! A/g
!       c_density_agg   = c_specific_agg * (xmax_c**3)*(1.0-eps_agg)*density_Fe3O4/(xmax_c**2) / 3.0 
      c_density_agg   = delC(4,l)/delT * Fconst * (xmax_c)*(1.0 - eps_agg) / 3.0

!       c_density_agg = -i_rxn_0(l)


      U = OSP(c0,cs)

      conc_cat = c0
      conc_an  = c0
      diff_cat = diff_agg
      diff_an  = diff_agg
      u_cat    = diff_cat/(Rigc*Temp)
      u_an     = diff_an/(Rigc*Temp)

!       kappa = Fconst**2*(z_cat**2*(diff_cat/(Rigc*Temp))*conc_cat + z_an**2 * (diff_an/(Rigc*Temp))*conc_an)
      kappa = Fconst**2*(z_cat**2*u_cat*conc_cat + z_an**2*u_an*conc_an)

! Butler-Volmer Equation and its partial derivatives

! ________ ELECTROCHEMICAL REACTION RATE [A/cm2]
      i_rxn_agg       = echem_rxn(c0, cs, Phi_1, Phi_2)  ![=] 

!       if ( (l.EQ.1).OR.(l.EQ.NJ) ) then
!         if ( (j.EQ.1).OR.(j.EQ.NJ_c) ) then
!             mA_per_crystal = spec_a_agg * i_rxn_agg/Fconst * 1000
!             gram_per_crystal = 4.0/3.0 * PI * (xmax_xtal)**3 * density_Fe3O4
!             mA_per_gram    = mA_per_crystal / gram_per_crystal
!             write(*,*) i_rxn_agg, mA_per_crystal, gram_per_crystal, mA_per_gram
!         end if
!       end if

! ________ NUMERICAL DERIVATIVE STEP SIZES
      step_c0    = 1.0d-6
      step_cs    = 1.0d-6
      step_theta = 1.0d-6
      step_phi_1 = 1.0d-6
      step_phi_2 = 1.0d-6



! _______ NUMERICAL DERIVATIVES WRT c0, ci, phi_1, phi_2
! *** NOTES *** For c0, ci, it is unphysical for those values to be below 0
! For cases, where those values are close to zero we use a simple forward stepping numerical derivative
      if (c0.LE.step_c0) then
        di_rxn_dc0    = ( echem_rxn(c0+step_c0,cs,Phi_1, Phi_2) - echem_rxn(c0,cs,Phi_1, Phi_2) )/(step_c0)
      else
        di_rxn_dc0    = ( echem_rxn(c0+step_c0,cs,Phi_1, Phi_2) - echem_rxn(c0-step_c0,cs,Phi_1, Phi_2) )/(2.0*step_c0)
      end if

      if (cs.LE.step_cs) then
        di_rxn_dcs    = ( echem_rxn(c0,cs+step_cs,Phi_1, Phi_2) - echem_rxn(c0,cs,Phi_1, Phi_2) )/(step_cs)
      else
        di_rxn_dcs    = ( echem_rxn(c0,cs+step_cs,Phi_1, Phi_2) - echem_rxn(c0,cs-step_cs,Phi_1, Phi_2) )/(2.0*step_cs)
      end if
      
      di_rxn_dphi_1 = ( echem_rxn(c0,cs,Phi_1+step_phi_1, Phi_2) - echem_rxn(c0,cs,Phi_1-step_phi_1, Phi_2) )/(2.0*step_phi_1)

      di_rxn_dphi_2 = ( echem_rxn(c0,cs,Phi_1, Phi_2+step_phi_2) - echem_rxn(c0,cs,Phi_1, Phi_2-step_phi_2) )/(2.0*step_phi_2)

!       if ((l.EQ.1).AND.(j.EQ.NJ_c)) then
!           write(*,*) c0, cprev(1,l), Phi_1, Phi_2, echem_rxn(c0,cs,Phi_1, Phi_2), spec_a_agg
!       end if

! Turn off reactions if the theta_beta fraction is too large
!       if ( theta_beta.GE.0.97 ) then
!         i_rxn_agg  = 0.0
!         di_rxn_dc0 = 0.0
!         di_rxn_dcs = 0.0
!         di_rxn_dphi_1 = 0.0
!         di_rxn_dphi_2 = 0.0
!       end if


!!! PHASE CHANGE !!!

    if (cs.LT.4.0*molar_volume_Fe3O4) then
      drxn_pchange_dcs    = 0.0
      drxn_pchange_dtheta = 0.0

      rxn_pchange         = 0.0

    else
      drxn_pchange_dcs    = ( Rxn_Beta(cs + step_cs, theta_beta) - Rxn_Beta(cs - step_cs, theta_beta) ) / (2.0*step_cs)
      drxn_pchange_dtheta = ( Rxn_Beta(cs, theta_beta + step_theta) - Rxn_Beta(cs, theta_beta - step_theta) ) / (2.0*step_theta)

      rxn_pchange      = Rxn_Beta(cs, theta_beta)

    end if


!**** NOTES ******
      ! It is important to benchmark the numerical derivatives and have a good way to determine
      ! the appropriate step size that will produce the most accurate results 
      ! In addition this method of multiplying a certain fraction does not work well when the variable
      ! is close to zero
      ! Use a conditional statement when close to zero to use additive steps_agg
      ! if ( ABS(variable).LE.1.0d-10 ) then
      !     step  = 1.0d-10
      !     df_dv = (f(v+step) - f(v-step))/(2*step)
      ! end if

      ! for the concentrations, they are in log() terms and log(-#) is undefined
      ! for these variables it might be more useful to instead do something like:
      ! if ( variable.LE.1.0d-10 ) then
      !     step  = 1.0d-10
      !     df_dv = (f(2*step) - f(step))/(step)
      ! end if

!  19    format(1(A5,1X),  2(A14,1X),  12(A14,1X))
!  20    format(1(F5.2,1X),2(ES14.5,1X),12(ES14.5,1X))
 21    format(20(ES14.5,1X)) ! all scientific notation

!************************************************************************************
!*********** x = 0 ********************* x = 0 ********************* x = 0 **********
!************************************************************************************

!______________________Boundary-condition at x = 0
! At the boundaries it is usually useful to set the volume to zero (delX_c(j) = 0)
!     dE_c(i,k)*dcdxE + rj_c(i,k)*c = g(i)

      if (j.eq.1) then 

        alphaE_c=delx_c(j)/(delx_c(j+1)+delx_c(j))
        betaE_c=2.0/(delx_c(j)+delx_c(j+1))

        do ic=1,N_c
          cE_c(ic)    = alphaE_c*cprev_c(l,ic,j+1) + (1.d0 - alphaE_c)*cprev_c(l,ic,j)
          dcdxE_c(ic) = betaE_c * (cprev_c(l,ic,j+1) - cprev_c(l,ic,j))
        end do

! *************** Electrolyte Concentration ****************
! Both Fixed Concentration and No Flux work

!         rj_c(1,1) = 1.0
!         smG_c(1)  = cbulk - rj_c(1,1)*cprev_c(l,1,j)
!  Fixed Concentration Condition
!     c0 = cbulk at x = 0
        
        dE_c(1,1) = 1.d0
        smG_c(1) = 0.0 - dE_c(1,1)*dcdxE_c(1)


! *************** Solid-State Concentration ****************
! for cs both dcs/dt = rxn and dcs/dx = 0 produce the same results
! because the solid-state concentration does not have any spatial gradients, 
! spatial boundary conditions are also not necessary

! dc/dt = rxn
! dc/dt = -spec_a*i_rxn/Fconst
! dc/dt - (di_rxn_dc * delC + di_rxn_dcs * delCs + di_rxn_dphi_1 * delPhi_1)/Fconst = i_rxn/Fconst

        dW_c(4,4) = 0.0
        dE_c(4,4) = 0.0
        fW_c(4,4) = 0.0
        fE_c(4,4) = 0.0
        
        rj_c(4,1) = -spec_a_agg*di_rxn_dc0/Fconst
        rj_c(4,4) = -spec_a_agg*di_rxn_dcs/Fconst - (1.0-eps_agg)/delT*(1.0 - theta_beta)
!         rj_c(4,4) = -spec_a_agg*di_rxn_dcs/Fconst - (1.0)/delT*(1.0 - theta_beta)

        rj_c(4,5) = -(1.0-eps_agg)/delT*(gamma_conc - cs)
!         rj_c(4,5) = -(1.0)/delT*(gamma - cs)
        
        rj_c(4,2) = -spec_a_agg*di_rxn_dphi_1/Fconst
        rj_c(4,3) = -spec_a_agg*di_rxn_dphi_2/Fconst

        smG_c(4)  = +spec_a_agg*i_rxn_agg/Fconst



        rj_c(5,1) = 0.0
        rj_c(5,4) = drxn_pchange_dcs/(gamma_conc)
        rj_c(5,5) = -(1.0)/delT
        rj_c(5,2) = 0.0
        rj_c(5,3) = 0.0

        smG_c(5)  = -rxn_pchange/(gamma_conc)

! *************** Solid-State Potential ****************
! At x = 0 this is where the electrode meets the separator
! At this point there cannot be any flow of current through the solid-state

! -(1-eps_agg) * simga * dV/dx = 0 ! no electron flux at the separator
        
        dE_c(2,2) = -(1.0 - eps_agg)*sigma
        smG_c(2)  =  0.0 - dE_c(2,2)*dcdxE_c(2)    

! *************** Solution Potential ****************    
! At x = 0 this is where the electrode meets the separator
! At this point there cannot be any flow of current through the solid-state
! Because there is no electronic current, all the current must be ionic
! **  -eps_agg * kappa * dPhi_2/dx = c_applied **
! However for some numerical reason this boundary condition is not working so instead
! we are assuming that Phi_2 = 0.0 at x = 0
! Note that this condition still leads to the condition -eps_agg * kappa * dPhi_2/dx = c_applied 
! being met

        dE_c(3,3) = 1.0
        smG_c(3)  = 0.0 - dE_c(3,3)*dcdxE_c(3)

        return
      end if

!*********************************************************************************
!*********** XMAX ********************* XMAX ********************* XMAX **********
!*********************************************************************************

!______________________Boundary-condition at x=xmax
! At the boundaries it is usually useful to set the volume to zero (delX_c(j) = 0)
! The general boundary condition equation looks like: dW_c(i,k)*dcdxW + rj_c(i,k)*delC = g(i)

      if (j.eq.NJ_c) then

        alphaW_c=delx_c(j-1)/(delx_c(j-1)+delx_c(j))
        betaW_c=2.0/(delx_c(j-1)+delx_c(j))

        do ic=1,N_c
          cW_c(ic)    = alphaW_c*cprev_c(l,ic,j) + (1.d0 - alphaW_c)*cprev_c(l,ic,j-1)
          dcdxW_c(ic) = betaW_c * (cprev_c(l,ic,j) - cprev_c(l,ic,j-1))
        end do

! *************** Electrolyte Concentration ****************
! at x = xmax this is where the electrode meets the current collector
! at this point there cannot be any flow of electrolyte through the solid electrode
! Therefore there is a no flux condition at xmax

!  No flux Condition
!     dc0/dx = 0 at x = xmax

!         dW_c(1,1) = 1.0
!         smG_c(1)  = 0.0 - dW_c(1,1)*dcdxW_c(1)

! Constant Concentration
        rj_c(1,1) = 1.0
        smG_c(1)  = c0_edge - rj_c(1,1)*cprev_c(l,1,j)

! *************** Solid-State Concentration ****************
! No Boundary conditions on cs
! dc/dt = rxn
! dc/dt = i_rxn/Fconst
! dc/dt - (di_rxn_dc * delC + di_rxn_dcs * delCs + di_rxn_dphi_1 * delPhi_1)/Fconst = i_rxn/Fconst

        dW_c(4,4) = 0.0
        dE_c(4,4) = 0.0
        fW_c(4,4) = 0.0
        fE_c(4,4) = 0.0
        
        rj_c(4,1) = -spec_a_agg*di_rxn_dc0/Fconst
        rj_c(4,4) = -spec_a_agg*di_rxn_dcs/Fconst - (1.0-eps_agg)/delT*(1.0 - theta_beta)
!         rj_c(4,4) = -spec_a_agg*di_rxn_dcs/Fconst - (1.0)/delT*(1.0 - theta_beta)

        rj_c(4,5) = -(1.0-eps_agg)/delT*(gamma_conc - cs)
!         rj_c(4,5) = -(1.0)/delT*(gamma - cs)
        
        rj_c(4,2) = -spec_a_agg*di_rxn_dphi_1/Fconst
        rj_c(4,3) = -spec_a_agg*di_rxn_dphi_2/Fconst

        smG_c(4)  = +spec_a_agg*i_rxn_agg/Fconst



        rj_c(5,1) = 0.0
        rj_c(5,4) = drxn_pchange_dcs/(gamma_conc)
        rj_c(5,5) = -(1.0)/delT
        rj_c(5,2) = 0.0
        rj_c(5,3) = 0.0

        smG_c(5)  = -rxn_pchange/(gamma_conc)

! *************** Solid-State Potential ****************
! at x = xmax this is where the electrode meets the current collector
! at this point there cannot be any flow of current in the electrolyte
! ie all the current must be in the solid-state (or electronic)

! -(1 - eps_agg) * sigma * dV/dx = c_density at x = xmax

! This boundary condition changes during the experiment:
! Discharge, recovery, charge

! *** NOTE ***
! c_density is defined as a positive (+) number is this code
! in the agglomerate code it was defined as a negative (-) number

!       if (time.LE.Discharge_Time) then                                     ! Lithiate
!         dW_c(2,2) = -(1.0-eps_agg)*sigma
!         smG_c(2)  = c_density_agg - dW_c(2,2)*dcdxW_c(2)
!         state = 'D'
! !         mAhg  = mAhg + 1000.0*c_specific*delT/3600.0

!       else if (time.LE.(Discharge_Time + Discharge_Relax_Time)) then        ! Recover
!         dW_c(2,2) = -(1.0-eps_agg)*sigma
!         smG_c(2)  = 0.0 - dW_c(2,2)*dcdxW_c(2)
!         state = 'R'
! !         mAhg  = mAhg 
!                                                                             ! Delithiate
!       else if ( (time.LE.(Discharge_Time+Discharge_Relax_Time+Charge_Time)) .AND. (cprev_c(l,2,NJ_c).LE.3.0) ) then ! Delithiate
!         dW_c(2,2) = -(1.0-eps_agg)*sigma
!         smG_c(2)  = -c_density_agg - dW_c(2,2)*dcdxW_c(2)
!         state = 'C'
! !         mAhg  = mAhg - 1000.0*c_specific*delT/3600.0

!       else if (Relax.EQ.1) then                                             ! Recover
!         dW_c(2,2) = -(1.0-eps_agg)*sigma
!         smG_c(2)  = 0.0 - dW_c(2,2)*dcdxW_c(2)
!         state = 'R'
! !         mAhg  = mAhg

!       else                                                                   ! Recover
!         dW_c(2,2) = -(1.0-eps_agg)*sigma
!         smG_c(2)  = 0.0 - dW_c(2,2)*dcdxW_c(2)
!         state = 'R'
! !         mAhg  = mAhg  

!       end if


      dW_c(2,2) = -(1.0-eps_agg)*sigma
      smG_c(2)  = c_density_agg - dW_c(2,2)*dcdxW_c(2)

! *************** Solution Potential ****************
! at x = xmax this is where the electrode meets the current collector
! at this point there cannot be any flow of current in the electrolyte
! Therefore a no flux condition is met for ionic current

! -eps_agg * kappa * dPhi_2/dx = 0 at x = xmax

!         dW_c(4,4) = -eps_agg*kappa
!         smG_c(4) = 0.0 - dW_c(4,4)*dcdxW_c(4)

! Phi_2 = 0 at x = xmax
        rj_c(3,3) = 1.0
        smG_c(3)  = Phi_2_elect - rj_c(3,3)*cprev_c(l,3,j)
           

        return
      end if


! !______________________Governing Equations
! Finite Difference Method (Control Volume Formulation)
! A material-balance for species i at node k

!       dc(i,k)/dt * delX_c(k) = Flux_west(i) - Flux_east(i) + R(i)*delX_c(k)

! where R(i) is the production rate per unit volume, Flux_west(i) is the flux of species i
! at the "west" face of node k, and Flux_east(i) is the flux of species i at the "east" face 
! of node k.

! In general the flux of species i can be affected by all of the variable and their gradients
      
!       Flux_west(i) = d_west(j,k) * dcdx_west(j) + f_west*c_west(j)
!       Flux_east(i) = d_east(j,k) * dcdx_east(j) + f_east*c_east(j)

! The general formula to apply to each node point can be written as:

!      g(i) = rj_c(i,k)*c + fW(i,k)*cW + dW_c(i,k)*dcdxW - fE(i,k)*cE - dE_c(i,k)*dcdxE

! The optimal interpolation formula to use to evaluate the variable c and their derivatives
! depend on the local Peclet number. The general form can be seen below
! For a low Peclet number, a linear interpolation is most appropriate
! Linear interpolation, alpha = 0.5, beta = 1/delX

! *** NOTE *** For large Peclet numbers this central-difference approximation causes an oscillatory
! behavior in the concentration. This can be eliminated by an upwind scheme, in which betaW_c and betaE_c
! remain the same, but 
! alphaE_c = 0, alphaW_c = 0 when flow is from west to east
! alphaE_c = 1, alphaW_c = 1 when flow is from east to west

        alphaW_c=delx_c(j-1)/(delx_c(j-1)+delx_c(j))
        alphaE_c=delx_c(j)/(delx_c(j+1)+delx_c(j))
        betaW_c=2.0/(delx_c(j-1)+delx_c(j))
        betaE_c=2.0/(delx_c(j)+delx_c(j+1))


        do ic=1,N_c
          cW_c(ic)    = alphaW_c*cprev_c(l,ic,j) + (1.d0 - alphaW_c)*cprev_c(l,ic,j-1)
          cE_c(ic)    = alphaE_c*cprev_c(l,ic,j+1) + (1.d0 - alphaE_c)*cprev_c(l,ic,j)
          dcdxW_c(ic) = betaW_c * (cprev_c(l,ic,j) - cprev_c(l,ic,j-1))
          dcdxE_c(ic) = betaE_c * (cprev_c(l,ic,j+1) - cprev_c(l,ic,j))
        end do

        rj_c(1,1) = 1.0
        smG_c(1)  = cbulk - rj_c(1,1)*cprev_c(l,1,j)
! *************** Electrolyte Concentration ****************
! dc/dt = D * d2c/dx2 + rxn
! delX * dc/dt = -D * dcdx_W - (-D * dcdx_E) - i_rxn/Fconst * delX
! delX * dc/dt + (di_rxn_dc * delC + di_rxn_dcs * delCs + di_rxn_dphi_1 * delPhi_1 + di_rxn_dphi_2 * delPhi_2)/Fconst * delX
!            = -D * dcdx_W - (-D * dcdx_E) - i_rxn/Fconst*delX

   ! GEOMETRY
        r_W = xx_c(j) - delx_c(j)/2.0
        r_E = xx_c(j) + delx_c(j)/2.0

        ! Rectangular
!         A_W = 1.0
!         A_E = 1.0
!         delV = r_E - r_W

!        Cylindrical
!         A_W = 2.0*PI*r_W
!         A_E = 2.0*PI*r_E
!         delV = PI*(r_E**2.0 - r_W**2.0)

!         ! Spherical
        A_W  = 4.0*PI*r_W**2.0
        A_E  = 4.0*PI*r_E**2.0
        delV = 4.0*PI/3.0*(r_E**3.0 - r_W**3.0)


        dW_c(1,1) = A_W*(-eps_agg*diff_cat)
        dE_c(1,1) = A_E*(-eps_agg*diff_cat)
        fW_c(1,1) = 0.d0
        fE_c(1,1) = 0.d0

        rj_c(1,1) = spec_a_agg*di_rxn_dc0/Fconst*delV - eps_agg/delT*delV
        rj_c(1,4) = spec_a_agg*di_rxn_dcs/Fconst*delV
        rj_c(1,2) = spec_a_agg*di_rxn_dphi_1/Fconst*delV
        rj_c(1,3) = spec_a_agg*di_rxn_dphi_2/Fconst*delV

        smG_c(1)  = -spec_a_agg*i_rxn_agg/Fconst*delV &
        &         - (fW_c(1,1)*cW_c(1) + dW_c(1,1)*dcdxW_c(1)) &
        &         + (fE_c(1,1)*cE_c(1) + dE_c(1,1)*dcdxE_c(1))


! *************** Solid-State Concentration ****************
! dc/dt = -rxn
! dc/dt = i_rxn/Fconst
! dc/dt - (di_rxn_dc * delC + di_rxn_dcs * delCs + di_rxn_dphi_1 * delPhi_1 + di_rxn_dphi_2 * delPhi_2)/Fconst = i_rxn/Fconst
        
        dW_c(4,4) = 0.0
        dE_c(4,4) = 0.0
        fW_c(4,4) = 0.0
        fE_c(4,4) = 0.0
        
        rj_c(4,1) = -spec_a_agg*di_rxn_dc0/Fconst
        rj_c(4,4) = -spec_a_agg*di_rxn_dcs/Fconst - (1.0-eps_agg)/delT*(1.0 - theta_beta)
!         rj_c(4,4) = -spec_a_agg*di_rxn_dcs/Fconst - (1.0)/delT*(1.0 - theta_beta)

        rj_c(4,5) = -(1.0-eps_agg)/delT*(gamma_conc - cs)
!         rj_c(4,5) = -(1.0)/delT*(gamma - cs)
        
        rj_c(4,2) = -spec_a_agg*di_rxn_dphi_1/Fconst
        rj_c(4,3) = -spec_a_agg*di_rxn_dphi_2/Fconst

        smG_c(4)  = +spec_a_agg*i_rxn_agg/Fconst



        rj_c(5,1) = 0.0
        rj_c(5,4) = drxn_pchange_dcs/(gamma_conc)
        rj_c(5,5) = -(1.0)/delT
        rj_c(5,2) = 0.0
        rj_c(5,3) = 0.0

        smG_c(5)  = -rxn_pchange/(gamma_conc)

! *************** Solid-State Potential ****************
! It is assumed that electrons travel very quickly so there is no accumulation of electrons
! (Kirchhoff's laws also assume no accumulation of current)
! dce/dt = sigma * d2V/dx2 + i_rxn 
!      0 = sigma * d2V/dx2 + i_rxn 
!  - delX*(di_rxn_dc * delC + di_rxn_dcs * delCs + di_rxn_dphi_1 * delPhi_1 + di_rxn_dphi_2 * delPhi_2) = -sigma * dcdx_W - (-sigma * dcdx_E) + delX*i_rxn

        dW_c(2,2) = A_W*(-(1.0-eps_agg)*sigma)
        dE_c(2,2) = A_E*(-(1.0-eps_agg)*sigma)
        fW_c(2,2) = 0.0
        fE_c(2,2) = 0.0

        rj_c(2,1) = -spec_a_agg*di_rxn_dc0*delV
        rj_c(2,4) = -spec_a_agg*di_rxn_dcs*delV
        rj_c(2,2) = -spec_a_agg*di_rxn_dphi_1*delV
        rj_c(2,3) = -spec_a_agg*di_rxn_dphi_2*delV

        smG_c(2)  = +spec_a_agg*i_rxn_agg*delV & 
        &         - (fW_c(2,2)*cW_c(2) + dW_c(2,2)*dcdxW_c(2)) &
        &         + (fE_c(2,2)*cE_c(2) + dE_c(2,2)*dcdxE_c(2))

! *************** Solution Potential **************** 
! Similar to what is assumed in the solid-state, in the electrolyte it is also assumed that ions do not accumulate 
! (also called conservation of charge) therefore dPhi_2/dt = 0 everywhere
! 0 = -eps_agg * kappa * d2Phi_2/dx2 - i_rxn
! + delX*(di_rxn_dc * delC + di_rxn_dcs * delCs + di_rxn_dphi_1 * delPhi_1 + di_rxn_dphi_2 * delPhi_2) = -esp * kappa * dcdx_W - (-eps_agg * kappa * dcdx_E) - delX*i_rxn       

!       kappa = Fconst**2*(z_cat**2*u_cat*conc_cat + z_an**2*u_an*conc_an)
!       kappa = Fconst**2*(z_cat**2*u_cat + z_an**2*u_an)*conc
        
        dW_c(3,1) = A_W*(-eps_agg*Fconst*(z_cat*diff_cat + z_an*diff_an))
        dE_c(3,1) = A_E*(-eps_agg*Fconst*(z_cat*diff_cat + z_an*diff_an))
        fW_c(3,1) = A_W*(-eps_agg*Fconst**2*(z_cat**2*u_cat + z_an**2*u_an)*dcdxW_c(3))
        fE_c(3,1) = A_E*(-eps_agg*Fconst**2*(z_cat**2*u_cat + z_an**2*u_an)*dcdxE_c(3))
        
        dW_c(3,3) = A_W*(-eps_agg*Fconst**2*(z_cat**2*u_cat + z_an**2*u_an)*cW_c(1))
        dE_c(3,3) = A_E*(-eps_agg*Fconst**2*(z_cat**2*u_cat + z_an**2*u_an)*cE_c(1))
        fW_c(3,3) = A_W*(0.0)
        fE_c(3,3) = A_E*(0.0)

        rj_c(3,1) = +spec_a_agg*di_rxn_dc0*delV
        rj_c(3,4) = +spec_a_agg*di_rxn_dcs*delV
        rj_c(3,2) = +spec_a_agg*di_rxn_dphi_1*delV
        rj_c(3,3) = +spec_a_agg*di_rxn_dphi_2*delV

        smG_c(3)  = -spec_a_agg*i_rxn_agg*delV &
        &         - (fW_c(3,3)*cW_c(3) + dW_c(3,3)*dcdxW_c(3)) &
        &         + (fE_c(3,3)*cE_c(3) + dE_c(3,3)*dcdxE_c(3))

      return

end subroutine fillmat_c



! Below this point in code are the fundamental subroutines. Do not edit anything below.
!************************************************************************************
!************************************************************************************
!************************************************************************************


!************************************ABDGXY******************************************

subroutine ABDGXY(h,j)
      use user_input
      use variables
      implicit double precision(a-h,o-z)

      if(j.eq.1) then
          do 1 ii=1,n
          do 10 kk=1,n
             X(ii,kk)=0.d0
             B(ii,kk)=rj(ii,kk) - (1.d0 - alphaE)*fE(ii,kk) + betaE*dE(ii,kk)
             D(ii,kk)= -alphaE*fE(ii,kk) - betaE*dE(ii,kk)

10         continue
             G(ii)=smG(ii)
1        continue
          return
      end if
      if (j.eq.NJ) then 
          do 2 ii=1,n
          do 20 kk=1,n
             Y(ii,kk)=0.d0
             A(ii,kk)=(1.d0 - alphaW)*fW(ii,kk) - betaW*dW(ii,kk)
             B(ii,kk)=rj(ii,kk) + betaW*dW(ii,kk) + alphaW*fW(ii,kk)
20         continue
             G(ii)=smG(ii)
2         continue
          return
      end if
      do 3 ii=1,n
      do 30 kk=1,n
             A(ii,kk)=(1.d0 - alphaW)*fW(ii,kk) - betaW*dW(ii,kk)
             B(ii,kk)=rj(ii,kk) + betaW*dW(ii,kk) + alphaW*fW(ii,kk) &
             &        - (1.d0 - alphaE)*fE(ii,kk) + betaE*dE(ii,kk)
             D(ii,kk)= -alphaE*fE(ii,kk) - betaE*dE(ii,kk)
30     continue
             G(ii)=smG(ii)
3     continue
      return
end subroutine ABDGXY

!***********************************MATINV*****************************************

SUBROUTINE MATINV(N,M,DETERM)     
use variables, only: A,B,delC,D,ID
implicit double precision (A-H,O-Z)    

      DETERM=1.0
      DO 1 I=1,N
1     ID(I)=0
      DO 18 NN=1,N
      BMAX=1.1
      DO 6 I=1,N
      IF (ID(I).NE.0) GOTO 6
      BNEXT=0.0
      BTRY=0.0
      DO 5 J=1,N
      IF (ID(J).NE.0) GOTO 5
      IF (DABS(B(I,J)).LE.BNEXT) GOTO 5
      BNEXT=DABS(B(I,J))
      IF (BNEXT.LE.BTRY) GOTO 5
      BNEXT=BTRY
      BTRY=DABS(B(I,J))
      JC=J
5     CONTINUE 
      IF (BNEXT.GE.BMAX*BTRY) GOTO 6
      BMAX=BNEXT/BTRY
      IROW=I
      JCOL=JC 
6     CONTINUE
      IF (ID(JC).EQ.0) GOTO 8
      DETERM=0.0
      RETURN
8     ID(JCOL)=1
      IF (JCOL.EQ.IROW) GOTO 12
9     DO 10 J=1,N
      SAVE=B(IROW,J)
      B(IROW,J)=B(JCOL,J)
10    B(JCOL,J)=SAVE
      DO 11 K=1,M
      SAVE=D(IROW,K)
      D(IROW,K)=D(JCOL,K)
11    D(JCOL,K)=SAVE
12    F=1.0/B(JCOL,JCOL)
      DO 13 J=1,N
13    B(JCOL,J)=B(JCOL,J)*F
      DO 14 K=1,M
14    D(JCOL,K)=D(JCOL,K)*F
      DO 18 I=1,N
      IF (I.EQ.JCOL) GOTO 18
      F=B(I,JCOL)
      DO 16 J=1,N
16    B(I,J)=B(I,J)-F*B(JCOL,J)
      DO 17 K=1,M
17    D(I,K)=D(I,K)-F*D(JCOL,K)
18    CONTINUE
      RETURN
      END

!*************************************BAND******************************************

SUBROUTINE BAND(J)
use variables, only: A,B,delC,D,G,X,Y,NP1,E
use user_input, only: N,NJ
implicit double precision (A-H,O-Z) 


101   FORMAT(15H DETERM=0 AT J=,I4)
      IF (J-2) 1,6,8
1     NP1=N+1
      DO 2 I=1,N
      D(I,2*N+1)=G(I)
      DO 2 L=1,N
      LPN=L+N
2     D(I,LPN)=X(I,L)
      CALL MATINV(N,2*N+1,DETERM)
      IF (DETERM) 4,3,4
3     PRINT 101,J
4     DO 5 K=1,N
      E(K,NP1,1)=D(K,2*N+1)
      DO 5 L=1,N
      E(K,L,1)=-D(K,L)
      LPN=L+N
5     X(K,L)=-D(K,LPN)
      RETURN
6     DO 7 I=1,N
      DO 7 K=1,N
      DO 7 L=1,N
7     D(I,K)=D(I,K)+A(I,L)*X(L,K)
8     IF (J-NJ) 11,9,9
9     DO 10 I=1,N
      DO 10 L=1,N
      G(I)=G(I)-Y(I,L)*E(L,NP1,J-2)
      DO 10 M=1,N
10    A(I,L)=A(I,L) + Y(I,M)*E(M,L,J-2)
11    DO 12 I=1,N
      D(I,NP1)=-G(I)
      DO 12 L=1,N
      D(I,NP1)=D(I,NP1)+A(I,L)*E(L,NP1,J-1)
      DO 12 K=1,N
12    B(I,K)=B(I,K) + A(I,L)*E(L,K,J-1)
      CALL MATINV(N,NP1,DETERM)
      IF (DETERM) 14,13,14
13    PRINT 101,J
14    DO 15 K=1,N
      DO 15 M=1,NP1
15    E(K,M,J)=-D(K,M)
      IF (J-NJ) 20,16,16
16    DO 17 K=1,N
17    delC(K,J)=E(K,NP1,J)
      DO 18 JJ=2,NJ
      M=NJ-JJ+1
      DO 18 K=1,N
      delc(K,M)=E(K,NP1,M)
      DO 18 L=1,N
18    delC(K,M)=delC(K,M) +E(K,L,M)*delC(L,M+1)  
      DO 19 L=1,N
      DO 19 K=1,N
19    delC(K,1)=delC(K,1)+X(K,L)*delC(L,3)
20    RETURN
      END

!************************************ Crystal Scale Operations ******************************************
!************************************ Crystal Scale Operations ******************************************
!************************************ Crystal Scale Operations ******************************************

!************************************ABDGXY_c******************************************

subroutine ABDGXY_c(h_c,j)
      use user_input
      use variables
      implicit double precision(a-h,o-z)

      if(j.eq.1) then
          do 1 ii=1,n_c
          do 10 kk=1,n_c
             X_c(ii,kk)=0.d0
             B_c(ii,kk)=rj_c(ii,kk) - (1.d0 - alphaE_c)*fE_c(ii,kk) + betaE_c*dE_c(ii,kk)
             D_c(ii,kk)= -alphaE_c*fE_c(ii,kk) - betaE_c*dE_c(ii,kk)

10         continue
             G_c(ii)=smG_c(ii)
1        continue
          return
      end if
      if (j.eq.NJ_c) then 
          do 2 ii=1,n_c
          do 20 kk=1,n_c
             Y_c(ii,kk)=0.d0
             A_c(ii,kk)=(1.d0 - alphaW_c)*fW_c(ii,kk) - betaW_c*dW_c(ii,kk)
             B_c(ii,kk)=rj_c(ii,kk) + betaW_c*dW_c(ii,kk) + alphaW_c*fW_c(ii,kk)
20         continue
             G_c(ii)=smG_c(ii)
2         continue
          return
      end if
      do 3 ii=1,n_c
      do 30 kk=1,n_c
             A_c(ii,kk)=(1.d0 - alphaW_c)*fW_c(ii,kk) - betaW_c*dW_c(ii,kk)
             B_c(ii,kk)=rj_c(ii,kk) + betaW_c*dW_c(ii,kk) + alphaW_c*fW_c(ii,kk) &
             &        - (1.d0 - alphaE_c)*fE_c(ii,kk) + betaE_c*dE_c(ii,kk)
             D_c(ii,kk)= -alphaE_c*fE_c(ii,kk) - betaE_c*dE_c(ii,kk)
30     continue
             G_c(ii)=smG_c(ii)
3     continue
      return
end subroutine ABDGXY_c

!***********************************MATINV_c*****************************************

SUBROUTINE MATINV_c(N_c,M,DETERM)     
 use variables, only: A_c,B_c,delC_c,D_c,ID_c ! A_c imported but not used
 implicit double precision (A-H,O-Z)    

      DETERM=1.0
      DO 1 I=1,N_c
1     ID_c(I)=0
      DO 18 NN=1,N_c
      BMAX=1.1
      DO 6 I=1,N_c
      IF (ID_c(I).NE.0) GOTO 6
      BNEXT=0.0
      BTRY=0.0
      DO 5 J=1,N_c
      IF (ID_c(J).NE.0) GOTO 5
      IF (DABS(B_c(I,J)).LE.BNEXT) GOTO 5
      BNEXT=DABS(B_c(I,J))
      IF (BNEXT.LE.BTRY) GOTO 5
      BNEXT=BTRY
      BTRY=DABS(B_c(I,J))
      JC=J
5     CONTINUE 
      IF (BNEXT.GE.BMAX*BTRY) GOTO 6
      BMAX=BNEXT/BTRY
      IROW=I
      JCOL=JC 
6     CONTINUE
      IF (ID_c(JC).EQ.0) GOTO 8
      DETERM=0.0
      RETURN
8     ID_c(JCOL)=1
      IF (JCOL.EQ.IROW) GOTO 12
9     DO 10 J=1,N_c
      SAVE=B_c(IROW,J)
      B_c(IROW,J)=B_c(JCOL,J)
10    B_c(JCOL,J)=SAVE
      DO 11 K=1,M
      SAVE=D_c(IROW,K)
      D_c(IROW,K)=D_c(JCOL,K)
11    D_c(JCOL,K)=SAVE
12    F=1.0/B_c(JCOL,JCOL)
      DO 13 J=1,N_c
13    B_c(JCOL,J)=B_c(JCOL,J)*F
      DO 14 K=1,M
14    D_c(JCOL,K)=D_c(JCOL,K)*F
      DO 18 I=1,N_c
      IF (I.EQ.JCOL) GOTO 18
      F=B_c(I,JCOL)
      DO 16 J=1,N_c
16    B_c(I,J)=B_c(I,J)-F*B_c(JCOL,J)
      DO 17 K=1,M
17    D_c(I,K)=D_c(I,K)-F*D_c(JCOL,K)
18    CONTINUE
      RETURN
      END

!*************************************BAND_c******************************************

SUBROUTINE BAND_c(J)
use variables, only: A_c,B_c,delC_c,D_c,G_c,X_c,Y_c,NP1_c,E_c
use user_input, only: N_c,NJ_c
implicit double precision (A-H,O-Z) 


101   FORMAT(15H DETERM=0 AT J=,I4)
      IF (J-2) 1,6,8
  1     NP1_c=N_c+1
        DO 2 I=1,N_c
        D_c(I,2*N_c+1)=G_c(I)
        DO 2 L=1,N_c
        LPN=L+N_c
  2     D_c(I,LPN)=X_c(I,L)
        CALL MATINV_c(N_c,2*N_c+1,DETERM)
        IF (DETERM) 4,3,4
  3     PRINT 101,J
  4     DO 5 K=1,N_c
        E_c(K,NP1_c,1)=D_c(K,2*N_c+1)
        DO 5 L=1,N_c
        E_c(K,L,1)=-D_c(K,L)
        LPN=L+N_c
  5     X_c(K,L)=-D_c(K,LPN)
        RETURN

  6     DO 7 I=1,N_c
        DO 7 K=1,N_c
        DO 7 L=1,N_c
  7     D_c(I,K)=D_c(I,K)+A_c(I,L)*X_c(L,K)


8     IF (J-NJ_c) 11,9,9
9     DO 10 I=1,N_c
      DO 10 L=1,N_c
      G_c(I)=G_c(I)-Y_c(I,L)*E_c(L,NP1_c,J-2)
      DO 10 M=1,N_c
10    A_c(I,L)=A_c(I,L) + Y_c(I,M)*E_c(M,L,J-2)
11    DO 12 I=1,N_c
      D_c(I,NP1_c)=-G_c(I)
      DO 12 L=1,N_c
      D_c(I,NP1_c)=D_c(I,NP1_c)+A_c(I,L)*E_c(L,NP1_c,J-1)
      DO 12 K=1,N_c
12    B_c(I,K)=B_c(I,K) + A_c(I,L)*E_c(L,K,J-1)
      CALL MATINV_c(N_c,NP1_c,DETERM)
      IF (DETERM) 14,13,14
13    PRINT 101,J
14    DO 15 K=1,N_c
      DO 15 M=1,NP1_c
15    E_c(K,M,J)=-D_c(K,M)
      IF (J-NJ_c) 20,16,16
16    DO 17 K=1,N_c
17    delC_c(K,J)=E_c(K,NP1_c,J)
      DO 18 JJ=2,NJ_c
      M=NJ_c-JJ+1
      DO 18 K=1,N_c
      delC_c(K,M)=E_c(K,NP1_c,M)
      DO 18 L=1,N_c
18    delC_c(K,M)=delC_c(K,M) +E_c(K,L,M)*delC_c(L,M+1)  
      DO 19 L=1,N_c
      DO 19 K=1,N_c
19    delC_c(K,1)=delC_c(K,1)+X_c(K,L)*delC_c(L,3)
20    RETURN
      
      END
