! ------------------------------------------------------------------------------
! . last update @20240521
! ------------------------------------------------------------------------------
module coordinate
! ------------------------------------------------------------------------------
  implicit none
  double precision,parameter::pi=3.14159265358979323846d0
  integer::irho
  integer,parameter::irho_min=0
  integer,parameter::irho_max=1024 ! 256
  double precision::rho(irho_min:irho_max)
  double precision,parameter::rho_min=0.1d0
  double precision,parameter::rho_max=1.0d0
  double precision,parameter::drho=(rho_max-rho_min)/dble(irho_max-irho_min)
  double precision,parameter::drho_sq=drho**2
  integer::n
  integer::ialpha
  integer,parameter::n_max=320 ! 40
  integer,parameter::ialpha_min=0
  integer,parameter::ialpha_max=1024 ! 128 ! : >3*n_max
  integer,parameter::numdata_alpha=ialpha_max-ialpha_min
  double precision,parameter::inumdata_alpha=1.0d0/dble(numdata_alpha)
  integer::l
  integer::iz
  integer,parameter::l_max=10
  integer,parameter::l_min=-l_max
  integer,parameter::iz_min=0
  integer,parameter::iz_max=32 ! : >3*l_max
  integer,parameter::numdata_z=iz_max-iz_min
  double precision::z(0:numdata_z-1)
  double precision,parameter::inumdata_z=1.0d0/dble(numdata_z)
  integer::itheta
  integer,parameter::itheta_min=0
  integer,parameter::itheta_max=4096 ! 512 ! : >3*m_res
  integer,parameter::numdata_theta=itheta_max-itheta_min
  double precision::theta(0:numdata_theta)
  double precision,parameter::inumdata_theta=1.0d0/dble(numdata_theta)
end module coordinate
! ------------------------------------------------------------------------------
module time
! ------------------------------------------------------------------------------
  implicit none
  double precision::t
  integer::it
  double precision,parameter::dt=0.005d0 ! 0.01d0 ! : <dt_crit
  integer,parameter::irestart=0
  integer,parameter::it_span=100000 ! 50000
  integer,parameter::it_start=it_span*irestart
  integer,parameter::it_end=it_span*(irestart+1)
  integer,parameter::it_skip=it_span/50
  integer,parameter::it_mskip=it_span/500
  character(len=100)::filename
  double precision::dt_rk(1:4)
end module time
! ------------------------------------------------------------------------------
module parameter
! ------------------------------------------------------------------------------
  implicit none
  double precision,parameter::eps=0.3d0
  double precision,parameter::beta=0.01d0
  double precision,parameter::mu=1.0d-5/2.0d0 ! 1.0d-5
  double precision,parameter::eta=1.0d-5/2.0d0 ! 1.0d-5
  double precision,parameter::chi=1.0d-6/2.0d0 ! 1.0d-5
end module parameter
! ------------------------------------------------------------------------------
module field
! ------------------------------------------------------------------------------
  use coordinate
  implicit none
  ! . equilibria
  double precision::U_eq(irho_min:irho_max)
  double precision::phi_eq(irho_min:irho_max)
  double precision::q(irho_min:irho_max)
  double precision::j_eq(irho_min:irho_max)
  double precision::p_eq(irho_min:irho_max)
  double precision::dU_eq(irho_min:irho_max)
  double precision::dphi_eq(irho_min:irho_max)
  double precision::dq(irho_min:irho_max)
  double precision::dj_eq(irho_min:irho_max)
  double precision::dp_eq(irho_min:irho_max)
  double precision::d2q(irho_min:irho_max)
  ! ..
  double precision,allocatable::l_shift(:,:)
  double precision,allocatable::k_theta(:,:)
  double precision,allocatable::k_para(:,:,:)
  double precision::fact_del_rho_bra_x(0:numdata_z-1)
  complex(kind(0d0)),allocatable::fact_bra_x(:,:,:)
  complex(kind(0d0)),allocatable::fact_del_rho_lap_perp(:,:,:)
  complex(kind(0d0)),allocatable::fact_lap_perp(:,:,:)
  complex(kind(0d0)),allocatable::ishifter(:,:,:)
  complex(kind(0d0)),allocatable::shifter(:,:,:)
  ! . spectra
  complex(kind(0d0))::U(irho_min:irho_max)
  complex(kind(0d0))::phi(irho_min:irho_max)
  complex(kind(0d0))::psi(irho_min:irho_max)
  complex(kind(0d0))::j(irho_min:irho_max)
  complex(kind(0d0))::p(irho_min:irho_max)
  complex(kind(0d0)),allocatable::U_nl(:,:,:)
  complex(kind(0d0)),allocatable::phi_nl(:,:,:)
  complex(kind(0d0)),allocatable::psi_nl(:,:,:)
  complex(kind(0d0)),allocatable::j_nl(:,:,:)
  complex(kind(0d0)),allocatable::p_nl(:,:,:)
  complex(kind(0d0)),allocatable::U_n(:,:,:)
  complex(kind(0d0)),allocatable::phi_n(:,:,:)
  complex(kind(0d0)),allocatable::psi_n(:,:,:)
  complex(kind(0d0)),allocatable::j_n(:,:,:)
  complex(kind(0d0)),allocatable::p_n(:,:,:)
  ! ..
  complex(kind(0d0))::U_old(irho_min:irho_max)
  complex(kind(0d0))::phi_old(irho_min:irho_max)
  complex(kind(0d0))::psi_old(irho_min:irho_max)
  complex(kind(0d0))::j_old(irho_min:irho_max)
  complex(kind(0d0))::p_old(irho_min:irho_max)
  complex(kind(0d0)),allocatable::U_old_nl(:,:,:)
  complex(kind(0d0)),allocatable::phi_old_nl(:,:,:)
  complex(kind(0d0)),allocatable::psi_old_nl(:,:,:)
  complex(kind(0d0)),allocatable::j_old_nl(:,:,:)
  complex(kind(0d0)),allocatable::p_old_nl(:,:,:)
  ! . linear terms
  complex(kind(0d0))::bra_phi_U_lin(irho_min:irho_max)
  complex(kind(0d0))::nab_para_j_lin(irho_min:irho_max)
  complex(kind(0d0))::bra_x_p(irho_min:irho_max)
  complex(kind(0d0))::lap_perp_U(irho_min:irho_max)
  complex(kind(0d0))::nab_para_phi_lin(irho_min:irho_max)
  complex(kind(0d0))::bra_phi_p_lin(irho_min:irho_max)
  complex(kind(0d0))::bra_x_phi(irho_min:irho_max)
  complex(kind(0d0))::lap_perp_p(irho_min:irho_max)
  complex(kind(0d0)),allocatable::bra_phi_U_lin_nl(:,:,:)
  complex(kind(0d0)),allocatable::nab_para_j_lin_nl(:,:,:)
  complex(kind(0d0)),allocatable::bra_x_p_nl(:,:,:)
  complex(kind(0d0)),allocatable::lap_perp_U_nl(:,:,:)
  complex(kind(0d0)),allocatable::nab_para_phi_lin_nl(:,:,:)
  complex(kind(0d0)),allocatable::bra_phi_p_lin_nl(:,:,:)
  complex(kind(0d0)),allocatable::bra_x_phi_nl(:,:,:)
  complex(kind(0d0)),allocatable::lap_perp_p_nl(:,:,:)
  complex(kind(0d0)),allocatable::bra_x_p_n(:,:,:)
  complex(kind(0d0)),allocatable::lap_perp_U_n(:,:,:)
  complex(kind(0d0)),allocatable::bra_x_phi_n(:,:,:)
  complex(kind(0d0)),allocatable::lap_perp_p_n(:,:,:)
  ! . nonlinear terms
  complex(kind(0d0))::bra_phi_U_non(irho_min:irho_max)
  complex(kind(0d0))::nab_para_j_non(irho_min:irho_max)
  complex(kind(0d0))::nab_para_phi_non(irho_min:irho_max)
  complex(kind(0d0))::bra_phi_p_non(irho_min:irho_max)
  complex(kind(0d0)),allocatable::bra_phi_U_non_nl(:,:,:)
  complex(kind(0d0)),allocatable::nab_para_j_non_nl(:,:,:)
  complex(kind(0d0)),allocatable::nab_para_phi_non_nl(:,:,:)
  complex(kind(0d0)),allocatable::bra_phi_p_non_nl(:,:,:)
  complex(kind(0d0)),allocatable::bra_phi_U_non_n(:,:,:)
  complex(kind(0d0)),allocatable::nab_para_j_non_n(:,:,:)
  complex(kind(0d0)),allocatable::nab_para_phi_non_n(:,:,:)
  complex(kind(0d0)),allocatable::bra_phi_p_non_n(:,:,:)
  ! . time integrals
  integer::irk
  integer::jrk
  complex(kind(0d0))::U_rk(irho_min:irho_max,1:4)
  complex(kind(0d0))::psi_rk(irho_min:irho_max,1:4)
  complex(kind(0d0))::p_rk(irho_min:irho_max,1:4)
  complex(kind(0d0)),allocatable::U_rk_nl(:,:,:,:)
  complex(kind(0d0)),allocatable::psi_rk_nl(:,:,:,:)
  complex(kind(0d0)),allocatable::p_rk_nl(:,:,:,:)
  ! ..
  complex(kind(0d0))::del_rho_U(irho_min:irho_max)
  complex(kind(0d0))::del_alpha_U(irho_min:irho_max)
  complex(kind(0d0))::del_rho_phi(irho_min:irho_max)
  complex(kind(0d0))::del_alpha_phi(irho_min:irho_max)
  complex(kind(0d0))::del_rho_psi(irho_min:irho_max)
  complex(kind(0d0))::del_alpha_psi(irho_min:irho_max)
  complex(kind(0d0))::del_rho_j(irho_min:irho_max)
  complex(kind(0d0))::del_alpha_j(irho_min:irho_max)
  complex(kind(0d0))::del_rho_p(irho_min:irho_max)
  complex(kind(0d0))::del_alpha_p(irho_min:irho_max)
  complex(kind(0d0)),allocatable::del_rho_U_n(:,:,:)
  complex(kind(0d0)),allocatable::del_alpha_U_n(:,:,:)
  complex(kind(0d0)),allocatable::del_rho_phi_n(:,:,:)
  complex(kind(0d0)),allocatable::del_alpha_phi_n(:,:,:)
  complex(kind(0d0)),allocatable::del_rho_psi_n(:,:,:)
  complex(kind(0d0)),allocatable::del_alpha_psi_n(:,:,:)
  complex(kind(0d0)),allocatable::del_rho_j_n(:,:,:)
  complex(kind(0d0)),allocatable::del_alpha_j_n(:,:,:)
  complex(kind(0d0)),allocatable::del_rho_p_n(:,:,:)
  complex(kind(0d0)),allocatable::del_alpha_p_n(:,:,:)
  complex(kind(0d0)),allocatable::del_rho_U_n_trans(:,:,:)
  complex(kind(0d0)),allocatable::del_alpha_U_n_trans(:,:,:)
  complex(kind(0d0)),allocatable::del_rho_phi_n_trans(:,:,:)
  complex(kind(0d0)),allocatable::del_alpha_phi_n_trans(:,:,:)
  complex(kind(0d0)),allocatable::del_rho_psi_n_trans(:,:,:)
  complex(kind(0d0)),allocatable::del_alpha_psi_n_trans(:,:,:)
  complex(kind(0d0)),allocatable::del_rho_j_n_trans(:,:,:)
  complex(kind(0d0)),allocatable::del_alpha_j_n_trans(:,:,:)
  complex(kind(0d0)),allocatable::del_rho_p_n_trans(:,:,:)
  complex(kind(0d0)),allocatable::del_alpha_p_n_trans(:,:,:)
  double precision::bra_phi_U_non_real(0:numdata_alpha-1)
  double precision::nab_para_j_non_real(0:numdata_alpha-1)
  double precision::nab_para_phi_non_real(0:numdata_alpha-1)
  double precision::bra_phi_p_non_real(0:numdata_alpha-1)
  complex(kind(0d0)),allocatable::bra_phi_U_non_n_trans(:,:,:)
  complex(kind(0d0)),allocatable::nab_para_j_non_n_trans(:,:,:)
  complex(kind(0d0)),allocatable::nab_para_phi_non_n_trans(:,:,:)
  complex(kind(0d0)),allocatable::bra_phi_p_non_n_trans(:,:,:)
  integer(kind=8)::plan_c2r
  integer(kind=8)::plan_r2c
  complex(kind(0d0))::f1(0:numdata_alpha/2)
  complex(kind(0d0))::f2(0:numdata_alpha/2)
  complex(kind(0d0))::f3(0:numdata_alpha/2)
  complex(kind(0d0))::f4(0:numdata_alpha/2)
  complex(kind(0d0))::f5(0:numdata_alpha/2)
  complex(kind(0d0))::f6(0:numdata_alpha/2)
  complex(kind(0d0))::f7(0:numdata_alpha/2)
  complex(kind(0d0))::f8(0:numdata_alpha/2)
  complex(kind(0d0))::f9(0:numdata_alpha/2)
  complex(kind(0d0))::f10(0:numdata_alpha/2)
  double precision::f1_idft(0:numdata_alpha-1)
  double precision::f2_idft(0:numdata_alpha-1)
  double precision::f3_idft(0:numdata_alpha-1)
  double precision::f4_idft(0:numdata_alpha-1)
  double precision::f5_idft(0:numdata_alpha-1)
  double precision::f6_idft(0:numdata_alpha-1)
  double precision::f7_idft(0:numdata_alpha-1)
  double precision::f8_idft(0:numdata_alpha-1)
  double precision::f9_idft(0:numdata_alpha-1)
  double precision::f10_idft(0:numdata_alpha-1)
  ! ..
  integer::ilu
  complex(kind(0d0))::coe_lu(1:3,irho_min+1:irho_max-1)
  complex(kind(0d0)),allocatable::coe_lu_n(:,:,:,:)
  complex(kind(0d0))::rhs_lu(irho_min+1:irho_max-1)
  ! ..
  integer(kind=8)::plan_back
  integer(kind=8)::plan_for
  complex(kind(0d0))::f(0:numdata_z-1)
  complex(kind(0d0))::f_dft(-numdata_z/2:numdata_z/2-1)
  ! . linear analysis
  integer::irho_peak
  double precision,allocatable::growth(:)
  double precision,allocatable::freq(:)
  double precision::growth_gather(0:n_max)
  double precision::freq_gather(0:n_max)
  complex(kind(0d0)),allocatable::p_old_n(:,:,:)
  double precision::p_amp(irho_min:irho_max)
  double precision::p_irho_peak
  ! . average
  double precision::U_avg(irho_min:irho_max)
  double precision::phi_avg(irho_min:irho_max)
  double precision::j_avg(irho_min:irho_max)
  double precision::p_avg(irho_min:irho_max)
  double precision::del_rho_U_avg(irho_min:irho_max)
  double precision::del_rho_phi_avg(irho_min:irho_max)
  double precision::del_rho_j_avg(irho_min:irho_max)
  double precision::del_rho_p_avg(irho_min:irho_max)
  ! . contours
  complex(kind(0d0)),allocatable::ishifter_cylin(:,:,:)
  complex(kind(0d0)),allocatable::U_n_cylin(:,:,:)
  complex(kind(0d0)),allocatable::phi_n_cylin(:,:,:)
  complex(kind(0d0)),allocatable::psi_n_cylin(:,:,:)
  complex(kind(0d0)),allocatable::j_n_cylin(:,:,:)
  complex(kind(0d0)),allocatable::p_n_cylin(:,:,:)
  double precision::U_real_cylin(irho_min:irho_max,0:numdata_theta)
  double precision::phi_real_cylin(irho_min:irho_max,0:numdata_theta)
  double precision::psi_real_cylin(irho_min:irho_max,0:numdata_theta)
  double precision::j_real_cylin(irho_min:irho_max,0:numdata_theta)
  double precision::p_real_cylin(irho_min:irho_max,0:numdata_theta)
  integer(kind=8)::plan_back_cylin
  complex(kind(0d0))::f_cylin(0:numdata_theta-1)
  ! . energy balances
  ! . kinetic energy
  complex(kind(0d0))::phi_cc(irho_min:irho_max)
  complex(kind(0d0))::phi_old_cc(irho_min:irho_max)
  double precision::dE_phi_nl_den(irho_min:irho_max)
  double precision::dE_phi_nl
  double precision,allocatable::dE_phi_n(:)
  double precision::dE_phi
  ! . linear energy sources
  double precision::S_phi_dphi_eq_nl_den(irho_min:irho_max)
  double precision::S_phi_q_nl_den(irho_min:irho_max)
  double precision::S_phi_dj_eq_nl_den(irho_min:irho_max)
  double precision::S_phi_x_p_nl_den(irho_min:irho_max)
  double precision::S_phi_mu_nl_den(irho_min:irho_max)
  double precision::S_phi_dphi_eq_nl
  double precision::S_phi_q_nl
  double precision::S_phi_dj_eq_nl
  double precision::S_phi_x_p_nl
  double precision::S_phi_mu_nl
  double precision,allocatable::S_phi_dphi_eq_n(:)
  double precision,allocatable::S_phi_q_n(:)
  double precision,allocatable::S_phi_dj_eq_n(:)
  double precision,allocatable::S_phi_x_p_n(:)
  double precision,allocatable::S_phi_mu_n(:)
  double precision::S_phi_dphi_eq
  double precision::S_phi_q
  double precision::S_phi_dj_eq
  double precision::S_phi_x_p
  double precision::S_phi_mu
  double precision::L_phi
  ! . nonlinear energy sources
  double precision::S_phi_phi_U_nl_den(irho_min:irho_max)
  double precision::S_phi_psi_j_nl_den(irho_min:irho_max)
  double precision::S_phi_phi_U_nl
  double precision::S_phi_psi_j_nl
  double precision,allocatable::S_phi_phi_U_n(:)
  double precision,allocatable::S_phi_psi_j_n(:)
  double precision::S_phi_phi_U
  double precision::S_phi_psi_j
  double precision::N_phi
  ! . magnetic energy
  complex(kind(0d0))::j_cc(irho_min:irho_max)
  complex(kind(0d0))::j_old_cc(irho_min:irho_max)
  double precision::dE_psi_nl_den(irho_min:irho_max)
  double precision::dE_psi_nl
  double precision,allocatable::dE_psi_n(:)
  double precision::dE_psi
  ! . linear energy sources
  double precision::S_j_q_nl_den(irho_min:irho_max)
  double precision::S_j_dphi_eq_nl_den(irho_min:irho_max)
  double precision::S_j_eta_nl_den(irho_min:irho_max)
  double precision::S_j_q_nl
  double precision::S_j_dphi_eq_nl
  double precision::S_j_eta_nl
  double precision,allocatable::S_j_q_n(:)
  double precision,allocatable::S_j_dphi_eq_n(:)
  double precision,allocatable::S_j_eta_n(:)
  double precision::S_j_q
  double precision::S_j_dphi_eq
  double precision::S_j_eta
  double precision::L_psi
  ! . nonlinear energy sources
  double precision::S_j_psi_phi_nl_den(irho_min:irho_max)
  double precision::S_j_psi_phi_nl
  double precision,allocatable::S_j_psi_phi_n(:)
  double precision::S_j_psi_phi
  double precision::N_psi
  ! . pressure energy
  double precision::dE_p_nl_den(irho_min:irho_max)
  double precision::dE_p_nl
  double precision,allocatable::dE_p_n(:)
  double precision::dE_p
  ! . linear energy sources
  complex(kind(0d0))::p_cc(irho_min:irho_max)
  double precision::S_p_dp_eq_nl_den(irho_min:irho_max)
  double precision::S_p_x_phi_nl_den(irho_min:irho_max)
  double precision::S_p_chi_nl_den(irho_min:irho_max)
  double precision::S_p_dp_eq_nl
  double precision::S_p_x_phi_nl
  double precision::S_p_chi_nl
  double precision,allocatable::S_p_dp_eq_n(:)
  double precision,allocatable::S_p_x_phi_n(:)
  double precision,allocatable::S_p_chi_n(:)
  double precision::S_p_dp_eq
  double precision::S_p_x_phi
  double precision::S_p_chi
  double precision::L_p
  ! . nonlinear energy sources
  double precision::S_p_phi_p_nl_den(irho_min:irho_max)
  double precision::S_p_phi_p_nl
  double precision,allocatable::S_p_phi_p_n(:)
  double precision::S_p_phi_p
  double precision::N_p
  ! . total energy
  double precision::dE_tot
  double precision::L_tot
  double precision::N_tot
  ! . energies
  ! . kinetic energy
  double precision::E_phi_nl_den(irho_min:irho_max)
  double precision::E_phi_nl
  double precision,allocatable::E_phi_n(:)
  double precision::E_phi
  double precision::E_phi_n_gather(0:n_max)
  ! . magnetic energy
  double precision::E_psi_nl_den(irho_min:irho_max)
  double precision::E_psi_nl
  double precision,allocatable::E_psi_n(:)
  double precision::E_psi
  double precision::E_psi_n_gather(0:n_max)
  ! . pressure energy
  double precision::E_p_nl_den(irho_min:irho_max)
  double precision::E_p_nl
  double precision,allocatable::E_p_n(:)
  double precision::E_p
  double precision::E_p_n_gather(0:n_max)
  ! . energy spectra
  ! . kinetic energy
  double precision::E_phi_old_nl_den(irho_min:irho_max)
  double precision::E_phi_old_nl
  double precision,allocatable::E_phi_old_n(:)
  double precision,allocatable::int_E_phi_n(:)
  double precision::int_E_phi_n_gather(0:n_max)
  ! . magnetic energy
  double precision::E_psi_old_nl_den(irho_min:irho_max)
  double precision::E_psi_old_nl
  double precision,allocatable::E_psi_old_n(:)
  double precision,allocatable::int_E_psi_n(:)
  double precision::int_E_psi_n_gather(0:n_max)
  ! . pressure energy
  double precision::E_p_old_nl_den(irho_min:irho_max)
  double precision::E_p_old_nl
  double precision,allocatable::E_p_old_n(:)
  double precision,allocatable::int_E_p_n(:)
  double precision::int_E_p_n_gather(0:n_max)
  ! . helicity balance
  complex(kind(0d0))::U_cc(irho_min:irho_max)
  complex(kind(0d0))::U_old_cc(irho_min:irho_max)
  double precision::dH_U_psi_nl_den(irho_min:irho_max)
  double precision::dH_U_psi_nl
  double precision,allocatable::dH_U_psi_n(:)
  double precision::dH_U_psi
  ! . linear helicity sources
  complex(kind(0d0))::psi_cc(irho_min:irho_max)
  double precision::S_psi_dU_eq_nl_den(irho_min:irho_max)
  double precision::S_U_psi_q_nl_den(irho_min:irho_max)
  double precision::S_psi_x_p_nl_den(irho_min:irho_max)
  double precision::S_psi_mu_nl_den(irho_min:irho_max)
  double precision::S_U_eta_nl_den(irho_min:irho_max)
  double precision::S_psi_dU_eq_nl
  double precision::S_U_psi_q_nl
  double precision::S_psi_x_p_nl
  double precision::S_psi_mu_nl
  double precision::S_U_eta_nl
  double precision,allocatable::S_psi_dU_eq_n(:)
  double precision,allocatable::S_U_psi_q_n(:)
  double precision,allocatable::S_psi_x_p_n(:)
  double precision,allocatable::S_psi_mu_n(:)
  double precision,allocatable::S_U_eta_n(:)
  double precision::S_psi_dU_eq
  double precision::S_U_psi_q
  double precision::S_psi_x_p
  double precision::S_psi_mu
  double precision::S_U_eta
  double precision::L_U_psi
  ! . nonlinear helicity sources
  double precision::S_psi_phi_U_nl_den(irho_min:irho_max)
  double precision::S_psi_psi_j_nl_den(irho_min:irho_max)
  double precision::S_U_psi_phi_nl_den(irho_min:irho_max)
  double precision::S_psi_phi_U_nl
  double precision::S_psi_psi_j_nl
  double precision::S_U_psi_phi_nl
  double precision,allocatable::S_psi_phi_U_n(:)
  double precision,allocatable::S_psi_psi_j_n(:)
  double precision,allocatable::S_U_psi_phi_n(:)
  double precision::S_psi_phi_U
  double precision::S_psi_psi_j
  double precision::S_U_psi_phi
  double precision::N_U_psi
  ! . helicity
  double precision::H_U_psi_nl_den(irho_min:irho_max)
  double precision::H_U_psi_nl
  double precision,allocatable::H_U_psi_n(:)
  double precision::H_U_psi
  double precision::H_U_psi_n_gather(0:n_max)
  ! . helicity spectrum
  double precision::H_U_psi_old_nl_den(irho_min:irho_max)
  double precision::H_U_psi_old_nl
  double precision,allocatable::H_U_psi_old_n(:)
  double precision,allocatable::int_H_U_psi_n(:)
  double precision::int_H_U_psi_n_gather(0:n_max)
end module field
! ------------------------------------------------------------------------------
module mpi_global
! ------------------------------------------------------------------------------
  use coordinate
  implicit none
  integer::myid
  integer::numprocs
  integer::ierr
  integer::iprocs
  integer::sn
  integer::en
  integer,allocatable::sn0(:)
  integer,allocatable::en0(:)
  integer::sirho
  integer::eirho
  integer,allocatable::sirho0(:)
  integer,allocatable::eirho0(:)
  ! . alltoallv
  integer::idata
  integer,allocatable::scounts_trans(:)
  integer,allocatable::rcounts_trans(:)
  integer,allocatable::sdispls_trans(:)
  integer,allocatable::rdispls_trans(:)
  complex(kind(0d0)),allocatable::a_send_trans(:)
  complex(kind(0d0)),allocatable::a_recv_trans(:)
  ! . allgatherv
  integer::scount_gather
  integer,allocatable::rcounts_gather(:)
  integer,allocatable::rdispls_gather(:)
  complex(kind(0d0)),allocatable::a_send_gather(:)
  complex(kind(0d0))::a_recv_gather(0:n_max)
  ! . allreduce
  double precision::a_send_op
  double precision::a_recv_op
end module mpi_global
! ------------------------------------------------------------------------------
program main
! ------------------------------------------------------------------------------
  use time
  use mpi_global
  implicit none

  call start_mpi
  call decompose_domains
  call allocate_variables
  call mpi_setups

  call coordinates
  call equilibria
  call time_steps
  call factorize
  call fftw_plans

  ! . time loop
  do it=it_start,it_end
    call present_time

    if(it==it_start)then
      if(irestart==0)then
        call initial_conditions
        goto 10
      else
        call load_restart
        goto 20
      end if
    end if

    call time_integrals

    if(it==it_end)then
      call save_restart
    end if

    10 continue

    if(mod(it,it_mskip)==0)then
      call energy_balances
      call energies
      call helicity_balance
      call helicity
    end if

    if(mod(it,it_skip)==0)then
      ! call linear_analysis
      call contours
      if(myid==0)then
        call average
      end if
    end if

    20 continue

    call energy_spectra
    call helicity_spectrum
  end do ! : it

  call end_mpi

end program main
! ------------------------------------------------------------------------------
subroutine start_mpi
! ------------------------------------------------------------------------------
  use mpi_global
  use mpi
  implicit none

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world,myid,ierr)
  call mpi_comm_size(mpi_comm_world,numprocs,ierr)

  if(myid==0)then
    write(6,*)'numprocs=',numprocs
  end if

  call mpi_barrier(mpi_comm_world,ierr)

end subroutine start_mpi
! ------------------------------------------------------------------------------
subroutine decompose_domains
! ------------------------------------------------------------------------------
  implicit none

  call decompose_domain_n
  call decompose_domain_rho

end subroutine decompose_domains
! ------------------------------------------------------------------------------
subroutine allocate_variables
! ------------------------------------------------------------------------------
  use coordinate
  use field
  use mpi_global
  implicit none

  ! . equilibria
  allocate(l_shift(irho_min:irho_max,sn:en))
  allocate(k_theta(irho_min:irho_max,sn:en))
  allocate(k_para(irho_min:irho_max,sn:en,l_min:l_max))
  allocate(fact_bra_x(irho_min:irho_max,sn:en,0:numdata_z-1))
  allocate(fact_del_rho_lap_perp(irho_min:irho_max,sn:en,0:numdata_z-1))
  allocate(fact_lap_perp(irho_min:irho_max,sn:en,0:numdata_z-1))
  allocate(ishifter(irho_min:irho_max,sn:en,0:numdata_z-1))
  allocate(shifter(irho_min:irho_max,sn:en,0:numdata_z-1))
  ! . spectra
  allocate(U_nl(irho_min:irho_max,sn:en,l_min:l_max))
  allocate(phi_nl(irho_min:irho_max,sn:en,l_min:l_max))
  allocate(psi_nl(irho_min:irho_max,sn:en,l_min:l_max))
  allocate(j_nl(irho_min:irho_max,sn:en,l_min:l_max))
  allocate(p_nl(irho_min:irho_max,sn:en,l_min:l_max))
  allocate(U_n(irho_min:irho_max,sn:en,0:numdata_z-1))
  allocate(phi_n(irho_min:irho_max,sn:en,0:numdata_z-1))
  allocate(psi_n(irho_min:irho_max,sn:en,0:numdata_z-1))
  allocate(j_n(irho_min:irho_max,sn:en,0:numdata_z-1))
  allocate(p_n(irho_min:irho_max,sn:en,0:numdata_z-1))
  ! ..
  allocate(U_old_nl(irho_min:irho_max,sn:en,l_min:l_max))
  allocate(phi_old_nl(irho_min:irho_max,sn:en,l_min:l_max))
  allocate(psi_old_nl(irho_min:irho_max,sn:en,l_min:l_max))
  allocate(j_old_nl(irho_min:irho_max,sn:en,l_min:l_max))
  allocate(p_old_nl(irho_min:irho_max,sn:en,l_min:l_max))
  ! . linear terms
  allocate(bra_phi_U_lin_nl(irho_min:irho_max,sn:en,l_min:l_max))
  allocate(nab_para_j_lin_nl(irho_min:irho_max,sn:en,l_min:l_max))
  allocate(bra_x_p_nl(irho_min:irho_max,sn:en,l_min:l_max))
  allocate(lap_perp_U_nl(irho_min:irho_max,sn:en,l_min:l_max))
  allocate(nab_para_phi_lin_nl(irho_min:irho_max,sn:en,l_min:l_max))
  allocate(bra_phi_p_lin_nl(irho_min:irho_max,sn:en,l_min:l_max))
  allocate(bra_x_phi_nl(irho_min:irho_max,sn:en,l_min:l_max))
  allocate(lap_perp_p_nl(irho_min:irho_max,sn:en,l_min:l_max))
  allocate(bra_x_p_n(irho_min:irho_max,sn:en,0:numdata_z-1))
  allocate(lap_perp_U_n(irho_min:irho_max,sn:en,0:numdata_z-1))
  allocate(bra_x_phi_n(irho_min:irho_max,sn:en,0:numdata_z-1))
  allocate(lap_perp_p_n(irho_min:irho_max,sn:en,0:numdata_z-1))
  ! . nonlinear terms
  allocate(bra_phi_U_non_nl(irho_min:irho_max,sn:en,l_min:l_max))
  allocate(nab_para_j_non_nl(irho_min:irho_max,sn:en,l_min:l_max))
  allocate(nab_para_phi_non_nl(irho_min:irho_max,sn:en,l_min:l_max))
  allocate(bra_phi_p_non_nl(irho_min:irho_max,sn:en,l_min:l_max))
  allocate(bra_phi_U_non_n(irho_min:irho_max,sn:en,0:numdata_z-1))
  allocate(nab_para_j_non_n(irho_min:irho_max,sn:en,0:numdata_z-1))
  allocate(nab_para_phi_non_n(irho_min:irho_max,sn:en,0:numdata_z-1))
  allocate(bra_phi_p_non_n(irho_min:irho_max,sn:en,0:numdata_z-1))
  ! . time integrals
  allocate(U_rk_nl(irho_min:irho_max,sn:en,l_min:l_max,1:4))
  allocate(psi_rk_nl(irho_min:irho_max,sn:en,l_min:l_max,1:4))
  allocate(p_rk_nl(irho_min:irho_max,sn:en,l_min:l_max,1:4))
  ! ..
  allocate(del_rho_U_n(irho_min:irho_max,sn:en,0:numdata_z-1))
  allocate(del_alpha_U_n(irho_min:irho_max,sn:en,0:numdata_z-1))
  allocate(del_rho_phi_n(irho_min:irho_max,sn:en,0:numdata_z-1))
  allocate(del_alpha_phi_n(irho_min:irho_max,sn:en,0:numdata_z-1))
  allocate(del_rho_psi_n(irho_min:irho_max,sn:en,0:numdata_z-1))
  allocate(del_alpha_psi_n(irho_min:irho_max,sn:en,0:numdata_z-1))
  allocate(del_rho_j_n(irho_min:irho_max,sn:en,0:numdata_z-1))
  allocate(del_alpha_j_n(irho_min:irho_max,sn:en,0:numdata_z-1))
  allocate(del_rho_p_n(irho_min:irho_max,sn:en,0:numdata_z-1))
  allocate(del_alpha_p_n(irho_min:irho_max,sn:en,0:numdata_z-1))
  allocate(del_rho_U_n_trans(0:n_max,sirho:eirho,0:numdata_z-1))
  allocate(del_alpha_U_n_trans(0:n_max,sirho:eirho,0:numdata_z-1))
  allocate(del_rho_phi_n_trans(0:n_max,sirho:eirho,0:numdata_z-1))
  allocate(del_alpha_phi_n_trans(0:n_max,sirho:eirho,0:numdata_z-1))
  allocate(del_rho_psi_n_trans(0:n_max,sirho:eirho,0:numdata_z-1))
  allocate(del_alpha_psi_n_trans(0:n_max,sirho:eirho,0:numdata_z-1))
  allocate(del_rho_j_n_trans(0:n_max,sirho:eirho,0:numdata_z-1))
  allocate(del_alpha_j_n_trans(0:n_max,sirho:eirho,0:numdata_z-1))
  allocate(del_rho_p_n_trans(0:n_max,sirho:eirho,0:numdata_z-1))
  allocate(del_alpha_p_n_trans(0:n_max,sirho:eirho,0:numdata_z-1))
  allocate(bra_phi_U_non_n_trans(0:n_max,sirho:eirho,0:numdata_z-1))
  allocate(nab_para_j_non_n_trans(0:n_max,sirho:eirho,0:numdata_z-1))
  allocate(nab_para_phi_non_n_trans(0:n_max,sirho:eirho,0:numdata_z-1))
  allocate(bra_phi_p_non_n_trans(0:n_max,sirho:eirho,0:numdata_z-1))
  ! ..
  allocate(coe_lu_n(1:3,irho_min+1:irho_max-1,sn:en,0:numdata_z-1))
  ! . linear analysis
  allocate(growth(sn:en))
  allocate(freq(sn:en))
  allocate(p_old_n(irho_min:irho_max,sn:en,0:numdata_z-1))
  ! . contours
  allocate(ishifter_cylin(irho_min:irho_max,sn:en,0:numdata_theta-1))
  allocate(U_n_cylin(irho_min:irho_max,sn:en,0:numdata_theta-1))
  allocate(phi_n_cylin(irho_min:irho_max,sn:en,0:numdata_theta-1))
  allocate(psi_n_cylin(irho_min:irho_max,sn:en,0:numdata_theta-1))
  allocate(j_n_cylin(irho_min:irho_max,sn:en,0:numdata_theta-1))
  allocate(p_n_cylin(irho_min:irho_max,sn:en,0:numdata_theta-1))
  ! . energy balances
  ! . kinetic energy
  allocate(dE_phi_n(sn:en))
  allocate(S_phi_dphi_eq_n(sn:en))
  allocate(S_phi_q_n(sn:en))
  allocate(S_phi_dj_eq_n(sn:en))
  allocate(S_phi_x_p_n(sn:en))
  allocate(S_phi_mu_n(sn:en))
  allocate(S_phi_phi_U_n(sn:en))
  allocate(S_phi_psi_j_n(sn:en))
  ! . magnetic energy
  allocate(dE_psi_n(sn:en))
  allocate(S_j_q_n(sn:en))
  allocate(S_j_dphi_eq_n(sn:en))
  allocate(S_j_eta_n(sn:en))
  allocate(S_j_psi_phi_n(sn:en))
  ! . pressure energy
  allocate(dE_p_n(sn:en))
  allocate(S_p_dp_eq_n(sn:en))
  allocate(S_p_x_phi_n(sn:en))
  allocate(S_p_chi_n(sn:en))
  allocate(S_p_phi_p_n(sn:en))
  ! . energies
  allocate(E_phi_n(sn:en))
  allocate(E_psi_n(sn:en))
  allocate(E_p_n(sn:en))
  ! . energy spectra
  allocate(E_phi_old_n(sn:en))
  allocate(int_E_phi_n(sn:en))
  allocate(E_psi_old_n(sn:en))
  allocate(int_E_psi_n(sn:en))
  allocate(E_p_old_n(sn:en))
  allocate(int_E_p_n(sn:en))
  ! . helicity balance
  allocate(dH_U_psi_n(sn:en))
  allocate(S_psi_dU_eq_n(sn:en))
  allocate(S_U_psi_q_n(sn:en))
  allocate(S_psi_x_p_n(sn:en))
  allocate(S_psi_mu_n(sn:en))
  allocate(S_U_eta_n(sn:en))
  allocate(S_psi_phi_U_n(sn:en))
  allocate(S_psi_psi_j_n(sn:en))
  allocate(S_U_psi_phi_n(sn:en))
  ! . helicity
  allocate(H_U_psi_n(sn:en))
  ! . helicity spectrum
  allocate(H_U_psi_old_n(sn:en))
  allocate(int_H_U_psi_n(sn:en))

end subroutine allocate_variables
! ------------------------------------------------------------------------------
subroutine mpi_setups
! ------------------------------------------------------------------------------
  implicit none

  call alltoallv_setups
  call allgatherv_setups

end subroutine mpi_setups
! ------------------------------------------------------------------------------
subroutine coordinates
! ------------------------------------------------------------------------------
  use coordinate
  implicit none

  do irho=irho_min,irho_max
    rho(irho)=rho_min+dble(irho-irho_min)*drho
  end do

  do iz=0,numdata_z-1
    z(iz)=-pi+2.0d0*pi*dble(iz)*inumdata_z
  end do

  do itheta=0,numdata_theta
    theta(itheta)=-pi+2.0d0*pi*dble(itheta)*inumdata_theta
  end do

end subroutine coordinates
! ------------------------------------------------------------------------------
subroutine equilibria
! ------------------------------------------------------------------------------
  use time
  use coordinate
  use parameter
  use field
  use mpi_global
  implicit none
  double precision,parameter::B_zeta_ir_min=1.0d0
  double precision::B_theta(irho_min:irho_max)
  double precision::dB_theta(irho_min:irho_max)
  double precision::d2B_theta(irho_min:irho_max)
  double precision::fact(irho_min:irho_max)
  double precision::B_theta_sq(irho_min:irho_max)
  integer::ipc
  double precision::B_theta_sq_pc(0:1)
  ! . linear analysis
  double precision::dp_eq_ir_peak

  ! . cylindrical MHD equilibrium
  do irho=irho_min,irho_max
    phi_eq(irho)=0.0d0
    q(irho)=2.0d0+2.0d0*rho(irho)**2
    p_eq(irho)=beta/(2.0d0*eps)*(1.0d0-rho(irho)**2/2.0d0)&
    *(1.0d0-tanh((rho(irho)-0.8d0)/0.05d0))
  end do

  call lap_rho_real(phi_eq,U_eq)
  call del_rho_real(U_eq,dU_eq)
  call del_rho_real(phi_eq,dphi_eq)
  call del_rho_real(q,dq)
  call del_rho_real(p_eq,dp_eq)

  do irho=irho_min,irho_max
    fact(irho)=(rho(irho)*q(irho)*dq(irho)-q(irho)**2&
    +(eps*rho(irho))**2)/rho(irho)**3
  end do

  B_theta_sq(irho_min)=(rho_min*B_zeta_ir_min/q(irho_min))**2

  do irho=irho_min+1,irho_max
    do ipc=0,1
      B_theta_sq_pc(ipc)=-2.0d0/((q(irho-1+ipc)/rho(irho-1+ipc))**2+eps**2)&
      *(eps*dp_eq(irho-1+ipc)/2.0d0+fact(irho-1+ipc)*B_theta_sq(irho-1+ipc))

      if(ipc==1)then
        B_theta_sq(irho)=B_theta_sq(irho-1)&
        +drho/2.0d0*(B_theta_sq_pc(0)+B_theta_sq_pc(1))
      else
        B_theta_sq(irho)=B_theta_sq(irho-1)+drho*B_theta_sq_pc(ipc)
      end if
    end do ! : ipc
  end do ! : irho

  do irho=irho_min,irho_max
    B_theta(irho)=sqrt(B_theta_sq(irho))
  end do

  call del_rho_real(B_theta,dB_theta)
  call del2_rho_real(B_theta,d2B_theta)

  do irho=irho_min,irho_max
    j_eq(irho)=(rho(irho)*dB_theta(irho)+B_theta(irho))/rho(irho)
    dj_eq(irho)=(rho(irho)**2*d2B_theta(irho)+rho(irho)*dB_theta(irho)&
    -B_theta(irho))/rho(irho)**2
  end do

  if(irestart==0)then
    if(myid==0)then
      open(unit=10,file='../data/equil_profile.dat',form='formatted')
      write(10,fmt='(6a16)')'rho','U_eq','phi_eq','q','j_eq','p_eq'
      do irho=irho_min,irho_max
        write(10,fmt='(6e16.4)')rho(irho),&
        U_eq(irho),phi_eq(irho),q(irho),j_eq(irho),p_eq(irho)
      end do
      close(10)

      open(unit=10,file='../data/equil_grad_profile.dat',form='formatted')
      write(10,fmt='(6a16)')'rho','dU_eq','dphi_eq','dq','dj_eq','dp_eq'
      do irho=irho_min,irho_max
        write(10,fmt='(6e16.4)')rho(irho),&
        dU_eq(irho),dphi_eq(irho),dq(irho),dj_eq(irho),dp_eq(irho)
      end do
      close(10)
    end if ! : myid==0
  end if ! : irestart==0

  call del2_rho_real(q,d2q)

  do n=sn,en
    do irho=irho_min,irho_max
      l_shift(irho,n)=dble(n)*q(irho)-floor(dble(n)*q(irho))
      k_theta(irho,n)=dble(n)*q(irho)/rho(irho)
    end do
  end do

  do l=l_min,l_max
    do n=sn,en
      do irho=irho_min,irho_max
        k_para(irho,n,l)=(dble(l)-l_shift(irho,n))/q(irho)
      end do
    end do
  end do

  do iz=0,numdata_z-1
    fact_del_rho_bra_x(iz)=sin(z(iz))

    do n=sn,en
      do irho=irho_min,irho_max
        fact_bra_x(irho,n,iz)=(0.0d0,1.0d0)*dble(n)&
        *(dq(irho)*z(iz)*sin(z(iz))+q(irho)/rho(irho)*cos(z(iz)))
        fact_del_rho_lap_perp(irho,n,iz)=1.0d0/rho(irho)&
        +2.0d0*(0.0d0,1.0d0)*dble(n)*dq(irho)*z(iz)
        fact_lap_perp(irho,n,iz)=-dble(n)&
        *(dble(n)*((dq(irho)*z(iz))**2+(q(irho)/rho(irho))**2)&
        -(0.0d0,1.0d0)*z(iz)*(d2q(irho)+dq(irho)/rho(irho)))
        ishifter(irho,n,iz)=exp(-(0.0d0,1.0d0)*l_shift(irho,n)*z(iz))
        shifter(irho,n,iz)=exp((0.0d0,1.0d0)*l_shift(irho,n)*z(iz))
      end do
    end do
  end do ! : iz

  irho_peak=irho_min
  dp_eq_ir_peak=0.0d0
  do irho=irho_min,irho_max
    if(dp_eq_ir_peak<abs(dp_eq(irho)))then
      irho_peak=irho
      dp_eq_ir_peak=abs(dp_eq(irho))
    end if
  end do

  do itheta=0,numdata_theta-1
    do n=sn,en
      do irho=irho_min,irho_max
        ishifter_cylin(irho,n,itheta)=exp((0.0d0,1.0d0)&
        *floor(dble(n)*q(irho))*theta(itheta))
      end do
    end do
  end do

end subroutine equilibria
! ------------------------------------------------------------------------------
subroutine time_steps
! ------------------------------------------------------------------------------
  use time
  implicit none

  dt_rk(1)=dt/2.0d0
  dt_rk(2)=dt/2.0d0
  dt_rk(3)=dt
  dt_rk(4)=dt/6.0d0

end subroutine time_steps
! ------------------------------------------------------------------------------
subroutine factorize
! ------------------------------------------------------------------------------
  use coordinate
  use field
  use mpi_global
  implicit none

  do iz=0,numdata_z-1
    do n=sn,en
      do irho=irho_min+1,irho_max-1
        coe_lu(1,irho)=1.0d0-drho/2.0d0*fact_del_rho_lap_perp(irho,n,iz)
        coe_lu(2,irho)=-2.0d0+drho_sq*fact_lap_perp(irho,n,iz)
        coe_lu(3,irho)=1.0d0+drho/2.0d0*fact_del_rho_lap_perp(irho,n,iz)
      end do

      do irho=irho_min+2,irho_max-1
        coe_lu(1,irho)=coe_lu(1,irho)/coe_lu(2,irho-1)
        coe_lu(2,irho)=coe_lu(2,irho)-coe_lu(1,irho)*coe_lu(3,irho-1)
      end do

      do irho=irho_min+1,irho_max-1
        do ilu=1,3
          coe_lu_n(ilu,irho,n,iz)=coe_lu(ilu,irho)
        end do
      end do
    end do ! : n
  end do ! : iz

end subroutine factorize
! ------------------------------------------------------------------------------
subroutine fftw_plans
! ------------------------------------------------------------------------------
  use coordinate
  use field
  implicit none
  include 'fftw3.f'

  call dfftw_plan_dft_1d(plan_back,numdata_z,f,f,fftw_backward,fftw_estimate)
  call dfftw_plan_dft_1d(plan_for,numdata_z,f,f,fftw_forward,fftw_estimate)
  call dfftw_plan_dft_c2r_1d(plan_c2r,numdata_alpha,f1,f1_idft,fftw_estimate)
  call dfftw_plan_dft_r2c_1d(plan_r2c,numdata_alpha,f1_idft,f1,fftw_estimate)
  call dfftw_plan_dft_1d(plan_back_cylin,numdata_theta,f_cylin,f_cylin,&
  fftw_backward,fftw_estimate)

end subroutine fftw_plans
! ------------------------------------------------------------------------------
subroutine present_time
! ------------------------------------------------------------------------------
  use time
  use mpi_global
  implicit none

  t=dt*dble(it)

  if(myid==0)then
    open(unit=10,file='../data/present_time.dat',form='formatted')
    write(10,*)'t=',t
    close(10)
  end if

end subroutine present_time
! ------------------------------------------------------------------------------
subroutine initial_conditions
! ------------------------------------------------------------------------------
  use coordinate
  use field
  use mpi_global
  implicit none
  integer::iseed
  integer::numseed
  integer,allocatable::seed(:)
  integer::irand
  double precision::U_rand1
  double precision::U_rand2
  double precision::j_rand1
  double precision::j_rand2
  double precision::p_rand1
  double precision::p_rand2

  call random_seed(size=numseed)
  allocate(seed(0:numseed-1))
  do iseed=0,numseed-1
    ! call system_clock(count=seed(iseed))
    seed(iseed)=numseed*(1+myid)+iseed
  end do
  call random_seed(put=seed)

  do l=l_min,l_max
    do n=sn,en
      do irho=irho_min,irho_max
        U(irho)=0.0d0
        j(irho)=0.0d0
        p(irho)=0.0d0
      end do

      if(n/=0)then
        do irand=1,10
          call random_number(U_rand1)
          call random_number(U_rand2)
          call random_number(j_rand1)
          call random_number(j_rand2)
          call random_number(p_rand1)
          call random_number(p_rand2)

          do irho=irho_min,irho_max
            U(irho)=U(irho)+1.0d-14/(1.0d0+dble(n)**2)&
            *sqrt(-2.0d0*log(U_rand1))*exp((0.0d0,1.0d0)*2.0d0*pi*U_rand2)&
            *sin(dble(irand)*pi*(rho(irho)-rho_min)/(rho_max-rho_min))
            j(irho)=j(irho)+1.0d-14/(1.0d0+dble(n)**2)&
            *sqrt(-2.0d0*log(j_rand1))*exp((0.0d0,1.0d0)*2.0d0*pi*j_rand2)&
            *sin(dble(irand)*pi*(rho(irho)-rho_min)/(rho_max-rho_min))
            p(irho)=p(irho)+1.0d-16/(1.0d0+dble(n)**2)&
            *sqrt(-2.0d0*log(p_rand1))*exp((0.0d0,1.0d0)*2.0d0*pi*p_rand2)&
            *sin(dble(irand)*pi*(rho(irho)-rho_min)/(rho_max-rho_min))
          end do
        end do ! : irand
      end if ! : n/=0
      ! else
      !   do irand=1,10
      !     call random_number(U_rand1)
      !     call random_number(U_rand2)
      !     call random_number(j_rand1)
      !     call random_number(j_rand2)
      !     call random_number(p_rand1)
      !     call random_number(p_rand2)

      !     do irho=irho_min,irho_max
      !       U(irho)=U(irho)+1.0d-10/(1.0d0+dble(n)**2)&
      !       *sqrt(-2.0d0*log(U_rand1))*cos(2.0d0*pi*U_rand2)&
      !       *sin(dble(irand)*pi*(rho(irho)-rho_min)/(rho_max-rho_min))
      !       j(irho)=j(irho)+1.0d-10/(1.0d0+dble(n)**2)&
      !       *sqrt(-2.0d0*log(j_rand1))*cos(2.0d0*pi*j_rand2)&
      !       *sin(dble(irand)*pi*(rho(irho)-rho_min)/(rho_max-rho_min))
      !       p(irho)=p(irho)+1.0d-12/(1.0d0+dble(n)**2)&
      !       *sqrt(-2.0d0*log(p_rand1))*cos(2.0d0*pi*p_rand2)&
      !       *sin(dble(irand)*pi*(rho(irho)-rho_min)/(rho_max-rho_min))
      !     end do
      !   end do ! : irand
      ! end if ! : n/=0

      U(irho_min)=0.0d0
      U(irho_max)=0.0d0
      p(irho_min)=0.0d0
      p(irho_max)=0.0d0

      do irho=irho_min,irho_max
        U_nl(irho,n,l)=U(irho)
        j_nl(irho,n,l)=j(irho)
        p_nl(irho,n,l)=p(irho)
      end do
    end do ! : n
  end do ! : l

  if(myid==0)then
    call conjugate_0l(U_nl)
    call conjugate_0l(j_nl)
    call conjugate_0l(p_nl)
  end if

  call idft_z(U_nl,U_n)
  call idft_z(j_nl,j_n)
  call idft_z(p_nl,p_n)

  if(myid==0)then
    call conjugate_0(U_n)
    call conjugate_0(j_n)
    call conjugate_0(p_n)
  end if

  do iz=0,numdata_z-1
    do n=sn,en
      do irho=irho_min,irho_max
        U(irho)=U_n(irho,n,iz)
        j(irho)=j_n(irho,n,iz)
      end do

      ! . definition of U
      call poisson_solver(U,phi)
      ! . Ampere's law
      call poisson_solver(-j,psi)

      do irho=irho_min,irho_max
        phi_n(irho,n,iz)=phi(irho)
        psi_n(irho,n,iz)=psi(irho)
      end do
    end do ! : n
  end do ! : iz

  if(myid==0)then
    call conjugate_0(phi_n)
    call conjugate_0(psi_n)
  end if

  call dft_z(phi_n,phi_nl)
  call dft_z(psi_n,psi_nl)

  if(myid==0)then
    call conjugate_0l(phi_nl)
    call conjugate_0l(psi_nl)
  end if

end subroutine initial_conditions
! ------------------------------------------------------------------------------
subroutine load_restart
! ------------------------------------------------------------------------------
  use time
  use coordinate
  use field
  use mpi_global
  implicit none

  write(filename,'("../data/restart_profile_",i3.3,".dat")')myid
  open(unit=10*(1+myid),file=filename,form='unformatted')
  do l=l_min,l_max
    do n=sn,en
      do irho=irho_min,irho_max
        read(10*(1+myid))U_nl(irho,n,l),psi_nl(irho,n,l),p_nl(irho,n,l)
      end do
    end do
  end do
  close(10*(1+myid))

  call phi_and_j

end subroutine load_restart
! ------------------------------------------------------------------------------
subroutine time_integrals
! ------------------------------------------------------------------------------
  use time
  use coordinate
  use parameter
  use field
  use mpi_global
  implicit none

  do l=l_min,l_max
    do n=sn,en
      do irho=irho_min,irho_max
        U_old_nl(irho,n,l)=U_nl(irho,n,l)
        psi_old_nl(irho,n,l)=psi_nl(irho,n,l)
        p_old_nl(irho,n,l)=p_nl(irho,n,l)

        phi_old_nl(irho,n,l)=phi_nl(irho,n,l)
        j_old_nl(irho,n,l)=j_nl(irho,n,l)
      end do ! : irho
    end do ! : n
  end do ! : l

  do irk=1,4
    call linear_terms
    call nonlinear_terms

    do l=l_min,l_max
      do n=sn,en
        do irho=irho_min,irho_max
          U_old(irho)=U_old_nl(irho,n,l)
          psi_old(irho)=psi_old_nl(irho,n,l)
          p_old(irho)=p_old_nl(irho,n,l)
          ! . linear terms
          bra_phi_U_lin(irho)=bra_phi_U_lin_nl(irho,n,l)
          nab_para_j_lin(irho)=nab_para_j_lin_nl(irho,n,l)
          bra_x_p(irho)=bra_x_p_nl(irho,n,l)
          lap_perp_U(irho)=lap_perp_U_nl(irho,n,l)
          nab_para_phi_lin(irho)=nab_para_phi_lin_nl(irho,n,l)
          j(irho)=j_nl(irho,n,l)
          bra_phi_p_lin(irho)=bra_phi_p_lin_nl(irho,n,l)
          bra_x_phi(irho)=bra_x_phi_nl(irho,n,l)
          lap_perp_p(irho)=lap_perp_p_nl(irho,n,l)
          ! . nonlinear terms
          bra_phi_U_non(irho)=bra_phi_U_non_nl(irho,n,l)
          nab_para_j_non(irho)=nab_para_j_non_nl(irho,n,l)
          nab_para_phi_non(irho)=nab_para_phi_non_nl(irho,n,l)
          bra_phi_p_non(irho)=bra_phi_p_non_nl(irho,n,l)

          ! . vorticity equation
          U_rk(irho,irk)=-bra_phi_U_lin(irho)+nab_para_j_lin(irho)&
          -bra_x_p(irho)+mu*lap_perp_U(irho)&
          -bra_phi_U_non(irho)+nab_para_j_non(irho)
          ! . generalized Ohm's law
          psi_rk(irho,irk)=-nab_para_phi_lin(irho)-eta*j(irho)&
          -nab_para_phi_non(irho)
          ! . pressure evolution equation
          p_rk(irho,irk)=-bra_phi_p_lin(irho)+2.0d0*beta*bra_x_phi(irho)&
          +chi*lap_perp_p(irho)&
          -bra_phi_p_non(irho)

          U_rk_nl(irho,n,l,irk)=U_rk(irho,irk)
          psi_rk_nl(irho,n,l,irk)=psi_rk(irho,irk)
          p_rk_nl(irho,n,l,irk)=p_rk(irho,irk)
        end do ! : irho

        if(irk==4)then
          do jrk=1,4
            do irho=irho_min,irho_max
              U_rk(irho,jrk)=U_rk_nl(irho,n,l,jrk)
              psi_rk(irho,jrk)=psi_rk_nl(irho,n,l,jrk)
              p_rk(irho,jrk)=p_rk_nl(irho,n,l,jrk)
            end do
          end do

          do irho=irho_min,irho_max
            U(irho)=U_old(irho)&
            +dt_rk(irk)*(U_rk(irho,1)+2.0d0*U_rk(irho,2)&
            +2.0d0*U_rk(irho,3)+U_rk(irho,4))
            psi(irho)=psi_old(irho)&
            +dt_rk(irk)*(psi_rk(irho,1)+2.0d0*psi_rk(irho,2)&
            +2.0d0*psi_rk(irho,3)+psi_rk(irho,4))
            p(irho)=p_old(irho)&
            +dt_rk(irk)*(p_rk(irho,1)+2.0d0*p_rk(irho,2)&
            +2.0d0*p_rk(irho,3)+p_rk(irho,4))
          end do
        else
          do irho=irho_min,irho_max
            U(irho)=U_old(irho)+dt_rk(irk)*U_rk(irho,irk)
            psi(irho)=psi_old(irho)+dt_rk(irk)*psi_rk(irho,irk)
            p(irho)=p_old(irho)+dt_rk(irk)*p_rk(irho,irk)
          end do
        end if ! : irk==4

        call boundary_conditions

        do irho=irho_min,irho_max
          U_nl(irho,n,l)=U(irho)
          psi_nl(irho,n,l)=psi(irho)
          p_nl(irho,n,l)=p(irho)
        end do
      end do ! : n
    end do ! : l

    call phi_and_j
  end do ! : irk

end subroutine time_integrals
! ------------------------------------------------------------------------------
subroutine save_restart
! ------------------------------------------------------------------------------
  use time
  use coordinate
  use field
  use mpi_global
  implicit none

  write(filename,'("../data/restart_profile_",i3.3,".dat")')myid
  open(unit=10*(1+myid),file=filename,form='unformatted')
  do l=l_min,l_max
    do n=sn,en
      do irho=irho_min,irho_max
        write(10*(1+myid))U_nl(irho,n,l),psi_nl(irho,n,l),p_nl(irho,n,l)
      end do
    end do
  end do
  close(10*(1+myid))

end subroutine save_restart
! ------------------------------------------------------------------------------
subroutine end_mpi
! ------------------------------------------------------------------------------
  use mpi_global
  use mpi
  implicit none

  call mpi_finalize(ierr)

end subroutine end_mpi
! ------------------------------------------------------------------------------
subroutine decompose_domain_n
! ------------------------------------------------------------------------------
  use coordinate
  use mpi_global
  use mpi
  implicit none
  integer::numdata
  integer::numfloor
  integer::numredistribute

  numdata=n_max-0+1
  numfloor=floor(dble(numdata)/dble(numprocs))
  numredistribute=numdata-numprocs*numfloor

  allocate(sn0(0:numprocs-1))
  allocate(en0(0:numprocs-1))

  sn0(0)=0
  do iprocs=1,numprocs-1
    sn0(iprocs)=0+iprocs*numfloor
  end do

  do iprocs=1,numredistribute
    sn0(iprocs)=sn0(iprocs)+iprocs
  end do
  do iprocs=numredistribute+1,numprocs-1
    sn0(iprocs)=sn0(iprocs)+numredistribute
  end do

  do iprocs=0,numprocs-2
    en0(iprocs)=sn0(iprocs+1)-1
  end do
  en0(numprocs-1)=n_max

  sn=sn0(myid)
  en=en0(myid)

  write(6,*)'myid=',myid,'sn=',sn,'en=',en,'en-sn+1=',en-sn+1

  call mpi_barrier(mpi_comm_world,ierr)

end subroutine decompose_domain_n
! ------------------------------------------------------------------------------
subroutine decompose_domain_rho
! ------------------------------------------------------------------------------
  use coordinate
  use mpi_global
  use mpi
  implicit none
  integer::numdata
  integer::numfloor
  integer::numredistribute

  numdata=irho_max-irho_min+1
  numfloor=floor(dble(numdata)/dble(numprocs))
  numredistribute=numdata-numprocs*numfloor

  allocate(sirho0(0:numprocs-1))
  allocate(eirho0(0:numprocs-1))

  sirho0(0)=irho_min
  do iprocs=1,numprocs-1
    sirho0(iprocs)=irho_min+iprocs*numfloor
  end do

  do iprocs=1,numredistribute
    sirho0(iprocs)=sirho0(iprocs)+iprocs
  end do
  do iprocs=numredistribute+1,numprocs-1
    sirho0(iprocs)=sirho0(iprocs)+numredistribute
  end do

  do iprocs=0,numprocs-2
    eirho0(iprocs)=sirho0(iprocs+1)-1
  end do
  eirho0(numprocs-1)=irho_max

  sirho=sirho0(myid)
  eirho=eirho0(myid)

  write(6,*)'myid=',myid,'sir=',sirho,'eir=',eirho,&
  'eir-sir+1=',eirho-sirho+1

  call mpi_barrier(mpi_comm_world,ierr)

end subroutine decompose_domain_rho
! ------------------------------------------------------------------------------
subroutine alltoallv_setups
! ------------------------------------------------------------------------------
  use coordinate
  use mpi_global
  implicit none
  integer::numdata_send
  integer::numdata_recv

  allocate(scounts_trans(0:numprocs-1))
  allocate(rcounts_trans(0:numprocs-1))

  do iprocs=0,numprocs-1
    scounts_trans(iprocs)=(en-sn+1)*(eirho0(iprocs)-sirho0(iprocs)+1)
    rcounts_trans(iprocs)=(en0(iprocs)-sn0(iprocs)+1)*(eirho-sirho+1)
  end do

  allocate(sdispls_trans(0:numprocs-1))
  allocate(rdispls_trans(0:numprocs-1))

  sdispls_trans(0)=0
  rdispls_trans(0)=0
  do iprocs=1,numprocs-1
    sdispls_trans(iprocs)=sdispls_trans(iprocs-1)+scounts_trans(iprocs-1)
    rdispls_trans(iprocs)=rdispls_trans(iprocs-1)+rcounts_trans(iprocs-1)
  end do

  numdata_send=(en-sn+1)*(irho_max-irho_min+1)
  numdata_recv=(n_max-0+1)*(eirho-sirho+1)

  allocate(a_send_trans(0:numdata_send-1))
  allocate(a_recv_trans(0:numdata_recv-1))

end subroutine alltoallv_setups
! ------------------------------------------------------------------------------
subroutine allgatherv_setups
! ------------------------------------------------------------------------------
  use coordinate
  use mpi_global
  implicit none

  scount_gather=en-sn+1

  allocate(rcounts_gather(0:numprocs-1))

  do iprocs=0,numprocs-1
    rcounts_gather(iprocs)=en0(iprocs)-sn0(iprocs)+1
  end do

  allocate(rdispls_gather(0:numprocs-1))

  rdispls_gather(0)=0
  do iprocs=1,numprocs-1
    rdispls_gather(iprocs)=rdispls_gather(iprocs-1)&
    +rcounts_gather(iprocs-1)
  end do

  allocate(a_send_gather(sn:en))

end subroutine allgatherv_setups
! ------------------------------------------------------------------------------
subroutine linear_terms
! ------------------------------------------------------------------------------
  use coordinate
  use field
  use mpi_global
  implicit none

  do l=l_min,l_max
    do n=sn,en
      do irho=irho_min,irho_max
        U(irho)=U_nl(irho,n,l)
        phi(irho)=phi_nl(irho,n,l)
        psi(irho)=psi_nl(irho,n,l)
        j(irho)=j_nl(irho,n,l)
        p(irho)=p_nl(irho,n,l)

        bra_phi_U_lin(irho)=(0.0d0,1.0d0)*k_theta(irho,n)&
        *(dphi_eq(irho)*U(irho)-dU_eq(irho)*phi(irho))
        nab_para_j_lin(irho)=(0.0d0,1.0d0)*(k_para(irho,n,l)*j(irho)&
        +k_theta(irho,n)*dj_eq(irho)*psi(irho))
        nab_para_phi_lin(irho)=(0.0d0,1.0d0)*(k_para(irho,n,l)*phi(irho)&
        +k_theta(irho,n)*dphi_eq(irho)*psi(irho))
        bra_phi_p_lin(irho)=(0.0d0,1.0d0)*k_theta(irho,n)&
        *(dphi_eq(irho)*p(irho)-dp_eq(irho)*phi(irho))

        bra_phi_U_lin_nl(irho,n,l)=bra_phi_U_lin(irho)
        nab_para_j_lin_nl(irho,n,l)=nab_para_j_lin(irho)
        nab_para_phi_lin_nl(irho,n,l)=nab_para_phi_lin(irho)
        bra_phi_p_lin_nl(irho,n,l)=bra_phi_p_lin(irho)
      end do ! : irho
    end do ! : n
  end do ! : l

  do iz=0,numdata_z-1
    do n=sn,en
      do irho=irho_min,irho_max
        U(irho)=U_n(irho,n,iz)
        phi(irho)=phi_n(irho,n,iz)
        p(irho)=p_n(irho,n,iz)
      end do

      call bra_x_spec(p,bra_x_p)
      call lap_perp_spec(U,lap_perp_U)
      call bra_x_spec(phi,bra_x_phi)
      call lap_perp_spec(p,lap_perp_p)

      do irho=irho_min,irho_max
        bra_x_p_n(irho,n,iz)=bra_x_p(irho)
        lap_perp_U_n(irho,n,iz)=lap_perp_U(irho)
        bra_x_phi_n(irho,n,iz)=bra_x_phi(irho)
        lap_perp_p_n(irho,n,iz)=lap_perp_p(irho)
      end do
    end do ! : n
  end do ! : iz

  call dft_z(bra_x_p_n,bra_x_p_nl)
  call dft_z(lap_perp_U_n,lap_perp_U_nl)
  call dft_z(bra_x_phi_n,bra_x_phi_nl)
  call dft_z(lap_perp_p_n,lap_perp_p_nl)

  ! if(myid==0)then
  !   n=0
  !   do l=l_min+1,l_max-1
  !     do irho=irho_min,irho_max
  !       bra_x_p_nl(irho,n,l)=bra_x_p_nl(irho,n,l)&
  !       -(0.0d0,1.0d0)/(2.0d0*rho(irho))&
  !       *(dble(l+1)*p_nl(irho,n,l+1)+dble(l-1)*p_nl(irho,n,l-1))
  !       bra_x_phi_nl(irho,n,l)=bra_x_phi_nl(irho,n,l)&
  !       -(0.0d0,1.0d0)/(2.0d0*rho(irho))&
  !       *(dble(l+1)*phi_nl(irho,n,l+1)+dble(l-1)*phi_nl(irho,n,l-1))
  !     end do
  !   end do
  !   l=l_min
  !   do irho=irho_min,irho_max
  !     bra_x_p_nl(irho,n,l)=bra_x_p_nl(irho,n,l)&
  !     -(0.0d0,1.0d0)/(2.0d0*rho(irho))*dble(l+1)*p_nl(irho,n,l+1)
  !     bra_x_phi_nl(irho,n,l)=bra_x_phi_nl(irho,n,l)&
  !     -(0.0d0,1.0d0)/(2.0d0*rho(irho))*dble(l+1)*phi_nl(irho,n,l+1)
  !   end do
  !   l=l_max
  !   do irho=irho_min,irho_max
  !     bra_x_p_nl(irho,n,l)=bra_x_p_nl(irho,n,l)&
  !     -(0.0d0,1.0d0)/(2.0d0*rho(irho))*dble(l-1)*p_nl(irho,n,l-1)
  !     bra_x_phi_nl(irho,n,l)=bra_x_phi_nl(irho,n,l)&
  !     -(0.0d0,1.0d0)/(2.0d0*rho(irho))*dble(l-1)*phi_nl(irho,n,l-1)
  !   end do
  ! end if

  do l=l_min+1,l_max-1
    do n=sn,en
      if(n<=1)then
        do irho=irho_min,irho_max
          bra_x_p_nl(irho,n,l)=bra_x_p_nl(irho,n,l)&
          -(0.0d0,1.0d0)/(2.0d0*rho(irho))&
          *((dble(l+1)-l_shift(irho,n))*p_nl(irho,n,l+1)&
          +(dble(l-1)-l_shift(irho,n))*p_nl(irho,n,l-1))
          bra_x_phi_nl(irho,n,l)=bra_x_phi_nl(irho,n,l)&
          -(0.0d0,1.0d0)/(2.0d0*rho(irho))&
          *((dble(l+1)-l_shift(irho,n))*phi_nl(irho,n,l+1)&
          +(dble(l-1)-l_shift(irho,n))*phi_nl(irho,n,l-1))
        end do
      end if
    end do
  end do
  l=l_min
  do n=sn,en
    if(n<=1)then
      do irho=irho_min,irho_max
        bra_x_p_nl(irho,n,l)=bra_x_p_nl(irho,n,l)&
        -(0.0d0,1.0d0)/(2.0d0*rho(irho))&
        *(dble(l+1)-l_shift(irho,n))*p_nl(irho,n,l+1)
        bra_x_phi_nl(irho,n,l)=bra_x_phi_nl(irho,n,l)&
        -(0.0d0,1.0d0)/(2.0d0*rho(irho))&
        *(dble(l+1)-l_shift(irho,n))*phi_nl(irho,n,l+1)
      end do
    end if
  end do
  l=l_max
  do n=sn,en
    if(n<=1)then
      do irho=irho_min,irho_max
        bra_x_p_nl(irho,n,l)=bra_x_p_nl(irho,n,l)&
        -(0.0d0,1.0d0)/(2.0d0*rho(irho))&
        *(dble(l-1)-l_shift(irho,n))*p_nl(irho,n,l-1)
        bra_x_phi_nl(irho,n,l)=bra_x_phi_nl(irho,n,l)&
        -(0.0d0,1.0d0)/(2.0d0*rho(irho))&
        *(dble(l-1)-l_shift(irho,n))*phi_nl(irho,n,l-1)
      end do
    end if
  end do

  ! do l=l_min+1,l_max-1
  !   do n=sn,en
  !     do irho=irho_min,irho_max
  !       bra_x_p_nl(irho,n,l)=bra_x_p_nl(irho,n,l)&
  !       -(0.0d0,1.0d0)/(2.0d0*rho(irho))&
  !       *((dble(l+1)-l_shift(irho,n))*p_nl(irho,n,l+1)&
  !       +(dble(l-1)-l_shift(irho,n))*p_nl(irho,n,l-1))
  !       bra_x_phi_nl(irho,n,l)=bra_x_phi_nl(irho,n,l)&
  !       -(0.0d0,1.0d0)/(2.0d0*rho(irho))&
  !       *((dble(l+1)-l_shift(irho,n))*phi_nl(irho,n,l+1)&
  !       +(dble(l-1)-l_shift(irho,n))*phi_nl(irho,n,l-1))
  !     end do
  !   end do
  ! end do
  ! l=l_min
  ! do n=sn,en
  !   do irho=irho_min,irho_max
  !     bra_x_p_nl(irho,n,l)=bra_x_p_nl(irho,n,l)&
  !     -(0.0d0,1.0d0)/(2.0d0*rho(irho))&
  !     *(dble(l+1)-l_shift(irho,n))*p_nl(irho,n,l+1)
  !     bra_x_phi_nl(irho,n,l)=bra_x_phi_nl(irho,n,l)&
  !     -(0.0d0,1.0d0)/(2.0d0*rho(irho))&
  !     *(dble(l+1)-l_shift(irho,n))*phi_nl(irho,n,l+1)
  !   end do
  ! end do
  ! l=l_max
  ! do n=sn,en
  !   do irho=irho_min,irho_max
  !     bra_x_p_nl(irho,n,l)=bra_x_p_nl(irho,n,l)&
  !     -(0.0d0,1.0d0)/(2.0d0*rho(irho))&
  !     *(dble(l-1)-l_shift(irho,n))*p_nl(irho,n,l-1)
  !     bra_x_phi_nl(irho,n,l)=bra_x_phi_nl(irho,n,l)&
  !     -(0.0d0,1.0d0)/(2.0d0*rho(irho))&
  !     *(dble(l-1)-l_shift(irho,n))*phi_nl(irho,n,l-1)
  !   end do
  ! end do

end subroutine linear_terms
! ------------------------------------------------------------------------------
subroutine nonlinear_terms
! ------------------------------------------------------------------------------
  use coordinate
  use field
  use mpi_global
  implicit none
  include 'fftw3.f'

  do iz=0,numdata_z-1
    do n=sn,en
      do irho=irho_min,irho_max
        U(irho)=U_n(irho,n,iz)
        phi(irho)=phi_n(irho,n,iz)
        psi(irho)=psi_n(irho,n,iz)
        j(irho)=j_n(irho,n,iz)
        p(irho)=p_n(irho,n,iz)
      end do

      call del_rho_spec(U,del_rho_U)
      call del_rho_spec(phi,del_rho_phi)
      call del_rho_spec(psi,del_rho_psi)
      call del_rho_spec(j,del_rho_j)
      call del_rho_spec(p,del_rho_p)

      do irho=irho_min,irho_max
        del_alpha_U(irho)=(0.0d0,1.0d0)*k_theta(irho,n)*U(irho)
        del_alpha_phi(irho)=(0.0d0,1.0d0)*k_theta(irho,n)*phi(irho)
        del_alpha_psi(irho)=(0.0d0,1.0d0)*k_theta(irho,n)*psi(irho)
        del_alpha_j(irho)=(0.0d0,1.0d0)*k_theta(irho,n)*j(irho)
        del_alpha_p(irho)=(0.0d0,1.0d0)*k_theta(irho,n)*p(irho)

        del_rho_U_n(irho,n,iz)=del_rho_U(irho)
        del_alpha_U_n(irho,n,iz)=del_alpha_U(irho)
        del_rho_phi_n(irho,n,iz)=del_rho_phi(irho)
        del_alpha_phi_n(irho,n,iz)=del_alpha_phi(irho)
        del_rho_psi_n(irho,n,iz)=del_rho_psi(irho)
        del_alpha_psi_n(irho,n,iz)=del_alpha_psi(irho)
        del_rho_j_n(irho,n,iz)=del_rho_j(irho)
        del_alpha_j_n(irho,n,iz)=del_alpha_j(irho)
        del_rho_p_n(irho,n,iz)=del_rho_p(irho)
        del_alpha_p_n(irho,n,iz)=del_alpha_p(irho)
      end do ! : irho
    end do ! : n
  end do ! : iz

  call transpose(del_rho_U_n,del_rho_U_n_trans)
  call transpose(del_alpha_U_n,del_alpha_U_n_trans)
  call transpose(del_rho_phi_n,del_rho_phi_n_trans)
  call transpose(del_alpha_phi_n,del_alpha_phi_n_trans)
  call transpose(del_rho_psi_n,del_rho_psi_n_trans)
  call transpose(del_alpha_psi_n,del_alpha_psi_n_trans)
  call transpose(del_rho_j_n,del_rho_j_n_trans)
  call transpose(del_alpha_j_n,del_alpha_j_n_trans)
  call transpose(del_rho_p_n,del_rho_p_n_trans)
  call transpose(del_alpha_p_n,del_alpha_p_n_trans)

  do iz=0,numdata_z-1
    do irho=sirho,eirho
      do ialpha=0,numdata_alpha/2
        f1(ialpha)=0.0d0
        f2(ialpha)=0.0d0
        f3(ialpha)=0.0d0
        f4(ialpha)=0.0d0
        f5(ialpha)=0.0d0
        f6(ialpha)=0.0d0
        f7(ialpha)=0.0d0
        f8(ialpha)=0.0d0
        f9(ialpha)=0.0d0
        f10(ialpha)=0.0d0
      end do

      do n=0,n_max
        f1(n)=del_rho_U_n_trans(n,irho,iz)
        f2(n)=del_alpha_U_n_trans(n,irho,iz)
        f3(n)=del_rho_phi_n_trans(n,irho,iz)
        f4(n)=del_alpha_phi_n_trans(n,irho,iz)
        f5(n)=del_rho_psi_n_trans(n,irho,iz)
        f6(n)=del_alpha_psi_n_trans(n,irho,iz)
        f7(n)=del_rho_j_n_trans(n,irho,iz)
        f8(n)=del_alpha_j_n_trans(n,irho,iz)
        f9(n)=del_rho_p_n_trans(n,irho,iz)
        f10(n)=del_alpha_p_n_trans(n,irho,iz)
      end do

      call dfftw_execute_dft_c2r(plan_c2r,f1,f1_idft)
      call dfftw_execute_dft_c2r(plan_c2r,f2,f2_idft)
      call dfftw_execute_dft_c2r(plan_c2r,f3,f3_idft)
      call dfftw_execute_dft_c2r(plan_c2r,f4,f4_idft)
      call dfftw_execute_dft_c2r(plan_c2r,f5,f5_idft)
      call dfftw_execute_dft_c2r(plan_c2r,f6,f6_idft)
      call dfftw_execute_dft_c2r(plan_c2r,f7,f7_idft)
      call dfftw_execute_dft_c2r(plan_c2r,f8,f8_idft)
      call dfftw_execute_dft_c2r(plan_c2r,f9,f9_idft)
      call dfftw_execute_dft_c2r(plan_c2r,f10,f10_idft)

      do ialpha=0,numdata_alpha-1
        bra_phi_U_non_real(ialpha)=f3_idft(ialpha)*f2_idft(ialpha)&
        -f4_idft(ialpha)*f1_idft(ialpha)
        nab_para_j_non_real(ialpha)=f6_idft(ialpha)*f7_idft(ialpha)&
        -f5_idft(ialpha)*f8_idft(ialpha)
        nab_para_phi_non_real(ialpha)=f6_idft(ialpha)*f3_idft(ialpha)&
        -f5_idft(ialpha)*f4_idft(ialpha)
        bra_phi_p_non_real(ialpha)=f3_idft(ialpha)*f10_idft(ialpha)&
        -f4_idft(ialpha)*f9_idft(ialpha)
      end do

      do ialpha=0,numdata_alpha-1
        f1_idft(ialpha)=bra_phi_U_non_real(ialpha)
        f2_idft(ialpha)=nab_para_j_non_real(ialpha)
        f3_idft(ialpha)=nab_para_phi_non_real(ialpha)
        f4_idft(ialpha)=bra_phi_p_non_real(ialpha)
      end do

      call dfftw_execute_dft_r2c(plan_r2c,f1_idft,f1)
      call dfftw_execute_dft_r2c(plan_r2c,f2_idft,f2)
      call dfftw_execute_dft_r2c(plan_r2c,f3_idft,f3)
      call dfftw_execute_dft_r2c(plan_r2c,f4_idft,f4)

      do n=0,n_max
        bra_phi_U_non_n_trans(n,irho,iz)=f1(n)*inumdata_alpha
        nab_para_j_non_n_trans(n,irho,iz)=f2(n)*inumdata_alpha
        nab_para_phi_non_n_trans(n,irho,iz)=f3(n)*inumdata_alpha
        bra_phi_p_non_n_trans(n,irho,iz)=f4(n)*inumdata_alpha
      end do
    end do ! : irho
  end do ! : iz

  call itranspose(bra_phi_U_non_n_trans,bra_phi_U_non_n)
  call itranspose(nab_para_j_non_n_trans,nab_para_j_non_n)
  call itranspose(nab_para_phi_non_n_trans,nab_para_phi_non_n)
  call itranspose(bra_phi_p_non_n_trans,bra_phi_p_non_n)

  call dft_z(bra_phi_U_non_n,bra_phi_U_non_nl)
  call dft_z(nab_para_j_non_n,nab_para_j_non_nl)
  call dft_z(nab_para_phi_non_n,nab_para_phi_non_nl)
  call dft_z(bra_phi_p_non_n,bra_phi_p_non_nl)

end subroutine nonlinear_terms
! ------------------------------------------------------------------------------
subroutine phi_and_j
! ------------------------------------------------------------------------------
  use coordinate
  use field
  use mpi_global
  implicit none

  if(myid==0)then
    call conjugate_0l(U_nl)
    call conjugate_0l(psi_nl)
    call conjugate_0l(p_nl)
  end if

  call idft_z(U_nl,U_n)
  call idft_z(psi_nl,psi_n)
  call idft_z(p_nl,p_n)

  if(myid==0)then
    call conjugate_0(U_n)
    call conjugate_0(psi_n)
    call conjugate_0(p_n)
  end if

  do iz=0,numdata_z-1
    do n=sn,en
      do irho=irho_min,irho_max
        U(irho)=U_n(irho,n,iz)
        psi(irho)=psi_n(irho,n,iz)
      end do

      ! . definition of U
      call poisson_solver(U,phi)
      ! . Ampere's law
      call lap_perp_spec(-psi,j)

      do irho=irho_min,irho_max
        phi_n(irho,n,iz)=phi(irho)
        j_n(irho,n,iz)=j(irho)
      end do
    end do ! : n
  end do ! : iz

  if(myid==0)then
    call conjugate_0(phi_n)
    call conjugate_0(j_n)
  end if

  call dft_z(phi_n,phi_nl)
  call dft_z(j_n,j_nl)

  if(myid==0)then
    call conjugate_0l(phi_nl)
    call conjugate_0l(j_nl)
  end if

end subroutine phi_and_j
! ------------------------------------------------------------------------------
subroutine boundary_conditions
! ------------------------------------------------------------------------------
  use coordinate
  use field
  implicit none

  U(irho_min)=0.0d0
  U(irho_max)=0.0d0
  psi(irho_min)=0.0d0
  psi(irho_max)=0.0d0
  p(irho_min)=0.0d0
  p(irho_max)=0.0d0

end subroutine boundary_conditions
! ------------------------------------------------------------------------------
subroutine del_rho_real(a,del_rho_a)
! ------------------------------------------------------------------------------
  use coordinate
  implicit none
  double precision::a(irho_min:irho_max)
  double precision::del_rho_a(irho_min:irho_max)

  do irho=irho_min+2,irho_max-2
    del_rho_a(irho)=(-a(irho+2)+8.0d0*a(irho+1)&
    -8.0d0*a(irho-1)+a(irho-2))/(12.0d0*drho)
  end do
  irho=irho_min+1
  del_rho_a(irho)=(a(irho+1)-a(irho-1))/(2.0d0*drho)
  irho=irho_max-1
  del_rho_a(irho)=(a(irho+1)-a(irho-1))/(2.0d0*drho)
  irho=irho_min
  del_rho_a(irho)=(-3.0d0*a(irho)+4.0d0*a(irho+1)&
  -a(irho+2))/(2.0d0*drho)
  irho=irho_max
  del_rho_a(irho)=(3.0d0*a(irho)-4.0d0*a(irho-1)&
  +a(irho-2))/(2.0d0*drho)

end subroutine del_rho_real
! ------------------------------------------------------------------------------
subroutine del2_rho_real(a,del2_rho_a)
! ------------------------------------------------------------------------------
  use coordinate
  implicit none
  double precision::a(irho_min:irho_max)
  double precision::del2_rho_a(irho_min:irho_max)

  do irho=irho_min+2,irho_max-2
    del2_rho_a(irho)=(-a(irho+2)+16.0d0*a(irho+1)-30.0d0*a(irho)&
    +16.0d0*a(irho-1)-a(irho-2))/(12.0d0*drho_sq)
  end do
  irho=irho_min+1
  del2_rho_a(irho)=(a(irho+1)-2.0d0*a(irho)+a(irho-1))/drho_sq
  irho=irho_max-1
  del2_rho_a(irho)=(a(irho+1)-2.0d0*a(irho)+a(irho-1))/drho_sq
  irho=irho_min
  del2_rho_a(irho)=(a(irho)-2.0d0*a(irho+1)+a(irho+2))/drho_sq
  irho=irho_max
  del2_rho_a(irho)=(a(irho)-2.0d0*a(irho-1)+a(irho-2))/drho_sq

end subroutine del2_rho_real
! ------------------------------------------------------------------------------
subroutine lap_rho_real(a,lap_rho_a)
! ------------------------------------------------------------------------------
  use coordinate
  implicit none
  double precision::a(irho_min:irho_max)
  double precision::lap_rho_a(irho_min:irho_max)
  double precision::del_rho_a(irho_min:irho_max)
  double precision::del2_rho_a(irho_min:irho_max)

  call del_rho_real(a,del_rho_a)
  call del2_rho_real(a,del2_rho_a)

  do irho=irho_min,irho_max
    lap_rho_a(irho)=del2_rho_a(irho)+del_rho_a(irho)/rho(irho)
  end do

end subroutine lap_rho_real
! ------------------------------------------------------------------------------
subroutine del_rho_spec(a,del_rho_a)
! ------------------------------------------------------------------------------
  use coordinate
  implicit none
  complex(kind(0d0))::a(irho_min:irho_max)
  complex(kind(0d0))::del_rho_a(irho_min:irho_max)

  do irho=irho_min+2,irho_max-2
    del_rho_a(irho)=(-a(irho+2)+8.0d0*a(irho+1)&
    -8.0d0*a(irho-1)+a(irho-2))/(12.0d0*drho)
  end do
  irho=irho_min+1
  del_rho_a(irho)=(a(irho+1)-a(irho-1))/(2.0d0*drho)
  irho=irho_max-1
  del_rho_a(irho)=(a(irho+1)-a(irho-1))/(2.0d0*drho)
  irho=irho_min
  del_rho_a(irho)=(-3.0d0*a(irho)+4.0d0*a(irho+1)&
  -a(irho+2))/(2.0d0*drho)
  irho=irho_max
  del_rho_a(irho)=(3.0d0*a(irho)-4.0d0*a(irho-1)&
  +a(irho-2))/(2.0d0*drho)

end subroutine del_rho_spec
! ------------------------------------------------------------------------------
subroutine del2_rho_spec(a,del2_rho_a)
! ------------------------------------------------------------------------------
  use coordinate
  implicit none
  complex(kind(0d0))::a(irho_min:irho_max)
  complex(kind(0d0))::del2_rho_a(irho_min:irho_max)

  do irho=irho_min+2,irho_max-2
    del2_rho_a(irho)=(-a(irho+2)+16.0d0*a(irho+1)-30.0d0*a(irho)&
    +16.0d0*a(irho-1)-a(irho-2))/(12.0d0*drho_sq)
  end do
  irho=irho_min+1
  del2_rho_a(irho)=(a(irho+1)-2.0d0*a(irho)+a(irho-1))/drho_sq
  irho=irho_max-1
  del2_rho_a(irho)=(a(irho+1)-2.0d0*a(irho)+a(irho-1))/drho_sq
  irho=irho_min
  del2_rho_a(irho)=(a(irho)-2.0d0*a(irho+1)+a(irho+2))/drho_sq
  irho=irho_max
  del2_rho_a(irho)=(a(irho)-2.0d0*a(irho-1)+a(irho-2))/drho_sq

end subroutine del2_rho_spec
! ------------------------------------------------------------------------------
subroutine bra_x_spec(a,bra_x_a)
! ------------------------------------------------------------------------------
  use coordinate
  use field
  implicit none
  complex(kind(0d0))::a(irho_min:irho_max)
  complex(kind(0d0))::bra_x_a(irho_min:irho_max)
  complex(kind(0d0))::del_rho_a(irho_min:irho_max)

  call del_rho_spec(a,del_rho_a)

  do irho=irho_min,irho_max
    bra_x_a(irho)=fact_del_rho_bra_x(iz)*del_rho_a(irho)&
    +fact_bra_x(irho,n,iz)*a(irho)
  end do

end subroutine bra_x_spec
! ------------------------------------------------------------------------------
subroutine lap_perp_spec(a,lap_perp_a)
! ------------------------------------------------------------------------------
  use coordinate
  use field
  implicit none
  complex(kind(0d0))::a(irho_min:irho_max)
  complex(kind(0d0))::lap_perp_a(irho_min:irho_max)
  complex(kind(0d0))::del_rho_a(irho_min:irho_max)
  complex(kind(0d0))::del2_rho_a(irho_min:irho_max)

  call del_rho_spec(a,del_rho_a)
  call del2_rho_spec(a,del2_rho_a)

  do irho=irho_min,irho_max
    lap_perp_a(irho)=del2_rho_a(irho)&
    +fact_del_rho_lap_perp(irho,n,iz)*del_rho_a(irho)&
    +fact_lap_perp(irho,n,iz)*a(irho)
  end do

end subroutine lap_perp_spec
! ------------------------------------------------------------------------------
subroutine poisson_solver(a,ilap_perp_a)
! ------------------------------------------------------------------------------
  use coordinate
  use field
  implicit none
  complex(kind(0d0))::a(irho_min:irho_max)
  complex(kind(0d0))::ilap_perp_a(irho_min:irho_max)

  do irho=irho_min+1,irho_max-1
    do ilu=1,3
      coe_lu(ilu,irho)=coe_lu_n(ilu,irho,n,iz)
    end do

    rhs_lu(irho)=drho_sq*a(irho)
  end do ! : irho

  do irho=irho_min+2,irho_max-1
    rhs_lu(irho)=rhs_lu(irho)-coe_lu(1,irho)*rhs_lu(irho-1)
  end do

  irho=irho_max-1
  rhs_lu(irho)=rhs_lu(irho)/coe_lu(2,irho)
  do irho=irho_max-2,irho_min+1,-1
    rhs_lu(irho)=(rhs_lu(irho)&
    -coe_lu(3,irho)*rhs_lu(irho+1))/coe_lu(2,irho)
  end do

  do irho=irho_min+1,irho_max-1
    ilap_perp_a(irho)=rhs_lu(irho)
  end do

  ilap_perp_a(irho_min)=0.0d0
  ilap_perp_a(irho_max)=0.0d0

end subroutine poisson_solver
! ------------------------------------------------------------------------------
subroutine transpose(a_n,a_n_trans)
! ------------------------------------------------------------------------------
  use coordinate
  use mpi_global
  use mpi
  implicit none
  complex(kind(0d0))::a_n(irho_min:irho_max,sn:en,0:numdata_z-1)
  complex(kind(0d0))::a_n_trans(0:n_max,sirho:eirho,0:numdata_z-1)

  do iz=0,numdata_z-1
    idata=0
    do iprocs=0,numprocs-1
      do n=sn,en
        do irho=sirho0(iprocs),eirho0(iprocs)
          a_send_trans(idata)=a_n(irho,n,iz)
          idata=idata+1
        end do
      end do
    end do

    call mpi_alltoallv(a_send_trans,scounts_trans,sdispls_trans,mpi_double_complex,&
    a_recv_trans,rcounts_trans,rdispls_trans,mpi_double_complex,&
    mpi_comm_world,ierr)

    idata=0
    do iprocs=0,numprocs-1
      do n=sn0(iprocs),en0(iprocs)
        do irho=sirho,eirho
          a_n_trans(n,irho,iz)=a_recv_trans(idata)
          idata=idata+1
        end do
      end do
    end do
  end do ! : iz

end subroutine transpose
! ------------------------------------------------------------------------------
subroutine itranspose(a_n_trans,a_n)
! ------------------------------------------------------------------------------
  use coordinate
  use mpi_global
  use mpi
  implicit none
  complex(kind(0d0))::a_n_trans(0:n_max,sirho:eirho,0:numdata_z-1)
  complex(kind(0d0))::a_n(irho_min:irho_max,sn:en,0:numdata_z-1)

  do iz=0,numdata_z-1
    idata=0
    do iprocs=0,numprocs-1
      do irho=sirho,eirho
        do n=sn0(iprocs),en0(iprocs)
          a_recv_trans(idata)=a_n_trans(n,irho,iz)
          idata=idata+1
        end do
      end do
    end do

    call mpi_alltoallv(a_recv_trans,rcounts_trans,rdispls_trans,mpi_double_complex,&
    a_send_trans,scounts_trans,sdispls_trans,mpi_double_complex,&
    mpi_comm_world,ierr)

    idata=0
    do iprocs=0,numprocs-1
      do irho=sirho0(iprocs),eirho0(iprocs)
        do n=sn,en
          a_n(irho,n,iz)=a_send_trans(idata)
          idata=idata+1
        end do
      end do
    end do
  end do ! : iz

end subroutine itranspose
! ------------------------------------------------------------------------------
subroutine conjugate_0l(a_nl)
! ------------------------------------------------------------------------------
  use coordinate
  use mpi_global
  implicit none
  complex(kind(0d0))::a_nl(irho_min:irho_max,sn:en,l_min:l_max)

  n=0
  do l=1,l_max
    do irho=irho_min,irho_max
      a_nl(irho,-n,-l)=dble(a_nl(irho,n,l))&
      +(0.0d0,1.0d0)*dble((0.0d0,1.0d0)*a_nl(irho,n,l))
    end do
  end do
  l=0
  do irho=irho_min,irho_max
    a_nl(irho,-n,-l)=dble(a_nl(irho,n,l))
  end do

end subroutine conjugate_0l
! ------------------------------------------------------------------------------
subroutine conjugate_0(a_n)
! ------------------------------------------------------------------------------
  use coordinate
  use mpi_global
  implicit none
  complex(kind(0d0))::a_n(irho_min:irho_max,sn:en,0:numdata_z-1)

  n=0
  do iz=0,numdata_z-1
    do irho=irho_min,irho_max
      a_n(irho,-n,iz)=dble(a_n(irho,n,iz))
    end do
  end do

end subroutine conjugate_0
! ------------------------------------------------------------------------------
subroutine idft_z(a_nl,a_n)
! ------------------------------------------------------------------------------
  use coordinate
  use field
  use mpi_global
  implicit none
  include 'fftw3.f'
  complex(kind(0d0))::a_nl(irho_min:irho_max,sn:en,l_min:l_max)
  complex(kind(0d0))::a_n(irho_min:irho_max,sn:en,0:numdata_z-1)

  do n=sn,en
    do irho=irho_min,irho_max
      do iz=0,numdata_z-1
        f(iz)=0.0d0
      end do

      do l=0,l_max
        f(l)=a_nl(irho,n,l)
      end do
      do l=l_min,-1
        f(l+numdata_z)=a_nl(irho,n,l)
      end do

      call dfftw_execute_dft(plan_back,f,f)

      do iz=0,numdata_z-1
        a_n(irho,n,iz)=f(iz)
      end do
    end do ! : irho
  end do ! : n

  do iz=0,numdata_z-1
    do n=sn,en
      do irho=irho_min,irho_max
        a_n(irho,n,iz)=a_n(irho,n,iz)&
        *ishifter(irho,n,iz)
      end do
    end do
  end do

end subroutine idft_z
! ------------------------------------------------------------------------------
subroutine dft_z(a_n,a_nl)
! ------------------------------------------------------------------------------
  use coordinate
  use field
  use mpi_global
  implicit none
  include 'fftw3.f'
  complex(kind(0d0))::a_n(irho_min:irho_max,sn:en,0:numdata_z-1)
  complex(kind(0d0))::a_nl(irho_min:irho_max,sn:en,l_min:l_max)
  complex(kind(0d0))::a_n_shifted(irho_min:irho_max,sn:en,0:numdata_z-1)

  do iz=0,numdata_z-1
    do n=sn,en
      do irho=irho_min,irho_max
        a_n_shifted(irho,n,iz)=a_n(irho,n,iz)&
        *shifter(irho,n,iz)
      end do
    end do
  end do

  do n=sn,en
    do irho=irho_min,irho_max
      do iz=0,numdata_z-1
        f(iz)=a_n_shifted(irho,n,iz)
      end do

      call dfftw_execute_dft(plan_for,f,f)

      do l=0,numdata_z/2-1
        f_dft(l)=f(l)*inumdata_z
      end do
      do l=numdata_z/2,numdata_z-1
        f_dft(l-numdata_z)=f(l)*inumdata_z
      end do

      do l=l_min,l_max
        a_nl(irho,n,l)=f_dft(l)
      end do
    end do ! : irho
  end do ! : n

end subroutine dft_z
! ------------------------------------------------------------------------------
subroutine idft_theta(a_nl,a_n_cylin)
! ------------------------------------------------------------------------------
  use coordinate
  use field
  use mpi_global
  implicit none
  include 'fftw3.f'
  complex(kind(0d0))::a_nl(irho_min:irho_max,sn:en,l_min:l_max)
  complex(kind(0d0))::a_n_cylin(irho_min:irho_max,sn:en,0:numdata_theta-1)

  do n=sn,en
    do irho=irho_min,irho_max
      do itheta=0,numdata_theta-1
        f_cylin(itheta)=0.0d0
      end do

      do l=0,l_max
        f_cylin(l)=a_nl(irho,n,l)
      end do
      do l=l_min,-1
        f_cylin(l+numdata_theta)=a_nl(irho,n,l)
      end do

      call dfftw_execute_dft(plan_back_cylin,f_cylin,f_cylin)

      do itheta=0,numdata_theta-1
        a_n_cylin(irho,n,itheta)=f_cylin(itheta)
      end do
    end do ! : irho
  end do ! : n

  do itheta=0,numdata_theta-1
    do n=sn,en
      do irho=irho_min,irho_max
        a_n_cylin(irho,n,itheta)=a_n_cylin(irho,n,itheta)&
        *ishifter_cylin(irho,n,itheta)
      end do
    end do
  end do

end subroutine idft_theta
! ------------------------------------------------------------------------------
subroutine linear_analysis
! ------------------------------------------------------------------------------
  use time
  use coordinate
  use field
  use mpi_global
  use mpi
  implicit none

  call idft_z(p_old_nl,p_old_n)

  if(it==it_start)then
    do n=sn,en
      growth(n)=0.0d0
      freq(n)=0.0d0
    end do
  else
    iz=numdata_z/2
    irho=irho_peak
    do n=sn,en
      p(irho)=p_n(irho,n,iz)
      p_old(irho)=p_old_n(irho,n,iz)

      growth(n)=dble(log(p(irho)/p_old(irho)))/dt
      freq(n)=dble((0.0d0,1.0d0)*log(p(irho)/p_old(irho)))/dt
    end do ! : n
  end if ! : it==it_start

  do n=sn,en
    a_send_gather(n)=growth(n)
  end do
  call mpi_allgatherv(a_send_gather,scount_gather,mpi_double_complex,&
  a_recv_gather,rcounts_gather,rdispls_gather,mpi_double_complex,&
  mpi_comm_world,ierr)
  do n=0,n_max
    growth_gather(n)=a_recv_gather(n)
  end do

  do n=sn,en
    a_send_gather(n)=freq(n)
  end do
  call mpi_allgatherv(a_send_gather,scount_gather,mpi_double_complex,&
  a_recv_gather,rcounts_gather,rdispls_gather,mpi_double_complex,&
  mpi_comm_world,ierr)
  do n=0,n_max
    freq_gather(n)=a_recv_gather(n)
  end do

  if(myid==0)then
    write(filename,'("../data/growth_freq_",i3.3,".dat")')it/it_skip
    write(6,*)filename
    open(unit=10,file=filename,form='formatted')
    write(10,fmt='(3a16)')'n','growth','freq'
    do n=0,n_max
      write(10,fmt='(3e16.4)')dble(n),growth_gather(n),freq_gather(n)
    end do
    close(10)
  end if

  iz=numdata_z/2
  do n=sn,en
    if(n==n_max/2.or.n==n_max/4)then
      do irho=irho_min,irho_max
        p_amp(irho)=abs(p_n(irho,n,iz))
      end do

      p_irho_peak=p_amp(irho_peak)
      do irho=irho_min,irho_max
        p(irho)=p_amp(irho)/p_irho_peak&
        *exp(-(0.0d0,1.0d0)*freq(n)*t)
      end do

      write(filename,'("../data/eigfunc_",i3.3,"_",i3.3,".dat")')n,it/it_skip
      write(6,*)filename
      open(unit=10*(1+myid),file=filename,form='formatted')
      write(10*(1+myid),fmt='(3a16)')'rho','Re(p_n)','Im(p_n)'
      do irho=irho_min,irho_max
        write(10*(1+myid),fmt='(3e16.4)')rho(irho),&
        dble(p(irho)),dble(-(0.0d0,1.0d0)*p(irho))
      end do
      close(10*(1+myid))
    end if ! : n==n_max/2.or.n==n_max/4
  end do ! : n

end subroutine linear_analysis
! ------------------------------------------------------------------------------
subroutine average
! ------------------------------------------------------------------------------
  use time
  use coordinate
  use field
  use mpi_global
  implicit none

  do irho=irho_min,irho_max
    U_avg(irho)=U_eq(irho)+dble(U_nl(irho,0,0))
    phi_avg(irho)=phi_eq(irho)+dble(phi_nl(irho,0,0))
    j_avg(irho)=j_eq(irho)+dble(j_nl(irho,0,0))
    p_avg(irho)=p_eq(irho)+dble(p_nl(irho,0,0))
  end do

  call del_rho_real(U_avg,del_rho_U_avg)
  call del_rho_real(phi_avg,del_rho_phi_avg)
  call del_rho_real(j_avg,del_rho_j_avg)
  call del_rho_real(p_avg,del_rho_p_avg)

  write(filename,'("../data/avg_profile_",i3.3,".dat")')it/it_skip
  write(6,*)filename
  open(unit=10,file=filename,form='formatted')
  write(10,fmt='(5a16)')'rho','U_avg','phi_avg','j_avg','p_avg'
  do irho=irho_min,irho_max
    write(10,fmt='(5e16.4)')rho(irho),&
    U_avg(irho),phi_avg(irho),j_avg(irho),p_avg(irho)
  end do
  close(10)

  write(filename,'("../data/avg_grad_profile_",i3.3,".dat")')it/it_skip
  write(6,*)filename
  open(unit=10,file=filename,form='formatted')
  write(10,fmt='(5a16)')'rho','del_rho_U_avg','del_rho_phi_avg',&
  'del_rho_j_avg','del_rho_p_avg'
  do irho=irho_min,irho_max
    write(10,fmt='(5e16.4)')rho(irho),del_rho_U_avg(irho),&
    del_rho_phi_avg(irho),del_rho_j_avg(irho),del_rho_p_avg(irho)
  end do
  close(10)

end subroutine average
! ------------------------------------------------------------------------------
subroutine contours
! ------------------------------------------------------------------------------
  use time
  use coordinate
  use field
  use mpi_global
  use mpi
  implicit none

  ! call idft_theta(U_nl,U_n_cylin)
  ! call idft_theta(phi_nl,phi_n_cylin)
  ! call idft_theta(psi_nl,psi_n_cylin)
  ! call idft_theta(j_nl,j_n_cylin)
  call idft_theta(p_nl,p_n_cylin)

  do itheta=0,numdata_theta-1
    do irho=irho_min,irho_max
      ! U_real_cylin(irho,itheta)=0.0d0
      ! phi_real_cylin(irho,itheta)=0.0d0
      ! psi_real_cylin(irho,itheta)=0.0d0
      ! j_real_cylin(irho,itheta)=0.0d0
      p_real_cylin(irho,itheta)=0.0d0
    end do

    do n=sn,en
      if(n==0)then
        do irho=irho_min,irho_max
          ! U_real_cylin(irho,itheta)=U_real_cylin(irho,itheta)&
          ! +dble(U_n_cylin(irho,n,itheta))
          ! phi_real_cylin(irho,itheta)=phi_real_cylin(irho,itheta)&
          ! +dble(phi_n_cylin(irho,n,itheta))
          ! psi_real_cylin(irho,itheta)=psi_real_cylin(irho,itheta)&
          ! +dble(psi_n_cylin(irho,n,itheta))
          ! j_real_cylin(irho,itheta)=j_real_cylin(irho,itheta)&
          ! +dble(j_n_cylin(irho,n,itheta))
          p_real_cylin(irho,itheta)=p_real_cylin(irho,itheta)&
          +dble(p_n_cylin(irho,n,itheta))
        end do
      else
        do irho=irho_min,irho_max
          ! U_real_cylin(irho,itheta)=U_real_cylin(irho,itheta)&
          ! +2.0d0*dble(U_n_cylin(irho,n,itheta))
          ! phi_real_cylin(irho,itheta)=phi_real_cylin(irho,itheta)&
          ! +2.0d0*dble(phi_n_cylin(irho,n,itheta))
          ! psi_real_cylin(irho,itheta)=psi_real_cylin(irho,itheta)&
          ! +2.0d0*dble(psi_n_cylin(irho,n,itheta))
          ! j_real_cylin(irho,itheta)=j_real_cylin(irho,itheta)&
          ! +2.0d0*dble(j_n_cylin(irho,n,itheta))
          p_real_cylin(irho,itheta)=p_real_cylin(irho,itheta)&
          +2.0d0*dble(p_n_cylin(irho,n,itheta))
        end do
      end if
    end do
  end do ! : itheta

  do itheta=0,numdata_theta-1
    do irho=irho_min,irho_max
      ! a_send_op=U_real_cylin(irho,itheta)
      ! call mpi_allreduce(a_send_op,a_recv_op,1,mpi_double_precision,mpi_sum,&
      ! mpi_comm_world,ierr)
      ! U_real_cylin(irho,itheta)=a_recv_op

      ! a_send_op=phi_real_cylin(irho,itheta)
      ! call mpi_allreduce(a_send_op,a_recv_op,1,mpi_double_precision,mpi_sum,&
      ! mpi_comm_world,ierr)
      ! phi_real_cylin(irho,itheta)=a_recv_op

      ! a_send_op=psi_real_cylin(irho,itheta)
      ! call mpi_allreduce(a_send_op,a_recv_op,1,mpi_double_precision,mpi_sum,&
      ! mpi_comm_world,ierr)
      ! psi_real_cylin(irho,itheta)=a_recv_op

      ! a_send_op=j_real_cylin(irho,itheta)
      ! call mpi_allreduce(a_send_op,a_recv_op,1,mpi_double_precision,mpi_sum,&
      ! mpi_comm_world,ierr)
      ! j_real_cylin(irho,itheta)=a_recv_op

      a_send_op=p_real_cylin(irho,itheta)
      call mpi_allreduce(a_send_op,a_recv_op,1,mpi_double_precision,mpi_sum,&
      mpi_comm_world,ierr)
      p_real_cylin(irho,itheta)=a_recv_op
    end do ! : irho
  end do ! : itheta

  do irho=irho_min,irho_max
    ! U_real_cylin(irho,numdata_theta)=U_real_cylin(irho,0)
    ! phi_real_cylin(irho,numdata_theta)=phi_real_cylin(irho,0)
    ! psi_real_cylin(irho,numdata_theta)=psi_real_cylin(irho,0)
    ! j_real_cylin(irho,numdata_theta)=j_real_cylin(irho,0)
    p_real_cylin(irho,numdata_theta)=p_real_cylin(irho,0)
  end do

  if(myid==0)then
    ! write(filename,'("../data/U_contour_",i3.3,".dat")')it/it_skip
    ! write(6,*)filename
    ! open(unit=10,file=filename,form='formatted')
    ! write(10,fmt='(3a16)')'x','Z','U_real'
    ! do itheta=0,numdata_theta
    !   do irho=irho_min,irho_max
    !     write(10,fmt='(3e16.4)')rho(irho)*cos(theta(itheta)),&
    !     rho(irho)*sin(theta(itheta)),U_real_cylin(irho,itheta)
    !   end do
    !   write(10,*)''
    ! end do
    ! close(10)

    ! write(filename,'("../data/phi_contour_",i3.3,".dat")')it/it_skip
    ! write(6,*)filename
    ! open(unit=10,file=filename,form='formatted')
    ! write(10,fmt='(3a16)')'x','Z','phi_real'
    ! do itheta=0,numdata_theta
    !   do irho=irho_min,irho_max
    !     write(10,fmt='(3e16.4)')rho(irho)*cos(theta(itheta)),&
    !     rho(irho)*sin(theta(itheta)),phi_real_cylin(irho,itheta)
    !   end do
    !   write(10,*)''
    ! end do
    ! close(10)

    ! write(filename,'("../data/psi_contour_",i3.3,".dat")')it/it_skip
    ! write(6,*)filename
    ! open(unit=10,file=filename,form='formatted')
    ! write(10,fmt='(3a16)')'x','Z','psi_real'
    ! do itheta=0,numdata_theta
    !   do irho=irho_min,irho_max
    !     write(10,fmt='(3e16.4)')rho(irho)*cos(theta(itheta)),&
    !     rho(irho)*sin(theta(itheta)),psi_real_cylin(irho,itheta)
    !   end do
    !   write(10,*)''
    ! end do
    ! close(10)

    ! write(filename,'("../data/j_contour_",i3.3,".dat")')it/it_skip
    ! write(6,*)filename
    ! open(unit=10,file=filename,form='formatted')
    ! write(10,fmt='(3a16)')'x','Z','j_real'
    ! do itheta=0,numdata_theta
    !   do irho=irho_min,irho_max
    !     write(10,fmt='(3e16.4)')rho(irho)*cos(theta(itheta)),&
    !     rho(irho)*sin(theta(itheta)),j_real_cylin(irho,itheta)
    !   end do
    !   write(10,*)''
    ! end do
    ! close(10)

    write(filename,'("../data/p_contour_",i3.3,".dat")')it/it_skip
    write(6,*)filename
    open(unit=10,file=filename,form='formatted')
    write(10,fmt='(3a16)')'x','Z','p_real'
    do itheta=0,numdata_theta
      do irho=irho_min,irho_max
        write(10,fmt='(3e16.4)')rho(irho)*cos(theta(itheta)),&
        rho(irho)*sin(theta(itheta)),p_real_cylin(irho,itheta)
      end do
      write(10,*)''
    end do
    close(10)
  end if ! : myid==0

end subroutine contours
! ------------------------------------------------------------------------------
subroutine energy_balances
! ------------------------------------------------------------------------------
  use time
  use coordinate
  use parameter
  use field
  use mpi_global
  use mpi
  implicit none

  ! . kinetic energy
  do n=sn,en
    dE_phi_n(n)=0.0d0
  end do

  do l=l_min,l_max
    do n=sn,en
      do irho=irho_min,irho_max
        U(irho)=U_nl(irho,n,l)
        phi(irho)=phi_nl(irho,n,l)
        U_old(irho)=U_old_nl(irho,n,l)
        phi_old(irho)=phi_old_nl(irho,n,l)
      end do

      call conjugate(phi,phi_cc)
      call conjugate(phi_old,phi_old_cc)

      do irho=irho_min,irho_max
        dE_phi_nl_den(irho)=-(dble(phi_cc(irho)*U(irho))&
        -dble(phi_old_cc(irho)*U_old(irho)))/(2.0d0*dt)
      end do

      dE_phi_nl=0.0d0
      do irho=irho_min+1,irho_max
        dE_phi_nl=dE_phi_nl+drho/2.0d0&
        *(dE_phi_nl_den(irho-1)*rho(irho-1)&
        +dE_phi_nl_den(irho)*rho(irho))
      end do

      dE_phi_n(n)=dE_phi_n(n)+dE_phi_nl
    end do ! : n
  end do ! : l

  dE_phi=0.0d0
  do n=sn,en
    if(n==0)then
      dE_phi=dE_phi+dE_phi_n(n)
    else
      dE_phi=dE_phi+2.0d0*dE_phi_n(n)
    end if
  end do

  ! . linear energy sources
  do n=sn,en
    S_phi_dphi_eq_n(n)=0.0d0
    S_phi_q_n(n)=0.0d0
    S_phi_dj_eq_n(n)=0.0d0
    S_phi_x_p_n(n)=0.0d0
    S_phi_mu_n(n)=0.0d0
  end do

  do l=l_min,l_max
    do n=sn,en
      do irho=irho_min,irho_max
        U(irho)=U_nl(irho,n,l)
        phi(irho)=phi_nl(irho,n,l)
        psi(irho)=psi_nl(irho,n,l)
        j(irho)=j_nl(irho,n,l)
        bra_x_p(irho)=bra_x_p_nl(irho,n,l)
        lap_perp_U(irho)=lap_perp_U_nl(irho,n,l)
      end do

      call conjugate(phi,phi_cc)

      do irho=irho_min,irho_max
        S_phi_dphi_eq_nl_den(irho)=dble((0.0d0,1.0d0)*k_theta(irho,n)&
        *dphi_eq(irho)*phi_cc(irho)*U(irho))
        S_phi_q_nl_den(irho)=-dble((0.0d0,1.0d0)*k_para(irho,n,l)&
        *phi_cc(irho)*j(irho))
        S_phi_dj_eq_nl_den(irho)=-dble((0.0d0,1.0d0)*k_theta(irho,n)&
        *dj_eq(irho)*phi_cc(irho)*psi(irho))
        S_phi_x_p_nl_den(irho)=dble(phi_cc(irho)*bra_x_p(irho))
        S_phi_mu_nl_den(irho)=-mu*dble(phi_cc(irho)*lap_perp_U(irho))
      end do

      S_phi_dphi_eq_nl=0.0d0
      S_phi_q_nl=0.0d0
      S_phi_dj_eq_nl=0.0d0
      S_phi_x_p_nl=0.0d0
      S_phi_mu_nl=0.0d0
      do irho=irho_min+1,irho_max
        S_phi_dphi_eq_nl=S_phi_dphi_eq_nl+drho/2.0d0&
        *(S_phi_dphi_eq_nl_den(irho-1)*rho(irho-1)&
        +S_phi_dphi_eq_nl_den(irho)*rho(irho))
        S_phi_q_nl=S_phi_q_nl+drho/2.0d0&
        *(S_phi_q_nl_den(irho-1)*rho(irho-1)&
        +S_phi_q_nl_den(irho)*rho(irho))
        S_phi_dj_eq_nl=S_phi_dj_eq_nl+drho/2.0d0&
        *(S_phi_dj_eq_nl_den(irho-1)*rho(irho-1)&
        +S_phi_dj_eq_nl_den(irho)*rho(irho))
        S_phi_x_p_nl=S_phi_x_p_nl+drho/2.0d0&
        *(S_phi_x_p_nl_den(irho-1)*rho(irho-1)&
        +S_phi_x_p_nl_den(irho)*rho(irho))
        S_phi_mu_nl=S_phi_mu_nl+drho/2.0d0&
        *(S_phi_mu_nl_den(irho-1)*rho(irho-1)&
        +S_phi_mu_nl_den(irho)*rho(irho))
      end do

      S_phi_dphi_eq_n(n)=S_phi_dphi_eq_n(n)+S_phi_dphi_eq_nl
      S_phi_q_n(n)=S_phi_q_n(n)+S_phi_q_nl
      S_phi_dj_eq_n(n)=S_phi_dj_eq_n(n)+S_phi_dj_eq_nl
      S_phi_x_p_n(n)=S_phi_x_p_n(n)+S_phi_x_p_nl
      S_phi_mu_n(n)=S_phi_mu_n(n)+S_phi_mu_nl
    end do ! : n
  end do ! : l

  S_phi_dphi_eq=0.0d0
  S_phi_q=0.0d0
  S_phi_dj_eq=0.0d0
  S_phi_x_p=0.0d0
  S_phi_mu=0.0d0
  do n=sn,en
    if(n==0)then
      S_phi_dphi_eq=S_phi_dphi_eq+S_phi_dphi_eq_n(n)
      S_phi_q=S_phi_q+S_phi_q_n(n)
      S_phi_dj_eq=S_phi_dj_eq+S_phi_dj_eq_n(n)
      S_phi_x_p=S_phi_x_p+S_phi_x_p_n(n)
      S_phi_mu=S_phi_mu+S_phi_mu_n(n)
    else
      S_phi_dphi_eq=S_phi_dphi_eq+2.0d0*S_phi_dphi_eq_n(n)
      S_phi_q=S_phi_q+2.0d0*S_phi_q_n(n)
      S_phi_dj_eq=S_phi_dj_eq+2.0d0*S_phi_dj_eq_n(n)
      S_phi_x_p=S_phi_x_p+2.0d0*S_phi_x_p_n(n)
      S_phi_mu=S_phi_mu+2.0d0*S_phi_mu_n(n)
    end if
  end do

  ! . nonlinear energy sources
  do n=sn,en
    S_phi_phi_U_n(n)=0.0d0
    S_phi_psi_j_n(n)=0.0d0
  end do

  do l=l_min,l_max
    do n=sn,en
      do irho=irho_min,irho_max
        phi(irho)=phi_nl(irho,n,l)
        bra_phi_U_non(irho)=bra_phi_U_non_nl(irho,n,l)
        nab_para_j_non(irho)=nab_para_j_non_nl(irho,n,l)
      end do

      call conjugate(phi,phi_cc)

      do irho=irho_min,irho_max
        S_phi_phi_U_nl_den(irho)=dble(phi_cc(irho)*bra_phi_U_non(irho))
        S_phi_psi_j_nl_den(irho)=-dble(phi_cc(irho)*nab_para_j_non(irho))
      end do

      S_phi_phi_U_nl=0.0d0
      S_phi_psi_j_nl=0.0d0
      do irho=irho_min+1,irho_max
        S_phi_phi_U_nl=S_phi_phi_U_nl+drho/2.0d0&
        *(S_phi_phi_U_nl_den(irho-1)*rho(irho-1)&
        +S_phi_phi_U_nl_den(irho)*rho(irho))
        S_phi_psi_j_nl=S_phi_psi_j_nl+drho/2.0d0&
        *(S_phi_psi_j_nl_den(irho-1)*rho(irho-1)&
        +S_phi_psi_j_nl_den(irho)*rho(irho))
      end do

      S_phi_phi_U_n(n)=S_phi_phi_U_n(n)+S_phi_phi_U_nl
      S_phi_psi_j_n(n)=S_phi_psi_j_n(n)+S_phi_psi_j_nl
    end do ! : n
  end do ! : l

  S_phi_phi_U=0.0d0
  S_phi_psi_j=0.0d0
  do n=sn,en
    if(n==0)then
      S_phi_phi_U=S_phi_phi_U+S_phi_phi_U_n(n)
      S_phi_psi_j=S_phi_psi_j+S_phi_psi_j_n(n)
    else
      S_phi_phi_U=S_phi_phi_U+2.0d0*S_phi_phi_U_n(n)
      S_phi_psi_j=S_phi_psi_j+2.0d0*S_phi_psi_j_n(n)
    end if
  end do

  a_send_op=dE_phi
  call mpi_allreduce(a_send_op,a_recv_op,1,mpi_double_precision,mpi_sum,&
  mpi_comm_world,ierr)
  dE_phi=a_recv_op

  a_send_op=S_phi_dphi_eq
  call mpi_allreduce(a_send_op,a_recv_op,1,mpi_double_precision,mpi_sum,&
  mpi_comm_world,ierr)
  S_phi_dphi_eq=a_recv_op

  a_send_op=S_phi_q
  call mpi_allreduce(a_send_op,a_recv_op,1,mpi_double_precision,mpi_sum,&
  mpi_comm_world,ierr)
  S_phi_q=a_recv_op

  a_send_op=S_phi_dj_eq
  call mpi_allreduce(a_send_op,a_recv_op,1,mpi_double_precision,mpi_sum,&
  mpi_comm_world,ierr)
  S_phi_dj_eq=a_recv_op

  a_send_op=S_phi_x_p
  call mpi_allreduce(a_send_op,a_recv_op,1,mpi_double_precision,mpi_sum,&
  mpi_comm_world,ierr)
  S_phi_x_p=a_recv_op

  a_send_op=S_phi_mu
  call mpi_allreduce(a_send_op,a_recv_op,1,mpi_double_precision,mpi_sum,&
  mpi_comm_world,ierr)
  S_phi_mu=a_recv_op

  a_send_op=S_phi_phi_U
  call mpi_allreduce(a_send_op,a_recv_op,1,mpi_double_precision,mpi_sum,&
  mpi_comm_world,ierr)
  S_phi_phi_U=a_recv_op

  a_send_op=S_phi_psi_j
  call mpi_allreduce(a_send_op,a_recv_op,1,mpi_double_precision,mpi_sum,&
  mpi_comm_world,ierr)
  S_phi_psi_j=a_recv_op

  L_phi=S_phi_dphi_eq+S_phi_q+S_phi_dj_eq+S_phi_x_p+S_phi_mu
  N_phi=S_phi_phi_U+S_phi_psi_j

  if(myid==0)then
    write(filename,'("../data/time_phi_energy_balance.dat")')
    if(it==it_start)then
      write(6,*)filename
      open(unit=10,file=filename,form='formatted')
      write(10,fmt='(4a16)')'t','dE_phi','L_phi','N_phi'
      dE_phi=0.0d0
      L_phi=0.0d0
      N_phi=0.0d0
    else
      open(unit=10,file=filename,form='formatted',position='append')
    end if
    write(10,fmt='(4e16.4)')t,dE_phi,L_phi,N_phi
    close(10)

    write(filename,'("../data/time_phi_energy_source_lin.dat")')
    if(it==it_start)then
      write(6,*)filename
      open(unit=10,file=filename,form='formatted')
      write(10,fmt='(6a16)')'t','S_phi_dphi_eq','S_phi_q','S_phi_dj_eq',&
      'S_phi_x_p','S_phi_mu'
      S_phi_dphi_eq=0.0d0
      S_phi_q=0.0d0
      S_phi_dj_eq=0.0d0
      S_phi_x_p=0.0d0
      S_phi_mu=0.0d0
    else
      open(unit=10,file=filename,form='formatted',position='append')
    end if
    write(10,fmt='(6e16.4)')t,S_phi_dphi_eq,S_phi_q,S_phi_dj_eq,&
    S_phi_x_p,S_phi_mu
    close(10)

    write(filename,'("../data/time_phi_energy_source_non.dat")')
    if(it==it_start)then
      write(6,*)filename
      open(unit=10,file=filename,form='formatted')
      write(10,fmt='(3a16)')'t','S_phi_phi_U','S_phi_psi_j'
      S_phi_phi_U=0.0d0
      S_phi_psi_j=0.0d0
    else
      open(unit=10,file=filename,form='formatted',position='append')
    end if
    write(10,fmt='(3e16.4)')t,S_phi_phi_U,S_phi_psi_j
    close(10)
  end if ! : myid==0

  ! . magnetic energy
  do n=sn,en
    dE_psi_n(n)=0.0d0
  end do

  do l=l_min,l_max
    do n=sn,en
      do irho=irho_min,irho_max
        psi(irho)=psi_nl(irho,n,l)
        j(irho)=j_nl(irho,n,l)
        psi_old(irho)=psi_old_nl(irho,n,l)
        j_old(irho)=j_old_nl(irho,n,l)
      end do

      call conjugate(j,j_cc)
      call conjugate(j_old,j_old_cc)

      do irho=irho_min,irho_max
        dE_psi_nl_den(irho)=(dble(j_cc(irho)*psi(irho))&
        -dble(j_old_cc(irho)*psi_old(irho)))/(2.0d0*dt)
      end do

      dE_psi_nl=0.0d0
      do irho=irho_min+1,irho_max
        dE_psi_nl=dE_psi_nl+drho/2.0d0&
        *(dE_psi_nl_den(irho-1)*rho(irho-1)&
        +dE_psi_nl_den(irho)*rho(irho))
      end do

      dE_psi_n(n)=dE_psi_n(n)+dE_psi_nl
    end do ! : n
  end do ! : l

  dE_psi=0.0d0
  do n=sn,en
    if(n==0)then
      dE_psi=dE_psi+dE_psi_n(n)
    else
      dE_psi=dE_psi+2.0d0*dE_psi_n(n)
    end if
  end do

  ! . linear energy sources
  do n=sn,en
    S_j_q_n(n)=0.0d0
    S_j_dphi_eq_n(n)=0.0d0
    S_j_eta_n(n)=0.0d0
  end do

  do l=l_min,l_max
    do n=sn,en
      do irho=irho_min,irho_max
        phi(irho)=phi_nl(irho,n,l)
        psi(irho)=psi_nl(irho,n,l)
        j(irho)=j_nl(irho,n,l)
      end do

      call conjugate(j,j_cc)

      do irho=irho_min,irho_max
        S_j_q_nl_den(irho)=-dble((0.0d0,1.0d0)*k_para(irho,n,l)&
        *j_cc(irho)*phi(irho))
        S_j_dphi_eq_nl_den(irho)=-dble((0.0d0,1.0d0)*k_theta(irho,n)&
        *dphi_eq(irho)*j_cc(irho)*psi(irho))
        S_j_eta_nl_den(irho)=-eta*abs(j(irho))**2
      end do

      S_j_q_nl=0.0d0
      S_j_dphi_eq_nl=0.0d0
      S_j_eta_nl=0.0d0
      do irho=irho_min+1,irho_max
        S_j_q_nl=S_j_q_nl+drho/2.0d0&
        *(S_j_q_nl_den(irho-1)*rho(irho-1)&
        +S_j_q_nl_den(irho)*rho(irho))
        S_j_dphi_eq_nl=S_j_dphi_eq_nl+drho/2.0d0&
        *(S_j_dphi_eq_nl_den(irho-1)*rho(irho-1)&
        +S_j_dphi_eq_nl_den(irho)*rho(irho))
        S_j_eta_nl=S_j_eta_nl+drho/2.0d0&
        *(S_j_eta_nl_den(irho-1)*rho(irho-1)&
        +S_j_eta_nl_den(irho)*rho(irho))
      end do

      S_j_q_n(n)=S_j_q_n(n)+S_j_q_nl
      S_j_dphi_eq_n(n)=S_j_dphi_eq_n(n)+S_j_dphi_eq_nl
      S_j_eta_n(n)=S_j_eta_n(n)+S_j_eta_nl
    end do ! : n
  end do ! : l

  S_j_q=0.0d0
  S_j_dphi_eq=0.0d0
  S_j_eta=0.0d0
  do n=sn,en
    if(n==0)then
      S_j_q=S_j_q+S_j_q_n(n)
      S_j_dphi_eq=S_j_dphi_eq+S_j_dphi_eq_n(n)
      S_j_eta=S_j_eta+S_j_eta_n(n)
    else
      S_j_q=S_j_q+2.0d0*S_j_q_n(n)
      S_j_dphi_eq=S_j_dphi_eq+2.0d0*S_j_dphi_eq_n(n)
      S_j_eta=S_j_eta+2.0d0*S_j_eta_n(n)
    end if
  end do

  ! . nonlinear energy sources
  do n=sn,en
    S_j_psi_phi_n(n)=0.0d0
  end do

  do l=l_min,l_max
    do n=sn,en
      do irho=irho_min,irho_max
        j(irho)=j_nl(irho,n,l)
        nab_para_phi_non(irho)=nab_para_phi_non_nl(irho,n,l)
      end do

      call conjugate(j,j_cc)

      do irho=irho_min,irho_max
        S_j_psi_phi_nl_den(irho)=-dble(j_cc(irho)*nab_para_phi_non(irho))
      end do

      S_j_psi_phi_nl=0.0d0
      do irho=irho_min+1,irho_max
        S_j_psi_phi_nl=S_j_psi_phi_nl+drho/2.0d0&
        *(S_j_psi_phi_nl_den(irho-1)*rho(irho-1)&
        +S_j_psi_phi_nl_den(irho)*rho(irho))
      end do

      S_j_psi_phi_n(n)=S_j_psi_phi_n(n)+S_j_psi_phi_nl
    end do ! : n
  end do ! : l

  S_j_psi_phi=0.0d0
  do n=sn,en
    if(n==0)then
      S_j_psi_phi=S_j_psi_phi+S_j_psi_phi_n(n)
    else
      S_j_psi_phi=S_j_psi_phi+2.0d0*S_j_psi_phi_n(n)
    end if
  end do

  a_send_op=dE_psi
  call mpi_allreduce(a_send_op,a_recv_op,1,mpi_double_precision,mpi_sum,&
  mpi_comm_world,ierr)
  dE_psi=a_recv_op

  a_send_op=S_j_q
  call mpi_allreduce(a_send_op,a_recv_op,1,mpi_double_precision,mpi_sum,&
  mpi_comm_world,ierr)
  S_j_q=a_recv_op

  a_send_op=S_j_dphi_eq
  call mpi_allreduce(a_send_op,a_recv_op,1,mpi_double_precision,mpi_sum,&
  mpi_comm_world,ierr)
  S_j_dphi_eq=a_recv_op

  a_send_op=S_j_eta
  call mpi_allreduce(a_send_op,a_recv_op,1,mpi_double_precision,mpi_sum,&
  mpi_comm_world,ierr)
  S_j_eta=a_recv_op

  a_send_op=S_j_psi_phi
  call mpi_allreduce(a_send_op,a_recv_op,1,mpi_double_precision,mpi_sum,&
  mpi_comm_world,ierr)
  S_j_psi_phi=a_recv_op

  L_psi=S_j_q+S_j_dphi_eq+S_j_eta
  N_psi=S_j_psi_phi

  if(myid==0)then
    write(filename,'("../data/time_psi_energy_balance.dat")')
    if(it==it_start)then
      write(6,*)filename
      open(unit=10,file=filename,form='formatted')
      write(10,fmt='(4a16)')'t','dE_psi','L_psi','N_psi'
      dE_psi=0.0d0
      L_psi=0.0d0
      N_psi=0.0d0
    else
      open(unit=10,file=filename,form='formatted',position='append')
    end if
    write(10,fmt='(4e16.4)')t,dE_psi,L_psi,N_psi
    close(10)

    write(filename,'("../data/time_psi_energy_source_lin.dat")')
    if(it==it_start)then
      write(6,*)filename
      open(unit=10,file=filename,form='formatted')
      write(10,fmt='(4a16)')'t','S_j_q','S_j_dphi_eq','S_j_eta'
      S_j_q=0.0d0
      S_j_dphi_eq=0.0d0
      S_j_eta=0.0d0
    else
      open(unit=10,file=filename,form='formatted',position='append')
    end if
    write(10,fmt='(4e16.4)')t,S_j_q,S_j_dphi_eq,S_j_eta
    close(10)

    write(filename,'("../data/time_psi_energy_source_non.dat")')
    if(it==it_start)then
      write(6,*)filename
      open(unit=10,file=filename,form='formatted')
      write(10,fmt='(2a16)')'t','S_j_psi_phi'
      S_j_psi_phi=0.0d0
    else
      open(unit=10,file=filename,form='formatted',position='append')
    end if
    write(10,fmt='(2e16.4)')t,S_j_psi_phi
    close(10)
  end if ! : myid==0

  ! . pressure energy
  do n=sn,en
    dE_p_n(n)=0.0d0
  end do

  do l=l_min,l_max
    do n=sn,en
      do irho=irho_min,irho_max
        p(irho)=p_nl(irho,n,l)
        p_old(irho)=p_old_nl(irho,n,l)
      end do

      do irho=irho_min,irho_max
        dE_p_nl_den(irho)=(abs(p(irho))**2&
        -abs(p_old(irho))**2)/(4.0d0*beta*dt)
      end do

      dE_p_nl=0.0d0
      do irho=irho_min+1,irho_max
        dE_p_nl=dE_p_nl+drho/2.0d0&
        *(dE_p_nl_den(irho-1)*rho(irho-1)&
        +dE_p_nl_den(irho)*rho(irho))
      end do

      dE_p_n(n)=dE_p_n(n)+dE_p_nl
    end do ! : n
  end do ! : l

  dE_p=0.0d0
  do n=sn,en
    if(n==0)then
      dE_p=dE_p+dE_p_n(n)
    else
      dE_p=dE_p+2.0d0*dE_p_n(n)
    end if
  end do

  ! . linear energy sources
  do n=sn,en
    S_p_dp_eq_n(n)=0.0d0
    S_p_x_phi_n(n)=0.0d0
    S_p_chi_n(n)=0.0d0
  end do

  do l=l_min,l_max
    do n=sn,en
      do irho=irho_min,irho_max
        phi(irho)=phi_nl(irho,n,l)
        p(irho)=p_nl(irho,n,l)
        bra_x_phi(irho)=bra_x_phi_nl(irho,n,l)
        lap_perp_p(irho)=lap_perp_p_nl(irho,n,l)
      end do

      call conjugate(p,p_cc)

      do irho=irho_min,irho_max
        S_p_dp_eq_nl_den(irho)=dble((0.0d0,1.0d0)*k_theta(irho,n)&
        *dp_eq(irho)*p_cc(irho)*phi(irho))/(2.0d0*beta)
        S_p_x_phi_nl_den(irho)=dble(p_cc(irho)*bra_x_phi(irho))
        S_p_chi_nl_den(irho)=chi/(2.0d0*beta)&
        *dble(p_cc(irho)*lap_perp_p(irho))
      end do

      S_p_dp_eq_nl=0.0d0
      S_p_x_phi_nl=0.0d0
      S_p_chi_nl=0.0d0
      do irho=irho_min+1,irho_max
        S_p_dp_eq_nl=S_p_dp_eq_nl+drho/2.0d0&
        *(S_p_dp_eq_nl_den(irho-1)*rho(irho-1)&
        +S_p_dp_eq_nl_den(irho)*rho(irho))
        S_p_x_phi_nl=S_p_x_phi_nl+drho/2.0d0&
        *(S_p_x_phi_nl_den(irho-1)*rho(irho-1)&
        +S_p_x_phi_nl_den(irho)*rho(irho))
        S_p_chi_nl=S_p_chi_nl+drho/2.0d0&
        *(S_p_chi_nl_den(irho-1)*rho(irho-1)&
        +S_p_chi_nl_den(irho)*rho(irho))
      end do

      S_p_dp_eq_n(n)=S_p_dp_eq_n(n)+S_p_dp_eq_nl
      S_p_x_phi_n(n)=S_p_x_phi_n(n)+S_p_x_phi_nl
      S_p_chi_n(n)=S_p_chi_n(n)+S_p_chi_nl
    end do ! : n
  end do ! : l

  S_p_dp_eq=0.0d0
  S_p_x_phi=0.0d0
  S_p_chi=0.0d0
  do n=sn,en
    if(n==0)then
      S_p_dp_eq=S_p_dp_eq+S_p_dp_eq_n(n)
      S_p_x_phi=S_p_x_phi+S_p_x_phi_n(n)
      S_p_chi=S_p_chi+S_p_chi_n(n)
    else
      S_p_dp_eq=S_p_dp_eq+2.0d0*S_p_dp_eq_n(n)
      S_p_x_phi=S_p_x_phi+2.0d0*S_p_x_phi_n(n)
      S_p_chi=S_p_chi+2.0d0*S_p_chi_n(n)
    end if
  end do

  ! . nonlinear energy sources
  do n=sn,en
    S_p_phi_p_n(n)=0.0d0
  end do

  do l=l_min,l_max
    do n=sn,en
      do irho=irho_min,irho_max
        p(irho)=p_nl(irho,n,l)
        bra_phi_p_non(irho)=bra_phi_p_non_nl(irho,n,l)
      end do

      call conjugate(p,p_cc)

      do irho=irho_min,irho_max
        S_p_phi_p_nl_den(irho)=-dble(p_cc(irho)*bra_phi_p_non(irho))/(2.0d0*beta)
      end do

      S_p_phi_p_nl=0.0d0
      do irho=irho_min+1,irho_max
        S_p_phi_p_nl=S_p_phi_p_nl+drho/2.0d0&
        *(S_p_phi_p_nl_den(irho-1)*rho(irho-1)&
        +S_p_phi_p_nl_den(irho)*rho(irho))
      end do

      S_p_phi_p_n(n)=S_p_phi_p_n(n)+S_p_phi_p_nl
    end do ! : n
  end do ! : l

  S_p_phi_p=0.0d0
  do n=sn,en
    if(n==0)then
      S_p_phi_p=S_p_phi_p+S_p_phi_p_n(n)
    else
      S_p_phi_p=S_p_phi_p+2.0d0*S_p_phi_p_n(n)
    end if
  end do

  a_send_op=dE_p
  call mpi_allreduce(a_send_op,a_recv_op,1,mpi_double_precision,mpi_sum,&
  mpi_comm_world,ierr)
  dE_p=a_recv_op

  a_send_op=S_p_dp_eq
  call mpi_allreduce(a_send_op,a_recv_op,1,mpi_double_precision,mpi_sum,&
  mpi_comm_world,ierr)
  S_p_dp_eq=a_recv_op

  a_send_op=S_p_x_phi
  call mpi_allreduce(a_send_op,a_recv_op,1,mpi_double_precision,mpi_sum,&
  mpi_comm_world,ierr)
  S_p_x_phi=a_recv_op

  a_send_op=S_p_chi
  call mpi_allreduce(a_send_op,a_recv_op,1,mpi_double_precision,mpi_sum,&
  mpi_comm_world,ierr)
  S_p_chi=a_recv_op

  a_send_op=S_p_phi_p
  call mpi_allreduce(a_send_op,a_recv_op,1,mpi_double_precision,mpi_sum,&
  mpi_comm_world,ierr)
  S_p_phi_p=a_recv_op

  L_p=S_p_dp_eq+S_p_x_phi+S_p_chi
  N_p=S_p_phi_p

  if(myid==0)then
    write(filename,'("../data/time_p_energy_balance.dat")')
    if(it==it_start)then
      write(6,*)filename
      open(unit=10,file=filename,form='formatted')
      write(10,fmt='(4a16)')'t','dE_p','L_p','N_p'
      dE_p=0.0d0
      L_p=0.0d0
      N_p=0.0d0
    else
      open(unit=10,file=filename,form='formatted',position='append')
    end if
    write(10,fmt='(4e16.4)')t,dE_p,L_p,N_p
    close(10)

    write(filename,'("../data/time_p_energy_source_lin.dat")')
    if(it==it_start)then
      write(6,*)filename
      open(unit=10,file=filename,form='formatted')
      write(10,fmt='(4a16)')'t','S_p_dp_eq','S_p_x_phi','S_p_chi'
      S_p_dp_eq=0.0d0
      S_p_x_phi=0.0d0
      S_p_chi=0.0d0
    else
      open(unit=10,file=filename,form='formatted',position='append')
    end if
    write(10,fmt='(4e16.4)')t,S_p_dp_eq,S_p_x_phi,S_p_chi
    close(10)

    write(filename,'("../data/time_p_energy_source_non.dat")')
    if(it==it_start)then
      write(6,*)filename
      open(unit=10,file=filename,form='formatted')
      write(10,fmt='(2a16)')'t','S_p_phi_p'
      S_p_phi_p=0.0d0
    else
      open(unit=10,file=filename,form='formatted',position='append')
    end if
    write(10,fmt='(2e16.4)')t,S_p_phi_p
    close(10)
  end if ! : myid==0

  ! . total energy
  dE_tot=dE_phi+dE_psi+dE_p
  L_tot=L_phi+L_psi+L_p
  N_tot=N_phi+N_psi+N_p

  if(myid==0)then
    write(filename,'("../data/time_tot_energy_balance.dat")')
    if(it==it_start)then
      write(6,*)filename
      open(unit=10,file=filename,form='formatted')
      write(10,fmt='(4a16)')'t','dE_tot','L_tot','N_tot'
    else
      open(unit=10,file=filename,form='formatted',position='append')
    end if
    write(10,fmt='(4e16.4)')t,dE_tot,L_tot,N_tot
    close(10)
  end if

end subroutine energy_balances
! ------------------------------------------------------------------------------
subroutine energies
! ------------------------------------------------------------------------------
  use time
  use coordinate
  use parameter
  use field
  use mpi_global
  use mpi
  implicit none

  ! . kinetic energy
  do n=sn,en
    E_phi_n(n)=0.0d0
  end do

  do l=l_min,l_max
    do n=sn,en
      do irho=irho_min,irho_max
        U(irho)=U_nl(irho,n,l)
        phi(irho)=phi_nl(irho,n,l)
      end do

      call conjugate(phi,phi_cc)

      do irho=irho_min,irho_max
        E_phi_nl_den(irho)=-dble(phi_cc(irho)*U(irho))/2.0d0
      end do

      E_phi_nl=0.0d0
      do irho=irho_min+1,irho_max
        E_phi_nl=E_phi_nl+drho/2.0d0&
        *(E_phi_nl_den(irho-1)*rho(irho-1)&
        +E_phi_nl_den(irho)*rho(irho))
      end do

      E_phi_n(n)=E_phi_n(n)+E_phi_nl
    end do ! : n
  end do ! : l

  E_phi=0.0d0
  do n=sn,en
    if(n==0)then
      E_phi=E_phi+E_phi_n(n)
    else
      E_phi=E_phi+2.0d0*E_phi_n(n)
    end if
  end do

  ! . magnetic energy
  do n=sn,en
    E_psi_n(n)=0.0d0
  end do

  do l=l_min,l_max
    do n=sn,en
      do irho=irho_min,irho_max
        psi(irho)=psi_nl(irho,n,l)
        j(irho)=j_nl(irho,n,l)
      end do

      call conjugate(j,j_cc)

      do irho=irho_min,irho_max
        E_psi_nl_den(irho)=dble(j_cc(irho)*psi(irho))/2.0d0
      end do

      E_psi_nl=0.0d0
      do irho=irho_min+1,irho_max
        E_psi_nl=E_psi_nl+drho/2.0d0&
        *(E_psi_nl_den(irho-1)*rho(irho-1)&
        +E_psi_nl_den(irho)*rho(irho))
      end do

      E_psi_n(n)=E_psi_n(n)+E_psi_nl
    end do ! : n
  end do ! : l

  E_psi=0.0d0
  do n=sn,en
    if(n==0)then
      E_psi=E_psi+E_psi_n(n)
    else
      E_psi=E_psi+2.0d0*E_psi_n(n)
    end if
  end do

  ! . pressure energy
  do n=sn,en
    E_p_n(n)=0.0d0
  end do

  do l=l_min,l_max
    do n=sn,en
      do irho=irho_min,irho_max
        p(irho)=p_nl(irho,n,l)
      end do

      do irho=irho_min,irho_max
        E_p_nl_den(irho)=abs(p(irho))**2/(4.0d0*beta)
      end do

      E_p_nl=0.0d0
      do irho=irho_min+1,irho_max
        E_p_nl=E_p_nl+drho/2.0d0&
        *(E_p_nl_den(irho-1)*rho(irho-1)&
        +E_p_nl_den(irho)*rho(irho))
      end do

      E_p_n(n)=E_p_n(n)+E_p_nl
    end do ! : n
  end do ! : l

  E_p=0.0d0
  do n=sn,en
    if(n==0)then
      E_p=E_p+E_p_n(n)
    else
      E_p=E_p+2.0d0*E_p_n(n)
    end if
  end do

  a_send_op=E_phi
  call mpi_allreduce(a_send_op,a_recv_op,1,mpi_double_precision,mpi_sum,&
  mpi_comm_world,ierr)
  E_phi=a_recv_op

  a_send_op=E_psi
  call mpi_allreduce(a_send_op,a_recv_op,1,mpi_double_precision,mpi_sum,&
  mpi_comm_world,ierr)
  E_psi=a_recv_op

  a_send_op=E_p
  call mpi_allreduce(a_send_op,a_recv_op,1,mpi_double_precision,mpi_sum,&
  mpi_comm_world,ierr)
  E_p=a_recv_op

  if(myid==0)then
    write(filename,'("../data/time_energy.dat")')
    if(it==it_start)then
      write(6,*)filename
      open(unit=10,file=filename,form='formatted')
      write(10,fmt='(4a16)')'t','E_phi','E_psi','E_p'
    else
      open(unit=10,file=filename,form='formatted',position='append')
    end if
    write(10,fmt='(4e16.4)')t,E_phi,E_psi,E_p
    close(10)
  end if

  do n=sn,en
    a_send_gather(n)=E_phi_n(n)
  end do
  call mpi_allgatherv(a_send_gather,scount_gather,mpi_double_complex,&
  a_recv_gather,rcounts_gather,rdispls_gather,mpi_double_complex,&
  mpi_comm_world,ierr)
  do n=0,n_max
    E_phi_n_gather(n)=a_recv_gather(n)
  end do

  do n=sn,en
    a_send_gather(n)=E_psi_n(n)
  end do
  call mpi_allgatherv(a_send_gather,scount_gather,mpi_double_complex,&
  a_recv_gather,rcounts_gather,rdispls_gather,mpi_double_complex,&
  mpi_comm_world,ierr)
  do n=0,n_max
    E_psi_n_gather(n)=a_recv_gather(n)
  end do

  do n=sn,en
    a_send_gather(n)=E_p_n(n)
  end do
  call mpi_allgatherv(a_send_gather,scount_gather,mpi_double_complex,&
  a_recv_gather,rcounts_gather,rdispls_gather,mpi_double_complex,&
  mpi_comm_world,ierr)
  do n=0,n_max
    E_p_n_gather(n)=a_recv_gather(n)
  end do

  if(myid==0)then
    write(filename,'("../data/time_phi_energy_spec.dat")')
    if(it==it_start)then
      write(6,*)filename
      open(unit=10,file=filename,form='formatted')
      write(10,fmt='(2a16)')'t','E_phi_n'
    else
      open(unit=10,file=filename,form='formatted',position='append')
    end if
    n=0
    write(10,fmt='(2e16.4)',advance='no')t,E_phi_n_gather(n)
    do n=1,n_max
      write(10,fmt='(e16.4)',advance='no')2.0d0*E_phi_n_gather(n)
    end do
    write(10,*)
    close(10)

    write(filename,'("../data/time_psi_energy_spec.dat")')
    if(it==it_start)then
      write(6,*)filename
      open(unit=10,file=filename,form='formatted')
      write(10,fmt='(2a16)')'t','E_psi_n'
    else
      open(unit=10,file=filename,form='formatted',position='append')
    end if
    n=0
    write(10,fmt='(2e16.4)',advance='no')t,E_psi_n_gather(n)
    do n=1,n_max
      write(10,fmt='(e16.4)',advance='no')2.0d0*E_psi_n_gather(n)
    end do
    write(10,*)
    close(10)

    write(filename,'("../data/time_p_energy_spec.dat")')
    if(it==it_start)then
      write(6,*)filename
      open(unit=10,file=filename,form='formatted')
      write(10,fmt='(2a16)')'t','E_p_n'
    else
      open(unit=10,file=filename,form='formatted',position='append')
    end if
    n=0
    write(10,fmt='(2e16.4)',advance='no')t,E_p_n_gather(n)
    do n=1,n_max
      write(10,fmt='(e16.4)',advance='no')2.0d0*E_p_n_gather(n)
    end do
    write(10,*)
    close(10)
  end if ! : myid==0

  if(mod(it,it_skip)==0)then
    if(myid==0)then
      write(filename,'("../data/energy_spec_",i3.3,".dat")')it/it_skip
      write(6,*)filename
      open(unit=10,file=filename,form='formatted')
      write(10,fmt='(4a16)')'n','E_phi_n','E_psi_n','E_p_n'
      n=0
      write(10,fmt='(4e16.4)')dble(n),E_phi_n_gather(n),&
      E_psi_n_gather(n),E_p_n_gather(n)
      do n=1,n_max
        write(10,fmt='(4e16.4)')dble(n),2.0d0*E_phi_n_gather(n),&
        2.0d0*E_psi_n_gather(n),2.0d0*E_p_n_gather(n)
      end do
      close(10)
    end if
  end if

end subroutine energies
! ------------------------------------------------------------------------------
subroutine energy_spectra
! ------------------------------------------------------------------------------
  use time
  use coordinate
  use parameter
  use field
  use mpi_global
  use mpi
  implicit none

  if(it==it_start)then
    if(irestart==0)then
      do n=sn,en
        int_E_phi_n(n)=0.0d0
        int_E_psi_n(n)=0.0d0
        int_E_p_n(n)=0.0d0
      end do
    else
      write(filename,'("../data/restart_energy_",i3.3,".dat")')myid
      open(unit=10*(1+myid),file=filename,form='unformatted')
      do n=sn,en
        read(10*(1+myid))int_E_phi_n(n),int_E_psi_n(n),int_E_p_n(n)
      end do
      close(10*(1+myid))
    end if
  else
    ! . kinetic energy
    do n=sn,en
      E_phi_n(n)=0.0d0
      E_phi_old_n(n)=0.0d0
    end do

    do l=l_min,l_max
      do n=sn,en
        do irho=irho_min,irho_max
          U(irho)=U_nl(irho,n,l)
          phi(irho)=phi_nl(irho,n,l)
          U_old(irho)=U_old_nl(irho,n,l)
          phi_old(irho)=phi_old_nl(irho,n,l)
        end do

        call conjugate(phi,phi_cc)
        call conjugate(phi_old,phi_old_cc)

        do irho=irho_min,irho_max
          E_phi_nl_den(irho)=-dble(phi_cc(irho)*U(irho))/2.0d0
          E_phi_old_nl_den(irho)=-dble(phi_old_cc(irho)*U_old(irho))/2.0d0
        end do

        E_phi_nl=0.0d0
        E_phi_old_nl=0.0d0
        do irho=irho_min+1,irho_max
          E_phi_nl=E_phi_nl+drho/2.0d0&
          *(E_phi_nl_den(irho-1)*rho(irho-1)&
          +E_phi_nl_den(irho)*rho(irho))
          E_phi_old_nl=E_phi_old_nl+drho/2.0d0&
          *(E_phi_old_nl_den(irho-1)*rho(irho-1)&
          +E_phi_old_nl_den(irho)*rho(irho))
        end do

        E_phi_n(n)=E_phi_n(n)+E_phi_nl
        E_phi_old_n(n)=E_phi_old_n(n)+E_phi_old_nl
      end do ! : n
    end do ! : l

    do n=sn,en
      int_E_phi_n(n)=int_E_phi_n(n)&
      +dt/2.0d0*(E_phi_n(n)+E_phi_old_n(n))
    end do

    ! . magnetic energy
    do n=sn,en
      E_psi_n(n)=0.0d0
      E_psi_old_n(n)=0.0d0
    end do

    do l=l_min,l_max
      do n=sn,en
        do irho=irho_min,irho_max
          psi(irho)=psi_nl(irho,n,l)
          j(irho)=j_nl(irho,n,l)
          psi_old(irho)=psi_old_nl(irho,n,l)
          j_old(irho)=j_old_nl(irho,n,l)
        end do

        call conjugate(j,j_cc)
        call conjugate(j_old,j_old_cc)

        do irho=irho_min,irho_max
          E_psi_nl_den(irho)=dble(j_cc(irho)*psi(irho))/2.0d0
          E_psi_old_nl_den(irho)=dble(j_old_cc(irho)*psi_old(irho))/2.0d0
        end do

        E_psi_nl=0.0d0
        E_psi_old_nl=0.0d0
        do irho=irho_min+1,irho_max
          E_psi_nl=E_psi_nl+drho/2.0d0&
          *(E_psi_nl_den(irho-1)*rho(irho-1)&
          +E_psi_nl_den(irho)*rho(irho))
          E_psi_old_nl=E_psi_old_nl+drho/2.0d0&
          *(E_psi_old_nl_den(irho-1)*rho(irho-1)&
          +E_psi_old_nl_den(irho)*rho(irho))
        end do

        E_psi_n(n)=E_psi_n(n)+E_psi_nl
        E_psi_old_n(n)=E_psi_old_n(n)+E_psi_old_nl
      end do ! : n
    end do ! : l

    do n=sn,en
      int_E_psi_n(n)=int_E_psi_n(n)&
      +dt/2.0d0*(E_psi_n(n)+E_psi_old_n(n))
    end do

    ! . pressure energy
    do n=sn,en
      E_p_n(n)=0.0d0
      E_p_old_n(n)=0.0d0
    end do

    do l=l_min,l_max
      do n=sn,en
        do irho=irho_min,irho_max
          p(irho)=p_nl(irho,n,l)
          p_old(irho)=p_old_nl(irho,n,l)
        end do

        do irho=irho_min,irho_max
          E_p_nl_den(irho)=abs(p(irho))**2/(4.0d0*beta)
          E_p_old_nl_den(irho)=abs(p_old(irho))**2/(4.0d0*beta)
        end do

        E_p_nl=0.0d0
        E_p_old_nl=0.0d0
        do irho=irho_min+1,irho_max
          E_p_nl=E_p_nl+drho/2.0d0&
          *(E_p_nl_den(irho-1)*rho(irho-1)&
          +E_p_nl_den(irho)*rho(irho))
          E_p_old_nl=E_p_old_nl+drho/2.0d0&
          *(E_p_old_nl_den(irho-1)*rho(irho-1)&
          +E_p_old_nl_den(irho)*rho(irho))
        end do

        E_p_n(n)=E_p_n(n)+E_p_nl
        E_p_old_n(n)=E_p_old_n(n)+E_p_old_nl
      end do ! : n
    end do ! : l

    do n=sn,en
      int_E_p_n(n)=int_E_p_n(n)&
      +dt/2.0d0*(E_p_n(n)+E_p_old_n(n))
    end do

    if(mod(it,it_skip)==0)then
      do n=sn,en
        a_send_gather(n)=int_E_phi_n(n)
      end do
      call mpi_allgatherv(a_send_gather,scount_gather,mpi_double_complex,&
      a_recv_gather,rcounts_gather,rdispls_gather,mpi_double_complex,&
      mpi_comm_world,ierr)
      do n=0,n_max
        int_E_phi_n_gather(n)=a_recv_gather(n)
      end do

      do n=sn,en
        a_send_gather(n)=int_E_psi_n(n)
      end do
      call mpi_allgatherv(a_send_gather,scount_gather,mpi_double_complex,&
      a_recv_gather,rcounts_gather,rdispls_gather,mpi_double_complex,&
      mpi_comm_world,ierr)
      do n=0,n_max
        int_E_psi_n_gather(n)=a_recv_gather(n)
      end do

      do n=sn,en
        a_send_gather(n)=int_E_p_n(n)
      end do
      call mpi_allgatherv(a_send_gather,scount_gather,mpi_double_complex,&
      a_recv_gather,rcounts_gather,rdispls_gather,mpi_double_complex,&
      mpi_comm_world,ierr)
      do n=0,n_max
        int_E_p_n_gather(n)=a_recv_gather(n)
      end do

      if(myid==0)then
        write(filename,'("../data/int_energy_spec_",i3.3,".dat")')it/it_skip
        write(6,*)filename
        open(unit=10,file=filename,form='formatted')
        write(10,fmt='(4a16)')'n','int_E_phi_n','int_E_psi_n','int_E_p_n'
        n=0
        write(10,fmt='(4e16.4)')dble(n),int_E_phi_n_gather(n),&
        int_E_psi_n_gather(n),int_E_p_n_gather(n)
        do n=1,n_max
          write(10,fmt='(4e16.4)')dble(n),2.0d0*int_E_phi_n_gather(n),&
          2.0d0*int_E_psi_n_gather(n),2.0d0*int_E_p_n_gather(n)
        end do
        close(10)
      end if
    end if ! : mod(it,it_skip)==0
  end if ! : it==it_start

  if(it==it_end)then
    write(filename,'("../data/restart_energy_",i3.3,".dat")')myid
    open(unit=10*(1+myid),file=filename,form='unformatted')
    do n=sn,en
      write(10*(1+myid))int_E_phi_n(n),int_E_psi_n(n),int_E_p_n(n)
    end do
    close(10*(1+myid))
  end if

end subroutine energy_spectra
! ------------------------------------------------------------------------------
subroutine helicity_balance
! ------------------------------------------------------------------------------
  use time
  use coordinate
  use parameter
  use field
  use mpi_global
  use mpi
  implicit none

  do n=sn,en
    dH_U_psi_n(n)=0.0d0
  end do

  do l=l_min,l_max
    do n=sn,en
      do irho=irho_min,irho_max
        U(irho)=U_nl(irho,n,l)
        psi(irho)=psi_nl(irho,n,l)
        U_old(irho)=U_old_nl(irho,n,l)
        psi_old(irho)=psi_old_nl(irho,n,l)
      end do

      call conjugate(U,U_cc)
      call conjugate(U_old,U_old_cc)

      do irho=irho_min,irho_max
        dH_U_psi_nl_den(irho)=(dble(U_cc(irho)*psi(irho))&
        -dble(U_old_cc(irho)*psi_old(irho)))/dt
      end do

      dH_U_psi_nl=0.0d0
      do irho=irho_min+1,irho_max
        dH_U_psi_nl=dH_U_psi_nl+drho/2.0d0&
        *(dH_U_psi_nl_den(irho-1)*rho(irho-1)&
        +dH_U_psi_nl_den(irho)*rho(irho))
      end do

      dH_U_psi_n(n)=dH_U_psi_n(n)+dH_U_psi_nl
    end do ! : n
  end do ! : l

  dH_U_psi=0.0d0
  do n=sn,en
    if(n==0)then
      dH_U_psi=dH_U_psi+dH_U_psi_n(n)
    else
      dH_U_psi=dH_U_psi+2.0d0*dH_U_psi_n(n)
    end if
  end do

  ! . linear helicity sources
  do n=sn,en
    S_psi_dU_eq_n(n)=0.0d0
    S_U_psi_q_n(n)=0.0d0
    S_psi_x_p_n(n)=0.0d0
    S_psi_mu_n(n)=0.0d0
    S_U_eta_n(n)=0.0d0
  end do

  do l=l_min,l_max
    do n=sn,en
      do irho=irho_min,irho_max
        U(irho)=U_nl(irho,n,l)
        phi(irho)=phi_nl(irho,n,l)
        psi(irho)=psi_nl(irho,n,l)
        j(irho)=j_nl(irho,n,l)
        bra_x_p(irho)=bra_x_p_nl(irho,n,l)
        lap_perp_U(irho)=lap_perp_U_nl(irho,n,l)
      end do

      call conjugate(U,U_cc)
      call conjugate(psi,psi_cc)

      do irho=irho_min,irho_max
        S_psi_dU_eq_nl_den(irho)=dble((0.0d0,1.0d0)*k_theta(irho,n)&
        *dU_eq(irho)*psi_cc(irho)*phi(irho))
        S_U_psi_q_nl_den(irho)=dble((0.0d0,1.0d0)*k_para(irho,n,l)&
        *(psi_cc(irho)*j(irho)-U_cc(irho)*phi(irho)))
        S_psi_x_p_nl_den(irho)=-dble(psi_cc(irho)*bra_x_p(irho))
        S_psi_mu_nl_den(irho)=mu*dble(psi_cc(irho)*lap_perp_U(irho))
        S_U_eta_nl_den(irho)=-eta*dble(U_cc(irho)*j(irho))
      end do

      S_psi_dU_eq_nl=0.0d0
      S_U_psi_q_nl=0.0d0
      S_psi_x_p_nl=0.0d0
      S_psi_mu_nl=0.0d0
      S_U_eta_nl=0.0d0
      do irho=irho_min+1,irho_max
        S_psi_dU_eq_nl=S_psi_dU_eq_nl+drho/2.0d0&
        *(S_psi_dU_eq_nl_den(irho-1)*rho(irho-1)&
        +S_psi_dU_eq_nl_den(irho)*rho(irho))
        S_U_psi_q_nl=S_U_psi_q_nl+drho/2.0d0&
        *(S_U_psi_q_nl_den(irho-1)*rho(irho-1)&
        +S_U_psi_q_nl_den(irho)*rho(irho))
        S_psi_x_p_nl=S_psi_x_p_nl+drho/2.0d0&
        *(S_psi_x_p_nl_den(irho-1)*rho(irho-1)&
        +S_psi_x_p_nl_den(irho)*rho(irho))
        S_psi_mu_nl=S_psi_mu_nl+drho/2.0d0&
        *(S_psi_mu_nl_den(irho-1)*rho(irho-1)&
        +S_psi_mu_nl_den(irho)*rho(irho))
        S_U_eta_nl=S_U_eta_nl+drho/2.0d0&
        *(S_U_eta_nl_den(irho-1)*rho(irho-1)&
        +S_U_eta_nl_den(irho)*rho(irho))
      end do

      S_psi_dU_eq_n(n)=S_psi_dU_eq_n(n)+S_psi_dU_eq_nl
      S_U_psi_q_n(n)=S_U_psi_q_n(n)+S_U_psi_q_nl
      S_psi_x_p_n(n)=S_psi_x_p_n(n)+S_psi_x_p_nl
      S_psi_mu_n(n)=S_psi_mu_n(n)+S_psi_mu_nl
      S_U_eta_n(n)=S_U_eta_n(n)+S_U_eta_nl
    end do ! : n
  end do ! : l

  S_psi_dU_eq=0.0d0
  S_U_psi_q=0.0d0
  S_psi_x_p=0.0d0
  S_psi_mu=0.0d0
  S_U_eta=0.0d0
  do n=sn,en
    if(n==0)then
      S_psi_dU_eq=S_psi_dU_eq+S_psi_dU_eq_n(n)
      S_U_psi_q=S_U_psi_q+S_U_psi_q_n(n)
      S_psi_x_p=S_psi_x_p+S_psi_x_p_n(n)
      S_psi_mu=S_psi_mu+S_psi_mu_n(n)
      S_U_eta=S_U_eta+S_U_eta_n(n)
    else
      S_psi_dU_eq=S_psi_dU_eq+2.0d0*S_psi_dU_eq_n(n)
      S_U_psi_q=S_U_psi_q+2.0d0*S_U_psi_q_n(n)
      S_psi_x_p=S_psi_x_p+2.0d0*S_psi_x_p_n(n)
      S_psi_mu=S_psi_mu+2.0d0*S_psi_mu_n(n)
      S_U_eta=S_U_eta+2.0d0*S_U_eta_n(n)
    end if
  end do

  ! . nonlinear helicity sources
  do n=sn,en
    S_psi_phi_U_n(n)=0.0d0
    S_psi_psi_j_n(n)=0.0d0
    S_U_psi_phi_n(n)=0.0d0
  end do

  do l=l_min,l_max
    do n=sn,en
      do irho=irho_min,irho_max
        U(irho)=U_nl(irho,n,l)
        psi(irho)=psi_nl(irho,n,l)
        bra_phi_U_non(irho)=bra_phi_U_non_nl(irho,n,l)
        nab_para_j_non(irho)=nab_para_j_non_nl(irho,n,l)
        nab_para_phi_non(irho)=nab_para_phi_non_nl(irho,n,l)
      end do

      call conjugate(U,U_cc)
      call conjugate(psi,psi_cc)

      do irho=irho_min,irho_max
        S_psi_phi_U_nl_den(irho)=-dble(psi_cc(irho)*bra_phi_U_non(irho))
        S_psi_psi_j_nl_den(irho)=dble(psi_cc(irho)*nab_para_j_non(irho))
        S_U_psi_phi_nl_den(irho)=-dble(U_cc(irho)*nab_para_phi_non(irho))
      end do

      S_psi_phi_U_nl=0.0d0
      S_psi_psi_j_nl=0.0d0
      S_U_psi_phi_nl=0.0d0
      do irho=irho_min+1,irho_max
        S_psi_phi_U_nl=S_psi_phi_U_nl+drho/2.0d0&
        *(S_psi_phi_U_nl_den(irho-1)*rho(irho-1)&
        +S_psi_phi_U_nl_den(irho)*rho(irho))
        S_psi_psi_j_nl=S_psi_psi_j_nl+drho/2.0d0&
        *(S_psi_psi_j_nl_den(irho-1)*rho(irho-1)&
        +S_psi_psi_j_nl_den(irho)*rho(irho))
        S_U_psi_phi_nl=S_U_psi_phi_nl+drho/2.0d0&
        *(S_U_psi_phi_nl_den(irho-1)*rho(irho-1)&
        +S_U_psi_phi_nl_den(irho)*rho(irho))
      end do

      S_psi_phi_U_n(n)=S_psi_phi_U_n(n)+S_psi_phi_U_nl
      S_psi_psi_j_n(n)=S_psi_psi_j_n(n)+S_psi_psi_j_nl
      S_U_psi_phi_n(n)=S_U_psi_phi_n(n)+S_U_psi_phi_nl
    end do ! : n
  end do ! : l

  S_psi_phi_U=0.0d0
  S_psi_psi_j=0.0d0
  S_U_psi_phi=0.0d0
  do n=sn,en
    if(n==0)then
      S_psi_phi_U=S_psi_phi_U+S_psi_phi_U_n(n)
      S_psi_psi_j=S_psi_psi_j+S_psi_psi_j_n(n)
      S_U_psi_phi=S_U_psi_phi+S_U_psi_phi_n(n)
    else
      S_psi_phi_U=S_psi_phi_U+2.0d0*S_psi_phi_U_n(n)
      S_psi_psi_j=S_psi_psi_j+2.0d0*S_psi_psi_j_n(n)
      S_U_psi_phi=S_U_psi_phi+2.0d0*S_U_psi_phi_n(n)
    end if
  end do

  a_send_op=dH_U_psi
  call mpi_allreduce(a_send_op,a_recv_op,1,mpi_double_precision,mpi_sum,&
  mpi_comm_world,ierr)
  dH_U_psi=a_recv_op

  a_send_op=S_psi_dU_eq
  call mpi_allreduce(a_send_op,a_recv_op,1,mpi_double_precision,mpi_sum,&
  mpi_comm_world,ierr)
  S_psi_dU_eq=a_recv_op

  a_send_op=S_U_psi_q
  call mpi_allreduce(a_send_op,a_recv_op,1,mpi_double_precision,mpi_sum,&
  mpi_comm_world,ierr)
  S_U_psi_q=a_recv_op

  a_send_op=S_psi_x_p
  call mpi_allreduce(a_send_op,a_recv_op,1,mpi_double_precision,mpi_sum,&
  mpi_comm_world,ierr)
  S_psi_x_p=a_recv_op

  a_send_op=S_psi_mu
  call mpi_allreduce(a_send_op,a_recv_op,1,mpi_double_precision,mpi_sum,&
  mpi_comm_world,ierr)
  S_psi_mu=a_recv_op

  a_send_op=S_U_eta
  call mpi_allreduce(a_send_op,a_recv_op,1,mpi_double_precision,mpi_sum,&
  mpi_comm_world,ierr)
  S_U_eta=a_recv_op

  a_send_op=S_psi_phi_U
  call mpi_allreduce(a_send_op,a_recv_op,1,mpi_double_precision,mpi_sum,&
  mpi_comm_world,ierr)
  S_psi_phi_U=a_recv_op

  a_send_op=S_psi_psi_j
  call mpi_allreduce(a_send_op,a_recv_op,1,mpi_double_precision,mpi_sum,&
  mpi_comm_world,ierr)
  S_psi_psi_j=a_recv_op

  a_send_op=S_U_psi_phi
  call mpi_allreduce(a_send_op,a_recv_op,1,mpi_double_precision,mpi_sum,&
  mpi_comm_world,ierr)
  S_U_psi_phi=a_recv_op

  L_U_psi=S_psi_dU_eq+S_U_psi_q+S_psi_x_p+S_psi_mu+S_U_eta
  N_U_psi=S_psi_phi_U+S_psi_psi_j+S_U_psi_phi

  if(myid==0)then
    write(filename,'("../data/time_helicity_balance.dat")')
    if(it==it_start)then
      write(6,*)filename
      open(unit=10,file=filename,form='formatted')
      write(10,fmt='(4a16)')'t','dH_U_psi','L_U_psi','N_U_psi'
      dH_U_psi=0.0d0
      L_U_psi=0.0d0
      N_U_psi=0.0d0
    else
      open(unit=10,file=filename,form='formatted',position='append')
    end if
    write(10,fmt='(4e16.4)')t,dH_U_psi,L_U_psi,N_U_psi
    close(10)

    write(filename,'("../data/time_helicity_source_lin.dat")')
    if(it==it_start)then
      write(6,*)filename
      open(unit=10,file=filename,form='formatted')
      write(10,fmt='(6a16)')'t','S_psi_dU_eq','S_U_psi_q','S_psi_x_p',&
      'S_psi_mu','S_U_eta'
      S_psi_dU_eq=0.0d0
      S_U_psi_q=0.0d0
      S_psi_x_p=0.0d0
      S_psi_mu=0.0d0
      S_U_eta=0.0d0
    else
      open(unit=10,file=filename,form='formatted',position='append')
    end if
    write(10,fmt='(6e16.4)')t,S_psi_dU_eq,S_U_psi_q,S_psi_x_p,&
    S_psi_mu,S_U_eta
    close(10)

    write(filename,'("../data/time_helicity_source_non.dat")')
    if(it==it_start)then
      write(6,*)filename
      open(unit=10,file=filename,form='formatted')
      write(10,fmt='(4a16)')'t','S_psi_phi_U','S_psi_psi_j','S_U_psi_phi'
      S_psi_phi_U=0.0d0
      S_psi_psi_j=0.0d0
      S_U_psi_phi=0.0d0
    else
      open(unit=10,file=filename,form='formatted',position='append')
    end if
    write(10,fmt='(4e16.4)')t,S_psi_phi_U,S_psi_psi_j,S_U_psi_phi
    close(10)
  end if ! : myid==0

end subroutine helicity_balance
! ------------------------------------------------------------------------------
subroutine helicity
! ------------------------------------------------------------------------------
  use time
  use coordinate
  use parameter
  use field
  use mpi_global
  use mpi
  implicit none

  do n=sn,en
    H_U_psi_n(n)=0.0d0
  end do

  do l=l_min,l_max
    do n=sn,en
      do irho=irho_min,irho_max
        U(irho)=U_nl(irho,n,l)
        psi(irho)=psi_nl(irho,n,l)
      end do

      call conjugate(U,U_cc)

      do irho=irho_min,irho_max
        H_U_psi_nl_den(irho)=dble(U_cc(irho)*psi(irho))
      end do

      H_U_psi_nl=0.0d0
      do irho=irho_min+1,irho_max
        H_U_psi_nl=H_U_psi_nl+drho/2.0d0&
        *(H_U_psi_nl_den(irho-1)*rho(irho-1)&
        +H_U_psi_nl_den(irho)*rho(irho))
      end do

      H_U_psi_n(n)=H_U_psi_n(n)+H_U_psi_nl
    end do ! : n
  end do ! : l

  H_U_psi=0.0d0
  do n=sn,en
    if(n==0)then
      H_U_psi=H_U_psi+H_U_psi_n(n)
    else
      H_U_psi=H_U_psi+2.0d0*H_U_psi_n(n)
    end if
  end do

  a_send_op=H_U_psi
  call mpi_allreduce(a_send_op,a_recv_op,1,mpi_double_precision,mpi_sum,&
  mpi_comm_world,ierr)
  H_U_psi=a_recv_op

  if(myid==0)then
    write(filename,'("../data/time_helicity.dat")')
    if(it==it_start)then
      write(6,*)filename
      open(unit=10,file=filename,form='formatted')
      write(10,fmt='(2a16)')'t','H_U_psi'
    else
      open(unit=10,file=filename,form='formatted',position='append')
    end if
    write(10,fmt='(2e16.4)')t,H_U_psi
    close(10)
  end if

  do n=sn,en
    a_send_gather(n)=H_U_psi_n(n)
  end do
  call mpi_allgatherv(a_send_gather,scount_gather,mpi_double_complex,&
  a_recv_gather,rcounts_gather,rdispls_gather,mpi_double_complex,&
  mpi_comm_world,ierr)
  do n=0,n_max
    H_U_psi_n_gather(n)=a_recv_gather(n)
  end do

  if(myid==0)then
    write(filename,'("../data/time_helicity_spec.dat")')
    if(it==it_start)then
      write(6,*)filename
      open(unit=10,file=filename,form='formatted')
      write(10,fmt='(2a16)')'t','H_U_psi_n'
    else
      open(unit=10,file=filename,form='formatted',position='append')
    end if
    n=0
    write(10,fmt='(2e16.4)',advance='no')t,H_U_psi_n_gather(n)
    do n=1,n_max
      write(10,fmt='(e16.4)',advance='no')2.0d0*H_U_psi_n_gather(n)
    end do
    write(10,*)
    close(10)
  end if

  if(mod(it,it_skip)==0)then
    if(myid==0)then
      write(filename,'("../data/helicity_spec_",i3.3,".dat")')it/it_skip
      write(6,*)filename
      open(unit=10,file=filename,form='formatted')
      write(10,fmt='(2a16)')'n','H_U_psi_n'
      n=0
      write(10,fmt='(2e16.4)')dble(n),H_U_psi_n_gather(n)
      do n=1,n_max
        write(10,fmt='(2e16.4)')dble(n),2.0d0*H_U_psi_n_gather(n)
      end do
      close(10)
    end if
  end if

end subroutine helicity
! ------------------------------------------------------------------------------
subroutine helicity_spectrum
! ------------------------------------------------------------------------------
  use time
  use coordinate
  use parameter
  use field
  use mpi_global
  use mpi
  implicit none

  if(it==it_start)then
    if(irestart==0)then
      do n=sn,en
        int_H_U_psi_n(n)=0.0d0
      end do
    else
      write(filename,'("../data/restart_helicity_",i3.3,".dat")')myid
      open(unit=10*(1+myid),file=filename,form='unformatted')
      do n=sn,en
        read(10*(1+myid))int_H_U_psi_n(n)
      end do
      close(10*(1+myid))
    end if
  else
    do n=sn,en
      H_U_psi_n(n)=0.0d0
      H_U_psi_old_n(n)=0.0d0
    end do

    do l=l_min,l_max
      do n=sn,en
        do irho=irho_min,irho_max
          U(irho)=U_nl(irho,n,l)
          psi(irho)=psi_nl(irho,n,l)
          U_old(irho)=U_old_nl(irho,n,l)
          psi_old(irho)=psi_old_nl(irho,n,l)
        end do

        call conjugate(U,U_cc)
        call conjugate(U_old,U_old_cc)

        do irho=irho_min,irho_max
          H_U_psi_nl_den(irho)=dble(U_cc(irho)*psi(irho))
          H_U_psi_old_nl_den(irho)=dble(U_old_cc(irho)*psi_old(irho))
        end do

        H_U_psi_nl=0.0d0
        H_U_psi_old_nl=0.0d0
        do irho=irho_min+1,irho_max
          H_U_psi_nl=H_U_psi_nl+drho/2.0d0&
          *(H_U_psi_nl_den(irho-1)*rho(irho-1)&
          +H_U_psi_nl_den(irho)*rho(irho))
          H_U_psi_old_nl=H_U_psi_old_nl+drho/2.0d0&
          *(H_U_psi_old_nl_den(irho-1)*rho(irho-1)&
          +H_U_psi_old_nl_den(irho)*rho(irho))
        end do

        H_U_psi_n(n)=H_U_psi_n(n)+H_U_psi_nl
        H_U_psi_old_n(n)=H_U_psi_old_n(n)+H_U_psi_old_nl
      end do ! : n
    end do ! : l

    do n=sn,en
      int_H_U_psi_n(n)=int_H_U_psi_n(n)&
      +dt/2.0d0*(H_U_psi_n(n)+H_U_psi_old_n(n))
    end do

    if(mod(it,it_skip)==0)then
      do n=sn,en
        a_send_gather(n)=int_H_U_psi_n(n)
      end do
      call mpi_allgatherv(a_send_gather,scount_gather,mpi_double_complex,&
      a_recv_gather,rcounts_gather,rdispls_gather,mpi_double_complex,&
      mpi_comm_world,ierr)
      do n=0,n_max
        int_H_U_psi_n_gather(n)=a_recv_gather(n)
      end do

      if(myid==0)then
        write(filename,'("../data/int_helicity_spec_",i3.3,".dat")')it/it_skip
        write(6,*)filename
        open(unit=10,file=filename,form='formatted')
        write(10,fmt='(2a16)')'n','int_H_U_psi_n'
        n=0
        write(10,fmt='(2e16.4)')dble(n),int_H_U_psi_n_gather(n)
        do n=1,n_max
          write(10,fmt='(2e16.4)')dble(n),2.0d0*int_H_U_psi_n_gather(n)
        end do
        close(10)
      end if
    end if ! : mod(it,it_skip)==0
  end if ! : it==it_start

  if(it==it_end)then
    write(filename,'("../data/restart_helicity_",i3.3,".dat")')myid
    open(unit=10*(1+myid),file=filename,form='unformatted')
    do n=sn,en
      write(10*(1+myid))int_H_U_psi_n(n)
    end do
    close(10*(1+myid))
  end if

end subroutine helicity_spectrum
! ------------------------------------------------------------------------------
subroutine conjugate(a,a_cc)
! ------------------------------------------------------------------------------
  use coordinate
  implicit none
  complex(kind(0d0))::a(irho_min:irho_max)
  complex(kind(0d0))::a_cc(irho_min:irho_max)

  do irho=irho_min,irho_max
    a_cc(irho)=dble(a(irho))&
    +(0.0d0,1.0d0)*dble((0.0d0,1.0d0)*a(irho))
  end do

end subroutine conjugate
