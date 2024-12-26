! ------------------------------------------------------------------------------
! . last update @20240530
! ------------------------------------------------------------------------------
module coordinate
! ------------------------------------------------------------------------------
  implicit none
  double precision,parameter::pi=3.14159265358979323846d0
  integer::ir
  integer,parameter::ir_min=0
  integer,parameter::ir_max=256
  double precision::r(ir_min:ir_max)
  double precision,parameter::r_min=0.1d0
  double precision,parameter::r_max=1.0d0
  double precision,parameter::dr=(r_max-r_min)/dble(ir_max-ir_min)
  double precision,parameter::dr_sq=dr**2
  integer::n
  integer::izeta
  integer,parameter::n_max=40
  integer,parameter::izeta_min=0
  integer,parameter::izeta_max=128 ! : >3*n_max
  integer,parameter::numdata_zeta=izeta_max-izeta_min
  double precision,parameter::inumdata_zeta=1.0d0/dble(numdata_zeta)
  integer::m
  integer::itheta
  integer,parameter::m_max=170 ! : >m_res
  integer,parameter::m_min=-m_max
  integer,parameter::itheta_min=0
  integer,parameter::itheta_max=512 ! : >3*m_max
  integer,parameter::numdata_theta=itheta_max-itheta_min
  double precision::theta(0:numdata_theta)
  double precision,parameter::inumdata_theta=1.0d0/dble(numdata_theta)
  integer::l
  integer,parameter::l_max=10
  integer,parameter::l_min=-l_max
end module coordinate
! ------------------------------------------------------------------------------
module time
! ------------------------------------------------------------------------------
  implicit none
  double precision::t
  integer::it
  double precision,parameter::dt=0.01d0 ! : <dt_crit
  integer,parameter::irestart=0
  integer,parameter::it_span=50000
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
  double precision,parameter::mu=1.0d-5
  double precision,parameter::eta=1.0d-5
  double precision,parameter::chi=1.0d-5
end module parameter
! ------------------------------------------------------------------------------
module field
! ------------------------------------------------------------------------------
  use coordinate
  implicit none
  ! . equilibria
  double precision::U_eq(ir_min:ir_max)
  double precision::phi_eq(ir_min:ir_max)
  double precision::q(ir_min:ir_max)
  double precision::j_eq(ir_min:ir_max)
  double precision::p_eq(ir_min:ir_max)
  double precision::dU_eq(ir_min:ir_max)
  double precision::dphi_eq(ir_min:ir_max)
  double precision::dq(ir_min:ir_max)
  double precision::dj_eq(ir_min:ir_max)
  double precision::dp_eq(ir_min:ir_max)
  ! ..
  double precision::k_theta(ir_min:ir_max,m_min:m_max)
  double precision,allocatable::k_para(:,:,:)
  double precision::inv_r(ir_min:ir_max)
  double precision::k_theta_sq(ir_min:ir_max,m_min:m_max)
  double precision::sin_theta(0:numdata_theta-1)
  double precision::cos_theta(0:numdata_theta-1)
  complex(kind(0d0)),allocatable::shifter(:,:,:)
  ! . spectra
  complex(kind(0d0))::U(ir_min:ir_max)
  complex(kind(0d0))::phi(ir_min:ir_max)
  complex(kind(0d0))::psi(ir_min:ir_max)
  complex(kind(0d0))::j(ir_min:ir_max)
  complex(kind(0d0))::p(ir_min:ir_max)
  complex(kind(0d0)),allocatable::U_nm(:,:,:)
  complex(kind(0d0)),allocatable::phi_nm(:,:,:)
  complex(kind(0d0)),allocatable::psi_nm(:,:,:)
  complex(kind(0d0)),allocatable::j_nm(:,:,:)
  complex(kind(0d0)),allocatable::p_nm(:,:,:)
  ! ..
  complex(kind(0d0))::U_old(ir_min:ir_max)
  complex(kind(0d0))::phi_old(ir_min:ir_max)
  complex(kind(0d0))::psi_old(ir_min:ir_max)
  complex(kind(0d0))::j_old(ir_min:ir_max)
  complex(kind(0d0))::p_old(ir_min:ir_max)
  complex(kind(0d0)),allocatable::U_old_nm(:,:,:)
  complex(kind(0d0)),allocatable::phi_old_nm(:,:,:)
  complex(kind(0d0)),allocatable::psi_old_nm(:,:,:)
  complex(kind(0d0)),allocatable::j_old_nm(:,:,:)
  complex(kind(0d0)),allocatable::p_old_nm(:,:,:)
  ! . linear terms
  complex(kind(0d0))::bra_phi_U_lin(ir_min:ir_max)
  complex(kind(0d0))::nab_para_j_lin(ir_min:ir_max)
  complex(kind(0d0))::bra_x_p(ir_min:ir_max)
  complex(kind(0d0))::lap_perp_U(ir_min:ir_max)
  complex(kind(0d0))::nab_para_phi_lin(ir_min:ir_max)
  complex(kind(0d0))::bra_phi_p_lin(ir_min:ir_max)
  complex(kind(0d0))::bra_x_phi(ir_min:ir_max)
  complex(kind(0d0))::lap_perp_p(ir_min:ir_max)
  complex(kind(0d0)),allocatable::bra_phi_U_lin_nm(:,:,:)
  complex(kind(0d0)),allocatable::nab_para_j_lin_nm(:,:,:)
  complex(kind(0d0)),allocatable::bra_x_p_nm(:,:,:)
  complex(kind(0d0)),allocatable::lap_perp_U_nm(:,:,:)
  complex(kind(0d0)),allocatable::nab_para_phi_lin_nm(:,:,:)
  complex(kind(0d0)),allocatable::bra_phi_p_lin_nm(:,:,:)
  complex(kind(0d0)),allocatable::bra_x_phi_nm(:,:,:)
  complex(kind(0d0)),allocatable::lap_perp_p_nm(:,:,:)
  ! . nonlinear terms
  complex(kind(0d0))::bra_phi_U_non(ir_min:ir_max)
  complex(kind(0d0))::nab_para_j_non(ir_min:ir_max)
  complex(kind(0d0))::nab_para_phi_non(ir_min:ir_max)
  complex(kind(0d0))::bra_phi_p_non(ir_min:ir_max)
  complex(kind(0d0)),allocatable::bra_phi_U_non_nm(:,:,:)
  complex(kind(0d0)),allocatable::nab_para_j_non_nm(:,:,:)
  complex(kind(0d0)),allocatable::nab_para_phi_non_nm(:,:,:)
  complex(kind(0d0)),allocatable::bra_phi_p_non_nm(:,:,:)
  ! . time integrals
  integer::irk
  integer::jrk
  complex(kind(0d0))::U_rk(ir_min:ir_max,1:4)
  complex(kind(0d0))::psi_rk(ir_min:ir_max,1:4)
  complex(kind(0d0))::p_rk(ir_min:ir_max,1:4)
  complex(kind(0d0)),allocatable::U_rk_nm(:,:,:,:)
  complex(kind(0d0)),allocatable::psi_rk_nm(:,:,:,:)
  complex(kind(0d0)),allocatable::p_rk_nm(:,:,:,:)
  ! ..
  complex(kind(0d0))::del_r_U(ir_min:ir_max)
  complex(kind(0d0))::del_theta_U(ir_min:ir_max)
  complex(kind(0d0))::del_r_phi(ir_min:ir_max)
  complex(kind(0d0))::del_theta_phi(ir_min:ir_max)
  complex(kind(0d0))::del_r_psi(ir_min:ir_max)
  complex(kind(0d0))::del_theta_psi(ir_min:ir_max)
  complex(kind(0d0))::del_r_j(ir_min:ir_max)
  complex(kind(0d0))::del_theta_j(ir_min:ir_max)
  complex(kind(0d0))::del_r_p(ir_min:ir_max)
  complex(kind(0d0))::del_theta_p(ir_min:ir_max)
  complex(kind(0d0)),allocatable::del_r_U_nm(:,:,:)
  complex(kind(0d0)),allocatable::del_theta_U_nm(:,:,:)
  complex(kind(0d0)),allocatable::del_r_phi_nm(:,:,:)
  complex(kind(0d0)),allocatable::del_theta_phi_nm(:,:,:)
  complex(kind(0d0)),allocatable::del_r_psi_nm(:,:,:)
  complex(kind(0d0)),allocatable::del_theta_psi_nm(:,:,:)
  complex(kind(0d0)),allocatable::del_r_j_nm(:,:,:)
  complex(kind(0d0)),allocatable::del_theta_j_nm(:,:,:)
  complex(kind(0d0)),allocatable::del_r_p_nm(:,:,:)
  complex(kind(0d0)),allocatable::del_theta_p_nm(:,:,:)
  complex(kind(0d0)),allocatable::del_r_U_n(:,:,:)
  complex(kind(0d0)),allocatable::del_theta_U_n(:,:,:)
  complex(kind(0d0)),allocatable::del_r_phi_n(:,:,:)
  complex(kind(0d0)),allocatable::del_theta_phi_n(:,:,:)
  complex(kind(0d0)),allocatable::del_r_psi_n(:,:,:)
  complex(kind(0d0)),allocatable::del_theta_psi_n(:,:,:)
  complex(kind(0d0)),allocatable::del_r_j_n(:,:,:)
  complex(kind(0d0)),allocatable::del_theta_j_n(:,:,:)
  complex(kind(0d0)),allocatable::del_r_p_n(:,:,:)
  complex(kind(0d0)),allocatable::del_theta_p_n(:,:,:)
  complex(kind(0d0)),allocatable::del_r_U_n_trans(:,:,:)
  complex(kind(0d0)),allocatable::del_theta_U_n_trans(:,:,:)
  complex(kind(0d0)),allocatable::del_r_phi_n_trans(:,:,:)
  complex(kind(0d0)),allocatable::del_theta_phi_n_trans(:,:,:)
  complex(kind(0d0)),allocatable::del_r_psi_n_trans(:,:,:)
  complex(kind(0d0)),allocatable::del_theta_psi_n_trans(:,:,:)
  complex(kind(0d0)),allocatable::del_r_j_n_trans(:,:,:)
  complex(kind(0d0)),allocatable::del_theta_j_n_trans(:,:,:)
  complex(kind(0d0)),allocatable::del_r_p_n_trans(:,:,:)
  complex(kind(0d0)),allocatable::del_theta_p_n_trans(:,:,:)
  double precision::bra_phi_U_non_real(0:numdata_zeta-1)
  double precision::nab_para_j_non_real(0:numdata_zeta-1)
  double precision::nab_para_phi_non_real(0:numdata_zeta-1)
  double precision::bra_phi_p_non_real(0:numdata_zeta-1)
  complex(kind(0d0)),allocatable::bra_phi_U_non_n_trans(:,:,:)
  complex(kind(0d0)),allocatable::nab_para_j_non_n_trans(:,:,:)
  complex(kind(0d0)),allocatable::nab_para_phi_non_n_trans(:,:,:)
  complex(kind(0d0)),allocatable::bra_phi_p_non_n_trans(:,:,:)
  complex(kind(0d0)),allocatable::bra_phi_U_non_n(:,:,:)
  complex(kind(0d0)),allocatable::nab_para_j_non_n(:,:,:)
  complex(kind(0d0)),allocatable::nab_para_phi_non_n(:,:,:)
  complex(kind(0d0)),allocatable::bra_phi_p_non_n(:,:,:)
  integer(kind=8)::plan_c2r
  integer(kind=8)::plan_r2c
  complex(kind(0d0))::f1(0:numdata_zeta/2)
  complex(kind(0d0))::f2(0:numdata_zeta/2)
  complex(kind(0d0))::f3(0:numdata_zeta/2)
  complex(kind(0d0))::f4(0:numdata_zeta/2)
  complex(kind(0d0))::f5(0:numdata_zeta/2)
  complex(kind(0d0))::f6(0:numdata_zeta/2)
  complex(kind(0d0))::f7(0:numdata_zeta/2)
  complex(kind(0d0))::f8(0:numdata_zeta/2)
  complex(kind(0d0))::f9(0:numdata_zeta/2)
  complex(kind(0d0))::f10(0:numdata_zeta/2)
  double precision::f1_idft(0:numdata_zeta-1)
  double precision::f2_idft(0:numdata_zeta-1)
  double precision::f3_idft(0:numdata_zeta-1)
  double precision::f4_idft(0:numdata_zeta-1)
  double precision::f5_idft(0:numdata_zeta-1)
  double precision::f6_idft(0:numdata_zeta-1)
  double precision::f7_idft(0:numdata_zeta-1)
  double precision::f8_idft(0:numdata_zeta-1)
  double precision::f9_idft(0:numdata_zeta-1)
  double precision::f10_idft(0:numdata_zeta-1)
  ! ..
  complex(kind(0d0)),allocatable::bra_x_p_n(:,:,:)
  complex(kind(0d0)),allocatable::bra_x_phi_n(:,:,:)
  ! ..
  integer::ilu
  double precision::coe_lu(1:3,ir_min+1:ir_max-1,m_min:m_max)
  complex(kind(0d0))::rhs_lu(ir_min+1:ir_max-1)
  ! ..
  integer(kind=8)::plan_back
  integer(kind=8)::plan_for
  complex(kind(0d0))::f(0:numdata_theta-1)
  complex(kind(0d0))::f_dft(-numdata_theta/2:numdata_theta/2-1)
  ! . linear analysis
  integer::ir_peak
  double precision,allocatable::growth(:)
  double precision,allocatable::freq(:)
  double precision::growth_gather(0:n_max)
  double precision::freq_gather(0:n_max)
  complex(kind(0d0)),allocatable::p_old_n(:,:,:)
  double precision::p_amp(ir_min:ir_max)
  double precision::p_ir_peak
  ! . average
  double precision::U_avg(ir_min:ir_max)
  double precision::phi_avg(ir_min:ir_max)
  double precision::j_avg(ir_min:ir_max)
  double precision::p_avg(ir_min:ir_max)
  double precision::del_r_U_avg(ir_min:ir_max)
  double precision::del_r_phi_avg(ir_min:ir_max)
  double precision::del_r_j_avg(ir_min:ir_max)
  double precision::del_r_p_avg(ir_min:ir_max)
  ! . contours
  complex(kind(0d0)),allocatable::U_n(:,:,:)
  complex(kind(0d0)),allocatable::phi_n(:,:,:)
  complex(kind(0d0)),allocatable::psi_n(:,:,:)
  complex(kind(0d0)),allocatable::j_n(:,:,:)
  complex(kind(0d0)),allocatable::p_n(:,:,:)
  double precision::U_real(ir_min:ir_max,0:numdata_theta)
  double precision::phi_real(ir_min:ir_max,0:numdata_theta)
  double precision::psi_real(ir_min:ir_max,0:numdata_theta)
  double precision::j_real(ir_min:ir_max,0:numdata_theta)
  double precision::p_real(ir_min:ir_max,0:numdata_theta)
  ! . energy balances
  ! . kinetic energy
  complex(kind(0d0))::phi_cc(ir_min:ir_max)
  complex(kind(0d0))::phi_old_cc(ir_min:ir_max)
  double precision::dE_phi_nm_den(ir_min:ir_max)
  double precision::dE_phi_nm
  double precision,allocatable::dE_phi_n(:)
  double precision::dE_phi
  ! . linear energy sources
  double precision::S_phi_dphi_eq_nm_den(ir_min:ir_max)
  double precision::S_phi_q_nm_den(ir_min:ir_max)
  double precision::S_phi_dj_eq_nm_den(ir_min:ir_max)
  double precision::S_phi_x_p_nm_den(ir_min:ir_max)
  double precision::S_phi_mu_nm_den(ir_min:ir_max)
  double precision::S_phi_dphi_eq_nm
  double precision::S_phi_q_nm
  double precision::S_phi_dj_eq_nm
  double precision::S_phi_x_p_nm
  double precision::S_phi_mu_nm
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
  double precision::S_phi_phi_U_nm_den(ir_min:ir_max)
  double precision::S_phi_psi_j_nm_den(ir_min:ir_max)
  double precision::S_phi_phi_U_nm
  double precision::S_phi_psi_j_nm
  double precision,allocatable::S_phi_phi_U_n(:)
  double precision,allocatable::S_phi_psi_j_n(:)
  double precision::S_phi_phi_U
  double precision::S_phi_psi_j
  double precision::N_phi
  ! . magnetic energy
  complex(kind(0d0))::j_cc(ir_min:ir_max)
  complex(kind(0d0))::j_old_cc(ir_min:ir_max)
  double precision::dE_psi_nm_den(ir_min:ir_max)
  double precision::dE_psi_nm
  double precision,allocatable::dE_psi_n(:)
  double precision::dE_psi
  ! . linear energy sources
  double precision::S_j_q_nm_den(ir_min:ir_max)
  double precision::S_j_dphi_eq_nm_den(ir_min:ir_max)
  double precision::S_j_eta_nm_den(ir_min:ir_max)
  double precision::S_j_q_nm
  double precision::S_j_dphi_eq_nm
  double precision::S_j_eta_nm
  double precision,allocatable::S_j_q_n(:)
  double precision,allocatable::S_j_dphi_eq_n(:)
  double precision,allocatable::S_j_eta_n(:)
  double precision::S_j_q
  double precision::S_j_dphi_eq
  double precision::S_j_eta
  double precision::L_psi
  ! . nonlinear energy sources
  double precision::S_j_psi_phi_nm_den(ir_min:ir_max)
  double precision::S_j_psi_phi_nm
  double precision,allocatable::S_j_psi_phi_n(:)
  double precision::S_j_psi_phi
  double precision::N_psi
  ! . pressure energy
  double precision::dE_p_nm_den(ir_min:ir_max)
  double precision::dE_p_nm
  double precision,allocatable::dE_p_n(:)
  double precision::dE_p
  ! . linear energy sources
  complex(kind(0d0))::p_cc(ir_min:ir_max)
  double precision::S_p_dp_eq_nm_den(ir_min:ir_max)
  double precision::S_p_x_phi_nm_den(ir_min:ir_max)
  double precision::S_p_chi_nm_den(ir_min:ir_max)
  double precision::S_p_dp_eq_nm
  double precision::S_p_x_phi_nm
  double precision::S_p_chi_nm
  double precision,allocatable::S_p_dp_eq_n(:)
  double precision,allocatable::S_p_x_phi_n(:)
  double precision,allocatable::S_p_chi_n(:)
  double precision::S_p_dp_eq
  double precision::S_p_x_phi
  double precision::S_p_chi
  double precision::L_p
  ! . nonlinear energy sources
  double precision::S_p_phi_p_nm_den(ir_min:ir_max)
  double precision::S_p_phi_p_nm
  double precision,allocatable::S_p_phi_p_n(:)
  double precision::S_p_phi_p
  double precision::N_p
  ! . total energy
  double precision::dE_tot
  double precision::L_tot
  double precision::N_tot
  ! . energies
  ! . kinetic energy
  double precision::E_phi_nm_den(ir_min:ir_max)
  double precision,allocatable::E_phi_nm(:,:)
  double precision,allocatable::E_phi_n(:)
  double precision::E_phi
  double precision::E_phi_n_gather(0:n_max)
  double precision::E_phi_nm_gather(0:n_max,m_min:m_max)
  ! . magnetic energy
  double precision::E_psi_nm_den(ir_min:ir_max)
  double precision,allocatable::E_psi_nm(:,:)
  double precision,allocatable::E_psi_n(:)
  double precision::E_psi
  double precision::E_psi_n_gather(0:n_max)
  double precision::E_psi_nm_gather(0:n_max,m_min:m_max)
  ! . pressure energy
  double precision::E_p_nm_den(ir_min:ir_max)
  double precision,allocatable::E_p_nm(:,:)
  double precision,allocatable::E_p_n(:)
  double precision::E_p
  double precision::E_p_n_gather(0:n_max)
  double precision::E_p_nm_gather(0:n_max,m_min:m_max)
  ! . energy spectra
  ! . kinetic energy
  double precision::E_phi_old_nm_den(ir_min:ir_max)
  double precision::E_phi_old_nm
  double precision,allocatable::E_phi_old_n(:)
  double precision,allocatable::int_E_phi_n(:)
  double precision::int_E_phi_n_gather(0:n_max)
  ! . magnetic energy
  double precision::E_psi_old_nm_den(ir_min:ir_max)
  double precision::E_psi_old_nm
  double precision,allocatable::E_psi_old_n(:)
  double precision,allocatable::int_E_psi_n(:)
  double precision::int_E_psi_n_gather(0:n_max)
  ! . pressure energy
  double precision::E_p_old_nm_den(ir_min:ir_max)
  double precision::E_p_old_nm
  double precision,allocatable::E_p_old_n(:)
  double precision,allocatable::int_E_p_n(:)
  double precision::int_E_p_n_gather(0:n_max)
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
  integer::sir
  integer::eir
  integer,allocatable::sir0(:)
  integer,allocatable::eir0(:)
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
  call decompose_domain_r

end subroutine decompose_domains
! ------------------------------------------------------------------------------
subroutine allocate_variables
! ------------------------------------------------------------------------------
  use coordinate
  use field
  use mpi_global
  implicit none

  ! . equilibria
  allocate(k_para(ir_min:ir_max,sn:en,m_min:m_max))
  allocate(shifter(ir_min:ir_max,sn:en,0:numdata_theta-1))
  ! . spectra
  allocate(U_nm(ir_min:ir_max,sn:en,m_min:m_max))
  allocate(phi_nm(ir_min:ir_max,sn:en,m_min:m_max))
  allocate(psi_nm(ir_min:ir_max,sn:en,m_min:m_max))
  allocate(j_nm(ir_min:ir_max,sn:en,m_min:m_max))
  allocate(p_nm(ir_min:ir_max,sn:en,m_min:m_max))
  ! ..
  allocate(U_old_nm(ir_min:ir_max,sn:en,m_min:m_max))
  allocate(phi_old_nm(ir_min:ir_max,sn:en,m_min:m_max))
  allocate(psi_old_nm(ir_min:ir_max,sn:en,m_min:m_max))
  allocate(j_old_nm(ir_min:ir_max,sn:en,m_min:m_max))
  allocate(p_old_nm(ir_min:ir_max,sn:en,m_min:m_max))
  ! . linear terms
  allocate(bra_phi_U_lin_nm(ir_min:ir_max,sn:en,m_min:m_max))
  allocate(nab_para_j_lin_nm(ir_min:ir_max,sn:en,m_min:m_max))
  allocate(bra_x_p_nm(ir_min:ir_max,sn:en,m_min:m_max))
  allocate(lap_perp_U_nm(ir_min:ir_max,sn:en,m_min:m_max))
  allocate(nab_para_phi_lin_nm(ir_min:ir_max,sn:en,m_min:m_max))
  allocate(bra_phi_p_lin_nm(ir_min:ir_max,sn:en,m_min:m_max))
  allocate(bra_x_phi_nm(ir_min:ir_max,sn:en,m_min:m_max))
  allocate(lap_perp_p_nm(ir_min:ir_max,sn:en,m_min:m_max))
  ! . nonlinear terms
  allocate(bra_phi_U_non_nm(ir_min:ir_max,sn:en,m_min:m_max))
  allocate(nab_para_j_non_nm(ir_min:ir_max,sn:en,m_min:m_max))
  allocate(nab_para_phi_non_nm(ir_min:ir_max,sn:en,m_min:m_max))
  allocate(bra_phi_p_non_nm(ir_min:ir_max,sn:en,m_min:m_max))
  ! . time integrals
  allocate(U_rk_nm(ir_min:ir_max,sn:en,m_min:m_max,1:4))
  allocate(psi_rk_nm(ir_min:ir_max,sn:en,m_min:m_max,1:4))
  allocate(p_rk_nm(ir_min:ir_max,sn:en,m_min:m_max,1:4))
  ! ..
  allocate(del_r_U_nm(ir_min:ir_max,sn:en,m_min:m_max))
  allocate(del_theta_U_nm(ir_min:ir_max,sn:en,m_min:m_max))
  allocate(del_r_phi_nm(ir_min:ir_max,sn:en,m_min:m_max))
  allocate(del_theta_phi_nm(ir_min:ir_max,sn:en,m_min:m_max))
  allocate(del_r_psi_nm(ir_min:ir_max,sn:en,m_min:m_max))
  allocate(del_theta_psi_nm(ir_min:ir_max,sn:en,m_min:m_max))
  allocate(del_r_j_nm(ir_min:ir_max,sn:en,m_min:m_max))
  allocate(del_theta_j_nm(ir_min:ir_max,sn:en,m_min:m_max))
  allocate(del_r_p_nm(ir_min:ir_max,sn:en,m_min:m_max))
  allocate(del_theta_p_nm(ir_min:ir_max,sn:en,m_min:m_max))
  allocate(del_r_U_n(ir_min:ir_max,sn:en,0:numdata_theta-1))
  allocate(del_theta_U_n(ir_min:ir_max,sn:en,0:numdata_theta-1))
  allocate(del_r_phi_n(ir_min:ir_max,sn:en,0:numdata_theta-1))
  allocate(del_theta_phi_n(ir_min:ir_max,sn:en,0:numdata_theta-1))
  allocate(del_r_psi_n(ir_min:ir_max,sn:en,0:numdata_theta-1))
  allocate(del_theta_psi_n(ir_min:ir_max,sn:en,0:numdata_theta-1))
  allocate(del_r_j_n(ir_min:ir_max,sn:en,0:numdata_theta-1))
  allocate(del_theta_j_n(ir_min:ir_max,sn:en,0:numdata_theta-1))
  allocate(del_r_p_n(ir_min:ir_max,sn:en,0:numdata_theta-1))
  allocate(del_theta_p_n(ir_min:ir_max,sn:en,0:numdata_theta-1))
  allocate(del_r_U_n_trans(0:n_max,sir:eir,0:numdata_theta-1))
  allocate(del_theta_U_n_trans(0:n_max,sir:eir,0:numdata_theta-1))
  allocate(del_r_phi_n_trans(0:n_max,sir:eir,0:numdata_theta-1))
  allocate(del_theta_phi_n_trans(0:n_max,sir:eir,0:numdata_theta-1))
  allocate(del_r_psi_n_trans(0:n_max,sir:eir,0:numdata_theta-1))
  allocate(del_theta_psi_n_trans(0:n_max,sir:eir,0:numdata_theta-1))
  allocate(del_r_j_n_trans(0:n_max,sir:eir,0:numdata_theta-1))
  allocate(del_theta_j_n_trans(0:n_max,sir:eir,0:numdata_theta-1))
  allocate(del_r_p_n_trans(0:n_max,sir:eir,0:numdata_theta-1))
  allocate(del_theta_p_n_trans(0:n_max,sir:eir,0:numdata_theta-1))
  allocate(bra_phi_U_non_n_trans(0:n_max,sir:eir,0:numdata_theta-1))
  allocate(nab_para_j_non_n_trans(0:n_max,sir:eir,0:numdata_theta-1))
  allocate(nab_para_phi_non_n_trans(0:n_max,sir:eir,0:numdata_theta-1))
  allocate(bra_phi_p_non_n_trans(0:n_max,sir:eir,0:numdata_theta-1))
  allocate(bra_phi_U_non_n(ir_min:ir_max,sn:en,0:numdata_theta-1))
  allocate(nab_para_j_non_n(ir_min:ir_max,sn:en,0:numdata_theta-1))
  allocate(nab_para_phi_non_n(ir_min:ir_max,sn:en,0:numdata_theta-1))
  allocate(bra_phi_p_non_n(ir_min:ir_max,sn:en,0:numdata_theta-1))
  ! ..
  allocate(bra_x_p_n(ir_min:ir_max,sn:en,0:numdata_theta-1))
  allocate(bra_x_phi_n(ir_min:ir_max,sn:en,0:numdata_theta-1))
  ! . linear analysis
  allocate(growth(sn:en))
  allocate(freq(sn:en))
  allocate(p_old_n(ir_min:ir_max,sn:en,0:numdata_theta-1))
  ! . contours
  allocate(U_n(ir_min:ir_max,sn:en,0:numdata_theta-1))
  allocate(phi_n(ir_min:ir_max,sn:en,0:numdata_theta-1))
  allocate(psi_n(ir_min:ir_max,sn:en,0:numdata_theta-1))
  allocate(j_n(ir_min:ir_max,sn:en,0:numdata_theta-1))
  allocate(p_n(ir_min:ir_max,sn:en,0:numdata_theta-1))
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
  allocate(E_phi_nm(sn:en,m_min:m_max))
  allocate(E_phi_n(sn:en))
  allocate(E_psi_nm(sn:en,m_min:m_max))
  allocate(E_psi_n(sn:en))
  allocate(E_p_nm(sn:en,m_min:m_max))
  allocate(E_p_n(sn:en))
  ! . energy spectra
  allocate(E_phi_old_n(sn:en))
  allocate(int_E_phi_n(sn:en))
  allocate(E_psi_old_n(sn:en))
  allocate(int_E_psi_n(sn:en))
  allocate(E_p_old_n(sn:en))
  allocate(int_E_p_n(sn:en))

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

  do ir=ir_min,ir_max
    r(ir)=r_min+dble(ir-ir_min)*dr
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
  double precision::B_theta(ir_min:ir_max)
  double precision::dB_theta(ir_min:ir_max)
  double precision::d2B_theta(ir_min:ir_max)
  double precision::fact(ir_min:ir_max)
  double precision::B_theta_sq(ir_min:ir_max)
  integer::ipc
  double precision::B_theta_sq_pc(0:1)
  ! . linear analysis
  double precision::dp_eq_ir_peak

  ! . cylindrical MHD equilibrium
  do ir=ir_min,ir_max
    phi_eq(ir)=0.0d0
    q(ir)=2.0d0+2.0d0*r(ir)**2
    p_eq(ir)=beta/(2.0d0*eps)*(1.0d0-r(ir)**2/2.0d0)&
    *(1.0d0-tanh((r(ir)-0.8d0)/0.05d0))
  end do

  call lap_r_real(phi_eq,U_eq)
  call del_r_real(U_eq,dU_eq)
  call del_r_real(phi_eq,dphi_eq)
  call del_r_real(q,dq)
  call del_r_real(p_eq,dp_eq)

  do ir=ir_min,ir_max
    fact(ir)=(r(ir)*q(ir)*dq(ir)-q(ir)**2&
    +(eps*r(ir))**2)/r(ir)**3
  end do

  B_theta_sq(ir_min)=(r_min*B_zeta_ir_min/q(ir_min))**2

  do ir=ir_min+1,ir_max
    do ipc=0,1
      B_theta_sq_pc(ipc)=-2.0d0/((q(ir-1+ipc)/r(ir-1+ipc))**2+eps**2)&
      *(eps*dp_eq(ir-1+ipc)/2.0d0+fact(ir-1+ipc)*B_theta_sq(ir-1+ipc))

      if(ipc==1)then
        B_theta_sq(ir)=B_theta_sq(ir-1)&
        +dr/2.0d0*(B_theta_sq_pc(0)+B_theta_sq_pc(1))
      else
        B_theta_sq(ir)=B_theta_sq(ir-1)+dr*B_theta_sq_pc(ipc)
      end if
    end do ! : ipc
  end do ! : ir

  do ir=ir_min,ir_max
    B_theta(ir)=sqrt(B_theta_sq(ir))
  end do

  call del_r_real(B_theta,dB_theta)
  call del2_r_real(B_theta,d2B_theta)

  do ir=ir_min,ir_max
    j_eq(ir)=(r(ir)*dB_theta(ir)+B_theta(ir))/r(ir)
    dj_eq(ir)=(r(ir)**2*d2B_theta(ir)+r(ir)*dB_theta(ir)&
    -B_theta(ir))/r(ir)**2
  end do

  if(irestart==0)then
    if(myid==0)then
      open(unit=10,file='../data/equil_profile.dat',form='formatted')
      write(10,fmt='(6a16)')'r','U_eq','phi_eq','q','j_eq','p_eq'
      do ir=ir_min,ir_max
        write(10,fmt='(6e16.4)')r(ir),&
        U_eq(ir),phi_eq(ir),q(ir),j_eq(ir),p_eq(ir)
      end do
      close(10)

      open(unit=10,file='../data/equil_grad_profile.dat',form='formatted')
      write(10,fmt='(6a16)')'r','dU_eq','dphi_eq','dq','dj_eq','dp_eq'
      do ir=ir_min,ir_max
        write(10,fmt='(6e16.4)')r(ir),&
        dU_eq(ir),dphi_eq(ir),dq(ir),dj_eq(ir),dp_eq(ir)
      end do
      close(10)
    end if ! : myid==0
  end if ! : irestart==0

  do m=m_min,m_max
    do ir=ir_min,ir_max
      k_theta(ir,m)=dble(m)/r(ir)
      k_theta_sq(ir,m)=k_theta(ir,m)**2
    end do

    do n=sn,en
      do ir=ir_min,ir_max
        k_para(ir,n,m)=dble(m)/q(ir)-dble(n)
      end do
    end do
  end do ! : m

  do ir=ir_min,ir_max
    inv_r(ir)=1.0d0/r(ir)
  end do

  do itheta=0,numdata_theta-1
    sin_theta(itheta)=sin(theta(itheta))
    cos_theta(itheta)=cos(theta(itheta))

    do n=sn,en
      do ir=ir_min,ir_max
        shifter(ir,n,itheta)=exp((0.0d0,1.0d0)&
        *floor(dble(n)*q(ir))*theta(itheta))
      end do
    end do
  end do ! : itheta

  ir_peak=ir_min
  dp_eq_ir_peak=0.0d0
  do ir=ir_min,ir_max
    if(dp_eq_ir_peak<abs(dp_eq(ir)))then
      ir_peak=ir
      dp_eq_ir_peak=abs(dp_eq(ir))
    end if
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
  implicit none

  do m=m_min,m_max
    do ir=ir_min+1,ir_max-1
      coe_lu(1,ir,m)=1.0d0-dr/2.0d0*inv_r(ir)
      coe_lu(2,ir,m)=-2.0d0-dr_sq*k_theta_sq(ir,m)
      coe_lu(3,ir,m)=1.0d0+dr/2.0d0*inv_r(ir)
    end do

    do ir=ir_min+2,ir_max-1
      coe_lu(1,ir,m)=coe_lu(1,ir,m)/coe_lu(2,ir-1,m)
      coe_lu(2,ir,m)=coe_lu(2,ir,m)-coe_lu(1,ir,m)*coe_lu(3,ir-1,m)
    end do
  end do ! : m

end subroutine factorize
! ------------------------------------------------------------------------------
subroutine fftw_plans
! ------------------------------------------------------------------------------
  use coordinate
  use field
  implicit none
  include 'fftw3.f'

  call dfftw_plan_dft_1d(plan_back,numdata_theta,f,f,fftw_backward,fftw_estimate)
  call dfftw_plan_dft_1d(plan_for,numdata_theta,f,f,fftw_forward,fftw_estimate)
  call dfftw_plan_dft_c2r_1d(plan_c2r,numdata_zeta,f1,f1_idft,fftw_estimate)
  call dfftw_plan_dft_r2c_1d(plan_r2c,numdata_zeta,f1_idft,f1,fftw_estimate)

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
  complex(kind(0d0))::U_nl(ir_min:ir_max,sn:en,l_min:l_max)
  complex(kind(0d0))::phi_nl(ir_min:ir_max,sn:en,l_min:l_max)
  complex(kind(0d0))::psi_nl(ir_min:ir_max,sn:en,l_min:l_max)
  complex(kind(0d0))::j_nl(ir_min:ir_max,sn:en,l_min:l_max)
  complex(kind(0d0))::p_nl(ir_min:ir_max,sn:en,l_min:l_max)

  call random_seed(size=numseed)
  allocate(seed(0:numseed-1))
  do iseed=0,numseed-1
    ! call system_clock(count=seed(iseed))
    seed(iseed)=numseed*(1+myid)+iseed
  end do
  call random_seed(put=seed)

  do l=l_min,l_max
    do n=sn,en
      do ir=ir_min,ir_max
        U(ir)=0.0d0
        j(ir)=0.0d0
        p(ir)=0.0d0
      end do

      if(n/=0)then
        do irand=1,10
          call random_number(U_rand1)
          call random_number(U_rand2)
          call random_number(j_rand1)
          call random_number(j_rand2)
          call random_number(p_rand1)
          call random_number(p_rand2)

          do ir=ir_min,ir_max
            U(ir)=U(ir)+1.0d-14/(1.0d0+dble(n)**2)&
            *sqrt(-2.0d0*log(U_rand1))*exp((0.0d0,1.0d0)*2.0d0*pi*U_rand2)&
            *sin(dble(irand)*pi*(r(ir)-r_min)/(r_max-r_min))
            j(ir)=j(ir)+1.0d-14/(1.0d0+dble(n)**2)&
            *sqrt(-2.0d0*log(j_rand1))*exp((0.0d0,1.0d0)*2.0d0*pi*j_rand2)&
            *sin(dble(irand)*pi*(r(ir)-r_min)/(r_max-r_min))
            p(ir)=p(ir)+1.0d-16/(1.0d0+dble(n)**2)&
            *sqrt(-2.0d0*log(p_rand1))*exp((0.0d0,1.0d0)*2.0d0*pi*p_rand2)&
            *sin(dble(irand)*pi*(r(ir)-r_min)/(r_max-r_min))
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

      !     do ir=ir_min,ir_max
      !       U(ir)=U(ir)+1.0d-10/(1.0d0+dble(n)**2)&
      !       *sqrt(-2.0d0*log(U_rand1))*cos(2.0d0*pi*U_rand2)&
      !       *sin(dble(irand)*pi*(r(ir)-r_min)/(r_max-r_min))
      !       j(ir)=j(ir)+1.0d-10/(1.0d0+dble(n)**2)&
      !       *sqrt(-2.0d0*log(j_rand1))*cos(2.0d0*pi*j_rand2)&
      !       *sin(dble(irand)*pi*(r(ir)-r_min)/(r_max-r_min))
      !       p(ir)=p(ir)+1.0d-12/(1.0d0+dble(n)**2)&
      !       *sqrt(-2.0d0*log(p_rand1))*cos(2.0d0*pi*p_rand2)&
      !       *sin(dble(irand)*pi*(r(ir)-r_min)/(r_max-r_min))
      !     end do
      !   end do ! : irand
      ! end if ! : n/=0

      U(ir_min)=0.0d0
      U(ir_max)=0.0d0
      p(ir_min)=0.0d0
      p(ir_max)=0.0d0

      do ir=ir_min,ir_max
        U_nl(ir,n,l)=U(ir)
        j_nl(ir,n,l)=j(ir)
        p_nl(ir,n,l)=p(ir)
      end do
    end do ! : n
  end do ! : l

  call convert(U_nl,U_nm)
  call convert(j_nl,j_nm)
  call convert(p_nl,p_nm)

  if(myid==0)then
    call conjugate_0(U_nm)
    call conjugate_0(j_nm)
    call conjugate_0(p_nm)
  end if

  do m=m_min,m_max
    do n=sn,en
      do ir=ir_min,ir_max
        U(ir)=U_nm(ir,n,m)
        j(ir)=j_nm(ir,n,m)
      end do

      ! . definition of U
      call poisson_solver(U,phi)
      ! . Ampere's law
      call poisson_solver(-j,psi)

      do ir=ir_min,ir_max
        phi_nm(ir,n,m)=phi(ir)
        psi_nm(ir,n,m)=psi(ir)
      end do
    end do ! : n
  end do ! : m

  if(myid==0)then
    call conjugate_0(phi_nm)
    call conjugate_0(psi_nm)
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
  do m=m_min,m_max
    do n=sn,en
      do ir=ir_min,ir_max
        read(10*(1+myid))U_nm(ir,n,m),psi_nm(ir,n,m),p_nm(ir,n,m)
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

  do m=m_min,m_max
    do n=sn,en
      do ir=ir_min,ir_max
        U_old_nm(ir,n,m)=U_nm(ir,n,m)
        psi_old_nm(ir,n,m)=psi_nm(ir,n,m)
        p_old_nm(ir,n,m)=p_nm(ir,n,m)

        phi_old_nm(ir,n,m)=phi_nm(ir,n,m)
        j_old_nm(ir,n,m)=j_nm(ir,n,m)
      end do ! : ir
    end do ! : n
  end do ! : m

  do irk=1,4
    call linear_terms
    call nonlinear_terms

    do m=m_min,m_max
      do n=sn,en
        do ir=ir_min,ir_max
          U_old(ir)=U_old_nm(ir,n,m)
          psi_old(ir)=psi_old_nm(ir,n,m)
          p_old(ir)=p_old_nm(ir,n,m)
          ! . linear terms
          bra_phi_U_lin(ir)=bra_phi_U_lin_nm(ir,n,m)
          nab_para_j_lin(ir)=nab_para_j_lin_nm(ir,n,m)
          bra_x_p(ir)=bra_x_p_nm(ir,n,m)
          lap_perp_U(ir)=lap_perp_U_nm(ir,n,m)
          nab_para_phi_lin(ir)=nab_para_phi_lin_nm(ir,n,m)
          j(ir)=j_nm(ir,n,m)
          bra_phi_p_lin(ir)=bra_phi_p_lin_nm(ir,n,m)
          bra_x_phi(ir)=bra_x_phi_nm(ir,n,m)
          lap_perp_p(ir)=lap_perp_p_nm(ir,n,m)
          ! . nonlinear terms
          bra_phi_U_non(ir)=bra_phi_U_non_nm(ir,n,m)
          nab_para_j_non(ir)=nab_para_j_non_nm(ir,n,m)
          nab_para_phi_non(ir)=nab_para_phi_non_nm(ir,n,m)
          bra_phi_p_non(ir)=bra_phi_p_non_nm(ir,n,m)

          ! . vorticity equation
          U_rk(ir,irk)=-bra_phi_U_lin(ir)+nab_para_j_lin(ir)&
          -bra_x_p(ir)+mu*lap_perp_U(ir)&
          -bra_phi_U_non(ir)+nab_para_j_non(ir)
          ! . generalized Ohm's law
          psi_rk(ir,irk)=-nab_para_phi_lin(ir)-eta*j(ir)&
          -nab_para_phi_non(ir)
          ! . pressure evolution equation
          p_rk(ir,irk)=-bra_phi_p_lin(ir)+2.0d0*beta*bra_x_phi(ir)&
          +chi*lap_perp_p(ir)&
          -bra_phi_p_non(ir)

          U_rk_nm(ir,n,m,irk)=U_rk(ir,irk)
          psi_rk_nm(ir,n,m,irk)=psi_rk(ir,irk)
          p_rk_nm(ir,n,m,irk)=p_rk(ir,irk)
        end do ! : ir

        if(irk==4)then
          do jrk=1,4
            do ir=ir_min,ir_max
              U_rk(ir,jrk)=U_rk_nm(ir,n,m,jrk)
              psi_rk(ir,jrk)=psi_rk_nm(ir,n,m,jrk)
              p_rk(ir,jrk)=p_rk_nm(ir,n,m,jrk)
            end do
          end do

          do ir=ir_min,ir_max
            U(ir)=U_old(ir)&
            +dt_rk(irk)*(U_rk(ir,1)+2.0d0*U_rk(ir,2)&
            +2.0d0*U_rk(ir,3)+U_rk(ir,4))
            psi(ir)=psi_old(ir)&
            +dt_rk(irk)*(psi_rk(ir,1)+2.0d0*psi_rk(ir,2)&
            +2.0d0*psi_rk(ir,3)+psi_rk(ir,4))
            p(ir)=p_old(ir)&
            +dt_rk(irk)*(p_rk(ir,1)+2.0d0*p_rk(ir,2)&
            +2.0d0*p_rk(ir,3)+p_rk(ir,4))
          end do
        else
          do ir=ir_min,ir_max
            U(ir)=U_old(ir)+dt_rk(irk)*U_rk(ir,irk)
            psi(ir)=psi_old(ir)+dt_rk(irk)*psi_rk(ir,irk)
            p(ir)=p_old(ir)+dt_rk(irk)*p_rk(ir,irk)
          end do
        end if ! : irk==4

        call boundary_conditions

        do ir=ir_min,ir_max
          U_nm(ir,n,m)=U(ir)
          psi_nm(ir,n,m)=psi(ir)
          p_nm(ir,n,m)=p(ir)
        end do
      end do ! : n
    end do ! : m

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
  do m=m_min,m_max
    do n=sn,en
      do ir=ir_min,ir_max
        write(10*(1+myid))U_nm(ir,n,m),psi_nm(ir,n,m),p_nm(ir,n,m)
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
subroutine decompose_domain_r
! ------------------------------------------------------------------------------
  use coordinate
  use mpi_global
  use mpi
  implicit none
  integer::numdata
  integer::numfloor
  integer::numredistribute

  numdata=ir_max-ir_min+1
  numfloor=floor(dble(numdata)/dble(numprocs))
  numredistribute=numdata-numprocs*numfloor

  allocate(sir0(0:numprocs-1))
  allocate(eir0(0:numprocs-1))

  sir0(0)=ir_min
  do iprocs=1,numprocs-1
    sir0(iprocs)=ir_min+iprocs*numfloor
  end do

  do iprocs=1,numredistribute
    sir0(iprocs)=sir0(iprocs)+iprocs
  end do
  do iprocs=numredistribute+1,numprocs-1
    sir0(iprocs)=sir0(iprocs)+numredistribute
  end do

  do iprocs=0,numprocs-2
    eir0(iprocs)=sir0(iprocs+1)-1
  end do
  eir0(numprocs-1)=ir_max

  sir=sir0(myid)
  eir=eir0(myid)

  write(6,*)'myid=',myid,'sir=',sir,'eir=',eir,&
  'eir-sir+1=',eir-sir+1

  call mpi_barrier(mpi_comm_world,ierr)

end subroutine decompose_domain_r
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
    scounts_trans(iprocs)=(en-sn+1)*(eir0(iprocs)-sir0(iprocs)+1)
    rcounts_trans(iprocs)=(en0(iprocs)-sn0(iprocs)+1)*(eir-sir+1)
  end do

  allocate(sdispls_trans(0:numprocs-1))
  allocate(rdispls_trans(0:numprocs-1))

  sdispls_trans(0)=0
  rdispls_trans(0)=0
  do iprocs=1,numprocs-1
    sdispls_trans(iprocs)=sdispls_trans(iprocs-1)+scounts_trans(iprocs-1)
    rdispls_trans(iprocs)=rdispls_trans(iprocs-1)+rcounts_trans(iprocs-1)
  end do

  numdata_send=(en-sn+1)*(ir_max-ir_min+1)
  numdata_recv=(n_max-0+1)*(eir-sir+1)

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

  do m=m_min,m_max
    do n=sn,en
      do ir=ir_min,ir_max
        U(ir)=U_nm(ir,n,m)
        phi(ir)=phi_nm(ir,n,m)
        psi(ir)=psi_nm(ir,n,m)
        j(ir)=j_nm(ir,n,m)
        p(ir)=p_nm(ir,n,m)

        bra_phi_U_lin(ir)=(0.0d0,1.0d0)*k_theta(ir,m)&
        *(dphi_eq(ir)*U(ir)-dU_eq(ir)*phi(ir))
        nab_para_j_lin(ir)=(0.0d0,1.0d0)*(k_para(ir,n,m)*j(ir)&
        +k_theta(ir,m)*dj_eq(ir)*psi(ir))
        nab_para_phi_lin(ir)=(0.0d0,1.0d0)*(k_para(ir,n,m)*phi(ir)&
        +k_theta(ir,m)*dphi_eq(ir)*psi(ir))
        bra_phi_p_lin(ir)=(0.0d0,1.0d0)*k_theta(ir,m)&
        *(dphi_eq(ir)*p(ir)-dp_eq(ir)*phi(ir))
      end do ! : ir

      call lap_perp_spec(U,lap_perp_U)
      call lap_perp_spec(p,lap_perp_p)

      do ir=ir_min,ir_max
        bra_phi_U_lin_nm(ir,n,m)=bra_phi_U_lin(ir)
        nab_para_j_lin_nm(ir,n,m)=nab_para_j_lin(ir)
        lap_perp_U_nm(ir,n,m)=lap_perp_U(ir)
        nab_para_phi_lin_nm(ir,n,m)=nab_para_phi_lin(ir)
        bra_phi_p_lin_nm(ir,n,m)=bra_phi_p_lin(ir)
        lap_perp_p_nm(ir,n,m)=lap_perp_p(ir)
      end do
    end do ! : n
  end do ! : m

  do m=m_min,m_max
    do n=sn,en
      do ir=ir_min,ir_max
        phi(ir)=phi_nm(ir,n,m)
        p(ir)=p_nm(ir,n,m)
      end do

      call del_r_spec(phi,del_r_phi)
      call del_r_spec(p,del_r_p)

      do ir=ir_min,ir_max
        del_theta_phi(ir)=(0.0d0,1.0d0)*k_theta(ir,m)*phi(ir)
        del_theta_p(ir)=(0.0d0,1.0d0)*k_theta(ir,m)*p(ir)

        del_r_phi_nm(ir,n,m)=del_r_phi(ir)
        del_theta_phi_nm(ir,n,m)=del_theta_phi(ir)
        del_r_p_nm(ir,n,m)=del_r_p(ir)
        del_theta_p_nm(ir,n,m)=del_theta_p(ir)
      end do ! : ir
    end do ! : n
  end do ! : m

  call idft_theta(del_r_phi_nm,del_r_phi_n)
  call idft_theta(del_theta_phi_nm,del_theta_phi_n)
  call idft_theta(del_r_p_nm,del_r_p_n)
  call idft_theta(del_theta_p_nm,del_theta_p_n)

  do itheta=0,numdata_theta-1
    do n=sn,en
      do ir=ir_min,ir_max
        del_r_phi(ir)=del_r_phi_n(ir,n,itheta)
        del_theta_phi(ir)=del_theta_phi_n(ir,n,itheta)
        del_r_p(ir)=del_r_p_n(ir,n,itheta)
        del_theta_p(ir)=del_theta_p_n(ir,n,itheta)

        bra_x_p(ir)=del_r_p(ir)*sin_theta(itheta)&
        +del_theta_p(ir)*cos_theta(itheta)
        bra_x_phi(ir)=del_r_phi(ir)*sin_theta(itheta)&
        +del_theta_phi(ir)*cos_theta(itheta)

        bra_x_p_n(ir,n,itheta)=bra_x_p(ir)
        bra_x_phi_n(ir,n,itheta)=bra_x_phi(ir)
      end do ! : ir
    end do ! : n
  end do ! : itheta

  call dft_theta(bra_x_p_n,bra_x_p_nm)
  call dft_theta(bra_x_phi_n,bra_x_phi_nm)

end subroutine linear_terms
! ------------------------------------------------------------------------------
subroutine nonlinear_terms
! ------------------------------------------------------------------------------
  use coordinate
  use field
  use mpi_global
  implicit none
  include 'fftw3.f'

  do m=m_min,m_max
    do n=sn,en
      do ir=ir_min,ir_max
        U(ir)=U_nm(ir,n,m)
        psi(ir)=psi_nm(ir,n,m)
        j(ir)=j_nm(ir,n,m)
      end do

      call del_r_spec(U,del_r_U)
      call del_r_spec(psi,del_r_psi)
      call del_r_spec(j,del_r_j)

      do ir=ir_min,ir_max
        del_theta_U(ir)=(0.0d0,1.0d0)*k_theta(ir,m)*U(ir)
        del_theta_psi(ir)=(0.0d0,1.0d0)*k_theta(ir,m)*psi(ir)
        del_theta_j(ir)=(0.0d0,1.0d0)*k_theta(ir,m)*j(ir)

        del_r_U_nm(ir,n,m)=del_r_U(ir)
        del_theta_U_nm(ir,n,m)=del_theta_U(ir)
        del_r_psi_nm(ir,n,m)=del_r_psi(ir)
        del_theta_psi_nm(ir,n,m)=del_theta_psi(ir)
        del_r_j_nm(ir,n,m)=del_r_j(ir)
        del_theta_j_nm(ir,n,m)=del_theta_j(ir)
      end do ! : ir
    end do ! : n
  end do ! : m

  call idft_theta(del_r_U_nm,del_r_U_n)
  call idft_theta(del_theta_U_nm,del_theta_U_n)
  call idft_theta(del_r_psi_nm,del_r_psi_n)
  call idft_theta(del_theta_psi_nm,del_theta_psi_n)
  call idft_theta(del_r_j_nm,del_r_j_n)
  call idft_theta(del_theta_j_nm,del_theta_j_n)

  call transpose(del_r_U_n,del_r_U_n_trans)
  call transpose(del_theta_U_n,del_theta_U_n_trans)
  call transpose(del_r_phi_n,del_r_phi_n_trans)
  call transpose(del_theta_phi_n,del_theta_phi_n_trans)
  call transpose(del_r_psi_n,del_r_psi_n_trans)
  call transpose(del_theta_psi_n,del_theta_psi_n_trans)
  call transpose(del_r_j_n,del_r_j_n_trans)
  call transpose(del_theta_j_n,del_theta_j_n_trans)
  call transpose(del_r_p_n,del_r_p_n_trans)
  call transpose(del_theta_p_n,del_theta_p_n_trans)

  do itheta=0,numdata_theta-1
    do ir=sir,eir
      do izeta=0,numdata_zeta/2
        f1(izeta)=0.0d0
        f2(izeta)=0.0d0
        f3(izeta)=0.0d0
        f4(izeta)=0.0d0
        f5(izeta)=0.0d0
        f6(izeta)=0.0d0
        f7(izeta)=0.0d0
        f8(izeta)=0.0d0
        f9(izeta)=0.0d0
        f10(izeta)=0.0d0
      end do

      do n=0,n_max
        f1(n)=dble(del_r_U_n_trans(n,ir,itheta))&
        +(0.0d0,1.0d0)*dble((0.0d0,1.0d0)*del_r_U_n_trans(n,ir,itheta))
        f2(n)=dble(del_theta_U_n_trans(n,ir,itheta))&
        +(0.0d0,1.0d0)*dble((0.0d0,1.0d0)*del_theta_U_n_trans(n,ir,itheta))
        f3(n)=dble(del_r_phi_n_trans(n,ir,itheta))&
        +(0.0d0,1.0d0)*dble((0.0d0,1.0d0)*del_r_phi_n_trans(n,ir,itheta))
        f4(n)=dble(del_theta_phi_n_trans(n,ir,itheta))&
        +(0.0d0,1.0d0)*dble((0.0d0,1.0d0)*del_theta_phi_n_trans(n,ir,itheta))
        f5(n)=dble(del_r_psi_n_trans(n,ir,itheta))&
        +(0.0d0,1.0d0)*dble((0.0d0,1.0d0)*del_r_psi_n_trans(n,ir,itheta))
        f6(n)=dble(del_theta_psi_n_trans(n,ir,itheta))&
        +(0.0d0,1.0d0)*dble((0.0d0,1.0d0)*del_theta_psi_n_trans(n,ir,itheta))
        f7(n)=dble(del_r_j_n_trans(n,ir,itheta))&
        +(0.0d0,1.0d0)*dble((0.0d0,1.0d0)*del_r_j_n_trans(n,ir,itheta))
        f8(n)=dble(del_theta_j_n_trans(n,ir,itheta))&
        +(0.0d0,1.0d0)*dble((0.0d0,1.0d0)*del_theta_j_n_trans(n,ir,itheta))
        f9(n)=dble(del_r_p_n_trans(n,ir,itheta))&
        +(0.0d0,1.0d0)*dble((0.0d0,1.0d0)*del_r_p_n_trans(n,ir,itheta))
        f10(n)=dble(del_theta_p_n_trans(n,ir,itheta))&
        +(0.0d0,1.0d0)*dble((0.0d0,1.0d0)*del_theta_p_n_trans(n,ir,itheta))
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

      do izeta=0,numdata_zeta-1
        bra_phi_U_non_real(izeta)=f3_idft(izeta)*f2_idft(izeta)&
        -f4_idft(izeta)*f1_idft(izeta)
        nab_para_j_non_real(izeta)=f6_idft(izeta)*f7_idft(izeta)&
        -f5_idft(izeta)*f8_idft(izeta)
        nab_para_phi_non_real(izeta)=f6_idft(izeta)*f3_idft(izeta)&
        -f5_idft(izeta)*f4_idft(izeta)
        bra_phi_p_non_real(izeta)=f3_idft(izeta)*f10_idft(izeta)&
        -f4_idft(izeta)*f9_idft(izeta)
      end do

      do izeta=0,numdata_zeta-1
        f1_idft(izeta)=bra_phi_U_non_real(izeta)
        f2_idft(izeta)=nab_para_j_non_real(izeta)
        f3_idft(izeta)=nab_para_phi_non_real(izeta)
        f4_idft(izeta)=bra_phi_p_non_real(izeta)
      end do

      call dfftw_execute_dft_r2c(plan_r2c,f1_idft,f1)
      call dfftw_execute_dft_r2c(plan_r2c,f2_idft,f2)
      call dfftw_execute_dft_r2c(plan_r2c,f3_idft,f3)
      call dfftw_execute_dft_r2c(plan_r2c,f4_idft,f4)

      do n=0,n_max
        bra_phi_U_non_n_trans(n,ir,itheta)=(dble(f1(n))&
        +(0.0d0,1.0d0)*dble((0.0d0,1.0d0)*f1(n)))*inumdata_zeta
        nab_para_j_non_n_trans(n,ir,itheta)=(dble(f2(n))&
        +(0.0d0,1.0d0)*dble((0.0d0,1.0d0)*f2(n)))*inumdata_zeta
        nab_para_phi_non_n_trans(n,ir,itheta)=(dble(f3(n))&
        +(0.0d0,1.0d0)*dble((0.0d0,1.0d0)*f3(n)))*inumdata_zeta
        bra_phi_p_non_n_trans(n,ir,itheta)=(dble(f4(n))&
        +(0.0d0,1.0d0)*dble((0.0d0,1.0d0)*f4(n)))*inumdata_zeta
      end do
    end do ! : ir
  end do ! : itheta

  call itranspose(bra_phi_U_non_n_trans,bra_phi_U_non_n)
  call itranspose(nab_para_j_non_n_trans,nab_para_j_non_n)
  call itranspose(nab_para_phi_non_n_trans,nab_para_phi_non_n)
  call itranspose(bra_phi_p_non_n_trans,bra_phi_p_non_n)

  call dft_theta(bra_phi_U_non_n,bra_phi_U_non_nm)
  call dft_theta(nab_para_j_non_n,nab_para_j_non_nm)
  call dft_theta(nab_para_phi_non_n,nab_para_phi_non_nm)
  call dft_theta(bra_phi_p_non_n,bra_phi_p_non_nm)

end subroutine nonlinear_terms
! ------------------------------------------------------------------------------
subroutine phi_and_j
! ------------------------------------------------------------------------------
  use coordinate
  use field
  use mpi_global
  implicit none

  if(myid==0)then
    call conjugate_0(U_nm)
    call conjugate_0(psi_nm)
    call conjugate_0(p_nm)
  end if

  do m=m_min,m_max
    do n=sn,en
      do ir=ir_min,ir_max
        U(ir)=U_nm(ir,n,m)
        psi(ir)=psi_nm(ir,n,m)
      end do

      ! . definition of U
      call poisson_solver(U,phi)
      ! . Ampere's law
      call lap_perp_spec(-psi,j)

      do ir=ir_min,ir_max
        phi_nm(ir,n,m)=phi(ir)
        j_nm(ir,n,m)=j(ir)
      end do
    end do ! : n
  end do ! : m

  if(myid==0)then
    call conjugate_0(phi_nm)
    call conjugate_0(j_nm)
  end if

end subroutine phi_and_j
! ------------------------------------------------------------------------------
subroutine boundary_conditions
! ------------------------------------------------------------------------------
  use coordinate
  use field
  implicit none

  U(ir_min)=0.0d0
  U(ir_max)=0.0d0
  psi(ir_min)=0.0d0
  psi(ir_max)=0.0d0
  p(ir_min)=0.0d0
  p(ir_max)=0.0d0

end subroutine boundary_conditions
! ------------------------------------------------------------------------------
subroutine del_r_real(a,del_r_a)
! ------------------------------------------------------------------------------
  use coordinate
  implicit none
  double precision::a(ir_min:ir_max)
  double precision::del_r_a(ir_min:ir_max)

  do ir=ir_min+2,ir_max-2
    del_r_a(ir)=(-a(ir+2)+8.0d0*a(ir+1)&
    -8.0d0*a(ir-1)+a(ir-2))/(12.0d0*dr)
  end do
  ir=ir_min+1
  del_r_a(ir)=(a(ir+1)-a(ir-1))/(2.0d0*dr)
  ir=ir_max-1
  del_r_a(ir)=(a(ir+1)-a(ir-1))/(2.0d0*dr)
  ir=ir_min
  del_r_a(ir)=(-3.0d0*a(ir)+4.0d0*a(ir+1)&
  -a(ir+2))/(2.0d0*dr)
  ir=ir_max
  del_r_a(ir)=(3.0d0*a(ir)-4.0d0*a(ir-1)&
  +a(ir-2))/(2.0d0*dr)

end subroutine del_r_real
! ------------------------------------------------------------------------------
subroutine del2_r_real(a,del2_r_a)
! ------------------------------------------------------------------------------
  use coordinate
  implicit none
  double precision::a(ir_min:ir_max)
  double precision::del2_r_a(ir_min:ir_max)

  do ir=ir_min+2,ir_max-2
    del2_r_a(ir)=(-a(ir+2)+16.0d0*a(ir+1)-30.0d0*a(ir)&
    +16.0d0*a(ir-1)-a(ir-2))/(12.0d0*dr_sq)
  end do
  ir=ir_min+1
  del2_r_a(ir)=(a(ir+1)-2.0d0*a(ir)+a(ir-1))/dr_sq
  ir=ir_max-1
  del2_r_a(ir)=(a(ir+1)-2.0d0*a(ir)+a(ir-1))/dr_sq
  ir=ir_min
  del2_r_a(ir)=(a(ir)-2.0d0*a(ir+1)+a(ir+2))/dr_sq
  ir=ir_max
  del2_r_a(ir)=(a(ir)-2.0d0*a(ir-1)+a(ir-2))/dr_sq

end subroutine del2_r_real
! ------------------------------------------------------------------------------
subroutine lap_r_real(a,lap_r_a)
! ------------------------------------------------------------------------------
  use coordinate
  implicit none
  double precision::a(ir_min:ir_max)
  double precision::lap_r_a(ir_min:ir_max)
  double precision::del_r_a(ir_min:ir_max)
  double precision::del2_r_a(ir_min:ir_max)

  call del_r_real(a,del_r_a)
  call del2_r_real(a,del2_r_a)

  do ir=ir_min,ir_max
    lap_r_a(ir)=del2_r_a(ir)+del_r_a(ir)/r(ir)
  end do

end subroutine lap_r_real
! ------------------------------------------------------------------------------
subroutine del_r_spec(a,del_r_a)
! ------------------------------------------------------------------------------
  use coordinate
  implicit none
  complex(kind(0d0))::a(ir_min:ir_max)
  complex(kind(0d0))::del_r_a(ir_min:ir_max)

  do ir=ir_min+2,ir_max-2
    del_r_a(ir)=(-a(ir+2)+8.0d0*a(ir+1)&
    -8.0d0*a(ir-1)+a(ir-2))/(12.0d0*dr)
  end do
  ir=ir_min+1
  del_r_a(ir)=(a(ir+1)-a(ir-1))/(2.0d0*dr)
  ir=ir_max-1
  del_r_a(ir)=(a(ir+1)-a(ir-1))/(2.0d0*dr)
  ir=ir_min
  del_r_a(ir)=(-3.0d0*a(ir)+4.0d0*a(ir+1)&
  -a(ir+2))/(2.0d0*dr)
  ir=ir_max
  del_r_a(ir)=(3.0d0*a(ir)-4.0d0*a(ir-1)&
  +a(ir-2))/(2.0d0*dr)

end subroutine del_r_spec
! ------------------------------------------------------------------------------
subroutine del2_r_spec(a,del2_r_a)
! ------------------------------------------------------------------------------
  use coordinate
  implicit none
  complex(kind(0d0))::a(ir_min:ir_max)
  complex(kind(0d0))::del2_r_a(ir_min:ir_max)

  do ir=ir_min+2,ir_max-2
    del2_r_a(ir)=(-a(ir+2)+16.0d0*a(ir+1)-30.0d0*a(ir)&
    +16.0d0*a(ir-1)-a(ir-2))/(12.0d0*dr_sq)
  end do
  ir=ir_min+1
  del2_r_a(ir)=(a(ir+1)-2.0d0*a(ir)+a(ir-1))/dr_sq
  ir=ir_max-1
  del2_r_a(ir)=(a(ir+1)-2.0d0*a(ir)+a(ir-1))/dr_sq
  ir=ir_min
  del2_r_a(ir)=(a(ir)-2.0d0*a(ir+1)+a(ir+2))/dr_sq
  ir=ir_max
  del2_r_a(ir)=(a(ir)-2.0d0*a(ir-1)+a(ir-2))/dr_sq

end subroutine del2_r_spec
! ------------------------------------------------------------------------------
subroutine lap_perp_spec(a,lap_perp_a)
! ------------------------------------------------------------------------------
  use coordinate
  use field
  implicit none
  complex(kind(0d0))::a(ir_min:ir_max)
  complex(kind(0d0))::lap_perp_a(ir_min:ir_max)
  complex(kind(0d0))::del_r_a(ir_min:ir_max)
  complex(kind(0d0))::del2_r_a(ir_min:ir_max)

  call del_r_spec(a,del_r_a)
  call del2_r_spec(a,del2_r_a)

  do ir=ir_min,ir_max
    lap_perp_a(ir)=del2_r_a(ir)+inv_r(ir)*del_r_a(ir)&
    -k_theta_sq(ir,m)*a(ir)
  end do

end subroutine lap_perp_spec
! ------------------------------------------------------------------------------
subroutine poisson_solver(a,ilap_perp_a)
! ------------------------------------------------------------------------------
  use coordinate
  use field
  implicit none
  complex(kind(0d0))::a(ir_min:ir_max)
  complex(kind(0d0))::ilap_perp_a(ir_min:ir_max)

  do ir=ir_min+1,ir_max-1
    rhs_lu(ir)=dr_sq*a(ir)
  end do

  do ir=ir_min+2,ir_max-1
    rhs_lu(ir)=rhs_lu(ir)-coe_lu(1,ir,m)*rhs_lu(ir-1)
  end do

  ir=ir_max-1
  rhs_lu(ir)=rhs_lu(ir)/coe_lu(2,ir,m)
  do ir=ir_max-2,ir_min+1,-1
    rhs_lu(ir)=(rhs_lu(ir)&
    -coe_lu(3,ir,m)*rhs_lu(ir+1))/coe_lu(2,ir,m)
  end do

  do ir=ir_min+1,ir_max-1
    ilap_perp_a(ir)=rhs_lu(ir)
  end do

  ilap_perp_a(ir_min)=0.0d0
  ilap_perp_a(ir_max)=0.0d0

end subroutine poisson_solver
! ------------------------------------------------------------------------------
subroutine transpose(a_n,a_n_trans)
! ------------------------------------------------------------------------------
  use coordinate
  use mpi_global
  use mpi
  implicit none
  complex(kind(0d0))::a_n(ir_min:ir_max,sn:en,0:numdata_theta-1)
  complex(kind(0d0))::a_n_trans(0:n_max,sir:eir,0:numdata_theta-1)

  do itheta=0,numdata_theta-1
    idata=0
    do iprocs=0,numprocs-1
      do n=sn,en
        do ir=sir0(iprocs),eir0(iprocs)
          a_send_trans(idata)=a_n(ir,n,itheta)
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
        do ir=sir,eir
          a_n_trans(n,ir,itheta)=a_recv_trans(idata)
          idata=idata+1
        end do
      end do
    end do
  end do ! : itheta

end subroutine transpose
! ------------------------------------------------------------------------------
subroutine itranspose(a_n_trans,a_n)
! ------------------------------------------------------------------------------
  use coordinate
  use mpi_global
  use mpi
  implicit none
  complex(kind(0d0))::a_n_trans(0:n_max,sir:eir,0:numdata_theta-1)
  complex(kind(0d0))::a_n(ir_min:ir_max,sn:en,0:numdata_theta-1)

  do itheta=0,numdata_theta-1
    idata=0
    do iprocs=0,numprocs-1
      do ir=sir,eir
        do n=sn0(iprocs),en0(iprocs)
          a_recv_trans(idata)=a_n_trans(n,ir,itheta)
          idata=idata+1
        end do
      end do
    end do

    call mpi_alltoallv(a_recv_trans,rcounts_trans,rdispls_trans,mpi_double_complex,&
    a_send_trans,scounts_trans,sdispls_trans,mpi_double_complex,&
    mpi_comm_world,ierr)

    idata=0
    do iprocs=0,numprocs-1
      do ir=sir0(iprocs),eir0(iprocs)
        do n=sn,en
          a_n(ir,n,itheta)=a_send_trans(idata)
          idata=idata+1
        end do
      end do
    end do
  end do ! : itheta

end subroutine itranspose
! ------------------------------------------------------------------------------
subroutine conjugate_0(a_nm)
! ------------------------------------------------------------------------------
  use coordinate
  use mpi_global
  implicit none
  complex(kind(0d0))::a_nm(ir_min:ir_max,sn:en,m_min:m_max)

  n=0
  do m=1,m_max
    do ir=ir_min,ir_max
      a_nm(ir,-n,-m)=dble(a_nm(ir,n,m))&
      +(0.0d0,1.0d0)*dble((0.0d0,1.0d0)*a_nm(ir,n,m))
    end do
  end do
  m=0
  do ir=ir_min,ir_max
    a_nm(ir,-n,-m)=dble(a_nm(ir,n,m))
  end do

end subroutine conjugate_0
! ------------------------------------------------------------------------------
subroutine idft_theta(a_nm,a_n)
! ------------------------------------------------------------------------------
  use coordinate
  use field
  use mpi_global
  implicit none
  include 'fftw3.f'
  complex(kind(0d0))::a_nm(ir_min:ir_max,sn:en,m_min:m_max)
  complex(kind(0d0))::a_n(ir_min:ir_max,sn:en,0:numdata_theta-1)

  do n=sn,en
    do ir=ir_min,ir_max
      do itheta=0,numdata_theta-1
        f(itheta)=0.0d0
      end do

      do m=0,m_max
        f(m)=a_nm(ir,n,m)
      end do
      do m=m_min,-1
        f(m+numdata_theta)=a_nm(ir,n,m)
      end do

      call dfftw_execute_dft(plan_back,f,f)

      do itheta=0,numdata_theta-1
        a_n(ir,n,itheta)=f(itheta)
      end do
    end do ! : ir
  end do ! : n

end subroutine idft_theta
! ------------------------------------------------------------------------------
subroutine dft_theta(a_n,a_nm)
! ------------------------------------------------------------------------------
  use coordinate
  use field
  use mpi_global
  implicit none
  include 'fftw3.f'
  complex(kind(0d0))::a_n(ir_min:ir_max,sn:en,0:numdata_theta-1)
  complex(kind(0d0))::a_nm(ir_min:ir_max,sn:en,m_min:m_max)

  do n=sn,en
    do ir=ir_min,ir_max
      do itheta=0,numdata_theta-1
        f(itheta)=a_n(ir,n,itheta)
      end do

      call dfftw_execute_dft(plan_for,f,f)

      do m=0,numdata_theta/2-1
        f_dft(m)=f(m)*inumdata_theta
      end do
      do m=numdata_theta/2,numdata_theta-1
        f_dft(m-numdata_theta)=f(m)*inumdata_theta
      end do

      do m=m_min,m_max
        a_nm(ir,n,m)=f_dft(m)
      end do
    end do ! : ir
  end do ! : n

end subroutine dft_theta
! ------------------------------------------------------------------------------
subroutine convert(a_nl,a_nm)
! ------------------------------------------------------------------------------
  use coordinate
  use field
  use mpi_global
  implicit none
  include 'fftw3.f'
  complex(kind(0d0))::a_nl(ir_min:ir_max,sn:en,l_min:l_max)
  complex(kind(0d0))::a_nm(ir_min:ir_max,sn:en,m_min:m_max)
  complex(kind(0d0))::a_n(ir_min:ir_max,sn:en,0:numdata_theta-1)

  do n=sn,en
    do ir=ir_min,ir_max
      do itheta=0,numdata_theta-1
        f(itheta)=0.0d0
      end do

      do l=0,l_max
        f(l)=a_nl(ir,n,l)
      end do
      do l=l_min,-1
        f(l+numdata_theta)=a_nl(ir,n,l)
      end do

      call dfftw_execute_dft(plan_back,f,f)

      do itheta=0,numdata_theta-1
        a_n(ir,n,itheta)=f(itheta)
      end do
    end do ! : ir
  end do ! : n

  do itheta=0,numdata_theta-1
    do n=sn,en
      do ir=ir_min,ir_max
        a_n(ir,n,itheta)=a_n(ir,n,itheta)&
        *shifter(ir,n,itheta)
      end do
    end do
  end do

  call dft_theta(a_n,a_nm)

end subroutine convert
! ------------------------------------------------------------------------------
subroutine linear_analysis
! ------------------------------------------------------------------------------
  use time
  use coordinate
  use field
  use mpi_global
  use mpi
  implicit none

  call idft_theta(p_nm,p_n)
  call idft_theta(p_old_nm,p_old_n)

  if(it==it_start)then
    do n=sn,en
      growth(n)=0.0d0
      freq(n)=0.0d0
    end do
  else
    itheta=numdata_theta/2
    ir=ir_peak
    do n=sn,en
      p(ir)=p_n(ir,n,itheta)
      p_old(ir)=p_old_n(ir,n,itheta)

      growth(n)=dble(log(p(ir)/p_old(ir)))/dt
      freq(n)=dble((0.0d0,1.0d0)*log(p(ir)/p_old(ir)))/dt
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

  itheta=numdata_theta/2
  do n=sn,en
    if(n==n_max/2.or.n==n_max/4)then
      do ir=ir_min,ir_max
        p_amp(ir)=abs(p_n(ir,n,itheta))
      end do

      p_ir_peak=p_amp(ir_peak)
      do ir=ir_min,ir_max
        p(ir)=p_amp(ir)/p_ir_peak&
        *exp(-(0.0d0,1.0d0)*freq(n)*t)
      end do

      write(filename,'("../data/eigfunc_",i3.3,"_",i3.3,".dat")')n,it/it_skip
      write(6,*)filename
      open(unit=10*(1+myid),file=filename,form='formatted')
      write(10*(1+myid),fmt='(3a16)')'r','Re(p_n)','Im(p_n)'
      do ir=ir_min,ir_max
        write(10*(1+myid),fmt='(3e16.4)')r(ir),&
        dble(p(ir)),dble(-(0.0d0,1.0d0)*p(ir))
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

  do ir=ir_min,ir_max
    U_avg(ir)=U_eq(ir)+dble(U_nm(ir,0,0))
    phi_avg(ir)=phi_eq(ir)+dble(phi_nm(ir,0,0))
    j_avg(ir)=j_eq(ir)+dble(j_nm(ir,0,0))
    p_avg(ir)=p_eq(ir)+dble(p_nm(ir,0,0))
  end do

  call del_r_real(U_avg,del_r_U_avg)
  call del_r_real(phi_avg,del_r_phi_avg)
  call del_r_real(j_avg,del_r_j_avg)
  call del_r_real(p_avg,del_r_p_avg)

  write(filename,'("../data/avg_profile_",i3.3,".dat")')it/it_skip
  write(6,*)filename
  open(unit=10,file=filename,form='formatted')
  write(10,fmt='(5a16)')'r','U_avg','phi_avg','j_avg','p_avg'
  do ir=ir_min,ir_max
    write(10,fmt='(5e16.4)')r(ir),&
    U_avg(ir),phi_avg(ir),j_avg(ir),p_avg(ir)
  end do
  close(10)

  write(filename,'("../data/avg_grad_profile_",i3.3,".dat")')it/it_skip
  write(6,*)filename
  open(unit=10,file=filename,form='formatted')
  write(10,fmt='(5a16)')'r','del_r_U_avg','del_r_phi_avg',&
  'del_r_j_avg','del_r_p_avg'
  do ir=ir_min,ir_max
    write(10,fmt='(5e16.4)')r(ir),&
    del_r_U_avg(ir),del_r_phi_avg(ir),del_r_j_avg(ir),del_r_p_avg(ir)
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

  ! call idft_theta(U_nm,U_n)
  call idft_theta(phi_nm,phi_n)
  ! call idft_theta(psi_nm,psi_n)
  ! call idft_theta(j_nm,j_n)
  call idft_theta(p_nm,p_n)

  do itheta=0,numdata_theta-1
    do ir=ir_min,ir_max
      ! U_real(ir,itheta)=0.0d0
      phi_real(ir,itheta)=0.0d0
      ! psi_real(ir,itheta)=0.0d0
      ! j_real(ir,itheta)=0.0d0
      p_real(ir,itheta)=0.0d0
    end do

    do n=sn,en
      if(n==0)then
        do ir=ir_min,ir_max
          ! U_real(ir,itheta)=U_real(ir,itheta)&
          ! +dble(U_n(ir,n,itheta))
          phi_real(ir,itheta)=phi_real(ir,itheta)&
          +dble(phi_n(ir,n,itheta))
          ! psi_real(ir,itheta)=psi_real(ir,itheta)&
          ! +dble(psi_n(ir,n,itheta))
          ! j_real(ir,itheta)=j_real(ir,itheta)&
          ! +dble(j_n(ir,n,itheta))
          p_real(ir,itheta)=p_real(ir,itheta)&
          +dble(p_n(ir,n,itheta))
        end do
      else
        do ir=ir_min,ir_max
          ! U_real(ir,itheta)=U_real(ir,itheta)&
          ! +2.0d0*dble(U_n(ir,n,itheta))
          phi_real(ir,itheta)=phi_real(ir,itheta)&
          +2.0d0*dble(phi_n(ir,n,itheta))
          ! psi_real(ir,itheta)=psi_real(ir,itheta)&
          ! +2.0d0*dble(psi_n(ir,n,itheta))
          ! j_real(ir,itheta)=j_real(ir,itheta)&
          ! +2.0d0*dble(j_n(ir,n,itheta))
          p_real(ir,itheta)=p_real(ir,itheta)&
          +2.0d0*dble(p_n(ir,n,itheta))
        end do
      end if
    end do
  end do ! : itheta

  do itheta=0,numdata_theta-1
    do ir=ir_min,ir_max
      ! a_send_op=U_real(ir,itheta)
      ! call mpi_allreduce(a_send_op,a_recv_op,1,mpi_double_precision,mpi_sum,&
      ! mpi_comm_world,ierr)
      ! U_real(ir,itheta)=a_recv_op

      a_send_op=phi_real(ir,itheta)
      call mpi_allreduce(a_send_op,a_recv_op,1,mpi_double_precision,mpi_sum,&
      mpi_comm_world,ierr)
      phi_real(ir,itheta)=a_recv_op

      ! a_send_op=psi_real(ir,itheta)
      ! call mpi_allreduce(a_send_op,a_recv_op,1,mpi_double_precision,mpi_sum,&
      ! mpi_comm_world,ierr)
      ! psi_real(ir,itheta)=a_recv_op

      ! a_send_op=j_real(ir,itheta)
      ! call mpi_allreduce(a_send_op,a_recv_op,1,mpi_double_precision,mpi_sum,&
      ! mpi_comm_world,ierr)
      ! j_real(ir,itheta)=a_recv_op

      a_send_op=p_real(ir,itheta)
      call mpi_allreduce(a_send_op,a_recv_op,1,mpi_double_precision,mpi_sum,&
      mpi_comm_world,ierr)
      p_real(ir,itheta)=a_recv_op
    end do ! : ir
  end do ! : itheta

  do ir=ir_min,ir_max
    ! U_real(ir,numdata_theta)=U_real(ir,0)
    phi_real(ir,numdata_theta)=phi_real(ir,0)
    ! psi_real(ir,numdata_theta)=psi_real(ir,0)
    ! j_real(ir,numdata_theta)=j_real(ir,0)
    p_real(ir,numdata_theta)=p_real(ir,0)
  end do

  if(myid==0)then
    ! write(filename,'("../data/U_contour_",i3.3,".dat")')it/it_skip
    ! write(6,*)filename
    ! open(unit=10,file=filename,form='formatted')
    ! write(10,fmt='(3a16)')'x','Z','U_real'
    ! do itheta=0,numdata_theta
    !   do ir=ir_min,ir_max
    !     write(10,fmt='(3e16.4)')r(ir)*cos(theta(itheta)),&
    !     r(ir)*sin(theta(itheta)),U_real(ir,itheta)
    !   end do
    !   write(10,*)''
    ! end do
    ! close(10)

    write(filename,'("../data/phi_contour_",i3.3,".dat")')it/it_skip
    write(6,*)filename
    open(unit=10,file=filename,form='formatted')
    write(10,fmt='(3a16)')'x','Z','phi_real'
    do itheta=0,numdata_theta
      do ir=ir_min,ir_max
        write(10,fmt='(3e16.4)')r(ir)*cos(theta(itheta)),&
        r(ir)*sin(theta(itheta)),phi_real(ir,itheta)
      end do
      write(10,*)''
    end do
    close(10)

    ! write(filename,'("../data/psi_contour_",i3.3,".dat")')it/it_skip
    ! write(6,*)filename
    ! open(unit=10,file=filename,form='formatted')
    ! write(10,fmt='(3a16)')'x','Z','psi_real'
    ! do itheta=0,numdata_theta
    !   do ir=ir_min,ir_max
    !     write(10,fmt='(3e16.4)')r(ir)*cos(theta(itheta)),&
    !     r(ir)*sin(theta(itheta)),psi_real(ir,itheta)
    !   end do
    !   write(10,*)''
    ! end do
    ! close(10)

    ! write(filename,'("../data/j_contour_",i3.3,".dat")')it/it_skip
    ! write(6,*)filename
    ! open(unit=10,file=filename,form='formatted')
    ! write(10,fmt='(3a16)')'x','Z','j_real'
    ! do itheta=0,numdata_theta
    !   do ir=ir_min,ir_max
    !     write(10,fmt='(3e16.4)')r(ir)*cos(theta(itheta)),&
    !     r(ir)*sin(theta(itheta)),j_real(ir,itheta)
    !   end do
    !   write(10,*)''
    ! end do
    ! close(10)

    write(filename,'("../data/p_contour_",i3.3,".dat")')it/it_skip
    write(6,*)filename
    open(unit=10,file=filename,form='formatted')
    write(10,fmt='(3a16)')'x','Z','p_real'
    do itheta=0,numdata_theta
      do ir=ir_min,ir_max
        write(10,fmt='(3e16.4)')r(ir)*cos(theta(itheta)),&
        r(ir)*sin(theta(itheta)),p_real(ir,itheta)
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

  do m=m_min,m_max
    do n=sn,en
      do ir=ir_min,ir_max
        U(ir)=U_nm(ir,n,m)
        phi(ir)=phi_nm(ir,n,m)
        U_old(ir)=U_old_nm(ir,n,m)
        phi_old(ir)=phi_old_nm(ir,n,m)
      end do

      call conjugate(phi,phi_cc)
      call conjugate(phi_old,phi_old_cc)

      do ir=ir_min,ir_max
        dE_phi_nm_den(ir)=-(dble(phi_cc(ir)*U(ir))&
        -dble(phi_old_cc(ir)*U_old(ir)))/(2.0d0*dt)
      end do

      dE_phi_nm=0.0d0
      do ir=ir_min+1,ir_max
        dE_phi_nm=dE_phi_nm+dr/2.0d0&
        *(dE_phi_nm_den(ir-1)*r(ir-1)&
        +dE_phi_nm_den(ir)*r(ir))
      end do

      dE_phi_n(n)=dE_phi_n(n)+dE_phi_nm
    end do ! : n
  end do ! : m

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

  do m=m_min,m_max
    do n=sn,en
      do ir=ir_min,ir_max
        U(ir)=U_nm(ir,n,m)
        phi(ir)=phi_nm(ir,n,m)
        psi(ir)=psi_nm(ir,n,m)
        j(ir)=j_nm(ir,n,m)
        bra_x_p(ir)=bra_x_p_nm(ir,n,m)
        lap_perp_U(ir)=lap_perp_U_nm(ir,n,m)
      end do

      call conjugate(phi,phi_cc)

      do ir=ir_min,ir_max
        S_phi_dphi_eq_nm_den(ir)=dble((0.0d0,1.0d0)*k_theta(ir,m)&
        *dphi_eq(ir)*phi_cc(ir)*U(ir))
        S_phi_q_nm_den(ir)=-dble((0.0d0,1.0d0)*k_para(ir,n,m)&
        *phi_cc(ir)*j(ir))
        S_phi_dj_eq_nm_den(ir)=-dble((0.0d0,1.0d0)*k_theta(ir,m)&
        *dj_eq(ir)*phi_cc(ir)*psi(ir))
        S_phi_x_p_nm_den(ir)=dble(phi_cc(ir)*bra_x_p(ir))
        S_phi_mu_nm_den(ir)=-mu*dble(phi_cc(ir)*lap_perp_U(ir))
      end do

      S_phi_dphi_eq_nm=0.0d0
      S_phi_q_nm=0.0d0
      S_phi_dj_eq_nm=0.0d0
      S_phi_x_p_nm=0.0d0
      S_phi_mu_nm=0.0d0
      do ir=ir_min+1,ir_max
        S_phi_dphi_eq_nm=S_phi_dphi_eq_nm+dr/2.0d0&
        *(S_phi_dphi_eq_nm_den(ir-1)*r(ir-1)&
        +S_phi_dphi_eq_nm_den(ir)*r(ir))
        S_phi_q_nm=S_phi_q_nm+dr/2.0d0&
        *(S_phi_q_nm_den(ir-1)*r(ir-1)&
        +S_phi_q_nm_den(ir)*r(ir))
        S_phi_dj_eq_nm=S_phi_dj_eq_nm+dr/2.0d0&
        *(S_phi_dj_eq_nm_den(ir-1)*r(ir-1)&
        +S_phi_dj_eq_nm_den(ir)*r(ir))
        S_phi_x_p_nm=S_phi_x_p_nm+dr/2.0d0&
        *(S_phi_x_p_nm_den(ir-1)*r(ir-1)&
        +S_phi_x_p_nm_den(ir)*r(ir))
        S_phi_mu_nm=S_phi_mu_nm+dr/2.0d0&
        *(S_phi_mu_nm_den(ir-1)*r(ir-1)&
        +S_phi_mu_nm_den(ir)*r(ir))
      end do

      S_phi_dphi_eq_n(n)=S_phi_dphi_eq_n(n)+S_phi_dphi_eq_nm
      S_phi_q_n(n)=S_phi_q_n(n)+S_phi_q_nm
      S_phi_dj_eq_n(n)=S_phi_dj_eq_n(n)+S_phi_dj_eq_nm
      S_phi_x_p_n(n)=S_phi_x_p_n(n)+S_phi_x_p_nm
      S_phi_mu_n(n)=S_phi_mu_n(n)+S_phi_mu_nm
    end do ! : n
  end do ! : m

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

  do m=m_min,m_max
    do n=sn,en
      do ir=ir_min,ir_max
        phi(ir)=phi_nm(ir,n,m)
        bra_phi_U_non(ir)=bra_phi_U_non_nm(ir,n,m)
        nab_para_j_non(ir)=nab_para_j_non_nm(ir,n,m)
      end do

      call conjugate(phi,phi_cc)

      do ir=ir_min,ir_max
        S_phi_phi_U_nm_den(ir)=dble(phi_cc(ir)*bra_phi_U_non(ir))
        S_phi_psi_j_nm_den(ir)=-dble(phi_cc(ir)*nab_para_j_non(ir))
      end do

      S_phi_phi_U_nm=0.0d0
      S_phi_psi_j_nm=0.0d0
      do ir=ir_min+1,ir_max
        S_phi_phi_U_nm=S_phi_phi_U_nm+dr/2.0d0&
        *(S_phi_phi_U_nm_den(ir-1)*r(ir-1)&
        +S_phi_phi_U_nm_den(ir)*r(ir))
        S_phi_psi_j_nm=S_phi_psi_j_nm+dr/2.0d0&
        *(S_phi_psi_j_nm_den(ir-1)*r(ir-1)&
        +S_phi_psi_j_nm_den(ir)*r(ir))
      end do

      S_phi_phi_U_n(n)=S_phi_phi_U_n(n)+S_phi_phi_U_nm
      S_phi_psi_j_n(n)=S_phi_psi_j_n(n)+S_phi_psi_j_nm
    end do ! : n
  end do ! : m

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

  do m=m_min,m_max
    do n=sn,en
      do ir=ir_min,ir_max
        psi(ir)=psi_nm(ir,n,m)
        j(ir)=j_nm(ir,n,m)
        psi_old(ir)=psi_old_nm(ir,n,m)
        j_old(ir)=j_old_nm(ir,n,m)
      end do

      call conjugate(j,j_cc)
      call conjugate(j_old,j_old_cc)

      do ir=ir_min,ir_max
        dE_psi_nm_den(ir)=(dble(j_cc(ir)*psi(ir))&
        -dble(j_old_cc(ir)*psi_old(ir)))/(2.0d0*dt)
      end do

      dE_psi_nm=0.0d0
      do ir=ir_min+1,ir_max
        dE_psi_nm=dE_psi_nm+dr/2.0d0&
        *(dE_psi_nm_den(ir-1)*r(ir-1)&
        +dE_psi_nm_den(ir)*r(ir))
      end do

      dE_psi_n(n)=dE_psi_n(n)+dE_psi_nm
    end do ! : n
  end do ! : m

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

  do m=m_min,m_max
    do n=sn,en
      do ir=ir_min,ir_max
        phi(ir)=phi_nm(ir,n,m)
        psi(ir)=psi_nm(ir,n,m)
        j(ir)=j_nm(ir,n,m)
      end do

      call conjugate(j,j_cc)

      do ir=ir_min,ir_max
        S_j_q_nm_den(ir)=-dble((0.0d0,1.0d0)*k_para(ir,n,m)&
        *j_cc(ir)*phi(ir))
        S_j_dphi_eq_nm_den(ir)=-dble((0.0d0,1.0d0)*k_theta(ir,m)&
        *dphi_eq(ir)*j_cc(ir)*psi(ir))
        S_j_eta_nm_den(ir)=-eta*abs(j(ir))**2
      end do

      S_j_q_nm=0.0d0
      S_j_dphi_eq_nm=0.0d0
      S_j_eta_nm=0.0d0
      do ir=ir_min+1,ir_max
        S_j_q_nm=S_j_q_nm+dr/2.0d0&
        *(S_j_q_nm_den(ir-1)*r(ir-1)&
        +S_j_q_nm_den(ir)*r(ir))
        S_j_dphi_eq_nm=S_j_dphi_eq_nm+dr/2.0d0&
        *(S_j_dphi_eq_nm_den(ir-1)*r(ir-1)&
        +S_j_dphi_eq_nm_den(ir)*r(ir))
        S_j_eta_nm=S_j_eta_nm+dr/2.0d0&
        *(S_j_eta_nm_den(ir-1)*r(ir-1)&
        +S_j_eta_nm_den(ir)*r(ir))
      end do

      S_j_q_n(n)=S_j_q_n(n)+S_j_q_nm
      S_j_dphi_eq_n(n)=S_j_dphi_eq_n(n)+S_j_dphi_eq_nm
      S_j_eta_n(n)=S_j_eta_n(n)+S_j_eta_nm
    end do ! : n
  end do ! : m

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

  do m=m_min,m_max
    do n=sn,en
      do ir=ir_min,ir_max
        j(ir)=j_nm(ir,n,m)
        nab_para_phi_non(ir)=nab_para_phi_non_nm(ir,n,m)
      end do

      call conjugate(j,j_cc)

      do ir=ir_min,ir_max
        S_j_psi_phi_nm_den(ir)=-dble(j_cc(ir)*nab_para_phi_non(ir))
      end do

      S_j_psi_phi_nm=0.0d0
      do ir=ir_min+1,ir_max
        S_j_psi_phi_nm=S_j_psi_phi_nm+dr/2.0d0&
        *(S_j_psi_phi_nm_den(ir-1)*r(ir-1)&
        +S_j_psi_phi_nm_den(ir)*r(ir))
      end do

      S_j_psi_phi_n(n)=S_j_psi_phi_n(n)+S_j_psi_phi_nm
    end do ! : n
  end do ! : m

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

  do m=m_min,m_max
    do n=sn,en
      do ir=ir_min,ir_max
        p(ir)=p_nm(ir,n,m)
        p_old(ir)=p_old_nm(ir,n,m)
      end do

      do ir=ir_min,ir_max
        dE_p_nm_den(ir)=(abs(p(ir))**2&
        -abs(p_old(ir))**2)/(4.0d0*beta*dt)
      end do

      dE_p_nm=0.0d0
      do ir=ir_min+1,ir_max
        dE_p_nm=dE_p_nm+dr/2.0d0&
        *(dE_p_nm_den(ir-1)*r(ir-1)&
        +dE_p_nm_den(ir)*r(ir))
      end do

      dE_p_n(n)=dE_p_n(n)+dE_p_nm
    end do ! : n
  end do ! : m

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

  do m=m_min,m_max
    do n=sn,en
      do ir=ir_min,ir_max
        phi(ir)=phi_nm(ir,n,m)
        p(ir)=p_nm(ir,n,m)
        bra_x_phi(ir)=bra_x_phi_nm(ir,n,m)
        lap_perp_p(ir)=lap_perp_p_nm(ir,n,m)
      end do

      call conjugate(p,p_cc)

      do ir=ir_min,ir_max
        S_p_dp_eq_nm_den(ir)=dble((0.0d0,1.0d0)*k_theta(ir,m)&
        *dp_eq(ir)*p_cc(ir)*phi(ir))/(2.0d0*beta)
        S_p_x_phi_nm_den(ir)=dble(p_cc(ir)*bra_x_phi(ir))
        S_p_chi_nm_den(ir)=chi/(2.0d0*beta)&
        *dble(p_cc(ir)*lap_perp_p(ir))
      end do

      S_p_dp_eq_nm=0.0d0
      S_p_x_phi_nm=0.0d0
      S_p_chi_nm=0.0d0
      do ir=ir_min+1,ir_max
        S_p_dp_eq_nm=S_p_dp_eq_nm+dr/2.0d0&
        *(S_p_dp_eq_nm_den(ir-1)*r(ir-1)&
        +S_p_dp_eq_nm_den(ir)*r(ir))
        S_p_x_phi_nm=S_p_x_phi_nm+dr/2.0d0&
        *(S_p_x_phi_nm_den(ir-1)*r(ir-1)&
        +S_p_x_phi_nm_den(ir)*r(ir))
        S_p_chi_nm=S_p_chi_nm+dr/2.0d0&
        *(S_p_chi_nm_den(ir-1)*r(ir-1)&
        +S_p_chi_nm_den(ir)*r(ir))
      end do

      S_p_dp_eq_n(n)=S_p_dp_eq_n(n)+S_p_dp_eq_nm
      S_p_x_phi_n(n)=S_p_x_phi_n(n)+S_p_x_phi_nm
      S_p_chi_n(n)=S_p_chi_n(n)+S_p_chi_nm
    end do ! : n
  end do ! : m

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

  do m=m_min,m_max
    do n=sn,en
      do ir=ir_min,ir_max
        p(ir)=p_nm(ir,n,m)
        bra_phi_p_non(ir)=bra_phi_p_non_nm(ir,n,m)
      end do

      call conjugate(p,p_cc)

      do ir=ir_min,ir_max
        S_p_phi_p_nm_den(ir)=-dble(p_cc(ir)*bra_phi_p_non(ir))/(2.0d0*beta)
      end do

      S_p_phi_p_nm=0.0d0
      do ir=ir_min+1,ir_max
        S_p_phi_p_nm=S_p_phi_p_nm+dr/2.0d0&
        *(S_p_phi_p_nm_den(ir-1)*r(ir-1)&
        +S_p_phi_p_nm_den(ir)*r(ir))
      end do

      S_p_phi_p_n(n)=S_p_phi_p_n(n)+S_p_phi_p_nm
    end do ! : n
  end do ! : m

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

  do m=m_min,m_max
    do n=sn,en
      do ir=ir_min,ir_max
        U(ir)=U_nm(ir,n,m)
        phi(ir)=phi_nm(ir,n,m)
      end do

      call conjugate(phi,phi_cc)

      do ir=ir_min,ir_max
        E_phi_nm_den(ir)=-dble(phi_cc(ir)*U(ir))/2.0d0
      end do

      E_phi_nm(n,m)=0.0d0
      do ir=ir_min+1,ir_max
        E_phi_nm(n,m)=E_phi_nm(n,m)+dr/2.0d0&
        *(E_phi_nm_den(ir-1)*r(ir-1)&
        +E_phi_nm_den(ir)*r(ir))
      end do

      E_phi_n(n)=E_phi_n(n)+E_phi_nm(n,m)
    end do ! : n
  end do ! : m

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

  do m=m_min,m_max
    do n=sn,en
      do ir=ir_min,ir_max
        psi(ir)=psi_nm(ir,n,m)
        j(ir)=j_nm(ir,n,m)
      end do

      call conjugate(j,j_cc)

      do ir=ir_min,ir_max
        E_psi_nm_den(ir)=dble(j_cc(ir)*psi(ir))/2.0d0
      end do

      E_psi_nm(n,m)=0.0d0
      do ir=ir_min+1,ir_max
        E_psi_nm(n,m)=E_psi_nm(n,m)+dr/2.0d0&
        *(E_psi_nm_den(ir-1)*r(ir-1)&
        +E_psi_nm_den(ir)*r(ir))
      end do

      E_psi_n(n)=E_psi_n(n)+E_psi_nm(n,m)
    end do ! : n
  end do ! : m

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

  do m=m_min,m_max
    do n=sn,en
      do ir=ir_min,ir_max
        p(ir)=p_nm(ir,n,m)
      end do

      do ir=ir_min,ir_max
        E_p_nm_den(ir)=abs(p(ir))**2/(4.0d0*beta)
      end do

      E_p_nm(n,m)=0.0d0
      do ir=ir_min+1,ir_max
        E_p_nm(n,m)=E_p_nm(n,m)+dr/2.0d0&
        *(E_p_nm_den(ir-1)*r(ir-1)&
        +E_p_nm_den(ir)*r(ir))
      end do

      E_p_n(n)=E_p_n(n)+E_p_nm(n,m)
    end do ! : n
  end do ! : m

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

    do m=m_min,m_max
      do n=sn,en
        a_send_gather(n)=E_phi_nm(n,m)
      end do
      call mpi_allgatherv(a_send_gather,scount_gather,mpi_double_complex,&
      a_recv_gather,rcounts_gather,rdispls_gather,mpi_double_complex,&
      mpi_comm_world,ierr)
      do n=0,n_max
        E_phi_nm_gather(n,m)=a_recv_gather(n)
      end do

      do n=sn,en
        a_send_gather(n)=E_psi_nm(n,m)
      end do
      call mpi_allgatherv(a_send_gather,scount_gather,mpi_double_complex,&
      a_recv_gather,rcounts_gather,rdispls_gather,mpi_double_complex,&
      mpi_comm_world,ierr)
      do n=0,n_max
        E_psi_nm_gather(n,m)=a_recv_gather(n)
      end do

      do n=sn,en
        a_send_gather(n)=E_p_nm(n,m)
      end do
      call mpi_allgatherv(a_send_gather,scount_gather,mpi_double_complex,&
      a_recv_gather,rcounts_gather,rdispls_gather,mpi_double_complex,&
      mpi_comm_world,ierr)
      do n=0,n_max
        E_p_nm_gather(n,m)=a_recv_gather(n)
      end do
    end do ! : m

    if(myid==0)then
      write(filename,'("../data/energy_spec_nm_",i3.3,".dat")')it/it_skip
      write(6,*)filename
      open(unit=10,file=filename,form='formatted')
      write(10,fmt='(5a16)')'n','m','E_phi_nm','E_psi_nm','E_p_nm'
      do m=m_min,m_max
        do n=n_max,1,-1
          write(10,fmt='(5e16.4)')-dble(n),dble(m),&
          E_phi_nm_gather(n,-m),E_psi_nm_gather(n,-m),E_p_nm_gather(n,-m)
        end do
        do n=0,n_max
          write(10,fmt='(5e16.4)')dble(n),dble(m),&
          E_phi_nm_gather(n,m),E_psi_nm_gather(n,m),E_p_nm_gather(n,m)
        end do
        write(10,*)''
      end do
      close(10)
    end if
  end if ! : mod(it,it_skip)==0

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

    do m=m_min,m_max
      do n=sn,en
        do ir=ir_min,ir_max
          U(ir)=U_nm(ir,n,m)
          phi(ir)=phi_nm(ir,n,m)
          U_old(ir)=U_old_nm(ir,n,m)
          phi_old(ir)=phi_old_nm(ir,n,m)
        end do

        call conjugate(phi,phi_cc)
        call conjugate(phi_old,phi_old_cc)

        do ir=ir_min,ir_max
          E_phi_nm_den(ir)=-dble(phi_cc(ir)*U(ir))/2.0d0
          E_phi_old_nm_den(ir)=-dble(phi_old_cc(ir)*U_old(ir))/2.0d0
        end do

        E_phi_nm(n,m)=0.0d0
        E_phi_old_nm=0.0d0
        do ir=ir_min+1,ir_max
          E_phi_nm(n,m)=E_phi_nm(n,m)+dr/2.0d0&
          *(E_phi_nm_den(ir-1)*r(ir-1)&
          +E_phi_nm_den(ir)*r(ir))
          E_phi_old_nm=E_phi_old_nm+dr/2.0d0&
          *(E_phi_old_nm_den(ir-1)*r(ir-1)&
          +E_phi_old_nm_den(ir)*r(ir))
        end do

        E_phi_n(n)=E_phi_n(n)+E_phi_nm(n,m)
        E_phi_old_n(n)=E_phi_old_n(n)+E_phi_old_nm
      end do ! : n
    end do ! : m

    do n=sn,en
      int_E_phi_n(n)=int_E_phi_n(n)&
      +dt/2.0d0*(E_phi_n(n)+E_phi_old_n(n))
    end do

    ! . magnetic energy
    do n=sn,en
      E_psi_n(n)=0.0d0
      E_psi_old_n(n)=0.0d0
    end do

    do m=m_min,m_max
      do n=sn,en
        do ir=ir_min,ir_max
          psi(ir)=psi_nm(ir,n,m)
          j(ir)=j_nm(ir,n,m)
          psi_old(ir)=psi_old_nm(ir,n,m)
          j_old(ir)=j_old_nm(ir,n,m)
        end do

        call conjugate(j,j_cc)
        call conjugate(j_old,j_old_cc)

        do ir=ir_min,ir_max
          E_psi_nm_den(ir)=dble(j_cc(ir)*psi(ir))/2.0d0
          E_psi_old_nm_den(ir)=dble(j_old_cc(ir)*psi_old(ir))/2.0d0
        end do

        E_psi_nm(n,m)=0.0d0
        E_psi_old_nm=0.0d0
        do ir=ir_min+1,ir_max
          E_psi_nm(n,m)=E_psi_nm(n,m)+dr/2.0d0&
          *(E_psi_nm_den(ir-1)*r(ir-1)&
          +E_psi_nm_den(ir)*r(ir))
          E_psi_old_nm=E_psi_old_nm+dr/2.0d0&
          *(E_psi_old_nm_den(ir-1)*r(ir-1)&
          +E_psi_old_nm_den(ir)*r(ir))
        end do

        E_psi_n(n)=E_psi_n(n)+E_psi_nm(n,m)
        E_psi_old_n(n)=E_psi_old_n(n)+E_psi_old_nm
      end do ! : n
    end do ! : m

    do n=sn,en
      int_E_psi_n(n)=int_E_psi_n(n)&
      +dt/2.0d0*(E_psi_n(n)+E_psi_old_n(n))
    end do

    ! . pressure energy
    do n=sn,en
      E_p_n(n)=0.0d0
      E_p_old_n(n)=0.0d0
    end do

    do m=m_min,m_max
      do n=sn,en
        do ir=ir_min,ir_max
          p(ir)=p_nm(ir,n,m)
          p_old(ir)=p_old_nm(ir,n,m)
        end do

        do ir=ir_min,ir_max
          E_p_nm_den(ir)=abs(p(ir))**2/(4.0d0*beta)
          E_p_old_nm_den(ir)=abs(p_old(ir))**2/(4.0d0*beta)
        end do

        E_p_nm(n,m)=0.0d0
        E_p_old_nm=0.0d0
        do ir=ir_min+1,ir_max
          E_p_nm(n,m)=E_p_nm(n,m)+dr/2.0d0&
          *(E_p_nm_den(ir-1)*r(ir-1)&
          +E_p_nm_den(ir)*r(ir))
          E_p_old_nm=E_p_old_nm+dr/2.0d0&
          *(E_p_old_nm_den(ir-1)*r(ir-1)&
          +E_p_old_nm_den(ir)*r(ir))
        end do

        E_p_n(n)=E_p_n(n)+E_p_nm(n,m)
        E_p_old_n(n)=E_p_old_n(n)+E_p_old_nm
      end do ! : n
    end do ! : m

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
subroutine conjugate(a,a_cc)
! ------------------------------------------------------------------------------
  use coordinate
  implicit none
  complex(kind(0d0))::a(ir_min:ir_max)
  complex(kind(0d0))::a_cc(ir_min:ir_max)

  do ir=ir_min,ir_max
    a_cc(ir)=dble(a(ir))&
    +(0.0d0,1.0d0)*dble((0.0d0,1.0d0)*a(ir))
  end do

end subroutine conjugate
