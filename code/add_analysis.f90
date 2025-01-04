! ------------------------------------------------------------------------------
! . last update @20240519
! ------------------------------------------------------------------------------
program main
  implicit none

  ! call energy_spec_nm_py
  ! call growth_rate
  call p_contour_py
  call avg_energy_spec

end program main
!-------------------------------------------------------------------------------
subroutine energy_spec_nm_py
!-------------------------------------------------------------------------------
  implicit none
  ! ..
  integer::it
  integer,parameter::it_start=0
  integer,parameter::it_end=50
  character::filename*100
  ! ..
  integer::n,m
  integer,parameter::n_max=40,m_max=170
  integer,parameter::n_min=-n_max,m_min=-m_max
  double precision::n_dummy,m_dummy
  ! ..
  double precision::E_phi_nm(n_min:n_max,m_min:m_max)
  double precision::E_psi_nm(n_min:n_max,m_min:m_max)
  double precision::E_p_nm(n_min:n_max,m_min:m_max)

  do it=it_start,it_end
    write(filename,'("../data/energy_spec_nm_",i3.3,".dat")')it
    open(10,file=filename,form='formatted')
    read(10,*)
    do m=m_min,m_max
      do n=n_min,n_max
        read(10,fmt='(5e16.4)')n_dummy,m_dummy,&
        E_phi_nm(n,m),E_psi_nm(n,m),E_p_nm(n,m)
      end do
      read(10,*)
    end do
    close(10)
    ! ..
    write(filename,'("./data/energy_spec_nm_py_",i3.3,".dat")')it
    write(6,*)filename
    open(10,file=filename,form='formatted')
    write(10,fmt='(5a16)')'n','m','E_phi_nm','E_psi_nm','E_p_nm'
    do m=m_min,m_max
      do n=n_min,n_max
        write(10,fmt='(5e16.4)')dble(m),dble(n),&
        E_phi_nm(n,m),E_psi_nm(n,m),E_p_nm(n,m)
      end do
      ! write(10,*)''
    end do
    close(10)
  end do

end subroutine energy_spec_nm_py
!-------------------------------------------------------------------------------
subroutine growth_rate
!-------------------------------------------------------------------------------
  implicit none
  ! ..
  integer::it
  integer,parameter::it_start=35
  integer,parameter::it_end=50
  character::filename*100
  ! ..
  integer::n
  integer,parameter::n_max=40
  double precision::n_dummy
  ! ..
  double precision::E_phi_n
  double precision::E_psi_n
  double precision::E_p_n_t_start(0:n_max)
  double precision::E_p_n_t_end(0:n_max)

  it=it_start
  write(filename,'("../data/energy_spec_",i3.3,".dat")')it
  open(10,file=filename,form='formatted')
  read(10,*)
  do n=0,n_max
    read(10,fmt='(4e16.4)')n_dummy,&
    E_phi_n,E_psi_n,E_p_n_t_start(n)
  end do
  close(10)
  it=it_end
  write(filename,'("../data/energy_spec_",i3.3,".dat")')it
  open(10,file=filename,form='formatted')
  read(10,*)
  do n=0,n_max
    read(10,fmt='(4e16.4)')n_dummy,&
    E_phi_n,E_psi_n,E_p_n_t_end(n)
  end do
  close(10)
  ! ..
  write(filename,'("./data/growth_rate_",i3.3,"_",i3.3,".dat")')it_start,it_end
  write(6,*)filename
  open(10,file=filename,form='formatted')
  write(10,fmt='(2a16)')'n','growth_rate'
  do n=0,n_max
    write(10,fmt='(2e16.4)')dble(n),&
    log(E_p_n_t_end(n)/E_p_n_t_start(n))&
    /(2.0d0*10.0d0*dble(it_end-it_start))
  end do
  close(10)

end subroutine growth_rate
!-------------------------------------------------------------------------------
subroutine p_contour_py
!-------------------------------------------------------------------------------
  implicit none
  ! ..
  integer::it
  integer,parameter::it_start=0
  integer,parameter::it_end=50
  character::filename*100
  ! ..
  integer::ir
  integer,parameter::ir_min=0
  integer,parameter::ir_max=1024
  double precision::r(ir_min:ir_max)
  double precision,parameter::r_min=0.1d0
  double precision,parameter::r_max=1.0d0
  double precision,parameter::dr=(r_max-r_min)/dble(ir_max-ir_min)
  ! ..
  double precision,parameter::pi=3.14159265358979323846d0
  integer::itheta
  integer,parameter::itheta_min=0
  integer,parameter::itheta_max=4096
  integer,parameter::numdata_theta=itheta_max-itheta_min
  double precision::theta(0:numdata_theta)
  double precision,parameter::inumdata_theta=1.0d0/dble(numdata_theta)
  ! ..
  double precision::x
  double precision::Z
  double precision::p_real(ir_min:ir_max,itheta_min:itheta_max)

  do ir=ir_min,ir_max
    r(ir)=r_min+dble(ir-ir_min)*dr
  end do
  do itheta=0,numdata_theta
    theta(itheta)=-pi+2.0d0*pi*dble(itheta)*inumdata_theta
  end do
  ! ..
  do it=it_start,it_end
    write(filename,'("../data/p_contour_",i3.3,".dat")')it
    open(10,file=filename,form='formatted')
    read(10,*)
    do itheta=0,numdata_theta
      do ir=ir_min,ir_max
        read(10,fmt='(3e16.4)')x,Z,p_real(ir,itheta)
      end do
      read(10,*)
    end do
    close(10)
    ! ..
    write(filename,'("./data/p_contour_py_",i3.3,".dat")')it
    write(6,*)filename
    open(10,file=filename,form='formatted')
    write(10,fmt='(3a16)')'r','theta','p_real'
    do itheta=0,numdata_theta
      do ir=ir_min,ir_max
        write(10,fmt='(3e16.4)')r(ir),theta(itheta),p_real(ir,itheta)
      end do
      !write(10,*)''
    end do
    close(10)
  end do

end subroutine p_contour_py
!-------------------------------------------------------------------------------
subroutine avg_energy_spec
!-------------------------------------------------------------------------------
  implicit none
  ! ..
  integer::it
  integer,parameter::it_start=35
  integer,parameter::it_end=50
  character::filename*100
  ! ..
  integer::n
  integer,parameter::n_max=320
  double precision::n_dummy
  ! ..
  double precision::int_E_phi_n_t_start(0:n_max)
  double precision::int_E_phi_n_t_end(0:n_max)
  double precision::int_E_psi_n_t_start(0:n_max)
  double precision::int_E_psi_n_t_end(0:n_max)
  double precision::int_E_p_n_t_start(0:n_max)
  double precision::int_E_p_n_t_end(0:n_max)

  it=it_start
  write(filename,'("../data/int_energy_spec_",i3.3,".dat")')it
  open(10,file=filename,form='formatted')
  read(10,*)
  do n=0,n_max
    read(10,fmt='(4e16.4)')n_dummy,&
    int_E_phi_n_t_start(n),int_E_psi_n_t_start(n),int_E_p_n_t_start(n)
  end do
  close(10)
  it=it_end
  write(filename,'("../data/int_energy_spec_",i3.3,".dat")')it
  open(10,file=filename,form='formatted')
  read(10,*)
  do n=0,n_max
    read(10,fmt='(4e16.4)')n_dummy,&
    int_E_phi_n_t_end(n),int_E_psi_n_t_end(n),int_E_p_n_t_end(n)
  end do
  close(10)
  ! ..
  write(filename,'("./data/avg_energy_spec_",i3.3,"_",i3.3,".dat")')it_start,it_end
  write(6,*)filename
  open(10,file=filename,form='formatted')
  write(10,fmt='(4a16)')'n','E_phi_n_avg','E_psi_n_avg','E_p_n_avg'
  do n=0,n_max
    write(10,fmt='(4e16.4)')dble(n),&
    (int_E_phi_n_t_end(n)-int_E_phi_n_t_start(n))&
    /(10.0d0*dble(it_end-it_start)),&
    (int_E_psi_n_t_end(n)-int_E_psi_n_t_start(n))&
    /(10.0d0*dble(it_end-it_start)),&
    (int_E_p_n_t_end(n)-int_E_p_n_t_start(n))&
    /(10.0d0*dble(it_end-it_start))
  end do
  close(10)

end subroutine avg_energy_spec