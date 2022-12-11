!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use mast_class,        only: mast
   use vfs_class,         only: vfs
   use matm_class,        only: matm
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   use string,            only: str_medium
   use mathtools,  only: Pi
   implicit none
   private
   
   !> Single two-phase flow solver, volume fraction solver, and material model set
   !> With corresponding time tracker
   type(mast),        public :: fs
   type(vfs),         public :: vf
   type(matm),        public :: matmod
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   type(ensight) :: ens_out
   type(event)   :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,cvgfile

   !> Case specific
   real(WP), dimension(3) :: center
   real(WP) :: r, int_dist
   
   public :: simulation_init,simulation_run,simulation_final
   
contains
    
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read, param_exists
      implicit none
      
      
      ! Initialize time tracker
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max steps',time%nmax)
         if (cfg%nz.eq.1) then
            time%tmax = 8.0_WP
         else
            time%tmax = 3.0_WP
         end if
         time%dt=time%dtmax
         time%itmax=1
      end block initialize_timetracker
      
      
      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use mms_geom, only: cube_refine_vol
         use vfs_class, only: r2p,lvira,VFhi,VFlo
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver with lvira reconstruction
         vf=vfs(cfg=cfg,reconstruction_method=lvira,name='VOF')
         ! Single-phase simulation; all gas
         vf%VF = 0.0_WP
         ! Initialize barycenters
         do k=vf%cfg%kmino_,vf%cfg%kmaxo_
            do j=vf%cfg%jmino_,vf%cfg%jmaxo_
               do i=vf%cfg%imino_,vf%cfg%imaxo_
                  vf%Lbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                  vf%Gbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
               end do
            end do
         end do
         ! Update the band
         call vf%update_band()
         ! Perform interface reconstruction from VOF field
         call vf%build_interface()
         ! Set interface at the boundaries
         call vf%set_full_bcond()
         ! Create discontinuous polygon mesh from IRL interface
         call vf%polygonalize_interface()
         ! Calculate distance from polygons
         call vf%distance_from_polygon()
         ! Calculate subcell phasic volumes
         call vf%subcell_vol()
         ! Calculate curvature
         call vf%get_curvature()
         ! Reset moments to guarantee compatibility with interface reconstruction
         call vf%reset_volume_moments()
      end block create_and_initialize_vof
      
      
      ! Create a compressible two-phase flow solver
      create_and_initialize_flow_solver: block
         use matm_class, only: none
         use ils_class,  only: gmres_smg,pcg_bbox
         use messager,   only: die
         integer :: i,j,k,n
         real(WP) :: gamm_g,dens,v0,d,int_dist,r
         ! Create material model class
         matmod=matm(cfg=cfg,name='Liquid-gas models')
         ! Get EOS parameters from input
         call param_read('Gas gamma',gamm_g)
         ! Register equations of state
         call matmod%register_idealgas('gas',gamm_g)
         ! Create flow solver
         fs=mast(cfg=cfg,name='Two-phase All-Mach',vf=vf)
         ! Register flow solver variables with material models
         call matmod%register_thermoflow_variables('liquid',fs%Lrho,fs%Ui,fs%Vi,fs%Wi,fs%LrhoE,fs%LP)
         call matmod%register_thermoflow_variables('gas'   ,fs%Grho,fs%Ui,fs%Vi,fs%Wi,fs%GrhoE,fs%GP)
         ! As a shocktube case, it is intended to be inviscid, no diffusion
         call matmod%register_diffusion_thermo_models(viscmodel_liquid=none,viscmodel_gas=none,hdffmodel_liquid=none,hdffmodel_gas=none)
         ! Surface tension should be set to 0
         fs%sigma = 0.0_WP
         ! Configure pressure solver
         call param_read('Pressure iteration',fs%psolv%maxit)
         call param_read('Pressure tolerance',fs%psolv%rcvg)
         ! Configure implicit momentum solver
         call param_read('Implicit iteration',fs%implicit%maxit)
         call param_read('Implicit tolerance',fs%implicit%rcvg)
         ! Setup the solver
         call fs%setup(pressure_ils=pcg_bbox,implicit_ils=gmres_smg)

         ! Initial conditions
         fs%Lrho = 1.0_WP; fs%LrhoE = 1.0_WP; 
         do k = fs%cfg%kmino_,fs%cfg%kmaxo_
            do j = fs%cfg%jmino_,fs%cfg%jmaxo_
               do i = fs%cfg%imino_,fs%cfg%imaxo
                  if (fs%cfg%nz.eq.1) then
                     fs%Ui(i,j,k) =-2.0_WP*sin(Pi*fs%cfg%xm(i))**2*sin(Pi*fs%cfg%ym(j))*cos(Pi*fs%cfg%ym(j))
                     fs%Vi(i,j,k) = 2.0_WP*sin(Pi*fs%cfg%ym(j))**2*sin(Pi*fs%cfg%xm(i))*cos(Pi*fs%cfg%xm(i))
                     fs%Wi(i,j,k) = 0.0_WP
                     fs%U(i,j,k) =-2.0_WP*sin(Pi*fs%cfg%x(i))**2*sin(Pi*fs%cfg%ym(j))*cos(Pi*fs%cfg%ym(j))
                     fs%V(i,j,k) = 2.0_WP*sin(Pi*fs%cfg%y(j))**2*sin(Pi*fs%cfg%xm(i))*cos(Pi*fs%cfg%xm(i))
                     fs%W(i,j,k) = 0.0_WP
                  else
                     fs%Ui(i,j,k) =+2.0_WP*sin(Pi*fs%cfg%xm(i))**2*sin(2.0_WP*Pi*fs%cfg%ym(j))*sin(2.0_WP*Pi*fs%cfg%zm(k))
                     fs%Vi(i,j,k) = -sin(2.0_WP*Pi*fs%cfg%xm(i))*sin(Pi*fs%cfg%ym(j))**2*sin(2.0_WP*Pi*fs%cfg%zm(k))
                     fs%Wi(i,j,k) = -sin(2.0_WP*Pi*fs%cfg%xm(i))*sin(2.0_WP*Pi*fs%cfg%ym(j))*sin(Pi*fs%cfg%zm(k))**2
                     fs%U(i,j,k) =+2.0_WP*sin(Pi*fs%cfg%x(i))**2*sin(2.0_WP*Pi*fs%cfg%ym(j))*sin(2.0_WP*Pi*fs%cfg%zm(k))
                     fs%V(i,j,k) = -sin(2.0_WP*Pi*fs%cfg%xm(i))*sin(Pi*fs%cfg%y(j))**2*sin(2.0_WP*Pi*fs%cfg%zm(k))
                     fs%W(i,j,k) = -sin(2.0_WP*Pi*fs%cfg%xm(i))*sin(2.0_WP*Pi*fs%cfg%ym(j))*sin(Pi*fs%cfg%z(k))**2
                  end if
               end do
            end do
         end do
         if (fs%cfg%nz.eq.1) then
            center = (/0.5_WP, 0.75_WP, fs%cfg%zm(fs%cfg%kmin_)/)
         else
            center = (/0.35_WP, 0.35_WP, 0.35_WP/)
         end if

         ! Default values
         r = 0.2_WP
         ! Get from input if present
         if (param_exists('Signal radius')) then
            call param_read('Signal radius', r)
         end if
         int_dist = 0.5_WP * r
         if (param_exists('Signal edge thickness')) then
            call param_read('Signal edge thickness', int_dist)
         end if
         do k=fs%cfg%kmino_,fs%cfg%kmaxo_
            do j=fs%cfg%jmino_,fs%cfg%jmaxo_
               do i=fs%cfg%imino_,fs%cfg%imaxo_
                  d = sqrt((fs%cfg%xm(i)-center(1))**2+(fs%cfg%ym(j)-center(2))**2+(fs%cfg%zm(k)-center(3))**2)
                  if (d.lt.r-int_dist) then
                     fs%Grho (i,j,k)=2.0_WP
                  else if (d.gt.r+int_dist) then
                     fs%Grho (i,j,k)=1.0_WP
                  else
                     fs%Grho (i,j,k)=1.5_WP-0.5_WP*sin((d-r)/(2.0_WP*int_dist)*Pi)
                  end if
               end do
            end do
         end do

         fs%RHO = fs%Grho

      end block create_and_initialize_flow_solver
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='SignalAdvect')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',fs%Ui,fs%Vi,fs%Wi)
         call ens_out%add_scalar('Density',fs%RHO)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         call vf%get_max()
         ! Create simulation monitor
         mfile=monitor(fs%cfg%amRoot,'simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%write()
         ! Create convergence monitor
         cvgfile=monitor(fs%cfg%amRoot,'cvg')
         call cvgfile%add_column(time%n,'Timestep number')
         call cvgfile%add_column(time%it,'Iteration')
         call cvgfile%add_column(time%t,'Time')
         call cvgfile%add_column(fs%impl_it_x,'Impl_x iteration')
         call cvgfile%add_column(fs%impl_rerr_x,'Impl_x error')
         call cvgfile%add_column(fs%impl_it_y,'Impl_y iteration')
         call cvgfile%add_column(fs%impl_rerr_y,'Impl_y error')
         call cvgfile%add_column(fs%implicit%it,'Impl_z iteration')
         call cvgfile%add_column(fs%implicit%rerr,'Impl_z error')
         call cvgfile%add_column(fs%psolv%it,'Pressure iteration')
         call cvgfile%add_column(fs%psolv%rerr,'Pressure error')
         ! Create CFL monitor
         cflfile=monitor(fs%cfg%amRoot,'cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(fs%CFLst,'STension CFL')
         call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
         call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
         call cflfile%write()
      end block create_monitor
      
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation - this mimicks NGA's old time integration for multiphase
   subroutine simulation_run
      use mast_class, only: central4th
      use mathtools,  only: Pi
      implicit none
      integer :: i,j,k
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         ! Modify time to get correct final time
         time%dt = min(time%dt,time%tmax-time%t)
         call time%increment()
         
         ! Reinitialize phase pressure by syncing it with conserved phase energy
         call fs%reinit_phase_pressure(vf,matmod)
         fs%Uiold=fs%Ui; fs%Viold=fs%Vi; fs%Wiold=fs%Wi
         fs%RHOold = fs%RHO
         ! Remember old flow variables (phase)
         fs%Grhoold = fs%Grho; fs%Lrhoold = fs%Lrho
         fs%GrhoEold=fs%GrhoE; fs%LrhoEold=fs%LrhoE
         fs%GPold   =   fs%GP; fs%LPold   =   fs%LP

         ! Remember old interface, including VF and barycenters
         call vf%copy_interface_to_old()
         
         ! Create in-cell reconstruction
         call fs%flow_reconstruct(vf,central4th)

         ! Zero pressure variables in predictor
         fs%P = 0.0_WP
         fs%Pjx = 0.0_WP; fs%Pjy = 0.0_WP; fs%Pjz = 0.0_WP
         fs%Hpjump = 0.0_WP

         ! Set flag as desired
         fs%sl_x = 1; fs%sl_y = 1; fs%sl_z = 1

         ! Set advective velocity for current time
         do k = fs%cfg%kmino_,fs%cfg%kmaxo_
            do j = fs%cfg%jmino_,fs%cfg%jmaxo_
               do i = fs%cfg%imino_,fs%cfg%imaxo
                  if (fs%cfg%nz.eq.1) then
                     fs%U(i,j,k) =-2.0_WP*sin(Pi*fs%cfg%x(i))**2*sin(Pi*fs%cfg%ym(j)) &
                                          *cos(Pi*fs%cfg%ym(j))*cos(Pi*(time%t + 0.5_WP*time%dt)/8.0_WP)
                     fs%V(i,j,k) = 2.0_WP*sin(Pi*fs%cfg%y(j))**2*sin(Pi*fs%cfg%xm(i)) &
                                          *cos(Pi*fs%cfg%xm(i))*cos(Pi*(time%t + 0.5_WP*time%dt)/8.0_WP) 
                     fs%W(i,j,k) = 0.0_WP
                  else
                     fs%U(i,j,k) =+2.0_WP*sin(Pi*fs%cfg%x(i))**2*sin(2.0_WP*Pi*fs%cfg%ym(j)) &
                                    *sin(2.0_WP*Pi*fs%cfg%zm(k))*cos(Pi*(time%t + 0.5_WP*time%dt)/3.0_WP)
                     fs%V(i,j,k) = -sin(2.0_WP*Pi*fs%cfg%xm(i))*sin(Pi*fs%cfg%y(j))**2 &
                                    *sin(2.0_WP*Pi*fs%cfg%zm(k))*cos(Pi*(time%t + 0.5_WP*time%dt)/3.0_WP)
                     fs%W(i,j,k) = -sin(2.0_WP*Pi*fs%cfg%xm(i))*sin(2.0_WP*Pi*fs%cfg%ym(j)) &
                                    *sin(Pi*fs%cfg%z(k))**2*cos(Pi*(time%t + 0.5_WP*time%dt)/3.0_WP)
                  end if
               end do
            end do
         end do
            
         ! Predictor step, involving advection and pressure terms
         call fs%advection_step(time%dt,vf,matmod)
            
         ! Record convergence monitor
         call cvgfile%write()

         ! Set velocity for visualization
         do k = fs%cfg%kmino_,fs%cfg%kmaxo_
            do j = fs%cfg%jmino_,fs%cfg%jmaxo_
               do i = fs%cfg%imino_,fs%cfg%imaxo
                  if (fs%cfg%nz.eq.1) then
                     fs%Ui(i,j,k) =-2.0_WP*sin(Pi*fs%cfg%xm(i))**2*sin(Pi*fs%cfg%ym(j)) &
                                          *cos(Pi*fs%cfg%ym(j))*cos(Pi*(time%t + 0.5_WP*time%dt)/8.0_WP)
                     fs%Vi(i,j,k) = 2.0_WP*sin(Pi*fs%cfg%ym(j))**2*sin(Pi*fs%cfg%xm(i)) &
                                          *cos(Pi*fs%cfg%xm(i))*cos(Pi*(time%t + 0.5_WP*time%dt)/8.0_WP) 
                     fs%Wi(i,j,k) = 0.0_WP
                  else
                     fs%Ui(i,j,k) =+2.0_WP*sin(Pi*fs%cfg%xm(i))**2*sin(2.0_WP*Pi*fs%cfg%ym(j)) &
                                    *sin(2.0_WP*Pi*fs%cfg%zm(k))*cos(Pi*(time%t + 0.5_WP*time%dt)/3.0_WP)
                     fs%Vi(i,j,k) = -sin(2.0_WP*Pi*fs%cfg%xm(i))*sin(Pi*fs%cfg%ym(j))**2 &
                                    *sin(2.0_WP*Pi*fs%cfg%zm(k))*cos(Pi*(time%t + 0.5_WP*time%dt)/3.0_WP)
                     fs%Wi(i,j,k) = -sin(2.0_WP*Pi*fs%cfg%xm(i))*sin(2.0_WP*Pi*fs%cfg%ym(j)) &
                                    *sin(Pi*fs%cfg%zm(k))**2*cos(Pi*(time%t + 0.5_WP*time%dt)/3.0_WP)
                  end if
               end do
            end do
         end do
         
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
         
         ! Perform and output monitoring
         call fs%get_max()
         call vf%get_max()
         call mfile%write()
         call cflfile%write()
         
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM
      use parallel,  only: MPI_REAL_WP
      implicit none
      real(WP) :: L2, d, buf1, ncell
      integer :: i,j,k,ierr
      
      ! Compare final density to initial density to calculate L2 error
      L2 = 0.0_WP
      ncell = 0.0_WP

      do k=fs%cfg%kmin_,fs%cfg%kmax_
         do j=fs%cfg%jmin_,fs%cfg%jmax_
            do i=fs%cfg%imin_,fs%cfg%imax_
               d = sqrt((fs%cfg%xm(i)-center(1))**2+(fs%cfg%ym(j)-center(2))**2+(fs%cfg%zm(k)-center(3))**2)
               if (d.lt.r-int_dist) then
                  L2 = L2 + (fs%Grho (i,j,k)-2.0_WP)**2
               else if (d.gt.r+int_dist) then
                  L2 = L2 + (fs%Grho (i,j,k)-1.0_WP)**2
               else
                  L2 = L2 + (fs%Grho (i,j,k)-(1.5_WP-0.5_WP*sin((d-r)/(2.0_WP*int_dist)*Pi)))**2
               end if
               ncell = ncell + 1.0_WP
            end do
         end do
      end do
      call MPI_ALLREDUCE(L2,buf1,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr);
      L2 = sqrt(buf1)
      call MPI_ALLREDUCE(ncell,buf1,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr);
      L2 = L2 / buf1
      if (fs%cfg%amRoot) print*,'L2 error per cell', L2
      
   end subroutine simulation_final
   
end module simulation
