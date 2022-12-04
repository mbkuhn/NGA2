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
   
   public :: simulation_init,simulation_run,simulation_final
   
contains
    
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max time',time%tmax)
         !call param_read('Max steps',time%nmax)
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
         use mathtools,  only: Pi
         use messager,   only: die
         integer :: i,j,k,n
         real(WP), dimension(3) :: xyz
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
         fs%Ui = 1.0_WP; fs%Vi = 0.0_WP; fs%Wi = 0.0_WP
         call fs%interp_vel_basic(vf,fs%Ui,fs%Vi,fs%Wi,fs%U,fs%V,fs%W)

         r = 0.2_WP
         int_dist = 0.5_WP * r
         do k=fs%cfg%kmino_,fs%cfg%kmaxo_
            do j=fs%cfg%jmino_,fs%cfg%jmaxo_
               do i=fs%cfg%imino_,fs%cfg%imaxo_
                  d = sqrt(fs%cfg%xm(i)**2+fs%cfg%ym(j)**2+fs%cfg%zm(k)**2)
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
      implicit none
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
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

         ! Set velocity for current time
         fs%Ui = 1.0_WP; fs%Vi = 0.0_WP; fs%Wi = 0.0_WP
         ! Interpolate for face velocity
         call fs%interp_vel_basic(vf,fs%Ui,fs%Vi,fs%Wi,fs%U,fs%V,fs%W)

         ! Perform sub-iterations
         do while (time%it.le.time%itmax)
            
            ! Predictor step, involving advection and pressure terms
            call fs%advection_step(time%dt,vf,matmod)
            
            ! Record convergence monitor
            call cvgfile%write()
            ! Increment sub-iteration counter
            time%it=time%it+1
            
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
      
      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! bcond
      ! timetracker
      
      ! Deallocate work arrays - none
      
   end subroutine simulation_final
   
end module simulation
