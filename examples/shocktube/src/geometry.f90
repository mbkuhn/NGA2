!> Various definitions and tools for initializing NGA2 config
module geometry
   use config_class, only: config
   use precision,    only: WP
   implicit none
   private
   
   !> Single config
   type(config), public :: cfg
   
   public :: geometry_init
   
contains
   
   
   !> Initialization of problem geometry
   subroutine geometry_init
      use sgrid_class, only: sgrid
      use param,       only: param_read
      implicit none
      type(sgrid) :: grid
      
      
      ! Create a grid from input params
      create_grid: block
         use sgrid_class, only: cartesian
         use messager,    only: die
         use string,      only: str_medium
         integer :: i,j,k,nx,ny,nz
         real(WP) :: Lx,Ly,Lz,x0
         real(WP), dimension(:), allocatable :: x,y,z
         character(len=str_medium) :: shtb_type
         
         ! Read in grid definition
         call param_read('nx',nx); allocate(x(nx+1))
         call param_read('Shocktube type',shtb_type)
         select case (trim(adjustl(shtb_type)))
         case ('sod')
            x0 = 0.0_WP
            Lx = 1.0_WP
         case ('shu-osher')
            x0 = -5.0_WP
            Lx = 10.0_WP
         case default
            call die("Unknown case for shocktube")
         end select

         ! 1D problem
         ny = 1; allocate(y(ny+1)); Ly = Lx/nx;
         nz = 1; allocate(z(nz+1)); Lz = Lx/nx;
         
         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=x0+real(i-1,WP)/real(nx,WP)*Lx
         end do
         do j=1,ny+1
            y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.false.,yper=.true.,zper=.true.,name='ShockTube')
         
      end block create_grid
      
      
      ! Create a config from that grid on our entire group
      create_cfg: block
         use parallel, only: group
         integer, dimension(3) :: partition
         
         ! Read in partition
         call param_read('Partition',partition,short='p')
         
         ! Create partitioned grid
         cfg=config(grp=group,decomp=partition,grid=grid)
         
      end block create_cfg
      
      
   end subroutine geometry_init
   
   
end module geometry
