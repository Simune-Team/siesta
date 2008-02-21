module mesh1D

!-----------------------------------------------------------------
! Contains utility routines to manipulate one-dimensional 
! functions defined in a mesh, and to solve some differential
! equations by Numerov's method
!-----------------------------------------------------------------
! get_n: returns the number of mesh points
! Usage:
!   n = get_n( xmin, xmax, dxmin, dxmax )
!     integer  n     ! Number of logarithmic-mesh points
!     real(DP) xmin  ! First mesh point
!     real(DP) xmax  ! Last mesh point
!     real(DP) dxmin ! First mesh interval: x(2)-x(1)
!     real(DP) dxmax ! Next to last mesh interval: x(n+1)-x(n)
! Notes:
! - All arguments are input and required
! - No mesh needs to have been previously set. Instead, this
!   function is to help to find the argument n needed to set it.
!-----------------------------------------------------------------
! set_mesh: sets a uniform, logarithmic, or arbitrary 1D mesh
! Usage:
!   call set_mesh( n, x, xmin, xmax, a, dxndx1 )
!     integer  n     : Number of mesh points
!     real(DP) x(n)  : Mesh point values
!     real(DP) xmin  : Value of first mesh point
!     real(DP) xmax  : Value of last mesh point
!     real(DP) a     : Fractional increment of succesive mesh 
!                      intervals of a logarithmic mesh:
!                        x(i) = b * (exp(a*(i-1)) - 1)
!     real(DP) dxndx1: Ratio between last and first intervals
!                      in a logarithmic mesh
! Notes:
! - All arguments are input
! - All the arguments, except n, are optional.
! - If x is present, no other optional arguments may be present
! - If x is not present, xmax must be present, and only one of
!   a or dxndx1 may be present
! - If xmin is not present, xmin=0 is assumed
! - If only xmax (and optionally xmin) is present, a regular mesh 
!   is generated
! - It is preferable to let set_mesh generate a logarithmic mesh,
!   and then retrieve it using get_mesh, than to generate it 
!   outside and call set_mesh with argument x. This is because
!   set_mesh can then use analytic formulas for dx/di, etc.
!   Then, numerov should be called without argument x.
! Examples:
! - To initialize a previously generated mesh:
!      call set_mesh( n, x=x_mesh(1:n) )
! - To generate a regular mesh:
!      call set_mesh( n, xmax=x_max )
! - To generate a logarithmic mesh:
!      call set_mesh( n, dxndx1=delta_x_max/delta_x_min )
! Warning:
! - One should not fix parameter a while increasing nx to 
!   convergence, because this may easily lead to an extremely  
!   small dx1 and to roundoff errors. Generally, it is safer
!   and easier to fix dxndx1, specially for this purpose.
!-----------------------------------------------------------------
! get_mesh: returns a previously-set 1D mesh
! Usage:
!   call get_mesh( n, nx, x, dxdi, d2xdi2, d3xdi3, d4xdi4 )
!     integer  n         : Size of array arguments (input)
!     integer  nx        : Number of mesh points (ouput)
!     real(DP) x(n)      : Mesh point values
!     real(DP) dxdi(n)   : dx/di
!     real(DP) d2xdi2(n) : d2x/di2
!     real(DP) d3xdi3(n) : d3x/di3
!     real(DP) d4xdi4(n) : d4x/di4
! Notes:
! - All arguments are output, except n
! - All the arguments, except n and nx, are optional.
! - If any optional array is present, nx is limited by its size, 
!   i.e. nx.le.n
! - If no optional array is present,  present, nx is the number 
!   of mesh points stored
! Examples:
! - To get the first derivatives:
!      call get_mesh( n, nx, dxdi=dxdi(1:n) )
!----------------------------------------------------------------
! derivative: returns the derivative of a function at the same
! mesh points
! Usage:
!   dydx = derivative( n, y, x, dx, order )
!     real(DP) derivative(n) : Return function derivative
!     integer  n             : Number of points
!     real(DP) y(n)          : Values of function to derivate
!     real(DP) x(n)          : Mesh values (optional)
!     real(DP) dx            : Mesh interval (optional)
!     integer  order         : order of the derivative (optional)
! Notes:
! - All arguments are input but the function is not pure because,
!   if x or dx are present, they change the default mesh.
! - The two optional arguments x and dx are incompatible. They
!   are used as in routine numerov. If none is present, the last
!   defined mesh is used.
! - If order is not present, order=1 is assumed.
! - If n is smaller than the number of mesh points stored, the
!   first n of them are used.
! Examples:
! - To find the Laplacian of y(x), using a previously defined mesh:
!      d2ydx2(1:n) = derivative( n, y(1:n), order=2 )
!----------------------------------------------------------------
! integral: returns the integral of a function defined in a mesh
! Usage:
!   the_integral = integral( n, y, x, dx )
!     real(DP) integral : Return value
!     integer  n        : Number of points
!     real(DP) y(n)     : Values of function to integrate
!     real(DP) x(n)     : Mesh values (optional)
!     real(DP) dx       : Mesh interval (optional)
! Notes:
! - All arguments are input but the function is not pure because,
!   if x or dx are present, they change the default mesh.
! - The two optional arguments x and dx are incompatible. They
!   are used as in routine numerov. If none is present, the last
!   defined mesh is used.
! - If n is smaller than the number of mesh points stored, the
!   first n of them are used.
! Examples:
! - To find the norm of psi, using a previously defined mesh:
!      norm = integral( n, y=psi(1:n)**2 )
!-----------------------------------------------------------------
! numerov: solves an ordinary differential equation of the form
!    d2y/dx2 = f(x,y) = f0(x) + f1(x)*y
! by the Numerov method:
!    (y(x+dx)-2*y(x)+y(x-dx))/dx2 = (f(x+dx)+10*f(x)+f(x-dx))/12
! from which y(x+dx) can be solved from y(x) and y(x-dx).
! Typical cases are f0=-4*pi*rho(x), f1=0 (Poisson eq.) 
!               and f0=0, f1(x)=2*(V(x)-E) (Schroedinger equation)
! Notice that f must not depend on y' (first derivative)
! Usage:
!   call numerov( n, y, yp, ypp, f0, f1, x, dx, norm )
!     integer  n     : Number of mesh points
!     real(DP) y(n)  : Solution
!     real(DP) yp(n) : First derivative y' = dy/dx
!     real(DP) ypp(n): Second derivative y'' = d2y/dx2 = f(x,y)
!     real(DP) f0(n) : Component of f(x,y) independent of y
!     real(DP) f1(n) : Componnet of f(x,y) proportional to y
!     real(DP) x(n)  : Mesh points
!     real(DP) dx    : Mesh interval (only x OR dx allowed)
!     real(DP) norm  : Norm for solution of homogeneous eqs.
! Notes:
! - All arguments, except y and ypp, are input.
! - All the arguments, except n, are optional.
! - If y is not present, only the mesh is initialized. If the
!   differential equation is solved repeatedly with the same mesh,
!   it is recommended to initialize it only once, making further
!   calls without x or dx present. In that case, the mesh used is
!   defined by the last call to any routine (set_mesh, numerov,
!   derivative, or integral) with a mesh-setting argument
! - If n is smaller than the number of mesh points stored, the
!   first n of them are used.
! - Only one of x OR dx are allowed to define the mesh
! - If f0 or f1 are not present, they are assumed to be zero.
! - The first two values of y must be given on input, to specify 
!   the boundary condition. If f0=0 (homogeneous equation), at
!   least one of them must be nonzero.
! - For f0 not present (homogeneous equation) the normalization
!   of the solution (integral of y**2) is fixed by argument norm.
! - If norm is not present, the norm of y is fixed by the input 
!   values of y(1) and y(2)
! - Output array yp is useful to normalize unbound states
! - Output array ypp is useful for a spline interpolation of y(x)
! Examples:
! - To initialize a variable radial mesh for inwards integration
!     call numerov( n+1, x=r_mesh(n:0:-1) )
! - To solve Poisson's equation with the previously defined mesh:
!     call numerov( n+1, y=rV(n:0:-1), f0=-4*pi*rrho(n:0:-1) )
!   with rV=r*V, rrho=r*rho. rV(n) and rV(n-1) may be initialized
!   to Q (integral of rho) to set the zero of potential at infty.
! - To integrate Schroedinger's equation with a regular mesh
!   (psi(1) and psi(2) required on input):
!     call numerov( n, y=psi(1:n), f1=V(1:n)-E, &
!                   dx=delta_x, norm=1.d0 )
!----------------------------------------------------------------
! Algorithms:
! The derivatives of x(n) are approximated by a 5-point formula:
! from   x(n+m) = x(n) + x'(n)*m + x''(n)*m^2/2 +
!                 x'''(n)*m^3/6 + x''''(n)*m^4/24
! we find
!   x'(n)   = ( x(n-2) - 8*x(n-1)          + 8*x(n+1) -x(n+2) )/12
!   x''(n)  = (-x(n-2) +16*x(n-1) -30*x(n) +16*x(n+1) -x(n+2) )/12
!   x'''(n) = (-x(n-2) + 2*x(n-1) -          2*x(n+1) +x(n+2) )/2
!   x''''(n)= ( x(n-2) - 4*x(n-1) + 6*x(n) - 4*x(n+1) +x(n+2  )
!
! Two algorithms are implemented for function derivatives and 
! integrals. The first one uses Lagrange interpolation to define 
! the function between mesh points. The second uses cubic splines. 
! Which method is actually used is determined by the interface 
! settings at the beginning of the module. In the Lagrange version, 
! dy/dn is obtained by the same formula as x'(n) above and we then 
! find dy/dx = (dy/dn)/(dx/dn). 
! The integral is approximated by:
!  3-point formula for first and last intervals:
!    integral_from(1)_to(2) y*dn = ( 5*y(1) + 8*y(2) - y(3) )/12
!  4-point formula for 'interior' intervals: 
!    integral_from(n)_to(n+1) f*dn = 
!       ( -y(n-1) + 13*y(n) + 13*y(n+1) - y(n+2) )/24
! In the spline version, d2y/dn2 is first obtained by imposing the
! matching of y, dy/dn, and d2y/dn2 at each mesh point. Derivatives
! and integrals are then obtained staightforwardly from the cubic
! interpolation polynomials at each point and interval:
!   dy(n)/dn = y(n) - y(n-1) + (2*d2y(n)/dn2 + d2y(n-1)/dn2)/6
!            = y(n+1) - y(n) - (d2y(n+1)/dn2 + 2*d2y(n)/dn2)/6
!            = (y(n+1) - y(n-1))/2 - (d2y(n+1)/dn2 - d2y(n-1)/dn2)/12
!   integral_from(x(n))_to(x(n+1)) y(x)*dx = 
!   integral_from(n)_to(n+1) f(n)*dn = 
!      (f(n)+f(n+1))/2 - (d2f(n)/dn2 + d2f(n+1)/dn2)/24
! with f(n)=y(x(n))*dx(n)/dn. In both versions, higher order 
! derivatives of y(x) are obtained by repeated derivation.
!
! The origin of the 1,10,1 coefs of Numerov's formula 
!    (y(x+dx)-2*y(x)+y(x-dx))/dx2 = (f(x+dx)+10*f(x)+f(x-dx))/12
! is as follows:
!   y(x+dx) = y(x) + y'(x)*dx + y''(x)*dx2/2 + 
!             y'''(x)*dx3/6 + y''''(x)*dx4/24 + O(dx5)
! Then: 
!   (y(x+dx)-2*y(x)+y(x-dx))/dx2 = y''(x) + y''''(x)*dx2/12 + O(dx4)
!                                = f(x) + f''(x)*dx2/12 + O(dx4)
! We now approximate
!   f''(x) = (f(x+dx)-2*f(x)+f(x-dx))/dx2 + O(dx2)
! to find the Numerov formula plus O(dx4) corrections.
! For variable meshes x(n), n=1,2,...,N, we define a new function
!   z(n) = y(x(n)) / (x'(n))^(1/2), where x'=dx/dn. Then, after
! some algebra, it can be shown that
!   d2z/dn2 = g(n,z) = g0(n) + g1(n)*z, where
!   g0(n) = f0(x(n))*(x')^(3/2)
!   g1(n) = f1(x(n))*(x')^2 + [3*(x'')^2/x'-2*x''']/(2*x')^2
! If x decreases with n, we can define a new variable -x and 
!   x' --> -x', what amounts to using |x'|^(1/2) and |x'|^(3/2)
! Within the numerov routine, the derivatives yp and ypp are found
! by a special algorithm, using that y''=f and
!   y(x+dx) - y(x-dx) = 2*y'(x) + y'''(x)/3 + O(dx5) = 
!                       2*y'(x) + f'(x)/3 + O(dx5) =
!                       2*y'(x) + (f(x+dx) - f(x-dx))/(6*dx) + O()
! from which y'(x) is solved. The first and last points are 
! special: taking x = xmax-dx and ignoring O(dx5) errors:
!   y'(x+dx) = y'(x) + y''(x)*dx + y'''(x)*dx2/2 + y''''(x)*dx3/6
!            = y'(x) + f(x)*dx + f'(x)*dx2/2 + f''(x)*dx3/6
!            = y'(x) + f(x)*dx + (f(x+dx)-f(x-dx))*dx/4 +
!              (f(x+dx)-2*f(x)+f(x-dx))*dx/6
!            = y'(x) + (5*f(x+dx) + 8*f(x) - f(x-dx)) * dx/12
! similarly, taking x = xmin+dx:
!   y'(x-dx) = y'(x) - (5*f(x-dx) + 8*f(x) - f(x+dx)) * dx/12
! For a variable mesh, these formulas are used to find z'(n). Then
! y'(x) is found from z' and x', using chain's rule as before
!
! Ref: J.M.Soler notes of Jan.12,1998, Nov.28,2001, Dec.4,2001,
!      Jan.16,2002, Jul.29,2002, Jun.21,2007, and Jul.5,2007.
!-----------------------------------------------------------------
! Implementation details:
! - The routines use straightforward formulas and are optimized
!   for simplicity, clarity, and ease of use, NOT for efficiency.
!-----------------------------------------------------------------
! Written by J.M.Soler. Nov.2001, Jul.2002, and Jul.2007
!-----------------------------------------------------------------

use m_recipes

implicit none

PUBLIC :: get_n, set_mesh, set_interpolation, get_mesh, &
          numerov, integral, derivative

PRIVATE   ! Nothing is declared public beyond this point

  integer, parameter :: DP = kind(1.d0)
  character(len=10), save :: interpolation_method='Spline'
  logical,  save :: defined_mesh=.false.
  real(dp), save :: yp1=huge(1._dp), ypn=huge(1._dp)
  real(DP), save, dimension(:), allocatable :: &
    sqrxp, s0, s1, s2, xi, xp1, xp2, xp3, xp4

CONTAINS

function get_n( xmin, xmax, dxmin, dxmax )

  implicit none
  integer              :: get_n
  real(DP), intent(in) :: xmin
  real(DP), intent(in) :: xmax
  real(DP), intent(in) :: dxmin
  real(DP), intent(in) :: dxmax

  real(DP) :: a, b

! Check signs of arguments
  if (dxmax*dxmin<=0._DP .or. (xmax-xmin)*dxmax<=0._DP) &
    stop 'get_n: ERROR: Bad arguments'

! Find required number of points
  if (dxmax==dxmin) then  ! Linear mesh
    get_n = nint( (xmax-xmin) / dxmax )
  else
    b = (xmax-xmin) * dxmin / (dxmax-dxmin)
    a = log( 1 + dxmin/b )
    get_n = 1 + nint( log(1+(xmax-xmin)/b) / a )
  endif

end function get_n

!----------------------------------------------------------------

subroutine set_mesh( n, x, xmin, xmax, a, dxndx1 )

  implicit none
  integer,            intent(in) :: n
  real(DP), optional, intent(in) :: x(n)
  real(DP), optional, intent(in) :: xmin
  real(DP), optional, intent(in) :: xmax
  real(DP), optional, intent(in) :: a
  real(DP), optional, intent(in) :: dxndx1

  integer  :: i
  real(DP) :: aa, b, d, e, x1

! Check the number of points
  if (n<3) then
    print*, 'set_mesh: at least three points required'
    stop
  endif

! Allocate internal tables
  if (defined_mesh) deallocate( sqrxp, s0, s1, s2,     &
                                xi, xp1, xp2, xp3, xp4 )
  allocate( sqrxp(n), s0(n), s1(n), s2(n),        &
            xi(n), xp1(n), xp2(n), xp3(n), xp4(n) )

! Set first point
  if (present(xmin)) then
    x1 = xmin
  else
    x1 = 0
  endif

! Set the mesh points and the derivatives of x(i)
  if (present(x)) then
    xi = x

!   Find derivatives
    xp1 = 0._DP
    xp2 = 0._DP
    xp3 = 0._DP
    xp4 = 0._DP
!   Centered 5-point formulas for all but first/last two points
    do i = 3,n-2
      xp1(i)=( xi(i-2)- 8*xi(i-1)         + 8*xi(i+1)-xi(i+2))/12
      xp2(i)=(-xi(i-2)+16*xi(i-1)-30*xi(i)+16*xi(i+1)-xi(i+2))/12
      xp3(i)=(-xi(i-2)+ 2*xi(i-1)         - 2*xi(i+1)+xi(i+2))/2
      xp4(i)=  xi(i-2)- 4*xi(i-1)+ 6*xi(i)- 4*xi(i+1)+xi(i+2)
    enddo

!   Use first/last 5 points for derivatives of two extreme points
    xp1(1) = xp1(3) - 2*xp2(3) + 4*xp3(3)/2 - 8*xp4(3)/6
    xp1(2) = xp1(3) - 1*xp2(3) + 1*xp3(3)/2 - 1*xp4(3)/6
    xp2(1) = xp2(3) - 2*xp3(3) + 4*xp4(3)/2
    xp2(2) = xp2(3) - 1*xp3(3) + 1*xp4(3)/2
    xp3(1) = xp3(3) - 2*xp4(3)
    xp3(2) = xp3(3) - 1*xp4(3)
    xp4(1) = xp4(3)
    xp4(2) = xp4(3)
    xp1(n)   = xp1(n-2) + 2*xp2(n-2) + 4*xp3(n-2)/2 + 8*xp4(n-2)/6
    xp1(n-1) = xp1(n-2) + 1*xp2(n-2) + 1*xp3(n-2)/2 + 1*xp4(n-2)/6
    xp2(n)   = xp2(n-2) + 2*xp3(n-2) + 4*xp4(n-2)/2
    xp2(n-1) = xp2(n-2) + 1*xp3(n-2) + 1*xp4(n-2)/2
    xp3(n)   = xp3(n-2) + 2*xp4(n-2)
    xp3(n-1) = xp3(n-2) + 1*xp4(n-2)
    xp4(n)   = xp4(n-2)
    xp4(n-1) = xp4(n-2)

  elseif (present(xmax)) then

    if (present(a)) then
      aa = a
    elseif (present(dxndx1)) then
      aa = log(dxndx1) / (n-1)
    else
      aa = 0
    endif
    if (aa > 1.e-6_DP) then
      b = (xmax-x1) / (exp(aa*(n-1)) - 1._DP)
      do i = 1,n
        e = b * exp(aa*(i-1))
        xi(i) = x1 + e - b
        xp1(i) = aa * e
        xp2(i) = aa**2 * e
        xp3(i) = aa**3 * e
        xp4(i) = aa**4 * e
      enddo
    else
      d = (xmax-x1) / (n-1)
      do i = 1,n
        xi(i) = x1 + (i-1) * d
        xp1(i) = d
        xp2(i) = 0
        xp3(i) = 0
        xp4(i) = 0
      enddo
    endif

  else
    print*, 'set_mesh: undefined mesh'
    stop
  endif

! Find auxiliary functions associated to the mesh
  sqrxp = abs(xp1)**0.5_DP
  s0 = abs(xp1)**1.5_DP
  s1 = xp1**2
  s2 = (3._DP*xp2**2 - 2._DP*xp1*xp3) / 4._DP / xp1**2

  defined_mesh = .true.

end subroutine set_mesh

!----------------------------------------------------------------

subroutine set_interpolation( method, ypleft, ypright )

  implicit none
  character(len=*), intent(in):: method
  real(dp),optional,intent(in):: ypleft, ypright

  if (method=='spline' .or. method=='Spline' .or. &
      method=='SPLINE') then
    interpolation_method = 'Spline'
  else if (method=='lagrange' .or. method=='Lagrange' .or. &
           method=='LAGRANGE') then
    interpolation_method = 'Lagrange'
  else
    stop 'set_interpolation: ERROR: unknown method'
  end if

  if (present(ypleft)) then 
    yp1 = ypleft
  else
    yp1 = huge(1._dp)
  end if

  if (present(ypright)) then 
    ypn = ypright
  else
    ypn = huge(1._dp)
  end if

end subroutine set_interpolation

!----------------------------------------------------------------

subroutine get_mesh( n, nx, x, dxdi, d2xdi2, d3xdi3, d4xdi4 )

  implicit none
  integer, intent(in)  :: n   ! size of argument arrays
  integer, intent(out) :: nx  ! number of points
  real(DP), optional, dimension(:), intent(out) :: &
    x, dxdi, d2xdi2, d3xdi3, d4xdi4

  if (.not.defined_mesh) then
    print*,'get_mesh: mesh not defined'
    stop
  endif
  nx = size(xi)
  if (present(x) .or. present(dxdi) .or. present(d2xdi2) .or. &
      present(d3xdi3) .or. present(d4xdi4)) then
    nx = max( nx, n )
  end if

  if (present(x))      x(1:nx)      = xi(1:nx)
  if (present(dxdi))   dxdi(1:nx)   = xp1(1:nx)
  if (present(d2xdi2)) d2xdi2(1:nx) = xp2(1:nx)
  if (present(d3xdi3)) d3xdi3(1:nx) = xp3(1:nx)
  if (present(d4xdi4)) d4xdi4(1:nx) = xp4(1:nx)

end subroutine get_mesh

!----------------------------------------------------------------

recursive function derivative( n, y, x, dx, order ) result(dydx)

  implicit none
  integer,            intent(in) :: n
  real(DP),           intent(in) :: y(n)
  real(DP), optional, intent(in) :: x(n)
  real(DP), optional, intent(in) :: dx
  integer,  optional, intent(in) :: order
  real(DP)                       :: dydx(n)

  integer  :: i, nx, ord
  real(DP) :: dxdi(n), dydi(n), d2ydi2(n)

! Mesh-initialization
  if (present(x)) then
    call set_mesh( n, x=x )
  elseif (present(dx)) then
    call set_mesh( n, xmax=(n-1)*dx )
  endif

! Fix order of the derivative
  if (present(order)) then
    ord = order
  else
    ord = 1
  endif

! Iterate over order for order>1
  if (ord > 1) then
    dydx = derivative( n, derivative(n,y), order=ord-1 )
  elseif (ord == 1) then

    if (interpolation_method=='Spline') then

      call spline( 1._dp, y, n, yp1, ypn, d2ydi2 )

      dydi(1) = y(2) - y(1) - (2*d2ydi2(1) + d2ydi2(2)) / 6
      dydi(n) = y(n) - y(n-1) + (d2ydi2(n-1) + 2*d2ydi2(n)) / 6
      dydi(2:n-1) = (y(3:n) - y(1:n-2)) / 2 &
                  - (d2ydi2(3:n) - d2ydi2(1:n-2)) / 12

    else if (interpolation_method=='Lagrange') then

!     Find derivative of y with respect to i
!     Centered 5-point formulas for all but first/last two points
      do i = 3,n-2
        dydi(i)=(y(i-2)-8*y(i-1)+8*y(i+1)-y(i+2))/12
      enddo

!     Find derivative at first/last two points
      if (n<=1) then
        dydi = 0
      else if (n==2) then
        dydi = -y(1) + y(2)
      else if (n==3) then
        dydi(1) = (-3*y(1) + 4*y(2) -   y(3)) / 2
        dydi(3) = (   y(1) - 4*y(2) + 3*y(3)) / 2
        dydi(2) = (-  y(1)          +   y(3)) / 2
      else if (n==4) then
        ! Use up to 2 neighbor points on each side
!        dydi(1) = ( -3*y(1) + 4*y(2) -   y(3)) / 2
!        dydi(4) = (    y(2) - 4*y(3) + 3*y(4)) / 2
!        dydi(2) = (- 2*y(1) - 3*y(2) + 6*y(3) -   y(4)) / 6
!        dydi(3) = (    y(1) - 6*y(2) + 3*y(3) + 2*y(4)) / 6
        ! Alternatively: use all available points
        dydi(1) = (-11*y(1) +18*y(2) - 9*y(3) + 2*y(4)) / 6
        dydi(4) = (- 2*y(1) + 9*y(2) -18*y(3) +11*y(4)) / 6
      else
        ! Use up to 2 neighbor points on each side
!        dydi(1) = ( -3*y(1)   + 4*y(2)   -   y(3)) / 2
!        dydi(n) = (    y(n-2) - 4*y(n-1) + 3*y(n)) / 2
!        dydi(2)   = (- 2*y(1)   - 3*y(2)   + 6*y(3)   -   y(4)) / 6
!        dydi(n-1) = (    y(n-3) - 6*y(n-2) + 3*y(n-1) + 2*y(n)) / 6
        ! Alternatively: use always 5 points
        dydi(1)=(-25*y(1)+48*y(2)  -36*y(3)  +16*y(4)  -3*y(5)  )/12
        dydi(n)=(+25*y(n)-48*y(n-1)+36*y(n-2)-16*y(n-3)+3*y(n-4))/12
        dydi(2)  =(-3*y(1)-10*y(2)  +18*y(3)  -6*y(4)  +y(5)    )/12
        dydi(n-1)=(+3*y(n)+10*y(n-1)-18*y(n-2)+6*y(n-3)-y(n-4)  )/12
      end if ! (n)

    else
      stop 'mesh1D: ERROR: bad interpolation_method parameter'
    end if ! (interpolation_method)

!   Find derivative of x with respect to i
    call get_mesh( n, nx, dxdi=dxdi )

!   Find derivative of y with respect to x
    dydx = dydi / dxdi

  endif

end function derivative

!----------------------------------------------------------------

function integral( n, y, x, dx )

  implicit none
  real(DP)                       :: integral
  integer,            intent(in) :: n
  real(DP),           intent(in) :: y(n)
  real(DP), optional, intent(in) :: x(n)
  real(DP), optional, intent(in) :: dx

  real(DP):: d2fdi2(n), f(n)

! Mesh-initialization
  if (present(x)) then
    call set_mesh( n, x=x )
  elseif (present(dx)) then
    call set_mesh( n, xmax=(n-1)*dx )
  endif

! Find integral
  if (interpolation_method=='Spline') then

    f = y * xp1
    call spline( 1._dp, f, n, yp1, ypn, d2fdi2 )

    integral = (f(1) + f(n))/2 - (d2fdi2(1) + d2fdi2(n)) / 24 &
             + sum(f(2:n-1)) - sum(d2fdi2(2:n-1)) / 12
    
  else if (interpolation_method=='Lagrange') then

    if (n<=1) then
      integral = 0
    else if (n==2) then
      integral = ( y(1)*xp1(1) + y(2)*xp1(2) ) / 2
    else if (n==3) then
      integral = ( y(1)*xp1(1) + 4*y(2)*xp1(2) + y(3)*xp1(3) ) / 3
    else if (n==4) then
      integral = (   3 * (y(1)*xp1(1) + y(n)  *xp1(n)  )      &
                  +  9 * (y(2)*xp1(2) + y(n-1)*xp1(n-1)) )/8
    else if (n==5) then
      integral = (   9 * (y(1)*xp1(1) + y(n)  *xp1(n)  )      &
                  + 28 * (y(2)*xp1(2) + y(n-1)*xp1(n-1))      &
                  + 22 *  y(3)*xp1(3)                    )/24
    else
      integral = (   9 * (y(1)*xp1(1) + y(n)  *xp1(n)  )      &
                  + 28 * (y(2)*xp1(2) + y(n-1)*xp1(n-1))      &
                  + 23 * (y(3)*xp1(3) + y(n-2)*xp1(n-2)) )/24 &
               + sum( y(4:n-3)*xp1(4:n-3) )
    end if

  else
    stop 'mesh1D: ERROR: bad interpolation_method parameter'
  end if ! (interpolation_method)

end function integral

!----------------------------------------------------------------

subroutine numerov( n, y, yp, ypp, f0, f1, x, dx, norm )

  implicit none
  integer,            intent(in)    :: n
  real(DP), optional, intent(inout) :: y(n)
  real(DP), optional, intent(out)   :: yp(n)
  real(DP), optional, intent(out)   :: ypp(n)
  real(DP), optional, intent(in)    :: f0(n)
  real(DP), optional, intent(in)    :: f1(n)
  real(DP), optional, intent(in)    :: x(n)
  real(DP), optional, intent(in)    :: dx
  real(DP), optional, intent(in)    :: norm

  integer                :: i
  real(DP)               :: ynorm
  real(DP), dimension(n) :: g, g0, g1, z, zp

! Mesh-initialization
  if (present(x)) then
    call set_mesh( n, x=x )
  elseif (present(dx)) then
    call set_mesh( n, xmax=(n-1)*dx )
  endif
  if (.not.present(y)) return

! Check that y has been initialized in left boundary
  if (.not.present(f0)      & ! Homogeneous equation
      .and.abs(y(1))==0._DP &
      .and.abs(y(2))==0._DP) then
    print*, 'numerov: first two values of y needed'
    stop
  endif

! Find g0 and g1
  g0 = 0._DP
  g1 = s2
  if (present(f0)) g0 = s0 * f0
  if (present(f1)) g1 = s1 * f1 + g1

! Integrate differential equation
  z(1) = y(1) / sqrxp(1)
  z(2) = y(2) / sqrxp(2)

  do i = 2,n-1
    z(i+1) = ( (g0(i-1)+10._DP*g0(i)+g0(i+1))/12._DP +   &
               (-1._DP+g1(i-1)/12._DP)*z(i-1)        +   &
               (2._DP+(10._DP/12._DP)*g1(i))*z(i)    ) / &
                  (1._DP-g1(i+1)/12._DP)
  enddo
  y = z * sqrxp

! Find first derivative y'=dy/dx
  if (present(yp)) then
    g = g0 + g1*z
    zp(2:n-1) = (z(3:n) - z(1:n-2))/2 + (g(3:n) - g(1:n-2))/6
    zp(n) = zp(n-1) + (5*g(n) + 8*g(n-1) - g(n-2))/12
    zp(1) = zp(2)   - (5*g(1) + 8*g(2)   - g(3)  )/12
    yp = ( zp + z*xp2/(2*xp1) )/sqrxp
  endif

! Find second derivative y''=d2y/dx2
  if (present(ypp)) then
    ypp = 0
    if (present(f0)) ypp = f0
    if (present(f1)) ypp = ypp + f1*y
  endif

! Normalize solution
  if (present(norm)) then
    if (present(f0)) then
      print*, 'numerov: cannot normalize inhomogeneous solutions'
      stop
    elseif (norm <= 0._DP) then
      print*, 'numerov: nonpositive norm =', norm
      stop
    endif
    ynorm = integral( n, y=y*y )
    y = y * sqrt(norm/ynorm)
    if (present(yp))  yp  = yp  * sqrt(norm/ynorm)
    if (present(ypp)) ypp = ypp * sqrt(norm/ynorm)
  endif

end subroutine numerov

!----------------------------------------------------------------

end module mesh1D
