module marsaglia
    !This module contains several very efficient and widely used random number 
    !generators by G. Marsaglia (Florida State University). 
    !All of the subroutines were taken from public sources and modified 
    !by Sergiy Bubin (April-June 2009)
    implicit none
    private
    public :: rkiss, ikiss, kiss_setseed, kiss_getseed, ziggurat_set, rndnor, rndexp, rndunitsph
    
    !Seeds for kiss generator
    integer(4),save :: ks_x=123456789
    integer(4),save :: ks_y=362436069
    integer(4),save :: ks_z=521288629
    integer(4),save :: ks_w=916191069
    
    !Parameters, tables, and seeds for random normals and random exponentials
    real(8),parameter :: m1=2147483648.0D0   
    real(8),parameter :: m2=2147483648.0D0
    real(8)           :: dn=3.442619855899D0
    real(8)           :: tn=3.442619855899D0
    real(8),parameter :: vn=0.00991256303526217D0
    real(8)           :: de=7.697117470131487D0
    real(8)           :: te=7.697117470131487D0
    real(8),parameter :: ve=0.003949659822581572D0
    integer(4)        :: ip
    integer(4)        :: kn(0:127)
    integer(4)        :: ke(0:255)
    integer(4)        :: hz
    real(8)           :: wn(0:127)
    real(8)           :: fn(0:127)
    real(8)           :: we(0:255)
    real(8)           :: fe(0:255)
    logical           :: init_done=.false.
    
    contains
    
    
    !The  KISS (Keep It Simple Stupid) random number generator of G. Marsaglia. 
    !It combines:
    !(1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
    !(2) A 3-shift shift-register generator, period 2^32-1,
    !(3) Two 16-bit multiply-with-carry generators, period 597273182964842497>2^59
    !Overall period>2^123;  Default seeds ks_x,ks_y,ks_z,ks_w. 
    !The real number generator is called rkiss(), the integer generator is called 
    !ikiss(). Both ikiss() and rkiss() share the same set of seeds.
    !Set your own seeds with 
    !call kiss_setseed(ix,iy,iz,iw).
    !Get current seeds with
    !call kiss_getseed(ix,iy,iz,iw)
    
    
    function ikiss()
    !Function ikiss generates a random integer. 
    integer(4)  ikiss
    ks_x = 69069 * ks_x + 1327217885
    ks_y = m (m (m (ks_y, 13), - 17), 5)
    ks_z = 18000 * iand (ks_z, 65535) + ishft (ks_z, - 16)
    ks_w = 30903 * iand (ks_w, 65535) + ishft (ks_w, - 16)
    ikiss = ks_x + ks_y + ishft (ks_z, 16) + ks_w
    contains
      function m(k, n)
      integer m, k, n
      m = ieor (k, ishft (k, n) )
      end function m
    end function ikiss
    
    
    function rkiss()
    !Function rkiss generates a random number between 0 and 1 
    !(numbers exactly equal to 0.0D0 and 1.0D0 do not happen)
    real(8)  rkiss 
    integer(4)  ikiss
    ks_x = 69069 * ks_x + 1327217885
    ks_y = m (m (m (ks_y, 13), - 17), 5)
    ks_z = 18000 * iand (ks_z, 65535) + ishft (ks_z, - 16)
    ks_w = 30903 * iand (ks_w, 65535) + ishft (ks_w, - 16)
    ikiss = ks_x + ks_y + ishft (ks_z, 16) + ks_w
    !rkiss = ikiss * 1/(2^32-1) + 2^31/(2^32-1)
    rkiss = ikiss*2.32830643D-10 + 0.500000000116415D0
    contains
      function m(k, n)
      integer m, k, n
      m = ieor (k, ishft (k, n) )
      end function m
    end function rkiss
    
    
    subroutine kiss_setseed(ix, iy, iz, iw)
    !Sets seeds for kiss generator
    integer(4) :: ix, iy, iz, iw
    ks_x = ix
    ks_y = iy
    ks_z = iz
    ks_w = iw
    end subroutine kiss_setseed
    
    
    subroutine kiss_getseed(ix, iy, iz, iw)
    !Gets current seeds of kiss generator 
    integer(4) :: ix, iy, iz, iw
    ix = ks_x
    iy = ks_y
    iz = ks_z
    iw = ks_w
    end subroutine kiss_getseed
    
    
    !Marsaglia & Tsang generator for random normals & random exponentials,
    !which uses the so called ziggurat algorithm.
    !Translated from C by Alan Miller. Modified by Sergiy Bubin.
    !Reference:
    !Marsaglia, G. & Tsang, W.W. (2000) `The ziggurat method for generating
    !random variables', J. Statist. Software, v5(8).
    
    
    subroutine ziggurat_set()
    !This subroutine initializes tables for the generators
    !of random normals and random exponentials. It is called once and 
    !automatically inside either rndnor or rndexp, when they are called 
    !for the first time.
    integer(4)             :: i
    real(8)                :: q
    
    !Tables for rndnor
    q = vn*exp(0.5D0*dn*dn)
    kn(0) = (dn/q)*m1
    kn(1) = 0
    wn(0) = q/m1
    wn(127) = dn/m1
    fn(0) = 1.0D0
    fn(127) = exp( -0.5D0*dn*dn )
    do  i = 126, 1, -1
      dn = sqrt( -2.0D0 * log( vn/dn + exp(-0.5D0*dn*dn ) ) )
      kn(i+1) = (dn/tn)*m1
      tn = dn
      fn(i) = exp(-0.5D0*dn*dn)
      wn(i) = dn/m1
    end do
    
    !Tables for rndexp
    q = ve*exp( de )
    ke(0) = (de/q)*m2
    ke(1) = 0
    we(0) = q/m2
    we(255) = de/m2
    fe(0) = 1.0D0
    fe(255) = exp( -de )
    do  i = 254, 1, -1
      de = -log( ve/de + exp( -de ) )
      ke(i+1) = m2 * (de/te)
      te = de
      fe(i) = exp( -de )
      we(i) = de/m2
    end do
    init_done = .true.
    return
    
    end subroutine ziggurat_set
    
    
    function rndnor() result(fn_val)
    !Function  rndnor generates random normals
    !Note that rndnor calls ikiss and rkiss. Hence, a call of rndnor changes
    !current seeds for those generators.
    real(8)             ::  fn_val
    real(8), parameter  ::  r = 3.442620D0
    real(8)             ::  x, y
    
    if ( .not. init_done ) call ziggurat_set()
    hz = ikiss()
    ip = iand( hz, 127 )
    if ( abs( hz ) < kn(ip) ) then
      fn_val = hz * wn(ip)
    else
      do
        if ( ip == 0 ) then
           do
             x = -0.2904764D0* log(rkiss())
             y = -log(rkiss())
             if( y+y >= x*x ) exit
           end do
           fn_val = r+x
           if ( hz <= 0 ) fn_val = -fn_val
           return
        end if
        x = hz * wn(ip)
        if ( fn(ip) + rkiss()*(fn(ip-1)-fn(ip)) < exp(-0.5D0*x*x) ) then
          fn_val = x
          return
        end if
        hz = ikiss()
        ip = iand( hz, 127 )
        if ( abs( hz ) < kn(ip) ) then
          fn_val = hz * wn(ip)
          return
        end if
      end do
    end if
    return
    
    end function rndnor
    
    
    function rndexp() result(fn_val)
    !Function  rndnor generates random exponentials
    !Note that rndexp calls ikiss and rkiss. Hence, a call of rndexp changes
    !current seeds for those generators.
    real(8)    ::  fn_val
    real(8)    ::  x
    integer(4) ::  jz
    
    if( .not. init_done ) call ziggurat_set()
    jz = ikiss()
    ip = iand( jz, 255 )
    if( abs( jz ) < ke(ip) ) then
      fn_val = abs(jz) * we(ip)
      return
    end if
    do
      if( ip == 0 ) then
        fn_val = 7.69711 - log( rkiss() )
        return
      end if
      x = abs( jz ) * we(ip)
      if( fe(ip) + rkiss()*(fe(ip-1) - fe(ip)) < exp( -x ) ) then
        fn_val = x
        return
      end if
      jz = ikiss()
      ip = iand( jz, 255 )
      if( abs( jz ) < ke(ip) ) then
        fn_val = abs( jz ) * we(ip)
        return
      end if
    end do
    return
    
    end function rndexp
    
    
    !Marsaglia's method of picking random points on a surface of a unit sphere
    !Reference:
    !Marsaglia, G. "Choosing a Point from the Surface of a Sphere." 
    !Ann. Math. Stat. 43, 645-646, 1972. 
    
    subroutine rndunitsph(x,y,z)
    !Subroutine rndunitsph returns in x,y,z a unit vector of random
    !direction. Note that rndunitsph calls rkiss. Hence, a call of rndunitsph 
    !changes current seeds for that generator.
    real(8)  x,y,z
    real(8)  x1,x2,t,r
    
    x1=2*(rkiss()-0.5D0)
    x2=2*(rkiss()-0.5D0)
    t=x1**2+x2**2
    do while (t>1.0D0)
      x1=2*(rkiss()-0.5D0)
      x2=2*(rkiss()-0.5D0)
      t=x1**2+x2**2
    enddo
    r=2.0D0*sqrt(1.0D0-t)
    x=x1*r
    y=x2*r
    z=1.0D0-2.0D0*t
    end subroutine rndunitsph
    
    
    end module marsaglia
    
    