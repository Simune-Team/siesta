c Common arrays for plrho routines
      integer maxx, maxy
      parameter ( maxx = 1240 )
      parameter ( maxy = 1024 )
      integer ixmin, ixmax, iymin, iymax,
     .        pixmap(0:maxx-1,0:maxy-1)
      common /compix/ ixmin, ixmax, iymin, iymax, pixmap

