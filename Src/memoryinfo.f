      module memoryinfo

C
C WordsInteger = number of words of memory for integer variables
C WordsSP      = number of words of memory for single precision variables
C WordsDP      = number of words of memory for double precision variables
C WordsLogical = number of words of memory for logical variables
C WordsSC      = number of words of memory for single complex variables
C WordsDC      = number of words of memory for double complex variables
C PeakMemory   = maximum amount of dynamic memory used
C PeakRoutine  = routine name where memory peak was reached
C ByteSize     = array of values specifying the size of a word for
C                each data type
C
      implicit none

      integer, save ::
     .  WordsInteger, WordsSP, WordsDP, WordsLogical, WordsSC, WordsDC,
     .  PeakMemory, CurrentMemory, ByteSize(6)

      character(len=30), save ::
     .  PeakRoutine

      data
     .  WordsInteger / 0 /,
     .  WordsSP / 0 /,
     .  WordsDP / 0 /,
     .  WordsLogical / 0 /,
     .  WordsSC / 0 /,
     .  WordsDC / 0 /,
     .  PeakMemory / 0 /,
     .  CurrentMemory / 0 /,
     .  ByteSize /4,4,8,4,8,16/

      end module memoryinfo

      module diagmemory
C
C  Stores the factor used to scale the default memory in rdiag/cdiag
C  By increasing this value it is possible to avoid failure to
C  converge eigenvalues.
C
C  real*8  MemoryFactor      : factor by which memory is scaled
C  logical TryMemoryIncrease : if .true. then diagonalisation will
C                            : run iteratively until enough memory
C                            : is obtained for convergence
C
      implicit none

      real*8, save ::
     .  MemoryFactor

      logical, save ::
     .  TryMemoryIncrease

      end module diagmemory
