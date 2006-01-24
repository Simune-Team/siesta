      module files
!
!     Contains the short system label, used to generate file names
!     slabel is currently set in reinit.
!
      integer, parameter, public                  :: label_length = 60
      character(len=label_length), save, public   :: slabel

      private

      end module files
