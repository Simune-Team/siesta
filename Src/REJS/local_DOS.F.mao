***************
*** 18,29 ****
        use atomlist,       only: datm, no_s, iaorb        
        use fdf,            only : fdf_block, fdf_convfac  
        use sys,            only: die                      
-       use m_mpi_utils,    only : broadcast               
        use files,          only : slabel     
        use Kpoint_grid
        use parallel,       only: IOnode                   
        use files,          only : label_length            
-       use parse
        use m_ntm
        use m_forces,       only: fa
        use m_eo
--- 18,27 ----
        use atomlist,       only: datm, no_s, iaorb        
        use fdf,            only : fdf_block, fdf_convfac  
        use sys,            only: die                      
        use files,          only : slabel     
        use Kpoint_grid
        use parallel,       only: IOnode                   
        use files,          only : label_length            
        use m_ntm
        use m_forces,       only: fa
        use m_eo
