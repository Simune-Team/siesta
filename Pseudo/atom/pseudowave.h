c
c     pswfnrd/pswfnru are the pseudowave functions
c
      double precision pswfnrd(lmax,nrmax), pswfnru(lmax,nrmax)
c
c     valence pswfs, only "major" set (down), with possibly
c     several shells with the same l
c
      integer nmax_shells
      parameter (nmax_shells = 10)

      integer nshells_stored
      double precision pswf(nmax_shells,nrmax)
      integer          n_pswf(nmax_shells)
      integer          l_pswf(nmax_shells)
c
      common /pseudowave/ pswfnrd, pswfnru
      common /pswf_data/ pswf
      common /pswf_qnumbers/ n_pswf, l_pswf, nshells_stored
      save /pseudowave/
      save /pswf_data/
      save /pswf_qnumbers/
c------
