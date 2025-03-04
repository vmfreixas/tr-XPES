    program test_read_movecs

    implicit none

    character*30 filename
    integer nbf
    integer iset ! 1 for closed shell
                 ! 2 for open shell
    double precision, dimension(:), allocatable :: occ
    double precision, dimension(:), allocatable :: evals
    double precision, dimension(:), allocatable :: g_vecs

    filename = 'HBQ_Neutral_janpa.movecs'
    iset = 1
    nbf = 519
    allocate(occ(nbf), evals(nbf))

    print*, filename

    call movecs_read(filename, iset) 

    print*, 'done'

    contains

!--------------------------------------------------------------------!
! movecs_read subrouitne adapted from NWChem ("src/ddscf/vectors.F") !
!--------------------------------------------------------------------!

      subroutine movecs_read(filename, iset)
      implicit none

      character*(*) filename
      integer iset              ! No. (1,2) of set of vectors to read
      double precision, allocatable, dimension(:) :: occ  ! Must be at least nbf long (not nmo)
      double precision, allocatable, dimension(:) :: evals ! Must be at least nbf long (not nmo)
      double precision, allocatable, dimension(:,:) :: g_vecs

      integer nsets             ! No. of sets of vectors
      integer nbf               ! No. of functions in basis
      integer nmo(2)            ! No. of vectors in each set
      integer ok, jset, i, j
      double precision, allocatable, dimension(:) :: l_vecs, k_vecs
      integer unitno
      parameter (unitno = 67)

      ok = 0
      open(unitno, status='old', form='unformatted', convert='little_endian', file=filename, err=1000)
!
!     Skip over uninteresting bits of the header
!
         read(unitno, err=1001, end=1001) ! convergence info
         read(unitno, err=1001, end=1001) ! scftype
         read(unitno, err=1001, end=1001) ! lentit
         read(unitno, err=1001, end=1001) ! title
         read(unitno, err=1001, end=1001) ! lenbas
         read(unitno, err=1001, end=1001) ! basis_name
         read(unitno, err=1001, end=1001) nsets
         read(unitno, err=1001, end=1001) nbf
         read(unitno, err=1001, end=1001) (nmo(i),i=1,nsets)

         print *, 'nsets: ', nsets
         print *, 'nbf: ', nbf
         print *, 'nmo: ', nmo

         allocate(l_vecs(nbf))
         allocate(k_vecs(nbf))
         allocate(occ(nbf))
         allocate(evals(nbf))
         allocate(g_vecs(nbf, nmo(1)))
!
!     Skip over unwanted sets
!
         do jset = 1, iset-1
            read(unitno, err=1001, end=1001)
            read(unitno, err=1001, end=1001)
            do i = 1, nmo(jset)
               read(unitno, err=1001, end=1001)
            enddo
         enddo
         read(unitno, err=1001, end=1001) (occ(j), j=1, nbf)
         read(unitno, err=1001, end=1001) (evals(j), j=1, nbf)
         do i = 1, nmo(iset)
            read(unitno, err=1001, end=1001) (k_vecs(j), j=1, nbf)
!            call sread(unitno, k_vecs, nbf)
            g_vecs(:,i) = k_vecs
         enddo
 9    close(unitno, err=1002)
      do i = 1, nmo(iset)
         print *, g_vecs(:,i)
      enddo
    goto 10

 1000 print *, ' movecs_read: failed to open ', trim(adjustl(filename))
    goto 10

 1001 print *, ' movecs_read: failing reading from ', trim(adjustl(filename))
      close(unitno, err=1002)
    goto 10
 1002 print *, ' movecs_read: failed to close', trim(adjustl(filename))

 10 continue

      end subroutine movecs_read

    end program test_read_movecs
