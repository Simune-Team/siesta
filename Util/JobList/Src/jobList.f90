module jobList

  implicit none

PUBLIC:: &
  runFile  ! reads and runs a job list file

PRIVATE    ! nothing is declared public beyond this point

  integer,parameter:: dp = kind(1.d0)
  integer,parameter:: maxLines = 1000
  integer,parameter:: maxWords = 100
  integer,parameter:: maxJobs = 1000
  integer,parameter:: unitIn = 42       ! I/O unit for .fdf output files
  integer,parameter:: unitOut = 43      ! I/O unit for .fdf output files

  ! Derived type to hold calculation results
  type resultsType
    private
    real(dp)        :: cell(3,3)   =0._dp
    real(dp)        :: stress(3,3) =0._dp
    real(dp)        :: volume      =0._dp
    real(dp)        :: energy      =0._dp
    real(dp)        :: pressure    =0._dp
    real(dp)        :: virial      =0._dp 
    real(dp),pointer:: coords(:,:) =>null()
    real(dp),pointer:: force(:,:)  =>null()
    integer, pointer:: za(:)       =>null()   ! atomic numbers
  end type resultsType

  integer,           save:: myUnit, totJobs, totLists, nResults
  character(len=120),save:: queue = 'siesta < $1.fdf > $1.out'
  character(len=120),save:: request(maxWords) = ' '

CONTAINS

!------------------------------------------------------------------------------

subroutine runFile( unit, nJobs, nLists )

  implicit none
  integer,         intent(in) :: unit   ! IO unit of datafile
  integer,optional,intent(out):: nJobs  ! total number of jobs
  integer,optional,intent(out):: nLists ! total number of job lists

  integer :: iLine, iWord, nc, nWords
  character(len=120):: line, listName, myDir, words(maxWords)

  ! Copy unit to module variable
  myUnit = unit

  ! Get current directory and add trailing '/' if necessary
  call getcwd(myDir)
  nc = len(trim(myDir))
  if (myDir(nc:nc) /= '/') myDir = myDir(1:nc) // '/'

  ! Initialize total number of jobs, lists, and requested results
  request(1) = 'energy'
  nResults = 1
  totJobs = 0
  totLists = 0

  ! Loop on lines of datafile
  do iLine = 1,maxLines

    ! Read and parse one line
    read(myUnit,'(a)',end=999) line
    call parser(line,words,nWords)

    if (nWords==0) then
      cycle
    elseif (trim(words(1))=='%queue') then
      line = adjustl(line)
      queue = adjustl(line(7:))              ! remove '%queue' from line
    elseif (trim(words(1))=='%result') then
      request = ' '
      request(1:nWords-1) = words(2:nWords)
      nResults = nWords-1
    elseif (trim(words(1))=='%list') then
      listName = trim(words(2))
      call runList(listName,myDir)
    endif

  end do ! iLine

999 continue  ! come here when end of file

  if (present(nJobs))  nJobs = totJobs
  if (present(nLists)) nLists = totLists

end subroutine runFile

!------------------------------------------------------------------------------

recursive subroutine runList( listName, dir )

  implicit none
  character(len=*),intent(in):: listName  ! name of job list
  character(len=*),intent(in):: dir       ! parent directory

  integer :: iJob, iLine, nWords, nJobs
  character(len=120):: jobDir(maxJobs), jobName(maxJobs), &
                       line, myDir, newList, words(maxWords)
  real(dp):: results(maxWords)
  logical :: finished

  ! Create a new directory for this list
  myDir = trim(dir) // trim(listName) // '/'
  call system('mkdir ' // trim(myDir))

  ! Copy input files to job directory
  call system('cp ' // trim(dir) //'siesta '// &
                       trim(dir) //'*.fdf ' // &
                       trim(dir) //'*.psf ' // trim(myDir) )

  ! Loop on lines of datafile
  finished = .false.
  nJobs = 0
  do iLine = 1,maxLines

    ! Read and parse one line, ignoring comments
    read(myUnit,'(a)',end=999) line
    call parser(line,words,nWords)

    if (nWords==0) then                  ! blank or comment line
      cycle
    elseif (trim(words(1))=='%queue') then
      line = adjustl(line)
      queue = adjustl(line(7:))          ! remove '%queue' from line
    elseif (trim(words(1))=='%result') then
      request = ' '
      request(1:nWords-1) = words(2:nWords)
      nResults = nWords-1
    elseif (words(1)=='%list') then      ! new sub-list
      newList = words(2)
      call runList(newList,myDir)
    elseif (words(1)=='%endlist') then   ! end of my list
      if (words(2)==listName) then
        finished = .true.
        exit ! iLine loop
      else
        print*,'runList ERROR: inconsistent %list = '//trim(listName)// &
               ', %endList = '//trim(words(2))
        stop
      endif
    else                                 ! new job
      nJobs = nJobs+1
      call runJob( words, nWords, myDir, jobDir(nJobs), jobName(nJobs) )
    endif

  end do ! iLine

  ! Write results summary file and return
  if (finished) then
    if (.false.) then  ! debug switch
    open(unitOut,file=trim(myDir)//trim(listName)//'.result')
    do iJob = 1,nJobs
      call readResult( jobDir(iJob), jobName(iJob), results=results )
      write(unitOut,'(10e15.6)') results(1:nResults) 
    end do
    close(unitOut)
    end if ! debug switch
    totJobs = totJobs+nJobs
    totLists = totLists+1
    return                               ! this is the normal return point
  else
    stop 'runList ERROR: parameter maxLines too small'
  end if ! finished

  ! Come here only in case of anomalous end of file condition
999 continue
  print*,'runList ERROR: unterminated %list = '//trim(listName)
  stop

end subroutine runList

!------------------------------------------------------------------------------

subroutine runJob( words, nWords, dir, jobDir, jobName )

  implicit none
  integer,         intent(in) :: nWords        ! number of words
  character(len=*),intent(in) :: words(nWords) ! words of the job line
  character(len=*),intent(in) :: dir           ! work directory
  character(len=*),intent(out):: jobDir        ! job directory
  character(len=*),intent(out):: jobName       ! job name

  integer:: ic, iLine, iWord, jWord, nc, nLines
  character(len=120):: line(nWords), queueJob, word

  jobName = ' '
  nLines = 0
  do iWord = 1,nWords
    nLines = nLines+1
    word = words(iWord)
    nc = len(trim(word))               ! number of characters in word
    if (word(nc-3:nc)=='.fdf') then    ! include new .fdf file
      line(nLines) = '%include '//trim(word)
      jobName = trim(jobName)//'_'//word(1:nc-4)
    else                               ! make a line with rest of words
      line(nLines) = word
      jobName = trim(jobName)//'_'//word
      do jWord = iWord+1,nWords
        line(nLines) = trim(line(nLines))//' '//words(jWord)
        jobName = trim(jobName)//'_'//words(jWord)
      end do
      exit ! do iWord loop
    end if
  end do ! iWord
  if (jobName(1:1)=='_') jobName=jobName(2:)

  ! Create a new directory for this job
  jobDir = trim(dir) // trim(jobName) // '/'
  call system('mkdir ' // trim(jobDir))

  ! Copy input files to job directory
  call system('cp siesta *.fdf *.psf '//trim(jobDir))

  ! Write .fdf file for job
  open(unitOut,file=trim(jobDir)//trim(jobName)//'.fdf')
  do iLine = nLines,1,-1               ! last lines (words) have priority
    write(unitOut,'(a)') trim(line(iLine))
  end do
  close(unitOut)

  ! Set job
  queueJob = queue
  nc = len(trim(queueJob))
  do ic = nc-1,1,-1
    if (queueJob(ic:ic+1)=='$1') queueJob(ic:) = trim(jobName)//queueJob(ic+2:)
  end do

  ! Submit job
  call chdir(trim(jobDir))
  call system(trim(queueJob))
  call chdir(trim(dir))

end subroutine runJob

!------------------------------------------------------------------------------

subroutine parser( line, words, nWords )

  implicit none
  character(len=*),intent(inout):: line      ! returns without comments
  character(len=*),intent(out)  :: words(:)  ! words in line, without comments
  integer,         intent(out)  :: nWords    ! number of words

  character(len=120):: myLine
  integer:: iWord, nc

  ! Remove comments
  nc = scan(line,'#')-1          ! last character before comment
  if (nc<0) nc=len(trim(line))   ! if no comment, last nonblank character
  line = line(1:nc)
  myLine = line

  ! Loop on words
  do iWord = 1,size(words)
    myLine = adjustl(myLine)     ! remove leading blanks
    nc = scan(myLine,' ')-1      ! last character before next blank
    if (nc<=0) then              ! no more words
      nWords = iWord-1
      return                     ! normal return point
    else
      words(iWord) = myLine(1:nc)
      myLine = myLine(nc+1:)
    endif
  end do ! iWord
  stop 'parser ERROR: size(words) too small'

end subroutine parser

!------------------------------------------------------------------------------

subroutine readResult( dir, name, result, results )

  implicit none
  character(len=*),          intent(in) :: dir         ! job directory
  character(len=*),          intent(in) :: name        ! job name
  type(resultsType),optional,intent(out):: result      ! job-results structure
  real(dp),         optional,intent(out):: results(:)  ! requested job results

  integer :: ia, ic, is, iostat, iResult, nAtoms, za
  real(dp):: c(3,3)
  character(len=120):: fileName
  type(resultsType) :: res

  ! Read final structure
  fileName = trim(dir)//trim(name)//'.XV'
  open(unitIn,file=trim(fileName),status='old',iostat=iostat)
  if (iostat==0) then
    do ic = 1,3
      read(unitIn,*) res%cell(:,ic)
    end do
    read(unitIn,*) nAtoms
    if (associated(res%za)) deallocate(res%za)
    if (associated(res%coords)) deallocate(res%coords)
    allocate( res%za(nAtoms), res%coords(3,nAtoms) )
    do ia = 1,nAtoms
      read(unitIn,*) is, res%za(ia), res%coords(:,ia)
    end do
  else
    res%cell = 0
    res%coords = 0
    res%za = 0
  end if
  close(unitIn)
  
  ! Read energy, force, and stress
  fileName = trim(dir)//'FORCE_STRESS'
  open(unitIn,file=trim(fileName),status='old',iostat=iostat)
  if (iostat==0) then
    read(unitIn,*) res%energy
    read(unitIn,*) res%stress
    read(unitIn,*) nAtoms
    if (associated(res%force)) deallocate(res%force)
    allocate(res%force(3,nAtoms))
    do ia = 1,nAtoms
      read(unitIn,*) is, res%za(ia), res%force(:,ia)
    end do
  else
    res%energy = 0
    res%stress = 0
    res%force = 0
  end if
  close(unitIn)

  ! Find cell volume, pressure, and virial
  c = res%cell
  res%volume = abs( c(1,1)*c(2,2)*c(3,3) - c(1,1)*c(2,3)*c(3,2) + &
                       c(1,2)*c(2,3)*c(3,1) - c(1,2)*c(2,1)*c(3,3) + &
                       c(1,3)*c(2,1)*c(3,2) - c(1,3)*c(2,2)*c(3,1) )
  res%pressure = -res%stress(1,1)-res%stress(2,2)-res%stress(3,3)
  res%virial = sum(res%coords*res%force)

  ! Copy results to output structure
  if (present(result)) then
    result = res
  end if

  ! Select requested (scalar) results
  if (present(results)) then
    do iResult = 1,nResults
      select case (trim(request(iResult)))
        case ('energy')
          results(iResult) = res%energy
        case ('volume')
          results(iResult) = res%volume
        case ('pressure')
          results(iResult) = res%pressure
        case ('virial')
          results(iResult) = res%virial
        case ('maxForce')
          results(iResult) = sqrt(maxval(sum(res%force**2,1)))
        case ('avgForce')
          results(iResult) = sqrt(sum(res%force**2)/nAtoms)
      end select
    end do
  end if

end subroutine readResult

END MODULE jobList
