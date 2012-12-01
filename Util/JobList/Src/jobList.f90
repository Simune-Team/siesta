module jobList

  implicit none

PUBLIC:: &
  countJobs,  &! count number of jobs and lists
  runJobs,    &! runs a job list read from a file unit
  getResults   ! collect requested results into summary files

PRIVATE    ! nothing is declared public beyond this point

  ! Internal parameters
  character(len=*),parameter:: dataSuffix = '.fdf'
  character(len=*),parameter:: endSuffix = '.EIG'
  character(len=*),parameter:: copyFiles  = 'cp -f *.fdf *.psf queue.sh'//' '
  character(len=*),parameter:: defaultQueue = &
                                       './siesta < $jobName.fdf > $jobName.out'
  character(len=*),parameter:: defaultRequest(1) = (/'energy'/)
  integer,parameter:: maxLines = 1000   ! Max lines in job-list file
  integer,parameter:: maxWords = 100    ! Max words in one line
  integer,parameter:: maxJobs = 1000    ! Max jobs in job list
  integer,parameter:: ll      = 500     ! Max characters per line
  integer,parameter:: wl      = 500     ! Max characters per word
  integer,parameter:: unitIn = 42       ! I/O unit for input files
  integer,parameter:: unitOut = 43      ! I/O unit for output files
  integer,parameter:: dp = kind(1.d0)

  ! Internal variables and arrays
  integer,           save:: myUnit, totJobs, totLists

  ! Derived type to hold calculation results
  type resultsType
    private
    real(dp)        :: cell(3,3)   =0._dp
    real(dp)        :: stress(3,3) =0._dp
    real(dp)        :: volume      =0._dp
    real(dp)        :: energy      =0._dp
    real(dp)        :: pressure    =0._dp
    real(dp)        :: virial      =0._dp 
    real(dp),pointer:: coord(:,:)  =>null()
    real(dp),pointer:: force(:,:)  =>null()
    integer, pointer:: za(:)       =>null()   ! atomic numbers
  end type resultsType

CONTAINS

!------------------------------------------------------------------------------

subroutine countJobs( unit, nJobs, nLists )

  implicit none
  integer,         intent(in) :: unit   ! IO unit of datafile
  integer,optional,intent(out):: nJobs  ! total number of jobs
  integer,optional,intent(out):: nLists ! total number of job lists

  integer :: nc
  character(len=wl):: myDir

  ! Copy input file unit for use across the module
  myUnit = unit

  ! Get current directory and add trailing '/' if necessary
  call getcwd(myDir)
  nc = len(trim(myDir))
  if (myDir(nc:nc) /= '/') myDir = myDir(1:nc) // '/'

  ! Initializations
  totJobs = 0             ! total number of jobs
  totLists = 0            ! total number of lists and sublists

  ! Scan the 'root' list (the whole file), specified by blank name
  call scanList(myDir,defaultQueue,defaultRequest,' ','count')

  ! Copy number of jobs and lists to output variables
  if (present(nJobs))  nJobs = totJobs
  if (present(nLists)) nLists = totLists

end subroutine countJobs

!------------------------------------------------------------------------------

subroutine runJobs( unit )

  implicit none
  integer,         intent(in) :: unit   ! IO unit of datafile

  integer :: nc
  character(len=wl):: myDir

  ! Copy input file unit for use across the module
  myUnit = unit

  ! Get current directory and add trailing '/' if necessary
  call getcwd(myDir)
  nc = len(trim(myDir))
  if (myDir(nc:nc) /= '/') myDir = myDir(1:nc) // '/'

  ! Initializations
  totJobs = 0             ! total number of jobs
  totLists = 0            ! total number of lists and sublists

  ! Scan the 'root' list (the whole file), specified by blank name
  call scanList(myDir,defaultQueue,defaultRequest,' ','run')

end subroutine runJobs

!------------------------------------------------------------------------------

subroutine getResults( unit )

  implicit none
  integer,         intent(in) :: unit   ! IO unit of datafile

  integer :: nc
  character(len=wl):: myDir

  ! Copy input file unit for use across the module
  myUnit = unit

  ! Get current directory and add trailing '/' if necessary
  call getcwd(myDir)
  nc = len(trim(myDir))
  if (myDir(nc:nc) /= '/') myDir = myDir(1:nc) // '/'

  ! Initializations
  totJobs = 0             ! total number of jobs
  totLists = 0            ! total number of lists and sublists

  ! Scan the 'root' list (the whole file), specified by blank name
  call scanList(myDir,defaultQueue,defaultRequest,' ','get')

end subroutine getResults

!------------------------------------------------------------------------------

recursive subroutine scanList( dir, queue, request, listName, task )

  implicit none
  character(len=*),intent(in) :: dir        ! parent directory
  character(len=*),intent(in) :: queue      ! queuing statement
  character(len=*),intent(in) :: request(:) ! requested results
  character(len=*),intent(in) :: listName   ! name of job list
  character(len=*),intent(in) :: task       ! ('count'|'run'|'get')

  integer :: iCase, iLine, iostat, nCases, nResults, nWords, nJobs
  character(len=wl):: caseDir(maxJobs), caseName(maxJobs), caseType(maxJobs), &
                       fileIn, fileOut, line, myDir, myRequest(maxWords), &
                       myQueue, newList, words(maxWords)
  real(dp):: results(maxWords)
  logical :: finished

  ! Set the list directory
  if (listName==' ') then  ! this is the 'root' list (the whole file)
    myDir = dir
  else                     ! this is a genuine list of jobs
    myDir = trim(dir) // trim(listName) // '/'
    ! Create a new directory for this list and copy files to it
    if (task=='run') then
      call system('mkdir -p ' // trim(myDir))
      call system(copyFiles // trim(myDir) )
      call chdir(trim(myDir))
    endif
  endif ! (listName==' ')

  ! Set my own copies of queue and request
  nResults = count(request/=' ')
  myRequest = ' '
  myRequest(1:nResults) = request(1:nResults)
  myQueue = queue

  ! Loop on lines of datafile
  finished = .false.       ! has the list terminated normally?
  nJobs = 0                ! number of jobs in list
  nCases = 0               ! number of 'cases' (jobs or sublists) in list
  do iLine = 1,maxLines

    ! Read one line of input file
    read(myUnit,'(a)',iostat=iostat) line

    ! End of file check
    if (iostat<0) then                      ! end of file
      if (listName==' ') then               ! this is the 'root' list
        finished = .true.
        exit ! iLine loop
      else
        print*,'runList ERROR: unterminated %list = '//trim(listName)
        stop
      endif ! (listName==' ')
    endif ! (iostat<0)

    ! Parse line, ignoring comments
    call parser(line,words,nWords)

    ! Act depending on line content
    if (nWords==0) then                     ! blank or comment line
      cycle
    elseif (trim(words(1))=='%queue') then  ! queue specification
      line = adjustl(line)
      myQueue = adjustl(line(7:))           ! remove '%queue' from line
    elseif (trim(words(1))=='%result') then ! result-request especification
      myRequest = ' '
      myRequest(1:nWords-1) = words(2:nWords) ! exclude 1st word (it is %result)
    elseif (words(1)=='%list') then         ! new sub-list
      nCases = nCases+1
      newList = words(2)
      caseType(nCases) = 'list'
      caseName(nCases) = newList
      caseDir(nCases) = trim(myDir) // trim(newList) // '/'
      call scanList(myDir,myQueue,myRequest,newList,task)
    elseif (words(1)=='%endlist') then      ! end of my list
      if (words(2)==listName) then
        finished = .true.
        exit ! iLine loop
      else
        print*,'runList ERROR: inconsistent %list = '//trim(listName)// &
               ', %endList = '//trim(words(2))
        stop
      endif
    else                                    ! new job
      nCases = nCases+1
      nJobs = nJobs+1
      caseType(nCases) = 'job'
      call nameJob( line, caseName(nCases) )
      caseDir(nCases) = trim(myDir) // trim(caseName(nCases)) // '/'
      if (task=='run') call runOneJob(myDir,myQueue,line)
    endif

  end do ! iLine
  if (.not.finished) stop 'runList ERROR: parameter maxLines too small'

  ! Write file of requested results
  if (task=='get') then
    if (listName==' ') then
      fileOut = trim(myDir)//'jobList.results'
    else
      fileOut = trim(myDir)//trim(listName)//'.results'
    endif
    call system('rm -f '//trim(fileOut))  ! clear output file, if it exists
    do iCase = 1,nCases
      if (caseType(iCase)=='list') then
        fileIn = trim(caseDir(iCase))//trim(caseName(iCase))//'.results'
        call system('cat '//trim(fileIn)//' >> '//trim(fileOut))
        call system('echo '' '' >> '//trim(fileOut))
      else
        call readResult( caseDir(iCase), myRequest, caseName(iCase), &
                         results=results )
        open(unitOut,file=trim(fileOut),position='append',status='unknown')
        write(unitOut,'(10e18.9)') results(1:nResults) 
        close(unitOut)
      endif
    end do
  endif ! (task=='get')

  totJobs = totJobs+nJobs
  totLists = totLists+1

end subroutine scanList

!------------------------------------------------------------------------------

subroutine nameJob( jobLine, jobName )

  implicit none
  character(len=*),intent(in) :: jobLine   ! line of job specification
  character(len=*),intent(out):: jobName   ! job name

  integer:: iWord, nc, nWords
  character(len=wl):: line, word, words(maxWords)

  ! Parse job line
  line = jobLine
  call parser(line,words,nWords)

  ! Set job name by concatenating words
  jobName = ' '
  do iWord = 1,nWords
    word = words(iWord)
    nc = len(trim(word))                         ! number of characters in word
    if (word(nc-3:nc)=='.fdf') word=word(1:nc-4) ! remove '.fdf' from word
    jobName = trim(jobName)//'_'//word           ! add word to job name
  end do ! iWord
  jobName=jobName(2:)                            ! remove leading '_'

end subroutine nameJob

!------------------------------------------------------------------------------

subroutine runOneJob( dir, queue, jobLine )

  implicit none
  character(len=*),intent(in) :: dir       ! work directory
  character(len=*),intent(in) :: queue     ! queuing statement
  character(len=*),intent(in) :: jobLine   ! line of job specification

  logical:: jobEnded
  integer:: ic, iLine, iWord, jc, jWord, nc, nLines, nWords
  character(len=wl):: jobDir, jobName, line(maxWords), &
                       queueJob, word, words(maxWords)

  ! Find job name
  call nameJob( jobLine, jobName )

  ! Create a new directory for this job and copy files to it
  jobDir = trim(dir) // trim(jobName) // '/'
  call system('mkdir -p ' // trim(jobDir))
  call system(copyFiles//trim(jobDir))

  ! Parse job line (using line(1) because parser argument is inout)
  line(1) = jobLine
  call parser(line(1),words,nWords)

  ! Convert job specifications to fdf format lines
  nLines = 0
  do iWord = 1,nWords
    nLines = nLines+1
    word = words(iWord)
    nc = len(trim(word))               ! number of characters in word
    if (word(nc-3:nc)=='.fdf') then    ! include new .fdf file
      line(nLines) = '%include '//trim(word)
    else                               ! make a line with rest of words
      line(nLines) = word
      do jWord = iWord+1,nWords
        line(nLines) = trim(line(nLines))//' '//words(jWord)
      end do
      exit ! do iWord loop
    end if
  end do ! iWord

  ! Write .fdf file for job
  open(unitOut,file=trim(jobDir)//trim(jobName)//'.fdf')
  write(unitOut,*) 'SystemLabel ',trim(jobName)   ! this is specific for siesta
  do iLine = nLines,1,-1               ! last lines (words) have priority
    write(unitOut,'(a)') trim(line(iLine))
  end do
  close(unitOut)

  ! Set job
  queueJob = queue
  nc = len(trim(queueJob))
  jc = len('$jobName')
  do ic = nc-jc+1,1,-1
    if (queueJob(ic:ic+jc-1)=='$jobName') &
      queueJob(ic:) = trim(jobName)//queueJob(ic+jc:)
  end do
!  print*, 'runOneJob: queueJob = ',trim(queueJob)

  ! Submit job only if directory does not contain terminated results already
  ! This is intended to re-run a whole list for failed jobs
  call chdir(trim(jobDir))
  inquire(file=trim(jobName)//endSuffix,exist=jobEnded)
  if (.not.jobEnded) call system(trim(queueJob))
  call chdir(trim(dir))

end subroutine runOneJob

!------------------------------------------------------------------------------

subroutine parser( line, words, nWords )

  implicit none
  character(len=*),intent(inout):: line      ! returns without comments
  character(len=*),intent(out)  :: words(:)  ! words in line, without comments
  integer,         intent(out)  :: nWords    ! number of words

  character(len=ll):: myLine
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

subroutine readResult( dir, request, name, result, results )

  implicit none
  character(len=*),          intent(in) :: dir         ! job directory
  character(len=*),          intent(in) :: request(:)  ! requested resuts
  character(len=*),          intent(in) :: name        ! job name
  type(resultsType),optional,intent(out):: result      ! job-results structure
  real(dp),         optional,intent(out):: results(:)  ! requested job results

  integer :: ia, ic, is, iostat, iResult, nAtoms, nResults, za
  real(dp):: c(3,3)
  character(len=wl):: fileName
  type(resultsType) :: res

  ! Read final structure
  fileName = trim(dir)//trim(name)//'.XV'
  open(unitIn,file=trim(fileName),status='old',iostat=iostat)
!  print'(a,i3,2x,a)','readResult: iostat,fileName=', iostat, trim(fileName)
  if (iostat==0) then
    do ic = 1,3
      read(unitIn,*) res%cell(:,ic)
    end do
    read(unitIn,*) nAtoms
    allocate( res%za(nAtoms), res%coord(3,nAtoms) )
    do ia = 1,nAtoms
      read(unitIn,*) is, res%za(ia), res%coord(:,ia)
    end do
  else
    nAtoms = 0
    allocate( res%za(nAtoms), res%coord(3,nAtoms) )
    res%cell = 0
    res%coord = 0
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
    deallocate(res%za)
    allocate(res%force(3,nAtoms),res%za(nAtoms))
    do ia = 1,nAtoms
      read(unitIn,*) is, res%za(ia), res%force(:,ia)
    end do
  else
    allocate( res%force(3,0) )
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
  res%virial = sum(res%coord*res%force)

  ! Copy results to output structure
  if (present(result)) then
    result = res
  end if

  ! Select requested (scalar) results
  if (present(results)) then
    nResults = count(request/=' ')
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
          results(iResult) = sqrt(sum(res%force**2)/max(1,nAtoms))
      end select
    end do
  end if

  ! Deallocate internal pointers
  if (associated(res%za)) deallocate(res%za)
  if (associated(res%coord)) deallocate(res%coord)
  if (associated(res%force)) deallocate(res%force)

end subroutine readResult

END MODULE jobList
