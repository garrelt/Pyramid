!>
!! \brief This module contains routines for calculating the ionization and temperature evolution of the entire grid (3D).
!!
!! Module for Capreole / C2-Ray (f90)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date:
!!
!! \b Version: 3D, no OpenMP

module evolve

  ! This module contains routines having to do with the calculation of
  ! the ionization evolution of the entire grid (3D).
     
  ! This version has been adapted for efficiency in order to be able
  ! to calculate large meshes.
    
  ! - evolve3D : step through grid
  ! - get_rate : take care of one grid point
  ! - evolution_of_cell: update entire grid

  use precision, only: dp
  use my_mpi
  use file_admin, only: logf,timefile,iterdump, results_dir
  use clocks, only: timestamp_wallclock
  use sizes, only: Ndim, mesh
  use grid, only: x_coordinate, y_coordinate, z_coordinate, cell_vol, cell_size
  use material, only: number_density_array, xHI_array, xHII_array, temperature_array
  use sourceprops, only: NumSrc, srcpos, NormFlux
  use radiation, only: NumFreqBnd
  use cosmology, only: time2zred
  use photonstatistics, only: state_before, calculate_photon_statistics, &
                              photon_loss, LLS_loss, report_photonstatistics, &
                              state_after, total_rates, total_ionizations, &
                              update_grandtotal_photonstatistics
  use c2ray_parameters, only: convergence_fraction, max_subbox, S_star_nominal



  implicit none

  save

  private

  public :: evolve3D, photoionization_array, heating_array, evolve_ini

  ! Periodic boundary conditions, has to be true for this version
  logical, parameter :: periodic_bc = .true.

  ! Minimum number of MPI processes for using the master-slave setup 
  integer, parameter ::  min_numproc_master_slave=10

  ! Variables connected to checking converge of global average ionization
  ! fraction
  ! Sum of intermediate ionization fraction xh_intermed(*,1)
  ! (used for convergence checking)
  real(kind=dp) :: current_sum_xHII 
  ! Previous value of this sum
  real(kind=dp) :: prev_sum_xHII
  ! Relative change in this sum (between iteration steps)
  real(kind=dp) :: delta_sum_xHII

  ! Grid variables

  !> H Photo-ionization rate on the entire grid
  real(kind=dp),dimension(:,:,:),allocatable :: photoionization_array
  real(kind=dp),dimension(:,:,:),allocatable :: heating_array

  !> Time-averaged H ionization fraction
  real(kind=dp),dimension(:,:,:),allocatable :: intermediate_average_xHI
  real(kind=dp),dimension(:,:,:),allocatable :: intermediate_average_xHII
  !> Intermediate result for H ionization fraction
  real(kind=dp),dimension(:,:,:),allocatable,public :: intermediate_end_xHI
  real(kind=dp),dimension(:,:,:),allocatable,public :: intermediate_end_xHII
  !> Intermediate result for Temperature
  real(kind=dp),dimension(:,:,:),allocatable,public :: intermediate_end_temperature
  real(kind=dp),dimension(:,:,:),allocatable,public :: intermediate_average_temperature
  !> H0 Column density (outgoing)
  real(kind=dp),dimension(:,:,:),allocatable :: coldensh_out
  !> Buffer for MPI communication
  real(kind=dp),dimension(:,:,:),allocatable :: buffer
  !> Photon loss from the grid
  real(kind=dp) :: photon_loss_all(1:NumFreqBnd)
  !> Photon loss from one source
  real(kind=dp),dimension(:,:),allocatable :: photon_loss_src_thread
  real(kind=dp) :: photon_loss_src(1:NumFreqBnd)

  ! mesh positions of end points for RT
  integer,dimension(Ndim) :: last_l !< mesh position of left end point for RT
  integer,dimension(Ndim) :: last_r !< mesh position of right end point for RT

  ! GM/121127: This variable should always be set. If not running OpenMP
  ! it should be equal to 1. We initialize it to 1 here.
  integer :: tn=1 !< thread number

contains

  ! =======================================================================

  !> Allocate the arrays needed for evolve
  subroutine evolve_ini ()
    
    allocate(photoionization_array(mesh(1),mesh(2),mesh(3)))
    allocate(heating_array(mesh(1),mesh(2),mesh(3)))
    allocate(intermediate_average_xHI(mesh(1),mesh(2),mesh(3)))
    allocate(intermediate_average_xHII(mesh(1),mesh(2),mesh(3)))
    allocate(intermediate_end_xHI(mesh(1),mesh(2),mesh(3))) 
    allocate(intermediate_end_xHII(mesh(1),mesh(2),mesh(3)))
    allocate(intermediate_end_temperature(mesh(1),mesh(2),mesh(3)))
    allocate(intermediate_average_temperature(mesh(1),mesh(2),mesh(3)))
    allocate(coldensh_out(mesh(1),mesh(2),mesh(3)))
    allocate(buffer(mesh(1),mesh(2),mesh(3)))
    allocate(photon_loss_src_thread(1:NumFreqBnd,nthreads))

  end subroutine evolve_ini

  ! ===========================================================================

  !> Evolve the entire grid over a time step dt
  subroutine evolve3D (time,dt,restart)

    ! Calculates the evolution of the hydrogen ionization state
     
    ! Author: Garrelt Mellema
     
    ! Date: 28-Feb-2008 (21-Aug-2006 (f77/OMP: 13-Jun-2005))

    ! Version: Multiple sources / Using average fractions to converge
    ! loop over sources
    
    ! History:
    ! 11-Jun-2004 (GM) : grid arrays now passed via common (in grid.h)
    !    and material arrays also (in material.h).
    ! 11-Jun-2004 (GM) : adapted for multiple sources.
    !  3-Jan-2005 (GM) : reintegrated with updated Ifront3D
    ! 20-May-2005 (GM) : split original eveolve0D into two routines
    ! 13-Jun-2005 (HM) : OpenMP version : Hugh Merz
    ! 21-Aug-2006 (GM) : MPI parallelization over the sources (static model).
    ! 28-Feb-2008 (GM) : Added master-slave model for distributing
    !                    over the processes. The program decides which
    !                    model to use.

    ! The time step
    real(kind=dp),intent(in) :: time !< time 
    real(kind=dp),intent(in) :: dt !< time step
    integer,intent(in) :: restart !< restart flag

    ! Loop variables
    integer :: iteration_counter  ! iteration counter

    ! Wall clock counting
    ! 8 bytes to beat the maxcount
    integer(kind=8) :: wallclock1
    integer(kind=8) :: wallclock2
    integer(kind=8) :: countspersec

    ! Flag variable (passed back from evolution_of_cell)
    integer :: convergence_failure

    ! Minimum number of cells which are allowed to be non-converged
    integer :: conv_criterion 

#ifdef MPI
    integer :: mympierror
#endif

    ! End of declarations

    ! Initialize wall clock counter (for dumps)
    call system_clock(wallclock1)

     ! Initial state (for photon statistics)
    call state_before (xHI_array, xHII_array)

    ! initialize average and intermediate results to initial values
    if (restart == 0) then

       intermediate_average_xHI(:,:,:) = xHI_array(:,:,:)
       intermediate_average_xHII(:,:,:) = xHII_array(:,:,:)	   
       intermediate_end_xHI(:,:,:) = xHI_array(:,:,:)
       intermediate_end_xHII(:,:,:) = xHII_array(:,:,:)
       intermediate_end_temperature(:,:,:) = temperature_array(:,:,:)
       intermediate_average_temperature(:,:,:) = temperature_array(:,:,:)
	   
       iteration_counter = 0 ! iteration starts at zero
       convergence_failure = mesh(1)*mesh(2)*mesh(3) ! initialize non-convergence 
       prev_sum_xHII = mesh(1)*mesh(2)*mesh(3) ! initialize non-convergence 
       delta_sum_xHII = 1.0 ! initialize non-convergence 
    else
       ! Reload xh_av,xh_intermed,photon_loss,iteration_counter
       call start_from_dump(restart,iteration_counter)
       call evolution_of_grid (convergence_failure,dt)
    endif

    ! Set the conv_criterion, if there are few sources we should make
    ! sure that things are converged around these sources.
    conv_criterion = min(int(convergence_fraction*mesh(1)*mesh(2)*mesh(3)),(NumSrc-1)/3)

    ! Report time
    if (rank == 0) write(timefile,"(A,F8.1)") &
         "Time before starting iteration: ", timestamp_wallclock ()

    ! Iterate to reach convergence for multiple sources
    do
       ! Update xh if converged and exit 
       current_sum_xHII = sum(intermediate_end_xHII(:,:,:))

       if (current_sum_xHII .gt. 0.0) then
          delta_sum_xHII = abs(current_sum_xHII-prev_sum_xHII)/current_sum_xHII
       else
          delta_sum_xHII = 1.0
       endif

       ! Convergence test
       if (convergence_failure .le. conv_criterion .or. delta_sum_xHII .lt. convergence_fraction) then

          xHI_array(:,:,:) = intermediate_end_xHI(:,:,:)
          xHII_array(:,:,:) = intermediate_end_xHII(:,:,:)
          temperature_array(:,:,:) = intermediate_end_temperature(:,:,:)
          ! Report
          if (rank == 0) then
             write(logf,*) "Multiple sources convergence reached"
             write(logf,*) "Test 1 values: ",convergence_failure, conv_criterion
             write(logf,*) "Test 2 values: ",delta_sum_xHII, &
                  convergence_fraction
          endif
          exit
       else
          if (iteration_counter > 100) then
             ! Complain about slow convergence
             if (rank == 0) write(logf,*) 'Multiple sources not converging'
             exit
          endif
       endif
       
       ! Save current value of mean ionization fraction
       prev_sum_xHII = current_sum_xHII

       ! Iteration loop counter
       iteration_counter = iteration_counter+1
       
       call total_rate_by_all_source (dt)

       if (rank == 0) then
          call system_clock(wallclock2,countspersec)
          ! Write iteration dump if more than 15 minutes have passed.
          ! system_clock starts counting at 0 when it reaches
          ! a max value. To catch this, test also for negative
          ! values of wallclock2-wallclock1
          write(logf,*) "Time and limit are: ", &
               wallclock2-wallclock1, 15.0*60.0*countspersec
          if (wallclock2-wallclock1 > 15*60*countspersec .or. &
               wallclock2-wallclock1 < 0 ) then
             call write_iteration_dump(iteration_counter)
             wallclock1=wallclock2
          endif
       endif

       call evolution_of_grid (convergence_failure,dt)

       ! Report time
       if (rank == 0) write(timefile,"(A,I3,A,F8.1)") &
            "Time after iteration ",iteration_counter," : ", timestamp_wallclock ()
    enddo

    ! Calculate photon statistics
    call calculate_photon_statistics (dt,xHI_array,xHII_array,intermediate_average_xHI,intermediate_average_xHII)
    call report_photonstatistics (dt)
    call update_grandtotal_photonstatistics (dt)

  end subroutine evolve3D

  ! ===========================================================================

  subroutine write_iteration_dump (iteration_counter)

    integer,intent(in) :: iteration_counter  ! iteration counter

    integer :: ndump=0
    
    character(len=20) :: iterfile

    ! Report time
    write(timefile,"(A,F8.1)") &
         "Time before writing iterdump: ", timestamp_wallclock ()

    ndump=ndump+1
    if (mod(ndump,2) == 0) then
       iterfile="iterdump2.bin"
    else
       iterfile="iterdump1.bin"
    endif

    open(unit=iterdump,file=iterfile,form="unformatted", &
         status="unknown")

    write(iterdump) iteration_counter,prev_sum_xHII
    write(iterdump) photon_loss_all
    write(iterdump) photoionization_array
    write(iterdump) heating_array
    write(iterdump) intermediate_average_xHI
    write(iterdump) intermediate_average_xHII
    write(iterdump) intermediate_end_xHI
    write(iterdump) intermediate_end_xHII
    write(iterdump) intermediate_average_temperature
    write(iterdump) intermediate_end_temperature

    close(iterdump)

    ! Report time
    write(timefile,"(A,F8.1)") "Time after writing iterdump: ", timestamp_wallclock ()

  end subroutine write_iteration_dump

  ! ===========================================================================

  subroutine start_from_dump(restart,iteration_counter)

    integer,intent(in) :: restart  ! restart flag
    integer,intent(out) :: iteration_counter  ! iteration counter

    character(len=20) :: iterfile

#ifdef MPI
    integer :: mympierror
#endif

    if (restart == 0) then
       if (rank == 0) &
            write(logf,*) "Warning: start_from_dump called incorrectly"
    else
       if (rank == 0) then

          ! Report time
          write(timefile,"(A,F8.1)") "Time before reading iterdump: ", timestamp_wallclock ()

          ! Set file to read (depending on restart flag)
          select case (restart)
          case (1) 
             iterfile="iterdump1.bin"
          case (2) 
             iterfile="iterdump2.bin"
          case (3) 
             iterfile="iterdump.bin"
          end select

          open(unit=iterdump,file=iterfile,form="unformatted",status="old")
       
          read(iterdump) iteration_counter,prev_sum_xHII
          read(iterdump) photon_loss_all
          read(iterdump) photoionization_array
		  read(iterdump) heating_array
          read(iterdump) intermediate_average_xHI
		  read(iterdump) intermediate_average_xHII
          read(iterdump) intermediate_end_xHI
		  read(iterdump) intermediate_end_xHII
		  read(iterdump) intermediate_average_temperature
		  read(iterdump) intermediate_end_temperature
          
          close(iterdump)
          write(logf,*) "Read iteration ",iteration_counter," from dump file"
          write(logf,*) 'photon loss counter: ',photon_loss_all
          write(logf,*) "Intermediate result for mean ionization fraction: ", &

          sum(intermediate_end_xHII(:,:,:))/real(mesh(1)*mesh(2)*mesh(3))

       endif
       
#ifdef MPI       
       ! Distribute the input parameters to the other nodes
       call MPI_BCAST(iteration_counter,1, &
            MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
       call MPI_BCAST(photon_loss_all,NumFreqBnd, &
            MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
       call MPI_BCAST(photoionization_array,mesh(1)*mesh(2)*mesh(3), &
            MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
	   call MPI_BCAST(heating_array,mesh(1)*mesh(2)*mesh(3), &
	        MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
       call MPI_BCAST(intermediate_average_xHI,mesh(1)*mesh(2)*mesh(3), &
            MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)       
       call MPI_BCAST(intermediate_average_xHII,mesh(1)*mesh(2)*mesh(3), &
            MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
	   call MPI_BCAST(intermediate_end_xHI,mesh(1)*mesh(2)*mesh(3), &
	        MPI_DOUBLE_PRECISION,0,&
	        MPI_COMM_NEW,mympierror)
	   call MPI_BCAST(intermediate_end_xHII,mesh(1)*mesh(2)*mesh(3), &
            MPI_DOUBLE_PRECISION,0,&
	        MPI_COMM_NEW,mympierror)
	   call MPI_BCAST(intermediate_average_temperature,mesh(1)*mesh(2)*mesh(3), &
	        MPI_DOUBLE_PRECISION,0,&
	        MPI_COMM_NEW,mympierror)
	   call MPI_BCAST(intermediate_end_temperature,mesh(1)*mesh(2)*mesh(3), &
            MPI_DOUBLE_PRECISION,0,&
	        MPI_COMM_NEW,mympierror)
#endif
       
       ! Report time
       write(timefile,"(A,F8.1)") &
            "Time after reading iterdump: ", timestamp_wallclock ()
       
    endif
  end subroutine start_from_dump
     
  ! ===========================================================================

  subroutine total_rate_by_all_source(dt)
    
    ! For random permutation of sources:
    use  m_ctrper, only: ctrper

    real(kind=dp),intent(in) :: dt  !< time step, passed on to evolve0D
    real(kind=dp) :: LLS_loss_all

#ifdef MPI
    integer :: mympierror
#endif

    if (rank == 0) write(logf,*) 'Doing all sources '
    ! reset global rates to zero for this iteration
    photoionization_array(:,:,:) = 0.0
    heating_array(:,:,:) = 0.0
    
    ! reset photon loss counters
    photon_loss(:) = 0.0
    LLS_loss = 0.0 ! make this a NumFreqBnd vector if needed later (GM/101129)

    ! Make a randomized list of sources :: call in serial
    ! disabled / GM110512
    !if ( rank == 0 ) call ctrper (SrcSeries(1:NumSrc),1.0)
    
#ifdef MPI
    ! Distribute the source list to the other nodes
    ! disabled / GM110512
    !call MPI_BCAST(SrcSeries,NumSrc,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
#endif
    
    ! Ray trace the whole grid for all sources.
    ! We can do this in two ways, depending on
    ! the number of processors. For many processors
    ! the master-slave setup should be more efficient.
    if (npr > min_numproc_master_slave) then
       call do_grid_master_slave (dt)
    else
       call do_grid_static (dt)
    endif
    
#ifdef MPI
    ! accumulate (sum) the MPI distributed photon losses
    call MPI_ALLREDUCE(photon_loss, photon_loss_all, NumFreqBnd, &
         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, mympierror)

    ! accumulate (sum) the MPI distributed photon losses
    call MPI_ALLREDUCE(LLS_loss, LLS_loss_all, 1, &
         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, mympierror)
    ! Put LLS_loss_all back in the LLS variable
    LLS_loss = LLS_loss_all

    ! accumulate (sum) the MPI distributed photoionization_array
    call MPI_ALLREDUCE(photoionization_array, buffer, mesh(1)*mesh(2)*mesh(3), &
         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, mympierror)
    
    ! Overwrite the processor local values with the accumulated value
    photoionization_array(:,:,:) = buffer(:,:,:)

    ! accumulate (sum) the MPI distributed heating_array
    call MPI_ALLREDUCE(heating_array, buffer, mesh(1)*mesh(2)*mesh(3), &
         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, mympierror)
    
    ! Overwrite the processor local values with the accumulated value
    heating_array(:,:,:) = buffer(:,:,:)

    photon_loss_all(:) = photon_loss(:)

#endif
    
  end subroutine total_rate_by_all_source

  ! ===========================================================================

  subroutine evolution_of_grid (convergence_failure,dt)

    ! Flag variable (passed back from evolution_of_cell)
    integer,intent(out) :: convergence_failure
    real(kind=dp),intent(in) :: dt  !< time step, passed on to evolve0D
	
    integer :: i,j,k  ! mesh position

    ! Mesh position of the cell being treated
    integer,dimension(Ndim) :: pos
	character(len=40) :: file1
    ! Report photon losses over grid boundary 
    if (rank == 0) write(logf,*) 'photon loss counter: ',photon_loss_all(:)
    
    ! Turn total photon loss into a mean per cell (used in evolution_of_cell)
    ! GM/110225: Possibly this should be 
    ! photon_loss_all(:)/(real(mesh(1))*real(mesh(2))
    ! Since this is the correct answer for an optically thin cell_volume
    ! Compromise: photon_loss_all(:)/(real(mesh(1))**3)*sum(xh_av(:,:,:,1))
    ! then if the box is fully ionized the optically thin 1/N^2 is used
    ! and if the box is more neutral something close to 1/N^3 is used.
    ! Not sure
    photon_loss(:)=photon_loss_all(:)/(real(mesh(1))*real(mesh(2))*real(mesh(3)))
 
    ! Report minimum value of xh_av(0) to check for zeros
    if (rank == 0) then

       write(logf,*) "min xh_av(0): ",minval(intermediate_average_xHI(:,:,:))

    endif
    
    ! Apply total photo-ionization rates from all sources (photoionization_array)
    convergence_failure = 0 ! will be used to check for convergence
    
    ! Loop through the entire mesh
    if (rank == 0) write(logf,*) 'Doing global '
    do k = 1,mesh(3)
       do j = 1,mesh(2)
          do i = 1,mesh(1)
             pos = (/ i,j,k /)
             call evolution_of_cell(dt,pos,convergence_failure)
          enddo
       enddo
    enddo
	
    ! Report on convergence and intermediate result
    if (rank == 0) then
       write(logf,*) "Number of non-converged points: ",convergence_failure
       write(logf,*) "Intermediate result for mean H ionization fraction: ", &
                     sum(intermediate_end_xHII(:,:,:))/real(mesh(1)*mesh(2)*mesh(3))			
    endif
    
    ! Report on photon conservation
    call calculate_photon_statistics (dt,intermediate_end_xHI,intermediate_end_xHII, &
                                      intermediate_average_xHI,intermediate_average_xHII)

    call report_photonstatistics (dt)
    
  end subroutine evolution_of_grid

  ! ===========================================================================
  
  !> Ray tracing the entire grid for all the sources using the
  !! master-slave model for distributing the sources over the
  !! MPI processes.
  subroutine do_grid_master_slave (dt)

    ! Ray tracing the entire grid for all the sources using the
    ! master-slave model for distributing the sources over the
    ! MPI processes.

    real(kind=dp),intent(in) :: dt  !< time step, passed on to evolve0D

	
    if (rank == 0) then
       call do_grid_master()
    else
       call do_grid_slave (dt)
    endif

  end subroutine do_grid_master_slave

  ! ===========================================================================

  !> The master task in the master-slave setup for distributing
  !! the ray-tracing over the sources over the MPI processes.
  subroutine do_grid_master ()

    ! The master task in the master-slave setup for distributing
    ! the ray-tracing over the sources over the MPI processes.
	
    integer :: source_index
    integer :: sources_done,whomax,who,answer
    ! counter for master-slave process
    integer,dimension(:),allocatable :: counts
#ifdef MPI
    integer :: mympierror
#endif

#ifdef MPI
    ! Source Loop - Master Slave with rank=0 as Master
    sources_done = 0
          
    source_index = 0
    
    ! Allocate counter for master-slave process
    if (.not.(allocated(counts))) allocate(counts(0:npr-1))

    ! send tasks to slaves 
    
    whomax = min(NumSrc,npr-1)
    do who = 1,whomax
       if (source_index <= NumSrc) then
          source_index=source_index+1
          call MPI_Send (source_index, 1, MPI_INTEGER, who, 1,  &
               MPI_COMM_NEW, mympierror)
       endif
    enddo
    
    do while (sources_done < NumSrc)
       
       ! wait for an answer from a slave. 
       
       call MPI_Recv (answer,     & ! address of receive buffer
            1,		   & ! number of items to receive
            MPI_INTEGER,	   & ! type of data
            MPI_ANY_SOURCE,  & ! can receive from any other
            1,		   & ! tag
            MPI_COMM_NEW,	   & ! communicator
            mympi_status,	   & ! status
            mympierror)
       
       who = mympi_status(MPI_SOURCE) ! find out who sent us the answer
       sources_done = sources_done+1 ! and the number of sources done
       
       ! put the slave on work again,
       ! but only if not all tasks have been sent.
       ! we use the value of num to detect this */
       if (source_index < NumSrc) then
          source_index = source_index+1
          call MPI_Send (source_index, 1, MPI_INTEGER, &
               who,		&	
               1,		&	
               MPI_COMM_NEW, &
               mympierror)
       endif
    enddo
    
    ! Now master sends a message to the slaves to signify that they 
    ! should end the calculations. We use a special tag for that:
    
    do who = 1,npr-1
       call MPI_Send (0, 1, MPI_INTEGER, &
            who,			  &
            2,			  & ! tag 
            MPI_COMM_NEW,	          &
            mympierror)
       
       ! the slave will send to master the number of calculations
       ! that have been performed. 
       ! We put this number in the counts array.
       
       call MPI_Recv (counts(who), & ! address of receive buffer
            1,                & ! number of items to receive
            MPI_INTEGER,      & ! type of data 
            who,              & ! receive from process who 
            7,                & ! tag 
            MPI_COMM_NEW,     & ! communicator 
            mympi_status,     & ! status
            mympierror)
    enddo
    
    write(logf,*) 'Mean number of sources per processor: ', real(NumSrc)/real(npr-1)
    write(logf,*) 'Counted mean number of sources per processor: ', real(sum(counts(1:npr-1)))/real(npr-1)
    write(logf,*) 'Minimum and maximum number of sources ', 'per processor: ', &
                  minval(counts(1:npr-1)),maxval(counts(1:npr-1))
    flush(logf)

#endif

  end subroutine do_grid_master

  ! ===========================================================================

  !> The slave task in the master-slave setup for distributing
  !! the ray-tracing over the sources over the MPI processes.
  subroutine do_grid_slave(dt)

    ! The slave task in the master-slave setup for distributing
    ! the ray-tracing over the sources over the MPI processes.

    real(kind=dp),intent(in) :: dt  !< time step, passed on to evolve0D
	
    integer :: local_count
    integer :: source_index
#ifdef MPI
    integer :: mympierror
#endif

#ifdef MPI
    local_count = 0
    call MPI_Recv (source_index,  & ! address of receive buffer
         1,    & ! number of items to receive
         MPI_INTEGER,  & ! type of data
         0,		  & ! can receive from master only
         MPI_ANY_TAG,  & ! can expect two values, so
         ! we use the wildcard MPI_ANY_TAG 
         ! here
         MPI_COMM_NEW, & ! communicator
         mympi_status, & ! status
         mympierror)
    
    ! if tag equals 2, then skip the calculations
    
    if (mympi_status(MPI_TAG) /= 2) then
       do 
#ifdef MPILOG
          ! Report
          write(logf,*) 'Processor ',rank,' received: ',source_index
          write(logf,*) ' that is source ',source_index !SrcSeries(source_index)
          write(logf,*) ' at:',srcpos(:,source_index)
          flush(logf)
#endif
          ! Do the source at hand
          call do_source(dt,source_index)
          
          ! Update local counter
          local_count = local_count+1
          
#ifdef MPILOG
          ! Report ionization fractions
          ! strange to report ionization fractiosn at this point
          !write(logf,*) sum(intermediate_end_xHII(:,:,:))/real(mesh(1)*mesh(2)*mesh(3))
          !write(logf,*) sum(intermediate_average_xHII(:,:,:))/real(mesh(1)*mesh(2)*mesh(3))
          !write(logf,*) local_count
#endif
          ! Send 'answer'
          call MPI_Send (local_count, 1,  & ! sending one int 
               MPI_INTEGER, 0, & ! to master
               1,              & ! tag
               MPI_COMM_NEW,   & ! communicator
               mympierror)
          
          ! Receive new source number
          call MPI_Recv (source_index,     & ! address of receive buffer
               1,            & ! number of items to receive
               MPI_INTEGER,  & ! type of data
               0,            & ! can receive from master only
               MPI_ANY_TAG,  & !  can expect two values, so
               !  we use the wildcard MPI_ANY_TAG 
               !  here
               MPI_COMM_NEW, & ! communicator
               mympi_status, & ! status
               mympierror)
          
          ! leave this loop if tag equals 2
          if (mympi_status(MPI_TAG) == 2) then
#ifdef MPILOG
             write(logf,*) 'Stop received'
             flush(logf)
#endif
             exit 
          endif
       enddo
    endif
    
    ! this is the point that is reached when a task is received with
    ! tag = 2
    
    ! send the number of calculations to master and return
    
#ifdef MPILOG
    ! Report
    write(logf,*) 'Processor ',rank,' did ',local_count,' sources'
    flush(logf)
#endif
    call MPI_Send (local_count,  &
         1,           & 
         MPI_INTEGER, & ! sending one int
         0,           & ! to master
         7,           & ! tag
         MPI_COMM_NEW,& ! communicator
         mympierror)
#endif

  end subroutine do_grid_slave

  ! ===========================================================================

  !> Does the ray-tracing over the sources by distributing
  !! the sources evenly over the available MPI processes-
  subroutine do_grid_static (dt)

    ! Does the ray-tracing over the sources by distributing
    ! the sources evenly over the available MPI processes-
    
    real(kind=dp),intent(in) :: dt  !< time step, passed on to evolve0D
	
    integer :: source_index

    ! Source Loop - distributed for the MPI nodes
    do source_index=1+rank,NumSrc,npr
       call do_source(dt,source_index)
    enddo

  end subroutine do_grid_static

  ! ===========================================================================
  
  !> Does the ray-tracing over the entire 3D grid for one source.
  !! The number of this source in the current list is source_index.
  subroutine do_source(dt,source_index)

    ! Does the ray-tracing over the entire 3D grid for one source.
    ! The number of this source in the current list is source_index.

    real(kind=dp),intent(in) :: dt  !< time step, passed on to evolve0D
    integer, intent(in) :: source_index !< number of the source being done

    integer :: i_axis,i_plane,i_quadrant
    integer :: i,j,k
    integer :: nbox
    integer :: nnt

    ! Mesh position of the cell being treated
    integer,dimension(Ndim) :: ray_traced_pos
    
    character(len=40) :: file1
		
    ! reset column densities for new source point
    ! coldensh_out is unique for each source point
    coldensh_out(:,:,:) = 0.0
    
    ! Find the mesh position for the end points of the loop
    ! We trace until we reach max_subbox (set in c2ray_parameters)
    ! or the end of the grid. In the periodic case the end of the
    ! grid is always mesh/2 away from the source. If the grid is
    ! even-sized we trave mesh/2 cells to the left and mesh/2-1
    ! cell to the right. If it is odd, it is mesh/2 in either direction.
    ! The mod(mesh,2) takes care of handling this.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! this thing need to be changed in adaptive method
    if (periodic_bc) then
       last_r(:) = srcpos(:,source_index)+mesh(:)/2-1+mod(mesh(:),2)
       last_l(:) = srcpos(:,source_index)-mesh(:)/2
    else
       last_r(:) = min(srcpos(:,source_index)+max_subbox,mesh(:))
       last_l(:) = max(srcpos(:,source_index)-max_subbox,1)
    endif

    ! Loop through grid in the order required by 
    ! short characteristics
    
    ! Transfer is done in a set of cubes of increasing size.
    ! If the HII region is small we do not waste time calculating
    ! column densities of parts of the grid where no radiation
    ! penetrates. To test whether the current subbox is large
    ! enough we use the photon_loss_src. If this is non-zero,
    ! photons are leaving this subbox and we need to do another
    ! one. We also stop once we have done the whole grid.
    !nbox=0 ! subbox counter
    !photon_loss_src(:)=NormFlux(source_index)*S_star_nominal !-1.0 ! to pass the first while test
    !last_r(:)=srcpos(:,source_index) ! to pass the first while test
    !last_l(:)=srcpos(:,source_index) ! to pass the first while test

     photon_loss_src(:) = 0.0 ! reset photon_loss_src to zero
     photon_loss_src_thread(:,:) = 0.0 ! reset photon_loss_src to zero

     ! First do source point (on first pass)

      ray_traced_pos(:) = srcpos(:,source_index)
      call get_rate(dt,ray_traced_pos,source_index)

    ! do independent areas of the mesh in parallel using OpenMP
     !$omp parallel default(shared) private(tn)
     !!!reduction(+:photon_loss_src)

     ! Find out your thread number
#ifdef MY_OPENMP
     tn = omp_get_thread_num()+1
#else
     tn = 1
#endif

          ! Then do the the axes
          !$omp do schedule(dynamic,1)
          do i_axis = 1,6
             call ray_trace_1D_axis(dt,source_index,i_axis)
          enddo
          !$omp end do
          
          ! Then the source planes
          !$omp do schedule (dynamic,1)
          do i_plane = 1,12
             call ray_trace_2D_plane(dt,source_index,i_plane)
          end do
          !$omp end do
          
          ! Then the quadrants
          !$omp do schedule (dynamic,1)
          do i_quadrant = 1,8
             call ray_trace_3D_quadrant(dt,source_index,i_quadrant)
          end do
          !$omp end do
          
          !$omp end parallel

          ! Collect photon losses for each thread
          do nnt=1,nthreads
             photon_loss_src(:) = photon_loss_src(:) + photon_loss_src_thread(:,nnt)
          enddo
		
    ! Record the final photon loss, this is the photon loss that leaves
    ! the grid.
    photon_loss(:) = photon_loss(:) + photon_loss_src(:)

  end subroutine do_source

  ! ===========================================================================

  ! Ray tracing for the axes going through the source point
  ! should be called after having done the source point
  subroutine ray_trace_1D_axis(dt,source_index,i_axis)

    real(kind=dp),intent(in) :: dt      ! passed on to get_rate
    integer,intent(in) :: source_index           ! current source
    integer,intent(in) :: i_axis        ! axis to do
    integer :: i,j,k
    integer,dimension(Ndim) :: ray_traced_pos ! mesh position

    select case (i_axis)
    case(1)
       ! sweep in +i direction
       ray_traced_pos(2:3)=srcpos(2:3,source_index)
       do i=srcpos(1,source_index)+1,last_r(1)
          ray_traced_pos(1)=i
          call get_rate(dt,ray_traced_pos,source_index) !# `positive' i
       enddo
    case(2)
       ! sweep in -i direction
       ray_traced_pos(2:3)=srcpos(2:3,source_index)
       do i=srcpos(1,source_index)-1,last_l(1),-1
          ray_traced_pos(1)=i
          call get_rate(dt,ray_traced_pos,source_index) !# `negative' i
       end do
    case(3)
       ! sweep in +j direction
       ray_traced_pos(1)=srcpos(1,source_index)
       ray_traced_pos(3)=srcpos(3,source_index)
       do j=srcpos(2,source_index)+1,last_r(2)
          ray_traced_pos(2)=j
          call get_rate(dt,ray_traced_pos,source_index) !# `positive' j
       end do
    case(4)
       ! sweep in -j direction
       ray_traced_pos(1)=srcpos(1,source_index)
       ray_traced_pos(3)=srcpos(3,source_index)
       do j=srcpos(2,source_index)-1,last_l(2),-1
          ray_traced_pos(2)=j
          call get_rate(dt,ray_traced_pos,source_index) !# `negative' j
       end do
    case(5)
       ! sweep in +k direction
       ray_traced_pos(1:2)=srcpos(1:2,source_index)
       do k=srcpos(3,source_index)+1,last_r(3)
          ray_traced_pos(3)=k
          call get_rate(dt,ray_traced_pos,source_index) !# `positive' k
       end do
    case(6)
       ! sweep in -k direction
       ray_traced_pos(1:2)=srcpos(1:2,source_index)
       do k=srcpos(3,source_index)-1,last_l(3),-1
          ray_traced_pos(3)=k
          call get_rate(dt,ray_traced_pos,source_index) !# `negative' k
       end do
    end select
    
  end subroutine ray_trace_1D_axis

  ! ===========================================================================

  !> Ray tracing for planes containing the source point
  !! should be called after ray_trace_1D_axis
  subroutine ray_trace_2D_plane(dt,source_index,i_plane)

    ! find column density for the axes going through the source point
    ! should be called after having done the source point
    
    real(kind=dp),intent(in) :: dt      ! passed on to get_rate
    integer,intent(in) :: source_index           ! current source
    integer,intent(in) :: i_plane        ! plane to do

    integer :: i,j,k
    integer,dimension(Ndim) :: ray_traced_pos ! mesh position

    select case (i_plane)
    case(1)
       ! sweep in +i,+j direction
       ray_traced_pos(3)=srcpos(3,source_index)
       do j=srcpos(2,source_index)+1,last_r(2)
          ray_traced_pos(2)=j
          do i=srcpos(1,source_index)+1,last_r(1)
             ray_traced_pos(1)=i
             call get_rate(dt,ray_traced_pos,source_index)
          enddo
       enddo
    case(2)
       ! sweep in +i,-j direction
       ray_traced_pos(3)=srcpos(3,source_index)
       do j=srcpos(2,source_index)-1,last_l(2),-1
          ray_traced_pos(2)=j
          do i=srcpos(1,source_index)+1,last_r(1)
             ray_traced_pos(1)=i
             call get_rate(dt,ray_traced_pos,source_index)
          enddo
       enddo
    case(3)
       ! sweep in -i,+j direction
       ray_traced_pos(3)=srcpos(3,source_index)
       do j=srcpos(2,source_index)+1,last_r(2)
          ray_traced_pos(2)=j
          do i=srcpos(1,source_index)-1,last_l(1),-1
             ray_traced_pos(1)=i
             call get_rate(dt,ray_traced_pos,source_index)
          enddo
       enddo
    case(4)
       ! sweep in -i,-j direction
       ray_traced_pos(3)=srcpos(3,source_index)
       do j=srcpos(2,source_index)-1,last_l(2),-1
          ray_traced_pos(2)=j
          do i=srcpos(1,source_index)-1,last_l(1),-1
             ray_traced_pos(1)=i
             call get_rate(dt,ray_traced_pos,source_index)
          enddo
       enddo
    case(5)
       ! sweep in +i,+k direction
       ray_traced_pos(2)=srcpos(2,source_index)
       do k=srcpos(3,source_index)+1,last_r(3)
          ray_traced_pos(3)=k
          do i=srcpos(1,source_index)+1,last_r(1)
             ray_traced_pos(1)=i
             call get_rate(dt,ray_traced_pos,source_index)
          enddo
       enddo
    case(6)
       ! sweep in -i,+k direction
       ray_traced_pos(2)=srcpos(2,source_index)
       do k=srcpos(3,source_index)+1,last_r(3)
          ray_traced_pos(3)=k
          do i=srcpos(1,source_index)-1,last_l(1),-1
             ray_traced_pos(1)=i
             call get_rate(dt,ray_traced_pos,source_index)
          enddo
       enddo
    case(7)
       ! sweep in -i,-k direction
       ray_traced_pos(2)=srcpos(2,source_index)
       do k=srcpos(3,source_index)-1,last_l(3),-1
          ray_traced_pos(3)=k
          do i=srcpos(1,source_index)-1,last_l(1),-1
             ray_traced_pos(1)=i
             call get_rate(dt,ray_traced_pos,source_index)
          enddo
       enddo
    case(8)
       ! sweep in +i,-k direction
       ray_traced_pos(2)=srcpos(2,source_index)
       do k=srcpos(3,source_index)-1,last_l(3),-1
          ray_traced_pos(3)=k
          do i=srcpos(1,source_index)+1,last_r(1)
             ray_traced_pos(1)=i
             call get_rate(dt,ray_traced_pos,source_index)
          enddo
       enddo
    case(9) 
       ! sweep in +j,+k direction
       ray_traced_pos(1)=srcpos(1,source_index)
       do k=srcpos(3,source_index)+1,last_r(3)
          ray_traced_pos(3)=k
          do j=srcpos(2,source_index)+1,last_r(2)
             ray_traced_pos(2)=j
             call get_rate(dt,ray_traced_pos,source_index)
          enddo
       enddo
    case(10) 
       ! sweep in -j,+k direction
       ray_traced_pos(1)=srcpos(1,source_index)
       do k=srcpos(3,source_index)+1,last_r(3)
          ray_traced_pos(3)=k
          do j=srcpos(2,source_index)-1,last_l(2),-1
             ray_traced_pos(2)=j
             call get_rate(dt,ray_traced_pos,source_index)
          enddo
       enddo
    case(11) 
       ! sweep in +j,-k direction
       ray_traced_pos(1)=srcpos(1,source_index)
       do k=srcpos(3,source_index)-1,last_l(3),-1
          ray_traced_pos(3)=k
          do j=srcpos(2,source_index)+1,last_r(2)
             ray_traced_pos(2)=j
             call get_rate(dt,ray_traced_pos,source_index)
          enddo
       enddo
    case(12) 
       ! sweep in -j,-k direction
       ray_traced_pos(1)=srcpos(1,source_index)
       do k=srcpos(3,source_index)-1,last_l(3),-1
          ray_traced_pos(3)=k
          do j=srcpos(2,source_index)-1,last_l(2),-1
             ray_traced_pos(2)=j
             call get_rate(dt,ray_traced_pos,source_index)
          enddo
       enddo
       
    end select
    
  end subroutine ray_trace_2D_plane

  ! ===========================================================================

  !> Ray tracing for the 8 octants 
  !! should be called after ray_trace_2D_plane
  subroutine ray_trace_3D_quadrant(dt,source_index,i_quadrant)

    ! find column density for a z-plane srcpos(3) by sweeping in x and y
    ! directions
    
    real(kind=dp),intent(in) :: dt     ! passed on to get_rate
    integer,intent(in) :: source_index           ! current source
    integer,intent(in) :: i_quadrant    ! which quadrant to do    

    integer :: i,j,k
    integer,dimension(Ndim) :: ray_traced_pos ! mesh position

    select case (i_quadrant)
    case (1)
       ! sweep in +i,+j,+k direction
       do k=srcpos(3,source_index)+1,last_r(3)
          ray_traced_pos(3) = k
          do j=srcpos(2,source_index)+1,last_r(2)
             ray_traced_pos(2) = j
             do i=srcpos(1,source_index)+1,last_r(1)
                ray_traced_pos(1) = i
                call get_rate(dt,ray_traced_pos,source_index)
             end do
          enddo
       enddo
    case (2)
       ! sweep in -i,+j,+k direction
       do k=srcpos(3,source_index)+1,last_r(3)
          ray_traced_pos(3) = k
          do j=srcpos(2,source_index)+1,last_r(2)
             ray_traced_pos(2) = j
             do i=srcpos(1,source_index)-1,last_l(1),-1
                ray_traced_pos(1) = i
                call get_rate(dt,ray_traced_pos,source_index) !# `negative' i
             end do
          end do
       enddo
    case (3)
       ! sweep in +i,-j,+k direction
       do k=srcpos(3,source_index)+1,last_r(3)
          ray_traced_pos(3) = k
          do j=srcpos(2,source_index)-1,last_l(2),-1
             ray_traced_pos(2) = j
             do i=srcpos(1,source_index)+1,last_r(1)
                ray_traced_pos(1) = i
                call get_rate(dt,ray_traced_pos,source_index) !# `negative' i
             end do
          end do
       enddo
    case(4)
       ! sweep in -i,-j,+k direction
       do k=srcpos(3,source_index)+1,last_r(3)
          ray_traced_pos(3) = k
          do j=srcpos(2,source_index)-1,last_l(2),-1
             ray_traced_pos(2) = j
             do i=srcpos(1,source_index)-1,last_l(1),-1
                ray_traced_pos(1) = i
                call get_rate(dt,ray_traced_pos,source_index) !# `negative' i
             end do
          end do
       enddo
    case (5)
       ! sweep in +i,+j,-k direction
       do k=srcpos(3,source_index)-1,last_l(3),-1
          ray_traced_pos(3) = k
          do j=srcpos(2,source_index)+1,last_r(2)
             ray_traced_pos(2) = j
             do i=srcpos(1,source_index)+1,last_r(1)
                ray_traced_pos(1) = i
                call get_rate(dt,ray_traced_pos,source_index) !# `positive' i
             end do
          enddo
       enddo
    case (6)
       ! sweep in -i,+j,-k direction
       do k=srcpos(3,source_index)-1,last_l(3),-1
          ray_traced_pos(3) = k
          do j=srcpos(2,source_index)+1,last_r(2)
             ray_traced_pos(2) = j
             do i=srcpos(1,source_index)-1,last_l(1),-1
                ray_traced_pos(1) = i
                call get_rate(dt,ray_traced_pos,source_index) !# `negative' i
             end do
          end do
       enddo
    case (7)
       ! sweep in +i,-j,-k direction
       do k=srcpos(3,source_index)-1,last_l(3),-1
          ray_traced_pos(3) = k
          do j=srcpos(2,source_index)-1,last_l(2),-1
             ray_traced_pos(2) = j
             do i=srcpos(1,source_index)+1,last_r(1)
                ray_traced_pos(1) = i
                call get_rate(dt,ray_traced_pos,source_index) !# `negative' i
             end do
          end do
       enddo
    case(8)
       ! sweep in -i,-j,-k direction
       do k=srcpos(3,source_index)-1,last_l(3),-1
          ray_traced_pos(3) = k
          do j=srcpos(2,source_index)-1,last_l(2),-1
             ray_traced_pos(2) = j
             do i=srcpos(1,source_index)-1,last_l(1),-1
                ray_traced_pos(1) = i
                call get_rate(dt,ray_traced_pos,source_index) !# `negative' i
             end do
          end do
       enddo
    end select

  end subroutine ray_trace_3D_quadrant

  !=======================================================================

  !> Calculates the photo-ionization rate for one cell due to one source
  !! and adds this contribution to the collective rate.
  subroutine get_rate(dt,ray_traced_pos,source_index)
    
    ! Calculates the photo-ionization rate for one cell due to one source
    ! and adds this contribution to the collective rate.
    
    ! Author: Garrelt Mellema
    
    ! Date: 01-Feb-2008 (21-Aug-2006, 20-May-2005, 5-Jan-2005, 02 Jun 2004)
    
    ! Version: multiple sources, fixed temperature
    
    ! Multiple sources
    ! We call this routine for every grid point and for every source (ns).
    ! The photo-ionization rates for each grid point are found and added
    ! to photoionization_array, but the ionization fractions are not updated.
    ! For the first pass (iteration_counter = 1) it makes sense to DO update the
    ! ionization fractions since this will increase convergence speed
    ! in the case of isolated sources.

    use tped, only: electrondens
    use doric_module, only: doric, coldens
    use radiation, only: photoion, photrates
    use material, only: clumping_point !,coldensh_LLS
    use c2ray_parameters, only: epsilon,convergence1,convergence2, &
         type_of_clumping, convergence_frac,use_LLS,type_of_LLS
    use mathconstants, only: pi
    use material, only: coldensh_LLS, LLS_point
    use photonstatistics, only: total_LLS_loss

    ! column density for stopping chemisty
    real(kind=dp),parameter :: max_coldensh=1.587e21_dp
	!real(kind=dp),parameter :: max_coldensh=2e19_dp 
    
    logical :: falsedummy ! always false, for tests
    parameter(falsedummy=.false.)

    real(kind=dp),intent(in) :: dt ! time step
    integer,dimension(Ndim),intent(in) :: ray_traced_pos ! cell position (for RT)
    integer,intent(in)      :: source_index ! source number 

    integer :: nx,nd,idim ! loop counters
    integer,dimension(Ndim) :: pos
    integer,dimension(Ndim) :: srcpos1
    real(kind=dp) :: dist2,path,cell_vol_ph
    real(kind=dp) :: xs,ys,zs
    real(kind=dp) :: coldensh_in
    real(kind=dp) :: number_density_atom
    real(kind=dp) :: convergence
    
    type(photrates) :: phi

    !write(*,*) ray_traced_pos
    ! set convergence tolerance
    convergence = convergence1

    ! coordinate transformation. from out of box frame back to box frame, but still global frame
    do idim = 1,Ndim
       pos(idim) = modulo(ray_traced_pos(idim)-1,mesh(idim))+1
    enddo

    ! If coldensh_out is zero, we have not yet done this point
    ! yet, so do it. Otherwise do nothing.
    if (coldensh_out(pos(1),pos(2),pos(3)) == 0.0) then
   
    ! Initialize local density and temperature
    number_density_atom = number_density_array(pos(1),pos(2),pos(3))

    ! Find the column density at the entrance point of the cell (short
    ! characteristics)
       
       if ( all( ray_traced_pos(:) == srcpos(:,source_index) ) ) then
          ! Do not call short_characteristic for the source point.
          ! Set coldensh and path by hand
          coldensh_in = 0.0
          path = 0.5*cell_size(1)
          
          ! Find the distance to the source (average?)
          !dist=0.5*dr(1) NOT NEEDED         ! this makes cell_vol=dx*dy*dz
          !cell_vol_ph=4.0/3.0*pi*dist**3
          cell_vol_ph = cell_size(1)*cell_size(2)*cell_size(3)
          
       else
          
          ! For all other points call short_characteristic to find the column density and path
          call short_characteristic(ray_traced_pos,srcpos(:,source_index),coldensh_in,path)
          path = path*cell_size(1)
 

		  

		  
          ! Find the distance to the source
          xs = cell_size(1)*real(ray_traced_pos(1)-srcpos(1,source_index))
          ys = cell_size(2)*real(ray_traced_pos(2)-srcpos(2,source_index))
          zs = cell_size(3)*real(ray_traced_pos(3)-srcpos(3,source_index))
          dist2 = xs*xs+ys*ys+zs*zs
          
          ! Find the cell_volume of the shell this cell is part of 
          ! (dilution factor).
          cell_vol_ph=4.0*pi*dist2*path

       endif
       
       ! Only ray trace and exit. Do not touch the ionization
       ! fractions. They are updated using photoionization_array in evolution_of_cell
          
       ! Add the (time averaged) column density of this cell
       ! to the total column density (for this source)
       ! and add the LLS column density to this.
       ! GM/110224: No! This messes up phi since phi is based
       !  upon the difference between the in and out column density.
       !  Instead add the LLS to coldensh_in, see above
       coldensh_out(pos(1),pos(2),pos(3))=coldensh_in + &
            coldens(path,intermediate_average_xHI(pos(1),pos(2),pos(3)),number_density_atom) !+ coldensh_LLS * path/cell_size(1)
	  
       ! Calculate (photon-conserving) photo-ionization rate
       if (coldensh_in < max_coldensh) then
          call photoion(phi,coldensh_in,coldensh_out(pos(1),pos(2),pos(3)), &
               cell_vol_ph,source_index)
          phi%h=phi%h/(intermediate_average_xHI(pos(1),pos(2),pos(3))*number_density_atom)

          call total_LLS_loss(phi%h_in*cell_vol/cell_vol_ph, coldensh_LLS * path/cell_size(1))
       else
          phi%h=0.0
          phi%h_out=0.0
       endif
       
       ! Add photo-ionization rate to the global array 
       ! (applied in evolution_of_cell)
       photoionization_array(pos(1),pos(2),pos(3)) = photoionization_array(pos(1),pos(2),pos(3)) + phi%h
       
       ! Add photo-ionization rate to the global array 
       ! (applied in evolution_of_cell)
       heating_array(pos(1),pos(2),pos(3)) = heating_array(pos(1),pos(2),pos(3)) + phi%hv_h
				 
       ! Photon statistics: register number of photons leaving the grid
       if ( (any(ray_traced_pos(:) == last_l(:))) .or. &
            (any(ray_traced_pos(:) == last_r(:))) ) then
          ! GM/121127: Make sure that tn is always set, even when we
          ! are not running OpenMP. In that case tn should be 1.
          photon_loss_src_thread(1,tn)=photon_loss_src_thread(1,tn) + &
               phi%h_out*cell_vol/cell_vol_ph
			   
       endif

    endif ! end of coldens test

  end subroutine get_rate

  ! =======================================================================

  !> Calculates the evolution of the hydrogen ionization state for
  !! one cell (mesh position pos) and multiple sources.
  subroutine evolution_of_cell(dt,pos,convergence_failure)

    ! Calculates the evolution of the hydrogen ionization state for
    ! one cell (pos) and multiple sources.

    ! Author: Garrelt Mellema

    ! Date: 11-Feb-2008 (20-May-2005, 5-Jan-2005, 02 Jun 2004)
    
    ! Version: Multiple sources (global update, no ray tracing)

    ! Multiple sources
    ! Global update: the collected rates are applied and the new ionization 
    ! fractions and temperatures are calculated.
    ! We check for convergence.
    
    use c2ray_parameters, only: convergence1,convergence2,type_of_clumping, convergence_frac

    real(kind=dp),intent(in) :: dt ! time step
    integer,dimension(Ndim),intent(in) :: pos ! position on mesh
    integer,intent(inout) :: convergence_failure ! convergence counter

    integer :: nx,nit ! loop counters

    real(kind=dp) :: current_end_xHI ! ionization fractions
	real(kind=dp) :: current_end_xHII
    real(kind=dp) :: current_begin_xHI
	real(kind=dp) :: current_begin_xHII
	real(kind=dp) :: current_avg_xHI
	real(kind=dp) :: current_avg_xHII
	real(kind=dp) :: current_begin_temper
	real(kind=dp) :: current_avg_temper
	real(kind=dp) :: current_end_temper
	real(kind=dp) :: prev_avg_xHI
	real(kind=dp) :: prev_end_temper
    real(kind=dp) :: number_density_atom ! local number density

    real(kind=dp) :: photoionization_rate ! local H photo-ionization rate (only non-zero when local=.false.!)
	real(kind=dp) :: heating_rate ! local H heating-ionization rate 

    real(kind=dp) :: convergence
    
    ! Set convergence tolerance
    convergence=convergence2

    ! Initialize local ionization states to global ones

    current_begin_xHI = xHI_array(pos(1),pos(2),pos(3))
    current_begin_xHII = xHII_array(pos(1),pos(2),pos(3))  
    current_avg_xHI = intermediate_average_xHI(pos(1),pos(2),pos(3))
    current_avg_xHII = intermediate_average_xHII(pos(1),pos(2),pos(3))
    
    ! Initialize local scalars for density and temperature
    current_begin_temper = temperature_array(pos(1),pos(2),pos(3))
    current_avg_temper = intermediate_average_temperature(pos(1),pos(2),pos(3))
    current_end_temper = temperature_array(pos(1),pos(2),pos(3))

    number_density_atom = number_density_array(pos(1),pos(2),pos(3))
	
    ! Use the collected photo-ionization rates
    photoionization_rate = photoionization_array(pos(1),pos(2),pos(3)) 
    heating_rate = heating_array(pos(1),pos(2),pos(3))
	
    call solve_ionization_equation (dt, number_density_atom, current_end_xHI, current_end_xHII, &
                                    current_avg_xHI, current_avg_xHII, current_begin_xHI, current_begin_xHII, &
                                    current_end_temper, current_avg_temper, current_begin_temper, &
                                    photoionization_rate, heating_rate, 0.0_dp, 0.0_dp, 0.0_dp, pos)

    ! Test for global convergence using the time-averaged neutral fraction.
    ! For low values of this number assume convergence

    prev_avg_xHI = intermediate_average_xHI(pos(1),pos(2),pos(3)) ! use previously calculated xh_av
    prev_end_temper = intermediate_end_temperature(pos(1),pos(2),pos(3))

    if (  (abs(current_avg_xHI-prev_avg_xHI) > convergence2 .and. &
          abs((current_avg_xHI-prev_avg_xHI)/current_avg_xHI) > convergence2 .and. &
          current_avg_xHI > convergence_frac) .or. &
		  (abs(current_end_temper-prev_end_temper)/prev_end_temper > convergence2	&	  
		  ) ) then 

       convergence_failure = convergence_failure+1
    endif

    ! Copy ion fractions to the global arrays.
    intermediate_end_xHI(pos(1),pos(2),pos(3))=current_end_xHI
    intermediate_end_xHII(pos(1),pos(2),pos(3))=current_end_xHII
    intermediate_average_xHI(pos(1),pos(2),pos(3))=current_avg_xHI
    intermediate_average_xHII(pos(1),pos(2),pos(3))=current_avg_xHII
    intermediate_end_temperature(pos(1),pos(2),pos(3))=current_end_temper
    intermediate_average_temperature(pos(1),pos(2),pos(3))=current_avg_temper
	
  end subroutine evolution_of_cell

  ! ===========================================================================

  subroutine solve_ionization_equation (dt, number_density_atom, end_xHI, end_xHII, avg_xHI, avg_xHII, &
                                        begin_xHI, begin_xHII, end_temper, avg_temper, begin_temper, &
                                        photoionization_rate, heating_rate, coldensh_in, path, cell_vol_ph, pos)

    use c2ray_parameters, only: convergence1,convergence2,type_of_clumping, & 
                                convergence_frac, add_photon_losses
    use tped, only: electrondens
    use doric_module, only: doric, coldens
    use material, only: clumping_point
    use radiation, only: photoion, photrates
    use thermalevolution, only: thermal
	
    real(kind=dp),intent(in) :: dt !< time step
    real(kind=dp),intent(in) :: number_density_atom
    real(kind=dp),intent(out) :: end_xHI
    real(kind=dp),intent(out) :: end_xHII
    real(kind=dp),intent(inout) :: avg_xHI
    real(kind=dp),intent(inout) :: avg_xHII
    real(kind=dp),intent(in) :: begin_xHI
    real(kind=dp),intent(in) :: begin_xHII
    real(kind=dp),intent(inout) :: end_temper
    real(kind=dp),intent(inout) :: avg_temper
    real(kind=dp),intent(in) :: begin_temper
    real(kind=dp),intent(in) :: photoionization_rate !< local photoionization rate
    real(kind=dp),intent(in) :: heating_rate !< local heating rate
    real(kind=dp),intent(in) :: coldensh_in
    real(kind=dp),intent(in) :: path
    real(kind=dp),intent(in) :: cell_vol_ph
    integer,dimension(Ndim),intent(in) :: pos !< position on mesh

    real(kind=dp) :: photoionization_rate_cell
    real(kind=dp) :: heating_rate_cell
    real(kind=dp) :: prev_avg_xHI
    real(kind=dp) :: prev_temper
    real(kind=dp) :: number_density_electron
    real(kind=dp) :: coldensh_cell

    integer :: nit

    type(photrates) :: phi

    ! Initialize yh to initial value
    end_xHI = begin_xHI
    end_xHII = begin_xHII

    ! Initialize local clumping (if type of clumping is appropriate)
    if (type_of_clumping == 5) call clumping_point (pos(1),pos(2),pos(3))
    
    nit=0
    do 
       nit=nit+1

       ! Save the values of yh_av found in the previous iteration
       prev_avg_xHI = avg_xHI
       prev_temper = end_temper
	   
       ! Copy ionic abundances back to initial values (doric assumes
       ! that it contains this)
       end_xHI = begin_xHI
       end_xHII = begin_xHII
              
       ! Calculate (mean) electron density
       number_density_electron = electrondens(number_density_atom,avg_xHI,avg_xHII)

       ! Find total photo-ionization rate 

       ! (direct plus photon losses)
       ! DO THIS HERE, yh_av is changing
       ! (if the cell is ionized, add a fraction of the lost photons)
       !if (xh_intermed(pos(1),pos(2),pos(3),1) > 0.5)
       photoionization_rate_cell = photoionization_rate
       heating_rate_cell=heating_rate 
       !if (add_photon_losses) photoionization_rate_cell=photoionization_rate_cell + & 
       !     photon_loss(1)/(cell_vol*yh_av(0)*number_density_atom)
       ! GM/110225: New approach to lost photons, taking into
       ! account optical depth of cells. See notes.
       if (add_photon_losses) then
          NormFlux(0)=sum(photon_loss(:))/S_star_nominal
          ! Calculate (time averaged) column density of cell
          coldensh_cell=coldens(cell_size(1),avg_xHI,number_density_atom)
          call photoion(phi,0.0d0,coldensh_cell,cell_vol,0)
          photoionization_rate_cell=photoionization_rate_cell + phi%h/(avg_xHI*number_density_atom)
       endif

       ! Calculate the new and mean ionization states
       call doric(dt,avg_temper,number_density_electron,number_density_atom, &
	              end_xHI,end_xHII,avg_xHI,avg_xHII,photoionization_rate_cell)

       ! Reset the final temperature to the beginning temperature
       end_temper = begin_temper 

       call thermal(dt,end_temper,avg_temper,number_density_electron,number_density_atom, &
                    end_xHI,end_xHII,avg_xHI,avg_xHII,begin_xHI,begin_xHII,heating_rate_cell)
	   
       ! Test for convergence on time-averaged neutral fraction
       ! For low values of this number assume convergence
       if ((abs((avg_xHI-prev_avg_xHI)/avg_xHI) < convergence2 &
             .or. (avg_xHI < convergence_frac)) &
			 .and. &
			 abs(end_temper-prev_temper)/end_temper .lt. convergence2) then
			  exit
       endif				 
				 
				  
				  
       ! Warn about non-convergence and terminate iteration
       if (nit > 100) then
          if (rank == 0) then
              write(logf,*) 'Convergence failing (global)'
              write(logf,*) 'xh: ',avg_xHI,prev_avg_xHI
          endif
          exit
       endif
    enddo

  end subroutine solve_ionization_equation

  ! ===========================================================================

  !> Finds the column density at pos as seen from the source point srcpos
  !! through interpolation. The interpolation
  !! depends on the orientation of the ray. The ray crosses either
  !! a z-plane, a y-plane or an x-plane.
  subroutine short_characteristic (ray_traced_pos,srcpos,cdensi,path)
    
    ! Author: Garrelt Mellema
    
    ! Date: 21-Mar-2006 (06-Aug-2004)
    
    ! History:
    ! Original routine written by Alex Raga, Garrelt Mellema, Jane Arthur
    ! and Wolfgang Steffen in 1999.
    ! This version: Modified for use with a grid based approach.
    ! Better handling of the diagonals.
    ! Fortran90
    
    ! does the interpolation to find the column density at pos
    ! as seen from the source point srcpos. the interpolation
    ! depends on the orientation of the ray. The ray crosses either
    ! a z-plane, a y-plane or an x-plane.
    
    integer,dimension(Ndim),intent(in) :: ray_traced_pos !< cell position (mesh)
    integer,dimension(Ndim),intent(in) :: srcpos !< source position (mesh)
    real(kind=dp),intent(out) :: cdensi !< column density to cell
    real(kind=dp),intent(out) :: path !< path length over cell

    real(kind=dp),parameter :: sqrt3=sqrt(3.0)
    real(kind=dp),parameter :: sqrt2=sqrt(2.0)

    integer :: i,j,k,i0,j0,k0

    integer :: idel,jdel,kdel
    integer :: idela,jdela,kdela
    integer :: im,jm,km
    integer :: ip,imp,jp,jmp,kp,kmp
    integer :: sgni,sgnj,sgnk
    real(kind=dp) :: alam,xc,yc,zc,dx,dy,dz,s1,s2,s3,s4
    real(kind=dp) :: c1,c2,c3,c4
    real(kind=dp) :: dxp,dyp,dzp
    real(kind=dp) :: w1,w2,w3,w4
    real(kind=dp) :: di,dj,dk

    ! map to local variables (should be pointers ;)
    i=ray_traced_pos(1)
    j=ray_traced_pos(2)
    k=ray_traced_pos(3)
    i0=srcpos(1)
    j0=srcpos(2)
    k0=srcpos(3)
    
    ! calculate the distance between the source point (i0,j0,k0) and 
    ! the destination point (i,j,k)
    idel = i-i0
    jdel = j-j0
    kdel = k-k0
    idela = abs(idel)
    jdela = abs(jdel)
    kdela = abs(kdel)
    
    ! Find coordinates of points closer to source
    sgni = sign(1,idel)
    sgnj = sign(1,jdel)
    sgnk = sign(1,kdel)

    im = i-sgni
    jm = j-sgnj
    km = k-sgnk
    di = real(idel)
    dj = real(jdel)
    dk = real(kdel)

    ! Z plane (bottom and top face) crossing
    ! we find the central (c) point (xc,xy) where the ray crosses 
    ! the z-plane below or above the destination (d) point, find the 
    ! column density there through interpolation, and add the contribution
    ! of the neutral material between the c-point and the destination
    ! point.
    
    if (kdela >= jdela.and.kdela >= idela) then
       
       ! alam is the parameter which expresses distance along the line s to d
       ! add 0.5 to get to the interface of the d cell.
       alam = (real(km-k0)+sgnk*0.5)/dk
              
       xc = alam*di+real(i0) ! x of crossing point on z-plane 
       yc = alam*dj+real(j0) ! y of crossing point on z-plane
       
       dx = 2.0*abs(xc-(real(im)+0.5*sgni)) ! distances from c-point to
       dy = 2.0*abs(yc-(real(jm)+0.5*sgnj)) ! the corners.
       
       s1 = (1.0-dx)*(1.0-dy)    ! interpolation weights of
       s2 = (1.0-dy)*dx         ! corner points to c-point
       s3 = (1.0-dx)*dy
       s4 = dx*dy
       
       ! Map to ray_traced_pos to mesh pos, assuming a periodic mesh
       ip = modulo(i-1,mesh(1))+1
       imp = modulo(im-1,mesh(1))+1
       jp = modulo(j-1,mesh(2))+1
       jmp = modulo(jm-1,mesh(2))+1
       kmp = modulo(km-1,mesh(3))+1
       c1 = coldensh_out(imp,jmp,kmp)    !# column densities at the
       c2 = coldensh_out(ip,jmp,kmp)     !# four corners
       c3 = coldensh_out(imp,jp,kmp)
       c4 = coldensh_out(ip,jp,kmp)
       
       ! extra weights for better fit to analytical solution
       w1 = s1*weight_function(c1)
       w2 = s2*weight_function(c2)
       w3 = s3*weight_function(c3)
       w4 = s4*weight_function(c4)
       ! column density at the crossing point
       cdensi = (c1*w1+c2*w2+c3*w3+c4*w4)/(w1+w2+w3+w4) 

       if (kdela == 1.and.(idela == 1.or.jdela == 1)) then
          if (idela == 1.and.jdela == 1) then
             cdensi = sqrt3*cdensi
          else
             cdensi = sqrt2*cdensi
          endif
       endif

!if (kdela == 1.and.(idela == 1.or.jdela == 1)) then
!if (idela == 1.and.jdela == 1) then
! cdensi = sqrt3*cdensi
!else if (idela == 1 .and. jdela == 0) then
! cdensi = sqrt2*cdensi
!else if (idela == 0 .and. jdela == 1) then
! cdensi = sqrt2*cdensi
!endif
!endif

       path = sqrt((di*di+dj*dj)/(dk*dk)+1.0) ! pathlength from c to d point  


       ! y plane (left and right face) crossing
       ! (similar approach as for the z plane, see comments there)
    elseif (jdela >= idela.and.jdela >= kdela) then
          
       alam = (real(jm-j0)+sgnj*0.5)/dj
       zc = alam*dk+real(k0)
       xc = alam*di+real(i0)
       dz = 2.0*abs(zc-(real(km)+0.5*sgnk))
       dx = 2.0*abs(xc-(real(im)+0.5*sgni))
       s1 = (1.0-dx)*(1.0-dz)
       s2 = (1.0-dz)*dx
       s3 = (1.0-dx)*dz
       s4 = dx*dz
       ip = modulo(i-1,mesh(1))+1
       imp = modulo(im-1,mesh(1))+1
       jmp = modulo(jm-1,mesh(2))+1
       kp = modulo(k-1,mesh(3))+1
       kmp = modulo(km-1,mesh(3))+1
       c1 = coldensh_out(imp,jmp,kmp)
       c2 = coldensh_out(ip,jmp,kmp)
       c3 = coldensh_out(imp,jmp,kp)
       c4 = coldensh_out(ip,jmp,kp)

       ! extra weights for better fit to analytical solution
       w1 = s1*weight_function(c1)
       w2 = s2*weight_function(c2)
       w3 = s3*weight_function(c3)
       w4 = s4*weight_function(c4)
       
       cdensi = (c1*w1+c2*w2+c3*w3+c4*w4)/(w1+w2+w3+w4)
       
       ! Take care of diagonals
       if (jdela == 1.and.(idela == 1.or.kdela == 1)) then
          if (idela == 1.and.kdela == 1) then
             cdensi = sqrt3*cdensi
          else
             cdensi = sqrt2*cdensi
          endif
       endif

!if (jdela == 1.and.(idela == 1.or.kdela == 1)) then
!if (idela == 1.and.kdela == 1) then
! cdensi = sqrt3*cdensi
!else if (idela == 1 .and. kdela == 0) then
! cdensi = sqrt2*cdensi
!else if (idela == 0 .and. kdela == 1) then
! cdensi = sqrt2*cdensi
!endif
!endif

       path = sqrt((di*di+dk*dk)/(dj*dj)+1.0)
       
       ! x plane (front and back face) crossing
       ! (similar approach as with z plane, see comments there)

    elseif(idela >= jdela.and.idela >= kdela) then

       alam = (real(im-i0)+sgni*0.5)/di
       zc = alam*dk+real(k0)
       yc = alam*dj+real(j0)
       dz = 2.0*abs(zc-(real(km)+0.5*sgnk))
       dy = 2.0*abs(yc-(real(jm)+0.5*sgnj))
       s1 = (1.0-dz)*(1.0-dy)
       s2 = (1.0-dz)*dy
       s3 = (1.0-dy)*dz
       s4 = dy*dz

       imp = modulo(im-1,mesh(1))+1
       jp = modulo(j-1,mesh(2))+1
       jmp = modulo(jm-1,mesh(2))+1
       kp = modulo(k-1,mesh(3))+1
       kmp = modulo(km-1,mesh(3))+1
       c1 = coldensh_out(imp,jmp,kmp)
       c2 = coldensh_out(imp,jp,kmp)
       c3 = coldensh_out(imp,jmp,kp)
       c4 = coldensh_out(imp,jp,kp)
       ! extra weights for better fit to analytical solution
       w1 = s1*weight_function(c1)
       w2 = s2*weight_function(c2)
       w3 = s3*weight_function(c3)
       w4 = s4*weight_function(c4)
       
       cdensi = (c1*w1+c2*w2+c3*w3+c4*w4)/(w1+w2+w3+w4)
       
       if ( idela == 1 .and. ( jdela == 1 .or. kdela == 1 ) ) then
          if ( jdela == 1 .and. kdela == 1 ) then
             cdensi = sqrt3*cdensi
          else
             cdensi = sqrt2*cdensi
          endif
       endif
       
!if (idela == 1.and.(jdela == 1.or.kdela == 1)) then
!if (jdela == 1.and.kdela == 1) then
! cdensi = sqrt3*cdensi
!else if (jdela == 1 .and. kdela == 0) then
! cdensi = sqrt2*cdensi
!else if (jdela == 0 .and. kdela == 1) then
! cdensi = sqrt2*cdensi
!endif
!endif

       path = sqrt(1.0+(dj*dj+dk*dk)/(di*di))
       
    end if
    
  end subroutine short_characteristic

  ! =========================================================================

  !> Weight function for interpolation in short_characteristic
  real(kind=dp) function weight_function (cd)

    use cgsphotoconstants, only: sigh

    real(kind=dp),intent(in) :: cd

    real(kind=dp),parameter :: minweight=1.0_dp/0.6_dp

    !weight_function=1.0
    ! weight_function=1.0/max(1.0d0,cd**0.54)
    ! weight_function=exp(-min(700.0,cd*0.15*6.3d-18))
    weight_function=1.0/max(0.6_dp,cd*sigh)

    ! weight_function=1.0/log(max(e_ln,cd))

  end function weight_function

  !----------------------------------------------------------------------------

  !> Produce intermediate output for a time frame
  subroutine output_intermediate(zred_now,iteration_counter)

    ! Simple output routine. Outputs intermediate result for ionization
    ! fraction for checking convergence

    ! Output format:
    ! SM3D

    real(kind=dp),intent(in) :: zred_now
    integer,intent(in) :: iteration_counter

    integer :: i,j,k
    character(len=6) :: zred_str
    character(len=3) :: iteration_counter_string
    character(len=40) :: file1

    ! Only produce output on rank 0
    if (rank == 0) then

       write(file1,"(f6.3)") zred_now
       write(iteration_counter_string,"(i3.3)") iteration_counter
       file1=trim(adjustl(results_dir))// &
            "xfrac3d_"//trim(adjustl(file1))//"_"//iteration_counter_string//".bin"
       open(unit=52,file=file1,form="unformatted",status="unknown")
       write(52) mesh(1),mesh(2),mesh(3)
       write(52) (((intermediate_end_xHII(i,j,k),i=1,mesh(1)),j=1,mesh(2)),k=1,mesh(3))
       close(52)

    endif

  end subroutine output_intermediate

end module evolve
