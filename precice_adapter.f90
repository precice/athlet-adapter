! Create ATHLET plugin interface module
#include "athlet/plugin.fmod"

!****************************************************************************
!     preCICE ATHLET adapter
!     Author: Gerasimos Chourdakis <chourdak@in.tum.de>
!     More information on https://github.com/precice/athlet-adapter
!****************************************************************************

module precice_adapter
  ! The adapter is an ATHLET plugin
   use athlet_plugin
   ! preCICE Fortran bindings (Fortran module based)
   use precice

   implicit none
   
   ! preCICE variables
   ! config                  : preCICE configuration file (e.g. precice-config.xml)
   ! participantName         : preCICE participant name (defined in preCICE config)
   ! meshName                : preCICE interface mesh name (defined in preCICE config)
   character(len=50)         :: config, participantName, meshName
   ! writeInitialData        : preCICE action to write initial data
   ! readItCheckp            : preCICE action to read a coupling iteration checkpoint
   ! writeItCheckp           : preCICE action to write a coupling iteration checkpoint
   character(len=50)         :: writeInitialData, readItCheckp, writeItCheckp
   ! dimensions              : Number of dimensions for fields (currently, preCICE only supports 2D and 3D)
   ! meshID                  : preCICE meshID (returned by the preCICE solver interface "constructor")
   ! vertexID                : Vertex ids defined on the coupling mesh. Currently, only 1 vertex.
   ! bool                    : A temporary auxiliary boolean variable
   ! ongoing                 : Is the coupling still ongoing? (used as boolean)
   integer                   :: dimensions=3, meshID, vertexID, bool, ongoing
   ! vertex                  : Coordinates of the interface mesh vertex
   real(8), dimension(:), allocatable :: vertex
   ! dtlimit                 : Maximum time step size allowed by preCICE
   real(8)                   :: dtlimit

   ! ATHLET plugin variables
   integer(kind=c_int)       :: c_accessorNameLength
   integer(kind=c_int)       :: c_configFileNameLength
   
contains

   subroutine adapterinitialize()
      ! Configure the adapter and initialize preCICE

      implicit none

      ! Placeholder / example: Access variables from other scopes
      type(HashMap_t),  pointer :: scope
      character(len=:), pointer :: output_id
      
      ! Get the command-line arguments scope
      type(HashMap_t),  pointer :: arg_scope
      
      scope       => getScope( state_scope, 'ccc' )
      output_id   => getCharPtr( scope, 'id' )
      write(*,*) "Output ID: ", output_id

      write (*,*) 'ATHLET preCICE adapter: Starting...'

      ! Bind pointer arg_scope to the athlet/global/arglist scope
      arg_scope => getScope( 'athlet', 'global', 'arglist' )
      ! Get the participantName
      participantName = char(string(get(arg_scope, '-participantName')))
      ! Get the meshName
      meshName = char(string(get(arg_scope, '-meshName')))
      ! Print the values
      write (*,*) 'Participant name:'
      write (*,*) participantName
      write (*,*) 'Mesh name:'
      write (*,*) meshName
      !
      ! We currently hard-code the configuration file.
      ! The simulation needs to be started from the same directory as precice-config.xml.
      config = 'precice-config.xml'

      ! Constants in Fortran have to be prefilled with blanks to be compatible with preCICE
      writeInitialData(1:50) ='                                                  '
      readItCheckp(1:50)     ='                                                  '
      writeItCheckp(1:50)    ='                                                  '

      write (*,*) 'Creating preCICE...'
      
      ! Create the preCICE solver interface
      ! participantName: preCICE participant name
      ! config: preCICE configuration file name
      ! 0: MPI rank
      ! 1: MPI communicator size
      ! 50: Maximum length for the participantName string
      ! 50: Maximum length for the config string
      call precicef_create(participantName, config, 0, 1, 50, 50)
      
      write (*,*) 'Created preCICE'
      
      ! Create the preCICE interface mesh
      ! We have only one vertex, we need 3 coordinates in 3D, 2 in 2D
      allocate(vertex(dimensions))
      ! Since this is currently a dummy adapter, initialize it to zero (arbitrary)
      vertex = 0
      ! Get the preCICE mesh id for the requested mesh name (50: maximum meshName string length)
      call precicef_get_mesh_id(meshName, meshID, 50)
      ! Add vertices to the preCICE mesh (here: only one vertex)
      call precicef_set_vertex(meshID, vertex, vertexID)
      ! We don't need this vertex anymore (we are not currently coupling anything)
      deallocate(vertex)
      
      ! Initialize preCICE
      ! This will return the maximum time step size that ATHLET is allowed to perform next
      call precicef_initialize(dtlimit)
      
      ! Do we need to write initial data? (specified in the config)
      call precicef_action_required(writeInitialData, bool, 50)
      if (bool.EQ.1) THEN
        ! This is still a dummy adapter, no actual data is being written
        write (*,*) 'DUMMY: Writing initial data'
      end if
      call precicef_initialize_data()

   end subroutine

   subroutine adapterexecute()
      ! Couple (if needed) at the end of every time loop iteration

      implicit none

      write (*,*) "ATHLET preCICE adapter: Executing..."

      ! Is the coupling still ongoing?
      call precicef_ongoing(ongoing)
      if (ongoing.NE.0) then

        ! Do we need to write a coupling iteration checkpoint? (only relevant in implicit coupling)
        call precicef_action_required(writeItCheckp, bool, 50)
        if (bool.EQ.1) THEN
          ! This is still a dummy adapter, no checkpoint is being written
          write (*,*) 'DUMMY: Writing iteration checkpoint'
          call precicef_mark_action_fulfilled(writeItCheckp, 50)
        end if

        ! Advance the coupling
        write (*,*) "ATHLET preCICE adapter: Calling advance..."
        call precicef_advance(dtlimit)
        ! After advancing, is the coupling still ongoing?
        call precicef_ongoing(ongoing)

        ! Do we need to go back to the last saved checkpoint?
        call precicef_action_required(readItCheckp, bool, 50)
        if (bool.EQ.1) THEN
          write (*,*) 'DUMMY: Reading iteration checkpoint'
          ! This is still a dummy adapter, no checkpoint is being read
          call precicef_mark_action_fulfilled(readItCheckp, 50)
        else
          write (*,*) 'DUMMY: Advancing in time'
        end if
      
      end if
      
      write (*,*) "ATHLET preCICE adapter: End of execute()"
      
   end subroutine

   subroutine adapterfinalize()
     ! Destroy the preCICE interface and solver connections
     
     write (*,*) 'ATHLET preCICE adapter: Finalizing...'
     call precicef_finalize()
   end subroutine

end module


!***********************************************************************
!     initialize_c has to be exported in every DLL for ATHLET
!***********************************************************************
!_PROC_EXPORT(initialize_c)
subroutine initialize_c()
   use precice_adapter

   call initialize_athlet_plugin()
   call adapterinitialize()
end subroutine


!***********************************************************************
!     execute_c has to be exported in every DLL for ATHLET
!***********************************************************************
!_PROC_EXPORT(execute_c)
subroutine execute_c()
  use precice_adapter

  call adapterexecute()
  
  ! This hooks adapterexecute() to "After completion of one FEBE time step and print of minor time step edit and last call to AFK."
  ok = connectCallback( hook_scope, 'ATRANS_ExtDataDone', adapterexecute )
  if (ok /= 1) then
     write(*,*) "Unable to set callback 'ATRANS_ExtDataDone'"
  end if

  ! This hooks adapterfinalize() to the end of ATHLET
  ok = connectCallback( hook_scope, 'ATHLET_End', adapterfinalize )
  if (ok /= 1) then
     write(*,*) "Unable to set callback 'ATHLET_End'"
  end if
end subroutine

