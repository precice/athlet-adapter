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
   ! interfaceLocation       : position of the interface on the ATHLET geometry (optiong: "left", "right")
   character(len=50)         :: config, participantName, meshName, interfaceLocation
   ! writeInitialData        : preCICE action to write initial data
   ! readItCheckp            : preCICE action to read a coupling iteration checkpoint
   ! writeItCheckp           : preCICE action to write a coupling iteration checkpoint
   character(len=50)         :: writeInitialData, readItCheckp, writeItCheckp
   ! dimensions              : Number of dimensions for fields (currently, preCICE only supports 2D and 3D)
   ! meshID                  : preCICE meshID (returned by the preCICE solver interface "constructor")
   ! vertexID                : Vertex ids defined on the coupling mesh. Currently, only 1 vertex.
   ! dataID_Pressure         : Data ID for Pressure
   ! bool                    : A temporary auxiliary boolean variable
   ! ongoing                 : Is the coupling still ongoing? (used as boolean)
   integer                   :: dimensions=3, meshID, vertexID, dataID_Pressure, bool, ongoing
   ! interfaceIndex          : ATHLET cell index for the interface
   integer                   :: interfaceIndex
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

      ! Command-line arguments scope
      type(HashMap_t),  pointer :: arg_scope

      ! Scopes to access ATHLET data
      type(HashMap_t),  pointer :: CDTF_scope
      real(8),           dimension(:), pointer :: PRESS

      ! Bind pointer arg_scope to the athlet/global/arglist scope
      arg_scope => getScope( "athlet", "global", "arglist" )

      ! Bind pointer CDTF_scope to the state_scope/cdtf
      CDTF_scope  => getScope( state_scope, "cdtf" )
      ! Pressure from CDTF_scope
      _refVar( CDTF_scope, "press",  PRESS )

      write(*,100) "ATHLET preCICE adapter: Starting..."

      ! Get the participantName
      participantName = char(string(get(arg_scope, "-participantName")))
      ! Get the meshName
      meshName = char(string(get(arg_scope, "-meshName")))
      ! Get the interfaceIndex
      interfaceLocation = char(string(get(arg_scope, "-interfaceLocation")))

      ! Print the values
      write(*,100) "Participant name: ", participantName
      write(*,100) "Mesh name: ", meshName
      write(*,100) "Interface location: ", interfaceLocation

      ! Set the interfaceIndex according to the interfaceLocation
      ! TODO: Make interfaceIndex more flexible, currently hard-coded
      if (interfaceLocation.EQ."inlet") then
        interfaceIndex = 0
      else if (interfaceLocation.EQ."outlet") then
        interfaceIndex = 20
      endif

      !
      ! We currently hard-code the configuration file.
      ! The simulation needs to be started from the same directory as precice-config.xml.
      config = "precice-config.xml"

      ! Constants in Fortran have to be prefilled with blanks to be compatible with preCICE
      writeInitialData(1:50) ="                                                  "
      readItCheckp(1:50)     ="                                                  "
      writeItCheckp(1:50)    ="                                                  "

      write(*,100) "Creating preCICE..."

      ! Create the preCICE solver interface
      ! participantName: preCICE participant name
      ! config: preCICE configuration file name
      ! 0: MPI rank
      ! 1: MPI communicator size
      ! 50: Maximum length for the participantName string
      ! 50: Maximum length for the config string
      call precicef_create(participantName, config, 0, 1, 50, 50)

      ! Create the preCICE interface mesh
      ! We have only one vertex, we need 3 coordinates in 3D, 2 in 2D
      allocate(vertex(dimensions))
      ! Since we are currently only supporting one vertex, the location does not matter: initialize it to zero (arbitrary)
      ! TODO: Define the vertex location according to user-defined coordinates
      vertex = 0
      ! Get the preCICE mesh id for the requested mesh name (50: maximum meshName string length)
      call precicef_get_mesh_id(meshName, meshID, 50)
      ! Add vertices to the preCICE mesh (here: only one vertex)
      call precicef_set_vertex(meshID, vertex, vertexID)
      ! We don't need this vertex anymore
      ! TODO: We only used it as a temporary point, normally we would use a vertex already defined by the solver
      deallocate(vertex)

      ! Initialize preCICE
      ! This will return the maximum time step size that ATHLET is allowed to perform next
      call precicef_initialize(dtlimit)

      ! TODO: Just trying arbitrary initial values for development puropses
      if (participantName.EQ."SolverOne") then
        write(*,100) "==== PRESS: ===="
        write(*,101) PRESS
        PRESS(interfaceIndex) = 12345
      endif

      ! Get the data id for Pressure (8: number of characters in "Pressure")
      call precicef_get_data_id("Pressure", meshID, dataID_Pressure, 8)

      ! Do we need to write initial data? (specified in the config)
      call precicef_action_required(writeInitialData, bool, 50)
      if (bool.EQ.1) then
        write(*,100) "Writing initial data"
        call precicef_write_sdata(dataID_Pressure, vertexID, PRESS(interfaceIndex))
      end if
      call precicef_initialize_data()

100   format("---[adapter] ", A)
101   format("---[adapter] ", 20(f10.1,1x))
102   format("---[adapter] ", I5)

   end subroutine

   subroutine adapterexecute()
      ! Couple (if needed) at the end of every time loop iteration

      implicit none

      ! Scopes to access ATHLET data
      type(HashMap_t),  pointer :: CDTF_scope
      real(8),           dimension(:), pointer :: PRESS

      CDTF_scope  => getScope( state_scope, "cdtf" )
      _refVar( CDTF_scope, "press",  PRESS  )

      write(*,100) "ATHLET preCICE adapter: Executing..."

      ! Is the coupling still ongoing?
      call precicef_ongoing(ongoing)
      if (ongoing.NE.0) then

        ! Do we need to write a coupling iteration checkpoint? (only relevant in implicit coupling)
        call precicef_action_required(writeItCheckp, bool, 50)
        if (bool.EQ.1) then
          ! TODO: Currently, no checkpoint is being written
          write(*,100) "Writing iteration checkpoint (! NOT IMPLEMENTED !)"
          call precicef_mark_action_fulfilled(writeItCheckp, 50)
        end if

        ! Write data to preCICE buffers
        if (participantName.EQ."SolverOne") then
          write(*,100) "Pressure before writing:"
          write(*,101) PRESS
          write(*,100) "Writing Pressure"
          call precicef_write_sdata(dataID_Pressure, vertexID, PRESS(interfaceIndex))
        endif

        ! Advance the coupling
        call precicef_advance(dtlimit)

        ! Read data from preCICE buffers
        if (participantName.EQ."SolverTwo") then
          write(*,100) "Pressure before reading:"
          write(*,101) PRESS
          write(*,100) "Reading Pressure"
          call precicef_read_sdata(dataID_Pressure, 0, PRESS(interfaceIndex))
          write(*,100) "Pressure after reading:"
          write(*,101) PRESS
        endif

        ! After advancing, is the coupling still ongoing?
        call precicef_ongoing(ongoing)

        ! Do we need to go back to the last saved checkpoint?
        call precicef_action_required(readItCheckp, bool, 50)
        if (bool.EQ.1) then
          write(*,100) "Reading iteration checkpoint (! NOT IMPLEMENTED !)"
          ! TODO: Currently, no checkpoint is being read
          call precicef_mark_action_fulfilled(readItCheckp, 50)
        end if

      end if

      write(*,100) "ATHLET preCICE adapter: End of execute()"

100   format("---[adapter] ", A)
101   format("---[adapter] ", 20(f10.1,1x))

   end subroutine

   subroutine adapterfinalize()
       ! Destroy the preCICE interface and solver connections

       write(*,100) "ATHLET preCICE adapter: Finalizing..."
       call precicef_finalize()
100   format("---[adapter] ", A)
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
100   format("---[adapter] ", A)
end subroutine


!***********************************************************************
!     execute_c has to be exported in every DLL for ATHLET
!***********************************************************************
!_PROC_EXPORT(execute_c)
subroutine execute_c()
      use precice_adapter

      call adapterexecute()

      ! This hooks adapterexecute() to "After completion of one FEBE time step and print of minor time step edit and last call to AFK."
      ok = connectCallback( hook_scope, "SOPLOT_WritePlotData", adapterexecute )
      if (ok /= 1) then
         write(*,100) "Unable to set callback 'SOPLOT_WritePlotData'"
      end if

      ! This hooks adapterfinalize() to the end of ATHLET
      ok = connectCallback( hook_scope, "ATHLET_End", adapterfinalize )
      if (ok /= 1) then
         write(*,100) "Unable to set callback 'ATHLET_End'"
      end if
100   format("---[adapter] ", A)
end subroutine
