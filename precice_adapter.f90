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
   ! dataID_TL               : Data ID for Temperature of liquid
   ! bool                    : A temporary auxiliary boolean variable
   ! ongoing                 : Is the coupling still ongoing? (used as boolean)
   integer                   :: dimensions=3, meshID, vertexID, dataID_Pressure, dataID_TL, bool, ongoing
   ! interfaceIndex          : ATHLET cell index for the interface
   integer                   :: interfaceIndex
   ! vertex                  : Coordinates of the interface mesh vertex
   real(8), dimension(:), allocatable :: vertex
   ! dt_limit                 : Maximum time step size allowed by preCICE
   real(8)                   :: dt_limit
   ! dt_athlet               : Time step size computed by ATHLET
   real(8), pointer          :: dt_athlet
   ! hmax_athlet             : Maximum time step size dictated to ATHLET
   real(8), pointer          :: hmax_athlet
   ! hmax_athlet_initial     : Maximum time step size defined by ATHLET, before the adapter modifies it
   real(8)                   :: hmax_athlet_initial
   ! time_end                : End time of ATHLET
   real(8), pointer          :: time_end
   ! time                    : Current time of ATHLET
   real(8), pointer          :: time

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
      type(HashMap_t),  pointer :: CDTF_scope, CDPR_scope, cca, ccf
      real(8),           dimension(:), pointer :: PRESS, TL

      ! Bind pointer arg_scope to the athlet/global/arglist scope
      arg_scope => getScope( "athlet", "global", "arglist" )

      ! Bind pointer CDTF_scope to the state_scope/cdtf
      CDTF_scope  => getScope( state_scope, "cdtf" )
      ! Bind pointer CDTF_scope to the state_scope/cdpr
      CDPR_scope  => getScope( state_scope, "cdpr" )
      ! Pressure from CDTF_scope
      _refVar( CDTF_scope, "press",  PRESS )
      ! Temperature of liquid from CDPR_scope
      _refVar( CDPR_scope, "tl",  TL )

      write(*,100) "===== ATHLET preCICE adapter: Starting... ====="

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
      if (interfaceLocation.EQ."inlet") then
        interfaceIndex = 1
        write(*,100, advance="no") "Assuming Inlet at cell index: "
        write(*,'(I5)') interfaceIndex
      else if (interfaceLocation.EQ."outlet") then
        interfaceIndex = UBOUND(PRESS,1)
        write(*,100, advance="no") "Assuming Outlet at cell index: "
        write(*,'(I5)') interfaceIndex
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
      ! We do not need this vertex anymore
      ! TODO: We only used it as a temporary point, normally we would use a vertex already defined by the solver
      deallocate(vertex)

      ! Get the time step size dt from ATHLET
      ccf => getScope( state_scope, "ccf" )
      _refVar( ccf, "hmax", hmax_athlet )

      ! Get the max allowed time step size hmax from ATHLET
      cca => getScope( state_scope, "cca" )
      _refVar( cca, "dt", dt_athlet )


      dt_limit = hmax_athlet

      write(*,103) dt_athlet, ": Current time step size of ATHLET (pointer to dt)"
      write(*,103) dt_limit, ": Next time step size limit before initialize(dt_limit) - not a pointer, preset as hmax"
      write(*,103) hmax_athlet, ": Next max allowed time step size of ATHLET (pointer to hmax)"
      ! Initialize preCICE
      ! This will return the maximum time step size that ATHLET is allowed to perform next
      call precicef_initialize(dt_limit)
      write(*,103) dt_athlet, ": Time step size of ATHLET (pointer to dt)"
      write(*,103) dt_limit, ": Time step size limit imposed by preCICE after initialize(dt_limit)"
      write(*,103) hmax_athlet, ": Current max allowed time step size of ATHLET (pointer to hmax)"

      ! Store the initial (user-defined) maximum time step size
      hmax_athlet_initial = hmax_athlet
      write(*,103) hmax_athlet_initial, ": Initial max allowed time step size of ATHLET (backup of hmax)"
      ! Adapt the time step size for the next iteration
      hmax_athlet = min(hmax_athlet_initial, dt_limit)
      write(*,103) hmax_athlet, ": --> Next max allowed time step size of ATHLET (pointer to hmax)"

      ! TODO: Just trying arbitrary initial values for development puropses
      if (participantName.EQ."SolverOne") then
        write(*,100) "==== PRESS: ===="
        write(*,101) PRESS
        ! PRESS(interfaceIndex) = 12345 ! Set only the interface node
        PRESS = 12345                   ! Set the complete array
        write(*,100) "==== PRESS: ===="
        write(*,101) PRESS

        write(*,100) "==== TL: ===="
        write(*,101) TL
        TL(interfaceIndex) = 123      ! Set only the interface node
        TL = 123                        ! Set the complete array
        write(*,100) "==== TL: ===="
        write(*,101) TL
      endif

      ! Get the data id for Pressure (8: number of characters in "Pressure")
      call precicef_get_data_id("Pressure", meshID, dataID_Pressure, 8)

      ! Get the data id for TL (2: number of characters in "TL")
      call precicef_get_data_id("TL", meshID, dataID_TL, 2)

      ! Do we need to write initial data? (specified in the config)
      call precicef_action_required(writeInitialData, bool, 50)
      if (bool.EQ.1) then
        write(*,100) "Writing initial data"
        call precicef_write_sdata(dataID_Pressure, vertexID, PRESS(interfaceIndex))
        call precicef_write_sdata(dataID_TL, vertexID, TL(interfaceIndex))
      end if
      call precicef_initialize_data()

      write(*,100) "================================================="
      write(*,*)

100   format("---[adapter] ", A, A)
101   format("---[adapter] ", 20(f10.1,1x))
102   format("---[adapter] ", I5)
103   format("---[adapter] ", f14.7, A)

   end subroutine

   subroutine adapterexecute()
      ! Couple (if needed) at the end of every time loop iteration

      implicit none

      ! Scopes to access ATHLET data
      type(HashMap_t),  pointer :: CDTF_scope, CDPR_scope, cca, ccf, cgcsm_scope
      real(8), dimension(:), pointer :: PRESS, TL

      CDTF_scope  => getScope( state_scope, "cdtf" )
      _refVar( CDTF_scope, "press",  PRESS  )

      CDPR_scope => getScope( state_scope, "cdpr" )
      _refVar( CDPR_scope, "tl",  TL )

      ! Get the time and current/max time step size from ATHLET
      ccf => getScope( state_scope, "ccf" )
      cca => getScope( state_scope, "cca" )
      _refVar( cca, "dt", dt_athlet )
      _refVar( ccf, "hmax", hmax_athlet)
      _refVar( cca, 't',  time )

      cgcsm_scope => getScope( state_scope, "cgcsm" )
      _refVar( cgcsm_scope, "te", time_end )

      write(*,100) "===== ATHLET preCICE adapter: Executing... ====="
      write(*,102) time, ": Current time"

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
        write(*,100) "Temperature of Liquid before writing:"
        write(*,101) TL
        write(*,100) "Writing Temperature of Liquid"
        call precicef_write_sdata(dataID_TL, vertexID, TL(interfaceIndex))
      endif

      write(*,102) dt_athlet, ": Time step size of ATHLET (pointer to dt)"
      write(*,102) dt_limit, ": Previous preCICE time step size limit - before advance(dt_limit)"
      write(*,102) hmax_athlet, ": Max allowed time step size of ATHLET (pointer to hmax)"

      ! We need to give the current time step size to preCICE and it will return
      ! the maximum time step size for the next iteration in the same variable. Therefore, we first
      ! copy the athlet time step size to the variable we give to preCICE.
      dt_limit = dt_athlet
      ! Advance the coupling
      call precicef_advance(dt_limit)
      write(*,102) dt_athlet, ": Current time step size of ATHLET (pointer to dt)"
      write(*,102) dt_limit, ": Next preCICE time step size limit - after advance(dt_limit)"
      write(*,102) hmax_athlet, ": Current max allowed time step size of ATHLET (pointer to hmax)"
      write(*,102) hmax_athlet_initial, ": Initial max allowed time step size of ATHLET (backup of hmax)"

      ! Adapt the max time step size for the next iteration
      hmax_athlet = min(hmax_athlet_initial, dt_limit)
      write(*,102) hmax_athlet, ": --> Next max allowed time step size of ATHLET (pointer to hmax)"

      ! Read data from preCICE buffers
      if (participantName.EQ."SolverTwo") then
        write(*,100) "Pressure before reading:"
        write(*,101) PRESS
        write(*,100) "TL before reading:"
        write(*,101) TL
        write(*,100) "Reading Pressure"
        call precicef_read_sdata(dataID_Pressure, vertexID, PRESS(interfaceIndex))
        write(*,100) "Temperature of Liquid before reading:"
        write(*,101) TL
        write(*,100) "Reading Temperature of Liquid"
        call precicef_read_sdata(dataID_TL, vertexID, TL(interfaceIndex))
        ! TODO Debugging
        ! PRESS(interfaceIndex+1) = PRESS(interfaceIndex)
        write(*,100) "Pressure after reading:"
        write(*,101) PRESS
        write(*,100) "TL after reading:"
        write(*,101) TL
      endif

      ! After advancing, is the coupling still ongoing?
      call precicef_ongoing(ongoing)
      if (ongoing .EQ. 0) then
        write(*,100) "Coupling is not ongoing. Will now stop the simulation."
        time_end = -1
      end if

      ! Do we need to go back to the last saved checkpoint?
      call precicef_action_required(readItCheckp, bool, 50)
      if (bool.EQ.1) then
        write(*,100) "Reading iteration checkpoint (! NOT IMPLEMENTED !)"
        ! TODO: Currently, no checkpoint is being read
        call precicef_mark_action_fulfilled(readItCheckp, 50)
      end if

      write(*,100) "================================================"
      write(*,*) ""

100   format("---[adapter] ", A)
101   format("---[adapter] ", 20(f10.1,1x))
102   format("---[adapter] ", f14.7, A)

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

      ! This hooks adapterinitialize() to "After steady state calculation"
      ok = connectCallback( hook_scope, "ATHLET_SteadyStateDone", adapterinitialize )
      if (ok /= 1) then
        write(*,100) "Unable to set callback 'ATHLET_SteadyStateDone'"
      end if

100   format("---[adapter] ", A)
end subroutine


!***********************************************************************
!     execute_c has to be exported in every DLL for ATHLET
!***********************************************************************
!_PROC_EXPORT(execute_c)
subroutine execute_c()
      use precice_adapter

      ! This hooks adapterexecute() to "TAfter completion of one FEBE time step and print of minor time step edit and last call to AFK."
      ok = connectCallback( hook_scope, "ATRANS_ExtDataDone", adapterexecute )
      if (ok /= 1) then
         write(*,100) "Unable to set callback 'ATRANS_ExtDataDone'"
      end if

      ! This hooks adapterfinalize() to the end of ATHLET
      ok = connectCallback( hook_scope, "ATHLET_End", adapterfinalize )
      if (ok /= 1) then
         write(*,100) "Unable to set callback 'ATHLET_End'"
      end if
100   format("---[adapter] ", A)
end subroutine
