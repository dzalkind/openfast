!**********************************************************************************************************************************
! The HydroDyn and HydroDyn_Types modules make up a template for creating user-defined calculations in the FAST Modularization 
! Framework. HydroDyns_Types will be auto-generated based on a description of the variables for the module.
!
! "HydroDyn" should be replaced with the name of your module. Example: HydroDyn
! "HydroDyn" (in HydroDyn_*) should be replaced with the module name or an abbreviation of it. Example: HD
!..................................................................................................................................
! LICENSING
! Copyright (C) 2013-2015  National Renewable Energy Laboratory
!
!    This file is part of HydroDyn.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!    
!**********************************************************************************************************************************
MODULE HydroDyn

   USE HydroDyn_Types   
   USE NWTC_Library
   USE NWTC_Num
   USE WAMIT
   USE WAMIT2
   USE HydroDyn_Input
   USE HydroDyn_Output
   USE Current
   USE Waves2
#ifdef USE_FIT
   USE FIT_MODULES
   USE FIT_Types
#endif      
   IMPLICIT NONE
   
   PRIVATE

  
   TYPE(ProgDesc), PARAMETER            :: HydroDyn_ProgDesc = ProgDesc( 'HydroDyn', '', '' )

    
   
   
      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: HydroDyn_Init                           ! Initialization routine
   PUBLIC :: HydroDyn_End                            ! Ending routine (includes clean up)
   
   PUBLIC :: HydroDyn_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating 
                                                    !   continuous states, and updating discrete states
   PUBLIC :: HydroDyn_CalcOutput                     ! Routine for computing outputs
   
   PUBLIC :: HydroDyn_CalcConstrStateResidual        ! Tight coupling routine for returning the constraint state residual
   PUBLIC :: HydroDyn_CalcContStateDeriv             ! Tight coupling routine for computing derivatives of continuous states
   !PUBLIC :: HydroDyn_UpdateDiscState                ! Tight coupling routine for updating discrete states
   ! CKA 3/9/2018 
   Public FHA_Force   
   
   CONTAINS
   
SUBROUTINE WvStretch_Init(WaveStMod, WtrDpth, NStepWave, NNodes,  &
                          NWaveElev, WaveElev, WaveKinzi, WaveTime, &
                          WaveVel0, WaveAcc0, WaveDynP0, &
                          WavePVel0, WavePAcc0, WavePDynP0, &
                          WaveVel , WaveAcc , WaveDynP , &
                          nodeInWater, ErrStat, ErrMsg )

 
   INTEGER,          INTENT(IN   )  :: WaveStMod
   REAL(SiKi),       INTENT(IN   )  :: WtrDpth
   INTEGER,          INTENT(IN   )  :: NStepWave
   INTEGER,          INTENT(IN   )  :: NNodes              !< TODO: WHY are there both NNodes and NWaveElev ??? GJH 2/1/2016
   INTEGER,          INTENT(IN   )  :: NWaveElev
   REAL(SiKi),       INTENT(IN   )  :: WaveElev(0:,:)
   REAL(SiKi),       INTENT(IN   )  :: WaveKinzi(:)
   REAL(SiKi),       INTENT(IN   )  :: WaveTime(0:)
   REAL(SiKi),       INTENT(IN   )  :: WaveVel0(0:,:,:)               !< Wave velocity in Global coordinate system at Z = 0.  Each point in this array has a corresponding entry (same index #) in the WaveVel array
   REAL(SiKi),       INTENT(IN   )  :: WaveAcc0(0:,:,:)
   REAL(SiKi),       INTENT(IN   )  :: WaveDynP0(0:,:)
   REAL(SiKi),       INTENT(IN   )  :: WavePVel0(0:,:,:)               !< Wave velocity in Global coordinate system at Z = 0.  Each point in this array has a corresponding entry (same index #) in the WaveVel array
   REAL(SiKi),       INTENT(IN   )  :: WavePAcc0(0:,:,:)
   REAL(SiKi),       INTENT(IN   )  :: WavePDynP0(0:,:)
   REAL(SiKi),       INTENT(INOUT)  :: WaveVel(0:,:,:)
   REAL(SiKi),       INTENT(INOUT)  :: WaveAcc(0:,:,:)
   REAL(SiKi),       INTENT(INOUT)  :: WaveDynP(0:,:)
   INTEGER(IntKi),   INTENT(INOUT)  :: nodeInWater(0:,:)
   INTEGER(IntKi),   INTENT(  OUT)  :: ErrStat             !< Error status of the operation
   CHARACTER(*),     INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! Local variables
   INTEGER(IntKi) ::  I, J                            !< Local loop counters
   REAL(SiKi) :: wavekinzloc ,WavePVel0loc
   
       ! Initialize ErrStat      
   ErrStat = ErrID_None         
   ErrMsg  = ""               
      
      
   DO I = 0,NStepWave-1       ! Loop through all time steps
       
      DO J = 1,NNodes
         
         SELECT CASE ( WaveStMod )  ! Which model are we using to extrapolate the incident wave kinematics to the instantaneous free surface?

            CASE ( 0 )                 ! None = no stretching.
               ! Since we have no stretching, the wave kinematics between the seabed and
               !   the mean sea level are left unchanged; below the seabed or above the
               !   mean sea level, the wave kinematics are zero:                     
               IF (   ( WaveKinzi(J) < -WtrDpth ) .OR. ( WaveKinzi(J) > 0.0          ) )  THEN   ! .TRUE. if the elevation of the point defined by WaveKinzi(J) lies below the seabed or above mean sea level (exclusive)

                  WaveDynP   (I,J  )  = 0.0
                  WaveVel    (I,J,:)  = 0.0
                  WaveAcc    (I,J,:)  = 0.0
                  nodeInWater(I,J  )  = 0
               ELSE   
                  nodeInWater(I,J  )  = 1
               END IF
            CASE ( 1 )                 ! Vertical stretching.


               ! Vertical stretching says that the wave kinematics above the mean sea level
               !   equal the wave kinematics at the mean sea level.  The wave kinematics
               !   below the mean sea level are left unchanged:
               IF (   ( WaveKinzi(J) < -WtrDpth ) .OR. ( WaveKinzi(J) > WaveElev(I,J) ) ) THEN   ! .TRUE. if the elevation of the point defined by WaveKinzi(J) lies below the seabed or above the instantaneous wave elevation (exclusive)

                  WaveDynP   (I,J  )  = 0.0
                  WaveVel    (I,J,:)  = 0.0
                  WaveAcc    (I,J,:)  = 0.0
                  nodeInWater(I,J  )  = 0
               ELSE 
                  nodeInWater(I,J  )  = 1
                  IF   ( WaveKinzi(J) >= 0.0_ReKi ) THEN
                     ! Set the wave kinematics to the kinematics at mean sea level for locations above MSL, but below the wave elevation.
                     WaveDynP   (I,J  )  = WaveDynP0  (I,J  )
                     WaveVel    (I,J,:)  = WaveVel0   (I,J,:)
                     WaveAcc    (I,J,:)  = WaveAcc0   (I,J,:)
                  END IF
                  ! Otherwise, do nothing because the kinematics have already be set correctly via the various Waves modules
               END IF
            



            CASE ( 2 )                 ! Extrapolation stretching.


            ! Extrapolation stretching uses a linear Taylor expansion of the wave
            !   kinematics (and their partial derivatives with respect to z) at the mean
            !   sea level to find the wave kinematics above the mean sea level.  The
            !   wave kinematics below the mean sea level are left unchanged:

              
               IF (   ( WaveKinzi(J) < -WtrDpth ) .OR. ( WaveKinzi(J) > WaveElev(I,J) ) ) THEN   ! .TRUE. if the elevation of the point defined by WaveKinzi(J) lies below the seabed or above the instantaneous wave elevation (exclusive)

                  WaveDynP   (I,J  )  = 0.0
                  WaveVel    (I,J,:)  = 0.0
                  WaveAcc    (I,J,:)  = 0.0
                  nodeInWater(I,J  )  = 0
               ELSE 
                  nodeInWater(I,J  )  = 1
                  wavekinzloc = WaveKinzi(J)
                  WavePVel0loc = WavePVel0   (I,J,1)
                  IF   ( WaveKinzi(J) >= 0.0_ReKi ) THEN
                     ! Set the wave kinematics to the kinematics at mean sea level for locations above MSL, but below the wave elevation.
                     WaveDynP   (I,J  )  = WaveDynP0  (I,J  ) + WaveKinzi(J)*WavePDynP0  (I,J  )
                     WaveVel    (I,J,:)  = WaveVel0   (I,J,:) + WaveKinzi(J)*WavePVel0   (I,J,:)
                     WaveAcc    (I,J,:)  = WaveAcc0   (I,J,:) + WaveKinzi(J)*WavePAcc0   (I,J,:)
                  END IF
                  ! Otherwise, do nothing because the kinematics have already be set correctly via the various Waves modules
               END IF


            CASE ( 3 )                 ! Wheeler stretching.


            ! Wheeler stretching says that wave kinematics calculated using Airy theory
            !   at the mean sea level should actually be applied at the instantaneous
            !   free surface and that Airy wave kinematics computed at locations between
            !   the seabed and the mean sea level should be shifted vertically to new
            !   locations in proportion to their elevation above the seabed.
            !
            ! Computing the wave kinematics with Wheeler stretching requires that first
            !   say that the wave kinematics we computed at the elevations defined by
            !   the WaveKinzi0Prime(:) array are actual applied at the elevations found
            !   by stretching the elevations in the WaveKinzi0Prime(:) array using the
            !   instantaneous wave elevation--these new elevations are stored in the
            !   WaveKinzi0St(:) array.  Next, we interpolate the wave kinematics
            !   computed without stretching to the desired elevations (defined in the
            !   WaveKinzi(:) array) using the WaveKinzi0St(:) array:

 
         ENDSELECT
      END DO                   ! J - All points where the incident wave kinematics will be computed
   END DO                      ! I - All time steps
   
   ! Set the ending timestep to the same as the first timestep
   WaveDynP (NStepWave,:  )  = WaveDynP (0,:  )
   WaveVel  (NStepWave,:,:)  = WaveVel  (0,:,:)
   WaveAcc  (NStepWave,:,:)  = WaveAcc  (0,:,:)
         
END SUBROUTINE WvStretch_Init
   
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the start of the simulation to perform initialization steps. 
!! The parameters are set here and not changed during the simulation.
!! The initial states and initial guess for the input are defined.
SUBROUTINE HydroDyn_Init( InitInp, u, p, x, xd, z, OtherState, y, m, Interval, InitOut, ErrStat, ErrMsg )
!..................................................................................................................................

      TYPE(HydroDyn_InitInputType),       INTENT(IN   )  :: InitInp     !< Input data for initialization routine. TODO: This does not follow the template due to the interface of HydroDyn_CopyInitInput()
      TYPE(HydroDyn_InputType),           INTENT(  OUT)  :: u           !< An initial guess for the input; input mesh must be defined
      TYPE(HydroDyn_ParameterType),       INTENT(  OUT)  :: p           !< Parameters      
      TYPE(HydroDyn_ContinuousStateType), INTENT(  OUT)  :: x           !< Initial continuous states
      TYPE(HydroDyn_DiscreteStateType),   INTENT(  OUT)  :: xd          !< Initial discrete states
      TYPE(HydroDyn_ConstraintStateType), INTENT(  OUT)  :: z           !< Initial guess of the constraint states
      TYPE(HydroDyn_OtherStateType),      INTENT(  OUT)  :: OtherState  !< Initial other states            
      TYPE(HydroDyn_OutputType),          INTENT(  OUT)  :: y           !< Initial system outputs (outputs are not calculated; 
                                                                        !!   only the output mesh is initialized)
      TYPE(HydroDyn_MiscVarType),         INTENT(  OUT)  :: m           !< Initial misc/optimization variables           
      REAL(DbKi),                         INTENT(INOUT)  :: Interval    !< Coupling interval in seconds: the rate that 
                                                                        !!   (1) HydroDyn_UpdateStates() is called in loose coupling &
                                                                        !!   (2) HydroDyn_UpdateDiscState() is called in tight coupling.
                                                                        !!   Input is the suggested time from the glue code; 
                                                                        !!   Output is the actual coupling interval that will be used 
                                                                        !!   by the glue code.
      TYPE(HydroDyn_InitOutputType),      INTENT(  OUT)  :: InitOut     !< Output for initialization routine
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      
         ! Local variables
         
      CHARACTER(1024)                        :: SummaryName                         ! name of the HydroDyn summary file   
      TYPE(HydroDyn_InitInputType)           :: InitLocal                           ! Local version of the initialization data, needed because the framework data (InitInp) is read-only
      TYPE(Waves_InitOutputType)             :: Waves_InitOut                       ! Initialization Outputs from the Waves module initialization
!      TYPE(Waves2_InitOutputType)            :: Waves2_InitOut                      ! Initialization Outputs from the Waves2 module initialization
      TYPE(Current_InitOutputType)           :: Current_InitOut                     ! Initialization Outputs from the Current module initialization
!      LOGICAL                                :: hasWAMITOuts                        ! Are there any WAMIT-related outputs
!      LOGICAL                                :: hasMorisonOuts                      ! Are there any Morison-related outputs
!      INTEGER                                :: numHydroOuts                        ! total number of WAMIT and Morison outputs
      INTEGER                                :: I, J                                ! Generic counters
      REAL(SiKi)                             :: WaveNmbr                            ! Wavenumber of the current frequency component (1/meter)
         ! These are dummy variables to satisfy the framework, but are not used 
         
      TYPE(Waves_InputType)                  :: Waves_u                             ! Waves module initial guess for the input; the input mesh is not defined because it is not used by the waves module
      TYPE(Waves_ParameterType)              :: Waves_p                             ! Waves module parameters
      TYPE(Waves_ContinuousStateType)        :: Waves_x                             ! Waves module initial continuous states
      TYPE(Waves_DiscreteStateType)          :: Waves_xd                            ! Waves module discrete states
      TYPE(Waves_ConstraintStateType)        :: Waves_z                             ! Waves module initial guess of the constraint states
      TYPE(Waves_OtherStateType)             :: WavesOtherState                     ! Waves module other states 
      TYPE(Waves_MiscVarType)                :: Waves_m                             ! Waves module misc/optimization data 
      TYPE(Waves_OutputType)                 :: Waves_y                             ! Waves module outputs   


      TYPE(Current_InputType)                :: Current_u                           ! Current module initial guess for the input; the input mesh is not defined because it is not used by the Current module
      TYPE(Current_ParameterType)            :: Current_p                           ! Current module parameters
      TYPE(Current_ContinuousStateType)      :: Current_x                           ! Current module initial continuous states
      TYPE(Current_DiscreteStateType)        :: Current_xd                          ! Current module discrete states
      TYPE(Current_ConstraintStateType)      :: Current_z                           ! Current module initial guess of the constraint states
      TYPE(Current_OtherStateType)           :: CurrentOtherState                   ! Current module other states 
      TYPE(Current_OutputType)               :: Current_y                           ! Current module outputs   
      TYPE(Current_MiscVarType)              :: Current_m                           ! Current module misc/optimization data 
      
 
#ifdef USE_FIT
         ! FIT - related data
      TYPE(FIT_InitInputType)                :: FITInitData 
      TYPE(FIT_InputType)                    :: FIT_u                             ! FIT module initial guess for the input; the input mesh is not defined because it is not used by the waves module
      TYPE(FIT_ParameterType)                :: FIT_p                             ! FIT module parameters
      TYPE(FIT_ContinuousStateType)          :: FIT_x                             ! FIT module initial continuous states
      TYPE(FIT_DiscreteStateType)            :: FIT_xd                            ! FIT module discrete states
      TYPE(FIT_ConstraintStateType)          :: FIT_z                             ! FIT module initial guess of the constraint states
      TYPE(FIT_OtherStateType)               :: FIT_OtherState                    ! FIT module other/optimization states 
      TYPE(FIT_OutputType)                   :: FIT_y                             ! FIT module outputs  
      TYPE(FIT_InitOutputType)               :: FIT_InitOut                       ! Initialization Outputs from the FIT module initialization
#endif

      Real(ReKi)                             :: Np      
      Real(ReKi)                             :: dftreal
      Real(ReKi)                             :: dftimag 
   
         ! Wave Stretching Data
      REAL(SiKi), ALLOCATABLE  :: tmpWaveKinzi(:    )
      INTEGER                  :: tmpNWaveElev
      REAL(SiKi), ALLOCATABLE  :: tmpWaveElevxi(:    )
      REAL(SiKi), ALLOCATABLE  :: tmpWaveElevyi(:    )
      REAL(SiKi), ALLOCATABLE  :: WaveElevSt  (:,:  ) 
      REAL(SiKi), ALLOCATABLE  :: WaveVel0    (:,:,:) 
      REAL(SiKi), ALLOCATABLE  :: WaveAcc0    (:,:,:)                              
      REAL(SiKi), ALLOCATABLE  :: WaveDynP0   (:,:  )  
      REAL(SiKi), ALLOCATABLE  :: WaveVel2S0  (:,:,:)
      REAL(SiKi), ALLOCATABLE  :: WaveAcc2S0  (:,:,:)                                   
      REAL(SiKi), ALLOCATABLE  :: WaveDynP2S0 (:,:  )   
      REAL(SiKi), ALLOCATABLE  :: WaveVel2D0  (:,:,:)    
      REAL(SiKi), ALLOCATABLE  :: WaveAcc2D0  (:,:,:)                              
      REAL(SiKi), ALLOCATABLE  :: WaveDynP2D0 (:,:  )                                     
                                       
                                                  
      INTEGER(IntKi)                         :: ErrStat2                            ! local error status
      CHARACTER(1024)                        :: ErrMsg2                             ! local error message
      CHARACTER(*), PARAMETER                :: RoutineName = 'HydroDyn_Init'
   

      
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      p%UnOutFile = -1 !bjj: this was being written to the screen when I had an error in my HD input file, so I'm going to initialize here.
      
#ifdef BETA_BUILD
   CALL DispBetaNotice( "This is a beta version of HydroDyn and is for testing purposes only."//NewLine//"This version includes user waves, WaveMod=6 and the ability to write example user waves." )
#endif
         ! Copy the initialization input data to a local version because the framework states InitInp should have INTENT (IN), but due to an issue the the
         ! copy routine, we needed to make it (INOUT), which means we actually don't need this local version!!  I'm leaving this with the idea that the
         ! copy routine will get modified so that InitInp can have INTENT (IN) again.  GJH 4-Apr-2013
         
      CALL HydroDyn_CopyInitInput( InitInp, InitLocal, MESH_NEWCOPY, ErrStat2, ErrMsg2 )   
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF
         
         ! CHANGE: CKA 3/16/2018 Add TMax to misc. data to be used in FHA Force
         m%TMax=InitInp%TMax
         ! END CHANGE
      
         ! Initialize the NWTC Subroutine Library
         
      CALL NWTC_Init(  )
     
        
         ! Display the module information

      CALL DispNVD( HydroDyn_ProgDesc )        
      
      
      
      
      IF ( InitInp%UseInputFile ) THEN
         
                  
         ! Parse all HydroDyn-related input files and populate the *_InitInputType derived types
         
         CALL HydroDynInput_GetInput( InitLocal, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) THEN
               CALL CleanUp()
               RETURN
            END IF
         
      END IF
           
      
         ! Start with the glue code's timestep.  This may be altered in the Input file processing, and we will check that afterwards.
                 
      InitLocal%DT  = Interval
      
      ! DZ Change: Copy Hull TMD filename to parameters
      m%TMDFile         = InitLocal%TMDFile
      m%TMDControlFile  = InitLocal%TMDControlFile
      m%OutRootName     = InitLocal%OutRootName
         ! Verify all the necessary initialization data. Do this at the HydroDynInput module-level 
         !   because the HydroDynInput module is also responsible for parsing all this 
         !   initialization data from a file


      CALL HydroDynInput_ProcessInitData( InitLocal, ErrStat2, ErrMsg2 )     
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF
      
      
        ! Since the Convolution Radiation module is currently the only module which requires knowledge of the time step size, 
        !  we will set Hydrodyn's time step to be that of the Convolution radiation module if it is being used.  Otherwise, we
        !  will set it to be equal to the glue-codes
      IF ((Initlocal%PotMod == 1) .AND. (Initlocal%WAMIT%RdtnMod == 1) ) THEN
         IF ( .NOT. EqualRealNos(Interval,InitLocal%WAMIT%Conv_Rdtn%RdtnDT) ) THEN
            CALL SetErrStat(ErrID_Fatal,'The value of RdtnDT is not equal to the glue code timestep.  This is not allowed in the current version of HydroDyn.',ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) THEN
               CALL CleanUp()
               RETURN
            END IF   
         END IF
         
         p%DT = InitLocal%WAMIT%Conv_Rdtn%RdtnDT
 
#ifdef USE_FIT
      ELSE IF (Initlocal%PotMod == 2) THEN
         ! This is the FIT potential flow model and the time step needs to be >= the driver timestep, and and integer multiple if larger
         ! We example WaveDT for this timestep size because FIT is tied to WaveDT
         IF ( ( .NOT. EqualRealNos(mod(real(Initlocal%Waves%WaveDT,ReKi), real(Interval,ReKi)) , 0.0_ReKi) ) .OR. Initlocal%Waves%WaveDT <= 0.0_DbKi ) THEn
            CALL SetErrStat(ErrID_Fatal,'The value of WaveDT is not greater than zero and an integer multiple of the glue code timestep.',ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) THEN
               CALL CleanUp()
               RETURN
            END IF 
         ELSE
            p%DT = Interval  !  If the above check is ok, then we can still set the module DT to the glue-code dt
            
         END IF
#endif

      ELSE
         
         p%DT = Interval
      END IF  
      
         ! Open a summary of the HydroDyn Initialization. Note: OutRootName must be set by the caller because there may not be an input file to obtain this rootname from.
         
      IF ( InitLocal%HDSum ) THEN 
         
         SummaryName = TRIM(InitLocal%OutRootName)//'.HD.sum'
         CALL HDOut_OpenSum( InitLocal%UnSum, SummaryName, HydroDyn_ProgDesc, ErrStat2, ErrMsg2 )    !this must be called before the Waves_Init() routine so that the appropriate wave data can be written to the summary file
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) THEN
               CALL CleanUp()
               RETURN
            END IF
      
      ELSE
         
         InitLocal%UnSum = -1
         
      END IF
      
         ! Copy Additional preload, stiffness, and damping to the parameters
      p%AddF0        = InitLocal%AddF0
      p%AddCLin      = InitLocal%AddCLin
      p%AddBLin      = InitLocal%AddBLin
      p%AddBQuad     = InitLocal%AddBQuad
      
      
         ! Set summary unit number in Waves, Radiation, and Morison initialization input data
         
      InitLocal%Waves%UnSum           = InitLocal%UnSum
      InitLocal%WAMIT%Conv_Rdtn%UnSum = InitLocal%UnSum
      InitLocal%Morison%UnSum         = InitLocal%UnSum      
    
      
         ! Now call each sub-module's *_Init subroutine
         ! to fully initialize each sub-module based on the necessary initialization data
      

         ! Initialize Current module
         
      CALL Current_Init(InitLocal%Current, Current_u, Current_p, Current_x, Current_xd, Current_z, CurrentOtherState, &
                                 Current_y, Current_m, Interval, Current_InitOut, ErrStat2, ErrMsg2 )   
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF

      ! Verify that Current_Init() did not request a different Interval!
      
      IF ( p%DT /= Interval ) THEN
         CALL SetErrStat(ErrID_Fatal,'Current Module attempted to change timestep interval, but this is not allowed.  Current Module must use the HydroDyn Interval.',ErrStat,ErrMsg,RoutineName)
         CALL CleanUp()
         RETURN
      END IF
      
      
         ! Move initialization output data from Current module into the initialization input data for the Waves module
                    
      IF (ALLOCATED(Current_InitOut%CurrVxi)) CALL Move_Alloc( Current_InitOut%CurrVxi, InitLocal%Waves%CurrVxi )
      IF (ALLOCATED(Current_InitOut%CurrVyi)) CALL Move_Alloc( Current_InitOut%CurrVyi, InitLocal%Waves%CurrVyi )
      
      InitLocal%Waves%PCurrVxiPz0   = Current_InitOut%PCurrVxiPz0
      InitLocal%Waves%PCurrVyiPz0   = Current_InitOut%PCurrVyiPz0
         

         ! Copy the WaveElevXY data in from the HydroDyn InitLocal (already a copy of InitInp)

      IF (ALLOCATED(InitLocal%WaveElevXY)) CALL MOVE_ALLOC(InitLocal%WaveElevXY, InitLocal%Waves%WaveElevXY)  
 
         ! Initialize Waves module
      
!==========================================================================
! Initialize Wave Stretching data for 1st Order Waves
!==========================================================================
      IF (InitLocal%Waves%WaveStMod > 0) THEN      
            ! Allocate the temporary storage array for the WvKinxi
         ALLOCATE ( tmpWaveKinzi(InitLocal%Waves%NWaveKin), STAT = ErrStat2 )
         IF ( ErrStat2 /= 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'Error allocating space for tmpWaveKinzi array.', ErrStat, ErrMsg, RoutineName)
            CALL CleanUp()
            RETURN
         END IF
            
            
         
         tmpWaveKinzi = InitLocal%Waves%WaveKinzi
         InitLocal%Waves%WaveKinzi = 0.0_ReKi         ! Force all zi coordinates to 0.0 for this version of the Waves initialization
         
         
            ! We will use the user-requested wave elevation arrays to compute the wave elevations for stretching at ALL node locations.
            ! We are going to store the user-requested wave elevation output locations so that we can restore them after we done.
         IF (InitLocal%Waves%NWaveElev > 0) THEN
            tmpNWaveElev = InitLocal%Waves%NWaveElev
            CALL MOVE_ALLOC( InitLocal%Waves%WaveElevxi, tmpWaveElevxi  )  ! (from, to)
            CALL MOVE_ALLOC( InitLocal%Waves%WaveElevyi, tmpWaveElevyi  ) 
         END IF
           
           
         ALLOCATE ( InitLocal%Waves%WaveElevxi(InitLocal%Waves%NWaveKin), STAT = ErrStat2 )
         IF ( ErrStat2 /= 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'Error allocating space for tmpWaveKinzi array.', ErrStat, ErrMsg, RoutineName)
            CALL CleanUp()
            RETURN
         END IF
         ALLOCATE ( InitLocal%Waves%WaveElevyi(InitLocal%Waves%NWaveKin), STAT = ErrStat2 )
         IF ( ErrStat2 /= 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'Error allocating space for tmpWaveKinzi array.', ErrStat, ErrMsg, RoutineName)
            CALL CleanUp()
            RETURN
         END IF    
         
         InitLocal%Waves%NWaveElev  = InitLocal%Waves%NWaveKin
         InitLocal%Waves%WaveElevxi = InitLocal%Waves%WaveKinxi
         InitLocal%Waves%WaveElevyi = InitLocal%Waves%WaveKinyi
         
         
         CALL Waves_Init(InitLocal%Waves, Waves_u, Waves_p, Waves_x, Waves_xd, Waves_z, WavesOtherState, &
                                    Waves_y, Waves_m, Interval, Waves_InitOut, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF
         
            ! Store the wave elevations coming out of the Waves_Init for use in the stretching calculations
         ALLOCATE ( WaveElevSt(0:Waves_InitOut%NStepWave,InitLocal%Waves%NWaveKin), STAT = ErrStat2 )
         IF ( ErrStat2 /= 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'Error allocating space for WaveElevSt array.', ErrStat, ErrMsg, RoutineName)
            CALL CleanUp()
            RETURN
         END IF    
         WaveElevSt = Waves_InitOut%WaveElev
         
         
            ! We need to reset the wave elevation arrays
         DEALLOCATE(InitLocal%Waves%WaveElevxi)
         DEALLOCATE(InitLocal%Waves%WaveElevyi)
         InitLocal%Waves%NWaveElev = tmpNWaveElev
         
         IF (InitLocal%Waves%NWaveElev > 0) THEN
            CALL MOVE_ALLOC( tmpWaveElevxi, InitLocal%Waves%WaveElevxi  )  ! (from, to)
            CALL MOVE_ALLOC( tmpWaveElevyi, InitLocal%Waves%WaveElevyi  ) 
         END IF
         
         ALLOCATE ( WaveDynP0 (0:Waves_InitOut%NStepWave,InitLocal%Waves%NWaveKin  ), STAT=ErrStat2 )
         IF (ErrStat2 /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveDynP0.', ErrStat, ErrMsg, RoutineName)

         ALLOCATE ( WaveVel0  (0:Waves_InitOut%NStepWave,InitLocal%Waves%NWaveKin,3), STAT=ErrStat2 )
         IF (ErrStat2 /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveVel0.',  ErrStat, ErrMsg, RoutineName)

         ALLOCATE ( WaveAcc0  (0:Waves_InitOut%NStepWave,InitLocal%Waves%NWaveKin,3), STAT=ErrStat2 )
         IF (ErrStat2 /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAcc0.',  ErrStat, ErrMsg, RoutineName)
              
         
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN         
         END IF
         
               ! Copy the init output arrays into the MSL versions
         WaveDynP0  =      Waves_InitOut%WaveDynP     
         WaveAcc0   =      Waves_InitOut%WaveAcc  
         WaveVel0   =      Waves_InitOut%WaveVel
         
         
         InitLocal%Waves%WaveKinzi =  tmpWaveKinzi
         
            ! Deallocate data which will be allocated again within the Waves_Init routine
         DEALLOCATE( Waves_InitOut%WaveDynP )
         DEALLOCATE( Waves_InitOut%WaveAcc )
         DEALLOCATE( Waves_InitOut%WaveVel )
         DEALLOCATE( Waves_InitOut%PWaveDynP0 )
         DEALLOCATE( Waves_InitOut%PWaveAcc0 )
         DEALLOCATE( Waves_InitOut%PWaveVel0 )
         DEALLOCATE( Waves_InitOut%WaveElevC0)   
         DEALLOCATE( Waves_InitOut%WaveDirArr)   
         DEALLOCATE( Waves_InitOut%WaveElev  )
         DEALLOCATE( Waves_InitOut%WaveTime  )
         DEALLOCATE( Waves_InitOut%NodeInWater  )
      END IF       
!==========================================================================     
          
      CALL Waves_Init(InitLocal%Waves, Waves_u, Waves_p, Waves_x, Waves_xd, Waves_z, WavesOtherState, &
                                 Waves_y, Waves_m, Interval, Waves_InitOut, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      END IF
      
      
      ! Verify that Waves_Init() did not request a different Interval!
      
      IF ( p%DT /= Interval ) THEN
         CALL SetErrStat(ErrID_Fatal,'Waves Module attempted to change timestep interval, but this is not allowed.  Waves Module must use the HydroDyn Interval.',ErrStat,ErrMsg,RoutineName)
         CALL CleanUp()
         RETURN
      END IF
     
         ! Copy the wave elevation time series corresponding to WaveElevXY to the output.

      IF (ALLOCATED(Waves_InitOut%WaveElevSeries)) CALL MOVE_ALLOC( Waves_InitOut%WaveElevSeries, InitOut%WaveElevSeries )
      IF (ALLOCATED(InitLocal%Waves%WaveElevXY)) CALL MOVE_ALLOC(InitLocal%Waves%WaveElevXY, InitLocal%WaveElevXY) ! move this back for waves2 later 

      
         ! Copy Waves initialization output into the initialization input type for the WAMIT module
      p%NWaveElev    = InitLocal%Waves%NWaveElev  
      p%NStepWave    = Waves_InitOut%NStepWave
      
      CALL MOVE_ALLOC( Waves_InitOut%WaveTime, p%WaveTime  ) 
      CALL MOVE_ALLOC( Waves_InitOut%WaveElev, p%WaveElev1 ) ! allocate p%WaveElev1, set p%WaveElev1 = Waves_InitOut%WaveElev, and deallocate Waves_InitOut%WaveElev
      
         ! Copy the first order wave elevation information to p%WaveElev1 so that we can output the total, first, and second order wave elevation separately
      ALLOCATE ( p%WaveElev   (0:p%NStepWave, p%NWaveElev ) , STAT=ErrStat2 )
      IF ( ErrStat2 /= 0 )  THEN
         CALL SetErrStat(ErrID_Fatal,'Error allocating memory for the WaveElev array.',ErrStat,ErrMsg,RoutineName)
         CALL CleanUp()
         RETURN         
      END IF
      p%WaveElev = p%WaveElev1



      m%LastIndWave = 1

      
      IF ( InitLocal%Waves%WaveMod /= 6 ) THEN
   
            !----------------------------------
            ! Initialize Waves2 module
            !----------------------------------
   
   
         IF (InitLocal%Waves2%WvDiffQTFF .OR. InitLocal%Waves2%WvSumQTFF ) THEN
               ! Set a few things from the Waves module output
            InitLocal%Waves2%NStepWave   = Waves_InitOut%NStepWave
            InitLocal%Waves2%NStepWave2  = Waves_InitOut%NStepWave2
            InitLocal%Waves2%WaveDOmega  = Waves_InitOut%WaveDOmega
                                                
               ! Copy the WaveElevXY data in from the HydroDyn InitLocal, already a copy from InitInp
            IF (ALLOCATED(InitLocal%WaveElevXY)) CALL MOVE_ALLOC(InitLocal%WaveElevXY, InitLocal%Waves2%WaveElevXY) 
   
               ! Temporarily move arrays to init input for Waves2 (save some space)
            CALL MOVE_ALLOC(p%WaveTime, InitLocal%Waves2%WaveTime) 
            CALL MOVE_ALLOC(Waves_InitOut%WaveElevC0, InitLocal%Waves2%WaveElevC0)
            CALL MOVE_ALLOC(Waves_InitOut%WaveDirArr, InitLocal%Waves2%WaveDirArr)
   
   !bjj: note that this doesn't get called if .not. (InitLocal%Waves2%WvDiffQTFF .OR. InitLocal%Waves2%WvSumQTFF), so p%waves2%* never get set
   ! however, they get queried later in the code!!!! I've set these parameters in an "else" statement, below

!==========================================================================
! Initialize Wave Stretching data for 2nd Order Waves
!==========================================================================
            IF (InitLocal%Waves%WaveStMod > 0) THEN      
                  ! Set the wave kinematics zi locations to zero to generate kinematics at MSL
               InitLocal%Waves2%WaveKinzi = 0
         
                  ! We will use the user-requested wave elevation arrays to compute the wave elevations for stretching at ALL node locations.
                  ! We are going to store the user-requested wave elevation output locations so that we can restore them after we done.
               IF (InitLocal%Waves2%NWaveElev > 0) THEN
                  tmpNWaveElev = InitLocal%Waves2%NWaveElev
                  CALL MOVE_ALLOC( InitLocal%Waves2%WaveElevxi, tmpWaveElevxi  )  ! (from, to)
                  CALL MOVE_ALLOC( InitLocal%Waves2%WaveElevyi, tmpWaveElevyi  ) 
               END IF
           
           
               ALLOCATE ( InitLocal%Waves2%WaveElevxi(InitLocal%Waves2%NWaveKin), STAT = ErrStat2 )
               IF ( ErrStat2 /= 0 ) THEN
                  CALL SetErrStat( ErrID_Fatal,'Error allocating space for WaveElevxi array.', ErrStat, ErrMsg, RoutineName)
                  CALL CleanUp()
                  RETURN
               END IF
               ALLOCATE ( InitLocal%Waves2%WaveElevyi(InitLocal%Waves2%NWaveKin), STAT = ErrStat2 )
               IF ( ErrStat2 /= 0 ) THEN
                  CALL SetErrStat( ErrID_Fatal,'Error allocating space for WaveElevyi array.', ErrStat, ErrMsg, RoutineName)
                  CALL CleanUp()
                  RETURN
               END IF    
         
               InitLocal%Waves2%NWaveElev  = InitLocal%Waves2%NWaveKin
               InitLocal%Waves2%WaveElevxi = InitLocal%Waves2%WaveKinxi
               InitLocal%Waves2%WaveElevyi = InitLocal%Waves2%WaveKinyi                        
                  
               CALL Waves2_Init(InitLocal%Waves2, m%u_Waves2, p%Waves2, x%Waves2, xd%Waves2, z%Waves2, OtherState%Waves2, &
                                          y%Waves2, m%Waves2, Interval, InitOut%Waves2, ErrStat2, ErrMsg2 )
                  CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
                  IF ( ErrStat >= AbortErrLev ) THEN
                     CALL CleanUp()
                     RETURN
                  END IF
            
  
                  ! Store the wave elevations coming out of the Waves_Init for use in the stretching calculations      
               WaveElevSt = WaveElevSt + p%Waves2%WaveElev2
         
                  ! We need to reset the wave elevation arrays
               DEALLOCATE(InitLocal%Waves2%WaveElevxi)
               DEALLOCATE(InitLocal%Waves2%WaveElevyi)
               InitLocal%Waves2%NWaveElev = tmpNWaveElev
         
               IF (InitLocal%Waves2%NWaveElev > 0) THEN
                  CALL MOVE_ALLOC( tmpWaveElevxi, InitLocal%Waves2%WaveElevxi  )  ! (from, to)
                  CALL MOVE_ALLOC( tmpWaveElevyi, InitLocal%Waves2%WaveElevyi  ) 
               END IF
                  
                  
               ALLOCATE ( WaveDynP2D0 (0:Waves_InitOut%NStepWave,InitLocal%Waves%NWaveKin  ), STAT=ErrStat2 )
               IF (ErrStat2 /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveDynP2D0.', ErrStat, ErrMsg, RoutineName)

               ALLOCATE ( WaveVel2D0  (0:Waves_InitOut%NStepWave,InitLocal%Waves%NWaveKin,3), STAT=ErrStat2 )
               IF (ErrStat2 /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveVel2D0.',  ErrStat, ErrMsg, RoutineName)

               ALLOCATE ( WaveAcc2D0  (0:Waves_InitOut%NStepWave,InitLocal%Waves%NWaveKin,3), STAT=ErrStat2 )
               IF (ErrStat2 /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAcc2D0.',  ErrStat, ErrMsg, RoutineName)
         
               ALLOCATE ( WaveDynP2S0 (0:Waves_InitOut%NStepWave,InitLocal%Waves%NWaveKin  ), STAT=ErrStat2 )
               IF (ErrStat2 /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveDynP2S0.', ErrStat, ErrMsg, RoutineName)

               ALLOCATE ( WaveVel2S0  (0:Waves_InitOut%NStepWave,InitLocal%Waves%NWaveKin,3), STAT=ErrStat2 )
               IF (ErrStat2 /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveVel2S0.',  ErrStat, ErrMsg, RoutineName)

               ALLOCATE ( WaveAcc2S0  (0:Waves_InitOut%NStepWave,InitLocal%Waves%NWaveKin,3), STAT=ErrStat2 )
               IF (ErrStat2 /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAcc2S0.',  ErrStat, ErrMsg, RoutineName)      
         
               IF ( ErrStat >= AbortErrLev ) THEN
                  CALL CleanUp()
                  RETURN         
               END IF

                     ! Copy the init output arrays into the MSL versions
               WaveDynP2D0  =      InitOut%Waves2%WaveDynP2D     
               WaveAcc2D0   =      InitOut%Waves2%WaveAcc2D  
               WaveVel2D0   =      InitOut%Waves2%WaveVel2D
               WaveDynP2S0  =      InitOut%Waves2%WaveDynP2S     
               WaveAcc2S0   =      InitOut%Waves2%WaveAcc2S  
               WaveVel2S0   =      InitOut%Waves2%WaveVel2S
         
                  ! Reset the wave kinematics zi locations 
               InitLocal%Waves2%WaveKinzi = InitLocal%Waves%WaveKinzi
         
                  ! Deallocate arrays which will be re-allocated in the next call to Waves2_Init
               DEALLOCATE ( p%Waves2%WaveElev2        )
               DEALLOCATE ( InitOut%Waves2%WaveVel2D  )
               DEALLOCATE ( InitOut%Waves2%WaveAcc2D  )
               DEALLOCATE ( InitOut%Waves2%WaveDynP2D )
               DEALLOCATE ( InitOut%Waves2%WaveVel2S  )
               DEALLOCATE ( InitOut%Waves2%WaveAcc2S  )
               DEALLOCATE ( InitOut%Waves2%WaveDynP2S )
               
            END IF       
!==========================================================================     
            
                               
            
            
            
            
            CALL Waves2_Init(InitLocal%Waves2, m%u_Waves2, p%Waves2, x%Waves2, xd%Waves2, z%Waves2, OtherState%Waves2, &
                                    y%Waves2, m%Waves2, Interval, InitOut%Waves2, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) THEN
               CALL CleanUp()
               RETURN
            END IF
   
               ! move arrays back
            CALL MOVE_ALLOC(InitLocal%Waves2%WaveTime, p%WaveTime) 
            CALL MOVE_ALLOC(InitLocal%Waves2%WaveElevC0, Waves_InitOut%WaveElevC0)
            CALL MOVE_ALLOC(InitLocal%Waves2%WaveDirArr, Waves_InitOut%WaveDirArr)
                  
            ! Verify that Waves2_Init() did not request a different Interval!
   
            IF ( p%DT /= Interval ) THEN
               CALL SetErrStat(ErrID_Fatal,'Waves2 Module attempted to change timestep interval, but this is not allowed. '// &
                                          ' Waves2 Module must use the HydroDyn Interval.',ErrStat,ErrMsg,RoutineName)
               CALL CleanUp()
               RETURN
            END IF
   
   
            ! If we calculated the wave elevation series data (for visualization purposes), add the second order corrections to the first order.
            IF (ALLOCATED(InitLocal%Waves2%WaveElevXY)) THEN
                  ! Make sure the sizes of the two resulting arrays are identical...
               IF ( SIZE(InitOut%WaveElevSeries,DIM=1) /= SIZE(InitOut%Waves2%WaveElevSeries2,DIM=1) .OR. &
                    SIZE(InitOut%WaveElevSeries,DIM=2) /= SIZE(InitOut%Waves2%WaveElevSeries2,DIM=2)) THEN
                  CALL SetErrStat(ErrID_Fatal,' WaveElevSeries arrays for first and second order wave elevations are of different sizes.',ErrStat,ErrMsg,RoutineName)
                  CALL CleanUp()
                  RETURN
               ELSE
                  DO J=1,SIZE(InitOut%WaveElevSeries,DIM=2)
                     DO I = 0,p%NStepWave
                        InitOut%WaveElevSeries(I,J)  =  InitOut%Waves2%WaveElevSeries2(I,J) + InitOut%WaveElevSeries(I,J)
                     ENDDO
                  ENDDO
               ENDIF
            ENDIF
   
            ! If we calculated wave elevations, it is now stored in p%WaveElev.  So we need to add the corrections.
            IF (p%Waves2%NWaveElev > 0 ) THEN
                  ! Make sure the sizes of the two resulting arrays are identical...
               IF ( SIZE(p%WaveElev,DIM=1) /= SIZE(p%Waves2%WaveElev2,DIM=1) .OR. &
                    SIZE(p%WaveElev,DIM=2) /= SIZE(p%Waves2%WaveElev2,DIM=2)) THEN
                  CALL SetErrStat(ErrID_Fatal,' WaveElev(NWaveElev) arrays for first and second order wave elevations are of different sizes.',ErrStat,ErrMsg,RoutineName)
                  CALL CleanUp()
                  RETURN
               ELSE
                  DO J=1,SIZE(p%Waves2%WaveElev2,DIM=2)
                     DO I = 0,p%NStepWave
                        p%WaveElev(I,J)  =  p%Waves2%WaveElev2(I,J) + p%WaveElev(I,J)
                     ENDDO
                  ENDDO
               ENDIF
            ENDIF
   
            ! The acceleration, velocity, and dynamic pressures will get added to the parts passed to the morrison module later...
   
         ELSE
                  ! these need to be set to zero since we don't have a UseWaves2 flag:
               p%Waves2%NWaveElev  = 0
               p%Waves2%WvDiffQTFF = .FALSE.
               p%Waves2%WvSumQTFF  = .FALSE.
               p%Waves2%NumOuts    = 0
               
         ENDIF
   
   
   
   
            ! Is there a WAMIT body? 
         
         IF ( InitLocal%PotMod == 1 ) THEN
            
               ! Copy Waves initialization output into the initialization input type for the WAMIT module
                  
            InitLocal%WAMIT%RhoXg        = Waves_InitOut%RhoXg
            InitLocal%WAMIT%NStepWave    = Waves_InitOut%NStepWave
            InitLocal%WAMIT%NStepWave2   = Waves_InitOut%NStepWave2
            InitLocal%WAMIT%WaveDirMin   = Waves_InitOut%WaveDirMin
            InitLocal%WAMIT%WaveDirMax   = Waves_InitOut%WaveDirMax
            InitLocal%WAMIT%WaveDOmega   = Waves_InitOut%WaveDOmega   
                        
               ! Temporarily move arrays to init input for WAMIT (save some space)
            CALL MOVE_ALLOC(p%WaveTime,               InitLocal%WAMIT%WaveTime) 
            CALL MOVE_ALLOC(Waves_InitOut%WaveElevC0, InitLocal%WAMIT%WaveElevC0) 
            CALL MOVE_ALLOC(Waves_InitOut%WaveDirArr, InitLocal%WAMIT%WaveDirArr) 
               
               !-----------------------------------------
               ! Initialize the WAMIT Calculations 
               !-----------------------------------------
              
            CALL WAMIT_Init(InitLocal%WAMIT, m%u_WAMIT, p%WAMIT, x%WAMIT, xd%WAMIT, z%WAMIT, OtherState%WAMIT, &
                                    y%WAMIT, m%WAMIT, Interval, InitOut%WAMIT, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) THEN
               CALL CleanUp()
               RETURN
            END IF
            
            ! Generate Summary file information for WAMIT module
                ! Compute the load contribution from hydrostatics:
            IF ( InitLocal%UnSum > 0 ) THEN
               
               WRITE( InitLocal%UnSum, '(A11)')          'WAMIT Model'
               WRITE( InitLocal%UnSum, '(A11)')          '-----------'
               WRITE( InitLocal%UnSum, '(A42,2X,ES15.6)') 'Displaced volume (m^3)                 :', p%WAMIT%PtfmVol0
               WRITE( InitLocal%UnSum, '(A42,2X,ES15.6)') 'X-offset of the center of buoyancy (m) :', p%WAMIT%PtfmCOBxt
               WRITE( InitLocal%UnSum, '(A42,2X,ES15.6)') 'Y-offset of the center of buoyancy (m) :', p%WAMIT%PtfmCOByt
               WRITE( InitLocal%UnSum,  '(/)' ) 
               WRITE( InitLocal%UnSum, '(A81)' ) 'Buoyancy loads from members modelled with WAMIT, summed about ( 0.0, 0.0, 0.0 )'
               WRITE( InitLocal%UnSum, '(18x,6(2X,A20))' ) ' BuoyFxi ', ' BuoyFyi ', ' BuoyFzi ', ' BuoyMxi ', ' BuoyMyi ', ' BuoyMzi '
               WRITE( InitLocal%UnSum, '(18x,6(2X,A20))' ) '   (N)   ', '   (N)   ', '   (N)   ', '  (N-m)  ', '  (N-m)  ', '  (N-m)  '
               WRITE( InitLocal%UnSum, '(A18,6(2X,ES20.6))') '  External:       ',0.0,0.0,p%WAMIT%RhoXg*p%WAMIT%PtfmVol0,p%WAMIT%RhoXg*p%WAMIT%PtfmVol0*p%WAMIT%PtfmCOByt, -p%WAMIT%RhoXg*p%WAMIT%PtfmVol0*p%WAMIT%PtfmCOBxt, 0.0   ! and the moment about Y due to the COB being offset from the WAMIT reference point
            
            END IF
            
            
               ! Verify that WAMIT_Init() did not request a different Interval!
         
            IF ( p%DT /= Interval ) THEN
               CALL SetErrStat(ErrID_Fatal,'WAMIT Module attempted to change timestep interval, but this is not allowed.  WAMIT Module must use the HydroDyn Interval.',ErrStat,ErrMsg,RoutineName)
               CALL CleanUp()
               RETURN
            END IF
   
               ! move arrays back
            CALL MOVE_ALLOC(InitLocal%WAMIT%WaveTime,               p%WaveTime  ) 
            CALL MOVE_ALLOC(InitLocal%WAMIT%WaveElevC0, Waves_InitOut%WaveElevC0) 
            CALL MOVE_ALLOC(InitLocal%WAMIT%WaveDirArr, Waves_InitOut%WaveDirArr) 
               
   
               !-----------------------------------------
               ! Initialize the WAMIT2 Calculations
               !-----------------------------------------
   
               ! Only call the WAMIT2_Init if one of the flags is set for a calculation
            IF ( InitLocal%WAMIT2%MnDriftF .OR. InitLocal%WAMIT2%NewmanAppF .OR. InitLocal%WAMIT2%DiffQTFF .OR. InitLocal%WAMIT2%SumQTFF ) THEN
   
               
               InitLocal%WAMIT2%RhoXg       = Waves_InitOut%RhoXg
               InitLocal%WAMIT2%NStepWave   = Waves_InitOut%NStepWave
               InitLocal%WAMIT2%NStepWave2  = Waves_InitOut%NStepWave2
               InitLocal%WAMIT2%WaveDirMin  = Waves_InitOut%WaveDirMin
               InitLocal%WAMIT2%WaveDirMax  = Waves_InitOut%WaveDirMax
               InitLocal%WAMIT2%WaveDOmega  = Waves_InitOut%WaveDOmega
               
                  ! Temporarily move arrays to init input for WAMIT2 (save some space)
               CALL MOVE_ALLOC(p%WaveTime, InitLocal%WAMIT2%WaveTime) 
               CALL MOVE_ALLOC(Waves_InitOut%WaveElevC0, InitLocal%WAMIT2%WaveElevC0) 
               CALL MOVE_ALLOC(Waves_InitOut%WaveDirArr, InitLocal%WAMIT2%WaveDirArr) 
               
               
               CALL WAMIT2_Init(InitLocal%WAMIT2, m%u_WAMIT2, p%WAMIT2, x%WAMIT2, xd%WAMIT2, z%WAMIT2, OtherState%WAMIT2, &
                                       y%WAMIT2, m%WAMIT2, Interval, InitOut%WAMIT2, ErrStat2, ErrMsg2 )
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
               IF ( ErrStat >= AbortErrLev ) THEN
                  CALL CleanUp()
                  RETURN
               END IF
   
                  ! move arrays back
               CALL MOVE_ALLOC(InitLocal%WAMIT2%WaveTime,               p%WaveTime  ) 
               CALL MOVE_ALLOC(InitLocal%WAMIT2%WaveElevC0, Waves_InitOut%WaveElevC0) 
               CALL MOVE_ALLOC(InitLocal%WAMIT2%WaveDirArr, Waves_InitOut%WaveDirArr) 
   
               
                  ! Verify that WAMIT2_Init() did not request a different Interval!
   
               IF ( p%DT /= Interval ) THEN
                  CALL SetErrStat(ErrID_Fatal,'WAMIT2 Module attempted to change timestep interval, but this is not allowed.  '// &
                                             'WAMIT2 Module must use the HydroDyn Interval.',ErrStat,ErrMsg,RoutineName)
                  CALL CleanUp()
                  RETURN
               END IF
               
            ELSE
               
               p%WAMIT2%NumOuts = 0  !This doesn't get initialized if we don't call WAMIT2_Init
   
            ENDIF

#ifdef USE_FIT 
         ELSE IF ( InitLocal%PotMod == 2  ) THEN  ! FIT 
            ! Set up the Initialization data for FIT
               ! General
            FITInitData%InputFile      = InitLocal%PotFile
            FITInitData%Gravity        = InitLocal%Gravity
            FITInitData%Rho            = InitLocal%Waves%WtrDens
            FITInitData%time_end       = InitLocal%TMax
            FITInitData%dtime          = InitLocal%Waves%WaveDT  ! Set the FIT module's timestep equal to the WaveDT timestep, this was checked earlier to make sure it is an integer muliple of the glue-code timestep!
               ! Waves
               ! Need to pre-process the incoming wave data to be compatible with FIT
            
            FITInitData%N_omega        = Waves_InitOut%NStepWave2
            FITInitData%Wave_angle     = Waves_InitOut%WaveDir
            
               ! allocate waves data arrays for FIT
            CALL AllocAry( FITInitData%Wave_amp, FITInitData%N_omega, "Wave_amp", ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            CALL AllocAry( FITInitData%Wave_omega, FITInitData%N_omega, "Wave_omega", ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            CALL AllocAry( FITInitData%Wave_number, FITInitData%N_omega, "Wave_number", ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            CALL AllocAry( FITInitData%Wave_phase, FITInitData%N_omega, "Wave_phase", ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               IF ( ErrStat >= AbortErrLev ) THEN
                  CALL Cleanup()
                  RETURN
               END IF
               
               ! Populate wave arrays
            Np = 2*(Waves_InitOut%WaveDOmega + 1)
            DO I = 1 , Waves_InitOut%NStepWave2
               
               dftreal        = Waves_InitOut%WaveElevC0( 1,ABS(I ) )
               dftimag        = Waves_InitOut%WaveElevC0( 2, ABS(I ) )*SIGN(1,I)
               FITInitData%Wave_amp   (I) = sqrt( dftreal**2 + dftimag**2 )  * 2.0 / Np
               FITInitData%Wave_omega (I) = I*Waves_InitOut%WaveDOmega
               FITInitData%Wave_number(I) = I*Waves_InitOut%WaveDOmega**2. / InitLocal%Gravity
               FITInitData%Wave_phase (I) = atan2( dftimag, dftreal ) 
              
            END DO         
         
  
              ! Output
            FITInitData%RootName       = trim(InitLocal%OutRootName)//'.FIT'
                              
      
            CALL FIT_Init(FITInitData, u%FIT, p%FIT, FIT_x, xd%FIT, FIT_z, OtherState%FIT, y%FIT, Interval, FIT_InitOut, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
#endif

         END IF
   
   
      END IF  ! Check for WaveMod = 6


         ! Are there Morison elements?
       
      IF ( InitLocal%Morison%NMembers > 0 ) THEN

         
                ! Copy Waves initialization output into the initialization input type for the Morison module                              
         
         InitLocal%Morison%NStepWave    = Waves_InitOut%NStepWave
         
         
            ! Temporarily move array to init input for Morison (save some space)
         CALL MOVE_ALLOC( p%WaveTime,               InitLocal%Morison%WaveTime )
         
            ! Permanently move these wave values to Morison init input (and note they are potentially modified by 2nd order stuff before being sent to Morison)
         CALL MOVE_ALLOC( Waves_InitOut%WaveAcc,   InitLocal%Morison%WaveAcc )            
         CALL MOVE_ALLOC( Waves_InitOut%WaveDynP,  InitLocal%Morison%WaveDynP )         
         CALL MOVE_ALLOC( Waves_InitOut%WaveVel,   InitLocal%Morison%WaveVel )         
         CALL MOVE_ALLOC( Waves_InitOut%nodeInWater,InitLocal%Morison%nodeInWater )  ! moved to Morison%p%nodeInWater in the init routine


               ! If we did some second order wave kinematics corrections to the acceleration, velocity or
               ! dynamic pressure using the Waves2 module, then we need to add these to the values that we
               ! will be passing into the Morrison module.

            ! Difference frequency results
         IF ( p%Waves2%WvDiffQTFF ) THEN

               ! Dynamic pressure -- difference frequency terms
            IF ( SIZE(InitLocal%Morison%WaveDynP,DIM=1) /= SIZE(InitOut%Waves2%WaveDynP2D,DIM=1) .OR. &
                 SIZE(InitLocal%Morison%WaveDynP,DIM=2) /= SIZE(InitOut%Waves2%WaveDynP2D,DIM=2)) THEN
               CALL SetErrStat(ErrID_Fatal, &
                  ' WaveDynP arrays for first and second order wave elevations are of different sizes.  '//NewLine// &
                  'Morrison: '// TRIM(Num2LStr(SIZE(InitLocal%Morison%WaveDynP,DIM=1)))//'x'//          &
                                 TRIM(Num2LStr(SIZE(InitLocal%Morison%WaveDynP,DIM=2)))//NewLine//      &
                  'Waves2:   '// TRIM(Num2LStr(SIZE(InitOut%Waves2%WaveDynP2D,DIM=1)))//'x'//            &
                                 TRIM(Num2LStr(SIZE(InitOut%Waves2%WaveDynP2D,DIM=2))),                  &
                  ErrStat,ErrMsg,RoutineName)
               CALL CleanUp()
               RETURN
            ELSE
               InitLocal%Morison%WaveDynP = InitLocal%Morison%WaveDynP + InitOut%Waves2%WaveDynP2D
               IF (InitLocal%Waves%WaveStMod > 0 ) WaveDynP0 = WaveDynP0 + WaveDynP2D0
            ENDIF

               ! Particle velocity -- difference frequency terms
            IF ( SIZE(InitLocal%Morison%WaveVel,DIM=1) /= SIZE(InitOut%Waves2%WaveVel2D,DIM=1) .OR. &
                 SIZE(InitLocal%Morison%WaveVel,DIM=2) /= SIZE(InitOut%Waves2%WaveVel2D,DIM=2) .OR. &
                 SIZE(InitLocal%Morison%WaveVel,DIM=3) /= SIZE(InitOut%Waves2%WaveVel2D,DIM=3)) THEN
               CALL SetErrStat(ErrID_Fatal, &
                  ' WaveVel arrays for first and second order wave elevations are of different sizes.',ErrStat,ErrMsg,RoutineName)
               CALL CleanUp()
               RETURN
            ELSE
               InitLocal%Morison%WaveVel = InitLocal%Morison%WaveVel + InitOut%Waves2%WaveVel2D
               IF (InitLocal%Waves%WaveStMod > 0 ) WaveVel0 = WaveVel0 + WaveVel2D0
            ENDIF


               ! Particle acceleration -- difference frequency terms
            IF ( SIZE(InitLocal%Morison%WaveAcc,DIM=1) /= SIZE(InitOut%Waves2%WaveAcc2D,DIM=1) .OR. &
                 SIZE(InitLocal%Morison%WaveAcc,DIM=2) /= SIZE(InitOut%Waves2%WaveAcc2D,DIM=2) .OR. &
                 SIZE(InitLocal%Morison%WaveAcc,DIM=3) /= SIZE(InitOut%Waves2%WaveAcc2D,DIM=3)) THEN
               CALL SetErrStat(ErrID_Fatal, &
                  ' WaveAcc arrays for first and second order wave elevations are of different sizes.',ErrStat,ErrMsg,RoutineName)
               CALL CleanUp()
               RETURN
            ELSE
               InitLocal%Morison%WaveAcc = InitLocal%Morison%WaveAcc + InitOut%Waves2%WaveAcc2D
               IF (InitLocal%Waves%WaveStMod > 0 ) WaveAcc0 = WaveAcc0 + WaveAcc2D0
            ENDIF

         ENDIF ! second order wave kinematics difference frequency results

            ! Sum frequency results
         IF ( p%Waves2%WvSumQTFF ) THEN

               ! Dynamic pressure -- sum frequency terms
            IF ( SIZE(InitLocal%Morison%WaveDynP,DIM=1) /= SIZE(InitOut%Waves2%WaveDynP2S,DIM=1) .OR. &
                 SIZE(InitLocal%Morison%WaveDynP,DIM=2) /= SIZE(InitOut%Waves2%WaveDynP2S,DIM=2)) THEN
               CALL SetErrStat(ErrID_Fatal, &
                  ' WaveDynP arrays for first and second order wave elevations are of different sizes.  '//NewLine// &
                  'Morrison: '// TRIM(Num2LStr(SIZE(InitLocal%Morison%WaveDynP,DIM=1)))//'x'//          &
                                 TRIM(Num2LStr(SIZE(InitLocal%Morison%WaveDynP,DIM=2)))//NewLine//      &
                  'Waves2:   '// TRIM(Num2LStr(SIZE(InitOut%Waves2%WaveDynP2D,DIM=1)))//'x'//            &
                                 TRIM(Num2LStr(SIZE(InitOut%Waves2%WaveDynP2D,DIM=2))),                  &
                  ErrStat,ErrMsg,RoutineName)
               CALL CleanUp()
               RETURN
            ELSE
               InitLocal%Morison%WaveDynP = InitLocal%Morison%WaveDynP + InitOut%Waves2%WaveDynP2S
               IF (InitLocal%Waves%WaveStMod > 0 ) WaveDynP0 = WaveDynP0 + WaveDynP2S0
            ENDIF

               ! Particle velocity -- sum frequency terms
            IF ( SIZE(InitLocal%Morison%WaveVel,DIM=1) /= SIZE(InitOut%Waves2%WaveVel2S,DIM=1) .OR. &
                 SIZE(InitLocal%Morison%WaveVel,DIM=2) /= SIZE(InitOut%Waves2%WaveVel2S,DIM=2) .OR. &
                 SIZE(InitLocal%Morison%WaveVel,DIM=3) /= SIZE(InitOut%Waves2%WaveVel2S,DIM=3)) THEN
               CALL SetErrStat(ErrID_Fatal, &
                  ' WaveVel arrays for first and second order wave elevations are of different sizes.',ErrStat,ErrMsg,RoutineName)
               CALL CleanUp()
               RETURN
            ELSE
               InitLocal%Morison%WaveVel = InitLocal%Morison%WaveVel + InitOut%Waves2%WaveVel2S
               IF (InitLocal%Waves%WaveStMod > 0 ) WaveVel0 = WaveVel0 + WaveVel2S0
            ENDIF

               ! Particle velocity -- sum frequency terms
            IF ( SIZE(InitLocal%Morison%WaveAcc,DIM=1) /= SIZE(InitOut%Waves2%WaveAcc2S,DIM=1) .OR. &
                 SIZE(InitLocal%Morison%WaveAcc,DIM=2) /= SIZE(InitOut%Waves2%WaveAcc2S,DIM=2) .OR. &
                 SIZE(InitLocal%Morison%WaveAcc,DIM=3) /= SIZE(InitOut%Waves2%WaveAcc2S,DIM=3)) THEN
               CALL SetErrStat(ErrID_Fatal, &
                  ' WaveAcc arrays for first and second order wave elevations are of different sizes.',ErrStat,ErrMsg,RoutineName)
               CALL CleanUp()
               RETURN
            ELSE
               InitLocal%Morison%WaveAcc = InitLocal%Morison%WaveAcc + InitOut%Waves2%WaveAcc2S
               IF (InitLocal%Waves%WaveStMod > 0 ) WaveAcc0 = WaveAcc0 + WaveAcc2S0
            ENDIF

         ENDIF ! second order wave kinematics sum frequency results

!==============================================================================
         ! TODO: 1/29/2016 GJH
         ! This is where we need to perform Wave Stretching, now that the wave kinematics have been combined.
         ! We will call a new subroutine to perform this work. 
         ! As an input, this code need the kinematics at the (X,Y,0) location which in a Z-line above/below all the nodes where kinematics are computed.
         ! This code will alter the kinematics for stretching AND alter the nodeInWater array based on the combined wave elevation information
         IF (InitLocal%Waves%WaveStMod > 0 ) THEN
            call WvStretch_Init( InitLocal%Waves%WaveStMod, InitLocal%Waves%WtrDpth, InitLocal%Morison%NStepWave, InitLocal%Morison%NNodes,  &
                              p%NWaveElev, WaveElevSt, InitLocal%Waves%WaveKinzi, InitLocal%Morison%WaveTime, &
                              WaveVel0, WaveAcc0, WaveDynP0, &
                              Waves_InitOut%PWaveVel0, Waves_InitOut%PWaveAcc0, Waves_InitOut%PWaveDynP0, &
                              InitLocal%Morison%WaveVel, InitLocal%Morison%WaveAcc, InitLocal%Morison%WaveDynP, &
                              InitLocal%Morison%nodeInWater, ErrStat, ErrMsg )  
            DEALLOCATE(WaveElevSt)
            DEALLOCATE(WaveVel0)
            DEALLOCATE(WaveAcc0)
            DEALLOCATE(WaveDynP0)
         END IF
!==============================================================================
         ! In this version, this can only be TRUE if the precomiler flag WRITE_WV_KIN set and WaveMod not equal to 5 or 6 and WvKinFile is a valid string  
         IF ( ( InitLocal%Waves%WaveMod == 5 .OR. InitLocal%Waves%WaveMod == 6 ) .AND.  InitLocal%Echo ) THEN
            call HDOut_WriteWvKinFiles( TRIM(InitLocal%Waves%WvKinFile)//'_ech', HydroDyn_ProgDesc, InitLocal%Morison%NStepWave, InitLocal%Morison%NNodes,  &
                                             p%NWaveElev, InitLocal%Morison%nodeInWater, p%WaveElev, InitLocal%Waves%WaveKinzi, InitLocal%Morison%WaveTime, &
                                        InitLocal%Morison%WaveVel, InitLocal%Morison%WaveAcc, InitLocal%Morison%WaveDynP, &
                                        ErrStat, ErrMsg )  
         ELSE IF (InitLocal%Waves%WriteWvKin ) THEN
            call HDOut_WriteWvKinFiles( TRIM(InitLocal%Waves%WvKinFile), HydroDyn_ProgDesc, InitLocal%Morison%NStepWave, InitLocal%Morison%NNodes,  &
                                             p%NWaveElev, InitLocal%Morison%nodeInWater, p%WaveElev, InitLocal%Waves%WaveKinzi, InitLocal%Morison%WaveTime, &
                                        InitLocal%Morison%WaveVel, InitLocal%Morison%WaveAcc, InitLocal%Morison%WaveDynP, &
                                        ErrStat, ErrMsg )  
         END IF





            ! Check the output switch to see if Morison is needing to send outputs back to HydroDyn via the WriteOutput array
            
         IF ( InitLocal%OutSwtch > 0 ) THEN
            InitLocal%Morison%OutSwtch     = 2  ! only HydroDyn or the Driver code will write outputs to the file, that's why we are forcing this to 2.
         ELSE
            InitLocal%Morison%OutSwtch     = 0
         END IF
        
            ! Initialize the Morison Element Calculations 
      
         CALL Morison_Init(InitLocal%Morison, u%Morison, p%Morison, x%Morison, xd%Morison, z%Morison, OtherState%Morison, &
                               y%Morison, m%Morison, Interval, InitOut%Morison, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF
         
            ! move array back
         CALL MOVE_ALLOC( InitLocal%Morison%WaveTime, p%WaveTime  )
         
         
         IF ( u%Morison%DistribMesh%Committed ) THEN
                  ! we need the translation displacement mesh for loads transfer:
            CALL MeshCopy ( SrcMesh  = u%Morison%DistribMesh            &
                    , DestMesh = m%MrsnDistribMesh_position   &
                    , CtrlCode = MESH_NEWCOPY        &
                    , IOS      = COMPONENT_INPUT     &
                    , TranslationDisp = .TRUE.       &
                    , ErrStat  = ErrStat2            &
                    , ErrMess  = ErrMsg2              )  ! automatically sets    DestMesh%RemapFlag = .TRUE.
                    
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'HydroDyn_Init:m%MrsnDistribMesh_position')
               IF ( ErrStat >= AbortErrLev ) THEN
                  CALL CleanUp()
                  RETURN
               END IF
            m%MrsnDistribMesh_position%TranslationDisp = 0.0  ! bjj: this is actually initialized in the ModMesh module, but I'll do it here anyway.
            
         END IF
         
         IF ( u%Morison%LumpedMesh%Committed ) THEN
                  ! we need the translation displacement mesh for loads transfer:
            CALL MeshCopy ( SrcMesh  = u%Morison%LumpedMesh           &
                    , DestMesh = m%MrsnLumpedMesh_position   &
                    , CtrlCode = MESH_NEWCOPY        &
                    , IOS      = COMPONENT_INPUT     &
                    , TranslationDisp = .TRUE.       &
                    , ErrStat  = ErrStat2            &
                    , ErrMess  = ErrMsg2             )  ! automatically sets    DestMesh%RemapFlag = .TRUE.
                    
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'HydroDyn_Init:m%MrsnLumpedMesh_position')
            IF ( ErrStat >= AbortErrLev ) THEN
               CALL CleanUp()
               RETURN
            END IF
            m%MrsnLumpedMesh_position%TranslationDisp = 0.0  ! bjj: this is actually initialized in the ModMesh module, but I'll do it here anyway.
            
         END IF
            ! Verify that Morison_Init() did not request a different Interval!
      
         IF ( p%DT /= Interval ) THEN
            CALL SetErrStat(ErrID_Fatal,'Morison Module attempted to change timestep interval, but this is not allowed.  Morison Module must use the HydroDyn Interval.',ErrStat,ErrMsg,RoutineName)
            CALL CleanUp()
            RETURN
         END IF
         
      END IF  ! ( InitLocal%Morison%NMembers > 0 )
    
!===============================================
      p%PotMod = InitLocal%Potmod      
      IF ( InitLocal%UnSum > 0 ) THEN
      
         IF (InitLocal%Waves%WaveMod /= 0 .AND. InitLocal%Waves%WaveMod /= 6)  THEN
               ! Write the header for this section
            WRITE( InitLocal%UnSum,  '(//)' )         
            WRITE( InitLocal%UnSum, '(1X,A15)' )   'Wave Kinematics'
            WRITE( InitLocal%UnSum,  '(/)' )
            WRITE( InitLocal%UnSum, '(1X,A10,2X,A14,2X,A14,2X,A14,2X,A19,2X,A19)' )  &
                     '    m   ', '    k    ', '   Omega[m]  ', '   Direction  ', 'REAL(DFT{WaveElev})','IMAG(DFT{WaveElev})'
            WRITE( InitLocal%UnSum, '(1X,A10,2X,A14,2X,A14,2X,A14,2X,A19,2X,A19)' )  &
                     '   (-)  ', '  (1/m)  ', '   (rad/s)   ', '     (deg)    ', '       (m)         ','       (m)         '

            ! Write the data
            DO I = -1*Waves_InitOut%NStepWave2+1,Waves_InitOut%NStepWave2
               WaveNmbr   = WaveNumber ( I*Waves_InitOut%WaveDOmega, InitLocal%Gravity, InitLocal%Waves%WtrDpth )
               WRITE( InitLocal%UnSum, '(1X,I10,2X,ES14.5,2X,ES14.5,2X,ES14.5,2X,ES14.5,7X,ES14.5)' ) I, WaveNmbr, I*Waves_InitOut%WaveDOmega, &
                      Waves_InitOut%WaveDirArr(ABS(I)),  Waves_InitOut%WaveElevC0( 1,ABS(I ) ) ,   Waves_InitOut%WaveElevC0( 2, ABS(I ) )*SIGN(1,I)
            END DO
         END IF
         
         
         IF ( InitLocal%PotMod == 1 .AND.  InitLocal%WAMIT%RdtnMod == 1) THEN
            ! Write the header for this section
            WRITE( InitLocal%UnSum,  '(//)' ) 
            WRITE( InitLocal%UnSum,  '(A)' ) 'Radiation memory effect kernel'
            WRITE( InitLocal%UnSum,  '(//)' ) 
            WRITE( InitLocal%UnSum, '(1X,A10,2X,A10,21(2X,A16))' )    '    n    ' , '     t    ', '   K11    ', '   K12    ', '    K13   ', '    K14    ', '    K15    ', '    K16    ', '    K22   ', '    K23   ', '    K24    ', '    K25    ', '    K26    ', '    K33    ', '    K34    ', '    K35    ',     'K36    ', '    K44    ', '    K45    ', '    K46    ', '    K55    ', '    K56    ', '    K66    '
            WRITE( InitLocal%UnSum, '(1X,A10,2X,A10,21(2X,A16))' )    '   (-)   ' , '    (s)   ', ' (kg/s^2) ', ' (kg/s^2) ', ' (kg/s^2) ', ' (kgm/s^2) ', ' (kgm/s^2) ', ' (kgm/s^2) ', ' (kg/s^2) ', ' (kg/s^2) ', ' (kgm/s^2) ', ' (kgm/s^2) ', ' (kgm/s^2) ', ' (kg/s^2)  ', ' (kgm/s^2) ', ' (kgm/s^2) ', ' (kgm/s^2) ', '(kgm^2/s^2)', '(kgm^2/s^2)', '(kgm^2/s^2)', '(kgm^2/s^2)', '(kgm^2/s^2)', '(kgm^2/s^2)'

               ! Write the data
            DO I = 0,p%WAMIT%Conv_Rdtn%NStepRdtn-1
   
               WRITE( InitLocal%UnSum, '(1X,I10,2X,E12.5,21(2X,ES16.5))' ) I, I*p%WAMIT%Conv_Rdtn%RdtnDT, p%WAMIT%Conv_Rdtn%RdtnKrnl(I,1,1), p%WAMIT%Conv_Rdtn%RdtnKrnl(I,1,2), p%WAMIT%Conv_Rdtn%RdtnKrnl(I,1,3), p%WAMIT%Conv_Rdtn%RdtnKrnl(I,1,4), p%WAMIT%Conv_Rdtn%RdtnKrnl(I,1,5), p%WAMIT%Conv_Rdtn%RdtnKrnl(I,1,6), p%WAMIT%Conv_Rdtn%RdtnKrnl(I,2,2), p%WAMIT%Conv_Rdtn%RdtnKrnl(I,2,3), p%WAMIT%Conv_Rdtn%RdtnKrnl(I,2,4), p%WAMIT%Conv_Rdtn%RdtnKrnl(I,2,5), p%WAMIT%Conv_Rdtn%RdtnKrnl(I,2,6), p%WAMIT%Conv_Rdtn%RdtnKrnl(I,3,3), p%WAMIT%Conv_Rdtn%RdtnKrnl(I,3,4), p%WAMIT%Conv_Rdtn%RdtnKrnl(I,3,5), p%WAMIT%Conv_Rdtn%RdtnKrnl(I,3,6), p%WAMIT%Conv_Rdtn%RdtnKrnl(I,4,4), p%WAMIT%Conv_Rdtn%RdtnKrnl(I,4,5), p%WAMIT%Conv_Rdtn%RdtnKrnl(I,4,6), p%WAMIT%Conv_Rdtn%RdtnKrnl(I,5,5), p%WAMIT%Conv_Rdtn%RdtnKrnl(I,5,6), p%WAMIT%Conv_Rdtn%RdtnKrnl(I,6,6)
      
            END DO
         END IF
         
      END IF

!==========================================
      
         ! Deallocate any remaining Waves Output data
      IF(ALLOCATED( Waves_InitOut%WaveElevC0 ))  DEALLOCATE( Waves_InitOut%WaveElevC0 )
      IF(ALLOCATED( Waves_InitOut%WaveAcc   ))  DEALLOCATE( Waves_InitOut%WaveAcc   )
      IF(ALLOCATED( Waves_InitOut%WaveDynP  ))  DEALLOCATE( Waves_InitOut%WaveDynP  )
      IF(ALLOCATED( Waves_InitOut%WaveTime   ))  DEALLOCATE( Waves_InitOut%WaveTime   )
      IF(ALLOCATED( Waves_InitOut%WaveVel   ))  DEALLOCATE( Waves_InitOut%WaveVel   )
      IF(ALLOCATED( Waves_InitOut%WaveElevC0 ))  DEALLOCATE( Waves_InitOut%WaveElevC0 )
      !IF(ALLOCATED( InitLocal%WAMIT%WaveElevC0 ))  DEALLOCATE( InitLocal%WAMIT%WaveElevC0)
      
         ! Close the summary file
      IF ( InitLocal%HDSum ) THEN
         CALL HDOut_CloseSum( InitLocal%UnSum, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) THEN
               CALL CleanUp()
               RETURN
            END IF
      END IF
      
      ! Define system output initializations (set up mesh) here:
      
      
          ! Create the input and output meshes associated with lumped load at the WAMIT reference point (WRP)
      
      CALL MeshCreate( BlankMesh        = u%Mesh            &
                     ,IOS               = COMPONENT_INPUT   &
                     ,Nnodes            = 1                 &
                     ,ErrStat           = ErrStat2          &
                     ,ErrMess           = ErrMsg2           &
                     ,TranslationDisp   = .TRUE.            &
                     ,Orientation       = .TRUE.            &
                     ,TranslationVel    = .TRUE.            &
                     ,RotationVel       = .TRUE.            &
                     ,TranslationAcc    = .TRUE.            &
                     ,RotationAcc       = .TRUE.)
         ! Create the node on the mesh
            
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF
         
      CALL MeshPositionNode (u%Mesh                                &
                              , 1                                  &
                              , (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi/)   &  
                              , ErrStat2                           &
                              , ErrMsg2                            )
      
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF
       
      
         ! Create the mesh element
      CALL MeshConstructElement (  u%Mesh              &
                                  , ELEMENT_POINT      &                         
                                  , ErrStat2           &
                                  , ErrMsg2            &
                                  , 1                  &
                                              )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF
      
      
      CALL MeshCommit ( u%Mesh   &
                      , ErrStat2            &
                      , ErrMsg2             )
   
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF
      

         
      CALL MeshCopy (   SrcMesh      = u%Mesh               &
                     ,DestMesh     = y%Mesh                 &
                     ,CtrlCode     = MESH_SIBLING           &
                     ,IOS          = COMPONENT_OUTPUT       &
                     ,ErrStat      = ErrStat2               &
                     ,ErrMess      = ErrMsg2                &
                     ,Force        = .TRUE.                 &
                     ,Moment       = .TRUE.                 )
     
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'HydroDyn_Init:y%Mesh')
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF      
      u%Mesh%RemapFlag  = .TRUE.
      y%Mesh%RemapFlag  = .TRUE.
     
     CALL MeshCopy (   SrcMesh     = y%Mesh                 &
                     ,DestMesh     = y%AllHdroOrigin        &
                     ,CtrlCode     = MESH_NEWCOPY           &
                     ,IOS          = COMPONENT_OUTPUT       &
                     ,ErrStat      = ErrStat2               &
                     ,ErrMess      = ErrMsg2                &
                     ,Force        = .TRUE.                 &
                     ,Moment       = .TRUE.                 )
     
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'HydroDyn_Init:y%AllHdroOrigin')
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF      
      y%AllHdroOrigin%RemapFlag  = .TRUE.
      
         ! we need the translation displacement mesh for loads transfer:
      CALL MeshCopy ( SrcMesh  = u%Mesh            &
                    , DestMesh = m%AllHdroOrigin_position   &
                    , CtrlCode = MESH_NEWCOPY        &
                    , IOS      = COMPONENT_INPUT     &
                    , TranslationDisp = .TRUE.       &
                    , ErrStat  = ErrStat2            &
                    , ErrMess  = ErrMsg2             )  ! automatically sets    DestMesh%RemapFlag = .TRUE.
                    
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF
      m%AllHdroOrigin_position%TranslationDisp = 0.0  ! bjj: this is actually initialized in the ModMesh module, but I'll do it here anyway.
      
     
         ! Create the Output file if requested
      
      p%OutSwtch      = InitLocal%OutSwtch 
      p%Delim         = ''
      !p%Morison%Delim = p%Delim  ! Need to set this from within Morison to follow framework
      !p%WAMIT%Delim   = p%Delim  ! Need to set this from within Morison to follow framework
      p%OutFmt        = InitLocal%OutFmt
      p%OutSFmt       = InitLocal%OutSFmt
      p%NumOuts       = InitLocal%NumOuts
      
      CALL HDOUT_Init( HydroDyn_ProgDesc, InitLocal, y,  p, m, InitOut, ErrStat2, ErrMsg2 )
      
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF
      
         ! Create some mesh mapping data
      CALL MeshCopy ( SrcMesh      = y%Mesh                 &
                     ,DestMesh     = m%y_mapped             &
                     ,CtrlCode     = MESH_NEWCOPY           &
                     ,IOS          = COMPONENT_OUTPUT       &
                     ,ErrStat      = ErrStat2               &
                     ,ErrMess      = ErrMsg2                &
                     ,Force        = .TRUE.                 &
                     ,Moment       = .TRUE.                 )
          
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      m%y_mapped%RemapFlag  = .TRUE.
 
      CALL MeshMapCreate( y%Mesh,                m%y_mapped, m%HD_MeshMap%HD_P_2_WRP_P, ErrStat2, ErrMsg2  );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      IF ( y%Morison%LumpedMesh%Committed ) THEN 
         CALL MeshMapCreate( y%Morison%LumpedMesh,  m%y_mapped, m%HD_MeshMap%M_P_2_WRP_P,  ErrStat2, ErrMsg2  );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      ENDIF
      IF ( y%Morison%DistribMesh%Committed ) THEN 
         CALL MeshMapCreate( y%Morison%DistribMesh, m%y_mapped, m%HD_MeshMap%M_L_2_WRP_P,  ErrStat2, ErrMsg2  );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      ENDIF
      
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF         
         
         ! Define initialization-routine output here:
         InitOut%Ver = HydroDyn_ProgDesc         
            ! These three come directly from processing the inputs, and so will exist even if not using Morison elements:
         InitOut%WtrDens = InitLocal%Morison%WtrDens
         InitOut%WtrDpth = InitLocal%Morison%WtrDpth
         InitOut%MSL2SWL = InitLocal%Morison%MSL2SWL
                                                                   
      IF ( InitInp%hasIce ) THEN
         IF ((InitLocal%Waves%WaveMod /= 0) .OR. (InitLocal%Current%CurrMod /= 0) ) THEN
            CALL SetErrStat(ErrID_Fatal,'Waves and Current must be turned off in HydroDyn when ice loading is computed. Set WaveMod=0 and CurrMod=0.',ErrStat,ErrMsg,RoutineName)
         END IF
      END IF
      
            
         ! Destroy the local initialization data
      CALL CleanUp()
         
CONTAINS
!................................
   SUBROUTINE CleanUp()
      
      CALL HydroDyn_DestroyInitInput( InitLocal,       ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL Waves_DestroyInitOutput(   Waves_InitOut,   ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) 
      CALL Current_DestroyInitOutput( Current_InitOut, ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) 
   
      
         ! These are dummy variables to satisfy the framework, but are not used again:
      
      CALL Waves_DestroyInput(       Waves_u,          ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      
      CALL Waves_DestroyParam(       Waves_p,          ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      
      CALL Waves_DestroyContState(   Waves_x,          ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      
      CALL Waves_DestroyDiscState(   Waves_xd,         ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      
      CALL Waves_DestroyConstrState( Waves_z,          ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      
      CALL Waves_DestroyOtherState(  WavesOtherState,  ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      
      CALL Waves_DestroyOutput(      Waves_y,          ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      

      
      CALL Current_DestroyInput(       Current_u,          ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      
      CALL Current_DestroyParam(       Current_p,          ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      
      CALL Current_DestroyContState(   Current_x,          ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      
      CALL Current_DestroyDiscState(   Current_xd,         ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      
      CALL Current_DestroyConstrState( Current_z,          ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      
      CALL Current_DestroyOtherState(  CurrentOtherState,  ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      
      CALL Current_DestroyOutput(      Current_y,          ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      
      CALL Current_DestroyMisc(        Current_m,          ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      
      
            
   END SUBROUTINE CleanUp
!................................
END SUBROUTINE HydroDyn_Init


!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the end of the simulation.
SUBROUTINE HydroDyn_End( u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )

      TYPE(HydroDyn_InputType),           INTENT(INOUT)  :: u           !< System inputs
      TYPE(HydroDyn_ParameterType),       INTENT(INOUT)  :: p           !< Parameters     
      TYPE(HydroDyn_ContinuousStateType), INTENT(INOUT)  :: x           !< Continuous states
      TYPE(HydroDyn_DiscreteStateType),   INTENT(INOUT)  :: xd          !< Discrete states
      TYPE(HydroDyn_ConstraintStateType), INTENT(INOUT)  :: z           !< Constraint states
      TYPE(HydroDyn_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Other/optimization states            
      TYPE(HydroDyn_OutputType),          INTENT(INOUT)  :: y           !< System outputs
      TYPE(HydroDyn_MiscVarType),         INTENT(INOUT)  :: m           !< Initial misc/optimization variables           
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None



         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      
         ! Place any last minute operations or calculations here:


            
         ! Write the HydroDyn-level output file data if the user requested module-level output
         ! and the current time has advanced since the last stored time step.
         
      IF ( p%OutSwtch == 1 .OR. p%OutSwtch == 3) THEN               
         CALL HDOut_WriteOutputs( m%LastOutTime, y, p, m%Decimate, ErrStat, ErrMsg )         
      END IF          
      
         ! Close files here:  
      CALL HDOut_CloseOutput( p, ErrStat, ErrMsg )           
          

         ! Destroy the input data:
         
      CALL HydroDyn_DestroyInput( u, ErrStat, ErrMsg )


         ! Destroy the parameter data:
      
      CALL HydroDyn_DestroyParam( p, ErrStat, ErrMsg )


         ! Destroy the state data:
         
      CALL HydroDyn_DestroyContState(   x,           ErrStat, ErrMsg )
      CALL HydroDyn_DestroyDiscState(   xd,          ErrStat, ErrMsg )
      CALL HydroDyn_DestroyConstrState( z,           ErrStat, ErrMsg )
      CALL HydroDyn_DestroyOtherState(  OtherState,  ErrStat, ErrMsg )
         
         ! Destroy misc variables:
      
      CALL HydroDyn_DestroyMisc( m, ErrStat, ErrMsg )

         ! Destroy the output data:
         
      CALL HydroDyn_DestroyOutput( y, ErrStat, ErrMsg )
      

END SUBROUTINE HydroDyn_End


!----------------------------------------------------------------------------------------------------------------------------------
!> Loose coupling routine for solving constraint states, integrating continuous states, and updating discrete states.
!! Continuous, constraint, and discrete states are updated to values at t + Interval.
SUBROUTINE HydroDyn_UpdateStates( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )

      REAL(DbKi),                         INTENT(IN   )  :: t               !< Current simulation time in seconds
      INTEGER(IntKi),                     INTENT(IN   )  :: n               !< Current step of the simulation: t = n*Interval
      TYPE(HydroDyn_InputType),           INTENT(INOUT ) :: Inputs(:)       !< Inputs at InputTimes
      REAL(DbKi),                         INTENT(IN   )  :: InputTimes(:)   !< Times in seconds associated with Inputs
      TYPE(HydroDyn_ParameterType),       INTENT(IN   )  :: p               !< Parameters
      TYPE(HydroDyn_ContinuousStateType), INTENT(INOUT)  :: x               !< Input: Continuous states at t;
                                                                            !!   Output: Continuous states at t + Interval
      TYPE(HydroDyn_DiscreteStateType),   INTENT(INOUT)  :: xd              !< Input: Discrete states at t;
                                                                            !!   Output: Discrete states at t + Interval
      TYPE(HydroDyn_ConstraintStateType), INTENT(INOUT)  :: z               !< Input: Constraint states at t;
                                                                            !!   Output: Constraint states at t + Interval
      TYPE(HydroDyn_OtherStateType),      INTENT(INOUT)  :: OtherState      !< Other states: Other states at t;
                                                                            !!   Output: Other states at t + Interval
      TYPE(HydroDyn_MiscVarType),         INTENT(INOUT)  :: m               !< Initial misc/optimization variables           
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat         !< Error status of the operation
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg          !< Error message if ErrStat /= ErrID_None

         ! Local variables
      INTEGER                                            :: I               ! Generic loop counter
      TYPE(HydroDyn_ContinuousStateType)                 :: dxdt            ! Continuous state derivatives at t
      TYPE(HydroDyn_DiscreteStateType)                   :: xd_t            ! Discrete states at t (copy)
      TYPE(HydroDyn_ConstraintStateType)                 :: z_Residual      ! Residual of the constraint state functions (Z)
      TYPE(HydroDyn_InputType)                           :: u               ! Instantaneous inputs
      INTEGER(IntKi)                                     :: ErrStat2        ! Error status of the operation (secondary error)
      CHARACTER(ErrMsgLen)                               :: ErrMsg2         ! Error message if ErrStat2 /= ErrID_None
      INTEGER                                            :: nTime           ! number of inputs 

      TYPE(WAMIT_InputType), ALLOCATABLE                 :: Inputs_WAMIT(:)  
      CHARACTER(*), PARAMETER                            :: RoutineName = 'HydroDyn_UpdateStates'
      
          ! Create dummy variables required by framework but which are not used by the module
            
#ifdef USE_FIT      
      TYPE(FIT_InputType), ALLOCATABLE                 :: Inputs_FIT(:) 
      TYPE(FIT_ConstraintStateType)      :: FIT_z              ! constraint states
      TYPE(FIT_ContinuousStateType)      :: FIT_x              ! Input: Continuous states at t;
#endif      
      
      REAL(ReKi)                         :: rotdisp(3)
         ! Initialize variables

      ErrStat   = ErrID_None           ! no error has occurred
      ErrMsg    = ""
      
      
         
         ! Return without doing any work if the we are not using a potential flow model
      IF ( p%PotMod == 0  ) RETURN
      
      ! Return without doing any work if the input mesh is not initialized (NOT USING WAMIT)
      !IF ( .NOT. Inputs(1)%WAMIT%Mesh%Initialized  ) RETURN
      
      nTime = size(Inputs)   
      
      
         ! Allocate array of WAMIT inputs
         ! TODO: We should avoid allocating this at each time step if we can!
         
!FIXME: Error handling appears to be broken here

   IF ( p%PotMod == 1 ) THEN
      ALLOCATE( Inputs_WAMIT(nTime), STAT = ErrStat2 )
      IF (ErrStat2 /=0) THEN
         CALL SetErrStat( ErrID_Fatal, 'Failed to allocate array Inputs_WAMIT.', ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF

         
         ! Loop over number of inputs and copy them into an array of WAMIT inputs
      
      DO I=1,nTime
                  
            ! Copy the inputs from the HD mesh into the WAMIT mesh         
         CALL MeshCopy( Inputs(I)%Mesh, Inputs_WAMIT(I)%Mesh, MESH_NEWCOPY, ErrStat2, ErrMsg2 )   
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )                  
         
      END DO
      
         
      IF (ErrStat < AbortErrLev) THEN    ! if there was an error copying the input meshes, we'll skip this step and then cleanup the temporary input meshes     
            ! Update the WAMIT module states
      
         CALL WAMIT_UpdateStates( t, n, Inputs_WAMIT, InputTimes, p%WAMIT, x%WAMIT, xd%WAMIT, z%WAMIT, OtherState%WAMIT, m%WAMIT, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )         
     
      END IF
      
         ! deallocate temporary inputs
      DO I=1,nTime
         CALL WAMIT_DestroyInput( Inputs_WAMIT(I), ErrStat2, ErrMsg2 )     
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )         
      END DO
      
      DEALLOCATE(Inputs_WAMIT)

#ifdef USE_FIT      
   ELSE IF ( p%PotMod == 2 ) THEN  ! FIT
      
      ALLOCATE( Inputs_FIT(nTime), STAT = ErrStat2 )
      IF (ErrStat2 /=0) THEN
         CALL SetErrStat( ErrID_Fatal, 'Failed to allocate array Inputs_FIT.', ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF

         
         ! Loop over number of inputs and copy them into an array of FIT inputs
      
      DO I=1,nTime
         
            ! Copy the inputs from the HD mesh into the FIT input variables
         
            ! Determine the rotational angles from the direction-cosine matrix
         rotdisp = GetSmllRotAngs ( Inputs(I)%Mesh%Orientation(:,:,1), ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )   
         Inputs_FIT(I)%roll     = rotdisp(1)
         Inputs_FIT(I)%pitch    = rotdisp(2)
         Inputs_FIT(I)%yaw      = rotdisp(3)
         Inputs_FIT(I)%si_t(:)  = Inputs(I)%Mesh%TranslationDisp(:,1)             
         Inputs_FIT(I)%vel_t(:) = Inputs(I)%Mesh%TranslationVel (:,1)  
      END DO
      
         
         
         ! Update the FIT module states
     
      CALL FIT_UpdateStates( t, n, Inputs_FIT, InputTimes, p%FIT, FIT_x, xd%FIT, FIT_z, OtherState%FIT, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )         
     

         ! deallocate temporary inputs
      DO I=1,nTime
         CALL FIT_DestroyInput( Inputs_FIT(I), ErrStat2, ErrMsg2 )     
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )         
      END DO
      
      DEALLOCATE(Inputs_FIT) 
#endif

   END IF
   
      
END SUBROUTINE HydroDyn_UpdateStates


!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for computing outputs, used in both loose and tight coupling.
SUBROUTINE HydroDyn_CalcOutput( Time, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )   
   
      REAL(DbKi),                         INTENT(IN   )  :: Time        !< Current simulation time in seconds
      REAL(DbKi)                                         :: F_FHA(6,1)       !! FHA force
      TYPE(HydroDyn_InputType),           INTENT(INOUT)  :: u           !< Inputs at Time (note that this is intent out because we're copying the u%mesh into m%u_wamit%mesh)
      TYPE(HydroDyn_ParameterType),       INTENT(IN   )  :: p           !< Parameters
      TYPE(HydroDyn_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at Time
      TYPE(HydroDyn_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at Time
      TYPE(HydroDyn_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at Time
      TYPE(HydroDyn_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states at Time
      TYPE(HydroDyn_OutputType),          INTENT(INOUT)  :: y           !< Outputs computed at Time (Input only so that mesh con-
                                                                        !!   nectivity information does not have to be recalculated)
      TYPE(HydroDyn_MiscVarType),         INTENT(INOUT)  :: m           !< Initial misc/optimization variables           
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat     !! Error status of the operation
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg      !! Error message if ErrStat /= ErrID_None

      INTEGER                                            :: I, J        ! Generic counters
      
      INTEGER(IntKi)                                     :: ErrStat2        ! Error status of the operation (secondary error)
      CHARACTER(ErrMsgLen)                               :: ErrMsg2         ! Error message if ErrStat2 /= ErrID_None

#ifdef USE_FIT       
      TYPE(FIT_ContinuousStateType)        :: FIT_x             ! Initial continuous states
      TYPE(FIT_ConstraintStateType)        :: FIT_z             ! Initial guess of the constraint states 
      TYPE(FIT_InputType)                  :: Inputs_FIT
#endif      
      REAL(ReKi)                           :: WaveElev (p%NWaveElev) ! Instantaneous total elevation of incident waves at each of the NWaveElev points where the incident wave elevations can be output (meters)
      REAL(ReKi)                           :: WaveElev1(p%NWaveElev)    ! Instantaneous first order elevation of incident waves at each of the NWaveElev points where the incident wave elevations can be output (meters)
      
      REAL(ReKi)                           :: q(6), qdot(6), qdotsq(6), qdotdot(6)
      REAL(ReKi)                           :: rotdisp(3)                              ! small angle rotational displacements
      REAL(ReKi)                           :: AllOuts(MaxHDOutputs)  
      
      
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      
         ! Compute outputs here:
         
         
         !-------------------------------------------------------------------
         ! Additional stiffness, damping forces.  These need to be placed on a point mesh which is located at the WAMIT reference point (WRP).
         ! This mesh will need to get mapped by the glue code for use by either ElastoDyn or SubDyn.
         !-------------------------------------------------------------------
!bjj: if these are false in the input file, the parameter verions of these variables don't get set:

         ! Deal with any output from the Waves2 module....
      IF (p%Waves2%WvDiffQTFF .OR. p%Waves2%WvSumQTFF ) THEN

            ! Waves2_CalcOutput is called only so that the wave elevations can be output (if requested).
         CALL Waves2_CalcOutput( Time, m%u_Waves2, p%Waves2, x%Waves2, xd%Waves2,  &
                                z%Waves2, OtherState%Waves2, y%Waves2, m%Waves2, ErrStat2, ErrMsg2 )

         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )                  
      END IF

!FIXME: Error handling appears to be broken here.

         ! Determine the rotational angles from the direction-cosine matrix
      rotdisp = GetSmllRotAngs ( u%Mesh%Orientation(:,:,1), ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )                  

      q         = reshape((/REAL(u%Mesh%TranslationDisp(:,1),ReKi),rotdisp(:)/),(/6/))
      qdot      = reshape((/u%Mesh%TranslationVel(:,1),u%Mesh%RotationVel(:,1)/),(/6/))
      qdotsq    = abs(qdot)*qdot
      qdotdot   = reshape((/u%Mesh%TranslationAcc(:,1),u%Mesh%RotationAcc(:,1)/),(/6/))
      
      
         ! Compute the load contirbution from user-supplied added stiffness and damping
         
      m%F_PtfmAdd = p%AddF0 - matmul(p%AddCLin, q) - matmul(p%AddBLin, qdot) - matmul(p%AddBQuad, qdotsq)
      
         ! Attach to the output point mesh
      y%Mesh%Force (:,1) = m%F_PtfmAdd(1:3)
      y%Mesh%Moment(:,1) = m%F_PtfmAdd(4:6)
    
         ! Compute the wave elevations at the requested output locations for this time.  Note that p%WaveElev has the second order added to it already.
         ! DZ: moved up so we can use to control TMD params
      DO I=1,p%NWaveElev   
         WaveElev1(I)   = InterpWrappedStpReal ( REAL(Time, SiKi), p%WaveTime(:), p%WaveElev1(:,I),          &
                                    m%LastIndWave, p%NStepWave + 1       )                      
         WaveElev(I)    = InterpWrappedStpReal ( REAL(Time, SiKi), p%WaveTime(:), p%WaveElev(:,I), &
                                    m%LastIndWave, p%NStepWave + 1       )

      END DO


         ! CKA 3/9/2018 - Start Change - Apply fluid harmonic absorber loading as external force
        CALL FHA_Force(Time, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg,F_FHA,q,qdot,qdotdot,WaveElev(1)) 
          !write(*,*) Time, F_FHA(3,1)       
        y%Mesh%Force (:,1) =y%Mesh%Force (:,1) + F_FHA(1:3,1)
        y%Mesh%Moment(:,1) =y%Mesh%Moment(:,1) + F_FHA(4:6,1)
         ! End Change
      
      IF ( p%PotMod == 1 ) THEN
         IF ( m%u_WAMIT%Mesh%Committed ) THEN  ! Make sure we are using WAMIT / there is a valid mesh
         
               ! Copy the inputs from the HD mesh into the WAMIT mesh
            CALL MeshCopy( u%Mesh, m%u_WAMIT%Mesh, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )   
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )                  
               IF ( ErrStat >= AbortErrLev ) RETURN
         
         
            CALL WAMIT_CalcOutput( Time, m%u_WAMIT, p%WAMIT, x%WAMIT, xd%WAMIT,  &
                                   z%WAMIT, OtherState%WAMIT, y%WAMIT, m%WAMIT, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )                  
         
               ! Add WAMIT forces to the HydroDyn output mesh
            y%Mesh%Force (:,1) = y%Mesh%Force (:,1) + y%WAMIT%Mesh%Force (:,1)
            y%Mesh%Moment(:,1) = y%Mesh%Moment(:,1) + y%WAMIT%Mesh%Moment(:,1)
         

               ! Copy the F_Waves1 information to the HydroDyn level so we can combine it with the 2nd order
            m%F_Waves   = m%WAMIT%F_Waves1

         
         END IF
        
#ifdef USE_FIT          
      ELSE IF ( p%PotMod ==2 ) THEN !FIT
         Inputs_FIT%roll     = rotdisp(1)
         Inputs_FIT%pitch    = rotdisp(2)
         Inputs_FIT%yaw      = rotdisp(3)
         Inputs_FIT%si_t(:)  = u%Mesh%TranslationDisp(:,1)             
         Inputs_FIT%vel_t(:) = u%Mesh%TranslationVel (:,1)  
         CALL FIT_CalcOutput( Time, Inputs_FIT, p%FIT, FIT_x, xd%FIT, FIT_z, OtherState%FIT, y%FIT, ErrStat2, ErrMsg2 ) 
         
            ! Add FIT forces to the HydroDyn output mesh
         y%Mesh%Force (:,1) = y%Mesh%Force (:,1) + y%FIT%F(:)
         y%Mesh%Moment(:,1) = y%Mesh%Moment(:,1) + y%FIT%M(:)
#endif  
         
      END IF
      

      IF ( m%u_WAMIT2%Mesh%Committed ) THEN  ! Make sure we are using WAMIT2 / there is a valid mesh

            ! Copy the inputs from the HD mesh into the WAMIT2 mesh
         CALL MeshCopy( u%Mesh, m%u_WAMIT2%Mesh, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )                  
            IF ( ErrStat >= AbortErrLev ) RETURN


         CALL WAMIT2_CalcOutput( Time, m%u_WAMIT2, p%WAMIT2, x%WAMIT2, xd%WAMIT2,  &
                                z%WAMIT2, OtherState%WAMIT2, y%WAMIT2, m%WAMIT2, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )                  

            ! Add WAMIT2 forces to the HydroDyn output mesh
         y%Mesh%Force (:,1) = y%Mesh%Force (:,1) + y%WAMIT2%Mesh%Force (:,1)
         y%Mesh%Moment(:,1) = y%Mesh%Moment(:,1) + y%WAMIT2%Mesh%Moment(:,1)

            ! Add the second order WAMIT forces to the first order WAMIT forces for the total (this is just to make the mesh match this misc var)
         m%F_Waves   =  m%F_Waves   +  m%WAMIT2%F_Waves2

      END IF



      IF ( u%Morison%LumpedMesh%Committed ) THEN  ! Make sure we are using Morison / there is a valid mesh
         CALL Morison_CalcOutput( Time, u%Morison, p%Morison, x%Morison, xd%Morison,  &
                                z%Morison, OtherState%Morison, y%Morison, m%Morison, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )                  
      END IF
      
      
         ! Integrate all the mesh loads onto the WAMIT reference Point (WRP) at (0,0,0)
      m%F_Hydro = CalcLoadsAtWRP( y, u, m%y_mapped, m%AllHdroOrigin_position, m%MrsnLumpedMesh_position, m%MrsnDistribMesh_position, m%HD_MeshMap, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )                  
      
      
      
      
      
          
      
         ! Write the HydroDyn-level output file data if the user requested module-level output
         ! and the current time has advanced since the last stored time step.
         
      IF ( (p%OutSwtch == 1 .OR. p%OutSwtch == 3) .AND. ( Time > m%LastOutTime ) ) THEN               
         CALL HDOut_WriteOutputs( m%LastOutTime, y, p, m%Decimate, ErrStat2, ErrMsg2 )         
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )                  
      END IF
      
      
         ! Map calculated results into the AllOuts Array
      CALL HDOut_MapOutputs( Time, y, p%NWaveElev, WaveElev, WaveElev1, m%F_PtfmAdd, m%F_Waves, m%F_Hydro, q, qdot, qdotdot, AllOuts, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )                  
      
      DO I = 1,p%NumOuts
            y%WriteOutput(I) = p%OutParam(I)%SignM * AllOuts( p%OutParam(I)%Indx )
      END DO    
      
         ! Aggregate the sub-module outputs 
         
      IF ( p%OutSwtch > 0) THEN
         
         J = p%NumOuts + 1        
         
         IF (ALLOCATED( p%WAMIT%OutParam ) .AND. p%WAMIT%NumOuts > 0) THEN
            DO I=1, p%WAMIT%NumOuts
               y%WriteOutput(J) = y%WAMIT%WriteOutput(I)
               J = J + 1
            END DO
         END IF
         
         IF (ALLOCATED( p%Waves2%OutParam ) .AND. p%Waves2%NumOuts > 0) THEN
            DO I=1, p%Waves2%NumOuts
               y%WriteOutput(J) = y%Waves2%WriteOutput(I)
               J = J + 1
            END DO
         END IF

         IF (ALLOCATED( p%WAMIT2%OutParam ) .AND. p%WAMIT2%NumOuts > 0) THEN
            DO I=1, p%WAMIT2%NumOuts
               y%WriteOutput(J) = y%WAMIT2%WriteOutput(I)
               J = J + 1
            END DO
         END IF

         IF (ALLOCATED( p%Morison%OutParam ) .AND. p%Morison%NumOuts > 0) THEN
            DO I=1, p%Morison%NumOuts
               y%WriteOutput(J) = y%Morison%WriteOutput(I)
               J = J + 1
            END DO
         END IF
         
      END IF
      
      m%LastOutTime   = Time
      
END SUBROUTINE HydroDyn_CalcOutput


!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for computing derivatives of continuous states
SUBROUTINE HydroDyn_CalcContStateDeriv( Time, u, p, x, xd, z, OtherState, m, dxdt, ErrStat, ErrMsg )  
   
      REAL(DbKi),                         INTENT(IN   )  :: Time        !< Current simulation time in seconds
      TYPE(HydroDyn_InputType),           INTENT(INOUT)  :: u           !< Inputs at Time (intent OUT only because we're copying the input mesh)
      TYPE(HydroDyn_ParameterType),       INTENT(IN   )  :: p           !< Parameters                             
      TYPE(HydroDyn_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at Time
      TYPE(HydroDyn_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at Time
      TYPE(HydroDyn_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at Time
      TYPE(HydroDyn_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states                    
      TYPE(HydroDyn_MiscVarType),         INTENT(INOUT)  :: m           !< Initial misc/optimization variables           
      TYPE(HydroDyn_ContinuousStateType), INTENT(  OUT)  :: dxdt        !< Continuous state derivatives at Time
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat     !< Error status of the operation     
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      CHARACTER(*), PARAMETER    :: RoutineName = 'HydroDyn_CalcContStateDeriv'
               
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      
         ! Compute the first time derivatives of the continuous states here:
      
   IF ( m%u_WAMIT%Mesh%Committed ) THEN  ! Make sure we are using WAMIT / there is a valid mesh
         
         ! Copy the inputs from the HD mesh into the WAMIT mesh
      CALL MeshCopy( u%Mesh, m%u_WAMIT%Mesh, MESH_UPDATECOPY, ErrStat, ErrMsg )   
         IF ( ErrStat >= AbortErrLev ) RETURN
      
      CALL WAMIT_CalcContStateDeriv( Time, m%u_WAMIT, p%WAMIT, x%WAMIT, xd%WAMIT, z%WAMIT, OtherState%WAMIT, m%WAMIT, dxdt%WAMIT, ErrStat, ErrMsg ) 

   END IF
   
END SUBROUTINE HydroDyn_CalcContStateDeriv


!----------------------------------------------------------------------------------------------------------------------------------
! Tight coupling routine for updating discrete states. Note that the WAMIT_UpdateDiscState violates the framework by having OtherStates
! be intent in/out. If/when this is fixed we can uncomment this routine.
!SUBROUTINE HydroDyn_UpdateDiscState( Time, n, u, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )   
!   
!   REAL(DbKi),                         INTENT(IN   )  :: Time        !< Current simulation time in seconds   
!   INTEGER(IntKi),                     INTENT(IN   )  :: n           !< Current step of the simulation: t = n*Interval
!   TYPE(HydroDyn_InputType),           INTENT(IN   )  :: u           !< Inputs at Time                       
!   TYPE(HydroDyn_ParameterType),       INTENT(IN   )  :: p           !< Parameters                                 
!   TYPE(HydroDyn_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at Time
!   TYPE(HydroDyn_DiscreteStateType),   INTENT(INOUT)  :: xd          !< Input: Discrete states at Time; 
!                                                                     !!   Output: Discrete states at Time + Interval
!   TYPE(HydroDyn_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at Time
!   TYPE(HydroDyn_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other/optimization states           
!   TYPE(HydroDyn_MiscVarType),         INTENT(INOUT)  :: m           !< Initial misc/optimization variables           
!   INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat     !< Error status of the operation
!   CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
!         
!   
!      ! Initialize ErrStat
!      
!   ErrStat = ErrID_None         
!   ErrMsg  = ""               
!      
!      
!         ! Update discrete states 
!         
!   IF ( m%u_WAMIT%Mesh%Committed ) THEN  ! Make sure we are using WAMIT / there is a valid mesh
!         
!         ! Copy the inputs from the HD mesh into the WAMIT mesh
!      CALL MeshCopy( u%Mesh, m%u_WAMIT%Mesh, MESH_UPDATECOPY, ErrStat, ErrMsg )   
!         IF ( ErrStat >= AbortErrLev ) RETURN
!      
!     CALL WAMIT_UpdateDiscState( Time, n, m%u_WAMIT, p%WAMIT, x%WAMIT, xd%WAMIT, z%WAMIT, OtherState%WAMIT, m%WAMIT, ErrStat, ErrMsg )          
!         IF ( ErrStat >= AbortErrLev ) RETURN
!  END IF
!
!END SUBROUTINE HydroDyn_UpdateDiscState


!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for solving for the residual of the constraint state equations
SUBROUTINE HydroDyn_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, m, z_residual, ErrStat, ErrMsg )   
   
   REAL(DbKi),                         INTENT(IN   )  :: Time        !< Current simulation time in seconds   
   TYPE(HydroDyn_InputType),           INTENT(INOUT)  :: u           !< Inputs at Time (intent OUT only because we're copying the input mesh)              
   TYPE(HydroDyn_ParameterType),       INTENT(IN   )  :: p           !< Parameters                           
   TYPE(HydroDyn_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at Time
   TYPE(HydroDyn_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at Time
   TYPE(HydroDyn_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at Time (possibly a guess)
   TYPE(HydroDyn_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other/optimization states                    
   TYPE(HydroDyn_MiscVarType),         INTENT(INOUT)  :: m           !< Initial misc/optimization variables           
   TYPE(HydroDyn_ConstraintStateType), INTENT(  OUT)  :: z_residual  !< Residual of the constraint state equations using  
                                                                     !!     the input values described above      
   INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

               
      ! Initialize ErrStat
         
   ErrStat = ErrID_None         
   ErrMsg  = ""               
      
      
         ! Solve for the constraint states here:
      
   IF ( m%u_WAMIT%Mesh%Committed ) THEN  ! Make sure we are using WAMIT / there is a valid mesh
         
         ! Copy the inputs from the HD mesh into the WAMIT mesh
      CALL MeshCopy( u%Mesh, m%u_WAMIT%Mesh, MESH_UPDATECOPY, ErrStat, ErrMsg )   
         IF ( ErrStat >= AbortErrLev ) RETURN
      
      call WAMIT_CalcConstrStateResidual( Time, m%u_WAMIT, p%WAMIT, x%WAMIT, xd%WAMIT, z%WAMIT, OtherState%WAMIT, m%WAMIT, z_residual%WAMIT, ErrStat, ErrMsg )

   END IF


END SUBROUTINE HydroDyn_CalcConstrStateResidual


!----------------------------------------------------------------------------------------------------------------------------------
FUNCTION CalcLoadsAtWRP( y, u, y_mapped, AllHdroOrigin_position, MrsnLumpedMesh_Postion, MrsnDistribMesh_Position, MeshMapData, ErrStat, ErrMsg )

   TYPE(HydroDyn_OutputType),  INTENT(INOUT)  :: y                   ! Hydrodyn outputs
   TYPE(HydroDyn_InputType),   INTENT(IN   )  :: u                   ! Hydrodyn inputs
   TYPE(MeshType),             INTENT(INOUT)  :: y_mapped            ! This is the mesh which data is mapped onto.  We pass it in to avoid allocating it at each call
   TYPE(MeshType),             INTENT(IN   )  :: AllHdroOrigin_position            ! This is the mesh which data is mapped onto.  We pass it in to avoid allocating it at each call
   TYPE(MeshType),             INTENT(IN   )  :: MrsnLumpedMesh_Postion            ! This is the mesh which data is mapped onto.  We pass it in to avoid allocating it at each call 
   TYPE(MeshType),             INTENT(IN   )  :: MrsnDistribMesh_Position            ! This is the mesh which data is mapped onto.  We pass it in to avoid allocating it at each call
   TYPE(HD_ModuleMapType),     INTENT(INOUT)  :: MeshMapData         ! Map  data structures 
   INTEGER(IntKi),             INTENT(  OUT)  :: ErrStat             ! Error status of the operation
   CHARACTER(*),               INTENT(  OUT)  :: ErrMsg              ! Error message if ErrStat /= ErrID_None                                                         
   REAL(ReKi)                                 :: CalcLoadsAtWRP(6)

      ! local variables
   INTEGER(IntKi)                                 :: ErrStat2                  ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                           :: ErrMsg2                   ! temporary Error message if ErrStat /= ErrID_None
   
   y%AllHdroOrigin%Force = 0.0
   y%AllHdroOrigin%Moment= 0.0
   
   IF ( y%Mesh%Committed  ) THEN

      ! Just transfer the loads because the meshes are at the same location (0,0,0)

      y%AllHdroOrigin%Force  =  y%Mesh%Force
      y%AllHdroOrigin%Moment =  y%Mesh%Moment

   END IF      
      
   IF ( y%Morison%LumpedMesh%Committed ) THEN 

         ! This is viscous drag associate with the WAMIT body and/or filled/flooded forces of the WAMIT body

      CALL Transfer_Point_to_Point( y%Morison%LumpedMesh, y_mapped, MeshMapData%M_P_2_WRP_P, ErrStat2, ErrMsg2, MrsnLumpedMesh_Postion, AllHdroOrigin_position )
         CALL CheckError( ErrStat2, ErrMsg2 )
         IF (ErrStat >= AbortErrLev) RETURN
            
      y%AllHdroOrigin%Force  = y%AllHdroOrigin%Force  + y_mapped%Force
      y%AllHdroOrigin%Moment = y%AllHdroOrigin%Moment + y_mapped%Moment

   END IF
   
   IF ( y%Morison%DistribMesh%Committed ) THEN 

      CALL Transfer_Line2_to_Point( y%Morison%DistribMesh, y_mapped, MeshMapData%M_L_2_WRP_P, ErrStat2, ErrMsg2,  MrsnDistribMesh_Position, AllHdroOrigin_position )
         CALL CheckError( ErrStat2, ErrMsg2 )
         IF (ErrStat >= AbortErrLev) RETURN
 
      y%AllHdroOrigin%Force  = y%AllHdroOrigin%Force  + y_mapped%Force
      y%AllHdroOrigin%Moment = y%AllHdroOrigin%Moment + y_mapped%Moment
         
   END IF
   
   CalcLoadsAtWRP(1:3) = y%AllHdroOrigin%Force(:,1)
   CalcLoadsAtWRP(4:6) = y%AllHdroOrigin%Moment(:,1)

CONTAINS   
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)

      INTEGER(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
      CHARACTER(1024)            :: ErrMsg3     ! The error message (ErrMsg)

      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................
      
      IF ( ErrID /= ErrID_None ) THEN

         IF ( LEN_TRIM(ErrMsg) > 0 ) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//' CalcLoadsAtWRP:'//TRIM(Msg)
         ErrStat = MAX(ErrStat, ErrID)
         
         !.........................................................................................................................
         ! Clean up if we're going to return on error: close files, deallocate local arrays
         !.........................................................................................................................
      END IF

   END SUBROUTINE CheckError
END FUNCTION CalcLoadsAtWRP
   
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE FHA_Force(Time, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg,F_FHA,q,qdot,qdotdot,w)   
USE NWTC_Library

    TYPE(HydroDyn_InputType),           INTENT(IN   )  :: u           !< Inputs at Time (note that this is intent out because we're copying the u%mesh into m%u_wamit%mesh)
    TYPE(HydroDyn_ParameterType),       INTENT(IN   )  :: p           !< Parameters
    TYPE(HydroDyn_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at Time
    TYPE(HydroDyn_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at Time
    TYPE(HydroDyn_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at Time
    TYPE(HydroDyn_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states at Time
    TYPE(HydroDyn_OutputType),          INTENT(IN   )  :: y           !< Outputs computed at Time (Input only so that mesh con-
    TYPE(HydroDyn_MiscVarType),         INTENT(INOUT)  :: m           !< Initial misc/optimization variables

    REAL(SiKi),                         INTENT(IN)     :: w           !< Disturbance (wave) input for controlling TMD frequency

    CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg      !! Error message if ErrStat /= ErrID_None

    INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat     !! Error status of the operation
    INTEGER                                            :: I,II,III        ! Loop counters
    INTEGER,SAVE                                       :: n=1           !! time step counter
    INTEGER(IntKi),SAVE                                :: n_FHA        !! number of FHAs in system
    INTEGER(IntKi),SAVE                                :: run_flag(1) ! Flag to set to either use the FHA routine or ignore it.  Routine will run if any of the masses in input file are greater than 0

    REAL(DbKi),                         INTENT(IN   )  :: Time        !< Current simulation time in seconds
    REAL(DbKi),SAVE                                    :: Time_last        !< Last step simulation time in seconds
    REAL(ReKi)                                         :: dt        !< time step
    REAL(DbKi),                         INTENT(  OUT)  :: F_FHA(6)       !! FHA force vector about the WAMIT ref pt.
    REAL(ReKi),                         INTENT(IN   )  :: q(6), qdot(6), qdotdot(6)  !! Platform Positions, Velocities, and Accelerations at WAMIT Ref Pt (m,rad)     
    REAL(DbKi)                                         :: Un(2,1) ! Results at time step n
    REAL(DbKi)                                         :: Up1(2,1) ! Results at time step n+1
    REAL(DbKi)                                         :: Um1(2,1) ! Results at time step n-1
    REAL(DbKi)                                         :: Uexp(2,1) ! Expected results at time step n+1  
    REAL(ReKi)                                         :: u_dummy(3)           !! dummy place holder for cross products               
    REAL(ReKi)                                         :: q_LF(3,1)   !! Platform positions at FHA in local ref frame
    REAL(ReKi)                                         :: qdot_LF(3,1)   !! Platform velocities at FHA in local ref frame
    REAL(ReKi)                                         :: qdotdot_LF(3,1)   !! Platform acceleration at FHA in local ref frame
    REAL(ReKi)                                         :: q_IF(3,1)   !! Platform positions at FHA in inertia frame
    REAL(ReKi)                                         :: qdot_IF(3,1)   !! Platform  velocities at FHA in inertia frame
    REAL(ReKi)                                         :: qdotdot_IF(3,1)   !! Platform  acceleration at FHA in inertia frame
    REAL(ReKi)                                         :: dx_LF(1)   !! Platform position at FHA in inertia frame along the orientation of the FHA's motion
    REAL(ReKi)                                         :: dxdot_LF(1)   !! Platform velocity at FHA in inertia frame along the orientation of the FHA's motion 
    REAL(ReKi)                                         :: dxdotdot_LF(1)   !! Platform acceleration at FHA in inertia frame along the orientation of the FHA's motion     
    REAL(ReKi)                                         :: T_IF_LF(3,3)    ! Eular angle transformation matrix to go from inertia frame to local frame of FHA
    REAL(ReKi)                                         :: T_LF_IF(3,3)    !  Eular angle transformation matrix to go from  local frame of FHA to inertia frame
    REAL(ReKi)                                         :: T_dq(3,3)    !  Eular angle transformation matrix to go from  fixed inertia frame to frame with system rotational displacements (w.r.t. the inertia frame)
    REAL(ReKi),SAVE                                    :: g(3,1) !=(\0;0;-9.807\)        ! Gravity vector in the local ref frame   
    REAL(ReKi),SAVE,ALLOCATABLE                        :: q_n(:,:), qdot_n(:,:), qdotdot_n(:,:)  !! Platform Positions, Velocities, and Accelerations at WAMIT Ref Pt (m,rad) from last time step      
    REAL(ReKi),SAVE,ALLOCATABLE                        :: u_FHA(:,:)           !! FHA  displacement at time n (m) 
    REAL(ReKi),SAVE,ALLOCATABLE                        :: up_FHA(:,:)           !! FHA  velocity at time n (m/s) 
    REAL(ReKi),SAVE,ALLOCATABLE                        :: upp_FHA(:,:)           !! FHA  acceleration at time n (m/s2)        
    REAL(ReKi),SAVE,ALLOCATABLE                        :: data_hold(:,:)           !! Output data     
    REAL(ReKi),SAVE,ALLOCATABLE                        :: mass(:) !             !! FHA mass (kg)
    REAL(ReKi),SAVE,ALLOCATABLE                        :: k(:)                  !! FHA stiffness (N/m)
    REAL(ReKi),SAVE,ALLOCATABLE                        :: c(:)                  !! FHA damping (N/(m/s))
    REAL(ReKi),SAVE,ALLOCATABLE                        :: cq(:)                  !! FHA nonlinear damping (N/(m/s)^2)
    REAL(ReKi),SAVE,ALLOCATABLE                        :: m_xyz(:,:)            !! position of masses wrt WAMIT output pt. (m)
    REAL(ReKi),SAVE,ALLOCATABLE                        :: u_bar(:,:)            !! Euler angles defining orientation of FHA w.r.t. inertia frame
    REAL(DbKi),SAVE,ALLOCATABLE                        :: F_LF(:,:,:)         !< Force rxn at FHA in local frame
    REAL(DbKi),SAVE,ALLOCATABLE                        :: F_IF(:,:,:)         !< Force rxn at FHA in inertia frame  
    REAL(DbKi),SAVE,ALLOCATABLE                        :: Time_n(:)         !< Array of time stamps over the simulation

    ! Sea State Est Vars
    REAL(DbKi),  DIMENSION(:), ALLOCATABLE, SAVE      :: control_per
    REAL(DbKi),  DIMENSION(:), ALLOCATABLE, SAVE      :: w_buff

    LOGICAL                                           :: reset, out_is_open
    INTEGER(IntKi)                                    :: status, n_buff, inst
    REAL(DbKi)                                        :: I0, minValue, maxValue, kp, ki, freq, damp, t_buff, Hs_est
    REAL(DbKi)                                        :: gamma, kappa, f0, dw, w_fll, w_hat, T_est, wn_control
    REAL(DbKi), SAVE                                  :: w_next
   !  CHARACTER(1024)                                   :: OL_InputFileName, control_filename
    REAL(DbKi), DIMENSION(2,2)                        :: A_fll, C_fll
    REAL(DbKi), DIMENSION(2,1)                        :: B_fll, dx_fll, x_fll, y_fll
    REAL(DbKi), DIMENSION(2,1), SAVE                  :: x_next
    REAL(DbKi), DIMENSION(1,1)                        :: u_fll, dw_temp

    ! TMD Control Vars
    REAL(DbKi),  DIMENSION(:,:), ALLOCATABLE, SAVE    :: Channels
    REAL(DbKi),  DIMENSION(:), ALLOCATABLE, SAVE      :: control_freq
    REAL(DbKi), SAVE                                  :: zeta, c_control, k_control  ! damping ratio, controlled damping, spring const.


    INTEGER                                            :: nvari        ! number of output variables
    CHARACTER(9000)                                    :: out_fmt
    CHARACTER(*), PARAMETER ::                                        FMT_OUT='aaaaaaa' 
    !CHARACTER(*), PARAMETER :: FMT_OUT_HEADER 
    !CHARACTER(*), PARAMETER :: blank
    !REAL(ReKi),ALLOCATABLE                         :: u0(:)           !! FHA  displacement at time 0 (m) 
    !REAL(ReKi),ALLOCATABLE                         :: up0(:)           !! FHA  velocity at time 0 (m/s)  
    
    ! Simulation time incriment  
    dt=p%DT
    !Update time stamp counter if we are in the next time stamp
    IF (Time/=Time_last) THEN
         n=n+1
    END IF    

    ! If first time step read inputs and set ICs for masses      
    IF (n == 1) THEN
        ! Read input file
        CALL Read_FHA_Input(run_flag,n_FHA,u_FHA,up_FHA,upp_FHA,Time_n,mass,k,c,cq,m_xyz,u_bar,q_n,qdot_n,qdotdot_n,F_LF,F_IF,m,p,data_hold)

        ! Read Control Input
         call read_input_ts(m%TMDControlFile,Channels,control_per)
         control_freq        = Channels(:,1)

         inquire(unit = 2000, opened=out_is_open)
         ! print *, TRIM(m%OutRootName)//"_TMD_ControlOut.txt"
         if (.NOT. out_is_open) THEN
            ! print *, TRIM(m%OutRootName)//"_TMD_ControlOut.txt"
            OPEN(unit = 2000,file = TRIM(m%OutRootName)//"_TMD_ControlOut.txt",STATUS='REPLACE')
            !                                         time , w, w_hat, T_est, Hs_est, wn_control, k_control, c_control
            WRITE(2000,'(99(a10,TR4:))') 'Time (s)', 'Wave (m)', 'Tp Est.', 'Hs Est.', 'w_n cntrl','stiffness','damping'
            ! WRITE(2000,*) 'Time (s)', 'Wave (m)', 'Tp Est.', 'Hs Est.', 'w_n cntrl','stiffness','damping'
            WRITE(2000,'(99(a10,TR4:))') '(s)', '(m)', '(s)', '(m)', '(rad/s)','(N/m)','(N-s/m)'
            ! WRITE(2000,*) '(s)', '(m)', '(s)', '(m)', '(rad/s)','(N/m)','(N-s/m)'
         end if

         ! Set initial params based on FHA Input
         zeta = c(1)/ (2 * sqrt(k(1) * mass(1)))
         c_control = c(1)
         k_control = k(1)

         
   
         ! Set buffer params
         t_buff = 200 ! sec
         n_buff = floor(t_buff/DT)
         call init_buffer(w_buff,n_buff)
         w_buff(:) = 0.0
        
        ! Inital positions and velocities along the vector which the FHA acts in the FHA local frame
        DO I=1,n_FHA  
            CALL Euler_Transformation(T_IF_LF,u_bar(1,I),u_bar(2,I),u_bar(3,I)) 
          !u_FHA(n,I) = u0(I)
          !up_FHA(n,I) =up0(I)
          
            u_dummy =  MATMUL(T_IF_LF,q(1:3)+CROSS_PRODUCT(q(4:6),m_xyz(1:3,I))) 
            u_FHA(n,I) = u_dummy(1)
            
            u_dummy =   MATMUL(T_IF_LF,qdot(1:3)+CROSS_PRODUCT(qdot(4:6),m_xyz(1:3,I)))
            up_FHA(n,I) = u_dummy(1)
            
            u_dummy =  MATMUL(T_IF_LF,qdotdot(1:3)+CROSS_PRODUCT(qdotdot(4:6),m_xyz(1:3,I))) 
        END DO
        
        ! System position, velocity, and acceleration at the WAMIT ref pt in the inertia frame
        DO III=1,6
            q_n(n,III)=q(III)
            qdot_n(n,III)=qdot(III)
        END DO
        
        ! Set gravity vector
        g(1,1)=0
        g(2,1)=0
        g(3,1)=-9.807
    END IF

    ! Sea State Estimator
    ! Init. vars
    if (n==1) THEN         

         ! Init State Space Vars
         x_fll          = reshape((/0,0/),(/2,1/))
         x_next         = reshape((/0,0/),(/2,1/))
         dx_fll         = reshape((/0,0/),(/2,1/))
         y_fll          = reshape((/0,0/),(/2,1/))

         w_fll    = 1
         w_next   = 1

         

      else
         reset = .FALSE.
         status  = 1
      end if

      ! fll params
    f0  = 1
    gamma = 0.16
    kappa = 1.4 
    freq     = 2 * 3.1415 / 300
    damp     = .707

   IF (Time/=Time_last) THEN

      ! Initialize Filter
      IF (Time < 2*DT) THEN
         reset    = .TRUE.
         I0       = w
         status   = 0
         
      END IF
      ! print *, time
       inst = 1
        ! Discrete State Space Rep.
        
        ! state update 
        w_fll        = max(w_next,1e-3)  ! force to be > 0
        x_fll       = x_next

        A_fll       = reshape((/-kappa*w_fll, 1d0, -w_fll*w_fll, 0d0/),(/2,2/))
        B_fll       = reshape((/kappa * w_fll, 0d0/),(/2,1/))
        C_fll       = reshape((/1d0,0d0,0d0,w_fll/),(/2,2/))

        u_fll       = reshape((/w/),(/1,1/))        
        dx_fll      = matmul(A_fll,x_fll) + matmul(B_fll,u_fll)
        y_fll       = matmul(C_fll,x_fll)
        dw         = -gamma * (w - x_fll(1,1)) * y_fll(2,1)

        ! Integrate: forward euler for now
        if (Time<m%TMax) THEN


            x_next      = x_fll + dx_fll * dt
            w_next       = w_fll + dw * dt

        END IF

        ! Filter
        w_hat = SecLPFilter(w_fll, DT, freq, damp, status, reset, inst)

        ! Est Period and saturate
        if (w_hat > 1e-8) THEN
            T_est   = 1.1 * 2 * 3.1415 / w_hat
        ELSE
            T_est   = 20.
        END IF

        IF (T_est > 20) THEN
            T_est   = 20.
        end if

        ! compute control from lookup table interpolation
        wn_control = interp1d(control_per,control_freq,T_est)

        k_control    = wn_control * wn_control * mass(1)
        c_control    = zeta * 2 * sqrt(k_control * mass(1))

        c(:)         = c_control
        k(:)         = k_control


        ! sig height est
        call pp_buffer(w_buff,w)
        Hs_est = 4 * std(w_buff)


        ! write output
      !   WRITE(198,*) w_fll, DT, freq, damp, status, reset, inst
      !   WRITE(199,'(99(F10.6,TR5:))') time , w, w_fll, dw, dt! , wn_control, Hs_est
        WRITE(2000,'(99(ES10.3E2,TR4:))') time , w, T_est, Hs_est, wn_control, k_control, c_control
      !   WRITE(200,*) time , w, T_est, Hs_est, wn_control, k_control, c_control
      !   WRITE(201,*) time , k(1), c(1), zeta, k_control, c_control



   END IF    
    
        ! Store system positions, velocities, and accelerations at WAMIT ref. pt
    DO I=1,6
       q_n(n,I)=q(I)
       qdot_n(n,I)=qdot(I)
       qdotdot_n(n,I)=qdotdot(I)
       
       ! Guess at q(n+1)
       q_n(n+1,I)=q_n(n,I)+qdot_n(n,I)*p%DT
       qdot_n(n+1,I)=qdot_n(n,I)+qdotdot_n(n,I)*p%DT
       
       ! Time stamp
       Time_n(n)=Time   
    END DO 
      

    IF (run_flag(1)==1) THEN
         DO I=1, n_FHA          
            IF (n<3) THEN ! We are using 1st order backward Euler
                Um1(1,1)=0
                Um1(2,1)=0
            ELSE ! We are using 2nd order backward Euler
                Um1(1,1)=u_FHA(n-1,I)
                Um1(2,1)=up_FHA(n-1,I)
            END IF

           ! Last incriment
            Un(1,1)=u_FHA(n,I)
            Un(2,1)=up_FHA(n,I)

            ! Inital guess at time n+1
            Up1(1,1)=Un(1,1)
            Up1(2,1)=Un(2,1)
            
            ! Calculate hull displacement and velocity at FHA location in FHA local ref frame
!             CALL Euler_Transformation(T_IF_LF,ATAN2(u_bar(3,I),u_bar(2,I)),ATAN2(u_bar(3,I),u_bar(1,I)),ATAN2(u_bar(2,I),u_bar(1,I)))
            CALL Euler_Transformation(T_IF_LF,u_bar(1,I),u_bar(2,I),u_bar(3,I)) 
            ! Position of the hull at the FHA location along the vector which the FHA acts in the FHA local frame
            q_IF(1:3,1)=q(1:3)+CROSS_PRODUCT(q(4:6),m_xyz(1:3,I))
            q_LF(1:3,1)=MATMUL(T_IF_LF,q_IF(1:3,1))
            dx_LF(1)=q_LF(1,1)
            
            ! Velocity of the hull at the FHA location along the vector which the FHA acts in the FHA local frame
            qdot_IF(1:3,1)=qdot(1:3)+CROSS_PRODUCT(qdot(4:6),m_xyz(1:3,I))
            qdot_LF(1:3,1)=MATMUL(T_IF_LF,qdot_IF(1:3,1))
            dxdot_LF(1)=qdot_LF(1,1)
            
            ! Reduce 
            CALL iterate_Up1(Un,Up1,Um1,k(I),c(I),cq(I),mass(I),dx_LF,dxdot_LF,m_xyz(1:3,I),dt,n) 
            
            ! Update position and velocity at n+1
            u_FHA(n+1,I)=Up1(1,1)
            up_FHA(n+1,I)=Up1(2,1)
        END DO
   END IF               

    ! Resolve forces to WAMIT ref point
    F_FHA(1)=0
    F_FHA(2)=0
    F_FHA(3)=0 
    F_FHA(4)=0
    F_FHA(5)=0
    F_FHA(6)=0     
    IF (run_flag(1)==1) THEN 
      
        DO I=1,n_FHA        
            ! Eular angle transformation matrix to go from inertia frame to local frame of FHA
            CALL Euler_Transformation(T_IF_LF,u_bar(1,I),u_bar(2,I),u_bar(3,I)) 
            ! Eular angle transformation matrix to go from local frame of FHA to inertia frame
            T_LF_IF= TRANSPOSE(T_IF_LF)
            
            ! Forces calculated about the local FHA frame
            ! Calculate external force due to relative displacement btw mass and undisplaced location in the hull 
            q_IF(1:3,1)=q(1:3)+CROSS_PRODUCT(q(4:6),m_xyz(1:3,I))
            q_LF(1:3,1)=MATMUL(T_IF_LF,q_IF(1:3,1))
            dx_LF(1)=q_LF(1,1)
            

            qdot_IF(1:3,1)=qdot(1:3)+CROSS_PRODUCT(qdot(4:6),m_xyz(1:3,I))
            qdot_LF(1:3,1)=MATMUL(T_IF_LF,qdot_IF(1:3,1))
            dxdot_LF(1)=qdot_LF(1,1)
            
            ! Acceleration of the hull at the FHA location along the vector which the FHA acts in the FHA local frame
            qdotdot_IF(1:3,1)=qdotdot(1:3)+CROSS_PRODUCT(qdotdot(4:6),m_xyz(1:3,I))
            qdotdot_LF(1:3,1)=MATMUL(T_IF_LF,qdotdot_IF(1:3,1))
            dxdotdot_LF(1)=qdotdot_LF(1,1)

            !F_LF(n,1,I)=(k(I)*(u_FHA(n,I)-dx_LF(1))+c(I)*(up_FHA(n,I)-dxdot_LF(1)))! Local frame
            F_LF(n,1,I)=k(I)*(u_FHA(n,I)-dx_LF(1))+c(I)*(up_FHA(n,I)-dxdot_LF(1))+cq(I)*(up_FHA(n,I)-dxdot_LF(1))*abs((up_FHA(n,I)-dxdot_LF(1)))! Local frame
            F_LF(n,2,I)=0
            F_LF(n,3,I)=0

            ! Forces due to gravity 
            F_LF(n,1:3,I)=F_LF(n,1:3,I)+MATMUL(T_IF_LF,g(1:3,1)*mass(I))
            
            ! Translational inertia force due to system acceleration along the FHA axes 2 and 3 (axis 1 is along the motion of the mass and is independent from the system)
            qdotdot_IF(1:3,1)=qdotdot(1:3)+CROSS_PRODUCT(qdotdot(4:6),m_xyz(1:3,I))
            qdotdot_LF(1:3,1)=MATMUL(T_IF_LF,qdotdot_IF(1:3,1))
            F_LF(n,2,I)=F_LF(n,2,I)+-qdotdot_LF(2,1)*mass(I)
            F_LF(n,3,I)=F_LF(n,3,I)+-qdotdot_LF(3,1)*mass(I)

   
            
            ! Forces calculated about the WAMIT ref point in the inertia frame
            F_IF(n,1:3,I)=MATMUL(T_LF_IF,F_LF(n,1:3,I))
            F_IF(n,4:6,I)=CROSS_PRODUCT(m_xyz(1:3,I),F_IF(n,1:3,I))

            ! Force at the WAMIT pt due to rotational stiffness
            F_IF(n,4,I)=F_IF(n,4,I)+-mass(I)*-9.807*m_xyz(3,I)*q(4)
            F_IF(n,5,I)=F_IF(n,5,I)+-mass(I)*-9.807*m_xyz(3,I)*q(5)  

            ! Store output data
            data_hold(n,(I-1)*3+1)=F_LF(n,1,I)
            data_hold(n,(I-1)*3+2)=u_FHA(n,I)
            data_hold(n,(I-1)*3+3)=up_FHA(n,I)

        END DO

        ! Sum all forces and moments due to all FHAs about the WAMIT ref pt
        DO II=1, 6
            DO I=1, n_FHA
                F_FHA(II)=F_FHA(II)+F_IF(n,II,I)
            END DO
        END DO
    END IF
      


    ! Update Time Last 
     Time_last=Time     
!
   IF (Time==m%TMax) THEN
       nvari=(n_FHA*3)
    
        
       
      OPEN(unit = 1,file = TRIM(m%OutRootName)//"_Hull_TMD_Output.txt",STATUS='REPLACE')
      WRITE(1,'(a10,TR5)',advance="no") 'Time (s)'
      
      DO I=1,n_FHA
          IF (I==n_FHA) THEN
            write(1,'(a10,TR5,a10,TR5,a10,TR5)') 'F (N)','x (m)','xp (m/s)'
          ELSE
            write(1,'(a10,TR5,a10,TR5,a10,TR5)',advance="no") 'F (N)','x (m)','xp (m/s)'
          ENDIF
      END DO
      
      out_fmt = '(F10.3,TR5,ES10.3E2,TR5)'

      out_fmt = '(F10.3,TR5'  !,ES10.3E2,TR5)'
      do I = 1,n_FHA*3
          out_fmt = trim(out_fmt) // ',ES10.3E2,TR5'
      end do
      out_fmt = trim(out_fmt) // ')'

      DO I=1,n-1
         write(1,out_fmt) Time_n(I), data_hold(I,:)
      END DO
      CLOSE (1)
      CLOSE (2000)
   END IF
   
   !

      
END SUBROUTINE FHA_Force
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Euler_Transformation(T,Y4,Y5,Y6)

    REAL(ReKi),INTENT(  OUT)                         :: T(3,3)
    REAL(ReKi),INTENT(IN   )                         :: Y4             
    REAL(ReKi),INTENT(IN   )                         :: Y5     
    REAL(ReKi),INTENT(IN   )                         :: Y6     

    T(1,1)=cos(Y5)*cos(Y6)
    T(1,2)=-cos(Y5)*sin(Y6)
    T(1,3)=sin(Y5)
    T(2,1)=(cos(Y4)*sin(Y6))+(cos(Y6)*sin(Y4)*sin(Y5))
    T(2,2)=(cos(Y4)*cos(Y6))-(sin(Y4)*sin(Y5)*sin(Y6))
    T(2,3)=-cos(Y5)*sin(Y4)
    T(3,1)=(sin(Y4)*sin(Y6))-(cos(Y4)*cos(Y6)*sin(Y5))
    T(3,2)=(cos(Y6)*sin(Y4))+(cos(Y4)*sin(Y5)*sin(Y6))
    T(3,3)=cos(Y4)*cos(Y5)

END SUBROUTINE Euler_Transformation
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Read_FHA_Input(run_flag,n_FHA,u_FHA,up_FHA,upp_FHA,Time_n,mass,k,c,cq,m_xyz,u_bar,q_n,qdot_n,qdotdot_n,F_LF,F_IF,m,p,data_hold)
    
    INTEGER(IntKi),INTENT(  OUT)                                 :: run_flag(1) ! Flag to set to either use the FHA routine or ignore it.  Routine will run if any of the masses in input file are greater than 0
    INTEGER(IntKi),INTENT(  OUT)                                 :: n_FHA        !! number of FHAs in system
    INTEGER(IntKi)                                               :: n_max        !! number of time steps over simulation
    INTEGER(IntKi)                                               :: III !Counter
    
    TYPE(HydroDyn_ParameterType),       INTENT(IN   )            :: p           !< Parameters
    TYPE(HydroDyn_MiscVarType),         INTENT(IN   )            :: m           !< Initial misc/optimization variables

    REAL(ReKi),ALLOCATABLE,INTENT(  OUT)                         :: u_FHA(:,:)           !! FHA  displacement at time n (m) 
    REAL(ReKi),ALLOCATABLE,INTENT(  OUT)                         :: up_FHA(:,:)           !! FHA  velocity at time n (m)  
    REAL(ReKi),ALLOCATABLE,INTENT(  OUT)                         :: upp_FHA(:,:)          !! FHA  acceleration at time n (m) 
    
    REAL(ReKi),ALLOCATABLE,INTENT(  OUT)                         :: mass(:)          !! FHA mass (kg)
    REAL(ReKi),ALLOCATABLE,INTENT(  OUT)                         :: k(:)           !! FHA stiffness (N/m)
    REAL(ReKi),ALLOCATABLE,INTENT(  OUT)                         :: c(:)           !! FHA damping (N/(m/s))
    REAL(ReKi),ALLOCATABLE,INTENT(  OUT)                         :: cq(:)           !! FHA nonlinear damping (N/(m/s))
    REAL(ReKi),ALLOCATABLE,INTENT(  OUT)                         :: m_xyz(:,:)            !! position of masses wrt WAMIT output pt. (m)
    REAL(ReKi),ALLOCATABLE,INTENT(  OUT)                         :: u_bar(:,:)            !! Euler angles defining orientation of FHA w.r.t. inertia frame
    REAL(ReKi),ALLOCATABLE,INTENT(  OUT)                         :: q_n(:,:), qdot_n(:,:), qdotdot_n(:,:)  !! Platform Positions, Velocities, and Accelerations at WAMIT Ref Pt (m,rad) from last time step    
    REAL(DbKi),ALLOCATABLE,INTENT(  OUT)                         :: F_LF(:,:,:)         !< Force rxn at FHA in local frame
    REAL(DbKi),ALLOCATABLE,INTENT(  OUT)                         :: F_IF(:,:,:)         !< Force rxn at FHA in inertia frame 
    REAL(DbKi),ALLOCATABLE,INTENT(  OUT)                         :: Time_n(:)           !! FHA stiffness (N/m) 
    REAL(ReKi),ALLOCATABLE,INTENT(  OUT)                         :: data_hold(:,:)          !! FHA mass (kg)
    !REAL(ReKi),ALLOCATABLE,INTENT(  OUT)                         :: u0(:)           !! FHA  displacement at time 0 (m) 
    !REAL(ReKi),ALLOCATABLE,INTENT(  OUT)                         :: up0(:)           !! FHA  velocity at time 0 (m/s)  

    n_max=CEILING(m%TMAX/p%DT)+2
!n_max=CEILING(10000/p%DT)*2
    ! Read input data
    !OPEN(unit = 1,file = "FHA_Input.dat")
    PRINT *, "Reading TMD File:"//trim(m%TMDFile)
    OPEN(unit = 1,file = trim(m%TMDFile))
    
    READ (1,*) 
    READ (1,*) n_FHA
    READ (1,*)
    
    ALLOCATE(u_FHA(n_max,n_FHA))
    ALLOCATE(up_FHA(n_max,n_FHA))
    ALLOCATE(upp_FHA(n_max,n_FHA))
    ALLOCATE(q_n(n_max,6))
    ALLOCATE(qdot_n(n_max,6))
    ALLOCATE(qdotdot_n(n_max,6))  
    ALLOCATE(F_LF(n_max,3,n_FHA))
    ALLOCATE(F_IF(n_max,6,n_FHA))
    ALLOCATE(Time_n(n_max))
    ALLOCATE(mass(n_FHA))
    ALLOCATE(k(n_FHA))
    ALLOCATE(c(n_FHA))
    ALLOCATE(cq(n_FHA))
    ALLOCATE(m_xyz(3,n_FHA))
    ALLOCATE(u_bar(3,n_FHA))
    ALLOCATE(data_hold(n_max,n_FHA*3))
    !
    !ALLOCATE(u0(n_FHA))
    !ALLOCATE(up0(n_FHA))

    DO III=1,n_FHA
        READ (1,*) mass(III),k(III),c(III), cq(III), m_xyz(1,III), m_xyz(2,III), m_xyz(3,III), u_bar(1,III), u_bar(2,III), u_bar(3,III) !, u0(III), up0(III)
    END DO

    CLOSE(1)
    
    ! Run flag, if any masses of FHAs from input file are greater than 0 then routine will run, else F_FHA is set to zeros(6,1)
    IF (sum(mass)>0) THEN
        run_flag(1)=1
    ELSE
        run_flag(1)=0
    END IF
        
        
!        DO III=1,n_FHA
!         WRITE (*,'(E8.3,2X,E8.3,2X,E8.3,2X,F8.3,2X,F8.3,2X,F8.3,2X,F8.3,2X,F8.3,2X,F8.3)') mass(III),k,c, m_xyz(1,III), m_xyz(2,III), m_xyz(3,III), u_bar(1,III), u_bar(2,III), u_bar(3,III)
!         !900 format (E5.2,1X,E5.2,1X,E5.2) !,F5.2,F5.2,F5.2,F5.2,F5.2,F5.2)
!        END DO

END SUBROUTINE Read_FHA_Input
!----------------------------------------------------------------------------------------------------------------------------------      
SUBROUTINE get_G(Up1,Un,Um1,dt,G,k,c,mass,F,n)   
 
      REAL(DbKi),                         INTENT(IN   )  :: Un(2,1) ! Results at time step n
      REAL(DbKi),                         INTENT(IN   )  :: Up1(2,1) ! Results at time step n+1 
      REAL(DbKi),                         INTENT(IN   )  :: Um1(2,1) ! Results at time step n-1 
      INTEGER,                            INTENT(IN   )  :: n        !< counter
      REAL(ReKi),                         INTENT(IN   )  :: dt        !< time incriment      
      REAL(DbKi),                         INTENT(IN   )  :: F(2,1)     !< Force vector    
      REAL(DbKi),                         INTENT(  OUT)  :: G(2,1)          
      REAL(ReKi),                         INTENT(IN   )  :: mass       !! FHA mass (kg)
      REAL(ReKi),                         INTENT(IN   )  :: k          !! FHA stiffness (N/m)
      REAL(ReKi),                         INTENT(IN   )  :: c          !! FHA damping (N/(m/s))
            
    IF (N<3) THEN ! Use 1st order accuracy
        G(1,1)=(Up1(1,1)-Un(1,1))/dt-F(1,1)
        G(2,1)=(Up1(2,1)-Un(2,1))/dt-F(2,1)
    ELSE ! Use 2nd order accuracy
        G(1,1)=(1.5*Up1(1,1)+-2*Un(1,1)+.5*Um1(1,1))/dt-F(1,1)
        G(2,1)=(1.5*Up1(2,1)+-2*Un(2,1)+.5*Um1(2,1))/dt-F(2,1)
    END IF
      
END SUBROUTINE get_G 
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE get_J(Up1,Un,Um1,dt,J,k,c,cq,mass,dx_LF,dxdot_LF,n)   

    REAL(DbKi),                         INTENT(IN   )  :: Un(2,1) ! Results at time step n
    REAL(DbKi),                         INTENT(IN   )  :: Up1(2,1) ! Results at time step n+1 
    REAL(DbKi),                         INTENT(IN   )  :: Um1(2,1) ! Results at time step n-1 
    INTEGER,                           INTENT(IN   )  :: n        !< counter
    REAL(ReKi),                         INTENT(IN   )   :: dt        !< time incriment
    REAL(ReKi),                         INTENT(IN   )  :: mass           !! FHA mass (kg)
    REAL(ReKi),                         INTENT(IN   )  :: k          !! FHA stiffness (N/m)
    REAL(ReKi),                         INTENT(IN   )  :: c            !! FHA damping (N/(m/s))
    REAL(ReKi),                         INTENT(IN   )  :: cq            !! FHA nonlinear damping (N/(m/s)^2)
    REAL(DbKi)                                         :: J(2,2)         !< Jacobian matrix   
    REAL(DbKi)                                         :: G1(2,1) 
    REAL(DbKi)                                         :: G2(2,1)
    REAL(DbKi)                                         :: F1(2,1) 
    REAL(DbKi)                                         :: F2(2,1)
    REAL(DbKi)                                         :: dU(2,1) 
    INTEGER                                            :: I        ! Generic counters
    REAL(ReKi),                         INTENT(IN   ) :: dx_LF(1)   !! Platform Position at FHA in inertia frame along the orientation of the FHA's motion
    REAL(ReKi),                         INTENT(IN   ):: dxdot_LF(1)   !! Platform velocity at FHA in inertia frame along the orientation of the FHA's motion

    DO I=1,2
        dU(1,1)=0
        dU(2,1)=0
        dU(I,1)=.00001

        CALL get_F(Up1,k,c,cq,mass,F1,dx_LF,dxdot_LF)
        CALL get_G(Up1+dU,Un,Um1,dt,G1,k,c,mass,F1,n)
        CALL get_F(Up1+dU,k,c,cq,mass,F2,dx_LF,dxdot_LF)
        CALL get_G(Up1,Un,Um1,dt,G2,k,c,mass,F2,n)

        J(1,I)=(G1(1,1)-G2(1,1))/dU(I,1)
        J(2,I)=(G1(2,1)-G2(2,1))/dU(I,1)
    END DO

      
END SUBROUTINE get_J 
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE get_F(U,k,c,cq,mass,F,dx_LF,dxdot_LF)   
 
    REAL(DbKi),                         INTENT(IN   )  :: U(2,1)       !< Predicted results for time step n      
    REAL(DbKi),                         INTENT(  OUT)  :: F(2,1)         !< Force vector          
    REAL(ReKi),                         INTENT(IN   )  :: mass           !! FHA mass (kg)
    REAL(ReKi),                         INTENT(IN   )  :: k          !! FHA stiffness (N/m)
    REAL(ReKi),                         INTENT(IN   )  :: c            !! FHA damping (N/(m/s))
    REAL(ReKi),                         INTENT(IN   )  :: cq            !! FHA nonlinear damping (N/(m/s)^2)
    REAL(DbKi)                                         :: M(2,2)           !! mass matrix
    REAL(DbKi)                                         :: Minv(2,2)           !! Inverse of mass matrix
    REAL(ReKi),                         INTENT(IN   ) :: dx_LF(1)   !! Platform Position at FHA in inertia frame along the orientation of the FHA's motion
    REAL(ReKi),                         INTENT(IN   ):: dxdot_LF(1)   !! Platform velocity at FHA in inertia frame along the orientation of the FHA's motion
    ! Assemble mass matrix M
    M(1,1)=1;
    M(1,2)=0;
    M(2,1)=0;
    M(2,2)=mass;

    ! Calc Inverse of matrix M
    Minv(1,1)=M(2,2)*(1/(M(1,1)*M(2,2)-M(1,2)*M(2,1)))
    Minv(1,2)=-M(1,2)*(1/(M(1,1)*M(2,2)-M(1,2)*M(2,1)))
    Minv(2,1)=-M(2,1)*(1/(M(1,1)*M(2,2)-M(1,2)*M(2,1)))
    Minv(2,2)=M(1,1)*(1/(M(1,1)*M(2,2)-M(1,2)*M(2,1)))

    F(1,1)=U(2,1)
    !F(2,1)=-(k*(U(1,1)-dx_LF(1))+c*(U(2,1)-dxdot_LF(1)))
 F(2,1)=-(k*(U(1,1)-dx_LF(1))+c*(U(2,1)-dxdot_LF(1))+cq*(U(2,1)-dxdot_LF(1))*abs((U(2,1)-dxdot_LF(1))))
    F(1,1)=(Minv(1,1)*F(1,1))+(Minv(1,2)*F(2,1))
    F(2,1)=(Minv(2,1)*F(1,1))+(Minv(2,2)*F(2,1))
END SUBROUTINE get_F 
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE iterate_Up1(Un,Up1,Um1,k,c,cq,mass,dx_LF,dxdot_LF,m_xyz,dt,n) 

    REAL(ReKi),                         INTENT(IN   )  :: dt
    INTEGER,                            INTENT(IN   )  :: n
    REAL(DbKi),                         INTENT(IN   )  :: Un(2,1)       !< Predicted results for time step n
    REAL(DbKi),                         INTENT(INOUT)  :: Up1(2,1)       !< Predicted results for time step n+1
    REAL(DbKi),                         INTENT(IN   )  :: Um1(2,1)       !< Predicted results for time step n-1 (currently not used)
    REAL(DbKi)                                         :: F(2,1)         !< Force vector          
    REAL(ReKi),                         INTENT(IN   )  :: mass           !! FHA mass (kg)
    REAL(ReKi),                         INTENT(IN   )  :: k          !! FHA stiffness (N/m)
    REAL(ReKi),                         INTENT(IN   )  :: c            !! FHA damping (N/(m/s))
    REAL(ReKi),                         INTENT(IN   )  :: cq            !! FHA nonlinear damping (N/(m/s)^2)
    REAL(ReKi),                         INTENT(IN   )  :: m_xyz(3,1)            !! position of masses wrt WAMIT output pt. (m)
    REAL(DbKi)                                         :: err ! Reduction loop error
    REAL(DbKi)                                         :: toll=1e-8 ! Reduction loop error tollerance for convergence 
    REAL(DbKi)                                         :: max_iter=200 ! Max iterations in reduction loop error tollerance for convergence  
    INTEGER                                            :: iter ! Reduction loop counter     
    REAL(ReKi)                                         :: u_dummy(3)           !! dummy place holder for cross products      
    REAL(DbKi)                                         :: G(2,1) ! Results at time step n
    REAL(DbKi)                                         :: G_new(2,1) ! Results at time step n+1
    REAL(DbKi)                                         :: J(2,2) !Jacobian matrix
    REAL(DbKi)                                         :: Jinv(2,2) ! Inverse of the Jacobian matrix
    REAL(ReKi),                         INTENT(IN   )  :: dx_LF(1)   !! Platform Position at FHA in inertia frame along the orientation of the FHA's motion
    REAL(ReKi),                         INTENT(IN   )  :: dxdot_LF(1)   !! Platform velocity at FHA in inertia frame along the orientation of the FHA's motion




    err=toll*2
    iter=1
    DO WHILE(err>toll.AND.iter<max_iter) ! Loop until convergence                            
        ! Call force vector
        CALL get_F(Up1,k,c,cq,mass,F,dx_LF,dxdot_LF)

        ! Get G vector
        CALL get_G(Up1,Un,Um1,dt,G,k,c,mass,F,n) 

        ! Call the Jacobian
        CALL get_J(Up1,Un,Um1,dt,J,k,c,cq,mass,dx_LF,dxdot_LF,n)  
        ! Invert Jacobian so we can mult.
        Jinv(1,1)=J(2,2)*(1/(J(1,1)*J(2,2)-J(1,2)*J(2,1)))
        Jinv(1,2)=-J(1,2)*(1/(J(1,1)*J(2,2)-J(1,2)*J(2,1)))
        Jinv(2,1)=-J(2,1)*(1/(J(1,1)*J(2,2)-J(1,2)*J(2,1)))
        Jinv(2,2)=J(1,1)*(1/(J(1,1)*J(2,2)-J(1,2)*J(2,1)))

        ! Update values of U at n+1
        Up1(1,1)=Up1(1,1)-((Jinv(1,1)*G(1,1))+(Jinv(1,2)*G(2,1)))
        Up1(2,1)=Up1(2,1)-((Jinv(2,1)*G(1,1))+(Jinv(2,2)*G(2,1)))

        ! Call force vector
        CALL get_F(Up1,k,c,cq,mass,F,dx_LF,dxdot_LF)

        ! Get G vector
        CALL get_G(Up1,Un,Um1,dt,G_new,k,c,mass,F,n)

        ! Calculate error and update counter 
        err=sqrt(((G(1,1)-G_new(1,1))*(G(1,1)-G_new(1,1)))+((G(2,1)-G_new(2,1))*(G(2,1)-G_new(2,1))))
        iter=iter+1
    END DO
                               
END SUBROUTINE iterate_Up1   
!----------------------------------------------------------------------------------------------------------------------------------
real function std(arr)
        REAL(DbKi), DIMENSION(:), INTENT(IN)       :: arr
        INTEGER(IntKi)                              :: i
        REAL(DbKi)                                 :: sum, sum_sq, var

        sum = 0.0
        sum_sq = 0.0

        do i = 1,size(arr)
            sum     = sum + arr(i)
            sum_sq  = sum_sq + arr(i) * arr(i)
        end do

        var     = (sum_sq - sum * sum / size(arr))/(size(arr)-1)
        std     = sqrt(var)

    end function std
!----------------------------------------------------------------------------------------------------------------------------------
    subroutine init_buffer(buff,n_buff)
      REAL(DbKi),ALLOCATABLE,INTENT(  OUT)                         :: buff(:)
      INTEGER(IntKi), INTENT(IN)                                   :: n_buff
      ALLOCATE(buff(n_buff))

   end subroutine init_buffer

    !----------------------------------------------------------------------------------------------------------------------------------
    subroutine pp_buffer(buffer,new_sig)
        ! push and pop buffer, new signal goes at end, all other signal indices decrease by 1
        
        REAL(DbKi), DIMENSION(:), INTENT(INOUT)            :: buffer
        REAL(ReKi), INTENT(IN)                             :: new_sig
        INTEGER(IntKi)                                      :: i

        do i=2,size(buffer)
            buffer(i-1) = buffer(i)
        end do
        buffer(size(buffer)) = new_sig



    end subroutine pp_buffer
!----------------------------------------------------------------------------------------------------------------------------------
    REAL FUNCTION SecLPFilter(InputSignal, DT, CornerFreq, Damp, iStatus, reset, inst)
    ! Discrete time Low-Pass Filter of the form:
    !                               Continuous Time Form:   H(s) = CornerFreq^2/(s^2 + 2*CornerFreq*Damp*s + CornerFreq^2)
    !                               Discrete Time From:     H(z) = (b2*z^2 + b1*z + b0) / (a2*z^2 + a1*z + a0)
        REAL(DbKi), INTENT(IN)         :: InputSignal
        REAL(ReKi), INTENT(IN)         :: DT                       ! time step [s]
        REAL(DbKi), INTENT(IN)         :: CornerFreq               ! corner frequency [rad/s]
        REAL(DbKi), INTENT(IN)         :: Damp                     ! Dampening constant
        INTEGER(IntKi), INTENT(IN)      :: iStatus                  ! A status flag set by the simulation as follows: 0 if this is the first call, 1 for all subsequent time steps, -1 if this is the final call at the end of the simulation.
        INTEGER(IntKi), INTENT(INOUT)   :: inst                     ! Instance number. Every instance of this function needs to have an unique instance number to ensure instances don't influence each other.
        LOGICAL(4), INTENT(IN)      :: reset                    ! Reset the filter to the input signal

        ! Local
        REAL(DbKi), DIMENSION(99), SAVE    :: a2                   ! Denominator coefficient 2
        REAL(DbKi), DIMENSION(99), SAVE    :: a1                   ! Denominator coefficient 1
        REAL(DbKi), DIMENSION(99), SAVE    :: a0                   ! Denominator coefficient 0
        REAL(DbKi), DIMENSION(99), SAVE    :: b2                   ! Numerator coefficient 2
        REAL(DbKi), DIMENSION(99), SAVE    :: b1                   ! Numerator coefficient 1
        REAL(DbKi), DIMENSION(99), SAVE    :: b0                   ! Numerator coefficient 0 
        REAL(DbKi), DIMENSION(99), SAVE    :: InputSignalLast1     ! Input signal the last time this filter was called. Supports 99 separate instances.
        REAL(DbKi), DIMENSION(99), SAVE    :: InputSignalLast2     ! Input signal the next to last time this filter was called. Supports 99 separate instances.
        REAL(DbKi), DIMENSION(99), SAVE    :: OutputSignalLast1    ! Output signal the last time this filter was called. Supports 99 separate instances.
        REAL(DbKi), DIMENSION(99), SAVE    :: OutputSignalLast2    ! Output signal the next to last time this filter was called. Supports 99 separate instances.

        ! Initialization
        IF ((iStatus == 0) .OR. reset )  THEN
            OutputSignalLast1(inst)  = InputSignal
            OutputSignalLast2(inst)  = InputSignal
            InputSignalLast1(inst)   = InputSignal
            InputSignalLast2(inst)   = InputSignal
            
            ! Coefficients
            a2(inst) = DT**2.0*CornerFreq**2.0 + 4.0 + 4.0*Damp*CornerFreq*DT
            a1(inst) = 2.0*DT**2.0*CornerFreq**2.0 - 8.0
            a0(inst) = DT**2.0*CornerFreq**2.0 + 4.0 - 4.0*Damp*CornerFreq*DT
            b2(inst) = DT**2.0*CornerFreq**2.0
            b1(inst) = 2.0*DT**2.0*CornerFreq**2.0
            b0(inst) = DT**2.0*CornerFreq**2.0
        ENDIF

        ! Filter
        SecLPFilter = 1.0/a2(inst) * (b2(inst)*InputSignal + b1(inst)*InputSignalLast1(inst) + b0(inst)*InputSignalLast2(inst) & 
        - a1(inst)*OutputSignalLast1(inst) - a0(inst)*OutputSignalLast2(inst))

        ! SecLPFilter = 1/(4+4*DT*Damp*CornerFreq+DT**2*CornerFreq**2) * ( (8-2*DT**2*CornerFreq**2)*OutputSignalLast1(inst) &
        !                 + (-4+4*DT*Damp*CornerFreq-DT**2*CornerFreq**2)*OutputSignalLast2(inst) + (DT**2*CornerFreq**2)*InputSignal &
        !                     + (2*DT**2*CornerFreq**2)*InputSignalLast1(inst) + (DT**2*CornerFreq**2)*InputSignalLast2(inst) )

        ! Save signals for next time step
        InputSignalLast2(inst)   = InputSignalLast1(inst)
        InputSignalLast1(inst)   = InputSignal
        OutputSignalLast2(inst)  = OutputSignalLast1(inst)
        OutputSignalLast1(inst)  = SecLPFilter

        inst = inst + 1

    END FUNCTION SecLPFilter
!----------------------------------------------------------------------------------------------------------------------------------
    REAL FUNCTION PIController(error, kp, ki, minValue, maxValue, DT, I0, reset, inst)
    ! PI controller, with output saturation

        IMPLICIT NONE
        ! Allocate Inputs
        REAL(DbKi), INTENT(IN)         :: error
        REAL(DbKi), INTENT(IN)         :: kp
        REAL(DbKi), INTENT(IN)         :: ki
        REAL(DbKi), INTENT(IN)         :: minValue
        REAL(DbKi), INTENT(IN)         :: maxValue
        REAL(DbKi), INTENT(IN)         :: DT
        INTEGER(IntKi), INTENT(INOUT)   :: inst
        REAL(DbKi), INTENT(IN)         :: I0
        LOGICAL, INTENT(IN)         :: reset     
        ! Allocate local variables
        INTEGER(IntKi)                      :: i                                            ! Counter for making arrays
        REAL(DbKi)                         :: PTerm                                        ! Proportional term
        REAL(DbKi), DIMENSION(99), SAVE    :: ITerm = (/ (real(9999.9), i = 1,99) /)       ! Integral term, current.
        REAL(DbKi), DIMENSION(99), SAVE    :: ITermLast = (/ (real(9999.9), i = 1,99) /)   ! Integral term, the last time this controller was called. Supports 99 separate instances.
        INTEGER(IntKi), DIMENSION(99), SAVE :: FirstCall = (/ (1, i=1,99) /)                ! First call of this function?
        
        ! Initialize persistent variables/arrays, and set inital condition for integrator term
        IF ((FirstCall(inst) == 1) .OR. reset) THEN
            ITerm(inst) = I0
            ITermLast(inst) = I0
            
            FirstCall(inst) = 0
            PIController = I0
        ELSE
            PTerm = kp*error
            ITerm(inst) = ITerm(inst) + DT*ki*error
            ITerm(inst) = saturate(ITerm(inst), minValue, maxValue)
            PIController = saturate(PTerm + ITerm(inst), minValue, maxValue)
        
            ITermLast(inst) = ITerm(inst)
        END IF
        inst = inst + 1
        
    END FUNCTION PIController
!----------------------------------------------------------------------------------------------------------------------------------
    REAL FUNCTION saturate(inputValue, minValue, maxValue)
    ! Saturates inputValue. Makes sure it is not smaller than minValue and not larger than maxValue

        IMPLICIT NONE

        REAL(DbKi), INTENT(IN)     :: inputValue
        REAL(DbKi), INTENT(IN)     :: minValue
        REAL(DbKi), INTENT(IN)     :: maxValue

        saturate = MIN(MAX(inputValue,minValue), maxValue)

    END FUNCTION saturate
!-------------------------------------------------------------------------------------------------------------------------------
    REAL FUNCTION interp1d(xData, yData, xq)
    ! interp1d 1-D interpolation (table lookup), xData should be monotonically increasing

        IMPLICIT NONE
        ! Inputs
        REAL(DbKi), DIMENSION(:), INTENT(IN)       :: xData        ! Provided x data (vector), to be interpolated
        REAL(DbKi), DIMENSION(:), INTENT(IN)       :: yData        ! Provided y data (vector), to be interpolated
        REAL(DbKi), INTENT(IN)                     :: xq           ! x-value for which the y value has to be interpolated
        INTEGER(IntKi)                              :: I            ! Iteration index
        
        ! Interpolate
        IF (xq <= MINVAL(xData)) THEN
            interp1d = yData(1)
        ELSEIF (xq >= MAXVAL(xData)) THEN
            interp1d = yData(SIZE(xData))
        ELSE
            DO I = 1, SIZE(xData)
                IF (xq <= xData(I)) THEN
                    interp1d = yData(I-1) + (yData(I) - yData(I-1))/(xData(I) - xData(I-1))*(xq - xData(I-1))
                    EXIT
                ELSE
                    CONTINUE
                END IF
            END DO
        END IF
        
    END FUNCTION interp1d
!----------------------------------------------------------------------------------------------------------------------------------

    subroutine read_input_ts(OL_InputFileName,Channels,Breakpoints)
        INTEGER(IntKi)                          :: accINFILE_size               ! size of DISCON input filename
        INTEGER(IntKi), PARAMETER               :: UnControllerParameters = 89
        INTEGER(IntKi)                          :: LoggingLevel
        CHARACTER(1024), INTENT(IN)         :: OL_InputFileName    ! DISCON input filename
        INTEGER(IntKi), PARAMETER   :: Unit_OL_Input           = 1009

        LOGICAL                 :: FileExists
        INTEGER                 :: IOS                                                 ! I/O status of OPEN.
        CHARACTER(1024)         :: Line              ! Temp variable for reading whole line from file
        INTEGER(IntKi)                                              :: NumComments
        INTEGER(IntKi)                                              :: NumDataLines
        INTEGER(IntKi), PARAMETER                                   :: NumChannels  = 1      ! Number of open loop channels being defined
        INTEGER(IntKi), PARAMETER                                   :: NumCols      = NumChannels + 1
        REAL(DbKi)                                                 :: TmpData(NumCols)  ! Temp variable for reading all columns from a line
        CHARACTER(15)                :: NumString



        REAL(DbKi),  DIMENSION(:,:), ALLOCATABLE, INTENT(OUT)     :: Channels    ! Rating at  Breakpoints
        REAL(DbKi),  DIMENSION(:), ALLOCATABLE, INTENT(OUT)     :: Breakpoints         ! Time index
        INTEGER(IntKi)                                              :: I,J

        ! OL_InputFileName = TRIM(OL_InputFileName)


        !-------------------------------------------------------------------------------------------------
        ! Read timeseries from input file, borrowed (read: copied) from (Open)FAST team...thanks!
        !-------------------------------------------------------------------------------------------------

        !-------------------------------------------------------------------------------------------------
        ! Open the file for reading
        !-------------------------------------------------------------------------------------------------

        INQUIRE (FILE = OL_InputFileName, EXIST = FileExists)

        IF ( .NOT. FileExists) THEN
            PRINT *, TRIM(OL_InputFileName)// ' does not exist, setting R = 1 for all time'
            ALLOCATE(Channels(2,1))
            ALLOCATE(Breakpoints(2))
            Channels(1,1) = 1;        Channels(2,1) = 1
            Breakpoints(1)      = 0;        Breakpoints(2)       = 90000;

        else
            print *, 'opening file'
            OPEN( Unit_OL_Input, FILE=TRIM(OL_InputFileName), STATUS='OLD', FORM='FORMATTED', IOSTAT=IOS, ACTION='READ' )

            IF (IOS /= 0) THEN
                PRINT *, 'Cannot open ' // TRIM(OL_InputFileName) // ', setting R = 1 for all time'
                ALLOCATE(Channels(2,1))
                ALLOCATE(Breakpoints(2))
                Channels(1,1) = 1;        Channels(2,1) = 1
                Breakpoints(1)      = 0;        Breakpoints(2)       = 90000;
                CLOSE(Unit_OL_Input)
            
            else
                ! Do all the stuff!

                !-------------------------------------------------------------------------------------------------
                ! Find the number of comment lines
                !-------------------------------------------------------------------------------------------------
                LINE = '!'                          ! Initialize the line for the DO WHILE LOOP
                NumComments = -1                    ! the last line we read is not a comment, so we'll initialize this to -1 instead of 0

                DO WHILE ( (INDEX( LINE, '!' ) > 0) .OR. (INDEX( LINE, '#' ) > 0) .OR. (INDEX( LINE, '%' ) > 0) ) ! Lines containing "!" are treated as comment lines
                NumComments = NumComments + 1
                
                READ(Unit_OL_Input,'( A )',IOSTAT=IOS) LINE

                ! NWTC_IO has some error catching here that we'll skip for now
            
                END DO !WHILE

                    !-------------------------------------------------------------------------------------------------
                ! Find the number of data lines
                !-------------------------------------------------------------------------------------------------
                NumDataLines = 0

                READ(LINE,*,IOSTAT=IOS) ( TmpData(I), I=1,NumCols ) ! this line was read when we were figuring out the comment lines; let's make sure it contains

                DO WHILE (IOS == 0)  ! read the rest of the file (until an error occurs)
                NumDataLines = NumDataLines + 1
                
                READ(Unit_OL_Input,*,IOSTAT=IOS) ( TmpData(I), I=1,NumCols )
                
                END DO !WHILE
            
            
                IF (NumDataLines < 1) THEN
                WRITE (NumString,'(I11)')  NumComments
                PRINT *, 'Error: '//TRIM(NumString)//' comment lines were found in the uniform wind file, '// &
                            'but the first data line does not contain the proper format.'
                CLOSE(Unit_OL_Input)
                END IF

                !-------------------------------------------------------------------------------------------------
                ! Allocate arrays for the uniform wind data
                !-------------------------------------------------------------------------------------------------
                ALLOCATE(Channels(NumDataLines,NumChannels))
                ALLOCATE(Breakpoints(NumDataLines))

                !-------------------------------------------------------------------------------------------------
                ! Rewind the file (to the beginning) and skip the comment lines
                !-------------------------------------------------------------------------------------------------

                REWIND( Unit_OL_Input )

                DO I=1,NumComments
                    READ(Unit_OL_Input,'( A )',IOSTAT=IOS) LINE
                END DO !I
            

                !-------------------------------------------------------------------------------------------------
                ! Read the data arrays
                !-------------------------------------------------------------------------------------------------
            
                DO I=1,NumDataLines
                
                    READ(Unit_OL_Input,*,IOSTAT=IOS) ( TmpData(J), J=1,NumCols )

                    if (IOS > 0) THEN
                        CLOSE(Unit_OL_Input)
                    endif

                    Breakpoints(I)          = TmpData(1)
                    Channels(I,:)           = TmpData(2:)
            
                END DO !I
                
                ! Rewind for next time
                REWIND( Unit_OL_Input )
            

            end if
        end if
    end subroutine read_input_ts

END MODULE HydroDyn
!**********************************************************************************************************************************
