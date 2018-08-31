#!/bin/tcsh
#____________________________ ATOMS 
setenv gnrm          1.0E-05
setenv Inflection  4
setenv NATOMS      10000                         # Number of atoms in the problem
setenv type1      Zr             # Still use the name of Si, but already modified codes for ZrCu system
setenv type2      Cu
#_____________________________ ART 
setenv EVENT_TYPE  NEW                  # Either 'NEW', 'REFINE_SADDLE' when further converging a saddle point
                                        # Or "REFINE_AND_RELAX", to refine at the saddle
			                            # and check the final minimum
setenv Temperature                     2.5   # Temperature in eV, if negative always reject the event
setenv Max_Number_Events               10  # Maximum number of events
setenv Type_of_Events                  local   # Initial move for events - global or local
setenv Radius_Initial_Deformation       3.692  #for Zr-center selection, 3.692 for Cu-center selection
setenv Central_Atom                           1
#setenv sym_break_dist                  0.001   # Breaks the symmetry of the crystal by randomly displacing
                                               # all atoms by this distance
setenv ENERGY_CALC                       LAM    # Still use the name of SWP, but already modified codes for EAM 


setenv Activation_MaxIter                150   # Maximum number of iteraction for reaching the saddle point
setenv Increment_Size                    0.1   # Overall scale for the increment moves
setenv Force_Threshold_Perp_Rel          0.025   # Threshold for perpendicular relaxation

setenv DT_MAX_FIRE                      0.2
#----------------------------- LOCAL FORCE
setenv LOCAL_FORCE                     .true.  # Use local forces
setenv INNER_REGION_RADIUS             15.0     # Radius of inner region          
setenv OUTER_REGION_WIDTH              6.0      #Thickness of outer region shell
#_____________________________ HARMONIC WELL

setenv Initial_Step_Size                0.5   # Size of initial displacement, in A
setenv Basin_Factor                      3.0   # Factor multiplying Increment_Size for leaving the basin
setenv Max_Perp_Moves_Basin                2   # Maximum number of perpendicular steps leaving basin
setenv Min_Number_KSteps                  0   # Min. number of ksteps before calling lanczos 
setenv Eigenvalue_Threshold              -1.0   # Eigenvalue threshold for leaving basin
setenv Max_Iter_Basin                     20   # Maximum number of iteraction for leaving the basin (kter)
#_____________________________ LANCZOS
setenv Lanczos_of_minimum            .False.    # Calculation of the Hessian for each minimum
setenv Number_Lanczos_Vectors             15   # Number of vectors included in lanczos procedure
setenv delta_disp_Lanczos              0.005   # Step of the numerical derivative of forces in lanczos (Ang)
#_____________________________ CONVERGENCE
setenv Exit_Force_Threshold              0.1   # Threshold for convergence at saddle point
setenv Prefactor_Push_Over_Saddle        0.15   # Fraction of displacement over the saddle
setenv Save_Conf_Int                  .False.   # Save the configuration at every step?
#_____________________________ INPUT 
setenv FILECOUNTER               filecounter   # File tracking  the file (event) number - facultative
setenv REFCONFIG                   refconfig   # Reference configuration (actual local minimum)
#_____________________________ OUTPUT 
setenv LOGFILE                      log.file   # General output for message
setenv EVENTSLIST                events.list   # list of events with success or failure
setenv Write_restart_file             .False.   # It is useful only for ab-initio
setenv RESTART                   restart.dat   # current data for restarting event
setenv Write_JMOL                     .False.   # Writes the configuration in jmol format.
###### RUN THE SIMULATION #################################################################
$HOME/artn/source/ARTn_exec
#$ART_EXE
