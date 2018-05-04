#!/bin/tcsh
#____________________________ ATOMS 
setenv gnrm          1.0E-05
setenv Inflection  4 
setenv NATOMS     2000                         # Number of atoms in the problem
setenv type1      Cu             # Still use the name of Si, but already modified codes for ZrCu system
setenv type2      Zr
#_____________________________ ART 
setenv EVENT_TYPE  NEW                  # Either 'NEW', 'REFINE_SADDLE' when further converging a saddle point
                                        # Or "REFINE_AND_RELAX", to refine at the saddle
			                            # and check the final minimum
setenv Temperature                      2.5   # Temperature in eV, if negative always reject the event
setenv Max_Number_Events               10   # Maximum number of events
setenv Type_of_Events                  local   # Initial move for events - global or local
setenv Radius_Initial_Deformation       4.108  #for Zr-center selection, 3.692 for Cu-center selection
#setenv Central_Atom        1  # change here everytime
#setenv sym_break_dist                  0.001   # Breaks the symmetry of the crystal by randomly displacing
                                               # all atoms by this distance
setenv ENERGY_CALC                       LAM    # Still use the name of SWP, but already modified codes for EAM 


setenv Activation_MaxIter                150   # Maximum number of iteraction for reaching the saddle point
setenv Increment_Size                    0.5   # Overall scale for the increment moves
setenv Force_Threshold_Perp_Rel          0.025   # Threshold for perpendicular relaxation

#_____________________________ HARMONIC WELL
setenv Initial_Step_Size                0.5   # Size of initial displacement, in A
setenv Basin_Factor                      2.4   # Factor multiplying Increment_Size for leaving the basin
setenv Max_Perp_Moves_Basin                5   # Maximum number of perpendicular steps leaving basin
setenv Min_Number_KSteps                  0   # Min. number of ksteps before calling lanczos 
setenv Eigenvalue_Threshold             -0.01   # Eigenvalue threshold for leaving basin
setenv Max_Iter_Basin                     20   # Maximum number of iteraction for leaving the basin (kter)
#_____________________________ LANCZOS
setenv Lanczos_of_minimum            .False.    # Calculation of the Hessian for each minimum
setenv Number_Lanczos_Vectors             15   # Number of vectors included in lanczos procedure
setenv delta_disp_Lanczos              0.005   # Step of the numerical derivative of forces in lanczos (Ang)
#_____________________________ CONVERGENCE
setenv Exit_Force_Threshold              0.05   # Threshold for convergence at saddle point
setenv Prefactor_Push_Over_Saddle        0.5   # Fraction of displacement over the saddle
setenv Save_Conf_Int                  .False.   # Save the configuration at every step?
#_____________________________ DIIS
setenv Iterative                      .False.   # Iterative use of Lanczos & DIIS
setenv Use_DIIS                       .False.   # Use DIIS for the final convergence to saddle
setenv DIIS_Force_Threshold              1.5   # Force threshold for call DIIS
setenv DIIS_Memory                         5   # Number of vectors kepts in memory for algorithm
setenv DIIS_Check_Eigenvector         .True.   # Check that the final state is indeed a saddle
setenv DIIS_Step_size                   0.03   # prefactor multiplying forces
setenv FACTOR_DIIS                       5.0   # times Increment_Size, max allowed diis step size
setenv MAX_DIIS                          150   # max diis iterations per call
#_____________________________ INPUT 
setenv FILECOUNTER               filecounter   # File tracking  the file (event) number - facultative
setenv REFCONFIG                   refconfig   # Reference configuration (actual local minimum)
#_____________________________ OUTPUT 
setenv LOGFILE                      log.file   # General output for message
setenv EVENTSLIST                events.list   # list of events with success or failure
setenv Write_restart_file             .False.   # It is useful only for ab-initio
setenv RESTART                   restart.dat   # current data for restarting event
setenv Write_JMOL                     .False.   # Writes the configuration in jmol format.
###### RUN THE SIMULATION #######################################################################
../../../../../source/ARTn_exec
