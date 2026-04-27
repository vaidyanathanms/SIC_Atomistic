!---------------To analyze properties of bulk-sei systems------------
!---------------Version 1: Apr-13-2026-------------------------------
!---------------Parameter File: params_statics.f90-------------------
!********************************************************************

PROGRAM ANALYSIS_MAIN

  USE STATICPARAMS

  IMPLICIT NONE

! Print headers

  PRINT *, "Static analysis of polymer system with SEI .."
  PRINT *, "Starting OMP Threads .."
!$OMP PARALLEL
  nproc = OMP_GET_NUM_THREADS()
!$OMP END PARALLEL
  PRINT *, "Number of threads: ", nproc

! Call functions
  CALL READ_ANA_IP_FILE()
!!$  CALL SORTALLARRAYS()
  CALL ALLOCATE_ANALYSIS_ARRAYS()
  CALL ANALYZE_TRAJECTORYFILE()
  CALL ALLOUTPUTS()
  IF(layerana_flag) THEN
     PRINT *, "Beginning layer-wise analysis..."
     CALL LAYERWISE_MAIN()
  END IF
  
!!$  PRINT *, "Beginning dynamical analysis..."
!!$  CALL DYNAMICS_MAIN()

  CALL DEALLOCATE_ARRAYS()

! Print completion
  PRINT *, "All Calculations Completed Succesfully :)"

END PROGRAM ANALYSIS_MAIN

!--------------------------------------------------------------------

SUBROUTINE READ_ANA_IP_FILE()

  USE STATICPARAMS

  IMPLICIT NONE
  
  INTEGER :: nargs,ierr,logflag,AllocateStatus,i,j,u
  INTEGER :: ntypes_per_group
  CHARACTER(256) :: dumchar,charline

  CALL DEFAULTVALUES()

  nargs = IARGC()
  IF(nargs .NE. 1) STOP "Input incorrect"

  logflag = 0

  CALL GETARG(nargs,ana_fname)

  OPEN(unit = anaread,file=trim(ana_fname),action="read",status="old"&
       &,iostat=ierr)
  
  IF(ierr /= 0) THEN

     PRINT *, trim(ana_fname), "not found"
     STOP

  END IF

  DO

     READ(anaread,*,iostat=ierr) dumchar

     IF(ierr .LT. 0) EXIT

     ! Read file and trajectory details
     IF(dumchar == 'datafile') THEN
        
        READ(anaread,*,iostat=ierr) data_fname
        CALL READ_DATAFILE()
        readdataflag = 1
        
     ELSEIF(dumchar == 'trajectory_file') THEN

        READ(anaread,*,iostat=ierr) traj_fname

     ELSEIF(dumchar == 'nframes') THEN

        READ(anaread,*,iostat=ierr) nframes

     ELSEIF(dumchar == 'skipfr') THEN

        READ(anaread,*,iostat=ierr) skipfr

     ELSEIF(dumchar == 'freqfr') THEN

        READ(anaread,*,iostat=ierr) freqfr
     
     ! Read group details
     ELSEIF(dumchar == 'num_groups') THEN

        IF(readdataflag /=1) STOP "Datafile should be read first ..."
        
        grpflag = 1
        READ(anaread,*,iostat=ierr) ngroups

        ALLOCATE(allgrp_typarr(0:ntotatoms,ngroups),stat&
             &=AllocateStatus)
        IF(AllocateStatus/=0) STOP "did not allocate allgrp_typarr"

        allgrp_typarr = 0 ! ZERO everything
        
        DO i = 1,ngroups
           
           READ(anaread,'(A)',iostat=ierr) charline
           IF(ierr /= 0) THEN
              PRINT *, "Error reading group size at group-",i
              STOP
           END IF
           
           READ(charline,*) allgrp_typarr(0,i), ntypes_per_group
           ALLOCATE(types_in_group(ntypes_per_group),stat&
                &=AllocateStatus)
           IF(AllocateStatus/=0) STOP "did not allocate types_in_group&
                &"

           READ(charline,*,iostat=ierr) u, ntypes_per_group,&
                & (types_in_group(j), j = 1,ntypes_per_group)
           IF(ierr /= 0) THEN
              PRINT *, "Mismatch in #of groups in group", i
              STOP
           END IF
           
           CALL FILL_GROUP_ARRAY(i,ntypes_per_group,types_in_group)
           DEALLOCATE(types_in_group)

        END DO

        CALL OUTPUT_ALL_GROUPS()
        
     ELSEIF(dumchar == 'ion_type') THEN
        
        READ(anaread,*,iostat=ierr) iontype
        
     ELSEIF(dumchar == 'cion_type') THEN

        READ(anaread,*,iostat=ierr) c_iontype

     ELSEIF(dumchar == 'pion_type') THEN

        READ(anaread,*,iostat=ierr) p_iontype       

     !Density profiles
     ELSEIF(dumchar == 'compute_dens') THEN

        IF(grpflag /= 1) STOP "Define groups before density..."
        densflag = 1

        READ(anaread,*,iostat=ierr) ndens_grps, major_axis, d_maxbin&
             &,densfreq
        ALLOCATE(densarray(0:d_maxbin,0:ndens_grps),stat&
             &=AllocateStatus)
        IF(AllocateStatus/=0) STOP "did not allocate densarray"

        DO i = 0,d_maxbin
           DO j = 0,ndens_grps
              densarray(i,j) = 0.0
           END DO
        END DO
        
        densarray(0,0) = 9999 ! dummy type for distance "r"
        READ(anaread,*,iostat=ierr) (densarray(0,j), j = 1,ndens_grps)
     !Layer-wise static properties
     ELSEIF(dumchar == 'compute_interfaces') THEN

        interfflag = 1;  layerana_flag = 1
        IF(densflag == 0) STOP "Set density groups before interfaces..&
             &"
        READ(anaread,*,iostat=ierr) interfgrp_a, interfgrp_b

     ELSEIF(dumchar == 'layer_groups') THEN
        layer_grpflag = 1
        READ(anaread,*,iostat=ierr) nlayer_groups, nmax_layers,&
             & epsinc, epspre, segper

        ALLOCATE(all_layergrp_typarr(0:ntotatoms,nlayer_groups),stat&
             &=AllocateStatus)
        IF(AllocateStatus/=0) STOP "did not allocate all_layergrp_typa&
             &rr"

        all_layergrp_typarr = 0 ! ZERO everything
        
        DO i = 1,nlayer_groups
           
           READ(anaread,'(A)',iostat=ierr) charline
           IF(ierr /= 0) THEN
              PRINT *, "Error reading layer group size at group-",i
              STOP
           END IF
           
           READ(charline,*) all_layergrp_typarr(0,i), ntypes_per_group
           ALLOCATE(types_in_group(ntypes_per_group),stat&
                &=AllocateStatus)
           IF(AllocateStatus/=0) STOP "did not allocate types_in_group&
                &"

           READ(charline,*,iostat=ierr) u, ntypes_per_group,&
                & (types_in_group(j), j = 1,ntypes_per_group)
           IF(ierr /= 0) THEN
              PRINT *, "Mismatch in #of groups in group", i
              STOP
           END IF
           
           CALL FILL_LAYERGROUP_ARRAY(i,ntypes_per_group&
                &,types_in_group)
           DEALLOCATE(types_in_group)

        END DO

     ELSEIF(dumchar == 'layer_rdf') THEN

        rdf2dflag = 1
        READ(anaread,*,iostat=ierr) rdf2dfreq,rdf2dmaxbin,npairs_2drdf
        
        ALLOCATE(pairs_2drdf_arr(npairs_2drdf,3),stat = AllocateStatus)
        IF(AllocateStatus/=0) STOP "did not allocate pairs_rdf"

        pairs_2drdf_arr = 0
        DO i = 1,npairs_2drdf

           READ(anaread,*,iostat=ierr) pairs_2drdf_arr(i,1),&
                & pairs_2drdf_arr(i,2)
           CALL COUNT_ATOMS_WITH_TYPE_I(pairs_2drdf_arr(i,2)&
                &,pairs_2drdf_arr(i,3))

        END DO

        
     !Global static properties
     ELSEIF(dumchar == 'compute_rdf') THEN

        rdfflag = 1
        READ(anaread,*,iostat=ierr) rdffreq, rmaxbin, rdomcut,npairs
        
        ALLOCATE(pairs_rdf(npairs,3),stat = AllocateStatus)
        IF(AllocateStatus/=0) STOP "did not allocate pairs_rdf"
      
        DO i = 1,npairs

           READ(anaread,*,iostat=ierr) pairs_rdf(i,1), pairs_rdf(i,2)

        END DO

     ELSEIF(dumchar == 'compute_boundpolrdf') THEN

        bfrdf_calc = 1
        READ(anaread,*,iostat=ierr) rcatpol_cut1

     ELSEIF(dumchar == 'compute_clust') THEN

        READ(anaread,*,iostat=ierr) clust_calc

     ELSEIF(dumchar == 'compute_rg') THEN

        rgcalc = 1
        READ(anaread,*,iostat=ierr) rgfreq,rgall,rgavg

     ELSEIF(dumchar == 'compute_catanneigh') THEN

        catan_neighcalc = 1
        READ(anaread,*,iostat=ierr) neighfreq,maxneighsize,rneigh_cut
        
     !Here onwards dynamic properties
     ELSEIF(dumchar == 'compute_iondiff') THEN

        READ(anaread,*,iostat=ierr) ion_diff, delta_t
        ion_dynflag = 1

     ELSEIF(dumchar == 'compute_piondiff') THEN !stored in polyionarray

        READ(anaread,*,iostat=ierr) pion_diff, delta_t
        pion_dynflag = 1

     ELSEIF(dumchar == 'compute_ciondiff') THEN

        READ(anaread,*,iostat=ierr) cion_diff, delta_t
        cion_dynflag = 1

     ELSEIF(dumchar == 'compute_catanrestime') THEN
        
        READ(anaread,*,iostat=ierr) rcatan_cut
        catan_autocfflag = 1
        ion_dynflag = 1; cion_dynflag = 1

     ELSEIF(dumchar == 'compute_catpolrestime') THEN

        READ(anaread,*,iostat=ierr) rcatpol_cut2
        catpol_autocfflag = 1; pion_dynflag = 1
        ion_dynflag = 1

     ! Read log filename
     ELSEIF(dumchar == 'log_file') THEN

        READ(anaread,*,iostat=ierr) log_fname
        logflag  = 1

     ELSE
        
        PRINT *, "unknown keyword: ", trim(dumchar)
        STOP

     END IF

  END DO

  IF(logflag == 0) log_fname = "log_"//trim(adjustl(traj_fname))
  OPEN(unit = logout,file=trim(log_fname),action="write",status="repla&
       &ce",iostat=ierr)

  PRINT *, "Analysis input file read finished .."

  CALL SANITY_CHECK_IONTYPES()

END SUBROUTINE READ_ANA_IP_FILE

!--------------------------------------------------------------------

SUBROUTINE DEFAULTVALUES()

  USE STATICPARAMS
  IMPLICIT NONE

  ! Frame, molecules and processor details
  nframes = 0; skipfr = 0; freqfr = 0; nfrcntr = 0
  nchains = 0; atperchain = 0

  ! Initialize flags
  grpflag = 0; densflag = 0; layerana_flag = 0
  interfflag = 0;
  rgall = 0; rgcalc = 0; rdfflag = 0
  ion_dynflag = 0; cion_dynflag = 0; pion_dynflag = 0
  ion_diff = 0; cion_diff = 0; pion_diff = 0
  bfrdf_calc = 0
  catan_autocfflag = 0; catpol_autocfflag = 0
  rdf2dflag = 0
  
  ! Initialize iontypes
  c_iontype = -1; p_iontype = -1; iontype = -1

  !Initialize system quantities
  ioncnt = 0; c_ioncnt = 0; p_ioncnt= 0

  ! Initialize distributions and frequencies
  rdffreq = 0; rgfreq = 1; densfreq = 1
  
  ! Initialize structural quantities
  rdomcut = 10.0;  rmaxbin = 100; rbinval = REAL(rdomcut)&
       &/REAL(rmaxbin)
  rcatan_cut = 0.0; rneigh_cut = 0.0

  ! Initialize density profile quantities
  d_maxbin   = 200; major_axis = 3

  ! Initialize layerwise property quantities
  num_mons_per_layer = 0; nmax_layers = 0
  epspre = 0; epsinit = 0; segper=0; epsinc = 0
  
  ! Interfacial profiles
  interfgrp_a = 0; interfgrp_b = 0; maxinterf = 4
  
  ! Initialize structural averages
  rvolavg = 0; rgavg = 0; dbinavg = 0.0; ddenavg = 0.0
  major_boxval = 0.0

  ! Initialize layer-wise structural averages
  rdf2dvolavg = 0.0; rdf2dfreq = 1; rdf2dmaxbin = 50
  rdf2dbinavg = 0.0
  
  ! Initialize dynamical quantities
  rcatpol_cut1 = 0.0; rcatpol_cut2 = 0.0
  
END SUBROUTINE DEFAULTVALUES

!--------------------------------------------------------------------

SUBROUTINE READ_DATAFILE()

  USE STATICPARAMS

  IMPLICIT NONE

  INTEGER :: i,j,k,ierr,u,AllocateStatus,imax
  INTEGER :: flag, cntr, nwords
  INTEGER :: aid,molid,atype,ix,iy,iz
  REAL    :: charge,rx,ry,rz
  REAL    :: xlo,xhi,ylo,yhi,zlo,zhi
  CHARACTER(256) :: rline,dumchar

  CALL COMPUTE_INIT_NLINES(imax)

  OPEN(unit=inpread,file = trim(data_fname),action =&
       & 'read', status='old',iostat=ierr) 
  
  IF(ierr .NE. 0) STOP "Data file not found"

  WRITE(logout,*) "Datafile used: ", trim(adjustl(data_fname))

  ntotatoms = 0;ntotbonds=0;ntotangls=0;ntotdihds=0;ntotimprs=0
  atomflag =0;velflag = 0;bondflag=0;anglflag=0;dihdflag=0;imprflag=0

  READ(inpread,*)
  READ(inpread,*) 

  DO i = 1,imax-2 !Two lines before "Masses" keyword
      
     READ(inpread,*) u, dumchar

     IF(dumchar == "atoms") THEN
        ntotatoms = u
     ELSEIF(dumchar == "bonds") THEN
        ntotbonds = u
     ELSEIF(dumchar == "angles") THEN
        ntotangls = u
     ELSEIF(dumchar == "dihedrals") THEN
        ntotdihds = u
     ELSEIF(dumchar == "impropers") THEN
        ntotdihds = u
     ELSEIF(dumchar == "atom" .OR. dumchar == "atomtypes") THEN
        ntotatomtypes = u
     ELSEIF(dumchar == "bond" .OR. dumchar == "bondtypes") THEN
        ntotbondtypes = u
     ELSEIF(dumchar == "angle" .OR. dumchar == "atomtypes") THEN
        ntotangltypes = u
     ELSEIF(dumchar == "dihedral" .OR. dumchar == "dihedraltypes") THEN
        ntotdihdtypes = u
     ELSEIF(dumchar == "improper" .OR. dumchar == "impropertypes") THEN
        ntotimprtypes = u       
     END IF
     
  END DO

  READ(inpread,*)
  READ(inpread,*) xlo, xhi
  READ(inpread,*) ylo, yhi
  READ(inpread,*) zlo, zhi
  
  box_xl = xhi - xlo
  box_yl = yhi - ylo
  box_zl = zhi - zlo

  PRINT *, "x-box  ", "y-box  ", "z-box  "
  PRINT *, box_xl, box_yl, box_zl

  PRINT *, "STATISTICS"
  PRINT *, "Number of atoms/atomtypes: " , ntotatoms,ntotatomtypes
  PRINT *, "Number of bonds/bondtypes: " , ntotbonds,ntotbondtypes
  PRINT *, "Number of angles/angletypes: " , ntotangls,ntotangltypes
  PRINT *, "Number of diheds/dihedtypes: " , ntotdihds,ntotdihdtypes
  flag = 0; cntr = 0

  CALL ALLOCATE_TOPO_ARRAYS()

  DO 

     READ(inpread,*,iostat=ierr) dumchar

     IF(ierr .LT. 0) EXIT

     !READ DATA HERE FOR CHARGES AND MOLID
     !READ EVERYTHING AND OVERWRITE LATER
     IF(trim(dumchar) == "Atoms") THEN
             
        atomflag = 1
        print *, "Reading ", trim(dumchar), " info"

        DO j = 1,ntotatoms

           READ(inpread,*) aid,molid,atype,charge,rx,ry,rz

           rx = rx - xlo
           ry = ry - ylo
           rz = rz - zlo

           aidvals(aid,1)     = aid
           aidvals(aid,2)     = molid
           aidvals(aid,3)     = atype
           charge_lmp(aid,1)  = charge
           rxyz_lmp(aid,1)    = rx
           rxyz_lmp(aid,2)    = ry
           rxyz_lmp(aid,3)    = rz

        END DO

     END IF

     IF(trim(dumchar) == "Masses") THEN

        ALLOCATE(masses(ntotatomtypes,1),stat = AllocateStatus)
        IF(AllocateStatus/=0) STOP "did not allocate masses"
         
        DO j = 1,ntotatomtypes

           READ(inpread,*) u, masses(u,1)

        END DO

     END IF

     IF(trim(dumchar) == "Velocities") THEN
             
        velflag = 1
        print *, "Reading ", trim(dumchar), " info"

        DO j = 1,ntotatoms
           
           READ(inpread,*) aid,vel_xyz(aid,2),vel_xyz(aid,3)&
                &,vel_xyz(aid,4)
           vel_xyz(aid,1) = aid

        END DO

        CALL CHECK_MOMENTUM(0)

     END IF

     IF(trim(dumchar) == "Bonds") THEN
             
        bondflag = 1
        print *, "Reading ", trim(dumchar), " info"

        DO j = 1,ntotbonds

           READ(inpread,*) bond_lmp(j,1),bond_lmp(j,2),bond_lmp(j,3)&
                &,bond_lmp(j,4)

        END DO
        
     END IF

     IF(trim(dumchar) == "Angles") THEN
             
        anglflag = 1
        print *, "Reading ", trim(dumchar), " info"

        DO j = 1,ntotangls

           READ(inpread,*) angl_lmp(j,1),angl_lmp(j,2),angl_lmp(j,3)&
                &,angl_lmp(j,4),angl_lmp(j,5)

        END DO

     END IF

     IF(trim(dumchar) == "Dihedrals") THEN
             
        dihdflag = 1
        print *, "Reading", trim(dumchar), "info"

        DO j = 1,ntotdihds

           READ(inpread,*) dihd_lmp(j,1),dihd_lmp(j,2),dihd_lmp(j,3)&
                &,dihd_lmp(j,4),dihd_lmp(j,5),dihd_lmp(j,6)

        END DO

     END IF
  
     IF(trim(dumchar) == "Impropers") THEN
             
        imprflag = 1
        print *, "Reading", trim(dumchar), "info"

        DO j = 1,ntotimprs

           READ(inpread,*) impr_lmp(j,1),impr_lmp(j,2),impr_lmp(j,3)&
                &,impr_lmp(j,4),impr_lmp(j,5),impr_lmp(j,6)

        END DO

     END IF

  END DO
  
  PRINT *, "Datafile read finished..."

END SUBROUTINE READ_DATAFILE

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_INIT_NLINES(imax)

  USE STATICPARAMS

  IMPLICIT NONE

  INTEGER, INTENT(OUT) :: imax
  INTEGER :: init, pos, ipos,u,nwords,lcnt,ierr, flag
  CHARACTER(LEN=120) :: charline
  

  OPEN(unit=inpread,file = trim(data_fname),action =&
       & 'read', status='old',iostat=ierr) 
  
  IF(ierr .NE. 0) STOP "Data file not found"
  
  lcnt = 0; flag = -1

  READ(inpread,*)

  DO 

     READ(inpread,'(A)',iostat=ierr) charline     

     IF(len_trim(charline) == 0) CYCLE

     lcnt = lcnt + 1
     pos = 1
     nwords = 0

     IF(INDEX(trim(charline),"Masses") > 0) THEN

        imax = lcnt-2
        flag = 1
        EXIT

     END IF

  END DO

  IF(flag == -1) THEN
     
     STOP "Masses keyword not found in datafile"
     
  END IF

  CLOSE(inpread)

END SUBROUTINE COMPUTE_INIT_NLINES

!--------------------------------------------------------------------

SUBROUTINE FILL_GROUP_ARRAY(group_id,ntypes_per_group,types_in_group_id)

  USE STATICPARAMS

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: group_id, ntypes_per_group
  INTEGER, INTENT(IN) :: types_in_group_id(1:ntypes_per_group)
  INTEGER :: i,j,cnt_atoms

  cnt_atoms = 0

  DO i = 1,ntotatoms

     IF(ANY(aidvals(i,3) == types_in_group_id)) THEN

        allgrp_typarr(i,group_id) = 1
        cnt_atoms = cnt_atoms + 1

     END IF

  END DO
  
  ! Write out data
  WRITE(logout,*) "*********Group ID: ",group_id, " ***************"
  WRITE(logout,*) "Group-ID/ntypes_in_group/types: ", group_id,&
       & ntypes_per_group,types_in_group
  WRITE(logout,*) "Number of atoms in group-", group_id, " is ",&
       & cnt_atoms
  WRITE(logout,*) "************************************************"
  PRINT *, "Number of atoms in group-", group_id, " is ", cnt_atoms

END SUBROUTINE FILL_GROUP_ARRAY

!--------------------------------------------------------------------

SUBROUTINE FILL_LAYERGROUP_ARRAY(group_id,ntypes_per_group&
     &,types_in_group_id)

  USE STATICPARAMS

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: group_id, ntypes_per_group
  INTEGER, INTENT(IN) :: types_in_group_id(1:ntypes_per_group)
  INTEGER :: i,j,cnt_atoms

  cnt_atoms = 0

  DO i = 1,ntotatoms

     IF(ANY(aidvals(i,3) == types_in_group_id)) THEN

        all_layergrp_typarr(i,group_id) = 1
        cnt_atoms = cnt_atoms + 1

     END IF

  END DO
  
  ! Write out data
  WRITE(logout,*) "*Layer-Group ID: ",group_id, " ***************"
  WRITE(logout,*) "Layer-Group-ID/ntypes_in_group/types: ", group_id,&
       & ntypes_per_group,types_in_group
  WRITE(logout,*) "Number of atoms in group-", group_id, " is ",&
       & cnt_atoms
  WRITE(logout,*) "************************************************"
  PRINT *, "Number of atoms in group-", group_id, " is ", cnt_atoms

END SUBROUTINE FILL_LAYERGROUP_ARRAY

!--------------------------------------------------------------------

SUBROUTINE OUTPUT_ALL_GROUPS()

  USE STATICPARAMS
  IMPLICIT NONE

  INTEGER :: i,j,ierr,cntr

  dum_fname = "allgroups.dat"
  OPEN(unit = dumwrite,file =trim(dum_fname),action="write"&
       &,status="replace",iostat=ierr)
  IF(ierr /= 0) THEN
     PRINT *, "Could not open", trim(dum_fname)
  END IF

  DO i = 1,ngroups

     WRITE(dumwrite,*) "GROUP-ID: ", allgrp_typarr(0,i)
     cntr = 0
     
     DO j = 1,ntotatoms

        IF(allgrp_typarr(j,i) == 1) THEN

           WRITE(dumwrite,'(I0,2X)',advance="no") aidvals(j,1)
           cntr = cntr + 1
           IF(MOD(cntr,8) == 0 .AND. cntr /= 0) WRITE(dumwrite,*)
           
        END IF

     END DO

     WRITE(dumwrite,*)

  END DO

  CLOSE(dumwrite)

END SUBROUTINE OUTPUT_ALL_GROUPS

!--------------------------------------------------------------------

SUBROUTINE COUNT_ATOMS_WITH_TYPE_I(inptype,outcnt)

  USE STATICPARAMS
  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: inptype
  INTEGER, INTENT(OUT) :: outcnt
  INTEGER :: atcnt

  outcnt = 0
  DO atcnt = 1,ntotatoms
     IF(aidvals(atcnt,3) == inptype) outcnt = outcnt + 1
  END DO

END SUBROUTINE COUNT_ATOMS_WITH_TYPE_I
!--------------------------------------------------------------------

SUBROUTINE ANALYZE_TRAJECTORYFILE()

  USE STATICPARAMS
  IMPLICIT NONE

  INTEGER :: i,j,aid,ierr,atchk,atype,jumpfr,jout
  REAL :: xlo,xhi,ylo,yhi,zlo,zhi

  atchk = 0
  OPEN(unit = 15,file =trim(traj_fname),action="read",status="old"&
       &,iostat=ierr)

  IF(ierr /= 0) STOP "trajectory file not found"

  PRINT *, "Trajectory file used: ",trim(adjustl(traj_fname))
  WRITE(logout,*) "Trajectory file used: "&
       &,trim(adjustl(traj_fname))

  
  PRINT *, "Analyzing trajectory file..."

  CALL STRUCT_INIT()
  CALL OPEN_STRUCT_OUTPUT_FILES()

  DO i = 1,skipfr

     DO j = 1,ntotatoms+9

        READ(15,*) 

     END DO

     IF(mod(i,100) == 0) PRINT *, "Skipped ", i, "frames"

  END DO

  DO i = 1,nframes

     nfrcntr = nfrcntr + 1
     IF(mod(i,100) == 0) PRINT *, "Processing ", i+1,"th frame"

     READ(15,*)
     READ(15,*) timestep

     READ(15,*) 
     READ(15,*) atchk

     READ(15,*) 
     READ(15,*) xlo, xhi
     READ(15,*) ylo, yhi
     READ(15,*) zlo, zhi

     READ(15,*)

     box_xl = xhi - xlo
     box_yl = yhi - ylo
     box_zl = zhi - zlo
     
     boxx_arr(i)  = box_xl
     boxy_arr(i)  = box_yl
     boxz_arr(i)  = box_zl

     IF(major_axis == 1) major_boxval = box_xl
     IF(major_axis == 2) major_boxval = box_yl
     IF(major_axis == 3) major_boxval = box_zl
     
     DO j = 1,atchk

        READ(15,*) aid,atype,rxyz_lmp(aid,1),rxyz_lmp(aid,2)&
             &,rxyz_lmp(aid,3)

        IF(atype .NE. aidvals(aid,3)) THEN

           PRINT *, "Incorrect atom ids"
           PRINT *, i,j,aid,atype,aidvals(aid,3)
           STOP

        END IF
        
     END DO
!$OMP PARALLEL FIRSTPRIVATE(i) PRIVATE(j)
!$OMP DO
     DO j = 1,atchk

        rxyz_lmp(j,1) = rxyz_lmp(j,1) - xlo
        rxyz_lmp(j,2) = rxyz_lmp(j,2) - ylo
        rxyz_lmp(j,3) = rxyz_lmp(j,3) - zlo

        IF(layerana_flag) THEN

           trx_lmp(j,i) = rxyz_lmp(j,1)
           try_lmp(j,i) = rxyz_lmp(j,2)
           trz_lmp(j,i) = rxyz_lmp(j,3)

        END IF

     END DO
!$OMP END DO        
!$OMP END PARALLEL
  
     IF(i == 1) PRINT *, "Beginning statics analysis..."
     CALL STRUCT_MAIN(nfrcntr)
     
     DO jumpfr = 1,freqfr

        READ(15,*)
        READ(15,*)        
        READ(15,*)
 
        READ(15,*) atchk

        DO j = 1,atchk+5

           READ(15,*) 

        END DO

     END DO

  END DO

  CLOSE(15)

  PRINT *, "Succesfully closed trajectory file ..."


END SUBROUTINE ANALYZE_TRAJECTORYFILE

!--------------------------------------------------------------------

SUBROUTINE STRUCT_INIT()

  USE STATICPARAMS
  IMPLICIT NONE

  INTEGER :: i,j,t1,t2,norm,acnt,fcnt,a1id,molid,flagch,flagpr,jmax
  INTEGER :: AllocateStatus

  IF(rdfflag) THEN

     rdfarray = 0.0
     rbinval = rdomcut/REAL(rmaxbin)

     DO i = 1, npairs

        t1 = 0; t2 = 0

        DO j = 1,ntotatoms

           IF(aidvals(j,3) == pairs_rdf(i,1)) t1 = t1+1
           IF(aidvals(j,3) == pairs_rdf(i,2)) t2 = t2+1
           
        END DO

        IF(pairs_rdf(i,1) == pairs_rdf(i,2)) THEN
           pairs_rdf(i,3) = t1*(t1-1) !g_AA(r)
        ELSE
           pairs_rdf(i,3) = t1*t2 !g_AB(r)
        END IF
 
     END DO

  END IF

END SUBROUTINE STRUCT_INIT

!--------------------------------------------------------------------

SUBROUTINE STRUCT_MAIN(tval)

  USE STATICPARAMS
  IMPLICIT NONE

  INTEGER, INTENT(IN):: tval
  INTEGER :: t1, t2
  INTEGER :: clock_rate, clock_max
  INTEGER :: i
  
!!$  IF(rgcalc .AND. mod(tval-1,rgfreq)==0) CALL COMPUTE_RADGYR(tval)

  IF(grpflag) THEN

     IF(tval == 1) THEN
        CALL SYSTEM_CLOCK(t1,clock_rate,clock_max)
        CALL COMPUTE_DENSPROFILES(tval)
        CALL SYSTEM_CLOCK(t2,clock_rate,clock_max)
        PRINT *, 'Elapsed real time for density analysis: ',REAL(t2&
             &-t1)/REAL(clock_rate), ' seconds'

     END IF
     IF(mod(tval-1,densfreq)==0) CALL COMPUTE_DENSPROFILES(tval)     

  END IF
        
  IF(rdfflag) THEN
     
     IF(tval == 1) THEN

        CALL SYSTEM_CLOCK(t1,clock_rate,clock_max)
        CALL COMPUTE_RDF(tval)
        CALL SYSTEM_CLOCK(t2,clock_rate,clock_max)
        PRINT *, 'Elapsed real time for RDF analysis: ',REAL(t2&
             &-t1)/REAL(clock_rate), ' seconds'

     END IF
     IF(mod(tval-1,rdffreq)==0) CALL COMPUTE_RDF(tval)     
     
  END IF

  IF(catan_neighcalc) THEN
     
     IF(tval == 1) THEN

        cat_an_neighavg = 0.0; an_cat_neighavg=0.0
        CALL SYSTEM_CLOCK(t1,clock_rate,clock_max)
        CALL CAT_AN_NEIGHS()
        CALL SYSTEM_CLOCK(t2,clock_rate,clock_max)
        PRINT *, 'Elapsed real time for neighbor analysis: ',REAL(t2&
             &-t1)/REAL(clock_rate), ' seconds'

     END IF
     
     IF(mod(tval,neighfreq) == 0) CALL CAT_AN_NEIGHS()
     
  END IF

  IF(bfrdf_calc) THEN
     
     IF(tval == 1) THEN
        rdf_p_fb = 0.0; rdf_p_ff=0.0; rdf_p_bb = 0.0        
        CALL SYSTEM_CLOCK(t1,clock_rate,clock_max)
        CALL  SORT_POLY_FREE_BOUND_COMPLEX(tval)
        CALL SYSTEM_CLOCK(t2,clock_rate,clock_max)
        PRINT *, 'Elapsed real time for bound/free RDF: ',REAL(t2&
             &-t1)/REAL(clock_rate), ' seconds'           
     END IF
     
     IF(mod(tval,rdffreq) == 0) CALL&
          & SORT_POLY_FREE_BOUND_COMPLEX(tval)
     
  END IF

  IF(clust_calc) THEN

     IF(tval == 1) THEN

        clust_avg = 0
        CALL SYSTEM_CLOCK(t1,clock_rate,clock_max)
        CALL CLUSTER_ANALYSIS(tval)
        CALL SYSTEM_CLOCK(t2,clock_rate,clock_max)
        PRINT *, 'Elapsed real time for cluster analysis= ',REAL(t2&
             &-t1)/REAL(clock_rate), ' seconds'           
     ELSE
        
        CALL CLUSTER_ANALYSIS(tval)
     
     END IF

  END IF


END SUBROUTINE STRUCT_MAIN

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_DENSPROFILES(tval)

  USE STATICPARAMS

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: tval
  INTEGER :: i,j,ibin,grtype,dval,grp_cntr,grp_id,grp_colid,atype
  REAL,DIMENSION(0:d_maxbin-1,1:ndens_grps) :: grpdens
  REAL :: dbinval,rval,ddenval,vnorm
  INTEGER :: flag

  dbinval  = major_boxval/d_maxbin
  dbinavg  = dbinavg + dbinval
  ddenval  = major_boxval/(box_xl*box_yl*box_zl*dbinval)
  ddenavg  = ddenavg + ddenval
  boxlzavg = boxlzavg + boxval

  IF(tval == 1) THEN
     
     DO i = 1,ndens_grps
        
        flag = -1
        
        DO j = 1,ngroups
           
           IF(INT(densarray(0,i)) == allgrp_typarr(0,j)) THEN
              
              flag = 0
              EXIT
              
           END IF
           
        END DO
        
        IF (flag == -1) THEN
           
           PRINT *, "All groups in densarray should be defined as grou&
                &ps"
           PRINT *, "Unknown group ID in densarray", INT(densarray(0,i))
           STOP
           
        END IF
        
     END DO
     
  END IF
  
!$OMP PARALLEL
  
!$OMP DO PRIVATE(grp_cntr,i)

  DO grp_cntr = 1,ndens_grps

     DO i = 0,d_maxbin-1
        
        grpdens(i,grp_cntr) = 0.0

     END DO

  END DO

!$OMP END DO

!$OMP DO PRIVATE(i,j,grp_cntr,grp_id,grp_colid,rval,atype,ibin) REDUCTION(+:grpdens)

  DO grp_cntr = 1,ndens_grps

     grp_id = densarray(0,grp_cntr)

     DO j = 1,ngroups

        IF(grp_id == allgrp_typarr(0,j)) THEN

           grp_colid = j

        END IF

     END DO
           
     DO i = 1,ntotatoms
     
        IF(allgrp_typarr(i,grp_colid)) THEN

           atype = aidvals(i,3)
           IF(major_axis == 1) rval = rxyz_lmp(i,1) 
           IF(major_axis == 2) rval = rxyz_lmp(i,2) 
           IF(major_axis == 3) rval = rxyz_lmp(i,3) 
           
           ibin = FLOOR(rval/dbinval)

           IF(ibin .LT. d_maxbin) THEN

              grpdens(ibin,grp_cntr) = grpdens(ibin,grp_cntr) + masses(atype,1)

           END IF

        END IF

     END DO

  END DO

!$OMP END DO

!$OMP DO PRIVATE(grp_cntr,i)
  DO grp_cntr = 1,ndens_grps

     DO i = 0,d_maxbin-1

        grpdens(i,grp_cntr) = grpdens(i,grp_cntr)*ddenval

     END DO

  END DO
!$OMP END DO  
  
!$OMP DO PRIVATE(grp_cntr,i)

  DO grp_cntr = 1,ndens_grps

     DO i = 0,d_maxbin-1
        !i+1 for densarray since i = 0 is the group type
        densarray(i+1,grp_cntr) = grpdens(i,grp_cntr) + densarray(i+1&
             &,grp_cntr)

     END DO

  END DO

!$OMP END DO
  
!$OMP END PARALLEL
  
END SUBROUTINE COMPUTE_DENSPROFILES

!--------------------------------------------------------------------

! pion_type is set to -1 in the beginning and is read if and only if
! pion_dynflag is activated
! This is obsolete with groups
SUBROUTINE SORTALLARRAYS()

  USE STATICPARAMS

  IMPLICIT NONE

  INTEGER :: i,j,a1type,cnt,AllocateStatus,ntotion_cnt,aid,molid
  INTEGER, DIMENSION(1:ntotatoms,2) :: dumsortarr,dumcionarr&
       &,dumpionarr

  dumsortarr = -1; dumcionarr = -1; dumpionarr = -1
  cnt = 0
  ntotion_cnt = 0

  DO i = 1,ntotatoms

     a1type = aidvals(i,3)

     IF(a1type == iontype) THEN
        ioncnt = ioncnt + 1
        dumsortarr(ioncnt,1) = i
        dumsortarr(ioncnt,2) = a1type

     ELSEIF(a1type == c_iontype) THEN
        c_ioncnt = c_ioncnt + 1
        dumcionarr(c_ioncnt,1) = i
        dumcionarr(c_ioncnt,2) = a1type

     ELSEIF(a1type == p_iontype) THEN
        p_ioncnt = p_ioncnt + 1
        dumpionarr(p_ioncnt,1) = i
        dumpionarr(p_ioncnt,2) = a1type

     END IF

  END DO
  
  ! Always identify ion and counter-ion types
  PRINT *, "Number of atoms of ion type: ", ioncnt
  PRINT *, "Number of atoms of cntion type: ", c_ioncnt
     
  ALLOCATE(ionarray(ioncnt,2),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate ionarray"
  ALLOCATE(counterarray(c_ioncnt,2),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate ionarray"
  
  ! Load ion array

  i = 0
     
  DO WHILE(dumsortarr(i+1,1) .NE. -1) 

     i = i + 1
     ionarray(i,1) = dumsortarr(i,1)
     ionarray(i,2) = dumsortarr(i,2)
     
  END DO

  IF(i .NE. ioncnt) THEN
     PRINT *, i, ioncnt
     STOP "Wrong total count in ionarray"
  END IF
  
  DO i = 1,ioncnt
     
     IF(ionarray(i,1) == -1 .OR. ionarray(i,2) == -1) THEN
        
        PRINT *, i,ionarray(i,1), ionarray(i,2)
        PRINT *, "Something wrong in assigning ionarray"
        STOP
        
     END IF
        
     IF(ionarray(i,2) .NE. iontype) THEN
           
        PRINT *, i,ionarray(i,1), ionarray(i,2)
        PRINT *, "Something wrong in ionarray type"
        STOP
        
     END IF
        
  END DO
         
  OPEN(unit = 93,file="iontypelist.txt",action="write",status="replace&
       &")
     
  WRITE(93,*) "Reference type/count: ", iontype, ioncnt
     
  DO i = 1,ioncnt
     WRITE(93,'(3(I0,1X))') i, ionarray(i,1), ionarray(i,2)
  END DO
     
  CLOSE(93)
     
  ! Load counterion array
     
  i = 0
  
  DO WHILE(dumcionarr(i+1,1) .NE. -1) 
     
     i = i + 1
     counterarray(i,1) = dumcionarr(i,1)
     counterarray(i,2) = dumcionarr(i,2)
     
  END DO
  
  IF(i .NE. c_ioncnt) THEN
     PRINT *, i, c_ioncnt
     STOP "Wrong total count in counterarray"
  END IF
     
  DO i = 1,c_ioncnt
     
     IF(counterarray(i,1) == -1 .OR. counterarray(i,2) == -1) THEN
        
        PRINT *, i,counterarray(i,1), counterarray(i,2)
        PRINT *, "Something wrong in assigning counterarray"
        STOP
        
     END IF
     
     IF(counterarray(i,2) .NE. c_iontype) THEN
        
        PRINT *, i,counterarray(i,1), counterarray(i,2)
        PRINT *, "Something wrong in counterionarray type"
        STOP
        
     END IF
        
  END DO
     
  
  OPEN(unit = 93,file="cntionlist.txt",action="write",status="repl&
       &ace")
  
  WRITE(93,*) "Reference type/count: ", c_iontype, c_ioncnt
     
  DO i = 1,c_ioncnt
     
     WRITE(93,'(3(I0,1X))') i, counterarray(i,1), counterarray(i,2)
     
  END DO
     
  CLOSE(93)
 

  ! Cluster calc requires to add iontype and c_iontype in same array
  IF (clust_calc) THEN

     ntotion_centers = ioncnt + c_ioncnt
     PRINT *, "Total number of ion centers: ", ntotion_centers
     cnt = 1

     ALLOCATE(allionids(ntotion_centers,2),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate allionids"
     ALLOCATE(clust_avg(ntotion_centers),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate clust_avg"

     allionids = 0 ! Initial allocation

     DO i = 1,ntotatoms

        a1type = aidvals(i,3)

        IF(a1type == iontype .OR. a1type == c_iontype) THEN
           
           allionids(cnt,1) = i
           allionids(cnt,2) = a1type
           cnt = cnt + 1

        END IF

     END DO

  ELSE

     ALLOCATE(allionids(1,2),stat = AllocateStatus)
     DEALLOCATE(allionids)
     ALLOCATE(clust_avg(1),stat = AllocateStatus)
     DEALLOCATE(clust_avg)

  END IF
     

  ! Polymerion array required only for diffusion systems
  IF (pion_dynflag .OR. bfrdf_calc) THEN

     PRINT *, "Number of atoms of polyion type: ",p_ioncnt

     ALLOCATE(polyionarray(p_ioncnt,2),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate polyionarray"

     i = 0

     DO WHILE(dumpionarr(i+1,1) .NE. -1) 

        i = i + 1
        polyionarray(i,1) = dumpionarr(i,1)
        polyionarray(i,2) = dumpionarr(i,2)
     
     END DO
     
     IF(i .NE. p_ioncnt) THEN
        PRINT *, i, p_ioncnt
        STOP "Wrong total count in counterarray"
     END IF
     
     DO i = 1,p_ioncnt
     
        IF(polyionarray(i,1) == -1 .OR. polyionarray(i,2) == -1) THEN
           
           PRINT *, i,polyionarray(i,1), polyionarray(i,2)
           PRINT *, "Something wrong in assigning polyionarray"
           STOP
           
        END IF
     
        IF(polyionarray(i,2) .NE. p_iontype) THEN
        
           PRINT *, i,polyionarray(i,1), polyionarray(i,2)
           PRINT *, "Something wrong in polyionarray"
           STOP
           
        END IF
     
     END DO
  
     
     OPEN(unit = 93,file="polyionlist.txt",action="write",status="repl&
          &ace")
  
     WRITE(93,*) "Reference type/count: ", p_iontype, p_ioncnt
     WRITE(93,*) "#  ","ID  ","molID  ","Type"
     DO i = 1,p_ioncnt
        aid = polyionarray(i,1)
        molid = aidvals(aid,2)
        WRITE(93,'(4(I0,1X))') i,polyionarray(i,1),molid, polyionarray(i,2)
        
     END DO
     
     CLOSE(93)

  ELSE

     ALLOCATE(polyionarray(1,2),stat = AllocateStatus)
     DEALLOCATE(polyionarray)

  END IF

END SUBROUTINE SORTALLARRAYS

!--------------------------------------------------------------------

SUBROUTINE SANITY_CHECK_IONTYPES()

  USE STATICPARAMS

  IMPLICIT NONE

  IF(ion_dynflag .OR. catan_neighcalc) THEN

     IF(iontype == -1) THEN
        
        PRINT *, "ion type undefined for neigh or diff calculation"
        STOP 

     END IF

  END IF

  IF(cion_dynflag .OR. catan_neighcalc) THEN

     IF(c_iontype == -1) THEN
        
        PRINT *, "counter-ion type undefined for neigh or diff calculation"
        STOP 

     END IF

  END IF

  IF(pion_dynflag) THEN

     IF(p_iontype == -1) THEN
        
        PRINT *, "polymer-ion type undefined for diff calculation"
        STOP 

     END IF

  END IF


END SUBROUTINE SANITY_CHECK_IONTYPES

!--------------------------------------------------------------------

SUBROUTINE MAP_REFTYPE(jin,atype,jout)
! Maps atomid into the corresponding place in array
  USE STATICPARAMS

  IMPLICIT NONE

  INTEGER :: i
  INTEGER, INTENT(IN):: jin,atype
  INTEGER, INTENT(OUT) :: jout

  jout = -1

  IF(atype == iontype) THEN

     DO i = 1,ioncnt
        
        IF(jin == ionarray(i,1)) THEN

           jout = i

           EXIT

        END IF

     END DO

  ELSEIF(atype == c_iontype) THEN

     DO i = 1,c_ioncnt
        
        IF(jin == counterarray(i,1)) THEN

           jout = i

           EXIT

        END IF

     END DO

  ELSEIF(atype == p_iontype) THEN

     DO i = 1,p_ioncnt
        
        IF(jin == polyionarray(i,1)) THEN
           
           jout = i

           EXIT
           
        END IF
        
     END DO

  END IF
  
  IF(jout == -1) THEN
     
     PRINT *, jin, atype
     STOP "Could not find a match"

  END IF


END SUBROUTINE MAP_REFTYPE

!--------------------------------------------------------------------

SUBROUTINE CAT_AN_NEIGHS()

  USE STATICPARAMS

  IMPLICIT NONE

  INTEGER :: i,j,a1id,a2id,neigh_cnt,tid
  INTEGER,DIMENSION(1:maxneighsize,0:nproc-1) :: cat_an_neigh_inst&
       &,an_cat_neigh_inst 
  REAL :: rxval, ryval, rzval, rval

  cat_an_neigh_inst = 0; an_cat_neigh_inst = 0

!$OMP PARALLEL PRIVATE(i,j,a1id,a2id,rxval,ryval,rzval,rval,neigh_cnt,tid)
!$OMP DO
  DO i = 1,ioncnt
     
     neigh_cnt = 0
     a1id = ionarray(i,1)
     tid = OMP_GET_THREAD_NUM()

     DO j = 1,c_ioncnt

        a2id = counterarray(j,1)
        
        rxval = rxyz_lmp(a1id,1) - rxyz_lmp(a2id,1) 
        ryval = rxyz_lmp(a1id,2) - rxyz_lmp(a2id,2) 
        rzval = rxyz_lmp(a1id,3) - rxyz_lmp(a2id,3) 
        
        rxval = rxval - box_xl*ANINT(rxval/box_xl)
        ryval = ryval - box_yl*ANINT(ryval/box_yl)
        rzval = rzval - box_zl*ANINT(rzval/box_zl)
        
        rval = sqrt(rxval**2 + ryval**2 + rzval**2)
        
        IF(rval .LT. rneigh_cut) THEN
           
           neigh_cnt = neigh_cnt + 1
           
        END IF

     END DO

     IF(neigh_cnt + 1 .GT. maxneighsize) THEN

        PRINT *, "Neighbor count exceeded max size"
        PRINT *, neigh_cnt, maxneighsize
        STOP

     END IF

     cat_an_neigh_inst(neigh_cnt+1,tid) = cat_an_neigh_inst(neigh_cnt&
          &+1,tid) + 1

  END DO
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL PRIVATE(i,j,a1id,a2id,rxval,ryval,rzval,rval,neigh_cnt,tid)
!$OMP DO
  DO i = 1,c_ioncnt
     
     neigh_cnt = 0
     a1id = counterarray(i,1)
     tid = OMP_GET_THREAD_NUM()

     DO j = 1,ioncnt

        a2id = ionarray(j,1)
        
        rxval = rxyz_lmp(a1id,1) - rxyz_lmp(a2id,1) 
        ryval = rxyz_lmp(a1id,2) - rxyz_lmp(a2id,2) 
        rzval = rxyz_lmp(a1id,3) - rxyz_lmp(a2id,3) 
        
        rxval = rxval - box_xl*ANINT(rxval/box_xl)
        ryval = ryval - box_yl*ANINT(ryval/box_yl)
        rzval = rzval - box_zl*ANINT(rzval/box_zl)
        
        rval = sqrt(rxval**2 + ryval**2 + rzval**2)
        
        IF(rval .LT. rneigh_cut) THEN
           
           neigh_cnt = neigh_cnt + 1
           
        END IF

     END DO

     IF(neigh_cnt + 1 .GT. maxneighsize) THEN

        PRINT *, "Neighbor count exceeded max size"
        PRINT *, neigh_cnt, maxneighsize
        STOP

     END IF

     an_cat_neigh_inst(neigh_cnt+1,tid) = an_cat_neigh_inst(neigh_cnt&
          &+1,tid) + 1

  END DO
!$OMP END DO


!$OMP DO 
  DO  i = 1,maxneighsize
     DO j = 0,nproc-1
        cat_an_neighavg(i) = cat_an_neighavg(i) + cat_an_neigh_inst(i&
             &,j)
        an_cat_neighavg(i) = an_cat_neighavg(i) + an_cat_neigh_inst(i&
             &,j)
     END DO
  END DO
!$OMP END DO

!$OMP END PARALLEL

END SUBROUTINE CAT_AN_NEIGHS

!--------------------------------------------------------------------
SUBROUTINE LAYERWISE_MAIN()

  USE STATICPARAMS
  IMPLICIT NONE

  INTEGER :: t1, t2, i
  INTEGER :: clock_rate, clock_max

  
  PRINT *, "Computing interfaces..."
  CALL SYSTEM_CLOCK(t1,clock_rate,clock_max)
  CALL COMPUTE_INTERFACES()
  CALL SYSTEM_CLOCK(t2,clock_rate,clock_max)
  PRINT *, 'Elapsed real time for interface analysis: ',REAL(t2-t1)&
       &/REAL(clock_rate), ' seconds'

  box_xl = boxx_arr(1)
  box_yl = boxy_arr(1)
  box_zl = boxz_arr(1)
  rvolval = box_xl*box_yl*box_zl

  ! Classify domain according to **number** densities
  CALL INITIAL_DOMAIN_CLASSIFICATION() 
  
  DO i = 1,nframes

     box_xl = boxx_arr(i)
     box_yl = boxy_arr(i)
     box_zl = boxz_arr(i)

     CALL COMPARTMENTALIZE_PARTICLES(i)

  END DO

  IF(rdf2dflag) CALL OUTPUT_LAYERRDF()
  
END SUBROUTINE LAYERWISE_MAIN

!--------------------------------------------------------------------
  
SUBROUTINE COMPUTE_INTERFACES()

  USE STATICPARAMS

  IMPLICIT NONE

!!$ To calculate the interface positions
  REAL*8  :: dum1, dum2, dum3, volnorm, duma, dumb
  INTEGER :: i,j,AllocateStatus,frnorm,maxdens,icount,iflag,ierr,jinit
  INTEGER :: interfcol_a, interfcol_b
  REAL, DIMENSION(d_maxbin*maxinterf) :: dummid_interarr,&
       & dumdens_interarr
  !Number of elements is arbitrary 

  dummid_interarr  = -1
  dumdens_interarr = -1

  dum_fname = "interf_"//trim(adjustl(traj_fname))
  OPEN(unit = dumwrite,file =trim(dum_fname),action="write"&
       &,status="replace",iostat=ierr)

  IF(ierr /= 0) THEN
     PRINT *, "Could not open", trim(dum_fname)
  END IF

  CALL MAP_GROUP_TO_COL(interfgrp_a,interfcol_a)
  CALL MAP_GROUP_TO_COL(interfgrp_b,interfcol_b)

  !Note densarray values start with j = 1; j=0 is 9999 or "r"
  !Using the density profiles
  icount = 0
  j = 1
  jinit = 2
  duma = densarray(j,interfcol_a)
  dumb = densarray(j,interfcol_b)
  
  !Left most end rho_a > rho_b
  IF(duma > dumb) THEN
     iflag = 1
  ELSEIF(dumb > duma) THEN !rho_a < rho_b
     iflag = 2
  ELSE !rho_a = rho_b
     icount = icount + 1
     dumdens_interarr(icount) = densarray(j,0)
     duma = densarray(j+1,interfcol_a)
     dumb = densarray(j+1,interfcol_b)
     jinit = 2
     IF(duma > dumb) iflag = 1
     IF(dumb > duma) iflag = 2
     IF(duma == dumb) STOP "Wrong density profiles"
  END IF
  
!Inside the box  
  DO j = jinit, d_maxbin

     IF(iflag == 1) THEN

        IF(densarray(j,interfcol_a) .LE. densarray(j,interfcol_b))&
             & THEN !Change sign

           icount = icount + 1
           dumdens_interarr(icount) = 0.5*(densarray(j-1,0)&
                &+densarray(j,0))
           iflag  = 2

        END IF

     ELSEIF(iflag == 2) THEN
        
        IF(densarray(j,interfcol_b) .LE. densarray(j,interfcol_a))&
             & THEN !Change sign
           icount = icount + 1
           dumdens_interarr(icount) = 0.5*(densarray(j-1,0)&
                &+densarray(j,0))
           iflag  = 1

        END IF

     END IF

  END DO

  maxdens = icount

  PRINT *, "Interfaces are calculated using densities of ",&
       & interfgrp_a, " and ", interfgrp_b
  PRINT *, "Number of interfaces obtained from density: ", maxdens
  PRINT *, dumdens_interarr(1:maxdens)

  WRITE(logout,*) "Number of interfaces from density: ", maxdens
  WRITE(logout,*) dumdens_interarr(1:maxdens)
  
! Need to determine which one needs to be used.
! Allocate to actual array. This can keep things clean

  maxinterf = maxdens

  IF(maxinterf .NE. 1 ) PRINT *, "More than one interface found.."
  
  ALLOCATE (interpos(maxinterf), stat = AllocateStatus)
  IF(AllocateStatus /=0 ) STOP "***Allocation inter not proper***"

  ALLOCATE (domtyp(maxinterf), stat = AllocateStatus)
  IF(AllocateStatus /=0 ) STOP "***Allocation domtyp not proper***"

  WRITE(logout,*) "Using density profiles to write interfaces..."

  DO i = 1, maxinterf

     IF(dumdens_interarr(i) == -1) THEN

        PRINT *, "Array updated incorrectly"
        STOP

     END IF

     interpos(i) = dumdens_interarr(i)

  END DO

  WRITE(dumwrite,*) maxinterf
  WRITE(logout,*) "Number of interfaces: ", maxinterf
  WRITE(logout,*) interpos

  DO i = 1, maxinterf

     WRITE(dumwrite,*) i, interpos(i)

  END DO

  CLOSE(dumwrite)

END SUBROUTINE COMPUTE_INTERFACES

!--------------------------------------------------------------------

SUBROUTINE MAP_GROUP_TO_COL(grp_type,col_val)

  USE STATICPARAMS
  IMPLICIT NONE
  
  INTEGER :: i
  INTEGER, INTENT(IN)  :: grp_type
  INTEGER, INTENT(OUT) :: col_val

  col_val = -1

  DO i = 1,SIZE(densarray(0,:))
        
     IF(grp_type == densarray(0,i)) THEN

        col_val = i
        EXIT
        
     END IF

  END DO

  IF(col_val == -1) THEN
     PRINT *, "Unknown type for calculating interfaces..",grp_type
     PRINT *, "Interface type should be a part of density group..."
     STOP
  END IF

END SUBROUTINE MAP_GROUP_TO_COL

!--------------------------------------------------------------------
SUBROUTINE INITIAL_DOMAIN_CLASSIFICATION(tval)

  ! Only done at t = 0 using the density profile
  USE STATICPARAMS
  IMPLICIT NONE

  INTEGER :: i,j, flagdom, AllocateStatus, flagnew
  REAL :: par_pos
  INTEGER, DIMENSION(1:maxinterf) :: domAdum, domBdum
  INTEGER, INTENT(IN) :: tval


  DO i = 1, maxinterf
        
     domAdum(i) = 0
     domBdum(i) = 0
     
  END DO

  
  DO i = 1,ntotatoms
     
     IF(major_axis == 1) THEN
        par_pos = trx_lmp(i,tval) - boxval*floor(trx_lmp(i,tval)&
             &/boxval)
     ELSE IF(major_axis == 2) THEN 
        par_pos = try_lmp(i,tval) - boxval*floor(try_lmp(i,tval)&
             &/boxval)
     ELSE
        par_pos = trz_lmp(i,tval) - boxval*floor(trz_lmp(i,tval)&
             &/boxval)
     END IF

     flagdom = 0
     IF(par_pos == 0.0000000) par_pos = 10**(-8)
     
     IF(par_pos .GE. 0.0 .AND. par_pos .LE. interpos(1)) THEN
        
        flagdom = 1
        IF(ANY(allgrp_typarr(interfgrp_a,:) == aidvals(i,3))) THEN

           domAdum(1) = domAdum(1) + 1

        ELSEIF(ANY(allgrp_typarr(interfgrp_b,:) == aidvals(i,3))) THEN

           domBdum(1) = domBdum(1) + 1

        END IF

     END IF
     
     IF(par_pos .GT. interpos(1) .AND. par_pos .LE. boxval) THEN
              
        IF(flagdom == 1) STOP "particle binned twice"

        IF(ANY(allgrp_typarr(interfgrp_a,:) == aidvals(i,3))) THEN
           domAdum(2) = domAdum(2) + 1
           
        ELSEIF(ANY(allgrp_typarr(interfgrp_b,:) == aidvals(i,3)))&
             & THEN              
           domBdum(2) = domBdum(2) + 1
           
        END IF
        
     END IF

  END DO

  !Nomenclature: domA=1,domB=2
  DO i = 1, maxinterf
     
     IF(domAdum(i) > domBdum(i)) THEN
        
        domtyp(i) = 1
        
     ELSEIF(domBdum(i) > domAdum(i)) THEN
        
        domtyp(i) = 2
        
     ELSE
        
        PRINT *,"Something wrong in distribution"
        PRINT *,"Number of particles in",i,"dom", domAdum(i),&
             & domBdum(i)
        
        WRITE(logout,*),"Something wrong in distribution"
        WRITE(logout,*),"Number of particles in",i,"dom", domAdum(i),&
             & domBdum(i)
        
        STOP
        
     END IF
     
  END DO
  
  WRITE(logout,*) "The domain types are "
  
  DO i = 1,maxinterf
     
     IF(domtyp(i) == 1) THEN
        
        WRITE(logout,*) i,domtyp(i), "Domain A"
        
     ELSE
        
        WRITE(logout,*) i,domtyp(i), "Domain B"
        
     END IF
     
  END DO
     
  ALLOCATE(widdoms(maxinterf),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate widdoms"
  
  widdoms(1) = interpos(1)
  widdoms(2) = boxval - interpos(maxinterf)
  
END SUBROUTINE INITIAL_DOMAIN_CLASSIFICATION

!--------------------------------------------------------------------

SUBROUTINE COMPARTMENTALIZE_PARTICLES(tval)

  USE STATICPARAMS

  IMPLICIT NONE

  INTEGER :: i,j,k,p,AllocateStatus,flagbin,a1id
  REAL    :: par_pos, domcheck
  REAL    :: segeps_domA,segeps_domB, eps1,eps2
  REAL    :: interin_domA, interout_domA,domwidth_domA
  REAL    :: interin_domB, interout_domB,domwidth_domB
  INTEGER :: interfdomA_col, interfdomB_col
  CHARACTER(LEN=3) :: epsnum
  INTEGER, INTENT(IN) :: tval
  INTEGER, ALLOCATABLE, DIMENSION(:)::dum_aid, dum_typ 

  !Nomenclature: Dom-A (Dom-B): left (right) of interface
  
  IF(tval == 1) THEN
     ALLOCATE(seg_typcnt(1:ntotatomtypes),stat=AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate seg_typcnt"
  END IF

  CALL MAP_GROUP_TO_COL(interfgrp_a,interfdomA_col)
  CALL MAP_GROUP_TO_COL(interfgrp_b,interfdomB_col)

  
  domcheck = 0.0; epsinit = 0.0
  ! Find the first value where DOMA starts to rise - to account for
  ! wall effects
  DO i = 1,d_maxbin ! threshold at 0.02
     IF(densarray(i,interfdomA_col) > 0.02) THEN
        domwidth_domA = interpos(1) - densarray(i,0)
        EXIT
     END IF
  END DO

  ! Since domB is asymmetric, need a separate width for domB
  DO i = d_maxbin-1,1,-1
     IF(densarray(i,interfdomB_col) > 0.02) THEN
        domwidth_domB = densarray(i,interfdomB_col) - interpos(1)
     END IF
  END DO

  segeps_domA  = segper*domwidth_domA
  segeps_domB  = segper*domwidth_domB
  
  IF(tval == 1) THEN
     WRITE(logout,*) "DomA domain width is ", domwidth_domA
     WRITE(logout,*) "DomA segmental width is ", segeps_domA
     WRITE(logout,*) "DomB domain width is ", domwidth_domB
     WRITE(logout,*) "DomB segmental width is ", segeps_domB
  END IF

  IF(tval == 1) THEN
     WRITE(logout,*) "Increment in DomA divisions: ", epsinc*segeps_domA
     WRITE(logout,*) "Width for each DomA segment: ", epspre*epsinc&
          &*segeps_domA
     WRITE(logout,*) "Increment in DomB divisions: ", epsinc&
          &*segeps_domB
     WRITE(logout,*) "Width for each DomB segment: ", epspre*epsinc&
          &*segeps_domB
  END IF
  
  ALLOCATE (dum_aid(ntotatoms), stat = AllocateStatus)
  IF(AllocateStatus /=0 ) STOP "Allocation dum_aid failed..."
  ALLOCATE (dum_typ(ntotatoms), stat = AllocateStatus)
  IF(AllocateStatus /=0 ) STOP "Allocation dum_typ failed..."

  !Identify the domain ID- 1(A), 2(B)
  ALLOCATE (seg_dtype(ntotatoms), stat = AllocateStatus)
  IF(AllocateStatus /=0 ) STOP "Allocation seg_dtype failed..."
  CALL ASSIGN_DOMAINID(tval)
  
  DO p = 1,nmax_layers

     dum_aid = -1 !initialize everything to -1
     dum_typ = -1 !initialize everything to -1
     seg_typcnt = 0 !initialize all counts to 0     
     eps1 = epsinit !takes monomers between eps1*segeps and eps2*segeps
     eps2 = eps1 + epspre*epsinc !prefactor is a fraction of epsinc
     
     ! Calculate upto 0.7 of domain to account for wall effects
     IF(eps2*segeps_domA > 0.7*domwidth_domA) THEN 
        
        WRITE(logout,*) "Domain between ", eps1*segeps_domA, "and",&
             & eps2*segeps_domA, "is greater than 70% of domain size"
        WRITE(logout,*) "Domain more than 70% of DOMA domain size"
        PRINT*,  "Domain width more than 70% of domain size at", tval
        STOP
        
     END IF
     
     ! Calculate upto 0.7 of domain to account for wall effects
     IF(eps2*segeps_domB > 0.7*domwidth_domB) THEN 
        
        WRITE(logout,*) "Domain between ", eps1*segeps_domB, "and",&
             & eps2*segeps_domB, "is greater than domB half domain siz&
             &e"
        WRITE(logout,*) "Domain more than domB half domain size"
        PRINT*,  "Domain width more than half domain size at", tval
        STOP
        
     END IF
     
     IF(tval == 1) THEN
        WRITE(epsnum,'(I3)') int(epsinit)
        dum_fname  = 'segcoord_'//trim(adjustl(traj_fname))&
             &//"_"//trim(adjustl(epsnum))
        
        OPEN(unit = dumwrite,file = dum_fname, status="replace",&
             & action="write")
        
     END IF
     ! For both DomB and DomA check respective box boundaries
     ! None of the segment boundary can be outside the box for 
     ! p p f boundary conditions in LAMMPS. So issue an error

     ! interf + right_boundary < z_box
     IF(((interpos(1) + eps1*segeps_domB) > boxval) .OR.&
          & ((interpos(1) + eps2*segeps_domB) > boxval)) THEN
        
        PRINT *, "ERROR: Right boundary exceeds box limits..."
        PRINT *, interpos(1), eps1*segeps_domB, eps2*segeps_domB,&
             & boxval
        EXIT
        
     ! interf - less_boundary > 0.0  
     ELSEIF(((interpos(1) - eps1*segeps_domA) < 0.0) .OR. &
          & ((interpos(1) - eps2*segeps_domA) < 0.0)) THEN
        
        PRINT *, "ERROR: Left boundary exceeds box limits..."
        PRINT *, interpos(1), eps1*segeps_domA, eps2*segeps_domA,&
             & boxval
        STOP 
        
     END IF
     
     ! Define boundaries for the sublayer
     ! Right of interface (z_out > z_in)
     interin_domB  = interpos(1) + eps1*segeps_domB
     interout_domB = interpos(1) + eps2*segeps_domB

     ! Left of interface (z_out < z_in)
     interin_domA  = interpos(1) - eps1*segeps_domA
     interout_domA = interpos(1) - eps2*segeps_domA
     
     !Initialize variables
     k = 0
     
     IF(tval == 1) WRITE(logout,*) "Box dimension: ", boxval

     !Loop through all atoms to decide the atoms in the sublayer
     DO i = 1, ntotatoms
        
        a1id = aidvals(i,1)
        IF(.NOT. ANY(all_layergrp_typarr(a1id,:) == 1)) THEN
           CYCLE
        END IF

        IF(major_axis == 1) THEN
           par_pos = trx_lmp(a1id,tval) - boxval*floor(trx_lmp(a1id&
                &,tval)/boxval)
        ELSE IF(major_axis == 2) THEN 
           par_pos = try_lmp(a1id,tval) - boxval*floor(try_lmp(a1id&
                &,tval)/boxval)
        ELSE
           par_pos = trz_lmp(a1id,tval) - boxval*floor(trz_lmp(a1id&
                &,tval)/boxval)
        END IF
        
        flagbin = 0

        ! Add to bin after checking either side of the interface
        IF((par_pos .GE. interin_domB .AND. par_pos .LT.&
             & interout_domB) .OR. (par_pos .GE. interout_domA&
             & .AND. par_pos .LT. interin_domA)) THEN
           
           flagbin = flagbin + 1
           k = k + 1
           dum_aid(k) = a1id
           dum_typ(k) = aidvals(a1id,3)
           seg_typcnt(aidvals(a1id,3)) = seg_typcnt(aidvals(a1id,3))+1
           j = maxinterf
           
           IF(tval == 1) WRITE(dumwrite,'(2(I0,1X),3(F14.8,1X))') i,&
                & aidvals(a1id,3),trx_lmp(a1id,tval)&
                &,try_lmp(a1id,tval),trz_lmp(a1id,tval)

           
        END IF
     
     END DO
  
     IF(tval == 1) THEN
        WRITE(logout,*) "Width considered for domA: ",(eps2-eps1)&
             &*segeps_domA
        WRITE(logout,*) "Width considered for domB: ",(eps2-eps1)&
             &*segeps_domB
        WRITE(logout,*) "Population of each type"
        DO i = 1,ntotatomtypes
           WRITE(logout,*) i,seg_typcnt(i)
        END DO
     END IF
          
     num_mons_per_layer = k
    
     ALLOCATE (seg_aid(num_mons_per_layer), stat = AllocateStatus)
     IF(AllocateStatus /=0 ) STOP "*** Allocation seg_aid not proper ***"
     
     ALLOCATE (seg_typ(num_mons_per_layer), stat = AllocateStatus)
     IF(AllocateStatus /=0 ) STOP "*** Allocation seg_typ not proper ***"

     IF(SUM(seg_typcnt) .NE. num_mons_per_layer) THEN
        PRINT *, "Total number of atoms in a layer does not match..."
        PRINT *, num_mons_per_layer
        DO i = 1,ntotatomtypes
           WRITE(logout,*) i,seg_typcnt(i)
        END DO
        STOP
        
     END IF
     
     DO i = 1,num_mons_per_layer
        
        IF(dum_aid(i) == -1 .OR. dum_typ(i) == -1) THEN
           
           PRINT *, "Dummy arrays for seg_mons updated incorrectly"
           
           STOP
           
        END IF
        
        seg_aid(i) = dum_aid(i)
        seg_typ(i) = dum_typ(i)
        
     END DO
     
     CALL LAYERWISE_ANALYSIS(tval,p,num_mons_per_layer)
     
     DEALLOCATE(seg_aid)
     DEALLOCATE(seg_typ)
     
     epsinit = epsinit + epsinc
     
  END DO
  
  IF(tval == 1) PRINT *, "Total number of sub-layers: ", p-1
  IF(p-1 .NE. nmax_layers) PRINT *, "Warning: nmax_layers exceed avail&
       &able sub-layers", nmax_layers, p-1
  DEALLOCATE(dum_aid)
  DEALLOCATE(dum_typ)
  DEALLOCATE(seg_dtype)
  
  WRITE(logout,*) "Segmental diffusion analysis complete..."
  WRITE(logout,*) "Segmental mode relaxation analysis complete..."
  
END SUBROUTINE COMPARTMENTALIZE_PARTICLES

!--------------------------------------------------------------------

SUBROUTINE ASSIGN_DOMAINID(tval)

  USE STATICPARAMS
  IMPLICIT NONE

  REAL :: par_pos
  INTEGER :: flagdom, domnum,i,j
  INTEGER, INTENT(IN) :: tval

  seg_dtype = 0
  IF(tval == 1) THEN
     dum_fname  = 'domdist_'//trim(adjustl(traj_fname))
     OPEN(unit = 18,file = dum_fname, status="replace", action = "writ&
          &e")
  END IF

!$OMP PARALLEL DO PRIVATE(i,par_pos,flagdom,domnum)
  DO i = 1, ntotatoms
     
     IF(major_axis == 1) THEN
        par_pos = trx_lmp(i,tval) - boxval*floor(trx_lmp(i,tval)&
             &/boxval)
     ELSE IF(major_axis == 2) THEN 
        par_pos = try_lmp(i,tval) - boxval*floor(try_lmp(i,tval)&
             &/boxval)
     ELSE
        par_pos = trz_lmp(i,tval) - boxval*floor(trz_lmp(i,tval)&
             &/boxval)
     END IF

     flagdom = 0
     domnum  = 0

     IF(par_pos .GE. 0.0 .AND. par_pos .LT. interpos(1)) THEN
        
        flagdom = 1
        domnum  = domtyp(1)
        
     ELSEIF(par_pos .GE. interpos(1) .AND. par_pos .LE. boxval) THEN
        
        flagdom = 1
        domnum  = domtyp(2)

     END IF

     IF(flagdom == 0 .OR. domnum == 0 .OR. domnum > 2) THEN
        
        PRINT *, "Particle Coordinate not assigned domain properly"
        PRINT *, "i,Coordinate,domnum,tval", i,par_pos,domnum,tval
        
        WRITE(logout,*),"Particle Coordinate not assinged domain properl&
             &y"
        WRITE(logout,*), "i,Coordinate,domnum,time", i,par_pos,domnum&
             &,tval
     
        STOP
        
     ELSE
        
        seg_dtype(i) = domnum
        
     END IF
     
  END DO
!$OMP END PARALLEL DO

  IF(tval == 1) THEN
     DO i = 1,ntotatoms
        WRITE(18,'(2(I0,1X),3(F14.8,1X))') i, seg_dtype(i)&
             &,trx_lmp(i,tval), try_lmp(i,tval), trz_lmp(i,tval)
        
     END DO

     CLOSE(18)
  END IF
  
END SUBROUTINE ASSIGN_DOMAINID

!--------------------------------------------------------------------

SUBROUTINE LAYERWISE_ANALYSIS(tval,ipos,num_mons)

  USE STATICPARAMS

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: tval,ipos,num_mons
  INTEGER :: t1,t2,clock_rate,clock_max

  IF(rdf2dflag .AND. tval == 1) THEN
     IF(ipos == 1) THEN
        rdf2darray = 0.0; rdf2dvolavg = 0.0
     END IF
     CALL SYSTEM_CLOCK(t1,clock_rate,clock_max)
     CALL RDF2D_LAYER(tval,ipos,num_mons)
     CALL SYSTEM_CLOCK(t2,clock_rate,clock_max)
     PRINT *, 'Elapsed real time for 2D RDF= ',REAL(t2-t1)/&
          & REAL(clock_rate)           
  ELSEIF(rdf2dflag .AND. mod(tval,rdffreq) == 0) THEN
     CALL RDF2D_LAYER(tval,ipos,num_mons)
  END IF

END SUBROUTINE LAYERWISE_ANALYSIS

!--------------------------------------------------------------------

SUBROUTINE RDF2D_LAYER(tval,ipos,num_mons)

  USE STATICPARAMS
  IMPLICIT NONE
  
  INTEGER :: i,j,a1type,a2type,a1id,a2id,ibin
  INTEGER :: ref1type, ref2type, paircnt, ref1cnt
  REAL :: rxval,ryval,rzval,rval,rboxval,normfac,rdf2dbinval
  INTEGER, INTENT(IN) :: tval, ipos, num_mons
  INTEGER,DIMENSION(0:rdf2dmaxbin-1,npairs_2drdf) :: dumrdfarray

  rvolval = box_xl*box_yl*box_zl
  IF(ipos == 1) rdf2dvolavg = rdf2dvolavg + rvolval
  rdf2dbinval = 0.5*REAL(MIN(box_xl,box_yl,box_zl))/REAL(rdf2dmaxbin)
  IF(ipos == 1) rdf2dbinavg = rdf2dbinavg + rdf2dbinval
  
  dumrdfarray = 0
  
!$OMP PARALLEL 
!$OMP DO PRIVATE(paircnt,i,j,ref1type,ref2type,a1type,a2type,&
!$OMP& a1id,a2id,rval,rxval,ryval,rzval,ibin) REDUCTION(+:dumrdfarray)

  DO paircnt = 1,npairs_2drdf

     ref1type = pairs_2drdf_arr(paircnt,1)
     ref2type = pairs_2drdf_arr(paircnt,2)

     IF(.NOT. ANY(seg_typ == ref1type)) CYCLE
     IF(.NOT. ANY(seg_typ == ref2type)) CYCLE
     
     DO i = 1,num_mons

        a1id   = seg_aid(i)     
        a1type = seg_typ(i)

        IF(a1type .NE. aidvals(a1id,3)) THEN
           PRINT *, "1", i, a1type, a1id, aidvals(a1id,3)
           STOP "ERROR: 2D-RDF Invalid i type found in 2DRDF"
        END IF

        IF(a1type .NE. ref1type) CYCLE
        
        DO j = 1,num_mons
        
           a2id   = seg_aid(j)
           a2type = seg_aid(j)

           IF(a2type .NE. aidvals(a2id,3)) THEN
              PRINT *, "2",j, a2type, a2id, aidvals(a2id,3)
              STOP "ERROR: 2D-RDF Invalid j type found in 2DRDF"
           END IF

           IF(a2type .NE. ref2type) CYCLE
           IF(a1id == a2id) CYCLE
           
           rxval = trx_lmp(a1id,tval) - trx_lmp(a2id,tval) 
           ryval = try_lmp(a1id,tval) - try_lmp(a2id,tval) 
           rzval = trz_lmp(a1id,tval) - trz_lmp(a2id,tval) 
           
           rxval = rxval - box_xl*ANINT(rxval/box_xl)
           ryval = ryval - box_yl*ANINT(ryval/box_yl)
           rzval = rzval - box_zl*ANINT(rzval/box_zl)
        
           rval = sqrt(rxval**2 + ryval**2 + rzval**2)
           ibin = FLOOR(rval/rdf2dbinval)
        
           IF(ibin .LT. rdf2dmaxbin) THEN
           
              dumrdfarray(ibin,paircnt) = dumrdfarray(ibin,paircnt) + 1
           
           END IF
           
        END DO
        
     END DO
     
  END DO
!$OMP END DO


!$OMP DO PRIVATE(i,j,ref1type,ref2type,ref1cnt,normfac)
  DO j = 1,npairs_2drdf

     ref1type = pairs_2drdf_arr(paircnt,1)
     ref2type = pairs_2drdf_arr(paircnt,2)
     normfac  = rvolval/(REAL(ref1cnt)*REAL(seg_typcnt(ref2type)))
     
     DO i = 0,rdf2dmaxbin-1

        rdf2darray(i,j,ipos) = rdf2darray(i,j,ipos) + REAL(dumrdfarray(i&
             &,j))*normfac

     END DO
     
  END DO
!$OMP END DO

!$OMP END PARALLEL

!!$  print *, rdf2darray(:,3,ipos), LiPcntr;pause;
!  PRINT *, tval,ipos,LiPcntr, LiOcntr, POcntr

END SUBROUTINE RDF2D_LAYER

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_RDF(iframe)

  USE STATICPARAMS
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iframe
  INTEGER :: i,j,a1type,a2type,ibin,a1id,a2id,paircnt,AllocateStatus
  REAL :: rxval,ryval,rzval,rval
  INTEGER :: a1ref,a2ref
  INTEGER,ALLOCATABLE, DIMENSION(:,:) :: dumrdfarray

  rvolval = box_xl*box_yl*box_zl
  rvolavg = rvolavg + rvolval 

  ALLOCATE(dumrdfarray(0:rmaxbin-1,npairs),stat=AllocateStatus)
  IF(AllocateStatus/=0) STOP "dumrdfarray not allocated"
  dumrdfarray = 0


!$OMP PARALLEL 
!$OMP DO PRIVATE(i,j,a1type,a2type,a1id,a2id,rval,rxval,ryval,rzval&
!$OMP& ,ibin,paircnt,a1ref,a2ref) REDUCTION(+:dumrdfarray)
  DO paircnt = 1,npairs

     a1ref = pairs_rdf(paircnt,1); a2ref = pairs_rdf(paircnt,2)
     
     DO i = 1,ntotatoms
        
        a1id   = aidvals(i,1)     
        a1type = aidvals(i,3)
        
        DO j = 1,ntotatoms

           a2id   = aidvals(j,1)        
           a2type = aidvals(j,3)

           ! Remove identical IDs when computing g_AA(r)
           IF(a1id == a2id .AND. a1ref == a2ref) CYCLE


           IF(a1type == a1ref .AND. a2type == a2ref) THEN

              rxval = rxyz_lmp(a1id,1) - rxyz_lmp(a2id,1) 
              ryval = rxyz_lmp(a1id,2) - rxyz_lmp(a2id,2) 
              rzval = rxyz_lmp(a1id,3) - rxyz_lmp(a2id,3) 
              
              rxval = rxval - box_xl*ANINT(rxval/box_xl)
              ryval = ryval - box_yl*ANINT(ryval/box_yl)
              rzval = rzval - box_zl*ANINT(rzval/box_zl)
           
              rval = sqrt(rxval**2 + ryval**2 + rzval**2)
              ibin = FLOOR(rval/rbinval)
           
              IF(ibin .LT. rmaxbin) THEN
                 
                 dumrdfarray(ibin,paircnt) = dumrdfarray(ibin&
                      &,paircnt) + 1
             
              END IF

           END IF
           
        END DO

     END DO

  END DO
!$OMP END DO

!$OMP DO PRIVATE(i,j)
  DO j = 1,npairs
     
     DO i = 0,rmaxbin-1

        rdfarray(i,j) = rdfarray(i,j) + REAL(dumrdfarray(i,j))&
             &*rvolval/(REAL(pairs_rdf(j,3)))
        
     END DO
     
  END DO
!$OMP END DO

!$OMP END PARALLEL

  DEALLOCATE(dumrdfarray)

END SUBROUTINE COMPUTE_RDF

!--------------------------------------------------------------------
!!$
!!$SUBROUTINE COMPUTE_RADGYR(iframe)
!!$
!!$  USE STATICPARAMS
!!$  IMPLICIT NONE
!!$
!!$  INTEGER :: i,j,molid,atype
!!$  REAL    :: rgsqavg, rgxxavg, rgyyavg, rgzzavg
!!$  REAL, DIMENSION(1:nchains) :: rgxx, rgyy, rgzz, rgsq
!!$  REAL, DIMENSION(1:nchains) :: rxcm, rycm, rzcm, totmass
!!$  INTEGER, INTENT(IN) :: iframe
!!$
!!$  rgxx = 0.0; rgyy =0.0; rgzz = 0.0; rgsq = 0.0
!!$  rgsqavg = 0.0; rgxxavg = 0.0; rgyyavg = 0.0; rgzzavg = 0.0
!!$  rxcm = 0.0; rycm = 0.0; rzcm = 0.0
!!$  totmass = 0.0
!!$
!!$  IF(iframe == 1) THEN
!!$     PRINT *, masses
!!$  END IF
!!$
!!$!$OMP PARALLEL 
!!$!$OMP DO PRIVATE(i,molid,atype) REDUCTION(+:totmass,rxcm,rycm,rzcm)
!!$
!!$  DO i = 1,ntotatoms
!!$
!!$     molid = aidvals(i,2)
!!$     atype = aidvals(i,3)
!!$
!!$     IF(ANY(polytyp_arr == atype)) THEN
!!$       
!!$        totmass(molid) = totmass(molid) + masses(atype,1)
!!$        rxcm(molid) = rxcm(molid)+ rxyz_lmp(i,1)*masses(atype,1)
!!$        rycm(molid) = rycm(molid)+ rxyz_lmp(i,2)*masses(atype,1)
!!$        rzcm(molid) = rzcm(molid)+ rxyz_lmp(i,3)*masses(atype,1)
!!$
!!$     END IF
!!$
!!$  END DO
!!$!$OMP END DO
!!$
!!$!$OMP DO PRIVATE(i)
!!$  DO i = 1, nchains
!!$     rxcm(i) = rxcm(i)/totmass(i)
!!$     rycm(i) = rycm(i)/totmass(i)
!!$     rzcm(i) = rzcm(i)/totmass(i)
!!$  END DO
!!$!$OMP END DO
!!$
!!$
!!$!$OMP DO PRIVATE(i,molid,atype) REDUCTION(+:rgxx, rgyy, rgzz, rgsq)
!!$  DO i = 1,ntotatoms
!!$
!!$     molid = aidvals(i,2)
!!$     atype = aidvals(i,3)
!!$
!!$     IF(ANY(polytyp_arr == atype)) THEN
!!$
!!$        rgxx(molid) = rgxx(molid) + masses(atype,1)*((rxyz_lmp(i,1)&
!!$             &-rxcm(molid))**2)
!!$        rgyy(molid) = rgyy(molid) + masses(atype,1)*((rxyz_lmp(i,2)&
!!$             &-rycm(molid))**2)
!!$        rgzz(molid) = rgzz(molid) + masses(atype,1)*((rxyz_lmp(i,3)&
!!$             &-rzcm(molid))**2)
!!$
!!$        rgsq(molid) = rgsq(molid) + masses(atype,1)*((rxyz_lmp(i,1)&
!!$             &-rxcm(molid))**2 + (rxyz_lmp(i,2)-rycm(molid))**2 +&
!!$             & (rxyz_lmp(i,3)-rzcm(molid))**2)
!!$
!!$
!!$     END IF
!!$
!!$  END DO
!!$!$OMP END DO
!!$
!!$!$OMP DO PRIVATE(i)
!!$  DO i = 1, nchains
!!$     rgxx(i) = rgxx(i)/totmass(i)
!!$     rgyy(i) = rgyy(i)/totmass(i)
!!$     rgzz(i) = rgzz(i)/totmass(i)
!!$     rgsq(i) = rgsq(i)/totmass(i)
!!$  END DO
!!$!$OMP END DO
!!$!$OMP END PARALLEL
!!$
!!$  IF(iframe == 1) THEN
!!$
!!$     OPEN(unit = 98,file ="totmasschk.txt",action="write",status="repl&
!!$          &ace")
!!$
!!$     DO i = 1,nchains
!!$
!!$        WRITE(98,'(I0,1X,8(F14.5,1X))') i, totmass(i),rxcm(i),&
!!$             & rycm(i), rzcm(i), rgxx(i), rgyy(i), rgzz(i), rgsq(i)
!!$
!!$     END DO
!!$
!!$     CLOSE(98)
!!$
!!$     OPEN(unit = 98,file ="molidchk.txt",action="write",status="repl&
!!$          &ace")
!!$
!!$     DO i = 1,ntotatoms
!!$
!!$        WRITE(98,'(I0,1X,I0)') i, aidvals(i,2)
!!$
!!$     END DO
!!$
!!$     CLOSE(98)
!!$
!!$  END IF
!!$
!!$
!!$  DO i = 1,nchains
!!$  
!!$     rgsqavg = rgsqavg + rgsq(i)
!!$     rgxxavg = rgxxavg + rgxx(i)
!!$     rgyyavg = rgyyavg + rgyy(i)
!!$     rgzzavg = rgzzavg + rgzz(i)
!!$     
!!$  END DO
!!$  
!!$  rgsqavg = rgsqavg/REAL(nchains)
!!$  rgxxavg = rgxxavg/REAL(nchains)
!!$  rgyyavg = rgyyavg/REAL(nchains)
!!$  rgzzavg = rgzzavg/REAL(nchains)
!!$  
!!$
!!$  IF(rgavg) THEN
!!$
!!$     WRITE(rgavgwrite,'(I0,1X,4(F14.6,1X))') timestep, sqrt(rgsqavg),&
!!$          & sqrt(rgxxavg), sqrt(rgyyavg), sqrt(rgzzavg)
!!$
!!$  END IF
!!$
!!$  IF(rgall) THEN
!!$     
!!$     WRITE(rgwrite,'(2(I0,1X),4(A3,1X))') timestep, nchains,"rgt","rg&
!!$          &x","rgy","rgz"
!!$
!!$     DO i = 1,nchains
!!$     
!!$        WRITE(rgwrite,'(I0,1X,4(F14.5,1X))') i,sqrt(rgsq(i))&
!!$             &,sqrt(rgxx(i)),sqrt(rgyy(i)),sqrt(rgzz(i))
!!$
!!$     END DO
!!$
!!$  END IF
!!$     
!!$END SUBROUTINE COMPUTE_RADGYR

!--------------------------------------------------------------------

SUBROUTINE CLUSTER_ANALYSIS(frnum)

  USE STATICPARAMS

  IMPLICIT NONE

!Ref Sevick et.al ., J Chem Phys 88 (2)

  INTEGER :: i,j,k,a2ptr,a1id,a2id,itype,jtype,jptr,idum,jflag,jcnt&
       &,iflag,jtot,jind,jprev
  INTEGER, DIMENSION(ntotion_centers,ntotion_centers) :: all_direct&
       &,catan_direct,all_neigh,catan_neigh
  INTEGER, DIMENSION(1:ntotion_centers) :: union_all,flag_catan,scnt&
       &,all_linked
  REAL :: rxval, ryval, rzval, rval
  INTEGER, INTENT(IN) :: frnum

!$OMP PARALLEL SHARED(catan_direct,catan_neigh)

!$OMP DO PRIVATE(i,j)
  DO i = 1,ntotion_centers

     scnt(i) = 0; all_linked(i)  = 0
     union_all(i) = -1; flag_catan(i) = -1

     DO j = 1,ntotion_centers

        IF(i .NE. j) THEN
           all_direct(i,j) = 0
           catan_direct(i,j) = 0
        END IF

        IF(i == j) THEN
           all_direct(i,j) = 1
           catan_direct(i,j) = 0
        END IF

        all_neigh(i,j) = 0
!        catan_neigh(i,j) = 0

     END DO

  END DO
!$OMP END DO

!Create Direct connectivity matrix
!allneigh - does not distinguish between Li and P neigh
!catan_neigh - neighbors with sequence cat-an-cat-an.. or an-cat-an-cat...

!$OMP DO PRIVATE(i,j,a1id,a2ptr,a2id,rxval,ryval,rzval,rval,itype&
!$OMP& ,jptr,jtype)  
  DO i = 1,ntotion_centers
     
     a1id = allionids(i,1)
     a2ptr = 1
     itype = aidvals(a1id,3)
     jptr  = 1
     all_neigh(i,i) = a1id
     
     DO j = 1,ntotion_centers
        
        a2id = allionids(j,1)
        
        rxval = rxyz_lmp(a1id,1) - rxyz_lmp(a2id,1) 
        ryval = rxyz_lmp(a1id,2) - rxyz_lmp(a2id,2) 
        rzval = rxyz_lmp(a1id,3) - rxyz_lmp(a2id,3) 
        
        rxval = rxval - box_xl*ANINT(rxval/box_xl)
        ryval = ryval - box_yl*ANINT(ryval/box_yl)
        rzval = rzval - box_zl*ANINT(rzval/box_zl)
        
        rval = sqrt(rxval**2 + ryval**2 + rzval**2)

        IF(rval .LT. rneigh_cut .AND. a1id .NE. a2id) THEN
           
           all_direct(i,j) = 1
           all_neigh(i,j)  = a2id
!           all_neigh(i,a2ptr) = a2id
!           a2ptr = a2ptr + 1
           
           jtype = aidvals(a2id,3)
           
           IF(itype .NE. jtype) THEN
              
              catan_direct(i,j) = 1
!              catan_neigh(i,j)  = a2id
!              catan_neigh(i,jptr+1) = a2id
              itype = jtype
!              jptr  = jptr + 1
              
           END IF
           
        END IF

     END DO
     
  END DO
!$OMP END DO  

!Check for symmetry
  IF(frnum == 1) THEN
!$OMP DO
     DO i = 1,ntotion_centers

        DO j = 1,ntotion_centers

           IF(all_direct(i,j) .NE. all_direct(j,i)) STOP "Unsymmetric&
                & all_direct"

          IF(all_neigh(i,j) .NE. 0) THEN

              IF(all_neigh(i,j) .NE. all_neigh(j,j) .OR. all_neigh(j&
                   &,i) .NE. all_neigh(i,i)) THEN

                 PRINT *, i,j,all_direct(i,j),all_direct(j,i)&
                      &,all_neigh(j,i),all_neigh(i,i)
                 STOP "Unsymmetric neighbor list"
                 
              END IF
              
           END IF

        END DO

     END DO
!$OMP END DO
  END IF

!$OMP END PARALLEL        


!Intersection

  DO i = 1,ntotion_centers-1 !Ref row counter

     iflag = 0
     idum  = i

     DO WHILE(iflag == 0 .AND. union_all(i) == -1)
     
        jflag = 0
        k    = 1 !Column counter
        j    = idum+1 !Other row counter
        
        DO WHILE(jflag == 0 .AND. k .LE. ntotion_centers) 
           
           IF((all_direct(i,k) == all_direct(j,k)).AND. all_direct(i&
                &,k)== 1) THEN
              
              jflag = 1
!!$              jprev = 0

              DO jcnt = 1,ntotion_centers
               
                 
!!$                 IF(all_direct(j,jcnt) == 1) jprev = 1

                 !Replace highest row by union of two rows

                 all_direct(j,jcnt) = all_direct(i,jcnt) .OR.&
                      & all_direct(j,jcnt)

!!$                 IF((all_direct(j,jcnt) == 1 .AND. all_direct(i,jcnt)&
!!$                      &==1) .AND. jprev == 0) THEN
!!$                    
!!$                    all_neigh(j,jcnt) = all_neigh(i,jcnt)
!!$                    jprev = 0 !Other condition is already
!!$                    ! incorporated before
!!$                 END IF
!!$                 
              END DO
              
              union_all(i) = 1 !One match implies the low ranked row
              ! is present in high ranked row
              
           ELSE
              
              k = k + 1
              
           END IF
           
        END DO
        
        IF(union_all(i) == 1) THEN
           
           iflag = 1
           
        ELSE
           
           idum  = idum + 1
           
        END IF
        
        IF(idum == ntotion_centers) iflag = 1

     END DO

  END DO
  
!Count
  jtot = 0
!$OMP PARALLEL PRIVATE(i,j,jind) 
!$OMP DO
  DO i = 1,ntotion_centers

     IF(union_all(i) == -1) THEN
        
        jind = 0

        DO j = 1,ntotion_centers

           IF(all_direct(i,j) == 1) jind = jind + 1

        END DO

        scnt(jind) = scnt(jind) + 1
        all_linked(i) = jind

     END IF
     
  END DO
!$OMP END DO

!$OMP DO

  DO i = 1,ntotion_centers

     clust_avg(i) = clust_avg(i) + scnt(i)

  END DO
!$OMP END DO

!$OMP END PARALLEL

  IF(frnum == 1) THEN
     OPEN(unit =90,file ="scnt.txt",action="write",status="replace")  
  END IF

  jtot = 0
  
  DO i = 1,ntotion_centers
     
     IF(frnum == 1) WRITE(90,*) i,scnt(i)
     jtot = jtot + all_linked(i)
     
  END DO
  
  IF(jtot .NE. ntotion_centers) THEN
     
     PRINT *, "Sum of ions not equal to total ion centers"
     PRINT *, jtot, ntotion_centers
     STOP
     
  END IF

  IF(frnum == 1) CLOSE(90)
  
  IF(frnum == 1) THEN

     OPEN(unit =90,file ="all_neigh.txt",action="write",status="replace")  
    
     DO i = 1,ntotion_centers
        
        IF(union_all(i) == -1) THEN
           
           WRITE(90,*) i,all_linked(i)
        
           DO j = 1,ntotion_centers
           
              IF(all_direct(i,j) == 1) WRITE(90,*) j,allionids(j,1),&
                   & allionids(j,2),all_direct(i,j)
           
           END DO

        END IF
        
     END DO

     CLOSE(90)

  END IF
  
END SUBROUTINE CLUSTER_ANALYSIS

!--------------------------------------------------------------------

SUBROUTINE SORT_POLY_FREE_BOUND_COMPLEX(tval)

  USE STATICPARAMS

  IMPLICIT NONE

  INTEGER :: i,j,a1id,a2id,boundflag,pcnt,catcnt,tid,cnt
  INTEGER :: pfree,pbound,pboundtot,pfreetot,AllocateStatus
  REAL :: rxval,ryval,rzval,rval
  INTEGER, DIMENSION(1:p_ioncnt,0:nproc-1) :: polbounddum,polfreedum
  INTEGER, INTENT(IN) :: tval

 pfree = 0; pbound = 0; pboundtot=0; pfreetot=0;catcnt = 1

!$OMP PARALLEL
!$OMP DO PRIVATE(i,j)
  DO i = 1,p_ioncnt
     
     DO j = 0,nproc-1
        
        polbounddum(i,j) = -1
        polfreedum(i,j)  = -1

     END DO

  END DO
!$OMP END DO

!$OMP DO PRIVATE(pcnt,tid,a1id,boundflag,a2id,rxval,ryval,rzval,rval)& 
!$OMP& REDUCTION(+:pboundtot,pfreetot) FIRSTPRIVATE(pfree,pbound,catcnt)
  DO pcnt = 1,p_ioncnt

     tid  = OMP_GET_THREAD_NUM()
     a1id = polyionarray(pcnt,1)
     catcnt = 1
     boundflag = 0
     
     DO WHILE(catcnt .LE. ioncnt)
     
        a2id  = ionarray(catcnt,1)

        rxval = rxyz_lmp(a1id,1) - rxyz_lmp(a2id,1) 
        ryval = rxyz_lmp(a1id,2) - rxyz_lmp(a2id,2) 
        rzval = rxyz_lmp(a1id,3) - rxyz_lmp(a2id,3) 
        
        rxval = rxval - box_xl*ANINT(rxval/box_xl)
        ryval = ryval - box_yl*ANINT(ryval/box_yl)
        rzval = rzval - box_zl*ANINT(rzval/box_zl)
        
        rval = sqrt(rxval**2 + ryval**2 + rzval**2)
        
        IF(rval .LT. rcatpol_cut1) THEN

           pbound = pbound + 1
           pboundtot = pboundtot + 1
           polbounddum(pbound,tid) = a1id
           catcnt = ioncnt + 1
           boundflag = 1

        ELSE
           
           catcnt = catcnt + 1

        END IF

     END DO
           
     IF(boundflag == 0) THEN

        pfree = pfree + 1
        pfreetot = pfreetot + 1
        polfreedum(pfree,tid) = a1id

     END IF

  END DO
!$OMP END DO
!$OMP END PARALLEL

  IF(pboundtot == 0 .AND. pfreetot == 0) THEN

     PRINT *, "Zero pboundtot/pfreetot"
     PRINT *, pboundtot, pfreetot, tval
     STOP

  END IF

  ALLOCATE(polboundarr(1:pboundtot),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate polboundarr"

  ALLOCATE(polfreearr(1:pfreetot),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate polfreearr"
  
  pfree = 0; pbound = 0
 
  DO i = 0,nproc-1
     
     cnt = 0

     DO WHILE(polfreedum(cnt+1,i) .NE. -1)

        cnt = cnt + 1
        pfree = pfree + 1
        polfreearr(pfree) = polfreedum(cnt,i)

     END DO
 
     cnt = 0

     DO WHILE(polbounddum(cnt+1,i) .NE. -1)
        
        cnt = cnt + 1
        pbound = pbound + 1
        polboundarr(pbound) = polbounddum(cnt,i)

     END DO

  END DO

  IF(pfree .NE. pfreetot .OR. pbound .NE. pboundtot) THEN

     PRINT *, "Unequal assignment in free and bound oxygens"
     PRINT *, pfree,pfreetot,pbound,pboundtot
     STOP

  END IF

  CALL FREE_BOUND_POLRDF(pfreetot,pboundtot)

  DEALLOCATE(polfreearr)
  DEALLOCATE(polboundarr)

END SUBROUTINE SORT_POLY_FREE_BOUND_COMPLEX

!--------------------------------------------------------------------

SUBROUTINE FREE_BOUND_POLRDF(pfreetot,pboundtot)

  USE STATICPARAMS

  IMPLICIT NONE

  INTEGER :: i,j,ptot,a1id,a2id,ibin
  INTEGER :: pbound_free, pbound_bound, pfree_free
  REAL    :: rval,rxval,ryval,rzval,rboxval
  INTEGER,INTENT(IN) :: pfreetot,pboundtot
  INTEGER,DIMENSION(0:rmaxbin-1) ::dumrdf_p_ff,dumrdf_p_bb,dumrdf_p_fb

  rvolval = box_xl*box_yl*box_zl
  IF(rdfflag == .false.) THEN !Already in RDF. So if it is true it is
     !already accounted in rdf computation
  
     rvolavg = rvolavg + rvolval
  
  END IF

  ptot = pfreetot + pboundtot

  pbound_free = 0; pbound_bound = 0; pfree_free = 0

!$OMP PARALLEL 

!$OMP DO PRIVATE(i)
  DO i = 0,rmaxbin-1
     
     dumrdf_p_fb(i) = 0
     dumrdf_p_ff(i) = 0
     dumrdf_p_bb(i) = 0
     
  END DO
!$OMP END DO

!$OMP DO REDUCTION(+:dumrdf_p_fb,pbound_free,pfree_free,pbound_bound)&
!$OMP& PRIVATE(i,j,a1id,a2id,rval,rxval,ryval,rzval,ibin)
  DO i = 1,pfreetot

     a1id = polfreearr(i)
     
     IF(aidvals(a1id,3) .NE. p_iontype) THEN
        PRINT *, a1id, aidvals(a1id,3), p_iontype
        STOP "Wrong polyion for a1id"
     END IF

     DO j = i,pboundtot

        a2id   = polboundarr(j)

        IF(aidvals(a2id,3) .NE. p_iontype) THEN
           PRINT *, a2id, aidvals(a2id,3), p_iontype
           STOP "Wrong polyion for a2id"
        END IF

        rxval = rxyz_lmp(a1id,1) - rxyz_lmp(a2id,1) 
        ryval = rxyz_lmp(a1id,2) - rxyz_lmp(a2id,2) 
        rzval = rxyz_lmp(a1id,3) - rxyz_lmp(a2id,3) 
        
        rxval = rxval - box_xl*ANINT(rxval/box_xl)
        ryval = ryval - box_yl*ANINT(ryval/box_yl)
        rzval = rzval - box_zl*ANINT(rzval/box_zl)
        
        rval = sqrt(rxval**2 + ryval**2 + rzval**2)
        ibin = FLOOR(rval/rbinval)

        IF(ibin .LT. rmaxbin) THEN

           dumrdf_p_fb(ibin) = dumrdf_p_fb(ibin) + 2
           pbound_free = pbound_free  + 1

        END IF

     END DO

  END DO
!$OMP END DO

!$OMP DO REDUCTION(+:dumrdf_p_ff,pbound_free,pfree_free,pbound_bound)&
!$OMP& PRIVATE(i,j,a1id,a2id,rval,rxval,ryval,rzval,ibin)
  DO i = 1,pfreetot

     a1id = polfreearr(i)
     
     IF(aidvals(a1id,3) .NE. p_iontype) THEN
        PRINT *, a1id, aidvals(a1id,3), p_iontype
        STOP "Wrong polyion for a1id"
     END IF

     DO j = i+1,pfreetot

        a2id   = polfreearr(j)

        IF(aidvals(a2id,3) .NE. p_iontype) THEN
           PRINT *, a2id, aidvals(a2id,3), p_iontype
           STOP "Wrong polyion for a2id"
        END IF

        rxval = rxyz_lmp(a1id,1) - rxyz_lmp(a2id,1) 
        ryval = rxyz_lmp(a1id,2) - rxyz_lmp(a2id,2) 
        rzval = rxyz_lmp(a1id,3) - rxyz_lmp(a2id,3) 
        
        rxval = rxval - box_xl*ANINT(rxval/box_xl)
        ryval = ryval - box_yl*ANINT(ryval/box_yl)
        rzval = rzval - box_zl*ANINT(rzval/box_zl)
        
        rval = sqrt(rxval**2 + ryval**2 + rzval**2)
        ibin = FLOOR(rval/rbinval)
        
        IF(ibin .LT. rmaxbin) THEN
           
           dumrdf_p_ff(ibin) = dumrdf_p_ff(ibin) + 2
           pfree_free = pfree_free + 1
           
        END IF

     END DO

  END DO
!$OMP END DO

!$OMP DO REDUCTION(+:dumrdf_p_bb,pbound_free,pfree_free,pbound_bound) &
!$OMP& PRIVATE(i,j,a1id,a2id,rval,rxval,ryval,rzval,ibin)
  DO i = 1,pboundtot

     a1id = polboundarr(i)

     IF(aidvals(a1id,3) .NE. p_iontype) THEN
        PRINT *, a1id, aidvals(a1id,3), p_iontype
        STOP "Wrong polyion for a1id"
     END IF

     DO j = i+1,pboundtot

        a2id   = polboundarr(j)

        IF(aidvals(a2id,3) .NE. p_iontype) THEN
           PRINT *, a2id, aidvals(a2id,3), p_iontype
           STOP "Wrong polyion for a2id"
        END IF

        rxval = rxyz_lmp(a1id,1) - rxyz_lmp(a2id,1) 
        ryval = rxyz_lmp(a1id,2) - rxyz_lmp(a2id,2) 
        rzval = rxyz_lmp(a1id,3) - rxyz_lmp(a2id,3) 
        
        rxval = rxval - box_xl*ANINT(rxval/box_xl)
        ryval = ryval - box_yl*ANINT(ryval/box_yl)
        rzval = rzval - box_zl*ANINT(rzval/box_zl)
        
        rval = sqrt(rxval**2 + ryval**2 + rzval**2)
        ibin = FLOOR(rval/rbinval)
        
        IF(ibin .LT. rmaxbin) THEN
           
           dumrdf_p_bb(ibin) = dumrdf_p_bb(ibin) + 2
           pbound_bound = pbound_bound + 1
           
        END IF

     END DO

  END DO
!$OMP END DO


!$OMP DO 
  DO i = 0,rmaxbin-1

     rdf_p_fb(i) = rdf_p_fb(i) + REAL(dumrdf_p_fb(i))*rvolval&
          &/(REAL(pboundtot)*REAL(pfreetot))
     rdf_p_ff(i) = rdf_p_ff(i) + REAL(dumrdf_p_ff(i))*rvolval&
          &/(REAL(pfreetot)*REAL(pfreetot))
     rdf_p_bb(i) = rdf_p_bb(i) + REAL(dumrdf_p_bb(i))*rvolval&
          &/(REAL(pboundtot)*REAL(pboundtot))

  END DO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE FREE_BOUND_POLRDF

!--------------------------------------------------------------------

SUBROUTINE DYNAMICS_MAIN()

  USE STATICPARAMS
  IMPLICIT NONE

  IF(ion_diff) THEN
     PRINT *, "Beginning ion diffusion calculation..."
     CALL DIFF_IONS()
     PRINT *, "Finished ion diffusion calculation..."
  END IF

  IF(cion_diff) THEN
     PRINT *, "Beginning counter-ion diffusion calculation..."
     CALL DIFF_COUNTERIONS()
     PRINT *, "Finished counter-ion diffusion calculation..."
  END IF

  IF(catan_autocfflag) THEN
     PRINT *, "Beginning cat-an residence time calculation..."
     CALL RESIDENCE_TIME_CATAN()
     PRINT *, "Finished cat-an residence time calculation..."
  END IF

  IF(catpol_autocfflag) THEN
     PRINT *, "Beginning cat-pol residence time calculation..."
     CALL RESIDENCE_TIME_CATPOL()
     PRINT *, "Beginning cat-pol residence time calculation..."
  END IF

END SUBROUTINE DYNAMICS_MAIN

!--------------------------------------------------------------------

SUBROUTINE CHECK_MOMENTUM(tval)

  USE STATICPARAMS
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: tval
  INTEGER :: i,aid,atype
  REAL :: xmom, ymom, zmom
  REAL, PARAMETER :: tol_mom = 1e-5

  xmom = 0; ymom = 0; zmom = 0

  DO i = 1, ntotatoms
     
     aid = vel_xyz(i,1); atype = aidvals(aid,3)
     xmom = xmom + masses(atype,1)*vel_xyz(i,2)
     ymom = ymom + masses(atype,1)*vel_xyz(i,3)
     zmom = zmom + masses(atype,1)*vel_xyz(i,4)

  END DO

  IF( (abs(xmom) .GT. tol_mom) .OR. (abs(ymom) .GT. tol_mom) .OR.&
       & (abs(zmom) .GT. tol_mom) ) THEN

     PRINT *, "WARNING: Net momentum not zero: ", tval,xmom,ymom,zmom

  ELSE

     IF(tval == 0) THEN

        PRINT *, "Momentum conserved at the beginning: ", xmom, ymom,&
             & zmom

     END IF

  END IF

END SUBROUTINE CHECK_MOMENTUM

!--------------------------------------------------------------------

SUBROUTINE DIFF_IONS()

  USE STATICPARAMS

  IMPLICIT NONE

  INTEGER :: i,j,tinc,ifin,tim,aid,atype,ierr,jout
  REAL    :: rxcm, rycm, rzcm
  REAL, DIMENSION(0:nframes-1) :: gxarr,gyarr,gzarr

  dum_fname = "iondiff_"//trim(adjustl(traj_fname))
  OPEN(unit = dumwrite,file = dum_fname, status="replace",action&
       &="write",iostat = ierr)

  IF(ierr /= 0) STOP "Ion diffusion file not found"

  WRITE(dumwrite,'(2(I0,1X),F14.8)') ioncnt, iontype, delta_t

! Ion Diffusion Analysis

  DO i = 0,nframes-1

     gxarr(i) = 0.0
     gyarr(i) = 0.0
     gzarr(i) = 0.0

  END DO

!$OMP PARALLEL DO PRIVATE(tinc,ifin,i,tim,j,rxcm,rycm,rzcm,aid)&
!$OMP&  REDUCTION(+:gxarr,gyarr,gzarr)
  DO tinc = 0, nframes-1

     ifin = nframes - tinc

     DO i = 1,ifin

        tim = i + tinc

        DO j = 1,ioncnt

           rxcm = itrx_lmp(j,tim) - itrx_lmp(j,i)
           rycm = itry_lmp(j,tim) - itry_lmp(j,i)
           rzcm = itrz_lmp(j,tim) - itrz_lmp(j,i)

           gxarr(tinc) = gxarr(tinc) + rxcm**2
           gyarr(tinc) = gyarr(tinc) + rycm**2
           gzarr(tinc) = gzarr(tinc) + rzcm**2


        END DO

     END DO

     gxarr(tinc) = gxarr(tinc)/(ifin*ioncnt)
     gyarr(tinc) = gyarr(tinc)/(ifin*ioncnt)
     gzarr(tinc) = gzarr(tinc)/(ifin*ioncnt)


  END DO
!$OMP END PARALLEL DO

  DO i = 0, nframes-1

     WRITE(dumwrite,"(I10,1X,3(F14.5,1X))") tarr_lmp(i+1), gxarr(i)&
          &,gyarr(i), gzarr(i)

  END DO

  CLOSE(dumwrite)

END SUBROUTINE DIFF_IONS

!--------------------------------------------------------------------

SUBROUTINE DIFF_COUNTERIONS()

  USE STATICPARAMS

  IMPLICIT NONE

  INTEGER :: i,j,tinc,ifin,tim,aid,atype,ierr,jout
  REAL    :: rxcm, rycm, rzcm
  REAL, DIMENSION(0:nframes-1) :: gxarr,gyarr,gzarr

  dum_fname = "countiondiff_"//trim(adjustl(traj_fname))
  OPEN(unit = dumwrite,file = dum_fname, status="replace",action&
       &="write",iostat = ierr)

  IF(ierr /= 0) STOP "Counter-ion diffusion file not found"

  WRITE(dumwrite,'(2(I0,1X),F14.8)') c_ioncnt, c_iontype, delta_t

! Ion Diffusion Analysis

  DO i = 0,nframes-1

     gxarr(i) = 0.0
     gyarr(i) = 0.0
     gzarr(i) = 0.0

  END DO

! To do shifted time average for segmental diffusion

!$OMP PARALLEL DO PRIVATE(tinc,ifin,i,tim,j,rxcm,rycm,rzcm,aid)&
!$OMP&  REDUCTION(+:gxarr,gyarr,gzarr)
  DO tinc = 0, nframes-1

     ifin = nframes - tinc
     
     DO i = 1,ifin
        
        tim = i + tinc
      
        DO j = 1,c_ioncnt

           rxcm = ctrx_lmp(j,tim) - ctrx_lmp(j,i)
           rycm = ctry_lmp(j,tim) - ctry_lmp(j,i)
           rzcm = ctrz_lmp(j,tim) - ctrz_lmp(j,i)

           gxarr(tinc) = gxarr(tinc) + rxcm**2
           gyarr(tinc) = gyarr(tinc) + rycm**2
           gzarr(tinc) = gzarr(tinc) + rzcm**2
                      
        END DO

     END DO
     
     gxarr(tinc) = gxarr(tinc)/(REAL(ifin*c_ioncnt))
     gyarr(tinc) = gyarr(tinc)/(REAL(ifin*c_ioncnt))
     gzarr(tinc) = gzarr(tinc)/(REAL(ifin*c_ioncnt))

     
  END DO
!$OMP END PARALLEL DO

  DO i = 0, nframes-1

     WRITE(dumwrite,"(I10,1X,3(F14.5,1X))") tarr_lmp(i+1), gxarr(i)&
          &,gyarr(i),gzarr(i)

  END DO

  CLOSE(dumwrite)

END SUBROUTINE DIFF_COUNTERIONS

!--------------------------------------------------------------------
!Ref: Borodin and Smith
!Macromolecules Vol: 39, No: 4, 1620-1629, 2006
!Here we look at the residence time on the polymer-ion
SUBROUTINE RESIDENCE_TIME_CATPOL()

  USE STATICPARAMS
  IMPLICIT NONE

  INTEGER :: i,j,tval,a1id,a2id,ierr,ifin,tinc,tim
  REAL :: rxval,ryval,rzval,rval
  INTEGER,DIMENSION(1:p_ioncnt,nframes) :: autocf
  REAL,DIMENSION(0:nframes-1) :: tplot_cf

  dum_fname = "autocorrpol_"//trim(adjustl(traj_fname))
  OPEN(unit = dumwrite,file = dum_fname, status="replace",action&
       &="write",iostat = ierr)

  IF(ierr /= 0) STOP "polymer residence time file not found"

!$OMP PARALLEL
!$OMP DO PRIVATE(i,j)
  DO j = 1,nframes

     DO i = 1,p_ioncnt

        autocf(i,j) = 0

     END DO
     
     tplot_cf(j-1) = 0

  END DO
!$OMP END DO

!$OMP DO PRIVATE(tval,i,j,a1id,a2id,rxval,ryval,rzval,rval)
  DO tval = 1,nframes 

     DO i = 1,p_ioncnt !populate autocorrelation fn array

        j = 1; a1id = polyionarray(i,1)

        IF(aidvals(a1id,3) .NE. p_iontype) THEN
              
           PRINT *, "Wrong poly-iontype atom type"
           PRINT *, tval, a1id, aidvals(a1id,3), p_iontype,&
                & polyionarray(i,1)
           STOP

        END IF
           
        DO WHILE(j .LE. ioncnt)

           a2id = ionarray(j,1)

           IF(aidvals(a2id,3) .NE. iontype) THEN
              
              PRINT *, "Wrong atom type"
              PRINT *, tval, a2id, aidvals(a2id,3), iontype,&
                   & ionarray(j,1)
              STOP

           END IF
           
           rxval = ptrx_lmp(i,tval) - itrx_lmp(j,tval) 
           ryval = ptry_lmp(i,tval) - itry_lmp(j,tval) 
           rzval = ptrz_lmp(i,tval) - itrz_lmp(j,tval) 
           
           rxval = rxval - box_xl*ANINT(rxval/box_xl)
           ryval = ryval - box_yl*ANINT(ryval/box_yl)
           rzval = rzval - box_zl*ANINT(rzval/box_zl)
           
           rval = sqrt(rxval**2 + ryval**2 + rzval**2)
           
           IF(rval .LT. rcatpol_cut2) THEN
              
              autocf(i,tval) = 1
              j = ioncnt+1
              
           ELSE

              j = j +1

           END IF

        END DO

     END DO

  END DO
!$OMP END DO

!$OMP DO PRIVATE(tinc,ifin,tim,i,j)

  DO tinc = 0, nframes-1 !compute spectral product

     ifin = nframes - tinc
     
     DO i = 1,ifin
        
        tim = i + tinc
      
        DO j = 1,p_ioncnt

           tplot_cf(tinc) = tplot_cf(tinc) + REAL(autocf(j,tim)&
                &*autocf(j,i))
           
        END DO

     END DO

     tplot_cf(tinc) = tplot_cf(tinc)/REAL(ifin*p_ioncnt)
     
  END DO
!$OMP END DO
!$OMP END PARALLEL

  DO tinc = 0, nframes-1

     WRITE(dumwrite,"(I0,1X,F14.6)") tinc, tplot_cf(tinc)

  END DO

  CLOSE(dumwrite)

END SUBROUTINE RESIDENCE_TIME_CATPOL

!--------------------------------------------------------------------

!Ref: Borodin and Smith
!Macromolecules Vol: 39, No: 4, 1620-1629, 2006
!Here we look at the residence time on the ANION
SUBROUTINE RESIDENCE_TIME_CATAN()

  USE STATICPARAMS
  IMPLICIT NONE

  INTEGER :: i,j,tval,a1id,a2id,ierr,ifin,tinc,tim
  REAL :: rxval,ryval,rzval,rval
  INTEGER,DIMENSION(1:c_ioncnt,nframes) :: autocf
  REAL,DIMENSION(0:nframes-1) :: tplot_cf

  dum_fname = "autocorrcion_"//trim(adjustl(traj_fname))
  OPEN(unit = dumwrite,file = dum_fname, status="replace",action&
       &="write",iostat = ierr)

  IF(ierr /= 0) STOP "c-ion residence time file not found"

!$OMP PARALLEL
!$OMP DO PRIVATE(i,j)
  DO j = 1,nframes

     DO i = 1,c_ioncnt

        autocf(i,j) = 0

     END DO
     
     tplot_cf(j-1) = 0

  END DO
!$OMP END DO

!$OMP DO PRIVATE(tval,i,j,a1id,a2id,rxval,ryval,rzval,rval)
  DO tval = 1,nframes 

     DO i = 1,c_ioncnt !populate autocorrelation fn array

        j = 1; a1id = counterarray(i,1)
        
        IF(aidvals(a1id,3) .NE. c_iontype) THEN
           
           PRINT *, "Wrong counter-ion type"
           PRINT *, tval, a1id, aidvals(a1id,3), c_iontype,&
                & counterarray(i,1)

        END IF

        DO WHILE(j .LE. ioncnt)

           a2id = ionarray(j,1)

           IF(aidvals(a2id,3) .NE. iontype) THEN
           
              PRINT *, "Wrong ion type"
              PRINT *, tval, a2id, aidvals(a2id,3), iontype&
                   &,ionarray(j,1)
              
           END IF
           
           rxval = ctrx_lmp(i,tval) - itrx_lmp(j,tval) 
           ryval = ctry_lmp(i,tval) - itry_lmp(j,tval) 
           rzval = ctrz_lmp(i,tval) - itrz_lmp(j,tval) 
           
           rxval = rxval - box_xl*ANINT(rxval/box_xl)
           ryval = ryval - box_yl*ANINT(ryval/box_yl)
           rzval = rzval - box_zl*ANINT(rzval/box_zl)
           
           rval = sqrt(rxval**2 + ryval**2 + rzval**2)
           
           IF(rval .LT. rcatan_cut) THEN
              
              autocf(i,tval) = 1
              j = ioncnt+1
              
           ELSE

              j = j +1

           END IF

        END DO

     END DO

  END DO
!$OMP END DO

!$OMP DO PRIVATE(tinc,ifin,tim,i,j)

  DO tinc = 0, nframes-1 !compute spectral product

     ifin = nframes - tinc
     
     DO i = 1,ifin
        
        tim = i + tinc
      
        DO j = 1,c_ioncnt

           tplot_cf(tinc) = tplot_cf(tinc) + REAL(autocf(j,tim)&
                &*autocf(j,i))
           
        END DO

     END DO

     tplot_cf(tinc) = tplot_cf(tinc)/REAL(ifin*c_ioncnt)
     
  END DO
!$OMP END DO
!$OMP END PARALLEL

  DO tinc = 0, nframes-1

     WRITE(dumwrite,"(I0,1X,F14.6)") tinc, tplot_cf(tinc)

  END DO

  CLOSE(dumwrite)

END SUBROUTINE RESIDENCE_TIME_CATAN

!--------------------------------------------------------------------

SUBROUTINE OPEN_STRUCT_OUTPUT_FILES()

  USE STATICPARAMS

  IMPLICIT NONE

  IF(rgcalc) THEN
     
     IF(rgavg) THEN
        dum_fname = "rgavg_"//trim(adjustl(traj_fname))
        OPEN(unit = rgavgwrite,file =trim(dum_fname),action="write"&
             &,status="replace")
        WRITE(rgavgwrite,'(A5,1X,4(A3,1X))') "tstep", "rgt","rgx","&
             &rgy","rgz"

     END IF

     IF(rgall) THEN
        dum_fname = "rgall_"//trim(adjustl(traj_fname))
        OPEN(unit = rgwrite,file =trim(dum_fname),action="write"&
             &,status="replace")
     END IF
 
  END IF

   
END SUBROUTINE OPEN_STRUCT_OUTPUT_FILES

!--------------------------------------------------------------------

SUBROUTINE ALLOUTPUTS()

  USE STATICPARAMS
  IMPLICIT NONE
  INTEGER :: i,ierr

  PRINT *, "Number of frames from start to end: ", nframes/(freqfr+1)
  PRINT *, "Frequency of Frames: ", freqfr + 1
  PRINT *, "Total number of Frames analyzed: ", nfrcntr

  WRITE(logout,*) "Number of frames from start to end: ", nframes&
       &/(freqfr+1)
  WRITE(logout,*) "Frequency of Frames: ", freqfr+1
  WRITE(logout,*) "Total number of Frames analyzed: ", nfrcntr

  IF(densflag) THEN
     PRINT *, "Writing densities"
     CALL OUTPUT_ALLDENS()
  END IF
  
  IF(rdfflag .or. bfrdf_calc) THEN
     PRINT *, "Writing RDFs .."
     CALL OUTPUT_ALLRDF()
  END IF

  IF(catan_neighcalc) THEN
     PRINT *, "Writing neighbors .."
     CALL OUTPUT_ALLNEIGHBORS()
  END IF

  IF(clust_calc) THEN

     dum_fname = "clust_"//trim(adjustl(traj_fname))
     OPEN(unit = dumwrite,file =trim(dum_fname),action="write"&
          &,status="replace",iostat=ierr)

     IF(ierr /= 0) PRINT *, "Unknown clust_filename"
     DO i = 1,ntotion_centers

        WRITE(dumwrite,'(I0,1X,F14.8,1X)') i, REAL(clust_avg(i))&
             &/REAL(nframes)

     END DO
     CLOSE(dumwrite)

  END IF

END SUBROUTINE ALLOUTPUTS
  
!--------------------------------------------------------------------

SUBROUTINE OUTPUT_ALLDENS()

  USE STATICPARAMS
  IMPLICIT NONE

  INTEGER :: i,j,ierr
  REAL, PARAMETER :: vconst = 4.0*pival/3.0
  REAL :: frnorm
  
  frnorm   = 1+INT(nfrcntr/densfreq)
  ddenavg  = ddenavg/REAL(frnorm)
  dbinavg  = dbinavg/REAL(frnorm)
  PRINT *, "Average inverse volume of box for densities", ddenavg
     
  dum_fname = "dens_"//trim(adjustl(traj_fname))
  OPEN(unit = dumwrite,file =trim(dum_fname),action="write"&
       &,status="replace",iostat=ierr)
     
  IF(ierr /= 0) THEN
     PRINT *, "Could not open", trim(dum_fname)
  END IF
  
  WRITE(dumwrite,'(4X,A,4X)',advance="no") "r"
  
  DO j = 1,ndens_grps
     
     WRITE(dumwrite,'(I0,4X)',advance="no") INT(densarray(0,j))
     
  END DO
  
  WRITE(dumwrite,*)

  DO i = 0,d_maxbin-1
           
     WRITE(dumwrite,'(F16.5,2X)',advance="no") 0.5*dbinavg&
          &*(REAL(2*i+1))
     densarray(i+1,0) = 0.5*dbinavg*(REAL(2*i+1))
     
     DO j = 1,ndens_grps
        
        densarray(i+1,j) = REAL(densarray(i+1,j))/REAL(frnorm)
        WRITE(dumwrite,'(F16.9,1X)',advance="no")densarray(i+1,j)

     END DO
        
     WRITE(dumwrite,*)
     
  END DO
  
  CLOSE(dumwrite)

  
END SUBROUTINE OUTPUT_ALLDENS

!--------------------------------------------------------------------

SUBROUTINE OUTPUT_LAYERRDF()

  USE STATICPARAMS
  IMPLICIT NONE
  INTEGER :: p,i,j,paircnt
  REAL, PARAMETER :: vconst = 4.0*pival/3.0
  REAL :: rlower,rupper,nideal,rdf2dfrnorm
  CHARACTER(LEN=3) :: intnum
  
  IF(rdf2dfreq == 1) THEN
     rdf2dfrnorm = 1
  ELSE
     rdf2dfrnorm = INT(nframes/rdf2dfreq)+1
  END IF

  rdf2dvolavg = rdf2dvolavg/REAL(rdf2dfrnorm)
  rdf2dbinavg = rdf2dbinavg/REAL(rdf2dfrnorm)
  PRINT *, "rdf2d-volavg/rdf2d-binavg: ", rdf2dvolavg,rdf2dbinavg
  
  DO p = 1,nmax_layers
     
     WRITE(intnum,'(I0)') p
     dum_fname  = "rdf2d_"//trim(intnum)//"_"&
          &//trim(adjustl(traj_fname))//".txt"
     
     OPEN(unit = dumwrite,file = trim(dum_fname), status="replace",&
          & action = "write")

     WRITE(dumwrite,'(5X,A1,2X)',advance="no") "r"

     DO i = 1,npairs_2drdf
        WRITE(dumwrite,'(I0,A1,I0,2X)',advance="no")&
             & pairs_2drdf_arr(i,1),"-",pairs_2drdf_arr(i,2)
     END DO
     WRITE(dumwrite,*)
     
     DO i = 0,rdf2dmaxbin-1
        
        rlower = real(i)*rdf2dbinavg
        rupper = rlower + rdf2dbinavg
        nideal = vconst*(rupper**3 - rlower**3)
        WRITE(dumwrite,'(F16.9,2X)',advance="no") 0.5*rdf2dbinavg&
             &*(REAL(2*i+1))

        DO paircnt = 1,npairs_2drdf
        
           WRITE(dumwrite,'(F16.9,2X)',advance="no") rdf2darray(i&
                &,paircnt,p)/(rdf2dfrnorm*nideal)

        END DO
        
        WRITE(dumwrite,*)
        
     END DO
     
     CLOSE(dumwrite)
     
  END DO
  
END SUBROUTINE OUTPUT_LAYERRDF

!--------------------------------------------------------------------

SUBROUTINE OUTPUT_ALLRDF()

  USE STATICPARAMS
  IMPLICIT NONE

  INTEGER :: i,j,ierr
  REAL, PARAMETER :: vconst = 4.0*pival/3.0
  REAL :: rlower,rupper,nideal,rdffrnorm,acrnorm
  
  IF(rdfflag .or. bfrdf_calc) THEN

     rdffrnorm = INT(nfrcntr/rdffreq)
     rvolavg = rvolavg/REAL(rdffrnorm)
     PRINT *, "Average volume of box", rvolavg
     
     IF(rdfflag) THEN
        dum_fname = "rdf_"//trim(adjustl(traj_fname))
        OPEN(unit = dumwrite,file =trim(dum_fname),action="write"&
             &,status="replace",iostat=ierr)
     
        IF(ierr /= 0) THEN
           PRINT *, "Could not open", trim(dum_fname)
        END IF
     
        WRITE(dumwrite,'(A,2X)',advance="no") "r"
     
        DO j = 1,npairs
        
           WRITE(dumwrite,'(2(I0,1X))',advance="no") pairs_rdf(j,1)&
                &,pairs_rdf(j,2)
        
        END DO
     
        WRITE(dumwrite,*)
     
        DO i = 0,rmaxbin-1
           
           rlower = real(i)*rbinval
           rupper = rlower + rbinval
           nideal = vconst*(rupper**3 - rlower**3)
        
           WRITE(dumwrite,'(F16.5,2X)',advance="no") 0.5*rbinval&
                &*(REAL(2*i+1))
        
           DO j = 1,npairs
           
              WRITE(dumwrite,'(F16.9,1X)',advance="no")rdfarray(i,j)&
                &/(rdffrnorm*nideal)
              
           END DO
        
           WRITE(dumwrite,*)
        
        END DO
     
        CLOSE(dumwrite)

     END IF

     IF(bfrdf_calc) THEN

        dum_fname = "freeboundrdf_"//trim(adjustl(traj_fname))
        OPEN(unit = dumwrite,file =trim(dum_fname),action="write"&
             &,status="replace",iostat=ierr)
     
        IF(ierr /= 0) THEN
           PRINT *, "Could not open", trim(dum_fname)
        END IF

        rlower = real(i)*rbinval
        rupper = rlower + rbinval
        nideal = vconst*(rupper**3 - rlower**3)
        
        DO i = 0,rmaxbin-1

           WRITE(dumwrite,'(4(F16.9,1X))') 0.5*rbinval*(REAL(2*i&
                &+1)),rdf_p_fb(i)/(rdffrnorm*nideal),rdf_p_ff(i)&
                &/(rdffrnorm*nideal),rdf_p_bb(i)/(rdffrnorm*nideal)

        END DO

        CLOSE(dumwrite)

     END IF

  END IF

END SUBROUTINE OUTPUT_ALLRDF

!--------------------------------------------------------------------

SUBROUTINE OUTPUT_ALLNEIGHBORS()

  USE STATICPARAMS

  IMPLICIT NONE

  INTEGER :: i,frnorm,ierr
  REAL :: totcat_an_neigh,totan_cat_neigh

  IF(neighfreq == 1) frnorm = nframes
  IF(neighfreq .NE. 1) frnorm = nframes/neighfreq + 1

  totcat_an_neigh = 0.0; totan_cat_neigh = 0.0

  IF(catan_neighcalc) THEN
!$OMP PARALLEL DO REDUCTION(+:totcat_an_neigh,totan_cat_neigh) PRIVATE(i)
  
     DO i = 1,maxneighsize

        totcat_an_neigh = totcat_an_neigh + REAL(cat_an_neighavg(i))
        totan_cat_neigh = totan_cat_neigh + REAL(an_cat_neighavg(i))

     END DO

!$OMP END PARALLEL DO

     IF(catan_neighcalc) THEN
        
        dum_fname = "catanneigh_"//trim(adjustl(traj_fname))
        OPEN(unit = dumwrite,file =trim(dum_fname),action="write"&
             &,status="replace")

        
     END IF

     
     DO i = 1,maxneighsize     

        WRITE(dumwrite,'(I0,1X,4(F14.8,1X))') i,&
             & REAL(cat_an_neighavg(i))/REAL(frnorm),100.0&
             &*REAL(cat_an_neighavg(i))/totcat_an_neigh&
             &,REAL(an_cat_neighavg(i))/REAL(frnorm),100.0&
             &*REAL(an_cat_neighavg(i))/totan_cat_neigh

     END DO

     CLOSE(dumwrite)

  END IF

END SUBROUTINE OUTPUT_ALLNEIGHBORS
     
!--------------------------------------------------------------------

SUBROUTINE ALLOCATE_TOPO_ARRAYS()

  USE STATICPARAMS
  IMPLICIT NONE

  INTEGER :: AllocateStatus

! Allocate LAMMPS structure

  ALLOCATE(aidvals(ntotatoms,3),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate aidvals"
  ALLOCATE(rxyz_lmp(ntotatoms,3),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate rxyz_lmp"
  ALLOCATE(charge_lmp(ntotatoms,1),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate charge_lmp"
  ALLOCATE(vel_xyz(ntotatoms,4),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate vel_xyz"

  IF(ntotbonds /= 0) THEN
     ALLOCATE(bond_lmp(ntotbonds,4),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate bond_lmp"
  ELSE
     PRINT *, "Warning: No bonds - Not correct for bonded systems"
     ALLOCATE(bond_lmp(1,1),stat = AllocateStatus)
     DEALLOCATE(bond_lmp)
  END IF
  
  IF(ntotangls /= 0) THEN
     ALLOCATE(angl_lmp(ntotangls,5),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate angl_lmp"
  ELSE
     ALLOCATE(angl_lmp(1,1),stat = AllocateStatus)
     DEALLOCATE(angl_lmp)
  END IF
     
  IF(ntotdihds /= 0) THEN
     ALLOCATE(dihd_lmp(ntotdihds,6),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate dihd_lmp"
  ELSE
     ALLOCATE(dihd_lmp(1,1),stat = AllocateStatus)
     DEALLOCATE(dihd_lmp)
  END IF
  
  IF(ntotimprs /= 0) THEN
     ALLOCATE(impr_lmp(ntotimprs,6),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate zlmp"
  ELSE
     ALLOCATE(impr_lmp(1,1),stat = AllocateStatus)
     DEALLOCATE(impr_lmp)
  END IF

! Allocate box details

  ALLOCATE(boxx_arr(nframes),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate boxx_arr"
  ALLOCATE(boxy_arr(nframes),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate boxy_arr"
  ALLOCATE(boxz_arr(nframes),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate boxz_arr"

  PRINT *, "Successfully allocated memory for topology"

END SUBROUTINE ALLOCATE_TOPO_ARRAYS

!--------------------------------------------------------------------

SUBROUTINE ALLOCATE_ANALYSIS_ARRAYS()

  USE STATICPARAMS
  IMPLICIT NONE

  INTEGER :: AllocateStatus

! Allocate for statics

  IF(rdfflag) THEN
     ALLOCATE(rdfarray(0:rmaxbin-1,npairs),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate rdfarray"
  ELSE
     ALLOCATE(rdfarray(1,1),stat = AllocateStatus)
     DEALLOCATE(rdfarray)
  END IF

  IF(bfrdf_calc) THEN
     ALLOCATE(rdf_p_fb(0:rmaxbin-1),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate rdf_p_fb"
     ALLOCATE(rdf_p_bb(0:rmaxbin-1),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate rdf_p_bb"
     ALLOCATE(rdf_p_ff(0:rmaxbin-1),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate rdf_p_ff"
  ELSE
     ALLOCATE(rdf_p_fb(1),stat = AllocateStatus)
     DEALLOCATE(rdf_p_fb)
     ALLOCATE(rdf_p_bb(1),stat = AllocateStatus)
     DEALLOCATE(rdf_p_bb)
     ALLOCATE(rdf_p_ff(1),stat = AllocateStatus)
     DEALLOCATE(rdf_p_ff)
  END IF

  IF(densflag == 0) THEN
     ALLOCATE(densarray(1,1),stat = AllocateStatus)
     DEALLOCATE(densarray)
  END IF

! Allocate for interfacial calcualtions

  IF(layerana_flag) THEN
     ALLOCATE(trx_lmp(ntotatoms,nframes),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate itrx_lmp"
     ALLOCATE(try_lmp(ntotatoms,nframes),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate itry_lmp"
     ALLOCATE(trz_lmp(ntotatoms,nframes),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate itrz_lmp"
  ELSE
     ALLOCATE(trx_lmp(1,nframes),stat = AllocateStatus)
     ALLOCATE(try_lmp(1,nframes),stat = AllocateStatus)
     ALLOCATE(trz_lmp(1,nframes),stat = AllocateStatus)
     DEALLOCATE(trx_lmp)
     DEALLOCATE(try_lmp)
     DEALLOCATE(trz_lmp)
  END IF
     
  IF(rdf2dflag) THEN
     ALLOCATE(rdf2darray(0:rdf2dmaxbin-1,npairs_2drdf,nmax_layers)&
          &,stat=AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate rdf2darray"
  ELSE
     ALLOCATE(rdf2darray(1,1,1))
     DEALLOCATE(rdf2darray)
  END IF
  
! Allocate for dynamics 

  IF(ion_dynflag .OR. cion_dynflag .OR. pion_dynflag) THEN
     ALLOCATE(tarr_lmp(nframes),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate tarr_lmp"
  ELSE
     ALLOCATE(tarr_lmp(1),stat = AllocateStatus)
     DEALLOCATE(tarr_lmp)
  END IF
  
  IF(ion_dynflag) THEN
     ALLOCATE(itrx_lmp(ioncnt,nframes),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate itrx_lmp"
     ALLOCATE(itry_lmp(ioncnt,nframes),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate itry_lmp"
     ALLOCATE(itrz_lmp(ioncnt,nframes),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate itrz_lmp"
  ELSE
     ALLOCATE(itrx_lmp(1,nframes),stat = AllocateStatus)
     ALLOCATE(itry_lmp(1,nframes),stat = AllocateStatus)
     ALLOCATE(itrz_lmp(1,nframes),stat = AllocateStatus)
     DEALLOCATE(itrx_lmp)
     DEALLOCATE(itry_lmp)
     DEALLOCATE(itrz_lmp)
  END IF

  IF(cion_dynflag) THEN
     ALLOCATE(ctrx_lmp(c_ioncnt,nframes),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate ctrx_lmp"
     ALLOCATE(ctry_lmp(c_ioncnt,nframes),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate ctry_lmp"
     ALLOCATE(ctrz_lmp(c_ioncnt,nframes),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate ctrz_lmp"
  ELSE
     ALLOCATE(ctrx_lmp(1,nframes),stat = AllocateStatus)
     ALLOCATE(ctry_lmp(1,nframes),stat = AllocateStatus)
     ALLOCATE(ctrz_lmp(1,nframes),stat = AllocateStatus)
     DEALLOCATE(ctrx_lmp)
     DEALLOCATE(ctry_lmp)
     DEALLOCATE(ctrz_lmp)
  END IF

  IF(pion_dynflag) THEN
     ALLOCATE(ptrx_lmp(p_ioncnt,nframes),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate ptrx_lmp"
     ALLOCATE(ptry_lmp(p_ioncnt,nframes),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate ptry_lmp"
     ALLOCATE(ptrz_lmp(p_ioncnt,nframes),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate ptrz_lmp"
  ELSE
     ALLOCATE(ptrx_lmp(1,nframes),stat = AllocateStatus)
     ALLOCATE(ptry_lmp(1,nframes),stat = AllocateStatus)
     ALLOCATE(ptrz_lmp(1,nframes),stat = AllocateStatus)
     DEALLOCATE(ptrx_lmp)
     DEALLOCATE(ptry_lmp)
     DEALLOCATE(ptrz_lmp)
  END IF
  
  IF(catan_neighcalc) THEN
     ALLOCATE(cat_an_neighavg(1:maxneighsize),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate cat_an_neighavg"
     ALLOCATE(an_cat_neighavg(1:maxneighsize),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate an_cat_neighavg"
  ELSE
     ALLOCATE(cat_an_neighavg(1),stat = AllocateStatus)
     ALLOCATE(an_cat_neighavg(1),stat = AllocateStatus)
     DEALLOCATE(cat_an_neighavg)
     DEALLOCATE(an_cat_neighavg)
  END IF

  PRINT *, "Successfully allocated memory for analyis"

END SUBROUTINE ALLOCATE_ANALYSIS_ARRAYS

!--------------------------------------------------------------------

SUBROUTINE DEALLOCATE_ARRAYS()

  USE STATICPARAMS

  IMPLICIT NONE

  !Global arrays
  DEALLOCATE(aidvals)
  DEALLOCATE(rxyz_lmp)
  DEALLOCATE(vel_xyz)
  DEALLOCATE(boxx_arr)
  DEALLOCATE(boxy_arr)
  DEALLOCATE(boxz_arr)

  !Topo arrays
  IF(ntotbonds /= 0) DEALLOCATE(bond_lmp)
  IF(ntotangls /= 0) DEALLOCATE(angl_lmp)
  IF(ntotdihds /= 0) DEALLOCATE(dihd_lmp)
  IF(ntotimprs /= 0) DEALLOCATE(impr_lmp)
  
  !Statics calculations arrays
  IF(rdfflag) DEALLOCATE(rdfarray)
  IF(densflag) DEALLOCATE(densarray)
  
  !Dynamic calculations arrays
  IF(ion_dynflag) THEN
     DEALLOCATE(itrx_lmp)
     DEALLOCATE(itry_lmp)
     DEALLOCATE(itrz_lmp)
  END IF

  IF(cion_dynflag) THEN
     DEALLOCATE(ctrx_lmp)
     DEALLOCATE(ctry_lmp)
     DEALLOCATE(ctrz_lmp)
  END IF

  IF(pion_dynflag) THEN
     DEALLOCATE(ptrx_lmp)
     DEALLOCATE(ptry_lmp)
     DEALLOCATE(ptrz_lmp)
  END IF

  !Dellocate groups
  IF(grpflag) THEN
     DEALLOCATE(types_in_group)
     DEALLOCATE(allgrp_typarr)
  END IF

  !Deallocate interface arrays
  IF(interfflag) THEN
     DEALLOCATE(interpos)
     DEALLOCATE(widdoms)
     DEALLOCATE(all_layergrp_typarr)
     DEALLOCATE(domtyp)
  END IF

  !Deallocate lmp_x,y,z as a function of time
  IF(layerana_flag) THEN
     DEALLOCATE(trx_lmp)
     DEALLOCATE(try_lmp)
     DEALLOCATE(trz_lmp)
  END IF
    
  !Layerwise calculation
  IF(rdf2dflag) THEN
     DEALLOCATE(rdf2darray)
     DEALLOCATE(pairs_2drdf_arr)
  END IF
  
  CLOSE(logout)
  
END SUBROUTINE DEALLOCATE_ARRAYS

!--------------------------------------------------------------------
