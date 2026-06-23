!---------------To analyze properties of bulk-sei systems------------
!---------------Version 2: Dev_Apr-27-2026---------------------------
!---------------Parameter File: params_statics.f90-------------------
!********************************************************************

PROGRAM ANALYSIS_MAIN

  USE SUBROUTINE_DEFS
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
  
  CALL DEALLOCATE_ARRAYS()

! Print completion
  PRINT *, "All Calculations Completed Succesfully :)"

END PROGRAM ANALYSIS_MAIN

!--------------------------------------------------------------------

SUBROUTINE READ_ANA_IP_FILE()

  USE SUBROUTINE_DEFS
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

  OPEN(unit=anaread,file=trim(ana_fname),action="read",status="old"&
       &,iostat=ierr)
  
  IF(ierr /= 0) THEN

     PRINT *, trim(ana_fname), "not found"
     STOP

  END IF

  DO

     dumchar = ''
     CALL READ_NEXT_KEYWORD(anaread, dumchar, ierr)

     IF(ierr < 0) EXIT
     IF(ierr > 0) STOP "Error reading input file"

     PRINT *, "Processing keyword: ", trim(dumchar)
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

        ALLOCATE(allgrp_atomarr(0:ntotatoms,ngroups),stat&
             &=AllocateStatus)
        IF(AllocateStatus/=0) STOP "did not allocate allgrp_atomarr"

        allgrp_atomarr = 0 ! ZERO everything

        ALLOCATE(allgrp_typarr(1:ntotatomtypes,ngroups),stat&
             &=AllocateStatus)
        IF(AllocateStatus/=0) STOP "did not allocate allgrp_typarr"

        allgrp_atomarr = 0 ! ZERO everything

        DO i = 1,ngroups
           
           READ(anaread,'(A)',iostat=ierr) charline
           IF(ierr /= 0) THEN
              PRINT *, "Error reading group size at group-",i
              STOP
           END IF
           
           READ(charline,*) allgrp_atomarr(0,i), ntypes_per_group
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

     ELSEIF(dumchar == 'layer_groups_interface') THEN
        IF(layer_grpflag_surf == 1) STOP "Conflicting options: Canno&
             &t bin based on both interface and bottom surface"
        layer_grpflag_interf = 1
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

     ELSEIF(dumchar == 'layer_groups_surface') THEN

        IF(layer_grpflag_interf == 1) STOP "Conflicting options: Canno&
             &t bin based on both interface and bottom surface"
        
        layer_grpflag_surf = 1
        READ(anaread,*,iostat=ierr) nlayer_groups, nmax_layers,&
             & zmin, zmax, deltaz_btw_layers

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

        IF(layer_grpflag_interf == 0 .AND. layer_grpflag_surf == 0)&
             & THEN
           PRINT *, "ERROR: Define domain subdivision (interface or su&
                &rface)"
           STOP
        END IF
           
        
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

     ELSEIF(dumchar == 'compute_ciondiff') THEN

        READ(anaread,*,iostat=ierr) cion_diff, delta_t
        cion_dynflag = 1

     ELSEIF(dumchar == 'compute_catanrestime') THEN
        
        READ(anaread,*,iostat=ierr) rcatan_cut
        catan_autocfflag = 1
        ion_dynflag = 1; cion_dynflag = 1

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

!!$  CALL SANITY_CHECK_IONTYPES()

END SUBROUTINE READ_ANA_IP_FILE

!--------------------------------------------------------------------

SUBROUTINE READ_NEXT_KEYWORD(iu, keyword, ierr)

  USE SUBROUTINE_DEFS
  USE STATICPARAMS
  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: iu
  CHARACTER(LEN=*), INTENT(OUT) :: keyword
  INTEGER, INTENT(out) :: ierr
  INTEGER :: ios
  character(len=256) :: dumline

  keyword = ''
  ierr = 0
  
  DO

     dumline = ''
     READ(iu,'(A256)',iostat=ios) dumline

      IF(ios < 0) THEN
         ierr = -1       ! EOF
         RETURN
      ELSEIF(ios > 0) then
         ierr = ios      ! Read error
         RETURN
      END IF
     
      ierr = 0
      dumline = ADJUSTL(dumline)
      
      ! Skip blank lines
      IF(len_trim(dumline) == 0) CYCLE
      
      ! Skip comments
      IF(dumline(1:1) == "#" .OR. dumline(1:1) == "!") CYCLE
      
      ! Extract first word from non-blank line
      READ(dumline, *, iostat=ios) keyword

      IF(ios /= 0) THEN
         ierr = ios
         RETURN
      END IF

      ierr = 0
      RETURN
     
  END DO
  
END SUBROUTINE READ_NEXT_KEYWORD

!--------------------------------------------------------------------

SUBROUTINE DEFAULTVALUES()

  USE SUBROUTINE_DEFS
  USE STATICPARAMS
  IMPLICIT NONE

  ! Frame, molecules and processor details
  nframes = 0; skipfr = 0; freqfr = 0; nfrcntr = 0
  nchains = 0; atperchain = 0

  ! Initialize flags
  grpflag = 0; densflag = 0; layerana_flag = 0
  interfflag = 0;
  rgall = 0; rgcalc = 0; rdfflag = 0
  ion_dynflag = 0; cion_dynflag = 0
  ion_diff = 0; cion_diff = 0
  bfrdf_calc = 0
  catan_autocfflag = 0; catpol_autocfflag = 0
  rdf2dflag = 0
  layer_grpflag_interf = 0; layer_grpflag_surf = 0
  
  ! Initialize iontypes
  c_iontype = -1; iontype = -1

  !Initialize system quantities
  ioncnt = 0; c_ioncnt = 0

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
  major_boxval = 0.0; boxlzavg = 0.0

  ! Initialize layer-wise structural averages
  rdf2dvolavg = 0.0; rdf2dfreq = 1; rdf2dmaxbin = 50
  rdf2dbinavg = 0.0
  
  ! Initialize dynamical quantities
  rcatpol_cut1 = 0.0; rcatpol_cut2 = 0.0
  
END SUBROUTINE DEFAULTVALUES

!--------------------------------------------------------------------

SUBROUTINE READ_DATAFILE()

  USE SUBROUTINE_DEFS
  USE STATICPARAMS
  IMPLICIT NONE

  INTEGER :: i,j,ierr,u,AllocateStatus,imax
  INTEGER :: flag, cntr
  INTEGER :: aid,molid,atype
  REAL    :: charge,rx,ry,rz
  REAL    :: xlo,xhi,ylo,yhi,zlo,zhi
  CHARACTER(256) :: dumchar

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

  USE SUBROUTINE_DEFS
  USE STATICPARAMS
  IMPLICIT NONE

  INTEGER, INTENT(OUT) :: imax
  INTEGER :: pos, nwords,lcnt,ierr, flag
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

SUBROUTINE FILL_GROUP_ARRAY(group_id,ntypes_per_group&
     &,types_in_group_id)

  USE STATICPARAMS

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: group_id, ntypes_per_group
  INTEGER, INTENT(IN) :: types_in_group_id(1:ntypes_per_group)
  INTEGER :: i,cnt_atoms

  cnt_atoms = 0

  DO i = 1,ntotatoms

     IF(ANY(aidvals(i,3) == types_in_group_id)) THEN

        allgrp_atomarr(i,group_id) = 1
        allgrp_typarr(aidvals(i,3),group_id) = 1
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
  INTEGER :: i,cnt_atoms

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
  WRITE(logout,*) "Number of atoms in layer group-", group_id, " is "&
       &, cnt_atoms
  WRITE(logout,*) "************************************************"
  PRINT *, "Number of atoms in layer group-", group_id, " is ",&
       & cnt_atoms

END SUBROUTINE FILL_LAYERGROUP_ARRAY

!--------------------------------------------------------------------

SUBROUTINE OUTPUT_ALL_GROUPS()

  USE STATICPARAMS
  IMPLICIT NONE

  INTEGER :: i,j,ierr,cntr

  dum_fname = "allgroups_id_types.dat"
  OPEN(unit = dumwrite,file =trim(dum_fname),action="write"&
       &,status="replace",iostat=ierr)
  IF(ierr /= 0) THEN
     PRINT *, "Could not open", trim(dum_fname)
  END IF

  DO i = 1,ngroups

     WRITE(dumwrite,*) "GROUP-ID: ", allgrp_atomarr(0,i)

     cntr = 0
     WRITE(dumwrite,*) "TYPES IN GROUP"
     DO j = 1,ntotatomtypes

        IF(allgrp_typarr(j,i) == 1) THEN
           
           WRITE(dumwrite,'(I0,2X)',advance="no") aidvals(j,1)
           cntr = cntr + 1
           IF(MOD(cntr,8) == 0 .AND. cntr /= 0) WRITE(dumwrite,*)

        END IF

     END DO

     WRITE(dumwrite,*)

     cntr = 0
     WRITE(dumwrite,*) "ATOM-IDS IN GROUP: "
     DO j = 1,ntotatoms

        IF(allgrp_atomarr(j,i) == 1) THEN

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

  USE SUBROUTINE_DEFS
  USE STATICPARAMS
  IMPLICIT NONE

  INTEGER :: i,j,aid,ierr,atchk,atype,jumpfr,AllocateStatus
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

  ! Allocate box details

  ALLOCATE(box_arr(nframes,3),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate box_arr"

  box_arr = 0.0
  
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

     box_arr(i,1)  = box_xl
     box_arr(i,2)  = box_yl
     box_arr(i,3)  = box_zl

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

  USE SUBROUTINE_DEFS
  USE STATICPARAMS
  IMPLICIT NONE

  INTEGER :: i,j,t1,t2

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

  USE SUBROUTINE_DEFS
  USE STATICPARAMS
  IMPLICIT NONE

  INTEGER, INTENT(IN):: tval
  INTEGER :: t1, t2
  INTEGER :: clock_rate, clock_max
  
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
        
!!$  IF(rdfflag) THEN
!!$     
!!$     IF(tval == 1) THEN
!!$
!!$        CALL SYSTEM_CLOCK(t1,clock_rate,clock_max)
!!$        CALL COMPUTE_RDF(tval)
!!$        CALL SYSTEM_CLOCK(t2,clock_rate,clock_max)
!!$        PRINT *, 'Elapsed real time for RDF analysis: ',REAL(t2&
!!$             &-t1)/REAL(clock_rate), ' seconds'
!!$
!!$     END IF
!!$     IF(mod(tval-1,rdffreq)==0) CALL COMPUTE_RDF(tval)     
!!$     
!!$  END IF
!!$
!!$  IF(catan_neighcalc) THEN
!!$     
!!$     IF(tval == 1) THEN
!!$
!!$        cat_an_neighavg = 0.0; an_cat_neighavg=0.0
!!$        CALL SYSTEM_CLOCK(t1,clock_rate,clock_max)
!!$        CALL CAT_AN_NEIGHS()
!!$        CALL SYSTEM_CLOCK(t2,clock_rate,clock_max)
!!$        PRINT *, 'Elapsed real time for neighbor analysis: ',REAL(t2&
!!$             &-t1)/REAL(clock_rate), ' seconds'
!!$
!!$     END IF
!!$     
!!$     IF(mod(tval,neighfreq) == 0) CALL CAT_AN_NEIGHS()
!!$     
!!$  END IF
!!$
!!$  IF(bfrdf_calc) THEN
!!$     
!!$     IF(tval == 1) THEN
!!$        rdf_p_fb = 0.0; rdf_p_ff=0.0; rdf_p_bb = 0.0        
!!$        CALL SYSTEM_CLOCK(t1,clock_rate,clock_max)
!!$        CALL  SORT_POLY_FREE_BOUND_COMPLEX(tval)
!!$        CALL SYSTEM_CLOCK(t2,clock_rate,clock_max)
!!$        PRINT *, 'Elapsed real time for bound/free RDF: ',REAL(t2&
!!$             &-t1)/REAL(clock_rate), ' seconds'           
!!$     END IF
!!$     
!!$     IF(mod(tval,rdffreq) == 0) CALL&
!!$          & SORT_POLY_FREE_BOUND_COMPLEX(tval)
!!$     
!!$  END IF
!!$
!!$  IF(clust_calc) THEN
!!$
!!$     IF(tval == 1) THEN
!!$
!!$        clust_avg = 0
!!$        CALL SYSTEM_CLOCK(t1,clock_rate,clock_max)
!!$        CALL CLUSTER_ANALYSIS(tval)
!!$        CALL SYSTEM_CLOCK(t2,clock_rate,clock_max)
!!$        PRINT *, 'Elapsed real time for cluster analysis= ',REAL(t2&
!!$             &-t1)/REAL(clock_rate), ' seconds'           
!!$     ELSE
!!$        
!!$        CALL CLUSTER_ANALYSIS(tval)
!!$     
!!$     END IF
!!$
!!$  END IF


END SUBROUTINE STRUCT_MAIN

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_DENSPROFILES(tval)

  USE SUBROUTINE_DEFS
  USE STATICPARAMS

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: tval
  INTEGER :: i,j,ibin,grp_cntr,grp_id,grp_colid,atype
  REAL,DIMENSION(0:d_maxbin-1,1:ndens_grps) :: grpdens
  REAL :: dbinval,rval,ddenval
  INTEGER :: flag

  dbinval  = major_boxval/d_maxbin
  dbinavg  = dbinavg + dbinval
  ddenval  = major_boxval/(box_xl*box_yl*box_zl*dbinval)
  ddenavg  = ddenavg + ddenval
  boxlzavg = boxlzavg + major_boxval

  IF(tval == 1) THEN
     
     DO i = 1,ndens_grps
        
        flag = -1
        
        DO j = 1,ngroups
           
           IF(INT(densarray(0,i)) == allgrp_atomarr(0,j)) THEN
              
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

        IF(grp_id == allgrp_atomarr(0,j)) THEN

           grp_colid = j
           EXIT

        END IF

     END DO
           
     DO i = 1,ntotatoms
     
        IF(allgrp_atomarr(i,grp_colid)) THEN

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

  INTEGER :: i,a1type,cnt,AllocateStatus,ntotion_cnt
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

  USE SUBROUTINE_DEFS
  USE STATICPARAMS
  IMPLICIT NONE

  INTEGER :: t1, t2, i
  INTEGER :: clock_rate, clock_max
  REAL :: rvolval
  
  PRINT *, "Computing interfaces..."
  CALL SYSTEM_CLOCK(t1,clock_rate,clock_max)
  CALL COMPUTE_INTERFACES()
  CALL SYSTEM_CLOCK(t2,clock_rate,clock_max)
  PRINT *, 'Elapsed real time for interface analysis: ',REAL(t2-t1)&
       &/REAL(clock_rate), ' seconds'

  box_xl = box_arr(1,1)
  box_yl = box_arr(1,2)
  box_zl = box_arr(1,3)  
  rvolval = box_xl*box_yl*box_zl

  
  ! Classify domain according to **number** densities
  CALL INITIAL_DOMAIN_CLASSIFICATION(1) 

  IF(layer_grpflag_interf == 1) THEN
     
     DO i = 1,nframes

        box_xl = box_arr(i,1)
        box_yl = box_arr(i,2)
        box_zl = box_arr(i,3)

        CALL COMPARTMENTALIZE_PARTICLES_INTERFACE(i)
        
     END DO
     
  ELSEIF(layer_grpflag_surf==1) THEN

     DO i = 1,nframes

        box_xl = box_arr(i,1)
        box_yl = box_arr(i,2)
        box_zl = box_arr(i,3)

        CALL COMPARTMENTALIZE_PARTICLES_BOTTOMSURF(i)

     END DO

  END IF

  IF(rdf2dflag) CALL OUTPUT_LAYERRDF()
  
END SUBROUTINE LAYERWISE_MAIN

!--------------------------------------------------------------------
  
SUBROUTINE COMPUTE_INTERFACES()

  USE SUBROUTINE_DEFS
  USE STATICPARAMS
  IMPLICIT NONE

!!$ To calculate the interface positions
  REAL*8  :: duma, dumb
  INTEGER :: i,j,AllocateStatus,maxdens,icount,iflag,ierr,jinit
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
  PRINT *, "Interfacial position(s) is/are: "
  PRINT *, dumdens_interarr(1:maxdens)

  WRITE(logout,*) "Number of interfaces from density: ", maxdens
  WRITE(logout,*) "Interfacial position(s) is/are: "
  WRITE(logout,*) dumdens_interarr(1:maxdens)
  
! Need to determine which one needs to be used.
! Allocate to actual array. This can keep things clean

  maxinterf = maxdens

  IF(maxinterf .NE. 1 ) PRINT *, "More than one interface found.."
  
  ALLOCATE (interpos(maxinterf), stat = AllocateStatus)
  IF(AllocateStatus /=0 ) STOP "***Allocation inter not proper***"

  ! For non-periodic, domtyp = n_interfaces + 1
  ALLOCATE (domtyp(maxinterf+1), stat = AllocateStatus)
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
  ! This is for non-periodic boundary conditions
  ! n_domains = n_interfaces + 1
  ! Only done at t = 0 using the density profile
  USE STATICPARAMS
  IMPLICIT NONE

  INTEGER :: i, flagdom, AllocateStatus
  REAL :: par_pos, boxval
  INTEGER, DIMENSION(1:maxinterf+1) :: domAdum, domBdum
  INTEGER, INTENT(IN) :: tval

  IF(major_axis == 1) boxval = box_arr(tval,1)
  IF(major_axis == 2) boxval = box_arr(tval,2)
  IF(major_axis == 3) boxval = box_arr(tval,3)

  DO i = 1, maxinterf+1
        
     domAdum(i) = 0
     domBdum(i) = 0
     domtyp(i)  = 0
     
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

        IF(allgrp_typarr(aidvals(i,3),interfgrp_a)) THEN

           domAdum(1) = domAdum(1) + 1

        ELSEIF(allgrp_typarr(aidvals(i,3),interfgrp_b)) THEN

           domBdum(1) = domBdum(1) + 1

        END IF

     END IF
     
     IF(par_pos > interpos(1) .AND. par_pos .LE. boxval) THEN

        IF(flagdom == 1) STOP "particle binned twice"

        IF(allgrp_typarr(aidvals(i,3),interfgrp_a)) THEN
           domAdum(2) = domAdum(2) + 1
           
        ELSEIF(allgrp_typarr(aidvals(i,3),interfgrp_b)) THEN
           
           domBdum(2) = domBdum(2) + 1
           
        END IF
        
     END IF

  END DO

  !Nomenclature: domA=1,domB=2
  DO i = 1, maxinterf+1 ! For non-periodic
     
     IF(domAdum(i) > domBdum(i)) THEN
        
        domtyp(i) = 1
        
     ELSEIF(domBdum(i) > domAdum(i)) THEN
        
        domtyp(i) = 2
        
     ELSE
        
        PRINT *,"Something wrong in distribution"
        PRINT *,"Number of particles in",i,"dom", domAdum(i),&
             & domBdum(i)
        PRINT *, "Interfacial groups: ", interfgrp_a, interfgrp_b
        
        WRITE(logout,*),"Something wrong in distribution"
        WRITE(logout,*),"Number of particles in",i,"dom", domAdum(i),&
             & domBdum(i)

        
        STOP
        
     END IF
     
  END DO
  
  WRITE(logout,*) "The domain types are "
  
  DO i = 1,maxinterf+1 !for non_periodic
     
     IF(domtyp(i) == 1) THEN

        PRINT *, i,domtyp(i),interfgrp_a, "Domain A"
        WRITE(logout,*) i,domtyp(i),interfgrp_a,"Domain A"
        
     ELSEIF(domtyp(i) == 2) THEN

        PRINT *, i,domtyp(i),interfgrp_b, "Domain B"
        WRITE(logout,*) i,domtyp(i),interfgrp_b,"Domain B"

     ELSE

        WRITE(logout,*) "Unknown domain type ",i,domtyp(i)
        STOP
        
     END IF
     
  END DO
     
  ALLOCATE(widdoms(maxinterf+1),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate widdoms"
  
  widdoms(1) = interpos(1)
  widdoms(2) = boxval - interpos(maxinterf)
  
END SUBROUTINE INITIAL_DOMAIN_CLASSIFICATION

!--------------------------------------------------------------------

SUBROUTINE COMPARTMENTALIZE_PARTICLES_INTERFACE(tval)

  USE SUBROUTINE_DEFS
  USE STATICPARAMS
  IMPLICIT NONE

  INTEGER :: i,j,k,p,AllocateStatus,flagbin,a1id
  REAL    :: par_pos, domcheck
  REAL    :: segeps_domA,segeps_domB, eps1,eps2
  REAL    :: interin_domA, interout_domA,domwidth_domA
  REAL    :: interin_domB, interout_domB,domwidth_domB
  REAL    :: segwid_domA, segwid_domB, segwid
  REAL    :: boxval
  INTEGER :: interfdomA_col, interfdomB_col
  CHARACTER(LEN=3) :: epsnum
  INTEGER, INTENT(IN) :: tval
  INTEGER, ALLOCATABLE, DIMENSION(:)::dum_aid, dum_typ 

  ! Find boxval
  IF(major_axis == 1) boxval = box_arr(tval,1)
  IF(major_axis == 2) boxval = box_arr(tval,2)
  IF(major_axis == 3) boxval = box_arr(tval,3)

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
        domwidth_domB = densarray(i,0) - interpos(1)
        EXIT
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
        WRITE(logout,*) "Width considered for segment in DomA: "&
             &,(eps2-eps1)*segeps_domA
        WRITE(logout,*) "Width considered for segment in DomB: "&
             &,(eps2-eps1)*segeps_domB
        WRITE(logout,*) "Population of each type"
        DO i = 1,ntotatomtypes
           WRITE(logout,*) i,seg_typcnt(i)
        END DO
     END IF
          
     num_mons_per_layer = k

     WRITE(logout,*) "Number of particles in layer# ", p , "at tval = &
          & ", tval, "is ", num_mons_per_layer

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

     segwid_domA = (eps2-eps1)*segeps_domA
     segwid_domA = (eps2-eps1)*segeps_domB
     segwid = MAX(segwid_domA,segwid_domB)
     
     CALL LAYERWISE_ANALYSIS(tval,p,num_mons_per_layer,segwid)

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
  DEALLOCATE(seg_typcnt)
  
!!$  WRITE(logout,*) "Segmental diffusion analysis complete..."
!!$  WRITE(logout,*) "Segmental mode relaxation analysis complete..."
  
END SUBROUTINE COMPARTMENTALIZE_PARTICLES_INTERFACE

!--------------------------------------------------------------------

SUBROUTINE COMPARTMENTALIZE_PARTICLES_BOTTOMSURF(tval)

  USE SUBROUTINE_DEFS
  USE STATICPARAMS

  IMPLICIT NONE

  INTEGER :: i,k,p,AllocateStatus,flagbin,a1id
  REAL    :: par_pos, layer_deltaz, delz_1, delz_2
  REAL    :: zinner, zouter, boxval
  INTEGER, INTENT(IN) :: tval
  INTEGER, ALLOCATABLE, DIMENSION(:)::dum_aid, dum_typ 

  ! Find boxval
  IF(major_axis == 1) boxval = box_arr(tval,1)
  IF(major_axis == 2) boxval = box_arr(tval,2)
  IF(major_axis == 3) boxval = box_arr(tval,3)
  
  ! delz = (zi-zmin)/n - (n-1)delta_l/n
  IF(tval == 1) THEN
     ALLOCATE(seg_typcnt(1:ntotatomtypes),stat=AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate seg_typcnt"
  END IF

  IF(zmin < 0.0 .OR. zmax > boxval) THEN
     PRINT *, zmin, zmax, 0.0, boxval
     PRINT *, "zmin/zmax exceeds box limits"
     STOP 
  END IF

  IF(maxinterf .NE. 1) STOP "layer_group_surf will work iff there is o&
       &nly 1 interface"
  
  ! based on zmin
  delz_1 = (1.0/REAL(nmax_layers))*(interpos(1) - zmin) -&
       & (REAL(nmax_layers - 1)/REAL(nmax_layers))*deltaz_btw_layers

  ! based on zmax
  delz_2 = (1.0/REAL(nmax_layers))*(zmax - interpos(1) ) -&
       & (REAL(nmax_layers - 1)/REAL(nmax_layers))*deltaz_btw_layers

  ALLOCATE (dum_aid(ntotatoms), stat = AllocateStatus)
  IF(AllocateStatus /=0 ) STOP "Allocation dum_aid failed..."
  ALLOCATE (dum_typ(ntotatoms), stat = AllocateStatus)
  IF(AllocateStatus /=0 ) STOP "Allocation dum_typ failed..."

  !Identify the domain ID- 1(A), 2(B)
  ALLOCATE (seg_dtype(ntotatoms), stat = AllocateStatus)
  IF(AllocateStatus /=0 ) STOP "Allocation seg_dtype failed..."
  CALL ASSIGN_DOMAINID(tval)

  layer_deltaz = MIN(delz_1, delz_2) ! Take the min of two
  zinner = zmin
  seg_typcnt = 0

  IF(tval == 1) PRINT *, "Initial layer thickness: ", layer_deltaz

  DO p = 1,2*nmax_layers !nmax_layers is per domain

     dum_aid = -1 !initialize everything to -1
     dum_typ = -1 !initialize everything to -1
     seg_typcnt = 0 !initialize all counts to 0     
     zouter = zinner + layer_deltaz

     !Sanity check
     IF(p == nmax_layers) THEN
        IF(ABS(zouter - interpos(1)) > 0.1) THEN
           PRINT *, "Something wrong with binning at tval: ", tval
           PRINT *, zinner, zouter, interpos(1), layer_deltaz, p
           STOP
        END IF
     END IF

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

        ! Add to bin if they are between zinner and zouter
        IF(par_pos .GE. zinner .AND. par_pos .LT. zouter) THEN
             
           flagbin = flagbin + 1
           k = k + 1
           dum_aid(k) = a1id
           dum_typ(k) = aidvals(a1id,3)
           seg_typcnt(aidvals(a1id,3)) = seg_typcnt(aidvals(a1id,3))+1
                
           IF(tval == 1) WRITE(dumwrite,'(2(I0,1X),3(F14.8,1X))') i,&
                & aidvals(a1id,3),trx_lmp(a1id,tval)&
                &,try_lmp(a1id,tval),trz_lmp(a1id,tval)

           
        END IF
     
     END DO

     IF(tval == 1) THEN
        WRITE(logout,*) zmin, zmax, nmax_layers, boxval
        WRITE(logout,*) "Width considered for segment: ", layer_deltaz
        WRITE(logout,*) "Population of each type"
        DO i = 1,ntotatomtypes
           WRITE(logout,*) i,seg_typcnt(i)
        END DO
     END IF
          
     num_mons_per_layer = k

     
     WRITE(logout,*) "Number of particles in layer# ", p ," coorespond&
          &ing to ", zinner, " < z < ", zouter, "at tval =  ", tval, "&
          &is ", num_mons_per_layer

     IF(num_mons_per_layer == 0) THEN
        PRINT *, "No particles found in layer between ", zinner, " < z&
             & < ",zouter
        CYCLE
     END IF
     
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

     CALL LAYERWISE_ANALYSIS(tval,p,num_mons_per_layer,layer_deltaz)

     DEALLOCATE(seg_aid)
     DEALLOCATE(seg_typ)
     
     zinner = zouter + deltaz_btw_layers
     
  END DO

  IF(tval == 1) PRINT *, "Total number of sub-layers: ", p-1
  DEALLOCATE(dum_aid)
  DEALLOCATE(dum_typ)
  DEALLOCATE(seg_dtype)

  
END SUBROUTINE COMPARTMENTALIZE_PARTICLES_BOTTOMSURF

!--------------------------------------------------------------------

SUBROUTINE ASSIGN_DOMAINID(tval)

  USE STATICPARAMS
  IMPLICIT NONE

  REAL :: par_pos
  INTEGER :: i, flagdom, domnum, boxval
  INTEGER, INTENT(IN) :: tval

  ! Find boxval
  IF(major_axis == 1) boxval = box_arr(tval,1)
  IF(major_axis == 2) boxval = box_arr(tval,2)
  IF(major_axis == 3) boxval = box_arr(tval,3)

  ! Zero domain types
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

SUBROUTINE LAYERWISE_ANALYSIS(tval,ipos,num_mons,segwidth)

  USE SUBROUTINE_DEFS
  USE STATICPARAMS

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: tval,ipos,num_mons
  REAL, INTENT(IN) :: segwidth
  INTEGER :: t1,t2,clock_rate,clock_max

  IF(rdf2dflag .AND. tval == 1) THEN
     IF(ipos == 1) THEN
        rdf2darray = 0.0
     END IF
     CALL SYSTEM_CLOCK(t1,clock_rate,clock_max)
     CALL RDF2D_LAYER(tval,ipos,num_mons,segwidth)
     CALL SYSTEM_CLOCK(t2,clock_rate,clock_max)
     PRINT *, 'Elapsed real time for 2D RDF= ',REAL(t2-t1)/&
          & REAL(clock_rate), " seconds"           
  ELSEIF(rdf2dflag .AND. mod(tval,rdf2dfreq) == 0) THEN
     CALL RDF2D_LAYER(tval,ipos,num_mons,segwidth)
  END IF

END SUBROUTINE LAYERWISE_ANALYSIS

!--------------------------------------------------------------------

SUBROUTINE RDF2D_LAYER(tval,ipos,num_mons,segwidth)

  USE STATICPARAMS
  IMPLICIT NONE
  
  INTEGER :: i,j,a1type,a2type,a1id,a2id,ibin
  INTEGER :: reftype, seltype, paircnt, selcnt
  REAL :: rxval,ryval,rzval,rval,normfac,rdf2dbinval,rvolval
  REAL, INTENT(IN) :: segwidth
  INTEGER, INTENT(IN) :: tval, ipos, num_mons
  INTEGER,DIMENSION(0:rdf2dmaxbin-1,npairs_2drdf) :: dumrdfarray

  rvolval = box_xl*box_yl*box_zl
  IF(ipos == 1) rdf2dvolavg = rdf2dvolavg + rvolval
  rdf2dbinval = 0.5*segwidth/REAL(rdf2dmaxbin)
  IF(ipos == 1) rdf2dbinavg = rdf2dbinavg + rdf2dbinval

!$OMP PARALLEL DO PRIVATE(i,j)
  DO j = 1,npairs_2drdf
     DO i = 0,rdf2dmaxbin-1
        dumrdfarray(i,j) = 0.0
     END DO
  END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL 
!$OMP DO PRIVATE(paircnt,i,j,reftype,seltype,a1type,a2type,&
!$OMP& a1id,a2id,rval,rxval,ryval,rzval,ibin) REDUCTION(+:dumrdfarray)

  DO paircnt = 1,npairs_2drdf

     reftype = pairs_2drdf_arr(paircnt,1)
     seltype = pairs_2drdf_arr(paircnt,2)

     IF(.NOT. ANY(seg_typ == reftype)) CYCLE
     IF(.NOT. ANY(seg_typ == seltype)) CYCLE
     
     DO i = 1,num_mons

        a1id   = seg_aid(i)     
        a1type = seg_typ(i)

        IF(a1type .NE. aidvals(a1id,3)) THEN
           PRINT *, "1", i, a1type, a1id, aidvals(a1id,3)
           STOP "ERROR: 2D-RDF Invalid i type found in 2DRDF"
        END IF

        IF(a1type .NE. reftype) CYCLE
        
        DO j = 1,num_mons
        
           a2id   = seg_aid(j)
           a2type = seg_typ(j)

           IF(a2type .NE. aidvals(a2id,3)) THEN
              PRINT *, "2",j, a2type, a2id, aidvals(a2id,3)
              STOP "ERROR: 2D-RDF Invalid j type found in 2DRDF"
           END IF

           IF(a2type .NE. seltype) CYCLE
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
!$OMP END PARALLEL

  
!$OMP PARALLEL DO PRIVATE(i,j,reftype,seltype,selcnt,normfac) 
  DO j = 1,npairs_2drdf

     reftype = pairs_2drdf_arr(j,1)
     seltype = pairs_2drdf_arr(j,2)
     selcnt  = pairs_2drdf_arr(j,3)
     normfac = rvolval/(REAL(selcnt)*REAL(seg_typcnt(reftype)))

     DO i = 0,rdf2dmaxbin-1

        rdf2darray(i,j,ipos) = rdf2darray(i,j,ipos) + REAL(dumrdfarray(i&
             &,j))*normfac

     END DO
     
  END DO
!$OMP END PARALLEL DO


END SUBROUTINE RDF2D_LAYER

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_RDF()

  USE STATICPARAMS
  IMPLICIT NONE

  INTEGER :: i,j,a1type,a2type,ibin,a1id,a2id,paircnt,AllocateStatus
  REAL :: rxval,ryval,rzval,rval,rvolval
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

  USE SUBROUTINE_DEFS
  USE STATICPARAMS
  IMPLICIT NONE

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

  INTEGER :: p,i,paircnt
  REAL :: zin, zout, zin_wrt_interf, zout_wrt_interf
  REAL, PARAMETER :: vconst = 4.0*pival/3.0
  REAL :: rlower,rupper,nideal,rdf2dfrnorm
  CHARACTER(LEN=3) :: intnum
  
  IF(rdf2dfreq == 1) THEN
     rdf2dfrnorm = nframes
  ELSE
     rdf2dfrnorm = INT(nframes/rdf2dfreq)+1
  END IF

  rdf2dvolavg = rdf2dvolavg/REAL(rdf2dfrnorm)
  rdf2dbinavg = rdf2dbinavg/REAL(rdf2dfrnorm)
  PRINT *, "rdf2d-volavg/rdf2d-binavg: ", rdf2dvolavg,rdf2dbinavg

  IF(layer_grpflag_interf) THEN
     DO p = 1,nmax_layers
     
        WRITE(intnum,'(I0)') p
        dum_fname  = "rdf2d_interfbased_"//trim(intnum)//"_"&
             &//trim(adjustl(traj_fname))//".txt"
        
        OPEN(unit = dumwrite,file = trim(dum_fname), status="replace"&
             &, action = "write")
        
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

  ELSEIF(layer_grpflag_surf) THEN

     zin = zmin
     DO p = 1,2*nmax_layers
     
        WRITE(intnum,'(I0)') p
        dum_fname  = "rdf2d_surfbased_"//trim(intnum)//"_"&
             &//trim(adjustl(traj_fname))//".txt"
        
        OPEN(unit = dumwrite,file = trim(dum_fname), status="replace"&
             &, action = "write")

        zout = zin + 2*rdf2dbinavg*rdf2dmaxbin !Factor of 2 because
        !rdfbin is defined based on half of the segment width

        IF(zin < interpos(1)) THEN
           zin_wrt_interf = interpos(1) - zout
        ELSE
           zin_wrt_interf  = zin - interpos(1)
        END IF
        zout_wrt_interf = zin_wrt_interf + 2*rdf2dbinavg*rdf2dmaxbin
           
        WRITE(dumwrite,'(F12.6,2X,A7,2X,F12.6)') zin, " < z < ", zout
        WRITE(dumwrite,'(F12.6,2X,A14,2X,F12.6)') zin_wrt_interf, " < &
             & |z-z_i| < ", zout_wrt_interf
        
        
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

        zin = zout + deltaz_btw_layers
        
     END DO

  END IF
  
END SUBROUTINE OUTPUT_LAYERRDF

!--------------------------------------------------------------------

SUBROUTINE OUTPUT_ALLRDF()

  USE STATICPARAMS
  IMPLICIT NONE

  INTEGER :: i,j,ierr
  REAL, PARAMETER :: vconst = 4.0*pival/3.0
  REAL :: rlower,rupper,nideal,rdffrnorm
  
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

  PRINT *, "Successfully allocated memory for topology"

END SUBROUTINE ALLOCATE_TOPO_ARRAYS

!--------------------------------------------------------------------

SUBROUTINE ALLOCATE_ANALYSIS_ARRAYS()

  USE SUBROUTINE_DEFS
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
     IF(AllocateStatus/=0) STOP "did not allocate trx_lmp"
     ALLOCATE(try_lmp(ntotatoms,nframes),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate try_lmp"
     ALLOCATE(trz_lmp(ntotatoms,nframes),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate trz_lmp"
  ELSE
     ALLOCATE(trx_lmp(1,nframes),stat = AllocateStatus)
     ALLOCATE(try_lmp(1,nframes),stat = AllocateStatus)
     ALLOCATE(trz_lmp(1,nframes),stat = AllocateStatus)
     DEALLOCATE(trx_lmp)
     DEALLOCATE(try_lmp)
     DEALLOCATE(trz_lmp)
  END IF
     
  IF(rdf2dflag) THEN
     IF(layer_grpflag_interf) THEN
        ALLOCATE(rdf2darray(0:rdf2dmaxbin-1,npairs_2drdf,nmax_layers)&
             &,stat=AllocateStatus)
        IF(AllocateStatus/=0) STOP "did not allocate rdf2darray"
     ELSEIF(layer_grpflag_surf) THEN
        ALLOCATE(rdf2darray(0:rdf2dmaxbin-1,npairs_2drdf,2&
             &*nmax_layers),stat=AllocateStatus)
        IF(AllocateStatus/=0) STOP "did not allocate rdf2darray"
     END IF
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

  USE SUBROUTINE_DEFS
  USE STATICPARAMS

  IMPLICIT NONE

  !Global arrays
  DEALLOCATE(aidvals)
  DEALLOCATE(rxyz_lmp)
  DEALLOCATE(vel_xyz)
  DEALLOCATE(box_arr)

  !Topo arrays
  IF(ntotbonds /= 0) DEALLOCATE(bond_lmp)
  IF(ntotangls /= 0) DEALLOCATE(angl_lmp)
  IF(ntotdihds /= 0) DEALLOCATE(dihd_lmp)
  IF(ntotimprs /= 0) DEALLOCATE(impr_lmp)
  
  !Statics calculations arrays
  IF(rdfflag) DEALLOCATE(rdfarray)
  IF(densflag) DEALLOCATE(densarray)

  !Dellocate groups
  IF(grpflag) THEN
     DEALLOCATE(allgrp_atomarr)
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
