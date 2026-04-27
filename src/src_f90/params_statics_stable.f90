!---------------To analyze properties of bulk-sei systems------------
!---------------Version 2: Apr-24-2026-------------------------------
!---------------Main File: analyze_statics.f90-----------------------
!********************************************************************

MODULE STATICPARAMS

  USE OMP_LIB
  IMPLICIT NONE

  !General flags
  INTEGER :: readdataflag
  
  !All bulk structure analysis flags
  INTEGER :: rgcalc, rgall, rgavg
  INTEGER :: catan_neighcalc
  INTEGER :: bfrdf_calc
  INTEGER :: clust_calc

  !All bulk dynamics analysis flags
  INTEGER :: ion_dynflag, cion_dynflag, pion_dynflag
  INTEGER :: ion_diff, cion_diff, pion_diff
  INTEGER :: catan_autocfflag, catpol_autocfflag

  !All layerwise analysis flags - structure
  INTEGER :: grpflag, densflag,rdfflag,interfflag
  INTEGER :: layerana_flag
  INTEGER :: layer_grpflag_interf, layer_grpflag_surf
  INTEGER :: layrdf_flag
  INTEGER :: rdf2dflag

  !Required Input Variables
  INTEGER :: initdist
  INTEGER :: nframes, skipfr, freqfr, nfrcntr
  INTEGER :: nchains, atperchain
  INTEGER :: nproc

  !Variables required for computing bulk static properties
  INTEGER :: rdffreq,rmaxbin,npairs,d_maxbin,densfreq
  INTEGER :: rgfreq, rdfpaircnt
  REAL    :: rvolavg,rdomcut,rbinval
  REAL    :: boxlzavg,ddenavg
  REAL    :: re2ave, re4ave, rg2ave, rg4ave, b2ave
  REAL    :: delta_t
  REAL    :: rneigh_cut,rcatan_cut,rcatpol_cut1,rcatpol_cut2

  !Variables required for computing layerwise static properties
  INTEGER :: major_axis
  INTEGER :: rdf2dfreq,rdf2dmaxbin,npairs_2drdf
  REAL    :: dbinavg, major_boxval
  REAL    :: epspre, epsinit, segper, epsinc
  REAL    :: rdf2dvolavg,rdf2dbinavg
  REAL    :: zmin, zmax, deltaz_btw_layers
  
  !Global group details
  INTEGER :: ngroups, ndens_grps
  INTEGER :: ioncnt,  iontype
  INTEGER :: c_ioncnt, c_iontype
  INTEGER :: maxneighsize, neighfreq
  INTEGER :: ntotion_centers

  !Layerwise group details
  INTEGER :: interfgrp_a, interfgrp_b
  INTEGER :: maxinterf, nmax_layers
  INTEGER :: num_mons_per_layer, nlayer_groups  
    
  !File names and unit numbers
  CHARACTER(LEN = 256) :: ana_fname,data_fname,traj_fname,log_fname
  CHARACTER(LEN = 256) :: rdf_fname, dum_fname
  INTEGER, PARAMETER :: anaread = 21,   logout = 3
  INTEGER, PARAMETER :: dumwrite = 50
  INTEGER, PARAMETER :: inpread = 100, rgwrite = 200,rgavgwrite = 300

  !Math Constants
  REAL*8, PARAMETER :: pival  = 3.14159265359
  REAL*8, PARAMETER :: pi2val = 2.0*pival

  !Global analysis variables and arrays
  INTEGER :: atomflag, velflag, bondflag, anglflag, dihdflag,imprflag
  INTEGER :: ntotatoms, ntotbonds, ntotangls,ntotdihds,ntotimprs
  INTEGER :: ntotatomtypes,ntotbondtypes,ntotangltypes,ntotdihdtypes&
       &,ntotimprtypes

  !Lammps trajectory file read details
  REAL :: box_xl,box_yl,box_zl, boxval
  INTEGER*8 :: timestep

  !Required Arrays - LAMMPS
  REAL*8, ALLOCATABLE, DIMENSION(:,:) :: rxyz_lmp, vel_xyz, charge_lmp&
       &,masses
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: bond_lmp, angl_lmp,&
       & dihd_lmp, impr_lmp,aidvals
  CHARACTER,ALLOCATABLE,DIMENSION(:,:) :: keywords
  REAL,ALLOCATABLE,DIMENSION(:,:) :: box_arr

  !General required arrays required for computing properties
  INTEGER,ALLOCATABLE,DIMENSION(:)   :: types_in_group
  INTEGER,ALLOCATABLE,DIMENSION(:,:) :: allgrp_typarr, allgrp_atomarr
  INTEGER,ALLOCATABLE,DIMENSION(:,:) :: ionarray,counterarray
  INTEGER,ALLOCATABLE,DIMENSION(:,:) :: allionids
  REAL,ALLOCATABLE,DIMENSION(:) :: clust_avg

  !Required Arrays - Bulk structural Quantities
  REAL,ALLOCATABLE,DIMENSION(:,:):: rdfarray, densarray
  REAL,ALLOCATABLE,DIMENSION(:) :: cat_an_neighavg,an_cat_neighavg
  REAL,ALLOCATABLE,DIMENSION(:) :: rdf_p_fb,rdf_p_bb,rdf_p_ff
  INTEGER,ALLOCATABLE,DIMENSION(:,:):: pairs_rdf
  INTEGER,ALLOCATABLE,DIMENSION(:) :: polboundarr,polfreearr

  !Required Arrays - Layerwise analysis
  REAL, ALLOCATABLE,DIMENSION(:)    :: interpos, widdoms
  INTEGER,ALLOCATABLE,DIMENSION(:,:):: all_layergrp_typarr
  INTEGER,ALLOCATABLE,DIMENSION(:)  :: domtyp
  INTEGER,ALLOCATABLE,DIMENSION(:)  :: seg_aid, seg_typ,seg_typcnt
  INTEGER,ALLOCATABLE,DIMENSION(:)  :: seg_dtype
  
  !Required Arrays - Layer-wise calculations
  REAL,ALLOCATABLE,DIMENSION(:,:,:) :: rdf2darray
  INTEGER,ALLOCATABLE,DIMENSION(:,:):: pairs_2drdf_arr


  !Required Arrays - Dynamic Quantities
  INTEGER*8,ALLOCATABLE,DIMENSION(:) :: tarr_lmp
  REAL*8,ALLOCATABLE,DIMENSION(:,:) :: trx_lmp,try_lmp,trz_lmp

END MODULE STATICPARAMS

!--------------------------------------------------------------------
