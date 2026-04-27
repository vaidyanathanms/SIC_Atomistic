!---------------To analyze properties of bulk-sei systems------------
!---------------Version 2: Dev_Apr-27-2026---------------------------
!---------------Main File: analyze_statics.f90-----------------------
!********************************************************************

MODULE SUBROUTINE_DEFS

  USE OMP_LIB
  IMPLICIT NONE

  ! Read all keywords
  INTERFACE
     SUBROUTINE READ_ANA_IP_FILE()
       IMPLICIT NONE
     END SUBROUTINE READ_ANA_IP_FILE
  END INTERFACE

  ! Parse keywords
  INTERFACE
     SUBROUTINE READ_NEXT_KEYWORD(iu, keyword, ierr)
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: iu
       CHARACTER(LEN=*), INTENT(OUT) :: keyword
       INTEGER, INTENT(out) :: ierr
     END SUBROUTINE READ_NEXT_KEYWORD
  END INTERFACE

  ! Set default values
  INTERFACE
     SUBROUTINE DEFAULTVALUES()
       IMPLICIT NONE
     END SUBROUTINE DEFAULTVALUES
  END INTERFACE

  ! Read LAMMPS datafile
  INTERFACE
     SUBROUTINE READ_DATAFILE()
       IMPLICIT NONE
     END SUBROUTINE READ_DATAFILE
  END INTERFACE

  ! Read # of lines in datafile before masses keyword
  INTERFACE
     SUBROUTINE COMPUTE_INIT_NLINES(imax)
       IMPLICIT NONE
       INTEGER, INTENT(OUT) :: imax
     END SUBROUTINE COMPUTE_INIT_NLINES
  END INTERFACE

  ! Generate arrays for types in each GLOBAL group
  INTERFACE
     SUBROUTINE FILL_GROUP_ARRAY(group_id,ntypes_per_group&
          &,types_in_group_id)
       IMPLICIT NONE
       INTEGER, INTENT(IN)::group_id, ntypes_per_group
       INTEGER, INTENT(IN)::types_in_group_id(1:ntypes_per_group)
     END SUBROUTINE FILL_GROUP_ARRAY
  END INTERFACE

  ! Generate arrays for types in each layer group
  INTERFACE
     SUBROUTINE FILL_LAYERGROUP_ARRAY(group_id,ntypes_per_group&
          &,types_in_group_id)
       IMPLICIT NONE
       INTEGER, INTENT(IN)::group_id, ntypes_per_group
       INTEGER, INTENT(IN)::types_in_group_id(1:ntypes_per_group)
     END SUBROUTINE FILL_LAYERGROUP_ARRAY
  END INTERFACE
  
  ! Output group details
  INTERFACE
     SUBROUTINE OUTPUT_ALL_GROUPS()
       IMPLICIT NONE
     END SUBROUTINE OUTPUT_ALL_GROUPS
  END INTERFACE

  ! Count number of atoms of i-th type
  INTERFACE
     SUBROUTINE COUNT_ATOMS_WITH_TYPE_I(inptype,outcnt)
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: inptype
       INTEGER, INTENT(OUT) :: outcnt
     END SUBROUTINE COUNT_ATOMS_WITH_TYPE_I
  END INTERFACE

  ! Analyze LAMMPS trajectory file
  INTERFACE
     SUBROUTINE ANALYZE_TRAJECTORYFILE()
       IMPLICIT NONE
     END SUBROUTINE ANALYZE_TRAJECTORYFILE
  END INTERFACE

  ! Initialize structural analysis details for RDF
  INTERFACE
     SUBROUTINE STRUCT_INIT()
       IMPLICIT NONE
     END SUBROUTINE STRUCT_INIT
  END INTERFACE

  ! Call to main unit of structural analysis
  INTERFACE
     SUBROUTINE STRUCT_MAIN(tval)
       IMPLICIT NONE
       INTEGER, INTENT(IN):: tval
     END SUBROUTINE STRUCT_MAIN
  END INTERFACE

  ! Computing density profile
  INTERFACE
     SUBROUTINE COMPUTE_DENSPROFILES(tval)
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: tval
     END SUBROUTINE COMPUTE_DENSPROFILES
  END INTERFACE

  ! Sort ions and counterion IDs and types
  INTERFACE
     SUBROUTINE SORTALLARRAYS()
       IMPLICIT NONE
     END SUBROUTINE SORTALLARRAYS
  END INTERFACE

  ! Sanity check for iontype definitions and counts
  INTERFACE
     SUBROUTINE SANITY_CHECK_IONTYPES()
       IMPLICIT NONE
     END SUBROUTINE SANITY_CHECK_IONTYPES
  END INTERFACE

  ! Used to map the atomtype to the sortedarray row
  INTERFACE
     SUBROUTINE MAP_REFTYPE(jin,atype,jout)
       IMPLICIT NONE
       INTEGER, INTENT(IN):: jin,atype
       INTEGER, INTENT(OUT) :: jout
     END SUBROUTINE MAP_REFTYPE
  END INTERFACE

  ! To compute # of anions around each cation
  INTERFACE
     SUBROUTINE CAT_AN_NEIGHS()
       IMPLICIT NONE
     END SUBROUTINE CAT_AN_NEIGHS
  END INTERFACE

  ! Call to main unit for layerwise analysis
  INTERFACE
     SUBROUTINE LAYERWISE_MAIN()
       IMPLICIT NONE
     END SUBROUTINE LAYERWISE_MAIN
  END INTERFACE

  ! Computing interface positions
  INTERFACE
     SUBROUTINE COMPUTE_INTERFACES
       IMPLICIT NONE
     END SUBROUTINE COMPUTE_INTERFACES
  END INTERFACE

  ! Mapping group ID to the column in grp_array
  INTERFACE
     SUBROUTINE MAP_GROUP_TO_COL(grp_type,col_val)
       IMPLICIT NONE
       INTEGER, INTENT(IN)  :: grp_type
       INTEGER, INTENT(OUT) :: col_val
     END SUBROUTINE MAP_GROUP_TO_COL
  END INTERFACE

  ! Initial domain classificaiton for nonperiodic BC
  INTERFACE
     SUBROUTINE INITIAL_DOMAIN_CLASSIFICATION(tval)
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: tval
     END SUBROUTINE INITIAL_DOMAIN_CLASSIFICATION
  END INTERFACE

  ! Bin particles based on INTERFACIAL position
  INTERFACE
     SUBROUTINE COMPARTMENTALIZE_PARTICLES_INTERFACE(tval)
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: tval
     END SUBROUTINE COMPARTMENTALIZE_PARTICLES_INTERFACE
  END INTERFACE

  ! Bin particles based on BOTTOM SURFACE
  INTERFACE
     SUBROUTINE COMPARTMENTALIZE_PARTICLES_BOTTOMSURF(tval)
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: tval
     END SUBROUTINE COMPARTMENTALIZE_PARTICLES_BOTTOMSURF
  END INTERFACE

  ! Assign domain IDs as either 1 or 2
  INTERFACE
     SUBROUTINE ASSIGN_DOMAINID(tval)
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: tval
     END SUBROUTINE ASSIGN_DOMAINID
  END INTERFACE

  ! Main call to doing layer-wise analysis
  INTERFACE
     SUBROUTINE LAYERWISE_ANALYSIS(tval,ipos,num_mons,segwidth)
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: tval,ipos,num_mons
       REAL, INTENT(IN) :: segwidth
     END SUBROUTINE LAYERWISE_ANALYSIS
  END INTERFACE

  ! 2D RDF analysis
  INTERFACE
     SUBROUTINE RDF2D_LAYER(tval,ipos,num_mons,segwidth)
       IMPLICIT NONE
       REAL, INTENT(IN) :: segwidth
       INTEGER, INTENT(IN) :: tval, ipos, num_mons
     END SUBROUTINE RDF2D_LAYER
  END INTERFACE

  ! Global RDF analysis - bulk
  INTERFACE
     SUBROUTINE COMPUTE_RDF()
       IMPLICIT NONE
     END SUBROUTINE COMPUTE_RDF
  END INTERFACE

  ! Main call to opening output files for structural analysis
  INTERFACE
     SUBROUTINE OPEN_STRUCT_OUTPUT_FILES()
       IMPLICIT NONE
     END SUBROUTINE OPEN_STRUCT_OUTPUT_FILES
  END INTERFACE

  ! Generic call to printing statistics
  INTERFACE
     SUBROUTINE ALLOUTPUTS()
       IMPLICIT NONE
     END SUBROUTINE ALLOUTPUTS
  END INTERFACE

  ! Ouput densities
  INTERFACE
     SUBROUTINE OUTPUT_ALLDENS()
       IMPLICIT NONE
     END SUBROUTINE OUTPUT_ALLDENS
  END INTERFACE

  ! Output layer-wise 2D RDF
  INTERFACE
     SUBROUTINE OUTPUT_LAYERRDF()
       IMPLICIT NONE
     END SUBROUTINE OUTPUT_LAYERRDF
  END INTERFACE

  ! Output global RDF
  INTERFACE
     SUBROUTINE OUTPUT_ALLRDF()
       IMPLICIT NONE
     END SUBROUTINE OUTPUT_ALLRDF
  END INTERFACE

  ! Allocate topology arrays (atom/bond/angle/dihedral/improper)
  INTERFACE
     SUBROUTINE ALLOCATE_TOPO_ARRAYS()
       IMPLICIT NONE
     END SUBROUTINE ALLOCATE_TOPO_ARRAYS
  END INTERFACE

  ! Allocate all analysis arrays
  INTERFACE
     SUBROUTINE ALLOCATE_ANALYSIS_ARRAYS()
       IMPLICIT NONE
     END SUBROUTINE ALLOCATE_ANALYSIS_ARRAYS
  END INTERFACE

  ! Global deallocation
  INTERFACE
     SUBROUTINE DEALLOCATE_ARRAYS()
       IMPLICIT NONE
     END SUBROUTINE DEALLOCATE_ARRAYS
  END INTERFACE
       
END MODULE SUBROUTINE_DEFS
