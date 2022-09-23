
                Chaneglog for Phenix distribution
                =================================

1.20.1 (January 2022)
=====================

- Add backwards compatibility for new solvent masking algorithm
- Bug fix SHELX HKLF format output
- Bug fix for phenix.dock_and_rebuild where no model is obtained
- Bug fix for map and model from phenix.douse not aligning
- Add --without-dials option to installation script

1.20 (December 2021)
====================

- New tools and methods
  - Phenix AlphaFold2 notebook: Run AlphaFold on Google Colab from Phenix GUI
  - phenix.process_predicted_model: Identify useful domains in AlphaFold model
  - phenix.dock_predicted_model: Dock domains of AlphaFold model into cryo-EM
  - phenix.rebuild_predicted_model: Rebuild AlphaFold model in cryo-EM map using docked domains
  - phenix.dock_and_rebuild : Process, dock and rebuild AlphaFold model with cryo-EM map
  - phenix.model_completion: Connect fragments and fill in gaps based on a map
  - phenix.rebuild_model: Rebuild a model using a map and keeping connectivity
  - phenix.replace_with_fragments_from_pdb: Rebuild a model using fragments from PDB
  - phenix.search_and_morph: SSM search PDB; morph to match target
  - phenix.fragment_search: Search for a fragment in PDB matching target
  - phenix.reverse_fragment: Reverse chain direction of a fragment
  - phenix.superpose_and_morph: SSM or least-squares superpose one model on another; optionally trim and morph to match
  - phenix.voyager.casp_rel_ellg: Calculate relative eLLG score for predicted model quality

- phenix.match_maps:
  - Bug fix (superposed map was not matching target map)

- phenix.real_space_refine:
  - Symmatry multiprocessing aware individual ADP and occupancy refinement
  - Multiple changes to improce runtime (for certain refinement strategies)
  - Make NQH flips symmetry aware

- phenix.superpose_pdbs:
  - Add feature to transform additional models with matrix found with moving model

- phenix.dock_in_map:
  - Allow splitting model into domains based on chain ID from phenix.process_predicted_model

- Restraints
  - GeoStd updated with 12k plus entity restraints files
  - cis-PRO default updated to EH99

- phenix.fetch_pdb, iotbx.cif_as_mtz:
  - Bug fix: Multiple datasets with different unit cells in a cif file now preserved
    as multiple crystals in mtz file.

1.19.2 (February 2021)
======================

- Avoid GUI crash when running validation after phenix.real_space_refine
- Other bug fixes

1.19.1 (January 2021)
=====================

- phenix.real_space_refine:
  - Fix bug in restoring previous jobs
  - Clarify NCS options

1.19 (December 2020)
=====================

- phenix.real_space_refine:
  - Improved rotamer fitting (use multiprocessing, in case of NCS constraints
    work on one copy only and propagate changes to all related copies, various
    performace improvements and bug fixes)
  - Improved map/restraints weight calculation
  - Morphing can now use multiprocessing (nproc)

- phenix.map_to_model and phenix.trace_and_build
  - Improved high-resolution model-building including detection of insertions
    and deletions

- New methods
  - phenix.local_resolution: calculates a local resolution map
  - phenix.local_aniso_sharpen: optimizes a map taking into account local
      resolution-dependence and anisotropy of the map and its errors

- New scripting tools
  - High-level scriptable Python tools are now available for map and model
     analyses and manipulation as well as for model-building

- Restraints
  - Adjusting the "positions" of atom names is (pseudo-)symmetric amino acid
    side chains is now the default
  - Improved restraints for ARG allows more flexibility of the CD atom

- Bug fixes
  - Chains with modified amino or nucleic acids can be aligned
  - Fixed handling of modified amino/nucleic acids in sequence alignment


1.18.2 (May 2020)
=================

- Fix weighting for iron sulfur clusters

1.18.1 (May 2020)
=================

- General bug fixes

- GUI
  - Add default location of Coot from CCP4 7.1 on macOS
  - Extra column for chiral outliers indicating probable cause

1.18 (April 2020)
=================

- Amber
  - phenix.refine GUI will automatically create Amber files if use_amber is
    checked and files not provided
  - AmberPrep GUI added for Amber file preparation

- Restraints
  - Engh & Huber restraints (2001) for cis-peptides now implemented
  - Restraints added for FES among others
  - Metal coordination library (MCL) is now the default for Zn2+ coordination
    as well as coordination of Iron-Sulfur clusters

- Rama-Z: A global Ramachandran score identifies protein structures with
  unlikely stereochemistry.
  - Command-line: phenix.rama_z, mmtbx.rama_z
  - GUI: validation reports

- Density modification for cryo-EM (ResolveCryoEM) now available in GUI
  - Includes model-based density modification with automatic model
    generation
  - Simplified GUI showing just the most important parameters
  - Optimized defaults
  - Additional documentation

- New tool for ordered solvent (water) picking into cryo-EM maps: phenix.douse
  (command line and GUI)

- Real-space refinement
  - Do not include H atoms into map target function. This improves the fit.
  - Add NQH flip option (enabled by default).

- Rotamer fitting in real- and reciprocal-space refinement
  - Multiple bug fixes and performance enhancements.
  - Enable multiprocessing (nproc).
  - Make NCS constraints aware (fit rotamers only in master copies).

- New map and model superposition tool: phenix.match_maps
  - Supports different unit cells and space groups of moving and fixed models.
  - Does not move model and map into an abstract frame of reference: moving
    model and map overlay with the fixed model and map.
  - Supports reflection data (crystallography) and real-space maps (cryo-EM).

- GUI
  - Added button to Project settings dialog for showing README.txt
  - Added FindProgram tool to find any phenix program with a text search and
    show how to access it in the GUI (if available)

1.17 (September 2019)
=====================

- phenix.reflection_file_converter
  - handling of SHELX output has been improved, with the output data type
    (intensities or amplitudes) matching the input instead of being converted
    to amplitudes: intensities are preferred for SHELX programs. The default
    filename extension is now .hkl, as expected by SHELX programs, but .shelx
    is still allowed for compatibility with existing scripts. Output formatted
    to preserve dynamic range of input and rescaled if necessary.

- Amber
  - Amber is now distributed with Phenix for immediate use in phenix.refine

- eLBOW
  - Optionally outputs files needed for Amber
  - Added support for QM package Orca

- dials.image_viewer is used for viewing diffraction images
  - Crashes may occur from unsupported file formats

- Bug fixes
  - Updated map smoothing
  - phenix.validation_cryoem fix inconsistency in clashscore values due to
    hydrogens being kept by default
  - Consistent behavior and logging of geometry_restraints.edits

1.16 (July 2019)
================

- NEW FEATURES

  - New GUI section for PDB Deposition
    - Updated tool for preparing mmCIF files for PDB deposition
      (mmtbx.prepare_pdb_deposition)
    - New tool to get PDB validation report (phenix.get_pdb_validation_report)
      from the PDB OneDep Validation web service.

  - Comprehensive Validation for cryo-EM now performs the same sequence check
    as Comprehensive Validation for X-ray/Neutron.

- Updated dependencies
  - Added amber_phenix
  - Added onedep_api

1.15 (February 2019)
====================

- NEW FEATURES

  - phenix.map_to_model version 2 uses the new trace_and_build algorithm to
      build 10x faster and with half as many chain breaks as version 1.

  - phenix.trace_and_build builds protein into density at resolutions
      as low as 4.5 A by working from high to low contour levels and
      choosing density that is non-branching.

  - phenix.fix_insertions_deletions improves low_resolution models by
      rebuilding at positions where the density is poor

  - phenix.sequence_from_map uses the new trace_and_build algorithm to
      improve identification of what sequence goes with a map.

  - phenix.refine_ca_model - optimizes a CA-only model by adding and
      subtracting CA atoms, matching CA positions to side-chain density and
      adjusting CA-CA distances

  - phenix.comparama - command-like tool to generate Kleywegt-like plots
    for models before and after refinement. Shows how residues moved during
    refinement.

  - eLBOW
    - New script to find unique instances of a ligand in the Protein Data
      Bank. Returns a list of ligands and resolutions via web services.
      Optionally, download and generate maps include a Polder OMIT map and
      load into Coot.

- Output model in PDB and mmCIF format by default by major tools:
  phenix.refine, phenix.real_space_refine, phenix.pdbtools.

- Updated mmCIF support
  - struct_conn loop
  - correct ligand restraints output

- Updated dependencies
  - Move to conda packages for all dependencies
  - Better compatibility with newer operating systems
  - More consistent builds across all platforms (linux, macOS, Windows)

- phenix.structure_search:
  - Improved alignment results in low similarity cases with direct SSM mapping.
  - Updated RCSB non-redundant structural library.

- phenix.ligand_identification:
  - Updated ligand collections for automated library compilation with EC,
    pfam, GO, CATH, SCOP, or InterPro Ids.

- Phaser-2.8.3:
  - bugfixes while we prepare for a new version, phaser_tng

- Validation
  - Improved reporting of cis-peptides for residues with altloc atoms

1.14 (August 2018)
=====================

- NEW FEATURES

  - Full set of cryo-EM tools now available in the GUI: map and model
       quality assessment, map sharpening, finding symmetry in maps,
       cutting out unique parts of maps, combining focused maps,
       docking models into maps, building models into maps,
       identifying sequence from maps, refinement, and validation.

  - New set of tools for interpreting a cryo-EM map if a template model
       is available:
    - phenix.map_symmetry finds the symmetry in your map
    - phenix.map_box with extract_unique=true extracts the unique
       part of the map
    - phenix.dock_in_map places your template model in the map

  - Phenix.map_symmetry:
    - Automatically identifies symmetry in map (C2, D7, helical, I ...)
    - Can find general symmetry operators in addition to point-group or helical
    - Creates Phenix .ncs_spec file with symmetry operators

  - Phenix.map_box:
    - New default is to keep the position of the cut-out density so it
        superimposes on original map (keep_origin=True)
    - New option extract_unique allows extraction of the unique portion
        of a map using supplied symmetry or symmetry operators

  - Phenix.dock_in_map:
    - Automatically identifies optimal position/orientation of a model
        supplied as a PDB file in a CCP4/MRC map.
    - Runs quickly and can use parallel processing

  - Phenix.sequence_from_map:
    - Guesses peptide sequences from a map
    - Identifies sequences that are compatible with a map
    - Carries out sequence assignment for main-chain that has been built

  - Phenix.combine_focused_maps:
    - Combines the best parts of unfocused and focused maps

  - phenix.cablam_idealization - Tool to fix Cablam outliers.

- eLBOW
  - better support for metals and metal clusters
  - added plugin for QM package Orca

- Phaser-2.8.2
  - bugfixes
    - Phassade substructure search when starting from seed substructure
    - fix crash in MR_ATOM
    - problems with cumulative intensity distribution for extremely weak data
  - improve computation and presentation of data information content

- GUI
  - New section in main window for cryo-EM tools
  - Separate validation GUI for cryo-EM structures
  - Added phenix.map_symmetry
  - Added phenix.dock_in_map
  - Added phenix.map_to_structure_factors
  - Added phenix.combine_focused_maps
  - Added phenix.sequence_from_map

- Updated dependencies
  - Python 2.7.14 -> 2.7.15
  - Numpy 1.8.2 -> 1.13.3

- Performance improvements
  - NCS search
  - Generation of secondary structure restraints
  - Clashscore calculation
  - AmberPrep is more robust
  - Restraints for ARG improved
  - Metal cluster DVT added to "superpose ideal ligand" list

- phenix.ligand_identification:
  - Added an option to generate ligand library based on sequence and structural
    homologs of the input pdb model.

1.13 (December 2017)
====================

- NEW FEATURES

  - Phenix.map_to_model:
    - Allows cases where symmetry that is neither point-group nor
       helical is present
    - Allows splitting runs into many small independent steps
    - Set defaults for quick=True to run somewhat more quickly

- GUI
  - Validation information added to phenix.real_space_refine results
  - phenix.real_space_refine GUI accepts input for reference models
  - Polygon stand-alone GUI deprecated; rolled into Comprehensive Validation GUI
    - Old Polygon jobs cannot be restored
  - Validation GUI
    - Stoplight colors (green/orange/red) added to validation summary
    - In rotamer outlier table, uncalculated chi angles now display blank
      instead of None
    - Validation table columns adjusted to reduce whitespace

- Updated dependencies
  - Python 2.7.13 -> 2.7.14
  - Boost 1.56 -> Boost 1.63
  - Numpy 1.8.1 -> 1.8.2

- phenix.structure_search:
  - Updated to November 2017 RCSB non-redundant library
    (representatives of 100% sequence-identity clusters.)
  - Changed to mmcif-based coordinate handling internally.

- phenix.ligand_identification:
  - Added search limit on ligand sizes (Non-H atoms)
  - Run LigandFit algorithm internally as default.

- Phaser-2.8.1
  - bugfixes
    - anisotropy correction of map coefficients after MR
    - detection of duplicate solutions
    - correct cell for PDB after refinement of cell for input data
    - packing places molecules near origin in compact arrangement

- Structure Comparison
  - Updated to include favoured, allowed and outliers for validations tabs
    such as rotamers
  - Add featurs such as ligands information & locations; water locations;
    cis/trans peptides; and histidine protonation
  - General increase of robustness and flexibility
  - Greater flexibility for model/data input including directory tree travesal
    and user defined scripts

- Restraints
  - Added better SF4 restraints
  - Added mechanism to improve SF4 linking

- Various performance enhancements

- Better handling of model mmCIF input/output

- New Phenix video tutorials:
  - Cryo-EM tools in Phenix - overview
  - How to run MolProbity - step by step (web interface and Phenix GUI)
  - MolProbity: All atom contacts tutorial

1.12 (July 2017)
================

- NEW TOOLS

  - phenix.secondary_structure_validation - tool to validate secondary structure
    section in .pdb and .mmcif files with option to filter bad annotations.

  - Structure Comparison: tool for parallel validation and analysis of
    near-identical protein structures. Results include the analysis of ligands,
    rotamers, Ramachandran angles, missing atoms, water locations and B-factors.

  - phenix.mtriage: command-line tool to validate cryo-EM models. GUI is
    available. Computes various map statistics including resolution estimates,
    FSC between half-maps, full and model maps. Outputs mask file used to
    compute FSC. Still under development.

  - phenix.real_space_diff_map: command-line tool to compute difference map.
    This is a real-space analogue of Fo-Fc map mainly designed for cryo-EM maps.
    See http://phenix-online.org/newsletter/CCN_2017_01.pdf for details.


- PHENIX VIDEO TUTORIALS

  Video tutorials about Phenix tools are now available on the Phenix Tutorial
  YouTube channel (www.youtube.com/c/phenixtutorials).

  The following topics are covered:
  - How to run phenix.refine
  - How to calculate polder maps
  - How to run phenix.real_space_refine
  - How to run PDB tools
  - How to change phenix.refine parameters via the GUI
  - Basic Phenix eLBOW tutorial
  - How to run Phaser-MR from the Phenix GUI
  - Explaining the atom selection syntax
  - How to run Structure Comparison


- NEW FEATURES

  - Phenix.auto_sharpen:  New map sharpening algorithms available based on
      (1) optimization of map clarity (adjusted surface area or kurtosis),
      (2) resolution-dependent map-model correlation,
      (3) resolution-dependent half-map correlation, and
      (4) local map sharpening

  - Phenix.map_to_model: Iterative map improvement using model from one
      cycle to optimize sharpening of the map for the next cycle

  - Phenix.map_to_object:  Now you can individually map each member of a
     group of chains onto your target

  - Phenix.chain_comparison: Now can compare a target to each of a group of
     models

  - Phenix.find_alt_orig_sym_mate: Now uses ssmalignment.stats to see what
    residue pairs were included by ssm in a hypothetical calculation of
    the ssm-rmsd. This is to avoid including atoms that even for an optimal
    structural alignmentwould have been too far apart that ssm would use them.
    Excluding these is necessary to avoid having the calculated rmsd blow up.

  - Now using angles for restraining hydrogen bonds in alpha helices and
    beta sheets in secondary structure restraints.

  - Phenix.polder: New option (compute_box=True) to compute a difference map
    excluding the bulk solvent mask within a user-defined box. The map can be
    used to inspect interesting difference density peaks.

  - phenix.map_model_cc (phenix.model_map_cc): various enhancements and bug
    fixes.

  - phenix.real_space_refine:
    - improved ADP refinement (still group)
    - improved NCS constraints definition
    - multi-model PDB file is not output by default anymore
    - (experimental) option to use masked histogram equalized map as refinement
      target. Using such map dramatically improves refined model geometry. This
      map is not suitable for ADP refinement.
    - show CCmask, CCvolume and CCpeaks that report on local and gloabal model
      to map fit.


- RESTRAINTS

  - Specialised restraints for ligands with differing chemistry in the same
    model allows using the same 3-letter code but different restraints for
    a specific instance - apply_cif_restraints

  - Detailed control of peptide planes and omega angles
    - apply_all_trans
    - apply_peptide_plane
    - apply_cis_trans_specification: use an atom selection (CA) to force a
      cis or trans omega angle

- REEL
  - Improvements to robustness and features

- GUI
  - New interfaces
    - phenix.mtriage
    - phenix.map_box
    - phenix.auto_sharpen
    - phenix.map_to_model
    - phenix.ncs_average (under "Other Tools")
    - phenix.sculpt_ensemble (alpha)
  - Updated interfaces
    - Comprehensive Validation shows real-space correlations from
      phenix.map_model_cc when map files are used
    - Structure Comparison

- Updated dependencies
  - Python 2.7.8 -> 2.7.13
  - biopython 1.66 -> 1.68
  - matplotlib 1.5.1 -> 2.0.0
  - certifi 2017.4.17

- Validation
  - improved speed of clashscore validation
  - bug fix to handling of cis-peptides in ramalyze
  - contours for ramalyze, rotalyze, and cablam now accessible at
      https://github.com/rlabduke/reference_data
  - validation GUI now shows altloc identifirs for most outliers
  - improved internal code accessibility of validation statistics

- Phaser-2.8.0
  - refinement of cell for data (electron micrcoscopy)
  - more analysis using eLLG for component addition
  - more analysis using eLLG for polyalanine residues for target eLLG given rms
  - more analysis using eLLG for single atoms more ordered than average
  - more analysis using eLLG for chains/gyre-and-gimble
  * Phassade: Fast Fourier Transform Anomalous Substructure Determination
  - SOLU TRIAL records RF-LLG for trial orientation
  - gyre and translation function report RF recorded for trial orientations
  - 'scons --ccp4' command line flag for ccp4 build (generates ccp4-style
    headers)
  * improved tNCS detection and analysis
  * gyre and gimble improved for output and python scripting
  * NMA outputs all perturbed files in MODEL delimited pdb file
    '<fileroot>.nma.pdb'
  * sceds outputs domains to chain delimited file '<fileroot>.nma.pdb'
  - sceds domain pdb file have TER separated domains
  - NMA option to output separate pdb files for each perturbation (XYZOUT NMA
    ALL OFF)
  - sceds option to output separate pdb files for each domain (XYZOUT NMA ALL
    OFF)
  - NMA has option to turn off output of input coordinates (NMA ORIG ON)
  - pdb chain identifies can be left unchanged from input (XYZOUT CHAIN COPY)
    (for Eleanor)
  - bugfixes
    - allow peaks and purge percent to override one another if not default
    - rms estimation developer options not being set correctly
    - throw out of memory error if Patterson map allocation fails

1.11.1 (October 2016)
=====================
- Minor bug fixes

1.11 (October 2016)
===================

- General
  - Improved geometry restraints:
    - Conformation-Dependent Library for omega added (omega_cdl=True)
  - Installation
    - Rosetta installation centralised for phenix.rosetta_refine,
      phenix.mr_rosetta and ERRASER
  - Improved NCS search procedure with simplified parameters. Provides
    status of user-supplied NCS groups during validation (refused/modified/ok)
  - Updated dependencies
    - biopython 1.64 -> 1.66
    - sphinx 1.2.2 -> 1.4.4
    - ipython 2.1.0 -> 3.2.1
    - pip 6.0.7 -> 8.0.2
  - Neutron scattering tables: support ions

- Amber refinement
  - Alpha release dev-2203

- GUI
  - New interfaces
    - phenix.map_comparison
    - phenix.polder
    - phenix.structure_search
    - phenix.real_space_refine
  - New selection editor
    - Unified interface for selecting atoms
      - Secondary structure annotations
      - NCS groups
      - TLS groups
      - Refinement strategy options
  - Unicode support
    - Non-ASCII characters are supported for most fields (e.g. file paths and
      job titles)
    - Using non-ASCII characters in projects/jobs will prevent earlier versions
      that do not have Unicode support from opening the project list correctly
  - Migrated validation after phenix.refine to use regular MolProbity for
    consistency (older versions of Phenix will not open new jobs)
  - phenix.real_space_correlation can now accept map files
  - phenix.molprobity can now accept map files
  - Updated dependencies for Linux
    - libpng 1.2.52 -> 1.5.26
    - freetype 2.4.2 -> 2.6.3
    - gettext 0.18.2 -> 0.19.7
    - glib 2.12.11 -> 2.46.2
    - expat 1.95.8 -> 2.1.0
    - fontconfig 2.3.95 -> 2.11.1
    - render 0.8 -> 0.11.1
    - xrender 0.8.3 -> 0.9.7
    - xft 2.1.2 -> 2.3.2
    - pixman 0.19.2 -> 0.34.0
    - cairo 1.8.10 -> 1.14.4
    - pango 1.16.1 -> 1.38.1
    - atk 1.9.1 -> 2.18.0
    - libtiff 3.6.1 -> 4.0.6
    - gtk+ 2.10.11 -> 2.24.29
    - matplotlib 1.3.1 -> 1.5.1

- Maps
  - phenix.polder calculates omit maps which exclude the bulk solvent around
    the omitted region. This way, weak densities possibly obscured by bulk
    solvent may become visible.
  - phenix.model_map: Given PDB file calculate model map and output as CCP4
    formatted binary map file.
  - phenix.mask: Given PDB file calculate bulk-solvent mask and output as CCP4
    formatted binary map file.
  - phenix.auto_sharpen: Optimizes resolution dependence of a map to improve
    clarity
  - phenix.segment_and_split_map: Carries out segmentation of a map

- Model-building
  - phenix.map_to_model builds models in low-resolution maps
    - builds any combination of RNA/DNA/protein
    - if NCS present, builds the asymmetric unit of NCS and expands
        to the entire map
  - phenix.segment_and_split_map breaks up a map into the asymmetric
     unit of NCS and further subdivides that into contiguous regions
     of density
  - phenix.chain_comparison compares CA or P atoms of two models and
     identifies how many residues match and how many are in the same direction

- phenix.real_space_refine:
  - Support output of refined model in mmCIF format
  - ADP refinement runs by default at the last macro-cycle. Several CPUs can be
    used to speed up refinement: use nproc=NUMBER_OF_CPUs parameter for this.

- phenix.model_idealization: tool to idealize model geometry while staying as
  close as possible to the starting model. Currently proteins only. Idealize
  covalent geometry and secondary structure, as well as eliminate rotamer and
  Ramachandran plot outliers, C-beta deviations.

- phenix.geometry_minimization
  - ability to use reference torsion restraints
  - ability to use NCS constraints

- phenix.phaser
   - SOLUTION HISTORY tracks solution through positions in
     RF/TF/PAK/RNP peak lists
   - selection by CHAIN and MODEL for PDB coordinate entry
   - automatic search number for single search ensemble
   - packing 'trace' molecule can be entered independently of coordinates
     and map
   - read TNCS/anisotropy binary files to avoid refinement, when running
     through scripts
   - write tNCS and anisotropy parameters to binary files (non-python interface)
   - default reading of I (or failing that, F) from mtz file (LABIN optional)
   - trace for ensembles from maps = hexgrid of 1000+/-100 points
   - trace for ensembles from coordinates above 1000 Calpha = hexgrid of
     1000+/-100 points
   - trace for ensembles from coordinates twixt 1000 atoms and
     1000 Calpha = Calpha atoms
   - trace for ensembles from coordinates under 1000 atoms = all atoms
   - packing by pairwise percent only, other packing modes obsoleted
   - packing test during FTF run by default with 50% pairwise packing cutoff
   - automatic tNCS NMOL determination in presence of commensurate modulation
   - added MODE GIMBLE, which splits ensembles by chain for rigid body
     refinement
   - support for unicode
   - solution coordinates placed nearest to input coordinates if possible

- phenix.reduce
  - Updated reduce_wwPDB_het_dict as of Aug 12, 2016
  - New script for updating het dict
  - Reduce no longer rotates methionine methyls by default. -DOROTMET flag
    reinstates old behavior

- phenix.ramalyze
  - Improved handling of residue connectivity
  - Summary statistics provided for altloc A specifically where multiple
    altlocs present

1.10.1 (October 2015)
=====================

- General
  - Using buildbot for automatic build control
  - New tools for planning and evaluating SAD data collection and SAD datasets
  - All command-line installers converted to Python
  - Removed dependency on Apple tools for Mac installer relocation
  - Improved geometry restraints:
    - Nucleic acids: basepair planarity, hydrogen bonding and stacking
      restraints
    - Angle and torsion restraints for disulfide bonds
    - Conformation-Dependent Library is now the default (cdl=True)
    - Automatic linking has been expanded with carbohydrate linking being the
      default in all applicable programs
  - Names and places of parameters controlling secondary structure, reference
    and NCS restraints were changed
  - NCS constraints ("strict" NCS) for refinement (phenix.refine and
    phenix.real_space_refine)
    - Applied to coordinates and ADPs
    - Refine NCS operators
  - Partial LINK records support:
    - only in phenix.refine, phenix.geometry_minimization &
      phenix.real_space_refine
    - output LINK records generated using automatic linking and applied .edits
      (LINK in input PDB file are ignored)
  - HELIX/SHEET records in input file will be preserved in output PDB file

- phenix.refine
  - Optionally use AFITT for geometry restraints of ligands
  - mmCIF output now includes _refln.pdbx_r_free_flag tag to explicitly record
    original R-free flags (instead of implicitly in status tag)
  - Updated Phenix GUI to handle reorganization of NCS implementation
  - main.secondary_structure_restraints, main.ncs and
    main.reference_model_restraints parameters were removed
  - Improved rotamer fitting (faster, less outliers)
  - Bug fixes

  - phenix.real_space_refine:
    - rigid-body refinement with possibility to specify rigid groups
    - NCS constraints are used by default
    - NCS groups are found automatically, also can be specified manually
    - support for custom restraints
    - support for ligand CIF files
    - map can have any origin (non-zero), refined model is not shifted to
      zero origin
    - if CCP4 map is used resolution must be specified using "resolution"
      keyword
    - Enable reference model restraints, including restraining to self
    - May bug fixes, many tests added

- phenix.model_idealization (new command, command-line only):
  - Idealize geometry of the input model. Only secondary structure elements
    of proteins are currently supported.

- phenix.secondary_structure_restraints (proteins and DNA/RNA)
  - Three methods are available to annotate protein secondary structure
  - Protein secondary structure definition may be outputted as HELIX/SHEET
    records in PDB format
  - New algorithm to search for nucleic acid base pairs and stacked pairs.

- phenix.omegalyze:
  - Identify non-trans peptide bonds in proteins. Three categories of peptide
    bond geometry: trans, cis, and twisted (>30deg away from planar).

- phenix.cablam_validate:
  - Validation of protein backbones now uses separate contours for
    trans-proline vs cis-proline.

- Graphical interface:
  - Integrated help browser now supported on Linux
  - improved display of Xtriage results, with summary of problems with data
  - support for 'within' syntax in graphical atom selections
  - Model Reconstruction GUI: new GUI for building complete asymmetric units
    and biological molecules using MTRIX and BIOMT records, respectively,
    in the PDB file
  - New interface for phenix.anomalous_signal
  - New interface for phenix.scale_and_merge

- Phaser:
  - now allows direct use of intensities instead of amplitudes for molecular
    replacement
    - new treatment of effect of intensity measurement error introduced
    - improved handling of weak data in normalization and anisotropy correction
  - new occupancy refinement mode to retain or prune residues along
    chain
    - highlights loops and domains that are likely to need major correction
    - automatically invoked if high TFZ-score solution rejected for bad packing
  - subgroups feature: Phaser lists possible subgroups of space group, and
    if twinning is suspected the space group can be set explicitly to one of
    its subgroups for a molecular replacement search.  In this case, the data
    will be expanded to the unique set for the subgroup.
  - new "Search" panel GUI making it less confusing to define multicomponent
    searches
  - replacing "Composition" with "ASU contents" on the GUI
  - phased translation function can be used, even as part of MR_AUTO search
  - magnification factor (cell scale) can be refined for electron density as a
    model: useful when cryoEM image reconstructions serve as models
  - SAD likelihood now reported as log-likelihood-gain: positive number relative
    to null hypothesis
  - Phaser now prints out suggestion for number of copies related by successive
    application of NCS translations (TNCS NMOL parameter), if suggested value > 2
  - fixed bug affecting xyz-values > 100 in PDB for substructure
  - various other bug fixes

- phenix.rosetta_refine:
  - more options for map type for real-space refinement, including composit
    omit and feature-enhanced maps
  - partial support for ligands (requires separate parameter files)

- phenix.molprobity:
  - incorporated RNA validation

- phenix.xtriage:
  - refactored for greater modularity and consistency between command line and
    Phenix GUI
  - simplified and clarified feedback, including new summary panel
  - many bug and stability fixes from PDB-wide testing
  - detection of elliptical truncation
  - analyze completeness of anomalous vs non-anomalous dataset

- phenix.plan_sad_experiment (new command, GUI and command-line):
  - Estimate the I/sigI for a dataset that will be necessary to solve a
    structure by SAD

- phenix.scale_and_merge (new command, GUI and command-line):
  - scale one or more SAD datasets, optimally using unmerged original index
    data, to create a scaled, merged dataset and to estimate the half-dataset
    anomalous correlation

- phenix.get_cc_ano (new command,command-line only):
  - Calculate the anomalous correlation between two datasets

- phenix.get_patterson_skew (new command, command-line only):
  - Get the skew of an anomalous difference Patterson and the ratio of
    anomalous differences to anomalous error.

- phenix.merging_statistics:
  - calculate CC(anom) for split half-datasets

- phenix.multistart_sa:
  - utility for running multi-start simulated annealing refinements with
    varying random seed - sample uncertainty due to stochastic effects
    - also performs map averaging
  - parallelizes across multicore systems or managed clusters

- phenix.find_alt_orig_sym_mate:
  - Created new alias, phenix.famos, which is an anagram of the acronym for find_alt_orig_sym_mate
    and is easier to remember for the commandline.

- phenix.subgroup_symmetry_datasets:
  Calculates all subgroup of the diffraction data, and outputs them on the
  conventional setting. Useful to explore options when multple twin laws are
  possible.

- phenix.structural_domain_search:
  Domain search program using an approach based on graph theory. Very quick and
  fairly accurate for simple problems.

- phenix.structure_search
  Quickly find homologous structures for a given model.  There is an option to
  compile and output a list of ligands found in structures of those homologs.

- phenix.map_model_cc (and its synonym phenix.model_map_cc):
  Compute model-map correlation coefficient (map CC or RSCC). Output overall,
  per whole model, per chain and per residue CC. Input map needs to be CCP4
  fromatted file.

- phenix.fem:
  - various bug fixes
  - make using "maximal synthesis" the default (before "minimal synthesis" was the
    default)
  - current version obsoletes all previous versions

- phenix.tls_analysis:
  Decomposition of TLS matrices into elemental rigid-body motions (axes and
  amplitudes of vibration and libration). Validate TLS matrices.

- phenix.tls_as_xyz:
  Explicit interpretation of refined TLS parameters by decomposing into
  elemental rigid-body motions and generating structural ensembles that
  are consistent with these motions

- phenix.chains_as_models:
  Convert PDB file with multiple chains into PDB file with multiple models
  (MODEL-ENDMDL).

- phenix.models_as_chains:
  Convert multi-model (MODEL-ENDMDL) file into one model file with each model
  assigned unique chain ID.

- phenix.helix_sheet_recs_as_pdb_files:
  Split PDB file into several PDB files with each file representing HELIX or
  SHEET record found in PDB file header.

- phenix.map_value_at_point:
  Output map values at specified (x,y,z) points in space.

- phenix.map_box:
  Extract box around selected atoms with map. Output PDB with atoms in the box
  and map files in CCP4, X-plor and MTZ formats.

- phenix.model_model_distances:
  Given two PDB files output distances per atom, residue, chain, model and
  overall. It is assumed (and asserted in the code!) that the amount and
  order of atoms in both files are identical.

- phenix.simple_ncs_from_pdb:
  Tool to find NCS groups in PDB model. Available in previous versions.
  Major code re-write to imporove re-usabiliy.


1.9 (May 2014)
==============

- General:
  - Re-organized and updated documentation
  - binary installers for Ubuntu (10.04 or later) available
  - Python built as shared library in Linux builds (facilitates Rosetta hybrid)
  - new commands: phenix.composite_omit_map, phenix.rosetta_refine,
    phenix.molprobity, phenix.sceds, phenix.fem, phenix.real_space_refine,
    and many smaller utilities

- GUI:
  - Rearranged main interface to accompany documentation changes
  - New Phenix releases may be downloaded directly from the GUI (requires
    previously registered academic email address or permanent user credentials)
  - 64-bit Mac installers updated to wxPython 3.0.0.0, which restores
    drag-and-drop functionality
  - added command to remove large data files (.geo, .kin, .ccp4, etc.)
    from a project directory to conserve disk space
  - bug fixes for Mac OS 10.9 ("Mavericks")

- phenix.refine:
  - similarity restraints for group isotropic B-factor refinement
    - set by group_b_iso.use_restraints (default=True)
    - data/restraints targets weight is determined automatically; can be
      adjusted by setting group_b_iso.restraints_weight (larger value makes
      restraints stronger)
  - various rotamer handling bug fixes
  - Automatic linking generation in phenix.refine and other programs that
    require restraints can be requested using the parameters link_all=True
    - handles any appropriate ligand-protein or ligand-nucleic acid covalent
      bonds, including sugars, amino acid modifications, and other prosthetic
      groups
    - behavior can be adjusted by setting maximum bond distances parameters
  - added "bfactor" and "occupancy" keywords to atom selection syntax (also
    applies to phenix.pdbtools)
  - Integration with DivCon suite (http://www.quantumbioinc.com) for QM
    refinement of ligands (requires separate license and installation)

- HySS (phenix.hyss):
  - Substructure search much more fully automated:
    - Search first with Patterson seeds and direct method
    - If no clear solution, try more seeds, Phaser completion, and varying
      resolution
    - Full parallelization makes the search rapid

- AutoSol:
  - Now much more powerful for weak SAD data:
    - Incorporates new automated Hyss search
    - Classic solvent flipping density modification after Resolve density
      modification
    - Tests several values of smoothing radius for solvent boundary
      identifcation for weak data
    - Iterates heavy-atom substructure search with Phaser completion including
      density after density modification and using model after model-building
    - Require solutions scored by correlation to be substantially better than
      altermatives before eliminating alternatives
  - Building with helices_strands_only set as default for data at lower than
    3 Angstroms

- phenix.rosetta_refine (new command):
  - Hybrid Rosetta/Phenix refinement of protein X-ray crystal structures
    for difficult cases at low resolution
  - reference: DiMaio et al. (2013) Nature Methods

- phenix.molprobity (new command):
  - Unified validation tool for command-line users; similar to existing
    comprehensive validation in the Phenix GUI
  - optional output of script for viewing results in Coot, or kinemage file

- phenix.composite_omit_map (new command):
  - Default method uses very fast iterative de-biasing procedure
  - Can also generate either a simple or composite omit map with refinement
    and optional simulated annealing
  - Will generate multiple map types at once
  - parallelizable on multi-core systems or managed clusters

- phenix.sceds (new command, formerly phaser.splitter):
  - implements normal mode-based domain finding method
  - reference: McCoy et al. (2013) Acta Cryst D69

- phenix.ensemble_refinement:
  - enabled use of experimental phases (MLHL target)

- phenix.find_alt_orig_sym_mate:
  - offset the produced MLAD value with log(0.9) to avoid
    negative values if identical structures were to be compared.
    Consequently older versions therefore yield MLAD values which are
    0.10536 smaller than the current version.

- AutoSol, AutoBuild, and LigandFit:
  - You can now specify the value of the test set flag for the free data with
    test_flag_value=xxxx

- phenix.morph_model:
  - fixed serious bug present in versions 1449-1626 that resulted in incorrect
    offset being applied

- phenix.ligand_identification:
  - Now you can use ligand cif files (in addition to .pdb files) to build your
    own custom ligand library.  Just specify "Ligand directory" (GUI) or use
    keyword "ligand_dir" (command line) to tell the program to use .pdb or .cif
    files in your own ligand directory.
  - Bug fix to properly handle mtz labels when using difference maps.

- phenix.real_space_refine (new command):
  - tool to refine a model against a map (especially EM maps)

- phenix.fem (new command):
  - Compute "feature-enhanced" map ("FEM") - this protocol modifies a 2mFo-DFc
    sigmaA-weighted map such that the resulting map has strengthened weak
    signal, if present, and a reduced amount of noise and model bias.
    - combines randomization and averaging, treatment of missing reflections,
      residual composite omit map, and histogram equalization
    - resulting maps are more interpretable and have less bias compared to the
      starting 2mFo-DFc map
    - very fast procedure, no refinement required

- phenix.map_to_structure_factors (new command):
  - convert map into structure factors (also intended for EM)

- phenix.angle (new command):
  - compute angle between two axes or two atom selections

- phenix.model_model_distances (new command):
  - Given two PDB files compute distances per atom, residue, chain, model and
    overall. It is assumed (and asserted in the code!) that the amount and
    order of atoms in both files are identical.

- phenix.fab_elbow_angle (new command):
  - Calculating Fragment Antigen-Binding (Fab) elbow angle

- phenix.f000 (new command):
  - Given PDB file estimate F(0,0,0); includes atomic model and bulk-solvent

1.8.4 (September 2013)
======================

- Fixes several bugs in 1.8.3

1.8.3 (September 2013)
======================

- General:
  - 64-bit Windows 7 installer
  - Mac/Linux installers partially converted to Python
  - Improved bulk solvent mask calculations (phenix.refine, phenix.maps, etc.):
    - analyzes mFo-DFc map to identify spurious negative density blobs in bulk
      solvent regions and applies correction
      - for map calculation, more aggressive treatment used
    - eliminates "holes" in void regions in protein structures
    - usually decreases R-factors by 0.5-2.0%, sometimes more

- Graphical interface:
  - added pause/resume function for most programs on Mac/Linux
  - more consistent layout of program GUIs; most programs run in subdirectory
  - new and improved eLBOW GUI with "wizard" style setup
  - added utilities for conversion to/from mmCIF
  - new GUIs: ERRASER, phenix.dynamics

- phenix.refine:
  - bug fixes in occupancy refinement, hydrogen addition
  - new local and global real-space refinement (RSR) methods:
    - Global RSR (refine all atoms, as before) runs only if:
        d_min > 1.5, r_work > 0.25, and not neutron data
      or "individual_sites_real_space" is the only selected strategy
    - Local RSR: residue-by-residue, replaces old fix_rotamers option
      - fits both main chain and sidechain
      - only applies to residues with rotamer outlier, poor CC, or clashes
      - new residue fit serves as torsion restrains for subsequent macro-cycle
      - normally not performed when explicit hydrogens present
  - added omit_selection keyword: sets atom occupancies for a selection to zero
    prior to refinement
  - optional harmonic restraints on starting positions - used to keep structure
    sensible during simulated annealing omit map calculation
  - experimental feature: placement of elemental ions (place_ions=True)
    - looks for Mg, Cl, Ca, and Zn by default; specify elements=Cl,Ni,... to
      give explicit list (recommended)
    - works best with anomalous data

- Phaser:
  - various speed enhancements
  - new normal modes analysis features for domain finding (phaser.splitter)
  - bug fixes

- HySS:
  - converted to use standard Phenix parameter syntax
  - Many new optional features:
    - Parallel pre-scoring of 2-site solutions
    - Parallel completion and scoring
    - Comparison with a fixed solution
    - Resolution cutoff for Phaser rescoring
    - Top-scoring solution written even if no matches found
    - All solutions written
    - Skip direct methods completion
    - LLG completion sigma cutoff
    - Scoring of input sites
    - Randomization of input sites

- phenix.ligand_pipeline:
  - numerous bug fixes
  - more flexible handling of input sequence files
  - optional Sculptor step prior to MR to fix residue identities
  - sidechain pruning enabled by default
  - interactive mode now includes ligand selection

- AutoSol:
  - overall_best_refine_data.mtz now contains uncorrected F+,sigF+,F-,sigF-
    data and uncorrected FP,SIGFP.  Previously only the FP,SIGFP was
    uncorrected, and F+,sigF+,F -,sigF- had been corrected for anisotropy and
    sharpened.  Note that output map files from phenix.autosol and
    phenix.autobuild are still sharpened and corrected for anisotropy.
  - For MIR, if quick=True then as soon as a good solution is found with HySS
    all other derivatives are solved by difference Fouriers.  This can speed
    up MIR structure solution. You can disable this feature with
    skip_hyss_other_derivs_if_quick=False.

- AutoSol and AutoBuild:
  - It is now possible to specify a different number of copies of different
    chains in the sequence file.  If the sequence file has n copies of chain A
    and m of chain B, it is now assumed that the stoichiometry is
    ncs_copies * (nA, mB).  Previously it had been assumed always that the
    stoichiometry was ncs_copies *(A, B).
  - For structure determination with extremely little phase information there is
    a new keyword extreme_dm=True (default is False). This allows density
    modification to use parameters optimized for very poor phasing if the FOM
    of phasing is less than fom_for_extreme_dm.

- AutoBuild
  - phenix.autobuild now requires that if a hires (or refinement) file is
    supplied, it must have a complete set of free R flags (test set) if the
    datafile also has a test set.  In this case the test set from the
    datafile is ignored and the one from the hires/refinement file is used.
    - You can use the Phenix GUI to extend the test set from your datafile to
      match the data in your hires/refinement file if necessary.

- phenix.parallel_autobuild (new command):
  - a tool for running phenix.autobuild many times in parallel, picking
    the best results, and iterating with the new maps as a starting point.
  - can make effective use of as many as 40 processors.

- phenix.guided_ligand_replacement:
  - Guided Ligand Replacement (GLR) can now search for a protein-ligand complex
    suitable for the guided fitting using web services on the PDB website or in
    a local directory.
  - GLR ligand input has been reworked to allow the use of the three-letter
    code relying on eLBOW to generate the restraints.

- REEL:
  - simplification and improvement of molecule building
  - support for cis/trans E/Z enantiomers
  - view .geo files

- MR-Rosetta:
  - You can now tell mr_rosetta to ignore long-running sub-processes and just
    go on with the finished jobs.

- PDBTools:
  - simple dynamics option moved to new command, phenix.dynamics
  - remove_alt_confs and truncate_to_poly_ala options made compatible with
    other model modifications

- phenix.dynamics (new command):
  - simple Cartesian dynamics, replaces equivalent function in PDBTools

Version 1.8.2 (February 2013)
=============================

- General:
  - phenix.refine, phenix.reduce, and MolProbity updated to use consistent,
    slightly revised, electron-cloud centers for hydrogen positions as default
    (for X-ray)
    - Hydrogen nuclear positions available in both systems (for neutron, etc);
      MolProbity vdW radii redetermined to match each choice.
    - See MolProbity website for details: http://molprobity.biochem.duke.edu
  - mmCIF format now supported in place of PDB format for input to many
    programs
  - new commands: phenix.ensemble_refinement, phenix.ligand_pipeline,
    phenix.sort_hetatms, phenix.cc_star, phenix.pdb_editor,
    phenix.maximum_entropy_map

- Graphical interface:
  - new GUIs: phenix.merging_statistics, phenix.pdb_editor
  - Table 1: now calculates merging statistics directly if unmerged intensities
    are specified

- phenix.refine:
  - better geometry when refining with explicit hydrogens
  - better handing of sidechain rotamers in torsion NCS
  - support for LLG map type - uses Phaser to compute residual anomalous LLG
    map
    - best when used with group_anomalous strategy, which will flatten the
      LLG map around existing anomalous scatterers with refined f'/f''.
    - alternatively, anomalous_residual map type (similar but less sensitive)
  - automatic linking in refinement and pdbtools for residues in the same chain
    can be turn on using automatic_linking.intra_chain=True; this includes
    carbohydrate linking to protein and between sugars, and covalently bound
    ligands.
  - added optional output of mmCIF format model and data files

- phenix.ensemble_refinement (new command):
  - New method for refining ensemble models
  - Combines MD simulation with X-ray restraints that simultaneously account
    for anisotopic and anharmonic atomic distributions.
  - For full description see: Burnley BT, Afonine PV, Adams PD, Gros P. 2012.
    Modelling dynamics in protein crystal structures by ensemble refinement.
    eLife 1:e00311. doi: 10.7554/elife.00311.

- phenix.ligand_pipeline (new command):
  - combines Xtriage, Phaser, eLBOW, phenix.refine, AutoBuild, and LigandFit
    to automatically solve protein/ligand complexes
  - optional integration with Coot (parameter interactive=True) allows
    semi-interactive operation where more rebuilding is required

- phenix.autobuild:
  - now uses torsion NCS parameterization when NCS is used in refinement
  - building starting from a very accurate but very small part of a model:
    You can now use the keyword rebuild_from_fragments=True to start
    rebuilding from fragments of a model. You might want to use this
    if you look for ideal helices using Phaser, then rebuild the
    resulting partial model, as in the Arcimboldo procedure.  The
    special feature of finding helices is that they can be very accurately
    placed in some cases. This really helps the subsequent rebuilding.
    If you have enough computer time, then run it several or even many
    times with different values of i_ran_seed.  Each time you'll get a
    slightly different result.
  - Base-pairing in RNA building -- phenix.autobuild will now try to
    guess which bases in a model are base-paired, and if there is no
    positive sequence match to the model, the bases that are base-paired
    will be chosen to be complementary. You can set the cutoff for base
    pairing with the keyword dist_cut_base.

- phenix.mr_rosetta
  - You can now give commands for Rosetta in mr_rosetta, including
    a command to specify where disulfide bonds are located
  - The Rosetta models are now identified by an ID number so that they have
    unique names
  - Default number of homology models to download is now 1 (was 5)
  - Default number of NCS copies (if ncs_copies=Auto) is number leading to
    solvent content closest to 50%; if ncs_copies=None then all plausible
    values of ncs_copies are tested.

- phenix.autobuild and phenix.autosol building:
  - waters are automatically named with the chain of the closest macromolecule
    if you set sort_hetatms=True. This is for the final model only.
  - you can supply a target position for your model
    with map_to_object=my_target.pdb.  Then at the very end of the run  your
    molecule will be placed as close to this as possible. The center of mass
    of the autobuild model will be superimposed on the center of mass of
    my_target.pdb using space group symmetry, taking any match closer than
    15 A within 3 unit cells of the original position. The new file will be
    overall_best_mapped.pdb.

- phenix.ready_set:
  - methyl rotations now automatic when adding explicit hydrogens (also in
    phenix.reduce)
  - more options adding deuteriums to a model have been introduced

- phenix.morph_model and phenix.autobuild morphing:
  - You can now specify that only main-chain and c-beta atoms are to be used
    in calculating the shifts for morphing. The keyword to do this is
    "morph_main=True".

- phenix.assign_sequence:
  - Now can make use of the order of segments from molecular replacement
    models in sequence assignment

- phenix.fit_loops and phenix.density_outside_model:
  - run in memory instead of writing many files

- phenix.xtriage:
  - include merging statistics if unmerged intensities are used as input

- phenix.cc_star (new command):
  - combined assessment of model and data quality, as outlined in Karplus &
    Diederichs (2012) Science 336:1030-3.
  - summary also output by phenix.model_vs_data if additional unmerged_data
    keyword is specified

- phenix.maximum_entropy_map (new command):
  - tool to compute maximum entropy map from map coefficients. The program
    reads input Fourier map coefficients and modifies corresponding synthesis
    using Maximum-Entropy Method (MEM). The MEM modified map is everywhere
    positive, smooth and is of higher resolution. The method uses is a
    modification of Gull & Daniell (1978) algorithm.

- phenix.sort_hetatms (new command):
  - rearranges heteroatoms (ligands, waters, etc.) in a model so they are
    paired with the nearest macromolecule chain, similar to the PDB's
    processing of models

- phenix.pdb_editor (new command):
  - graphical editor for PDB files, based on tree view of model hierarchy
  - perform common operations such as chain renaming, renumbering, manipulation
    of atomic properties, add/remove atoms

- phenix.maps:
  - added wavelength parameter and support for anomalous residual and Phaser
    LLG map coefficients

Version 1.8.1 (September 2012)
==============================

- General:
  - primarily a bug fix release
  - new commands: phenix.erraser, phenix.merging_statistics
  - added geometry restraints for over 400 non-standard amino acid entries

- Graphical interface:
  - new GUI: phenix.find_alt_orig_sym_mate

- phenix.refine:
  - overall B-cart now applied to individual atomic B-factors (restores
    pre-1.8 behavior)
  - improvements to torsion NCS rotamer correction
  - real-space refinement (including rotamer fitting) for twinned structures

- phenix.erraser (new command):
  - New method from Fang-Chieh Chou and Rhiju Das (Stanford) combining
    Rosetta nucleic acid rebuilding with phenix.refine (requires separate
    installation of Rosetta 3.4 or newer)

- phenix.real_space_correlation:
  - simplified command-line behavior

- phaser (and AutoMR):
  - updated to version 2.5.2
  - Scattering from HETATM records that represent modified amino acids
    (e.g. MSE/SeMet) included so that maps no longer show holes for these
    HETATM records in the model; occupancies will not be reset.
  - Packing test was incorrectly applied to reoriented molecules where
    there was PTNCS

- phenix.merging_statistics (new command):
  - calculate statistics for unmerged data in a variety of formats,
    including <I/sigmaI>, R-merge, R-meas, R-pim, CC1/2
  - References:
    Diederichs K & Karplus PA (1997) Nature Structural Biology 4:269-275
    Weiss MS (2001) J Appl Cryst 34:130-135.

- phenix.kinemage:
  - improvements to handling of alternate conformations.

- phenix.combine_models:
  - Combine two models, taking the best parts of each to create a new model

- phenix.autobuild
  - Unassigned chains are now given the segid of "UNK" and are normally given
    chain ID="U"

- phenix.autosol
  - Scale of derivative in SIR (for RIP phasing) can be set with keyword
    "derscale"

- eLBOW:
  - New option, --multiple-planes, to generate planes around double bonds
  - Saturated rings now have dihedral restraints to maintain chair conformation

- REEL:
  - View relationships between ligands using PubChem fingerprint and tanimoto
    score.

Version 1.8 (June 2012)
=======================

- General:
  - New bulk-solvent and overall anisotropic scaling procedure:
    - approximately 100x faster than previous method
    - novel improved bulk-solvent scale model (no k_sol and B_sol)
    - affects phenix.refine, phenix.maps, validation, and any other programs
      which call these
    - more details:
      http://www.phenix-online.org/presentations/PhenixSantaFe2012_PA.pdf
  - Experimental Windows installer (32-bit only)
  - Improved support for anomalous difference arrays (F SIGF DANO SIGDANO ISYM)
    - propagated to output MTZ files with original data types
  - Overall Molprobity score now calculated in phenix.model_vs_data and the
    validation/refinement GUIs
  - New commands: phenix.den_refine, phenix.find_peaks_holes

- Graphical interface:
  - new GUIs: phenix.multi_crystal_average, phenix.morph_model,
    phenix.table_one, phenix.find_peaks_holes
  - improved layout for AutoBuild and reflection file editor GUIs
  - editor for custom bond/angle/plane restraints in phenix.refine
  - project groups: manage defaults for similar projects
  - phenix.table_one: generate table of crystallographic statistics for
    publication, as text or RTF document
    - data processing log files can be parsed for merging/scaling statistics
    - capable of handling multiple structures
  - validation: now using Top8000 database for Ramachandran statistics
    - added Ile/Val as separate class, and split proline into cis- and trans-

- phenix.refine:
  - change in default NCS restraints behavior:
    - use torsion-angle parameterization
    - restraint on NCS-related B-factors turned off
  - automatic matching of chains to reference model; allows multiple chains
    to be restrained to a single reference chain
  - bulk solvent and scaling speedup
    - old procedure is optionally available; in phenix.refine:
      "bulk_bulk_solvent_and_scale.mode=slow"
    - currently new procedure is not used in twin refinement
  - support for multiple reference model files
  - new atom selection feature for adp.individual.isotropic and anisotropic:
    "optional" keyword to allow for empty atom selections
    - typical example: adp.individual.isotropic="optional element H"
    - this will work both in the presence and absence of hydrogen atoms.
  - new improved filling of missing Fobs: use simple density modification.
  - bug fixes:
    - using ionic radii for atoms with charge column set, or recognized as
      elemental ions
    - hydrogen atoms now included in residual map calculation
    - limit sum of occupancies of hydrogen atoms on alternative conformers
    - more stable individual and group B-factor refinement (especially at
      low resolution)
  - slightly faster weight optimization
  - DEN refinement moved to phenix.den_refine

- phenix.den_refine (new command):
  - independent program for running DEN refinement (references:
    Schroder et al. (2010) Nature 464:1218-22.
    Schroder et al. (2007) Structure 15:1630-41.)
  - defaults for annealing modified for sensible compatibility with DEN
  - optimization turned on by default
  - weighting scheme improved for optimal DEN network influence

- phaser (version 2.5):
  - TNCS support enabled by default for both MR and SAD

- AutoSol:
  - TNCS support for SAD phasing

- phenix.pdbtools:
  - segID modifications: set to chain ID, or clear field
  - ability to set atomic charge field (used for X-ray scattering factors and
    VDW radii in nonbonded restraints)

- phenix.find_peaks_holes (new command):
  - standalone map peak search program: highlight difference map features
  - includes analysis of anomalous maps, flag suspicious "water" molecules

- phenix.maps:
  - simplified execution on command line: "phenix.maps model.pdb data.mtz"
    - generates default map coefficients (2mFo-DFc, mFo-DFc, anomalous if
      available) and 2mFo-DFc CCP4-format map.

- phenix.fetch_pdb:
  - added --maps option: generates 2mFo-DFc and mFo-DFc map coefficients

- phenix.get_cc_mtz_pdb:
  - added atom_selection keyword: compute CC for a subset of atoms

- phenix.real_space_correlation:
  - simplified execution on command line:
      "phenix.real_space_correlation model.pdb data.mtz"
    - calculates CC of model to 2mFo-DFc map
  - added atom_selection keyword: compute CC for a single selection (instead of
    per-atom or per-residue statistics)

- phenix_regression.wizards.list [characters-to-search]
  - Lists all the wizard regression tests (or all that contain
    characters-to-search if given). Convenient for finding a test that will
    serve as a template for your procedure and for making sure that everything
    should work.

Version 1.7.3 (December 2011)
=============================

- General:
  - new commands: phenix.ncs_average, phenix.model_vs_sequence
  - support for CNS v1.3 reflection files
  - numerous bug fixes

- Graphical interface:
  - new GUI: MR-Rosetta
  - phenix.model_vs_sequence incorporated into validation output
  - viewer for diffraction images (also available as phenix.image_viewer)
  - editable heavy-atom sites in HySS and AutoSol

- phenix.refine:
  - reduced memory footprint for map-related calculations
  - improved handling of parameter files in GUI
  - RNA geometry target now includes pucker & base type-specific Chi 1
    parameters
  - TLS group definitions from PDB file header will be automatically used if
    no groups are provided by the user
  - bug fix in clashscore calculation for large structures with duplicate
    chain IDs plus SEGIDs
  - torsion NCS restraints use SEGID if present

- phaser (version 2.4):
  - support for pseudo-translational non-crystallographic symmetry (tNCS) in
    molecular replacement
    - detects presence of tNCS from inspection of native Patterson
    - characterizes tNCS parameters (translation, small rotational difference
      between copies)
    - accounts for effect of tNCS in molecular replacement calculations
    - note: support currently limited to one tNCS operator relating two
      molecules (or sets of molecules)
  - FAST search method turned on by default, replacing FULL search method
  - MR possible with a model containing a single atom
  - B-factor refinement turned on by default in molecular replacement
  - Bugfix for occasional and irreproducible segmentation faults for 32 bit
    binary running on 64 bit Linux
  - Packing function accepts cases where an ensemble with internal point group
    symmetry is on a special position, and pdb file output deletes atoms
    overlapped due to the symmetry of special position

- phenix.ready_set/phenix.reduce:
  - improved handling of alternate conformations

- AutoSol:
  - Fixed problems with clean_up=True deleting overall_best.pdb
  - Allow data file names of any length
  - Use NCS in scoring by default if ncs copies > 2
  - Correct anisotropy and sharpen by default

- AutoBuild:
  - automatic anisotropy correction and sharpening for all maps

- phenix.multi_crystal_average:
  - automatic anisotropy correction and sharpening
  - automatic filling of reflections to highest resolution of any dataset
  - masking heavy-atom sites from averaging
  - inputs grouped by crystal

- phenix.model_vs_sequence (new command):
  - Verifies sequences of chains in PDB file against a sequence file
  - also part of the validation GUI

- phenix.cif_as_mtz
 - support for files containing multiple datasets
 - extraction of Hendrickson-Lattman coefficients and anomalous arrays
 - new map_to_asu and remove_systematic_absences options

- phaser.MRage (version 0.1.0):
  - new command-line switches for common actions, e.g. setting the verbosity,
    getting PHIL parameters, etc
  - define component by sequence, but calculate mw from available models if
    sequence is omitted
  - calculate number of copies to find from Matthews coefficient
  - new input scope: "model_collection": PDBs listed will be superposed with
    phenix.ensembler and optionally trimmed
  - new input scope: "search": ability to perform homology search (currently
    only BLAST is supported, either local installation or NCBI service)
  - alignments automatically generated for "template" scope if target sequence
    is known
  - various options for queueing systems, e.g. qslot
  - option to automatically write out solutions
  - fully parallel space group exploration

- phaser.MRage.solutions:
  - list solutions in XML and phaser script format

- phenix.sculptor, phenix.ensembler:
  - new command-line switches for common actions, e.g. setting the verbosity,
    getting PHIL parameters, etc

- phenix.maps:
  - added atom selection parameter for simple omit maps (omit.selection)

- phenix.ligand_identification:
  - Implemented methods to generate custom ligand library based on parent
    protein's SCOP or CATH terms, Pfam accession numbers, GO accession numbers,
    or InterPro ID
  - Improved ligand geometry after real-space refinement

- phenix.find_alt_orig_sym_mate (new command):
  - For different molecular replacement solutions from the same dataset, finds
    the copy of moving_pdb closest to reference_pdb with respect to all
    symmetry operations and alternative origin shifts permitted by the
    spacegroup in moving_pdb

- phenix.ncs_average (new command, mostly for development purposes):
  - Carries out NCS averaging of a map without any density modification
    (simple NCS averaging).
  - Useful for creating accurate maps for cases with high NCS symmetry and
    for evaluating what an NCS density modified map should be expected to
    look like.

Version 1.7.2 (September 2011)
==============================

- General:
  - new commands: phenix.cut_out_density, phenix.data_viewer,
    phenix.morph_model
  - default map output format changed to CCP4
  - Native support for CIF-format reflection files (e.g. from PDB) in many
    commands and GUI

- Graphical interface:
  - phenix.data_viewer added under Reflection Tools
  - new GUIs: HySS, phenix.cut_out_density, phenix.reciprocal_space_arrays
  - many bug fixes for 64-bit Mac build

- phenix.refine:
  - torsion-space NCS (ncs.type=torsion)
    - automatic parameterization, similar to reference model restraints and
      using top-out potential to allow genuine deviations between NCS groups
    - incorporates optional rotamer correction before refinement
  - new Xray/restraints weight optimization. Keywords optimize_wxc and
    optimize_wxu replaced with optimize_xyz_weight and optimize_adp_weight,
    respectively.
  - grid searches parallelized: use 'nproc' parameter to set number of CPUs
    to use
    - this affects bulk solvent mask optimization, and XYZ and ADP weight
      optimization (in addition to TLS identification)
    - can reduce runtime by 75% or more depending on settings (see manual for
      details)
  - option to use Conformation-Dependent Library restraints (Berkholz et al.
    [2009] Structure 17:1316)
  - improved handling of hydrogen atoms in refinement: results in lower
    clashscores at low resolution
  - increased nonbonded restraints weight when hydrogens are absent: improves
    clashscore and Ramachandran statistics
  - rigid-body refinement mode defaults to grouping by chain
  - printout of validation statistics during coordinate refinement
  - bug fix in phase-combined map calculation results in significantly improved
    maps
  - bug fix in X-ray/ADP restraints weight calculation results in smaller
    B-factors deviations for bonded atoms

- phenix.cut_out_density (new command):
  - carve out a region from map coefficients (MTZ file)
  - density may be specified as a sphere, a box, or masked around a PDB file

- phenix.data_viewer (new command):
  - visualization of reciprocal-space reflection data (similar to 'hklview' in
    CCP4i)
  - 3D OpenGL view of all data, or 2D view of planes (pseudo-precession
    photograph)

- phenix.morph_model (new command):
  - Use after molecular replacement to improve model for rebuilding
  - Morphs model to match a prime-and-switch map

- Phaser (and AutoMR wizard):
 - various bug fixes and minor enhancements
 - significant speedup of fast rotation function

- phenix.reel:
  - better integration with main GUI, especially eLBOW

- phenix.fobs_minus_fobs_map:
  - option to disable strict unit cell isomorphism check

- phenix.ramalyze:
  - optional output of PNG images of Ramachandran plots; add --plot to command
    line arguments

- phenix.fetch_pdb:
  - new flags: --all fetches PDB, CIF, and FASTA files simultaneously; --mtz
    runs phenix.cif_as_mtz on downloaded CIF file

Version 1.7.1 (April 2011)
==========================

- General:
  - Graphical installers for Macintosh
  - Binary installer for 64-bit Mac Intel (OS 10.6 only)
  - new commands: phenix.french_wilson, phenix.reciprocal_space_arrays,
    phaser.brunett
  - many bug fixes

- Graphical interface:
  - New GUIs: phenix.phase_and_build, phenix.mtz2map, phenix.french_wilson
  - Structure comparison GUI:
    - Support for homologous sequences: uses MUSCLE to map residue numbering
      to reference structure
    - Insertion codes skipped when extracting sequence
    - PyMOL plugin: automatically visualize regions of difference
    - improved parallelization
  - Plugin for running multiple sequence alignment using MUSCLE (Edgar B.
    [2004] NAR 32:1792)

- phenix.refine:
  - output MTZ file with original data, amplitudes used in refinement, and
    map coefficients (replaces *_map_coeffs.mtz)
  - automatically run French&Wilson treatment if using intensities
  - secondary structure restraints not included in RMS(bonds) calculation
  - optional automatic rigid-body grouping by polymer chain
    (rigid_groups_mode=one_group_per_chain)
  - "top-out" potential with smooth cutoff for reference model restraints
  - automatic determination of Saenger class for nucleic acid base pairs (if
    not provided by user)
  - ability to run phenix.find_tls_groups internally to obtain atom selections
    - option: tls.find_automatically=True
  - automatic non-bonded distance correction for H-bonded atoms

- phenix.mr_rosetta:
  - Changed to allow "*" in sequence file
  - Fixed bug in which identity was set to 100% when running from hhr files
    because mr_model_preparation changes the sequence. Now identity is
    determined by comparison of original model sequence and sequence file.
  - Added web site with data used in Nature Rosetta+Phenix paper at
    http://www.phenix-online.org/phenix_data/terwilliger/rosetta_2011/
  - Added keyword include_solvation_energy=False to allow use of mr_rosetta
    for membrane proteins.

- phaser.brunett (new command):
  - automated molecular replacement using incremental exploration
  - parallel option using local CPUs or queuing systems (SGE, LSF and PBS are
    supported)
  - quick SSM-superposition-based LLG evaluation for alternative models for the
    same component
  - "puzzle"-assembly of solutions for a more complete solution
  - novel space group determination algorithm that retains possible space
    groups until they are significantly outperformed by others
  - handles custom models (used without modification), templates (a PDB
    template and an alignment, Sculptor is run to create models from it) and
    homology searches (PDB and alignment obtained from homology search,
    currently HHPRED and BlastXML; data downloaded from EBI)

- phenix.french_wilson (new command):
  - Implements French & Wilson correction of negative and weak intensities (see
    Acta Cryst A34:517 [1978])
  - outputs modified amplitudes in MTZ format
  - used internally by default in Phaser GUI and phenix.refine when intensities
    are provided

- eLBOW:
  - Geometries can be generated using CSD Mogul substructure searches of
    experimental geometries. Estimated standard deviations of the bonds
    and angles are also included in the restraints file.

- phenix.reciprocal_space_arrays (new command):
  - this tool takes a PDB model and a data file and outputs an MTZ file
    containing a number of arrays that are typically used in various
    crystallographic calculations, such as R-factors, maps, refinement targets,     etc. The output arrays include:
    - Fmodel, Fcalc, Fmask, Fbulk, FOM, resolution, Hendrickson-Lattman
      coefficients (input original, computed from model and combined),
      FB_CART (overall anisotropic scale), input Fobs and Rfree-flags and more.

- phaser/AutoMR:
  - bug fix in residue number input, now allows four-digit IDs

- phenix.multi_crystal_average:
  - Added resolve_size keyword

- phenix.get_cc_mtz_mtz:
  - Removed excessive printout for clarity

- phenix.autosol:
 - Added the text "+/-2SD" to printout for uncertainties for clarification

- phenix.find_ncs:
 - Fixed minor bug in which ncs_domain_pdb was not defined

- KiNG (phenix.king):
 - Improved opening multiple electron density maps: tiered windows, with
   automatic deciphering of file name to preset sigma levels and colors.
 - Added -phenix command line option to preset map colors to be more like Coot.

Version 1.7 (January 2011)
==========================

- General:
  - 64-bit Snow Leopard supported when installing from source
  - New commands: phenix.find_tls_groups, phenix.mr_rosetta

- Graphical interface:
  - Simplified AutoMR interface for single-component searches
  - Multi-criterion validation kinemage, opens in KiNG
  - Validation of RNA, isotropic ADPs, occupancies
  - Improved project management interface
    - when upgrading from a previous version of Phenix, will automatically
      import old projects to new format
  - Integrated help browser on Mac
    - contextual help buttons embedded in program windows
  - phenix.find_tls_groups available as part of phenix.refine GUI
  - improved and expanded documentation for wizard and Phaser GUIs

- phenix.refine:
  - Ramachandran plot phi/psi restraints for low-resolution refinement
    - depending on weight, can eliminate Ramachandran outliers
    - command-line switch: ramachandran_restraints=True
  - new option: set_zero_occupancies_if_in_poor_density
    - option to set occupancies to zero for residues where density is poor
      or absent (map CC is less than defined value).
    - this allows avoiding gaps in chains by including residues that miss
      density.  The occupancies are reset autmatically and updated every
      macro-cycle.
  - new option: use_statistical_model_for_missing_atoms
    - supposed to help refining very incomplete models and produce much
      improved maps. Similar to protocol described in Blanc et al. Acta
      Cryst D60:2210-21 (2004).
    - still very very developmental version but can be used if needed.
  - new real-space refinement protocol (mmtbx.lockit)
    - fast and advanced real-space refinement that can be done as part of
      refinement. Typically it doesn't do any improvement for already good
      models, but can reduce R-factors by up to 10% for initially poor models.
  - reference model methods now allow for sequence alignment between the
    working model and reference model, allowing for insertions, deletions,
    and mutations
  - bug fix in occupancy refinement: in some very rare cases the sum of
    occupancies was not equal to 1 exactly after group-constrained occupancy
    refinement.
  - improved N/Q/H-flips: residues that clash in both orientations are flagged,
    pointing the user to a potential building error that should be evaluated by
    hand.
  - improved secondary-structure/base-pair restraints (see below)
  - automatic support for D-peptides when using standard PDB residue names:
     DAL DAR DAS DCY DGL DGN DHI DIL DLE DLY
     DPN DPR DSG DSN DTH DTR DTY DVA MED

- phenix.mr_rosetta (new command):
  - Integration of Phenix molecular replacement and model-building with
    Rosetta model-building for difficult structures
  - requires Rosetta version 3.2 (expected January 2011)

- phenix.find_tls_groups (new command):
  - identifies suitable atom selections for TLS refinement
  - similar to TLSMD, but uses cross-validation to yield one unique solution

- phenix.model_vs_data:
  - updated database (to use in POLYGON, phenix.r_factor_statistics, etc..)
  - SIGMA-A plot vs resolution
  - separate ADP statistics for macromolecule, solvent and ligands

- phenix.guided_ligand_replacement:
  - Works with different crystal settings (Nat the space group and ASU
    and number of copies of the protein is now handled.)
  - Real space refinement used lockit and can include some of the
    surrounding protein

- eLBOW
  - When enumerating chiral isomers, the quantum chemistry optimisation
    selected is performed on each isomer.
  - Enumerating chiral isomers can be done in parallel (most helpful
    when doing QM)
  - Improved the dihedral specifications in planes including rotatable
    conjugated branches and peptide planes in rings.

- phenix.secondary_structure_restraints:
  - expanded RNA base pair restraints to include non-Watson-Crick pairings,
    including entire Saenger classification
    - Saenger class now used by default to specify bonding pattern
  - many bug fixes

- phaser:
  - option to specify fixed ensemble (partial known structure) in MR searches
  - new fast search method that finds multiple copies much more efficiently
  - various minor enhancements and bug fixes

- phenix.ensembler:
  - new sequence alignment mode: ssm (using code from Eugene Krissinel).
  - new structure clustering to help classifying similar structures (either
    for molecular replacement or for conformational analysis).
  - several small improvements, including faster superposition code and
    improvements in weight iteration.

- phenix.sculptor:
  - new mainchain deletion mode: instead of specifying sequence similarity
    threshold, specify target completeness and it works out threshold
    automatically (more widely applicable defaults).
  - possiblity to input target sequence index in an alignment file in case
    target sequence is not the top one.
  - input target sequence and let sculptor do an alignment with the protein
    chain.
  - faster accessible solution area calculation.
  - fast alignment search (for large multi-sequence alignment files).
  - internal optimisations in structure manipulation for faster results.

- mr_model_preparation:
  - new program to prepare files for MR. Run with homology search results from
    HHPred. For each hit (max number can be limited), download structure from
    the EBI, extract alignment from HHPred and runs sculptor on all of them.

- phenix.kinemage:
  - Multi-criterion kinemage generation now handles custom restraints
    (.cif), useful for validation of ligands and other non-standard
    geometries.  Available in the GUI and from the command line.

- phenix.rna_validate (new command) :
  - RNA model validation is now available on the command line. Validation
    includes pucker outliers, bond length and angle deviations, and suite
    outliers, with output written to stdout.

- phenix.autobuild:
  - bug fix: failure to apply space group correctly when space_group=xxx is
    specified on the command line and input file is an MTZ file

Version 1.6.4 (July 2010)
=========================

- Bug fix release
- Multi-criterion plot in validation

Version 1.6.3 (July 2010)
=========================

- General:
  - SOLVE converted to C++; FORTRAN compiler no longer required
  - CCP4 map output supported for all applications that write XPLOR maps
  - New commands: phenix.superpose_maps, phenix.apply_ncs,
    phenix.phase_and_build, phenix.build_one_model, phenix.assign_sequence

- Graphical interface:
  - Supports execution of jobs on Sun Grid Engine queueing system
  - Maps for PyMOL now output in (smaller, faster) CCP4 format
  - new GUIs for phenix.find_ncs, phenix.superpose_maps
  - support for execution on Sun Grid Engine queueing system
  - display of sequence and secondary structure in wizards and validation

- phenix.refine:
  - improved map parameters (identical to phenix.maps)
  - new custom_nonbonded_symmetry_exclusions atom selection parameter
  - reference model refinement now supports selections for both the
    reference and working models (refinement.reference_model.reference_group)
    - this method is most useful in cases such as the reference model having
      one chain and the working model having multiple copies of the same
      chain, all of which should be restrained to the single chain in the
      reference model.

- phenix.autosol:
  - much faster and better model-building after density modification

- phenix.find_ncs:
  - now supports identification of NCS operators from density map

- phenix.fit_loops:
  - two main algorithms now used: fit short gaps using a loop library derived
    from high-resolution structures in the PDB, or build loops directly by
    iterative extension with tripeptides
  - the loop library approach (specified with loop_lib=True) is very fast,
    and is currently applicable for short gaps of up to 3 residues
  - the iterative extension approach is slower, but can be used for longer
    gaps (typically up to 10-15 residues)

- phenix.superpose_maps:
  - transforms maps to follow a molecular superposition
  - output in CCP4 or XPLOR format

- phenix.phase_and_build:
  - carries out an iterative process of building a model as rapidly as
    possible and uses this model in density modification to improve the map
  - up to 10x faster as AutoBuild, but model quality is nearly as good

- phenix.build_one_model:
  - quickly build a single model from a map and sequence file, or extend an
    existing model.
  - example of use (for building a new model):
      phenix.build_one_model my_map.mtz my_sequence.dat

- phenix.assign_sequence:
  - carry out an improved sequence assignment of a model that you have already
    built
  - once the sequence has been assigned, this method will use the sequence
    and proximity to identify chains that should be connected, and it will
    connect those that have the appropriate relationships using the new loop
    libraries available in phenix.fit_loops

- phenix.apply_ncs:
  - apply NCS operators from a my_ncs.ncs_spec file (from phenix.find_ncs) to
    a single copy of your protein to create all the NCS-related copies
  - example of use:
      phenix.apply_ncs my_ncs.ncs_spec my_model_one_ncs_copy.pdb

Version 1.6.2 (June 2010)
=========================

- Disk image-based installation for Macs

- Graphical interface:
  - new GUIs: eLBOW, phenix.cif_as_mtz
  - structure comparison tool: evaluate multiple structures of a single
    protein, and identify differences in geometry; superposes models and maps
    for viewing in Coot and PyMOL
  - rearranged phenix.refine, validation, and wizard GUIs for wider screens
  - many bug fixes

- phenix.refine:
  - automatic flipping of incorrect N/Q/H rotamers
  - custom planarity restraints
  - hydrogen-bond restraints for Watson-Crick base pairs
  - reflections with Fobs=0 used in refinement
  - automatic optimization of mask parameters
  - reference model restraints: used to steer refinement of working model
    - advantageous in cases where the working data set is low resolution, but
      there is a know related structure solved at higher resolution.
    - the higher resolution reference model is used to generate a set of
      dihedral restraints that are applied to each matching dihedral in the
      working model.
    - to use (also available in GUI):
        main.reference_model_restraints=True
        reference_model.file=my_reference.pdb
  - Alternate ideal angles for torsion definitions: CIF definitions in geostd
    now have an 'alt_value_angle' tag for torsion angle definitions.
    'alt_value_angle' allows for non-periodic specification of target angle
    values, and can take a single value or a set of comma-separated values.
  - phenix.refine can now write a geometry restraints file *after* refinement
    for comparison with the standard .geo file that is written *before*
    refinement.  To turn on this option, use output.write_final_geo_file=True.
    After refinement, a <run name>_final.geo file will be written.

- New program - KiNG:
  - The Java kinemage viewer KiNG (Protein Science 2009;18:2403-2409) can
    now be launched via 'phenix.king'.  KiNG requires a Java version of
    1.5.1 or later to be installed on your system.

- phenix.hyss:
  - Phaser rescoring in HySS: Trial solutions after real/reciprocal space
    recycling are evaluated using the SAD likelihood target in phaser (only
    applicable to SAD). LLG scores are then used for identifying top solutions
    and terminating the search if a clear solution has been found.

- phenix.maps:
  - B-factor sharpening for improving low-resolution maps

- phenix.mtz2map:
  - new command for converting map coefficients to XPLOR or CCP4 format maps

- phenix.ligand_identification:
  - Incorporated real-space refinement in ligand fitting
    - keyword: refine_ligand=True (default)
  - Added intermolecular non-bonded interaction terms in ligand scoring and
    ranking
    - keyword: non_bonded=True (default)

- more robust build of Fortran programs (solve, resolve_pattern) with
  static ifort libraries

Version 1.6.1 (March 2010)
==========================

- New commands:
  - phenix.fmodel: structure factors from PDB file, as real or complex numbers
  - phenix.maps: generate any map coefficients or XPLOR maps

- Graphical interface:
  - numerous bug fixes
  - support for "detached" processes run outside GUI
  - new GUIs for phenix.superpose_pdbs, phenix.get_cc_mtz_mtz,
    phenix.get_cc_mtz_pdb, phenix.fmodel
  - Coot interface for running phenix.sculptor
  - general maps GUI now uses phenix.maps
  - validation GUI consolidated, works with or without reflections
  - alpha-test features available via preferences setting
  - 20x faster loading of PDB files in selection editor

- phenix.refine
  - much faster rotamer fixing
  - initial implementation of secondary structure restraints (hydrogen bonds)
    - automatic assignment of sec. str. using KSDSSP (from UCSF CGL)
  - proper treatment of charges on metal ions for scattering factors

Version 1.6 (January 2010)
==========================

- C++ version of RESOLVE

- Graphical interface:
  - new interfaces for Phaser (MR and SAD), AutoMR, ligand identification,
    ensembler and sculptor tools, NCS identification
  - reflection file editor:
    - more options for manipulating data arrays
    - assign R-free flags in thin shells
    - allow changing of symmetry by expanding or merging reflections
  - maps GUI:
    - allow custom map types
    - more options for creating XPLOR-format maps
  - improved startup time for main GUI and several apps
  - preferences setting for font sizes
  - expanded documentation
  - numerous bug fixes

- phenix.refine:
  - support for residues located on symmetry elements:
    - Atoms with occupancy == 1 are subject to automatic handling of the
      symmetry.
    - Atoms with occupancy != 1 are assumed to take the special position
      factor into account in the occupancy. Automatic handling of the
      symmetry is turned off.
    - If all atoms in a residue have occupancies != 1, all nonbonded symmetry
      interactions between the atoms in the residue are disabled.
  - Important bug fix affecting handling of branched oligosaccharides:
    - apply_cif_link can now be applied multiple times to the same residue
  - RNA sugar-pucker specific restraints, slightly revised restraints
    for all RNA residues, automatic chain links to non-standard RNA/DNA
  - Revised PRO dihedral angle restraints compatible with "endo" rotamer
  - Torsion-angle simulated-annealing bug fixes and small enhancements
  - fix_rotamers bug fixes

- Wizards:
  - optional real-space refinement of ligands in LigandFit

- New commands:
  - iotbx.reflection_file_editor: command-line version of graphical editor
  - New experimental standalone real-space refinement utility: mmtbx.lockit
    Reads a pdb file and map Fourier coefficients from an mtz file.
    (No documentation yet.)

Version 1.5 (September 2009)
============================

- General:
  - Improved installer and application launchers for Mac
    - Installation to /usr/local is no longer required
  - Nightly builds of installer available on web
  - Extensive additions to FAQs pages
  - new commands: phenix.fit_loops, phenix.fobs_minus_fobs_map,
    phenix.clashscore, phenix.grow_density, plus many graphical programs

- Graphical interface:
  - new main GUI (old GUI still used for AutoMR)
  - two-way communication with Coot, automatic updating of model/maps
  - phenix.reflection_file_editor:
    - Combine data from multiple reflection files, extend R-free sets
    - Any input format supported, output as MTZ
  - phenix.create_maps:
    - Generate standard likelihood-weighted maps, e.g 2mFo-DFc, mFo-DFc,
      anomalous, etc.

- LABELIT added to Phenix distribution (http://cci.lbl.gov/labelit)
  - autoindexing program designed to correct commonly experienced problems
  - detects pseudotranslation in diffraction patterns
  - can interface with MOSFLM for integration of data

- Phaser:
  - updated to latest version of Phaser (2.2)
  - reduced memory requirements
  - for anisotropic data, map coefficients sharpened to falloff of strongest
    direction
  - SAD phasing features
    - improvements to SAD log-likelihood-gradient substructure completion
    - automated reassignment of atom types, when more than one type of
      scatterer
    - improved sensitivity of LLG completion protocol
  - MR features
    - facility to compare MR solutions with prior ("template") solution, e.g.
      from another program
    - refine relative B-factors of different components of MR solution
    - recognize internal symmetry in MR model, use to identify equivalent
      solutions
    - updated atom names consistent with PDBv3
  - phenix.sculptor:
    - new tool to automate selection and editing of molecular replacement (MR)
      models
  - phenix.ensembler:
    - new multiple superposition tool to automate construction of ensembles
      for MR

- phenix.refine:
  - Torsion-angle dynamics for simulated annealing
  - Real-space refinement using "simple" and difference density map targets
  - Allow custom constrained occupancy groups in occupancy refinement
  - Make use of electron density correlation and electron density value at
    atomic center, to make water picking more reliable
  - Option to compute Average Kick Maps
  - Enhancements for neutron and joint XN refinement: initial implementation
    of map-based DOD placement and local fit, local real-space search to
    correctly orient rotationally ambiguous X-H/D positions
  - Automatic finding of best residue rotamers based on Rotamer library,
    combined with real-space refinement (in development)
  - Use better (faster, less memory expensive, with one bug fixed) code for
    mask calculation
  - Complete re-do of the code that determines bulk-solvent ksol&Bsol,
    anisotropic scale matrix and twin-fraction

- Wizards:
  - new GUIs for AutoSol, AutoBuild, and LigandFit
  - Simplified and general format for wizard parameters files
  - Any combination of MAD/SAD/MIR datasets in AutoSol supported through
    parameters files

- phenix.elbow:
  - restraints and geometry files can be output for all chiral combinations
  - faster generation of ligand topology (3-5x)
  - more robust graph matching
  - extensive testing of MOL, SDF input
  - save a 2D respresentation
  - load chemical components using key
  - send output directly to REEL

- phenix.reel:
  - R/S nomenclature for chiral centres displayed
  - switch between planar and tetrahedral geometry
  - ability to compare geometries with ideal restraints
  - load chemical components using key

- phenix.ready_set:
  - more robust for the edge cases
  - pass a ligand CIF file to use instead of eLBOW generated file
  - add exchangeable deuteriums (including HIS)
  - replace all hydrogens with deuteriums
  - allow better metal coordination via more options
  - add hydrogene to water
  - converts PDB LINKS to edits
  - can remove all waters

- phenix.ligand_identification:
  - command line version of old GUI ligand identification strategy
  - allow users to input biological information of macromolecule to generate
    custom-made libraries (e.g. function=kinase; or EC=1.1.1.1 etc.)
  - allow custom library as list of 3-letter codes
    (e.g. ligand_list="ATP CTP APN FMN")
  - supports multi-thread computing

- Validation:
  - phenix.clashscore: all-atom contact analysis using Probe
  - phenix.real_space_correlation:
    - general tool to compute correlation between two maps
    - can work with ensembles of structures
  - phenix.validate_model (GUI only):
    - Shows all geometry outliers, with links to Coot/PyMOL
    - Includes Ramachandran, rotamer, C-beta, and bad contact analyses
  - phenix.validate (GUI only):
    - Extends phenix.validate_model with results from phenix.model_vs_data
    - Includes real-space correlation plots and outliers
  - phenix.polygon (GUI only):
    - Graphically ompare model statistics to similar structures in the PDB
    - reference: Acta Cryst. D65, 297-300 (2009)

- phenix.model_vs_data:
  - automatically detects and reports twinning (using Xtriage)
  - reports Molprobity-style set of numbers: Ramachandran statistics,
    clash-scores, CB-outliers
  - various bug fixes and enhancements.

- phenix.tls: tool for manipulating TLS records in PDB file headers

- phenix.grow_density:
  - a tool for local density improvement, as originally described in Acta
    Cryst. (1997). D53, 540-543 (in development)

- phenix.xtriage:
  - GUI rewritten with embedded graphs, tables, and text summaries
  - twin law selection integrated with phenix.refine GUI

- phenix.fit_loops: fast building of missing loops with RESOLVE

- phenix.find_helices_strands:
  - Rapid chain tracing (20 residues/sec) for evaluating map quality

- phenix.multi_crystal_average:
  - Multi-crystal averaging

- phenix.fobs_minus_fobs_map:
  - Utility for creating isomorphous difference maps (also in GUI)

Version 1.4 (December 2008)
===========================

- General:
  - Command-line validation programs: phenix.ramalyze, phenix.rotalyze, and
    phenix.cbetadev
  - PyMOL now included as general linux and Mac binaries

- phenix.refine:
  - Missing F(obs) filled in with F(calc) in default maps
  - Phase-combined maps created when experimental phases used
  - Selection of constrained occupancy groups
  - Graphical interface: phenix.refine_gui

- phenix.xtriage:
  - Graphical interface: phenix.xtriage_gui

Version 1.3 (August 2008)
=========================

- General:
  - Release version 2008
  - Changes to installer to check for csh/tcsh
  - changes in svn_checkout command
  - incorporated the PDB chemical components dictionary to aid with the
    interpretation and processing of ligands
  - infrastructure for OpenMP compilation added

- Wizards:
  - AutoSol:
    - Scoring now uses Bayesian estimation of correlation coefficient
    - Building changed to helices and stands by default (faster)
    - Building will continue with more thorough algorithm if it
      fails to build with the fast algorithm
    - When testing both hands in building the best model from the first
      hand will not be overwritten with the current model from the
      second hand unless it is better
  - AutoBuild:
    - Parallelized for running on SMP architectures using threads

- phenix.refine
  - complete rewrite of PDB handling - new model hierarchy, faster structure
    reading, and direct structure writing.
  - constrained group occupancy refinement applied to alternate conformations
    by default
  - algorithms implemented for hydrogen and deuterium picking in neutron/X-ray
    maps
  - numerous changes to better handle joint refinement of X-ray/neutron data
  - map coefficients written for an anomalous difference Fourier map if
    anomalous data is used

- phenix.xtriage
  - analysis of systematic absences

- phenix.elbow
  - extended to more element types
  - bug fix to deal with close atoms recognized as bonds
  - Numeric SVD code used for larger structures (addresses bug in tntbx)

- phenix.reel
  - more model manipulation features added

- Phaser:
  - minor bug fixes

- phenix.reduce
  - minor bug fixes

Version 1.3b
============

- General:
  - support for 2 letter chain IDs.
  - support for Hybrid36 atom and residue numbers.
  - Mac OSX port (PowerPC and Intel)
  - GUI modifications to support the Mac
  - Conversion to latest wxPython (2.8.4.0) and associated libraries
  - Installer modifications to support Mac-OSX port
  - Removed C++ library test on linux platforms
  - New tools: phenix.xmanip, phenix.pdbtools, phenix.reel,
               phenix.find_helices_strands,
  - New document generation system
  - Full documentation

- Wizards:
  - Simple command-line versions of all Wizards
  - AutoBuild:
    - Automated nucleic acid model building
    - Simple-omit, sa-omit, and iterative-build-omit maps
    - Automatic treatment of ligands when rebuilding models
    - Simple commands for building multiple models compatible with data
    - Simple command for rebuilding to produce one very good model
    - Prime-and-switch maps
    - Iterative model-building of strands and helices at low resolution
  - AutoMR
    - now using Phaser-2.1.2:
    - Default packing criteria relaxed in AUTO_MR mode
    - Composition can be estimated from solvent content
    - Default composition corresponds to 50% solvent
    - If AUTO_MR search finds some, but not all, components, the partial
      solution is output
  - AutoSol
    - SAD phasing starting from MR partial model using Phaser
  - Ligand identification
    - Improved ligand scoring and ranking algorithm.
    - Implemented additional analysis of ranking results
    - Improved task display options
    - Task now accepts user-provided ligand library

- phenix.find_helices_strands:
  - rapid fitting of secondary structure to an electron density map

- phenix.probe
  - better support for alternate names.
  - numerous bug fixes.

- phenix.reduce
  - updated het dictionary.
  - numerous bug fixes.

- phenix.refine:
  - Updated to handle the new atom names from the remediated
    PDB-V3. Currently protein and RNA/DNA only. For ligands with new names,
    use elbow.builder to generate CIF files with matching restraints.
  - Writing of a PDB deposition header
  - Constraints for occupancy refinement added (used in refinement of atoms
    in alternative conformations);
  - Default behavior change: now phenix.refine always refines occupancies
    for atoms in alternative conformations or for atoms with partial
    occupancies;
  - Water picking procedure improved: bug fix to remove waters located too
    far from the molecule, H-bond criteria added;
  - Significant improvement in handling of hydrogens in refinement;
  - Bulk-solvent and scaling procedure made even more robust (B_cart
    optimization is done a few times in-between grid search steps in
    ksol and Bsol);
  - Documentationa reworked
  - Vdw_1_4_factor changed from 2/3 to 0.8 (to improve
    geometry and MolProbity clash scores in particular)
  - Obsolete remove_hydrogens option removed
  - Site symmetry enforced for ADP shaking (weight calculation) in
    anisotropic ADP refinement (no more eventual crashes due to corresponding
    assertion error)
  - Complete set of tools for refinement at ultra-high resolution (higher
    than 0.9A; modeling of bond electron density with interatomic scatterers,
    IAS)
  - New support for custom angle restraints
  - Determination of nonbonded distances distinguishes between refinement
    with and without hydrogens (based on remove_hydrogens=True/False);
  - Numerous enhancements for joint Xray+Neutron refinement (clearer handling
    of input data files, better targets weight calculation protocol)
  - MZ (multiple zone) rigid body refinement protocol is further optimized to
    decrease the runtime, accommodate a provision for multiple refinable
    rigid bodies.
  - A set of examples is added.
  - Bug fixes:
    - Proper handling of chain breaks, to avoid incorrect angle,
      dihedral, chirality and planarity restraints. This bug
      was likely to lead to bad restraints for structures with
      alternative conformations.
    - Avoid Wilson B crash for no data if low resolution omitted
    - Proper handling of negative residue numbers
    - Avoid "Segmentation fault" on Macs in twin refinement
    - Twin refinement (fix for wrong weight for X-ray/ADP term => refinement
      unstable)
    - Anisotropic ADP refinement for a model with hydrogens
      (riding refinement of H)
    - Fix inconsistency in group B refinement due to deep_copy method
      problem in twin_f_model.py

- phenix.pdbtools (new command-line tool for performing various manipulations
  on a model in a PDB file and general model statistics):
  - shaking of coordinates (random coordinate shifts)
  - rotation-translation shift of coordinates
  - shaking of occupancies
  - shaking of ADP
  - shifting of ADP (addition of a constant value)
  - scaling of ADP (multiplication by a constant value)
  - setting ADP to a given value
  - conversion to isotropic ADP
  - conversion to anisotropic ADP
  - removal of selected parts of a model
  - generates complete statistics for model geometry and Atomic Displacement
    Parameters.
  - tools for structure factors calculation from PDB model.
  - A documentation page added.

- phenix.cif_as_mtz (tool for conversion of reflection files in CIF format
  into MTZ format). Main features include:
  - distinguishes between Iobs and Fobs;
  - distinguishes between Xray and Neutron data sets;
  - analysis and picks up valid R-free flags;
  - simple sanity check of data.
  - A documentation page added.

- phenix.elbow
  - Perform simple eLBOW jobs in COOT
  - Automatically perform the appropriate calculations on all the unknown
    ligands in a PDB file and combine the CIF results into one file
  - Covalent bonding of ligands to macromolecules (phase I implementation)
  - Print the sequence of PDB file (elbow.print_sequence)
  - User control of addition and writing of Hydrogens
  - Output of TRIPOS mol2 files
  - User control over naming of output files
  - elbow.join_cif_files has re-ordered arguments (target as last argument)
  - elbow.link_edits will generate "edits" from PDB LINK records for input
    to phenix.refine
  - elbow.metal_coordination will generate phenix.refine "edits"
    for metal coordination spheres including angles
  - improved chiral centre handling
  - access to all MSDChem SMILES strings
  - converted to PDBv3.0

- phenix.reel
  - alpha release of restraints editor
  - read, edit, optimise and write CIF files
  - run eLBOW from a GUI window

- phenix.xtriage
  - now handles unmerged data

Version 1.26.1b
===============

- Conversion to latest wxPython (2.8.1.1) and associated libraries
- Withdrawn support for SGI-Irix and Alpha-Tru64 (GUI libraries are
  no longer supportable on these platforms)

Version 1.26b
=============

- General
  - Reduce (from Jane and Dave Richardson, Duke University) integrated into
    PHENIX, available as phenix.reduce

- Wizards:
  - improved model building in AutoSol
  - improved model building in AutoBuild
  - Automated detection and application of NCS
  - Automated loop building, crossovers between chains in different models
    of a structure, and side-chain optimization
  - SA-composite omit maps
  - Simple omit maps
  - Iterative-build composite omit maps
  - Omit around atoms in PDB file
  - Model rebuilding in place

- Phaser:
  - Updated to version 2.1
  - Molecular replacement changes:
    - Clash test is more forgiving
    - Partial solution is output if full solution fails
  - SAD phasing includes partial structure option to start from molecular
    replacement solution

- Textal:
  - Minor bug fixes

- Resolve ligand identification:
  - Improved ligand ranking algorithm

- phenix.xtriage:
  - If pseudo translation peaks implicate a higher symmetry and/or smaller
    unit cell, the appropriate unit cell and spacegroup are listed.
  - Optional detwining and treatment of data, including protocols for
    anisotropy correction, outliers detection and twinning or detwinning
    of data
  - xtriage code reorganized in such a way that it can be easily used in other
    scripts and the wizards

- phenix.refine:
  - Bug fixes (major and minor)
  - Testing and bug fixes in TLS refinement
  - Final f_model exported as proper MTZ or CNS file (was ad-hoc before)
  - Group_anomalous refinement (for anomalous atoms such as Se)
  - Refinment of hemihedrally twinned data using a lsq target function
  - Output of 'detwinned' and gradient maps for twinned data
  - Support for refinement with hydrogens (not fully tested)
  - Simultaneous use of NCS restraints and automatic water picking
  - Individual anisotropic ADP restraints (similar to the SIMU restraint in
    SHELX)
  - Flexible model parameterization (refinement can use rigid body, isoADP,
    anisoADP, etc simultaneously)
  - Simultaneous refinement of structures against X-ray and neutron data
  - Bug fix in processing of multiple apply_cif_link blocks applied to the
    same residue
  - Changes to allow user specified covalent bonds: apply_cif_modification,
    apply_cif_link
  - Fixes to binning to deal with sparsely distributed free-R set (e.g. thin
    shells)
  - Individual occupancy refinement
  - Group occupancy refinement
  - Automatic NCS detection
  - Reorganization of phenix.refine parameters

- eLBOW
  - Perform simple eLBOW jobs in COOT
  - Automatically perform the appropriate calculations on all the unknown
    ligands in a PDB file and combine the CIF results into one file
  - Covalent bonding of ligands to macromolecules (phase I implementation)
  - Print the sequence of PDB file
  - Added classes to monitor peptide linkages, nitrogen planarity and acid
    groups to maintain geometry during coordinate calculations
  - Greatly improved molecule matching procedures
  - User control of addition and writing of Hydrogens
  - User control of bond cutoff for auto bonding
  - User control of number optimisation steps and convergence tolerance
  - Output of TRIPOS mol2 files
  - User control over naming of output files
  - User control of total memory usage
  - More tests added
  - AM1 optimisation includes peptide linkage correction and more efficient
    search procedures
  - elbow.join_cif_files has re-ordered arguments (target as last argument)

- phenix.explore_metric_symmetry:
  - Renamed from iotbx.ehms
  - Generation of point group graphs given a unit cell and space group
    combination
  - Comparison of unit cells and their sublattices

- phenix.twin_map_utils:
  - Refinement of twin fractions and map parameters for all possible twin
    laws given a fixed model and a given dataset

- mmtbx.remove_outliers
  - Outlier removal using Wilson or model based statistics

- PHENIX GUI:
  - Minor bug fixes

- Installer:
  - New solve/resolve CVS installation
  - Updates to libtbx package installation
  - Improved developer tools for updating PHENIX components (cvs_checkout,
    svn_checkout)
  - Backwards compatibility checks in phenix_env to detect an installation
    from an earlier OS and use it if necessary

Version 1.25b
=============

- Minor bug fixes

Version 1.24b
=============

- Wizards:
  - LigandFit Wizard made more robust and tested against 9000
    ligands in PDB
  - IterativeBuild and ModelBuild Wizards combined into single
    AutoBuild Wizard
  - AutoBuild Wizard
    - Introduced "rebuild_in_place" for rebuilding MR models
    - Full-omit map and multiple-model generation available
  - AutoMR Wizard
    - Automatic rebuilding of MR solutions using AutoBuild
      rebuild_in_place
  - AutoSol Wizard decision-making improved using model completeness
    as a quality measure
  - All Wizards can be fully run with scripts
  - Sample scripts available for all Wizards with data in $PHENIX/examples

- Phaser:
  - Integrated with AutoMr Wizard
  - Sphericity restraint on overall anisotropy
  - Improved pruning of list of plausible solutions
  - Various bug fixes for unusual space groups and non-Linux architectures

- Textal:
  - created a single integrated command-line script: textal.build
    options: build backbone only (C-alpha chains)            [--capra_only]
             build complete model (backbone and side-chains) [default]
             build side-chains for user-supplied C-alpha's   [--input_model]
  - reads reflection files or XPLOR maps for building
    options for region of space to build:
      - runs FINDMOL to automatically create map centered over a contiguous,
        symmetry-unique molecular region using symmetry ops and clustering
      - build in ASU [--asu]
      - build in user-supplied XPLOR map [--input_map]
  - will do sequence alignment to improve amino acid identity if
    amino-acid sequence file is provided [--sequence]
  - user can provide coodinates of Selenium sites (if available) to potentially
    improve accuracy of sequence alignment by identifying methionine
    residues [--se_sites]
  - now runs simulated annealing (phenix.refine) automatically as final
    step to estimate R-factor (as an indicator of model quality)
  - new unified Textal task/strategy in Phenix GUI
    - graphical front-end to textal.build script; same options, all on
      one dialog box

- eLBOW
  - Documentation available in installation and on-line
    - computer resources page
  - Uses bond length defaults from HF/6-31G(d,p) quantum calculations
  - Uses GDIIS to improve geometry optimisation convergence and DIIS in
    quantum convergence
  - Restart mechanism begins from previous best geometry
  - Expert/Novice modes
  - Can use GAMESS for geometry optimisation
  - Provides estimates of the computer resources needed for a calculation
  - Provides suggestions upon completion of a run
  - eLBOW will extract a specified ligand from a PDB file containing both
    protein and ligand
  - eLBOW will dynamically update geometry visualised in PyMOL
  - eLBOW can construct a short protein sequence
  - eLBOW can use a file as a template to name the ligand atoms
  - Database of SMILES and PDB files included
  - Input can be piped to eLBOW

- iotbx.ehms:
  - Minor updates to functionality and output
  - Centring type is inputable rather than just a spacegroup

- mmtbx.xtriage:
  - Maximum likelihood based estimation of twin fraction
  - Maximum likelihood based estimation of twin fraction
    while taking into account possible NCS parallel to
    the twin axis.
  - PHIL Interface rationalisation:
    - Detwinned data can be written out
    - various parameters can be set by the user if desired
  - Help/Manual updated (mmtbx.xtriage --h)
  - various bug fixes

- phenix.refine:
  - NCS restraints (coordinates and ADPs)
  - TLS refinement
    - Recalculation of TLS parameters from overall ADPs values
    - Parameter modification to remove positive definite atomic aniso-Us
    - TLS parameters written in output PDB file
    - Options:
      - TLS alone
      - TLS plus ADPs
      - TLS plus ADPs and coordinates
  - Group B-factor refinement
  - Rigid body refinement
  - Simple structure factor calculations (with or without bulk solvent
    and scaling)
  - Output MTZ coefficients for use in COOT
  - Increased checking of input PDB files and reflection data
  - Updated documentation

- PHENIX GUI:
  - color indicators added to show when a strategy/wizard icon has been
    double-clicked
  - strategies and tasks rationalized
  - improved default map and coordinate viewing options when displaying
    results in PyMOL

- Installer:
  - updated to 0.99rc6 version of PyMol
  - updated to more recent versions of freetype, fontconfig and pango on linux.
    This fixes a memory allocation problem on Fedora Core 5.
  - added --alias_mtype option to installer to force installation of
    a different machine type (binaries only). May be necessary on some
    SuSE systems.
  - now using freeglut and latest 8.4.x versions of tcl/tk
  - solve/resolve source code now included in the distribution
  - examples directory restructured and more examples/tests added

Version 1.23a
=============

- New tool: iotbx.ehms
  Explore Higher Metric Symmetry. Allows one to easely compute/list
  all possible point and spacegroups allowed by a certain point or
  spacegroup. A graph is constructed and is plotted if graphiv is
  installed.

- mmtbx.xtriage updates:
  - Data strength analysis
  - Low resolution completeness reporting
  - Wilson plot sanity analyses
  - Outlier detection based on Wilson statistics
  - Ice ring detection
  - Anomalous difference analyses
  - More elaborate reporting on type of twin laws
  - If calculated data is supplied next to experimental data,
    R.vs.R statistics is reported (Lebedev et al, Acta Cryst.
   (2006). D62, 83-95)
  - Missing rotational symmetry analyses. If symmetry is too low,
    new spacegroup is suggested on basis of R-value.
  - various bug fixes

- New tool: mmtbx.fest
  - Simple delta F and Fa construction from 2-lambda MAD, SAD SIR or RIP
    data.
  - Various scalings protocol can be chosen. Most (if not all)
    critical parameters can be set/modified by the user.

- New tool: eLBOW

Version 1.22a
=============

- Updated version of SOLVE/RESOLVE (2.11)
- Generation of environment file for sh/bash

Version 1.21a
=============

- Updated wizard infrastructure:
  - Wizards can call other wizards
  - Bug fixes and changes to AutoSol and AutoMR
  - New loop-fitting algorithm implemented in RESOLVE
  - New Iterative, TextalBuild and ResolveBuild Wizards
- New automatic ligand geometry and restraint building (ELBOW):
  - semi-empirical method for determination of coordinates
  - generation of coordinates from SMILES
  - generation of monomer library (cif) files for refinement
- Updated phenix.refine:
  - improved performance (1.5 to 2 times faster)
  - automated water picking
  - local-neighborhood individual ADP restraints
  - improved bulk-solvent/anisotropic scaling calculations
  - new automated weight generation scheme
  - refinement against neutron diffraction data
- New data analysis tool: mmtbx.xtriage
  - automated detection of anisotropy, twinning, and
    non-crystallographic translational symmetry
  - maximum likelihood determination of anisotropic and absolute
    scaling parameters
  - automatic determination of possible twin laws
- Infrastructure:
  - latest versions of thrid-party GUI libraries (glib, gtk+,
    wxWidgets and wxPython)
  - removed BLT from installation
  - linux installation modified to account for new x86_64 architectures
    and new linux 2.6 kernel
  - SGI installer

Version 1.2a
============

- Updated wizards

Version 1.1a
============

- Release version based on 1.08a
- Minor bug fixes for TEXTAL, SOLVE/RESOLVE strategies
- Minor bug fixes for PHENIX GUI
- Additional documentation - HySS and phenix.refine

Version 1.08a
=============

- Structure refinement (phenix.refine) including:
  - LBFGS minimization
  - Cartesian dynamics simulated annealing
  - Maximum likelihood with amplitudes
  - Maximum likelihood with amplitudes and phases
  - Restrained individual B-factor refinement
  - CCP4 monomer library for geometry restraints generation
  - Automatic input PDB interpretation
  - Automatic input data interpretation
- Latest version of PHASER (1.3)
  - Bug fixes from Phaser-1.2
  - Improved automated MR
  - Final phasing and refinement to full resolution
  - Improved packing
    - RNA/DNA
    - close-packed oligomers
    - only most homologous model in ensemble used for analysis
  - Reduced memory requirements
- Latest version of TEXTAL
  - Builds models from structure factors
    - users no longer required to prepare a centered density map
    - will automatically adjust resolution as necessary, and can
      now build models for datasets at 3.0A resolution or higher
  - Automated detection, centering, and tracing of protein molecules
    using the FINDMOL algorithm
  - Simplex optimized model building for better CA positions and side-chain
    identities
  - Task re-organization into hierarchical groups:
    - New high-level tasks to simplify model building
    - utilities and primitive tasks (for experts)
  - Invocation of many Textal routines from command line
- Wizards
  - Updated AutoSol and AutoBuild wizards
  - New AutoMR wizard for automated molecular replacement
  - New AutoLig wizard for automated ligand fitting
- Latest solve/resolve version - 2.09
- Updated solve/resolve strategies
- Infrastructure
  - PHENIX installation co-exists with other cctbx-based packages
  - PHENIX co-exists with Linux GNOME environment
  - Now using Python 2.4.1
  - Now using version 2.0+ of GTK+
  - Now using latest versions of wxGTK (wxWidgets) and wxPython
  - Now using PyMOL 0.97

Version 1.07a
=============

- Added Wizard infrastructure to PHENIX
- Wizards for automatic structure solution and model building
- Improved use of UNIX domain sockets
- Fixes in PHASER code for SGI and Alpha platforms
- Updated PHASER strategies
- Updated TEXTAL strategies
- New TEXTAL stitch algorithm for gap closure
- Updated SOLVE/RESOLVE strategies
- New strategies added for running phase improvement with RESOLVE
- Updates to the PHENIX GUI infrastructure
- Splash-page documentation added
- SOLVE/RESOLVE incorporated in the PHENIX distribution
- Test for appropriate version of make for source installation
- Fixes to glib configure script for Alpha
- Fixes to wxGTK configure and source code for Alpha
- Installer compiles Python with C++ compiler if required
  for wxPython compilation
- RedHat 7.3 removed from supported platforms
- Additional command line options for HySS
- Improved the killing of strategies
- Added ability to run strategies in queuing systems
- Sun Grid Engine currently supported for queuing (this may work
  with other queuing systems, but is not tested)

Version 1.06a
=============

- Minor bug fixes in GUI
- Added option to use rsh instead of ssh to run remote jobs
- New detection of network status on startup
- Use UNIX domain sockets if network is not functioning
- Updated tasks and strategies for PHASER
- Improved map reading speed in TEXTAL tasks

Version 1.05a
=============

- New tasks and strategies for PHASER
- New tasks and strategies for TEXTAL
- New map database for TEXTAL
- Bitmaps for task and strategy canvas in PHENIX GUI

Version 1.04a
=============

- New tasks and strategies for SOLVE/RESOLVE
- Improvements to hostname resolution to deal with incorrectly
  configured systems

Version 1.03a
=============

- New installation scripts
- End-user installation from complete binary bundle
- Added SGI Irix support (remote job execution with rsh only)
- Added Redhat 9.0 support - to deal with NPL problems
- Generic Linux installation for platforms other than Redhat 7.3
- Bug fixes in GUI
- Improved support for remote job execution
- Moved to hostname only resolution, no use of explicit IPs
- Improvements to HySS performance

Version 1.02a
=============

- Fixed bug with making connections between tasks in GUI
- Fixed file ownership problems for binary bundles in installation
- Fixed problem with absolute links in binary bundle (all linked files removed)
- Changed install script back to using csh rather than tcsh

Version 1.01a
=============

- Numerous bug fixes in GUI
- Added text to better document PHASER and TEXTAL tasks
- Fixed file permission problems for TEXTAL PERL scripts
- Improvements to HySS procedure

Version 1.0a
============

- First alpha test release
