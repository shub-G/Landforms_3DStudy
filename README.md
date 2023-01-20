# Landforms_3DStudy
This repository contains the source code for the 3D modelling study of Landforms evolution along continental margins over one glacial cycle, that was part of the following article, "Characteristics and controls of erosional landforms formed by offshore groundwater flow in siliclastic continental margins", submitted to GRL.

## Instructions for installation
* System requirements: 
  * Ububtu LTS 18.04 or higher
  * packages: git, auto-conf, cmake (>=3.13), clang, g++, gcc, gfortran, superlu (dev), lapack, libblas, libboost, metis (dev), parmetis (dev), gmsh, paraview, openmpi 
* Execute following commands in terminal (in given order):
  * mkdir dune_2_8
  * cd mkdir dune_2_8
  * mkdir source
  * cd mkdir source
  * cp /downloads/installer.sh .
  * chmod 755 installer.sh
  * ./installer.sh dune
  * cd dune
  * ./buildmobules.sh
  * cd ../..
  * chmod 755 newproject.sh
  * ./newproject.sh LandslideLandscape
    * On prompt for "2) Which modules should this module depend on?", enter: dune-common dune-geometry dune-uggrid dune-grid dune-localfunctions dune-istl dune-typetree dune-functions dune-alugrid dune-pdelab
    * On prompt for "3) Project/Module version?", enter 1
    * On promt for "4) Maintainer's email address?", enter your email address.
  * chmod 755 buildproject.sh
  * ./buildproject.sh LandslideLandscape
  * cd LandslideLandscape/src/
  * rm -rf LandslideLandscape.cc
  * rm -rf CMakeLists.txt
  * cp \_all_source_files_in_repo\_ .
  * cd ../..
  * chmod 755 compile.sh
  * ./compile.sh LandslideLandscape

## To run the simulations:
* Execute following commands in terminal (in given order):
  * cd \_HOME\_/dune_2_8/LandslideLandscape/release-build/src
  * ./main \_input-file\_  
    * In this case: default_inputs

## Files included in this repo:
* installation files:
  * installation/installer.sh
  * installation/newproject.sh
  * installation/buildproject.sh
  * installation/compile.sh
* source files 
  * duneincludes.hh
  * main.cc
    * Implement the correct MODEL_PATH, INPUT_PATH, and OUTPUT_PATH before compiling.
  * driver.hh
  * lop_elasticity.hh
  * lop_flow.hh
  * postprocess.hh
* problem-specific files
  * parameters.hh
  * default_inputs.ini
  * sample_numerical_runs.tar.gz
  * numerical_runs_readme.txt
    * contains instructions for reading the numerical_runs
