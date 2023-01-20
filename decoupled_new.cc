/*
 * decoupled_new.cc
 *
 *  Created on: May 18, 2022
 *      Author: sgupta
 */


#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include<math.h>
#include<iostream>
#include<fstream>
#include<vector>
#include<map>
#include<string>
#include<stdlib.h>
#include<time.h>
#include<exception>
#include<chrono>

/***************************************************
 * SELECT DIMENSION AND PROBLEM
 ***************************************************/
#define DIMENSION 3
#define PROBLEM_CONTINENTALMARGIN02

/***************************************************
 * INCLUDE SOURCES AND OPERATORS
 ***************************************************/
#include"landslidelandscape/duneincludes.hh"
#ifdef PROBLEM_CONTINENTALMARGIN01
#include"landslidelandscape/decoupled_new/problem_continentalmargin01/parameters_1.hh"
#elif defined(PROBLEM_CONTINENTALMARGIN02)
#include"landslidelandscape/decoupled_new/problem_continentalmargin02/parameters.hh"
#elif defined(PROBLEM_CONTINENTALMARGIN03)
#include"landslidelandscape/decoupled_new/problem_continentalmargin03/parameters.hh"
#endif
#include"landslidelandscape/decoupled_new/operators/postprocess.hh"
#include"landslidelandscape/decoupled_new/operators/lop_elasticity.hh"
#include"landslidelandscape/decoupled_new/operators/lop_flow.hh"
#include"landslidelandscape/decoupled_new/driver.hh"


int main(int argc, char** argv)
{
	try{
		//----------------------------------------------------------------------------------
		//––––––––––––––––––––––––––––––––––––––––––
	    // Maybe initialize MPI
	    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
	    if(helper.rank()==0){
	    std::cout << "Hello World! This is MODEL:DECOUPLED_NEW of PROJECT:LANDSLIDELANDSCAPE." << std::endl;
	    }
	    if(Dune::MPIHelper::isFake){
	      std::cout<< "This is a sequential program." << std::endl;
	    }
	    else {
	    	if(helper.rank()==0){
	    		std::cout<<"I am rank "<<helper.rank()<<" of "<<helper.size()<<" processes!"<<std::endl;
	    	}
	    }

	    //––––––––––––––––––––––––––––––––––––––––––
		// INPUTS
	    if (argc!=2)
	    {
	    	if(helper.rank()==0){
	    		std::cout << "usage: ./decoupled_new <input-file> " << std::endl;
	    	}
	        return 1;
	    }

		//––––––––––––––––––––––––––––––––––––––––––
		// DUNE MODEL PATH
	    std::string MODEL_PATH = "/home/sgupta/dune_2_8/LandslideLandscape/src/landslidelandscape/decoupled_new/";
	    // PROBLEM PATH
	    std::string PROBLEM_PATH;
		#ifdef PROBLEM_CONTINENTALMARGIN01
	    PROBLEM_PATH = "problem_continentalmargin01/";
		#elif defined(PROBLEM_CONTINENTALMARGIN02)
		PROBLEM_PATH = "problem_continentalmargin02/";
		#elif defined(PROBLEM_CONTINENTALMARGIN03)
		PROBLEM_PATH = "problem_continentalmargin03/";
		#endif
	    // INPUT PATH
	    std::string INPUT_PATH = "/home/sgupta/dune_2_8/LandslideLandscape/src/inputs/decoupled_new/" + PROBLEM_PATH ;
	    // OUTPUT PATH
	    std::string OUTPUT_PATH = "/home/sgupta/dune_2_8/LandslideLandscape/src/outputs/decoupled_new/" + PROBLEM_PATH ;
	    // INI-FILE FOR USER-DEFINED INPUTS
	    char INI_FILE[100];
	    sscanf(argv[1],"%99s", INI_FILE);
	    std::string input_file = INPUT_PATH;
	    input_file += INI_FILE;
	    input_file += ".ini";
        if(helper.rank()==0){
	    std::cout<< "input file: " << input_file << std::endl ;
        }

        //––––––––––––––––––––––––––––––––––––––––––
        // PARAMETER TREE
	    Dune::ParameterTree ptree;
	    Dune::ParameterTreeParser ptreeparser;
	    ptreeparser.readINITree(input_file,ptree);
	    ptreeparser.readOptions(argc,argv,ptree);

	    //––––––––––––––––––––––––––––––––––––––––––
		// MESH
		/* READ MESH
		 * AND GENERATE GRID VIEW
		 * WITH UG MESH
		 */
		using GridType = Dune::UGGrid<DIMENSION>;
		GridType grid_type;

		const std::string grid_name = ptree.get("mesh.name",(std::string)"continental_margin_default");
		auto grid_file = MODEL_PATH + PROBLEM_PATH ;
		grid_file += grid_name;
		grid_file += ".msh";
		Dune::GmshReader<GridType> gmshreader;
		std::shared_ptr<GridType> grid(gmshreader.read(grid_file,true,false));
		grid->loadBalance();

		using GV = GridType::LeafGridView;
		GV gv = grid->leafGridView();
	    using ES = Dune::PDELab::NonOverlappingEntitySet<GV>;
	    ES es(gv);

		//––––––––––––––––––––––––––––––––––––––––––
    	// MATERIAL PROPERTIES, NUMERICAL AND TEST PARAMETERS, CONSTITUTIVE RELATIONSHIPS
    	using Params = Parameters<ES,Dune::ParameterTree>;
    	Params params(es,ptree);

    	//––––––––––––––––––––––––––––––––––––––––––
		// DRIVER
		driver( es,
				ptree,
    			params,
				OUTPUT_PATH,
				helper );
    	//----------------------------------------------------------------------------------
	}
	catch (Dune::Exception &e){
		std::cerr << "Dune reported error: " << e << std::endl;
	}
	catch (...){
		std::cerr << "Unknown exception thrown!" << std::endl;
	}
}
