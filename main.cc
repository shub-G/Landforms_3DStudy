/*
 * main.cc
 *
 *  	CREATED ON: May 18, 2022
 *      AUTHOR: Shubhangi Gupta
 *	AFFIL: GEOMAR HElmholtz Center for Ocean Research Kiel, Germany
 *	CONTACT: sgupta@geomar.de
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
 * SELECT DIMENSION
 ***************************************************/
#define DIMENSION 3

/***************************************************
 * INCLUDE SOURCES AND OPERATORS
 ***************************************************/
#include"duneincludes.hh"
#include"parameters.hh"
#include"postprocess.hh"
#include"lop_elasticity.hh"
#include"lop_flow.hh"
#include"driver.hh"


int main(int argc, char** argv)
{
	try{
		//----------------------------------------------------------------------------------
		//––––––––––––––––––––––––––––––––––––––––––
	    	// Maybe initialize MPI
	    	Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
	    	if(helper.rank()==0){
	    	std::cout << "Hello World! This is PROJECT:LANDSLIDELANDSCAPE." << std::endl;
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
	    		std::cout << "usage: ./main <input-file> " << std::endl;
	    		}
	        	return 1;
	    	}

		//––––––––––––––––––––––––––––––––––––––––––
		// DUNE MODEL PATH
	    	std::string MODEL_PATH = __PATH__TO__SRC__FOLDER__;
	    	// INPUT PATH
	    	std::string INPUT_PATH = MODEL_PATH + "/inputs/" ;
	    	// OUTPUT PATH
	    	std::string OUTPUT_PATH = MODEL_PATH + "/outpts/" ;
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
