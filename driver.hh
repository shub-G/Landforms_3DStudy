/*
 * driver.hh
 *
 *  Created on: May 18, 2022
 *      Author: sgupta
 */

#ifndef LANDSLIDELANDSCAPE_DECOUPLED_NEW_DRIVER_HH_
#define LANDSLIDELANDSCAPE_DECOUPLED_NEW_DRIVER_HH_

template<typename ES,typename PTree,typename Parameters>
void driver( const ES& es,
			 const PTree& ptree,
			 Parameters& parameter,
			 std::string output_path,
			 Dune::MPIHelper& helper){

	std::string pathName = output_path;
	auto pathExt = ptree.get("output.path_name",(std::string)"test0");
	pathName += pathExt;
	pathName += "/";
	auto fileName = ptree.get("output.file_name",(std::string)"test");

	double time=0.;
	double dt = ptree.get("time_stepping.dt0",(double) 1000.0) ; /*years*/
	dt *= 86400.0*365.0;
	double dtstart = dt; //s
	double op_interval = ptree.get("time_stepping.output_interval",(double)1000.); //years
	op_interval *= 86400.0*365.0;
	double end_time = ptree.get("time_stepping.end_time",(double)120000.); //years
	end_time *= 86400.0*365.0;
	bool isAdaptive = ptree.get("time_stepping.adaptivity.flag",(bool)false);
	double dt_min = ptree.get("time_stepping.adaptivity.dt_min",(double)1.e-9);//years
	dt_min *= 86400.0*365.0;
	double dt_max = ptree.get("time_stepping.adaptivity.dt_max",(double)10.);//years
	dt_max *= 86400.0*365.0;
	int maxAllowableIterations = ptree.get("time_stepping.adaptivity.max_newton_steps",(int)6);
	int minAllowableIterations = ptree.get("time_stepping.adaptivity.min_newton_steps",(int)4);

	// prepare vtk writer
	using GV = typename ES::Traits::GridView;
	int subsampling = 1;
	using VTKWRITER = Dune::SubsamplingVTKWriter<GV>;
	VTKWRITER vtkwriter(es.gridView(),Dune::refinementIntervals(subsampling));
	using VTKSEQUENCEWRITER = Dune::VTKSequenceWriter<GV>;
	VTKSEQUENCEWRITER vtkSequenceWriter( std::make_shared<VTKWRITER>(vtkwriter),fileName,pathName,"");

	using Coord = typename GV::Grid::ctype;
	const int dim = GV::dimension;

	using NOCON0 = Dune::PDELab::ConformingDirichletConstraints;
	using CON0 = Dune::PDELab::ConformingDirichletConstraints;
	using VBE0 = Dune::PDELab::ISTL::VectorBackend<> ;
	using MBE = Dune::PDELab::ISTL::BCRSMatrixBackend<>;

	//-------------------------
	// FEM SPACE, GFS0, GFS
	//-------------------------
	// GOVERNING SYSTEM: FLOW (pressure and porosity)
	const int degree_f = 1;/*k-P1fem*/
	using FEM_F0 = Dune::PDELab::PkLocalFiniteElementMap<ES, Coord, double, degree_f>;
	FEM_F0 fem_f0(es);
	using GFS_F0 = Dune::PDELab::GridFunctionSpace<ES,FEM_F0,CON0,VBE0> ;
	GFS_F0 gfs_f0(es,fem_f0);
	using GFS_F = Dune::PDELab::PowerGridFunctionSpace< GFS_F0,
													    FlowVariables::num,
													    Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::none>,
													    Dune::PDELab::LexicographicOrderingTag >;
	GFS_F gfs_f(gfs_f0);

	// GOVERNING SYSTEM: ELASTICITY
	const int degree_u = 1;/*k-P1fem*/
	using FEM_U0 = Dune::PDELab::PkLocalFiniteElementMap<ES, Coord, double, degree_u>;
	FEM_U0 fem_u0(es);
	using GFS_U0 = Dune::PDELab::GridFunctionSpace<ES,FEM_U0,CON0,VBE0> ;
	GFS_U0 gfs_u0(es,fem_u0);
	using GFS_U = Dune::PDELab::PowerGridFunctionSpace< GFS_U0,
														DIMENSION,
														Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::none>,
														Dune::PDELab::LexicographicOrderingTag >;
	GFS_U gfs_u(gfs_u0);

	//
	// INVARIANTS OF STRESS, YIELD SURFACE
	auto gt = Dune::GeometryTypes::simplex(dim);
	using FEM_I = Dune::PDELab::P0LocalFiniteElementMap<Coord,double,dim>;
	FEM_I fem_i(gt);
	using GFS_I = Dune::PDELab::GridFunctionSpace<ES,FEM_I,NOCON0,VBE0> ;
	GFS_I gfs_i(es,fem_i);

	//-------------------
	// VECTOR CONTAINERS
	//-------------------
	// MAIN
	using F = Dune::PDELab::Backend::Vector<GFS_F,double>;
	F f_new( gfs_f,0.0 );
	F f_old( gfs_f,0.0 );
	F f_ini( gfs_f,0.0 );
	using U = Dune::PDELab::Backend::Vector<GFS_U,double>;
	U u_new( gfs_u,0.0 );
	// INVARIANTS OF STRESS, YIELD SURFACE
	using I = Dune::PDELab::Backend::Vector<GFS_I,double>;
	I meanstress( gfs_i,0.0 );
	I shearstress( gfs_i,0.0 );
	I yieldsurface( gfs_i,0.0 );

	//----------------------
	// INITIAL CONDITIONS
	//----------------------
	using Initial = InitialConditions<ES,PTree,Parameters>;
	Initial initial(es,ptree,parameter);
	auto pw_ic_local = [&](const auto& e, const auto& x){return initial.pressure(e,x);};
	auto pw_ic = Dune::PDELab::makeGridFunctionFromCallable(es,pw_ic_local);
	auto por_ic_local = [&](const auto& e, const auto& x){return initial.porosity(e,x);};
	auto por_ic = Dune::PDELab::makeGridFunctionFromCallable(es,por_ic_local);
	using  IC_F = Dune::PDELab::CompositeGridFunction<decltype(pw_ic),
													  decltype(por_ic)>;
	IC_F ic_f( pw_ic,por_ic );
	Dune::PDELab::interpolate( ic_f, gfs_f, f_ini );
	f_old=f_ini;
	f_new=f_ini;

	auto u0_ic_local = [&](const auto& e, const auto& x){return initial.displacement_u0(e,x);};
	auto u0_ic = Dune::PDELab::makeGridFunctionFromCallable(es,u0_ic_local);
	auto u1_ic_local = [&](const auto& e, const auto& x){return initial.displacement_u1(e,x);};
	auto u1_ic = Dune::PDELab::makeGridFunctionFromCallable(es,u1_ic_local);
#if DIMENSION==3
	auto u2_ic_local = [&](const auto& e, const auto& x){return initial.displacement_u2(e,x);};
	auto u2_ic = Dune::PDELab::makeGridFunctionFromCallable(es,u2_ic_local);
	using  IC_U = Dune::PDELab::CompositeGridFunction<decltype(u0_ic),
													  decltype(u1_ic),
													  decltype(u2_ic)>;
	IC_U ic_u( u0_ic,u1_ic,u2_ic );
#elif DIMENSION==2
	using  IC_U = Dune::PDELab::CompositeGridFunction<decltype(u0_ic),
													  decltype(u1_ic)>;
	IC_U ic_u( u0_ic,u1_ic );
#endif

	// 	Initialize the solution at t=0 with the given initial values
	Dune::PDELab::interpolate( ic_u, gfs_u, u_new );

	//----------------------
	// BOUNDARY CONDITIONS
	//----------------------
	using BCT_F0 = BoundaryTypesFLOW<ES,PTree,Parameters> ;
	BCT_F0 bct_pw( es,ptree,parameter,FlowVariables::pw,&time,&dt );
	BCT_F0 bct_por(es,ptree,parameter,FlowVariables::porosity,&time,&dt );
	using BCT_F = Dune::PDELab::CompositeConstraintsParameters<BCT_F0,BCT_F0>;
	BCT_F bct_f( bct_pw, bct_por );
	using C_F = typename GFS_F::template ConstraintsContainer<double>::Type;
	C_F c_f;
	c_f.clear();
	Dune::PDELab::constraints( bct_f, gfs_f, c_f, false);  // to artificial boundaries
	if(helper.rank()==0){
	std::cout << "constrained dofs (FLOW) = " << c_f.size() << " of " << gfs_f.globalSize() << std::endl;
	}

	using BCV_F = NeumannBoundaryValuesFLOW<ES,PTree,Parameters>;
	BCV_F bcv_f( es,ptree,parameter ) ;

	using BCVEXT_F0 = DirichletBoundaryValuesFLOW<ES,PTree,Parameters> ;
	BCVEXT_F0 bcvext_pw( es,ptree,parameter,FlowVariables::pw,&time,&dt ) ;
	BCVEXT_F0 bcvext_por(es,ptree,parameter,FlowVariables::porosity,&time,&dt ) ;
	typedef Dune::PDELab::CompositeGridFunction< BCVEXT_F0,
												 BCVEXT_F0> BCVEXT_F;
	BCVEXT_F bcvext_f( bcvext_pw, bcvext_por );
	Dune::PDELab::interpolate( bcvext_f, gfs_f, f_new );

	using BCT_U0 = BoundaryTypesELASTICITY<ES,PTree,Parameters> ;
	BCT_U0 bct_u0( es,ptree,parameter,0,&time,&dt );
	BCT_U0 bct_u1( es,ptree,parameter,1,&time,&dt );
#if DIMENSION==3
	BCT_U0 bct_u2( es,ptree,parameter,2,&time,&dt );
	using BCT_U = Dune::PDELab::CompositeConstraintsParameters<BCT_U0,BCT_U0,BCT_U0>;
	BCT_U bct_u( bct_u0, bct_u1, bct_u2 );
#elif DIMENSION==2
	using BCT_U = Dune::PDELab::CompositeConstraintsParameters<BCT_U0,BCT_U0>;
	BCT_U bct_u( bct_u0, bct_u1 );
#endif
	using C_U = typename GFS_U::template ConstraintsContainer<double>::Type;
	C_U c_u;
	c_u.clear();
	Dune::PDELab::constraints( bct_u, gfs_u, c_u, false);  // to artificial boundaries
	if(helper.rank()==0){
	std::cout << "constrained dofs (ELASTICITY) = " << c_u.size() << " of " << gfs_u.globalSize() << std::endl;
	}

	using BCV_U = NeumannBoundaryValuesELASTICITY<ES,PTree,Parameters>;
	BCV_U bcv_u( es,ptree,parameter ) ;

	using BCVEXT_U0 = DirichletBoundaryValuesELASTICITY<ES,PTree,Parameters> ;
	BCVEXT_U0 bcvext_u0( es,ptree,parameter,0,&time,&dt ) ;
	BCVEXT_U0 bcvext_u1( es,ptree,parameter,1,&time,&dt ) ;
#if DIMENSION==3
	BCVEXT_U0 bcvext_u2( es,ptree,parameter,2,&time,&dt ) ;
	typedef Dune::PDELab::CompositeGridFunction< BCVEXT_U0,
												 BCVEXT_U0,
												 BCVEXT_U0> BCVEXT_U;
	BCVEXT_U bcvext_u( bcvext_u0, bcvext_u1, bcvext_u2 );
#elif DIMENSION==2
	typedef Dune::PDELab::CompositeGridFunction< BCVEXT_U0,
												 BCVEXT_U0> BCVEXT_U;
	BCVEXT_U bcvext_u( bcvext_u0, bcvext_u1 );
#endif
	Dune::PDELab::interpolate( bcvext_u, gfs_u, u_new );

	//-------------------
	// LOCAL OPERATOR:
	// FLOW
	//-------------------
	unsigned int intorder = ptree.get("cg_parameters.intorder",(int)4);
	using LOP_F = LocalOperatorFLOW<ES,GFS_F,F,GFS_I,I,Parameters,BCT_F,BCV_F>;
	LOP_F lop_f(es,
				gfs_f,&f_old,
				gfs_i,&yieldsurface,
				parameter,
				bct_f,bcv_f,
				&time,&dt,
				intorder);

	//-------------------
	// LOCAL OPERATOR:
	// ELASTICITY
	//-------------------
	using LOP_U = LocalOperatorELASTICITY<ES,GFS_F,F,Parameters,BCT_U,BCV_U>;
	LOP_U lop_u(es,
				gfs_f,&f_new,&f_ini,
				parameter,
				bct_u,bcv_u,
				&time,&dt,
				intorder);

	//-------------------
	// GRID OPERATOR
	// FLOW
	//-------------------
	MBE mbe_f(20); // Maximal number of nonzeroes per row
	using GOLOP_F = Dune::PDELab::GridOperator<GFS_F,GFS_F,LOP_F,MBE,double,double,double,C_F,C_F>;
	GOLOP_F golop_f(gfs_f,c_f,gfs_f,c_f,lop_f,mbe_f);
	typename GOLOP_F::Traits::Jacobian jac_f(golop_f);
	if(helper.rank()==0){
		std::cout << jac_f.patternStatistics() << std::endl;
	}

	//-------------------
	// GRID OPERATOR
	//-------------------
	MBE mbe_u(20); // Maximal number of nonzeroes per row
	using GOLOP_U = Dune::PDELab::GridOperator<GFS_U,GFS_U,LOP_U,MBE,double,double,double,C_U,C_U>;
	GOLOP_U golop_u(gfs_u,c_u,gfs_u,c_u,lop_u,mbe_u);
	typename GOLOP_U::Traits::Jacobian jac_u(golop_u);
	if(helper.rank()==0){
		std::cout << jac_u.patternStatistics() << std::endl;
	}

	//-------------------
	// LINEAR SOLVER
	// FLOW
	//-------------------
	using LS_F = Dune::PDELab::ISTLBackend_NOVLP_BCGS_AMG_SSOR<GOLOP_F>;
	LS_F ls_f(golop_f,1000,1,true,true);
	Dune::Amg::Parameters ls_params_f = ls_f.parameters();
	ls_params_f.setCoarsenTarget(1000000);// max DoF at coarsest level
	ls_f.setParameters(ls_params_f);

	//-------------------
	// LINEAR SOLVER
	// ELASTICITY
	//-------------------
	using LS_U = Dune::PDELab::ISTLBackend_NOVLP_BCGS_AMG_SSOR<GOLOP_U>;
	LS_U ls_u(golop_u,1000,1,true,true);
	Dune::Amg::Parameters ls_params_u = ls_u.parameters();
	ls_params_u.setCoarsenTarget(1000000);// max DoF at coarsest level
	ls_u.setParameters(ls_params_u);

	//-------------------
	// NON-LINEAR SOLVER
	// FLOW
	//-------------------
	using PDESOLVER_F = Dune::PDELab::Newton<GOLOP_F,LS_F,F>;
	PDESOLVER_F pdesolver_f(golop_f,ls_f);
	bool line_search_flag_f = ptree.get("newton.flow.line_search",(bool)true);
	if(!line_search_flag_f) pdesolver_f.setLineSearchStrategy(PDESOLVER_F::Strategy::noLineSearch);
	pdesolver_f.setReassembleThreshold(0.0);
	pdesolver_f.setVerbosityLevel(ptree.get("newton.flow.verbosity",(int)2));
	pdesolver_f.setReduction(ptree.get("newton.flow.reduction",(double)1e-6));
	pdesolver_f.setMinLinearReduction(ptree.get("newton.flow.min_lin_reduction",(double)1e-9));
	pdesolver_f.setMaxIterations(ptree.get("newton.flow.max_iterations",(int)15));
	pdesolver_f.setForceIteration(ptree.get("newton.flow.force_iteration",(bool)false));
	pdesolver_f.setAbsoluteLimit(ptree.get("newton.flow.abs_error",(double)1.e-6));

	//-------------------
	// NON-LINEAR SOLVER
	// ELASTICITY
	//-------------------
	using PDESOLVER_U = Dune::PDELab::Newton<GOLOP_U,LS_U,U>;
	PDESOLVER_U pdesolver_u(golop_u,ls_u);
	bool line_search_flag_u = ptree.get("newton.elasticity.line_search",(bool)true);
	if(!line_search_flag_u) pdesolver_u.setLineSearchStrategy(PDESOLVER_U::Strategy::noLineSearch);
	pdesolver_u.setReassembleThreshold(0.0);
	pdesolver_u.setVerbosityLevel(ptree.get("newton.elasticity.verbosity",(int)2));
	pdesolver_u.setReduction(ptree.get("newton.elasticity.reduction",(double)1e-6));
	pdesolver_u.setMinLinearReduction(ptree.get("newton.elasticity.min_lin_reduction",(double)1e-9));
	pdesolver_u.setMaxIterations(ptree.get("newton.elasticity.max_iterations",(int)15));
	pdesolver_u.setForceIteration(ptree.get("newton.elasticity.force_iteration",(bool)false));
	pdesolver_u.setAbsoluteLimit(ptree.get("newton.elasticity.abs_error",(double)1.e-6));

	//-------------------
	// POSTPROCESS
	//-------------------
	PostProcess postprocess;
	postprocess.evaluate_meanstress(es,parameter,gfs_u,u_new,gfs_i,&meanstress);
	postprocess.evaluate_shearstress(es,parameter,gfs_u,u_new,gfs_i,&shearstress);
	postprocess.evaluate_yieldsurface(es,parameter,gfs_u,u_new,gfs_i,&yieldsurface);

	//-------------------
	// UPDATE F_INI and U_INI
	//-------------------
//	pdesolver_u.apply(u_new);
//	postprocess.update_initial_conditions(es,parameter,gfs,u_new,u_old,&u_ini);
	//f_ini=f_old;
	//f_new=f_ini;

	//-------------------
	// VTK OUTPUT
	//-------------------
	// Add names to the components for VTK output
	using namespace Dune::TypeTree::Indices;
	gfs_f.child(FlowVariables::pw).name("pressure");
	gfs_f.child(FlowVariables::porosity).name("porosity");
	Dune::PDELab::addSolutionToVTKWriter(vtkSequenceWriter,gfs_f,f_new);

	for(int i=0; i<DIMENSION; i++){
		gfs_u.child(i).name("u_"+std::to_string(i));
	}
	Dune::PDELab::addSolutionToVTKWriter(vtkSequenceWriter,gfs_u,u_new);

	gfs_i.name("stress_p");
	Dune::PDELab::addSolutionToVTKWriter(vtkSequenceWriter,gfs_i,meanstress);

	gfs_i.name("stress_q");
	Dune::PDELab::addSolutionToVTKWriter(vtkSequenceWriter,gfs_i,shearstress);

	gfs_i.name("yield_surface");
	Dune::PDELab::addSolutionToVTKWriter(vtkSequenceWriter,gfs_i,yieldsurface);

	vtkSequenceWriter.write(time,Dune::VTK::appendedraw);

	Dune::PDELab::constraints( bct_f, gfs_f, c_f, false);  // to artificial boundaries
	Dune::PDELab::interpolate( bcvext_f, gfs_f , f_new );

	Dune::PDELab::constraints( bct_u, gfs_u, c_u, false);  // to artificial boundaries
	Dune::PDELab::interpolate( bcvext_u, gfs_u , u_new );

	/****************************************************/
	// TIME-LOOP
	/****************************************************/
	int opcount = 1;
	double timecount = time;
	double dtLast = dtstart;
	int dtFlag = 0;
	bool exceptionCaught = false;
	int newton_iterations = 0;
	int max_iter=ptree.get("FP.max_iterations",(int)5);

	while( time <= end_time ){

		if( exceptionCaught==false ){
			dt = std::max(dt,dt_min);
		}

		if(helper.rank()==0){
		std::cout<< "_____________________________________________________" <<std::endl;
		std::cout<< " current opcount = " << opcount - 1 << std::endl;
		}
		try{
			
			for( int iter=0; iter<max_iter; iter++){
				std::cout<< "+++++++ ITER: " << iter << " ++++++++++" << std::endl;
				// FLOW
				std::cout<< "SOLVING FLOW PROBLEM" << std::endl;
				pdesolver_f.apply(f_new);
				int newton_iterations_f = pdesolver_f.result().iterations;
				
				// ELASTICITY
				std::cout<< "SOLVING ELASTICITY PROBLEM" << std::endl;
				pdesolver_u.apply(u_new);
				int newton_iterations_u = pdesolver_u.result().iterations;
				
				if(iter==0) newton_iterations = std::max(newton_iterations_f,newton_iterations_u);

				// YIELDSURFACE
				std::cout<< "COMPUTING YIELD SURFACE" << std::endl;
				postprocess.evaluate_yieldsurface(es,parameter,gfs_u,u_new,gfs_i,&yieldsurface);
				
				if(newton_iterations_f==0) iter=max_iter;
				
				std::cout<< "------------" << std::endl;
				exceptionCaught = false;
			}

		}catch ( Dune::Exception &e ) {
			exceptionCaught = true;
			if( dt > 1e-15 ){
				if(helper.rank()==0){
				std::cout << "Catched Error, Dune reported error: " << e << std::endl;
				}
				f_new = f_old;

				Dune::PDELab::constraints( bct_f, gfs_f, c_f, false);  // to artificial boundaries
				Dune::PDELab::interpolate( bcvext_f, gfs_f , f_new );

				Dune::PDELab::constraints( bct_u, gfs_u, c_u, false);  // to artificial boundaries
				Dune::PDELab::interpolate( bcvext_u, gfs_u , u_new );

				newton_iterations = 0;
				dt *= 0.5;
				continue;
			}
			else
			{
				if(helper.rank()==0){
					std::cout << "ABORTING, due to DUNE error: " << e << std::endl;
				}
				exit(0);
			}
		}
		if(helper.rank()==0){
		std::cout<<"DONE"<<std::endl;
		std::cout<<"_____________________________________________________"<<std::endl;
		}

		/***********************************************
		 * OUTPUT
		 ***********************************************/
		/* At op_interval
		 *
		 */
		if( time+dt == op_interval*opcount )
		{
			// EVALUATE SECONDARY VARS
			postprocess.evaluate_meanstress(es,parameter,gfs_u,u_new,gfs_i,&meanstress);
			postprocess.evaluate_shearstress(es,parameter,gfs_u,u_new,gfs_i,&shearstress);
//			postprocess.evaluate_yieldsurface(es,parameter,gfs_u,u_new,gfs_i,&yieldsurface);

			// WRITE OUTPUT
			vtkSequenceWriter.write(time,Dune::VTK::appendedraw);

			if(helper.rank()==0){
			std::cout<< " ******************************************************************* " << std::endl;
			std::cout<< " OUTPUT WRITTEN " << opcount << std::endl;
			std::cout<< " ******************************************************************* " << std::endl;
			std::cout<< std::flush;
			}

			timecount = time;
			opcount = opcount+1;
		}

		/***********************************************/
		// PREPARE FOR NEXT TIME INTEGRATION
		/***********************************************/
		//1. ASSIGN THE 'NEXT' VALUE TO 'OLD' VARIABLE
		f_old = f_new;
		//2. ADVANCE TIME:
		time += dt;
		//3. UPDATE DIRICHLET CONSTRAINTS
		Dune::PDELab::constraints( bct_f, gfs_f, c_f, false);  // to artificial boundaries
		Dune::PDELab::interpolate( bcvext_f, gfs_f , f_new );
		Dune::PDELab::constraints( bct_u, gfs_u, c_u, false);  // to artificial boundaries
		Dune::PDELab::interpolate( bcvext_u, gfs_u , u_new );

		if(helper.rank()==0){
		std::cout<<" "<< std::endl;
		std::cout<< " time = " << time ;
		std::cout<< std::flush;
		}
		
		if( isAdaptive ){
			if(newton_iterations>maxAllowableIterations){
				dt=std::max(dt*0.75,dt_min);
			}
			else if(newton_iterations<=minAllowableIterations){
				dt=std::min(dt*1.25,dt_max);
			}
		}
		else{
			dt = dtstart;
		}

		if(helper.rank()==0){
		std::cout << " , time+dt = " << (time + dt)
				  << " , opTime = "  << op_interval * opcount ;
		std::cout<< std::flush;
		}

		if( time + dt  > op_interval * opcount){
			dtLast = dt;
			dt 	 = op_interval * opcount - time ;

			if(helper.rank()==0){
			std::cout<< " , because timeNext > opNext , dt set to : " << dt << std::endl;
			std::cout<< std::flush;
			}

			dtFlag = 0;
			exceptionCaught = true;
		}
		dtFlag += 1;

		if(helper.rank()==0){
		std::cout<< " , dt  : " << dt << std::endl;
		std::cout<<" "<< std::endl;
		std::cout << " READY FOR NEXT ITERATION. " << std::endl;
		std::cout<< std::flush;
		}

	}//END: time-loop

}

#endif /* LANDSLIDELANDSCAPE_DECOUPLED_NEW_DRIVER_HH_ */
