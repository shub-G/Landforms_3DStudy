/*
 * postprocess.hh
 *
 *  Created on: May 18, 2022
 *      Author: sgupta
 */

#ifndef LANDSLIDELANDSCAPE_DECOUPLED_NEW_OPERATORS_POSTPROCESS_HH_
#define LANDSLIDELANDSCAPE_DECOUPLED_NEW_OPERATORS_POSTPROCESS_HH_

class PostProcess{
private:
	unsigned int intorder = 4;

public:

	template<class GV, class Parameters, class GFSU, class U, class GFSI, class I>
	void evaluate_meanstress(const GV& gv,
							 const Parameters& param,
							 GFSU gfsu,
							 U unew,
							 GFSI gfsi,
							 I *stress_p){

	    typedef Dune::PDELab::LocalFunctionSpace< GFSU > LFSU;
	    LFSU lfsu(gfsu);
	    typedef Dune::PDELab::LFSIndexCache<LFSU> LFSCacheU;
		LFSCacheU lfs_cache_u(lfsu);
		typedef typename U::template LocalView<LFSCacheU> VectorViewU;
		VectorViewU u_view( unew );

	    typedef Dune::PDELab::LocalFunctionSpace< GFSI > LFSI;
	    LFSI lfsi(gfsi);
	    typedef Dune::PDELab::LFSIndexCache<LFSI> LFSCacheI;
		LFSCacheI lfs_cache_i(lfsi);
		typedef typename I::template LocalView<LFSCacheI> VectorViewI;
		VectorViewI i_view( *stress_p );

		typedef typename GV::Traits::template Codim<0>::Iterator LeafIterator;
		typedef typename GV::Grid::template Codim<0>::Entity Element;
		typedef typename GV::IntersectionIterator IntersectionIterator;
		typedef typename IntersectionIterator::Intersection Intersection;
	    typedef typename GV::IndexSet IndexSet;
	    static const int dim = GV::dimension;

	    typedef typename LFSU::template Child<0>::Type::Traits::FiniteElementType::
	    		Traits::LocalBasisType::Traits::RangeType RangeType;
		typedef typename LFSU::template Child<0>::Type::Traits::FiniteElementType::
				Traits::LocalBasisType::Traits::JacobianType JacobianType;
		typedef typename LFSU::Traits::SizeType SizeType;

		// Iterate over each element
		LeafIterator beginElem = gv.template begin< 0 >();
		LeafIterator endElem = gv.template end< 0 >();

		for ( LeafIterator self = beginElem; self!= endElem; ++self )
		{
			/***********************************************/
			// Get Element geometry
			typedef typename LeafIterator::Entity::Geometry ElementGeometry;
			const auto& geo = (*self).geometry();
			auto ref_el = referenceElement(geo);
			auto center_local = ref_el.position(0,0);
			auto center_global = geo.global(center_local);

	        /***********************************************/
	        // save local coeffs in map
	        lfsu.bind(*self);
	        lfsi.bind(*self);

			lfs_cache_u.update();
			lfs_cache_i.update();

			u_view.bind(lfs_cache_u);
			i_view.bind(lfs_cache_i);

		    std::vector<double> vl_unew(lfsu.size());
		    u_view.read( vl_unew );
	        std::vector<double> vl_stress_p(lfsi.size());
	        for(int i = 0. ; i < lfsi.size() ; i++){
	        	vl_stress_p[i] = 0.;
	        }

		    // transform gradients from reference element to real element
		    typename Element::Geometry::JacobianInverseTransposed jac = (*self).geometry().jacobianInverseTransposed( center_local );

			std::vector< std::vector<RangeType> > phi_u( dim );
			std::vector< std::vector<JacobianType> > js_u( dim );
			std::vector< std::vector< Dune::FieldVector<double,dim> > > gradphi_u( dim );
			std::vector< Dune::FieldVector<double,dim > > grad_u( dim );
			for( int d=0; d<dim; d++ ){

				phi_u[d] = std::vector<RangeType> ( lfsu.child(d).size() );
				lfsu.child(d).finiteElement().localBasis().evaluateFunction( center_local,phi_u[d] );

				js_u[d] = std::vector<JacobianType> ( lfsu.child(d).size() );
				lfsu.child(d).finiteElement().localBasis().evaluateJacobian( center_local,js_u[d] );

				gradphi_u[d] = std::vector<Dune::FieldVector<double,dim>> ( lfsu.child(d).size() );
				for (SizeType i=0; i < lfsu.child(d).size(); i++){
					jac.mv(js_u[d][i][0],gradphi_u[d][i]);
				}

				grad_u[d] =  Dune::FieldVector<double,dim >  (0.0);
				for (SizeType i=0; i<lfsu.child(d).size(); i++){
					grad_u[d].axpy(vl_unew[lfsu.child(d).localIndex(i)],gradphi_u[d][i]);
				}
			}

			//evaluate stress
			// stress = 2 mu strain + lambda tr(strain) I
			// LP[0] = lambda
			// LP[1] = mu
			auto LP = param.LameParameters(center_global);
			Dune::FieldMatrix<double,dim,dim> strain(0.);
			Dune::FieldMatrix<double,dim,dim> Id(0.);
			for(int i=0; i<dim; i++ ){
				Id[i][i] = 1.0;
				for(int j=0; j<dim; j++ ){
					strain[i][j] = (1./2.)*( grad_u[i][j] + grad_u[j][i] );
				}
			}
			double tr_strain = 0.0;
			for(int i=0; i<dim; i++ ) tr_strain += strain[i][i];
			Dune::FieldMatrix<double,dim,dim> stress(0.);
			for(int i=0; i<dim; i++ ){
				for(int j=0; j<dim; j++ ){
					stress[i][j] = 2.*LP[1]*strain[i][j] + LP[0]*tr_strain*Id[i][j];
				}
			}

			//stress vector
			std::vector<double> stress_vector(StressIndices::num,0.0);
			#if DIMENSION==3
			stress_vector[StressIndices::S00] = stress[0][0];
			stress_vector[StressIndices::S11] = stress[1][1];
			stress_vector[StressIndices::S22] = stress[2][2];
			stress_vector[StressIndices::S01] = stress[0][1];
			stress_vector[StressIndices::S12] = stress[1][2];
			stress_vector[StressIndices::S20] = stress[2][0];
			#elif DIMENSION==2
			stress_vector[StressIndices::S00] = stress[0][0];
			stress_vector[StressIndices::S11] = stress[1][1];
			stress_vector[StressIndices::S01] = stress[0][1];
			#endif

			// evaluate meanstress
			vl_stress_p[lfsi.localIndex(0)] = param.Stress_p(stress_vector);

			i_view.write( vl_stress_p );
			i_view.commit();
			i_view.unbind();

		}//END: iterate over elements

	}//END: evaluate_meanstress



	template<class GV, class Parameters, class GFSU, class U, class GFSI, class I>
	void evaluate_shearstress(const GV& gv,
							  const Parameters& param,
							  GFSU gfsu,
							  U unew,
							  GFSI gfsi,
							  I *stress_q){

	    typedef Dune::PDELab::LocalFunctionSpace< GFSU > LFSU;
	    LFSU lfsu(gfsu);
	    typedef Dune::PDELab::LFSIndexCache<LFSU> LFSCacheU;
		LFSCacheU lfs_cache_u(lfsu);
		typedef typename U::template LocalView<LFSCacheU> VectorViewU;
		VectorViewU u_view( unew );

	    typedef Dune::PDELab::LocalFunctionSpace< GFSI > LFSI;
	    LFSI lfsi(gfsi);
	    typedef Dune::PDELab::LFSIndexCache<LFSI> LFSCacheI;
		LFSCacheI lfs_cache_i(lfsi);
		typedef typename I::template LocalView<LFSCacheI> VectorViewI;
		VectorViewI i_view( *stress_q );

		typedef typename GV::Traits::template Codim<0>::Iterator LeafIterator;
		typedef typename GV::Grid::template Codim<0>::Entity Element;
		typedef typename GV::IntersectionIterator IntersectionIterator;
		typedef typename IntersectionIterator::Intersection Intersection;
	    typedef typename GV::IndexSet IndexSet;
	    static const int dim = GV::dimension;

	    typedef typename LFSU::template Child<0>::Type::Traits::FiniteElementType::
	    		Traits::LocalBasisType::Traits::RangeType RangeType;
		typedef typename LFSU::template Child<0>::Type::Traits::FiniteElementType::
				Traits::LocalBasisType::Traits::JacobianType JacobianType;
		typedef typename LFSU::Traits::SizeType SizeType;

		// Iterate over each element
		LeafIterator beginElem = gv.template begin< 0 >();
		LeafIterator endElem = gv.template end< 0 >();

		for ( LeafIterator self = beginElem; self!= endElem; ++self )
		{
			/***********************************************/
			// Get Element geometry
			typedef typename LeafIterator::Entity::Geometry ElementGeometry;
			const auto& geo = (*self).geometry();
//			auto center_local  = geo.local( geo.center());
//			auto center_global = geo.global(x);//geo.global(geo.center());
			auto ref_el = referenceElement(geo);
			auto center_local = ref_el.position(0,0);
			auto center_global = geo.global(center_local);

	        /***********************************************/
	        // save local coeffs in map
	        lfsu.bind(*self);
	        lfsi.bind(*self);

			lfs_cache_u.update();
			lfs_cache_i.update();

			u_view.bind(lfs_cache_u);
			i_view.bind(lfs_cache_i);

		    std::vector<double> vl_unew(lfsu.size());
		    u_view.read( vl_unew );
	        std::vector<double> vl_stress_q(lfsi.size());
	        for(int i = 0. ; i < lfsi.size() ; i++){
	        	vl_stress_q[i] = 0.;
	        }

		    // transform gradients from reference element to real element
		    typename Element::Geometry::JacobianInverseTransposed jac = (*self).geometry().jacobianInverseTransposed( center_local );

			std::vector< std::vector<RangeType> > phi_u( dim );
			std::vector< std::vector<JacobianType> > js_u( dim );
			std::vector< std::vector< Dune::FieldVector<double,dim> > > gradphi_u( dim );
			std::vector< Dune::FieldVector<double,dim > > grad_u( dim );
			for( int d=0; d<dim; d++ ){

				phi_u[d] = std::vector<RangeType> ( lfsu.child(d).size() );
				lfsu.child(d).finiteElement().localBasis().evaluateFunction( center_local,phi_u[d] );

				js_u[d] = std::vector<JacobianType> ( lfsu.child(d).size() );
				lfsu.child(d).finiteElement().localBasis().evaluateJacobian( center_local,js_u[d] );

				gradphi_u[d] = std::vector<Dune::FieldVector<double,dim>> ( lfsu.child(d).size() );
				for (SizeType i=0; i < lfsu.child(d).size(); i++){
					jac.mv(js_u[d][i][0],gradphi_u[d][i]);
				}

				grad_u[d] =  Dune::FieldVector<double,dim >  (0.0);
				for (SizeType i=0; i<lfsu.child(d).size(); i++){
					grad_u[d].axpy(vl_unew[lfsu.child(d).localIndex(i)],gradphi_u[d][i]);
				}
			}

			//evaluate stress
			// stress = 2 mu strain + lambda tr(strain) I
			// LP[0] = lambda
			// LP[1] = mu
			auto LP = param.LameParameters(center_global);
			Dune::FieldMatrix<double,dim,dim> strain(0.);
			Dune::FieldMatrix<double,dim,dim> Id(0.);
			for(int i=0; i<dim; i++ ){
				Id[i][i] = 1.0;
				for(int j=0; j<dim; j++ ){
					strain[i][j] = (1./2.)*( grad_u[i][j] + grad_u[j][i] );
				}
			}
			double tr_strain = 0.0;
			for(int i=0; i<dim; i++ ) tr_strain += strain[i][i];
			Dune::FieldMatrix<double,dim,dim> stress(0.);
			for(int i=0; i<dim; i++ ){
				for(int j=0; j<dim; j++ ){
					stress[i][j] = 2.*LP[1]*strain[i][j] + LP[0]*tr_strain*Id[i][j];
				}
			}

			//stress vector
			std::vector<double> stress_vector(StressIndices::num,0.0);
			#if DIMENSION==3
			stress_vector[StressIndices::S00] = stress[0][0];
			stress_vector[StressIndices::S11] = stress[1][1];
			stress_vector[StressIndices::S22] = stress[2][2];
			stress_vector[StressIndices::S01] = stress[0][1];
			stress_vector[StressIndices::S12] = stress[1][2];
			stress_vector[StressIndices::S20] = stress[2][0];
			#elif DIMENSION==2
			stress_vector[StressIndices::S00] = stress[0][0];
			stress_vector[StressIndices::S11] = stress[1][1];
			stress_vector[StressIndices::S01] = stress[0][1];
			#endif

			// evaluate shearstress
			vl_stress_q[lfsi.localIndex(0)] = param.Stress_q(stress_vector);

			i_view.write( vl_stress_q );
			i_view.commit();
			i_view.unbind();

		}//END: iterate over elements

	}//END: evaluate_shearstress



	template<class GV, class Parameters, class GFSU, class U, class GFSI, class I>
	void evaluate_yieldsurface(const GV& gv,
							   const Parameters& param,
							   GFSU gfsu,
							   U unew,
							   GFSI gfsi,
							   I *F){

	    typedef Dune::PDELab::LocalFunctionSpace< GFSU > LFSU;
	    LFSU lfsu(gfsu);
	    typedef Dune::PDELab::LFSIndexCache<LFSU> LFSCacheU;
		LFSCacheU lfs_cache_u(lfsu);
		typedef typename U::template LocalView<LFSCacheU> VectorViewU;
		VectorViewU u_view( unew );

	    typedef Dune::PDELab::LocalFunctionSpace< GFSI > LFSI;
	    LFSI lfsi(gfsi);
	    typedef Dune::PDELab::LFSIndexCache<LFSI> LFSCacheI;
		LFSCacheI lfs_cache_i(lfsi);
		typedef typename I::template LocalView<LFSCacheI> VectorViewI;
		VectorViewI i_view( *F );

		typedef typename GV::Traits::template Codim<0>::Iterator LeafIterator;
		typedef typename GV::Grid::template Codim<0>::Entity Element;
		typedef typename GV::IntersectionIterator IntersectionIterator;
		typedef typename IntersectionIterator::Intersection Intersection;
	    typedef typename GV::IndexSet IndexSet;
	    static const int dim = GV::dimension;

	    typedef typename LFSU::template Child<0>::Type::Traits::FiniteElementType::
	    		Traits::LocalBasisType::Traits::RangeType RangeType;
		typedef typename LFSU::template Child<0>::Type::Traits::FiniteElementType::
				Traits::LocalBasisType::Traits::JacobianType JacobianType;
		typedef typename LFSU::Traits::SizeType SizeType;

		// Iterate over each element
		LeafIterator beginElem = gv.template begin< 0 >();
		LeafIterator endElem = gv.template end< 0 >();

		for ( LeafIterator self = beginElem; self!= endElem; ++self )
		{
			/***********************************************/
			// Get Element geometry
			typedef typename LeafIterator::Entity::Geometry ElementGeometry;
			const auto& geo = (*self).geometry();
//			auto center_local  = geo.local( geo.center());
//			auto center_global = geo.global(x);//geo.global(geo.center());
			auto ref_el = referenceElement(geo);
			auto center_local = ref_el.position(0,0);
			auto center_global = geo.global(center_local);

	        /***********************************************/
	        // save local coeffs in map
	        lfsu.bind(*self);
	        lfsi.bind(*self);

			lfs_cache_u.update();
			lfs_cache_i.update();

			u_view.bind(lfs_cache_u);
			i_view.bind(lfs_cache_i);

		    std::vector<double> vl_unew(lfsu.size());
		    u_view.read( vl_unew );
	        std::vector<double> vl_F(lfsi.size());
	        for(int i = 0. ; i < lfsi.size() ; i++){
	        	vl_F[i] = 0.;
	        }

		    // transform gradients from reference element to real element
		    typename Element::Geometry::JacobianInverseTransposed jac = (*self).geometry().jacobianInverseTransposed( center_local );

			std::vector< std::vector<RangeType> > phi_u( dim );
			std::vector< std::vector<JacobianType> > js_u( dim );
			std::vector< std::vector< Dune::FieldVector<double,dim> > > gradphi_u( dim );
			std::vector< Dune::FieldVector<double,dim > > grad_u( dim );
			for( int d=0; d<dim; d++ ){

				phi_u[d] = std::vector<RangeType> ( lfsu.child(d).size() );
				lfsu.child(d).finiteElement().localBasis().evaluateFunction( center_local,phi_u[d] );

				js_u[d] = std::vector<JacobianType> ( lfsu.child(d).size() );
				lfsu.child(d).finiteElement().localBasis().evaluateJacobian( center_local,js_u[d] );

				gradphi_u[d] = std::vector<Dune::FieldVector<double,dim>> ( lfsu.child(d).size() );
				for (SizeType i=0; i < lfsu.child(d).size(); i++){
					jac.mv(js_u[d][i][0],gradphi_u[d][i]);
				}

				grad_u[d] =  Dune::FieldVector<double,dim >  (0.0);
				for (SizeType i=0; i<lfsu.child(d).size(); i++){
					grad_u[d].axpy(vl_unew[lfsu.child(d).localIndex(i)],gradphi_u[d][i]);
				}
			}

			//evaluate stress
			// stress = 2 mu strain + lambda tr(strain) I
			// LP[0] = lambda
			// LP[1] = mu
			auto LP = param.LameParameters(center_global);
			Dune::FieldMatrix<double,dim,dim> strain(0.);
			Dune::FieldMatrix<double,dim,dim> Id(0.);
			for(int i=0; i<dim; i++ ){
				Id[i][i] = 1.0;
				for(int j=0; j<dim; j++ ){
					strain[i][j] = (1./2.)*( grad_u[i][j] + grad_u[j][i] );
				}
			}
			double tr_strain = 0.0;
			for(int i=0; i<dim; i++ ) tr_strain += strain[i][i];
			Dune::FieldMatrix<double,dim,dim> stress(0.);
			for(int i=0; i<dim; i++ ){
				for(int j=0; j<dim; j++ ){
					stress[i][j] = 2.*LP[1]*strain[i][j] + LP[0]*tr_strain*Id[i][j];
				}
			}

			//stress vector
			std::vector<double> stress_vector(StressIndices::num,0.0);
			#if DIMENSION==3
			stress_vector[StressIndices::S00] = stress[0][0];
			stress_vector[StressIndices::S11] = stress[1][1];
			stress_vector[StressIndices::S22] = stress[2][2];
			stress_vector[StressIndices::S01] = stress[0][1];
			stress_vector[StressIndices::S12] = stress[1][2];
			stress_vector[StressIndices::S20] = stress[2][0];
			#elif DIMENSION==2
			stress_vector[StressIndices::S00] = stress[0][0];
			stress_vector[StressIndices::S11] = stress[1][1];
			stress_vector[StressIndices::S01] = stress[0][1];
			#endif

			// evaluate yieldsurface
			vl_F[lfsi.localIndex(0)] = param.YieldSurface(center_global,stress_vector);

			i_view.write( vl_F );
			i_view.commit();
			i_view.unbind();

		}//END: iterate over elements

	}//END: evaluate_yieldsurface

};

#endif /* LANDSLIDELANDSCAPE_DECOUPLED_NEW_OPERATORS_POSTPROCESS_HH_ */
