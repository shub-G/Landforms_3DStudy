/*
 * lop_elasticity.hh
 *
 *  Created on: May 18, 2022
 *      Author: sgupta
 */

#ifndef LANDSLIDELANDSCAPE_DECOUPLED_NEW_OPERATORS_LOP_ELASTICITY_HH_
#define LANDSLIDELANDSCAPE_DECOUPLED_NEW_OPERATORS_LOP_ELASTICITY_HH_

template < class GV,
		   class GFS_F, class F,
		   class Params,
		   class BCT, class BCV>
class LocalOperatorELASTICITY :
	//public Dune::PDELab::NumericalJacobianApplyVolume	<LocalOperatorELASTICITY<GV,GFS_F,F,Params,BCT,BCV> >,
	//public Dune::PDELab::NumericalJacobianVolume		<LocalOperatorELASTICITY<GV,GFS_F,F,Params,BCT,BCV>>,
	//public Dune::PDELab::NumericalJacobianApplyBoundary	<LocalOperatorELASTICITY<GV,GFS_F,F,Params,BCT,BCV> >,
	//public Dune::PDELab::NumericalJacobianBoundary	<LocalOperatorELASTICITY<GV,GFS_F,F,Params,BCT,BCV> >,
	public Dune::PDELab::FullVolumePattern,
	public Dune::PDELab::LocalOperatorDefaultFlags,
	public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
private:
	    const GV& gv;
	    GFS_F gfs_f;
	    F *fvars_new;
	    F *fvars_ini;
		const Params& params;
		const BCT& bct;
		const BCV& bcv;
		double *time;
		double *dt;
		unsigned int intorder ;

public:
  		// pattern assembly flags
  		enum { doPatternVolume 	 = true  };

  		// residual assembly flags
  		enum { doAlphaVolume  	= true 	};
  		enum { doLambdaBoundary = true 	};

  		typedef typename GV::IndexSet IndexSet;

  		typedef Dune::PDELab::LocalFunctionSpace<GFS_F> LFS_F;
  		typedef Dune::PDELab::LFSIndexCache<LFS_F> LFSCache_F;
  		typedef typename F::template LocalView<LFSCache_F> VectorView_F;

  		// constructor stores parameters
  		LocalOperatorELASTICITY(const GV& gv_,
  								GFS_F gfs_f_,
								F *fvars_new_,
								F *fvars_ini_,
								const Params& params_,
								const BCT& bct_,
								const BCV& bcv_,
								double *time_,
								double *dt_,
								unsigned int intorder_ = 4  )
  		:  gv(gv_),
  		   gfs_f(gfs_f_),
		   fvars_new(fvars_new_),
		   fvars_ini(fvars_ini_),
		   params( params_ ),
		   bct( bct_ ),
  		   bcv( bcv_ ),
		   time(time_),
		   dt(dt_),
		   intorder( intorder_ )
  		{}


        template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
        void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, M & mat) const
        {
        	// Define types
        	using LFSU_SUB = typename LFSU::template Child<0>::Type;
        	using RF = typename M::value_type;
        	using JacobianType = typename LFSU_SUB::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType;
        	using size_type = typename LFSU_SUB::Traits::SizeType;

        	// dimensions
        	const int dim = EG::Entity::dimension;

        	// Reference to cell
        	const auto& cell = eg.entity();

        	// get geometry
        	auto geo = eg.geometry();

        	// Transformation
        	typename EG::Geometry::JacobianInverseTransposed jac;

			// loop over quadrature points
			for (const auto& qp : quadratureRule(geo,intorder)){

				auto ip_global = geo.global(qp.position());

				// transform gradient to real element
				jac = geo.jacobianInverseTransposed(qp.position());

				// material parameters
				auto LP = params.LameParameters(ip_global);

				// geometric weight
				auto factor = qp.weight() * geo.integrationElement(qp.position());

				for(int d=0; d<dim; ++d){
					// evaluate gradient of shape functions
					std::vector<JacobianType> jsu(lfsu.child(d).size());
					std::vector<JacobianType> jsv(lfsv.child(d).size());
					lfsu.child(d).finiteElement().localBasis().evaluateJacobian(qp.position(),jsu);
					lfsv.child(d).finiteElement().localBasis().evaluateJacobian(qp.position(),jsv);

					std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsu.child(d).size());
					for (size_type i=0; i<lfsu.child(d).size(); i++)
						jac.umv(jsu[i][0],gradphi[i]);

					std::vector<Dune::FieldVector<RF,dim> > gradpsi(lfsv.child(d).size());
					for (size_type i=0; i<lfsv.child(d).size(); i++)
						jac.umv(jsv[i][0],gradpsi[i]);

					for (size_type i=0; i<lfsu.child(d).size(); i++){

						for (int k=0; k<dim; k++){

							for (size_type j=0; j<lfsv.child(k).size(); j++){
								// integrate 2.*\mu grad_u:grad_v + \lambda div_u I:grad_v

								mat.accumulate(lfsv.child(k),j,lfsu.child(k),i,
								2.*LP[1] * gradphi[i][d] * gradpsi[j][d]
								*factor);

								mat.accumulate(lfsv.child(k),j,lfsu.child(d),i,
								LP[0] * gradphi[i][d] * gradpsi[j][k]
								*factor);

							}
						}
					}// integrate -> mat.accumulate
				}//END: loop over each dim
			}//END: loop over uad-points
        }//END: jacobian_volume
       

		// volume integral depending on test and ansatz functions
		template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
		void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
		{
			typedef typename LFSU::template Child<0>::Type::Traits::FiniteElementType::
			  Traits::LocalBasisType::Traits::DomainFieldType DF;
			typedef typename LFSU::template Child<0>::Type::Traits::FiniteElementType::
			  Traits::LocalBasisType::Traits::RangeFieldType RF;
			typedef typename LFSU::template Child<0>::Type::Traits::FiniteElementType::
			  Traits::LocalBasisType::Traits::JacobianType JacobianType;
			typedef typename LFSU::template Child<0>::Type::Traits::FiniteElementType::
			  Traits::LocalBasisType::Traits::RangeType RangeType;
			typedef typename LFSU::Traits::SizeType SizeType;

			// dimensions
			const int dim = EG::Entity::dimension;
			const int order = std::max(	lfsu.child(0).finiteElement().localBasis().order(),
										lfsv.child(0).finiteElement().localBasis().order());

			// Get cell
			const auto& cell = eg.entity();

			// Get geometry
			auto geo = eg.geometry();

			// Transformation matrix
			typename EG::Geometry::JacobianInverseTransposed jac;

			// evaluate diffusion tensor at cell center, assume it is constant over elements
			auto ref_el = referenceElement(geo);
			auto localcenter = ref_el.position(0,0);

		    // Flow variables: pressure and porosity
			LFS_F lfs_f(gfs_f);
			LFSCache_F lfscache_f(lfs_f);
			VectorView_F fvars_ini_view((*fvars_ini));
			VectorView_F fvars_new_view((*fvars_new));
			lfs_f.bind(cell);
			lfscache_f.update();
			fvars_ini_view.bind(lfscache_f);
		    std::vector<RF> fl_ini(lfs_f.size());
		    fvars_ini_view.read( fl_ini );
		    fvars_new_view.bind(lfscache_f);
		    std::vector<RF> fl_new(lfs_f.size());
		    fvars_new_view.read( fl_new );

		    double t_new = (*time)+(*dt);

			// loop over quadrature points
			for (const auto& ip : quadratureRule(geo,intorder))
			{
				auto ip_global = geo.global(ip.position());
				auto ip_local = geo.local(ip_global);

		        // transform gradients from reference element to real element
		        jac = eg.geometry().jacobianInverseTransposed(ip.position());

		        // integration element
		        RF factor = ip.weight()*eg.geometry().integrationElement(ip.position());

		        //------------------------------------------------------------------------
				// new and ini pressure
				std::vector<RangeType> phi_pw(lfs_f.child(FlowVariables::pw).size(),0.);
				lfs_f.child(FlowVariables::pw).finiteElement().localBasis().evaluateFunction(ip_local,phi_pw);
				RF pw_new = 0., pw_ini=0.;
				for (SizeType i=0;  i<lfs_f.child(FlowVariables::pw).size(); i++){
					pw_new += fl_new[lfs_f.child(FlowVariables::pw).localIndex(i)] * phi_pw[i];
					pw_ini += fl_ini[lfs_f.child(FlowVariables::pw).localIndex(i)] * phi_pw[i];
				}
				// new and ini porosity
				std::vector<RangeType> phi_por(lfs_f.child(FlowVariables::porosity).size(),0.);
				lfs_f.child(FlowVariables::porosity).finiteElement().localBasis().evaluateFunction(ip_local,phi_por);
				RF por_new = 0., por_ini=0.;
				for (SizeType i=0;  i<lfs_f.child(FlowVariables::porosity).size(); i++){
					por_new += fl_new[lfs_f.child(FlowVariables::porosity).localIndex(i)] * phi_por[i];
					por_ini += fl_ini[lfs_f.child(FlowVariables::porosity).localIndex(i)] * phi_por[i];
				}
				//std::cout<< por_ini << '\t' << por_new << '\t' << pw_ini << '\t' << pw_new << std::endl;

				// grad pw_new and grad pw_ini
		        std::vector<JacobianType> js_pw( lfs_f.child(FlowVariables::pw).size(),0. );
		        lfs_f.child(FlowVariables::pw).finiteElement().localBasis().evaluateJacobian(ip_local,js_pw);

		        std::vector<Dune::FieldVector<RF,dim>> gradphi_pw(lfs_f.child(FlowVariables::pw).size());
				for (SizeType i=0; i < lfs_f.child(FlowVariables::pw).size(); i++){
					jac.mv(js_pw[i][0],gradphi_pw[i]);
				}

		        Dune::FieldVector<RF,dim> grad_pw_new(0.0);
		        Dune::FieldVector<RF,dim> grad_pw_ini(0.0);
		        for (SizeType i=0; i<lfs_f.child(FlowVariables::pw).size(); i++){
					grad_pw_new.axpy( fl_new[lfs_f.child(FlowVariables::pw).localIndex(i)],gradphi_pw[i] );
					grad_pw_ini.axpy( fl_ini[lfs_f.child(FlowVariables::pw).localIndex(i)],gradphi_pw[i] );
		        }
		        //------------------------------------------------------------------------

		        // properties
		        auto LP = params.LameParameters(ip_global);
		        auto g = params.GravityVector();
		        auto rhos = params.SedimentDensity();
		        auto rhoeff_new = std::max((1.-por_new)*rhos,0.);
		        auto rhoeff_ini = std::max((1.-por_ini)*rhos,0.);
		        auto alpha = 1.; //Biot-Willis constant

				for( SizeType d=0; d<dim; d++ ){

					// evaluate basis functions
			        std::vector<RangeType> phi( lfsu.child(d).size(),0. );
			        std::vector<RangeType> psi( lfsv.child(d).size(),0. );
			        lfsu.child(d).finiteElement().localBasis().evaluateFunction(ip_local,phi);
			        lfsv.child(d).finiteElement().localBasis().evaluateFunction(ip_local,psi);

			        // evaluate gradient of basis functions on reference element
			        std::vector<JacobianType> jsu( lfsu.child(d).size(),0. );
			        std::vector<JacobianType> jsv( lfsv.child(d).size(),0. );
			        lfsu.child(d).finiteElement().localBasis().evaluateJacobian(ip_local,jsu);
			        lfsv.child(d).finiteElement().localBasis().evaluateJacobian(ip_local,jsv);

			        std::vector<Dune::FieldVector<RF,dim>> gradphi(lfsu.child(d).size());
					for (SizeType i=0; i < lfsu.child(d).size(); i++){
						jac.mv(jsu[i][0],gradphi[i]);
					}

			        std::vector<Dune::FieldVector<RF,dim>> gradpsi(lfsv.child(d).size());
					for (SizeType i=0; i < lfsv.child(d).size(); i++){
						jac.mv(jsv[i][0],gradpsi[i]);
					}

			        // compute gradient of dth component of u
			        Dune::FieldVector<RF,dim> gradu( 0. );
			        for (SizeType i=0; i<lfsu.child(d).size(); i++){
			        	gradu.axpy(x(lfsu.child(d),i),gradphi[i]);
			        }

		    		for (SizeType i=0; i<lfsv.child(d).size(); i++){
		                for (int k=0; k<dim; k++){

		                	r.accumulate( lfsv.child(d) , i , 2. * LP[1] * gradu[k] * gradpsi[i][k] * factor );

		                	r.accumulate( lfsv.child(k) , i , LP[0] * gradu[d] * gradpsi[i][k] * factor );

		                }
						r.accumulate( lfsv.child(d) , i , - ( alpha*grad_pw_new[d] - alpha*grad_pw_ini[d] ) * psi[i] * factor );
						r.accumulate( lfsv.child(d) , i , + ( rhoeff_new*g[d] - rhoeff_ini*g[d] ) * psi[i] * factor );
		    		}
				}

			}//End Quadrature Rule

		}//END:alpha_volume

		template<typename IG, typename LFSV, typename R>
		void lambda_boundary (const IG& ig, const LFSV& lfsv, R& r) const
		{
			double t_new = (*time)+(*dt);

			// domain and range field type (assume both components have same RF)
			typedef typename LFSV::template Child<0>::Type::Traits::FiniteElementType::
			  Traits::LocalBasisType::Traits::DomainFieldType DF;
			typedef typename LFSV::template Child<0>::Type::Traits::FiniteElementType::
			  Traits::LocalBasisType::Traits::RangeFieldType RF;
			typedef typename LFSV::template Child<0>::Type::Traits::FiniteElementType::
			  Traits::LocalBasisType::Traits::JacobianType JacobianType;
			typedef typename LFSV::template Child<0>::Type::Traits::FiniteElementType::
			  Traits::LocalBasisType::Traits::RangeType RangeType;
			typedef typename LFSV::Traits::SizeType SizeType;

			// dimensions
			const unsigned int dim = IG::Entity::dimension;

			// References to inside and outside cells
			const auto& cell_inside = ig.inside();

			// Get geometries
			auto geo = ig.geometry();
			auto geo_inside = cell_inside.geometry();

			// Get geometry of intersection in local coordinates of cell_inside and cell_outside
			auto geo_in_inside = ig.geometryInInside();

			// loop over quadrature points
			for (const auto& ip : quadratureRule(geo,intorder))
			{
				// position of quadrature point in local coordinates of elements
				auto local = geo_in_inside.global(ip.position());
				auto normal = ig.unitOuterNormal(ip.position());
				auto global = geo.global(ip.position());

				// evaluate basis functions at integration point
				std::vector< std::vector< RangeType> > psi( dim );
				for( SizeType d=0; d<dim; d++ ){
					psi[d] = std::vector<RangeType> ( lfsv.child(d).size() );
		          	lfsv.child(d).finiteElement().localBasis().evaluateFunction( local,psi[d]);
				}

				// integrate j
				RF factor = ip.weight()*ig.geometry().integrationElement(ip.position());

		        // RESIDUALS
		        RF tmp = 0.0 ;
		        bool BCTisDirichlet = true; //default
		        typename BCV::Traits::RangeType f;

		        for( SizeType d=0; d<dim; d++ ){

					if(d==0)
						BCTisDirichlet = bct.template child<0>().isDirichlet(ig,ip.position() );
					else if(d==1)
						BCTisDirichlet = bct.template child<1>().isDirichlet(ig,ip.position() );
#if DIMENSION==3
					else if(d==2)
						BCTisDirichlet = bct.template child<2>().isDirichlet(ig,ip.position() );
#endif

					if( BCTisDirichlet == false ){
						f = bcv.evaluate( ig , ip.position() , t_new, d ) ;
						for (SizeType i=0; i<lfsv.child(d).size(); ++i){
							tmp = -f * psi[d][i];
							r.accumulate( lfsv.child(d) , i , tmp * factor );
						}
					}
		        }//END: d^th residual

			}//END: quadrature rule

		}//END:lambda_boundary
};

#endif /* LANDSLIDELANDSCAPE_DECOUPLED_NEW_OPERATORS_LOP_ELASTICITY_HH_ */
