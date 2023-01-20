/*
 * lop_flow.hh
 *
 *  Created on: May 18, 2022
 *      Author: sgupta
 */

#ifndef LANDSLIDELANDSCAPE_DECOUPLED_NEW_OPERATORS_LOP_FLOW_HH_
#define LANDSLIDELANDSCAPE_DECOUPLED_NEW_OPERATORS_LOP_FLOW_HH_

template <class GV,
		  class GFS_F, class F,
		  class GFS_YS, class YS,
		  class Params,
		  class BCT, class BCV>
class LocalOperatorFLOW :
  public Dune::PDELab::NumericalJacobianApplyVolume		< LocalOperatorFLOW<GV,GFS_F,F,GFS_YS,YS,Params,BCT,BCV> >,
  public Dune::PDELab::NumericalJacobianVolume			< LocalOperatorFLOW<GV,GFS_F,F,GFS_YS,YS,Params,BCT,BCV> >,
  public Dune::PDELab::NumericalJacobianApplySkeleton	< LocalOperatorFLOW<GV,GFS_F,F,GFS_YS,YS,Params,BCT,BCV> >,
  public Dune::PDELab::NumericalJacobianSkeleton		< LocalOperatorFLOW<GV,GFS_F,F,GFS_YS,YS,Params,BCT,BCV> >,
  public Dune::PDELab::NumericalJacobianApplyBoundary	< LocalOperatorFLOW<GV,GFS_F,F,GFS_YS,YS,Params,BCT,BCV> >,
  public Dune::PDELab::NumericalJacobianBoundary		< LocalOperatorFLOW<GV,GFS_F,F,GFS_YS,YS,Params,BCT,BCV> >,
  public Dune::PDELab::FullSkeletonPattern,                     // matrix entries skeleton
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags,
  public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{

private:
    const GV& gv;
    GFS_F gfs_f;
    F *f_old;
    GFS_YS gfs_ys;
    YS *ys;
	const Params& param;
	const BCT& bct;
	const BCV& bcv;
	double *time;
	double *dt;
	unsigned int intorder ;

public:
	// pattern assembly flags
	  enum { doPatternVolume	= true };
	  enum { doPatternSkeleton	= false };

	  // residual assembly flags
	  enum { doAlphaVolume  	= true };
	  enum { doAlphaSkeleton	= false };
	  enum { doLambdaBoundary	= true };

	  typedef Dune::PDELab::LocalFunctionSpace<GFS_F> LFS_F;
	  typedef Dune::PDELab::LFSIndexCache<LFS_F> LFSCache_F;
	  typedef typename F::template LocalView<LFSCache_F> VectorView_F;

	  typedef Dune::PDELab::LocalFunctionSpace<GFS_YS> LFS_YS;
	  typedef Dune::PDELab::LFSIndexCache<LFS_YS> LFSCacheYS;
	  typedef typename YS::template LocalView<LFSCacheYS> VectorViewYS;

	  // constructor stores parameters
	  LocalOperatorFLOW( const GV&	gv_,
			  	  	  	 GFS_F gfs_f_,
						 F *f_old_,
					 	 GFS_YS gfs_ys_,
						 YS *ys_,
						 const Params& param_,
						 const BCT& bct_,
						 const BCV& bcv_,
						 double	*time_	,
						 double *dt_	,
						 unsigned int intorder_ = 4)
	  : gv(gv_),
		gfs_f(gfs_f_),
		f_old(f_old_),
		gfs_ys(gfs_ys_),
		ys(ys_),
		param(param_),
		bct( bct_ ),
		bcv( bcv_ ),
		time(time_),
		dt( dt_ ),
		intorder( intorder_ )
	  {}

	  // volume integral depending on test and ansatz functions
	  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
	  {
		  using FEM = typename LFSU::template Child<0>::Type;
		  using LocalBasis = typename FEM::Traits::FiniteElementType::Traits::LocalBasisType;
		  using Cache = Dune::PDELab::LocalBasisCache<LocalBasis>; // a cache for local basis evaluations
		  Cache cache;

		  // subspaces
		  // pressure
		  const auto& lfsv_pw = lfsv.template child<FlowVariables::pw>();
		  const auto& lfsu_pw = lfsu.template child<FlowVariables::pw>();
		  // porosity
		  const auto& lfsv_por = lfsv.template child<FlowVariables::porosity>();
		  const auto& lfsu_por = lfsu.template child<FlowVariables::porosity>();

		  // define types
		  using RF = typename LFSU::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType;
		  using RangeType = typename LFSU::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType;
		  using JacobianType = typename LFSU::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType;
		  using size_type = typename LFSU::template Child<0>::Type::Traits::SizeType;

		  // dimensions
		  const int dim = EG::Entity::dimension;

		  // Reference to cell
		  const auto& cell = eg.entity();

		  // Get geometry
		  auto geo = eg.geometry();

		  // Transformation matrix
		  typename EG::Geometry::JacobianInverseTransposed jac;

		  LFS_F lfs_f(gfs_f);
		  LFSCache_F lfscache_f(lfs_f);
		  VectorView_F f_old_view((*f_old));
		  lfs_f.bind(cell);
		  lfscache_f.update();
		  f_old_view.bind(lfscache_f);
		  std::vector<RF> fl_old(lfs_f.size());
		  f_old_view.read( fl_old );

		  LFS_YS lfs_ys(gfs_ys);
		  LFSCacheYS lfscache_ys(lfs_ys);
		  VectorViewYS ys_view((*ys));
		  lfs_ys.bind(cell);
		  lfscache_ys.update();
		  ys_view.bind(lfscache_ys);
		  std::vector<RF> ysl(lfs_ys.size());
		  ys_view.read( ysl );

		  // Initialize vectors outside for loop
		  std::vector<Dune::FieldVector<RF,dim>> gradphi_pw(lfsu_pw.size());
		  std::vector<Dune::FieldVector<RF,dim>> gradpsi_pw(lfsv_pw.size());
		  Dune::FieldVector<RF,dim> grad_pw(0.0);
		  Dune::FieldVector<RF,dim> Kgrad_pw(0.0);
		  std::vector<Dune::FieldVector<RF,dim>> gradphi_por(lfsu_por.size());
		  std::vector<Dune::FieldVector<RF,dim>> gradpsi_por(lfsv_por.size());
  		  Dune::FieldVector<RF,dim> grad_por(0.0);
		  Dune::FieldVector<RF,dim> Cgrad_por(0.0);

		  // loop over quadrature points
		  for (const auto &ip : quadratureRule(geo,intorder))
		  {
			  double t_new = (*time)+(*dt);

			  // position of quadrature point in local coordinates of elements
			  auto ip_global = geo.global(ip.position());
			  auto ip_local = geo.local(ip_global);

			  // transform gradients from reference element to real element
			  jac = eg.geometry().jacobianInverseTransposed(ip.position());

			  // integration factor
			  RF factor = ip.weight()*eg.geometry().integrationElement(ip.position());

			  // evaluate basis functions
			  auto &phi_pw  = cache.evaluateFunction(ip.position(), lfsu_pw.finiteElement().localBasis() );
			  auto &psi_pw  = cache.evaluateFunction(ip.position(), lfsv_pw.finiteElement().localBasis() );
			  auto &phi_por = cache.evaluateFunction(ip.position(), lfsu_por.finiteElement().localBasis());
			  auto &psi_por = cache.evaluateFunction(ip.position(), lfsv_por.finiteElement().localBasis());

			  // evaluate pressure
			  RF pw = 0.0;
			  for (size_type i = 0; i < lfsu_pw.size(); i++)
				  pw += x(lfsu_pw, i) * phi_pw[i];
			  // evaluate porosity
			  RF porosity = 0.0;
			  for (size_type i = 0; i < lfsu_por.size(); i++)
				  porosity += x(lfsu_por, i) * phi_por[i];

			  // old pressure
			  std::vector<RangeType> phi_pw_old(lfs_f.child(FlowVariables::pw).size(),0.);
			  lfs_f.child(FlowVariables::pw).finiteElement().localBasis().evaluateFunction(ip.position(),phi_pw_old);
			  RF pw_old = 0.;
			  for (size_type i=0;  i<lfs_f.child(FlowVariables::pw).size(); i++){
				  pw_old += fl_old[lfs_f.child(FlowVariables::pw).localIndex(i)] * phi_pw_old[i];
			  }
			  // old porosity
			  std::vector<RangeType> phi_por_old(lfs_f.child(FlowVariables::porosity).size(),0.);
			  lfs_f.child(FlowVariables::porosity).finiteElement().localBasis().evaluateFunction(ip.position(),phi_por_old);
			  RF por_old = 0.;
			  for (size_type i=0;  i<lfs_f.child(FlowVariables::porosity).size(); i++){
				  por_old += fl_old[lfs_f.child(FlowVariables::porosity).localIndex(i)] * phi_por_old[i];
			  }

			  // time derivatives
			  auto dt_pw  = (pw-pw_old)/(*dt);
			  auto dt_por = (porosity-por_old)/(*dt);

			  // yield-surface
			  std::vector<RangeType> phi_ys(lfs_ys.size(),0.);
			  	  lfs_ys.finiteElement().localBasis().evaluateFunction(ip.position(),phi_ys);
			  RF yieldfunction=0.;
			  for (size_type i = 0; i < lfs_ys.size(); i++){
				  yieldfunction += ysl[lfs_ys.localIndex(i)] * phi_ys[i];
			  }

			  // PROPERTIES
			  // phase densities [kg/m^3]
			  auto rhow = param.WaterDensity();
			  auto rhos = param.SedimentDensity();
			  auto rhoeff = std::max((1.-porosity)*rhos,0.);
			  //Hydraulic mobility (K/mu) [m^2/Pa.s]
			  auto K = param.PermeabilityTensor(ip_global,porosity);
			  K *= 1./param.WaterViscosity();
			  //source term [1/s]
			  auto qw = param.Sources(ip_global,t_new,pw);
			  // phase compressibilities [1/Pa]
			  auto Cs = param.SedimentCompressibility(ip_global);
			  auto Cw = param.WaterCompressibility();
			  auto Ceff = (1.-porosity)*Cs + porosity*Cw;
			  //Sedimentation rate [m/s]
			  auto vs = param.SedimentationRate(ip_global,t_new);
			  // landslide source
			  auto r_qs = param.LandslideRate(ip_global,yieldfunction,porosity);
			  //gravity
			  auto g = param.GravityVector();
			  // sediment diffusion/creep coefficient
			  auto C = param.CreepTensor(ip_global,porosity);

			  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			  // PRESSURE EQUATION

			  // evaluate gradient of pressure basis functions
			  auto &jsu_pw = cache.evaluateJacobian(ip.position(), lfsu_pw.finiteElement().localBasis());
			  auto &jsv_pw = cache.evaluateJacobian(ip.position(), lfsv_pw.finiteElement().localBasis());
			  for (size_type i = 0; i < lfsu_pw.size(); i++)
				  jac.mv(jsu_pw[i][0], gradphi_pw[i]);
			  for (size_type i = 0; i < lfsv_pw.size(); i++)
				  jac.mv(jsv_pw[i][0], gradpsi_pw[i]);

			  // evaluate gradients
			  grad_pw = 0.0;
			  for (size_type i = 0; i < lfsu_pw.size(); i++)
				  grad_pw.axpy(x(lfsu_pw, i), gradphi_pw[i]);
			  auto Fw = g;
			  Fw *= rhow;
			  Fw += grad_pw;
			  K.mv(Fw, Kgrad_pw);
			  // magnitude of vw
			  double vw_mag_sq = 0.;
			  for(int i=0; i<dim; i++){
				  vw_mag_sq += ( (-1.)*Kgrad_pw[i] )*( (-1.)*Kgrad_pw[i] );
			  }
			  auto vw_mag = sqrt(vw_mag_sq);
			  // erosion rate
			  auto r_er = param.ErosionRate(ip_global,vw_mag,yieldfunction,porosity);

			  // integrate
			  for (size_type i = 0; i < lfsv_pw.size(); i++) {
				  auto term = Ceff * dt_pw * psi_pw[i]
							+ (Kgrad_pw * gradpsi_pw[i])
							- ((1.-porosity)*Cs*rhos + porosity*Cw*rhow)*(vs*g) * psi_pw[i]
							//- vs * gradpsi_pw[i]
							- qw * psi_pw[i]
							- r_qs * psi_pw[i]
							- r_er * psi_pw[i]
							;
//				  std::cout<< "volume: " << term << std::endl;
				  r.accumulate(lfsv_pw, i, term*factor);
			  }

			  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			  // POROSITY EQUATION

			  // evaluate gradient of porosity basis functions
			  auto &jsu_por = cache.evaluateJacobian(ip.position(), lfsu_por.finiteElement().localBasis());
			  auto &jsv_por = cache.evaluateJacobian(ip.position(), lfsv_por.finiteElement().localBasis());
			  for (size_type i = 0; i < lfsu_por.size(); i++)
				  jac.mv(jsu_por[i][0], gradphi_por[i]);
			  for (size_type i = 0; i < lfsv_por.size(); i++)
				  jac.mv(jsv_por[i][0], gradpsi_por[i]);
			  grad_por = 0.0;
			  for (size_type i = 0; i < lfsu_por.size(); i++)
				  grad_por.axpy(x(lfsu_por, i), gradphi_por[i]);
			  C.mv(grad_por, Cgrad_por);
				  
			  // integrate
			  for (size_type i = 0; i < lfsv_por.size(); i++) {
				  auto term = (-1.) * dt_por * psi_por[i]
							+ (1.-porosity) * Cs * dt_pw * psi_por[i]
							- (Cgrad_por * gradpsi_por[i])
							- ((1.-porosity)*Cs*rhos)*(vs*g) * psi_por[i]
							//- (1.-porosity) * (vs * gradpsi_por[i])
							- r_qs * psi_por[i]
							- r_er * psi_por[i]
						    ;
				  r.accumulate(lfsv_por, i, term*factor);
			  }

			  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		  }//END::loop_over_quad_points

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

				// integrate j
				RF factor = ip.weight()*ig.geometry().integrationElement(ip.position());

		        RF tmp = 0.0 ;
		        typename BCV::Traits::RangeType f;
		        bool BCTisDirichlet = true; //default

		        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		        // PRESSURE EQN
				std::vector< RangeType> psi_pw(lfsv.child(FlowVariables::pw).size(),0.0);
				lfsv.child(FlowVariables::pw).finiteElement().localBasis().evaluateFunction( local,psi_pw);

		        BCTisDirichlet = bct.template child<FlowVariables::pw>().isDirichlet(ig,ip.position() );
				if( BCTisDirichlet == false ){
					f = bcv.evaluate( ig , ip.position() , t_new, FlowVariables::pw ) ;
					for (SizeType i=0; i<lfsv.child(FlowVariables::pw).size(); ++i){
						tmp = -f * psi_pw[i];
						r.accumulate( lfsv.child(FlowVariables::pw) , i , tmp * factor );
					}
				}//END::pw residual

				//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		        // POROSITY EQN
		        	tmp = 0.;
				std::vector< RangeType> psi_por(lfsv.child(FlowVariables::porosity).size(),0.0);
				lfsv.child(FlowVariables::porosity).finiteElement().localBasis().evaluateFunction( local,psi_por);

		        BCTisDirichlet = bct.template child<FlowVariables::porosity>().isDirichlet(ig,ip.position() );
				if( BCTisDirichlet == false ){
					f = bcv.evaluate( ig , ip.position() , t_new, FlowVariables::porosity ) ;
					for (SizeType i=0; i<lfsv.child(FlowVariables::porosity).size(); ++i){
						tmp = -f * psi_por[i];
						r.accumulate( lfsv.child(FlowVariables::porosity) , i , tmp * factor );
					}
				}//END::porosity residual

				//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

			}//END: quadrature rule

		}//END:lambda_boundary


};

#endif /* LANDSLIDELANDSCAPE_DECOUPLED_NEW_OPERATORS_LOP_FLOW_HH_ */
