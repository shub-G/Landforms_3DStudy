/*
 * parameters.hh
 *
 *  Created on: June 8, 2022
 *      Author: sgupta
 */

#ifndef LANDSLIDELANDSCAPE_DECOUPLED_NEW_PROBLEM_CONTINENTALMARGIN02_PARAMETERS_HH_
#define LANDSLIDELANDSCAPE_DECOUPLED_NEW_PROBLEM_CONTINENTALMARGIN02_PARAMETERS_HH_

struct StressIndices
{
#if DIMENSION==2
  enum Type { num=3, S00=0, S11=1, S01=2 };
#else
  enum Type { num=6, S00=0, S11=1, S22=2, S01=3, S12=4, S20=5 };
#endif
};

struct FlowVariables {
	enum Type { num = 2, pw=0, porosity=1 };
};

template<typename GV,typename PTree>
class Parameters
{
private:
	  const GV& gv;
	  const PTree& ptree;
	  const static int dim = GV::dimension;
	  double eps = 1.0e-6;
	  using FieldVector = typename Dune::FieldVector<double,dim>;
	  using FieldMatrix = typename Dune::FieldMatrix<double,dim,dim>;

	  bool gravity_flag;
	  double gravity_magnitude;
	  double theta;

	  bool isGraded;
	  double dl1;
	  double n;
	  bool isActiveLens;

	  double sed_por0_x0;
	  double sed_por0_xL;
	  double lens_por0;

	  double sed_K0_x0;
	  double sed_K0_xL;
	  double lens_K0;

	  double sed_KF_x0;
	  double sed_KF_xL;
	  double lens_KF;

	  double sed_A0_x0;
	  double sed_A0_xL;
	  double lens_A0;

	  double sed_beta_x0;
	  double sed_beta_xL;
	  double lens_beta;

	  double sed_beta_c_x0;
	  double sed_beta_c_xL;
	  double lens_beta_c;

	  double sed_E_x0, sed_nu_x0, sed_alpha_x0, sed_c_x0;
	  double sed_E_xL, sed_nu_xL, sed_alpha_xL, sed_c_xL;
	  double lens_E, lens_nu, lens_alpha, lens_c;

	  double sed_R_x0, sed_R_xL;
	  double lens_R;

	  double sed_er_rate_x0, sed_F_cr_x0, sed_vw_cr_x0;
	  double sed_er_rate_xL, sed_F_cr_xL, sed_vw_cr_xL;
	  double lens_er_rate, lens_F_cr, lens_vw_cr;
	  double er_exponent;
	  double beta_er;
	  
  	  double sed_C_x0;
	  double sed_C_xL;
	  double lens_C;

	  double L1, L2, L3, LT;
	  double H1, H2, H3, HT;
	  double DL, dL, dc;
	  double m1, m2;
	  double HB, LB;
	  double theta1, theta2, theta3;
#if DIMENSION==3
	  double Z;
	  double ZL;
#endif

	  double HSL0, Tm, T, HSL0_offset;

	  double Cw, lens_Cs, sed_Cs_x0, sed_Cs_xL;

	  double burial_rate;

public:
	//! construct from grid view
	Parameters (const GV& gv_, const PTree& ptree_)
	: gv( gv_ ) ,
	  ptree(ptree_)
	{
		gravity_flag = ptree.get("gravity.flag",(bool)true);
		gravity_magnitude = 9.81; /*m/s^2*/
		theta = ptree.get("gravity.theta",(double)0.0);//degrees

		isGraded = ptree.get("sediment.is_graded",(bool)false);
		dl1 = ptree.get("sediment.shelf_break.offset",(double)0.0);//km
		dl1 *= 1000.; //m
		n = ptree.get("sediment.shelf_break.n",(double)1.0);//-

		sed_por0_x0 = ptree.get("sediment.onshore.porosity",(double)0.1);/*-*/
		sed_K0_x0 	= ptree.get("sediment.onshore.permeability",(double)1.e-16);/*m^2*/
		sed_KF_x0 	= ptree.get("sediment.onshore.Kzz_over_Kxx",(double)1.0);/*-*/
		sed_A0_x0 	= ptree.get("sediment.onshore.porosity_scaling_factor",(double)0.0);
		sed_beta_x0 = ptree.get("sediment.onshore.decay_index",(double)3.2e-3); //1/m
		sed_beta_c_x0 = ptree.get("sediment.onshore.decay_index_cohesion",(double)2.0e-3); //1/m
		sed_Cs_x0 	= ptree.get("sediment.onshore.compressibility",(double)1.e-8);//1/Pa
		sed_E_x0 	= ptree.get("sediment.onshore.youngs_modulus",(double)1.e10);//Pa
		sed_nu_x0 	= ptree.get("sediment.onshore.poissons_ratio",(double)0.15);//-
		sed_alpha_x0= ptree.get("sediment.onshore.friction_coefficient",(double)0.5);//-
		sed_c_x0 	= ptree.get("sediment.onshore.cohesive_strength",(double)1.e6);//Pa
		sed_R_x0 	= ptree.get("sediment.onshore.landslide_rate",(double)0.0);
		sed_er_rate_x0 	= ptree.get("sediment.onshore.erosion.rate",(double)0.0);//1/m
		sed_F_cr_x0 	= ptree.get("sediment.onshore.erosion.critical_F",(double)1.e6);//Pa
		sed_vw_cr_x0 	= ptree.get("sediment.onshore.erosion.critical_vw",(double)0.0);//m/s
		sed_C_x0 	= ptree.get("sediment.onshore.creep_coefficient",(double)0.0);//TODO

		sed_por0_xL = ptree.get("sediment.offshore.porosity",(double)0.1);/*-*/
		sed_K0_xL	= ptree.get("sediment.offshore.permeability",(double)1.e-16);/*m^2*/
		sed_KF_xL 	= ptree.get("sediment.offshore.Kzz_over_Kxx",(double)1.0);/*-*/
		sed_A0_xL 	= ptree.get("sediment.offshore.porosity_scaling_factor",(double)0.0);
		sed_beta_xL = ptree.get("sediment.offshore.decay_index",(double)3.2e-3); //1/m
		sed_beta_c_xL = ptree.get("sediment.offshore.decay_index_cohesion",(double)1.0e-3); //1/m
		sed_Cs_xL 	= ptree.get("sediment.offshore.compressibility",(double)1.e-8);//1/Pa
		sed_E_xL 	= ptree.get("sediment.offshore.youngs_modulus",(double)1.e10);//Pa
		sed_nu_xL 	= ptree.get("sediment.offshore.poissons_ratio",(double)0.15);//-
		sed_alpha_xL= ptree.get("sediment.offshore.friction_coefficient",(double)0.5);//-
		sed_c_xL 	= ptree.get("sediment.offshore.cohesive_strength",(double)1.e6);//Pa
		sed_R_xL 	= ptree.get("sediment.offshore.landslide_rate",(double)0.0);
		sed_er_rate_xL 	= ptree.get("sediment.offshore.erosion.rate",(double)0.0);//1/m
		sed_F_cr_xL 	= ptree.get("sediment.offshore.erosion.critical_F",(double)1.e6);//Pa
		sed_vw_cr_xL 	= ptree.get("sediment.offshore.erosion.critical_vw",(double)0.0);//m/s
		sed_C_xL 	= ptree.get("sediment.offshore.creep_coefficient",(double)0.0);//TODO

		isActiveLens = ptree.get("lens.is_active",(bool)false);
		lens_por0 = ptree.get("lens.porosity",(double)0.5);/*-*/
		lens_K0 = ptree.get("lens.permeability",(double)1.e-13);/*m^2*/
		lens_KF = ptree.get("lens.Kzz_over_Kxx",(double)1.0);/*-*/
		lens_A0 = ptree.get("lens.porosity_scaling_factor",(double)0.0);
		lens_beta = ptree.get("lens.decay_index",(double)1.5e-3); //1/m
		lens_beta_c = ptree.get("lens.decay_index_cohesion",(double)2.0e-3); //1/m
		lens_Cs = ptree.get("lens.compressibility",(double)1.e-6);//1/Pa
		lens_E = ptree.get("lens.youngs_modulus",(double)1.e10);//Pa
		lens_nu = ptree.get("lens.poissons_ratio",(double)0.15);//-
		lens_alpha = ptree.get("lens.friction_coefficient",(double)0.5);//-
		lens_c = ptree.get("lens.cohesive_strength",(double)1.e6);//Pa
		lens_R = ptree.get("lens.landslide_rate",(double)0.0);
		lens_er_rate 	= ptree.get("lens.erosion.rate",(double)0.0);//1/m
		lens_F_cr 	= ptree.get("lens.erosion.critical_F",(double)1.e6);//Pa
		lens_vw_cr 	= ptree.get("lens.erosion.critical_vw",(double)0.0);//m/s
		lens_C 	= ptree.get("lens.creep_coefficient",(double)0.0);//TODO

		er_exponent = ptree.get("sediment.erosion_law_exponent",(double)1.0);
		beta_er = ptree.get("sediment.erosion_rate_decay",(double)0.0);

		Cw = ptree.get("water.compressibility",(double)1.e-10);//1/Pa

		double scaleX=1.0;
		double scaleY=1.0;
		L1= ptree.get("domain.L1",(double)50.0);//km
		L1 *= 1000.0/scaleX;
		double hT1 = ptree.get("domain.dH1",(double)0.120);//km
		double LT1 = ptree.get("domain.dL1",(double)150.0);//km
		H1=L1*((hT1*1000.0*scaleY)/(LT1*1000*scaleX));
		L2=ptree.get("domain.L2",(double)30.0)*1000/scaleX;
		H2=ptree.get("domain.H2",(double)2.0)*1000/scaleY;
		theta1=std::atan(H1/L1);
		theta2=std::atan(H2/L2);
		HT=H1+H2;
		LT=L1+L2;
		DL=ptree.get("domain.DL",(double)0.300)*1000/scaleY;
		dL=ptree.get("domain.dL",(double)0.025)*1000/scaleY;
		m1=H1/L1;
		m2=H2/L2;
		dc=ptree.get("domain.dc",(double)2.0)*1000/scaleX;
		double DB=ptree.get("domain.H2",(double)0.500)*1000*scaleY;
		HB=HT-(DL+dL+DB);
		LB=(115./100.)*L1;
#if DIMENSION==3
		Z=ptree.get("domain.DZ",(double)2.0)*1000;
		ZL=ptree.get("domain.DZL",(double)0.5)*1000;
#endif

		HSL0 = ptree.get("sealevel.HSL0",(double)120.0); //m
		Tm	 = ptree.get("sealevel.Tm",(double)100000.0);//years
		Tm  *= 86400.0*365.0; //s
		T	 = ptree.get("sealevel.T",(double)120000.0);//years
		T   *= 86400.0*365.0; //s
		HSL0_offset = hT1*1.e3-30.0-HSL0;

		burial_rate = ptree.get("burial.rate",(double)0.);//m/years
		burial_rate *= 1./(86400.0*365.0);

		std::cout<< sed_KF_x0 << '\t' << sed_KF_xL << std::endl;
	}

	/*****************************************************************
	 * SEAFLOOR
	 */
	double InitialSeafloorHeight(FieldVector globalpos /*m*/) const {
		double h=HT;
		if( globalpos[0]<=L1 ) h=HT-m1*globalpos[0];
		else if( globalpos[0]>L1 and globalpos[0]<=L1+L2 ) h=(HT-H1)-m2*(globalpos[0]-L1);
		else h = HT-H1-H2;
		return h;/*m*/
	}

	double SeaLevel(FieldVector globalpos /*m*/, double time )const{
		double h= HSL0;
		if( time < Tm ) h = HSL0 * (1. - time/Tm);
		else if( time >= Tm and time < T ) h = HSL0 * ( (time-Tm)/(T-Tm) );
		else h = HSL0;

		double hSL = HT-H1+h+HSL0_offset;
		return hSL;
	}

	/*****************************************************************
	 * BOUNDARY SEGMENTS
	 */
	bool isSeafloorBoundary(FieldVector globalpos /*m*/)const{
		double HSF = InitialSeafloorHeight(globalpos);
		if(globalpos[1] > HSF-eps ) return true;
		else return false;
	}

	bool isBedrockBoundary(FieldVector globalpos /*m*/)const{
		if( globalpos[0]<LB+eps and globalpos[1]<HB+eps ) return true;
		else return false;
	}

	bool isBottomBoundary(FieldVector globalpos/*m*/)const{
		if(globalpos[0]>LB and globalpos[1]<0.+eps) return true;
		else return false;
	}

	bool isUpstreamBoundary(FieldVector globalpos /*m*/)const{
		if( globalpos[0]<0.+eps and globalpos[1]>HB) return true;
		else return false;
	}

	bool isDownstreamBoundary(FieldVector globalpos /*m*/)const{
		if( globalpos[0]>LT-eps) return true;
		else return false;
	}
#if DIMENSION==3
	bool isFrontBoundary(FieldVector globalpos /*m*/)const{
		if( globalpos[2]<0.+eps) return true;
		else return false;
	}

	bool isBackBoundary(FieldVector globalpos /*m*/)const{
		if( globalpos[2]>Z-eps) return true;
		else return false;
	}
#endif
	/*****************************************************************
	 * INTERIOR ZONES
	 */
	bool isSandLens(FieldVector globalpos /*m*/)const{
		if( isActiveLens
			and globalpos[1] 	< (HT-DL)-m1*globalpos[0] + eps
			and globalpos[1]	> (HT-DL-dL)-m1*globalpos[0] - eps
			and globalpos[0] 	< (HT-H1-globalpos[1])/m2 + (L1-dc) + eps
#if DIMENSION==3
			and globalpos[2] > Z/2.-ZL/2.-eps and globalpos[2] < Z/2.+ZL/2.+eps
#endif
		) return true;
		else return false;
	}

	bool isClaySediment(FieldVector globalpos /*m*/)const{
		if( isSandLens(globalpos) ) return false;
		else return true;
	}

	/*****************************************************************
	 * PROPERTIES
	 */
	double WaterDensity() const {
		return 1000.0; //kg/m^3
	}

	double WaterViscosity() const {
		return 0.001;// Pa.s
	}

	double WaterCompressibility() const {
		return Cw;// 1/Pa
	}

	double SedimentDensity() const {
		return 2500.0; //kg/m^3
	}

	double SedimentCompressibility(FieldVector globalpos) const {
		double Zmax = InitialSeafloorHeight(globalpos);
		double Fz_x0 = std::exp(-sed_beta_c_x0*(Zmax-globalpos[1]));
		double Fz_xL = std::exp(-sed_beta_c_xL*(Zmax-globalpos[1]));
		double Fz_lens = std::exp(-lens_beta_c*(Zmax-globalpos[1]));

		double Cs=0.;
		if( isSandLens(globalpos) ) Cs=lens_Cs*Fz_lens;
		else{
			Cs=sed_Cs_xL*Fz_xL;
			// if(isGraded) Cs += (sed_Cs_x0-sed_Cs_xL)*(1.-std::min(globalpos[0],L1+dl1)/(L1+dl1));
			if(isGraded){
				double x0=0., y0=sed_Cs_x0*Fz_x0;
				double x1=L1+dl1, y1=n*sed_Cs_xL*Fz_xL;
				double x2=LT, y2=sed_Cs_xL*Fz_xL;
				double a = ( x0*(y2-y1) + x1*(y0-y2) + x2*(y1-y0) )/(( x0-x2)*(x1-x0)*(x2-x1) );
				double b = (y2-y0)/(x2-x0) - a*(x2+x0);
				double c = y2 - a*x2*x2 - b*x2;
				Cs = a*globalpos[0]*globalpos[0] + b*globalpos[0] + c;
			}
		}
		return Cs;
	}

	FieldVector GravityVector() const {
		FieldVector gravity( 0. );
		double g = 0.;
		if(gravity_flag) g = gravity_magnitude;
		gravity[1] = g;
		gravity[0] = g*std::sin(theta*M_PI/180.0);
		return gravity; //m/s^2
	}

	double InitialPorosity(FieldVector globalpos/*m*/)const{
		double por=0.;
		if( isSandLens(globalpos) ) por=lens_por0;
		else{
			por=sed_por0_xL;
			// if(isGraded) por += (sed_por0_x0-sed_por0_xL)*(1.-std::min(globalpos[0],L1+dl1)/(L1+dl1));
			if(isGraded){
				double x0=0., y0=sed_por0_x0;
				double x1=L1+dl1, y1=n*sed_por0_xL;
				double x2=LT, y2=sed_por0_xL;
				double a = ( x0*(y2-y1) + x1*(y0-y2) + x2*(y1-y0) )/(( x0-x2)*(x1-x0)*(x2-x1) );
				double b = (y2-y0)/(x2-x0) - a*(x2+x0);
				double c = y2 - a*x2*x2 - b*x2;
				por = a*globalpos[0]*globalpos[0] + b*globalpos[0] + c;
			}
		}
		return por;
	}

	FieldMatrix
	PermeabilityTensor(FieldVector globalpos/*m*/, double porosity) const {

		double K0, Fz, Fpor, KF;
		double Zmax = InitialSeafloorHeight(globalpos);
		double porosity0 = InitialPorosity(globalpos);

		double Fz_x0 = std::exp(-sed_beta_x0*(Zmax-globalpos[1]));
		double Fpor_x0=std::exp( sed_A0_x0*(porosity-porosity0)/(1.-porosity0));
		double Fz_xL = std::exp(-sed_beta_xL*(Zmax-globalpos[1]));
		double Fpor_xL=std::exp( sed_A0_xL*(porosity-porosity0)/(1.-porosity0));
		double KF_x0=sed_KF_x0;
		double KF_xL=sed_KF_xL;

		if( isSandLens(globalpos) ){
			Fz = std::exp(-lens_beta*(Zmax-globalpos[1]));
			Fpor=std::exp( lens_A0*(porosity-porosity0)/(1.-porosity0));
			KF=lens_KF;
			K0=lens_K0*Fz*Fpor;
		}else{
			Fz = Fz_xL;
			Fpor=Fpor_xL;
			K0=sed_K0_xL*Fz_xL*Fpor_xL;
			// if(isGraded) {
			// 	K0 		+= (sed_K0_x0*Fz_x0*Fpor_x0-sed_K0_xL*Fz_xL*Fpor_xL)*(1.-std::min(globalpos[0],L1+dl1)/(L1+dl1));
			// }
			if(isGraded){
				double x0=0., y0=sed_K0_x0;
				double x1=L1+dl1, y1=n*sed_K0_xL;
				double x2=LT, y2=sed_K0_xL;
				double a = ( x0*(y2-y1) + x1*(y0-y2) + x2*(y1-y0) )/(( x0-x2)*(x1-x0)*(x2-x1) );
				double b = (y2-y0)/(x2-x0) - a*(x2+x0);
				double c = y2 - a*x2*x2 - b*x2;
				K0 = a*globalpos[0]*globalpos[0] + b*globalpos[0] + c;
				// K0 *= Fz_xL*Fpor_xL;

				y0=Fz_x0;
				y1=n*Fz_xL;
				y2=Fz_xL;
				a = ( x0*(y2-y1) + x1*(y0-y2) + x2*(y1-y0) )/(( x0-x2)*(x1-x0)*(x2-x1) );
				b = (y2-y0)/(x2-x0) - a*(x2+x0);
				c = y2 - a*x2*x2 - b*x2;
				Fz = a*globalpos[0]*globalpos[0] + b*globalpos[0] + c;

				y0=Fpor_x0;
				y1=n*Fpor_xL;
				y2=Fpor_xL;
				a = ( x0*(y2-y1) + x1*(y0-y2) + x2*(y1-y0) )/(( x0-x2)*(x1-x0)*(x2-x1) );
				b = (y2-y0)/(x2-x0) - a*(x2+x0);
				c = y2 - a*x2*x2 - b*x2;
				Fpor = a*globalpos[0]*globalpos[0] + b*globalpos[0] + c;

				y0=KF_x0;
				y1=n*KF_xL;
				y2=KF_xL;
				a = ( x0*(y2-y1) + x1*(y0-y2) + x2*(y1-y0) )/(( x0-x2)*(x1-x0)*(x2-x1) );
				b = (y2-y0)/(x2-x0) - a*(x2+x0);
				c = y2 - a*x2*x2 - b*x2;
				KF = a*globalpos[0]*globalpos[0] + b*globalpos[0] + c;

				K0 *= Fz*Fpor;
			}
		}

		FieldMatrix K(0.);
//		for( int d=0; d<dim; d++) K[d][d] = K0*Fz;
//		K[0][0] = K0;
//		K[1][1] = KF*K0;
//		K[2][2] = K0;
		double theta=0.;
		if( globalpos[0]<=L1) theta=theta1;
		else theta=theta2;
		K[0][0] = K0*std::cos(theta);		K[0][1] = K0*KF*std::sin(theta);
		K[1][0] = K0*(-1.)*std::sin(theta);	K[1][1] = K0*KF*std::cos(theta);
		K[2][2] = K0;

		return K; //m^2
	}
	
	
	FieldMatrix
	CreepTensor(FieldVector globalpos/*m*/, double porosity) const {

		double C0;

		if( isSandLens(globalpos) ){
			C0=lens_C;
		}else{
			C0=sed_C_xL;
			// if(isGraded) {
			// 	K0 		+= (sed_K0_x0*Fz_x0*Fpor_x0-sed_K0_xL*Fz_xL*Fpor_xL)*(1.-std::min(globalpos[0],L1+dl1)/(L1+dl1));
			// }
			if(isGraded){
				double x0=0., y0=sed_C_x0;
				double x1=L1+dl1, y1=n*sed_C_xL;
				double x2=LT, y2=sed_C_xL;
				double a = ( x0*(y2-y1) + x1*(y0-y2) + x2*(y1-y0) )/(( x0-x2)*(x1-x0)*(x2-x1) );
				double b = (y2-y0)/(x2-x0) - a*(x2+x0);
				double c = y2 - a*x2*x2 - b*x2;
				C0 = a*globalpos[0]*globalpos[0] + b*globalpos[0] + c;
			}
		}

		FieldMatrix C(0.);
		for( int d=0; d<dim; d++) C[d][d] = C0*(1-porosity);
		
		C[1][1] *= 10.;

		return C; //TODO
	}
	

	FieldVector SedimentationRate(FieldVector globalpos /*m*/, double time /*s*/) const {
		FieldVector vs( 0. );
		vs[1] = burial_rate;
		vs[0] = 0.;
		return vs; //m/s
	}

	double Sources(FieldVector globalpos /*m*/, double time /*s*/, double P /*Pa*/)const{
		return 0.; /*1/s*/
	}

	double LandslideRate(FieldVector globalpos /*m*/, double F /*Pa*/, double porosity)const{
		double r=0.;
		if( F>1.e5) {
			if( isSandLens(globalpos) ) r = lens_R;
			else{
				r=sed_R_xL;
				// if(isGraded) r += (sed_R_x0-sed_R_xL)*(1.-std::min(globalpos[0],L1+dl1)/(L1+dl1));
				if(isGraded){
				double x0=0., y0=sed_R_x0;
				double x1=L1+dl1, y1=n*sed_R_xL;
				double x2=LT, y2=sed_R_xL;
				double a = ( x0*(y2-y1) + x1*(y0-y2) + x2*(y1-y0) )/(( x0-x2)*(x1-x0)*(x2-x1) );
				double b = (y2-y0)/(x2-x0) - a*(x2+x0);
				double c = y2 - a*x2*x2 - b*x2;
				r = a*globalpos[0]*globalpos[0] + b*globalpos[0] + c;
				}
			}
		}
		return -1. * r * (1.-porosity);
	}

	double ErosionRate (FieldVector globalpos /*m*/, double vw_mag /*m/s*/, double F /*Pa*/, double porosity) const {
		double qer=0., erosion_rate=0., F_cr=0., vw_cr=0.;
		double Zmax = InitialSeafloorHeight(globalpos);
		double Fz = std::exp(-beta_er*(Zmax-globalpos[1]));
		
		if( isSandLens(globalpos) ){
			erosion_rate = lens_er_rate;
			F_cr = lens_F_cr;
			vw_cr = lens_vw_cr;
		}
		else{
			erosion_rate = sed_er_rate_xL;
			F_cr = sed_F_cr_xL;
			vw_cr = sed_vw_cr_xL;
			// if(isGraded) r += (sed_R_x0-sed_R_xL)*(1.-std::min(globalpos[0],L1+dl1)/(L1+dl1));
			if(isGraded){
				double x0=0., y0=sed_er_rate_x0;
				double x1=L1+dl1, y1=n*sed_er_rate_xL;
				double x2=LT, y2=sed_er_rate_xL;
				double a = ( x0*(y2-y1) + x1*(y0-y2) + x2*(y1-y0) )/(( x0-x2)*(x1-x0)*(x2-x1) );
				double b = (y2-y0)/(x2-x0) - a*(x2+x0);
				double c = y2 - a*x2*x2 - b*x2;
				erosion_rate = a*globalpos[0]*globalpos[0] + b*globalpos[0] + c;

				y0=sed_F_cr_x0;
				y1=n*sed_F_cr_xL;
				y2=sed_F_cr_xL;
				a = ( x0*(y2-y1) + x1*(y0-y2) + x2*(y1-y0) )/(( x0-x2)*(x1-x0)*(x2-x1) );
				b = (y2-y0)/(x2-x0) - a*(x2+x0);
				c = y2 - a*x2*x2 - b*x2;
				F_cr = a*globalpos[0]*globalpos[0] + b*globalpos[0] + c;

				y0=sed_vw_cr_x0;
				y1=n*sed_vw_cr_xL;
				y2=sed_vw_cr_xL;
				a = ( x0*(y2-y1) + x1*(y0-y2) + x2*(y1-y0) )/(( x0-x2)*(x1-x0)*(x2-x1) );
				b = (y2-y0)/(x2-x0) - a*(x2+x0);
				c = y2 - a*x2*x2 - b*x2;
				vw_cr = a*globalpos[0]*globalpos[0] + b*globalpos[0] + c;
			}
		}

		if(F>F_cr and vw_mag>vw_cr and globalpos[0]>10.*1000.){
			double F0 = 1.e5;
			//qer = (-1.) * (F/(F_cr+1.e-15)) * erosion_rate /*1/m*/
			//			* std::pow(std::max((1.-porosity),0.),er_exponent)
			//			* (vw_mag-vw_cr)/*m/s*/;
			qer = (-1.) * erosion_rate * Fz /*1/m*/
						* std::pow(std::max((1.-porosity),0.),er_exponent)
						* (vw_mag-vw_cr)/*m/s*/;
		}
		return qer; //1/s
	}

	double YoungsModulus (FieldVector globalpos) const {
		double E=0.;
		if( isSandLens(globalpos) ) E= lens_E;
		else {
			E=sed_E_xL;
			// if(isGraded) E += (sed_E_x0-sed_E_xL)*(1.-std::min(globalpos[0],L1+dl1)/(L1+dl1));
			if(isGraded){
				double x0=0., y0=sed_E_x0;
				double x1=L1+dl1, y1=n*sed_E_xL;
				double x2=LT, y2=sed_E_xL;
				double a = ( x0*(y2-y1) + x1*(y0-y2) + x2*(y1-y0) )/(( x0-x2)*(x1-x0)*(x2-x1) );
				double b = (y2-y0)/(x2-x0) - a*(x2+x0);
				double c = y2 - a*x2*x2 - b*x2;
				E = a*globalpos[0]*globalpos[0] + b*globalpos[0] + c;
			}
		}
		return E;
	}

	double PoissonRatio(FieldVector globalpos) const {
		double nu=0.;
		if( isSandLens(globalpos) ) nu= lens_nu;
		else {
			nu=sed_nu_xL;
			// if(isGraded) nu += (sed_nu_x0-sed_nu_xL)*(1.-std::min(globalpos[0],L1+dl1)/(L1+dl1));
			if(isGraded){
				double x0=0., y0=sed_nu_x0;
				double x1=L1+dl1, y1=n*sed_nu_xL;
				double x2=LT, y2=sed_nu_xL;
				double a = ( x0*(y2-y1) + x1*(y0-y2) + x2*(y1-y0) )/(( x0-x2)*(x1-x0)*(x2-x1) );
				double b = (y2-y0)/(x2-x0) - a*(x2+x0);
				double c = y2 - a*x2*x2 - b*x2;
				nu = a*globalpos[0]*globalpos[0] + b*globalpos[0] + c;
			}
		}
		return nu;
	}

	std::vector<double> LameParameters(FieldVector globalpos) const {
	  // sigma = 2 mu epsilon + lambda tr(epsilon) I
	  // LP[0] = lambda
	  // LP[1] = mu

	  double nu = PoissonRatio(globalpos);
	  double E = YoungsModulus(globalpos);

	  if( nu == 0.5 ){
		  std::cerr<< "error in: "<< __FILE__ << '\t'
				   << "Poisson's ratio is 1/2." << std::endl;
		  exit(1);
	  }

	  std::vector<double> LP(2,0.);
	  LP[1] = E/(2.*(1+nu)); //mu
	  LP[0] = E*nu/((1+nu)*(1.-2*nu)); //lambda

	  return LP;
	}

	std::vector< std::vector<double> >
	ElasticStiffness(FieldVector globalpos) const {

	  auto LP = LameParameters(globalpos);

	  std::vector< std::vector<double> > D(StressIndices::num);
	  for( int i=0;i<StressIndices::num;i++ ){
		  D[i] = std::vector<double>(StressIndices::num,0.);
	  }
	  for( int i=0;i<DIMENSION;i++ ){
		  D[i][i] = 2*LP[1];
		  for( int j=0;j<DIMENSION;j++ ){
			  D[i][j] += 2.*LP[0];
		  }
	  }
	  for( int i=DIMENSION;i<StressIndices::num;i++ ){
		  D[i][i] = LP[1];
	  }
	  return D;
	}

	double Invariant_I1 ( std::vector<double> stress ) const {
		double I1 = 0.;
		for(int i=0; i<DIMENSION; i++){
			I1 += stress[i];
		}
		return I1;
	}

	// Mean Stress p
	double Stress_p ( std::vector<double> stress ) const {
		auto p = (1./DIMENSION)*Invariant_I1(stress);
		return p;
	}

	std::vector< std::vector<double> >
	P () const {
		std::vector< std::vector<double> > matP(StressIndices::num);
		for( int i=0;i<StressIndices::num;i++ ){
			matP[i] = std::vector<double>(StressIndices::num,0.);
		}

		for( int i=0; i<DIMENSION; i++ ){
			matP[i][i] = 1.;
			for( int j=0; j<DIMENSION; j++ ){
				matP[i][j] += -1./DIMENSION;
			}
		}

		for( int i=DIMENSION; i<StressIndices::num; i++ ){
			matP[i][i] = 2.;
		}

		return matP;
	}

	double Invariant_J2 ( std::vector<double> stress ) const {
		//J2 = (1./2) * var^T . matP . var

		auto matP = P();

		std::vector<double> tmp(stress.size(),0.);
		double J2 = 0.;
		for(int i=0;i<matP.size();i++){
			tmp[i] = 0.;
			for(int j=0;j<matP[i].size();j++){
				tmp[i] += matP[i][j] * stress[j];
			}
			J2 += (1./2.) * stress[i] * tmp[i];
		}

		return (J2+1.e-20);
	}

	// Dev stress q
	double Stress_q ( std::vector<double> stress ) const {
		auto q = std::sqrt( DIMENSION * Invariant_J2(stress) );
		return q;
	}

	double FrictionCoefficient(FieldVector globalpos) const{
		double alpha=0.;
		if( isSandLens(globalpos) ) alpha= lens_alpha;
		else {
			alpha=sed_alpha_xL;
			// if(isGraded) alpha += (sed_alpha_x0-sed_alpha_xL)*(1.-std::min(globalpos[0],L1+dl1)/(L1+dl1));
			if(isGraded){
				double x0=0., y0=sed_alpha_x0;
				double x1=L1+dl1, y1=n*sed_alpha_xL;
				double x2=LT, y2=sed_alpha_xL;
				double a = ( x0*(y2-y1) + x1*(y0-y2) + x2*(y1-y0) )/(( x0-x2)*(x1-x0)*(x2-x1) );
				double b = (y2-y0)/(x2-x0) - a*(x2+x0);
				double c = y2 - a*x2*x2 - b*x2;
				alpha = a*globalpos[0]*globalpos[0] + b*globalpos[0] + c;
			}
		}
		return alpha;
	}

	double CohesiveStrength(FieldVector globalpos) const{
		double Zmax = InitialSeafloorHeight(globalpos);
		double Fz_x0 = std::exp(-sed_beta_c_x0*(Zmax-globalpos[1]));
		double Fz_xL = std::exp(-sed_beta_c_xL*(Zmax-globalpos[1]));
		double Fz_lens = std::exp(-lens_beta_c*(Zmax-globalpos[1]));

		double cc=0.,Fz=0.;
		if( isSandLens(globalpos) ) {
			Fz=Fz_lens;
			cc=lens_c/Fz;
		}
		else {
			Fz=Fz_xL;
			cc=sed_c_xL/Fz_xL;
			// if(isGraded) {
			// 	cc  += (sed_c_x0/Fz_x0-sed_c_xL/Fz_xL)*(1.-std::min(globalpos[0],L1+dl1)/(L1+dl1));
			// }
			if(isGraded){
				double x0=0., y0=sed_c_x0/Fz_x0;
				double x1=L1+dl1, y1=n*sed_c_xL/Fz_xL;
				double x2=LT, y2=sed_c_xL/Fz_xL;
				double a = ( x0*(y2-y1) + x1*(y0-y2) + x2*(y1-y0) )/(( x0-x2)*(x1-x0)*(x2-x1) );
				double b = (y2-y0)/(x2-x0) - a*(x2+x0);
				double c = y2 - a*x2*x2 - b*x2;
				cc = a*globalpos[0]*globalpos[0] + b*globalpos[0] + c;
			}
		}
		return cc;
	}

	double YieldSurface( FieldVector globalpos, std::vector<double> stress )const{
		auto p = Stress_p(stress);
		auto q = Stress_q(stress);
		auto c = CohesiveStrength(globalpos);
		auto alpha = FrictionCoefficient(globalpos);

		auto F = q + alpha*p - c;

		return F;
	}

	/*****************************************************************/

};


struct ConvectionDiffusionBoundaryConditions
{
  enum Type { Dirichlet=1, Neumann=-1, Outflow=-2 };
};

template<typename GV,typename PTree,typename Parameters>
class BoundaryTypesFLOW : public Dune::PDELab::DirichletConstraintsParameters
{
private :
	const GV& gv ;
	const PTree& ptree;
	const Parameters& parameter;
	int var_id; /*PV IDs -> pw, por*/
	double *t;
	double *dt;
	double eps = 1.e-6;

public:
	  //! construct from grid view
	BoundaryTypesFLOW( const GV& gv_,
				   	   const PTree& ptree_,
					   const Parameters& parameter_,
					   int var_id_, double *t_, double *dt_ )
	: gv ( gv_ ),
	  ptree(ptree_),
	  parameter(parameter_),
	  var_id(var_id_), t(t_), dt(dt_)
	{}

	using Traits = Dune::PDELab::BoundaryGridFunctionTraits< GV,int,1,Dune::FieldVector<int,1> >;

	template<typename I>
	bool isDirichlet( const I& i, const typename Traits::DomainType& xlocal ) const {

		const int dim = Traits::GridViewType::Grid::dimension;
		typedef typename Traits::GridViewType::Grid::ctype ctype;
		auto x = i.geometry().global(xlocal);

		bool isDirichlet = false;

		if(var_id==FlowVariables::pw ){
			double hw = parameter.SeaLevel(x,((*t)+(*dt))) - parameter.InitialSeafloorHeight(x);
			if( parameter.isSeafloorBoundary(x) and hw>0.0 ){
				isDirichlet = true;
			}
		}
		
		if(var_id==FlowVariables::porosity ){
			if( parameter.isUpstreamBoundary(x) ) isDirichlet=true;
		}

		return isDirichlet;
	}

	// ! get a reference to the gridview
	inline const GV& getGridView () { return gv ; }

};


template<typename GV,typename PTree,typename Parameters>
class BoundaryTypesELASTICITY : public Dune::PDELab::DirichletConstraintsParameters
{
private :
	const GV& gv ;
	const PTree& ptree;
	const Parameters& parameter;
	int var_id; /*PV IDs -> ux, uy, uz*/
	double *t;
	double *dt;
	double eps = 1.e-6;

public:
	  //! construct from grid view
	BoundaryTypesELASTICITY( const GV& gv_,
				   const PTree& ptree_,
				   const Parameters& parameter_,
				   int var_id_, double *t_, double *dt_ )
	: gv ( gv_ ),
	  ptree(ptree_),
	  parameter(parameter_),
	  var_id(var_id_), t(t_), dt(dt_)
	{}

	using Traits = Dune::PDELab::BoundaryGridFunctionTraits< GV,int,1,Dune::FieldVector<int,1> >;

	template<typename I>
	bool isDirichlet( const I& i, const typename Traits::DomainType& xlocal ) const {

		const int dim = Traits::GridViewType::Grid::dimension;
		typedef typename Traits::GridViewType::Grid::ctype ctype;
		auto x = i.geometry().global(xlocal);

		bool isDirichlet = false;

		if( var_id==0 ){
			if( 	!parameter.isSeafloorBoundary(x)
				//and !parameter.isUpstreamBoundary(x)
				and !parameter.isDownstreamBoundary(x)
#if DIMENSION==3
				and !parameter.isFrontBoundary(x)
				and !parameter.isBackBoundary(x)
#endif
				//and !parameter.isBottomBoundary(x)
				//parameter.isBedrockBoundary(x)
			){
				isDirichlet = true;
			}
		}else if(var_id==1 ){
			if( 	!parameter.isSeafloorBoundary(x)
				and !parameter.isUpstreamBoundary(x)
				and !parameter.isDownstreamBoundary(x)
#if DIMENSION==3
				and !parameter.isFrontBoundary(x)
				and !parameter.isBackBoundary(x)
#endif
				//and !parameter.isBottomBoundary(x)
				//parameter.isBedrockBoundary(x)
			){
				isDirichlet = true;
			}
		}
#if DIMENSION==3
		else if(var_id==2 ){//Y-AXIS
					if( 	!parameter.isSeafloorBoundary(x)
						and !parameter.isUpstreamBoundary(x)
						and !parameter.isDownstreamBoundary(x)
						//and !parameter.isFrontBoundary(x)
						//and !parameter.isBackBoundary(x)
						//and !parameter.isBottomBoundary(x)
						//parameter.isBedrockBoundary(x)
					){
						isDirichlet = true;
					}
				}
#endif

		return isDirichlet;
	}

	// ! get a reference to the gridview
	inline const GV& getGridView () { return gv ; }

};


template<typename GV,typename PTree,typename Parameters>
class NeumannBoundaryValuesFLOW
		: public Dune::PDELab::BoundaryGridFunctionBase< Dune::PDELab::BoundaryGridFunctionTraits<GV,double,1,Dune::FieldVector<double,1>>,
		  	  	  	  	  	  	  	  	  	  	  	  	 NeumannBoundaryValuesFLOW<GV,PTree,Parameters>>
{
public :

	using Traits = Dune::PDELab::BoundaryGridFunctionTraits< GV,double,1,Dune::FieldVector<double,1>> ;

	// ! construct from gridview
	NeumannBoundaryValuesFLOW ( const GV& gv_,
				   	 	 		const PTree& ptree_,
								const Parameters& parameter_)
	: gv ( gv_ ),
	  ptree(ptree_),
	  parameter(parameter_)
	{
		F_up = ptree.get("flux.upstream",(double)1.0)*(1./(86400.0*365.0));//m/s
		F_down = ptree.get("flux.downstream",(double)0.0)*(1./(86400.0*365.0));//m/s
		F_rain = ptree.get("flux.rainfall",(double)0.0)*(1./(86400.0*365.0));//m/s
	}

	//! Neumann boundary condition
	template<typename I, typename X>
	double evaluate (const I& i, const X& x, double& t, int tag) const {
		auto xglobal = i.geometry().global(x);

		double value=0.0;
		if( tag==FlowVariables::pw){
			double hw = parameter.SeaLevel(xglobal,t) - parameter.InitialSeafloorHeight(xglobal);
			if( parameter.isUpstreamBoundary(xglobal) ) return F_up; //m/s
			else if( parameter.isDownstreamBoundary(xglobal) ) return F_down; //m/s
			else if( parameter.isSeafloorBoundary(xglobal) and hw<=0. ) return F_rain;
		}

		return value; //m/s
	}

private :
	const GV& gv ;
	const PTree& ptree;
	const Parameters& parameter;
	double F_up, F_down, F_rain;
} ;


template<typename GV,typename PTree,typename Parameters>
class NeumannBoundaryValuesELASTICITY
		: public Dune::PDELab::BoundaryGridFunctionBase< Dune::PDELab::BoundaryGridFunctionTraits<GV,double,1,Dune::FieldVector<double,1>>,
		  NeumannBoundaryValuesELASTICITY<GV,PTree,Parameters>>
{
public :

	using Traits = Dune::PDELab::BoundaryGridFunctionTraits< GV,double,1,Dune::FieldVector<double,1>> ;

	// ! construct from gridview
	NeumannBoundaryValuesELASTICITY ( const GV& gv_,
				   	 	 			  const PTree& ptree_,
									  const Parameters& parameter_)
	: gv ( gv_ ),
	  ptree(ptree_),
	  parameter(parameter_)
	{}

	//! Neumann boundary condition
	template<typename I, typename X>
	double evaluate (const I& i, const X& x, double& t, int tag) const {
		auto xglobal = i.geometry().global(x);

		double value=0.0;
		if( tag==0 ){
			value=0.;
		}else if( tag==1 ){
			value=0.;
		}
#if DIMENSION==3
		else if( tag==2 ){
			value=0.;
		}
#endif
		return value;
	}

private :
	const GV& gv ;
	const PTree& ptree;
	const Parameters& parameter;
} ;


/** \brief A function that defines Dirichlet boundary conditions AND
    its extension to the interior */
template<typename GV, class PTree, class Parameters>
class DirichletBoundaryValuesFLOW
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,double,1,Dune::FieldVector<double,1> >,
	DirichletBoundaryValuesFLOW<GV,PTree,Parameters> >
{
	const GV& gv ;
	const PTree& ptree;
	const Parameters& parameter;
	double *t;
	double *dt;
	int var_id;

public:
  typedef Dune::PDELab::GridFunctionTraits<GV,double,1,Dune::FieldVector<double,1> > Traits;
  typedef typename GV::IntersectionIterator IntersectionIterator;
  typedef typename IntersectionIterator::Intersection Intersection;

  //! construct from grid view
  DirichletBoundaryValuesFLOW (const GV& gv_,
						   	   const PTree& ptree_,
							   const Parameters& parameter_,
							   int var_id_,
							   double *t_,
							   double *dt_ )
  : gv ( gv_ ),
	ptree(ptree_),
	parameter(parameter_),
	var_id(var_id_),
	t( t_ ),
	dt( dt_ )
  {}

  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const {
//	  y=0.;

	  const int dim = Traits::GridViewType::Grid::dimension;
	  typedef typename Traits::GridViewType::Grid::ctype ctype;
	  auto x = e.geometry().global(xlocal);

	  double t_new = (*t)+(*dt);

	  if(var_id==FlowVariables::pw ){
		  double hw = parameter.SeaLevel(x,((*t)+(*dt))) - parameter.InitialSeafloorHeight(x);
		  if( parameter.isSeafloorBoundary(x) and hw>0.0 ){
			  y = parameter.WaterDensity()*9.81*std::max(hw,0.0);
		  }
		  else{
				auto Zmax = parameter.SeaLevel(x,((*t)+(*dt)));
				y=parameter.WaterDensity()*9.81*(Zmax-x[1]);
		  }
	  }
	  
	  if(var_id==FlowVariables::porosity ){
	  	y=parameter.InitialPorosity(x);
	  }
	  return;
  }

  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
};



/** \brief A function that defines Dirichlet boundary conditions AND
    its extension to the interior */
template<typename GV, class PTree, class Parameters>
class DirichletBoundaryValuesELASTICITY
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,double,1,Dune::FieldVector<double,1> >,
	DirichletBoundaryValuesELASTICITY<GV,PTree,Parameters> >
{
	const GV& gv ;
	const PTree& ptree;
	const Parameters& parameter;
	double *t;
	double *dt;
	int var_id;

public:
  typedef Dune::PDELab::GridFunctionTraits<GV,double,1,Dune::FieldVector<double,1> > Traits;
  typedef typename GV::IntersectionIterator IntersectionIterator;
  typedef typename IntersectionIterator::Intersection Intersection;

  //! construct from grid view
  DirichletBoundaryValuesELASTICITY (const GV& gv_,
						   	   	     const PTree& ptree_,
									 const Parameters& parameter_,
									 int var_id_,
									 double *t_,
									 double *dt_ )
  : gv ( gv_ ),
	ptree(ptree_),
	parameter(parameter_),
	var_id(var_id_),
	t( t_ ),
	dt( dt_ )
  {}

  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const {
//	  y=0.;

	  const int dim = Traits::GridViewType::Grid::dimension;
	  typedef typename Traits::GridViewType::Grid::ctype ctype;
	  auto x = e.geometry().global(xlocal);

	  double t_new = (*t)+(*dt);

	  if( var_id==0 ){
		  if( 	!parameter.isSeafloorBoundary(x)
				//and !parameter.isUpstreamBoundary(x)
				and !parameter.isDownstreamBoundary(x)
#if DIMENSION==3
				and !parameter.isFrontBoundary(x)
				and !parameter.isBackBoundary(x)
#endif
				//and !parameter.isBottomBoundary(x)
				//parameter.isBedrockBoundary(x)
		  ){
			  y=0.;
		  }
	  }else if(var_id==1 ){
		  if( 	!parameter.isSeafloorBoundary(x)
				  and !parameter.isUpstreamBoundary(x)
				  and !parameter.isDownstreamBoundary(x)
#if DIMENSION==3
				  and !parameter.isFrontBoundary(x)
				  and !parameter.isBackBoundary(x)
#endif
				  //and !parameter.isBottomBoundary(x)
				  //parameter.isBedrockBoundary(x)
		  ){
			  y=0.;
		  }
	  }
#if DIMENSION==3
	  else if(var_id==2 ){
		  if( 	!parameter.isSeafloorBoundary(x)
				  and !parameter.isUpstreamBoundary(x)
				  and !parameter.isDownstreamBoundary(x)
				  //and !parameter.isFrontBoundary(x)
				  //and !parameter.isBackBoundary(x)
				  //and !parameter.isBottomBoundary(x)
				  //parameter.isBedrockBoundary(x)
		  ){
			  y=0.;
		  }
	  }
#endif
	  return;
  }

  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
};

template<typename GV,typename PTree,typename Parameters>
class InitialConditions
{
private:
	const GV& gv;
	const PTree& ptree;
	const Parameters& param;
	const static int dim = GV::dimension;

public:
	//! construct from grid view
	InitialConditions ( const GV& gv_ ,
				 	 	const PTree& ptree_,
						const Parameters& param_)
	: gv( gv_ ) ,
	  ptree(ptree_),
	  param(param_)
	{}

	template<typename E, typename X>
	double pressure (const E& e, const X& x) const {
		auto xglobal = e.geometry().global(x);

		auto Zmax = param.SeaLevel(xglobal,0.0);
		double pw0=param.WaterDensity()*param.GravityVector()[1]*(Zmax-xglobal[1]);
//		double hw = param.SeaLevel(xglobal,0.0) - param.InitialSeafloorHeight(xglobal);
//		if( param.isSeafloorBoundary(xglobal) and hw>0.0 ){
//			  pw0 = param.WaterDensity()*9.81*std::max(hw,0.0);
//		}
		return pw0; //--
	}

	template<typename E, typename X>
	double porosity (const E& e, const X& x) const {
		auto xglobal = e.geometry().global(x);

		double por0=param.InitialPorosity(xglobal);
		return por0; //--
	}

	template<typename E, typename X>
	double displacement_u0 (const E& e, const X& x ) const {
		auto xglobal = e.geometry().global(x);
		double u=0.;
		return u; //
	}

	template<typename E, typename X>
	double displacement_u1 (const E& e, const X& x ) const {
		auto xglobal = e.geometry().global(x);
		double u=0.;
		return u; //
	}

#if DIMENSION==3
	template<typename E, typename X>
	double displacement_u2 (const E& e, const X& x ) const {
		auto xglobal = e.geometry().global(x);
		double u=0.;
		return u; //
	}
#endif

};

#endif /* LANDSLIDELANDSCAPE_DECOUPLED_NEW_PROBLEM_CONTINENTALMARGIN02_PARAMETERS_HH_ */
