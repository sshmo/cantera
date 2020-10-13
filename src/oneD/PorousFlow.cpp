//! @file PorousFlow.cpp

#include "cantera/oneD/PorousFlow.h"
#include "cantera/oneD/StFlow.h"
#include "cantera/base/ctml.h"
#include "cantera/transport/TransportBase.h"
#include "cantera/numerics/funcs.h"
#include "cantera/numerics/polyfit.h"
#include "cantera/oneD/OneDim.h"

using namespace std;

namespace Cantera
{

PorousFlow::PorousFlow(ThermoPhase* ph, size_t nsp, size_t points):
    StFlow(ph, nsp, points),
    pore1(0.835),pore2(0.87),
    diam1(0.00029),diam2(0.00152),
    scond1(1.3),scond2(1.771),  
	Omega1(0.8),Omega2(0.8), 
    srho(510),sCp(824),
    m_zmid(0.035),m_dzmid(0.002),
    m_adapt(0.1), m_porea(0.1), m_poreb(0.1), 
	m_porec(0.1), m_pored(0.1), m_diama(0.1), m_diamb(0.1), 
	m_diamc(0.1), m_diamd(0.1)
        {
	   Tw.resize(points);
	   dq.resize(points);
	   hconv.resize(points);
        }

void PorousFlow::setupGrid(size_t n, const doublereal* z)
{
    vector_fp TwTmp = Tw;
    vector_fp dqTmp = dq;
    Tw.resize(n);
    dq.resize(n);

    size_t j = 0;
    for (size_t i=0;i<n;i++)
    {
      if (z[i] <= m_z[0] )
      {
        Tw[i]=TwTmp[0];
        dq[i]=dqTmp[0];
      }
      else if (z[i] >= m_z[m_points-1] )
      {
        Tw[i]=TwTmp[m_points-1];
        dq[i]=dqTmp[m_points-1];
      }
      else 
      {
        while ( z[i] > m_z[j+1] )
        {
          j++;
          if ( j+1 > m_points-1 ) 
          {
            throw 10;
          }
        }
        double tmp = (z[i]-m_z[j])/(m_z[j+1]-m_z[j]);
        Tw[i] = (1.0-tmp)*TwTmp[j] + tmp*TwTmp[j+1];
        dq[i] = (1.0-tmp)*dqTmp[j] + tmp*dqTmp[j+1];
      }
    }
	StFlow::setupGrid(n,z);

    ifstream in("Properties.txt"); //Read in the solid properties
    double proper;
    in>>proper;
    pore1=proper;
    in>>proper;
    pore2=proper;
    in>>proper;
    diam1=proper;
    in>>proper;
    diam2=proper;
    in>>proper;
    Omega1=proper;
    in>>proper;
    Omega2=proper;
    in>>proper;
    srho=proper;
    in>>proper;
    sCp=proper;
    in>>proper;
    m_zmid=proper;
    in>>proper;
    m_dzmid=proper;
    in.close();
    
}

void PorousFlow::eval(size_t jg, doublereal* xg, doublereal* rg, integer* diagg, doublereal rdt)
{
    // if evaluating a Jacobian, and the global point is outside
    // the domain of influence for this domain, then skip
    // evaluating the residual
    if (jg != npos && (jg + 1 < firstPoint() || jg > lastPoint() + 1)) {
        return;
    }

    // start of local part of global arrays
    doublereal* x = xg + loc();
    doublereal* rsd = rg + loc();
    integer* diag = diagg + loc();

    size_t jmin, jmax;

    if (jg == npos) {      // evaluate all points
        jmin = 0;
        jmax = m_points - 1;
    } else {          // evaluate points for Jacobian
        size_t jpt = (jg == 0) ? 0 : jg - firstPoint();
        jmin = std::max<size_t>(jpt, 1) - 1;
        jmax = std::min(jpt+1,m_points-1);
    }

    // properties are computed for grid points from j0 to j1
    size_t j0 = std::max<size_t>(jmin, 1) - 1;
    size_t j1 = std::min(jmax+1,m_points-1);

    size_t j, k;
    m_dovisc = 1;

    //-----------------------------------------------------
    //              update properties
    //-----------------------------------------------------

    updateThermo(x, j0, j1);
    // update transport properties only if a Jacobian is not
    // being evaluated
    if (jg == npos) {
        updateTransport(x, j0, j1);
    }

    // update the species diffusive mass fluxes whether or not a
    // Jacobian is being evaluated
    updateDiffFluxes(x, j0, j1);


    //----------------------------------------------------
    // evaluate the residual equations at all required
    // grid points
    //----------------------------------------------------

    doublereal sum, sum2, dtdzj;

    doublereal lam, visc, Re; //Defining new variables.
    int length=m_points;
    hconv.resize(length);  

    //initialize property vectors
    //
    //
    pore.resize(length);
    diam.resize(length);
    scond.resize(length);
    //vector<double> scond(length);
    vector<double> Omega(length);
    vector<double> Cmult(length);
    vector<double> mpow(length);
    vector<double> RK(length);
   
    for (int i=0; i<=length-1;i++)
    {
      if (z(i)<m_zmid-m_dzmid)
      {
        pore[i]=pore1;
        diam[i]=diam1;
      }
      else if (z(i)>m_zmid+m_dzmid)
      {
        pore[i]=pore2;
        diam[i]=diam2;
      }
      else
      {
        pore[i]=(((pore2-pore1)/(2*m_dzmid))*(z(i)-(m_zmid-m_dzmid) ))+pore1;
        diam[i]=(((diam2-diam1)/(2*m_dzmid))*(z(i)-(m_zmid-m_dzmid) ))+diam1;
	
       }
       RK[i]=(3*(1-pore[i])/diam[i]);   //extinction coefficient, PSZ, Hsu and Howell(1992)
       Cmult[i]=-400*diam[i]+0.687;	// Nusselt number coefficients
       mpow[i]=443.7*diam[i]+0.361;
       scond[i]=0.188-17.5*diam[i];    //solid phase thermal conductivity, PSZ, Hsu and Howell(1992) 
    }
    for (int i=0; i<=length-1;i++)
    {
       if (z(i)<m_zmid)
       {
          Omega[i]=Omega1;		//scattering albedo/extinction
       }
       else
       {
          Omega[i]=Omega2;
       }
    }
    
   int solidenergy=0;
    //loop over gas energy vector. If it is going to be solved then find hv
    for(j=jmin;j<=jmax;j++)
    {
       //std::cout << "gas energy" << m_do_energy[j] << std::endl;
       solidenergy+=m_do_energy[j];
    }
    solidenergy=1;
//    std::cout << "solid energy" << solidenergy << std::endl;
    if (solidenergy!=0)
    {
       for (j = jmin; j <= jmax; j++)
       {
          lam=m_tcon[j];	//Gas phase thermal conductivity 
          visc=m_visc[j];
            
          Re= (rho_u(x,j)*pore[j]*diam[j])/visc;
          doublereal nusselt = Cmult[j]*pow(Re,mpow[j]);
          hconv[j] = (lam * nusselt)/pow(diam[j],2);
          
       }
       //m_dosolid = true;
      //std::cout << "dosolid" << Domain1D::container().dosolid << std::endl;   

     //solid(x,hconv,scond,RK,Omega,srho,sCp,rdt);
       if (container().dosolid==1 )
       {
          solid(x,hconv,scond,RK,Omega,srho,sCp,rdt);
          (*m_container).dosolid=0;
       }
    }


    for (j = jmin; j <= jmax; j++) {


        //----------------------------------------------
        //         left boundary
        //----------------------------------------------

        if (j == 0) {

            // these may be modified by a boundary object

            // Continuity. This propagates information right-to-left,
            // since rho_u at point 0 is dependent on rho_u at point 1,
            // but not on mdot from the inlet.
            rsd[index(c_offset_U,0)] =
                -(rho_u(x,1) - rho_u(x,0))/m_dz[0]
                -(density(1)*V(x,1) + density(0)*V(x,0));

            // the inlet (or other) object connected to this one
            // will modify these equations by subtracting its values
            // for V, T, and mdot. As a result, these residual equations
            // will force the solution variables to the values for
            // the boundary object
            rsd[index(c_offset_V,0)] = V(x,0);
            if (doEnergy(0)) {
                rsd[index(c_offset_T,0)] = T(x,0);
            } else {
                rsd[index(c_offset_T,0)] = T(x,0) - T_fixed(0);
            }
            rsd[index(c_offset_L,0)] = -rho_u(x,0);

            // The default boundary condition for species is zero
            // flux. However, the boundary object may modify
            // this.
            sum = 0.0;
            for (k = 0; k < m_nsp; k++) {
                sum += Y(x,k,0);	//MODIFIED
		rsd[index(c_offset_Y + k, 0)] =
                    -(m_flux(k,0) + rho_u(x,0)* Y(x,k,0));
            }
            rsd[index(c_offset_Y + leftExcessSpecies(), 0)] = 1.0 - sum;

            // set residual of poisson's equ to zero
            rsd[index(c_offset_E, 0)] = x[index(c_offset_E, j)];
        } else if (j == m_points - 1) {
            evalRightBoundary(x, rsd, diag, rdt);
            // set residual of poisson's equ to zero
            rsd[index(c_offset_E, j)] = x[index(c_offset_E, j)];
        } else { // interior points
            //evalContinuity(j, x, rsd, diag, rdt,pore);
            rsd[index(c_offset_U,j)] =
               -(rho_u(x,j+1)*pore[j+1] - rho_u(x,j)*pore[j])/m_dz[j] //added porosity
               -(density(j+1)*V(x,j+1) + density(j)*V(x,j));

            diag[index(c_offset_U,j)] = 0.0;
            
            // set residual of poisson's equ to zero
            rsd[index(c_offset_E, j)] = x[index(c_offset_E, j)];
            //------------------------------------------------
            //    Radial momentum equation
            //
            //    \rho dV/dt + \rho u dV/dz + \rho V^2
            //       = d(\mu dV/dz)/dz - lambda
            //
            //-------------------------------------------------
            rsd[index(c_offset_V,j)]
            = (shear(x,j) - lambda(x,j) - rho_u(x,j)*dVdz(x,j)
               - m_rho[j]*V(x,j)*V(x,j))/m_rho[j]
              - rdt*(V(x,j) - V_prev(j));
            diag[index(c_offset_V, j)] = 1;

            //-------------------------------------------------
            //    Species equations
            //
            //   \rho dY_k/dt + \rho u dY_k/dz + dJ_k/dz
            //   = M_k\omega_k
            //
            //-------------------------------------------------
            getWdot(x,j);

            doublereal convec, diffus;
            for (k = 0; k < m_nsp; k++) {
                convec = rho_u(x,j)*dYdz(x,k,j)*pore[j]; //added porosity
                diffus = 2.0*(m_flux(k,j)*pore[j] - m_flux(k,j-1)*pore[j-1]) //added porosity, m_flux is  mass flux of species k in units of kg m-3 s-1
                         /(z(j+1) - z(j-1));
                rsd[index(c_offset_Y + k, j)]
                = (m_wt[k]*(wdot(k,j)*pore[j]) //added porosity
                   - convec - diffus)/(m_rho[j]*pore[j]) //added porosity
                  - rdt*(Y(x,k,j) - Y_prev(k,j));
                diag[index(c_offset_Y + k, j)] = 1;
            }

            //-----------------------------------------------
            //    energy equation
            //
            //    \rho c_p dT/dt + \rho c_p u dT/dz
            //    = d(k dT/dz)/dz
            //      - sum_k(\omega_k h_k_ref)
            //      - sum_k(J_k c_p_k / M_k) dT/dz
            //-----------------------------------------------

            if (m_do_energy[j]) {

                setGas(x,j);

                // heat release term
                const vector_fp& h_RT = m_thermo->enthalpy_RT_ref();
                const vector_fp& cp_R = m_thermo->cp_R_ref();

                sum = 0.0; //chemical source term sum(wdot*h_RT)*R*T*porosity ( GasConstant is R)
                sum2 = 0.0; // diffusive molar flux
                doublereal flxk;
                for (k = 0; k < m_nsp; k++) {
                    flxk = 0.5*(m_flux(k,j-1) + m_flux(k,j));
                    sum += wdot(k,j)*h_RT[k];
                    sum2 += flxk*cp_R[k]/m_wt[k];
                }
                sum *= GasConstant * T(x,j);
                dtdzj = dTdz(x,j);
                sum2 *= GasConstant * dtdzj;
                rsd[index(c_offset_T, j)]   =
                    - m_cp[j]*rho_u(x,j)*dtdzj
                    - divHeatFlux(x,j) - sum - sum2;
                rsd[index(c_offset_T, j)] =  rsd[index(c_offset_T, j)] 
                                         -(hconv[j]*(T(x,j)-Tw[j]))/pore[j]; //added convective term
                rsd[index(c_offset_T, j)] /= (m_rho[j]*m_cp[j]);

                rsd[index(c_offset_T, j)] -= rdt*(T(x,j) - T_prev(j));
                diag[index(c_offset_T, j)] = 1;

		//Porosity related modifcations
               // sum *= GasConstant * T(x,j) * pore[j]; //added porosity 
	       // dtdzj = (T(x, j)*pore[j] - T(x, j-1)*pore[j-1])/m_dz[j-1]; //added porosity
               // sum2 *= GasConstant * dtdzj;

	       // doublereal c1 = m_tcon[j-1]*(T(x,j)*pore[j] - T(x,j-1)*pore[j]);
	       // doublereal c2 = m_tcon[j]*(T(x,j+1)*pore[j+1] - T(x,j)*pore[j]);
               // doublereal divHeatFlux_pore = -2.0*(c2/(z(j+1) - z(j)) - c1/(z(j) - z(j-1)))/(z(j+1) - z(j-1));
   
               //  rsd[index(c_offset_T, j)]   =
               //     - m_cp[j]*rho_u(x,j)*pore[j]*dTdz(x,j)//added porosity, advective term
               //    - divHeatFlux_pore - sum - sum2; //divHeatFlux is conductive heat flux/ replaced with porosity term in temperature gradient
               // 
               // rsd[index(c_offset_T, j)] -= (hconv[j]*(T(x,j)-Tw[j])); //added convective term
               // rsd[index(c_offset_T, j)] /= (m_rho[j]*m_cp[j]*pore[j]); //added porosity

               // rsd[index(c_offset_T, j)] -= rdt*(T(x,j) - T_prev(j)); //transient term
               // diag[index(c_offset_T, j)] = 1;
            } else {
                // residual equations if the energy equation is disabled
                rsd[index(c_offset_T, j)] = T(x,j) - T_fixed(j);
                diag[index(c_offset_T, j)] = 0;
            }

            rsd[index(c_offset_L, j)] = lambda(x,j) - lambda(x,j-1);
            diag[index(c_offset_L, j)] = 0;
        }
    }
}

void PorousFlow::restore(const XML_Node& dom, doublereal* soln, int loglevel)
{
	StFlow::restore(dom,soln,loglevel);
	vector_fp x;
	if (dom.hasChild("Solid")) {
        XML_Node& ref = dom.child("Solid");
        
        pore1 = getFloat(ref, "pore1" );
        pore2 = getFloat(ref, "pore2" );
        diam1 = getFloat(ref, "diam1" );
        diam2 = getFloat(ref, "diam2" );
        scond1 = getFloat(ref, "scond1");
	scond2 =  getFloat(ref, "scond2");
	Omega1= getFloat(ref, "Omega1");
        Omega2= getFloat(ref, "Omega2");
        srho  = getFloat(ref, "rho"   );
        sCp   = getFloat(ref, "Cp"   );
        
        m_zmid = getFloat(ref, "zmid"  );
        m_dzmid= getFloat(ref, "dzmid" );
        
	getFloatArray(ref, x, false, "", "Tsolid");
	Tw.resize(nPoints());
	if (x.size() == nPoints()) {
            for (size_t i = 0; i < x.size(); i++) {
                Tw[i] = x[i];
            }
        } else if (!x.empty()) {
            throw CanteraError("PorousFlow::restore", "Tw is of length" +
                               to_string(x.size()) + "but should be length" +
                               to_string(nPoints()));
        }
	getFloatArray(ref, x, false, "", "Radiation");
	dq.resize(nPoints());
	if (x.size() == nPoints()) {
            for (size_t i = 0; i < x.size(); i++) {
                dq[i] = x[i];
            }
        } else if (!x.empty()) {
            throw CanteraError("PorousFlow::restore", "dq is of length" +
                               to_string(x.size()) + "but should be length" +
                               to_string(nPoints()));
        }
        getFloatArray(ref, x, false, "", "Porosity");
	pore.resize(nPoints());
	if (x.size() == nPoints()) {
            for (size_t i = 0; i < x.size(); i++) {
                pore[i] = x[i];
            }
        } else if (!x.empty()) {
            throw CanteraError("PorousFlow::restore", "Porosity is of length" +
                               to_string(x.size()) + "but should be length" +
                               to_string(nPoints()));
        }
        getFloatArray(ref, x, false, "", "Diameter");
		diam.resize(nPoints());
		if (x.size() == nPoints()) {
            for (size_t i = 0; i < x.size(); i++) {
                diam[i] = x[i];
            }
        } else if (!x.empty()) {
            throw CanteraError("PorousFlow::restore", "diam is of length" +
                               to_string(x.size()) + "but should be length" +
                               to_string(nPoints()));
        }

        getFloatArray(ref, x, false, "", "SolidConductivity");
		scond.resize(nPoints());
		if (x.size() == nPoints()) {
            for (size_t i = 0; i < x.size(); i++) {
                scond[i] = x[i];
            }
        } else if (!x.empty()) {
            throw CanteraError("PorousFlow::restore", "scond is of length" +
                               to_string(x.size()) + "but should be length" +
                               to_string(nPoints()));
        }

	getFloatArray(ref, x, false, "", "Hconv");
	hconv.resize(nPoints());
	if (x.size() == nPoints()) {
            for (size_t i = 0; i < x.size()-1; i++) { 
                hconv[i] = x[i];
            }
        } else if (!x.empty()) {
            throw CanteraError("PorousFlow::restore", "hconv is of length" +
                               to_string(x.size()) + "but should be length" +
                               to_string(nPoints()));
        }

    }
}

XML_Node& PorousFlow::save(XML_Node& o, const doublereal* const sol)
{
    XML_Node& flow = StFlow::save(o, sol);

	vector_fp values(nPoints());
	XML_Node& solid = flow.addChild("Solid");
    
    // addFloat(solid, "pore1",  pore1  );
    // addFloat(solid, "pore2",  pore2  );
    // addFloat(solid, "diam1",  diam1  );
    // addFloat(solid, "diam2",  diam2  );
    // addFloat(solid, "scond1",  scond1  );
    // addFloat(solid, "scond2",  scond2  );
    // addFloat(solid, "Omega1", Omega1 );
    // addFloat(solid, "Omega2", Omega2 );
    // addFloat(solid, "rho",    srho   );
    // addFloat(solid, "Cp",    sCp    );
    // addFloat(solid, "zmid" ,  m_zmid );
    // addFloat(solid, "dzmid",  m_dzmid);

    for (size_t i = 0; i < nPoints(); i++) {
        values[i] = Tw[i];
    }
    addNamedFloatArray(solid, "Tsolid", nPoints(), &values[0]);
    

    for (size_t i = 0; i < nPoints(); i++) {
        values[i] = dq[i];
    }
    addNamedFloatArray(solid, "Radiation", nPoints(), &values[0]);

    for (size_t i = 0; i < nPoints(); i++) {
        values[i] = pore[i];
    }
    addNamedFloatArray(solid, "Porosity", nPoints(), &values[0]);
    
    for (size_t i = 0; i < nPoints(); i++) {
        values[i] = diam[i];
    }
    addNamedFloatArray(solid, "Diameter", nPoints(), &values[0]);

    for (size_t i = 0; i < nPoints(); i++) {
        values[i] = scond[i];
    }
    addNamedFloatArray(solid, "SolidConductivity", nPoints(), &values[0]);

    for (size_t i = 0; i < nPoints()-1; i++) { 
        values[i] = hconv[i];
    }
    addNamedFloatArray(solid, "Hconv", nPoints(), &values[0]);

    return flow;
}

void PorousFlow::solid(doublereal* x, vector<double> &hconv, vector<double>& scond, vector<double>& RK, vector<double>&Omega,double & srho,double & sCp, double rdt) 
//Solid solver
{
   int length=m_points; //
   Twprev = Tw;
   
      //Start of Conduction Radiation Stuff
   //
   //Vector Iinitialization
   //
   vector<double> edia(length);
   vector<double> fdia(length);
   vector<double> gdia(length);
   vector<double> rhs(length);
   vector<double> dqnew(length);
   double sigma=5.67e-8;
   double change1=1;

   //Vector Population
   for(int i=0;i<=length-1;i++)
   {
      dq[i]=0;
   }
   //double T0=300;
   //double T1=300;
   int count1=0;
   int fail1=0;
   while (change1>0.000001)
   {
      count1=count1+1;
      for(int i=0;i<=length-1;i++)
      {
         if (i==0)
         {
            edia[i]=0;
            fdia[i]=1;
            gdia[i]=-1;
            //gdia[i]=0;
            //Tw[i]  =T(x,0);
            rhs[i]=0;
         }
         else if (i==length-1)
         {
            edia[i]=-1;
            //edia[i]=0;
            fdia[i]=1;
            gdia[i]=0;
            rhs[i]=0;
            //Tw[i]  =1000.0;
            //Tw[i]  =0.0;
         }
         else
         {
            //std::cout << pore[i] << std::endl;
	    edia[i]=(2*scond[i])/((z(i)-z(i-1))*(z(i+1)-z(i-1)));
            //fdia[i]=-(2*scond[i])/((z(i+1)-z(i))*(z(i+1)-z(i-1)))-(2*scond[i])/((z(i)-z(i-1))*(z(i+1)-z(i-1)))-hconv[i]-srho*(1-pore[i])*sCp*rdt; //added porosity
            fdia[i]=-(2*scond[i])/((z(i+1)-z(i))*(z(i+1)-z(i-1)))-(2*scond[i])/((z(i)-z(i-1))*(z(i+1)-z(i-1)))-hconv[i]-srho*sCp*rdt;
            gdia[i]=(2*scond[i])/((z(i+1)-z(i))*(z(i+1)-z(i-1)));
            //Tw[i]=-hconv[i]*T(x,i)+dq[i]-srho*sCp*rdt*Twprev[i];
            rhs[i]=-hconv[i]*T(x,i)+dq[i]-srho*sCp*rdt*Twprev[i]; 
         }
      }

     
      //Decomposition
      for(int i=1;i<=length-1;i++)
      {
         edia[i]=edia[i]/fdia[i-1];
         fdia[i]=fdia[i]-edia[i]*gdia[i-1];
      }

      //Forward Substitution
      for(int i=1;i<=length-1;i++)
      {
         rhs[i]=rhs[i]-edia[i]*rhs[i-1];
      }

      //Back Substitution
      Tw[length-1]=rhs[length-1]/fdia[length-1];
      for(int i=length-2;i>=0;i--)
      {
         Tw[i]=(rhs[i]-gdia[i]*Tw[i+1])/fdia[i];
      }
      //T0=Tw[0];
      //T1=Tw[length-1];

      //Radiation Time
      //Vector Initialization
      vector<double> qplus(length);
      vector<double> qpnew(length);
      vector<double> qminus(length);
      vector<double> qmnew(length);
      double change2=1;

      //Vector Population
      double temp2 = T(x,0);
      //double temp2 = 300;
      for(int i=0;i<=length-1;i++)
      {
         if (i==0)
         {
            //double temp=Tw[i];
            //qplus[i]=sigma*pow(temp,4);
            //qpnew[i]=sigma*pow(temp,4);
            qplus[i]=sigma*pow(temp2,4);
            qpnew[i]=sigma*pow(temp2,4);
            qminus[i]=0;
            qmnew[i]=0;
         }
         else if (i==length-1)
         {
            //double temp=Tw[i];
            qplus[i]=0;
            qpnew[i]=0;
            //qminus[i]=sigma*pow(temp,4);
            //qmnew[i]=sigma*pow(temp,4);
            qminus[i]=sigma*pow(temp2,4);
            qmnew[i]=sigma*pow(temp2,4);
            //qminus[i]=0.0;
            //qmnew[i]=0.0;
         }
         else
         {
            qplus[i]=0;
            qpnew[i]=0;
            qminus[i]=0;
            qmnew[i]=0;
         }
      }
      int count=0;
      int fail=0;
      //S2 method
      while (change2>0.000001)
      {
         count=count+1;
         for(int i=1;i<=length-1;i++)
         {
            double temp=Tw[i];
            qpnew[i]=(qpnew[i-1]+RK[i]*(z(i)-z(i-1))*Omega[i]*qminus[i]+
                  2*RK[i]*(z(i)-z(i-1))*(1-Omega[i])*sigma*pow(temp,4))/
                  (1+(z(i)-z(i-1))*RK[i]*(2-Omega[i]) );
         }
         for(int i=length-2;i>=0;i--)
         {
            double temp=Tw[i];
            qmnew[i]=(qmnew[i+1]+RK[i]*(z(i+1)-z(i))*Omega[i]*qpnew[i]+
                  2*RK[i]*(z(i+1)-z(i))*(1-Omega[i])*sigma*pow(temp,4))/
                  (1+(z(i+1)-z(i))*RK[i]*(2-Omega[i]));
         }
         double norm1=0;
         double norm2=0;
         for(int i=0;i<=length-1;i++)
         {
            norm1+=(qpnew[i]-qplus[i])*(qpnew[i]-qplus[i]);
            norm2+=(qmnew[i]-qminus[i])*(qmnew[i]-qminus[i]);
            qplus[i]=qpnew[i];
            qminus[i]=qmnew[i];
         }
         norm1=sqrt(norm1);
         norm2=sqrt(norm2);
         if (count>100)
         {
            change2=0;
            fail=1;
         }
         else
         {
            change2=max(norm1,norm2);
         }
      }
      if (fail==1)
      {
         for(int i=0;i<=length-1;i++)
         {
            dqnew[i]=dq[i];
         }
         writelog("Rad Stall");
      }
      else
      {
         for(int i=0;i<=length-1;i++)
         {
            double temp=Tw[i];
            dqnew[i]=4*RK[i]*(1-Omega[i])*(sigma*pow(temp,4)-0.5*qplus[i]-0.5*qminus[i]); 
         }
      }
      double norm=0;
      double a=0.1;
      for (int i=0;i<=length-1;i++)
      {
         norm+=(dqnew[i]-dq[i])*(dqnew[i]-dq[i]);
         dq[i]=a*dqnew[i]+(1-a)*dq[i];
      }
      if (count1>400)
      {
         fail1=1;
         change1=0;
      }
      else
      {
         change1=sqrt(norm);
      }
   }
   if (fail1==1)
   {
      for (int i=0;i<=length-1;i++)
      {
         Tw[i]=Twprev[i];
      }
      writelog("Rad not Converged");
   }

   if ( m_refiner ) {
      refiner().setExtraVar(Tw.data());
   }


}

}// namespace