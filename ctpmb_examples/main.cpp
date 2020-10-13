#include <iostream>
#include <vector>
#include <fstream>
#include <time.h>

#include "cantera/oneD/Sim1D.h"
#include "cantera/oneD/Boundary1D.h"
#include "cantera/oneD/StFlow.h"
#include "cantera/oneD/PorousFlow.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/transport.h"

using namespace std;
using namespace Cantera;
using fmt::print;


void Combust()
{
	//New Parameters to allow for correlations
	//Define Material Properties
	double pore1 = 0.69;//0.835;
	double pore2 = 0.7;//0.870;
	double diam1 = 0.0015;//.00029;
	double diam2 = 0.0015;//0.00152
	double Omega1 = 0.8;
	double Omega2 = 0.8;
	double srho = 510;
	double sCp = 824;
	double m_zmid =0.035;
	double m_dzmid =0.002;
	//Export Proerties to fill
	ofstream fid("Properties.txt");
	fid << pore1 << endl;
	fid << pore2 << endl;
	fid << diam1 << endl;
	fid << diam2 << endl;
	fid << Omega1 << endl;
	fid << Omega2 << endl;
	fid << srho << endl;
	fid << sCp << endl;
    fid << m_zmid << endl;
	fid << m_dzmid << endl;
	fid.close();


    //auto sol = newSolution("drm19.yaml", "drm19", "None");
	auto sol = newSolution("gri30.yaml", "gri30", "None");
    auto gas = sol->thermo();
	//Define Fluid Properties
	double Pressure = OneAtm;
	double Tburner = 300;
	double uin = 0.45; //0.45
	double phi = 2.5;

    size_t nsp = gas->nSpecies();
    vector_fp x(nsp, 0.0);

    gas->setEquivalenceRatio(phi, "CH4", "O2:0.9,N2:0.1"); //"CH4", "O2:0.21,N2:0.78, AR:0.01"
    gas->setState_TP(Tburner, Pressure);
    gas->getMoleFractions(x.data());

	double rho_in = gas->density();

    vector_fp yin(nsp);
    gas->getMassFractions(&yin[0]);

    gas->equilibrate("HP");
	std::cout << gas->report() << std::endl;

	vector_fp yout(nsp);
	gas->getMassFractions(&yout[0]);

    double rho_out = gas->density();
    double Tad = gas->temperature();
    print("phi = {}, Tad = {}\n", phi, Tad);

    //=============  build each domain ========================


    //-------- step 1: create the flow -------------

	//Create the flow object
    PorousFlow flow(gas);
    flow.setAxisymmetricFlow();

	// create an initial grid
//	int nn = 1;
//	int nz = 300*nn+1;
	int nz = 6;
	double lz = 0.2;//0.1
	vector_fp z(nz);

	double dz = lz/((double)(nz-1));
//	double dz1 = 0.0006/nn;
//	double dz2 = 0.00005/nn;
//	double dz3 = 0.00041/nn;

	for (int iz = 0; iz < nz; iz++) {
        z[iz] = ((double)iz) * dz;
//		if (iz <= 50*nn)
//		{
//			z[iz] = ((double)iz) * dz1;
//		}
//		else if (iz <= 250*nn)
//		{
//			z[iz] = .03 + ((double)(iz - 50*nn)) * dz2;
//		}
//		else
//		{
//			z[iz] = .04 + ((double)(iz - 250*nn)) * dz3;
//		}
	}

	flow.setupGrid(nz, &z[0]);

    // specify the objects to use to compute kinetic rates and
    // transport properties

	std::unique_ptr<Transport> trmix(newTransportMgr("Mix", sol->thermo().get()));
	std::unique_ptr<Transport> trmulti(newTransportMgr("Multi", sol->thermo().get()));

	flow.setTransport(*trmix);
	flow.setKinetics(*sol->kinetics());
	flow.setPressure(Pressure);

	//------- step 2: create the inlet  -----------------------

	Inlet1D inlet;

	inlet.setMoleFractions(x.data());
	double mdot=uin*rho_in;
	inlet.setMdot(mdot);
	inlet.setTemperature(Tburner);


	//------- step 3: create the outlet  ---------------------

	Outlet1D outlet;

	//=================== create the container and insert the domains =====

	std::vector<Domain1D*> domains { &inlet, &flow, &outlet };
	Sim1D flame(domains);

	 //----------- Supply initial guess----------------------

//	//Build initial guess
//	vector_fp locs;
//	vector_fp value;
//	double z1 = 0.55;
//	double z2 = 0.62;
//	double uout = inlet.mdot() / rho_out;
//	//Velocity Profile
//	locs.resize(2);
//	value.resize(2);
//	locs[0] = 0;
//	locs[1] = 1;
//	value[0] = uin;
//	value[1] = uout;
//	flame.setInitialGuess("u", locs, value);
//	//Species Profiles
//	locs.resize(3);
//	value.resize(3);
//	locs[0] = 0;
//	locs[1] = z1;
//	locs[2] = 1;
//	for (int i = 0; i < nsp; i++)
//	{
//		value[0] = yin[i];
//		value[1] = yout[i];
//		value[2] = yout[i];
//		flame.setInitialGuess(gas->speciesName(i), locs, value);
//	}
//	//Temperature Profile
//	locs.resize(4);
//	value.resize(4);
//	locs[0] = 0;
//	locs[1] = z1;
//	locs[2] = z2;
//	locs[3] = 1;
//	value[0] = Tburner;
//	value[1] = Tburner;
//	value[2] = 2000;
//	value[3] = Tad;
//	flame.setInitialGuess("T", locs, value);

    vector_fp locs{0.0, 0.3, 0.7, 1.0};
	vector_fp value;

	double uout = inlet.mdot()/rho_out;
	value = {uin, uin, uout, uout};
	flame.setInitialGuess("velocity",locs,value);
	value = {Tburner, Tburner, Tad, Tad};
	flame.setInitialGuess("T",locs,value);

	for (size_t i=0; i<nsp; i++) {
		value = {yin[i], yin[i], yout[i], yout[i]};
		flame.setInitialGuess(gas->speciesName(i),locs,value);
	}

	inlet.setMoleFractions(x.data());
    inlet.setMdot(mdot);
    inlet.setTemperature(Tburner);

	flame.showSolution();


	double rtolSS = 1.0e-4;
	double atolSS = 1.0e-9;
	double rtolTS = 1.0e-4;
	double atolTS = 1.0e-9;
	flow.setSteadyTolerances(rtolSS, atolSS);
	flow.setTransientTolerances(rtolTS, atolTS);

	double SSJacAge = 5;
	double TSJacAge = 10;
	flame.setJacAge(SSJacAge, TSJacAge);

	int flowdomain = 1;
	double ratio = 4;
	double slope = 0.4;
	double curve = 0.4;
	//double prune = -0.001;

	flame.setRefineCriteria(flowdomain,ratio,slope,curve);//,prune);

	// Solve freely propagating flame

	// Linearly interpolate to find location where this temperature would
	// exist. The temperature at this location will then be fixed for
	// remainder of calculation.
	flame.setFixedTemperature(0.5 * (Tburner + Tad));
	int loglevel = 1;
	bool refine_grid = true;
	//flame.solve(loglevel,refine_grid);

	//now switch to multicomponent transport
	//flow.setTransport(*trmulti);
	//flame.solve(loglevel, refine_grid);

	// now enable Soret diffusion
	//flow.enableSoret(true);
	//flame.solve(loglevel, refine_grid);

	// now enable Energy equation
    //refine_grid = true;
    flow.solveEnergyEqn();
	flame.solve(loglevel, refine_grid);

	flame.save("gradienttest.xml", "run", "solution without energy equation");
	flame.writeStats();

    // open fod stream and save the results in results.txt
    ofstream fod("results.txt");

    // get t_wall and save it as Tsolid
    fod << "Tsolid" << ",";
    for (size_t i = 0; i < flow.nPoints(); i++) {
        fod << flow.getTw(i) << ",";
    }
    fod << endl;

    // get Hvonv and save it as h_conv
    fod << "h_conv" << ",";
    for (size_t i = 0; i < flow.nPoints(); i++) {
        fod << flow.getHconv(i) << ",";
    }
    fod << endl;

    // get Tprev and save it as Tgas
    fod << "Tgas" << ",";
    for (size_t i = 0; i < flow.nPoints(); i++) {
        fod << flow.getT(i) << ",";
    }
    fod << endl;

    // get m_z and save it as Z
    fod << "Z" << ",";
    for (size_t i = 0; i < flow.nPoints(); i++) {
        fod << flow.getZ(i) << ",";
    }
    fod << endl;

    for (size_t k = 0; k < nsp; k++) {

            fod << gas->speciesName(k) << ",";
            for (size_t i = 0; i < flow.nPoints(); i++) {
                fod << flow.getY(k,i) << ",";
            }
            fod << endl;
    }

    fod.close();

}
int main() {

	Combust();
   return 1;
}

