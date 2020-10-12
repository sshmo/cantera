//! @file PorousFlow.h


#ifndef CT_POROUSFLOW_H
#define CT_POROUSFLOW_H

#include "cantera/oneD/StFlow.h"

namespace Cantera
{
    /**
 * A class for flow through porous material.
 * @ingroup onedim
 */

class PorousFlow : public StFlow
{
public:
    PorousFlow(ThermoPhase* ph = 0, size_t nsp = 1, size_t points = 1);

    //! Delegating constructor
    PorousFlow(shared_ptr<ThermoPhase> th, size_t nsp = 1, size_t points = 1) :
    PorousFlow(th.get(), nsp, points) {
    }

    virtual void setupGrid(size_t n, const doublereal* z);
    //! initialize the solid solver as well as the radiant flux vector
    virtual void eval(size_t j, doublereal* x, doublereal* r,
                      integer* mask, doublereal rdt);
    void solid(doublereal* x, vector_fp & hconv, vector_fp & scond,
               vector_fp & RK, vector_fp & Omega, double & srho, 
               double & sCp, double rdt);
    virtual XML_Node& save(XML_Node& o, const doublereal* const sol);
	
    virtual void restore(const XML_Node& dom, doublereal* soln,
                         int loglevel);
						 
    //! initialize the solid properties
    double pore1;
    double pore2;
    double diam1;
    double diam2;
    double scond1;
    double scond2;
    double Omega1;
    double Omega2;
    double srho;
    double sCp;
    double m_zmid;
    double m_dzmid;
    double m_porea;
    double m_poreb;
    double m_porec;
    double m_pored;
    double m_diama;
    double m_diamb;
    double m_diamc;
    double m_diamd;


    int geometry; 
    vector_fp dq;
    doublereal getTw  (const int & i) { return Tw[i];   }
    doublereal getDq  (const int & i) { return dq[i];   }
    doublereal getPore(const int & i) { return pore[i]; }
    doublereal getDiam(const int & i) { return diam[i]; }
    doublereal getScond(const int & i) { return scond[i]; }
    doublereal getHconv(const int & i) {return hconv[i]; } 
    
    doublereal getZ (const int & i) { return m_z[i];   }
    doublereal getT (const int & i) { return T_prev(i);   }
    doublereal getY (const int & k, const int & i) { return Y_prev(k, i);   }

    virtual std::string flowType() {
        return "Porous Stagnation";
    }
private:
    // porous burner
    vector_fp Tw;
    vector_fp pore;
    vector_fp diam;
    vector_fp scond;
    vector_fp Twprev;
    vector_fp Twprev1;
    vector_fp zprev;
    vector_fp hconv;
    int m_adapt;
};

}

#endif