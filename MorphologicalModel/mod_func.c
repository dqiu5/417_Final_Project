#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _A2hh_k_reg();
extern void _A2hh_na_reg();
extern void _AII_KA_reg();
extern void _AII_KM_reg();
extern void _AII_Na_reg();
extern void _IinjLT_reg();
extern void _ipulse1_reg();
extern void _ipulse2_reg();
extern void _ipulse3_reg();

void modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," A2hh_k.mod");
fprintf(stderr," A2hh_na.mod");
fprintf(stderr," AII_KA.mod");
fprintf(stderr," AII_KM.mod");
fprintf(stderr," AII_Na.mod");
fprintf(stderr," IinjLT.mod");
fprintf(stderr," ipulse1.mod");
fprintf(stderr," ipulse2.mod");
fprintf(stderr," ipulse3.mod");
fprintf(stderr, "\n");
    }
_A2hh_k_reg();
_A2hh_na_reg();
_AII_KA_reg();
_AII_KM_reg();
_AII_Na_reg();
_IinjLT_reg();
_ipulse1_reg();
_ipulse2_reg();
_ipulse3_reg();
}
