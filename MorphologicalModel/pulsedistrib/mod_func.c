#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _ipulse1_reg();
extern void _ipulse2_reg();
extern void _ipulse3_reg();

void modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," ipulse1.mod");
fprintf(stderr," ipulse2.mod");
fprintf(stderr," ipulse3.mod");
fprintf(stderr, "\n");
    }
_ipulse1_reg();
_ipulse2_reg();
_ipulse3_reg();
}
