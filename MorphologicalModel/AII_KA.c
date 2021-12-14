/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__AIIka
#define _nrn_initial _nrn_initial__AIIka
#define nrn_cur _nrn_cur__AIIka
#define _nrn_current _nrn_current__AIIka
#define nrn_jacob _nrn_jacob__AIIka
#define nrn_state _nrn_state__AIIka
#define _net_receive _net_receive__AIIka 
#define rates rates__AIIka 
#define states states__AIIka 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define gkabar _p[0]
#define mtau _p[1]
#define minf _p[2]
#define hinf _p[3]
#define htau1 _p[4]
#define htau2 _p[5]
#define c _p[6]
#define m _p[7]
#define h1 _p[8]
#define h2 _p[9]
#define ek _p[10]
#define Dm _p[11]
#define Dh1 _p[12]
#define Dh2 _p[13]
#define ik _p[14]
#define _g _p[15]
#define _ion_ek	*_ppvar[0]._pval
#define _ion_ik	*_ppvar[1]._pval
#define _ion_dikdv	*_ppvar[2]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 /* external NEURON variables */
 /* declaration of user functions */
 static void _hoc_find_htau2(void);
 static void _hoc_rates(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_AIIka", _hoc_setdata,
 "find_htau2_AIIka", _hoc_find_htau2,
 "rates_AIIka", _hoc_rates,
 0, 0
};
#define find_htau2 find_htau2_AIIka
 extern double find_htau2( double );
 /* declare global and static user variables */
#define kha kha_AIIka
 double kha = 2;
#define kma kma_AIIka
 double kma = 7;
#define vhalf_h vhalf_h_AIIka
 double vhalf_h = -40.5;
#define vhalf_m vhalf_m_AIIka
 double vhalf_m = -10;
#define v_init v_init_AIIka
 double v_init = -62;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "v_init_AIIka", "mV",
 "kma_AIIka", "mV",
 "kha_AIIka", "mV",
 "vhalf_m_AIIka", "mV",
 "vhalf_h_AIIka", "mV",
 "gkabar_AIIka", "mho/cm2",
 "mtau_AIIka", "ms",
 0,0
};
 static double delta_t = 1;
 static double h20 = 0;
 static double h10 = 0;
 static double maxerr = 1e-05;
 static double m0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "v_init_AIIka", &v_init_AIIka,
 "kma_AIIka", &kma_AIIka,
 "kha_AIIka", &kha_AIIka,
 "vhalf_m_AIIka", &vhalf_m_AIIka,
 "vhalf_h_AIIka", &vhalf_h_AIIka,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[3]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"AIIka",
 "gkabar_AIIka",
 "mtau_AIIka",
 0,
 "minf_AIIka",
 "hinf_AIIka",
 "htau1_AIIka",
 "htau2_AIIka",
 "c_AIIka",
 0,
 "m_AIIka",
 "h1_AIIka",
 "h2_AIIka",
 0,
 0};
 static Symbol* _k_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 16, _prop);
 	/*initialize range parameters*/
 	gkabar = 0.08;
 	mtau = 1;
 	_prop->param = _p;
 	_prop->param_size = 16;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_k_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ek */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dikdv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _AII_KA_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("k", -10000.);
 	_k_sym = hoc_lookup("k_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 16, 4);
  hoc_register_dparam_semantics(_mechtype, 0, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 AIIka C:/Users/daoda/Documents/EECS417/Project/Final_Code/AII_KA.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "AII A-type K channel for oscillating Vm in degenerating retina";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rates(double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static double *_temp1;
 static int _slist1[3], _dlist1[3];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   rates ( _threadargscomma_ v ) ;
   Dm = ( minf - m ) / mtau ;
   Dh1 = ( hinf - h1 ) / htau1 ;
   Dh2 = ( hinf - h2 ) / htau2 ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 rates ( _threadargscomma_ v ) ;
 Dm = Dm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / mtau )) ;
 Dh1 = Dh1  / (1. - dt*( ( ( ( - 1.0 ) ) ) / htau1 )) ;
 Dh2 = Dh2  / (1. - dt*( ( ( ( - 1.0 ) ) ) / htau2 )) ;
  return 0;
}
 /*END CVODE*/
 
static int states () {_reset=0;
 {
   rates ( _threadargscomma_ v ) ;
   Dm = ( minf - m ) / mtau ;
   Dh1 = ( hinf - h1 ) / htau1 ;
   Dh2 = ( hinf - h2 ) / htau2 ;
   }
 return _reset;}
 
double find_htau2 (  double _lv ) {
   double _lfind_htau2;
 double _ltmp_tau ;
 _ltmp_tau = 26.0 + ( _lv + 17.0 ) * ( _lv + 17.0 ) / 4.0 ;
   if ( _ltmp_tau < 100.0 ) {
     _lfind_htau2 = _ltmp_tau ;
     }
   else {
     _lfind_htau2 = 100.0 ;
     }
   
return _lfind_htau2;
 }
 
static void _hoc_find_htau2(void) {
  double _r;
   _r =  find_htau2 (  *getarg(1) );
 hoc_retpushx(_r);
}
 
static int  rates (  double _lv ) {
   c = 1.0 / ( 1.0 + exp ( - ( _lv + 45.0 ) / 15.0 ) ) ;
   minf = 1.0 / ( 1.0 + exp ( - ( _lv - vhalf_m ) / kma ) ) ;
   hinf = ( 1.0 - 0.83 ) + 0.83 / ( 1.0 + exp ( ( _lv - vhalf_h ) / kha ) ) ;
   htau1 = 25.0 - 20.0 / ( 1.0 + exp ( - ( _lv + 35.0 ) / 6.0 ) ) ;
   htau2 = find_htau2 ( _threadargscomma_ _lv ) ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
   _r = 1.;
 rates (  *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 3;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ek = _ion_ek;
     _ode_spec1 ();
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 3; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 ();
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ek = _ion_ek;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_k_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_k_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_k_sym, _ppvar, 2, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  h2 = h20;
  h1 = h10;
  m = m0;
 {
   m = 0.00059 ;
   h1 = 0.1791 ;
   h2 = 0.1791 ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  ek = _ion_ek;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   ik = gkabar * m * ( c * h1 + ( 1.0 - c ) * h2 ) * ( v - ek ) ;
   }
 _current += ik;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
  ek = _ion_ek;
 _g = _nrn_current(_v + .001);
 	{ double _dik;
  _dik = ik;
 _rhs = _nrn_current(_v);
  _ion_dikdv += (_dik - ik)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ik += ik ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
double _dtsav = dt;
if (secondorder) { dt *= 0.5; }
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
  ek = _ion_ek;
 { error =  adrunge(_ninits, 3, _slist1, _dlist1, _p, &t, dt, states, &_temp1, maxerr);
 if(error){fprintf(stderr,"at line 57 in file AII_KA.mod:\n	SOLVE states\n"); nrn_complain(_p); abort_run(error);}
    if (secondorder) {
    int _i;
    for (_i = 0; _i < 3; ++_i) {
      _p[_slist1[_i]] += dt*_p[_dlist1[_i]];
    }}
  states();
 } }}
 dt = _dtsav;
}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(m) - _p;  _dlist1[0] = &(Dm) - _p;
 _slist1[1] = &(h1) - _p;  _dlist1[1] = &(Dh1) - _p;
 _slist1[2] = &(h2) - _p;  _dlist1[2] = &(Dh2) - _p;
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "AII_KA.mod";
static const char* nmodl_file_text = 
  "TITLE AII A-type K channel for oscillating Vm in degenerating retina\n"
  ":\n"
  ": Based on 3-compartmental model from Riecke et al, 2014, J Neurophysiol\n"
  ": Written by Kyle Loizos, September 2016\n"
  ":\n"
  "\n"
  "INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}\n"
  "\n"
  "NEURON {\n"
  "	SUFFIX AIIka\n"
  "	USEION k READ ek WRITE ik\n"
  "	RANGE gkabar\n"
  "	RANGE htau1, htau2, c\n"
  "	RANGE minf, hinf, mtau\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(mA) = (milliamp)\n"
  "	(mV) = (millivolt)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	gkabar = 0.08 (mho/cm2)\n"
  "	ek = -77 (mV)\n"
  "	v_init = -62 (mV)\n"
  "	kma = 7 (mV)\n"
  "	kha = 2 (mV)\n"
  "	vhalf_m = -10 (mV)\n"
  "	vhalf_h = -40.5 (mV)\n"
  "	mtau = 1 (ms)\n"
  "	\n"
  "	dt (ms)\n"
  "	v (mV)\n"
  "}\n"
  "\n"
  "STATE {\n"
  "	m h1 h2\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	ik	(mA/cm2)\n"
  "	minf\n"
  "	hinf\n"
  "	htau1\n"
  "	htau2\n"
  "	c\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  ": The initial values were determined at a resting value of -62 mV (minf and hinf eqns)\n"
  "	m = 0.00059\n"
  "	h1 = 0.1791\n"
  "	h2 = 0.1791\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE states\n"
  "	ik = gkabar * m * (c*h1 + (1-c)*h2) * (v - ek)\n"
  "}\n"
  "\n"
  "DERIVATIVE states {\n"
  "	rates(v)\n"
  "	m' = (minf - m)/mtau\n"
  "	h1' = (hinf - h1)/htau1\n"
  "	h2' = (hinf - h2)/htau2\n"
  "}\n"
  "\n"
  "FUNCTION find_htau2(v(mV)) { LOCAL tmp_tau\n"
  "	tmp_tau = 26 + (v+17)*(v+17)/4\n"
  "\n"
  "	if (tmp_tau < 100) {\n"
  "		find_htau2 = tmp_tau\n"
  "	} else {\n"
  "		find_htau2 = 100\n"
  "	}\n"
  "}\n"
  "\n"
  "PROCEDURE rates(v(mV)) {\n"
  "	c = 1/(1+exp(-(v+45)/15))\n"
  "	minf = 1/(1+exp(-(v-vhalf_m)/kma))\n"
  "	hinf = (1-0.83) + 0.83/(1+exp((v-vhalf_h)/kha))\n"
  "	htau1 = 25 - 20/(1+exp(-(v+35)/6))\n"
  "	htau2 = find_htau2(v)\n"
  "}\n"
  "\n"
  "UNITSON\n"
  ;
#endif
