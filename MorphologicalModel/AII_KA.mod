TITLE AII A-type K channel for oscillating Vm in degenerating retina
:
: Based on 3-compartmental model from Riecke et al, 2014, J Neurophysiol
: Written by Kyle Loizos, September 2016
:

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX AIIka
	USEION k READ ek WRITE ik
	RANGE gkabar
	RANGE htau1, htau2, c
	RANGE minf, hinf, mtau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gkabar = 0.08 (mho/cm2)
	ek = -77 (mV)
	v_init = -62 (mV)
	kma = 7 (mV)
	kha = 2 (mV)
	vhalf_m = -10 (mV)
	vhalf_h = -40.5 (mV)
	mtau = 1 (ms)
	
	dt (ms)
	v (mV)
}

STATE {
	m h1 h2
}

ASSIGNED {
	ik	(mA/cm2)
	minf
	hinf
	htau1
	htau2
	c
}

INITIAL {
: The initial values were determined at a resting value of -62 mV (minf and hinf eqns)
	m = 0.00059
	h1 = 0.1791
	h2 = 0.1791
}

BREAKPOINT {
	SOLVE states
	ik = gkabar * m * (c*h1 + (1-c)*h2) * (v - ek)
}

DERIVATIVE states {
	rates(v)
	m' = (minf - m)/mtau
	h1' = (hinf - h1)/htau1
	h2' = (hinf - h2)/htau2
}

FUNCTION find_htau2(v(mV)) { LOCAL tmp_tau
	tmp_tau = 26 + (v+17)*(v+17)/4

	if (tmp_tau < 100) {
		find_htau2 = tmp_tau
	} else {
		find_htau2 = 100
	}
}

PROCEDURE rates(v(mV)) {
	c = 1/(1+exp(-(v+45)/15))
	minf = 1/(1+exp(-(v-vhalf_m)/kma))
	hinf = (1-0.83) + 0.83/(1+exp((v-vhalf_h)/kha))
	htau1 = 25 - 20/(1+exp(-(v+35)/6))
	htau2 = find_htau2(v)
}

UNITSON