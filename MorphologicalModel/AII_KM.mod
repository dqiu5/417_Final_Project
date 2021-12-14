TITLE AII M-type K channel for oscillating Vm in degenerating retina
:
: Based on 3-compartmental model from Riecke et al, 2014, J Neurophysiol
: Written by Kyle Loizos, September 2016
:

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX AIIkm
	USEION k READ ek WRITE ik
	RANGE gkmbar
	RANGE minf
	RANGE mtau	
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gkmbar = 0.03 (mho/cm2)
	ek = -77 (mV)
	v_init = -62 (mV)
	kmm  = 4 (mV)
	vhalf_m = -40 (mV)
	mtau = 50 (ms)
	
	dt (ms)
	v (mV)
}

STATE {
	m
}

ASSIGNED {
	ik	(mA/cm2)
	minf
}

INITIAL {
: The initial values were determined at a resting value of -62 mV (minf and hinf eqns)
	m = 0.0041
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = gkmbar * m * (v - ek)
}

DERIVATIVE states {
	rates(v)
	m' = (minf - m)/mtau
}

PROCEDURE rates(v(mV)) {
	minf = 1/(1+exp(-(v-vhalf_m)/kmm))
	
}

UNITSON