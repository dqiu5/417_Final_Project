TITLE AII Na channel for oscillating Vm in degenerating retina
:
: Based on 3-compartmental model from Riecke et al, 2014, J Neurophysiol
: Written by Kyle Loizos, September 2016
:

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX AIIna
	USEION na READ ena WRITE ina
	RANGE gnabar
	RANGE minf, hinf
	RANGE mtau, htau	
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gnabar = 0.2 (mho/cm2)
	ena = 50 (mV)
	v_init = -62 (mV)
	kmna  = 5 (mV)
	khna = 2 (mV)
	vhalf_m = -48 (mV)
	vhalf_h = -49.5 (mV)
	mtau = 0.01 (ms)
	htau = 0.5 (ms)
	
	dt (ms)
	v (mV)
}

STATE {
	m h
}

ASSIGNED {
	ina	(mA/cm2)
	minf
	hinf
}

INITIAL {
: The initial values were determined at a resting value of -62 mV (minf and hinf eqns)
	m = 0.0573
	h = 0.9981
}

BREAKPOINT {
	SOLVE states
	ina = gnabar * m*m*m * h * (v - ena)
}

DERIVATIVE states {
	rates(v)
	m' = (minf - m)/mtau
	h' = (hinf - h)/htau
}

PROCEDURE rates(v(mV)) {
	minf = 1/(1+exp(-(v-vhalf_m)/kmna))
	hinf = 1/(1+exp((v-vhalf_h)/khna))
	
}

UNITSON