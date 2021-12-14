clc;
clear;
close all;

%Project Plotting Voltage vs Time

data = importdata('trial_custom_stimulus4.txt');
Vext = data(:,2);
xpos = data(:,1);
plot(xpos, Vext);
