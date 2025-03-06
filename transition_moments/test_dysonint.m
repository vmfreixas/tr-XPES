%% TestBottom_up
% This is a script for test. 

clear;
clc;

r=1;
% load('test2017.mat');
load('/Users/victormanuelfreixaslemus/Desktop/Projects/Photoelectron_spectroscopy/tr_XPES_code/tr-XPES/transition_moments/testdyson.mat');
Options.BasisSet = testcase(r).Basis;

% out = calcBottom_up(testcase(r).xyz,testcase(r).Elements,testcase(r).TotalCharge,Options);
% out = calcdysonint_up(testcase(r).xyz,testcase(r).Elements,testcase(r).TotalCharge,Options);
%out = calcdysonint_up_rot(testcase(r).xyz,testcase(r).Elements,testcase(r).TotalCharge,Options);
%out = calcdysonint_up_rot3(testcase(r).xyz,testcase(r).Elements,testcase(r).TotalCharge,Options);
%out = calcdysonint_up_rot4(testcase(r).xyz,testcase(r).Elements,testcase(r).TotalCharge,Options);
%out = calcdysonint_up_rot74(testcase(r).xyz,testcase(r).Elements,testcase(r).TotalCharge,Options);
out = calcdysonint_up_rot74pes(testcase(r).xyz,testcase(r).Elements,testcase(r).TotalCharge,Options);
%out = calcdysonint_up_rot74part2(testcase(r).xyz,testcase(r).Elements,testcase(r).TotalCharge,Options);
%out = calcdysonint_up_rot74part3(testcase(r).xyz,testcase(r).Elements,testcase(r).TotalCharge,Options);
%out = calcdysonint_up_rot38part2(testcase(r).xyz,testcase(r).Elements,testcase(r).TotalCharge,Options);
%out = calcdysonint_up_rot38part3(testcase(r).xyz,testcase(r).Elements,testcase(r).TotalCharge,Options);
%out = calcdysonint_up_rot38part4(testcase(r).xyz,testcase(r).Elements,testcase(r).TotalCharge,Options);
%out = calcdysonint_up_rotR(testcase(r).xyz,testcase(r).Elements,testcase(r).TotalCharge,Options);

% this part is pyrrole 
% writematrix(real(out.dysonx),'dysonorbint/kzdysonxreal.txt');
% writematrix(imag(out.dysonx),'dysonorbint/kzdysonximag.txt');
% 
% writematrix(real(out.dysony),'dysonorbint/kzdysonyreal.txt');
% writematrix(imag(out.dysony),'dysonorbint/kzdysonyimag.txt');
% 
% writematrix(real(out.dysonz),'dysonorbint/kzdysonzreal.txt');
% writematrix(imag(out.dysonz),'dysonorbint/kzdysonzimag.txt');

% writematrix(real(out.dysonx),'dysonorbint/kyzdysonxreal.txt');
% writematrix(imag(out.dysonx),'dysonorbint/kyzdysonximag.txt');
% 
% writematrix(real(out.dysony),'dysonorbint/kyzdysonyreal.txt');
% writematrix(imag(out.dysony),'dysonorbint/kyzdysonyimag.txt');
% 
% writematrix(real(out.dysonz),'dysonorbint/kyzdysonzreal.txt');
% writematrix(imag(out.dysonz),'dysonorbint/kyzdysonzimag.txt');
% 
% writematrix(real(out.dysonx),'dysonorbint/kxydysonxreal.txt');
% writematrix(imag(out.dysonx),'dysonorbint/kxydysonximag.txt');
% 
% writematrix(real(out.dysony),'dysonorbint/kxydysonyreal.txt');
% writematrix(imag(out.dysony),'dysonorbint/kxydysonyimag.txt');
% 
% writematrix(real(out.dysonz),'dysonorbint/kxydysonzreal.txt');
% writematrix(imag(out.dysonz),'dysonorbint/kxydysonzimag.txt');


%this part is MTSO turn this part off for rotational average. turn this on
%for oriented
% writematrix(real(out.dysonx),'dysonorbint/MTSOdysonxreal.txt');
% writematrix(imag(out.dysonx),'dysonorbint/MTSOdysonximag.txt');
% 
% writematrix(real(out.dysony),'dysonorbint/MTSOdysonyreal.txt');
% writematrix(imag(out.dysony),'dysonorbint/MTSOdysonyimag.txt');
% 
% writematrix(real(out.dysonz),'dysonorbint/MTSOdysonzreal.txt');
% writematrix(imag(out.dysonz),'dysonorbint/MTSOdysonzimag.txt');


