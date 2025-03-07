function test_dysonint(inputFileName, r, integralDir, e0, ef, de)

%% TestBottom_up
% This is a script for test. 

%clear;
%clc;

%r=1;
%load('/Users/victormanuelfreixaslemus/Desktop/Projects/Photoelectron_spectroscopy/tr_XPES_code/tr-XPES/transition_moments/testdyson.mat');
load(inputFileName)
Options.BasisSet = testcase(r).Basis;

out = calcdysonint_up_rot74pes(testcase(r).xyz,testcase(r).Elements,testcase(r).TotalCharge,Options,integralDir, e0, ef, de);

%
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

end
