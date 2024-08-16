   % subjID: e.g. S001
clear all;
clc;

SubjID   = 'slh';
sess_num = 2;
run_num  = 5;         

%%  
offset = [0,0];
CurrDir = pwd;
SetupRand;
% set_test_gamma;
HideCursor; 
parameters;
 
%%
% EyelinkExample;
% sample;
% prac;
test;
% training;