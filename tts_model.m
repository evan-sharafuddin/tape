
clear
clc
close all

%%% Generate the system dynamics for the linear tape transport system
%   dynamics

%%% Create tape parameters
% i -- takeup/machine reel
% o -- supply/file reel

% time/LPOS varying parameters
syms J_i R_i J_o R_o K_T D_T 

% constant parameters
syms mu beta_i beta_o

% other parameters
syms alpha beta

A = [  0                      , 1                                    , 0                      , 0                         ; 
      -(1+mu)*R_i^2*K_T / J_i , ( -(1+mu)*R_i^2*D_T - beta_i ) / J_i , (1+mu)*R_i^2*K_T / J_i , (1+mu)*R_i^2*D_T / J_i    ;
       0                      , 0                                    , 0                      , 1                         ;
       R_o^2*K_T / J_o        , R_o^2*D_T / J_o                      ,-R_o^2*K_T / J_o        , ( -R_o^2 - beta_o ) / J_o ];