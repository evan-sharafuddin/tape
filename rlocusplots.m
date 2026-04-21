
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
syms mu beta_i beta_o K_i K_o

% other parameters
syms alpha beta

A = [  0                      , 1                                    , 0                      , 0                              ; 
      -(1+mu)*R_i^2*K_T / J_i , ( -(1+mu)*R_i^2*D_T - beta_i ) / J_i , (1+mu)*R_i^2*K_T / J_i , (1+mu)*R_i^2*D_T / J_i         ;
       0                      , 0                                    , 0                      , 1                              ;
       R_o^2*K_T / J_o        , R_o^2*D_T / J_o                      ,-R_o^2*K_T / J_o        , ( -R_o^2*D_T - beta_o ) / J_o ];

syms g_v
syms s
K_v_i = g_v*J_i/R_i/K_i;
K_v_o = g_v*J_o/R_o/K_o;

B = [ 0                     , 0                     ;
      R_i*K_i / J_i * K_v_i , 0                     ;
      0                     , 0                     ;
      0                     , R_o*K_o / J_o * K_v_o ];

C = [ 0   , 1     , 0    , 0    ;
      0   , 0     , 0    , 1    ;
      K_T , D_T   , -K_T , -D_T ];

Hs_den = simplify( det(s*eye(4)-A) );
Hs_num = simplify( C * adjoint( s*eye(4)-A ) * B );
Hs = Hs_num / Hs_den;

sysvals = load("tts_system.mat");

% symbolic variable names (as symbolic variables, not strings)
symvec = [J_i, R_i, J_o, R_o, K_T, D_T, ...
          mu, beta_i, beta_o, K_i, K_o, ...
          alpha, beta];

% corresponding numerical values from struct
symval = [sysvals.J_i, sysvals.R_i, sysvals.J_o, sysvals.R_o, ...
          sysvals.K_T, sysvals.D_T, ...
          sysvals.mu, sysvals.beta_i, sysvals.beta_o, ...
          sysvals.K_i, sysvals.K_o, ...
          sysvals.alpha, sysvals.beta];

Hs = subs(Hs, symvec, symval);

% when look at analytic solution, it is clear that g_v scales each term in
% the matrix. Therefore, we can just subs in 1 for it and use rlocfind to
% determine where g_v=40 gets us

Hs = subs(Hs, g_v, 1);

syms s

% assume Hs is already defined as your 2x2 symbolic matrix

% preallocate
H_tf = tf(zeros(2));

for i = 1:2
    for j = 1:2
        % get numerator and denominator
        [num, den] = numden(Hs(i,j));
        
        % expand to polynomials in s
        num = expand(num);
        den = expand(den);
        
        % convert to coefficient vectors
        num_coeff = sym2poly(num);
        den_coeff = sym2poly(den);
        
        % build transfer function
        H_tf(i,j) = tf(num_coeff, den_coeff);
    end
end

H_tf = zpk(H_tf);

figure

for i = 1:2
    for j = 1:2
        subplot(2,2,(i-1)*2 + j)
        
        G = H_tf(i,j);
        
        % Root locus
        rlocus(G)
        hold on
        
        % Poles at K = 40
        p = rlocus(G, 40);
        plot(real(p), imag(p), 'rx', 'MarkerSize', 10, 'LineWidth', 2)
        
        title(sprintf('H(%d,%d)', i, j))
        % grid on
    end
end
% figure, rlocfind(H_tf(1,1))
% figure, rlocfind(H_tf(2,1))
% figure, rlocfind(H_tf(1,2))
% figure, rlocfind(H_tf(2,2))
