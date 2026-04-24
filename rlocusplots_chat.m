% clear; clc; close all;
clear
clc

syms J_i R_i J_o R_o K_T D_T
syms mu beta_i beta_o K_i K_o
syms alpha beta
syms g_v g_T s

F_NAIVE_GAIN = 1;

A = [  0                      , 1                                    , 0                      , 0                              ;
      -(1+mu)*R_i^2*K_T / J_i , ( -(1+mu)*R_i^2*D_T - beta_i ) / J_i , (1+mu)*R_i^2*K_T / J_i , (1+mu)*R_i^2*D_T / J_i         ;
       0                      , 0                                    , 0                      , 1                              ;
       R_o^2*K_T / J_o        , R_o^2*D_T / J_o                      ,-R_o^2*K_T / J_o        , ( -R_o^2*D_T - beta_o ) / J_o ];

K_v_i = g_v*J_i/R_i/K_i;
K_v_o = g_v*J_o/R_o/K_o;
% K_T_tf = g_T * 1 / ( R_i*K_i/J_i + R_o*K_o/J_o );

if F_NAIVE_GAIN
    B = [ 0                     , 0                     ;
          R_i*K_i / J_i * g_v , 0                     ;
          0                     , 0                     ;
          0                     , R_o*K_o / J_o * g_v ];
else
    B = [ 0                     , 0                 ;
      R_i*K_i / J_i * K_v_i , 0                     ;
      0                     , 0                     ;
      0                     , R_o*K_o / J_o * K_v_o ];
    % B = [ 0                     , 0                     , 0        ;
    %       R_i*K_i / J_i * K_v_i , 0                     , R_i*K_i / J_i * K_T_tf      ;
    %       0                     , 0                     , 0       ;
    %       0                     , R_o*K_o / J_o * K_v_o , R_o*K_o / J_o * K_T_tf         ];
end

C = [ 0   , 1     , 0    , 0    ;
      0   , 0     , 0    , 1    ;
      K_T , D_T   , -K_T , -D_T ];

Hs_den = simplify(det(s*eye(4)-A));
Hs_num = simplify(C * adjoint(s*eye(4)-A) * B);
Hs_sym = Hs_num / Hs_den;

%%% Load BOTH systems
files = {"tts_system_bot.mat", "tts_system_eot.mat"};
colors = {'b','r'}; % BOT = blue, EOT = red

% --- Precompute BOTH systems ONCE ---
H_all = cell(1,2);

for f = 1:2
    sysvals = load(files{f});

    symvec = [J_i, 
        R_i, 
        J_o, 
        R_o, 
        K_T, 
        D_T, 
        mu, 
        beta_i, 
        beta_o, 
        K_i, 
        K_o, 
        alpha, 
        beta];

    symval = [sysvals.J_i,
        sysvals.R_i, 
        sysvals.J_o, 
        sysvals.R_o, 
        sysvals.K_T, 
        sysvals.D_T, 
        sysvals.mu, 
        sysvals.beta_i, 
        sysvals.beta_o, 
        sysvals.K_i, 
        sysvals.K_o, 
        sysvals.alpha, 
        sysvals.beta];

    Hs = subs(Hs_sym, symvec, symval);
    Hs = subs(Hs, g_v, 1);
    % Hs = subs(Hs, g_T, 1);

    H_tf = tf(zeros(2));

    for i = 1:2
        for j = 1:2
            [num, den] = numden(Hs(i,j));
            num = sym2poly(expand(num));
            den = sym2poly(expand(den));
            H_tf(i,j) = tf(num, den);
        end
    end

    H_all{f} = zpk(H_tf);
end


% --- Plot ---
figure

vara = ["u_i","u_o"];
varb = ["v_i","v_o", "T"];

t = tiledlayout(2,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';

for i = 1:2
    for j = 1:2
        nexttile

        G1 = H_all{1}(i,j); % BOT
        G2 = H_all{2}(i,j); % EOT

        h = rlocusplot(G1,'r', G2,'b');

        setoptions(h, ...
            'XLim', [-20 5], ...
            'YLim', [-75 75], ...
            'Grid','off')

        hold on
        
        % BOT: J_i K_i/R_i = 0.2
        % EOT: J_i K_i/R_i = 0.2337
        if F_NAIVE_GAIN
            p1 = rlocus(G1,  40*0.2337);
            p2 = rlocus(G2,  40*0.2337);
            p11 = rlocus(G1, 20*0.2337);
            p22 = rlocus(G2, 20*0.2337);
        else
            p1 = rlocus(G1, 40);
            p2 = rlocus(G2, 40);
            p11 = rlocus(G1, 20);
            p22 = rlocus(G2, 20);
        end

        plot(real(p11), imag(p11), 'rx', 'MarkerSize',10,'LineWidth',2)
        plot(real(p22), imag(p22), 'bx', 'MarkerSize',10,'LineWidth',2)

        title(sprintf('%s \\rightarrow %s', vara(j), varb(i)), ...
            'FontWeight','bold','FontSize',16)

        xlabel('Re','FontWeight','bold','FontSize',14)
        ylabel('Im','FontWeight','bold','FontSize',14)

        ax = findall(gcf, 'Type', 'axes');
        set(ax, 'FontSize', 14, 'LineWidth', 1.2, 'FontWeight', 'bold');
    end
end

% set(gcf, 'Units', 'pixels');
% set(gcf, 'Position', [100 100 1e3 1e3]);

legend({'BOT','EOT','BOT','EOT'}, 'FontSize', 12)
