%% Sets latex as default text interpreter
set(groot, 'DefaultAxesFontSize', 12);
set(groot, 'DefaultTextInterpreter',          'latex');
set(groot, 'DefaultLegendInterpreter',        'latex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');

%% Figure 3 a) and b)
load("./Resultados/m_0_kerr.mat")
analytics_a = 2.*sqrt(lambda_c);
analytics_b = sqrt(2).*lambda_c.^(-1/2);

figure('Name','Figure 3a','NumberTitle','off')
    hold on;
    An = mean(An);
    A = An(:,1,:);B = An(:,2,:);
    A = permute(A, [3 2 1]);
    B = permute(B, [3 2 1]);
    plot(lambda_c,A,"rd")
    plot(lambda_c,analytics_a,"k")
    plot(lambda_c,B,"b*")
    plot(lambda_c,analytics_b,"k--")
    xlabel("$\lambda$");
    legend("CS$-a$","$2\sqrt{2}$",...
           "CS$-b$","$\sqrt{2}\lambda^{-1/2}$",...
           'Position',[0.4 0.7 0.1 0.2])
clear
figure('Name','Figure 3b','NumberTitle','off')
load("./Resultados/m_0_kerr_sech.mat")
analytics_a = 2.17.*sqrt(lambda_c);
analytics_b = 0.78.*lambda_c.^(-1/2);
    hold on;
    An = mean(An);
    A = An(:,1,:);B = An(:,2,:);
    A = permute(A, [3 2 1]);
    B = permute(B, [3 2 1]);
    plot(lambda_c,A,"rd")
    plot(lambda_c,analytics_a,"k")
    plot(lambda_c,B,"b*")
    plot(lambda_c,analytics_b,"k--")
    xlabel("$\lambda$");
    legend("CS$-a$","$2\sqrt{\lambda}$",...
           "CS$-b$","$\sqrt{2}\lambda^{-1/2}$",...
           'Position',[0.4 0.7 0.1 0.2])
%% figure 4 a) and b)  
load("./Resultados/m_1_sat_m_fig_2_s_05.mat")
analytics_a = 2.308.*lambda_c+0.005773;
analytics_b = 1.875.*lambda_c.^(-1/2)+0.2207;

figure('Name','Figure 4a','NumberTitle','off')
    hold on;
    A = z(:,1);
    B = z(:,2);
    plot(lambda_c,A,"rd")
    plot(lambda_c,analytics_a,"k")
    plot(lambda_c,B,"b*")
    plot(lambda_c,analytics_b,"k--")
    xlabel("$\lambda$");
    legend("CS$-a$","$2.308\sqrt{\lambda}+0.005773$",...
           "CS$-b$","$1.875\lambda^{-1/2}+0.2207$",...
           'Position',[0.4 0.7 0.1 0.2])
f
clear
load("./Resultados/m_1_sat_m_fig_2_s_2.mat")
analytics_a = 2.460.*lambda_c-0.08048;
analytics_b = 1.8075.*lambda_c.^(-1/2)+0.1525.*lambda_c+0.4424;

figure('Name','Figure 4b','NumberTitle','off')
    hold on;
    A = z(:,1);
    B = z(:,2);
    plot(lambda_c,A,"rd")
    plot(lambda_c,analytics_a,"k")
    plot(lambda_c,B,"b*")
    plot(lambda_c,analytics_b,"k--")
    xlabel("$\lambda$");
    legend("CS$-a$","$2.460\sqrt{\lambda}-0.08048$",...
           "CS$-b$","$1.8075\lambda^{-1/2}+0.1525\lambda+0.4423$",...
           'Position',[0.4 0.7 0.1 0.2])
       