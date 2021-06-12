x=[2,4,6,8,10];
hold on;
plot(x,T,'sqk')
plot(x,58.38.*x-63.09,'k')

legend('Average Time','58.38$\times$NVP - 63.09','location', 'northwest')

annotation('textbox', [0.2, 0.7, 0.1, 0.1], 'String', "$R^2=0.9909$")

xlabel("Number of Variational Parameters (NVP)")
ylabel("Time [s]")