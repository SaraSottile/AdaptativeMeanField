%Complete graph
c = [128, 64, 32, 16, 8, 4, 2, 1, 1/2, 1/4, 1/8, 1/16, 1/32, 1/64, 1/128];
y_infty = [0.5747, 0.5691, 0.5593, 0.5436, 0.5208, 0.4905, 0.4524, 0.4064, 0.3538, 0.2974, 0.2416, 0.1903, 0.1460, 0.1097, 0.0810];
y_peak = [0.6348, 0.6327, 0.6289, 0.6225, 0.6121, 0.5956, 0.5716, 0.5414, 0.4998, 0.4447, 0.3803, 0.3128, 0.2485, 0.2, 0.2];

figure; set(gcf,'position',[200 200 400 300]); hold all;
plot(c, y_infty, 'o-','Linewidth',2,'DisplayName','steady-state','Color', [0.4660, 0.6740, 0.1880]);
plot(c, y_peak, '<-','Linewidth',2,'DisplayName','peak','Color', [0.6350, 0.0780, 0.1840]);
box on;
set(gca, 'XScale', 'log');
legend('Interpreter','latex');
ylim([0 1]);
xlabel('parameter $c$');
ylabel('prevalence');

%Cycle graph
c = [128, 64, 32, 16, 8, 4, 2, 1, 1/2, 1/4, 1/8, 1/16, 1/32, 1/64, 1/128];
y_infty = [0.1762, 0.1918, 0.2096, 0.2282, 0.2449, 0.2567, 0.2604, 0.2541, 0.2376, 0.2134, 0.1851, 0.1561, 0.1287, 0.1041, 0.0829];
y_peak = [0.3184, 0.3261, 0.3360, 0.3475, 0.3660, 0.3950, 0.4068, 0.3995, 0.3879, 0.3653, 0.3373, 0.3134, 0.2873, 0.2593, 0.2298];

figure; set(gcf,'position',[200 200 400 300]); hold all;
plot(c, y_infty, 'o-','Linewidth',2,'DisplayName','steady-state','Color', [0.4660, 0.6740, 0.1880]);
plot(c, y_peak, '<-','Linewidth',2,'DisplayName','peak','Color', [0.6350, 0.0780, 0.1840]);
box on;
set(gca, 'XScale', 'log');
legend('Interpreter','latex');
ylim([0 1]);
xlabel('parameter $c$');
ylabel('prevalence');