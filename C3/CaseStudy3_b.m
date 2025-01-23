clearvars;set(0,'defaulttextinterpreter','latex'); format compact;% close all;%clc
rng(2); % setting seed=1 for reproducibility

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
n = 50;                   % number of groups
A =  load('graph4.mat');% Adjacency matrix
A = table2array(struct2table(A)) + eye(n);% add diagonal to A
beta = 0.8*rand(n,n);     % Infection rate \beta_{ij}
delta = 1*ones(n,1);      % Curing rate \delta_i
zeta = 2*rand(n,n);       % Link-breaking rate \zeta_{ij}
xi = 1*rand(n,n);         % Link-creation rate \xi_{ij}
vector_c = linspace(0, 1, 11);

% initial prevalence
y_init = zeros(n,1); y_init(1) = 0.2;   % initial prevalence
z_init = 1*zeros(n,n);                   % initial links

dt = 0.01;    % Numerical integration time step
tmax = 50;    % Time range

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise simulations
vector_y_peak = zeros(length(vector_c), 1);
vector_y_ss = zeros(length(vector_c), 1);
vector_t = 0:dt:tmax;
Nt = length(vector_t);
beta_matrix = beta .* A;

for k=1:length(vector_c)

    % Functional responses (link-breaking fbr and link-creation fcr)
    % DONT FORGET: ADD ELEMENTWISE OPERATOR .
        c = vector_c(k)
    fbr_in = @(y,yglobal) (y.^2);   % Link-breaking within a community
    fcr_in = @(y,yglobal) (1-y.*y);          % Link-creation within a community
    fbr_out = @(y1,y2,yglobal) c*y1.*y2 + (1-c)*yglobal.^2;   % Link-breaking between communities
    fcr_out = @(y1,y2,yglobal) (1-y1.*y2);   % Link-creation between communities
    
    % Initialise simulations
    tic;
    R0 = calculateR0(n, delta, beta_matrix, zeta, xi, fbr_in, fcr_in, fbr_out, fcr_out)
    y = zeros(Nt, n);
    z = zeros(Nt, n, n);
    y(1, :) = y_init;
    z(1, :, :) = z_init;

    % Perform simulations (simple Forward Euler)
    for t=2:Nt
        for i=1:n
            y(t,i) = y(t-1,i) + dt * f1(i, y(t-1,:), squeeze(z(t-1,:,:)), delta, beta_matrix);
            for j=1:n
                if beta_matrix(i,j)>0 || i==j
                    z(t,i,j) = z(t-1,i,j) + dt * f2(i, j, y(t-1,:), squeeze(z(t-1,:,:)), zeta, xi, fbr_in, fcr_in, fbr_out, fcr_out);
                end
            end
        end
    end
    
    % Compute average prevalence and link density
    total_link_weight = sum(A, 'all');
    y_average = mean(y, 2);
    z_average = sum(z, [2,3]) / total_link_weight;
    z_in_average = zeros(Nt, 1);
    z_out_average = zeros(Nt, 1);
    for t = 1:Nt
        d = squeeze(z(t,:,:));
        z_in_average(t) = sum(diag(d)) / n;
        z_out_average(t) = sum(d - diag(diag(d)), 'all') / (total_link_weight - n);
    end
    y_ss = y_average(end)
    y_peak = max(y, [], 'all')
    
    vector_y_ss(k) = y_ss;
    vector_y_peak(k) = y_peak;

end

%% Plot results
% Prevalence over time
figure; set(gcf,'position',[200 200 400 300]); hold all;
plot(vector_t, y(:,1), 'b', 'DisplayName','individual prevalences $y_i$');
for i=2:n
    plot(vector_t, y(:,2:end), 'b','HandleVisibility','off');
end
plot(vector_t, y_average, 'k','Linewidth',2,'DisplayName','average prevalence $\bar{y}$');
legend('Interpreter','Latex');
box on;
ylim([0 1]);
xlabel('time $t$');
ylabel('prevalence');

% Link density over time
figure; set(gcf,'position',[200 200 400 300]); hold all;
plot(vector_t, z(:,1,1), 'r', 'DisplayName','internal link densities $z_{ii}$');
plot(vector_t, z(:,1,1), 'b', 'DisplayName','external link densities $z_{ij}$');
for i=1:n
    for j=1:n
        if i~=j
            plot(vector_t, z(:,i,j), 'b','HandleVisibility','off');
        else
            plot(vector_t, z(:,i,j), 'r','HandleVisibility','off');
        end
    end
end
plot(vector_t, z_average, 'k','Linewidth',2,'DisplayName','average link density');
legend('Interpreter','Latex');
box on;
ylim([0 1]);
xlabel('time $t$');
ylabel('link density');

% Link densities
figure; set(gcf,'position',[200 200 400 300]); hold all;
plot(vector_t, y_average, 'Linewidth',2,'DisplayName','prevalence');
plot(vector_t, z_average, 'Linewidth',2,'DisplayName','average link density');
plot(vector_t, z_in_average, '--', 'Linewidth',2,'DisplayName','average IN link density','Color', [0.8500, 0.3250, 0.0980]);
plot(vector_t, z_out_average, ':', 'Linewidth',2,'DisplayName','average OUT link density','Color', [0.8500, 0.3250, 0.0980]);
% legend('Interpreter','Latex');
box on;
ylim([0 1]);
xlabel('time $t$');
ylabel('fraction');

% Plot solutions in (z,y) plane
y0 = y_average(1);
z0 = z_average(1);
yend = y_average(end);
zend = z_average(end);
figure; set(gcf,'position',[200 200 400 300]); hold all;
plot(z_average, y_average, 'b');    % Show average trajectory
plot(z0, y0, 'b*');                 % show average starting point
plot(zend, yend, 'k.', 'Markersize', 22); % Show average end point
% plot(1, 0, 'r.','Markersize',22);   % Show DFE (TO BE DONE!)
hold off;
xlim([0 1]);
ylim([0 1]);
box on;
xlabel('average link density $\bar{z}$');
ylabel('average prevalence $\bar{y}$');

% Plot for various c values
figure; set(gcf,'position',[200 200 400 300]); hold all;
plot(vector_c, vector_y_peak, '.-', 'DisplayName', 'peak');
plot(vector_c, vector_y_ss, '.-', 'DisplayName','steady-state');
legend('Interpreter','Latex');
xlim([0, 1]);
ylim([0, 1]);
box on;
xlabel('c');
ylabel('prevalence');

% Peak prevalence (bordeaux red)
figure; set(gcf,'position',[200 200 400 300]); hold all;
plot(vector_c, vector_y_peak, '.-', 'DisplayName','peak', 'Color', [0.6350, 0.0780, 0.1840]);
hold off;
xlim([0, 1]);
box on;
xlabel('c');
ylabel('prevalence');

% Steady-state prevalence (green)
figure; set(gcf,'position',[200 200 400 300]); hold all;
plot(vector_c, vector_y_ss, '.-', 'DisplayName','steady-state', 'Color', [0.4660, 0.6740, 0.1880]);
hold off;
xlim([0, 1]);
box on;
xlabel('c');
ylabel('prevalence');

toc;


% NIMFA equation (1a)
function output = f1(i, y, z, delta, beta)
    output = - delta(i)*y(i) + (1-y(i))*(beta(i,:).*z(i,:))*y';
end
% Network-changing equation (1b)
function output = f2(i, j, y, z, zeta, xi, fbr_in, fcr_in, fbr_out, fcr_out)
    if (i==j)
        output = - zeta(i,j) * z(i,j) * fbr_in(y(i), mean(y)) + xi(i,j) * (1-z(i,j)) * fcr_in(y(i), mean(y));
    else
        output = - zeta(i,j) * z(i,j) * fbr_out(y(i), y(j), mean(y)) + xi(i,j) * (1-z(i,j)) * fcr_out(y(i), y(j), mean(y));
    end
end
% Calculate the basic reproduction number R0 based on Eq. (3), (4)
function R0 = calculateR0(n, delta, beta, zeta, xi, fbr_in, fcr_in, fbr_out, fcr_out)
    z_DFE = zeros(n, n);
    for i=1:n
        for j=1:n
            if (i==j)
                z_DFE(i,j) = xi(i,j) * fcr_in(0,0) / (zeta(i,j) * fbr_in(0,0) + xi(i,j)* fcr_in(0,0));
            else
                z_DFE(i,j) = xi(i,j) * fcr_out(0,0,0) / (zeta(i,j) * fbr_out(0,0,0) + xi(i,j)* fcr_out(0,0,0));
            end
        end
    end
    M = z_DFE .* beta;
    V_inv = diag(1./ delta);
    F = M * V_inv;
    R0 = max(eig(F), [], 'ComparisonMethod', 'real');
end