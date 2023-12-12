%% Load variables
load('Decision_Rules/ap11.mat')
load('Decision_Rules/ap12.mat')
load('Decision_Rules/c11.mat')
load('Decision_Rules/c12.mat')
load('Decision_Rules/faz11.mat')
load('Decision_Rules/faz12.mat')
load('Decision_Rules/p11.mat')
load('Decision_Rules/p12.mat')
load('Decision_Rules/v11.mat')
load('Decision_Rules/v12.mat')
load('Decision_Rules/agrid.mat')

%%  Other Variables of interest
amax = 5; 
amin = 0;
I = 1000;              
J = 2;  

%% Total Plot Section

% Value Functions
figure(1)

subplot(2,2,[3, 4])
hold on
plot(agrid, v11(:,1) ,'LineWidth', 2, 'color', "#ff9696")
plot(agrid, v11(:,2) ,'LineWidth', 2, 'color', "#FF0000")
plot(agrid, v12(:,1)  ,'LineWidth', 2, 'color', "#8293ff")
plot(agrid, v12(:,2) ,'LineWidth', 2, 'color', "#0023ff")
xlabel('$a_t$', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('$v^j(a_t,\varepsilon_t)$', 'Interpreter', 'latex', 'FontSize',16)
hl = legend('$J=1$ \, $\varepsilon_t = \varepsilon^L$', '$J=1$ \, $\varepsilon_t =\varepsilon^H$', '$J = 2$ \, $\varepsilon_t =\varepsilon^L$', '$J = 2$ \, $\varepsilon_t =\varepsilon^L$');
set(hl,'Interpreter','latex','Location','southeast')

ax = gca;
set(gca,'TickLabelInterpreter', 'latex');
    set(gca,...
            'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
         'FontName','cmr10',...
        'FontSize',14,...
        'Box','off'); 
hold off



% Consumption Functions
subplot(2,2,1)
hold on
plot(agrid, c11(:,1),'LineWidth',2,'color',"#ff9696")
plot(agrid, c11(:,2),'LineWidth',2,'color',"#FF0000")
plot(agrid, c12(:,1),'LineWidth',2,'color',"#8293ff")
plot(agrid, c12(:,2),'LineWidth',2,'color',"#0023ff")
xlabel('$a_t$', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('${c^j}_t(a_t,\varepsilon_t)$', 'Interpreter', 'latex', 'FontSize',16)
ax = gca;
set(gca,'TickLabelInterpreter', 'latex');
    set(gca,...
            'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
         'FontName','cmr10',...
        'FontSize',14,...
        'Box','off'); 
hold off


% Asset Functions
subplot(2,2,2)
hold on
plot(agrid, ap11(:,1),'DisplayName','J = 1, Low State' ,'LineWidth',2,'color',"#ff9696")
plot(agrid, ap11(:,2),'DisplayName','J = 1, High State','LineWidth',2,'color',"#FF0000")
plot(agrid, ap12(:,1),'DisplayName','J = 2, Low State' ,'LineWidth',2,'color',"#8293ff")
plot(agrid, ap12(:,2),'DisplayName','J = 2, High State','LineWidth',2,'color',"#0023ff")
xlabel('$a_t$', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('${a^j}_{t+1}(a_t,\varepsilon_t)$', 'Interpreter', 'latex', 'FontSize',16)
ax = gca;
set(gca,'TickLabelInterpreter', 'latex');
    set(gca,...
            'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
         'FontName','cmr10',...
        'FontSize',14,...
        'Box','off'); 
hold off

% Occupational Probabilities
figure(2)
hold on
plot(agrid, p11(:,1),'DisplayName','J = 1, Low State' ,'LineWidth',2,'color',"#ff9696")
plot(agrid, p11(:,2),'DisplayName','J = 1, High State','LineWidth',2,'color',"#FF0000")
plot(agrid, p12(:,1),'DisplayName','J = 2, Low State' ,'LineWidth',2,'color',"#8293ff")
plot(agrid, p12(:,2),'DisplayName','J = 2, High State','LineWidth',2,'color',"#0023ff")
yline(0.5,'--','DisplayName','50% Probability')
xlabel('$a_t$', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('$\theta$', 'Interpreter', 'latex', 'FontSize',16)
ax = gca;
set(gca,'TickLabelInterpreter', 'latex');
    set(gca,...
            'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
         'FontName','cmr10',...
        'FontSize',14,...
        'Box','off'); 
hold off
legend

% % Wealth distribution for occupation 1
figure(5)
title('Occupation 1')
hold on
fa = sum(faz11,2);
abar = linspace(amin,amax,1000); 
fabar = interp1(agrid,fa,abar); 
b = bar(abar,fabar,'hist'); 
%ylim([0,0.05]);
%xlim([-0.2,agrid(I)]);
xlabel('$a$','Interpreter','latex','FontSize',14); 
ylabel('Wealth distribution','Interpreter','latex','FontSize',14); 
hold off
% 
% % Wealth distribution for occupation 2
figure(6)
title('Occupation 2')
hold on
fa = sum(faz12,2);
abar = linspace(amin,amax,1000); 
fabar = interp1(agrid,fa,abar); 
b = bar(abar,fabar,'hist'); 
%ylim([0,0.05]);
%xlim([-0.2,agrid(I)]);
xlabel('$a$','Interpreter','latex','FontSize',14); 
ylabel('Wealth distribution','Interpreter','latex','FontSize',14); 
hold off
% 
% % Total wealth distribution
figure(7)
title('Total Wealth Distribution')
hold on
fa = sum(faz11.*p11 + faz12.*p12, 2);
abar = linspace(amin,amax,1000); 
fabar = interp1(agrid,fa,abar); 
b = bar(abar,fabar,'hist'); 
%ylim([0,0.05]);
%xlim([-0.2,agrid(I)]);
xlabel('$a$','Interpreter','latex','FontSize',14); 
ylabel('Wealth distribution','Interpreter','latex','FontSize',14); 
hold off

%% DECISION RULE Plots state 1

% Value function example with a taste shock. Set phi to 0.8 to get 
% a discontinuity
phi = 0.8;
shockstate = 1;
V_h = [v11(:,shockstate) + 0, v12(:,shockstate) + phi];
[V_1, argmax_1] = max(V_h,[],2);

% Recompute C
C_1 = zeros(I,1);
for i=1:I
    if argmax_1(i) == 1
        C_1(i) = c11(i,shockstate);
    else 
        C_1(i) = c12(i,shockstate);
    end
end

% Recompute AP
AP_1 = zeros(I,1);
for i=1:I
    if argmax_1(i) == 1
        AP_1(i) = ap11(i,shockstate);
    else 
        AP_1(i) = ap12(i,shockstate);
    end
end

% Value function plot
figure(3)
subplot(2,2,1)
hold on
% plot(agrid, v11(:,shockstate),'--','DisplayName','J = 1','LineWidth',2,'color'    ,"#FF0000")
% plot(agrid, v12(:,shockstate) + phi,'--','DisplayName','J = 2' ,'LineWidth',2,'color',"#0023ff")
plot(agrid, V_1,'DisplayName','$V$','LineWidth',2,'color',"k")
xlabel('$a_t$','Interpreter','latex','FontSize',14); 
ylabel('$V(a_t,\varepsilon_t = \varepsilon^L)$','Interpreter','latex','FontSize',14)
ax = gca;
set(gca,'TickLabelInterpreter', 'latex');
    set(gca,...
            'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
         'FontName','cmr10',...
        'FontSize',14,...
        'Box','off'); 
hold off

% Consumption Function
subplot(2,2,2)
hold on
% plot(agrid, c11(:,shockstate),'--','DisplayName','J = 1',   'LineWidth',2,'color',"#FF0000")
% plot(agrid, c12(:,shockstate),'--','DisplayName','J = 2' ,'LineWidth',2,'color'  ,"#0023ff")
plot(agrid, C_1,'DisplayName','$c_t$','LineWidth',2,'color',"k")
xlabel('$a_t$','Interpreter','latex','FontSize',14); 
ylabel('$c(a_t,\varepsilon_t = \varepsilon^L)$','Interpreter','latex','FontSize',14)
ax = gca;
set(gca,'TickLabelInterpreter', 'latex');
    set(gca,...
            'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
         'FontName','cmr10',...
        'FontSize',14,...
        'Box','off');
hold off

% Asset Function
subplot(2,2,3)
hold on
plot(agrid, agrid,'--','DisplayName','45 line','LineWidth',2)
%plot(agrid, ap11(:,shockstate),'--','DisplayName','J = 1',   'LineWidth',2,'color',"#FF0000")
%plot(agrid, ap12(:,shockstate),'--','DisplayName','J = 2' ,'LineWidth',2,'color'  ,"#0023ff")
plot(agrid, AP_1,'DisplayName','$a_{t+1}$','LineWidth',2,'color',"k")
xlabel('$a_t$','Interpreter','latex','FontSize',14); 
ylabel('$a_{t+1}(a_t,\varepsilon_t = \varepsilon^L)$','Interpreter','latex','FontSize',14)
ax = gca;
set(gca,'TickLabelInterpreter', 'latex');
    set(gca,...
            'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
         'FontName','cmr10',...
        'FontSize',14,...
        'Box','off');
hold off
hl = legend('show','Location','southeast');
set(hl, 'Interpreter','latex')

% Occupation Decision Rule
subplot(2,2,4)
hold on
plot(agrid, argmax_1,'LineWidth',2,'color',"k")
xlabel('$a_t$','Interpreter','latex','FontSize',14); 
ylabel('$o$','Interpreter','latex','FontSize',14)
ax = gca;
set(gca,'TickLabelInterpreter', 'latex');
    set(gca,...
            'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
         'FontName','cmr10',...
        'FontSize',14,...
        'Box','off');
hold off

%% SHock state 2
% Value function example with a taste shock. Set phi to 0.8 to get 
% a discontinuity
phi = 0.8;
shockstate = 2;
V_h = [v11(:,shockstate) + 0, v12(:,shockstate) + phi];
[V_2, argmax_2] = max(V_h,[],2);

% Recompute C
C_2 = zeros(I,1);
for i=1:I
    if argmax_2(i) == 1
        C_2(i) = c11(i,shockstate);
    else 
        C_2(i) = c12(i,shockstate);
    end
end

% Recompute AP
AP_2 = zeros(I,1);
for i=1:I
    if argmax_2(i) == 1
        AP_2(i) = ap11(i,shockstate);
    else 
        AP_2(i) = ap12(i,shockstate);
    end
end

% Value function plot
figure(4)
subplot(2,2,1)
hold on
% plot(agrid, v11(:,shockstate),'--','DisplayName','J = 1','LineWidth',2,'color'    ,"#FF0000")
% plot(agrid, v12(:,shockstate) + phi,'--','DisplayName','J = 2' ,'LineWidth',2,'color',"#0023ff")
plot(agrid, V_2,'DisplayName','$V$','LineWidth',2,'color',"k")
xlabel('$a_t$','Interpreter','latex','FontSize',14); 
ylabel('$V(a_t,\varepsilon_t = \varepsilon^H)$','Interpreter','latex','FontSize',14)
ax = gca;
set(gca,'TickLabelInterpreter', 'latex');
    set(gca,...
            'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
         'FontName','cmr10',...
        'FontSize',14,...
        'Box','off'); 
hold off
hl = legend('show');
set(hl, 'Interpreter','latex','Location','southeast')

% Consumption Function
subplot(2,2,2)
hold on
% plot(agrid, c11(:,shockstate),'--','DisplayName','J = 1',   'LineWidth',2,'color',"#FF0000")
% plot(agrid, c12(:,shockstate),'--','DisplayName','J = 2' ,'LineWidth',2,'color'  ,"#0023ff")
plot(agrid, C_2,'DisplayName','$c_t$','LineWidth',2,'color',"k")
xlabel('$a_t$','Interpreter','latex','FontSize',14); 
ylabel('$c(a_t,\varepsilon_t = \varepsilon^H)$','Interpreter','latex','FontSize',14)
ax = gca;
set(gca,'TickLabelInterpreter', 'latex');
    set(gca,...
            'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
         'FontName','cmr10',...
        'FontSize',14,...
        'Box','off');
hold off
hl = legend('show','Location','southeast');
set(hl, 'Interpreter','latex')

% Asset Function
subplot(2,2,3)
hold on
plot(agrid, agrid,'--','DisplayName','45 line','LineWidth',2)
% plot(agrid, ap11(:,shockstate),'--','DisplayName','J = 1',   'LineWidth',2,'color',"#FF0000")
% plot(agrid, ap12(:,shockstate),'--','DisplayName','J = 2' ,'LineWidth',2,'color'  ,"#0023ff")
plot(agrid, AP_2,'DisplayName','$a_{t+1}$','LineWidth',2,'color',"k")
xlabel('$a_t$','Interpreter','latex','FontSize',14); 
ylabel('$a_{t+1}(a_t,\varepsilon_t = \varepsilon^H)$','Interpreter','latex','FontSize',14)
ax = gca;
set(gca,'TickLabelInterpreter', 'latex');
    set(gca,...
            'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
         'FontName','cmr10',...
        'FontSize',14,...
        'Box','off');
hold off
hl = legend('show','Location','southeast');
set(hl, 'Interpreter','latex')

% Occupation Decision Rule
subplot(2,2,4)
hold on
plot(agrid, argmax_2,'LineWidth',2,'color',"k")
xlabel('$a_t$','Interpreter','latex','FontSize',14); 
ylabel('$o$','Interpreter','latex','FontSize',14)
ax = gca;
set(gca,'TickLabelInterpreter', 'latex');
    set(gca,...
            'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
         'FontName','cmr10',...
        'FontSize',14,...
        'Box','off');
hold off

%% Wealth dynamics and Occupational Choice
figure(5)
subplot(1,2,1)
hold on
plot(agrid, agrid,'--','DisplayName','45 line','LineWidth',2)
%plot(agrid, ap11(:,shockstate),'--','DisplayName','J = 1',   'LineWidth',2,'color',"#FF0000")
%plot(agrid, ap12(:,shockstate),'--','DisplayName','J = 2' ,'LineWidth',2,'color'  ,"#0023ff")
plot(agrid, AP_1,'DisplayName','$a_{t+1}$','LineWidth',2,'color',"k")
xlabel('$a_t$','Interpreter','latex','FontSize',14); 
ylabel('$a_{t+1}(a_t,\varepsilon_t = \varepsilon^L)$','Interpreter','latex','FontSize',14)
ax = gca;
set(gca,'TickLabelInterpreter', 'latex');
    set(gca,...
            'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
         'FontName','cmr10',...
        'FontSize',14,...
        'Box','off');
hold off
hl = legend('show','Location','southeast');
set(hl, 'Interpreter','latex')

subplot(1,2,2)
hold on
plot(agrid, agrid,'--','DisplayName','45 line','LineWidth',2)
% plot(agrid, ap11(:,shockstate),'--','DisplayName','J = 1',   'LineWidth',2,'color',"#FF0000")
% plot(agrid, ap12(:,shockstate),'--','DisplayName','J = 2' ,'LineWidth',2,'color'  ,"#0023ff")
plot(agrid, AP_2,'DisplayName','$a_{t+1}$','LineWidth',2,'color',"k")
xlabel('$a_t$','Interpreter','latex','FontSize',14); 
ylabel('$a_{t+1}(a_t,\varepsilon_t = \varepsilon^H)$','Interpreter','latex','FontSize',14)
ax = gca;
set(gca,'TickLabelInterpreter', 'latex');
    set(gca,...
            'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
         'FontName','cmr10',...
        'FontSize',14,...
        'Box','off');
hold off
hl = legend('show','Location','southeast');
set(hl, 'Interpreter','latex')

%% Aggregated Wealth Distributions
% to make better plot I aggregate in 10 bins the wealth distribution

figure(14)
title('Occupation 1')
hold on
fa = sum(faz11, 2);
abar = linspace(amin,amax,10); 
fabar = interp1(agrid, fa, abar); 
b = bar(abar,fabar,'hist'); 
ylim([0,0.05]);
xlim([-0.2,agrid(I)]);
xlabel('$a$','Interpreter','latex','FontSize',14); 
ylabel('Wealth distribution','Interpreter','latex','FontSize',14); 
hold off

%% Numerical derivative of c11

lhs = (c11(2:1000, 1).^-2).*diff(c11(:, 1));
rhs = (c12(2:1000, 1).^-2).*diff(c12(:, 1));
varo = (p11(2:1000, 1).*(1 - p11(2:1000, 1)));
dtda = varo.*(lhs - rhs);

lhs(lhs < 10^-03) = NaN;
rhs(rhs < 10^-03) = NaN;

lhs = smoothdata(lhs);
rhs = smoothdata(rhs);

figure(15)
hold on
xlabel('$a_t$','Interpreter','latex')
plot(agrid(2:1000),lhs,'DisplayName',"$u'(c^1)\partial c1$")
plot(agrid(2:1000),rhs,'DisplayName',"$u'(c^2)\partial c2$")
ax = gca;
set(gca,'TickLabelInterpreter', 'latex');
    set(gca,...
            'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
         'FontName','cmr10',...
        'FontSize',14,...
        'Box','off');
hold off
hl = legend('show','Location','northeast');
set(hl, 'Interpreter','latex')


%%
distribution = sum(faz11.*p11 + faz12.*p12, 2);
m = mean(distribution);
s = std(distribution);
pdf = lognpdf(agrid, 2.2941, 8.3937);
figure(100)
hold on
plot(agrid, pdf)
grid on

%%
figure(1)
hold on
plot(agrid, argmax_1,'LineWidth',2,'color',"k")
xlabel('$a_t$','Interpreter','latex','FontSize',14); 
ylabel('$o$','Interpreter','latex','FontSize',14)
ax = gca;
ax.YTick = [1,2];
set(gca,'TickLabelInterpreter', 'latex');
    set(gca,...
            'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
         'FontName','cmr10',...
        'FontSize',14,...
        'Box','off');
hold off

figure(2)
hold on
plot(agrid, argmax_2,'LineWidth',2,'color',"k")
xlabel('$a_t$','Interpreter','latex','FontSize',14); 
ylabel('$o$','Interpreter','latex','FontSize',14)
ax = gca;
ax.YTick = [1,2];
set(gca,'TickLabelInterpreter', 'latex');
    set(gca,...
            'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
         'FontName','cmr10',...
        'FontSize',14,...
        'Box','off');
hold off
