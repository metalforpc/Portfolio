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
xlabel('$a$', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('$v(a,\varepsilon, j)$', 'Interpreter', 'latex', 'FontSize',16)
hl = legend('$J=1$ \, $\varepsilon = \varepsilon^L$', '$J=1$ \, $\varepsilon =\varepsilon^H$', '$J = 2$ \, $\varepsilon =\varepsilon^L$', '$J = 2$ \, $\varepsilon =\varepsilon^L$');
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
xlabel('$a$', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('${c}(a,\varepsilon)$', 'Interpreter', 'latex', 'FontSize',16)
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
xlabel('$a$', 'Interpreter', 'latex', 'FontSize', 16)
ylabel("${a'}(a, \varepsilon)$", 'Interpreter', 'latex', 'FontSize',16)
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
xlabel('$a$', 'Interpreter', 'latex', 'FontSize', 16)
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
xlabel('$a$','Interpreter','latex','FontSize',14); 
ylabel('$V(a,\varepsilon = \varepsilon^L)$','Interpreter','latex','FontSize',14)
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
plot(agrid, C_1,'DisplayName','$c$','LineWidth',2,'color',"k")
xlabel('$a$','Interpreter','latex','FontSize',14); 
ylabel('$c(a,\varepsilon = \varepsilon^L)$','Interpreter','latex','FontSize',14)
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
plot(agrid, agrid,'--','DisplayName','45° line','LineWidth',2)
%plot(agrid, ap11(:,shockstate),'--','DisplayName','J = 1',   'LineWidth',2,'color',"#FF0000")
%plot(agrid, ap12(:,shockstate),'--','DisplayName','J = 2' ,'LineWidth',2,'color'  ,"#0023ff")
plot(agrid, AP_1,'DisplayName',"$a'$",'LineWidth',2,'color',"k")
xlabel('$a$','Interpreter','latex','FontSize',14); 
ylabel("$a'(a,\varepsilon = \varepsilon^L)$",'Interpreter','latex','FontSize',14)
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
xlabel('$a$','Interpreter','latex','FontSize',14); 
ylabel('$j$','Interpreter','latex','FontSize',14)
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
xlabel('$a$','Interpreter','latex','FontSize',14); 
ylabel('$V(a,\varepsilon = \varepsilon^H)$','Interpreter','latex','FontSize',14)
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
xlabel('$a$','Interpreter','latex','FontSize',14); 
ylabel('$c(a,\varepsilon = \varepsilon^H)$','Interpreter','latex','FontSize',14)
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
plot(agrid, agrid,'--','DisplayName','45° line','LineWidth',2)
% plot(agrid, ap11(:,shockstate),'--','DisplayName','J = 1',   'LineWidth',2,'color',"#FF0000")
% plot(agrid, ap12(:,shockstate),'--','DisplayName','J = 2' ,'LineWidth',2,'color'  ,"#0023ff")
plot(agrid, AP_2,'DisplayName',"$a'$",'LineWidth',2,'color',"k")
xlabel('$a$','Interpreter','latex','FontSize',14); 
ylabel("$a'(a,\varepsilon = \varepsilon^H)$",'Interpreter','latex','FontSize',14)
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
xlabel('$a$','Interpreter','latex','FontSize',14); 
ylabel('$j$','Interpreter','latex','FontSize',14)
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
plot(agrid, agrid,'--','DisplayName','45° line','LineWidth',2)
%plot(agrid, ap11(:,shockstate),'--','DisplayName','J = 1',   'LineWidth',2,'color',"#FF0000")
%plot(agrid, ap12(:,shockstate),'--','DisplayName','J = 2' ,'LineWidth',2,'color'  ,"#0023ff")
plot(agrid, AP_1,'DisplayName',"$a'$",'LineWidth',2,'color',"k")
xlabel('$a$','Interpreter','latex','FontSize',14); 
ylabel("$a'(a,\varepsilon = \varepsilon^L)$",'Interpreter','latex','FontSize',14)
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
plot(agrid, agrid,'--','DisplayName','45° line','LineWidth',2)
% plot(agrid, ap11(:,shockstate),'--','DisplayName','J = 1',   'LineWidth',2,'color',"#FF0000")
% plot(agrid, ap12(:,shockstate),'--','DisplayName','J = 2' ,'LineWidth',2,'color'  ,"#0023ff")
plot(agrid, AP_2,'DisplayName',"$a'$",'LineWidth',2,'color',"k")
xlabel('$a$','Interpreter','latex','FontSize',14); 
ylabel("$a'(a,\varepsilon = \varepsilon^H)$",'Interpreter','latex','FontSize',14)
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

