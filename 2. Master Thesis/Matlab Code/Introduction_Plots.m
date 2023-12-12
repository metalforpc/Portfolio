%% Parameters
gamma11 = 1.5;
gamma12 = 0.5;
gamma21 = 0.5;
gamma22 = 1.5;

%% Constant Labor Supply
L1S = @(M11) M11*gamma11;
L2S = @(M22) M22*gamma22;

figure(1)
subplot(1,2,1)
hold on

plot(0:.1:1, repelem(gamma11,11))

hold off

subplot(1,2,2)
hold on
plot(0:.1:1, repelem(gamma11,11))
hold off



%% Shift in supply due to type 1.
l1S = @(M11, M21) M11*gamma11 + M21*gamma21;
l2S = @(M12, M22) M12*gamma12 + M22*gamma22;

M11 = 0:.1:1;
M12 = 1 - M11;

M21 = 0.5;
M22 = 1 - M21;


L1S = l1S(M11, M21);
L2S = l2S(M12, M22);

subplot(1,2,1)
hold on 
xlabel('$M_{1, 2}$', 'Interpreter', 'latex')
ylabel('${L_1}^S$', 'Interpreter', 'latex')
plot(M12, L1S,'red')
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

subplot(1,2,2)
hold on
xlabel('$M_{1, 2}$', 'Interpreter', 'latex')
ylabel('${L_2}^S$', 'Interpreter', 'latex')
plot(M12, L2S,'blue')
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

%% L1 an L2 with different combinations of Mass
l1S = @(M11, M21) M11*gamma11 + M21*gamma21;
l2S = @(M12, M22) M12*gamma12 + M22*gamma22;

M11 = 0.5;
M12 = 1-M11;
M21 = 0.5;
M22 = 1 - M21;

figure(2)
hold on
for i=0.01:.1:0.99
    for j = 0.01:.1:0.99
        M12 = 1 - i;
        M22 = 1 - j;
        L1 = l1S(i, j);
        L2 = l2S(M12, M22);
        scatter(L1,L2,'blue');
    end
end
xlabel('${L_1}^S$', 'Interpreter', 'latex')
ylabel('${L_2}^S$', 'Interpreter', 'latex')
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

%% Grid for wage movements

[M11grid, M12grid, M21grid, M22grid] = ndgrid((0.01:.1:0.99), (0.01:.1:0.99), (0.01:.1:0.99), (0.01:.1:0.99));
m11 = M11grid(:);
m12 = M12grid(:);
m21 = M21grid(:);
m22 = M22grid(:);

mass_domain = [m11, m12, m21, m22];
mass_domain(:) = mass_domain(randperm(numel(mass_domain)));

L1S = mass_domain(:,1)*gamma11 + mass_domain(:,3)*gamma21;
L2S = mass_domain(:,2)*gamma12 + mass_domain(:,4)*gamma22;

omega = 0.33;
eta = 0.33;

w1 = [];
w2 = [];

for i=1:10000
    w1(i) = omega*(L1S(i)^(1-omega))*(L2S(i)^eta);
    w2(i) = eta*(L1S(i)^(omega))*(L2S(i)^(1 - eta));
end


%% PLOT
figure(1)
subplot(1,2,1)
hold on 
plot(w1(1:10),'DisplayName','$w_1$')
plot(w2(1:10),'DisplayName','$w_2$')
ax = gca;
set(gca,'TickLabelInterpreter', 'latex');
    set(gca,...
            'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
         'FontName','cmr10',...
        'FontSize',14,...
        'Box','off'); 
hl = legend('show');
set(hl, 'Interpreter','latex','Location','southeast')
hold off

subplot(1,2,2)
hold on 
plot(L1S(1:10),'DisplayName','${L_1}^S$')
plot(L2S(1:10),'DisplayName','${L_2}^S$')
ax = gca;
set(gca,'TickLabelInterpreter', 'latex');
    set(gca,...
            'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
         'FontName','cmr10',...
        'FontSize',14,...
        'Box','off'); 
hl = legend('show');
set(hl, 'Interpreter','latex','Location','southeast')
hold off