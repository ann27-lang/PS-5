%% Problem 3

%% phase-portrait plots for the gene toggle switch

close all;
clear all;

%part b and c
nval = [1 2];
alpha = 10;
coef = ["n = 1", 'n = 2'];
vs = [-0.5+sqrt(41)/2, 2];
us = [-0.5+sqrt(41)/2, 2];
for z = 1:2
    n = nval(z);
    phaseplot3(n,alpha);
    title([coef(z) 'Phase portrait plot'])
    xlabel('u concentration')
    ylabel('v concentration')
    legend('phase vectors','u nullcline','v nullcline')
    
    %part e
    fu = -1;
    fv = (-vs(z)^(nval(z)-1)*alpha*nval(z))/(1+vs(z)^nval(z))^2;
    gu = (-us(z)^(nval(z)-1)*alpha*nval(z))/(1+us(z)^nval(z))^2;
    gv = -1;
    
    J = [fu fv; gu gv];
    eigval(z,:) = eig(J);
end

%% part e bonus attempt (:
variable = linspace(0,15,40);

n = 2; % CHANGE WHEN NEED NEXT ACRIT VALUE

%Conditions
a1 = variable;		% synthesis rate parameter

%Initialize
a1xmulti = [];
a1xmono1 = [];
a1xmono2 = [];
x1ymulti = [];
x2ymulti = [];
x1ymono1 = [];
x2ymono1 = [];
x1ymono2 = [];
x2ymono2 = [];

%Plots
for z = 1:length(variable)
    alpha = a1(z);
    [x1_null,x2_null,x1_test, x2_test] = phaseplot3(n,alpha);
    [a1c, a2c] = polyxpoly(x1_test,x2_null,x1_null,x2_test);
    a1c = a1c';
    a2c = a2c';
    if z>=6 && z<=20
        a1xmulti = [a1xmulti repelem(alpha,length(a1c))];
        x1ymulti = [x1ymulti a1c];
        x2ymulti = [x2ymulti a2c];
    elseif  z<=5
        a1xmono1 = [a1xmono1 repelem(alpha,length(a1c))];
        x1ymono1 = [x1ymono1 a1c];
        x2ymono1 = [x2ymono1 a2c];
        if z == 5
            val5 = alpha
        end
    else
        a1xmono2 = [a1xmono2 repelem(alpha,length(a1c))];
        x1ymono2 = [x1ymono2 a1c];
        x2ymono2 = [x2ymono2 a2c];
        if z == 21
            val21 = alpha
        end
    end
end

figure
hold on
plot(a1xmono1,x1ymono1,'-r',a1xmulti,x1ymulti,'*r',a1xmono2,x1ymono2,'-r')
xlabel('alpha synthesis rate')
ylabel('x1 Protein concentration')
legend('monostable','multistable')
title('alpha synthesis rate v x1 Protein concentrations')
hold off

figure
hold on
plot(a1xmono1,x2ymono1,'-b',a1xmulti,x2ymulti,'*b',a1xmono2,x2ymono2,'-b')
xlabel('alpha synthesis rate')
ylabel('x2 Protein concentration')
legend('monostable','multistable')
title('alpha synthesis rate v x2 Protein concentrations')
hold off

function [u_null,v_null,u_test,v_test] = phaseplot3(n,alpha)

%Test R1 and R2 concentrations
size_test = 85; %number of test points in grid
u_test = linspace(0,10,size_test);
v_test = linspace(0,10,size_test);

%Calculate phase vectors
for i = 1:size_test
    for j = 1:size_test
        du(j,i) = alpha/(1+v_test(j)^n) - u_test(i);
        dv(j,i) = alpha/(1+u_test(i)^n) - v_test(j);
    end
end

%Calculate null clines
for i = 1:size_test
    u_null(i) = alpha/(1+v_test(i)^n);
    v_null(i) = alpha/(1+u_test(i)^n);
end
if alpha == 10
    figure
    quiver(u_test,v_test,du,dv,1.5,'k')
    hold on
    plot(u_test,v_null,'-r')
    plot(u_null,v_test,'-b')
    xlim([0 10])
    ylim([0 10])
    hold off
end
end