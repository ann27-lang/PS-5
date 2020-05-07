%CHEME 5440
%Problem 2
%part e
%Parameters
cao = 1;
cro = 10;
xlimca = [0 200 10];
ylimcr = [0 100 10];
xinit = [cao cro];

ror = 1;
rr = 100;
roa = 100;
ra = 5000;
da = 30;

f = @(t,c) [-da*c(1) + (roa+ra*c(1)^2)/(1+c(1)^2+c(2)^2); ...
    -c(2) + (ror+rr*c(1)^2)/(1+c(1)^2)];
[t,cf] = ode45(f,[0 50],[1 10]);

figure
plot(t,cf(:,1))
hold on
plot(t,cf(:,2))
title("time versus concentration")
legend("Species A","Species B")
xlabel("time")
ylabel("concentration")
ylim([0 200])
