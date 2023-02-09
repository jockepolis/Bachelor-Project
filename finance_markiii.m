%% The starting values and the exact value of the integral

clear all; close all; clc;
fprintf('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*\n');
fprintf('Script that integrates the s-dimensional Rosenbrock function\n');
fprintf('over the unit-square in every dim: [0,1]^s .\n');
fprintf('By Axel Englund & Joakim Svensson April 2022\n');
fprintf('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*\n');
D = input('Number of dimensions to integrate over: '); 
a = 1; b = 100;
x=linspace(0,1,100);
N=1000;
P = sobolset(D);
Int_exact=(D-1)*(3*b+15*a^2-15*a+5)/15;
err_MC = zeros(N,1);
err_QMC = zeros(N,1);
%% MC and QMC algorithms with the errors and a fitted polynomial

for n=1:N
    f_MC=zeros(N,1);
    for i=1:n
        x_MC = rand(D,1).*(x(end)-x(1)) + x(1);
        f_MC_elements = f_MC + Rosenbrock3(x_MC,a,b);
        f_MC(i) = f_MC_elements;
    end
    Sobol_coordinates = net(P,n);
    
    f_QMC=zeros(N,1);
    for i=1:n
        x_QMC = Sobol_coordinates(i,:);
        f_QMC_elements = f_QMC + Rosenbrock3(x_QMC,a,b);
        f_QMC(i) = f_QMC_elements;
    end
    
    err_MC(i)=var(f_MC);
    err_QMC(i)=var(f_QMC);
    if mod(i,N/100)==0
        fprintf(['',num2str(i/N *100),' %% done \n'])
    end
end

A =[log(1:1:N)', ones(N,1)];
fit_QMC = inv(A'*A)*A'*log(err_QMC);
fit_MC = inv(A'*A)*A'*log(err_MC);
%% Plots

figure(1)

loglog(1:1:N,err_MC,'k.',1:1:N,err_QMC,'r.',1:1:N,1./sqrt(1:1:N),'k-',1:1:N,1./(1:1:N),'r-')
legend('Monte-Carlo Method','Quasi-Monte-Carlo Method','$1/\sqrt{N}$',...
    '$1/N$','Location','Best','Interpreter','Latex');
xlabel('Number of evaluation points','Interpreter','Latex');
ylabel('Error','Interpreter','Latex');
title([num2str(D),' Dimensions'],'Interpreter','Latex');
grid on;
set(gcf,'color','w');

figure(2)

loglog(1:1:N,err_MC,'k.',1:1:N,err_QMC,'r.',1:1:N,exp(fit_MC(2)).*(1:1:N).^(fit_MC(1)),'k-',...
    1:1:N,exp(fit_QMC(2)).*(1:1:N).^(fit_QMC(1)),'r-')
legend('Monte-Carlo Method','Quasi-Monte-Carlo Method',...
    ['$\mathcal{O}\left(',num2str(exp(fit_MC(2))),'\cdot N^{',num2str(fit_MC(1)),'}\right)_{MC}$'],...
    ['$\mathcal{O}\left(',num2str(exp(fit_QMC(2))),'\cdot N^{',num2str(fit_QMC(1)),'}\right)_{QMC}$'],...
    'Location','Best','Interpreter','Latex');
xlabel('Number of evaluation points','Interpreter','Latex');
ylabel('Error','Interpreter','Latex');
title([num2str(D),' Dimensions'],'Interpreter','Latex');
grid on;
set(gcf,'color','w');

%% Rosenbrock's function

function y = Rosenbrock3(x,a,b)
    if nargin<3
        a = 1;
        b = 100;
    end
    N=length(x);
    y=0;
    for i=1:N-1
        y = y + (b*(x(i+1) - x(i)^2)^2 + (a-x(i))^2);
    end
end