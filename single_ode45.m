%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rudimentary script to solve an initial value problem consisting of the
% following initial value problem y'=f(t,y) on the interval a<=t<=b
%
% y'(t) = t(1-y^2)  y(a) = y_0
%
% Look for "*** Your entry i ***" where i=1,2,3,4,5,6 to make appropriate
% changes to solve your specific initial value problem.
%
% Running the script as is outputs two solutions: one decreasing to the
% equilibrium solution (y=1) and one increasing to it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Housekeeping

clc;        % clear the command window
clear;      % clear all variables
close all;  % close all figures

    % ODE45 options

options = odeset('AbsTol',1e-10,'RelTol',1e-7);

    % Independent variable

a = 0;      % left endpoint of t interval  *** Your entry 1 ***
b = 5;      % right endpoint of t interval *** Your entry 2 ***
h = 0.001;  % stepsize                     *** Your entry 3 ***
tt = a:h:b;

    % Initial values y(a) (create solutions for 
    % multiple initial values if wanted)

y1 = 1.5;    % initial value #1     *** Your entry 4 ***
y2 = -0.55;    % initial value #2     *** Your entry 5 ***
y0 = [y1 y2];

    % Solve the equation using ode45

[t,soln] = ode45(@(t,y) rhs(t,y),tt,y0,options);

    % Plot solution(s)

figure(1);
hold on
for i = 1:length(y0)
    plot(t,soln(:,i));
end
xlabel('$t$','Interpreter','latex','FontSize',20);
ylabel('$y$','Interpreter','latex','FontSize',20,'Rotation',0)

    % Right hand side of the differential equation

function [yprime] = rhs(t,y)   % do not change this line
    % *** Your entry 6 ***  (start entering your f(t,y) here)
    yprime = t.*(1-y.^2);
    % *** Your entry 6 ***  (stop entering your f(t,y) here)
end                        % do not change this line
