%% APPM 2360 Project 1 - Tribble Population Modeling
%  Model:  dy/dt = f(y) = a*y - b*y^2 - H(y),   H(y) = p*y^3/(y^3 + q)
%  t in days, y in hundreds of tribbles.
%
%  Sections:
%    Q4a - check of the analytic logistic solution (p = 0) against ode45
%    Q6  - plots of f(y) and equilibrium solutions (fzero)
%    Q7  - direction fields with overlaid solutions (quiver + ode45)
%    Q8  - stability table and basins of attraction
%    Q9  - numbers that support the stocking recommendation
%    Q10 - seasonal hunting H(t,y) = p*y^3/(q+y^3)*|sin(pi*t/365)|
%
%  Run the whole script (or one section at a time with "Run Section").
%  Figures are also saved as PNG files in the current folder.

clear; close all; clc;

%% Parameters
a     = 0.75;               % 1/day
p     = 1.5;                % hundreds of tribbles / day
q     = 1.25;               % (hundreds of tribbles)^3
bvals = [0.005 0.05 0.10];  % 1/(hundreds of tribbles * day)
names = {'brown (b = 0.005)', 'white (b = 0.05)', 'grey (b = 0.10)'};

H    = @(y)   p*y.^3 ./ (y.^3 + q);          % hunting term
f    = @(y,b) a*y - b*y.^2 - H(y);           % right-hand side of (1)
dfdy = @(y,b) (f(y+1e-6,b) - f(y-1e-6,b)) / 2e-6;   % numerical f'(y)

odeopts = odeset('RelTol',1e-8,'AbsTol',1e-10);

%% Q4a - Analytic logistic solution (H = 0) vs. numerical solution
%  y(t) = (a/b) / (1 + (a/(b*y0) - 1) * exp(-a*t)),  y -> a/b as t -> inf
b  = 0.05;  y0 = 2;  tspan = [0 20];
yLogistic = @(t) (a/b) ./ (1 + (a/(b*y0) - 1)*exp(-a*t));
[tn, yn]  = ode45(@(t,y) a*y - b*y.^2, tspan, y0, odeopts);

figure('Name','Q4a');
plot(tn, yn, 'b-', 'LineWidth', 2); hold on;
tt = linspace(tspan(1), tspan(2), 40);
plot(tt, yLogistic(tt), 'ro', 'MarkerSize', 5);
plot(tspan, [a/b a/b], 'k--');
xlabel('t (days)'); ylabel('y (hundreds of tribbles)');
title(sprintf('Logistic model, no hunting (b = %.2f, y_0 = %g)', b, y0));
legend('ode45', 'analytic solution', 'carrying capacity a/b', 'Location', 'southeast');
grid on;
saveas(gcf, 'Q4a_logistic_check.png');
fprintf('Q4a: max |analytic - ode45| = %.2e\n\n', max(abs(yn - yLogistic(tn))));

%% Q6 - Plot f(y) and find the equilibrium solutions
eqs  = cell(1,3);    % equilibria for each b
stab = cell(1,3);    % 'stable' / 'unstable'

figure('Name','Q6','Position',[100 100 1200 650]);
for k = 1:3
    b = bvals(k);

    % --- find equilibria: y = 0 plus every sign change of f on a fine grid
    yGrid = linspace(1e-6, a/b + 5, 200000);   % f < 0 for all y > a/b
    fGrid = f(yGrid, b);
    idx   = find(fGrid(1:end-1).*fGrid(2:end) < 0);
    roots_k = 0;
    for i = idx
        roots_k(end+1) = fzero(@(y) f(y,b), [yGrid(i) yGrid(i+1)]); %#ok<SAGROW>
    end
    eqs{k} = roots_k;

    % --- classify by the sign of f'(y*)
    s = cell(size(roots_k));
    for j = 1:numel(roots_k)
        if dfdy(roots_k(j), b) < 0, s{j} = 'stable'; else, s{j} = 'unstable'; end
    end
    stab{k} = s;

    % --- full-range plot
    yMax = 1.1*a/b;
    yy = linspace(0, yMax, 4000);
    subplot(2,3,k);
    plot(yy, f(yy,b), 'b-', 'LineWidth', 1.5); hold on;
    plot([0 yMax], [0 0], 'k-');
    plot(roots_k, zeros(size(roots_k)), 'ro', 'MarkerFaceColor', 'r');
    xlabel('y (hundreds)'); ylabel('f(y)  (hundreds/day)');
    title(sprintf('f(y), b = %.3g', b)); grid on; xlim([0 yMax]);

    % --- zoom near y = 0 so small equilibria are not missed
    subplot(2,3,k+3);
    yz = linspace(0, 4, 2000);
    plot(yz, f(yz,b), 'b-', 'LineWidth', 1.5); hold on;
    plot([0 4], [0 0], 'k-');
    rz = roots_k(roots_k <= 4);
    plot(rz, zeros(size(rz)), 'ro', 'MarkerFaceColor', 'r');
    xlabel('y (hundreds)'); ylabel('f(y)  (hundreds/day)');
    title(sprintf('Zoom 0 \\leq y \\leq 4, b = %.3g', b)); grid on;
end
saveas(gcf, 'Q6_f_of_y.png');

%% Q8(a) - Table of equilibria and their stability (printed to Command Window)
fprintf('Q6/Q8(a): Equilibrium solutions of dy/dt = ay - by^2 - py^3/(y^3+q)\n');
fprintf('%-8s %-14s %-12s %-12s %-10s\n', 'b', 'y* (hundreds)', 'y* (tribbles)', 'f''(y*)', 'stability');
fprintf('%s\n', repmat('-', 1, 62));
for k = 1:3
    for j = 1:numel(eqs{k})
        fprintf('%-8.3f %-14.4f %-12.0f %-12.4f %-10s\n', bvals(k), eqs{k}(j), ...
                100*eqs{k}(j), dfdy(eqs{k}(j), bvals(k)), stab{k}{j});
    end
    fprintf('\n');
end

%% Q8(b) - Basins of attraction of the stable equilibria (y0 >= 0)
%  For y' = f(y), the basin of a stable y* runs between the neighbouring
%  unstable equilibria (or to +infinity if there is none above it).
fprintf('Q8(b): Basins of attraction (initial populations y0, any t0)\n');
for k = 1:3
    e = eqs{k};  s = stab{k};
    for j = 1:numel(e)
        if strcmp(s{j}, 'stable')
            lower = 0;  upper = Inf;
            below = find(strcmp(s(1:j-1), 'unstable'), 1, 'last');
            above = find(strcmp(s(j+1:end), 'unstable'), 1, 'first');
            if ~isempty(below), lower = e(below); end
            if ~isempty(above), upper = e(j+above); end
            fprintf('  b = %.3f: y* = %8.4f  <-  %.4f < y0 < %s\n', bvals(k), e(j), ...
                    lower, num2str(upper, '%.4f'));
        end
    end
end
fprintf('\n');

%% Q7 - Direction fields with overlaid solutions
figure('Name','Q7','Position',[100 100 1300 420]);
tEnd = [80 30 30];          % days shown (b = 0.005 needs longer: slow passage near y ~ 1-1.5)
for k = 1:3
    b = bvals(k);  e = eqs{k};
    yTop = 1.2*max(e);  if yTop < 3, yTop = 3; end

    subplot(1,3,k); hold on;

    % direction field (arrows normalised to equal length)
    [T, Y] = meshgrid(linspace(0, tEnd(k), 25), linspace(0, yTop, 25));
    dY = f(Y, b);  dT = ones(size(dY));
    L  = sqrt(dT.^2 + (dY*tEnd(k)/yTop).^2);        % scale to axis aspect
    quiver(T, Y, dT./L, dY./L, 0.5, 'Color', [0.6 0.6 0.6]);

    % equilibrium solutions
    for j = 1:numel(e)
        if strcmp(stab{k}{j}, 'stable'), ls = 'r-'; else, ls = 'r--'; end
        plot([0 tEnd(k)], [e(j) e(j)], ls, 'LineWidth', 1.5);
    end

    % initial conditions: just above/below every equilibrium + spread across the range
    ics = [e(:)' + 0.05*max(1,e(:)'), e(e>0) - 0.05*max(1,e(e>0)), ...
           linspace(0.1, yTop, 8)];
    ics = unique(ics(ics > 0 & ics <= yTop));
    for y0 = ics
        [t, y] = ode45(@(t,y) f(y,b), [0 tEnd(k)], y0, odeopts);
        plot(t, y, 'b-', 'LineWidth', 1);
    end

    xlim([0 tEnd(k)]); ylim([0 yTop]);
    xlabel('t (days)'); ylabel('y (hundreds of tribbles)');
    title(sprintf('Direction field, b = %.3g', b)); box on;
end
saveas(gcf, 'Q7_direction_fields.png');

%% Q8(c)/Q9 - Supporting numbers for the long-term behaviour & stocking choice
fprintf('Q8(c)/Q9: Long-term populations\n');
for k = 1:3
    st = eqs{k}(strcmp(stab{k}, 'stable'));
    fprintf('  %-18s stable populations: %s tribbles\n', names{k}, ...
            mat2str(round(100*st)));
end
e2 = eqs{2};
fprintf(['  White tribbles need y0 > %.4f (about %d tribbles) to reach the\n' ...
         '  high equilibrium of %.0f tribbles; otherwise they settle at %.0f.\n\n'], ...
         e2(3), ceil(100*e2(3)), 100*e2(4), 100*e2(2));

%% Q10 - Seasonal hunting  H(t,y) = p*y^3/(q+y^3) * |sin(pi*t/365)|
Hs = @(t,y) p*y.^3 ./ (q + y.^3) .* abs(sin(pi*t/365));
years = 5;  tspan = [0 365*years];
opts10 = odeset(odeopts, 'MaxStep', 1);   % resolve the seasonal forcing

figure('Name','Q10','Position',[100 100 1300 420]);
for k = 1:3
    b = bvals(k);
    subplot(1,3,k); hold on;
    ics = unique([0.2 0.5 1 1.5 2 2.5 linspace(0.1, 1.1*a/b, 6)]);
    for y0 = ics
        [t, y] = ode45(@(t,y) a*y - b*y.^2 - Hs(t,y), tspan, y0, opts10);
        plot(t/365, y, 'LineWidth', 1);
    end
    % reference: autonomous equilibria from Q6
    for j = 1:numel(eqs{k})
        plot([0 years], [eqs{k}(j) eqs{k}(j)], 'k:');
    end
    xlabel('t (years)'); ylabel('y (hundreds of tribbles)');
    title(sprintf('Seasonal hunting, b = %.3g', b)); grid on; box on;
    xlim([0 years]);
end
saveas(gcf, 'Q10_seasonal_hunting.png');

% Zoom on one year for b = 0.05 to show the seasonal oscillation
figure('Name','Q10 detail');
b = 0.05;
subplot(2,1,1); hold on;
for y0 = [0.5 1 1.5 2 2.5 5 10 15]
    [t, y] = ode45(@(t,y) a*y - b*y.^2 - Hs(t,y), tspan, y0, opts10);
    plot(t, y, 'LineWidth', 1);
end
ylabel('y (hundreds of tribbles)');
title('Seasonal hunting, b = 0.05'); grid on; xlim(tspan);
subplot(2,1,2);
tt = linspace(tspan(1), tspan(2), 2000);
plot(tt, abs(sin(pi*tt/365)), 'k-');
xlabel('t (days)'); ylabel('|sin(\pi t/365)|');
title('Hunting intensity (0 = midwinter, 1 = midsummer)'); grid on; xlim(tspan);
saveas(gcf, 'Q10_seasonal_b005_detail.png');

fprintf('Done. Figures saved as PNG files in: %s\n', pwd);