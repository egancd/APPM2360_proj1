clc; clear; close all;

%% Parameters
a     = 0.75;               % 1/day
p     = 1.5;                % hundreds of tribbles / day
q     = 1.25;               % (hundreds of tribbles)^3
bvals = [0.005 0.05 0.10];  % 1/(hundreds of tribbles * day)
names = {'brown (b = 0.005)', 'white (b = 0.05)', 'grey (b = 0.10)'};

H    = @(y)   p*y.^3 ./ (y.^3 + q);                  % hunting term
f    = @(y,b) a*y - b*y.^2 - H(y);                   % right-hand side of (1)
dfdy = @(y,b) (f(y+1e-6,b) - f(y-1e-6,b)) / 2e-6;    % numerical f'(y)

options = odeset('AbsTol',1e-10,'RelTol',1e-7);

%  y(t) = (a/b) / (1 + (a/(b*y0) - 1) * exp(-a*t)),  y -> a/b as t -> inf
b  = 0.05;
t0 = 0;  tf = 20;  h = 0.01;  tt = t0:h:tf;          % as in single_ode45.m
y0 = [2 20];                                         % one below, one above a/b
yLogistic = @(t,y0) (a/b) ./ (1 + (a./(b*y0) - 1).*exp(-a*t));

[t, soln] = ode45(@(t,y) a*y - b*y.^2, tt, y0, options);

figure(1); hold on;
for i = 1:length(y0)
    plot(t, soln(:,i), 'b-', 'LineWidth', 2);
    plot(t(1:100:end), yLogistic(t(1:100:end), y0(i)), 'ro', 'MarkerSize', 5);
end
plot([t0 tf], [a/b a/b], 'k--');
xlabel('t (days)'); ylabel('y (hundreds of tribbles)');
title(sprintf('Logistic model, no hunting (b = %.2f)', b));
legend('ode45', 'analytic solution', 'Location', 'east'); grid on;
saveas(gcf, 'Q4a_logistic_check.png');
err = max(max(abs(soln - [yLogistic(t,y0(1)) yLogistic(t,y0(2))])));
fprintf('Q4a: max |analytic - ode45| = %.2e\n\n', err);

eqs  = cell(1,3);    % equilibria for each b
stab = cell(1,3);    % 'stable' / 'unstable'

figure(2); set(gcf, 'Position', [100 100 1200 650]);
for k = 1:3
    b = bvals(k);

    % find equilibria: y = 0 plus every sign change of f on a fine grid
    yGrid = linspace(1e-6, a/b + 5, 200000);   % f < 0 for all y > a/b
    fGrid = f(yGrid, b);
    idx   = find(fGrid(1:end-1).*fGrid(2:end) < 0);
    roots_k = 0;
    for i = idx
        roots_k(end+1) = fzero(@(y) f(y,b), [yGrid(i) yGrid(i+1)]); %#ok<SAGROW>
    end
    eqs{k} = roots_k;

    % classify by the sign of f'(y*)
    s = cell(size(roots_k));
    for j = 1:numel(roots_k)
        if dfdy(roots_k(j), b) < 0, s{j} = 'stable'; else, s{j} = 'unstable'; end
    end
    stab{k} = s;

    % full-range plot
    yMax = 1.1*a/b;
    yy = linspace(0, yMax, 4000);
    subplot(2,3,k);
    plot(yy, f(yy,b), 'b-', 'LineWidth', 1.5); hold on;
    plot([0 yMax], [0 0], 'k-');
    plot(roots_k, zeros(size(roots_k)), 'ro', 'MarkerFaceColor', 'r');
    xlabel('y (hundreds)'); ylabel('f(y)  (hundreds/day)');
    title(sprintf('f(y), b = %.3g', b)); grid on; xlim([0 yMax]);

    % zoom near y = 0 so small equilibria are not missed
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

tEnd = [80 30 30];      % days shown (b = 0.005 needs longer: slow passage near y ~ 1-1.5)
h    = 0.01;            % output step for ode45 (as in single_ode45.m)

figure(3); set(gcf, 'Position', [100 100 1300 420]);
for k = 1:3
    b = bvals(k);  e = eqs{k};
    yTop = 1.2*max(e);  if yTop < 3, yTop = 3; end

    subplot(1,3,k); hold on;

    % dirfield.m evaluates func2str(f) inside its own workspace, where a, b,
    % p, q do not exist, so the parameter values are written into the
    % function as numbers with str2func.
    fdir = str2func(sprintf('@(t,y) %.12g*y - %.12g*y.^2 - %.12g*y.^3./(y.^3 + %.12g)', ...
                            a, b, p, q));
    dirfield(fdir, linspace(0, tEnd(k), 25), linspace(0, yTop, 25), 0.5);

    for j = 1:numel(e)
        if strcmp(stab{k}{j}, 'stable'), ls = 'k-'; else, ls = 'k--'; end
        plot([0 tEnd(k)], [e(j) e(j)], ls, 'LineWidth', 1.5);
    end

    y0 = [e + 0.05*max(1,e), e(e>0) - 0.05*max(1,e(e>0)), linspace(0.1, yTop, 8)];
    y0 = unique(y0(y0 > 0 & y0 <= yTop));

    tt = 0:h:tEnd(k);
    [t, soln] = ode45(@(t,y) f(y,b), tt, y0, options);
    for i = 1:length(y0)
        plot(t, soln(:,i), 'b-', 'LineWidth', 1);
    end

    axis([0 tEnd(k) 0 yTop]);
    xlabel('t (days)'); ylabel('y (hundreds of tribbles)');
    title(sprintf('Direction field, b = %.3g', b)); box on;
end
saveas(gcf, 'Q7_direction_fields.png');

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

Hs = @(t,y) p*y.^3 ./ (q + y.^3) .* abs(sin(pi*t/365));
years = 5;
h  = 0.5;
tt = 0:h:365*years;
options10 = odeset(options, 'MaxStep', 1);

figure(4); set(gcf, 'Position', [100 100 1300 750]);
for k = 1:3
    b  = bvals(k);
    y0 = unique([0.2 0.5 1 1.5 2 2.5 linspace(0.1, 1.1*a/b, 6)]);
    [t, soln] = ode45(@(t,y) a*y - b*y.^2 - Hs(t,y), tt, y0, options10);

    % top row: 5 years
    subplot(2,3,k); hold on;
    for i = 1:length(y0)
        plot(t/365, soln(:,i), 'LineWidth', 1);
    end
    for j = 1:numel(eqs{k})                 % constant-hunting equilibria (Q6)
        plot([0 years], [eqs{k}(j) eqs{k}(j)], 'k:');
    end
    xlabel('t (years)'); ylabel('y (hundreds of tribbles)');
    title(sprintf('Seasonal hunting, b = %.3g', b)); grid on; box on;
    xlim([0 years]);

    subplot(2,3,k+3); hold on;
    early = t <= 30;
    for i = 1:length(y0)
        plot(t(early), soln(early,i), 'LineWidth', 1);
    end
    xlabel('t (days)'); ylabel('y (hundreds of tribbles)');
    title(sprintf('First 30 days, b = %.3g', b)); grid on; box on;
end
saveas(gcf, 'Q10_seasonal_hunting.png');

figure(5);
b  = 0.05;
y0 = [0.5 1 1.5 2 2.5 5 10 15];
[t, soln] = ode45(@(t,y) a*y - b*y.^2 - Hs(t,y), tt, y0, options10);
subplot(2,1,1); hold on;
for i = 1:length(y0)
    plot(t, soln(:,i), 'LineWidth', 1);
end
ylabel('y (hundreds of tribbles)');
title('Seasonal hunting, b = 0.05'); grid on; xlim([tt(1) tt(end)]);
subplot(2,1,2);
plot(t, abs(sin(pi*t/365)), 'k-');
xlabel('t (days)'); ylabel('|sin(\pi t/365)|');
title('Hunting intensity (0 = midwinter, 1 = midsummer)'); grid on;
xlim([tt(1) tt(end)]);
saveas(gcf, 'Q10_seasonal_b005_detail.png');

fprintf('Done. Figures saved as PNG files in: %s\n', pwd);