clear all;
close all;

% Causal bounds over arms
h = [0.6506, 0.7527];
l = [0.5905, 0.4839];
% Expected rewards of arms
miu = [0.6201, 0.6948];

T = 5000;
K = 100;
N = length(h);

e = 1e-15;

opts = zeros(1, T);
vals = zeros(1, T);
for k = 1 : K
    counts = zeros(N, 2); 
    ysum = 0;
    k
    for t = 1 : T
        x = 0;
        if t <= N
            x = t;
        else
            value = zeros(1, N);
            
            % Finding the kl-UCB bound subject to causal bound constraints.
            for i = 1 : N
                theta = (counts(i, 1)) / (sum(counts(i, :)));
                d = (log(t) + 3 * log(log(t)))/double(sum(counts(i,:)));
                % d is f(t)/Nx(t) in line 4 of algorithm 2

                ll = theta;
                uu = h(i);

                while uu - ll > e
                    m = (ll + uu) / 2;
                    if fun1(m, theta) > d %kl() should be fun1 m is mu* theta is mu_x
                        uu = m;
                    else
                        ll = m;
                    end
                end
                value(i) = (ll + uu) / 2;
            end
            
            [val, x] = max(value);
        end
        
        y = binornd(1, miu(x));

        counts(x, 1) = counts(x, 1) + y;
        counts(x, 2) = counts(x, 2) + 1 - y;
        
        [vopt, opt] = max(miu);
        
        ysum = ysum + miu(opt) - miu(x);

        opts(t) = opts(t) + (x == opt);
        vals(t) = vals(t) + ysum;
    end
end

function [ val ] = fun1(v, theta)
    if theta == 0
        val = log(1 / (1 - v));
    elseif theta == 1
        val = log(1 / v);
    else
        val = theta * log(theta / v) + (1 - theta) * log((1 - theta)/(1 - v));
    end %theta is mu_x and v is mu*
end

for k = 1 : K
    counts = zeros(N, 2);
    ysum = 0;
    k
    for t = 1 : T
        x = 0;
        if t <= N
            x = t;
        else
            value = zeros(1, N);
            
            % Finding the kl-UCB bound subject to causal bound constraints.
            for i = 1 : N
                theta = (counts(i, 1)) / (sum(counts(i, :)));
                d = (log(t) + 3 * log(log(t)))/double(sum(counts(i,:)));

                ll = theta;
                uu = h(i);
                
                while uu - ll > e
                    m = (ll + uu) / 2;
                    if fun1(m, theta) > d
                        uu = m;
                    else
                        ll = m;
                    end
                end
                value(i) = (ll + uu) / 2;
            end
            
            [val, x] = max(value);
        end
        
        y = binornd(1, miu(x));

        counts(x, 1) = counts(x, 1) + y;
        counts(x, 2) = counts(x, 2) + 1 - y;
        
        [vopt, opt] = max(miu);
        
        ysum = ysum + miu(opt) - miu(x);

        opts(t) = opts(t) + (x == opt);
        vals(t) = vals(t) + ysum;
    end
end

opts2 = zeros(1, T);
vals2 = zeros(1, T);
for k = 1 : K
    counts = zeros(N, 2);
    ysum = 0;
    k
    for t = 1 : T
        x = 0;
        if t <= N
            x = t;
        else
            value = zeros(1, N);

            % Finding the kl-UCB bound
            for i = 1 : N
                theta = (counts(i, 1)) / (sum(counts(i, :)));
                d = (log(t) + 3 * log(log(t)))/double(sum(counts(i,:)));

                ll = theta;
                uu = 1.0;
                
                while uu - ll > e
                    m = (ll + uu) / 2;
                    if fun1(m, theta) > d
                        uu = m;
                    else
                        ll = m;
                    end
                end
                value(i) = (ll + uu) / 2;
            end
            
            [val, x] = max(value);
        end
        
        y = binornd(1, miu(x));

        counts(x, 1) = counts(x, 1) + y;
        counts(x, 2) = counts(x, 2) + 1 - y;
        
        [vopt, opt] = max(miu);
        
        ysum = ysum + miu(opt) - miu(x);

        opts2(t) = opts2(t) + (x == opt);
        vals2(t) = vals2(t) + ysum;
    end
end

% Plot the cumulative regrets
h = figure(1);
set(h, 'Position', [100, 100, 1200, 1000]);

t1 = [1:T];

opts = opts ./ K;
plot(t1, opts(t1), 'LineStyle', '-', 'LineWidth', 2);
hold on;
opts2 = opts2 ./ K;
plot(t1, opts2(t1), 'LineStyle', '-', 'LineWidth', 2);
hold on;


axis square
set(gca,'FontSize',40);
title('Probability of Optimal Action');
ylabel('Probability');
xlabel('Trials');
l = legend({'b-klucb', 'klucb'}, 'Location', 'SouthEast');
set(l, 'FontSize', 50);

regretAxis = [0, T, 0, 50];
h = figure(2);
set(h, 'Position', [100, 100, 1200, 1000]);

t2 = [[1:500:T],T];

vals = vals ./ K;
plot(t2, vals(t2), 'LineStyle', '-', 'Marker', 'o', 'LineWidth', 5);
hold on;
vals2 = vals2 ./ K;
plot(t2, vals2(t2), 'LineStyle', '-.', 'Marker', '>', 'LineWidth', 5);
hold on;

axis square
axis(regretAxis);
set(gca,'FontSize',40);
title('Regret');
ylabel('Cumulative Regret');
xlabel('Trials');

l = legend({'b-klucb', 'klucb'}, 'Location', 'SouthEast');
set(l, 'FontSize', 50);

