clc,clear;close all;

mu = 10000;
d_q = 5*10^-3;                      % 既存手法の許容遅延 [s]
d_ii = 5*10^-3;                     % オフロードしないジョブの許容遅延 [s]
d_ij = 5*10^-3;                     % オフロードするジョブの許容遅延 [s]
tui = 2*10^-3;                      % ユーザとクラウドレットの伝搬遅延 [s]
tij = [0 0.75/10^3; 0.75/10^3 0];   % クラウドレット間の伝搬遅延 [s]

n = 2;

omega1 = 5*10^2;                    % ジョブの処理ごとの報酬
omega2 = 1*10^6;                    % ジョブのオフロードの罰とオフロードされたジョブの報酬
omega3 = 5*10^8;                    % dqを超えた時の罰と範囲内の時の報酬
zeta = 300;                         % ジョブごとの処理コスト
eta = 700;                          % クラウドレットの稼働コスト

num = 101;
varphix = linspace(0, 1, num);
varphiy = linspace(0, 1, num);
utility1_matrix = zeros(num);
utility2_matrix = zeros(num);
utility1_matrix_p = zeros(num);
utility2_matrix_p = zeros(num);
latency1_matrix = zeros(num);
latency2_matrix = zeros(num);
utilization1_matrix = zeros(num);
utilization2_matrix = zeros(num);
x_opt = [];
y_opt = [];
x_opt_p = [];
y_opt_p = [];
lambda_p = [];
c1c2dif = [];
c1c2dif_p = [];
u1 = []; u1_p = [];
u2 = []; u2_p = [];

syms varphi [n, n] positive;

for k=1:90
    lambda = [k*100, 9000];
    %利得関数の定義

    %utility1 = @(x, y) omega1*lambda(1)/mu + omega2*y*lambda(2)/mu - omega2*x*lambda(1)/mu - (zeta*((1-x)*lambda(1) + y*lambda(2))/mu + eta) - omega3*lambda(1)/mu*(tui - d_q + (1-x)/(mu - (1-x)*lambda(1) - y*lambda(2)) + x*(tij(1,2) + 1/(mu - x*lambda(1) - (1-y)*lambda(2))));
    %utility2 = @(x, y) omega1*lambda(2)/mu + omega2*x*lambda(1)/mu - omega2*y*lambda(2)/mu - (zeta*((1-y)*lambda(2) + x*lambda(1))/mu + eta) - omega3*lambda(2)/mu*(tui - d_q + (1-y)/(mu - (1-y)*lambda(2) - x*lambda(1)) + y*(tij(2,1) + 1/(mu - y*lambda(2) - (1-x)*lambda(1))));
    %utility1_prop = @(x, y) omega1*lambda(1)/mu + omega2*y*lambda(2)/mu - omega2*x*lambda(1)/mu - (zeta*((1-x)*lambda(1) + y*lambda(2))/mu + eta) - omega3*0.5*lambda(1)/mu*((tui - d_ii + 1/(mu - (1-x)*lambda(1) - y*lambda(2))) + (tui - d_ij + (tij(1,2) + 1/(mu - x*lambda(1) - (1-y)*lambda(2)))));
    %utility2_prop = @(x, y) omega1*lambda(2)/mu + omega2*x*lambda(1)/mu - omega2*y*lambda(2)/mu - (zeta*((1-y)*lambda(2) + x*lambda(1))/mu + eta) - omega3*0.5*lambda(2)/mu*((tui - d_ii + 1/(mu - (1-y)*lambda(2) - x*lambda(1))) + (tui - d_ij + (tij(2,1) + 1/(mu - y*lambda(2) - (1-x)*lambda(1)))));
    utility1 = @(x, y) - (tui - d_q + (1-x)/(mu - (1-x)*lambda(1) - y*lambda(2)) + x*(tij(1,2) + 1/(mu - x*lambda(1) - (1-y)*lambda(2))));
    utility2 = @(x, y) - (tui - d_q + (1-y)/(mu - (1-y)*lambda(2) - x*lambda(1)) + y*(tij(2,1) + 1/(mu - y*lambda(2) - (1-x)*lambda(1))));
    utility1_prop = @(x, y) - ((tui - d_ii + 1/(mu - (1-x)*lambda(1) - y*lambda(2))) + (tui - d_ij + (tij(1,2) + 1/(mu - x*lambda(1) - (1-y)*lambda(2)))));
    utility2_prop = @(x, y) - ((tui - d_ii + 1/(mu - (1-y)*lambda(2) - x*lambda(1))) + (tui - d_ij + (tij(2,1) + 1/(mu - y*lambda(2) - (1-x)*lambda(1)))));

    utility_min = [utility1(0,0), utility2(0,0)];
    utility_prop_min = [utility1_prop(0,0), utility2_prop(0,0)];
    latency1 = @(x, y) 1/(mu - lambda(1)*(1-x) - lambda(2)*y);
    latency2 = @(x, y) 1/(mu - lambda(2)*(1-y) - lambda(1)*x);
    utilization1 = @(x, y) lambda(1)*(1-x) + lambda(2)*y;
    utilization2 = @(x, y) lambda(2)*(1-y) + lambda(1)*x;
    for i=1:num
        for j=1:num
            utility1_matrix(i,j) = utility1(varphix(i), varphiy(j));
            utility2_matrix(i,j) = utility2(varphix(i), varphiy(j));
            utility1_matrix_p(i,j) = utility1_prop(varphix(i), varphiy(j));
            utility2_matrix_p(i,j) = utility2_prop(varphix(i), varphiy(j));
            latency1_matrix(i,j) = latency1(varphix(i), varphiy(j));
            latency2_matrix(i,j) = latency2(varphix(i), varphiy(j));
            utilization1_matrix(i,j) = utilization1(varphix(i), varphiy(j));
            utilization2_matrix(i,j) = utilization2(varphix(i), varphiy(j));

            if lambda(1)*(1-varphix(i))+lambda(2)*varphiy(j) > mu
                utility1_matrix(i,j) = 0;
                utility2_matrix(i,j) = 0;
                utility1_matrix_p(i,j) = 0;
                utility2_matrix_p(i,j) = 0;
            end
            if lambda(1)*varphix(i)+lambda(2)*(1-varphiy(j)) > mu
                utility1_matrix(i,j) = 0;
                utility2_matrix(i,j) = 0;
                utility1_matrix_p(i,j) = 0;
                utility2_matrix_p(i,j) = 0;
            end
            if utility1_matrix(i,j) < utility_min(1)
                utility1_matrix(i,j) = 0;
            end
            if utility2_matrix(i,j) < utility_min(2)
                utility2_matrix(i,j) = 0;
            end
            if utility1_matrix_p(i,j) < utility_prop_min(1)
                utility1_matrix_p(i,j) = 0;
            end
            if utility2_matrix_p(i,j) < utility_prop_min(2)
                utility2_matrix_p(i,j) = 0;
            end
        end
    end

    temp_y = 0;
    while 1
        temp = utility1_matrix(:,temp_y+1);
        [max_utility1,temp_x] = max(temp);
        temp = utility2_matrix(temp_x,:);
        [max_utility2,temp_y1] = max(temp);
        temp = utility1_matrix(:,temp_y1);
        [max_utility1_1, temp_x1] = max(temp);

        temp_y = temp_y + 1;
        if max_utility1_1 == max_utility1
            break;
        end
        if temp_y > num-1
            temp_y = 1;
            temp_x = 1;
            break;
        end
    end
    temp_y_p = 0;
    while 1
        temp = utility1_matrix_p(:, temp_y_p+1);
        [max_utility1_p, temp_x_p] = max(temp);
        temp = utility2_matrix_p(temp_x_p, :);
        [max_utility2_p, temp_y1_p] = max(temp);
        temp = utility1_matrix_p(:, temp_y1_p);
        [max_utility1_1_p, temp_x1_p] = max(temp);

        temp_y_p = temp_y_p + 1;

        if max_utility1_1_p == max_utility1_p
            break;
        end
        if temp_y_p > num-1
            temp_y_p = 1;
            temp_x_p = 1;
            break;
        end
    end
    
    x_opt = [x_opt (temp_x-1)*0.01];
    y_opt = [y_opt (temp_y-1)*0.01];
    x_opt_p = [x_opt_p (temp_x_p-1)*0.01];
    y_opt_p = [y_opt_p (temp_y_p-1)*0.01];
    lambda_p = [lambda_p lambda(1)];
    c1c2dif = [c1c2dif abs((lambda(2)*(1-(temp_y-1)*0.01)+lambda(1)*(temp_x-1)*0.01)/mu - (lambda(1)*(1-(temp_x-1)*0.01)+lambda(2)*(temp_y-1)*0.01)/mu)];
    c1c2dif_p = [c1c2dif_p abs((lambda(2)*(1-(temp_y_p-1)*0.01)+lambda(1)*(temp_x_p-1)*0.01)/mu - (lambda(1)*(1-(temp_x_p-1)*0.01)+lambda(2)*(temp_y_p-1)*0.01)/mu)];
    u1 = [u1 utility1_matrix(temp_x, temp_y)];
    u1_p = [u1_p utility1_matrix(temp_x_p, temp_y_p)];
    u2 = [u2 utility2_matrix_p(temp_x, temp_y)];
    u2_p = [u2_p utility2_matrix_p(temp_x_p, temp_y_p)];
    utility1_matrix = [];
    utility2_matrix = [];
    utility1_matrix_p = [];
    utility2_matrix_p = [];
end

f1 = figure; f2 = figure; f3 = figure; f4 = figure; f5 = figure;
figure(f1);
hold on
scatter(lambda_p/mu, x_opt(1:90), 100, 'filled', 'blue');
scatter(lambda_p/mu, x_opt_p(1:90), 100, '^', "MarkerEdgeColor", "#FF0000", "MarkerFaceColor", "#FF0000");xlabel('Initial utilization of cloudlet 1'); ylabel('Optimal offload fraction \phi_{12}^*');
l5 = legend('Conventional method [8]', 'Proposed method');
l5.FontSize = 10;
ylim([0 1]);
hold off
figure(f2);
hold on
scatter(lambda_p/mu, y_opt(1:90), 100, 'filled', 'blue');
scatter(lambda_p/mu, y_opt_p(1:90), 100, '^', "MarkerEdgeColor", "#FF0000", "MarkerFaceColor", "#FF0000");
xlabel('Initial utilization of cloudlet 1'); ylabel('Optimal offload fraction \phi_{21}^*');
l6 = legend('Conventional method [8]', 'Proposed method');
l6.FontSize = 10;
ylim([0 1]);
hold off
figure(f3);
hold on
scatter(lambda_p/mu, c1c2dif(1:90), 100, 'filled', 'blue');
scatter(lambda_p/mu, c1c2dif_p(1:90), 100, '^', "MarkerEdgeColor", "#FF0000", "MarkerFaceColor", "#FF0000");
xlabel('Initial utilization of cloudlet 1'); ylabel('Difference in utilization after offloading');
l6 = legend('Conventional method [8]', 'Proposed method');
l6.FontSize = 10;
hold off
figure(f4);
hold on
scatter(lambda_p/mu, u1(1:90), 100, 'filled', 'blue');
scatter(lambda_p/mu, u1_p(1:90), 100, '^', "MarkerEdgeColor", "#FF0000", "MarkerFaceColor", "#FF0000");
xlabel('Initial utilization of cloudlet 1'); ylabel('Utility value of cloudlet 1');
l6 = legend('Conventional method [8]', 'Proposed method');
l6.FontSize = 10;
hold off
figure(f5);
hold on
scatter(lambda_p/mu, u2(1:90), 100, 'filled', 'blue');
scatter(lambda_p/mu, u2_p(1:90), 100, '^', "MarkerEdgeColor", "#FF0000", "MarkerFaceColor", "#FF0000");
xlabel('Initial utilization of cloudlet 1'); ylabel('Utility value of cloudlet 2');
l6 = legend('Conventional method [8]', 'Proposed method');
l6.FontSize = 10;
hold off