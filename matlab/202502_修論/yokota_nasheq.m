clc,clear;close all;

num = 101;
maxlambda = 101;
initial = 1;
%% パラメータ設定
N = 2; % ユーザ数
lambda = [750, 0]; % 各ユーザのジョブ到着率
mu = 1000;b=1/mu; % 各ユーザの処理能力
d_rq = 5*10^-3; % 許容遅延
omega1 = 1;
omega2 = 1;
omega3 = 100;

%% 変数定義
%利得関数
utility1 = zeros(num); % クラウドレット１の利得行列
utility2 = zeros(num); % クラウドレット2の利得行列
utility1_0 = 0;
utility2_0 = 0;
latency_rq1 = zeros(num); % 許容遅延の項
latency_rq2 = zeros(num); % 許容遅延の項

%平均遅延
E_nooffload1 = zeros(num); % クラウドレット１のオフロードしないジョブの平均遅延
E_offload1 = zeros(num); % クラウドレット１のオフロードするジョブの平均遅延
E_nooffload2 = zeros(num); % クラウドレット２のオフロードしないジョブの平均遅延
E_offload2 = zeros(num); % クラウドレット２のオフロードするジョブの平均遅延
E_nopriority1 = zeros(num); % クラウドレット１のFCFSにおける遅延
E_nopriority2 = zeros(num); %クラウドレット２のFCFSにおける遅延

%到着率
lambda_offload1 = zeros(num); % クラウドレット１がオフロードするジョブの到着率
lambda_nooffload1 = zeros(num); % クラウドレット１がオフロードしないジョブの到着率
lambda_offload2 = zeros(num); % クラウドレット２がオフロードするジョブの到着率
lambda_nooffload2 = zeros(num); % クラウドレット２がオフロードしないジョブの到着率
lambda1 = zeros(num); %クラウドレット１のオフロード後の到着率
lambda2 = zeros(num); %クラウドレット２のオフロード後の到着率

%利用率
rho1 = zeros(num); %クラウドレット１の利用率
rho2 = zeros(num); %クラウドレット２の利用率

%オフロード割合
varphi = linspace(0, 1, num);
varphi_1 = 0;
varphi_2 = 0;

%% 利得の計算
% オフロード優先
lambda_p = [];
avarphi1_p = [];
avarphi2_p = [];
autility1_p = [];
autility2_p = [];
autility_mean_p = [];
aE_offload1_p = [];
aE_nooffload1_p = [];
aE_offload2_p = [];
aE_nooffload2_p = [];

% オフロード非優先
bvarphi1_p = [];
bvarphi2_p = [];
butility1_p = [];
butility2_p = [];
butility_mean_p = [];
utility1_0_p = [];
utility2_0_p = [];
bE_offload1_p = [];
bE_nooffload1_p = [];
bE_offload2_p = [];
bE_nooffload2_p = [];

% FCFS
cvarphi1_p = [];
cvarphi2_p = [];
cutility1_p = [];
cutility2_p = [];
cutility_mean_p = [];
ctility1_0_p = [];
ctility2_0_p = [];
cE_nopriority1_p = [];
cE_nopriority2_p = [];

%% オフロード優先
for n=initial:maxlambda
    lambda(2) = 10*(n-1);
    lambda_p = [lambda_p lambda(2)];

    for i=1:num
        varphi_1 = varphi(i);
        for j=1:num
            varphi_2 = varphi(j);

            lambda_offload1(i,j) = varphi_1*lambda(1);
            lambda_nooffload1(i,j) = (1-varphi_1)*lambda(1);
            lambda_offload2(i,j) = varphi_2*lambda(2);
            lambda_nooffload2(i,j) = (1-varphi_2)*lambda(2);
            lambda1(i,j) = (1-varphi_1)*lambda(1) + varphi_2*lambda(2);
            lambda2(i,j) = varphi_1*lambda(1) + (1-varphi_2)*lambda(2);
            rho1(i,j) = (1-varphi_1)*lambda(1)*b + varphi_2*lambda(2)*b;
            rho2(i,j) = varphi_1*lambda(1)*b + (1-varphi_2)*lambda(2)*b;
            if rho1(i,j) >= 1 || rho2(i,j) >= 1
                lambda1(i,j) = 0;
                lambda2(i,j) = 0;
                lambda_offload1(i,j) = 0;lambda_nooffload1(i,j) = 0;
                lambda_offload2(i,j) = 0;lambda_nooffload2(i,j) = 0;
                utility1(i,j) = -1;
                utility2(i,j) = -1;
                continue;
            end
            if varphi_1 == 0
                E_offload1(i,j) = 0;
                E_nooffload2(i,j) = 1/(1/b - lambda2(i,j));
                latency_rq1(i,j) = 0;
            else
                E_offload1(i,j) = (varphi_1*lambda(1)*b*b + (1-varphi_2)*lambda(2)*b*b)/((1-varphi_1*lambda(1)*b)) + b + 2.0*10^-3;
                E_nooffload2(i,j) = (varphi_1*lambda(1)*b*b + (1-varphi_2)*lambda(2)*b*b)/((1-varphi_1*lambda(1)*b)*(1-varphi_1*lambda(1)*b-(1-varphi_2)*lambda(2)*b)) + b;
                latency_rq1(i,j) = 1;
            end
            if varphi_2 == 0
                E_offload2(i,j) = 0;
                E_nooffload1(i,j) = 1/(1/b - lambda1(i,j));
                latency_rq2(i,j) = 0;
            else
                E_offload2(i,j) = ((1-varphi_1)*lambda(1)*b*b + varphi_2*lambda(2)*b*b)/((1-varphi_2*lambda(2)*b)) + b + 2.0*10^-3;
                E_nooffload1(i,j) = ((1-varphi_1)*lambda(1)*b*b + varphi_2*lambda(2)*b*b)/((1-varphi_2*lambda(2)*b)*(1-varphi_2*lambda(2)*b-(1-varphi_1)*lambda(1)*b)) + b;
                latency_rq2(i,j) = 1;
            end

            utility1(i,j) = (-omega1*lambda_offload1(i,j)+omega1*lambda_offload2(i,j)+omega2*lambda1(i,j))*b + omega3*(latency_rq1(i,j)*(d_rq-E_offload1(i,j))+(d_rq-E_nooffload1(i,j)));
            utility2(i,j) = (-omega1*lambda_offload2(i,j)+omega1*lambda_offload1(i,j)+omega2*lambda2(i,j))*b + omega3*(latency_rq2(i,j)*(d_rq-E_offload2(i,j))+(d_rq-E_nooffload2(i,j)));

        end
    end
    [pureNashEquilibria, payoffs]= findPureNashEquilibriaWithPayoffs(utility1, utility2);
    if isempty(payoffs)
        pureNashEquilibria = [0,0];
        payoffs = [0,0];
    end
    rowMeans = mean(payoffs,2);
    [truePayoffs, trueIndex] = max(rowMeans);
    pureNashEquilibria(1) = pureNashEquilibria(trueIndex, 1);
    pureNashEquilibria(2) = pureNashEquilibria(trueIndex, 2);
    payoffs(1) = payoffs(trueIndex, 1);
    payoffs(2) = payoffs(trueIndex, 2);

    if pureNashEquilibria(1) == 0
        aE_offload1_p = [aE_offload1_p 0];
        aE_nooffload1_p = [aE_nooffload1_p 0];
        aE_offload2_p = [aE_offload2_p 0];
        aE_nooffload2_p = [aE_nooffload2_p 0];
    else
        aE_offload1_p = [aE_offload1_p E_offload1(pureNashEquilibria(1),pureNashEquilibria(2))];
        aE_nooffload1_p = [aE_nooffload1_p E_nooffload1(pureNashEquilibria(1),pureNashEquilibria(2))];
        aE_offload2_p = [aE_offload2_p E_offload2(pureNashEquilibria(1),pureNashEquilibria(2))];
        aE_nooffload2_p = [aE_nooffload2_p E_nooffload2(pureNashEquilibria(1),pureNashEquilibria(2))];
    end
    utility1_0 = 100*(d_rq-(1/(mu-lambda(1))));
    utility2_0 = 100*(d_rq-(1/(mu-lambda(2))));


    avarphi1_p = [avarphi1_p pureNashEquilibria(1,1)];
    avarphi2_p = [avarphi2_p pureNashEquilibria(1,2)];
    autility1_p = [autility1_p payoffs(1,1)];
    autility2_p = [autility2_p payoffs(1,2)];
    autility_mean_p = [autility_mean_p mean(payoffs(1,:))];
    utility1_0_p = [utility1_0_p utility1_0];
    utility2_0_p = [utility2_0_p utility2_0];

    E_offload1 = zeros(num);
    E_offload2 = zeros(num);
    E_nooffload1 = zeros(num);
    E_nooffload2 = zeros(num);
end

%% オフロード非優先
for n=initial:maxlambda
    lambda(2) = 10*(n-1);

    for i=1:num
        varphi_1 = varphi(i);
        for j=1:num
            varphi_2 = varphi(j);
            lambda_offload1(i,j) = varphi_1*lambda(1);
            lambda_nooffload1(i,j) = (1-varphi_1)*lambda(1);
            lambda_offload2(i,j) = varphi_2*lambda(2);
            lambda_nooffload2(i,j) = (1-varphi_2)*lambda(2);
            lambda1(i,j) = (1-varphi_1)*lambda(1) + varphi_2*lambda(2);
            lambda2(i,j) = varphi_1*lambda(1) + (1-varphi_2)*lambda(2);
            rho1(i,j) = (1-varphi_1)*lambda(1)*b + varphi_2*lambda(2)*b;
            rho2(i,j) = varphi_1*lambda(1)*b + (1-varphi_2)*lambda(2)*b;
            if rho1(i,j) >= 1 || rho2(i,j) >= 1
                lambda1(i,j) = 0;
                lambda2(i,j) = 0;
                lambda_offload1(i,j) = 0;lambda_nooffload1(i,j) = 0;
                lambda_offload2(i,j) = 0;lambda_nooffload2(i,j) = 0;
                utility1(i,j) = -1;
                utility2(i,j) = -1;
                continue;
            end
            if varphi_1 == 0
                E_offload1(i,j) = 0;
                E_nooffload2(i,j) = 1/(1/b - lambda2(i,j));
                latency_rq1(i,j) = 0;
            else
                E_offload1(i,j) = (varphi_1*lambda(1)*b*b + (1-varphi_2)*lambda(2)*b*b)/((1-(1-varphi_2)*lambda(2)*b)*(1-(1-varphi_2)*lambda(2)*b-varphi_1*lambda(1)*b)) + b + 2.0*10^-3;
                E_nooffload2(i,j) = (varphi_1*lambda(1)*b*b + (1-varphi_2)*lambda(2)*b*b)/((1-(1-varphi_2)*lambda(2)*b)) + b;
                latency_rq1(i,j) = 1;
            end
            if varphi_2 == 0
                E_offload2(i,j) = 0;
                E_nooffload1(i,j) = 1/(1/b - lambda1(i,j));
                latency_rq2(i,j) = 0;
            else
                E_offload2(i,j) = ((1-varphi_1)*lambda(1)*b*b + varphi_2*lambda(2)*b*b)/((1-(1-varphi_1)*lambda(1)*b)*(1-(1-varphi_1)*lambda(1)*b-varphi_2*lambda(2)*b)) + b + 2.0*10^-3;
                E_nooffload1(i,j) = ((1-varphi_1)*lambda(1)*b*b + varphi_2*lambda(2)*b*b)/((1-(1-varphi_1)*lambda(1)*b)) + b;
                latency_rq2(i,j) = 1;
            end

            utility1(i,j) = (-omega1*lambda_offload1(i,j)+omega1*lambda_offload2(i,j)+omega2*lambda1(i,j))*b + omega3*(latency_rq1(i,j)*(d_rq-E_offload1(i,j))+(d_rq-E_nooffload1(i,j)));
            utility2(i,j) = (-omega1*lambda_offload2(i,j)+omega1*lambda_offload1(i,j)+omega2*lambda2(i,j))*b + omega3*(latency_rq2(i,j)*(d_rq-E_offload2(i,j))+(d_rq-E_nooffload2(i,j)));
        end
    end
    [pureNashEquilibria, payoffs]= findPureNashEquilibriaWithPayoffs(utility1, utility2);
    if isempty(payoffs)
        pureNashEquilibria = [0,0];
        payoffs = [0,0];
    end
    rowMeans = mean(payoffs,2);
    [truePayoffs, trueIndex] = max(rowMeans);
    pureNashEquilibria(1) = pureNashEquilibria(trueIndex, 1);
    pureNashEquilibria(2) = pureNashEquilibria(trueIndex, 2);
    payoffs(1) = payoffs(trueIndex, 1);
    payoffs(2) = payoffs(trueIndex, 2);
    if pureNashEquilibria(1) == 0
        bE_offload1_p = [bE_offload1_p 0];
        bE_nooffload1_p = [bE_nooffload1_p 0];
        bE_offload2_p = [bE_offload2_p 0];
        bE_nooffload2_p = [bE_nooffload2_p 0];
    else
        bE_offload1_p = [bE_offload1_p E_offload1(pureNashEquilibria(1),pureNashEquilibria(2))];
        bE_nooffload1_p = [bE_nooffload1_p E_nooffload1(pureNashEquilibria(1),pureNashEquilibria(2))];
        bE_offload2_p = [bE_offload2_p E_offload2(pureNashEquilibria(1),pureNashEquilibria(2))];
        bE_nooffload2_p = [bE_nooffload2_p E_nooffload2(pureNashEquilibria(1),pureNashEquilibria(2))];
    end

    bvarphi1_p = [bvarphi1_p pureNashEquilibria(1,1)];
    bvarphi2_p = [bvarphi2_p pureNashEquilibria(1,2)];
    butility1_p = [butility1_p payoffs(1,1)];
    butility2_p = [butility2_p payoffs(1,2)];
    butility_mean_p = [butility_mean_p mean(payoffs(1,:))];

    E_offload1 = zeros(num);
    E_offload2 = zeros(num);
    E_nooffload1 = zeros(num);
    E_nooffload2 = zeros(num);
end

%% FCFS
for n=initial:maxlambda
    lambda(2) = 10*(n-1);

    for i=1:num
        varphi_1 = varphi(i);
        for j=1:num
            varphi_2 = varphi(j);
            lambda_offload1(i,j) = varphi_1*lambda(1);
            lambda_nooffload1(i,j) = (1-varphi_1)*lambda(1);
            lambda_offload2(i,j) = varphi_2*lambda(2);
            lambda_nooffload2(i,j) = (1-varphi_2)*lambda(2);
            lambda1(i,j) = (1-varphi_1)*lambda(1) + varphi_2*lambda(2);
            lambda2(i,j) = varphi_1*lambda(1) + (1-varphi_2)*lambda(2);
            rho1(i,j) = (1-varphi_1)*lambda(1)*b + varphi_2*lambda(2)*b;
            rho2(i,j) = varphi_1*lambda(1)*b + (1-varphi_2)*lambda(2)*b;
            if rho1(i,j) >= 1 || rho2(i,j) >= 1
                lambda1(i,j) = 0;
                lambda2(i,j) = 0;
                lambda_offload1(i,j) = 0;lambda_nooffload1(i,j) = 0;
                lambda_offload2(i,j) = 0;lambda_nooffload2(i,j) = 0;
                utility1(i,j) = -1;
                utility2(i,j) = -1;
                continue;
            end
            if varphi_1 == 0
                latency_rq1(i,j) = 0;
            else
                latency_rq1(i,j) = 1;
            end
            if varphi_2 == 0
                latency_rq2(i,j) = 0;
            else
                latency_rq2(i,j) = 1;
            end
            E_nopriority1(i,j) =  1/(1/b - lambda1(i,j));
            E_nopriority2(i,j) =  1/(1/b - lambda2(i,j));
            E_offload1(i,j) = 0;
            E_offload2(i,j) = 0;
            utility1(i,j) = (-omega1*lambda_offload1(i,j)+omega1*lambda_offload2(i,j)+omega2*lambda1(i,j))*b + omega3*((d_rq-E_nopriority1(i,j))+latency_rq2(i,j)*(d_rq-E_nopriority2(i,j)));
            utility2(i,j) = (-omega1*lambda_offload2(i,j)+omega1*lambda_offload1(i,j)+omega2*lambda2(i,j))*b + omega3*((d_rq-E_nopriority2(i,j))+latency_rq1(i,j)*(d_rq-E_nopriority1(i,j)));
        end
    end
    [pureNashEquilibria, payoffs]= findPureNashEquilibriaWithPayoffs(utility1, utility2);
    if isempty(payoffs)
        pureNashEquilibria = [0,0];
        payoffs = [0,0];
    end
    rowMeans = mean(payoffs,2);
    [truePayoffs, trueIndex] = max(rowMeans);
    pureNashEquilibria(1) = pureNashEquilibria(trueIndex, 1);
    pureNashEquilibria(2) = pureNashEquilibria(trueIndex, 2);
    payoffs(1) = payoffs(trueIndex, 1);
    payoffs(2) = payoffs(trueIndex, 2);
    if pureNashEquilibria(1) == 0
        cE_nopriority1_p = [cE_nopriority1_p 0];
        cE_nopriority2_p = [cE_nopriority2_p 0];
    else
        cE_nopriority1_p = [cE_nopriority1_p E_nopriority1(pureNashEquilibria(1),pureNashEquilibria(2))];
        cE_nopriority2_p = [cE_nopriority2_p E_nopriority2(pureNashEquilibria(1),pureNashEquilibria(2))];
    end

    cvarphi1_p = [cvarphi1_p pureNashEquilibria(1,1)];
    cvarphi2_p = [cvarphi2_p pureNashEquilibria(1,2)];
    cutility1_p = [cutility1_p payoffs(1,1)];
    cutility2_p = [cutility2_p payoffs(1,2)];
    cutility_mean_p = [cutility_mean_p mean(payoffs(1,:))];

    E_nopriority1 = zeros(num);
    E_nopriority2 = zeros(num);
end

%% グラフ表示
f1 = figure;
f2 = figure;
f3 = figure;

figure(f1);
hold on
grid on
scatter(lambda_p/mu, avarphi1_p(1:maxlambda-initial+1)*0.01,75, 'x', 'red');
scatter(lambda_p/mu, avarphi2_p(1:maxlambda-initial+1)*0.01, 75, 'filled', 'red');
scatter(lambda_p/mu, bvarphi1_p(1:maxlambda-initial+1)*0.01,75, 'x', 'blue');
scatter(lambda_p/mu, bvarphi2_p(1:maxlambda-initial+1)*0.01, 75, 'filled', 'blue');
scatter(lambda_p/mu, cvarphi1_p(1:maxlambda-initial+1)*0.01,75, 'x', 'green');
scatter(lambda_p/mu, cvarphi2_p(1:maxlambda-initial+1)*0.01, 75, 'filled', 'green');
ylim([0 1]);
xlabel('Initial Utilization \rho_2 of Cloudlet 2'); ylabel('Nash Equilibrium strategy');
l1 = legend('OP: Offload fraction of cloudlet 1', 'OP: Offload fraction of cloudlet 2', 'ONP: Offload fraction of cloudlet 1', 'ONP: Offload fraction of cloudlet 2','FCFS: Offload fraction of cloudlet 1', 'FCFS: Offload fraction of cloudlet 1', 'FCFS: Offload fraction of cloudlet 2');
%l1.FontSize = 10;
xlim([0.75 1]); %0.77
hold off

figure(f2);
hold on
%grid on
scatter(lambda_p/mu, autility_mean_p(1:maxlambda-initial+1), 20, 'filled', 'red');
%scatter(lambda_p/mu, autility1_p(1:maxlambda-initial+1), 40, 'x', 'red');
%scatter(lambda_p/mu, autility2_p(1:maxlambda-initial+1), 40, '^', 'red');

scatter(lambda_p/mu, butility_mean_p(1:maxlambda-initial+1), 20, 'filled', 'blue');
%scatter(lambda_p/mu, butility1_p(1:maxlambda-initial+1), 40, 'x', 'blue');
%scatter(lambda_p/mu, butility2_p(1:maxlambda-initial+1), 40, '^', 'blue');
    
scatter(lambda_p/mu, cutility_mean_p(1:maxlambda-initial+1), 20, 'filled', 'green');
%scatter(lambda_p/mu, cutility1_p(1:maxlambda-initial+1), 40, 'x', 'green');
%scatter(lambda_p/mu, cutility2_p(1:maxlambda-initial+1), 40, '^', 'green');
ylim([0 1.0]);
xlim([0.75 1]); %0.77
xlabel('Initial Utilization \rho_2 of Cloudlet 2'); ylabel('Average utility $\bar{U}$', 'Interpreter', 'latex');
l1 = legend('OP', 'ONP', 'FCFS');
l1.FontSize = 10;
hold off

figure(f3);
hold on
grid on
scatter(lambda_p/mu, aE_nooffload1_p(1:maxlambda-initial+1), 40, 'filled', 'red');
scatter(lambda_p/mu, aE_offload1_p(1:maxlambda-initial+1), 40, 'x', 'red');
scatter(lambda_p/mu, aE_nooffload2_p(1:maxlambda-initial+1), 40, '^', 'red');
scatter(lambda_p/mu, aE_offload2_p(1:maxlambda-initial+1), 40, '*', 'red');

scatter(lambda_p/mu, bE_nooffload1_p(1:maxlambda-initial+1), 40, 'filled', 'blue');
scatter(lambda_p/mu, bE_offload1_p(1:maxlambda-initial+1), 40, 'x', 'blue');
scatter(lambda_p/mu, bE_nooffload2_p(1:maxlambda-initial+1), 40, '^', 'blue');
scatter(lambda_p/mu, bE_offload2_p(1:maxlambda-initial+1), 40, '*', 'blue');

scatter(lambda_p/mu, cE_nopriority1_p(1:maxlambda-initial+1), 40, 'filled', 'green');
scatter(lambda_p/mu, cE_nopriority2_p(1:maxlambda-initial+1), 40, 'x', 'green');
ylim([-inf inf]);
xlabel('Initial Utilization \rho_2 of Cloudlet 2'); ylabel('Average Latency');
l1 = legend('OP: Not offloaded jobs at cloudlet 1', 'OP: Offloaded jobs at cloudlet 1', 'OP: Not offloaded jobs at cloudlet 2', 'OP: Offloaded jobs at cloudlet 2','ONP: Not offloaded jobs at cloudlet 1', 'ONP: Offloaded jobs at cloudlet 1', 'ONP: Not offloaded jobs at cloudlet 2', 'ONP: Offloaded jobs at cloudlet 2','FCFS: Jobs at cloudlet 1', 'FCFS: Jobs at cloudlet 2');

l1.FontSize = 10;
xlim([0.75 1]); %0.77
hold off
%% ナッシュ均衡探索

function [pureNashEquilibria, payoffs] = findPureNashEquilibriaWithPayoffs(U1, U2)
    % findPureNashEquilibriaWithPayoffs - 純粋戦略ナッシュ均衡とその利得を探索
    %
    % 使用法:
    %   [pureNashEquilibria, payoffs] = findPureNashEquilibriaWithPayoffs(U1, U2)
    %
    % 入力:
    %   U1 - プレイヤー1の利得行列 (m x n 行列)
    %   U2 - プレイヤー2の利得行列 (m x n 行列)
    %
    % 出力:
    %   pureNashEquilibria - 純粋戦略ナッシュ均衡のリスト (k x 2 行列)
    %                        各行は [プレイヤー1の戦略, プレイヤー2の戦略] を表す
    %                        均衡が存在しない場合は空行列 []
    %   payoffs - ナッシュ均衡における利得のリスト (k x 2 行列)
    %             各行は [プレイヤー1の利得, プレイヤー2の利得] を表す

    % 利得行列のサイズを取得
    [numStrategies1, numStrategies2] = size(U1);

    % 結果のリストを初期化
    pureNashEquilibria = [];
    payoffs = [];

    % 各戦略の組み合わせを評価
    for i = 1:numStrategies1
        for j = 1:numStrategies2
            % プレイヤー1の最適反応を確認
            if U1(i, j) >= max(U1(:, j))
                % プレイヤー2の最適反応を確認
                if U2(i, j) >= max(U2(i, :))
                    % 純粋戦略ナッシュ均衡を追加
                    pureNashEquilibria = [pureNashEquilibria; i, j];
                    payoffs = [payoffs; U1(i, j), U2(i, j)];
                end
            end
        end
    end

    % 結果の表示
    % if isempty(pureNashEquilibria)
    %     disp('純粋戦略ナッシュ均衡は存在しません。');
    % else
    %     disp('純粋戦略ナッシュ均衡:');
    %     for k = 1:size(pureNashEquilibria, 1)
    %         fprintf('戦略: プレイヤー1 = %d, プレイヤー2 = %d | 利得: プレイヤー1 = %d, プレイヤー2 = %d\n', ...
    %             pureNashEquilibria(k, 1), pureNashEquilibria(k, 2), payoffs(k, 1), payoffs(k, 2));
    %     end
    % end
end
