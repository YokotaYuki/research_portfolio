clc,clear;close all;

n = 2;                              % 競争するクラウドレットの数
w = 0.1;                            % ステップサイズ
t = 0; t_p = 0;                     % 反復回数
tolerance = 10^-3;                  % NEを決めるための許容値
maxiter = 10^3;                     % maximum no. of iterations;

mu_ii = 10000;                      % サービス率 [jobs/s]
lambda = [0 0];                     % 到着率の定義 [jobs/s]
lambda_p = [];                      % 出力用の箱

d_q = 5*10^-3;                      % 既存手法の許容遅延 [s]
d_ii = 5*10^-3;                     % オフロードしないジョブの許容遅延 [s]
d_ij = 5*10^-3;                     % オフロードするジョブの許容遅延 [s]
tui = 2*10^-3;                      % ユーザとクラウドレットの伝搬遅延 [s]
tij = [0 0.75/10^3; 0.75/10^3 0];   % クラウドレット間の伝搬遅延 [s]

omega1 = 5*10^2;                    % ジョブの処理ごとの報酬
omega2 = 1*10^6;                    % ジョブのオフロードの罰とオフロードされたジョブの報酬
omega3 = 5*10^8;                    % dqを超えた時の罰と範囲内の時の報酬
zeta = 300;                         % ジョブごとの処理コスト
eta = 700;                          % クラウドレットの稼働コスト

%出力用の各値の箱の定義
latency1_dif = []; latency2_dif = [];latency1_dif_p = []; latency2_dif_p = [];
latency11 = []; latency12 = []; latency21 = []; latency22 = [];
latency11_prop = []; latency12_prop = []; latency21_prop = []; latency22_prop = [];

x = []; y = []; x_prop = []; y_prop = []; x_s = []; y_s = []; x_prop_s = []; y_prop_s = [];
x_p = []; y_p = [];
xx_p = []; yy_p = [];
c1 = []; c2 = []; c1_prop = []; c2_prop = [];
c1c2dif = []; c1c2dif_prop = [];
temp12 = 0;temp21 = 0;
temp12_p = 0;temp21_p = 0;
varphi12 = 0;varphi21 = 0;
varphi12_p = 0; varphi21_p = 0;
utility_p = []; utility_pprop = [];
latency_var1 = [];latency_var2 = [];
latency_var1_prop = [];latency_var2_prop = [];

maxlambda = 90;                      % 到着率の最大値

utility1_values = zeros(1,maxlambda+1);
utility2_values = zeros(1,maxlambda+1);
utility1_p_values = zeros(1,maxlambda+1);
utility2_p_values = zeros(1,maxlambda+1);

utility1_comparison = zeros(1, maxlambda+1);
utility2_comparison = zeros(1,maxlambda+1);
utility1_p_comparison = zeros(1,maxlambda+1);
utility2_p_comparison = zeros(1,maxlambda+1);

utility_sum = zeros(1, maxlambda+1);
utility_p_sum = zeros(1, maxlambda+1);

keta = 1001;

utility1_matrix = zeros(keta);
utility2_matrix = zeros(keta);
utility1_matrix_prop = zeros(keta);
utility2_matrix_prop = zeros(keta);


%アルゴリズムの開始
for k=0:+1:maxlambda
    lambda = [k*100, 9000];
    utility1 = @(x, y) omega1*lambda(1)/mu_ii + omega2*y*lambda(2)/mu_ii - omega2*x*lambda(1)/mu_ii - omega3*lambda(1)/mu_ii*(tui - d_q + (1-x)/(mu_ii - (1-x)*lambda(1) - y*lambda(2)) + x*(tij(1,2) + 1/(mu_ii - x*lambda(1) - (1-y)*lambda(2))));
    utility2 = @(x, y) omega1*lambda(2)/mu_ii + omega2*x*lambda(1)/mu_ii - omega2*y*lambda(2)/mu_ii - omega3*lambda(2)/mu_ii*(tui - d_q + (1-y)/(mu_ii - (1-y)*lambda(2) - x*lambda(1)) + y*(tij(2,1) + 1/(mu_ii - y*lambda(2) - (1-x)*lambda(1))));
    utility1_prop = @(x, y) omega1*lambda(1)/mu_ii + omega2*y*lambda(2)/mu_ii - omega2*x*lambda(1)/mu_ii - omega3*lambda(1)/mu_ii*((tui - d_ii + 1/(mu_ii - (1-x)*lambda(1) - y*lambda(2))) + (tui - d_ij + (tij(1,2) + 1/(mu_ii - x*lambda(1) - (1-y)*lambda(2)))));
    utility2_prop = @(x, y) omega1*lambda(2)/mu_ii + omega2*x*lambda(1)/mu_ii - omega2*y*lambda(2)/mu_ii - omega3*lambda(2)/mu_ii*((tui - d_ii + 1/(mu_ii - (1-y)*lambda(2) - x*lambda(1))) + (tui - d_ij + (tij(2,1) + 1/(mu_ii - y*lambda(2) - (1-x)*lambda(1)))));
    utility_min = [utility1(0,0), utility2(0,0)];
    utility_prop_min = [utility1_prop(0,0), utility2_prop(0,0)];
    
    latency = @(x, y)[(tui + 1/(mu_ii - (1-x)*lambda(1) - y*lambda(2))), (tui + tij(1,2) + 1/(mu_ii - x*lambda(1) - (1-y)*lambda(2)));...
                      (tui + tij(2,1) + 1/(mu_ii - y*lambda(2) - (1-x)*lambda(1))), (tui + 1/(mu_ii - (1-y)*lambda(2) - x*lambda(1)))];
%    if lambda(1) == 8800
%        disp(lambda);
%    end
    for i=1:+1:keta
        %利得行列の要素を戦略と対応
        xx(1) = 0.0 + 0.001*(i-1); %クラウドレット1の戦略

        for j=1:+1:keta
            %利得行列の要素を戦略と対応        
            xx(2) = 0.0 + 0.001*(j-1); %クラウドレット2の戦略
            %条件を満たす要素のみ計算、満たさないものはスキップ
            if mu_ii - lambda(1)*(1-xx(1)) - lambda(2)*xx(2) < 0
                continue;
            end
            if mu_ii - lambda(2)*(1-xx(2)) - lambda(1)*xx(1) < 0
                continue;
            end
            %利得行列の計算
            utility1_matrix(i,j) = utility1(xx(1), xx(2));
            utility2_matrix(i,j) = utility2(xx(1), xx(2));
            if utility1_matrix(i,j) < utility_min(1)
                utility1_matrix(i,j) = utility_min(1);
            end
            if utility2_matrix(i,j) < utility_min(2)
                utility2_matrix(i,j) = utility_min(2);
            end
            % if (2*(lambda(1)^2)*lambda(2)^2)/(mu_ii - (1-xx(1))*lambda(1) - xx(2)*lambda(2))^3 + (2*(lambda(1)^2)*lambda(2)^2)/(mu_ii - xx(1)*lambda(1) - (1-xx(2))*lambda(2))^3 < 2*(lambda(1)*lambda(2))^3 / (mu_ii^2 - ((1-xx(1))^2)*lambda(1)^2 - 2*mu_ii*(1-xx(1))*lambda(1)*lambda(2) - (xx(2)^2)*lambda(2)^2 + 2*(1-xx(1))*xx(2)*lambda(2) + (xx(1)^2)*lambda(1)^2 - 2*mu_ii*xx(1)*lambda(1)*(1-xx(2)) - (1-xx(2))^2*lambda(2)^2)^(3/2)
            %     continue;
            % end
            % if 2*lambda(1)^2 / (mu_ii - (1-xx(1))*lambda(1) - xx(2)*lambda(2))^3 + 2*lambda(1)^2 / (mu_ii - xx(1)*lambda(1) - (1-xx(2))*lambda(2))^3 < 2*lambda(1)^3 / (mu_ii^2 - ((1-xx(1))^2)*lambda(1)^2 - 2*mu_ii*(1-xx(1))*lambda(1)*lambda(2) - (xx(2)^2)*lambda(2)^2 + 2*(1-xx(1))*lambda(1)*xx(2)*lambda(2) + (xx(1)^2)*lambda(1)^2 - 2*mu_ii*xx(1)*lambda(1)*(1-xx(2)) - ((1-xx(2))^2)*lambda(2)^2)^(3/2)
            %     continue;
            % end
            % if 2*lambda(2)^2 / (mu_ii - (1-xx(1))*lambda(1) - xx(2)*lambda(2))^3 + 2*lambda(2)^2 / (mu_ii - xx(1)*lambda(1) - (1-xx(2))*lambda(2))^3 < 2*lambda(2)^3 / (mu_ii^2 - ((1-xx(1))^2)*lambda(1)^2 - 2*mu_ii*(1-xx(1))*lambda(1)*lambda(2) - (xx(2)^2)*lambda(1)^2 + 2*(1-xx(1))*lambda(1)*xx(2)*lambda(2) + (xx(1)^2)*lambda(1)^2 - 2*mu_ii*xx(1)*lambda(1)*(1-xx(2)) - ((1-xx(2))^2)*lambda(2)^2)^(3/2)
            %     continue;
            % end
            %利得行列の計算
            utility1_matrix_prop(i,j) = utility1_prop(xx(1), xx(2));
            utility2_matrix_prop(i,j) = utility2_prop(xx(1), xx(2));
            if utility1_matrix_prop(i,j) < utility_prop_min(1)
                utility1_matrix_prop(i,j) = utility_prop_min(1);
            end
            if utility2_matrix_prop(i,j) < utility_prop_min(2)
                utility2_matrix_prop(i,j) = utility_prop_min(2);
            end
        end
    end
   %ナッシュ均衡の導出
   %ランダムな要素を選択
   r = randi([1 101], 1);
   while(1)
       [M2_1, I2_1] = max(utility2_matrix(r,:));
       [M1_1, I1_1] = max(utility1_matrix(:,I2_1));
       [M2_2, I2_2] = max(utility2_matrix(I1_1,:));
       [M1_2, I1_2] = max(utility1_matrix(:,I2_2));
       if I1_2 == I1_1
           break;
       end
       r = I1_2; %randi([1 1001], 1);
       t = t + 1;
       if t == 10000
           disp(I1_2);
           disp(I2_2);
           disp(r);
           t = 0;
           break;
       end
   end

   r = randi([1 101], 1);
   while(1)
       [M2_1_prop, I2_1_prop] = max(utility2_matrix_prop(r,:));
       [M1_1_prop, I1_1_prop] = max(utility1_matrix_prop(:,I2_1_prop));
       [M2_2_prop, I2_2_prop] = max(utility2_matrix_prop(I1_1_prop,:));
       [M1_2_prop, I1_2_prop] = max(utility1_matrix_prop(:,I2_2_prop));
       if I1_2_prop == I1_1_prop
           break;
       end
       r = I1_2_prop;
       t_p = t_p + 1;
       if t_p == 10000
           disp(I1_2_prop);
           disp(I2_2_prop);
           disp(mu_ii);
           t_p = 0;
           break;
       end
   end
    varphi12 = 0.0 + 0.001*(I1_2); varphi21 = 0.0 + 0.001*(I2_2);
    varphi12_p = 0.0 + 0.001*(I1_2_prop); varphi21_p = 0.0 + 0.001*(I2_2_prop);
    utility1_values(k+1) = M1_2;
    utility2_values(k+1) = M2_2;
    utility1_p_values(k+1) = M1_2_prop;
    utility2_p_values(k+1) = M2_2_prop;
    utility1_comparison(k+1) = utility1_matrix_prop(I1_2, I2_2);
    utility2_comparison(k+1) = utility2_matrix_prop(I1_2, I2_2);
    utility1_p_comparison(k+1) = M1_2_prop;
    utility2_p_comparison(k+1) = M2_2_prop;
    utility_sum(k+1) = utility1_matrix_prop(I1_2, I2_2) + utility2_matrix_prop(I1_2, I2_2);
    utility_p_sum(k+1) = M1_2_prop + M2_2_prop;

    %導出したオフロード割合の格納とその時の遅延時間と処理数の計算
    x_s = [x_s varphi12]; y_s = [y_s varphi21];
    x_prop_s = [x_prop_s varphi12_p]; y_prop_s = [y_prop_s varphi21_p];
    c1 = [c1 lambda(1)*(1-varphi12)+varphi21*lambda(2)]; c2 = [c2 lambda(2)*(1-varphi21)+lambda(1)*varphi12];
    c1_prop = [c1_prop lambda(1)*(1-varphi12_p)+varphi21_p*lambda(2)]; c2_prop = [c2_prop lambda(2)*(1-varphi21_p)+lambda(1)*varphi12_p];
    c1c2dif = [c1c2dif abs((lambda(1)*(1-varphi12)+varphi21*lambda(2))/mu_ii-(lambda(2)*(1-varphi21)+lambda(1)*varphi12)/mu_ii)];
    c1c2dif_prop = [c1c2dif_prop abs((lambda(1)*(1-varphi12_p)+varphi21_p*lambda(2))/mu_ii-(lambda(2)*(1-varphi21_p)+lambda(1)*varphi12_p)/mu_ii)];
    lambda_p = [lambda_p lambda(1)];
    latency1 = latency(varphi12,varphi21);latency1_p = latency(varphi12_p,varphi21_p);
    if lambda(1)==0
        latency1 = [0,0;latency1(2,1), latency1(2,2)];
        latency1_p = [0,0;latency1_p(2,1), latency1_p(2,2)];
    end
    latency11 = [latency11 latency1(1,1)]; latency12 = [latency12 latency1(1,2)]; latency21 = [latency21 latency1(2,1)]; latency22 = [latency22 latency1(2,2)];
    latency11_prop = [latency11_prop latency1_p(1,1)];latency12_prop = [latency12_prop latency1_p(1,2)];latency21_prop = [latency21_prop latency1_p(2,1)];latency22_prop = [latency22_prop latency1_p(2,2)];
    latency1_dif = [latency1_dif abs(latency1(1,1)-latency1(1,2))]; latency2_dif = [latency2_dif abs(latency1(2,1)-latency1(2,2))];
    latency1_dif_p = [latency1_dif_p abs(latency1_p(1,1)-latency1_p(1,2))]; latency2_dif_p = [latency2_dif_p abs(latency1_p(2,1)-latency1_p(2,2))];

% 箱の中身を空にして到着率を変化させて繰り返し
    varphi12 = 0; varphi21 = 0; varphi12_p = 0; varphi21_p = 0;
    utility1_matrix = zeros(keta);
    utility2_matrix = zeros(keta);
    utility1_matrix_prop = zeros(keta);
    utility2_matrix_prop = zeros(keta);
end

% グラフの定義
f1 = figure;
f2 = figure;
f3 = figure;
f4 = figure;
f5 = figure;
f6 = figure;
f7 = figure;
f8 = figure;
f9 = figure;
f10 = figure;
f11 = figure;
f12 = figure;
f13 = figure;
f14 = figure;
% 結果の出力
figure(f1);
hold on
scatter(lambda_p,latency11(1:maxlambda+1),100, '+', 'blue');
scatter(lambda_p,latency12(1:maxlambda+1),100, 'x', 'blue');
scatter(lambda_p,latency11_prop(1:maxlambda+1),100, 'o', 'red');
scatter(lambda_p,latency12_prop(1:maxlambda+1),100, 's', "MarkerEdgeColor", "#FF0000");  
xlabel('Average arrival rate of Cloudlet 1 \lambda_{1} [jobs/s]'); ylabel('Latency of Cloudlet 1 [s]');
l1 = legend('[Mondal+, OJ-COMS, 2020]: T_{11}', '[Mondal+, 2020]: T_{12}', 'Proposed method: T_{11}', 'Proposed method: T_{12}');
l1.FontSize = 10;
hold off
figure(f2);
hold on
scatter(lambda_p,latency21(1:maxlambda+1),100, '+', 'blue');
scatter(lambda_p,latency22(1:maxlambda+1),100, 'x', 'blue');
scatter(lambda_p,latency21_prop(1:maxlambda+1),100, 'o', 'red');
scatter(lambda_p,latency22_prop(1:maxlambda+1),100, 's', "MarkerEdgeColor", "#FF0000");
xlabel('Average arrival rate of Cloudlet 1 \lambda_{1} [jobs/s]'); ylabel('Latency of Cloudlet 2 [s]');
l2 = legend('[Mondal+, OJ-COMS, 2020]: T_{21}', '[Mondal+, 2020]: T_{22}', 'Proposed method: T_{21}', 'Proposed method: T_{22}');
l2.FontSize = 10;
hold off
figure(f3);
hold on
scatter(lambda_p(1:maxlambda+1),latency1_dif(1:maxlambda+1),100, 'filled', 'blue');
scatter(lambda_p(1:maxlambda+1),latency1_dif_p(1:maxlambda+1),100, '^', "MarkerEdgeColor", "#FF0000", "MarkerFaceColor", "#FF0000");
xlabel('Average arrival rate of Cloudlet 1 \lambda_{1} [jobs/s]'); ylabel('Difference between latency |T_{11}-T_{12}|');
l3 = legend('[Mondal+, OJ-COMS, 2020]', 'Proposed method');
l3.FontSize = 10;
hold off
figure(f4);
hold on
scatter(lambda_p(1:maxlambda+1),latency2_dif(1:maxlambda+1),100, 's', "MarkerEdgeColor", "#005AFF", "MarkerFaceColor", "#005AFF");
scatter(lambda_p(1:maxlambda+1),latency2_dif_p(1:maxlambda+1),100, 'v', "MarkerEdgeColor", "#FF4B00", "MarkerFaceColor", "#FF4B00");
xlabel('Average arrival rate of Cloudlet 1 \lambda_{1} [jobs/s]'); ylabel('Difference between latency |T_{21}-T_{22}|');
l4 = legend('[Mondal+, OJ-COMS, 2020]', 'Proposed method');
l4.FontSize = 10;
hold off
figure(f5);
hold on
scatter(lambda_p,x_s(1:maxlambda+1),100, 'filled', 'blue');
scatter(lambda_p,x_prop_s(1:maxlambda+1),100, '^', "MarkerEdgeColor", "#FF0000", "MarkerFaceColor", "#FF0000");
xlabel('Average arrival rate of Cloudlet 1 \lambda_{1} [jobs/s]'); ylabel('Offload fraction \phi_{12}');
l5 = legend('[Mondal+, OJ-COMS, 2020]', 'Proposed method');
l5.FontSize = 10;
hold off
figure(f6);
hold on
scatter(lambda_p,y_s(1:maxlambda+1),100, 's', "MarkerEdgeColor", "#005AFF", "MarkerFaceColor", "#005AFF");
scatter(lambda_p,y_prop_s(1:maxlambda+1),100, 'v', "MarkerEdgeColor", "#FF4B00", "MarkerFaceColor", "#FF4B00");
xlabel('Average arrival rate of Cloudlet 1 \lambda_{1} [jobs/s]'); ylabel('Offload fraction \phi_{21}');
l6 = legend('[Mondal+, OJ-COMS, 2020]', 'Proposed method');
l6.FontSize = 10;
hold off
figure(f7);
hold on
scatter(lambda_p,c1(1:maxlambda+1),100, 's', "MarkerEdgeColor", "#005AFF", "MarkerFaceColor", "#005AFF");
scatter(lambda_p,c1_prop(1:maxlambda+1),100, '+', 'red');
xlabel('Average arrival rate of Cloudlet 1 \lambda_{1} [jobs/s]'); ylabel('Actual arriving jobs at Cloudlet 1');
l7 = legend('[Mondal+, OJ-COMS, 2020]', 'Proposed method');
l7.FontSize = 10;
hold off
figure(f8);
hold on
scatter(lambda_p,c2(1:maxlambda+1),100, 'v', "MarkerEdgeColor", "#005AFF", "MarkerFaceColor", "#005AFF");
scatter(lambda_p,c2_prop(1:maxlambda+1),100, 'x', 'red');
xlabel('Average arrival rate of Cloudlet 1 \lambda_{1} [jobs/s]'); ylabel('Actual arriving jobs at Cloudlet 2');
l8 = legend('[Mondal+, OJ-COMS, 2020]', 'Proposed method');
l8.FontSize = 10;
hold off
figure(f9);
hold on
scatter(lambda_p/10000, c1c2dif(1:maxlambda+1),100, 'filled', 'blue');
scatter(lambda_p/10000, c1c2dif_prop(1:maxlambda+1),100, '^', "MarkerEdgeColor", "#FF0000", "MarkerFaceColor", "#FF0000");
xlabel('Initial utilization of cloudlet 1 \lambda_{1} [jobs/s]'); ylabel('Difference in utilization after offloading');
l9 = legend('Conventional method [9]', 'Proposed method');
l9.FontSize = 10;
hold off
figure(f10);
hold on
scatter(lambda_p, utility1_values(1:maxlambda+1),100, 'filled', 'blue');
scatter(lambda_p, utility1_p_values(1:maxlambda+1),100, '^', "MarkerEdgeColor", "#FF0000", "MarkerFaceColor", "#FF0000");
xlabel('Average arrival rate of Cloudlet 1 \lambda_{1} [jobs/s]'); ylabel('Utility of cloudlet 1');
l9 = legend('[Mondal+, OJ-COMS, 2020]', 'Proposed method');
l9.FontSize = 10;
hold off
figure(f11);
hold on
scatter(lambda_p, utility2_values(1:maxlambda+1),100, 'filled', 'blue');
scatter(lambda_p, utility2_p_values(1:maxlambda+1),100, '^', "MarkerEdgeColor", "#FF0000", "MarkerFaceColor", "#FF0000");
xlabel('Average arrival rate of Cloudlet 1 \lambda_{1} [jobs/s]'); ylabel('Utility of cloudlet2');
l9 = legend('[Mondal+, OJ-COMS, 2020]', 'Proposed method');
l9.FontSize = 10;
hold off
figure(f12);
hold on
scatter(lambda_p, utility1_comparison(1:maxlambda+1),100, 'filled', 'blue');
scatter(lambda_p, utility1_p_comparison(1:maxlambda+1),100, '^', "MarkerEdgeColor", "#FF0000", "MarkerFaceColor", "#FF0000");
xlabel('Average arrival rate of Cloudlet 1 \lambda_{1} [jobs/s]'); ylabel('Utility difference of cloudlet 1');
l9 = legend('Conventional method', 'Proposed method');
l9.FontSize = 10;
hold off
figure(f13);
hold on
scatter(lambda_p, utility2_comparison(1:maxlambda+1),100, 'filled', 'blue');
scatter(lambda_p, utility2_p_comparison(1:maxlambda+1),100, '^', "MarkerEdgeColor", "#FF0000", "MarkerFaceColor", "#FF0000");
xlabel('Average arrival rate of Cloudlet 1 \lambda_{1} [jobs/s]'); ylabel('Utility difference of cloudlet2');
l9 = legend('Conventional method', 'Proposed method');
l9.FontSize = 10;
hold off
figure(f14);
hold on
scatter(lambda_p, utility_sum(1:maxlambda+1),100, 'filled', 'blue');
scatter(lambda_p, utility_p_sum(1:maxlambda+1),100, '^', "MarkerEdgeColor", "#FF0000", "MarkerFaceColor", "#FF0000");
xlabel('Average arrival rate of Cloudlet 1 \lambda_{1} [jobs/s]'); ylabel('Total Utility');
l9 = legend('Conventional method', 'Proposed method');
l9.FontSize = 10;
hold off