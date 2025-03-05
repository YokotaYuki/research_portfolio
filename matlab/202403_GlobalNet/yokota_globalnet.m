clc,clear;close all;

n = 2;
m = 2;
maxlambda = 10;
latency_priority11 = [];
latency_priority12 = [];
latency_priority21 = [];
latency_priority22 = [];
latency_priority1 = [];
latency_priority2 = [];
latency_nooffloading = [];
latency_nopriority = [];
varphi1_p = [];
varphi2_p = [];
varphi11_p = [];
varphi22_p = [];

lambda = [6000, 0;
          3000, 2000];
lambda_offload = zeros(n);
lambda_i = [sum(lambda(1,:)), sum(lambda(2,:))];
lambda_p = [];
b = [0.6*10^-4, 1.2*10^-4];
b_i(1) = (lambda(1,1)*b(1)+lambda(1,2)*b(2))/lambda_i(1);
b_i(2) = (lambda(2,1)*b(1)+lambda(2,2)*b(2))/lambda_i(2);
rho = zeros(n);

varphi = zeros(n);
varphi_x = linspace(0, 1, 11);
varphi_y = linspace(0, 1, 11);

E_wp1 = zeros(n); %オフロードしたジョブの優先度が低い場合の遅延
E_wp2 = zeros(n); %オフロードしたジョブの優先度が高い場合の遅延
E_wnp = zeros(n); %優先度をつけない場合の遅延
latency1_11 = zeros(11); latency2_11 = zeros(11);
latency1_12 = zeros(11); latency2_12 = zeros(11);
latency1_21 = zeros(11); latency2_21 = zeros(11);
latency1_22 = zeros(11); latency2_22 = zeros(11);
latency1_sum = zeros(11); latency2_sum = zeros(11);
latency_nopriority_sum = zeros(11);
latency_c1 = zeros(11);
latency_c2 = zeros(11);
rho1 = zeros(11);
rho2 = zeros(11);

for i=1:maxlambda
    %lambda(1,1) = 2000;
    %lambda(2,1) = 2000;
    lambda(1,2) = 1000*(i-1);
    %lambda(2,2) = 100*(i-1);
    lambda_i = [sum(lambda(1,:)), sum(lambda(2,:))];

    lambda_p = [lambda_p lambda(1,2)];
    %全てのオフロード組み合わせごとに平均遅延時間を求める
    for l=1:11
        varphi(1,1) = 1-varphi_x(l);
        varphi(1,2) = varphi_x(l);
        for k=1:11
            varphi(2,1) = varphi_y(k);
            varphi(2,2) = 1-varphi_y(k);
    
            % オフロードしないジョブの到着率
            % オフロードするジョブの到着率
            
            % オフロードしないジョブの到着率
            % オフロードするジョブの到着率
            lambda_offload(1,1) = lambda_i(1)*varphi(1,1);
            lambda_offload(1,2) = lambda_i(1)*varphi(1,2);
            lambda_offload(2,1) = lambda_i(2)*varphi(2,1);
            lambda_offload(2,2) = lambda_i(2)*varphi(2,2);
    
            % オフロードしないジョブの利用率
            % オフロードするジョブの利用率
            rho(1,1) = lambda_offload(1,1)*b_i(1);
            rho(1,2) = lambda_offload(1,2)*b_i(1);
            rho(2,1) = lambda_offload(2,1)*b_i(2);
            rho(2,2) = lambda_offload(2,2)*b_i(2);
    
            E_wp1(1,1) = 0.5*(lambda_offload(1,1)*b_i(1)^2 + lambda_offload(2,1)*b_i(2)^2)/(1-rho(1,1));
            E_wp1(1,2) = 0.5*(lambda_offload(1,2)*b_i(1)^2 + lambda_offload(2,2)*b_i(2)^2)/((1-rho(1,2))*(1-rho(1,2)-rho(2,2)));
            E_wp1(2,1) = 0.5*(lambda_offload(1,1)*b_i(1)^2 + lambda_offload(2,1)*b_i(2)^2)/((1-rho(2,1))*(1-rho(2,1)-rho(1,1)));
            E_wp1(2,2) = 0.5*(lambda_offload(1,2)*b_i(1)^2 + lambda_offload(2,2)*b_i(2)^2)/(1-rho(2,2));
            E_wp2(1,1) = 0.5*(lambda_offload(1,1)*b_i(1)^2 + lambda_offload(2,1)*b_i(2)^2)/((1-rho(1,1))*(1-rho(1,1)-rho(2,1)));
            E_wp2(1,2) = 0.5*(lambda_offload(1,2)*b_i(1)^2 + lambda_offload(2,2)*b_i(2)^2)/(1-rho(1,2));
            E_wp2(2,1) = 0.5*(lambda_offload(1,1)*b_i(1)^2 + lambda_offload(2,1)*b_i(2)^2)/(1-rho(2,1));
            E_wp2(2,2) = 0.5*(lambda_offload(1,2)*b_i(1)^2 + lambda_offload(2,2)*b_i(2)^2)/((1-rho(2,2))*(1-rho(1,2)-rho(2,2)));        
            E_wnp(1,1) = 0.5*(lambda_offload(1,1)+lambda_offload(2,1))*((b_i(1)*lambda_offload(1,1)+b_i(2)*lambda_offload(2,1))/(lambda_offload(1,1)+lambda_offload(2,1)))^2/(1-(lambda_offload(1,1)+lambda_offload(2,1))*(b_i(1)*lambda_offload(1,1)+b_i(2)*lambda_offload(2,1))/(lambda_offload(1,1)+lambda_offload(2,1)));
            E_wnp(2,2) = 0.5*(lambda_offload(2,2)+lambda_offload(1,2))*((b_i(2)*lambda_offload(2,2)+b_i(1)*lambda_offload(1,2))/(lambda_offload(2,2)+lambda_offload(1,2)))^2/(1-(lambda_offload(2,2)+lambda_offload(1,2))*(b_i(2)*lambda_offload(2,2)+b_i(1)*lambda_offload(1,2))/(lambda_offload(2,2)+lambda_offload(1,2)));
            
            latency1_11(l,k) = E_wp1(1,1); latency2_11(l,k) = E_wp2(1,1);
            latency1_12(l,k) = E_wp1(1,2); latency2_12(l,k) = E_wp2(1,2);
            latency1_21(l,k) = E_wp1(2,1); latency2_21(l,k) = E_wp2(2,1);
            latency1_22(l,k) = E_wp1(2,2); latency2_22(l,k) = E_wp2(2,2);
            latency1_sum(l,k) = varphi(1,1)*E_wp1(1,1)+varphi(1,2)*E_wp1(1,2)+varphi(2,1)*E_wp1(2,1)+varphi(2,2)*E_wp1(2,2);
            latency2_sum(l,k) = varphi(1,1)*E_wp2(1,1)+varphi(1,2)*E_wp2(1,2)+varphi(2,1)*E_wp2(2,1)+varphi(2,2)*E_wp2(2,2);
            latency_c1(l,k) = E_wnp(1,1);
            latency_c2(l,k) = E_wnp(2,2);
            latency_nopriority_sum(l,k) = varphi(1,1)*E_wnp(1,1)+varphi(2,1)*E_wnp(1,1)+varphi(2,2)*E_wnp(2,2)+varphi(1,2)*E_wnp(2,2);
            if latency1_11(l,k) <= 0
                latency1_11(l,k) = 1;
            end
            if latency1_12(l,k) <= 0
                latency1_12(l,k) = 1;
            end
            if latency1_21(l,k) <= 0
                latency1_21(l,k) = 1;
            end
            if latency1_22(l,k) <= 0
                latency1_22(l,k) = 1;
            end
            if latency2_11(l,k) <= 0
                latency2_11(l,k) = 1;
            end
            if latency2_12(l,k) <= 0
                latency2_12(l,k) = 1;
            end
            if latency2_21(l,k) <= 0
                latency2_21(l,k) = 1;
            end
            if latency2_22(l,k) <= 0
                latency2_22(l,k) = 1;
            end
            if latency_c1(l,k) <= 0
                latency_c1(l,k) = 1;
            end
            if latency_c2(l,k) <= 0
                latency_c2(l,k) = 1;
            end
            if latency1_sum(l,k) <= 0
                latency1_sum(l,k) = 1;
            end
            if latency2_sum(l,k) <= 0
                latency2_sum(l,k) = 1;
            end
            if latency_nopriority_sum(l,k) <= 0
                latency_nopriority_sum(l,k) = 1;
            end
            rho1(l,k) = rho(1);
            rho2(l,k) = rho(2);
        end
    end
    [m1, i1] = min(latency1_11+latency1_22);
    [m2, i2] = min(m1);
    [m11, i11] = min(latency_nopriority_sum);
    [m22, i22] = min(m11);
    varphi1_p = [varphi1_p 0.1*i1-0.1];
    varphi2_p = [varphi2_p 0.1*i2-0.1];
    varphi11_p = [varphi11_p 0.1*i11-0.1];
    varphi22_p = [varphi22_p 0.1*i22-0.1];
    latency_priority11 = [latency_priority11 m2];
    latency_priority12 = [latency_priority12 min(min(latency1_12+latency1_21))];
    latency_priority21 = [latency_priority21 min(min(latency2_11+latency2_22))];
    latency_priority22 = [latency_priority22 min(min(latency2_12+latency2_21))];
    latency_nooffloading = [latency_nooffloading latency_c1(1,1)+latency_c2(1,1)];
    latency_nopriority = [latency_nopriority m22];
    latency_priority1 = [latency_priority1 min(min(latency1_sum))];
    latency_priority2 = [latency_priority2 min(min(latency2_sum))];
end

f1 = figure;f2 = figure; f3=figure;

 figure(f1);
 hold on
 grid on
 scatter(lambda_p, latency_nooffloading(1:maxlambda),75, 'filled', 'blue');
 scatter(lambda_p, latency_priority1(1:maxlambda),75, 'filled', 'red');
 scatter(lambda_p, latency_nopriority(1:maxlambda),75, 'filled', 'yellow');
 scatter(lambda_p, latency_priority11(1:maxlambda),75, 'x', 'green');
 scatter(lambda_p, latency_priority12(1:maxlambda),75, '+', 'magenta');
 ylim([0 inf]);
 xlabel('Number of prioritized jobs'); ylabel('Average latency');
 l1 = legend('Without offloading', 'With priority','Without priority', 'High priority', 'Low priority');
 l1.FontSize = 10;
 hold off

figure(f2);
hold on
scatter(lambda_p, latency_nooffloading(1:maxlambda),100, 'filled', 'blue');
scatter(lambda_p, latency_priority2(1:maxlambda),100, 'filled', 'red');
scatter(lambda_p, latency_priority22(1:maxlambda),100, '+', 'magenta');
scatter(lambda_p, latency_priority21(1:maxlambda),100, 'x', 'green');

xlabel('Arrival rate of \lambda_{12}'); ylabel('Average latency');
l1 = legend('Total latency of unprioritized jobs', 'Total latency of prioritized jobs', 'High priority', 'Low priority');
l1.FontSize = 10;
hold off

figure(f3);
hold on
scatter(lambda_p, varphi1_p(1:maxlambda),70, 'filled', 'blue');
scatter(lambda_p, varphi2_p(1:maxlambda),70, 'x', 'red');
scatter(lambda_p, varphi11_p(1:maxlambda),70, 'filled', 'green');
scatter(lambda_p, varphi22_p(1:maxlambda),70, '+', 'magenta');
ylim([0 1]);
xlabel('Arrival rate of \lambda_{12}'); ylabel('Offload fraction');
l1 = legend('\phi_{12}^{priority*}', '\phi_{21}^{priority*}', '\phi_{12}^{offload*}', '\phi_{21}^{offload*}');
l1.FontSize = 10;
hold offx