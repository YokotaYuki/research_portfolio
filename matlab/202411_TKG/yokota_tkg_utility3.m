clc,clear;close all;

maxlambda = 11;
latency_nooffload_p = [];
latency_offload_p = [];
latency_mm1_1_p = [];
latency_mm1_2_p = [];

lambda = [400, 0];
lambda_p = [];
b = 1.0*10^-3;

lambda_offload1 = zeros(11);
lambda_nooffload1 = zeros(11);
lambda_offload2 = zeros(11);
lambda_nooffload2 = zeros(11);
lambda1 = zeros(11);
lambda2 = zeros(11);
rho1 = zeros(11);
rho2 = zeros(11);
E_nopriority1 = zeros(11);
E_nopriority2 = zeros(11);
E_offload1 = zeros(11);
E_offload2 = zeros(11);

utility1 = zeros(11);
utility2 = zeros(11);

E_mm1_1 = 0;
E_mm1_2 = 0;


varphi = linspace(0, 1, 11);
varphi_1 = 0;
varphi_2 = 0;

d_rq = 5*10^-3;

for n=1:maxlambda
    lambda(2) = 100*(n-1);
    lambda_p = [lambda_p lambda(2)];
    
    E_mm1_1 = 1/(1/b - lambda(1));
    E_mm1_2 = 1/(1/b - lambda(2));
    latency_nooffload_p = [latency_nooffload_p E_nopriority1];
    latency_mm1_1_p = [latency_mm1_1_p E_mm1_1];
    latency_mm1_2_p = [latency_mm1_2_p E_mm1_2];
end

lambda = [500, 900];
for i=1:11
    varphi_1 = varphi(i);
    for j=1:11
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
        E_nopriority1(i,j) =  1/(1/b - lambda1(i,j));
        E_nopriority2(i,j) =  1/(1/b - lambda2(i,j));
        E_offload1(i,j) = 0;
        E_offload2(i,j) = 0;
        utility1(i,j) = (-lambda_offload1(i,j) + lambda_offload2(i,j) + 0.1*lambda1(i,j))*b + 100*(lambda1(i,j)*b*(d_rq-E_nopriority1(i,j)));
        utility2(i,j) = (-lambda_offload2(i,j) + lambda_offload1(i,j) + 0.1*lambda2(i,j))*b + 100*(lambda2(i,j)*b*(d_rq-E_nopriority2(i,j)));
    end
end

utility1_latency_offload = d_rq - E_offload1;
utility1_latency_nooffload = d_rq - E_nopriority1;
utility2_latency_offload = d_rq - E_offload2;
utility2_latency_nooffload = d_rq - E_nopriority2;




f1 = figure;
f2 = figure;
f3 = figure;
f4 = figure;
f5 = figure;

 figure(f1);
 hold on
 grid on
 scatter(lambda_p, latency_mm1_1_p(1:maxlambda),75, 'filled', 'blue');
 scatter(lambda_p, latency_mm1_2_p(1:maxlambda), 75, 'filled', 'red'); 
 ylim([0 inf]);
 xlabel('Number of prioritized jobs'); ylabel('Average latency');
 l1 = legend('cloudlet1', 'cloudlet2');
 l1.FontSize = 10;
 hold off

figure(f2);
hold on
grid on
for i=1:11
    scatter3(varphi(i), varphi, utility1_latency_offload(i,:), 100, 'filled','blue');
    scatter3(varphi(i), varphi, utility1_latency_nooffload(i,:), 100, 'filled', 'red');
end
zlim([-3*10^-3 inf]);
hold off
figure(f3);
hold on
grid on
for i=1:11
    scatter3(varphi(i), varphi, utility2_latency_offload(i,:), 100, 'filled','blue');
    scatter3(varphi(i), varphi, utility2_latency_nooffload(i,:), 100, 'filled', 'red');
end
zlim([-3*10^-3 inf]);
hold off
figure(f4);
hold on
grid on
for i=1:11
    scatter3(varphi(i), varphi, lambda1(i,:), 100, 'filled','blue');
    scatter3(varphi(i), varphi, lambda2(i,:), 100, 'filled', 'red');
end
hold off
figure(f5);
hold on
grid on
for i=1:11
    scatter3(varphi(i), varphi, utility1(i,:), 100, 'filled','blue');
    scatter3(varphi(i), varphi, utility2(i,:), 100, 'filled', 'red');
end
hold off