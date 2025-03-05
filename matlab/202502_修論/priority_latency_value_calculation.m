clc;
clear;
close all;

% 定数の定義
b = 1.0 * 10^-3; % サービス時間
x_range = linspace(0, 1, 500); % λb の範囲（0 から 0.9 未満）
varphi = 0.3;

% レイテンシ関数の定義
E_offload = @(x) (x .* b) ./ (1 - varphi .* x) + b;
E_nooffload = @(x) (x .* b) ./ ((1 - varphi .* x) .* (1 - varphi .* x - (1-varphi) .* x)) + b;
E_mm1 = @(x) 1 ./ (1 / b - x / b);

% レイテンシ計算
latency_offload_p = E_offload(x_range);
latency_nooffload_p = E_nooffload(x_range);
latency_mm1 = E_mm1(x_range);

% グラフの描画
f1 = figure;
hold on;
%grid on;

% 線でプロット
plot(x_range, latency_offload_p, '-r', 'LineWidth', 2, 'DisplayName', 'T_{21,1}');
plot(x_range, latency_nooffload_p, '-b', 'LineWidth', 2, 'DisplayName', 'T_{11,2}');
plot(x_range, latency_mm1, '-g', 'LineWidth', 2, 'DisplayName', 'T_{11}');
yline(0.005, '--k', 'D_{rq} = 0.005', 'LineWidth', 1, 'LabelHorizontalAlignment', 'left', 'FontSize', 10, 'DisplayName', 'D_{rq}');
% 軸ラベルと凡例
xlabel('Utilization \rho_1', 'FontSize', 12);
ylabel('Average Latency', 'FontSize', 12);
ylim([0 0.03]);
legend('Location', 'best', 'FontSize', 10);

% 終了処理
hold off;
