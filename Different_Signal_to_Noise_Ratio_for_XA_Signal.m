% 1. 生成原始信号 Run the single_sensor_collect_data.m
t_sc = 1000000; % 时间单位转换
t = kgrid.t_array * t_sc;  % 时间轴
s = sensor_data.p(1, :);  % 原始信号（正弦波）
SNR_dB = [10, 5, 0];  % 信噪比 (dB)
noisy_signals = cell(length(SNR_dB), 1);  % 存储加噪信号

% 2. 定义不同SNR并添加高斯白噪声
for i = 1:length(SNR_dB)
    % 计算噪声标准差
    signal_power = var(s);  % 信号功率（方差）
    noise_power = signal_power / (10^(SNR_dB(i)/10));  % 噪声功率
    sigma = sqrt(noise_power);  % 噪声标准差
    
    % 生成高斯白噪声
    noise = sigma * randn(size(s));  % 均值为0，标准差为sigma
    
    % 叠加噪声
    noisy_signals{i} = s + noise;
end

% 3. 绘制原始信号和不同SNR的加噪信号
figure;
subplot(length(SNR_dB)+1, 1, 1);
plot(t, s, 'Color', [91 171 84]/255, 'LineWidth', 1.5);
title('Initial XA signal without noise');
xlabel('Time (μs)');
ylabel('Pressure (Pa)');
grid on;

for i = 1:length(SNR_dB)
    subplot(length(SNR_dB)+1, 1, i+1);
    plot(t, noisy_signals{i}, 'Color', [136 95 156]/255, 'LineWidth', 1);
    title(['Noisy signal (SNR = ', num2str(SNR_dB(i)), ' dB)']);
    xlabel('Time (μs)');
    ylabel('Pressure (Pa)');
    grid on;
end

% 4.计算实际SNR（验证）
for i = 1:length(SNR_dB)
    signal_power_actual = var(s);  % 信号功率
    noise_power_actual = var(noisy_signals{i} - s);  % 噪声功率
    SNR_actual = 10*log10(signal_power_actual / noise_power_actual);  % 实际SNR

    disp(['设定 SNR = ', num2str(SNR_dB(i)), ' dB | 实际 SNR = ', num2str(SNR_actual), ' dB']);
end
