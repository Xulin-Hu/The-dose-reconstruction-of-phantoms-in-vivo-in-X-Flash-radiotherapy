
% ʱ���źŴ洢�ڱ��� signal ��
signal = sensor_data.p(1, :); % �������������ʱ���ź�
sample_rate = 1 / kgrid.dt; % ���ò���Ƶ��
    
% ���и���Ҷ�任
freq_signal = fft(signal);

% ��ȡƵ����
n = length(signal);
f = (0:n-1) * (sample_rate / n); % Ƶ������

% ���ȹ�һ��
Y_magnitude = abs(freq_signal/n); 

% ����Ƶ��
figure;
plot(f, Y_magnitude, 'LineWidth', 1,'Color', 'b');
xlabel('Frequency (Hz)','FontSize',10), ylabel('Magnitude','FontSize',10);
title('Frequency Spectrum');
xlim([0 sample_rate/2]); % ֻ��ʾ��Ƶ�ʲ���
grid on;