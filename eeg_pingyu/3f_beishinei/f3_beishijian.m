%%
clear;clc;
%% Specify Basic information of different groups
% 1. 定义三组数据所在的路径 (请修改为你实际的Group 3路径)
group1_dir = 'path_to_your_data';
group2_dir = 'path_to_your_data';
group3_dir = 'path_to_your_data'; % %% 修改：添加第三组路径

% 2. 筛选三组数据
group1_files = dir([group1_dir, filesep, '*.set']);
group2_files = dir([group2_dir, filesep, '*.set']);
group3_files = dir([group3_dir, filesep, '*.set']); % %% 修改：添加第三组文件列表

% 定义数据的采样率信息
Fs = 250;
% 数据的长度
L = 500;
NFFT = 2^nextpow2(L);
% 制作频率列表
f = linspace(0,Fs/2,NFFT/2+1);

%% Load and perform FFT transform on data of different groups
% ================= Group 1 =================
for i=1:length(group1_files)
    subj_fn = group1_files(i).name;
    EEG = pop_loadset([group1_dir, filesep, subj_fn]);
    for ii=1:size(EEG.data,1)
        for jj=1:size(EEG.data,3)
            y = squeeze(EEG.data(ii,:,jj));
            temp = fft(y,NFFT);
            Y(jj,:) = abs(temp(1:NFFT/2+1))*2/L;    
        end
        group1_FFT_power(i,ii,:) = squeeze(mean(Y,1)); clear Y;
    end
end

% ================= Group 2 =================
for i=1:length(group2_files)
    subj_fn = group2_files(i).name;
    EEG = pop_loadset([group2_dir, filesep, subj_fn]);
    for ii=1:size(EEG.data,1)
        for jj=1:size(EEG.data,3)
            y = squeeze(EEG.data(ii,:,jj));
            temp = fft(y,NFFT);
            Y(jj,:) = abs(temp(1:NFFT/2+1))*2/L;
        end
        group2_FFT_power(i,ii,:) = squeeze(mean(Y,1)); clear Y;
    end
end

% ================= Group 3 (%% 修改：新增读取循环) =================
for i=1:length(group3_files)
    subj_fn = group3_files(i).name;
    EEG = pop_loadset([group3_dir, filesep, subj_fn]);
    for ii=1:size(EEG.data,1)
        for jj=1:size(EEG.data,3)
            y = squeeze(EEG.data(ii,:,jj));
            temp = fft(y,NFFT);
            Y(jj,:) = abs(temp(1:NFFT/2+1))*2/L;
        end
        group3_FFT_power(i,ii,:) = squeeze(mean(Y,1)); clear Y;
    end
end

%% plot power and scalp maps (画频谱图和地形图)
% 假设我们感兴趣的是30Hz一下的部分
f_idx = find(f <= 30);
f_plot = f(f_idx);

% --- 画频谱图 (Spectrum) ---
figure;
Cz = 28; % 定义电极
% 绘制三组的线
plot(f_plot, squeeze(mean(group1_FFT_power(:,Cz, f_idx),1)),'r','linewidth',1.5); hold on;
plot(f_plot, squeeze(mean(group2_FFT_power(:,Cz, f_idx),1)),'g','linewidth',1.5);
plot(f_plot, squeeze(mean(group3_FFT_power(:,Cz, f_idx),1)),'b','linewidth',1.5); % %% 修改：添加第三组，蓝色
legend('Group 1', 'Group 2', 'Group 3'); % 添加图例
title(['Group level FFT at Channel ' num2str(Cz)], 'fontsize', 16);
xlabel('Frequency (Hz)');
ylabel('Amplitude (\muV)');

% --- 画地形图 (Topoplot) ---
% 定义Alpha频段
alpha_idx = find((f >=8) & (f<=12));

% 计算各组Alpha频段均值
group1_alpha_mag_avg = squeeze(mean(mean(group1_FFT_power(:,:,alpha_idx),3),1));
group2_alpha_mag_avg = squeeze(mean(mean(group2_FFT_power(:,:,alpha_idx),3),1));
group3_alpha_mag_avg = squeeze(mean(mean(group3_FFT_power(:,:,alpha_idx),3),1)); % %% 修改

% 确定Colorbar范围，保证三张图颜色标准一致
c_lim = [0 max([max(group1_alpha_mag_avg), max(group2_alpha_mag_avg), max(group3_alpha_mag_avg)])];

figure;
subplot(131); topoplot(group1_alpha_mag_avg, EEG.chanlocs,'maplimits',c_lim); title('Group 1', 'fontsize', 14);
subplot(132); topoplot(group2_alpha_mag_avg, EEG.chanlocs,'maplimits',c_lim); title('Group 2', 'fontsize', 14);
subplot(133); topoplot(group3_alpha_mag_avg, EEG.chanlocs,'maplimits',c_lim); title('Group 3', 'fontsize', 14); % %% 修改
colorbar; % 在最后一个图旁显示颜色条

%% Specify (ROI) channels and compare across frequency points (ANOVA)
% 选定电极，对比各频率点上的组间差异
Pz = 46;

% 提取数据: Sub * Freq
g1_data = squeeze(group1_FFT_power(:,Pz,:));
g2_data = squeeze(group2_FFT_power(:,Pz,:));
g3_data = squeeze(group3_FFT_power(:,Pz,:)); % %% 修改

f_idx = find(f<=30);
f_band = f(f_idx);

% 预分配数组
pvals = zeros(1, length(f_band));
fvals = zeros(1, length(f_band));

% 获取各组人数 (防止三组人数不一致)
n1 = size(g1_data, 1);
n2 = size(g2_data, 1);
n3 = size(g3_data, 1);

% 制作分组标签 (用于anova1函数)
group_labels = [ones(n1,1); 2*ones(n2,1); 3*ones(n3,1)];

for i=1:length(f_band)
    % 提取当前频率点的数据
    d1 = g1_data(:,i);
    d2 = g2_data(:,i);
    d3 = g3_data(:,i);
    
    % 拼合数据: (n1+n2+n3) * 1
    data_all = [d1; d2; d3];
    
    % %% 修改：使用单因素方差分析 anova1
    % 'off' 表示不显示ANOVA自带的箱线图
    [p, tbl, stats] = anova1(data_all, group_labels, 'off');
    
    pvals(i) = p;
    fvals(i) = tbl{2,5}; % 提取F值 (F-statistic)
end

figure;
subplot(2,1,1);
plot(f_band, pvals, 'k', 'LineWidth', 1.5);
hold on; 
yline(0.05, 'r--'); % 画出0.05显著性水平线
ylim([0 0.1]); % 重点看0到0.1的范围
xlim([0 30]);
title(['ANOVA P Values at Channel ' num2str(Pz)]);
xlabel('Frequency (Hz)');
ylabel('p-value');

subplot(2,1,2);
plot(f_band, fvals, 'b', 'LineWidth', 1.5);
xlim([0 30]);
title('F Statistics');
xlabel('Frequency (Hz)');
ylabel('F-value');

%% Specify (ROI) frequency band and compare across channels (ANOVA Topoplot)
clear pvals fvals;

% 1. 选定频段 (Alpha)
f_ROI = find((f>=8) & (f<=12));

% 2. 提取数据: 沿着频率做平均 -> Sub * Ch
g1_ch = squeeze(mean(group1_FFT_power(:,:,f_ROI),3));
g2_ch = squeeze(mean(group2_FFT_power(:,:,f_ROI),3));
g3_ch = squeeze(mean(group3_FFT_power(:,:,f_ROI),3)); % %% 修改

% 获取人数
n1 = size(g1_ch, 1);
n2 = size(g2_ch, 1);
n3 = size(g3_ch, 1);
group_labels = [ones(n1,1); 2*ones(n2,1); 3*ones(n3,1)];

pvals = zeros(1, EEG.nbchan);
fvals = zeros(1, EEG.nbchan);

% 对于每个通道进行ANOVA
for i=1:EEG.nbchan
    d1 = g1_ch(:,i);
    d2 = g2_ch(:,i);
    d3 = g3_ch(:,i);
    
    data_all = [d1; d2; d3];
    
    % %% 修改：ANOVA
    [p, tbl, stats] = anova1(data_all, group_labels, 'off');
    
    pvals(i) = p;
    fvals(i) = tbl{2,5}; % F值
end

figure;
% 画F值的地形图 (反映组间差异的大小)
subplot(121);
topoplot(fvals, EEG.chanlocs); 
title('Alpha Band: ANOVA F-values','fontsize',12);
colorbar;

% 画P值的地形图 (通常比较难看，因为不显著的地方P值很大，建议反向画 1-p 或者只画显著点)
% 这里演示一种方法：画出显著(p<0.05)的通道位置
subplot(122);
topoplot(fvals, EEG.chanlocs, 'pmask', pvals < 0.05); 
title('Significant Channels (p<0.05)','fontsize',12);

%% 提取SPSS数据 (单因素三水平)
% 提取所有被试在某个电极、特定频段的均值，用于导出到Excel/SPSS做后续分析

f_ROI = find((f>=8) & (f<=12)); % Alpha频段

% 提取 Sub * Ch
g1_avg = squeeze(mean(group1_FFT_power(:,:,f_ROI),3));
g2_avg = squeeze(mean(group2_FFT_power(:,:,f_ROI),3));
g3_avg = squeeze(mean(group3_FFT_power(:,:,f_ROI),3));

target_ch = 28; % 选择想提取的电极 (例如 Cz)

% 提取该电极的数据 (列向量)
data_spss_1 = g1_avg(:, target_ch);
data_spss_2 = g2_avg(:, target_ch);
data_spss_3 = g3_avg(:, target_ch);

% 打印提示
fprintf('已提取第 %d 通道的数据。\n', target_ch);
fprintf('Group 1 样本数: %d\n', length(data_spss_1));
fprintf('Group 2 样本数: %d\n', length(data_spss_2));
fprintf('Group 3 样本数: %d\n', length(data_spss_3));

% 如果想保存为 CSV 供 SPSS 使用：
% 由于各组人数可能不同，不能直接拼成矩阵，通常保存为 "长格式" (Long Format)
% 第一列：数据，第二列：组别标签
all_data = [data_spss_1; data_spss_2; data_spss_3];
all_labels = [ones(length(data_spss_1),1); 2*ones(length(data_spss_2),1); 3*ones(length(data_spss_3),1)];

% 生成一个 table 并保存
T = table(all_data, all_labels, 'VariableNames', {'Power', 'Group'});
writetable(T, 'EEG_Alpha_Data_For_SPSS.csv');
disp('数据已保存为 EEG_Alpha_Data_For_SPSS.csv');