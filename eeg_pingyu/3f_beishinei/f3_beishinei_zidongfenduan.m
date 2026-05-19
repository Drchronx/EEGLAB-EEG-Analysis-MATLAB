%%
clear;clc;
%% 1. 基础设置
group1_dir = 'path_to_your_data'; % 条件1
group2_dir = 'path_to_your_data'; % 条件2
group3_dir = 'path_to_your_data'; % 条件3

group1_files = dir([group1_dir, filesep, '*.set']);
group2_files = dir([group2_dir, filesep, '*.set']);
group3_files = dir([group3_dir, filesep, '*.set']);

% === 安全检查：确保三个文件夹文件数量一致 ===
if length(group1_files) ~= length(group2_files) || length(group1_files) ~= length(group3_files)
    error('错误：组内设计要求三个文件夹的文件数量必须严格一致！请检查文件。');
end

% 定义数据的采样率信息
Fs = 500; 
% 数据的长度 (2秒 * 500Hz = 1000点)
L = 1000; 

NFFT = 2^nextpow2(L);
f = linspace(0,Fs/2,NFFT/2+1);

%% 2. 读取数据并分段 (Load Data & Epoching)
n_subj = length(group1_files); % 被试数量

% ==================== 读取条件 1 ====================
for i=1:n_subj
    subj_fn = group1_files(i).name;
    EEG = pop_loadset([group1_dir, filesep, subj_fn]);
    
    % %% 修改：新增2秒分段代码 %%
    % 如果数据是连续的(2维)，则切成2秒一段
    if length(size(EEG.data)) == 2
        % recurrence: 重复间隔2秒; limits: [0 2]表示取0到2秒的数据; rmbase: NaN表示不去基线(频域分析通常不需要)
        EEG = eeg_regepochs(EEG, 'recurrence', 2, 'limits', [0 2], 'rmbase', NaN);
        EEG = eeg_checkset(EEG); % 检查数据结构完整性
    end
    
    for ii=1:size(EEG.data,1)
        for jj=1:size(EEG.data,3)
            y = squeeze(EEG.data(ii,:,jj));
            temp = fft(y,NFFT);
            Y(jj,:) = abs(temp(1:NFFT/2+1))*2/L;    
        end
        group1_FFT_power(i,ii,:) = squeeze(mean(Y,1)); clear Y;
    end
end

% ==================== 读取条件 2 ====================
for i=1:n_subj
    subj_fn = group2_files(i).name; 
    EEG = pop_loadset([group2_dir, filesep, subj_fn]);
    
    % %% 修改：新增2秒分段代码 %%
    if length(size(EEG.data)) == 2
        EEG = eeg_regepochs(EEG, 'recurrence', 2, 'limits', [0 2], 'rmbase', NaN);
        EEG = eeg_checkset(EEG);
    end
    
    for ii=1:size(EEG.data,1)
        for jj=1:size(EEG.data,3)
            y = squeeze(EEG.data(ii,:,jj));
            temp = fft(y,NFFT);
            Y(jj,:) = abs(temp(1:NFFT/2+1))*2/L;
        end
        group2_FFT_power(i,ii,:) = squeeze(mean(Y,1)); clear Y;
    end
end

% ==================== 读取条件 3 ====================
for i=1:n_subj
    subj_fn = group3_files(i).name;
    EEG = pop_loadset([group3_dir, filesep, subj_fn]);
    
    % %% 修改：新增2秒分段代码 %%
    if length(size(EEG.data)) == 2
        EEG = eeg_regepochs(EEG, 'recurrence', 2, 'limits', [0 2], 'rmbase', NaN);
        EEG = eeg_checkset(EEG);
    end
    
    for ii=1:size(EEG.data,1)
        for jj=1:size(EEG.data,3)
            y = squeeze(EEG.data(ii,:,jj));
            temp = fft(y,NFFT);
            Y(jj,:) = abs(temp(1:NFFT/2+1))*2/L;
        end
        group3_FFT_power(i,ii,:) = squeeze(mean(Y,1)); clear Y;
    end
end

%% 3. 画频谱图 (Spectrum)
f_idx = find(f <= 30);
f_plot = f(f_idx);
Cz = 12; % 定义电极

figure;
% 使用 mean 计算每种条件的平均值
plot(f_plot, squeeze(mean(group1_FFT_power(:,Cz, f_idx),1)),'r','linewidth',1.5); hold on;
plot(f_plot, squeeze(mean(group2_FFT_power(:,Cz, f_idx),1)),'g','linewidth',1.5);
plot(f_plot, squeeze(mean(group3_FFT_power(:,Cz, f_idx),1)),'b','linewidth',1.5);
legend('Condition 1', 'Condition 2', 'Condition 3');
title(['Grand Average FFT at Channel ' num2str(Cz)], 'fontsize', 16);
xlabel('Frequency (Hz)'); ylabel('Amplitude (\muV)');

%% 4. 画地形图 (Topoplot)
alpha_idx = find((f >=4) & (f<=8));

group1_alpha_mag_avg = squeeze(mean(mean(group1_FFT_power(:,:,alpha_idx),3),1));
group2_alpha_mag_avg = squeeze(mean(mean(group2_FFT_power(:,:,alpha_idx),3),1));
group3_alpha_mag_avg = squeeze(mean(mean(group3_FFT_power(:,:,alpha_idx),3),1));

c_lim = [0 max([max(group1_alpha_mag_avg), max(group2_alpha_mag_avg), max(group3_alpha_mag_avg)])];

figure;
subplot(131); topoplot(group1_alpha_mag_avg, EEG.chanlocs,'maplimits',c_lim); title('Cond 1', 'fontsize', 14);
subplot(132); topoplot(group2_alpha_mag_avg, EEG.chanlocs,'maplimits',c_lim); title('Cond 2', 'fontsize', 14);
subplot(133); topoplot(group3_alpha_mag_avg, EEG.chanlocs,'maplimits',c_lim); title('Cond 3', 'fontsize', 14);
colorbar;

%% 5. 单因素重复测量方差分析 (Repeated Measures ANOVA)
% === 修改开始：自动查找电极索引，不再硬编码数字 ===
target_elec_name = 'Cz'; % 在这里输入你想要的电极名字，比如 'Pz', 'Cz', 'Fz'

% 在 chanlocs 中查找该名字对应的索引
Pz = find(strcmpi({EEG.chanlocs.labels}, target_elec_name));

% 安全检查：如果没找到该电极
if isempty(Pz)
    disp('目前的电极列表如下：');
    disp({EEG.chanlocs.labels});
    error(['错误：在当前数据中找不到电极 [' target_elec_name ']！请检查上方列表，确认名字拼写是否正确（如 PZ vs Pz）。']);
else
    fprintf('已找到电极 %s，其当前索引为第 %d 通道。\n', target_elec_name, Pz);
end
% === 修改结束 ===

% 提取数据: Sub * Freq
g1_data = squeeze(group1_FFT_power(:,Pz,:));
g2_data = squeeze(group2_FFT_power(:,Pz,:));
g3_data = squeeze(group3_FFT_power(:,Pz,:));

f_idx = find(f<=30);
f_band = f(f_idx);

pvals = zeros(1, length(f_band));
fvals = zeros(1, length(f_band));

Conditions = [1 2 3]';
varNames = {'C1', 'C2', 'C3'};

fprintf('正在计算频率点的重复测量方差分析，请稍候...\n');

for i=1:length(f_band)
    % 准备当前频率点的数据表 (宽格式)
    d1 = g1_data(:,i);
    d2 = g2_data(:,i);
    d3 = g3_data(:,i);
    
    t = table(d1, d2, d3, 'VariableNames', varNames);
    
    % 创建重复测量模型
    rm = fitrm(t, 'C1-C3 ~ 1', 'WithinDesign', table(Conditions));
    
    % 进行方差分析
    ranovatbl = ranova(rm);
    
    % 提取结果
    pvals(i) = ranovatbl.pValue(1);
    fvals(i) = ranovatbl.F(1);
end

figure;
subplot(2,1,1);
plot(f_band, pvals, 'k', 'LineWidth', 1.5);
hold on; yline(0.05, 'r--');
ylim([0 0.1]); xlim([0 30]);
title(['RM-ANOVA P Values at ' target_elec_name]); xlabel('Frequency (Hz)'); ylabel('p-value');

subplot(2,1,2);
plot(f_band, fvals, 'b', 'LineWidth', 1.5);
xlim([0 30]); title('F Statistics'); xlabel('Frequency (Hz)'); ylabel('F-value');
%% 6. 全脑统计地形图 (Repeated Measures)
f_ROI = find((f>=8) & (f<=12)); % Alpha

g1_ch = squeeze(mean(group1_FFT_power(:,:,f_ROI),3));
g2_ch = squeeze(mean(group2_FFT_power(:,:,f_ROI),3));
g3_ch = squeeze(mean(group3_FFT_power(:,:,f_ROI),3));

pvals = zeros(1, EEG.nbchan);
fvals = zeros(1, EEG.nbchan);

fprintf('正在计算全脑通道的重复测量方差分析...\n');

for i=1:EEG.nbchan
    d1 = g1_ch(:,i);
    d2 = g2_ch(:,i);
    d3 = g3_ch(:,i);
    t = table(d1, d2, d3, 'VariableNames', varNames);
    rm = fitrm(t, 'C1-C3 ~ 1', 'WithinDesign', table(Conditions));
    ranovatbl = ranova(rm);
    pvals(i) = ranovatbl.pValue(1);
    fvals(i) = ranovatbl.F(1);
end

figure;
subplot(121); topoplot(fvals, EEG.chanlocs); title('RM-ANOVA F-values','fontsize',12); colorbar;
subplot(122); topoplot(fvals, EEG.chanlocs, 'pmask', pvals < 0.05); title('Significant (p<0.05)','fontsize',12);

%% 7. 提取 SPSS 数据 (宽格式 Wide Format)
% 组内设计在 SPSS 里需要宽格式：每个被试一行，Condition 1, 2, 3 分别占三列

% === 修改开始：输入你想提取的电极名称 ===
target_elec_name_spss = 'Cz'; % 在这里修改，例如 'Cz', 'Pz', 'Fz'
% === 修改结束 ===

% 1. 自动查找电极索引
target_ch = find(strcmpi({EEG.chanlocs.labels}, target_elec_name_spss));

% 安全检查
if isempty(target_ch)
    error(['错误：无法为SPSS提取找到电极 [' target_elec_name_spss ']！请检查拼写。']);
else
    fprintf('正在提取 %s (第 %d 通道) 的数据供 SPSS 使用...\n', target_elec_name_spss, target_ch);
end

% 2. 定义感兴趣的频段 (例如 Alpha 波 8-12Hz)
f_ROI_idx = find((f>=4) & (f<=8)); 

% 3. 提取数据
% 逻辑：取指定通道 -> 对频段内的频率点取平均 -> 得到 (被试x1) 的向量
col_cond1 = squeeze(mean(group1_FFT_power(:, target_ch, f_ROI_idx), 3));
col_cond2 = squeeze(mean(group2_FFT_power(:, target_ch, f_ROI_idx), 3));
col_cond3 = squeeze(mean(group3_FFT_power(:, target_ch, f_ROI_idx), 3));

% 4. 创建表格
% 格式：SubjID | Cond1 | Cond2 | Cond3
SubjID = (1:length(col_cond1))';

% 动态生成列名，带上电极名字，防止混淆
varName1 = ['Cond1_', target_elec_name_spss];
varName2 = ['Cond2_', target_elec_name_spss];
varName3 = ['Cond3_', target_elec_name_spss];

T = table(SubjID, col_cond1, col_cond2, col_cond3, ...
    'VariableNames', {'SubjectID', varName1, varName2, varName3});

% 5. 保存文件 (文件名也带上电极名字)
save_filename = ['EEG_Alpha_WithinSubject_' target_elec_name_spss '.csv'];
writetable(T, save_filename);

fprintf('SPSS 数据已成功保存为: %s\n', save_filename);
disp('请在 SPSS 中使用 Analyze -> General Linear Model -> Repeated Measures');