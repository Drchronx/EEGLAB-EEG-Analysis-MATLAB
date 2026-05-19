%数据导入并分段
clear all; clc; close all
data_path = 'path_to_your_data';%预处理完数据路径
cd(data_path)
files = dir('*.set');
fn = {files.name};
Cond = {  '1', '2', '3', '4'}; %% condition name  mark名称
for i = 1:length(fn)
     setname = fn{i}; %% filename of set file
     EEG = pop_loadset('filename',setname,'filepath',data_path); 
     EEG = eeg_checkset( EEG );
    for j = 1:length(Cond)%对各Mark重新分段
        EEG_new = pop_epoch( EEG, Cond(j), [-0.2  1], 'newname', 'Merged datasets pruned with ICA', 'epochinfo', 'yes'); %% epoch by conditions, input to EEG_new
        EEG_new = eeg_checkset( EEG_new );
        EEG_new = pop_rmbase( EEG_new, [-200     0]); %% baseline correction for EEG_new
        EEG_new = eeg_checkset( EEG_new );
        EEG_avg(i,j,:,:) = squeeze(mean(EEG_new.data,3));  %% average across trials for EEG_new, EEG_avg dimension: subj*cond*channel*time
    end 
end
EEG.times = EEG.times(:, 401:1000);%更新一下时间数据

%%  多因素找显著时间窗

data_test = squeeze(EEG_avg(:,:,28,:)); %% 选电极找显著, data_test: subj*cond*time
%对于每个时间点
for i = 1:size(data_test,3)
    %提取所有被试 所有条件 第i个时间点的数据
    %data_anova sub * cond
    data_anova = squeeze(data_test(:,:,i)); %% select the data at time point i
    %在第i个时间点下进行重复测量方差分析
    %注意 在使用anova_rm函数时 要保证函数在当前路径下 或者set path中
    %此函数只能做单因素重复测量方差分析
    %要求输入的数据组织形式是 被试为行 条件为列，off关闭弹窗
    %输出变量p中的第一个值是条件的主效应
    %table是方差分析的表
    sub = 27; % 被试数量
    F1_con = 2; %F1自变量水平数
    F2_con = 2; %F2自变量水平数
    total_con = F1_con * F2_con ;
    S = repmat([1:sub], 1, total_con).';
    % F1和F2 为自变量的水平
    F1 = sort(repmat([1:F1_con], 1, sub*F2_con)).';
    F2 = repmat(sort(repmat([1:F2_con], 1, sub)), 1, F1_con).';
    %FACTNAMES为变量名称
    FACTNAMES = {'一条件', '二条件'};
    
    %重复测量方差分析
    tab = rm_anova2(data_anova , S, F1, F2, FACTNAMES);

    %输出p值
    p2 =tab{3, 6}%第二条件主效应
    p1 =tab{2, 6}%第一条件主效应
    p =tab{4, 6};%交互作用
    %汇总每次统计下的p值
    P_anovaz1(i) = p1(1); %% save the data from ANOVA
    P_anovaz2(i) = p2(1); 
    P_anova(i) = p(1); 
end

mean_data = squeeze(mean(data_test,1)); %%可提取为画图数据 可以直接MATLAB画图，也可以复制黏贴出来放origion画图  dimension: cond*time
figure; 
subplot(211);plot(EEG.times, mean_data,'linewidth', 1.5); %% waveform for different condition 
set(gca,'YDir','reverse');
axis([-200 1000 -35 25]);
subplot(212);plot(EEG.times,P_anovaz1); axis([-200 1000 0 0.05]); %

figure; 
subplot(211);plot(EEG.times, mean_data,'linewidth', 1.5); %% waveform for different condition 
set(gca,'YDir','reverse');
axis([-200 1000 -35 25]);
subplot(212);plot(EEG.times,P_anovaz2); axis([-200 1000 0 0.05]); %

figure; 
subplot(211);plot(EEG.times, mean_data,'linewidth', 1.5); %% waveform for different condition 
set(gca,'YDir','reverse');
axis([-200 1000 -35 25]);
subplot(212);plot(EEG.times,P_anova); axis([-200 1000 0 0.05]); %


%% ====================================================================================
%   附加功能一：全新2x2布局专业绘图 (最终修正版)
% ====================================================================================
% 描述: 此模块使用您代码前面计算出的结果变量
% (mean_data, P_anovaz1, P_anovaz2, P_anova)
% 来生成一个全新的、更清晰的2x2布局图表。

fprintf('\nPart (附加一): 正在生成2x2布局的结果图表...\n');

% --- 绘图配置 ---
% 您可以在此修改绘图的参数
cfg_plot = struct();
cfg_plot.roi_name = 'CZ'; % 请确保这与您前面分析的电极一致
cfg_plot.factor1_name = '因素A'; % 您可以修改为真实的因素名
cfg_plot.factor2_name = '因素B';
cfg_plot.p_value_threshold = 0.05;

% --- 开始绘图 ---
figure('Name', 'ROI时序分析结果(2x2布局)', 'NumberTitle', 'off', 'Position', [100, 100, 1400, 900]);
sgtitle(sprintf('ROI [ %s ] 的 2x2 ANOVA 结果', cfg_plot.roi_name), 'FontSize', 16, 'FontWeight', 'bold');

% --- 1. 左上角: ERP波形图 (4条线) ---
subplot(2, 2, 1);
hold on;
% mean_data 的维度是: 条件 x 时间
plot(EEG.times, mean_data', 'LineWidth', 2);
hold off;
xline(0, 'k--'); yline(0, 'k--'); set(gca, 'YDir', 'reverse'); grid on;
title('总平均ERP波形', 'FontSize', 14);
xlabel('时间 (ms)'); ylabel('幅值 (\muV)');
legend({'Cond 1', 'Cond 2', 'Cond 3', 'Cond 4'}, 'Location', 'best');
xlim(EEG.times([1, end]));


% --- 2. 右上角: 因素1主效应P值图 ---
subplot(2, 2, 2);
hold on;
plot(EEG.times, P_anovaz1, 'k-', 'LineWidth', 2);
sig_points_1 = find(P_anovaz1 < cfg_plot.p_value_threshold);
plot(EEG.times(sig_points_1), P_anovaz1(sig_points_1), 'r.', 'MarkerSize', 8);
hold off;
yline(cfg_plot.p_value_threshold, 'r--', 'Label', sprintf('p = %.2f', cfg_plot.p_value_threshold));
set(gca, 'YDir', 'reverse'); grid on;
title(sprintf('P值: %s 主效应', cfg_plot.factor1_name), 'FontSize', 14);
xlabel('时间 (ms)'); ylabel('p-value');
xlim(EEG.times([1, end])); ylim([0, 1]);


% --- 3. 左下角: 因素2主效应P值图 ---
subplot(2, 2, 3);
hold on;
plot(EEG.times, P_anovaz2, 'k-', 'LineWidth', 2);
sig_points_2 = find(P_anovaz2 < cfg_plot.p_value_threshold);
plot(EEG.times(sig_points_2), P_anovaz2(sig_points_2), 'r.', 'MarkerSize', 8);
hold off;
yline(cfg_plot.p_value_threshold, 'r--', 'Label', sprintf('p = %.2f', cfg_plot.p_value_threshold));
set(gca, 'YDir', 'reverse'); grid on;
title(sprintf('P值: %s 主效应', cfg_plot.factor2_name), 'FontSize', 14);
xlabel('时间 (ms)'); ylabel('p-value');
xlim(EEG.times([1, end])); ylim([0, 1]);


% --- 4. 右下角: 交互效应P值图 ---
subplot(2, 2, 4);
hold on;
plot(EEG.times, P_anova, 'k-', 'LineWidth', 2);
sig_points_inter = find(P_anova < cfg_plot.p_value_threshold);
plot(EEG.times(sig_points_inter), P_anova(sig_points_inter), 'r.', 'MarkerSize', 8);
hold off;
yline(cfg_plot.p_value_threshold, 'r--', 'Label', sprintf('p = %.2f', cfg_plot.p_value_threshold));
set(gca, 'YDir', 'reverse'); grid on;
title(sprintf('P值: %s x %s 交互效应', cfg_plot.factor1_name, cfg_plot.factor2_name), 'FontSize', 14);
xlabel('时间 (ms)'); ylabel('p-value');
xlim(EEG.times([1, end])); ylim([0, 1]);

%% ====================================================================================
%% ====================================================================================
%% ====================================================================================
%   附加功能二：科学的异常值检测 (最终修正版)
% ====================================================================================
% 描述: 此模块使用您代码前面加载的 EEG_avg 和电极信息，
% 在您指定的时间窗和电极ROI上，进行客观的异常值检测。

% --- 异常值检测配置 ---
cfg_outlier = struct();
cfg_outlier.enabled = true; % 设置为 true 来开启此功能
cfg_outlier.time_window = [500, 600]; % <--- 您想检查异常值的时间窗 (ms)
cfg_outlier.roi_names = {'CZ', 'C1', 'C2', 'PZ', 'P1', 'P2', 'CPZ', 'CP1', 'CP2'}; % <--- 您想检查的电极ROI
cfg_outlier.sd_threshold = 2; % 定义异常值的标准差阈值

% --- 开始异常值检测 ---
if cfg_outlier.enabled
    fprintf('\nPart (附加二): 开始进行探索性的异常值检测...\n');
    fprintf('!!! 警告: 请谨慎解读此部分结果，并遵守科研规范。!!!\n');

    % 准备数据
    outlier_ch_indices = find(ismember({EEG.chanlocs.labels}, cfg_outlier.roi_names));
    outlier_time_indices = EEG.times >= cfg_outlier.time_window(1) & EEG.times <= cfg_outlier.time_window(2);
    if isempty(outlier_ch_indices), error('未找到异常值检测的ROI电极。'); end
    
    fprintf('在ROI ( %s ) 和时间窗 [ %d, %d ] ms 内进行分析...\n', strjoin(cfg_outlier.roi_names, ', '), cfg_outlier.time_window(1), cfg_outlier.time_window(2));
    
    % EEG_avg 维度: 被试 x 条件 x 通道 x 时间
    data_for_outlier = squeeze(mean(mean(EEG_avg(:, :, outlier_ch_indices, outlier_time_indices), 3), 4));

    % 对于2x2设计，最关键的效应是交互效应
    % 交互效应的差值 = (Cond1 - Cond2) - (Cond3 - Cond4)
    interaction_effect_diff = (data_for_outlier(:,1) - data_for_outlier(:,2)) - (data_for_outlier(:,3) - data_for_outlier(:,4));

    % 基于标准差识别异常值
    mean_diff = mean(interaction_effect_diff);
    std_diff = std(interaction_effect_diff);
    z_scores = (interaction_effect_diff - mean_diff) / std_diff;
    outlier_indices = find(abs(z_scores) > cfg_outlier.sd_threshold);

    if isempty(outlier_indices)
        fprintf('\n====== 在 %.1f 个标准差的阈值下，未发现异常值被试。 ======\n', cfg_outlier.sd_threshold);
    else
        fprintf('\n====== 已识别出 %d 名可能的异常值被试 (阈值 > %.1f SD) ======\n', length(outlier_indices), cfg_outlier.sd_threshold);
        outlier_report = table(fn(outlier_indices)', interaction_effect_diff(outlier_indices), z_scores(outlier_indices), 'VariableNames', {'FileName', 'InteractionEffect', 'Z_Score'});
        disp('被识别为异常值的被试信息:');
        disp(outlier_report);
        
        % --- 模拟剔除异常值后的ANOVA结果 ---
        fprintf('\n--- 正在模拟剔除以上被试后的ANOVA结果 (在相同的ROI和时间窗内) ---\n');
        
        num_subjects_total = size(EEG_avg, 1);
        subjects_to_keep = setdiff(1:num_subjects_total, outlier_indices);
        data_for_anova_clean = data_for_outlier(subjects_to_keep, :);
        
        fprintf('剔除后剩余 %d 名被试。\n', size(data_for_anova_clean, 1));
        
        % --- [代码修正] ---
        % 之前这里创建的table使用了默认列名(Var1, Var2)，可能导致fitrm出错。
        % 我们现在明确地为它命名为'Cond1','Cond2'等，以确保变量名能被正确识别。
        num_conds_clean = size(data_for_anova_clean, 2);
        var_names = sprintfc('Cond%d', 1:num_conds_clean);
        data_table_clean = array2table(data_for_anova_clean, 'VariableNames', var_names);

        % 相应地，更新模型公式以使用新的列名
        model_formula = sprintf('Cond1-Cond%d ~ 1', num_conds_clean);
        % --- [修正结束] ---
        
        % 使用标准的fitrm进行模拟检验
        f1_name = 'FactorA'; f2_name = 'FactorB';
        [f1, f2] = ndgrid([1, 2], [1, 2]);
        within_design = table(categorical(f1(:)), categorical(f2(:)), 'VariableNames', {f1_name, f2_name});
        within_model_str = sprintf('%s*%s', f1_name, f2_name);
        
        rm_model_clean = fitrm(data_table_clean, model_formula, 'WithinDesign', within_design);
        ranova_results_clean = ranova(rm_model_clean, 'WithinModel', within_model_str);
        
        disp('模拟剔除异常值被试后的ANOVA结果表:');
        disp(ranova_results_clean);
    end
end
fprintf('\n====== 附加功能模块运行完毕。 ======\n');

%% point-by-point repeated measures of ANOVA across channels
%找到感兴趣时间范围内的采样点的位置信息
test_idx = find((EEG.times>=500)&(EEG.times<=600)); %% define the intervals
%sub * con * ch
data_test = squeeze(mean(EEG_avg(:,:,:,test_idx),4)); %% select the data in [197 217]ms, subj*cond*channel
%对于每个通道
for i = 1:size(data_test,3)
    data_anova = squeeze(data_test(:,:,i)); %% select the data at channel i
    %在第i个时间点下进行重复测量方差分析
    %注意 在使用anova_rm函数时 要保证函数在当前路径下 或者set path中
    %此函数只能做单因素重复测量方差分析
    %要求输入的数据组织形式是 被试为行 条件为列，off关闭弹窗
    %输出变量p中的第一个值是条件的主效应
    %table是方差分析的表
    sub = 28; % 被试数量
    F1_con = 2; %F1自变量水平数
    F2_con = 2; %F2自变量水平数
    total_con = F1_con * F2_con ;
    S = repmat([1:sub], 1, total_con).';
    % F1和F2 为自变量的水平
    F1 = sort(repmat([1:F1_con], 1, sub*F2_con)).';
    F2 = repmat(sort(repmat([1:F2_con], 1, sub)), 1, F1_con).';
    %FACTNAMES为变量名称
    FACTNAMES = {'一条件', '二条件'};
    
    %重复测量方差分析
    tab = rm_anova2(data_anova , S, F1, F2, FACTNAMES);

    %输出p值
    %p =tab{3, 6}第二条件主效应
    %p =tab{2, 6}第一条件主效应
    %输出p值
    p2 =tab{3, 6}%第二条件主效应
    p1 =tab{2, 6}%第一条件主效应
    p =tab{4, 6}%交互作用
    %汇总每次统计下的p值
    P_anovaz11(i) = p1(1); %% save the data from ANOVA
    P_anovaz22(i) = p2(1); 
    %汇总每次统计下的p值
    P_anova2(i) = p(1); %% save the data from ANOVA
    
end

figure;
for i = 1:4
    subplot(1,5,i); 
    topoplot(squeeze(mean(data_test(:,i,:),1)),EEG.chanlocs,'maplimits',[-20 20]); 
end
subplot(1,5,5); topoplot( P_anovaz11,EEG.chanlocs,'maplimits',[0 0.05]);

figure;
for i = 1:4
    subplot(1,5,i); 
    topoplot(squeeze(mean(data_test(:,i,:),1)),EEG.chanlocs,'maplimits',[-20 20]); 
end
subplot(1,5,5); topoplot( P_anovaz22,EEG.chanlocs,'maplimits',[0 0.05]);

figure; 
for i = 1:4
    subplot(1,5,i); 
    topoplot(squeeze(mean(data_test(:,i,:),1)),EEG.chanlocs,'maplimits',[-20 20]); 
end
subplot(1,5,5); topoplot( P_anova2,EEG.chanlocs,'maplimits',[0 0.05]);
%subplot(1,5,5); topoplot( F_anova2,EEG.chanlocs); 

%% 提取脑电统计的数据，选择通道，时间在上环节提取，这些数据放spss里
clear all; clc; close all
data_path = 'path_to_your_data';
cd(data_path)
files = dir('*.set');
fn = {files.name};
%Cond = {  'S 21'  'S 22'  'S 23'  'S 24'  };
Cond = {  '1'  '2'  '3'  '4' }; %% condition name
for i = 1:length(fn)
     setname = fn{i}; %% filename of set file
     EEG = pop_loadset('filename',setname,'filepath',data_path); 
     EEG = eeg_checkset( EEG );
    for j = 1:length(Cond)
        EEG_new = pop_epoch( EEG, Cond(j), [-0.2  1], 'newname', 'Merged datasets pruned with ICA', 'epochinfo', 'yes'); %% epoch by conditions, input to EEG_new
        EEG_new = eeg_checkset( EEG_new );
        EEG_new = pop_rmbase( EEG_new, [-200     0]); %% baseline correction for EEG_new
        EEG_new = eeg_checkset( EEG_new );
        EEG_avg(i,j,:,:) = squeeze(mean(EEG_new.data,3));  %% average across trials for EEG_new, EEG_avg dimension: subj*cond*channel*time
    end 
end
EEG.times = EEG.times(:, 401:1000);%更新一下时间数据方便后面画图

data_test = squeeze(EEG_avg(:,:,28,:)); %% 选电极找显著, data_test: subj*cond*time
mean_data = squeeze(mean(data_test,1)); %%可提取为画图数据 可以直接MATLAB画图，也可以复制黏贴出来放origion画图  dimension: cond*time



test_idx = find((EEG.times>=80)&(EEG.times<=100)); %
data_test = squeeze(mean(EEG_avg(:,:,:,test_idx),4));
data_djtongji= squeeze(data_test(:,:,29));%选择想
% define the intervals
%sub * con * ch提取电极



%% ====================================================================================

%% ====================================================================================
%  Part 8 (修正版): 批量导出 SPSS 统计数据
%  前置要求：请确保前面的数据读取代码已运行，且工作区有 'EEG_avg' 和 'fn' 变量
% ====================================================================================

% --- 1. 关键变量桥接与初始化 (修复报错的核心) ---

% 检查是否存在 EEG_avg (你的读取代码生成的变量)
if exist('EEG_avg', 'var')
    all_subjects_4D_data = EEG_avg; % 统一变量名
else
    error('错误: 找不到变量 EEG_avg。请先运行数据读取部分的代码！');
end

% 自动获取条件数量 (根据数据维度)
% EEG_avg 维度: sub * cond * channel * time
num_conditions = size(all_subjects_4D_data, 2); 

% 创建 subject_info 表格 (使用文件名作为 ID)
if exist('fn', 'var')
    SubjID = fn'; % 将文件名转置为列向量
else
    % 如果没有文件名，就用数字 ID 生成
    num_subs = size(all_subjects_4D_data, 1);
    SubjID = arrayfun(@(x) sprintf('Sub_%02d', x), (1:num_subs)', 'UniformOutput', false);
end
subject_info = table(SubjID); %这就定义了报错的那个变量

% --- 2. 参数配置 ---
cfg = []; % 初始化结构体，防止报错
%cfg.spss_channels = {'CPZ', 'P1', 'P2', 'CPZ', 'CP1', 'CP2'}; 
cfg.spss_channels = {'PZ'}; 
cfg.spss_time_window = [400, 500]; 
cfg.spss_output_filename = 'SPSS_Data_Export.csv'; 

fprintf('\nPart 8: 正在为 SPSS 导出 %d 个电极的数据 (时间窗: %d-%d ms)...\n', ...
    length(cfg.spss_channels), cfg.spss_time_window(1), cfg.spss_time_window(2));

% --- 3. 开始处理 ---
% 找到时间窗索引 (注意：这里需要确保 EEG 结构体存在且包含 times)
if ~exist('EEG', 'var') || ~isfield(EEG, 'times')
     error('错误: 找不到 EEG.times 变量。请确保数据已加载。');
end
spss_time_idx = EEG.times >= cfg.spss_time_window(1) & EEG.times <= cfg.spss_time_window(2);

% 初始化输出表格
output_table_spss = subject_info;

% 循环每一个电极
for i = 1:length(cfg.spss_channels)
    curr_chan = cfg.spss_channels{i};
    
    % 找电极索引 (注意：需要 chanlocs 变量)
    % 如果你的工作区没有 chanlocs，这里会报错。通常 EEG.chanlocs 存在。
    if exist('chanlocs', 'var')
        current_locs = chanlocs;
    elseif isfield(EEG, 'chanlocs')
        current_locs = EEG.chanlocs;
    else
        error('错误: 找不到电极定位信息 (chanlocs)。');
    end
    
    ch_idx = find(strcmp({current_locs.labels}, curr_chan));
    
    if isempty(ch_idx)
        fprintf('  警告: 数据中找不到电极 [ %s ]，跳过。\n', curr_chan);
        continue;
    end
    
    % 提取该电极、该时间窗的数据
    % 维度: 被试 x 条件 x 1 x 时间 -> squeeze -> 被试 x 条件
    chan_data_avg = squeeze(mean(all_subjects_4D_data(:, :, ch_idx, spss_time_idx), 4));
    
    % 将每个条件作为一列加入表格
    for c = 1:num_conditions
        % 构建列名: 电极_Cond1, 电极_Cond2 ...
        col_name = sprintf('%s_Cond%d', curr_chan, c);
        % 注意：chan_data_avg 可能是 (被试x条件)，直接取列
        output_table_spss.(col_name) = chan_data_avg(:, c);
    end
    
    fprintf('  -> 已添加电极 [ %s ]\n', curr_chan);
end

% --- 4. 保存 ---
try
    writetable(output_table_spss, cfg.spss_output_filename);
    fprintf('====== SPSS 数据已保存至: %s ======\n', cfg.spss_output_filename);
catch ME
    fprintf('错误: 保存失败，文件可能被占用。\n');
end
%% Part 10 (修正版): 批量导出 Origin/Excel 画图数据 (总平均波形)
% ====================================================================================
% 功能：导出指定电极的“总平均”波形数据 (所有被试平均)。
% 格式：输出 CSV 文件。
%       第 1 行：时间轴 (Time Points)
%       第 2 行：条件 1 的波形数值
%       第 3 行：条件 2 的波形数值 ... 以此类推
% 用法：直接拖入 Origin，或者在 Excel 中打开复制。
% ====================================================================================

% --- 1. 变量准备 (防报错) ---
if exist('EEG_avg', 'var')
    plot_data_source = EEG_avg; % 使用读取好的数据
else
    error('错误: 找不到 EEG_avg 变量，请先运行数据读取代码。');
end

if exist('EEG', 'var') && isfield(EEG, 'times')
    curr_times = EEG.times; % 获取时间轴
else
    error('错误: 找不到 EEG.times 时间轴数据。');
end

% 确保 chanlocs 存在
if exist('chanlocs', 'var')
    curr_locs = chanlocs;
elseif isfield(EEG, 'chanlocs')
    curr_locs = EEG.chanlocs;
else
    error('错误: 找不到电极定位信息 (chanlocs)。');
end

% --- 2. 参数配置 ---
% 在这里填入所有需要画图的电极
cfg = [];
cfg.plot_channels = {'CZ', 'PZ', 'CPZ', 'FZ'}; 

% 导出文件的前缀
cfg.plot_filename_prefix = 'Origin_Plot_GrandAvg'; 

fprintf('\nPart 10: 正在导出 Origin 画图数据 (%d 个电极)...\n', length(cfg.plot_channels));

% --- 3. 循环处理每个电极 ---
for i = 1:length(cfg.plot_channels)
    curr_chan = cfg.plot_channels{i};
    
    % 找电极索引
    ch_idx = find(strcmp({curr_locs.labels}, curr_chan));
    
    if isempty(ch_idx)
        fprintf('  警告: 找不到电极 [ %s ]，跳过。\n', curr_chan);
        continue; 
    end
    
    % --- 核心计算 ---
    % 数据维度: 被试(Sub) x 条件(Cond) x 通道(Chan) x 时间(Time)
    
    % 1. 提取该通道的所有数据
    % 结果: 被试 x 条件 x 1 x 时间
    chan_data_all = plot_data_source(:, :, ch_idx, :);
    
    % 2. 对“被试”维度求平均 (Grand Average)
    % 结果: 1 x 条件 x 1 x 时间
    chan_data_avg = mean(chan_data_all, 1);
    
    % 3. 压缩维度 (Squeeze)
    % 结果: 条件(Rows) x 时间(Cols)  (注意: squeeze后维度取决于数据，通常是 Cond x Time)
    wave_data = squeeze(chan_data_avg);
    
    % 检查维度方向，确保行是条件，列是时间
    % 如果条件数很少(比如4)，时间点很多(比如600)，Matlab通常会弄成 4x600。
    % 如果反了，就转置一下。
    if size(wave_data, 2) ~= length(curr_times)
        wave_data = wave_data'; 
    end
    
    % --- 4. 拼接输出矩阵 ---
    % 第1行是时间，下面几行是数据
    output_matrix = [curr_times; wave_data];
    
    % --- 5. 保存文件 ---
    fname = sprintf('%s_%s.csv', cfg.plot_filename_prefix, curr_chan);
    
    try
        writematrix(output_matrix, fname);
        fprintf('  -> [ %s ] 导出成功: %s\n', curr_chan, fname);
    catch
        fprintf('  -> [ %s ] 保存失败 (文件可能被占用)。\n', curr_chan);
    end
end

fprintf('====== 所有画图数据导出完成 ======\n');
fprintf('提示: CSV文件中，第1行为时间轴(X轴)，第2行起为各条件(Cond1, Cond2...)的波形(Y轴)。\n');
%% Part 11: [新功能] 便捷绘制各条件地形图 (Topography)
%% Part 11 (修正版): 绘制指定时间窗的地形图
% ====================================================================================
% 功能：指定一个时间窗，快速画出所有组、所有条件的电压地形图。
%% Part 11 (最终修正版): 绘制指定时间窗的地形图 [Image of EEG topographic map example]
% ====================================================================================
% 功能：指定一个时间窗，快速画出所有组、所有条件的电压地形图。
% 自动修复：会自动检测缺失的 chanlocs 和数据变量。
% ====================================================================================

% --- 0. 变量自动修复 (关键步骤) ---

% 1. 修复数据源 (all_subjects_4D_data)
if exist('EEG_avg', 'var') && ~exist('all_subjects_4D_data', 'var')
    all_subjects_4D_data = EEG_avg; 
    fprintf('系统提示: 已自动加载 EEG_avg 作为绘图数据。\n');
elseif ~exist('all_subjects_4D_data', 'var')
    error('严重错误: 找不到数据变量 (EEG_avg 或 all_subjects_4D_data)。请重新运行数据读取部分。');
end

% 2. 修复电极定位 (chanlocs)
if exist('chanlocs', 'var')
    % 变量已存在，无需操作
    fprintf('系统提示: 检测到电极位置信息 (chanlocs)。\n');
elseif exist('EEG', 'var') && isfield(EEG, 'chanlocs')
    chanlocs = EEG.chanlocs; % 从 EEG 结构体中提取
    fprintf('系统提示: 已从 EEG 结构体中提取 chanlocs。\n');
else
    % 如果连 EEG 变量都没有，尝试加载一个示例数据集来获取位置（通常不建议，除非你确实没办法）
    error('严重错误: 找不到 chanlocs (电极位置)。请确保 EEG 变量已加载 (例如运行 EEG = pop_loadset...)。');
end

% --- 1. 参数配置 ---
cfg = [];
cfg.topo_plot_window = [350, 450];    % 时间窗 (ms)
cfg.topo_plot_limits = [-5, 5];       % 颜色范围 (uV)

% 分组信息默认处理
if ~exist('subject_info', 'var')
    % 如果没有分组表，假定所有人是一组
    num_subs = size(all_subjects_4D_data, 1);
    subject_info = table(ones(num_subs,1), 'VariableNames', {'Group'});
    cfg.groups(1).name = 'All Subjects';
else
    cfg.groups(1).name = 'Group 1'; % 请根据实际修改
    cfg.groups(2).name = 'Group 2'; 
end

fprintf('\nPart 11: 正在生成时间窗 [ %d, %d ] ms 的地形图...\n', cfg.topo_plot_window(1), cfg.topo_plot_window(2));

% --- 2. 数据准备 ---
% 检查时间轴
if ~exist('EEG', 'var') || ~isfield(EEG, 'times')
    error('错误: 找不到时间轴 EEG.times。');
end

% 找到时间索引
t_idx = EEG.times >= cfg.topo_plot_window(1) & EEG.times <= cfg.topo_plot_window(2);

% 如果时间窗没有选到任何点（比如范围写错了），报错
if sum(t_idx) == 0
    error('错误: 选择的时间窗 [%d, %d] 超出了数据时间范围 (%.1f 到 %.1f ms)。', ...
        cfg.topo_plot_window(1), cfg.topo_plot_window(2), EEG.times(1), EEG.times(end));
end

% 计算该时间窗内的平均值 (压缩时间维)
% 结果维度: 被试 x 条件 x 通道
data_topo_avg_time = squeeze(mean(all_subjects_4D_data(:, :, :, t_idx), 4));

% --- 3. 绘图循环 ---
figure('Name', '地形图概览', 'Color', 'w', 'Position', [100, 100, 1200, 500]);
sgtitle(sprintf('地形图: %d-%d ms', cfg.topo_plot_window(1), cfg.topo_plot_window(2)), ...
    'FontSize', 16, 'FontWeight', 'bold');

num_conds = size(all_subjects_4D_data, 2); 
unique_groups = unique(subject_info.Group);
num_groups = length(unique_groups);

plot_counter = 1;

for g = 1:num_groups
    group_val = unique_groups(g);
    
    % 获取组名
    if g <= length(cfg.groups)
        curr_group_name = cfg.groups(g).name;
    else
        curr_group_name = sprintf('Group %d', group_val);
    end
    
    % 找到该组被试索引
    group_indices = find(subject_info.Group == group_val);
    if isempty(group_indices); continue; end
    
    for c = 1:num_conds
        % 取出数据: 该组、该条件
        temp_data = data_topo_avg_time(group_indices, c, :);
        
        % 对被试维度求平均 -> 得到 (1 x 通道)
        if length(group_indices) > 1
            topo_data_final = squeeze(mean(temp_data, 1));
        else
            topo_data_final = squeeze(temp_data); % 只有一个被试直接降维
        end
        
        % 绘图位置
        subplot(num_groups, num_conds, plot_counter);
        
        % *** 核心绘图函数 ***
        try
            topoplot(topo_data_final, chanlocs, ...
                'maplimits', cfg.topo_plot_limits, ...
                'style', 'map', ...
                'electrodes', 'on', ...
                'headrad', 'rim', ...
                'shading', 'interp');
        catch ME
            % 捕获 topoplot 内部错误（通常是电极名称不匹配）
            warning('Topoplot 出错: %s', ME.message);
            title('Error Plotting');
            plot_counter = plot_counter + 1;
            continue;
        end
            
        title(sprintf('%s - Cond %d', curr_group_name, c), 'FontSize', 11);
        
        % 最后一个图加色标
        if c == num_conds
            cbar = colorbar;
            cbar.Label.String = 'uV';
        end
        
        plot_counter = plot_counter + 1;
    end
end

fprintf('====== 地形图绘制完成 ======\n');
%% 画地形图
N2_interval=find((EEG.times>=250)&(EEG.times<=350)); %% N2 interval
P2_interval=find((EEG.times>=400)&(EEG.times<=600)); %% P2 interval
%EEG_avg 4D sub * con * ch * t
%提取所有被试 所有条件 所有通道 感兴趣时间范围内的数据
%沿着时间做平均 沿着被试做平均
%-> con * ch
N2_amplitude=squeeze(mean(mean(EEG_avg(:,:,:,N2_interval),4),1)); %% N2 amplitude for each subject, condition, and channels
P2_amplitude=squeeze(mean(mean(EEG_avg(:,:,:,P2_interval),4),1)); %% P2 amplitude for each subject, condition, and channels
%将两个条件的数据进行纵向拼接
all_amp = [N2_amplitude;P2_amplitude];
%将数据取绝对值
abs_all_amp = abs(all_amp);
%找最大值,把并向下取整
max_val = floor(0.7*max(abs_all_amp(:)));
%生成一个空画布
figure; %% divide the panel into 4 rows and 2 colums
%对于每个条件
for i = 1:4
    N2_data = N2_amplitude(i,:); %% average across subjects
    %绘制2行四列的第i张
    subplot(2,4,i); 
    topoplot(N2_data,EEG.chanlocs,'maplimits',[-1*max_val max_val]); title(['L',num2str(i),' N2 Amplitude'],'fontsize',16); %% plot N2 scalp map (group-level)
    P2_data = P2_amplitude(i,:); %% average across subjets
    %绘制两行四列的i+4张
    subplot(2,4,i+4); 
    topoplot(P2_data,EEG.chanlocs,'maplimits',[-1*max_val max_val]); title(['L',num2str(i),' P2 Amplitude'],'fontsize',16); %% plot P2 scamp map (group-level)
end