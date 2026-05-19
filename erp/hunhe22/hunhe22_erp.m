

%  此脚本提供了一个完整的EEG数据分析流程，用于混合实验设计。
%  1. Part 1-4: 执行逐时间点的混合设计ANOVA，找到显著的效应和时间窗。
%  2. Part 5:   可视化时间序列的P值和指定ROI的ERP波形。
%  3. Part 6-7: 在指定时间窗内，执行逐通道的ANOVA，并通过地形图可视化结果。
%  4. Part 8:   提取指定ROI和时间窗的平均幅值，以供发表或在其他软件中使用。
% ====================================================================================

clear all; clc; close all;

%% Part 1: 实验配置 (Experiment Configuration)
% ====================================================================================
% --- 请在此处修改您的实验参数 ---
fprintf('Part 1: 正在配置实验参数...\n');

cfg = struct();

% --- 实验组设计与路径配置 ---
cfg.groups(1).name = 'zy';
cfg.groups(1).data_dir = 'path_to_your_data';
cfg.groups(1).conditions = {'1', '2'};

cfg.groups(2).name = 'bzy';
cfg.groups(2).data_dir = 'path_to_your_data';
cfg.groups(2).conditions = {'3', '4'};

% --- 组内因素定义 ---
cfg.within_factors.names = {'rw'};
cfg.within_factors.levels = [1, 2];

% --- 分析参数 ---
cfg.p_value_threshold = 0.05;      % 显著性水平阈值
cfg.epoch_time = [-0.2, 1.0];      % 分段的时间窗口
cfg.baseline_time = [-200, 0];     % 基线校正的时间窗口 (ms)






%% Part 2: 统一数据加载与准备 (加载所有通道)
% ====================================================================================
fprintf('\nPart 2: 正在加载所有被试的完整通道数据...\n');

all_subjects_4D_data = []; % 存储4D数据: 被试 x 条件 x 通道 x 时间
subject_info = table();
subject_counter = 0;
chanlocs = [];
eeg_times = [];

for g = 1:length(cfg.groups)
    group_name = cfg.groups(g).name;
    data_dir = cfg.groups(g).data_dir;
    conditions = cfg.groups(g).conditions;
    file_list = dir(fullfile(data_dir, '*.set'));
    
    fprintf('--- 正在处理组: %s (%d 个文件) ---\n', group_name, length(file_list));
    
    for s = 1:length(file_list)
        fprintf('  加载被试: %s\n', file_list(s).name);
        EEG = pop_loadset('filename', file_list(s).name, 'filepath', data_dir);
        
        if isempty(chanlocs)
            chanlocs = EEG.chanlocs;
        end

        subject_all_cond_data = [];
        is_subject_valid = true;
        
        for c = 1:length(conditions)
            try
                EEG_epoch = pop_epoch(EEG, {conditions{c}}, cfg.epoch_time, 'epochinfo', 'yes');
                if EEG_epoch.trials == 0
                    fprintf('  警告: 被试 %s 缺少条件 %s 的有效试次，将跳过该被试。\n', file_list(s).name, conditions{c});
                    is_subject_valid = false;
                    break;
                end
                EEG_epoch = pop_rmbase(EEG_epoch, cfg.baseline_time);
                
                % 平均trials，保留所有通道 (通道 x 时间)
                avg_data = squeeze(mean(EEG_epoch.data, 3));
                subject_all_cond_data(c, :, :) = avg_data;
                
                if isempty(eeg_times)
                    eeg_times = EEG_epoch.times;
                end
            catch ME
                fprintf('  警告: 被试 %s 在处理条件 %s 时出错，将跳过该被试。\n', file_list(s).name, ME.message);
                is_subject_valid = false;
                break;
            end
        end
        
        if is_subject_valid
            subject_counter = subject_counter + 1;
            all_subjects_4D_data(subject_counter, :, :, :) = subject_all_cond_data;
            subject_info.SubjectID(subject_counter) = subject_counter;
            subject_info.Group(subject_counter) = g;
            subject_info.GroupName{subject_counter} = cfg.groups(g).name;
        end
    end
end

fprintf('\n====== 数据加载完成! ======\n');
[total_subjects, num_conditions, num_channels, num_timepoints] = size(all_subjects_4D_data);
fprintf('数据维度: %d 被试 x %d 条件 x %d 通道 x %d 时间点\n', total_subjects, num_conditions, num_channels, num_timepoints);
if total_subjects == 0; error('致命错误: 未加载任何有效被试数据。'); end
disp('被试信息概览:');
disp(subject_info);


%% Part 3: 数据完整性校验 (无需修改)
% ====================================================================================
% (此部分代码与之前相同，保持不变)
fprintf('\nPart 3: 正在进行数据完整性校验...\n');
for g = 1:length(cfg.groups)
    group_indices = find(subject_info.Group == g);
    if length(group_indices) < 2
        fprintf('警告: 组 %d (%s) 少于2名被试，无法进行组间比较。\n', g, cfg.groups(g).name);
        continue;
    end
    group_data_sample = squeeze(all_subjects_4D_data(group_indices, 1, 1, :));
    if all(var(group_data_sample, 0, 1) < 1e-12)
        error('致命错误: 第 %d 组 (%s) 的所有被试数据完全相同!', g, cfg.groups(g).name);
    end
end
fprintf('====== 数据校验通过。 ======\n');

%% Part 4: 逐时间点混合设计ANOVA (在指定ROI上)
% ====================================================================================
% --- Part 4 & 5: ROI时序分析参数 ---
cfg.roi_channel_name = 'CZ';       % 您主要感兴趣的ROI电极，用于逐时间点分析
fprintf('\nPart 4: 在ROI电极 [ %s ] 上执行逐时间点ANOVA...\n', cfg.roi_channel_name);

% 从4D数据中提取出ROI的数据
roi_channel_index = find(strcmp({chanlocs.labels}, cfg.roi_channel_name));
if isempty(roi_channel_index)
    error('错误: 在电极列表中未找到指定的ROI电极 "%s"', cfg.roi_channel_name);
end
all_subjects_roi_data = squeeze(all_subjects_4D_data(:, :, roi_channel_index, :));

% --- 准备ANOVA模型 ---
within_factor_name = cfg.within_factors.names{1};
within_design = table(categorical(cfg.within_factors.levels'), 'VariableNames', {within_factor_name});
model_formula = sprintf('Cond1-Cond%d ~ Group', num_conditions);

% --- 预分配结果 ---
effect_names = {'Group', within_factor_name, ['Group_' within_factor_name]};
results_roi.p_values = array2table(nan(num_timepoints, length(effect_names)), 'VariableNames', effect_names);
results_roi.f_values = array2table(nan(num_timepoints, length(effect_names)), 'VariableNames', effect_names);

% --- 逐时间点循环 (与之前修正版相同) ---
for t = 1:num_timepoints
    data_t = squeeze(all_subjects_roi_data(:, :, t));
    data_table = array2table(data_t, 'VariableNames', sprintfc('Cond%d', 1:num_conditions));
    data_table.Group = categorical(subject_info.Group);
    try
        rm_model = fitrm(data_table, model_formula, 'WithinDesign', within_design);
        ranova_results = ranova(rm_model, 'WithinModel', within_factor_name);
        
        if any(strcmp(ranova_results.Properties.RowNames, 'Group'))
            results_roi.p_values.Group(t) = ranova_results{'Group', 'pValue'};
            results_roi.f_values.Group(t) = ranova_results{'Group', 'F'};
        end
        emotion_main_effect_name = ['(Intercept):' within_factor_name];
        if any(strcmp(ranova_results.Properties.RowNames, emotion_main_effect_name))
            results_roi.p_values.(within_factor_name)(t) = ranova_results{emotion_main_effect_name, 'pValue'};
            results_roi.f_values.(within_factor_name)(t) = ranova_results{emotion_main_effect_name, 'F'};
        end
        interaction_effect_name = ['Group:' within_factor_name];
        if any(strcmp(ranova_results.Properties.RowNames, interaction_effect_name))
            results_roi.p_values.(['Group_' within_factor_name])(t) = ranova_results{interaction_effect_name, 'pValue'};
            results_roi.f_values.(['Group_' within_factor_name])(t) = ranova_results{interaction_effect_name, 'F'};
        end
    catch
        continue; % 如果单个点出错，跳过
    end
end
fprintf('====== ROI逐时间点分析完成! ======\n');


%% Part 5: 可视化 - ROI波形图与分离的P值曲线 (V2 - 2x2布局)
% ====================================================================================
fprintf('\nPart 5: 正在生成2x2布局的ROI分析结果图表...\n');

% --- 准备波形图数据 ---
% 1. 找到每个组的被试索引
group1_indices = find(subject_info.Group == 1);
group2_indices = find(subject_info.Group == 2);

% 2. 分别计算每个组在每个条件下的平均ERP波形
%    all_subjects_roi_data 维度: 被试 x 条件 x 时间
g1_c1_erp = squeeze(mean(all_subjects_roi_data(group1_indices, 1, :), 1));
g1_c2_erp = squeeze(mean(all_subjects_roi_data(group1_indices, 2, :), 1));
g2_c1_erp = squeeze(mean(all_subjects_roi_data(group2_indices, 1, :), 1));
g2_c2_erp = squeeze(mean(all_subjects_roi_data(group2_indices, 2, :), 1));

% --- 开始绘图 ---
figure('Name', 'ROI时序分析结果(2x2布局)', 'NumberTitle', 'off', 'Position', [100, 100, 1400, 900]);
sgtitle(sprintf('ROI [ %s ] 的详细分析结果', cfg.roi_channel_name), 'FontSize', 16, 'FontWeight', 'bold');

% --- 1. 左上角: ERP波形图 (4条线) ---
subplot(2, 2, 1);
hold on;
plot(eeg_times, g1_c1_erp, 'LineWidth', 2, 'LineStyle', '-');
plot(eeg_times, g2_c1_erp, 'LineWidth', 2, 'LineStyle', '-');
plot(eeg_times, g1_c2_erp, 'LineWidth', 2, 'LineStyle', '--');
plot(eeg_times, g2_c2_erp, 'LineWidth', 2, 'LineStyle', '--');
hold off;
xline(0, 'k--');
yline(0, 'k--');
set(gca, 'YDir', 'reverse', 'FontSize', 10);
grid on;
title('总平均ERP波形', 'FontSize', 14);
xlabel('时间 (ms)');
ylabel('幅值 (\muV)');
legend({
    sprintf('%s - Cond 1', cfg.groups(1).name), ...
    sprintf('%s - Cond 1', cfg.groups(2).name), ...
    sprintf('%s - Cond 2', cfg.groups(1).name), ...
    sprintf('%s - Cond 2', cfg.groups(2).name)
}, 'Location', 'best');
xlim(cfg.epoch_time * 1000);

% --- 提取各效应的P值，方便后续绘图 ---
p_group = results_roi.p_values.Group;
p_emotion = results_roi.p_values.(within_factor_name);
p_interaction = results_roi.p_values.(['Group_' within_factor_name]);

% --- 2. 右上角: 组别主效应P值图 ---
subplot(2, 2, 2);
hold on;
plot(eeg_times, p_group, 'k-', 'LineWidth', 2);
sig_points = find(p_group < cfg.p_value_threshold);
plot(eeg_times(sig_points), p_group(sig_points), 'r.', 'MarkerSize', 8);
hold off;
yline(cfg.p_value_threshold, 'r--', 'LineWidth', 1.5, 'Label', sprintf('p = %.2f', cfg.p_value_threshold));
set(gca, 'YDir', 'reverse', 'FontSize', 10);
grid on;
title('P值: 组别 (Group) 主效应', 'FontSize', 14);
xlabel('时间 (ms)');
ylabel('p-value');
xlim(cfg.epoch_time * 1000);
ylim([0, 1]);

% --- 3. 左下角: 情绪主效应P值图 ---
subplot(2, 2, 3);
hold on;
plot(eeg_times, p_emotion, 'k-', 'LineWidth', 2);
sig_points = find(p_emotion < cfg.p_value_threshold);
plot(eeg_times(sig_points), p_emotion(sig_points), 'r.', 'MarkerSize', 8);
hold off;
yline(cfg.p_value_threshold, 'r--', 'LineWidth', 1.5, 'Label', sprintf('p = %.2f', cfg.p_value_threshold));
set(gca, 'YDir', 'reverse', 'FontSize', 10);
grid on;
title(sprintf('P值: %s 主效应', within_factor_name), 'FontSize', 14);
xlabel('时间 (ms)');
ylabel('p-value');
xlim(cfg.epoch_time * 1000);
ylim([0, 1]);

% --- 4. 右下角: 交互效应P值图 ---
subplot(2, 2, 4);
hold on;
plot(eeg_times, p_interaction, 'k-', 'LineWidth', 2);
sig_points = find(p_interaction < cfg.p_value_threshold);
plot(eeg_times(sig_points), p_interaction(sig_points), 'r.', 'MarkerSize', 8);
hold off;
yline(cfg.p_value_threshold, 'r--', 'LineWidth', 1.5, 'Label', sprintf('p = %.2f', cfg.p_value_threshold));
set(gca, 'YDir', 'reverse', 'FontSize', 10);
grid on;
title(sprintf('P值: Group x %s 交互效应', within_factor_name), 'FontSize', 14);
xlabel('时间 (ms)');
ylabel('p-value');
xlim(cfg.epoch_time * 1000);
ylim([0, 1]);

%% Part 6: 事后分析 - 逐通道混合ANOVA (在指定时间窗内)
% ====================================================================================
% --- Part 6 & 7: 事后分析与地形图参数 ---
% !! 根据Part 5的P值图结果，在这里定义您感兴趣的时间窗 (单位: ms) !!
cfg.posthoc_time_window = [170, 210];
cfg.topo_map_limits = [-3, 3];     % 地形图的电压范围 (uV)


fprintf('\nPart 6: 在时间窗 [ %d, %d ] ms 内执行逐通道ANOVA...\n', cfg.posthoc_time_window(1), cfg.posthoc_time_window(2));

% 找到时间窗对应的采样点
time_idx = eeg_times >= cfg.posthoc_time_window(1) & eeg_times <= cfg.posthoc_time_window(2);
% 在该时间窗内对数据进行平均
data_time_avg = squeeze(mean(all_subjects_4D_data(:, :, :, time_idx), 4)); % 维度: 被试 x 条件 x 通道

% --- 预分配结果 ---
results_channels.p_values = array2table(nan(num_channels, length(effect_names)), 'VariableNames', effect_names);
results_channels.f_values = array2table(nan(num_channels, length(effect_names)), 'VariableNames', effect_names);

% --- 逐通道循环 ---
for ch = 1:num_channels
    fprintf('  处理通道 %d / %d: %s\n', ch, num_channels, chanlocs(ch).labels);
    data_ch = squeeze(data_time_avg(:, :, ch));
    data_table = array2table(data_ch, 'VariableNames', sprintfc('Cond%d', 1:num_conditions));
    data_table.Group = categorical(subject_info.Group);
    
    try
        rm_model = fitrm(data_table, model_formula, 'WithinDesign', within_design);
        ranova_results = ranova(rm_model, 'WithinModel', within_factor_name);
        
        if any(strcmp(ranova_results.Properties.RowNames, 'Group'))
            results_channels.p_values.Group(ch) = ranova_results{'Group', 'pValue'};
            results_channels.f_values.Group(ch) = ranova_results{'Group', 'F'};
        end
        emotion_main_effect_name = ['(Intercept):' within_factor_name];
        if any(strcmp(ranova_results.Properties.RowNames, emotion_main_effect_name))
            results_channels.p_values.(within_factor_name)(ch) = ranova_results{emotion_main_effect_name, 'pValue'};
            results_channels.f_values.(within_factor_name)(ch) = ranova_results{emotion_main_effect_name, 'F'};
        end
        interaction_effect_name = ['Group:' within_factor_name];
        if any(strcmp(ranova_results.Properties.RowNames, interaction_effect_name))
            results_channels.p_values.(['Group_' within_factor_name])(ch) = ranova_results{interaction_effect_name, 'pValue'};
            results_channels.f_values.(['Group_' within_factor_name])(ch) = ranova_results{interaction_effect_name, 'F'};
        end
    catch
        continue;
    end
end
fprintf('====== 逐通道分析完成! ======\n');


%% Part 7: 可视化 - 电压地形图 & P值地形图 (V3 - 按组和条件分别绘制)
% ====================================================================================
fprintf('\nPart 7: 正在生成地形图...\n');

% --- [新功能] 按组和条件分别计算平均电压 ---
% 1. 找到每个组的被试索引
group1_indices = find(subject_info.Group == 1);
group2_indices = find(subject_info.Group == 2);

% 2. 分别计算每个组在每个条件下的平均数据
%    data_time_avg 的维度是: 被试 x 条件 x 通道
g1_c1_topo = squeeze(mean(data_time_avg(group1_indices, 1, :), 1));
g1_c2_topo = squeeze(mean(data_time_avg(group1_indices, 2, :), 1));
g2_c1_topo = squeeze(mean(data_time_avg(group2_indices, 1, :), 1));
g2_c2_topo = squeeze(mean(data_time_avg(group2_indices, 2, :), 1));

% 将要绘制的数据和标题放入元胞数组(cell array)，方便后面循环绘图
topo_data_to_plot = {g1_c1_topo, g2_c1_topo, g1_c2_topo, g2_c2_topo};
topo_titles = {
    sprintf('%s - Cond 1', cfg.groups(1).name), ...
    sprintf('%s - Cond 1', cfg.groups(2).name), ...
    sprintf('%s - Cond 2', cfg.groups(1).name), ...
    sprintf('%s - Cond 2', cfg.groups(2).name)
};
% --- [新功能结束] ---

% 从分析结果中获取效应的数量
effect_names = results_channels.p_values.Properties.VariableNames;
num_effects = length(effect_names);

figure('Name', '地形图分析', 'NumberTitle', 'off', 'Position', [100, 100, 1600, 800]); % 稍微加宽画布以容纳所有图
sgtitle(sprintf('时间窗 [ %d, %d ] ms 的地形图结果', cfg.posthoc_time_window(1), cfg.posthoc_time_window(2)), 'FontSize', 16, 'FontWeight', 'bold');

% 确定子图的列数，以美观地容纳所有电压图和P值图
num_cols = max(4, num_effects); % 列数取4（电压图数量）和效应数量中的较大值

% --- 绘制电压地形图 (现在是4个独立的图) ---
% 第一行：电压地形图
for i = 1:4
    subplot(2, num_cols, i);
    topoplot(topo_data_to_plot{i}, chanlocs, 'maplimits', cfg.topo_map_limits, 'style', 'map', 'electrodes', 'off');
    title(topo_titles{i}, 'FontSize', 12);
    colorbar;
end

% --- 绘制P值地形图 (保持不变，绘制在第二行) ---
% 第二行：P值地形图
for i = 1:num_effects
    current_effect = effect_names{i};
    p_vals_ch = results_channels.p_values.(current_effect);
    
    subplot(2, num_cols, num_cols + i);
    % 找到并标记显著的电极
    sig_channels = find(p_vals_ch < cfg.p_value_threshold);
    topoplot(p_vals_ch, chanlocs, 'maplimits', [0, 0.1], 'style', 'map', 'emarker2', {sig_channels, '*', 'w', 10, 1}, 'electrodes', 'off');
    title(sprintf('P值图: %s', strrep(current_effect, '_', ':')), 'FontSize', 12);
    colorbar;
end



%% Part 8: (修改版) 多通道数据提取 (用于SPSS或其他软件)
% ====================================================================================
% 功能：同时导出多个电极在指定时间窗内的平均幅值。
% 输出：CSV文件，每个电极的每个条件都会生成独立的一列。
% ====================================================================================

% --- Part 8: 数据提取参数 ---
% 【修改处】在这里用大括号 {'A', 'B'} 填入您想导出的所有电极
cfg.extract_channels = {'C1', 'CZ', 'C2', 'P1', 'P2', 'PZ'}; 

% 提取的时间窗 (ms)
cfg.extract_time_window = [600, 800]; 

% 输出文件名
cfg.extract_output_filename = 'Multi_Channel_ERP_Data.csv'; 

fprintf('\nPart 8: 正在提取 %d 个电极的平均幅值...\n', length(cfg.extract_channels));

% 1. 找到时间窗索引
extract_time_idx = eeg_times >= cfg.extract_time_window(1) & eeg_times <= cfg.extract_time_window(2);

% 2. 初始化输出表格 (先放入被试信息)
output_table = subject_info;

% 3. 循环处理列表中的每一个电极
for i = 1:length(cfg.extract_channels)
    curr_chan_name = cfg.extract_channels{i};
    
    % 查找该电极在数据中的索引
    ch_idx = find(strcmp({chanlocs.labels}, curr_chan_name));
    
    if isempty(ch_idx)
        fprintf('  警告: 未找到电极 "%s"，跳过该电极。\n', curr_chan_name);
        continue; 
    end
    
    % 提取数据: 
    % 原始维度: 被试 x 条件 x 通道 x 时间
    % 第一步: 选定特定通道 -> 被试 x 条件 x 1 x 时间
    % 第二步: 选定时间窗并求平均 -> 被试 x 条件 x 1
    chan_data = squeeze(mean(all_subjects_4D_data(:, :, ch_idx, extract_time_idx), 4)); 
    
    % 将该电极的数据添加到表格中
    for c = 1:num_conditions
        % 构建列名，例如: P1_Cond1, P1_Cond2
        col_name = sprintf('%s_Cond%d', curr_chan_name, c);
        
        % 添加列 (Matlab table会自动匹配行数)
        output_table.(col_name) = chan_data(:, c);
    end
    
    fprintf('  -> 已添加电极 [ %s ] 的数据\n', curr_chan_name);
end

% 4. 保存为CSV文件
try
    writetable(output_table, cfg.extract_output_filename);
    fprintf('====== 数据已成功提取并保存至 [ %s ] ======\n', cfg.extract_output_filename);
    
    % 显示前几行预览
    disp('数据表预览 (前5行):');
    disp(head(output_table, 5));
catch ME
    fprintf('错误: 保存文件失败。请确保文件未被Excel打开。\n错误信息: %s\n', ME.message);
end

fprintf('\n====== 所有分析流程已完成! ======\n');
% ====================================================================================
%%
% --- Part 9: 异常值检测参数 ---
cfg.outlier_roi_names = {'CZ', 'C1', 'C2', 'FZ', 'F1', 'F2', 'FCZ', 'FC1', 'FC2'}; % 您怀疑有问题的电极ROI
cfg.outlier_time_window = [250, 350];      % 您怀疑有问题的时间窗 (ms)
cfg.outlier_sd_threshold = 2;           % 定义异常值的标准差阈值 (通常在2到3之间)

fprintf('\nPart 9: 开始进行探索性的异常值检测...\n');
fprintf('!!! 警告: 请谨慎解读此部分结果，并遵守科研规范。!!!\n');

% --- 步骤1: 准备数据 ---
% 找到我们感兴趣的电极和时间点
outlier_ch_indices = find(ismember({chanlocs.labels}, cfg.outlier_roi_names));
outlier_time_indices = eeg_times >= cfg.outlier_time_window(1) & eeg_times <= cfg.outlier_time_window(2);

if isempty(outlier_ch_indices)
    error('错误: 在Part 9的配置中，找不到任何指定的ROI电极。');
end

fprintf('在ROI ( %s ) 和时间窗 [ %d, %d ] ms 内进行分析...\n', strjoin(cfg.outlier_roi_names, ', '), cfg.outlier_time_window(1), cfg.outlier_time_window(2));

% 提取数据 -> 平均电极 -> 平均时间
% 结果维度: 被试 x 条件
data_for_outlier = squeeze(mean(mean(all_subjects_4D_data(:, :, outlier_ch_indices, outlier_time_indices), 3), 4));

% --- 步骤2: 计算关键效应值 ---
% 我们以组内因素(Emotion)的差值作为衡量标准。
% 一个被试如果在这个差值上表现得与群体格格不入，就可能是异常值。
effect_diff = data_for_outlier(:, 1) - data_for_outlier(:, 2);

% --- 步骤3: 基于标准差识别异常值 ---
mean_diff = mean(effect_diff);
std_diff = std(effect_diff);
z_scores = (effect_diff - mean_diff) / std_diff; % 计算每个被试的Z-score

% 找到Z-score的绝对值超过阈值的被试
outlier_indices = find(abs(z_scores) > cfg.outlier_sd_threshold);

% --- 步骤4: 报告识别出的异常值 ---
if isempty(outlier_indices)
    fprintf('\n====== 在 %.1f 个标准差的阈值下，未发现异常值被试。 ======\n', cfg.outlier_sd_threshold);
else
    fprintf('\n====== 已识别出 %d 名可能的异常值被试 (阈值 > %.1f SD) ======\n', length(outlier_indices), cfg.outlier_sd_threshold);
    
    outlier_report = table();
    outlier_report.SubjectID = subject_info.SubjectID(outlier_indices);
    outlier_report.GroupName = subject_info.GroupName(outlier_indices);
    outlier_report.EffectDiff = effect_diff(outlier_indices);
    outlier_report.Z_Score = z_scores(outlier_indices);
    
    disp('被识别为异常值的被试信息:');
    disp(outlier_report);
    
    % --- 步骤5: 模拟剔除异常值后的ANOVA结果 ---
    fprintf('\n--- 正在模拟剔除以上被试后的ANOVA结果 ---\n');
    
    % 创建一个剔除了异常值的新数据表和信息表
    subjects_to_keep = true(total_subjects, 1);
    subjects_to_keep(outlier_indices) = false;
    
    % 同样在ROI和时间窗内进行检验
    data_for_anova_clean = data_for_outlier(subjects_to_keep, :);
    subject_info_clean = subject_info(subjects_to_keep, :);
    
    % 准备ANOVA数据表
    data_table_clean = array2table(data_for_anova_clean, 'VariableNames', sprintfc('Cond%d', 1:num_conditions));
    data_table_clean.Group = categorical(subject_info_clean.Group);

    % 执行混合设计ANOVA
    fprintf('剔除后剩余 %d 名被试。\n', height(data_table_clean));
    rm_model_clean = fitrm(data_table_clean, model_formula, 'WithinDesign', within_design);
    ranova_results_clean = ranova(rm_model_clean, 'WithinModel', within_factor_name);
    
    disp('模拟剔除异常值被试后的ANOVA结果表:');
    disp(ranova_results_clean);
end

fprintf('\n====== 探索性分析结束。 ======\n');

%% Part 10: [升级版] 批量导出多电极波形数据 (用于Origin/Excel作图)
% ====================================================================================
% 功能：批量导出指定列表中所有电极的平均波形数据。
% 格式：每个“组-电极”对生成一个CSV文件。
%       第1行：时间点 (Time)
%       第2行：条件1平均波幅
%       第3行：条件2平均波幅...
% ====================================================================================

% --- 参数配置 ---
% 【修改处】在这里填入所有想要导出的电极名称，用单引号括起来，逗号分隔
cfg.export_wave_channels = {'CZ', 'FZ', 'FCZ', 'CPZ'}; 

cfg.export_wave_filename_prefix = 'Waveform'; % 导出文件的前缀名

fprintf('\nPart 10: 正在批量导出 %d 个电极的波形数据...\n', length(cfg.export_wave_channels));

% --- 第一层循环：遍历您列表中的每一个电极 ---
for ch = 1:length(cfg.export_wave_channels)
    
    target_chan_name = cfg.export_wave_channels{ch};
    
    % 1. 找到当前电极在数据中的索引
    wave_ch_idx = find(strcmp({chanlocs.labels}, target_chan_name));
    
    if isempty(wave_ch_idx)
        fprintf('  警告: 数据中找不到电极 [ %s ]，跳过该电极。\n', target_chan_name);
        continue; % 如果找不到，跳过这个，继续下一个
    end
    
    fprintf('  正在处理电极: %s ...\n', target_chan_name);

    % --- 第二层循环：遍历每一个组 (Group) ---
    for g = 1:length(cfg.groups)
        group_name = cfg.groups(g).name;
        group_indices = find(subject_info.Group == g);
        
        if isempty(group_indices)
            continue;
        end
        
        % 提取该组数据: 
        % 原始维度: all_subjects_4D_data (被试 x 条件 x 通道 x 时间)
        % 取该组被试 -> 取特定通道 -> 对被试维(dim1)求平均 -> 压缩维度
        % 结果维度: 条件 x 时间
        group_avg_wave = squeeze(mean(all_subjects_4D_data(group_indices, :, wave_ch_idx, :), 1));
        
        % --- 构建输出矩阵 ---
        % 第1行放时间轴 (eeg_times)
        % 后续行放波形数据
        output_matrix = [eeg_times; group_avg_wave];
        
        % --- 生成文件名 ---
        % 格式: 前缀_组名_电极名.csv (例如: Waveform_DI_CZ.csv)
        fname = sprintf('%s_%s_%s.csv', cfg.export_wave_filename_prefix, group_name, target_chan_name);
        
        % --- 写入文件 ---
        try
            writematrix(output_matrix, fname);
            % fprintf('    -> 已保存: %s\n', fname);
        catch ME
            fprintf('    错误: 保存文件 %s 失败，请检查文件是否被打开。\n', fname);
        end
    end
end

fprintf('====== 批量波形导出完成! (共生成 %d 个文件) ======\n', length(cfg.export_wave_channels) * length(cfg.groups));


%% Part 11: [新功能] 便捷绘制各条件地形图 (Topography)
% ====================================================================================
% 功能：指定一个时间窗，快速画出所有组、所有条件的电压地形图
% ====================================================================================

% --- 参数配置 ---
cfg.topo_plot_window = [600, 800];    % 【在此处修改】您想看哪个时间段的地形图 (ms)
cfg.topo_plot_limits = [-5, 5];       % 【在此处修改】颜色范围 (uV)，保持一致以便对比

fprintf('\nPart 11: 正在生成时间窗 [ %d, %d ] ms 的地形图...\n', cfg.topo_plot_window(1), cfg.topo_plot_window(2));

% 1. 找到时间索引并取平均
t_idx = eeg_times >= cfg.topo_plot_window(1) & eeg_times <= cfg.topo_plot_window(2);
% data_topo_all: 被试 x 条件 x 通道 (压缩了时间维)
data_topo_all = squeeze(mean(all_subjects_4D_data(:, :, :, t_idx), 4));

% 2. 准备画布
figure('Name', '各条件地形图概览', 'Color', 'w', 'Position', [100, 100, 1000, 600]);
sgtitle(sprintf('Time: %d - %d ms', cfg.topo_plot_window(1), cfg.topo_plot_window(2)), 'FontSize', 16, 'FontWeight', 'bold');

plot_counter = 1;
num_groups = length(cfg.groups);
num_conds = num_conditions;

% 3. 双重循环：组 x 条件
for g = 1:num_groups
    group_indices = find(subject_info.Group == g);
    g_name = cfg.groups(g).name;
    
    for c = 1:num_conds
        % 计算该组、该条件的平均地形数据
        % mean over subjects
        topo_data = squeeze(mean(data_topo_all(group_indices, c, :), 1));
        
        % 绘图
        subplot(num_groups, num_conds, plot_counter);
        
        topoplot(topo_data, chanlocs, ...
            'maplimits', cfg.topo_plot_limits, ...
            'style', 'map', ...
            'electrodes', 'on', ...  % 显示电极点，方便定位
            'headrad', 'rim', ...
            'shading', 'interp');
            
        title(sprintf('%s - Cond %d', g_name, c), 'FontSize', 12);
        colorbar; % 显示色标
        
        plot_counter = plot_counter + 1;
    end
end

fprintf('====== 地形图绘制完成 ======\n');