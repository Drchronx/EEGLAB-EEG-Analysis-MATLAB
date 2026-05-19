% ====================================================================================
%  单因素 (Single Factor: One Group x Multiple Conditions) STFT 分析脚本
% ====================================================================================
clear all; clc; close all;

%% Part 1: 实验配置 (Experiment Configuration)
% ====================================================================================
fprintf('Part 1: 正在配置实验参数...\n');

cfg = struct();

% --- 数据路径配置 (单组) ---
cfg.data_dir = 'path_to_your_data'; % 修改为您的数据路径
cfg.conditions = {'1', '2', '3'};  % 修改为您的条件名称/Trigger号

% --- STFT 核心参数 ---
cfg.epoch_time = [-1, 2];          % 分段的时间窗口 (s)
cfg.rmbase_time = [-1000, 0];      % 预处理时的基线校正 (ms)
cfg.stft_freqs = 1:1:30;           % 频率轴 (Hz)
cfg.winsize = 0.400;               % 时间窗大小 (s)
cfg.baseline_window = [-0.8, -0.2];% STFT后的基线校正时间窗 (s)

% --- 通道选择 (内存优化) ---
% 如果内存不足，请改为 {'CZ', 'PZ'} 等具体电极
cfg.analyze_channels = 'all'; 

% --- 统计参数 ---
cfg.p_value_threshold = 0.05;


%% Part 2: 数据加载与 STFT 计算
% ====================================================================================
fprintf('\nPart 2: 正在加载数据并执行 STFT...\n');

P_BC = []; % 存储: 被试 x 条件 x 通道 x 频率 x 时间
subject_info = table();
subject_counter = 0;
chanlocs = [];
stft_times = [];
stft_freqs = cfg.stft_freqs;

file_list = dir(fullfile(cfg.data_dir, '*.set'));
fprintf('检测到 %d 个数据文件。\n', length(file_list));

for s = 1:length(file_list)
    fprintf('  处理被试 (%d/%d): %s ... ', s, length(file_list), file_list(s).name);
    
    try
        % 1. 加载数据
        EEG = pop_loadset('filename', file_list(s).name, 'filepath', cfg.data_dir);
        if isempty(chanlocs), chanlocs = EEG.chanlocs; end
        
        % 确定通道索引
        if ischar(cfg.analyze_channels) && strcmp(cfg.analyze_channels, 'all')
            ch_indices = 1:EEG.nbchan;
        else
            ch_indices = find(ismember({EEG.chanlocs.labels}, cfg.analyze_channels));
        end
        
        subj_cond_data = []; % 临时: 条件 x 通道 x 频率 x 时间
        is_subject_valid = true;
        
        for c = 1:length(cfg.conditions)
            % 2. 分段
            EEG_new = pop_epoch(EEG, {cfg.conditions{c}}, cfg.epoch_time, 'newname', 'STFT', 'epochinfo', 'yes', 'verbose', 'off');
            
            if EEG_new.trials == 0
                fprintf('[条件 %s 无试次] ', cfg.conditions{c}); 
                is_subject_valid = false; break;
            end
            
            % 3. 时域基线校正
            EEG_new = pop_rmbase(EEG_new, cfg.rmbase_time);
            
            % 准备参数
            Fs = EEG_new.srate;
            xtimes = EEG_new.times / 1000;
            if isempty(stft_times), stft_times = xtimes; end
            
            % 4. 逐通道 STFT
            for ch_idx = 1:length(ch_indices)
                real_ch = ch_indices(ch_idx);
                x = squeeze(EEG_new.data(real_ch, :, :)); % Time x Trials
                
                % 调用 sub_stft
                [~, P, ~, ~] = sub_stft(x, xtimes, xtimes, cfg.stft_freqs, Fs, cfg.winsize);
                % P: Freq x Time x Trials -> 平均 Trials
                pow_avg = squeeze(mean(P, 3)); 
                
                subj_cond_data(c, ch_idx, :, :) = pow_avg;
            end
        end
        
        if is_subject_valid
            subject_counter = subject_counter + 1;
            
            % 5. 频域基线校正 (dB转换: 10*log10(Data/Baseline)) 或 减法 (Data - Baseline)
            % 这里沿用您之前的【减法】逻辑
            bl_idx = find(stft_times >= cfg.baseline_window(1) & stft_times <= cfg.baseline_window(2));
            
            % 计算基线均值 (对时间维平均)
            baseline_vals = mean(subj_cond_data(:, :, :, bl_idx), 4);
            % 扩展维度以匹配
            baseline_expanded = repmat(baseline_vals, [1, 1, 1, size(subj_cond_data, 4)]);
            
            % 存入大矩阵
            P_BC(subject_counter, :, :, :, :) = subj_cond_data - baseline_expanded;
            
            subject_info.SubjectID(subject_counter) = subject_counter;
            subject_info.FileName{subject_counter} = file_list(s).name;
            fprintf('完成。\n');
        else
            fprintf('\n');
        end
        
    catch ME
        fprintf('  [错误: %s]\n', ME.message);
    end
end

fprintf('\n====== 数据计算完成! ======\n');
[n_sub, n_cond, n_chan, n_freq, n_time] = size(P_BC);
fprintf('数据维度: %d被试 x %d条件 x %d通道 x %d频率 x %d时间\n', n_sub, n_cond, n_chan, n_freq, n_time);


%% Part 3: 保存处理后的数据
% ====================================================================================
save_filename = 'SingleFactor_STFT_Data.mat';
fprintf('\n正在保存数据至: %s ...\n', save_filename);
save(save_filename, 'P_BC', 'subject_info', 'chanlocs', 'stft_times', 'stft_freqs', 'cfg', '-v7.3');
fprintf('保存成功! 下次可直接从 Part 4 开始。\n');


%% Part 4: 统计分析 (单因素重复测量 ANOVA)
% ====================================================================================
% --- 自动加载 ---
if ~exist('P_BC', 'var')
    load('SingleFactor_STFT_Data.mat');
end

cfg.roi_channel_name = 'CZ';    % 【修改】ROI电极
cfg.stat_freq_range  = [4 8];   % 【修改】统计频率窗
cfg.stat_time_range  = [0.5 1.0]; % 【修改】统计时间窗

fprintf('\nPart 4: 正在进行单因素统计 (ROI: %s)...\n', cfg.roi_channel_name);

% 1. 找电极索引
if ischar(cfg.analyze_channels) && strcmp(cfg.analyze_channels, 'all')
    roi_idx = find(strcmp({chanlocs.labels}, cfg.roi_channel_name));
else
    roi_idx = find(strcmp(cfg.analyze_channels, cfg.roi_channel_name));
end

% 2. 提取数据 (被试 x 条件 x 频率 x 时间)
roi_data = squeeze(P_BC(:, :, roi_idx, :, :));

% 3. 逐像素统计
p_values_grid = nan(n_freq, n_time);

% 准备 ANOVA 模型
% 因子名: Condition
factor_names = {'Condition'};
within_design = table(categorical(1:n_cond)', 'VariableNames', factor_names);
model_formula = sprintf('C1-C%d ~ 1', n_cond); % "~ 1" 表示只有截距(单组)

wb = waitbar(0, '正在计算统计...');
for f = 1:n_freq
    waitbar(f/n_freq, wb);
    for t = 1:n_time
        % 提取: 被试 x 条件
        dat = squeeze(roi_data(:, :, f, t));
        t_tab = array2table(dat, 'VariableNames', sprintfc('C%d', 1:n_cond));
        
        try
            rm = fitrm(t_tab, model_formula, 'WithinDesign', within_design);
            ra = ranova(rm);
            % 提取 'Condition' 的 P值
            % ranova 输出行名通常包含 '(Intercept):Condition' 或 'Condition'
            p_idx = find(contains(ra.Row, 'Condition'));
            if ~isempty(p_idx)
                p_values_grid(f, t) = ra.pValue(p_idx(1));
            end
        catch
        end
    end
end
close(wb);
fprintf('统计完成。\n');

%% Part 5: 可视化
% ====================================================================================
fprintf('\nPart 5: 绘图...\n');

% 1. 计算每个条件的平均图 (被试平均)
% roi_data: Subj x Cond x Freq x Time
cond_avgs = squeeze(mean(roi_data, 1)); % -> Cond x Freq x Time

% 确定色标
clim_range = [min(cond_avgs(:)), max(cond_avgs(:))];
% clim_range = [-2 2]; % 手动设置

figure('Color', 'w', 'Position', [100, 100, 1200, 600]);
sgtitle(sprintf('Single Factor Analysis - %s', cfg.roi_channel_name));

% 布局计算: 如果条件数<=3，排一行；否则排两行
num_plots = n_cond + 1; % 条件图 + 1个P值图
cols = min(num_plots, 4);
rows = ceil(num_plots / cols);

% 2. 绘制每个条件
for c = 1:n_cond
    subplot(rows, cols, c);
    imagesc(stft_times, stft_freqs, squeeze(cond_avgs(c,:,:)));
    axis xy;
    clim(clim_range);
    colorbar;
    title(sprintf('Cond: %s', cfg.conditions{c}));
    xlabel('Time (s)'); ylabel('Freq (Hz)');
    xline(0, '--k');
    
    % 画感兴趣的框
    rect_pos = [cfg.stat_time_range(1), cfg.stat_freq_range(1), ...
                diff(cfg.stat_time_range), diff(cfg.stat_freq_range)];
    hold on; rectangle('Position', rect_pos, 'EdgeColor', 'k'); hold off;
end

% 3. 绘制 P值图
subplot(rows, cols, n_cond + 1);
p_plot = p_values_grid;
p_plot(p_plot > cfg.p_value_threshold) = NaN; % 遮蔽不显著

imagesc(stft_times, stft_freqs, p_plot);
axis xy;
colormap(gca, flipud(hot));
clim([0 0.05]);
colorbar;
title('Main Effect of Condition (p<0.05)');
xlabel('Time (s)'); ylabel('Freq (Hz)');


%% Part 6: 地形图
% ====================================================================================
cfg.topo_t = [0.6 0.8]; 
cfg.topo_f = [4 7];
cfg.topo_lim = [-0.5 0.5];

fprintf('\nPart 6: 绘制地形图...\n');
t_idx = stft_times >= cfg.topo_t(1) & stft_times <= cfg.topo_t(2);
f_idx = stft_freqs >= cfg.topo_f(1) & stft_freqs <= cfg.topo_f(2);

% 提取: 被试 x 条件 x 通道
topo_raw = squeeze(mean(mean(P_BC(:, :, :, f_idx, t_idx), 4), 5));
% 平均被试 -> 条件 x 通道
topo_avg = squeeze(mean(topo_raw, 1));

figure('Color', 'w', 'Position', [100, 100, 800, 400]);
for c = 1:n_cond
    subplot(1, n_cond, c);
    topoplot(topo_avg(c, :), chanlocs, 'maplimits', cfg.topo_lim, 'style', 'map');
    title(sprintf('Cond: %s', cfg.conditions{c}));
    colorbar;
end

%% Part 7: 数据导出 (SPSS/Origin)
% ====================================================================================
cfg.extract_chans = {'CZ', 'PZ'}; % 要导出的电极
cfg.extract_fname = 'SingleFactor_Data_Export.csv';

fprintf('\nPart 7: 导出数据至 %s...\n', cfg.extract_fname);

out_table = subject_info;

for i = 1:length(cfg.extract_chans)
    ch_name = cfg.extract_chans{i};
    
    % 找索引
    if ischar(cfg.analyze_channels) && strcmp(cfg.analyze_channels, 'all')
         real_idx = find(strcmp({chanlocs.labels}, ch_name));
    else
         real_idx = find(strcmp(cfg.analyze_channels, ch_name));
    end
    
    % 提取数据: 被试 x 条件
    val = squeeze(mean(mean(P_BC(:, :, real_idx, f_idx, t_idx), 4), 5));
    
    for c = 1:n_cond
        col_name = sprintf('%s_%s', ch_name, cfg.conditions{c});
        out_table.(col_name) = val(:, c);
    end
end

writetable(out_table, cfg.extract_fname);
fprintf('完成。\n');