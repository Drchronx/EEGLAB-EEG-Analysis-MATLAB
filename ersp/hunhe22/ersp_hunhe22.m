% ====================================================================================
%  混合设计 STFT (短时傅里叶) 分析脚本 - 适配 Group x Condition
% ====================================================================================
% 注意：请确保 sub_stft 函数在您的路径中，或者使用脚本末尾提供的替代函数。

clear all; clc; close all;

%% Part 1: 实验配置 (Experiment Configuration)
% ====================================================================================
fprintf('Part 1: 正在配置 STFT 分析参数...\n');

cfg = struct();

% --- 实验组设计与路径配置 ---
cfg.groups(1).name = 'cz';
cfg.groups(1).data_dir = 'path_to_your_data';
cfg.groups(1).conditions = {'1', '2'};

cfg.groups(2).name = 'gj';
cfg.groups(2).data_dir = 'path_to_your_data';
cfg.groups(2).conditions = {'5', '6'};

% --- 组内因素定义 ---
cfg.within_factors.names = {'rw'};
cfg.within_factors.levels = [1, 2];

% --- STFT 核心参数 (按您提供的参数设置) ---
cfg.epoch_time = [-1, 2];          % 分段的时间窗口 (s)
cfg.rmbase_time = [-1000, 0];      % 预处理时的基线校正 (ms)
cfg.stft_freqs = 1:1:30;           % 频率轴 (Hz)
cfg.winsize = 0.400;               % 时间窗大小 (s)
cfg.baseline_window = [-0.8, -0.2];% STFT后的基线校正时间窗 (s)

% --- 通道选择 (内存优化) ---
% !!! 警告: STFT 全通道数据非常占内存。
% 建议先用 {'CZ', 'PZ'} 测试跑通，内存够大再换 'all'
cfg.analyze_channels = 'all'; % 或者 {'CZ', 'PZ', 'P1', 'P2'}

% --- 统计参数 ---
cfg.p_value_threshold = 0.05;

%% Part 2: 数据加载与 STFT 计算
% ====================================================================================
fprintf('\nPart 2: 正在加载数据并执行 STFT (Fast Fourier)...\n');

% 初始化存储矩阵
% 维度预估: 我们稍后根据第一个被试的数据来初始化大矩阵以节省内存重新分配的时间
P_BC = []; % 最终存储基线校正后的数据
subject_info = table();
subject_counter = 0;
chanlocs = [];
stft_times = [];
stft_freqs = cfg.stft_freqs;

for g = 1:length(cfg.groups)
    group_name = cfg.groups(g).name;
    data_dir = cfg.groups(g).data_dir;
    conditions = cfg.groups(g).conditions;
    file_list = dir(fullfile(data_dir, '*.set'));
    
    fprintf('--- 正在处理组: %s (%d 个文件) ---\n', group_name, length(file_list));
    
    for s = 1:length(file_list)
        fprintf('  处理被试: %s ... ', file_list(s).name);
        
        try
            % 1. 加载数据
            EEG = pop_loadset('filename', file_list(s).name, 'filepath', data_dir);
            EEG = eeg_checkset(EEG);
            
            if isempty(chanlocs), chanlocs = EEG.chanlocs; end
            
            % 确定要处理的通道
            if ischar(cfg.analyze_channels) && strcmp(cfg.analyze_channels, 'all')
                ch_indices = 1:EEG.nbchan;
            else
                ch_indices = find(ismember({EEG.chanlocs.labels}, cfg.analyze_channels));
            end
            
            is_subject_valid = true;
            subj_cond_data = []; % 临时存储该被试所有条件的数据
            
            for c = 1:length(conditions)
                % 2. 重新分段 (按您的代码)
                EEG_new = pop_epoch(EEG, {conditions{c}}, cfg.epoch_time, 'newname', 'STFT_Epochs', 'epochinfo', 'yes', 'verbose', 'off');
                EEG_new = eeg_checkset(EEG_new);
                
                if EEG_new.trials == 0
                    fprintf('[无试次] '); is_subject_valid = false; break;
                end
                
                % 3. 时域基线校正 (按您的代码)
                EEG_new = pop_rmbase(EEG_new, cfg.rmbase_time);
                EEG_new = eeg_checkset(EEG_new);
                
                % 准备 STFT 参数
                Fs = EEG_new.srate;
                xtimes = EEG_new.times / 1000; % ms -> s
                t_points = xtimes; % 在每个时间点计算
                
                % 4. 逐通道 STFT
                cond_data_temp = [];
                for ch_idx = 1:length(ch_indices)
                    real_ch = ch_indices(ch_idx);
                    x = squeeze(EEG_new.data(real_ch, :, :)); % Time x Trials
                    
                    % --- 调用 sub_stft ---
                    % 确保 sub_stft 在路径中，否则使用后面的替代函数
                    try
                        [~, P, ~, ~] = sub_stft(x, xtimes, t_points, cfg.stft_freqs, Fs, cfg.winsize);
                        % P 的维度通常是: 频率 x 时间 x 试次 (需确认 sub_stft 输出)
                        % 您的代码中: P_data(..., squeeze(mean(P,3)))
                        % 假设 P 是 Freq x Time x Trials
                        
                        pow_avg = squeeze(mean(P, 3)); % 平均试次 -> Freq x Time
                        
                    catch ME
                        if strcmp(ME.identifier, 'MATLAB:UndefinedFunction')
                             error('找不到 sub_stft 函数。请将该函数放入文件夹，或使用本脚本末尾提供的替代方案。');
                        else
                             rethrow(ME);
                        end
                    end
                    
                    cond_data_temp(ch_idx, :, :) = pow_avg; % Chan x Freq x Time
                end
                
                subj_cond_data(c, :, :, :) = cond_data_temp;
                
                if isempty(stft_times)
                    % 尝试修正时间轴 (sub_stft可能会稍微改变时间点，这取决于具体实现)
                    % 这里暂且沿用 EEG_new.times，如果 sub_stft 返回了 U 或 F，请在此处更新
                    stft_times = xtimes; 
                end
            end
            
            if is_subject_valid
                subject_counter = subject_counter + 1;
                
                % --- 立即进行基线校正 (节省内存策略) ---
                % 您的逻辑: P_BC = temp_data - mean(baseline)
                % 数据维度: Cond x Chan x Freq x Time
                
                % 1. 找到基线索引
                bl_idx = find(stft_times >= cfg.baseline_window(1) & stft_times <= cfg.baseline_window(2));
                
                % 2. 执行校正 (减法)
                % 计算基线均值: 对时间维(dim 4)的基线段求平均
                baseline_vals = mean(subj_cond_data(:, :, :, bl_idx), 4); 
                
                % 扩展基线矩阵以匹配数据维度进行减法
                baseline_expanded = repmat(baseline_vals, [1, 1, 1, size(subj_cond_data, 4)]);
                
                % 存储结果 (P_BC)
                % 最终大矩阵维度: 被试 x 条件 x 通道 x 频率 x 时间
                P_BC(subject_counter, :, :, :, :) = subj_cond_data - baseline_expanded;
                
                % 记录信息
                subject_info.SubjectID(subject_counter) = subject_counter;
                subject_info.Group(subject_counter) = g;
                subject_info.GroupName{subject_counter} = cfg.groups(g).name;
                
                fprintf('完成 (已校正)。\n');
            else
                fprintf('\n');
            end
            
        catch ME
            fprintf('  [错误: %s]\n', ME.message);
        end
    end
end

fprintf('\n====== 所有数据处理完成! ======\n');
[n_sub, n_cond, n_chan, n_freq, n_time] = size(P_BC);
fprintf('最终数据维度 (P_BC): %d被试 x %d条件 x %d通道 x %d频率 x %d时间\n', n_sub, n_cond, n_chan, n_freq, n_time);
%% --- [新增] 保存计算好的 ERSP 数据 ---
% 数据处理通常很耗时，保存下来可以避免重复计算
save_filename = 'STFT_Processed_Data.mat';
fprintf('\n正在将处理好的数据保存至: %s ...\n', save_filename);

% 保存所有后续分析需要的变量
% '-v7.3' 是为了支持超过2GB的大文件 (ERSP数据通常很大)
save(save_filename, ...
    'P_BC', ...           % 最核心的数据 (已基线校正)
    'subject_info', ...   % 被试分组信息
    'chanlocs', ...       % 电极定位信息
    'stft_times', ...     % 时间轴
    'stft_freqs', ...     % 频率轴
    'cfg', ...            % 配置信息
    '-v7.3'); 

fprintf('====== 数据保存成功! 下次可直接从 Part 4 开始运行 ======\n');
%% 
%%



%% Part 4: 统计分析 (ROI 混合设计 ANOVA)
% ====================================================================================

% --- [新增] 智能加载数据 ---
% 检查工作区是否有 P_BC 变量，如果没有，说明是新开的 MATLAB，需要加载
if ~exist('P_BC', 'var')
    fprintf('检测到工作区无数据，正在加载上次保存的 "STFT_Processed_Data.mat" ...\n');
    try
        load('STFT_Processed_Data.mat');
        fprintf('数据加载成功！\n');
    catch
        error('错误: 找不到文件 "STFT_Processed_Data.mat"。\n请先运行 Part 1-3 生成并保存数据，或者检查文件名是否正确。');
    end
else
    fprintf('检测到数据已在工作区，直接进行统计分析...\n');
end

% --- 下面接您原来的参数设置 ---
cfg.roi_channel_name = 'CZ';    % 选择分析电极
% ... (后续代码保持不变)
% ====================================================================================
%% Part 4: 统计分析 (ROI 混合设计 ANOVA)
% ====================================================================================
cfg.roi_channel_name = 'CZ';    % 选择分析电极
% cfg.roi_channel_name = '19';  % 如果您的电极名是数字字符串
cfg.stat_freq_range  = [4 8];   % 统计频率窗 (Hz)
cfg.stat_time_range  = [0.5 1.0]; % 统计时间窗 (s) - 请根据您的需求修改

fprintf('\nPart 4: 正在对通道 [ %s ] 进行混合设计 ANOVA...\n', cfg.roi_channel_name);

% 1. 找到电极索引
if ischar(cfg.analyze_channels) && strcmp(cfg.analyze_channels, 'all')
    roi_idx = find(strcmp({chanlocs.labels}, cfg.roi_channel_name));
    % 如果找不到且是数字，尝试转换
    if isempty(roi_idx) && ~isnan(str2double(cfg.roi_channel_name))
         roi_idx = str2double(cfg.roi_channel_name);
    end
else
    roi_idx = find(strcmp(cfg.analyze_channels, cfg.roi_channel_name));
end

if isempty(roi_idx), error('未找到ROI电极，请检查名称。'); end

% 2. 提取数据 (被试 x 条件 x 频率 x 时间)
roi_data = squeeze(P_BC(:, :, roi_idx, :, :));

% 3. 逐像素统计 (Group x Condition 交互作用)
p_values_grid = nan(n_freq, n_time);

within_factor = cfg.within_factors.names{1};
within_design = table(categorical(cfg.within_factors.levels'), 'VariableNames', {within_factor});
model_formula = sprintf('C1-C%d ~ Group', n_cond);

fprintf('  正在计算逐点 ANOVA (这可能需要几分钟)...\n');
wb = waitbar(0, '计算中...');
for f = 1:n_freq
    waitbar(f/n_freq, wb);
    for t = 1:n_time
        % 提取: 被试 x 条件
        dat = squeeze(roi_data(:, :, f, t));
        t_tab = array2table(dat, 'VariableNames', sprintfc('C%d', 1:n_cond));
        t_tab.Group = categorical(subject_info.Group);
        
        try
            rm = fitrm(t_tab, model_formula, 'WithinDesign', within_design);
            ra = ranova(rm, 'WithinModel', within_factor);
            
            % 提取交互作用 P值
            inter_term = ['Group:' within_factor];
            if any(strcmp(ra.Row, inter_term))
                p_values_grid(f, t) = ra{inter_term, 'pValue'};
            end
        catch
        end
    end
end
close(wb);
fprintf('  统计完成。\n');


%% --- [新增] 保存中间结果 ---
% 这里的 '-v7.3' 是为了支持大文件保存（ERSP数据通常很大）
save_filename = 'ERSP_Analysis_Results.mat';
fprintf('\n正在保存处理后的数据和统计结果至: %s ...\n', save_filename);

% 保存画图和后续导出所需的所有关键变量
save(save_filename, ...
    'P_BC', ...           % 完整的基线校正后数据 (大矩阵)
    'subject_info', ...   % 被试信息
    'chanlocs', ...       % 电极定位
    'stft_times', ...     % 时间轴
    'stft_freqs', ...     % 频率轴
    'p_values_grid', ...  % 统计结果 (P值矩阵)
    'roi_data', ...       % 提取出的 ROI 数据 (方便直接画图)
    'cfg', ...            % 配置参数
    '-v7.3');             % 启用大文件支持

fprintf('====== 数据保存成功! 您下次可以直接加载该文件进行画图 ======\n');
%% Part 5: 可视化 (2x2 STFT 图 + 显著性框)
%% Part 5: 可视化 (2x2 STFT 图 + 显著性框)
% ====================================================================================
% ====================================================================================
fprintf('\nPart 5: 正在绘图...\n');

% 1. 计算组平均
g1_idx = find(subject_info.Group == 1);
g2_idx = find(subject_info.Group == 2);

% 数据: 条件 x 频率 x 时间 (对被试平均)
avg_g1 = squeeze(mean(roi_data(g1_idx, :, :, :), 1));
avg_g2 = squeeze(mean(roi_data(g2_idx, :, :, :), 1));

% 准备绘图数据
plot_data = {squeeze(avg_g1(1,:,:)), squeeze(avg_g2(1,:,:)), ...
             squeeze(avg_g1(2,:,:)), squeeze(avg_g2(2,:,:))};
titles = {sprintf('%s - C1', cfg.groups(1).name), sprintf('%s - C1', cfg.groups(2).name), ...
          sprintf('%s - C2', cfg.groups(1).name), sprintf('%s - C2', cfg.groups(2).name)};

% 2. 确定色标范围 (根据您的代码，可以是 min max)
all_vals = [avg_g1(:); avg_g2(:)];
clim_range = [min(all_vals), max(all_vals)];
% 或者手动: clim_range = [-2 2];

figure('Color', 'w', 'Position', [100, 100, 1000, 800]);
% 注意: sgtitle 是 R2018b 引入的。如果报错 'sgtitle' 无法识别，请改用 suptitle (需下载) 或手动用 annotation
sgtitle(sprintf('STFT Power (Baseline Corrected) - Ch %s', cfg.roi_channel_name));

% 3. 绘制 4 个条件
for i = 1:4
    subplot(3, 2, i);
    imagesc(stft_times, stft_freqs, plot_data{i});
    axis xy; % 频率轴向上
    
    % --- 修改处 1: 将 clim 替换为 caxis ---
    caxis(clim_range); 
    % ------------------------------------
    
    colorbar;
    title(titles{i});
    xlabel('Time (s)'); ylabel('Freq (Hz)');
    
    % 画0时刻线
    xline(0, '--k');
    
    % 画感兴趣的框 (ROI Box)
    % 矩形: [x, y, w, h]
    rect_pos = [cfg.stat_time_range(1), cfg.stat_freq_range(1), ...
                diff(cfg.stat_time_range), diff(cfg.stat_freq_range)];
    hold on;
    rectangle('Position', rect_pos, 'EdgeColor', 'k', 'LineWidth', 2, 'LineStyle', '-');
    hold off;
end

% 4. 绘制 P值图
subplot(3, 2, [5 6]);
p_plot = p_values_grid;
p_plot(p_plot > cfg.p_value_threshold) = NaN; % 不显著设为透明(白色)

imagesc(stft_times, stft_freqs, p_plot);
axis xy;
colormap(gca, flipud(hot)); % 显著的颜色深

% --- 修改处 2: 将 clim 替换为 caxis ---
caxis([0 0.05]); 
% ------------------------------------

colorbar;
title(sprintf('Interaction Effect (p < %.2f)', cfg.p_value_threshold));
xlabel('Time (s)'); ylabel('Freq (Hz)');
xline(0, '--k');
%% Part 11: 地形图 (特定频段和时间窗)
% ====================================================================================
% 按您的要求: 0.6-0.8s, 4-7Hz (theta)
cfg.topo_t = [0.6 0.8]; 
cfg.topo_f = [4 7];
cfg.topo_lim = [-0.5 0.5]; % 根据您的数据调整

fprintf('\nPart 11: 绘制地形图 [%.1f-%.1f s, %d-%d Hz]...\n', ...
    cfg.topo_t(1), cfg.topo_t(2), cfg.topo_f(1), cfg.topo_f(2));

% 1. 索引
t_idx = stft_times >= cfg.topo_t(1) & stft_times <= cfg.topo_t(2);
f_idx = stft_freqs >= cfg.topo_f(1) & stft_freqs <= cfg.topo_f(2);

% 2. 提取并平均 (被试 x 条件 x 通道)
topo_raw = squeeze(mean(mean(P_BC(:, :, :, f_idx, t_idx), 4), 5));

figure('Color', 'w', 'Position', [100, 100, 800, 600]);
sgtitle('Topography Analysis');

cnt = 1;
for g = 1:length(cfg.groups)
    g_idx = find(subject_info.Group == g);
    for c = 1:length(conditions)
        % 组平均
        topo_dat = squeeze(mean(topo_raw(g_idx, c, :), 1));
        
        subplot(2, 2, cnt);
        topoplot(topo_dat, chanlocs, 'maplimits', cfg.topo_lim, 'style', 'map', 'electrodes', 'off');
        title(sprintf('%s - C%d', cfg.groups(g).name, c));
        colorbar;
        cnt = cnt + 1;
    end
end

fprintf('====== 分析结束 ======\n');


%% 附录: sub_stft 替代函数
% 如果您没有 sub_stft.m，请取消以下注释并将其保存为 sub_stft.m 或放在脚本最后
% (注意: MATLAB 脚本中函数必须放在最后)

function [S, P, F, T] = sub_stft(x, xtimes, t_req, f_req, Fs, winsize_sec)
    % 一个基于 spectrogram 的简易实现，模拟 sub_stft 的行为
    % 输入:
    %   x: 信号 (Time x Trials) 或 (Time x 1)
    %   winsize_sec: 窗口大小 (秒)
    
    % 参数转换
    window = round(winsize_sec * Fs);
    noverlap = round(window * 0.9); % 90% 重叠，平滑度高
    nfft = max(256, 2^nextpow2(window)); 
    
    % 如果 x 是多试次的，我们需要循环或重塑。
    % 这里假设 sub_stft 输出 P 为: Freq x Time x Trials
    
    [n_samples, n_trials] = size(x);
    P = [];
    
    for i = 1:n_trials
        [s, f, t, p] = spectrogram(x(:,i), window, noverlap, f_req, Fs);
        % s: Freq x Time
        if i == 1
            F = f;
            % 调整时间轴以匹配原始 xtimes (spectrogram的时间是从0开始的相对时间)
            % 这里做一个简单的映射，这只是近似
            T = t + xtimes(1); 
        end
        P(:, :, i) = abs(s).^2; % Power
    end
    
    S = []; % Complex spectrum (optional)
end