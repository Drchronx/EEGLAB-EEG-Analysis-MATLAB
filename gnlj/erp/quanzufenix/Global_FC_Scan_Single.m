% =========================================================================
%  全脑全局连接时频扫描 (Global Connectivity TF Scanner)
%  版本: 单因素双水平 (1 Factor, 2 Levels)
%  功能: 支持【组内设计 (Paired)】和【组间设计 (Independent)】切换
%  目的: 快速找出有显著差异的时间和频率窗口 (ROI)
% =========================================================================

clc; clear all; close all;

%% 1. 选择实验设计类型
design_type = questdlg('请选择您的实验设计类型:', ...
    'Design Type', ...
    '组内设计 (Paired, 同一批人)', ...
    '组间设计 (Independent, 两组人)', ...
    '组内设计 (Paired, 同一批人)');

if isempty(design_type), return; end

is_paired = strcmp(design_type, '组内设计 (Paired, 同一批人)');

if is_paired
    disp('>>> 模式：组内设计 (Paired t-test) <<<');
    titles_dlg = {'条件 1 (Cond 1)', '条件 2 (Cond 2)'};
else
    disp('>>> 模式：组间设计 (Independent t-test) <<<');
    titles_dlg = {'组别 1 (Group 1)', '组别 2 (Group 2)'};
end

%% 2. 参数设置
prompt = {'时间轴 (ms) [Start Step End]:', '频率轴 (Hz) [Start Step End]:', ...
          '基线范围 (ms) [Min Max]:', ...
          '条件/组别名称 1:', '条件/组别名称 2:'};
def = {'[-1000 4 1996]', '[1 1 30]', '[-800 -200]', 'A', 'B'};
ans_set = inputdlg(prompt, 'Scan Parameters', 1, def);
if isempty(ans_set), return; end

time_axis = str2num(ans_set{1});
f_axis    = str2num(ans_set{2});
baseline  = str2num(ans_set{3});
Name_1    = ans_set{4};
Name_2    = ans_set{5};

%% 3. 载入数据
% 容器: [Freq, Time, Subj]
Global_Data1 = []; 
Global_Data2 = [];
nSubj1 = 0; nSubj2 = 0;

% --- 载入第一个数据 (Cond1 或 Group1) ---
fprintf('正在载入: %s ...\n', titles_dlg{1});
[fn1, pn1] = uigetfile('*.mat', ['Select: ' titles_dlg{1}]);
if isequal(fn1,0), return; end
tmp = load(fullfile(pn1, fn1)); vars = fieldnames(tmp);
data_5d = tmp.(vars{1}); % [F, T, Ch, Ch, Subj]

[nF, nT, ~, ~, nS] = size(data_5d);
nSubj1 = nS;

% 全局平均 + 基线校正
g_avg = squeeze(mean(mean(data_5d, 3, 'omitnan'), 4, 'omitnan')); 
base_idx = time_axis >= baseline(1) & time_axis <= baseline(2);
base_val = mean(g_avg(:, base_idx, :), 2);
Global_Data1 = g_avg - repmat(base_val, [1, nT, 1]);
clear data_5d g_avg

% --- 载入第二个数据 (Cond2 或 Group2) ---
fprintf('正在载入: %s ...\n', titles_dlg{2});
[fn2, pn2] = uigetfile('*.mat', ['Select: ' titles_dlg{2}]);
if isequal(fn2,0), return; end
tmp = load(fullfile(pn2, fn2)); vars = fieldnames(tmp);
data_5d = tmp.(vars{1});

[nF2, nT2, ~, ~, nS2] = size(data_5d);
nSubj2 = nS2;

if nF~=nF2 || nT~=nT2
    error('两个文件的时间点或频率点数量不一致！');
end

% 组内设计必须检查人数一致
if is_paired && (nSubj1 ~= nSubj2)
    error('错误：组内设计要求两个条件的文件被试数量必须完全相同！(N1=%d, N2=%d)', nSubj1, nSubj2);
end

% 全局平均 + 基线校正
g_avg = squeeze(mean(mean(data_5d, 3, 'omitnan'), 4, 'omitnan')); 
base_val = mean(g_avg(:, base_idx, :), 2);
Global_Data2 = g_avg - repmat(base_val, [1, nT, 1]);
clear data_5d g_avg

fprintf('数据就绪。比较: %s vs %s\n', Name_1, Name_2);

%% 4. 快速扫描 (逐点 t-test)
fprintf('正在扫描显著性...\n');
P_Map = ones(nF, nT);
h = waitbar(0, 'Scanning...');

for f = 1:nF
    waitbar(f/nF, h);
    for t = 1:nT
        % 提取当前点数据
        Y1 = squeeze(Global_Data1(f, t, :)); % [Subj1 x 1]
        Y2 = squeeze(Global_Data2(f, t, :)); % [Subj2 x 1]
        
        if is_paired
            % 组内: 配对 t 检验
            [~, p] = ttest(Y1, Y2);
        else
            % 组间: 独立样本 t 检验
            [~, p] = ttest2(Y1, Y2);
        end
        P_Map(f, t) = p;
    end
end
close(h);

%% 5. 绘图 (FDR 校正)
% FDR 校正
[~, Mask] = fdr(P_Map, 0.05);

figure('Color','w', 'Position', [100, 100, 600, 500], 'Name', 'Single Factor Scout');

% 绘制 P 值 (Masked)
p_plot = P_Map;
p_plot(Mask == 0) = NaN; % 只显示显著区域

if all(isnan(p_plot(:)))
    text(mean(time_axis), mean(f_axis), 'No Significant Area', 'Horiz','center', 'FontSize',14);
    xlim([min(time_axis) max(time_axis)]); ylim([min(f_axis) max(f_axis)]);
else
    imagesc(time_axis, f_axis, p_plot); 
    axis xy; colorbar;
    caxis([0 0.05]); % 锁定 0 到 0.05
    colormap(gca, flipud(parula));
end

if is_paired
    design_label = 'Paired';
else
    design_label = 'Independent';
end
title_str = sprintf('Significant Area: %s vs %s (%s)', Name_1, Name_2, design_label);
title(title_str, 'FontSize', 12, 'FontWeight','bold');
xlabel('Time (ms)'); ylabel('Freq (Hz)');

msgbox({'扫描完成！'; ''; '请根据图中的色块确定 ROI (Time & Freq)'; '然后填入全脑分析代码。'});
