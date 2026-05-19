% =========================================================================
%  静息态全脑全局连接时频扫描 (Resting Global Connectivity TF Scanner)
%  版本: 单因素双水平 (Single Factor, 2 Levels)
%  功能: 
%    1. 支持【组内 (Paired)】和【组间 (Independent)】一键切换
%    2. 计算全脑全局平均连接 (Global Mean FC)
%    3. 快速扫描显著时频窗口 (ROI)
% =========================================================================

clc; clear all; close all;

%% 1. 选择实验设计类型
design_type = questdlg('请选择实验设计类型:', ...
    'Design Type', ...
    '组内 (Paired, 同一批人)', ...
    '组间 (Independent, 两组人)', ...
    '组内 (Paired, 同一批人)');

if isempty(design_type), return; end
is_paired = strcmp(design_type, '组内 (Paired, 同一批人)');

if is_paired
    disp('>>> 模式：组内设计 (Paired t-test) <<<');
    titles_dlg = {'条件 1 (Cond 1)', '条件 2 (Cond 2)'};
else
    disp('>>> 模式：组间设计 (Independent t-test) <<<');
    titles_dlg = {'组别 1 (Group 1)', '组别 2 (Group 2)'};
end

%% 2. 参数设置
prompt = {'时间轴 (ms/Seg) [Start Step End]:', '频率轴 (Hz) [Start Step End]:', ...
          '基线范围 (ms) [Min Max] (无基线填0):', ...
          '名称 1:', '名称 2:'};
% 静息态通常是分段的，假设时间轴是分段编号或时间
def = {'[0 2 1000]', '[1 1 30]', '0', 'A', 'B'};
ans_set = inputdlg(prompt, 'Scan Params', 1, def);
if isempty(ans_set), return; end

time_axis = str2num(ans_set{1});
f_axis    = str2num(ans_set{2});
baseline  = str2num(ans_set{3});
Name_1    = ans_set{4};
Name_2    = ans_set{5};

%% 3. 载入数据 (并计算全局平均)
% 容器: [Freq, Time, Subj]
Global_Data1 = []; 
Global_Data2 = [];
nSubj1 = 0; nSubj2 = 0;

% --- 载入数据 1 ---
fprintf('正在载入: %s ...\n', titles_dlg{1});
[fn1, pn1] = uigetfile('*.mat', ['Select: ' titles_dlg{1}]);
if isequal(fn1,0), return; end
tmp = load(fullfile(pn1, fn1)); vars = fieldnames(tmp);
data_5d = tmp.(vars{1}); % [F, T, Ch, Ch, Subj]

[nF, nT, ~, ~, nS] = size(data_5d);
nSubj1 = nS;

% !!! 核心：对通道维(3,4)求平均 !!!
g_avg = squeeze(mean(mean(data_5d, 3, 'omitnan'), 4, 'omitnan')); 

% 基线校正 (如果有)
if length(baseline) > 1
    base_idx = time_axis >= baseline(1) & time_axis <= baseline(2);
    base_val = mean(g_avg(:, base_idx, :), 2);
    g_avg = g_avg - repmat(base_val, [1, nT, 1]);
end
Global_Data1 = g_avg;
clear data_5d g_avg

% --- 载入数据 2 ---
fprintf('正在载入: %s ...\n', titles_dlg{2});
[fn2, pn2] = uigetfile('*.mat', ['Select: ' titles_dlg{2}]);
if isequal(fn2,0), return; end
tmp = load(fullfile(pn2, fn2)); vars = fieldnames(tmp);
data_5d = tmp.(vars{1});

[nF2, nT2, ~, ~, nS2] = size(data_5d);
nSubj2 = nS2;

if nF~=nF2 || nT~=nT2, error('两个文件的时间/频率维度不一致！'); end
if is_paired && (nSubj1 ~= nSubj2)
    error('错误：组内设计要求被试数量必须一致！(N1=%d, N2=%d)', nSubj1, nSubj2);
end

g_avg = squeeze(mean(mean(data_5d, 3, 'omitnan'), 4, 'omitnan')); 
if length(baseline) > 1
    base_val = mean(g_avg(:, base_idx, :), 2);
    g_avg = g_avg - repmat(base_val, [1, nT, 1]);
end
Global_Data2 = g_avg;
clear data_5d g_avg

%% 4. 快速扫描 (逐点 t-test)
fprintf('正在扫描显著性...\n');
P_Map = ones(nF, nT);
h = waitbar(0, 'Scanning...');

for f = 1:nF
    waitbar(f/nF, h);
    for t = 1:nT
        Y1 = squeeze(Global_Data1(f, t, :)); 
        Y2 = squeeze(Global_Data2(f, t, :)); 
        
        if is_paired
            [~, p] = ttest(Y1, Y2);
        else
            [~, p] = ttest2(Y1, Y2);
        end
        P_Map(f, t) = p;
    end
end
close(h);

%% 5. 绘图 (FDR 校正)
[~, Mask] = fdr(P_Map, 0.05);

figure('Color','w', 'Position', [100, 100, 600, 500], 'Name', 'Single Factor Scout');

% 绘制 P 值
p_plot = P_Map;
p_plot(Mask == 0) = NaN; % 只显示显著区域

if all(isnan(p_plot(:)))
    text(mean(time_axis), mean(f_axis), 'No Significant Area', 'Horiz','center', 'FontSize',14);
    xlim([min(time_axis) max(time_axis)]); ylim([min(f_axis) max(f_axis)]);
else
    imagesc(time_axis, f_axis, p_plot); 
    axis xy; colorbar;
    caxis([0 0.05]); 
    colormap(gca, flipud(parula));
end

if is_paired
    design_label = 'Paired';
else
    design_label = 'Independent';
end
title_str = sprintf('Significant ROI: %s vs %s (%s)', Name_1, Name_2, design_label);
title(title_str, 'FontSize', 12, 'FontWeight','bold');
xlabel('Time/Segment'); ylabel('Freq (Hz)');

msgbox({'扫描完成！'; ''; '请根据图中的色块确定显著 ROI (Time & Freq)'; '然后运行全脑分析代码。'});
