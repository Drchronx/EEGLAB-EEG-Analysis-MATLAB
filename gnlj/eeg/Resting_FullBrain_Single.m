% =========================================================================
%  静息态全脑连接组学分析 (Resting Whole-brain Connectome Analysis)
%  版本: 单因素双水平 (1 Factor, 2 Levels)
%  功能: 
%    1. 支持【组内】和【组间】切换
%    2. 提取指定时频 ROI 平均值
%    3. 全脑逐边统计 (Mass-univariate) + FDR 校正
%    4. 导出 BrainNet Viewer (.node / .edge)
% =========================================================================

clc; clear all; close all;

%% 1. 模式选择与参数
design_type = questdlg('请选择实验设计类型:', 'Design Type', ...
    '组内 (Paired)', '组间 (Independent)', '组内 (Paired)');
if isempty(design_type), return; end
is_paired = strcmp(design_type, '组内 (Paired)');

prompt = { ...
    '输入显著频率范围 (Hz) [Min Max]:', ...
    '输入显著时间范围 (ms/Seg) [Min Max]:', ...
    '时间轴 (ms/Seg) [Start Step End]:', ...
    '频率轴 (Hz) [Start Step End]:', ...
    '名称 1 (e.g. Pre):', '名称 2 (e.g. Post):' ...
};
def = {'[4 7]', '[0 1000]', '[0 2 2000]', '[1 1 30]', 'Cond1', 'Cond2'};
ans_set = inputdlg(prompt, 'FullBrain Params', 1, def);
if isempty(ans_set), return; end

f_roi = str2num(ans_set{1});
t_roi = str2num(ans_set{2});
time_axis = str2num(ans_set{3});
f_axis = str2num(ans_set{4});
Name_1 = ans_set{5}; Name_2 = ans_set{6};

% 载入电极
disp('>>> [重要] 请选择一个 *预处理后* (已剔除电极) 的 .set 文件 <<<');
[fn, pn] = uigetfile('*.set', 'Load .set file');
if isequal(fn,0), return; end
EEG = pop_loadset('filename', fn, 'filepath', pn);
chanlocs = EEG.chanlocs;
nChans = length(chanlocs);
nPairs = nChans * (nChans - 1) / 2; 

fprintf('检测到 %d 个电极，共 %d 条连接。\n', nChans, nPairs);

%% 2. 载入数据 (提取 ROI 均值)
pairs_idx = nchoosek(1:nChans, 2); 
Data_1 = []; Data_2 = [];
nSubj1 = 0; nSubj2 = 0;

if is_paired
    file_titles = {'Condition 1', 'Condition 2'};
    stat_label = 'Paired';
else
    file_titles = {'Group 1', 'Group 2'};
    stat_label = 'Independent';
end

% --- 读取文件 1 ---
fprintf('正在读取: %s ...\n', file_titles{1});
[fn, pn] = uigetfile('*.mat', ['Select: ' file_titles{1}]);
if isequal(fn,0), return; end
tmp = load(fullfile(pn, fn)); vars = fieldnames(tmp);
data_5d = tmp.(vars{1}); % [Freq, Time, Ch, Ch, Subj]
[~, ~, ~, ~, nSubj1] = size(data_5d);

% 提取 ROI 并拉直
t_idx = time_axis >= t_roi(1) & time_axis <= t_roi(2);
f_idx = f_axis >= f_roi(1) & f_axis <= f_roi(2);
temp_avg = squeeze(mean(mean(data_5d(f_idx, t_idx, :, :, :), 1, 'omitnan'), 2, 'omitnan'));

pairs_data = zeros(nSubj1, nPairs);
for s = 1:nSubj1
    curr = temp_avg(:,:,s);
    for p = 1:nPairs
        pairs_data(s, p) = curr(pairs_idx(p,1), pairs_idx(p,2));
    end
end
Data_1 = pairs_data; 
clear data_5d temp_avg pairs_data

% --- 读取文件 2 ---
fprintf('正在读取: %s ...\n', file_titles{2});
[fn, pn] = uigetfile('*.mat', ['Select: ' file_titles{2}]);
if isequal(fn,0), return; end
tmp = load(fullfile(pn, fn)); vars = fieldnames(tmp);
data_5d = tmp.(vars{1});
[~, ~, ~, ~, nSubj2] = size(data_5d);

if is_paired && nSubj1 ~= nSubj2
    error('组内设计被试数不一致！');
end

temp_avg = squeeze(mean(mean(data_5d(f_idx, t_idx, :, :, :), 1, 'omitnan'), 2, 'omitnan'));
pairs_data = zeros(nSubj2, nPairs);
for s = 1:nSubj2
    curr = temp_avg(:,:,s);
    for p = 1:nPairs
        pairs_data(s, p) = curr(pairs_idx(p,1), pairs_idx(p,2));
    end
end
Data_2 = pairs_data;
clear data_5d temp_avg pairs_data

%% 3. 全脑逐边统计 (Mass-univariate)
fprintf('开始全脑统计 (%s)...\n', stat_label);
h_bar = waitbar(0, 'Calculating...');

p_values = ones(1, nPairs);
t_values = zeros(1, nPairs); % 存一下 t 值

for p = 1:nPairs
    if mod(p, 500) == 0, waitbar(p/nPairs, h_bar); end
    
    y1 = Data_1(:, p);
    y2 = Data_2(:, p);
    
    if is_paired
        [~, p_val, ~, stats] = ttest(y1, y2);
    else
        [~, p_val, ~, stats] = ttest2(y1, y2);
    end
    p_values(p) = p_val;
    t_values(p) = stats.tstat;
end
close(h_bar);

%% 4. FDR 校正
fprintf('\n--- 统计结果 ---\n');
[~, mask] = fdr(p_values, 0.05);
n_sig = sum(mask);
fprintf('显著连接数 (FDR < 0.05): %d 条\n', n_sig);

if n_sig == 0
    warning('未发现显著连接。');
end

%% 5. 导出 BrainNet Viewer
choice_exp = questdlg('是否导出 BrainNet Viewer 3D 文件?', 'Export 3D', 'Yes', 'No', 'Yes');

if strcmp(choice_exp, 'Yes')
    save_dir = uigetdir(pwd, '选择保存结果的文件夹');
    if isequal(save_dir,0), return; end
    
    % --- A. Node 文件 ---
    node_file = fullfile(save_dir, 'Electrodes.node');
    fid = fopen(node_file, 'w');
    scale = 85; 
    for i = 1:nChans
        if isfield(chanlocs, 'X') && ~isempty(chanlocs(i).X)
            x=chanlocs(i).X; y=chanlocs(i).Y; z=chanlocs(i).Z;
        else
            [y,x,z] = pol2cart(chanlocs(i).theta*pi/180, chanlocs(i).radius);
            z = cos(chanlocs(i).radius * pi/2); 
        end
        vec = [x, y, z]; if norm(vec)>0, vec=vec/norm(vec)*scale; end
        % 颜色1 (红色)，大小2
        fprintf(fid, '%.4f %.4f %.4f 1 2 %s\n', vec(1), vec(2), vec(3), chanlocs(i).labels);
    end
    fclose(fid);
    
    % --- B. Edge 文件 ---
    if n_sig > 0
        edge_mat = zeros(nChans, nChans);
        sig_idx = find(mask == 1);
        
        for s = 1:length(sig_idx)
            idx = sig_idx(s);
            ch1 = pairs_idx(idx, 1);
            ch2 = pairs_idx(idx, 2);
            
            % 存入 t 值，便于用颜色区分增强/减弱
            val = t_values(idx); 
            edge_mat(ch1, ch2) = val;
            edge_mat(ch2, ch1) = val;
        end
        
        fname = fullfile(save_dir, [Name_1 '_vs_' Name_2 '_FDR.edge']);
        dlmwrite(fname, edge_mat, 'delimiter', '\t');
        
        msgbox({'导出完成！'; ''; ...
                '请使用 BrainNet Viewer 加载 .node 和 .edge 文件。'; ...
                '提示：Edge 文件保存的是 t 值。'; ...
                '建议在 BrainNet Viewer 中将 Edge Color 设置为 Jet，以区分正负差异。'});
    else
        msgbox('无显著结果，仅导出了 Node 文件。');
    end
end
