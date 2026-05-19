% =========================================================================
%  单因素功能连接 ROI 分析 (Single-Factor ROI FC Analysis)
%  适用: 针对【特定电极对】进行时频统计分析
%  功能: 
%    1. 支持 组内 (Paired) 和 组间 (Independent) 切换
%    2. 逐时频点 t 检验 + FDR 校正
%    3. 绘制时频图 (TFR)
%    4. 提取数据供 SPSS 使用
%    5. 导出 BrainNet Viewer 3D 图
% =========================================================================

clc; clear all; close all;

%% 1. 实验设计与参数设置
% 1.1 选择设计类型
design_type = questdlg('请选择实验设计类型:', 'Design Selection', ...
    '组内 (Paired, 同一批人)', '组间 (Independent, 两组人)', '组内 (Paired, 同一批人)');
if isempty(design_type), return; end
is_paired = strcmp(design_type, '组内 (Paired, 同一批人)');

% 1.2 输入参数
prompt = { ...
    '输入分析电极对 (e.g. [6 5]):', ...
    '时间轴 (ms) [Start Step End]:', ...
    '频率轴 (Hz) [Start Step End]:', ...
    '基线范围 (ms) [Min Max]:', ...
    '条件/组别 1 名称:', '条件/组别 2 名称:' ...
};
def = {'[6 5]', '[-1000 4 1996]', '[1 1 30]', '[-800 -200]', 'CondA', 'CondB'};
ans_set = inputdlg(prompt, 'ROI FC Parameters', 1, def);
if isempty(ans_set), return; end

chan_pair = str2num(ans_set{1});
time_axis = str2num(ans_set{2});
f_axis    = str2num(ans_set{3});
baseline  = str2num(ans_set{4});
Name_1    = ans_set{5};
Name_2    = ans_set{6};

%% 2. 载入数据 (自动提取 ROI)
% 容器
Data_1 = []; % [Freq, Time, Subj]
Data_2 = []; 
nSubj1 = 0; nSubj2 = 0;

if is_paired
    titles = {'Condition 1', 'Condition 2'};
else
    titles = {'Group 1', 'Group 2'};
end

% --- 载入数据集 1 ---
fprintf('正在读取 %s ...\n', titles{1});
[fn1, pn1] = uigetfile('*.mat', ['Select: ' titles{1}]);
if isequal(fn1,0), return; end
tmp = load(fullfile(pn1, fn1)); vars = fieldnames(tmp);
data_5d = tmp.(vars{1}); % [Freq, Time, Ch, Ch, Subj]

[nF, nT, ~, ~, nS1] = size(data_5d);
nSubj1 = nS1;

% 提取 ROI 并基线校正
roi_raw = squeeze(data_5d(:,:,chan_pair(1), chan_pair(2), :));
base_idx = time_axis >= baseline(1) & time_axis <= baseline(2);
base_mean = mean(roi_raw(:, base_idx, :), 2);
Data_1 = roi_raw - repmat(base_mean, [1, nT, 1]);
clear data_5d roi_raw

% --- 载入数据集 2 ---
fprintf('正在读取 %s ...\n', titles{2});
[fn2, pn2] = uigetfile('*.mat', ['Select: ' titles{2}]);
if isequal(fn2,0), return; end
tmp = load(fullfile(pn2, fn2)); vars = fieldnames(tmp);
data_5d = tmp.(vars{1});

[nF2, nT2, ~, ~, nS2] = size(data_5d);
nSubj2 = nS2;

% 检查
if nF~=nF2 || nT~=nT2, error('两个文件的时间/频率维度不一致！'); end
if is_paired && nS1~=nS2, error('组内设计要求被试数量必须一致！'); end

% 提取 ROI 并基线校正
roi_raw = squeeze(data_5d(:,:,chan_pair(1), chan_pair(2), :));
base_mean = mean(roi_raw(:, base_idx, :), 2);
Data_2 = roi_raw - repmat(base_mean, [1, nT, 1]);
clear data_5d roi_raw

fprintf('数据准备就绪。N1=%d, N2=%d\n', nSubj1, nSubj2);

%% 3. 逐点统计分析 (t-test)
fprintf('正在进行时频统计...\n');
P_Map = ones(nF, nT);
T_Map = zeros(nF, nT); % 记录 t 值，用于看方向

h = waitbar(0, 'Calculating Statistics...');
for f = 1:nF
    waitbar(f/nF, h);
    for t = 1:nT
        y1 = squeeze(Data_1(f, t, :));
        y2 = squeeze(Data_2(f, t, :));
        
        if is_paired
            [~, p, ~, stats] = ttest(y1, y2);
        else
            [~, p, ~, stats] = ttest2(y1, y2);
        end
        P_Map(f, t) = p;
        T_Map(f, t) = stats.tstat;
    end
end
close(h);

%% 4. FDR 校正与绘图
[~, Mask] = fdr(P_Map, 0.05);

% --- 绘图 ---
figure('Color', 'w', 'Position', [100, 100, 1000, 700], 'Name', 'Single Factor ROI Analysis');

% 1. Condition 1 Mean
subplot(2,2,1);
imagesc(time_axis, f_axis, mean(Data_1, 3)); axis xy; colorbar;
title(['Mean: ' Name_1]); xlabel('Time'); ylabel('Freq');

% 2. Condition 2 Mean
subplot(2,2,2);
imagesc(time_axis, f_axis, mean(Data_2, 3)); axis xy; colorbar;
title(['Mean: ' Name_2]); xlabel('Time'); ylabel('Freq');

% 3. Uncorrected P
subplot(2,2,3);
p_plot = P_Map; p_plot(P_Map >= 0.05) = NaN;
imagesc(time_axis, f_axis, p_plot); axis xy; colorbar;
colormap(gca, flipud(parula)); caxis([0 0.05]);
title('Uncorrected p < 0.05');

% 4. FDR Corrected
subplot(2,2,4);
p_fdr = P_Map; p_fdr(Mask == 0) = NaN;
if all(isnan(p_fdr(:)))
    text(mean(time_axis), mean(f_axis), 'No Sig. Area', 'Horiz', 'center');
    xlim([min(time_axis) max(time_axis)]); ylim([min(f_axis) max(f_axis)]);
else
    imagesc(time_axis, f_axis, p_fdr); axis xy; colorbar;
    caxis([0 0.05]);
end
title('FDR Corrected p < 0.05');

%% 5. 提取数据导出 Excel (SPSS)
choice_spss = questdlg('是否提取显著区域数据导出到 Excel?', 'Export SPSS', 'Yes', 'No', 'Yes');
if strcmp(choice_spss, 'Yes')
    prompt_r = {'Time Range [min max]:', 'Freq Range [min max]:'};
    ans_r = inputdlg(prompt_r, 'ROI Definition', 1, {'[200 400]', '[4 8]'});
    if ~isempty(ans_r)
        tr = str2num(ans_r{1}); fr = str2num(ans_r{2});
        t_idx = time_axis>=tr(1) & time_axis<=tr(2);
        f_idx = f_axis>=fr(1) & f_axis<=fr(2);
        
        % 提取均值
        val1 = squeeze(mean(mean(Data_1(f_idx, t_idx, :), 1), 2));
        val2 = squeeze(mean(mean(Data_2(f_idx, t_idx, :), 1), 2));
        
        if is_paired
            % 宽格式: ID, Cond1, Cond2
            T = table((1:nSubj1)', val1, val2, 'VariableNames', {'SubID', Name_1, Name_2});
        else
            % 长格式: ID, Group, Value
            g1_col = ones(nSubj1, 1);
            g2_col = ones(nSubj2, 1) * 2;
            T = table([(1:nSubj1)'; (1:nSubj2)'], [g1_col; g2_col], [val1; val2], ...
                'VariableNames', {'SubID', 'Group', 'Value'});
        end
        
        [fn, pn] = uiputfile('*.xlsx', 'Save SPSS Data');
        if fn~=0, writetable(T, fullfile(pn, fn)); end
    end
end

%% 6. 导出 BrainNet Viewer 3D 文件 (Node + Edge)
choice_3d = questdlg('是否导出 3D 脑图文件 (.node/.edge)?', '3D Export', 'Yes', 'No', 'Yes');
if strcmp(choice_3d, 'Yes')
    disp('>>> 请选择一个 .set 文件以获取电极坐标 <<<');
    [fn, pn] = uigetfile('*.set', 'Load .set file');
    if isequal(fn,0), return; end
    EEG = pop_loadset('filename', fn, 'filepath', pn);
    chanlocs = EEG.chanlocs;
    nChans = length(chanlocs);
    
    save_dir = uigetdir(pwd, '选择保存文件夹');
    if isequal(save_dir, 0), return; end
    
    % --- A. 生成 Node (高亮 ROI 电极) ---
    node_file = fullfile(save_dir, 'ROI_Electrodes.node');
    fid = fopen(node_file, 'w');
    scale = 85;
    for i = 1:nChans
        if isfield(chanlocs, 'X') && ~isempty(chanlocs(i).X)
            x=chanlocs(i).X; y=chanlocs(i).Y; z=chanlocs(i).Z;
        else
            [y,x,z] = pol2cart(chanlocs(i).theta*pi/180, chanlocs(i).radius);
            z = cos(chanlocs(i).radius*pi/2);
        end
        vec = [x, y, z]; if norm(vec)>0, vec=vec/norm(vec)*scale; end
        
        % 颜色逻辑: 目标电极=红色(1), 其他=灰色(2)
        if i == chan_pair(1) || i == chan_pair(2)
            col=1; sz=5;
        else
            col=2; sz=2;
        end
        fprintf(fid, '%.4f %.4f %.4f %d %d %s\n', vec(1), vec(2), vec(3), col, sz, chanlocs(i).labels);
    end
    fclose(fid);
    
    % --- B. 生成 Edge (差异强度) ---
    % 计算全时频段的平均差异，或者只计算 ROI 差异
    % 这里我们导出 ROI 区域内的 t 值，这样颜色可以表示差异方向
    if exist('t_idx', 'var')
        % 如果刚才选了 ROI，就用 ROI 内的平均 t 值
        roi_t = mean(mean(T_Map(f_idx, t_idx)));
    else
        % 否则用全图最大 t 值代表
        roi_t = max(T_Map(:));
    end
    
    edge_mat = zeros(nChans, nChans);
    edge_mat(chan_pair(1), chan_pair(2)) = roi_t;
    edge_mat(chan_pair(2), chan_pair(1)) = roi_t;
    
    edge_file = fullfile(save_dir, 'ROI_Diff.edge');
    dlmwrite(edge_file, edge_mat, 'delimiter', '\t');
    
    msgbox({'导出成功！'; ''; 'Node文件: ROI_Electrodes.node'; 'Edge文件: ROI_Diff.edge (值为 t-value)'; ...
            '请使用 Draw3D_Brain.m 脚本打开。'});
end
