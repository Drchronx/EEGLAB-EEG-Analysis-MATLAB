% =========================================================================
%  静息态单因素功能连接 ROI 全程分析 (Resting-State Single-Factor ROI FC)
%  功能: 
%    1. 支持 组内 (Paired) 和 组间 (Independent) 模式切换
%    2. 自动提取特定电极对 (ROI)
%    3. 寻找显著的时间-频率窗口 (T-test + FDR)
%    4. 绘制 2D 时频图
%    5. 导出 SPSS 数据
%    6. 导出 BrainNet Viewer 3D 绘图文件
% =========================================================================

clc; clear all; close all;

%% ================= 1. 实验设置与数据载入 =================
% 1.1 选择设计类型
design_type = questdlg('请选择实验设计类型:', 'Design Selection', ...
    '组内 (Paired, 同一批人)', '组间 (Independent, 两组人)', '组间 (Independent, 两组人)');
if isempty(design_type), return; end
is_paired = strcmp(design_type, '组内 (Paired, 同一批人)');

% 1.2 参数设置
% 注意：静息态通常没有"刺激前基线"，这里的基线通常指全段平均或特定时间段
prompt = { ...
    '输入分析电极对 (e.g. [6 5]):', ...
    '时间轴 (ms/s) [Start Step End]:', ...
    '频率轴 (Hz) [Start Step End]:', ...
    '是否进行基线校正? (1=是, 0=否):', ...
    '基线范围 (ms/s) [Min Max] (如不校正可忽略):', ...
    '条件/组别 1 名称:', '条件/组别 2 名称:' ...
};
% 默认值示例 (假设是分段后的静息态数据)
def = {'[6 5]', '[0 4 2000]', '[1 1 30]', '0', '[0 200]', 'GroupA', 'GroupB'};
ans_set = inputdlg(prompt, 'Resting FC Parameters', 1, def);
if isempty(ans_set), return; end

chan_pair   = str2num(ans_set{1});
time_axis   = str2num(ans_set{2});
f_axis      = str2num(ans_set{3});
do_baseline = str2num(ans_set{4});
baseline    = str2num(ans_set{5});
Name_1      = ans_set{6};
Name_2      = ans_set{7};

% 1.3 载入电极信息 (用于3D画图)
disp('>>> [重要] 请选择一个 *预处理后* (已剔除电极) 的 .set 文件 <<<');
[fn_set, pn_set] = uigetfile('*.set', 'Load EEG .set file');
if isequal(fn_set,0), return; end
EEG = pop_loadset('filename', fn_set, 'filepath', pn_set);
chanlocs = EEG.chanlocs;
nChans_Total = length(chanlocs);

% 1.4 载入连接矩阵数据
% 容器: [Freq, Time, Subj]
Data_1 = []; Data_2 = [];
nSubj1 = 0; nSubj2 = 0;

if is_paired
    titles = {'Condition 1', 'Condition 2'};
else
    titles = {'Group 1', 'Group 2'};
end

% --- 读取数据 1 ---
fprintf('正在读取 %s ...\n', titles{1});
[fn1, pn1] = uigetfile('*.mat', ['Select: ' titles{1}]);
if isequal(fn1,0), return; end
tmp = load(fullfile(pn1, fn1)); vars = fieldnames(tmp);
data_5d = tmp.(vars{1}); % 假设格式: [Freq, Time, Ch, Ch, Subj]

[nF, nT, nCh, ~, nS1] = size(data_5d);
if nCh ~= nChans_Total, error('数据通道数与 .set 文件不匹配！'); end
nSubj1 = nS1;

% 提取 ROI
roi_raw = squeeze(data_5d(:,:,chan_pair(1), chan_pair(2), :)); % [Freq, Time, Subj]

% 基线校正 (可选)
if do_baseline
    base_idx = time_axis >= baseline(1) & time_axis <= baseline(2);
    base_mean = mean(roi_raw(:, base_idx, :), 2);
    Data_1 = roi_raw - repmat(base_mean, [1, nT, 1]);
else
    Data_1 = roi_raw;
end
clear data_5d roi_raw

% --- 读取数据 2 ---
fprintf('正在读取 %s ...\n', titles{2});
[fn2, pn2] = uigetfile('*.mat', ['Select: ' titles{2}]);
if isequal(fn2,0), return; end
tmp = load(fullfile(pn2, fn2)); vars = fieldnames(tmp);
data_5d = tmp.(vars{1});

[nF2, nT2, ~, ~, nS2] = size(data_5d);
nSubj2 = nS2;

if nF~=nF2 || nT~=nT2, error('两个文件的时间/频率维度不一致！'); end
if is_paired && nS1~=nS2, error('组内设计要求被试数量必须一致！'); end

% 提取 ROI
roi_raw = squeeze(data_5d(:,:,chan_pair(1), chan_pair(2), :));
if do_baseline
    base_mean = mean(roi_raw(:, base_idx, :), 2);
    Data_2 = roi_raw - repmat(base_mean, [1, nT, 1]);
else
    Data_2 = roi_raw;
end
clear data_5d roi_raw

fprintf('数据准备就绪。N1=%d, N2=%d\n', nSubj1, nSubj2);

%% ================= 2. 统计分析 (找显著窗口) =================
fprintf('正在进行逐点统计检验...\n');
P_Map = ones(nF, nT);
T_Map = zeros(nF, nT); 

h_bar = waitbar(0, 'Calculating Stats...');
for f = 1:nF
    waitbar(f/nF, h_bar);
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
close(h_bar);

% FDR 校正
[~, Mask] = fdr(P_Map, 0.05);

%% ================= 3. 绘制 2D 时频图 =================
figure('Color', 'w', 'Position', [100, 100, 1000, 700], 'Name', 'Resting State ROI Analysis');

% 1. Group/Cond 1 Mean
subplot(2,2,1);
imagesc(time_axis, f_axis, mean(Data_1, 3)); axis xy; colorbar;
title(['Mean: ' Name_1]); xlabel('Time/Segment'); ylabel('Freq (Hz)');

% 2. Group/Cond 2 Mean
subplot(2,2,2);
imagesc(time_axis, f_axis, mean(Data_2, 3)); axis xy; colorbar;
title(['Mean: ' Name_2]); xlabel('Time/Segment'); ylabel('Freq (Hz)');

% 3. Uncorrected P
subplot(2,2,3);
p_plot = P_Map; p_plot(P_Map >= 0.05) = NaN;
imagesc(time_axis, f_axis, p_plot); axis xy; colorbar;
colormap(gca, flipud(parula)); caxis([0 0.05]);
title('Uncorrected p < 0.05'); xlabel('Time'); ylabel('Freq');

% 4. FDR Corrected (显著时间窗看这里)
subplot(2,2,4);
p_fdr = P_Map; p_fdr(Mask == 0) = NaN;
if all(isnan(p_fdr(:)))
    text(mean(time_axis), mean(f_axis), 'No Sig. Area', 'Horiz', 'center');
    xlim([min(time_axis) max(time_axis)]); ylim([min(f_axis) max(f_axis)]);
else
    imagesc(time_axis, f_axis, p_fdr); axis xy; colorbar;
    caxis([0 0.05]);
end
title('FDR Corrected (Significant Windows)'); xlabel('Time'); ylabel('Freq');

msgbox('请观察图4 (右下角)，确定显著的时间和频率范围，并在下一步输入。');

%% ================= 4. 提取数据导出 Excel (SPSS) =================
choice_spss = questdlg('是否提取显著区域数据导出到 Excel?', 'Export SPSS', 'Yes', 'No', 'Yes');
mean_val1 = 0; mean_val2 = 0; % 用于3D画图的强度

if strcmp(choice_spss, 'Yes')
    prompt_r = {'输入显著时间范围 [min max]:', '输入显著频率范围 [min max]:'};
    ans_r = inputdlg(prompt_r, 'ROI Definition', 1, ...
        {sprintf('[%d %d]', time_axis(1), time_axis(end)), '[4 8]'});
    
    if ~isempty(ans_r)
        tr = str2num(ans_r{1}); fr = str2num(ans_r{2});
        t_idx = time_axis>=tr(1) & time_axis<=tr(2);
        f_idx = f_axis>=fr(1) & f_axis<=fr(2);
        
        % 提取 ROI 内的均值 (先平均Freq和Time，保留Subj)
        val1 = squeeze(mean(mean(Data_1(f_idx, t_idx, :), 1), 2));
        val2 = squeeze(mean(mean(Data_2(f_idx, t_idx, :), 1), 2));
        
        % 计算总均值 (用于3D画图的颜色强度)
        mean_val1 = mean(val1);
        mean_val2 = mean(val2);
        
        if is_paired
            T = table((1:nSubj1)', val1, val2, 'VariableNames', {'SubID', Name_1, Name_2});
        else
            g1_col = ones(nSubj1, 1);
            g2_col = ones(nSubj2, 1) * 2;
            T = table([(1:nSubj1)'; (1:nSubj2)'], [g1_col; g2_col], [val1; val2], ...
                'VariableNames', {'SubID', 'Group', 'Value'});
        end
        
        [fn, pn] = uiputfile('*.xlsx', 'Save SPSS Data');
        if fn~=0, writetable(T, fullfile(pn, fn)); end
    end
else
    % 如果不导出SPSS，默认取全时段平均用于画图
    mean_val1 = mean(Data_1(:));
    mean_val2 = mean(Data_2(:));
end

%% ================= 5. 导出 BrainNet Viewer 3D 文件 =================
% 这里的逻辑：只画这一条线。
% 生成两个 .edge 文件：
% 1. Condition 1 的强度
% 2. Condition 2 的强度
% (或者你可以选择导出差异值 t-value)

choice_3d = questdlg('是否导出 3D 脑图文件 (.node/.edge)?', '3D Export', 'Yes', 'No', 'Yes');
if strcmp(choice_3d, 'Yes')
    save_dir = uigetdir(pwd, '选择保存文件夹');
    if isequal(save_dir, 0), return; end
    
    % --- A. 生成 Node (高亮 ROI 电极) ---
    node_file = fullfile(save_dir, 'ROI_Electrodes.node');
    fid = fopen(node_file, 'w');
    scale = 85;
    for i = 1:nChans_Total
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
    
    % --- B. 生成 Edge 文件 ---
    % 我们导出 3 个文件：Cond1, Cond2, 和 Diff(t-value)
    
    % 1. Cond 1
    edge_1 = zeros(nChans_Total);
    edge_1(chan_pair(1), chan_pair(2)) = mean_val1;
    edge_1(chan_pair(2), chan_pair(1)) = mean_val1;
    dlmwrite(fullfile(save_dir, [Name_1 '.edge']), edge_1, 'delimiter', '\t');
    
    % 2. Cond 2
    edge_2 = zeros(nChans_Total);
    edge_2(chan_pair(1), chan_pair(2)) = mean_val2;
    edge_2(chan_pair(2), chan_pair(1)) = mean_val2;
    dlmwrite(fullfile(save_dir, [Name_2 '.edge']), edge_2, 'delimiter', '\t');
    
    % 3. Difference (T-value) -> 如果刚才选了 ROI，就用 ROI 内的 t 值
    if exist('f_idx', 'var')
        t_val = mean(mean(T_Map(f_idx, t_idx)));
    else
        t_val = max(T_Map(:));
    end
    edge_t = zeros(nChans_Total);
    edge_t(chan_pair(1), chan_pair(2)) = t_val;
    edge_t(chan_pair(2), chan_pair(1)) = t_val;
    dlmwrite(fullfile(save_dir, 'Diff_T_Value.edge'), edge_t, 'delimiter', '\t');
    
    msgbox({'导出成功！'; ''; ...
            ['1. Node文件: ROI_Electrodes.node']; ...
            ['2. Edge文件: ' Name_1 '.edge (显示条件1强度)']; ...
            ['3. Edge文件: ' Name_2 '.edge (显示条件2强度)']; ...
            ['4. Edge文件: Diff_T_Value.edge (显示差异显著性)']; ...
            ''; ...
            '请使用 Draw3D_Brain.m 脚本打开。'});
end
