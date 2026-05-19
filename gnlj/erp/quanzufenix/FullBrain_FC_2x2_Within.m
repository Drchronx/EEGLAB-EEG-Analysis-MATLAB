% =========================================================================
%  全脑连接组学分析 (Whole-brain Connectome Analysis)
%  实验设计: 2x2 组内设计 (Within-Subject Design)
%  统计方法: 逐边重复测量方差分析 (Mass-univariate RM-ANOVA) + FDR 校正
%  输出: BrainNet Viewer (.node 和 .edge) 文件
% =========================================================================

clc; clear all; close all;

%% ================= 1. 参数设置与电极载入 =================
% 1.1 设置参数
prompt = { ...
    '输入感兴趣频率范围 (Hz) [Min Max]:', ...
    '输入感兴趣时间范围 (ms) [Min Max]:', ...
    '时间轴 (ms) [Start Step End] (需与计算FC时一致):', ...
    '频率轴 (Hz) [Start Step End] (需与计算FC时一致):', ...
    '因素A名称 (e.g. Angle):', ...
    '因素B名称 (e.g. Distance):' ...
};
% 默认值根据你的预处理猜测，请按实际修改
def = {'[4 7]', '[200 600]', '[-1000 4 1996]', '[1 1 30]', 'FactorA', 'FactorB'};
ans_set = inputdlg(prompt, '全脑分析参数', 1, def);
if isempty(ans_set), return; end

f_roi_range = str2num(ans_set{1});
t_roi_range = str2num(ans_set{2});
time_axis   = str2num(ans_set{3});
f_axis      = str2num(ans_set{4});
Name_A      = ans_set{5};
Name_B      = ans_set{6};

% 1.2 载入电极信息
disp('>>> [重要] 请选择一个 *预处理后* (已剔除电极) 的 .set 文件 <<<');
[fn, pn] = uigetfile('*.set', 'Load processed .set file');
if isequal(fn,0), return; end
EEG = pop_loadset('filename', fn, 'filepath', pn);
chanlocs = EEG.chanlocs;
nChans = length(chanlocs);
nPairs = nChans * (nChans - 1) / 2; % 下三角连接数

fprintf('检测到 %d 个电极，共 %d 条连接。\n', nChans, nPairs);

%% ================= 2. 载入功能连接数据 =================
% 逻辑：你需要有 4 个 .mat 文件，每个文件包含所有被试在该条件下的 5D 矩阵
% 文件顺序必须固定：A1B1, A1B2, A2B1, A2B2
% 假设文件内变量为 plv/pli/wpli [Freq, Time, Chan, Chan, Subj]

titles = {
    ['1. Cond: ' Name_A '1 & ' Name_B '1'], ...
    ['2. Cond: ' Name_A '1 & ' Name_B '2'], ...
    ['3. Cond: ' Name_A '2 & ' Name_B '1'], ...
    ['4. Cond: ' Name_A '2 & ' Name_B '2']
};

% 生成电极对索引
pairs_idx = nchoosek(1:nChans, 2); 

% 初始化数据容器
% 结构: Cell array {4} -> 每个里面是 [Subj x nPairs]
Data_Conditions = cell(1, 4); 
nSubj = 0;

for c = 1:4
    fprintf('正在读取条件 %d: %s ...\n', c, titles{c});
    [fn, pn] = uigetfile('*.mat', ['选择 ' titles{c}]);
    if isequal(fn,0), error('用户取消'); end
    
    tmp = load(fullfile(pn, fn));
    vars = fieldnames(tmp);
    % 自动寻找 5D 矩阵变量
    data_5d = tmp.(vars{1}); 
    
    [nF, nT, nCh, ~, nS] = size(data_5d);
    
    % 检查维度
    if nCh ~= nChans
        error('数据矩阵的通道数 (%d) 与 .set 文件的电极数 (%d) 不一致！请检查是否选对了 .set 文件。', nCh, nChans);
    end
    if c == 1
        nSubj = nS;
    elseif nS ~= nSubj
        error('不同条件下的被试数量不一致！');
    end
    
    % --- 核心步骤：提取 ROI 平均值 ---
    % 1. 找到 Time 和 Freq 的索引
    t_idx = time_axis >= t_roi_range(1) & time_axis <= t_roi_range(2);
    f_idx = f_axis >= f_roi_range(1) & f_axis <= f_roi_range(2);
    
    % 2. 切片并平均: [Freq, Time, Ch, Ch, Subj] -> [Ch, Ch, Subj]
    % 注意 mean 的维度变化
    temp_avg = squeeze(mean(mean(data_5d(f_idx, t_idx, :, :, :), 1, 'omitnan'), 2, 'omitnan'));
    
    % 3. 拉直为 [Subj x nPairs]
    pairs_data = zeros(nS, nPairs);
    for s = 1:nS
        curr_mat = temp_avg(:,:,s);
        for p = 1:nPairs
            ch1 = pairs_idx(p,1); 
            ch2 = pairs_idx(p,2);
            pairs_data(s, p) = curr_mat(ch1, ch2);
        end
    end
    
    Data_Conditions{c} = pairs_data;
    clear data_5d temp_avg
end

fprintf('数据载入完毕。N = %d\n', nSubj);

%% ================= 3. 全脑 2x2 组内 ANOVA =================
fprintf('开始全脑统计 (2x2 Within-Subject ANOVA)...\n');
h_bar = waitbar(0, '正在计算...');

% 初始化 P 值
p_A   = ones(1, nPairs);
p_B   = ones(1, nPairs);
p_Int = ones(1, nPairs);

% 预备矩阵计算，加速循环
% Y 结构: [Subj, 4] -> A1B1, A1B2, A2B1, A2B2
Y_all = zeros(nSubj, 4); 

for p = 1:nPairs
    if mod(p, 100) == 0, waitbar(p/nPairs, h_bar); end
    
    % 提取当前连接的所有条件数据
    for c = 1:4
        Y_all(:, c) = Data_Conditions{c}(:, p);
    end
    
    % 调用快速 ANOVA 函数 (在脚本末尾定义)
    [p1, p2, p3] = fast_rm_anova_2x2(Y_all);
    
    p_A(p)   = p1;
    p_B(p)   = p2;
    p_Int(p) = p3;
end
close(h_bar);

%% ================= 4. FDR 校正与结果报告 =================
fprintf('\n--- 统计结果 ---\n');

% FDR 校正 (Benjamini-Hochberg)
alpha = 0.05;
[~, mask_A]   = fdr(p_A, alpha);
[~, mask_B]   = fdr(p_B, alpha);
[~, mask_Int] = fdr(p_Int, alpha);

fprintf('主效应 %s 显著连接数: %d\n', Name_A, sum(mask_A));
fprintf('主效应 %s 显著连接数: %d\n', Name_B, sum(mask_B));
fprintf('交互作用显著连接数: %d\n', sum(mask_Int));

if sum(mask_A)+sum(mask_B)+sum(mask_Int) == 0
    warning('未发现任何经过 FDR 校正的显著连接。可能需要放宽阈值或改用 NBS 方法。');
end

%% ================= 5. 导出 BrainNet Viewer 文件 =================
choice_exp = questdlg('是否导出 BrainNet Viewer 文件?', 'Export 3D', 'Yes', 'No', 'Yes');

if strcmp(choice_exp, 'Yes')
    save_dir = uigetdir(pwd, '选择保存结果的文件夹');
    if isequal(save_dir,0), return; end
    
    % ---------------- A. 生成 Node 文件 ----------------
    % 这里的关键是处理被删掉电极后的坐标
    % BrainNet Viewer 需要 MNI 坐标。如果 chanlocs 是标准坐标，直接用。
    % 如果是头皮坐标，我们做一个简单的球形投影。
    
    node_file = fullfile(save_dir, 'Electrodes_Filtered.node');
    fid = fopen(node_file, 'w');
    
    scale = 85; % 缩放系数
    for i = 1:nChans
        if isfield(chanlocs, 'X') && ~isempty(chanlocs(i).X)
            x=chanlocs(i).X; y=chanlocs(i).Y; z=chanlocs(i).Z;
        else
            % 极坐标转笛卡尔 (近似)
            [y,x,z] = pol2cart(chanlocs(i).theta*pi/180, chanlocs(i).radius);
            z = cos(chanlocs(i).radius * pi/2); % 假定 z 轴
        end
        % 归一化并放大
        vec = [x, y, z]; 
        if norm(vec)>0, vec = vec/norm(vec)*scale; end
        
        % 格式: X Y Z Color Size Label
        fprintf(fid, '%.4f %.4f %.4f 1 2 %s\n', vec(1), vec(2), vec(3), chanlocs(i).labels);
    end
    fclose(fid);
    
    % ---------------- B. 生成 Edge 文件 ----------------
    masks = {mask_A, mask_B, mask_Int};
    names = {['Main_' Name_A], ['Main_' Name_B], 'Interaction'};
    
    for k = 1:3
        curr_mask = masks{k};
        if sum(curr_mask) > 0
            edge_mat = zeros(nChans, nChans);
            sig_idx = find(curr_mask == 1);
            
            for s = 1:length(sig_idx)
                idx = sig_idx(s);
                ch1 = pairs_idx(idx, 1);
                ch2 = pairs_idx(idx, 2);
                % 赋值为 1 (二值化网络) 或者 可以赋值为 F值/T值
                edge_mat(ch1, ch2) = 1;
                edge_mat(ch2, ch1) = 1;
            end
            
            fname = fullfile(save_dir, [names{k} '_FDR05.edge']);
            dlmwrite(fname, edge_mat, 'delimiter', '\t');
            fprintf('已导出: %s\n', fname);
        end
    end
    msgbox('导出完成！请使用 BrainNet Viewer 加载生成的 .node 和 .edge 文件。');
end

% =========================================================================
%  子函数: 快速 2x2 组内 ANOVA (Vectorized for Speed)
%  输入 Y: [Subj x 4] matrix (A1B1, A1B2, A2B1, A2B2)
% =========================================================================
function [pA, pB, pInt] = fast_rm_anova_2x2(Y)
    % 基于对比 (Contrasts) 的快速计算法，等价于 RM-ANOVA
    % A1B1 (1), A1B2 (2), A2B1 (3), A2B2 (4)
    
    [n, ~] = size(Y);
    
    % 1. 主效应 A (A1 vs A2): (Cond1+Cond2) - (Cond3+Cond4)
    % 实际上比较的是 (A1B1+A1B2)/2 vs (A2B1+A2B2)/2
    contrast_A = (Y(:,1) + Y(:,2)) - (Y(:,3) + Y(:,4));
    [~, pA] = ttest(contrast_A, 0); % 单样本 t 检验 (对比值是否显著不为0)
    
    % 2. 主效应 B (B1 vs B2): (Cond1+Cond3) - (Cond2+Cond4)
    contrast_B = (Y(:,1) + Y(:,3)) - (Y(:,2) + Y(:,4));
    [~, pB] = ttest(contrast_B, 0);
    
    % 3. 交互作用 (Difference of Differences)
    % (A1B1 - A1B2) - (A2B1 - A2B2)
    contrast_Int = (Y(:,1) - Y(:,2)) - (Y(:,3) - Y(:,4));
    [~, pInt] = ttest(contrast_Int, 0);
    
    % 注意: ttest 返回的 p 值与 F 检验的 p 值在 df=1 时是相同的 (F = t^2)
end