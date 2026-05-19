% =========================================================================
%  2x2 混合设计功能连接分析 (Mixed Design ANOVA for FC)
%  设计: 1个组间因子 (Between) x 1个组内因子 (Within)
%  适用于: 比如 (病人组 vs 对照组) x (条件A vs 条件B)
% =========================================================================

clc; clear all; close all;

%% ================= 1. 参数设置 =================
% 弹出窗口获取参数
prompt = { ...
    '输入电极对 (e.g., [6 5]):', ...
    '时间轴 (ms, e.g., -1000:4:1996):', ...
    '频率轴 (Hz, e.g., 1:30):', ...
    '基线范围 (ms, e.g., [-800 -200]):', ...
    '组间因子名称 (Between, e.g., Group):', ...
    '组内因子名称 (Within, e.g., Condition):' ...
};
definput = {'[6 5]', '-1000:4:1996', '1:30', '[-800 -200]', 'Group', 'Condition'};
answer = inputdlg(prompt, '混合设计参数设置', [1 60], definput);

if isempty(answer), return; end

chan_pair = str2num(answer{1});
time_axis = str2num(answer{2});
f_axis    = str2num(answer{3});
baseline  = str2num(answer{4});
Name_Btwn = answer{5}; % e.g., Group
Name_Wthn = answer{6}; % e.g., Emotion

%% ================= 2. 载入数据 (关键步骤) =================
% 逻辑：
% Group 1 的人 和 Group 2 的人是不同的。
% 必须分别载入：
% 1. G1_Cond1
% 2. G1_Cond2
% 3. G2_Cond1
% 4. G2_Cond2

disp('>>> 注意：混合设计数据载入顺序 <<<');
disp(['1. ' Name_Btwn ' Level 1 (e.g. Controls) - ' Name_Wthn ' Level 1']);
disp(['2. ' Name_Btwn ' Level 1 (e.g. Controls) - ' Name_Wthn ' Level 2']);
disp(['3. ' Name_Btwn ' Level 2 (e.g. Patients) - ' Name_Wthn ' Level 1']);
disp(['4. ' Name_Btwn ' Level 2 (e.g. Patients) - ' Name_Wthn ' Level 2']);

titles = {
    ['1. Select: ' Name_Btwn '1 & ' Name_Wthn '1'], ...
    ['2. Select: ' Name_Btwn '1 & ' Name_Wthn '2'], ...
    ['3. Select: ' Name_Btwn '2 & ' Name_Wthn '1'], ...
    ['4. Select: ' Name_Btwn '2 & ' Name_Wthn '2']
};

% 容器初始化
Data_G1 = []; % 将存储 Group 1 的数据 [Freq, Time, Subj_G1, 2]
Data_G2 = []; % 将存储 Group 2 的数据 [Freq, Time, Subj_G2, 2]

% --- 载入 Group 1 的两个条件 ---
for c = 1:2
    fprintf('正在载入 Group 1 - Condition %d ...\n', c);
    [fn, pn] = uigetfile('*.mat', titles{c});
    if isequal(fn,0), error('用户取消'); end
    
    tmp = load(fullfile(pn, fn));
    vars = fieldnames(tmp);
    raw_5d = tmp.(vars{1}); % [Freq, Time, Ch, Ch, Subj]
    
    % 维度检查与提取
    [nF, nT, ~, ~, nS] = size(raw_5d);
    if c == 1
        nSubj_G1 = nS;
        Data_G1 = zeros(nF, nT, nSubj_G1, 2);
        % 检查轴
        if nF~=length(f_axis) || nT~=length(time_axis)
             error('Group1 数据维度与定义的时间/频率轴不符！');
        end
    else
        if nS ~= nSubj_G1, error('Group 1 的两个条件被试数量不一致！'); end
    end
    
    % 提取 ROI 并基线校正
    roi_data = squeeze(raw_5d(:,:,chan_pair(1), chan_pair(2),:));
    base_idx = time_axis >= baseline(1) & time_axis <= baseline(2);
    base_mean = mean(roi_data(:, base_idx, :), 2);
    Data_G1(:,:,:,c) = roi_data - repmat(base_mean, [1, nT, 1]);
    
    clear raw_5d roi_data tmp
end

% --- 载入 Group 2 的两个条件 ---
for c = 1:2
    fprintf('正在载入 Group 2 - Condition %d ...\n', c);
    [fn, pn] = uigetfile('*.mat', titles{c+2}); % titles 3 & 4
    if isequal(fn,0), error('用户取消'); end
    
    tmp = load(fullfile(pn, fn));
    vars = fieldnames(tmp);
    raw_5d = tmp.(vars{1});
    
    [nF, nT, ~, ~, nS] = size(raw_5d);
    if c == 1
        nSubj_G2 = nS; % Group 2 可以和 Group 1 人数不同
        Data_G2 = zeros(nF, nT, nSubj_G2, 2);
    else
        if nS ~= nSubj_G2, error('Group 2 的两个条件被试数量不一致！'); end
    end
    
    % 提取 ROI 并基线校正
    roi_data = squeeze(raw_5d(:,:,chan_pair(1), chan_pair(2),:));
    base_idx = time_axis >= baseline(1) & time_axis <= baseline(2);
    base_mean = mean(roi_data(:, base_idx, :), 2);
    Data_G2(:,:,:,c) = roi_data - repmat(base_mean, [1, nT, 1]);
    
    clear raw_5d roi_data tmp
end

fprintf('数据载入完成。\n Group 1 (N=%d)\n Group 2 (N=%d)\n', nSubj_G1, nSubj_G2);

%% ================= 3. 逐点混合方差分析 (Mixed ANOVA) =================
%% ================= 3. 逐点混合方差分析 (Mixed ANOVA) =================
% 修改说明：已将计算公式内置到循环中，不再调用外部函数，彻底解决 "无法识别函数" 报错。

fprintf('开始计算 2x2 Mixed ANOVA (Group x Condition)...\n');

% 获取维度信息
[nF, nT, nSubj_G1, ~] = size(Data_G1);
[~, ~, nSubj_G2, ~]   = size(Data_G2);

% 初始化 P 值矩阵
P_Btwn = ones(nF, nT); % 组间主效应 (Group)
P_Wthn = ones(nF, nT); % 组内主效应 (Condition)
P_Int  = ones(nF, nT); % 交互作用 (Group x Condition)

h_bar = waitbar(0, 'Running Mixed ANOVA (Inline Calculation)...');

% --- 预计算自由度 (Degrees of Freedom) ---
% 自由度只跟人数有关，不需要在循环里重复算
df_grp = 1;
df_subj_within_grp = (nSubj_G1 - 1) + (nSubj_G2 - 1); % Error term for Between
df_cond = 1;
df_int = 1;
df_error_within = df_subj_within_grp; % Error term for Within

% 总人数
N_total = nSubj_G1 + nSubj_G2;

for f = 1:nF
    waitbar(f/nF, h_bar, sprintf('Freq %d / %d', f, nF));
    
    for t = 1:nT
        % 1. 提取当前时频点的数据
        % Y1: Group 1 数据 [nSubj_G1 x 2] -> [Cond1, Cond2]
        Y1 = squeeze(Data_G1(f, t, :, :)); 
        % Y2: Group 2 数据 [nSubj_G2 x 2] -> [Cond1, Cond2]
        Y2 = squeeze(Data_G2(f, t, :, :));
        
        % 处理单被试可能导致维度变成向量的情况
        if nSubj_G1 == 1, Y1 = Y1(:)'; end
        if nSubj_G2 == 1, Y2 = Y2(:)'; end

        % 2. 计算各种均值
        % 组均值 (Group Means)
        m_g1 = mean(Y1(:)); 
        m_g2 = mean(Y2(:));
        
        % 总均值 (Grand Mean)
        grand_mean = (sum(Y1(:)) + sum(Y2(:))) / (nSubj_G1*2 + nSubj_G2*2);
        
        % 条件均值 (Condition Means)
        % 注意：如果是由于组人数不等，不能简单 mean([Y1; Y2])，需要加权平均或者由 GLM 决定
        % 这里采用 Type III SS 的思路 (Unweighted means for factors)
        % Cond 1 total mean
        m_c1 = (mean(Y1(:,1)) * nSubj_G1 + mean(Y2(:,1)) * nSubj_G2) / N_total;
        m_c2 = (mean(Y1(:,2)) * nSubj_G1 + mean(Y2(:,2)) * nSubj_G2) / N_total;
        
        % 单元格均值 (Cell Means)
        cell_11 = mean(Y1(:,1)); cell_12 = mean(Y1(:,2));
        cell_21 = mean(Y2(:,1)); cell_22 = mean(Y2(:,2));
        
        % 被试均值 (Subject Means - across conditions)
        subj_m_g1 = mean(Y1, 2);
        subj_m_g2 = mean(Y2, 2);

        % 3. 计算平方和 (Sum of Squares)
        
        % --- SS Between (Group Effect) ---
        % SS_Group = 2 * n1 * (m_g1 - GM)^2 + ...
        SS_grp = 2 * nSubj_G1 * (m_g1 - grand_mean)^2 + ...
                 2 * nSubj_G2 * (m_g2 - grand_mean)^2;
             
        % --- SS Error Between (Subject within Group) ---
        SS_subj_wg = 2 * sum((subj_m_g1 - m_g1).^2) + ...
                     2 * sum((subj_m_g2 - m_g2).^2);
                 
        % --- SS Within (Condition Effect) ---
        SS_cond = N_total * (m_c1 - grand_mean)^2 + ...
                  N_total * (m_c2 - grand_mean)^2;
        
        % --- SS Interaction ---
        % SS_cells = sum(n * (cell_mean - grand_mean)^2)
        SS_cells = nSubj_G1*(cell_11-grand_mean)^2 + nSubj_G1*(cell_12-grand_mean)^2 + ...
                   nSubj_G2*(cell_21-grand_mean)^2 + nSubj_G2*(cell_22-grand_mean)^2;
        SS_int = SS_cells - SS_grp - SS_cond;
        
        % --- SS Error Within ---
        % 计算残差平方和
        SS_w_g1 = sum(sum((Y1 - repmat(subj_m_g1,1,2)).^2));
        SS_w_g2 = sum(sum((Y2 - repmat(subj_m_g2,1,2)).^2));
        SS_error_w = (SS_w_g1 + SS_w_g2) - SS_cond - SS_int;

        % 4. 计算 F 值和 P 值
        
        % Mean Squares
        MS_grp = SS_grp / df_grp;
        MS_subj_wg = SS_subj_wg / df_subj_within_grp;
        MS_cond = SS_cond / df_cond;
        MS_int = SS_int / df_int;
        MS_error_w = SS_error_w / df_error_within;
        
        % F Statistics
        % 避免分母为0的情况
        if MS_subj_wg < 1e-10, MS_subj_wg = 1e-10; end
        if MS_error_w < 1e-10, MS_error_w = 1e-10; end
        
        F_grp  = MS_grp / MS_subj_wg;
        F_cond = MS_cond / MS_error_w;
        F_int  = MS_int / MS_error_w;
        
        % P Values (使用 fcdf)
        % P = 1 - cdf
        P_Btwn(f,t) = 1 - fcdf(F_grp, df_grp, df_subj_within_grp);
        P_Wthn(f,t) = 1 - fcdf(F_cond, df_cond, df_error_within);
        P_Int(f,t)  = 1 - fcdf(F_int, df_int, df_error_within);
    end
end
close(h_bar);

% 防止数值误差出现负数P值
P_Btwn(P_Btwn < 0) = 0; P_Wthn(P_Wthn < 0) = 0; P_Int(P_Int < 0) = 0;

% FDR 校正
[~, Mask_Btwn] = fdr(P_Btwn, 0.05);
[~, Mask_Wthn] = fdr(P_Wthn, 0.05);
[~, Mask_Int]  = fdr(P_Int, 0.05);

disp('混合方差分析计算完成。');

%% ================= 4. 基础绘图 (TF Maps) =================
% 定义画图结构体
Plots = {
    struct('p', P_Btwn, 'mask', Mask_Btwn, 'title', ['Main Effect: ' Name_Btwn]), ...
    struct('p', P_Wthn, 'mask', Mask_Wthn, 'title', ['Main Effect: ' Name_Wthn]), ...
    struct('p', P_Int,  'mask', Mask_Int,  'title', ['Interaction: ' Name_Btwn ' x ' Name_Wthn])
};

for i = 1:3
    figure('Color','w', 'Name', Plots{i}.title);
    
    % 上图：Uncorrected P
    subplot(2,1,1);
    p_raw = Plots{i}.p;
    p_raw(p_raw >= 0.05) = NaN;
    imagesc(time_axis, f_axis, p_raw); axis xy; colorbar;
    colormap(gca, flipud(parula)); caxis([0 0.05]);
    title([Plots{i}.title ' (Uncorrected p<0.05)']);
    
    % 下图：FDR Corrected
    subplot(2,1,2);
    p_fdr = Plots{i}.p;
    p_fdr(Plots{i}.mask == 0) = NaN;
    if all(isnan(p_fdr(:)))
        text(mean(time_axis), mean(f_axis), 'No Significant Area (FDR)', 'Horiz','center');
        xlim([min(time_axis) max(time_axis)]); ylim([min(f_axis) max(f_axis)]);
    else
        imagesc(time_axis, f_axis, p_fdr); axis xy; colorbar;
        caxis([0 0.05]);
    end
    title([Plots{i}.title ' (FDR Corrected)']);
end



%% ================= 4. 高级绘图 (美化版 - 仿3D质感) =================
choice_plot = questdlg('是否进行高级绘图 (美化版)?', 'Advanced Plot', 'Yes', 'No', 'Yes');

if strcmp(choice_plot, 'Yes')
    % 1. 获取电极坐标
    if ~exist('chanlocs', 'var') || isempty(chanlocs)
        [fn, pn] = uigetfile('*.set', 'Load .set file for locations');
        if isequal(fn,0), return; end
        EEG_tmp = pop_loadset('filename', fn, 'filepath', pn);
        chanlocs = EEG_tmp.chanlocs;
    end
    
    % 2. 计算绘图数据 (均值)
    % 确保前面已经运行过ANOVA并有了 f_idx, t_idx
    if ~exist('t_idx','var')
        % 如果未定义，临时手动输入
        prompt_roi = {'Time Range [min max]:', 'Freq Range [min max]:'};
        ans_roi = inputdlg(prompt_roi, 'Plot ROI', 1, {'[200 400]', '[4 8]'});
        if isempty(ans_roi), return; end
        t_r = str2num(ans_roi{1}); f_r = str2num(ans_roi{2});
        t_idx = time_axis >= t_r(1) & time_axis <= t_r(2);
        f_idx = f_axis >= f_r(1) & f_axis <= f_r(2);
    end
    
    % 计算四个条件的均值
    % Data_G1: [Freq, Time, Subj, 2]
    m_g1c1 = mean(mean(mean(Data_G1(f_idx, t_idx, :, 1), 1), 2), 3);
    m_g1c2 = mean(mean(mean(Data_G1(f_idx, t_idx, :, 2), 1), 2), 3);
    m_g2c1 = mean(mean(mean(Data_G2(f_idx, t_idx, :, 1), 1), 2), 3);
    m_g2c2 = mean(mean(mean(Data_G2(f_idx, t_idx, :, 2), 1), 2), 3);
    means = [m_g1c1, m_g1c2, m_g2c1, m_g2c2];
    
    titles_plot = {
        [Name_Btwn '1 - ' Name_Wthn '1'], [Name_Btwn '1 - ' Name_Wthn '2'], ...
        [Name_Btwn '2 - ' Name_Wthn '1'], [Name_Btwn '2 - ' Name_Wthn '2']
    };
    
    % 3. 准备连线坐标
    target_pair = chan_pair; % [Ch1, Ch2]
    
    % --- 开始美化绘图 ---
    figure('Color', 'w', 'Position', [100, 100, 900, 750], 'Name', 'Publication Quality Topo');
    
    % 颜色设置
    max_val = max(abs(means)); if max_val==0, max_val=1; end
    c_lim = [-max_val*1.2, max_val*1.2];
    
    % 使用更高级的配色方案 (Matlab自带的 parula 或者 jet)
    % 建议：如果要仿照你发的图，可以用红色系表示强连接
    cmap = colormap(parula(256)); 
    
    for i = 1:4
        subplot(2,2,i);
        
        % 关键美化参数：
        % 'style', 'map': 填色图
        % 'electrodes', 'off': 隐藏杂乱的电极点
        % 'headrad', 'rim': 绘制头部轮廓
        % 'shading', 'interp': 极度平滑的插值，去除锯齿
        % 'numcontour', 0: 去除等高线，更干净
        % 'whitebk', 'on': 白色背景
        topoplot([], chanlocs, 'style', 'blank', ...
            'electrodes', 'off', ... 
            'headrad', 'rim', ...
            'shading', 'interp', ...
            'whitebk', 'on', ...
            'noplot', 'on'); % 先画个空架子获取坐标
            
        hold on;
        % 重新画轮廓
        topoplot([], chanlocs, 'style', 'blank', 'electrodes', 'off', 'headrad', 0.5);
        
        % 获取连线坐标
        % 这里使用更稳健的方法获取坐标
        [~, Grid, plotrad, xmesh, ymesh] = topoplot([], chanlocs, 'style', 'blank', 'noplot', 'on');
        x = xmesh(1,:); y = ymesh(:,1); 
        % 找到对应电极的最近坐标
        % (简化处理：直接用chanlocs的坐标转2D)
        [y_ch, x_ch] = pol2cart([chanlocs.theta]*pi/180, [chanlocs.radius]);
        x1 = x_ch(target_pair(1)); y1 = y_ch(target_pair(1));
        x2 = x_ch(target_pair(2)); y2 = y_ch(target_pair(2));

        % 根据数值设定连线颜色和粗细
        val = means(i);
        % 归一化颜色
        c_idx = round( (val - c_lim(1)) / (c_lim(2) - c_lim(1)) * 255 ) + 1;
        c_idx = max(1, min(256, c_idx));
        line_color = cmap(c_idx, :);
        
        % 画出像“光束”一样的连接线
        % 技巧：画一条粗的半透明线做光晕，再画一条细实线
        plot([x1 x2], [y1 y2], 'Color', [line_color, 0.3], 'LineWidth', 6); % 光晕
        
        if mod(i,2)==1, ls='-'; else, ls='--'; end
        plot([x1 x2], [y1 y2], 'Color', line_color, 'LineWidth', 2.5, 'LineStyle', ls); % 核心线
        
        % 画两个端点（像你图里的高亮区域）
        scatter([x1 x2], [y1 y2], 80, line_color, 'filled', 'MarkerEdgeColor', 'k');
        
        axis equal; axis off;
        title(titles_plot{i}, 'FontSize', 14, 'FontWeight', 'bold');
        
        % 在下方标数值
        text(0, -0.6, sprintf('FC: %.3f', val), 'Horiz','center', 'FontSize', 12);
    end
    
    % 美化 Colorbar
    cb = colorbar('Position', [0.92, 0.15, 0.02, 0.7]);
    ylabel(cb, 'Connectivity Strength', 'FontSize', 12, 'FontWeight', 'bold');
    caxis(c_lim);
    colormap(cmap);
    
    sgtitle(['2x2 Mixed Design Connectivity'], 'FontSize', 16, 'FontWeight', 'bold');
    disp('美化绘图完成。');
end


%% ================= 5. 数据提取 (SPSS 格式) =================
choice = questdlg('是否提取ROI数据导出Excel (SPSS)?', 'Export', 'Yes', 'No', 'Yes');
if strcmp(choice, 'Yes')
    prompt_roi = {'Time Range [min max]:', 'Freq Range [min max]:'};
    def_roi = {'[200 400]', '[4 8]'};
    ans_roi = inputdlg(prompt_roi, 'ROI Definition', 1, def_roi);
    
    if ~isempty(ans_roi)
        t_r = str2num(ans_roi{1}); f_r = str2num(ans_roi{2});
        t_idx = time_axis >= t_r(1) & time_axis <= t_r(2);
        f_idx = f_axis >= f_r(1) & f_axis <= f_r(2);
        
        % 提取 Group 1
        % 先求ROI均值: [1, 1, Subj, 2] -> [Subj, 2]
        temp1 = squeeze(mean(mean(Data_G1(f_idx, t_idx, :, :), 1, 'omitnan'), 2, 'omitnan'));
        if nSubj_G1 == 1, temp1 = temp1(:)'; end
        
        % 提取 Group 2
        temp2 = squeeze(mean(mean(Data_G2(f_idx, t_idx, :, :), 1, 'omitnan'), 2, 'omitnan'));
        if nSubj_G2 == 1, temp2 = temp2(:)'; end
        
        % 构造表格
        % 格式: [ID, Group, Cond1, Cond2]
        % Group 1: Label=1, Group 2: Label=2
        T1 = table((1:nSubj_G1)', ones(nSubj_G1,1), temp1(:,1), temp1(:,2), ...
            'VariableNames', {'SubID', 'Group', 'Cond1', 'Cond2'});
        T2 = table((1:nSubj_G2)', ones(nSubj_G2,1)*2, temp2(:,1), temp2(:,2), ...
            'VariableNames', {'SubID', 'Group', 'Cond1', 'Cond2'});
        
        T_Final = [T1; T2];
        
        [fn, pn] = uiputfile('*.xlsx', 'Save ROI Data');
        if fn~=0
            writetable(T_Final, fullfile(pn, fn));
            msgbox('SPSS数据已保存。Group 1=1, Group 2=2.');
        end
    end
end

%% ================= 6. 高级绘图 (网络拓扑 + 柱状图) =================
% 整合了你之前提供的绘图代码
choice_plot = questdlg('是否进行高级绘图 (拓扑图+柱状图)?', 'Advanced Plot', 'Yes', 'No', 'Yes');

if strcmp(choice_plot, 'Yes')
    % 1. 获取电极坐标
    disp('选择 .set 文件以获取电极位置...');
    [fn, pn] = uigetfile('*.set', 'Load .set file');
    if isequal(fn,0), warning('无坐标文件，跳过拓扑图'); chanlocs=[]; else
        EEG_tmp = pop_loadset('filename', fn, 'filepath', pn);
        chanlocs = EEG_tmp.chanlocs;
    end
    
    % 2. 再次确认ROI (用于画图)
    if ~exist('t_idx','var') % 如果第5步没做
        prompt_roi = {'Time Range [min max]:', 'Freq Range [min max]:'};
        ans_roi = inputdlg(prompt_roi, 'Plot ROI', 1, {'[200 400]', '[4 8]'});
        t_r = str2num(ans_roi{1}); f_r = str2num(ans_roi{2});
        t_idx = time_axis >= t_r(1) & time_axis <= t_r(2);
        f_idx = f_axis >= f_r(1) & f_axis <= f_r(2);
    end
    
    % 3. 计算用于绘图的均值和标准误
    % 需要 4 个值: G1C1, G1C2, G2C1, G2C2
    
    % Group 1
    g1_raw = Data_G1(f_idx, t_idx, :, :); % [ROI, Subj, 2]
    % 平均频时 -> [1, 1, Subj, 2] -> squeeze -> [Subj, 2]
    g1_vals = squeeze(mean(mean(g1_raw, 1, 'omitnan'), 2, 'omitnan')); 
    if nSubj_G1==1, g1_vals = g1_vals(:)'; end
    
    % Group 2
    g2_raw = Data_G2(f_idx, t_idx, :, :);
    g2_vals = squeeze(mean(mean(g2_raw, 1, 'omitnan'), 2, 'omitnan'));
    if nSubj_G2==1, g2_vals = g2_vals(:)'; end
    
    % 均值
    m_g1 = mean(g1_vals, 1); % [Mean_C1, Mean_C2]
    m_g2 = mean(g2_vals, 1);
    means = [m_g1, m_g2]; % Order: G1C1, G1C2, G2C1, G2C2
    
    % 标准误
    se_g1 = std(g1_vals, 0, 1) / sqrt(nSubj_G1);
    se_g2 = std(g2_vals, 0, 1) / sqrt(nSubj_G2);
    ses = [se_g1, se_g2];
    
    titles_plot = {
        [Name_Btwn '1 - ' Name_Wthn '1'], [Name_Btwn '1 - ' Name_Wthn '2'], ...
        [Name_Btwn '2 - ' Name_Wthn '1'], [Name_Btwn '2 - ' Name_Wthn '2']
    };

    % --- A. 2x2 拓扑图 ---
    if ~isempty(chanlocs)
        % 获取坐标 (自动适配)
        fig_hid = figure('Visible','off'); 
        topoplot([], chanlocs, 'style', 'blank', 'electrodes', 'on');
        h_pts = findobj(gca, 'Type','line', 'Marker','.');
        if isempty(h_pts), h_pts = findobj(gca, 'Type','scatter'); end
        try
            x_all = h_pts.XData; y_all = h_pts.YData;
            x1 = x_all(chan_pair(1)); y1 = y_all(chan_pair(1));
            x2 = x_all(chan_pair(2)); y2 = y_all(chan_pair(2));
        catch
            warning('坐标提取失败，尝试极坐标转换');
            [y, x] = pol2cart([chanlocs.theta]*pi/180, [chanlocs.radius]);
            x1=x(chan_pair(1)); y1=y(chan_pair(1));
            x2=x(chan_pair(2)); y2=y(chan_pair(2));
        end
        close(fig_hid);
        
        % 绘图参数
        max_v = max(abs(means)); if max_v==0, max_v=1; end
        c_lim = [-max_v*1.2, max_v*1.2];
        cmap = jet(256);
        
        figure('Color','w', 'Position', [100, 100, 800, 700], 'Name', 'Mixed Design Network');
        for i = 1:4
            subplot(2,2,i);
            topoplot([], chanlocs, 'style', 'blank', 'electrodes', 'on', 'headrad', 'rim');
            hold on;
            
            % 颜色
            val = means(i);
            n_v = (val - c_lim(1)) / (c_lim(2) - c_lim(1));
            n_v = max(0, min(1, n_v));
            col = cmap(round(n_v*255)+1, :);
            
            % 线型 (区分组内条件)
            % i=1,3 是 Cond1; i=2,4 是 Cond2
            if mod(i,2)==1, lst='-'; lw=3.5; else, lst='--'; lw=2; end
            
            plot([x1 x2], [y1 y2], 'Color', col, 'LineWidth', lw, 'LineStyle', lst);
            title(titles_plot{i}, 'FontSize', 12, 'FontWeight', 'bold');
            text(0, -0.7, sprintf('%.3f', val), 'Horiz', 'center', 'FontWeight','bold');
            axis equal off; hold off;
        end
        colormap(cmap); caxis(c_lim);
        cb = colorbar('Position', [0.93, 0.15, 0.02, 0.7]);
        ylabel(cb, 'FC Strength');
    end
    
    % --- B. 交互作用柱状图 ---
    figure('Color','w', 'Position', [200,200,500,600], 'Name', 'Bar Plot');
    % Data Format for Bar: Row=Groups, Col=Conditions
    bar_d = [means(1), means(2); means(3), means(4)];
    err_d = [ses(1), ses(2); ses(3), ses(4)];
    
    b = bar(bar_d, 'grouped'); hold on;
    b(1).FaceColor = [0.85, 0.33, 0.10]; % Cond 1
    b(2).FaceColor = [0.00, 0.45, 0.74]; % Cond 2
    
    % Errorbars
    [ngr, nba] = size(bar_d);
    x_pos = nan(nba, ngr);
    for i=1:nba, x_pos(i,:) = b(i).XEndPoints; end
    errorbar(x_pos', bar_d, err_d, 'k', 'linestyle', 'none', 'LineWidth', 1.5);
    
    set(gca, 'XTickLabel', {[Name_Btwn ' 1'], [Name_Btwn ' 2']}, 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Connectivity Strength');
    legend({[Name_Wthn ' 1'], [Name_Wthn ' 2']}, 'Location', 'best');
    title(['Interaction: ' Name_Btwn ' x ' Name_Wthn]);
    box off; hold off;
    
    disp('高级绘图完成。');
end

% =========================================================================
%% ================= 7. 导出至 BrainNet Viewer (3D 脑网络可视化) =================
choice_bnv = questdlg('是否导出数据到 BrainNet Viewer 进行 3D 绘图?', 'BrainNet Viewer Export', 'Yes', 'No', 'Yes');

if strcmp(choice_bnv, 'Yes')
    disp('>>> 开始生成 BrainNet Viewer 文件...');
    
    % 1. 检查必要信息
    if ~exist('chanlocs', 'var') || isempty(chanlocs)
        [fn, pn] = uigetfile('*.set', '请选择带电极坐标的 .set 文件');
        if isequal(fn,0), return; end
        EEG_tmp = pop_loadset('filename', fn, 'filepath', pn);
        chanlocs = EEG_tmp.chanlocs;
    end
    
    % 确保有计算好的均值数据 (G1C1, G1C2, G2C1, G2C2)
    % 如果没有，这里给个默认值测试用，或者基于上一部分的 means
    if ~exist('means', 'var')
        warning('未找到 means 变量，将使用随机数据进行演示导出！');
        means = rand(1, 4); 
    end
    
    % 2. 创建保存文件夹
    save_dir = uigetdir(pwd, '选择保存 .node 和 .edge 文件的文件夹');
    if isequal(save_dir, 0), return; end
    
    % ---------------------------------------------------------
    % A. 生成 .node 文件 (节点信息：坐标、颜色、大小、标签)
    % ---------------------------------------------------------
    nChans = length(chanlocs);
    node_data = zeros(nChans, 6); % [X, Y, Z, Color, Size, Label(text)]
    
    % 尝试获取坐标并转换为 MNI 近似空间
    % EEG 坐标通常在单位球上，BrainNet Viewer 需要 MNI (mm, range approx -90 to 90)
    % 我们将坐标放大以便它包围住大脑模型
    scale_factor = 85; 
    
    for i = 1:nChans
        % 获取/转换坐标
        if isfield(chanlocs, 'X') && ~isempty(chanlocs(i).X)
            x = chanlocs(i).X; y = chanlocs(i).Y; z = chanlocs(i).Z;
        elseif isfield(chanlocs, 'theta')
            % 极坐标转笛卡尔
            [y, x, z] = pol2cart(chanlocs(i).theta * pi/180, chanlocs(i).radius);
            % pol2cart Z 是0，我们需要把 radius 投影到球面上
            % 简单的球体投影:
            phi = (chanlocs(i).radius) * pi/2; % 假设 max radius = 1 对应 90度
            x = sin(phi) * cos(chanlocs(i).theta * pi/180);
            y = sin(phi) * sin(chanlocs(i).theta * pi/180);
            z = cos(phi);
        else
            x=0; y=0; z=0;
        end
        
        % 归一化并放大到 MNI 尺寸
        vec = [x, y, z];
        if norm(vec) > 0
            vec = vec / norm(vec) * scale_factor;
        end
        
        % 存入矩阵
        node_data(i, 1) = vec(1); % X
        node_data(i, 2) = vec(2); % Y
        node_data(i, 3) = vec(3); % Z
        
        % 设置颜色模块 (例如: 左半球=1, 右半球=2, 中线=3)
        if vec(1) < -5, col = 1;      % Left
        elseif vec(1) > 5, col = 2;   % Right
        else, col = 3; end            % Mid
        node_data(i, 4) = col; 
        
        % 设置节点大小 (默认小一点，关键电极可以大一点)
        node_data(i, 5) = 2; 
    end
    
    % 将感兴趣的两个电极标记为大红色
    ch1 = chan_pair(1); ch2 = chan_pair(2);
    node_data(ch1, 4) = 4; node_data(ch1, 5) = 5; % Color 4 (Red), Size 5
    node_data(ch2, 4) = 4; node_data(ch2, 5) = 5;
    
    % 写入 .node 文件
    node_file = fullfile(save_dir, 'EEG_Electrodes.node');
    fid = fopen(node_file, 'w');
    for i = 1:nChans
        fprintf(fid, '%.4f %.4f %.4f %d %d %s\n', ...
            node_data(i,1), node_data(i,2), node_data(i,3), ...
            node_data(i,4), node_data(i,5), chanlocs(i).labels);
    end
    fclose(fid);
    fprintf('已生成节点文件: %s\n', node_file);
    
    % ---------------------------------------------------------
    % B. 生成 .edge 文件 (连接矩阵)
    % ---------------------------------------------------------
    % 我们有4个条件，生成4个 edge 文件
    file_suffixes = {'_G1_Cond1', '_G1_Cond2', '_G2_Cond1', '_G2_Cond2'};
    
    for k = 1:4
        edge_matrix = zeros(nChans, nChans);
        
        % 填入我们的 FC 强度
        val = means(k);
        
        % BrainNet Viewer 读取下三角或全矩阵
        % 将感兴趣的电极对连线赋值
        edge_matrix(ch1, ch2) = val;
        edge_matrix(ch2, ch1) = val;
        
        % 保存 .edge 文件
        edge_fname = ['Connectivity' file_suffixes{k} '.edge'];
        dlmwrite(fullfile(save_dir, edge_fname), edge_matrix, 'delimiter', '\t');
        fprintf('已生成连边文件: %s (Strength: %.4f)\n', edge_fname, val);
    end
    
    msgbox({'BrainNet Viewer 文件生成完毕！'; ...
            ''; ...
            '使用步骤:'; ...
            '1. 打开 BrainNet Viewer'; ...
            '2. Load File -> Surface: 选择 "ICBM152.nv" (软件自带)'; ...
            '3. Load File -> Node: 选择刚才生成的 "EEG_Electrodes.node"'; ...
            '4. Load File -> Edge: 选择任意一个生成的 ".edge" 文件'; ...
            '5. 点击 OK，即可看到 3D 效果！'}, ...
            '导出成功');
end