clc; clear all; close all

%% 1. 指定相关信息 (保留原有交互风格)
% -------------------------------------------------------------------------
% 弹出对话框输入参数
prompt = { ...
    'The channel pair to be test (e.g., [6 5] for lower triangle)', ...
    'The time axis of TFR (in ms, e.g., -1000:4:1996)', ...
    'The frequency axis of TFR (in Hz, e.g., 1:30)', ...
    'The baseline limits (in ms, e.g., [-800 -200])', ...
    'Name of Factor A (e.g., Group/Condition)', ...
    'Name of Factor B (e.g., Emotion/Stimuli)' ...
};
dlgtitle = 'Parameters for 2x2 ANOVA FC Analysis';
dims = [1 60];
definput = {'[6 5]', '-1000:4:1996', '1:30', '[-800 -200]', 'FactorA', 'FactorB'};
answer = inputdlg(prompt, dlgtitle, dims, definput);

if isempty(answer), return; end % 如果点击取消则停止

% 解析输入参数
channel_pair = str2num(answer{1});
time_axis = str2num(answer{2});
f_axis = str2num(answer{3});
baseline = str2num(answer{4});
name_Fac1 = answer{5};
name_Fac2 = answer{6};

%% 2. 载入数据 (依次选择4个条件的文件)
% -------------------------------------------------------------------------
% 逻辑：2x2 设计需要4个文件。顺序至关重要！
% Order: 
% 1. FactorA_Level1 & FactorB_Level1
% 2. FactorA_Level1 & FactorB_Level2
% 3. FactorA_Level2 & FactorB_Level1
% 4. FactorA_Level2 & FactorB_Level2

titles = {
    ['1. Select Mat File: ' name_Fac1 ' (L1) & ' name_Fac2 ' (L1)'], ...
    ['2. Select Mat File: ' name_Fac1 ' (L1) & ' name_Fac2 ' (L2)'], ...
    ['3. Select Mat File: ' name_Fac1 ' (L2) & ' name_Fac2 ' (L1)'], ...
    ['4. Select Mat File: ' name_Fac1 ' (L2) & ' name_Fac2 ' (L2)']
};

% 初始化大矩阵 [Freq, Time, Subj, 4]
% 先读取第一个文件来确定被试数量和维度
disp('>>> 请按照弹出窗口提示，依次选择4个条件的 .mat 文件 <<<');
[filename1, pathname1] = uigetfile('*.mat', titles{1});
if isequal(filename1,0), return; end

% 预读取第一个文件，获取维度信息
file_path = fullfile(pathname1, filename1);
temp_struct = load(file_path);
vars = fieldnames(temp_struct);
data_temp = temp_struct.(vars{1}); % 假设 mat 里只有一个主变量 (如 wpli)

[nFreq, nTime, nCh, ~, nSubj] = size(data_temp);

% 检查维度一致性
if nFreq ~= length(f_axis) || nTime ~= length(time_axis)
    error('错误：输入的时间或频率轴长度与数据文件维度不匹配！请检查输入参数。');
end

% 初始化存储矩阵
All_Data_ROI = zeros(nFreq, nTime, nSubj, 4);

% 循环读取4个文件并处理
file_list_paths = {file_path, '', '', ''}; % 存储路径以便循环

% 填充第一个文件的数据
fprintf('正在处理条件 1: %s ...\n', filename1);
% 提取 ROI [Freq, Time, Subj]
roi_data = squeeze(data_temp(:,:,channel_pair(1), channel_pair(2),:));
% 基线校正
base_idx = time_axis >= baseline(1) & time_axis <= baseline(2);
base_mean = mean(roi_data(:, base_idx, :), 2);
All_Data_ROI(:,:,:,1) = roi_data - repmat(base_mean, [1, nTime, 1]);
clear data_temp temp_struct roi_data base_mean

% 读取剩余3个文件
for c = 2:4
    [filename, pathname] = uigetfile('*.mat', titles{c});
    if isequal(filename,0), error('用户取消操作'); end
    
    fprintf('正在处理条件 %d: %s ...\n', c, filename);
    file_path = fullfile(pathname, filename);
    temp_struct = load(file_path);
    vars = fieldnames(temp_struct);
    data_raw = temp_struct.(vars{1}); % 获取 5D 矩阵
    
    % 检查被试数
    if size(data_raw, 5) ~= nSubj
        error('条件 %d 的被试数量与条件 1 不一致！', c);
    end
    
    % 提取 ROI
    roi_data = squeeze(data_raw(:,:,channel_pair(1), channel_pair(2),:));
    
    % 基线校正
    base_mean = mean(roi_data(:, base_idx, :), 2);
    All_Data_ROI(:,:,:,c) = roi_data - repmat(base_mean, [1, nTime, 1]);
    
    clear data_raw temp_struct roi_data
end

%% 3. 对逐个时频点进行 2x2 ANOVA，并使用 FDR 矫正
% -------------------------------------------------------------------------
fprintf('开始进行方差分析计算，请稍候...\n');

% 准备 ANOVA 参数 (对应 rm_anova2)
sub = 25;
F1_con = 2; F2_con = 2; total_con = 4;
FACTNAMES = {name_Fac1, name_Fac2};

% 构造因子向量
S = repmat([1:sub], 1, total_con).';
F1 = sort(repmat([1:F1_con], 1, sub*F2_con)).'; 
F2 = repmat(sort(repmat([1:F2_con], 1, sub)), 1, F1_con).';

% 初始化 P 值矩阵
P_MainA = ones(nFreq, nTime);
P_MainB = ones(nFreq, nTime);
P_Inter = ones(nFreq, nTime);

h_wait = waitbar(0, 'Running ANOVA...');

for f = 1:nFreq
    waitbar(f/nFreq, h_wait, sprintf('Frequency %d / %d', f, nFreq));
    for t = 1:nTime
        % 提取数据: [Subj, 4] -> 转为长向量 [Subj*4, 1]
        data_slice = squeeze(All_Data_ROI(f, t, :, :));
        Y = data_slice(:);
        
        try
            tab = rm_anova2(Y, S, F1, F2, FACTNAMES);
            % tab{2,6}: Factor A Main Effect
            % tab{3,6}: Factor B Main Effect
            % tab{4,6}: Interaction
            p_vec = [tab{2,6}, tab{3,6}, tab{4,6}]; 
            
            P_MainA(f,t) = p_vec(1);
            P_MainB(f,t) = p_vec(2);
            P_Inter(f,t) = p_vec(3);
        catch
            % Catch errors (e.g. constant data)
        end
    end
end
close(h_wait);

% FDR 校正
[~, Mask_MainA] = fdr(P_MainA, 0.05);
[~, Mask_MainB] = fdr(P_MainB, 0.05);
[~, Mask_Inter] = fdr(P_Inter, 0.05);

%% 4. 绘图 (直接使用循环绘图，解决子函数报错问题)
% -------------------------------------------------------------------------

% 将需要画的三组数据整理到 Cell 数组中，方便循环处理
plot_data = {
    struct('p', P_MainA, 'mask', Mask_MainA, 'title', ['Main Effect: ' name_Fac1]), ...
    struct('p', P_MainB, 'mask', Mask_MainB, 'title', ['Main Effect: ' name_Fac2]), ...
    struct('p', P_Inter, 'mask', Mask_Inter, 'title', ['Interaction: ' name_Fac1 ' x ' name_Fac2])
};

for i = 1:3
    % 获取当前循环的数据
    curr_p = plot_data{i}.p;
    curr_mask = plot_data{i}.mask;
    curr_title = plot_data{i}.title;
    
    figure('Color', 'w', 'Name', curr_title); % 新建一个窗口
    
    % --- 上图：原始 P 值 (未校正，只显示 p < 0.05) ---
    subplot(2,1,1); 
    p_plot_raw = curr_p;
    p_plot_raw(curr_p >= 0.05) = NaN; % 不显著的设为 NaN (透明)
    
    imagesc(time_axis, f_axis, p_plot_raw);
    axis xy; 
    colorbar;
    title([curr_title ' (Uncorrected p < 0.05)']);
    xlabel('Time (ms)'); ylabel('Frequency (Hz)');
    
    % 设置颜色映射：翻转 parula，让 p 值越小(越显著)颜色越暖/深
    colormap(gca, flipud(parula)); 
    caxis([0 0.05]); % 锁定颜色范围 0 到 0.05
    
    % --- 下图：FDR 校正后的显著区 ---
    subplot(2,1,2); 
    p_plot_fdr = curr_p;
    p_plot_fdr(curr_mask == 0) = NaN; % 只保留通过 FDR 校正的点
    
    if all(isnan(p_plot_fdr(:)))
        % 如果没有显著区域，在图中间写字
        text(mean(time_axis), mean(f_axis), 'No Significant Area (FDR)', ...
            'HorizontalAlignment', 'center', 'FontSize', 12);
        xlim([min(time_axis) max(time_axis)]); 
        ylim([min(f_axis) max(f_axis)]);
    else
        imagesc(time_axis, f_axis, p_plot_fdr);
        axis xy; 
        colorbar;
        caxis([0 0.05]);
    end
    
    title([curr_title ' (FDR Corrected)']);
    xlabel('Time (ms)'); ylabel('Frequency (Hz)');
    colormap(gca, flipud(parula));
end

disp('分析完成！已生成 3 张结果图。');


%% 5. 提取感兴趣区域 (ROI) 数据导出为 Excel (供 SPSS 使用)
% -------------------------------------------------------------------------
% 询问用户是否需要提取数据
choice = questdlg('分析完成。是否提取特定时频区域的平均值，保存为 Excel 供 SPSS 分析?', ...
    '导出数据', 'Yes', 'No', 'Yes');

if strcmp(choice, 'Yes')
    % 1. 弹出对话框输入感兴趣的范围 (ROI)
    % 建议参考刚刚生成的图，找到显著的区域范围
    prompt_roi = { ...
        'Enter Time Range (ms) [min max]:', ...
        'Enter Frequency Range (Hz) [min max]:'};
    dlg_title_roi = 'Define ROI for SPSS';
    def_roi = {'[200 400]', '[4 8]'}; % 默认值，可修改
    ans_roi = inputdlg(prompt_roi, dlg_title_roi, 1, def_roi);

    if ~isempty(ans_roi)
        t_range = str2num(ans_roi{1});
        f_range = str2num(ans_roi{2});

        % 2. 找到对应的时间和频率索引
        t_idx = time_axis >= t_range(1) & time_axis <= t_range(2);
        f_idx = f_axis >= f_range(1) & f_axis <= f_range(2);
        
        % 检查是否有有效的数据点
        if sum(t_idx) == 0 || sum(f_idx) == 0
            errordlg('选定的范围内没有数据点，请检查时间或频率范围输入是否正确。');
            return;
        end

        fprintf('正在提取 ROI 数据: Time [%d %d] ms, Freq [%d %d] Hz ...\n', ...
            t_range(1), t_range(2), f_range(1), f_range(2));

        % 3. 提取数据并计算平均值
        % All_Data_ROI 维度: [Freq, Time, Subj, 4]
        roi_block = All_Data_ROI(f_idx, t_idx, :, :);
        
        % 第一步：对频率维(dim 1)求平均
        temp_mean = mean(roi_block, 1, 'omitnan'); 
        % 第二步：对时间维(dim 2)求平均 -> 得到 [1, 1, Subj, 4]
        temp_mean = mean(temp_mean, 2, 'omitnan');
        
        % 压缩维度 -> 得到 [Subj, 4]
        Data_SPSS = squeeze(temp_mean);
        
        % 处理单被试的特殊情况 (防止维度转置错误)
        if nSubj == 1
            Data_SPSS = Data_SPSS(:)'; % 强制转为行向量 1x4
        end

        % 4. 准备保存为表格
        % 定义列名 (对应 SPSS 里的变量名)
        % 顺序严格对应你文件载入的顺序：
        % 1: A1B1, 2: A1B2, 3: A2B1, 4: A2B2
        VarNames = {'A1B1', 'A1B2', 'A2B1', 'A2B2'};
        
        % 创建 Table
        T = array2table(Data_SPSS, 'VariableNames', VarNames);
        
        % 增加一列被试编号
        SubjectID = (1:nSubj)';
        T = addvars(T, SubjectID, 'Before', 1);

        % 5. 保存文件
        [file_spss, path_spss] = uiputfile('*.xlsx', '保存 SPSS 数据文件');
        if file_spss ~= 0
            full_save_path = fullfile(path_spss, file_spss);
            writetable(T, full_save_path);
            
            msgbox({['数据已成功保存至: ' file_spss]; ...
                    ''; ...
                    'SPSS 变量对应顺序:'; ...
                    'A1B1: FactorA(1) - FactorB(1)'; ...
                    'A1B2: FactorA(1) - FactorB(2)'; ...
                    'A2B1: FactorA(2) - FactorB(1)'; ...
                    'A2B2: FactorA(2) - FactorB(2)'}, ...
                    '导出成功');
        end
    end
end


%% =========================================================================
%% =========================================================================
%  PART 4: 准备数据进行精美绘图 (2x2 网络图 + 柱状图)
% =========================================================================

% 1. 询问用户是否进行高级绘图
choice_plot = questdlg('ANOVA 分析完成。是否进行高级绘图 (2x2 网络拓扑图 + 柱状图)?', ...
    '高级绘图', 'Yes', 'No', 'Yes');

if strcmp(choice_plot, 'Yes')
    
    % ---------------------------------------------------------------------
    % A. 获取电极位置信息 (用于画头皮图)
    % ---------------------------------------------------------------------
    disp('>>> 为了绘制拓扑图，请选择一个包含电极位置信息(chanlocs)的 .set 文件 <<<');
    [name_set, path_set] = uigetfile('*.set', '选择 .set 文件');
    if isequal(name_set, 0)
        warning('未选择 .set 文件，无法绘制拓扑图，仅显示柱状图。');
        chanlocs = [];
    else
        EEG_temp = pop_loadset('filename', name_set, 'filepath', path_set);
        chanlocs = EEG_temp.chanlocs;
        clear EEG_temp;
    end
    
    % ---------------------------------------------------------------------
    % B. 确定 ROI 范围 (用于提取数值画图)
    % ---------------------------------------------------------------------
    % 默认使用之前输入的时间/频率轴，但通常我们需要一个特定的显著区域来画拓扑图
    prompt_roi = { ...
        'Enter Time Range for Plot (ms) [min max]:', ...
        'Enter Frequency Range for Plot (Hz) [min max]:'};
    dlg_title_roi = 'Define ROI for Plotting';
    
    % 尝试根据之前的输入给个默认值
    def_roi = {sprintf('[%d %d]', time_axis(1), time_axis(end)), ...
               sprintf('[%d %d]', f_axis(1), f_axis(end))}; 
    
    ans_roi = inputdlg(prompt_roi, dlg_title_roi, 1, def_roi);
    
    if isempty(ans_roi), return; end
    
    t_roi = str2num(ans_roi{1});
    f_roi = str2num(ans_roi{2});
    
    % ---------------------------------------------------------------------
    % C. 提取数据 (计算均值和标准误)
    % ---------------------------------------------------------------------
    % 找到索引
    t_idx = time_axis >= t_roi(1) & time_axis <= t_roi(2);
    f_idx = f_axis >= f_roi(1) & f_axis <= f_roi(2);
    
    % 初始化
    mean_values = zeros(1, 4);
    se_values   = zeros(1, 4);
    num_sub = nSubj;
    
    % 循环提取4个条件的数据
    for i = 1:4
        % All_Data_ROI 维度: [Freq, Time, Subj, 4]
        % 提取 ROI 内的数据 -> [Freq_pts, Time_pts, Subj]
        temp_block = All_Data_ROI(f_idx, t_idx, :, i);
        
        % 平均频率和时间 -> [1, 1, Subj]
        sub_vals = squeeze(mean(mean(temp_block, 1, 'omitnan'), 2, 'omitnan'));
        
        % 存均值和标准误
        mean_values(i) = mean(sub_vals);
        se_values(i)   = std(sub_vals) / sqrt(num_sub);
    end
    
    % ---------------------------------------------------------------------
    % D. 绘图1：2x2 脑网络图 (基于您提供的逻辑)
    % ---------------------------------------------------------------------
    if ~isempty(chanlocs)
        % 定义绘图用的标题
        plot_titles = {
            [name_Fac1 '(L1) - ' name_Fac2 '(L1)'], ...
            [name_Fac1 '(L1) - ' name_Fac2 '(L2)'], ...
            [name_Fac1 '(L2) - ' name_Fac2 '(L1)'], ...
            [name_Fac1 '(L2) - ' name_Fac2 '(L2)']
        };
        
        % 1. 获取电极坐标
        % 注意：channel_pair 是您最开始输入的 [ch1 ch2]
        ch1_idx = channel_pair(1);
        ch2_idx = channel_pair(2);
        
        % 为了获取准确坐标，画个隐藏图
        fig_temp = figure('Visible', 'off'); 
        topoplot([], chanlocs, 'style', 'blank', 'electrodes', 'on', 'headrad', 'rim');
        h_points = findobj(gca, 'Type', 'line', 'Marker', '.');
        if isempty(h_points), h_points = findobj(gca, 'Type', 'scatter'); end
        
        try
            x_all = h_points.XData; y_all = h_points.YData;
            % 注意：topoplot画点的顺序可能和chanlocs一致，但也可能倒序，这里假设一致
            % 如果坐标错乱，建议使用 pol2cart 方法作为备用
            x1 = x_all(ch1_idx); y1 = y_all(ch1_idx);
            x2 = x_all(ch2_idx); y2 = y_all(ch2_idx);
        catch
            % 备用方案：极坐标转笛卡尔
            th = [chanlocs.theta]; rd = [chanlocs.radius];
            [y_try, x_try] = pol2cart(th * pi / 180, rd); 
            x1 = y_try(ch1_idx); y1 = x_try(ch1_idx); 
            x2 = y_try(ch2_idx); y2 = x_try(ch2_idx);
        end
        close(fig_temp);
        
        % 2. 设定颜色参数
        max_val = max(abs(mean_values));
        if max_val == 0, max_val = 1; end
        c_lim = [-max_val*1.2, max_val*1.2]; 
        n_colors = 256;
        cmap = jet(n_colors); 
        
        figure('Color', 'w', 'Position', [100, 100, 800, 700], 'Name', '2x2 Network Topoplot');
        
        for i = 1:4
            ax = subplot(2, 2, i);
            
            % 画底图
            topoplot([], chanlocs, 'style', 'blank', 'electrodes', 'on', ...
                     'headrad', 'rim', 'shading', 'interp', 'whitebk', 'on');
            hold on;
            
            % 计算线条颜色
            val = mean_values(i);
            norm_val = (val - c_lim(1)) / (c_lim(2) - c_lim(1));
            norm_val = max(0, min(1, norm_val));
            color_idx = round(norm_val * (n_colors - 1)) + 1;
            this_rgb = cmap(color_idx, :);
            
            % 设定线型 (区分 Factor B)
            % 假设 1,3 是 FactorB Level 1; 2,4 是 FactorB Level 2
            if mod(i, 2) == 1 
                l_style = '-'; l_width = 3.5; % 实线
            else
                l_style = '--'; l_width = 2;  % 虚线
            end
            
            % 画线
            plot([x1 x2], [y1 y2], 'Color', this_rgb, 'LineWidth', l_width, 'LineStyle', l_style);
            
            axis equal; axis tight; axis off;
            title(plot_titles{i}, 'FontSize', 12, 'FontWeight', 'bold');
            text(0, -0.7, sprintf('%.3f', val), 'HorizontalAlignment', 'center', ...
                 'FontSize', 11, 'FontWeight', 'bold');
            hold off;
        end
        
        % Colorbar
        hp = colormap(cmap);
        cb = colorbar('Position', [0.93, 0.15, 0.02, 0.7]); 
        ylabel(cb, 'FC Strength', 'FontSize', 12);
        caxis(c_lim);
    end
    
    % ---------------------------------------------------------------------
    % E. 绘图2：分组柱状图 (Bar Chart)
    % ---------------------------------------------------------------------
    figure('Color', 'w', 'Position', [150, 150, 500, 600], 'Name', 'Interaction Bar Plot');
    
    % 整理数据: 行=FactorA(Level 1,2), 列=FactorB(Level 1,2)
    % Data Order in mean_values: A1B1, A1B2, A2B1, A2B2
    bar_data = [mean_values(1), mean_values(2); mean_values(3), mean_values(4)];
    err_data = [se_values(1), se_values(2); se_values(3), se_values(4)];
    
    % 绘图
    b = bar(bar_data, 'grouped');
    hold on;
    
    % 设置颜色
    b(1).FaceColor = [0.85, 0.33, 0.10]; % B Level 1 (橙)
    b(2).FaceColor = [0.00, 0.45, 0.74]; % B Level 2 (蓝)
    
    % 绘制误差棒
    [ngroups, nbars] = size(bar_data);
    x_pos = nan(nbars, ngroups);
    for i = 1:nbars
        x_pos(i,:) = b(i).XEndPoints;
    end
    errorbar(x_pos', bar_data, err_data, 'k', 'linestyle', 'none', 'LineWidth', 1.5);
    
    % 美化
    set(gca, 'XTickLabel', {[name_Fac1 ' L1'], [name_Fac1 ' L2']}, 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Connectivity Strength', 'FontSize', 12);
    legend({[name_Fac2 ' L1'], [name_Fac2 ' L2']}, 'Location', 'best');
    title(['Interaction: ' name_Fac1 ' x ' name_Fac2], 'FontSize', 14);
    
    box off;
    hold off;
    
    disp('高级绘图完成！');
end


%% =========================================================================
%  PART 5: 导出 BrainNet Viewer 3D 绘图文件 (针对特定电极对)
% =========================================================================
% 逻辑：生成 .node 文件(突出显示ROI电极) 和 4个 .edge 文件(仅包含该连线)

choice_bnv = questdlg('是否导出数据到 BrainNet Viewer 进行 3D 可视化?', ...
    '3D Export', 'Yes', 'No', 'Yes');

if strcmp(choice_bnv, 'Yes')
    disp('>>> 开始生成 3D 脑图文件...');
    
    % 1. 检查必要变量
    if ~exist('chanlocs', 'var') || isempty(chanlocs)
        disp('>>> 请选择一个 *预处理后* (已剔除电极) 的 .set 文件以获取坐标 <<<');
        [name_set, path_set] = uigetfile('*.set', '选择 .set 文件');
        if isequal(name_set, 0), return; end
        EEG_temp = pop_loadset('filename', name_set, 'filepath', path_set);
        chanlocs = EEG_temp.chanlocs;
    end
    
    if ~exist('mean_values', 'var')
        errordlg('未找到均值数据(mean_values)，请先运行上方的高级绘图部分(Part 4)计算数值。');
        return;
    end

    % 2. 选择保存路径
    save_dir = uigetdir(pwd, '选择保存 .node 和 .edge 文件的文件夹');
    if isequal(save_dir, 0), return; end
    
    % ---------------------------------------------------------------------
    % A. 生成 .node 文件 (高亮显示选中的电极对)
    % ---------------------------------------------------------------------
    nChans = length(chanlocs);
    node_data = zeros(nChans, 6); % [X, Y, Z, Color, Size, Label]
    
    % 坐标转换参数
    scale_factor = 85; % 放大以适配 BrainNet 模板
    
    % 获取目标电极索引
    ch1 = channel_pair(1);
    ch2 = channel_pair(2);
    
    for i = 1:nChans
        % 坐标转换: 极坐标 -> 笛卡尔 -> 投影到球体
        if isfield(chanlocs, 'X') && ~isempty(chanlocs(i).X)
            x = chanlocs(i).X; y = chanlocs(i).Y; z = chanlocs(i).Z;
        else
            [y, x, z] = pol2cart(chanlocs(i).theta * pi/180, chanlocs(i).radius);
            z = cos(chanlocs(i).radius * pi/2); % 简单球面投影
        end
        
        % 归一化并放大
        vec = [x, y, z];
        if norm(vec) > 0, vec = vec / norm(vec) * scale_factor; end
        
        node_data(i, 1:3) = vec;
        
        % --- 关键：设置颜色和大小 ---
        if i == ch1 || i == ch2
            node_data(i, 4) = 1; % 颜色 1 (红色，需在BNV中设置)
            node_data(i, 5) = 5; % 大小 5 (大球)
        else
            node_data(i, 4) = 2; % 颜色 2 (灰色/蓝色)
            node_data(i, 5) = 2; % 大小 2 (小球)
        end
    end
    
    % 写入 .node 文件
    node_filename = fullfile(save_dir, 'ROI_Electrodes.node');
    fid = fopen(node_filename, 'w');
    for i = 1:nChans
        % 最后一列写 Label，如果 Label 是空的则写 number
        lbl = chanlocs(i).labels;
        if isempty(lbl), lbl = num2str(i); end
        
        fprintf(fid, '%.4f %.4f %.4f %d %d %s\n', ...
            node_data(i,1), node_data(i,2), node_data(i,3), ...
            node_data(i,4), node_data(i,5), lbl);
    end
    fclose(fid);
    
    % ---------------------------------------------------------------------
    % B. 生成 4 个 .edge 文件 (对应 4 个条件)
    % ---------------------------------------------------------------------
    % 定义文件名后缀
    suffixes = {
        ['_' name_Fac1 '_L1_' name_Fac2 '_L1'], ...
        ['_' name_Fac1 '_L1_' name_Fac2 '_L2'], ...
        ['_' name_Fac1 '_L2_' name_Fac2 '_L1'], ...
        ['_' name_Fac1 '_L2_' name_Fac2 '_L2']
    };

    for k = 1:4
        edge_mat = zeros(nChans, nChans);
        
        % 填入连接强度 (mean_values 来自 Part 4)
        val = mean_values(k);
        
        % BrainNet Viewer 矩阵是对称的
        edge_mat(ch1, ch2) = val;
        edge_mat(ch2, ch1) = val;
        
        % 保存
        edge_fname = fullfile(save_dir, ['FC_Strength' suffixes{k} '.edge']);
        dlmwrite(edge_fname, edge_mat, 'delimiter', '\t');
        
        fprintf('已生成: %s (强度: %.4f)\n', ['FC_Strength' suffixes{k} '.edge'], val);
    end
    
    % ---------------------------------------------------------------------
    % C. 完成提示
    % ---------------------------------------------------------------------
    msgbox({'BrainNet Viewer 文件生成成功！'; ...
            ''; ...
            '请打开 BrainNet Viewer:'; ...
            '1. Surface: Load "ICBM152.nv"'; ...
            '2. Node: Load "ROI_Electrodes.node"'; ...
            '3. Edge: Load 任意一个生成的 .edge 文件'; ...
            ''; ...
            '提示: 在 Options -> Node 中设置 Color 1 为红色，Color 2 为灰色，'; ...
            '即可突出显示您分析的电极对。'}, ...
            '完成');
end