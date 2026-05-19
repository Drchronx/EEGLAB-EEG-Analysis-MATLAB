% =========================================================================
%  全脑连接组学分析 (Whole-brain Connectome Analysis)
%  版本: 2x2 混合设计 (Mixed Design: 1 Between + 1 Within)
%  适用: 两组不同被试 (Group 1 vs Group 2)，每人做两个条件 (Cond 1 vs Cond 2)
%  统计: 逐边混合方差分析 (Mass-univariate Mixed ANOVA) + FDR 校正
%  输出: BrainNet Viewer (.node 和 .edge)
% =========================================================================

clc; clear all; close all;

%% ================= 1. 参数设置与电极载入 =================
prompt = { ...
    '输入感兴趣频率范围 (Hz) [Min Max]:', ...
    '输入感兴趣时间范围 (ms) [Min Max]:', ...
    '时间轴 (ms) [Start Step End]:', ...
    '频率轴 (Hz) [Start Step End]:', ...
    '组间因素名 (e.g. Group):', ...
    '组内因素名 (e.g. Condition):' ...
};
% 默认参数 (请根据您的预处理修改)
def = {'[4 7]', '[300 500]', '[-1000 4 1996]', '[1 1 30]', 'Group', 'Condition'};
ans_set = inputdlg(prompt, '全脑混合设计参数', 1, def);
if isempty(ans_set), return; end

f_roi_range = str2num(ans_set{1});
t_roi_range = str2num(ans_set{2});
time_axis   = str2num(ans_set{3});
f_axis      = str2num(ans_set{4});
Name_Btwn   = ans_set{5};
Name_Wthn   = ans_set{6};

% 1.2 载入电极信息
disp('>>> [重要] 请选择一个 *预处理后* (已剔除电极) 的 .set 文件 <<<');
[fn, pn] = uigetfile('*.set', 'Load processed .set file');
if isequal(fn,0), return; end
EEG = pop_loadset('filename', fn, 'filepath', pn);
chanlocs = EEG.chanlocs;
nChans = length(chanlocs);
nPairs = nChans * (nChans - 1) / 2; 

fprintf('检测到 %d 个电极，共 %d 条连接。\n', nChans, nPairs);

%% ================= 2. 载入数据 (4个文件) =================
% 顺序约定: 
% 1. Group 1 - Cond 1
% 2. Group 1 - Cond 2
% 3. Group 2 - Cond 1
% 4. Group 2 - Cond 2

titles = {
    ['1. ' Name_Btwn '1 & ' Name_Wthn '1'], ...
    ['2. ' Name_Btwn '1 & ' Name_Wthn '2'], ...
    ['3. ' Name_Btwn '2 & ' Name_Wthn '1'], ...
    ['4. ' Name_Btwn '2 & ' Name_Wthn '2']
};

% 生成电极对索引
pairs_idx = nchoosek(1:nChans, 2); 

% 容器: 存储处理后的连边数据 [Subj x nPairs]
Data_G1_C1 = []; Data_G1_C2 = [];
Data_G2_C1 = []; Data_G2_C2 = [];
nSubj_G1 = 0; nSubj_G2 = 0;

for c = 1:4
    fprintf('正在读取文件 %d: %s ...\n', c, titles{c});
    [fn, pn] = uigetfile('*.mat', ['Select: ' titles{c}]);
    if isequal(fn,0), error('用户取消'); end
    
    tmp = load(fullfile(pn, fn)); vars = fieldnames(tmp);
    data_5d = tmp.(vars{1}); % [Freq, Time, Ch, Ch, Subj]
    
    [~, ~, nCh, ~, nS] = size(data_5d);
    
    if nCh ~= nChans
        error('数据通道数与 .set 文件不符！');
    end
    
    % 检查每组内被试数一致性
    if c == 1
        nSubj_G1 = nS;
    elseif c == 2 && nS ~= nSubj_G1
        error('Group 1 的两个条件被试数不一致！');
    elseif c == 3
        nSubj_G2 = nS; % Group 2 可以和 Group 1 人数不同
    elseif c == 4 && nS ~= nSubj_G2
        error('Group 2 的两个条件被试数不一致！');
    end
    
    % --- 提取 ROI 平均值 ---
    t_idx = time_axis >= t_roi_range(1) & time_axis <= t_roi_range(2);
    f_idx = f_axis >= f_roi_range(1) & f_axis <= f_roi_range(2);
    
    % 平均 Freq/Time -> [Ch, Ch, Subj]
    temp_avg = squeeze(mean(mean(data_5d(f_idx, t_idx, :, :, :), 1, 'omitnan'), 2, 'omitnan'));
    
    % 拉直为 [Subj x nPairs]
    pairs_data = zeros(nS, nPairs);
    for s = 1:nS
        curr_mat = temp_avg(:,:,s);
        for p = 1:nPairs
            ch1 = pairs_idx(p,1); ch2 = pairs_idx(p,2);
            pairs_data(s, p) = curr_mat(ch1, ch2);
        end
    end
    
    % 存入对应变量
    if c==1, Data_G1_C1 = pairs_data;
    elseif c==2, Data_G1_C2 = pairs_data;
    elseif c==3, Data_G2_C1 = pairs_data;
    elseif c==4, Data_G2_C2 = pairs_data;
    end
    
    clear data_5d temp_avg pairs_data
end

fprintf('数据载入完毕。Group 1 (N=%d), Group 2 (N=%d)\n', nSubj_G1, nSubj_G2);

%% ================= 3. 全脑 2x2 混合方差分析 =================
fprintf('开始全脑统计 (2x2 Mixed ANOVA)...\n');
h_bar = waitbar(0, '正在计算全脑连接...');

p_Btwn = ones(1, nPairs);
p_Wthn = ones(1, nPairs);
p_Int  = ones(1, nPairs);

% 预计算自由度
df_grp = 1; 
df_err_b = (nSubj_G1 - 1) + (nSubj_G2 - 1);
df_cond = 1; 
df_int = 1;
df_err_w = df_err_b; 

for p = 1:nPairs
    if mod(p, 500) == 0, waitbar(p/nPairs, h_bar); end
    
    % 提取当前连线的数据
    % Group 1: [N1 x 2] matrix
    Y1 = [Data_G1_C1(:, p), Data_G1_C2(:, p)];
    % Group 2: [N2 x 2] matrix
    Y2 = [Data_G2_C1(:, p), Data_G2_C2(:, p)];
    
    % --- 快速混合 ANOVA 计算 (基于 SS 公式) ---
    
    % 1. Means
    m_g1 = mean(Y1(:)); 
    m_g2 = mean(Y2(:));
    grand_mean = (sum(Y1(:)) + sum(Y2(:))) / (nSubj_G1*2 + nSubj_G2*2);
    
    % Type III SS Logic approximations for speed
    
    % --- SS Between (Group) ---
    % SS_grp = sum n_i * (mean_i - GM)^2
    SS_grp = 2*nSubj_G1*(m_g1 - grand_mean)^2 + 2*nSubj_G2*(m_g2 - grand_mean)^2;
    
    subj_m_g1 = mean(Y1, 2); 
    subj_m_g2 = mean(Y2, 2);
    SS_err_b = 2*sum((subj_m_g1 - m_g1).^2) + 2*sum((subj_m_g2 - m_g2).^2);
    
    % --- SS Within (Condition) ---
    m_c1 = (sum(Y1(:,1)) + sum(Y2(:,1))) / (nSubj_G1 + nSubj_G2);
    m_c2 = (sum(Y1(:,2)) + sum(Y2(:,2))) / (nSubj_G1 + nSubj_G2);
    SS_cond = (nSubj_G1+nSubj_G2) * ((m_c1 - grand_mean)^2 + (m_c2 - grand_mean)^2);
    
    % --- SS Interaction ---
    % Cell means
    c11 = mean(Y1(:,1)); c12 = mean(Y1(:,2));
    c21 = mean(Y2(:,1)); c22 = mean(Y2(:,2));
    
    SS_cells = nSubj_G1*(c11-grand_mean)^2 + nSubj_G1*(c12-grand_mean)^2 + ...
               nSubj_G2*(c21-grand_mean)^2 + nSubj_G2*(c22-grand_mean)^2;
    
    SS_int = SS_cells - SS_grp - SS_cond;
    
    % --- SS Error Within ---
    % Residuals
    SS_w_g1 = sum(sum((Y1 - repmat(subj_m_g1,1,2)).^2));
    SS_w_g2 = sum(sum((Y2 - repmat(subj_m_g2,1,2)).^2));
    SS_err_w = (SS_w_g1 + SS_w_g2) - SS_cond - SS_int; 
    
    % Stats
    F_grp = (SS_grp/df_grp) / (SS_err_b/df_err_b);
    F_cond = (SS_cond/df_cond) / (SS_err_w/df_err_w);
    F_int = (SS_int/df_int) / (SS_err_w/df_err_w);
    
    p_Btwn(p) = 1 - fcdf(F_grp, df_grp, df_err_b);
    p_Wthn(p) = 1 - fcdf(F_cond, df_cond, df_err_w);
    p_Int(p)  = 1 - fcdf(F_int, df_int, df_err_w);
end
close(h_bar);

%% ================= 4. FDR 校正 =================
fprintf('\n--- 统计结果 (FDR q < 0.05) ---\n');
alpha = 0.05;

[~, mask_Btwn] = fdr(p_Btwn, alpha);
[~, mask_Wthn] = fdr(p_Wthn, alpha);
[~, mask_Int]  = fdr(p_Int, alpha);

fprintf('组间主效应 (%s) 显著连接: %d\n', Name_Btwn, sum(mask_Btwn));
fprintf('组内主效应 (%s) 显著连接: %d\n', Name_Wthn, sum(mask_Wthn));
fprintf('交互作用显著连接: %d\n', sum(mask_Int));

if sum(mask_Btwn)+sum(mask_Wthn)+sum(mask_Int) == 0
    warning('未发现 FDR 校正后的显著连接。');
end

%% ================= 5. 导出 BrainNet Viewer 文件 (3D) =================
choice_exp = questdlg('是否导出 BrainNet Viewer 文件?', 'Export 3D', 'Yes', 'No', 'Yes');

if strcmp(choice_exp, 'Yes')
    save_dir = uigetdir(pwd, '选择保存结果的文件夹');
    if isequal(save_dir,0), return; end
    
    % --- A. 生成 Node 文件 ---
    node_file = fullfile(save_dir, 'Electrodes_Mixed.node');
    fid = fopen(node_file, 'w');
    scale = 85; 
    for i = 1:nChans
        if isfield(chanlocs, 'X') && ~isempty(chanlocs(i).X)
            x=chanlocs(i).X; y=chanlocs(i).Y; z=chanlocs(i).Z;
        else
            [y,x,z] = pol2cart(chanlocs(i).theta*pi/180, chanlocs(i).radius);
            z = cos(chanlocs(i).radius * pi/2); 
        end
        vec = [x, y, z]; 
        if norm(vec)>0, vec = vec/norm(vec)*scale; end
        fprintf(fid, '%.4f %.4f %.4f 1 2 %s\n', vec(1), vec(2), vec(3), chanlocs(i).labels);
    end
    fclose(fid);
    
    % --- B. 生成 Edge 文件 ---
    masks = {mask_Btwn, mask_Wthn, mask_Int};
    names = {['Main_' Name_Btwn], ['Main_' Name_Wthn], 'Interaction'};
    
    for k = 1:3
        curr_mask = masks{k};
        if sum(curr_mask) > 0
            edge_mat = zeros(nChans, nChans);
            sig_idx = find(curr_mask == 1);
            
            for s = 1:length(sig_idx)
                idx = sig_idx(s);
                ch1 = pairs_idx(idx, 1);
                ch2 = pairs_idx(idx, 2);
                edge_mat(ch1, ch2) = 1;
                edge_mat(ch2, ch1) = 1;
            end
            
            fname = fullfile(save_dir, [names{k} '_FDR05.edge']);
            dlmwrite(fname, edge_mat, 'delimiter', '\t');
            fprintf('已导出: %s\n', fname);
        end
    end
    msgbox('导出完成！请使用 BrainNet Viewer 加载生成的文件。');
end