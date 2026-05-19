% =========================================================================
%  全脑全局连接时频扫描 (Global Connectivity TF Scanner)
%  版本: 2x2 混合设计 (Mixed Design: 1 Between + 1 Within)
%  目的: 快速找出有显著差异的时间和频率窗口 (ROI)
%  原理: 全脑连接平均 -> 全局指标 -> 逐点跑 2x2 Mixed ANOVA
% =========================================================================

clc; clear all; close all;

%% 1. 参数设置
prompt = {'时间轴 (ms) [Start Step End]:', '频率轴 (Hz) [Start Step End]:', ...
          '基线范围 (ms) [Min Max]:', ...
          '组间因素名 (e.g. Group):', '组内因素名 (e.g. Condition):'};
def = {'[-1000 4 1996]', '[1 1 30]', '[-800 -200]', 'Group', 'Condition'};
ans_set = inputdlg(prompt, 'Global Scan Params (Mixed)', 1, def);
if isempty(ans_set), return; end

time_axis = str2num(ans_set{1});
f_axis    = str2num(ans_set{2});
baseline  = str2num(ans_set{3});
Name_Btwn = ans_set{4};
Name_Wthn = ans_set{5};

%% 2. 载入数据 (Group 1 & Group 2)
% 顺序约定: 
% 1. Group 1 - Cond 1
% 2. Group 1 - Cond 2
% 3. Group 2 - Cond 1
% 4. Group 2 - Cond 2

disp('>>> 请依次选择 4 个文件 (包含 5D 连接矩阵) <<<');
titles = {
    ['1. ' Name_Btwn '1 & ' Name_Wthn '1'], ...
    ['2. ' Name_Btwn '1 & ' Name_Wthn '2'], ...
    ['3. ' Name_Btwn '2 & ' Name_Wthn '1'], ...
    ['4. ' Name_Btwn '2 & ' Name_Wthn '2']
};

% 容器
Global_G1 = []; % [Freq, Time, Subj_G1, 2]
Global_G2 = []; % [Freq, Time, Subj_G2, 2]
nSubj_G1 = 0;
nSubj_G2 = 0;

% --- 载入 Group 1 ---
for c = 1:2
    fprintf('Loading Group 1 - %s ...\n', titles{c});
    [fn, pn] = uigetfile('*.mat', ['Select: ' titles{c}]);
    if isequal(fn,0), return; end
    
    tmp = load(fullfile(pn, fn)); vars = fieldnames(tmp);
    data_5d = tmp.(vars{1}); % [F, T, Ch, Ch, Subj]
    
    [nF, nT, ~, ~, nS] = size(data_5d);
    
    if c == 1
        nSubj_G1 = nS;
        Global_G1 = zeros(nF, nT, nSubj_G1, 2);
    elseif nS ~= nSubj_G1
        error('Group 1 的两个条件被试数不一致！');
    end
    
    % !!! 核心：全脑平均 !!!
    global_avg = squeeze(mean(mean(data_5d, 3, 'omitnan'), 4, 'omitnan')); 
    
    % 基线校正
    base_idx = time_axis >= baseline(1) & time_axis <= baseline(2);
    base_val = mean(global_avg(:, base_idx, :), 2);
    global_avg = global_avg - repmat(base_val, [1, nT, 1]);
    
    Global_G1(:,:,:,c) = global_avg;
    clear data_5d global_avg
end

% --- 载入 Group 2 ---
for c = 1:2
    fprintf('Loading Group 2 - %s ...\n', titles{c+2});
    [fn, pn] = uigetfile('*.mat', ['Select: ' titles{c+2}]);
    if isequal(fn,0), return; end
    
    tmp = load(fullfile(pn, fn)); vars = fieldnames(tmp);
    data_5d = tmp.(vars{1});
    
    [nF, nT, ~, ~, nS] = size(data_5d);
    
    if c == 1
        nSubj_G2 = nS;
        Global_G2 = zeros(nF, nT, nSubj_G2, 2);
    elseif nS ~= nSubj_G2
        error('Group 2 的两个条件被试数不一致！');
    end
    
    global_avg = squeeze(mean(mean(data_5d, 3, 'omitnan'), 4, 'omitnan'));
    
    base_idx = time_axis >= baseline(1) & time_axis <= baseline(2);
    base_val = mean(global_avg(:, base_idx, :), 2);
    global_avg = global_avg - repmat(base_val, [1, nT, 1]);
    
    Global_G2(:,:,:,c) = global_avg;
    clear data_5d global_avg
end

fprintf('数据载入完成。Group 1 (N=%d), Group 2 (N=%d)\n', nSubj_G1, nSubj_G2);

%% 3. 快速扫描 (逐点 Mixed ANOVA)
fprintf('正在扫描全时频显著性 (Global Mean Connectivity)...\n');

P_Btwn = ones(nF, nT); % Group Effect
P_Wthn = ones(nF, nT); % Condition Effect
P_Int  = ones(nF, nT); % Interaction

h = waitbar(0, 'Scanning (Mixed ANOVA)...');

% 自由度预计算
df_grp = 1; 
df_err_b = (nSubj_G1 - 1) + (nSubj_G2 - 1);
df_cond = 1; 
df_int = 1;
df_err_w = df_err_b; % 对于 2x2，组内误差自由度与组间误差自由度数值一样

for f = 1:nF
    waitbar(f/nF, h);
    for t = 1:nT
        % 提取数据 [Subj x 2]
        Y1 = squeeze(Global_G1(f, t, :, :)); 
        Y2 = squeeze(Global_G2(f, t, :, :));
        if nSubj_G1==1, Y1=Y1(:)'; end
        if nSubj_G2==1, Y2=Y2(:)'; end
        
        % --- 内嵌 2x2 Mixed ANOVA 计算 (Type III SS logic for unequal N) ---
        
        % 1. 计算均值
        m_g1 = mean(mean(Y1));     % Group 1 mean
        m_g2 = mean(mean(Y2));     % Group 2 mean
        grand_mean = (m_g1 + m_g2) / 2; % Unweighted Grand Mean
        
        % 2. SS Between (Group)
        % SS_A = n * sum((mean_group - grand_mean)^2) -> but n is unequal
        % For scanning, we use weighted sums for SS but unweighted for means usually
        % Let's use standard GLM approach logic simplified:
        SS_grp = 2*nSubj_G1*(m_g1 - grand_mean)^2 + 2*nSubj_G2*(m_g2 - grand_mean)^2;
        
        % SS Error Between (Subject within Group)
        subj_m_g1 = mean(Y1, 2); 
        subj_m_g2 = mean(Y2, 2);
        SS_err_b = 2 * sum((subj_m_g1 - m_g1).^2) + 2 * sum((subj_m_g2 - m_g2).^2);
        
        % 3. SS Within (Condition)
        m_c1_g1 = mean(Y1(:,1)); m_c2_g1 = mean(Y1(:,2));
        m_c1_g2 = mean(Y2(:,1)); m_c2_g2 = mean(Y2(:,2));
        
        m_c1 = (m_c1_g1 + m_c1_g2) / 2; % Unweighted Condition 1 Mean
        m_c2 = (m_c2_g1 + m_c2_g2) / 2; 
        
        % SS Condition calculated on the deviations
        % Simplified SS_cond = N_total * variance of cond means? No, use cell means.
        % Using the contrast method for speed and robustness (Logic of t-test on difference)
        
        % --- Method 2: Contrast Method (Equivalent to ANOVA F) ---
        % Much faster and safer for Unequal N
        
        % Interaction: Difference of Differences
        % D1 = C1-C2 for G1; D2 = C1-C2 for G2
        d_g1 = Y1(:,1) - Y1(:,2);
        d_g2 = Y2(:,1) - Y2(:,2);
        [~, p_int_val] = ttest2(d_g1, d_g2); % Independent t-test on difference scores
        
        % Within (Condition):
        % Test if the grand mean of difference scores is different from 0
        % This is complex with unequal N in ttest2.
        % Back to ANOVA Sums for Condition:
        % SS_cond = N * (m_c1 - GM)^2 ... let's stick to F statistic construction
        
        % Re-calculation of SS for Condition (Unweighted)
        SS_cond = (nSubj_G1+nSubj_G2) * ((m_c1 - grand_mean)^2 + (m_c2 - grand_mean)^2); 
        
        % SS Interaction
        % SS_cells for interaction
        SS_int = nSubj_G1*(m_c1_g1 - m_g1 - m_c1 + grand_mean)^2 + ...
                 nSubj_G1*(m_c2_g1 - m_g1 - m_c2 + grand_mean)^2 + ...
                 nSubj_G2*(m_c1_g2 - m_g2 - m_c1 + grand_mean)^2 + ...
                 nSubj_G2*(m_c2_g2 - m_g2 - m_c2 + grand_mean)^2;
        
        % SS Error Within
        % Sum of variances of difference scores / 2 roughly
        % Or: SS_total_within - SS_cond - SS_int (per group)
        ss_w_g1 = sum(sum((Y1 - repmat(subj_m_g1,1,2)).^2));
        ss_w_g2 = sum(sum((Y2 - repmat(subj_m_g2,1,2)).^2));
        SS_err_w = ss_w_g1 + ss_w_g2 - SS_int; % Note: SS_cond is orthogonal usually
        % Let's use the robust residual calculation
        SS_err_w = sum(sum((Y1 - repmat([m_c1_g1 m_c2_g1],nSubj_G1,1)).^2)) + ...
                   sum(sum((Y2 - repmat([m_c1_g2 m_c2_g2],nSubj_G2,1)).^2));
                   
        % F Statistics
        MS_grp = SS_grp / df_grp;
        MS_err_b = SS_err_b / df_err_b;
        MS_cond = SS_cond / df_cond;
        MS_int = SS_int / df_int;
        MS_err_w = SS_err_w / df_err_w;
        
        F_grp = MS_grp / MS_err_b;
        F_cond = MS_cond / MS_err_w;
        F_int = MS_int / MS_err_w;
        
        P_Btwn(f,t) = 1 - fcdf(F_grp, 1, df_err_b);
        P_Wthn(f,t) = 1 - fcdf(F_cond, 1, df_err_w);
        P_Int(f,t)  = 1 - fcdf(F_int, 1, df_err_w);
        
    end
end
close(h);

%% 4. 绘图 (FDR 校正后的显著区)
% FDR 校正
[~, Mask_Btwn] = fdr(P_Btwn, 0.05);
[~, Mask_Wthn] = fdr(P_Wthn, 0.05);
[~, Mask_Int]  = fdr(P_Int,  0.05);

plot_titles = {['Main Effect: ' Name_Btwn], ['Main Effect: ' Name_Wthn], 'Interaction'};
masks = {Mask_Btwn, Mask_Wthn, Mask_Int};
pvals = {P_Btwn, P_Wthn, P_Int};

figure('Color','w', 'Position', [100, 100, 1200, 400], 'Name', 'Global Connectivity Scout (Mixed)');
for i = 1:3
    subplot(1,3,i);
    
    p_plot = pvals{i};
    p_plot(masks{i} == 0) = NaN; % 只显示过 FDR 的区域
    
    if all(isnan(p_plot(:)))
        text(mean(time_axis), mean(f_axis), 'No Sig.', 'Horiz','center', 'FontSize',14);
        xlim([min(time_axis) max(time_axis)]); ylim([min(f_axis) max(f_axis)]);
    else
        imagesc(time_axis, f_axis, p_plot); 
        axis xy; colorbar;
        caxis([0 0.05]); 
        colormap(gca, flipud(parula));
    end
    
    title(plot_titles{i}, 'FontSize', 12, 'FontWeight','bold');
    xlabel('Time (ms)'); ylabel('Freq (Hz)');
end

msgbox({'扫描完成！(2x2 混合设计)'; ''; '请根据图中的色块确定 ROI'; '然后将数值填入全脑网络分析代码。'});