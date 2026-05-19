% =========================================================================
%  全脑全局连接时频扫描 (Global Connectivity TF Scanner)
%  目的: 快速找出有显著差异的时间和频率窗口 (ROI)
%  方法: 将所有电极对平均 -> 得到全局指标 -> 跑 2x2 Mixed ANOVA
% =========================================================================

clc; clear all; close all;

%% 1. 参数设置
prompt = {'时间轴 (ms) [Start Step End]:', '频率轴 (Hz) [Start Step End]:', ...
          '基线范围 (ms) [Min Max]:', ...
          '组间因素 (Group):', '组内因素 (Condition):'};
def = {'[-1000 4 1996]', '[1 1 30]', '[-800 -200]', 'Group', 'Condition'};
ans_set = inputdlg(prompt, 'Global Scan Params', 1, def);
if isempty(ans_set), return; end

time_axis = str2num(ans_set{1});
f_axis    = str2num(ans_set{2});
baseline  = str2num(ans_set{3});
Name_Btwn = ans_set{4};
Name_Wthn = ans_set{5};

%% 2. 载入数据 (Group 1 & Group 2)
% 逻辑同混合设计：G1C1, G1C2, G2C1, G2C2
disp('>>> 请依次选择 4 个文件 (包含 5D 连接矩阵) <<<');
titles = {'G1_Cond1', 'G1_Cond2', 'G2_Cond1', 'G2_Cond2'};

% 容器: 存储 [Freq, Time, Subj] 的全局平均数据
Global_G1 = []; 
Global_G2 = [];

% --- Group 1 ---
for c = 1:2
    fprintf('Loading Group 1 - %s ...\n', titles{c});
    [fn, pn] = uigetfile('*.mat', titles{c});
    tmp = load(fullfile(pn, fn)); vars = fieldnames(tmp);
    data_5d = tmp.(vars{1}); % [F, T, Ch, Ch, Subj]
    
    [nF, nT, ~, ~, nS] = size(data_5d);
    
    % !!! 核心步骤：对通道维(3,4)求平均，得到全局指标 !!!
    % mean(mean(data, 3), 4) -> [F, T, 1, 1, Subj]
    global_avg = squeeze(mean(mean(data_5d, 3, 'omitnan'), 4, 'omitnan')); 
    
    % 基线校正
    base_idx = time_axis >= baseline(1) & time_axis <= baseline(2);
    base_val = mean(global_avg(:, base_idx, :), 2);
    global_avg = global_avg - repmat(base_val, [1, nT, 1]);
    
    if c == 1
        nSubj_G1 = nS;
        Global_G1 = zeros(nF, nT, nSubj_G1, 2);
    end
    Global_G1(:,:,:,c) = global_avg;
    clear data_5d global_avg
end

% --- Group 2 ---
for c = 1:2
    fprintf('Loading Group 2 - %s ...\n', titles{c+2});
    [fn, pn] = uigetfile('*.mat', titles{c+2});
    tmp = load(fullfile(pn, fn)); vars = fieldnames(tmp);
    data_5d = tmp.(vars{1});
    
    [nF, nT, ~, ~, nS] = size(data_5d);
    
    global_avg = squeeze(mean(mean(data_5d, 3, 'omitnan'), 4, 'omitnan'));
    
    base_idx = time_axis >= baseline(1) & time_axis <= baseline(2);
    base_val = mean(global_avg(:, base_idx, :), 2);
    global_avg = global_avg - repmat(base_val, [1, nT, 1]);
    
    if c == 1
        nSubj_G2 = nS;
        Global_G2 = zeros(nF, nT, nSubj_G2, 2);
    end
    Global_G2(:,:,:,c) = global_avg;
    clear data_5d global_avg
end

%% 3. 快速扫描 (逐点 Mixed ANOVA)
fprintf('正在扫描全时频显著性 (Global Mean Connectivity)...\n');

P_Btwn = ones(nF, nT);
P_Wthn = ones(nF, nT);
P_Int  = ones(nF, nT);

h = waitbar(0, 'Scanning...');

% 预计算参数
df_grp = 1; 
df_err_b = (nSubj_G1 - 1) + (nSubj_G2 - 1);
df_cond = 1; 
df_err_w = df_err_b;

for f = 1:nF
    waitbar(f/nF, h);
    for t = 1:nT
        Y1 = squeeze(Global_G1(f, t, :, :)); % [Subj_G1, 2]
        Y2 = squeeze(Global_G2(f, t, :, :)); % [Subj_G2, 2]
        if nSubj_G1==1, Y1=Y1(:)'; end
        if nSubj_G2==1, Y2=Y2(:)'; end
        
        % 直接内嵌 ANOVA 计算 (避免函数调用错误)
        % 1. Means
        m_g1 = mean(Y1(:)); m_g2 = mean(Y2(:));
        gm = mean([Y1(:); Y2(:)]); % Grand Mean (Simplified for equal weight)
        % Note: For precise Mixed ANOVA with unequal N, we use Type III SS logic
        % Here simplified for quick scanning:
        
        % SS Between
        ss_grp = 2*nSubj_G1*(m_g1-gm)^2 + 2*nSubj_G2*(m_g2-gm)^2;
        
        subj_m_g1 = mean(Y1,2); subj_m_g2 = mean(Y2,2);
        ss_err_b = 2*sum((subj_m_g1-m_g1).^2) + 2*sum((subj_m_g2-m_g2).^2);
        
        % SS Within
        m_c1 = (sum(Y1(:,1)) + sum(Y2(:,1))) / (nSubj_G1+nSubj_G2);
        m_c2 = (sum(Y1(:,2)) + sum(Y2(:,2))) / (nSubj_G1+nSubj_G2);
        ss_cond = (nSubj_G1+nSubj_G2)*(m_c1-gm)^2 + (nSubj_G1+nSubj_G2)*(m_c2-gm)^2;
        
        % SS Interaction
        % cell means
        c11=mean(Y1(:,1)); c12=mean(Y1(:,2));
        c21=mean(Y2(:,1)); c22=mean(Y2(:,2));
        ss_cells = nSubj_G1*(c11-gm)^2 + nSubj_G1*(c12-gm)^2 + ...
                   nSubj_G2*(c21-gm)^2 + nSubj_G2*(c22-gm)^2;
        ss_int = ss_cells - ss_grp - ss_cond;
        
        % SS Error Within
        ss_w_g1 = sum(sum((Y1 - repmat(subj_m_g1,1,2)).^2));
        ss_w_g2 = sum(sum((Y2 - repmat(subj_m_g2,1,2)).^2));
        ss_err_w = (ss_w_g1 + ss_w_g2) - ss_cond - ss_int;
        
        % F & P
        f_grp = (ss_grp/1) / (ss_err_b/df_err_b);
        f_cond = (ss_cond/1) / (ss_err_w/df_err_w);
        f_int = (ss_int/1) / (ss_err_w/df_err_w);
        
        P_Btwn(f,t) = 1 - fcdf(f_grp, 1, df_err_b);
        P_Wthn(f,t) = 1 - fcdf(f_cond, 1, df_err_w);
        P_Int(f,t)  = 1 - fcdf(f_int, 1, df_err_w);
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

figure('Color','w', 'Position', [100, 100, 1200, 400], 'Name', 'Global Connectivity Scout');
for i = 1:3
    subplot(1,3,i);
    % 绘制 P 值，非显著区域透明
    p_plot = pvals{i};
    p_plot(masks{i} == 0) = NaN; % 只显示过 FDR 的区域
    
    if all(isnan(p_plot(:)))
        text(mean(time_axis), mean(f_axis), 'No Sig.', 'Horiz','center', 'FontSize',14);
        xlim([min(time_axis) max(time_axis)]); ylim([min(f_axis) max(f_axis)]);
    else
        imagesc(time_axis, f_axis, p_plot); 
        axis xy; colorbar;
        caxis([0 0.05]); % 锁定 0 到 0.05
        colormap(gca, flipud(parula));
    end
    title(plot_titles{i}, 'FontSize', 12, 'FontWeight','bold');
    xlabel('Time (ms)'); ylabel('Freq (Hz)');
end

msgbox('扫描完成！请根据图中的显著色块确定 ROI，然后填入全脑分析代码。');