% =========================================================================
%  全脑全局连接时频扫描 (Global Connectivity TF Scanner)
%  版本: 2x2 纯组内设计 (Repeated Measures)
%  目的: 快速找出有显著差异的时间和频率窗口 (ROI)
%  方法: 将所有电极对平均 -> 得到全局指标 -> 跑 2x2 RM-ANOVA (快速对比法)
% =========================================================================

clc; clear all; close all;

%% 1. 参数设置
prompt = {'时间轴 (ms) [Start Step End]:', '频率轴 (Hz) [Start Step End]:', ...
          '基线范围 (ms) [Min Max]:', ...
          '因素A名称 (e.g. Load):', '因素B名称 (e.g. Stimulus):'};
def = {'[-1000 4 1996]', '[1 1 30]', '[-800 -200]', 'Load', 'Stimulus'};
ans_set = inputdlg(prompt, 'Global Scan Params (Within)', 1, def);
if isempty(ans_set), return; end

time_axis = str2num(ans_set{1});
f_axis    = str2num(ans_set{2});
baseline  = str2num(ans_set{3});
Name_A    = ans_set{4};
Name_B    = ans_set{5};

%% 2. 载入数据 (4个条件，同一批被试)
% 顺序约定: 
% 1: A1 B1
% 2: A1 B2
% 3: A2 B1
% 4: A2 B2

disp('>>> 请依次选择 4 个文件 (包含 5D 连接矩阵) <<<');
titles = {
    ['1. ' Name_A '1 & ' Name_B '1'], ...
    ['2. ' Name_A '1 & ' Name_B '2'], ...
    ['3. ' Name_A '2 & ' Name_B '1'], ...
    ['4. ' Name_A '2 & ' Name_B '2']
};

% 容器: 存储 [Freq, Time, Subj, 4] 的全局平均数据
Global_Data = []; 
nSubj = 0;

for c = 1:4
    fprintf('Loading Condition %d: %s ...\n', c, titles{c});
    [fn, pn] = uigetfile('*.mat', ['Select: ' titles{c}]);
    if isequal(fn,0), return; end
    
    tmp = load(fullfile(pn, fn)); 
    vars = fieldnames(tmp);
    data_5d = tmp.(vars{1}); % [F, T, Ch, Ch, Subj]
    
    [nF, nT, ~, ~, nS] = size(data_5d);
    
    % 检查被试数一致性
    if c == 1
        nSubj = nS;
        Global_Data = zeros(nF, nT, nSubj, 4);
    elseif nS ~= nSubj
        error('错误：条件 %d 的被试数 (%d) 与条件 1 (%d) 不一致！', c, nS, nSubj);
    end
    
    % !!! 核心步骤：对通道维(3,4)求平均，得到全局指标 !!!
    % mean(mean(data, 3), 4) -> [F, T, 1, 1, Subj]
    % 这一步把全脑所有连接压缩成一个值，代表“大脑整体连接强度”
    global_avg = squeeze(mean(mean(data_5d, 3, 'omitnan'), 4, 'omitnan')); 
    
    % 基线校正
    base_idx = time_axis >= baseline(1) & time_axis <= baseline(2);
    base_val = mean(global_avg(:, base_idx, :), 2);
    global_avg = global_avg - repmat(base_val, [1, nT, 1]);
    
    % 存入矩阵
    Global_Data(:,:,:,c) = global_avg;
    clear data_5d global_avg
end

fprintf('数据载入完成。N = %d\n', nSubj);

%% 3. 快速扫描 (逐点 2x2 RM-ANOVA)
fprintf('正在扫描全时频显著性 (基于快速对比法)...\n');

P_MainA = ones(nF, nT);
P_MainB = ones(nF, nT);
P_Int   = ones(nF, nT);

h = waitbar(0, 'Scanning...');

for f = 1:nF
    waitbar(f/nF, h);
    for t = 1:nT
        % 提取当前点所有被试、所有条件的数据 [Subj x 4]
        Y = squeeze(Global_Data(f, t, :, :)); 
        % 如果只有1个被试，Y会变成向量，需转置确保是行向量 (虽然1个被试没法跑统计)
        if nSubj == 1, Y = Y(:)'; end
        
        y1 = Y(:,1); % A1B1
        y2 = Y(:,2); % A1B2
        y3 = Y(:,3); % A2B1
        y4 = Y(:,4); % A2B2
        
        % --- 统计核心：利用 t 检验计算 ANOVA 效应 ---
        
        % 1. 主效应 A (A1 vs A2)
        % 比较 (A1B1+A1B2) - (A2B1+A2B2) 是否 != 0
        diff_A = (y1 + y2) - (y3 + y4);
        [~, P_MainA(f,t)] = ttest(diff_A);
        
        % 2. 主效应 B (B1 vs B2)
        % 比较 (A1B1+A2B1) - (A1B2+A2B2) 是否 != 0
        diff_B = (y1 + y3) - (y2 + y4);
        [~, P_MainB(f,t)] = ttest(diff_B);
        
        % 3. 交互作用 (Difference of Differences)
        % 比较 (A1B1-A1B2) - (A2B1-A2B2) 是否 != 0
        diff_Int = (y1 - y2) - (y3 - y4);
        [~, P_Int(f,t)] = ttest(diff_Int);
    end
end
close(h);

%% 4. 绘图 (FDR 校正后的显著区)
% FDR 校正
[~, Mask_A]   = fdr(P_MainA, 0.05);
[~, Mask_B]   = fdr(P_MainB, 0.05);
[~, Mask_Int] = fdr(P_Int,   0.05);

plot_titles = {['Main Effect: ' Name_A], ['Main Effect: ' Name_B], 'Interaction'};
masks = {Mask_A, Mask_B, Mask_Int};
pvals = {P_MainA, P_MainB, P_Int};

figure('Color','w', 'Position', [100, 100, 1200, 400], 'Name', 'Global Connectivity Scout (Within)');
for i = 1:3
    subplot(1,3,i);
    
    % 准备绘图数据
    p_plot = pvals{i};
    p_plot(masks{i} == 0) = NaN; % 只显示过 FDR 的区域
    
    if all(isnan(p_plot(:)))
        text(mean(time_axis), mean(f_axis), 'No Sig.', 'Horiz','center', 'FontSize',14);
        xlim([min(time_axis) max(time_axis)]); ylim([min(f_axis) max(f_axis)]);
    else
        imagesc(time_axis, f_axis, p_plot); 
        axis xy; colorbar;
        caxis([0 0.05]); % 锁定 P值颜色范围
        colormap(gca, flipud(parula));
    end
    
    title(plot_titles{i}, 'FontSize', 12, 'FontWeight','bold');
    xlabel('Time (ms)'); ylabel('Freq (Hz)');
end

msgbox({'扫描完成！'; ''; '1. 观察图中的色块区域'; '2. 确定显著的时间段和频率段 (ROI)'; '3. 将这些数值填入全脑网络分析代码中。'});