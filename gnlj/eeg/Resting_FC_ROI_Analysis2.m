clc;
clear;
close all;

%% =========================================================================
%  Resting / TF FC ROI Analysis (Universal Version)
%  兼容两类输入：
%  1) [Ch, Ch, Subj]                    -> 静息态频段平均 FC
%  2) [Freq, Time, Ch, Ch, Subj]        -> 时频 FC
%  3) [Freq, Time, Ch, Ch]              -> 单被试时频 FC
%
%  功能：
%    1. 支持组内/组间设计
%    2. 自动识别数据维度
%    3. ROI 电极对提取
%    4. 时频数据：逐点 t 检验 + FDR + 时频图 + ROI 导出 + BrainNet 导出
%    5. 频段平均数据：ROI 电极对统计 + Excel 导出 + BrainNet 导出
%% =========================================================================

%% ================= 1. 基本设置 =================
design_type = questdlg('请选择实验设计类型:', 'Design Selection', ...
    '组内 (Paired, 同一批人)', '组间 (Independent, 两组人)', '组间 (Independent, 两组人)');
if isempty(design_type)
    return;
end
is_paired = strcmp(design_type, '组内 (Paired, 同一批人)');

prompt = { ...
    '输入分析电极对 (e.g. [6 5]):', ...
    '时间轴参数 [Start Step End]，若为频段平均数据可填 []:', ...
    '频率轴参数 [Start Step End]，若为频段平均数据可填 []:', ...
    '是否进行基线校正? (1=是, 0=否；频段平均数据一般填0):', ...
    '基线范围 [Min Max] (不用可填 []):', ...
    '条件/组别 1 名称:', ...
    '条件/组别 2 名称:' ...
};

def = {'[6 5]', '[0 4 2000]', '[1 1 30]', '0', '[]', 'GroupA', 'GroupB'};
ans_set = inputdlg(prompt, 'FC ROI Parameters', 1, def);
if isempty(ans_set)
    return;
end

chan_pair   = str2num(ans_set{1}); %#ok<ST2NM>
time_par    = str2num(ans_set{2}); %#ok<ST2NM>
f_par       = str2num(ans_set{3}); %#ok<ST2NM>
do_baseline = str2double(ans_set{4});
baseline    = str2num(ans_set{5}); %#ok<ST2NM>
Name_1      = strtrim(ans_set{6});
Name_2      = strtrim(ans_set{7});

if numel(chan_pair) ~= 2
    error('电极对输入应为如 [6 5]。');
end

if isempty(time_par)
    time_axis = [];
else
    if numel(time_par) ~= 3
        error('时间轴参数应为 [Start Step End]，或者填 []。');
    end
    time_axis = time_par(1):time_par(2):time_par(3);
end

if isempty(f_par)
    f_axis = [];
else
    if numel(f_par) ~= 3
        error('频率轴参数应为 [Start Step End]，或者填 []。');
    end
    f_axis = f_par(1):f_par(2):f_par(3);
end

if isempty(baseline)
    do_baseline = 0;
end

%% ================= 2. 载入电极信息 =================
disp('>>> [重要] 请选择一个 *预处理后* (已剔除电极) 的 .set 文件 <<<');
[fn_set, pn_set] = uigetfile('*.set', 'Load EEG .set file');
if isequal(fn_set, 0)
    return;
end

EEG = pop_loadset('filename', fn_set, 'filepath', pn_set);
chanlocs = EEG.chanlocs;
nChans_Total = length(chanlocs);

if any(chan_pair < 1) || any(chan_pair > nChans_Total)
    error('电极编号超出范围，总电极数 = %d', nChans_Total);
end

%% ================= 3. 读取两组数据 =================
if is_paired
    titles = {'Condition 1', 'Condition 2'};
else
    titles = {'Group 1', 'Group 2'};
end

fprintf('正在读取 %s ...\n', titles{1});
[fn1, pn1] = uigetfile('*.mat', ['Select: ' titles{1}]);
if isequal(fn1, 0), return; end
tmp1 = load(fullfile(pn1, fn1));
[data1, varname1] = find_fc_matrix_from_mat(tmp1);
fprintf('第一个文件读取变量: %s\n', varname1);
disp(['变量维度: ', mat2str(size(data1))]);

fprintf('正在读取 %s ...\n', titles{2});
[fn2, pn2] = uigetfile('*.mat', ['Select: ' titles{2}]);
if isequal(fn2, 0), return; end
tmp2 = load(fullfile(pn2, fn2));
[data2, varname2] = find_fc_matrix_from_mat(tmp2);
fprintf('第二个文件读取变量: %s\n', varname2);
disp(['变量维度: ', mat2str(size(data2))]);

%% ================= 4. 判断数据类型 =================
type1 = detect_fc_data_type(data1);
type2 = detect_fc_data_type(data2);

if ~strcmp(type1, type2)
    error('两个文件的数据类型不一致：第一个=%s, 第二个=%s', type1, type2);
end

fprintf('识别到的数据类型为：%s\n', type1);

%% ================= 5. 分析 =================
switch type1

    case 'static_fc'   % [Ch, Ch, Subj]
        run_static_fc_analysis(data1, data2, chan_pair, chanlocs, nChans_Total, ...
            is_paired, Name_1, Name_2);

    case 'tf_fc'       % [Freq, Time, Ch, Ch, Subj]
        run_tf_fc_analysis(data1, data2, chan_pair, chanlocs, nChans_Total, ...
            is_paired, Name_1, Name_2, time_axis, f_axis, do_baseline, baseline);

    otherwise
        error('暂不支持该数据格式。');
end

disp('分析完成。');

%% =========================================================================
%% 子函数区
%% =========================================================================

function [data_out, var_name] = find_fc_matrix_from_mat(tmp)
    vars = fieldnames(tmp);
    data_out = [];
    var_name = '';

    best_score = -inf;

    for i = 1:length(vars)
        x = tmp.(vars{i});
        if ~isnumeric(x)
            continue;
        end

        sx = size(x);
        nd = ndims(x);

        score = -inf;

        % 偏好最像 FC 矩阵的变量
        if nd == 5
            if sx(3) == sx(4)
                score = 1000 + prod(double(sx));
            else
                score = 700 + prod(double(sx));
            end
        elseif nd == 4
            if sx(3) == sx(4)
                score = 900 + prod(double(sx));
            else
                score = 600 + prod(double(sx));
            end
        elseif nd == 3
            if sx(1) == sx(2)
                score = 800 + prod(double(sx));
            else
                score = 500 + prod(double(sx));
            end
        else
            score = prod(double(sx));
        end

        if score > best_score
            best_score = score;
            data_out = x;
            var_name = vars{i};
        end
    end

    if isempty(data_out)
        error('.mat 文件中未找到合适的数值型 FC 变量。');
    end
end

function dtype = detect_fc_data_type(x)
    sx = size(x);
    nd = ndims(x);

    if nd == 5 && sx(3) == sx(4)
        dtype = 'tf_fc';
    elseif nd == 4 && sx(3) == sx(4)
        % 单被试时频数据，补成5维
        dtype = 'tf_fc';
    elseif nd == 3 && sx(1) == sx(2)
        dtype = 'static_fc';
    else
        error('无法识别数据格式，size = %s', mat2str(size(x)));
    end
end

function run_static_fc_analysis(A_fc, B_fc, chan_pair, chanlocs, nChans_Total, is_paired, Name_1, Name_2)

    % 兼容 [Ch, Ch] -> 单被试补成 [Ch, Ch, 1]
    if ndims(A_fc) == 2
        A_fc = reshape(A_fc, size(A_fc,1), size(A_fc,2), 1);
    end
    if ndims(B_fc) == 2
        B_fc = reshape(B_fc, size(B_fc,1), size(B_fc,2), 1);
    end

    [ch1, ch2, nS1] = size(A_fc);
    [ch3, ch4, nS2] = size(B_fc);

    if ch1 ~= ch2 || ch3 ~= ch4
        error('静息态 FC 数据应为 [Ch, Ch, Subj]。');
    end
    if ch1 ~= nChans_Total || ch3 ~= nChans_Total
        error('FC 数据通道数与 .set 文件通道数不一致。');
    end
    if is_paired && nS1 ~= nS2
        error('组内设计要求两组被试数一致。');
    end

    % ROI 电极对取值
    val1 = squeeze(A_fc(chan_pair(1), chan_pair(2), :));
    val2 = squeeze(B_fc(chan_pair(1), chan_pair(2), :));
    val1 = val1(:);
    val2 = val2(:);

    % 统计
    if is_paired
        [~, p, ~, stats] = ttest(val1, val2);
    else
        [~, p, ~, stats] = ttest2(val1, val2);
    end

    fprintf('\n================ ROI统计结果 ================\n');
    fprintf('ROI电极对: [%d %d]\n', chan_pair(1), chan_pair(2));
    fprintf('%s: Mean = %.6f, SD = %.6f, N = %d\n', Name_1, mean(val1,'omitnan'), std(val1,'omitnan'), length(val1));
    fprintf('%s: Mean = %.6f, SD = %.6f, N = %d\n', Name_2, mean(val2,'omitnan'), std(val2,'omitnan'), length(val2));
    fprintf('t = %.6f, p = %.6f\n', stats.tstat, p);
    fprintf('=============================================\n');

    % 画柱状/散点
    figure('Color','w','Position',[200 150 600 500],'Name','Static FC ROI');
    subplot(1,1,1);
    hold on;
    bar(1, mean(val1,'omitnan'));
    bar(2, mean(val2,'omitnan'));
    scatter(ones(size(val1))*1, val1, 35, 'filled', 'jitter','on', 'jitterAmount',0.08);
    scatter(ones(size(val2))*2, val2, 35, 'filled', 'jitter','on', 'jitterAmount',0.08);
    errorbar([1 2], [mean(val1,'omitnan') mean(val2,'omitnan')], ...
        [std(val1,'omitnan')/sqrt(length(val1)) std(val2,'omitnan')/sqrt(length(val2))], 'k', 'LineStyle','none');
    set(gca,'XTick',[1 2],'XTickLabel',{Name_1, Name_2});
    ylabel('FC Value');
    title(sprintf('ROI [%d %d]: t = %.3f, p = %.4f', chan_pair(1), chan_pair(2), stats.tstat, p), 'Interpreter','none');
    box off;

    % 导出 Excel
    choice_excel = questdlg('是否导出 ROI 数据到 Excel?', 'Export Excel', 'Yes', 'No', 'Yes');
    if strcmp(choice_excel, 'Yes')
        if is_paired
            T = table((1:length(val1))', val1, val2, 'VariableNames', {'SubID', Name_1, Name_2});
        else
            T = table([(1:length(val1))'; (1:length(val2))'], ...
                      [ones(length(val1),1); ones(length(val2),1)*2], ...
                      [val1; val2], ...
                      'VariableNames', {'SubID', 'Group', 'Value'});
        end
        [fn, pn] = uiputfile('*.xlsx', 'Save ROI Data');
        if ~isequal(fn,0)
            writetable(T, fullfile(pn, fn));
        end
    end

    % 导出 BrainNet
    choice_3d = questdlg('是否导出 BrainNet Viewer 3D 文件?', '3D Export', 'Yes', 'No', 'Yes');
    if strcmp(choice_3d, 'Yes')
        save_dir = uigetdir(pwd, '选择保存文件夹');
        if isequal(save_dir,0)
            return;
        end

        mean_val1 = mean(val1, 'omitnan');
        mean_val2 = mean(val2, 'omitnan');
        t_val = stats.tstat;

        export_brainnet_files(save_dir, chanlocs, nChans_Total, chan_pair, Name_1, Name_2, mean_val1, mean_val2, t_val);
    end
end

function run_tf_fc_analysis(Data_1, Data_2, chan_pair, chanlocs, nChans_Total, is_paired, Name_1, Name_2, time_axis, f_axis, do_baseline, baseline)

    % 兼容单被试 [F,T,Ch,Ch]
    if ndims(Data_1) == 4
        Data_1 = reshape(Data_1, size(Data_1,1), size(Data_1,2), size(Data_1,3), size(Data_1,4), 1);
    end
    if ndims(Data_2) == 4
        Data_2 = reshape(Data_2, size(Data_2,1), size(Data_2,2), size(Data_2,3), size(Data_2,4), 1);
    end

    [nF, nT, nCh1, nCh2, nS1] = size(Data_1);
    [nF2, nT2, nCh3, nCh4, nS2] = size(Data_2);

    if nCh1 ~= nCh2 || nCh3 ~= nCh4
        error('时频 FC 数据应满足第3/4维相等。');
    end
    if nCh1 ~= nChans_Total || nCh3 ~= nChans_Total
        error('时频 FC 数据通道数与 .set 文件不匹配。');
    end
    if nF ~= nF2 || nT ~= nT2
        error('两个文件的频率/时间维度不一致。');
    end
    if is_paired && nS1 ~= nS2
        error('组内设计要求两组被试数一致。');
    end

    if isempty(time_axis) || length(time_axis) ~= nT
        error('时频数据必须提供正确的时间轴参数，长度需等于 %d。', nT);
    end
    if isempty(f_axis) || length(f_axis) ~= nF
        error('时频数据必须提供正确的频率轴参数，长度需等于 %d。', nF);
    end

    % 提取 ROI
    roi1 = squeeze(Data_1(:, :, chan_pair(1), chan_pair(2), :)); % [F,T,S]
    roi2 = squeeze(Data_2(:, :, chan_pair(1), chan_pair(2), :));

    if ndims(roi1) == 2
        roi1 = reshape(roi1, [nF, nT, 1]);
    end
    if ndims(roi2) == 2
        roi2 = reshape(roi2, [nF, nT, 1]);
    end

    % 基线校正
    if do_baseline
        if isempty(baseline) || numel(baseline) ~= 2
            error('基线范围应为 [min max]。');
        end
        base_idx = time_axis >= baseline(1) & time_axis <= baseline(2);
        if ~any(base_idx)
            error('基线范围未落入时间轴范围内。');
        end

        base_mean1 = mean(roi1(:, base_idx, :), 2);
        base_mean2 = mean(roi2(:, base_idx, :), 2);

        roi1 = roi1 - repmat(base_mean1, [1 nT 1]);
        roi2 = roi2 - repmat(base_mean2, [1 nT 1]);
    end

    fprintf('数据准备就绪。N1=%d, N2=%d\n', size(roi1,3), size(roi2,3));

    %% 逐点统计
    fprintf('正在进行逐点统计检验...\n');
    P_Map = ones(nF, nT);
    T_Map = zeros(nF, nT);

    h_bar = waitbar(0, 'Calculating Stats...');
    for f = 1:nF
        waitbar(f/nF, h_bar);
        for t = 1:nT
            y1 = squeeze(roi1(f, t, :));
            y2 = squeeze(roi2(f, t, :));

            y1 = y1(~isnan(y1));
            y2 = y2(~isnan(y2));

            if isempty(y1) || isempty(y2)
                P_Map(f,t) = NaN;
                T_Map(f,t) = NaN;
                continue;
            end

            if is_paired
                if length(y1) ~= length(y2)
                    P_Map(f,t) = NaN;
                    T_Map(f,t) = NaN;
                else
                    [~, p, ~, stats] = ttest(y1, y2);
                    P_Map(f,t) = p;
                    T_Map(f,t) = stats.tstat;
                end
            else
                [~, p, ~, stats] = ttest2(y1, y2);
                P_Map(f,t) = p;
                T_Map(f,t) = stats.tstat;
            end
        end
    end
    close(h_bar);

    Mask = fdr_bh_mask(P_Map, 0.05);

    %% 绘图
    figure('Color', 'w', 'Position', [100, 100, 1000, 700], 'Name', 'TF ROI Analysis');

    subplot(2,2,1);
    imagesc(time_axis, f_axis, mean(roi1, 3, 'omitnan'));
    axis xy; colorbar;
    title(['Mean: ' Name_1], 'Interpreter','none');
    xlabel('Time'); ylabel('Freq (Hz)');

    subplot(2,2,2);
    imagesc(time_axis, f_axis, mean(roi2, 3, 'omitnan'));
    axis xy; colorbar;
    title(['Mean: ' Name_2], 'Interpreter','none');
    xlabel('Time'); ylabel('Freq (Hz)');

    subplot(2,2,3);
    p_plot = P_Map;
    p_plot(P_Map >= 0.05) = NaN;
    imagesc(time_axis, f_axis, p_plot);
    axis xy; colorbar;
    caxis([0 0.05]);
    title('Uncorrected p < 0.05');
    xlabel('Time'); ylabel('Freq');

    subplot(2,2,4);
    p_fdr = P_Map;
    p_fdr(~Mask) = NaN;
    if all(isnan(p_fdr(:)))
        imagesc(time_axis, f_axis, nan(size(P_Map)));
        axis xy; colorbar;
        text(mean(time_axis), mean(f_axis), 'No Sig. Area', ...
            'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
    else
        imagesc(time_axis, f_axis, p_fdr);
        axis xy; colorbar;
        caxis([0 0.05]);
    end
    title('FDR Corrected (Significant Windows)');
    xlabel('Time'); ylabel('Freq');

    msgbox('请观察右下角图，输入显著时间和频率范围用于提取。');

    %% 导出 ROI 数据
    choice_spss = questdlg('是否提取显著区域数据导出到 Excel?', 'Export SPSS', 'Yes', 'No', 'Yes');

    mean_val1 = mean(roi1(:), 'omitnan');
    mean_val2 = mean(roi2(:), 'omitnan');
    t_val = max(abs(T_Map(:)), [], 'omitnan');

    if strcmp(choice_spss, 'Yes')
        prompt_r = {'输入显著时间范围 [min max]:', '输入显著频率范围 [min max]:'};
        ans_r = inputdlg(prompt_r, 'ROI Definition', 1, ...
            {sprintf('[%g %g]', time_axis(1), time_axis(end)), sprintf('[%g %g]', f_axis(1), f_axis(end))});

        if ~isempty(ans_r)
            tr = str2num(ans_r{1}); %#ok<ST2NM>
            fr = str2num(ans_r{2}); %#ok<ST2NM>

            t_idx = time_axis >= tr(1) & time_axis <= tr(2);
            f_idx = f_axis >= fr(1) & f_axis <= fr(2);

            if ~any(t_idx)
                error('输入的时间范围未覆盖任何时间点。');
            end
            if ~any(f_idx)
                error('输入的频率范围未覆盖任何频点。');
            end

            val1 = squeeze(mean(mean(roi1(f_idx, t_idx, :), 1, 'omitnan'), 2, 'omitnan'));
            val2 = squeeze(mean(mean(roi2(f_idx, t_idx, :), 1, 'omitnan'), 2, 'omitnan'));
            val1 = val1(:);
            val2 = val2(:);

            mean_val1 = mean(val1, 'omitnan');
            mean_val2 = mean(val2, 'omitnan');

            tmp_t = T_Map(f_idx, t_idx);
            t_val = mean(tmp_t(:), 'omitnan');

            if is_paired
                T = table((1:length(val1))', val1, val2, 'VariableNames', {'SubID', Name_1, Name_2});
            else
                T = table([(1:length(val1))'; (1:length(val2))'], ...
                          [ones(length(val1),1); ones(length(val2),1)*2], ...
                          [val1; val2], ...
                          'VariableNames', {'SubID', 'Group', 'Value'});
            end

            [fn, pn] = uiputfile('*.xlsx', 'Save SPSS Data');
            if ~isequal(fn,0)
                writetable(T, fullfile(pn, fn));
            end
        end
    end

    %% 导出 BrainNet
    choice_3d = questdlg('是否导出 3D 脑图文件 (.node/.edge)?', '3D Export', 'Yes', 'No', 'Yes');
    if strcmp(choice_3d, 'Yes')
        save_dir = uigetdir(pwd, '选择保存文件夹');
        if isequal(save_dir,0)
            return;
        end
        export_brainnet_files(save_dir, chanlocs, nChans_Total, chan_pair, Name_1, Name_2, mean_val1, mean_val2, t_val);
    end
end

function export_brainnet_files(save_dir, chanlocs, nChans_Total, chan_pair, Name_1, Name_2, mean_val1, mean_val2, t_val)

    node_file = fullfile(save_dir, 'ROI_Electrodes.node');
    fid = fopen(node_file, 'w');
    if fid == -1
        error('无法创建 node 文件。');
    end

    scale = 85;
    for i = 1:nChans_Total
        if isfield(chanlocs, 'X') && ~isempty(chanlocs(i).X) && ...
           isfield(chanlocs, 'Y') && ~isempty(chanlocs(i).Y) && ...
           isfield(chanlocs, 'Z') && ~isempty(chanlocs(i).Z)

            x = chanlocs(i).X;
            y = chanlocs(i).Y;
            z = chanlocs(i).Z;
        else
            if isfield(chanlocs, 'theta') && isfield(chanlocs, 'radius') && ...
               ~isempty(chanlocs(i).theta) && ~isempty(chanlocs(i).radius)

                theta_rad = chanlocs(i).theta * pi / 180;
                r = chanlocs(i).radius;
                [x, y] = pol2cart(theta_rad, r);
                z = cos(r * pi / 2);
            else
                x = 0; y = 0; z = 0;
            end
        end

        vec = [x, y, z];
        if norm(vec) > 0
            vec = vec / norm(vec) * scale;
        end

        if i == chan_pair(1) || i == chan_pair(2)
            col = 1; sz = 5;
        else
            col = 2; sz = 2;
        end

        if isfield(chanlocs(i), 'labels') && ~isempty(chanlocs(i).labels)
            label_str = chanlocs(i).labels;
        else
            label_str = sprintf('Ch%d', i);
        end

        fprintf(fid, '%.4f\t%.4f\t%.4f\t%d\t%d\t%s\n', vec(1), vec(2), vec(3), col, sz, label_str);
    end
    fclose(fid);

    edge_1 = zeros(nChans_Total);
    edge_1(chan_pair(1), chan_pair(2)) = mean_val1;
    edge_1(chan_pair(2), chan_pair(1)) = mean_val1;
    dlmwrite(fullfile(save_dir, [Name_1 '.edge']), edge_1, 'delimiter', '\t');

    edge_2 = zeros(nChans_Total);
    edge_2(chan_pair(1), chan_pair(2)) = mean_val2;
    edge_2(chan_pair(2), chan_pair(1)) = mean_val2;
    dlmwrite(fullfile(save_dir, [Name_2 '.edge']), edge_2, 'delimiter', '\t');

    edge_t = zeros(nChans_Total);
    edge_t(chan_pair(1), chan_pair(2)) = t_val;
    edge_t(chan_pair(2), chan_pair(1)) = t_val;
    dlmwrite(fullfile(save_dir, 'Diff_T_Value.edge'), edge_t, 'delimiter', '\t');

    msgbox({ ...
        '导出成功！'; ...
        ' '; ...
        '1. ROI_Electrodes.node'; ...
        ['2. ' Name_1 '.edge']; ...
        ['3. ' Name_2 '.edge']; ...
        '4. Diff_T_Value.edge' ...
        });
end

function Mask = fdr_bh_mask(P_Map, q)
    p = P_Map(:);
    valid_idx = ~isnan(p);
    p_valid = p(valid_idx);

    Mask = false(size(P_Map));
    if isempty(p_valid)
        return;
    end

    [p_sorted, ~] = sort(p_valid(:), 'ascend');
    m = numel(p_sorted);
    thresh = (1:m)' / m * q;
    below = p_sorted <= thresh;

    if ~any(below)
        return;
    end

    k = find(below, 1, 'last');
    p_crit = p_sorted(k);

    tmp_mask = false(size(p));
    tmp_mask(valid_idx) = p(valid_idx) <= p_crit;
    Mask = reshape(tmp_mask, size(P_Map));
end