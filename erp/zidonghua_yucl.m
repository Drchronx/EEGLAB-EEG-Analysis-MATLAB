clc; clear; close all;

%% ============================================================
%  EEGLAB task EEG one-click preprocessing pipeline
%  只需要改这一段 CONFIG，下面代码不要反复改路径。
%  流程:
%    1. raw .dap -> _step1
%    2. 自动坏导插值 + 坏段剔除 -> _preica
%    3. ICA -> _ica
%    4. ICLabel + 最终坏段筛查 -> _done
%% ============================================================

cfg = struct();

% ---------- 路径 ----------
cfg.group_dir = 'J:\脑电实验\处理数据\ai信息框架个人数据\demo';
cfg.raw_pattern = '*.dap';

cfg.step1_dirname = '_step1';
cfg.preica_dirname = '_preica';
cfg.ica_dirname = '_ica';
cfg.done_dirname = 'done';

% ---------- 是否运行各阶段 ----------
cfg.run_step1 = true;
cfg.run_preica_qc = true;
cfg.run_ica = true;
cfg.run_iclabel_done = true;
cfg.overwrite_existing = true;

% ---------- EEGLAB ----------
cfg.start_eeglab_nogui = true;
cfg.curry_locations = 'False';

% ---------- Step 1: 基础预处理 ----------
cfg.srate = 500;
cfg.remove_channels = {'HEO','VEO','EKG','EMG','TRIGGER','CB1','CB2'};

% 可以写通道编号，也可以写通道名。
% 重要: 如果用编号，必须确认删掉无用电极之后编号仍然对应你的参考电极。
cfg.ref_channels = [33 43];
% cfg.ref_channels = {'M1','M2'};

cfg.highpass_hz = 0.1;
cfg.lowpass_hz = 40;
cfg.notch_band_hz = [48 52];

cfg.epoch_markers = {'1','2','3','4','5','6','7','8','9','10','11','12', ...
                     '20','21','22','23','24','25'};
cfg.epoch_window_sec = [-1 2];
cfg.baseline_ms = [-1000 0];

% ---------- Step 2: _preica 质控 ----------
cfg.badchan_threshold = 3;             % 3较严格；误判多可改4或5
cfg.badchan_measure = 'kurt';          % 'kurt' 或 'prob'
cfg.preica_threshold_candidates_uV = [300 400];
cfg.preica_max_reject_rate = 0.50;
cfg.preica_min_remaining_trials = 20;

% ---------- Step 3: ICA ----------
cfg.ica_type = 'runica';
cfg.ica_extended = 1;
cfg.ica_interrupt = 'off';             % 一键跑建议off；想手动中断可改'on'

% ---------- Step 4: ICLabel + 最终筛查 ----------
% ICLabel行顺序:
% Brain, Muscle, Eye, Heart, Line Noise, Channel Noise, Other
cfg.icflag_thresholds = [ ...
    NaN NaN; ...
    0.3 1; ...
    0.3 1; ...
    0.3 1; ...
    0.3 1; ...
    NaN NaN; ...
    NaN NaN];

cfg.final_thresh_uV = 150;             % ICA后最终筛查: ±150 uV
cfg.final_max_reject_rate = 0.50;
cfg.final_min_remaining_trials = 30;

%% ============================================================
%  One-click run
%% ============================================================

if cfg.start_eeglab_nogui
    if ~exist('ALLCOM', 'var')
        eeglab nogui;
    end
end

cfg.step1_dir = fullfile(cfg.group_dir, cfg.step1_dirname);
cfg.preica_dir = fullfile(cfg.group_dir, cfg.preica_dirname);
cfg.ica_dir = fullfile(cfg.group_dir, cfg.ica_dirname);
cfg.done_dir = fullfile(cfg.group_dir, cfg.done_dirname);

ensure_dir(cfg.step1_dir);
ensure_dir(cfg.preica_dir);
ensure_dir(cfg.ica_dir);
ensure_dir(cfg.done_dir);

fprintf('\n==============================\n');
fprintf('EEGLAB 一键预处理开始\n');
fprintf('数据根目录: %s\n', cfg.group_dir);
fprintf('==============================\n');

if cfg.run_step1
    run_step1(cfg);
end

if cfg.run_preica_qc
    run_preica_qc(cfg);
end

if cfg.run_ica
    run_ica_stage(cfg);
end

if cfg.run_iclabel_done
    run_iclabel_done(cfg);
end

fprintf('\n==============================\n');
fprintf('全部流程完成。\n');
fprintf('最终好被试数据目录: %s\n', cfg.done_dir);
fprintf('==============================\n');

%% ============================================================
%  Local functions
%% ============================================================

function ensure_dir(pathstr)
    if ~exist(pathstr, 'dir')
        mkdir(pathstr);
    end
end

function run_step1(cfg)
    fprintf('\n========== Step 1: raw -> _step1 ==========\n');
    files = dir(fullfile(cfg.group_dir, cfg.raw_pattern));
    if isempty(files)
        warning('未找到原始文件: %s', fullfile(cfg.group_dir, cfg.raw_pattern));
        return;
    end

    for i = 1:length(files)
        subj_fn = files(i).name;
        save_name = [subj_fn(1:end-4), '.set'];
        save_path = fullfile(cfg.step1_dir, save_name);

        if exist(save_path, 'file') && ~cfg.overwrite_existing
            fprintf('[跳过] %s 已存在。\n', save_name);
            continue;
        end

        fprintf('\n[%d/%d] Step1处理: %s\n', i, length(files), subj_fn);
        try
            EEG = loadcurry(fullfile(cfg.group_dir, subj_fn), ...
                'CurryLocations', cfg.curry_locations);

            EEG = pop_resample(EEG, cfg.srate);

            rm_chs = existing_channels(EEG, cfg.remove_channels);
            if ~isempty(rm_chs)
                EEG = pop_select(EEG, 'rmchannel', rm_chs);
            end

            ref_chs = resolve_ref_channels(EEG, cfg.ref_channels);
            if ~isempty(ref_chs)
                EEG = pop_reref(EEG, ref_chs);
            else
                warning('%s 未找到参考电极，跳过重参考。', subj_fn);
            end

            EEG = pop_eegfiltnew(EEG, 'locutoff', cfg.highpass_hz, 'plotfreqz', 0);
            EEG = pop_eegfiltnew(EEG, 'hicutoff', cfg.lowpass_hz, 'plotfreqz', 0);
            EEG = pop_eegfiltnew(EEG, ...
                'locutoff', cfg.notch_band_hz(1), ...
                'hicutoff', cfg.notch_band_hz(2), ...
                'revfilt', 1, 'plotfreqz', 0);

            EEG = pop_epoch(EEG, cfg.epoch_markers, cfg.epoch_window_sec, ...
                'newname', 'Curry resampled epochs', 'epochinfo', 'yes');
            EEG = pop_rmbase(EEG, cfg.baseline_ms, []);
            EEG = eeg_checkset(EEG);

            EEG = pop_saveset(EEG, 'filename', save_name, 'filepath', cfg.step1_dir);
            fprintf('[完成] 保存: %s\n', save_path);
        catch ME
            fprintf(2, '[错误] Step1失败 %s: %s\n', subj_fn, ME.message);
        end
    end
end

function run_preica_qc(cfg)
    fprintf('\n========== Step 2: _step1 -> _preica ==========\n');
    files = dir(fullfile(cfg.step1_dir, '*.set'));
    if isempty(files)
        warning('未找到_step1数据: %s', cfg.step1_dir);
        return;
    end

    good_rows = {};
    bad_rows = {};

    good_varnames = {'subject','original_trials','threshold_uV', ...
        'rejected_trials','reject_rate_percent','remaining_trials', ...
        'bad_channel_count','bad_channel_names','median_absmax_uV', ...
        'p95_absmax_uV','max_absmax_uV'};
    bad_varnames = {'subject','original_trials','best_threshold_uV', ...
        'rejected_trials_at_best_threshold','reject_rate_percent', ...
        'remaining_if_rejected','bad_channel_count','bad_channel_names', ...
        'median_absmax_uV','p95_absmax_uV','max_absmax_uV','reason'};

    for i = 1:length(files)
        subj_fn = files(i).name;
        fprintf('\n[%d/%d] PreICA质控: %s\n', i, length(files), subj_fn);

        try
            EEG = pop_loadset('filename', subj_fn, 'filepath', cfg.step1_dir);
            EEG = eeg_checkset(EEG);
            original_trials = EEG.trials;

            [median_absmax, p95_absmax, max_absmax] = epoch_absmax_stats(EEG);
            fprintf('数据幅值: median=%.1f uV, p95=%.1f uV, max=%.1f uV\n', ...
                median_absmax, p95_absmax, max_absmax);

            EEG_original = EEG;
            [EEG, bad_chs] = pop_rejchan(EEG, ...
                'elec', 1:EEG.nbchan, ...
                'threshold', cfg.badchan_threshold, ...
                'norm', 'on', ...
                'measure', cfg.badchan_measure);

            bad_chs = bad_chs(:)';
            bad_ch_names = channel_names_from_indices(EEG_original, bad_chs);
            fprintf('自动检测坏导数量: %d\n', length(bad_chs));
            if ~isempty(bad_ch_names)
                fprintf('坏导: %s\n', bad_ch_names);
            end

            if ~isempty(bad_chs)
                EEG = pop_interp(EEG, EEG_original.chanlocs, 'spherical');
                EEG = eeg_checkset(EEG);
                fprintf('已完成坏导插值。\n');
            end

            [chosen_thresh, chosen_rej, chosen_rate, chosen_n_rej, ...
                best_thresh, best_rate, best_n_rej] = choose_epoch_threshold( ...
                    EEG, cfg.preica_threshold_candidates_uV, cfg.preica_max_reject_rate);

            if isempty(chosen_thresh)
                reason = sprintf('bad_epoch_rate_too_high_even_at_%duV', best_thresh);
                remaining_if_rejected = original_trials - best_n_rej;
                fprintf('[剔除] %s，不保存到_preica。\n', reason);

                bad_rows(end+1, :) = {subj_fn, original_trials, best_thresh, ...
                    best_n_rej, best_rate * 100, remaining_if_rejected, ...
                    length(bad_chs), bad_ch_names, median_absmax, ...
                    p95_absmax, max_absmax, reason};
                continue;
            end

            if chosen_n_rej > 0
                EEG = pop_rejepoch(EEG, chosen_rej(:)', 0);
                EEG = eeg_checkset(EEG);
            end

            remaining_trials = EEG.trials;
            fprintf('采用阈值 ±%d uV，删除坏段 %d，剩余 trials: %d\n', ...
                chosen_thresh, chosen_n_rej, remaining_trials);

            if remaining_trials < cfg.preica_min_remaining_trials
                reason = sprintf('too_few_remaining_trials_less_than_%d', ...
                    cfg.preica_min_remaining_trials);
                fprintf('[剔除] %s，不保存到_preica。\n', reason);

                bad_rows(end+1, :) = {subj_fn, original_trials, chosen_thresh, ...
                    chosen_n_rej, chosen_rate * 100, remaining_trials, ...
                    length(bad_chs), bad_ch_names, median_absmax, ...
                    p95_absmax, max_absmax, reason};
                continue;
            end

            EEG = pop_saveset(EEG, 'filename', subj_fn, 'filepath', cfg.preica_dir);
            fprintf('[好被试] 保存到_preica: %s\n', fullfile(cfg.preica_dir, subj_fn));

            good_rows(end+1, :) = {subj_fn, original_trials, chosen_thresh, ...
                chosen_n_rej, chosen_rate * 100, remaining_trials, ...
                length(bad_chs), bad_ch_names, median_absmax, ...
                p95_absmax, max_absmax};
        catch ME
            fprintf(2, '[错误] PreICA失败 %s: %s\n', subj_fn, ME.message);
            bad_rows(end+1, :) = {subj_fn, NaN, NaN, NaN, NaN, NaN, ...
                NaN, '', NaN, NaN, NaN, ['error_', ME.message]};
        end
    end

    write_log_table(fullfile(cfg.group_dir, 'good_subjects_preica_log.csv'), ...
        good_rows, good_varnames);
    write_log_table(fullfile(cfg.group_dir, 'bad_subjects_preica_log.csv'), ...
        bad_rows, bad_varnames);

    fprintf('PreICA好被试数量: %d\n', size(good_rows, 1));
    fprintf('PreICA坏被试数量: %d\n', size(bad_rows, 1));
end

function run_ica_stage(cfg)
    fprintf('\n========== Step 3: _preica -> _ica ==========\n');
    files = dir(fullfile(cfg.preica_dir, '*.set'));
    if isempty(files)
        warning('未找到_preica数据: %s', cfg.preica_dir);
        return;
    end

    for i = 1:length(files)
        subj_fn = files(i).name;
        save_path = fullfile(cfg.ica_dir, subj_fn);

        if exist(save_path, 'file') && ~cfg.overwrite_existing
            fprintf('[跳过] %s 已存在。\n', subj_fn);
            continue;
        end

        fprintf('\n[%d/%d] ICA: %s\n', i, length(files), subj_fn);
        try
            EEG = pop_loadset('filename', subj_fn, 'filepath', cfg.preica_dir);
            EEG = eeg_checkset(EEG);
            EEG = pop_runica(EEG, 'icatype', cfg.ica_type, ...
                'extended', cfg.ica_extended, 'interrupt', cfg.ica_interrupt);
            EEG = eeg_checkset(EEG);
            EEG = pop_saveset(EEG, 'filename', subj_fn, 'filepath', cfg.ica_dir);
            fprintf('[完成] ICA保存: %s\n', save_path);
        catch ME
            fprintf(2, '[错误] ICA失败 %s: %s\n', subj_fn, ME.message);
        end
    end
end

function run_iclabel_done(cfg)
    fprintf('\n========== Step 4: _ica -> _done ==========\n');
    if exist('pop_iclabel', 'file') ~= 2
        error('未找到 ICLabel 插件函数 pop_iclabel，请先在 EEGLAB 中安装 ICLabel。');
    end

    files = dir(fullfile(cfg.ica_dir, '*.set'));
    if isempty(files)
        warning('未找到_ica数据: %s', cfg.ica_dir);
        return;
    end

    good_rows = {};
    bad_rows = {};

    good_varnames = {'subject','original_trials','removed_ic_count', ...
        'removed_ic_indices','threshold_uV','rejected_trials', ...
        'reject_rate_percent','remaining_trials'};
    bad_varnames = {'subject','original_trials','removed_ic_count', ...
        'removed_ic_indices','threshold_uV','rejected_trials', ...
        'reject_rate_percent','remaining_trials','reason'};

    for i = 1:length(files)
        subj_fn = files(i).name;
        fprintf('\n[%d/%d] ICLabel + done终筛: %s\n', i, length(files), subj_fn);

        try
            EEG = pop_loadset('filename', subj_fn, 'filepath', cfg.ica_dir);
            EEG = eeg_checkset(EEG);
            original_trials = EEG.trials;

            EEG = pop_iclabel(EEG, 'default');
            EEG = pop_icflag(EEG, cfg.icflag_thresholds);
            ic_to_remove = find(EEG.reject.gcompreject);
            removed_ic_count = length(ic_to_remove);
            removed_ic_indices = num2str(ic_to_remove);

            fprintf('ICLabel标记伪迹成分数量: %d\n', removed_ic_count);
            if removed_ic_count > 0
                fprintf('删除IC成分: %s\n', removed_ic_indices);
                EEG = pop_subcomp(EEG, ic_to_remove, 0);
            else
                EEG = pop_subcomp(EEG, [], 0);
            end
            EEG = eeg_checkset(EEG);

            EEG_tmp = pop_eegthresh(EEG, ...
                1, 1:EEG.nbchan, ...
                -cfg.final_thresh_uV, cfg.final_thresh_uV, ...
                EEG.xmin, EEG.xmax, ...
                0, 0);

            if isfield(EEG_tmp.reject, 'rejthresh') && ~isempty(EEG_tmp.reject.rejthresh)
                rej = logical(EEG_tmp.reject.rejthresh);
            else
                rej = false(1, EEG_tmp.trials);
            end

            n_rej = sum(rej);
            reject_rate = n_rej / EEG_tmp.trials;
            remaining_trials = EEG_tmp.trials - n_rej;

            fprintf('ICA后 ±%d uV终筛: 标记坏段 %d/%d (%.1f%%)\n', ...
                cfg.final_thresh_uV, n_rej, EEG_tmp.trials, reject_rate * 100);

            if reject_rate > cfg.final_max_reject_rate
                reason = sprintf('final_bad_epoch_rate_over_%.0f_percent', ...
                    cfg.final_max_reject_rate * 100);
                fprintf('[剔除] %s，不保存到_done。\n', reason);
                bad_rows(end+1, :) = {subj_fn, original_trials, ...
                    removed_ic_count, removed_ic_indices, cfg.final_thresh_uV, ...
                    n_rej, reject_rate * 100, remaining_trials, reason};
                continue;
            end

            if remaining_trials < cfg.final_min_remaining_trials
                reason = sprintf('remaining_trials_less_than_%d', ...
                    cfg.final_min_remaining_trials);
                fprintf('[剔除] %s，不保存到_done。\n', reason);
                bad_rows(end+1, :) = {subj_fn, original_trials, ...
                    removed_ic_count, removed_ic_indices, cfg.final_thresh_uV, ...
                    n_rej, reject_rate * 100, remaining_trials, reason};
                continue;
            end

            if n_rej > 0
                EEG = pop_rejepoch(EEG, rej(:)', 0);
                EEG = eeg_checkset(EEG);
            end
            remaining_trials = EEG.trials;

            EEG = pop_saveset(EEG, 'filename', subj_fn, 'filepath', cfg.done_dir);
            fprintf('[好被试] 保存到_done: %s\n', fullfile(cfg.done_dir, subj_fn));

            good_rows(end+1, :) = {subj_fn, original_trials, ...
                removed_ic_count, removed_ic_indices, cfg.final_thresh_uV, ...
                n_rej, reject_rate * 100, remaining_trials};
        catch ME
            fprintf(2, '[错误] ICLabel/done失败 %s: %s\n', subj_fn, ME.message);
            bad_rows(end+1, :) = {subj_fn, NaN, NaN, '', ...
                cfg.final_thresh_uV, NaN, NaN, NaN, ['error_', ME.message]};
        end
    end

    write_log_table(fullfile(cfg.group_dir, 'done_good_subjects_log.csv'), ...
        good_rows, good_varnames);
    write_log_table(fullfile(cfg.group_dir, 'done_bad_subjects_log.csv'), ...
        bad_rows, bad_varnames);

    fprintf('Done好被试数量: %d\n', size(good_rows, 1));
    fprintf('Done坏被试数量: %d\n', size(bad_rows, 1));
    fprintf('Done数据目录: %s\n', cfg.done_dir);
end

function chs = existing_channels(EEG, wanted)
    chs = {};
    if isempty(wanted) || ~isfield(EEG, 'chanlocs') || isempty(EEG.chanlocs)
        return;
    end
    labels = {EEG.chanlocs.labels};
    labels_upper = upper(labels);
    for i = 1:length(wanted)
        if any(strcmpi(wanted{i}, labels_upper)) || any(strcmpi(wanted{i}, labels))
            chs{end+1} = wanted{i}; %#ok<AGROW>
        end
    end
end

function ref = resolve_ref_channels(EEG, ref_cfg)
    if isnumeric(ref_cfg)
        ref = ref_cfg;
        return;
    end
    ref = [];
    if isempty(ref_cfg)
        return;
    end
    labels = {EEG.chanlocs.labels};
    labels_upper = upper(labels);
    for i = 1:length(ref_cfg)
        idx = find(strcmp(labels_upper, upper(ref_cfg{i})), 1);
        if isempty(idx)
            warning('未找到参考电极: %s', ref_cfg{i});
        else
            ref(end+1) = idx; %#ok<AGROW>
        end
    end
end

function names = channel_names_from_indices(EEG, idxs)
    names = '';
    if isempty(idxs)
        return;
    end
    try
        labels = {EEG.chanlocs(idxs).labels};
        names = strjoin(labels, ',');
    catch
        names = num2str(idxs);
    end
end

function [median_absmax, p95_absmax, max_absmax] = epoch_absmax_stats(EEG)
    epoch_absmax = squeeze(max(max(abs(EEG.data), [], 1), [], 2));
    epoch_absmax = epoch_absmax(:);
    sorted_absmax = sort(epoch_absmax);
    p95_idx = max(1, round(0.95 * length(sorted_absmax)));
    median_absmax = median(epoch_absmax);
    p95_absmax = sorted_absmax(p95_idx);
    max_absmax = max(epoch_absmax);
end

function [chosen_thresh, chosen_rej, chosen_rate, chosen_n_rej, ...
          best_thresh, best_rate, best_n_rej] = choose_epoch_threshold( ...
          EEG, threshold_candidates, max_reject_rate)

    chosen_thresh = [];
    chosen_rej = [];
    chosen_rate = [];
    chosen_n_rej = [];

    best_thresh = NaN;
    best_rate = Inf;
    best_n_rej = NaN;

    for t = 1:length(threshold_candidates)
        thresh_uV = threshold_candidates(t);
        EEG_tmp = pop_eegthresh(EEG, ...
            1, 1:EEG.nbchan, ...
            -thresh_uV, thresh_uV, ...
            EEG.xmin, EEG.xmax, ...
            0, 0);

        if isfield(EEG_tmp.reject, 'rejthresh') && ~isempty(EEG_tmp.reject.rejthresh)
            rej = logical(EEG_tmp.reject.rejthresh);
        else
            rej = false(1, EEG_tmp.trials);
        end

        n_rej = sum(rej);
        rej_rate = n_rej / EEG_tmp.trials;

        fprintf('阈值 ±%d uV: 标记坏段 %d/%d (%.1f%%)\n', ...
            thresh_uV, n_rej, EEG_tmp.trials, rej_rate * 100);

        if rej_rate < best_rate
            best_rate = rej_rate;
            best_thresh = thresh_uV;
            best_n_rej = n_rej;
        end

        if rej_rate <= max_reject_rate
            chosen_thresh = thresh_uV;
            chosen_rej = rej;
            chosen_rate = rej_rate;
            chosen_n_rej = n_rej;
            break;
        end
    end
end

function write_log_table(filepath, rows, varnames)
    if isempty(rows)
        tbl = cell2table(cell(0, length(varnames)), 'VariableNames', varnames);
    else
        tbl = cell2table(rows, 'VariableNames', varnames);
    end
    writetable(tbl, filepath);
    fprintf('日志已保存: %s\n', filepath);
end
