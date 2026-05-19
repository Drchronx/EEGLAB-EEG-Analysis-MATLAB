clc;
clear;
close all;

%% ================= 1. 选择 EEG .set 文件 =================
[fn_set, pn_set] = uigetfile('*.set', '请选择一个包含电极坐标的 EEG .set 文件');
if isequal(fn_set, 0)
    error('未选择 .set 文件。');
end

EEG = pop_loadset('filename', fn_set, 'filepath', pn_set);
chanlocs = EEG.chanlocs;
nChans = length(chanlocs);

%% ================= 2. 找到 FCz / F5 / F6 =================
labels_all = cell(nChans,1);
for i = 1:nChans
    labels_all{i} = strtrim(chanlocs(i).labels);
end

idx_FCz = find(strcmpi(labels_all, 'FCz'), 1);
idx_F5  = find(strcmpi(labels_all, 'F5'), 1);
idx_F6  = find(strcmpi(labels_all, 'F6'), 1);

if isempty(idx_FCz)
    error('未找到电极 FCz');
end
if isempty(idx_F5)
    error('未找到电极 F5');
end
if isempty(idx_F6)
    error('未找到电极 F6');
end

fprintf('找到电极:\n');
fprintf('FCz = %d\n', idx_FCz);
fprintf('F5  = %d\n', idx_F5);
fprintf('F6  = %d\n', idx_F6);

%% ================= 3. 选择输出文件夹 =================
save_dir = uigetdir(pwd, '选择输出 BrainNet 文件的文件夹');
if isequal(save_dir, 0)
    error('未选择保存文件夹。');
end

%% ================= 4. 生成 node 文件 =================
% BrainNet .node 格式:
% x y z color size label

node_file = fullfile(save_dir, 'FCz_F5_F6.node');
fid = fopen(node_file, 'w');
if fid == -1
    error('无法创建 node 文件。');
end

scale = 85;  % 坐标缩放

for i = 1:nChans
    % 优先使用 X/Y/Z
    if isfield(chanlocs, 'X') && ~isempty(chanlocs(i).X) && ...
       isfield(chanlocs, 'Y') && ~isempty(chanlocs(i).Y) && ...
       isfield(chanlocs, 'Z') && ~isempty(chanlocs(i).Z)

        x = chanlocs(i).X;
        y = chanlocs(i).Y;
        z = chanlocs(i).Z;
    else
        % 没有 XYZ 时尝试 theta/radius 近似转换
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

    % 高亮 FCz / F5 / F6
    if i == idx_FCz
        col = 1;   % 红色组
        sz  = 7;
    elseif i == idx_F5
        col = 2;   % 第二种颜色组
        sz  = 6;
    elseif i == idx_F6
        col = 3;   % 第三种颜色组
        sz  = 6;
    else
        col = 4;   % 其他电极
        sz  = 2;
    end

    fprintf(fid, '%.4f\t%.4f\t%.4f\t%d\t%d\t%s\n', ...
        vec(1), vec(2), vec(3), col, sz, chanlocs(i).labels);
end

fclose(fid);

%% ================= 5. 生成 edge 文件 =================
% 你可以修改下面的连接强度
w_FCz_F5 = 1;
w_FCz_F6 = 1;

% --- 1) FCz-F5 ---
edge_1 = zeros(nChans);
edge_1(idx_FCz, idx_F5) = w_FCz_F5;
edge_1(idx_F5, idx_FCz) = w_FCz_F5;
dlmwrite(fullfile(save_dir, 'FCz_F5.edge'), edge_1, 'delimiter', '\t');

% --- 2) FCz-F6 ---
edge_2 = zeros(nChans);
edge_2(idx_FCz, idx_F6) = w_FCz_F6;
edge_2(idx_F6, idx_FCz) = w_FCz_F6;
dlmwrite(fullfile(save_dir, 'FCz_F6.edge'), edge_2, 'delimiter', '\t');

% --- 3) FCz-F5 和 FCz-F6 同时显示 ---
edge_3 = zeros(nChans);
edge_3(idx_FCz, idx_F5) = w_FCz_F5;
edge_3(idx_F5, idx_FCz) = w_FCz_F5;
edge_3(idx_FCz, idx_F6) = w_FCz_F6;
edge_3(idx_F6, idx_FCz) = w_FCz_F6;
dlmwrite(fullfile(save_dir, 'FCz_F5_and_F6.edge'), edge_3, 'delimiter', '\t');

%% ================= 6. 提示 =================
msgbox({ ...
    '导出完成！'; ...
    ' '; ...
    ['Node 文件: ' fullfile(save_dir, 'FCz_F5_F6.node')]; ...
    ['Edge 文件: ' fullfile(save_dir, 'FCz_F5.edge')]; ...
    ['Edge 文件: ' fullfile(save_dir, 'FCz_F6.edge')]; ...
    ['Edge 文件: ' fullfile(save_dir, 'FCz_F5_and_F6.edge')]; ...
    ' '; ...
    '在 BrainNet Viewer 中加载 .node 和对应 .edge 后点击 Draw 即可。' ...
    });