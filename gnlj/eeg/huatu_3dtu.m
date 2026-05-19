clc;
clear;
close all;

%% 1. 检查 BrainNet Viewer
if exist('BrainNet_MapCfg', 'file') ~= 2
    error('未找到 BrainNet_MapCfg，请先把 BrainNet Viewer 加入 MATLAB 路径。');
end

%% 2. 选择 node / edge / 脑模板
[node_name, node_path] = uigetfile('*.node', '请选择 node 文件');
if isequal(node_name, 0)
    error('未选择 node 文件。');
end
node_file = fullfile(node_path, node_name);

[edge_name, edge_path] = uigetfile('*.edge', '请选择 edge 文件');
if isequal(edge_name, 0)
    error('未选择 edge 文件。');
end
edge_file = fullfile(edge_path, edge_name);

% 自动找 BrainNet 脑模板
brainnet_root = fileparts(which('BrainNet_MapCfg'));
surf_file = '';

cand = { ...
    fullfile(brainnet_root, 'BrainMesh_ICBM152.nv'), ...
    fullfile(brainnet_root, 'BrainMesh_Ch2.nv'), ...
    fullfile(brainnet_root, 'Data', 'BrainMesh_ICBM152.nv'), ...
    fullfile(brainnet_root, 'Data', 'BrainMesh_Ch2.nv') ...
    };

for i = 1:numel(cand)
    if exist(cand{i}, 'file')
        surf_file = cand{i};
        break;
    end
end

if isempty(surf_file)
    [surf_name, surf_path] = uigetfile('*.nv', '请选择 BrainNet 的脑模板 .nv 文件');
    if isequal(surf_name, 0)
        error('未选择脑模板 .nv 文件。');
    end
    surf_file = fullfile(surf_path, surf_name);
end

%% 3. 选择输出图片
[out_name, out_path] = uiputfile('*.png', '保存输出图片为');
if isequal(out_name, 0)
    error('未设置输出图片文件名。');
end
out_img = fullfile(out_path, out_name);

%% 4. 直接调用 BrainNet Viewer
% 第4个参数传空字符串，使用 BrainNet 默认配置
BrainNet_MapCfg(surf_file, node_file, edge_file, '', out_img);

disp('绘图完成。');
disp(out_img);