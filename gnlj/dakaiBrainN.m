%% === 强制启动 BrainNet Viewer ===
% 1. 设置 BrainNet Viewer 的模板路径
% (通常在 BrainNet Viewer 安装目录的 Data\SurfTemplate 下)
% 请找到 BrainMesh_ICBM152_Smoothed.nv 文件的位置
surface_file = 'path_to_BrainNetViewer_surface_file.nv'; 

% 如果 MATLAB 找不到这个文件，请提供完整路径，例如：
% surface_file = 'path_to_BrainNetViewer_surface_file.nv';

% 2. 设置你刚才生成的 Node 和 Edge 文件路径
% (假设你刚才导出在 'path_to_your_data' 文件夹下)
node_file = 'path_to_your_node_file.node';  % 修改为你的实际路径
edge_file = 'path_to_your_edge_file.edge';    % 修改为你的实际路径

% 3. 一键打开（绕过 GUI 报错）全脑用
% BrainNet_MapCfg(surface_file, node_file, edge_file);

%% 3. 电极对用  设置可视化参数 (Option) - 进阶技巧！
% 如果您已经保存过一个好看的配置文件（设置好了透明度0.3、红线等），
% 可以把那个 .mat 文件的路径写在这里。
% 如果是第一次运行，留空即可，脚本会用默认设置打开。
option_file = '-'; 

% 示例：如果您保存了配置文件，就把下面这行取消注释
% option_file = 'path_to_your_option_file.mat';

%% 4. 启动绘图 (核心命令)
% 语法: BrainNet_MapCfg(Surf, Node, Edge, Volume/Mapping, Config)
% 注意：第四个参数我们传 '-', 代表 Mapping 为空
fprintf('正在生成 3D 脑图...\n');
fprintf('Surface: %s\n', surface_file);
fprintf('Node:    %s\n', node_file);
fprintf('Edge:    %s\n', edge_file);

try
    BrainNet_MapCfg(surface_file, node_file, edge_file, '-', option_file);
    fprintf('>>> 绘图成功！请在弹出的窗口中调整美化。\n');
catch ME
    fprintf('>>> 报错了！错误信息:\n%s\n', ME.message);
end