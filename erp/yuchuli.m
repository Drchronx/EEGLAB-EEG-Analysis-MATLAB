%%
%% 基于 MATLAB 和 EEGLAB 的脑电任务态数据处理示例流程
% 本脚本作为公开发布的非商业 EEG 分析工具的一部分发布。使用前请根据自己的数据路径、
% 事件标记、通道设置和实验设计参数进行配置。
% 包含三部分
%1.bp设备脑电数据预处理批处理代码 
%2.erp时域分析（单因素双因素重复测量 找数据显著时间窗口  找数据显著电极 提取spss统计数据 提取波形图绘图数据 绘制地形图波形图 ）
%3.时频分析（单因素双因素重复测量 找数据显著频率段和时间窗口  找数据显著电极 提取spss统计数据 提取能量图谱图和能量地形图）

%%  第一部分，基于curry设备预处理批处理，这一部分因为我电脑处理器原因用的是2020版MATLAB编写
%% 数据格式格式转化
clc;clear;
%% Specify Basic information of different groups     将cdt文件转换成set文件，降采样为500，去除无用电极，保存。
group_dir = 'path_to_your_data';     % 此处路径需要设置为自己的文件目录
group_files = dir([group_dir, filesep, '*.dap']);  %filesep是\的意思
for i=1:length(group_files)
    subj_fn = group_files(i).name;
    EEG = loadcurry(strcat(group_dir, filesep, subj_fn), 'CurryLocations', 'False');    %导入原始数据
    EEG = pop_resample( EEG, 500);   %降采样
    EEG = pop_select( EEG, 'rmchannel',{'HEO','VEO','EKG','EMG','TRIGGER','CB1','CB2'});
    EEG = pop_reref( EEG, [33 43] );
    EEG = pop_eegfiltnew(EEG, 'locutoff',0.1,'plotfreqz',1);
    EEG = pop_eegfiltnew(EEG, 'hicutoff',40,'plotfreqz',1);
    EEG = pop_eegfiltnew(EEG, 'locutoff',48,'hicutoff',52,'revfilt',1,'plotfreqz',1); 
    EEG = pop_epoch( EEG, {   '1'  '2'  '3'  '4'  }, [-1  2], 'newname', 'Neuroscan Curry file resampled epochs', 'epochinfo', 'yes');
    EEG = pop_rmbase( EEG, [-1000 0] ,[]);

    EEG = pop_saveset( EEG, 'filename',strcat(group_files(i).name(1:end-4), '.set'), 'filepath',strcat(group_dir, filesep, '_step1'));   %注意需要在运行代码之前，文件目录下建一个_resam_remch的文件夹，以下雷同
end

%%文件夹_step1用于存放我们经过降采样、带通滤波、陷波滤波、去除无关电极之后的数据；
%文件夹_preica用于存放我们经过检查有无坏段坏导之后的数据；
%文件夹_ica用于存放我们跑完ICA之后的数据；
%文件夹_rm_ica用于存放ICLabel自动去完伪迹以及全脑平均重参考之后的数据；
%母文件夹demo下存放原始采集的数据。
%% 第一部分，基于curry设备预处理批处理
%% 数据格式格式转化
clc; clear; close all; % close all 可以关闭之前残留的图窗

%% 初始化
% 启动 EEGLAB 但不显示图形界面，节省资源
if ~exist('ALLCOM', 'var')
    eeglab nogui;
end

%% Specify Basic information of different groups
group_dir = 'path_to_your_data';      % 此处路径需要设置为自己的文件目录
save_dir = [group_dir, filesep, '_step1']; % 定义保存路径

% 自动检查文件夹是否存在，不存在则创建（不用手动建了）
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

group_files = dir([group_dir, filesep, '*.dap']);  % filesep是\的意思

%% 循环处理
for i = 1:length(group_files)
    subj_fn = group_files(i).name;
    
    % 1. 导入原始数据
    % 注意：loadcurry有时在不同版本表现不同，如果报错请检查路径是否拼接正确
    EEG = loadcurry(fullfile(group_dir, subj_fn), 'CurryLocations', 'False'); 
    
    % 2. 降采样
    EEG = pop_resample( EEG, 500); 
    
    % 3. 去除无用电极
    EEG = pop_select( EEG, 'rmchannel',{'HEO','VEO','EKG','EMG','TRIGGER','CB1','CB2'});
    
    % 4. 重参考 (M1 M2 一般是 33 43，请确保你的通道号对应正确)
    EEG = pop_reref( EEG, [33 43] );
    
    % 5. 滤波 (关键修改：plotfreqz 改为 0)
    % 高通 0.1Hz
    EEG = pop_eegfiltnew(EEG, 'locutoff', 0.1, 'plotfreqz', 0); 
    % 低通 40Hz
    EEG = pop_eegfiltnew(EEG, 'hicutoff', 40, 'plotfreqz', 0);
    % 陷波 50Hz (48-52Hz)
    EEG = pop_eegfiltnew(EEG, 'locutoff', 48, 'hicutoff', 52, 'revfilt', 1, 'plotfreqz', 0); 
    
    % 6. 分段 (Epoch)
    EEG = pop_epoch( EEG, {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'}, [-1  2], ...
        'newname', 'Neuroscan Curry file resampled epochs', 'epochinfo', 'yes');
    
    % 7. 基线校正
    EEG = pop_rmbase( EEG, [-1000 0] ,[]);
    
    % 8. 保存数据
    % 自动保存到 _step1 文件夹
    save_name = [subj_fn(1:end-4), '.set'];
    EEG = pop_saveset( EEG, 'filename', save_name, 'filepath', save_dir); 
    
    fprintf('文件 %s 处理完成 (%d / %d)\n', save_name, i, length(group_files));
end

disp('所有数据预处理完成！');

%% %% 


%% 


%% 
%% 标准化 Mark 合并工具 (从 Step1 到 Step2)
clc; clear;

% 输入文件夹 (上一步处理好的 _step1)
input_dir = 'path_to_your_data'; 

% 输出文件夹 (自动新建 _step2)
output_dir = 'path_to_your_data';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% 获取所有 .set 文件
files = dir([input_dir, filesep, '*.set']);

% 2. 定义合并规则 (在这里修改最方便！)
% 规则A：将 4, 5, 8, 9 合并为 11 (示例：4个合并成1个)
Target_Mark_A = '11'; 
Source_Marks_A = {'4', '5', '8', '9'}; 

% 规则B：将 6, 7 合并为 12
Target_Mark_B = '12';
Source_Marks_B = {'6', '7'};

% 3. 循环处理
for i = 1:length(files)
    filename = files(i).name;
    fprintf('正在处理第 %d/%d 个被试: %s ...\n', i, length(files), filename);
    
    % --- 读取数据 ---
    EEG = pop_loadset('filename', filename, 'filepath', input_dir);
    
    % --- 核心：修改 Mark ---
    % 遍历所有 Event
    for e = 1:length(EEG.event)
        curr_type = EEG.event(e).type; % 获取当前Mark
        % 确保 curr_type 是字符串格式 (有些时候导入后是数字，转一下保险)
        if isnumeric(curr_type)
            curr_type = num2str(curr_type);
        end
        % 应用规则 A (使用 ismember 函数，只要属于集合中的任意一个，就修改)
        if ismember(curr_type, Source_Marks_A)
            EEG.event(e).type = Target_Mark_A;
        % 应用规则 B
        elseif ismember(curr_type, Source_Marks_B)
            EEG.event(e).type = Target_Mark_B;
        end
    end
    % --- 检查并同步 ---
    % 这一步至关重要，确保 epoch 信息和 event 信息一致
    EEG = eeg_checkset(EEG, 'eventconsistency');
    % --- 保存数据 ---
    % 保持原文件名，存入 _step2
    EEG = pop_saveset(EEG, 'filename', filename, 'filepath', output_dir);
end

disp('>>> 所有数据 Mark 合并完成，已保存至 _step2 文件夹。');

%% 菜单操作：打开eeglab导入_step1的文件，看有无坏导坏段，若有坏导则插值；若有坏段则删掉。完了之后保存数据到_preica文件夹，文件名称保持不变。

%%文件导入：打开eeglab，点击File->load existing dataset；
%文件保存：点击File->save current dataset as。

%%   运行 ICA
clc;clear;
group1_dir = 'path_to_your_data';     % 此处路径需要设置为自己的文件目录
group1_dir1 = 'path_to_your_data';     % 此处路径需要设置为自己的文件目录
group1_files = dir([group1_dir1, filesep, '*.set']);  %filesep是\的意思
for i=1:length(group1_files)
    subj_fn = group1_files(i).name;
    EEG = pop_loadset('filename',strcat(subj_fn(1:end-4), '.set'), 'filepath', strcat(group1_dir, filesep, '_preica')); %导入数据
    EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'interrupt','on');   % 跑ICA
    EEG = pop_saveset( EEG, 'filename',strcat(group1_files(i).name(1:end-4), '.set'), 'filepath',strcat(group1_dir, filesep, '_ica'));  %保存数据
    
end
%% 使用ICLabel自动去除ICA成分
clc;clear;
group1_dir ='path_to_your_data';        % 此处路径需要设置为自己的文件目录
group1_dir2 = 'path_to_your_data'; 
group1_files = dir([group1_dir2, filesep, '*.set']);  %filesep是\的意思
for i=1:length(group1_files)
    subj_fn = group1_files(i).name;
    EEG = pop_loadset('filename',strcat(subj_fn(1:end-4), '.set'), 'filepath', group1_dir2);
    EEG = pop_iclabel(EEG, 'default');
    EEG = pop_icflag(EEG, [NaN NaN;0.3 1;0.3 1;0.3 1;0.3 1;NaN NaN;NaN NaN]); % 标记伪迹成分。这里可以自定义设定阈值，依次为Brain, Muscle, Eye, Heart, Line Noise, Channel Noise, Other.
    EEG = pop_subcomp( EEG, [], 0)   %去除上述伪迹成分
    %EEG = pop_reref( EEG, []);    %全脑平均重参考
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',strcat(group1_files(i).name(1:end-4), '.set'), 'filepath',strcat(group1_dir, filesep, '_rm_ica')); 
end

%% 手动浏览数据 确保预处理之后的数据是干净的
