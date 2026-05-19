%%%%%%%%%%
%% eeg数据
clc;clear;
%% Specify Basic information of different groups     将cdt文件转换成set文件，降采样为500，去除无用电极，保存。
group_dir = 'path_to_your_data';     % 此处路径需要设置为自己的文件目录
group_files = dir([group_dir, filesep, '*.dap']);  
eventTypeNumber = 15; % 例如21, 22, 23等

% 您想要提取的事件类型的数值
eventTypeNumber2 = 20; % 例如21, 22, 23等

%filesep是\的意思
for i=1:length(group_files)
    subj_fn = group_files(i).name;
    EEG = loadcurry(strcat(group_dir, filesep, subj_fn), 'CurryLocations', 'False');    %导入原始数据
    EEG = pop_resample( EEG, 250);   %降采样
    EEG = pop_select( EEG, 'rmchannel',{'HEO','VEO','EKG','EMG','TRIGGER','CB1','CB2'});
    %EEG = pop_reref( EEG, []); %全脑参考
    EEG = pop_reref( EEG, [33 43] );%双耳乳突参考
     event = EEG.event([EEG.event.type] == eventTypeNumber);
% 提取该事件的latency值
    latencyValue = event.latency;
% 显示latency值
    %disp(['Type ' num2str(eventTypeNumber) ' 的latency值是: ' num2str(latencyValue)]);

% 假设您已经加载了EEG数据集
% 假设您已经加载了EEG数据集

% 查找事件类型为21的事件
    event = EEG.event([EEG.event.type] == eventTypeNumber2);
% 提取该事件的latency值
    latencyValue2 = event.latency;
% 显示latency值
  % disp(['Type ' num2str(eventTypeNumber2) ' 的latency值是: ' num2str(latencyValue2)]);

    EEG = pop_select( EEG,'point',[latencyValue latencyValue2] );
    
    EEG = pop_eegfiltnew(EEG, 'locutoff',0.1,'plotfreqz',1);
    EEG = pop_eegfiltnew(EEG, 'hicutoff',40,'plotfreqz',1);
    EEG = pop_eegfiltnew(EEG, 'locutoff',48,'hicutoff',52,'revfilt',1,'plotfreqz',1); 
    %静息态分段 分成两秒一段
    EEG = eeg_regepochs(EEG, 'recurrence', 2, 'limits',[0 2], 'rmbase',NaN);
    %假设分成m秒一段，
    %EEG = eeg_regepochs(EEG, 'recurrence', m, 'limits',[0 m], 'rmbase',NaN);
    %重新绘制eeglab视窗
    eeglab redraw
    
    EEG = pop_saveset( EEG, 'filename',strcat(group_files(i).name(1:end-4), '.set'), 'filepath',strcat(group_dir, filesep, '_step1'));   %注意需要在运行代码之前，文件目录下建一个_resam_remch的文件夹，以下雷同
end%%
%%
%%


%%
%% 这部分是选择多段实验一起的eeg 可以先整体做set然后再分段
clc;clear;
%% Specify Basic information of different groups     将cdt文件转换成set文件，降采样为500，去除无用电极，保存。
group_dir = 'path_to_your_data';     % 此处路径需要设置为自己的文件目录
group_files = dir([group_dir, filesep, '*.dap']);  
eventTypeNumber = 25; % 例如21, 22, 23等

% 您想要提取的事件类型的数值
eventTypeNumber2 = 30; % 例如21, 22, 23等

%filesep是\的意思
for i=1:length(group_files)
    subj_fn = group_files(i).name;
    EEG = loadcurry(strcat(group_dir, filesep, subj_fn), 'CurryLocations', 'False');    %导入原始数据
    EEG = pop_resample( EEG, 250);   %降采样
    EEG = pop_select( EEG, 'rmchannel',{'HEO','VEO','EKG','EMG','TRIGGER','CB1','CB2'});
    %EEG = pop_reref( EEG, []); %全脑参考
    EEG = pop_reref( EEG, [33 43] );%双耳乳突参考
    event = EEG.event([EEG.event.type] == eventTypeNumber);
% 提取该事件的latency值
    latencyValue = event.latency;
% 显示latency值
    %disp(['Type ' num2str(eventTypeNumber) ' 的latency值是: ' num2str(latencyValue)]);

% 查找事件类型为21的事件
    event = EEG.event([EEG.event.type] == eventTypeNumber2);
% 提取该事件的latency值
    latencyValue2 = event.latency;
% 显示latency值
  % disp(['Type ' num2str(eventTypeNumber2) ' 的latency值是: ' num2str(latencyValue2)]);

    EEG = pop_select( EEG,'point',[latencyValue latencyValue2] );
    
    EEG = pop_eegfiltnew(EEG, 'locutoff',0.1,'plotfreqz',1);
    EEG = pop_eegfiltnew(EEG, 'hicutoff',40,'plotfreqz',1);
    EEG = pop_eegfiltnew(EEG, 'locutoff',48,'hicutoff',52,'revfilt',1,'plotfreqz',1); 
    %静息态分段 分成两秒一段
    %EEG = eeg_regepochs(EEG, 'recurrence', 2, 'limits',[0 2], 'rmbase',NaN);
    %假设分成m秒一段，
    %EEG = eeg_regepochs(EEG, 'recurrence', m, 'limits',[0 m], 'rmbase',NaN);
    %重新绘制eeglab视窗
    eeglab redraw
    
    EEG = pop_saveset( EEG, 'filename',strcat(group_files(i).name(1:end-4), '.set'), 'filepath',strcat(group_dir, filesep, '_step1'));   %注意需要在运行代码之前，文件目录下建一个_resam_remch的文件夹，以下雷同
end
%%


%% 上面是整个预处理 后面进行对选择的mark进行分段提取
clc;clear;
group1_dir = 'path_to_your_data';     % 此处路径需要设置为自己的文件目录
group1_dir1 = 'path_to_your_data';     % 此处路径需要设置为自己的文件目录
group1_files = dir([group1_dir1, filesep, '*.set']);  %filesep是\的意思
eventTypeNumber = 25; % 例如21, 22, 23等

% 您想要提取的事件类型的数值
eventTypeNumber2 = 27; % 例如21, 

for i=1:length(group1_files)
    subj_fn = group1_files(i).name;
    EEG = pop_loadset('filename',strcat(subj_fn(1:end-4), '.set'), 'filepath', strcat(group1_dir, filesep, '_step1')); %导入数据
    
    eventType1 = num2str(eventTypeNumber);
    eventType2 = num2str(eventTypeNumber2);
    
    % 找到对应type的索引
    idx1 = find(strcmp({EEG.event.type}, eventType1));
    idx2 = find(strcmp({EEG.event.type}, eventType2));
    latencyValue = EEG.event(idx1(1)).latency;
    latencyValue2 = EEG.event(idx2(1)).latency;
    
    %event = EEG.event([EEG.event.type] == eventTypeNumber);
% 提取该事件的latency值
    %latencyValue = event.latency;
% 显示latency值
    %disp(['Type ' num2str(eventTypeNumber) ' 的latency值是: ' num2str(latencyValue)]);

% 假设您已经加载了EEG数据集
% 假设您已经加载了EEG数据集

% 查找事件类型为21的事件
    %event = EEG.event([EEG.event.type] == eventTypeNumber2);
% 提取该事件的latency值
    %latencyValue2 = event.latency;
% 显示latency值
  % disp(['Type ' num2str(eventTypeNumber2) ' 的latency值是: ' num2str(latencyValue2)]);

    EEG = pop_select( EEG,'point',[latencyValue latencyValue2] );
    %静息态分段 分成两秒一段
    EEG = eeg_regepochs(EEG, 'recurrence', 2, 'limits',[0 2], 'rmbase',NaN);   % 跑ICA
    EEG = pop_saveset( EEG, 'filename',strcat(group1_files(i).name(1:end-4), '.set'), 'filepath',strcat(group1_dir, filesep, '_step2'));  %文件夹可多任务断保存数据
end



%%  第四部分，基于curry设备EEG数据预处理批处理，这一部分因为我电脑处理器原因用的是2020版MATLAB编写
%% 数据格式格式转化
clc;clear;
%% Specify Basic information of different groups     将cdt文件转换成set文件，降采样为500，去除无用电极，保存。
group_dir = 'path_to_your_data';     % 此处路径需要设置为自己的文件目录
group_files = dir([group_dir, filesep, '*.dap']);  %filesep是\的意思
for i=1:length(group_files)
    subj_fn = group_files(i).name;
    EEG = loadcurry(strcat(group_dir, filesep, subj_fn), 'CurryLocations', 'False');    %导入原始数据
    EEG = pop_resample( EEG, 250);   %降采样
    EEG = pop_select( EEG, 'rmchannel',{'HEO','VEO','EKG','EMG','TRIGGER','CB1','CB2'});
    EEG = pop_reref( EEG, []); %全脑参考
    %EEG = pop_reref( EEG, [33 43] );%双耳乳突参考
    EEG = pop_eegfiltnew(EEG, 'locutoff',0.1,'plotfreqz',1);
    EEG = pop_eegfiltnew(EEG, 'hicutoff',40,'plotfreqz',1);
    EEG = pop_eegfiltnew(EEG, 'locutoff',48,'hicutoff',52,'revfilt',1,'plotfreqz',1); 
    %静息态分段 分成两秒一段
    EEG = eeg_regepochs(EEG, 'recurrence', 2, 'limits',[0 2], 'rmbase',NaN);
    %假设分成m秒一段，
    %EEG = eeg_regepochs(EEG, 'recurrence', m, 'limits',[0 m], 'rmbase',NaN);
    %重新绘制eeglab视窗
    eeglab redraw
    
    EEG = pop_saveset( EEG, 'filename',strcat(group_files(i).name(1:end-4), '.set'), 'filepath',strcat(group_dir, filesep, '_step1'));   %注意需要在运行代码之前，文件目录下建一个_resam_remch的文件夹，以下雷同
end

%%文件夹_step1用于存放我们经过降采样、带通滤波、陷波滤波、去除无关电极之后的数据；
%文件夹_preica用于存放我们经过检查有无坏段坏导之后的数据；
%文件夹_ica用于存放我们跑完ICA之后的数据；
%文件夹_rm_ica用于存放ICLabel自动去完伪迹以及全脑平均重参考之后的数据；
%母文件夹demo下存放原始采集的数据。
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
    EEG = pop_icflag(EEG, [NaN NaN;0.5 1;0.5 1;NaN NaN;NaN NaN;NaN NaN;NaN NaN]); % 标记伪迹成分。这里可以自定义设定阈值，依次为Brain, Muscle, Eye, Heart, Line Noise, Channel Noise, Other.
    EEG = pop_subcomp( EEG, [], 0)   %去除上述伪迹成分
    %EEG = pop_reref( EEG, []);    %全脑平均重参考
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',strcat(group1_files(i).name(1:end-4), '.set'), 'filepath',strcat(group1_dir, filesep, '_rm_ica')); 
end


%% 
clear; clc; close all;

%% 1. 基础参数设置
% 定义路径
group1_dir = 'path_to_your_data';
group2_dir = 'path_to_your_data';

group1_files = dir([group1_dir, filesep, '*.set']);
group2_files = dir([group2_dir, filesep, '*.set']);

% 定义采样率和长度 (根据之前的预处理设定)
Fs = 250;
L = 500; 

% FFT 参数
NFFT = 2^nextpow2(L);
f = linspace(0, Fs/2, NFFT/2+1);

%% 2. 循环读取数据并进行 FFT
% === Group 1 ===
disp('正在读取 Group 1 数据...');
for i = 1:length(group1_files)
    subj_fn = group1_files(i).name;
    EEG = pop_loadset([group1_dir, filesep, subj_fn]);
    
    % 维度: Channel * Time * Epochs
    for ii = 1:size(EEG.data, 1) % Channel
        for jj = 1:size(EEG.data, 3) % Epochs
            y = squeeze(EEG.data(ii, :, jj));
            temp = fft(y, NFFT);
            % 计算单边频谱
            Y(jj, :) = abs(temp(1:NFFT/2+1)) * 2 / L;
        end
        % 对分段取平均 -> Sub * Ch * Freq
        group1_FFT_power(i, ii, :) = squeeze(mean(Y, 1)); 
        clear Y;
    end
end

% === Group 2 ===
disp('正在读取 Group 2 数据...');
for i = 1:length(group2_files)
    subj_fn = group2_files(i).name;
    EEG = pop_loadset([group2_dir, filesep, subj_fn]);
    
    for ii = 1:size(EEG.data, 1)
        for jj = 1:size(EEG.data, 3)
            y = squeeze(EEG.data(ii, :, jj));
            temp = fft(y, NFFT);
            Y(jj, :) = abs(temp(1:NFFT/2+1)) * 2 / L;
        end
        group2_FFT_power(i, ii, :) = squeeze(mean(Y, 1)); 
        clear Y;
    end
end

%% 3. 高级绘图：带标准误阴影的频谱图 (Publication Quality)
% 定义感兴趣的电极和频率范围
target_ch_idx = 28; % 假设是 Cz (第28通道)
target_ch_name = EEG.chanlocs(target_ch_idx).labels; % 自动获取电极名称
f_idx = find(f <= 30);
f_plot = f(f_idx);

% 提取绘图数据 (Sub * Freq)
g1_plot_data = squeeze(group1_FFT_power(:, target_ch_idx, f_idx));
g2_plot_data = squeeze(group2_FFT_power(:, target_ch_idx, f_idx));

% 计算均值 (Mean) 和 标准误 (SEM)
mean1 = mean(g1_plot_data, 1);
sem1 = std(g1_plot_data, 0, 1) / sqrt(size(g1_plot_data, 1));

mean2 = mean(g2_plot_data, 1);
sem2 = std(g2_plot_data, 0, 1) / sqrt(size(g2_plot_data, 1));

figure('Color','w'); % 白色背景
% 绘制 Group 1 阴影
fill([f_plot, fliplr(f_plot)], [mean1-sem1, fliplr(mean1+sem1)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% 绘制 Group 2 阴影
fill([f_plot, fliplr(f_plot)], [mean2-sem2, fliplr(mean2+sem2)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% 绘制均值线
h1 = plot(f_plot, mean1, 'r', 'LineWidth', 2);
h2 = plot(f_plot, mean2, 'b', 'LineWidth', 2);

% 美化
title(['Power Spectrum at ' target_ch_name], 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Frequency (Hz)', 'FontSize', 12);
ylabel('Power (\muV)', 'FontSize', 12);
legend([h1, h2], {'Group 1', 'Group 2'}, 'Location', 'northeast');
xlim([0 30]); box off;
grid on;

%% 4. 高级绘图：地形图 + 统计显著性标记
% 定义感兴趣频段 (Alpha)
freq_band = [8 12];
f_ROI = find((f >= freq_band(1)) & (f <= freq_band(2)));

% 提取该频段均值 -> Sub * Ch
g1_topo = squeeze(mean(group1_FFT_power(:, :, f_ROI), 3));
g2_topo = squeeze(mean(group2_FFT_power(:, :, f_ROI), 3));

% 组平均 (Ch * 1)
g1_mean_topo = mean(g1_topo, 1);
g2_mean_topo = mean(g2_topo, 1);

% === 统计检验 (独立样本 T 检验) ===
pvals = zeros(1, EEG.nbchan);
tvals = zeros(1, EEG.nbchan);
for ch = 1:EEG.nbchan
    [~, p, ~, stats] = ttest2(g1_topo(:, ch), g2_topo(:, ch));
    pvals(ch) = p;
    tvals(ch) = stats.tstat;
end

% 统一 Colorbar 范围 (美观)
clim_max = max([g1_mean_topo, g2_mean_topo]);
clim_min = min([g1_mean_topo, g2_mean_topo]);

figure('Color','w', 'Position', [100, 100, 1000, 300]);
% 图1: Group 1
subplot(1,3,1); 
topoplot(g1_mean_topo, EEG.chanlocs, 'maplimits', [clim_min, clim_max]); 
title('Group 1 (Alpha)', 'FontSize', 12); colorbar;

% 图2: Group 2
subplot(1,3,2); 
topoplot(g2_mean_topo, EEG.chanlocs, 'maplimits', [clim_min, clim_max]); 
title('Group 2 (Alpha)', 'FontSize', 12); colorbar;

% 图3: T-values (差异图) + 显著点
subplot(1,3,3); 
% 'pmask' 参数会自动在 p<0.05 的位置画 'x' 或 '*'
topoplot(tvals, EEG.chanlocs, 'maplimits', [-3 3], 'pmask', pvals < 0.05); 
title('T-Test (p<0.05 masked)', 'FontSize', 12); 
colorbar;
colormap(subplot(1,3,3), jet); % T值图通常用 jet 颜色

%% 5. 提取数据供 SPSS 使用
% 提取指定电极、指定频段的数据
target_elec = 28; % Cz
g1_spss = mean(group1_FFT_power(:, target_elec, f_ROI), 3); % 均值
g2_spss = mean(group2_FFT_power(:, target_elec, f_ROI), 3);

% 构建长格式数据 (Long Format)
% 第一列: Power, 第二列: Group Label (1 or 2)
data_all = [g1_spss; g2_spss];
group_labels = [ones(length(g1_spss), 1); 2*ones(length(g2_spss), 1)];

% 保存为 CSV
T = table(data_all, group_labels, 'VariableNames', {'AlphaPower', 'Group'});
writetable(T, 'EEG_Data_For_SPSS_Independent.csv');

fprintf('分析完成！\n');
fprintf('1. 频谱图已生成 (带误差阴影)。\n');
fprintf('2. 地形图已生成 (包含 Group 1, Group 2 和 T检验显著性图)。\n');
fprintf('3. SPSS 数据已保存为 EEG_Data_For_SPSS_Independent.csv\n');