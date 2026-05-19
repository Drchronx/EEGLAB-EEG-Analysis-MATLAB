
%第三部分时频分析，可以结合时域结果一起补充撰写，可解释数据多一维度，有时候有新发现新惊喜。
clear all; clc;
data_path = 'path_to_your_data';%定义预处理完的数据路径，我这就不该了按自己数据把路径更改
cd(data_path)
files = dir('*.set');
fn = {files.name};
%定义条件信息
Cond={  '1'  '2'  '3'  '4'  };

%% time-frequency analysis for multiple conditions

% 对于所有被试 
for i=1:length(fn)
    setname = fn{i}; %% filename of set file
    EEG = pop_loadset('filename',setname,'filepath',data_path); 
    EEG= eeg_checkset( EEG );
    for j=1:length(Cond)
        %已重新分段的方式提取j条件的数据
        EEG_new = pop_epoch( EEG, Cond(j), [-1  2], 'newname', 'Merged datasets pruned with ICA   epochs epochs', 'epochinfo', 'yes'); 
        EEG_new = eeg_checkset( EEG_new );
        EEG_new = pop_rmbase( EEG_new, [-1000     0]); 
        EEG_new = eeg_checkset( EEG_new );
        %对每个通道
        for nchan=1:size(EEG_new.data,1)
            %提取第nchan通道下的所有时间分段的数据
            %时间 *分段
            x = squeeze(EEG_new.data(nchan,:,:)); 
            %定义时间轴信息
            xtimes=EEG_new.times/1000; 
            %定义在什么时间点加窗
            t=EEG_new.times/1000;
            %定义频率轴
            f=1:1:30; 
            %定义采样率
            Fs = EEG.srate;
            %定义时间窗
            winsize = 0.400; 
            %进行STFT
            %
            [S, P, F, U] = sub_stft(x, xtimes, t, f, Fs, winsize); 
            % 被试 * 条件 * 通道 * 频率 * 时间点
            %进行个体水平叠加平均（分段）
            %P_DATA 5维  被试*通道*频率*时间
            P_data(i,j,nchan,:,:)=squeeze(mean(P,3)); %%P_data (without baseline correciton):  subj*cond*chan*f*time
        end
    end    
end
%%
%save('P_data.mat','P_data')
%load('P_data.mat');
%% baseline correction 
%p_data 5维 被试、条件、通道、频率、时间
% 对于所有被试、条件、通道、频率，都做基线校正
%定义基线范围
t_pre_idx=find((t>=-0.8)&(t<=-0.2));
for i=1:size(P_data,1)
    %对每个被试
    for j=1:size(P_data,2)
        %对于每个通道
        for ii=1:size(P_data,3)
            %对于每个频率点
            for jj=1:size(P_data,4)
                %对于时间序列
                temp_data=squeeze(P_data(i,j,ii,jj,:));
                %进行基线校正
                P_BC(i,j,ii,jj,:)=temp_data-mean(temp_data(t_pre_idx));
            end
        end
    end
end
%%
save('P_BC.mat','P_BC')
%load('P_BC.mat');

%% %%单因素找显著


% 提取感兴趣的通道（例如第19通道），得到 data_test: 被试 x 条件(2) x 频率 x 时间
data_test = squeeze(P_BC(:,:,19,:,:)); 

% 初始化结果矩阵，避免动态增长变量导致变慢
P_ttest = zeros(size(data_test,3), size(data_test,4));
T_ttest = zeros(size(data_test,3), size(data_test,4));

% 对于每个频率点
for i = 1:size(data_test,3)
    % 对于每个时间点
    for j = 1:size(data_test,4)
        % 分别提取两个条件ij时频点下的数据
        % 【修改点】：这里必须对应你的 Cond 定义，只有 1 和 2
        data_1 = squeeze(data_test(:, 1, i, j)); % 提取第1个条件 (对应 Cond{'1'})
        data_2 = squeeze(data_test(:, 2, i, j)); % 提取第2个条件 (对应 Cond{'2'})
        
        % 配对样本t检验
        [h, p, ci, stats] = ttest(data_1, data_2); 
        P_ttest(i,j) = p; 
        T_ttest(i,j) = stats.tstat; 
    end
end




%% 方案一：绘制 T 值图并叠加显著性轮廓 (推荐)
figure;
% 1. 绘制背景：使用 T 值 (T_ttest) 而不是 P 值
% T值能显示方向性：正值(暖色)表示 Cond1 > Cond2，负值(冷色)表示 Cond1 < Cond2
imagesc(xtimes, f, T_ttest); 
hold on;
axis xy; % 翻转Y轴，让低频在下，高频在上

% 2. 叠加轮廓线：将 P < 0.05 的区域圈出来
% contour 的参数 [0.05 0.05] 表示只画数值为 0.05 的等高线
[M, c] = contour(xtimes, f, P_ttest, [0.05 0.05], 'LineColor', 'k', 'LineWidth', 2);

% 3. 美化设置
% 设置色标范围：为了让0点对应白色/绿色，通常取绝对值的最大值做对称范围
max_t = max(abs(T_ttest(:))); 
caxis([-max_t max_t]); 

% 使用红蓝配色 (如果没有 redblue 函数，可以用 jet)
colormap(jet); 
colorbar;
title('Time-Frequency Differences (T-values with p<0.05 contour)');
xlabel('Time (s)');
ylabel('Frequency (Hz)');

% 加一条 0时刻的竖线
xline(0, '--k', 'LineWidth', 1.5);



%%
%%
%%
%% 重新导入的情况下1. 数据准备 (Data Preparation)
% 定义时间轴 (已知的)
xtimes = -1000:2:1998; 

% 定义频率轴 (【必须修改】这里假设是 1-30Hz)
% ------------------------------------------------------
start_freq = 1;   % 分析的起始频率
end_freq = 30;    % 分析的结束频率
% ------------------------------------------------------
num_freqs = size(T_ttest, 1); % 自动获取行数
f = linspace(start_freq, end_freq, num_freqs); % 生成 f 向量

%% 2. 维度检查 (安全检查)
% 确保 X轴 和 Y轴 的长度分别对应 数据的 列数 和 行数
if length(xtimes) ~= size(T_ttest, 2) || length(f) ~= size(T_ttest, 1)
    error('维度不匹配！\n f 的长度是 %d, 数据行数是 %d \n xtimes 的长度是 %d, 数据列数是 %d',...
        length(f), size(T_ttest, 1), length(xtimes), size(T_ttest, 2));
end

%% 3. 绘图 (Plotting)
figure;
% 绘制 T 值背景
imagesc(xtimes, f, T_ttest); 
hold on;
axis xy; % 翻转Y轴，低频在下

% 叠加 P < 0.05 轮廓
% 注意：P_ttest 必须和 T_ttest 维度完全一样
[M, c] = contour(xtimes, f, P_ttest, [0.05 0.05], 'LineColor', 'k', 'LineWidth', 2);

% 设置颜色范围 (对称显示)
max_t = max(abs(T_ttest(:)));
if max_t == 0, max_t = 1; end
caxis([-max_t max_t]); 

% 配色与标注
colormap(jet); 
colorbar;
title('Time-Frequency T-values (p<0.05 outlined)');
xlabel('Time (ms)');
ylabel('Frequency (Hz)');

% 0时刻竖线
xline(0, '--k', 'LineWidth', 1.5);
%% 
%% 
%%



%% 提取数据  补充的1. 补全缺失的变量 t 和 f
%% 

% 定义时间 (毫秒)，并转化为 秒 (s) 以匹配你的 ROI [0.5 1.0]
% 你的数据是 -1000ms 到 1998ms，间隔 2ms
t_ms = -1000:2:1998; 
t = t_ms / 1000;  % 【关键】转换成秒！也就是 -1.0s 到 1.998s

% 定义频率 f (假设是1-30Hz，需根据你 data_test 的第3维度修改)
% 请务必确认 size(data_test, 3) 是多少，如果是 30，就写 1:30
f = 1:30;   % <--- 如果你的 f 还没定义，请取消注释并修改它

%% 2. 你的提取代码 (现在可以运行了)

data_test = squeeze(P_BC(:,:,20,:,:)); 
% 此时 data_test 应该是：被试(N) × 条件(2) × 频率(F) × 时间(T)

ROI1_t = [0.20 0.30]; % 时间窗：0.5秒 到 1.0秒
ROI1_f = [4 7];       % 频率窗：Theta波 4-8Hz

% 查找对应的索引
ROI1_t_idx = find((t >= ROI1_t(1)) & (t <= ROI1_t(2)));
ROI1_f_idx = find((f >= ROI1_f(1)) & (f <= ROI1_f(2)));

% 检查一下是否找到了索引 (防止空值)
if isempty(ROI1_t_idx) || isempty(ROI1_f_idx)
    error('未找到对应的时间或频率索引！请检查 t 和 f 的范围是否正确。');
end

% 提取并平均
% dim3 是频率，dim4 是时间 -> 两次 mean 变成一个数值
data_test1 = squeeze(mean(mean(data_test(:,:,ROI1_f_idx, ROI1_t_idx), 3), 4));

% 结果 data_test1 应该是一个 [被试 × 条件] 的矩阵
disp('数据提取完成，data_test1 尺寸为：');
disp(size(data_test1));



%% 画图
%% 0. 准备工作：选择通道和提取数据
% 假设你关注的是第19通道 (和前面统计检验保持一致)
chan_idx = 19; 

% 计算群组平均 (Grand Average)
% P_BC 维度: 被试 x 条件 x 通道 x 频率 x 时间
% 对第1维(被试)求平均，squeeze后得到: 频率 x 时间
GrandAvg_Cond1 = squeeze(mean(P_BC(:, 1, chan_idx, :, :), 1)); 
GrandAvg_Cond2 = squeeze(mean(P_BC(:, 2, chan_idx, :, :), 1)); 

%% 1. 计算统一的色标范围 (关键步骤)
% 找出两个图中数据的最大值和最小值，确保两张图颜色代表的数值一致
all_values = [GrandAvg_Cond1(:); GrandAvg_Cond2(:)];
% 方案A：自动范围 (从最小到最大)
c_lims = [min(all_values) max(all_values)]; 

% 方案B：对称范围 (推荐，如果你的数据是关于0对称的，比如ERD/ERS)
% limit = max(abs(all_values));
% c_lims = [-limit limit]; 

%% 2. 绘图
figure('Color', 'w', 'Position', [100, 100, 1200, 400]); % 设置一个宽一点的画布

% --- 绘制条件 1 ---
subplot(1, 2, 1);
imagesc(xtimes, f, GrandAvg_Cond1); 
axis xy;                % 翻转Y轴，让低频在下
caxis(c_lims);          % 【重要】锁定色标范围
title(['Condition 1 (Channel ' num2str(chan_idx) ')']);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
xline(0, '--k', 'LineWidth', 1.5); % 0时刻竖线
colorbar;               % 显示色条
colormap(jet);          % 设定色系

% --- 绘制条件 2 ---
subplot(1, 2, 2);
imagesc(xtimes, f, GrandAvg_Cond2); 
axis xy;
caxis(c_lims);          % 【重要】使用和左图一样的范围
title(['Condition 2 (Channel ' num2str(chan_idx) ')']);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
xline(0, '--k', 'LineWidth', 1.5);
colorbar;
colormap(jet);

%% 3. (可选) 保存图片
% saveas(gcf, 'TF_Comparison_Cond1_vs_Cond2.png');
