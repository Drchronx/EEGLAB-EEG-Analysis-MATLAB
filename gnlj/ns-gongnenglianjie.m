%%  ����̬��������



%%��һ���֣�����curry�豸Ԥ�������������һ������Ϊ�ҵ��Դ�����ԭ���õ���2020��MATLAB��д
%% ���ݸ�ʽ��ʽת��
clc;clear;
%% Specify Basic information of different groups     ��cdt�ļ�ת����set�ļ���������Ϊ500��ȥ�����õ缫�����档
group_dir = 'path_to_your_data';     % �˴�·����Ҫ����Ϊ�Լ����ļ�Ŀ¼
group_files = dir([group_dir, filesep, '*.dap']);  %filesep��\����˼
for i=1:length(group_files)
    subj_fn = group_files(i).name;
    EEG = loadcurry(strcat(group_dir, filesep, subj_fn), 'CurryLocations', 'False');    %����ԭʼ����
    EEG = pop_resample( EEG, 500);   %������
    EEG = pop_select( EEG, 'rmchannel',{'HEO','VEO','EKG','EMG','TRIGGER','CB1','CB2'});
    EEG = pop_reref( EEG, [33 43] );
    EEG = pop_epoch( EEG, {  '1'  '2'  '3'  '4' }, [-1  2], 'newname', 'Neuroscan Curry file resampled epochs', 'epochinfo', 'yes');
    EEG = pop_rmbase( EEG, [-1000 0] ,[]);
    EEG = pop_eegfiltnew(EEG, 'locutoff',0.1,'plotfreqz',1);
    EEG = pop_eegfiltnew(EEG, 'hicutoff',40,'plotfreqz',1);
    EEG = pop_eegfiltnew(EEG, 'locutoff',48,'hicutoff',52,'revfilt',1,'plotfreqz',1); 
    EEG = pop_saveset( EEG, 'filename',strcat(group_files(i).name(1:end-4), '.set'), 'filepath',strcat(group_dir, filesep, '_step1'));   %ע����Ҫ�����д���֮ǰ���ļ�Ŀ¼�½�һ��_resam_remch���ļ��У�������ͬ
end

%%�ļ���_step1���ڴ�����Ǿ�������������ͨ�˲����ݲ��˲���ȥ���޹ص缫֮������ݣ�
%�ļ���_preica���ڴ�����Ǿ���������޻��λ���֮������ݣ�
%�ļ���_ica���ڴ����������ICA֮������ݣ�
%�ļ���_rm_ica���ڴ��ICLabel�Զ�ȥ��α���Լ�ȫ��ƽ���زο�֮������ݣ�
%ĸ�ļ���demo�´��ԭʼ�ɼ������ݡ�
%% �˵���������eeglab����_step1���ļ��������޻������Σ����л������ֵ�����л�����ɾ�������֮�󱣴����ݵ�_preica�ļ��У��ļ����Ʊ��ֲ��䡣

%%�ļ����룺��eeglab�����File->load existing dataset��
%�ļ����棺���File->save current dataset as��

%%   ���� ICA
clc;clear;
group1_dir = 'path_to_your_data';     % �˴�·����Ҫ����Ϊ�Լ����ļ�Ŀ¼
group1_dir1 = 'path_to_your_data';     % �˴�·����Ҫ����Ϊ�Լ����ļ�Ŀ¼
group1_files = dir([group1_dir1, filesep, '*.set']);  %filesep��\����˼
for i=1:length(group1_files)
    subj_fn = group1_files(i).name;
    EEG = pop_loadset('filename',strcat(subj_fn(1:end-4), '.set'), 'filepath', strcat(group1_dir, filesep, '_preica')); %��������
    EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'interrupt','on');   % ��ICA
    EEG = pop_saveset( EEG, 'filename',strcat(group1_files(i).name(1:end-4), '.set'), 'filepath',strcat(group1_dir, filesep, '_ica'));  %��������
    
end
%% ʹ��ICLabel�Զ�ȥ��ICA�ɷ�
clc;clear;
group1_dir ='path_to_your_data';        % �˴�·����Ҫ����Ϊ�Լ����ļ�Ŀ¼
group1_dir2 = 'path_to_your_data'; 
group1_files = dir([group1_dir2, filesep, '*.set']);  %filesep��\����˼
for i=1:length(group1_files)
    subj_fn = group1_files(i).name;
    EEG = pop_loadset('filename',strcat(subj_fn(1:end-4), '.set'), 'filepath', group1_dir2);
    EEG = pop_iclabel(EEG, 'default');
    EEG = pop_icflag(EEG, [NaN NaN;0.5 1;0.5 1;NaN NaN;NaN NaN;NaN NaN;NaN NaN]); % ���α���ɷ֡���������Զ����趨��ֵ������ΪBrain, Muscle, Eye, Heart, Line Noise, Channel Noise, Other.
    EEG = pop_subcomp( EEG, [], 0)   %ȥ������α���ɷ�
    %EEG = pop_reref( EEG, []);    %ȫ��ƽ���زο�
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',strcat(group1_files(i).name(1:end-4), '.set'), 'filepath',strcat(group1_dir, filesep, '_rm_ica')); 
end

%% �ֶ�������� ȷ��Ԥ����֮��������Ǹɾ���

%%
%%
%%
%%  �ڶ��������ݵ��������з��ļ���

clear;clc
%�����������ڵ�·��
data_path = 'path_to_your_data';
%�������ݵı���·��
save_path = 'path_to_your_data';
%���������ڵ�·������Ϊ����·��
cd(data_path)
%ɸѡ��ǰ·�������е�vhdr��β���ļ�
files = dir('*.set');
%��ȡ�ļ���
fn = {files.name};
%����ÿ������
for i = 1:length(fn)
    %����set��ʽ������
    EEG = pop_loadset('filename',fn{i},'filepath',data_path);
    EEG = eeg_checkset( EEG ); 
    %EEG = pop_select( EEG,'nochannel',{'AF3' 'AF4' 'F7' 'F8' 'FT7' 'FT8' 'T7' 'T8' 'TP7' 'TP8' 'P7' 'P8' 'PO7' 'PO8'});
    EEG = eeg_checkset( EEG );
    EEG = pop_resample( EEG, 250); 
    EEG = pop_epoch( EEG, {  '4'  }, [-1  2], 'newname', 'Neuroscan Curry file resampled epochs', 'epochinfo', 'yes');
    EEG = pop_rmbase( EEG, [-1000 0] ,[]);
    %����set��ʽ������
    %EEG.subject = fn{1,i}(1:4);
    %EEG.group = '0';
    %EEG.condition = '22';
    %original_subject_name = fn{i}(1:4); 
    % ����EEG.subjectΪ�µ�����
    %EEG.subject = [original_subject_name '_' EEG.group '_' EEG.condition];
    EEG = pop_saveset( EEG, 'filename',fn{i},'filepath',save_path);
end



%%
%%
%% ������ csdת���������ɲ�����
%compute_CSDdata.m����ÿ�����ݵ�CSD����ѡ����
%���д˽ű�ǰ����Ҫ��CSDtoolbox���ؽ�MATLAB·�������д˽ű��󣬻ᵯ�������Ի�����ֱ��������Ի���ָ������ת����EEG�������ڵ��ļ��С��͡�ת�������ɵ�CSD�������ڵ�·����



%%
%% ���Ĳ������ݽ��м���

%%PS_computation.m
% �ýű��ɼ����¼�������ݵ�PLV/PLI/wPLI��
%���д˽ű�ǰ������ÿ������������𣩸�һ���ļ��С�ÿ���ļ��д洢��Ӧ����������µ����б��Ե����ݣ���ʾ��������ʾ
%�˽ű�����PLV/PLI/wPLI�ǻ���STFT�ģ������Ҫsub_stft.m����Ҫȷ���ú�����MATLAB��current folder��������·����

clc;clear all; close all;

%% ָ�������Ϣ
Data_Dir = uigetdir([],'Path of the EEG datasets'); % ָ���Ե��������ڵ�·��
Output_Dir = uigetdir([],'Path to store the measures'); % ָ�����ɵ�ָ�걣���·��

f_roi = inputdlg('the frequencies of interest');% STFT������Ƶ�ʵ�
f_roi = str2num(f_roi{1});

winsize = inputdlg('the winsize of STFT (in seconds)'); % STFT������ʱ�䴰
winsize = str2num(winsize{1});

prefix = inputdlg('the prefix of computed measures'); % ��Ҫ�����ָ�������ǰ׺���������������������
prefix = prefix{1};
%% ����PLV/PLI/wPLI
Dir_Data = dir(fullfile(Data_Dir,'*.set'));  
filenames = {Dir_Data.name};

for subj = 1:length(filenames)
    subj
    EEG = pop_loadset('filename',filenames{1, subj},'filepath',Data_Dir);
    EEG = eeg_checkset( EEG );
    % STFT
    Fs = EEG.srate;
    xtimes = EEG.times/1000;
    t = EEG.times/1000;   
    
    for nch = 1:EEG.nbchan
       [S, P, F, U] = sub_stft(squeeze(EEG.data(nch,:,:)), xtimes, t, f_roi, Fs, winsize);   %��ĳ�����Ե��������ͨ������STFT 
       S_subject(:,:,:,nch) = S; % ���ñ���STFT��Ľ����źŸ�ֵ��S_subject ά����Ƶ��*ʱ��*�Դ�*ͨ��
       clear P S F U
    end    
    %
    angle_subject = angle(S_subject); % ��ȡ��λ   
    for x = 1:EEG.nbchan
        for y = 1:EEG.nbchan
            if x >= y % ��PLV PLI��wPLI��Ϊ����ָ�꣬������ֻ�Ǽ���缫*�缫�����Ӿ����������
                % PLV & PLI 
                tempx = squeeze(angle_subject(:,:,:,x));
                tempy = squeeze(angle_subject(:,:,:,y));
                relative_phase = tempx - tempy;
                plv(:,:,x,y,subj) = abs(sum(exp(1i*relative_phase),3)/size(relative_phase,3));
                pli(:,:,x,y,subj) = abs(mean(sign((abs(relative_phase)- pi).*relative_phase),3));
                % wPLI               
                tempx2 = squeeze(S_subject(:,:,:,x));
                tempy2 = squeeze(S_subject(:,:,:,y));
                crossspec = tempx2.* conj(tempy2);
                crossspec_imag = imag(crossspec);
                wpli(:,:,x,y,subj) = abs(mean(crossspec_imag,3))./mean(abs(crossspec_imag),3);
           end
        end
    end    
    clear angle_subject tempx tempy tempx2 tempy2 relative_phase S_subject crossspec crossspec_imag
    
    waitbar(subj/length(filenames))   
end

wpli(isnan(wpli)) = 0;

save(strcat(Output_Dir,'\',prefix,'_plv.mat'), 'plv', '-v7.3');
save(strcat(Output_Dir,'\',prefix,'_pli.mat'), 'pli', '-v7.3');
save(strcat(Output_Dir,'\',prefix,'_wpli.mat'), 'wpli', '-v7.3');



%%
%% ���岿�Խ������ͳ�Ʒ���  ֻ�Ƽ������Ƚ�

%��ѡ�����缫������ʱƵ��

%��ֱ��ʹ��tf_fc_test.m

%Ҳ��ʹ���������
clc;clear all;close all

%% ָ�������Ϣ
channel_pair = inputdlg('The channel pair to be test (e.g., [6 5])');
channel_pair = str2num(channel_pair{1});

time_axis = inputdlg('The time axis of TFR (in ms)');
time_axis = str2num(time_axis{1});

f_axis = inputdlg('The frequency axis of TFR (in Hz)');
f_axis = str2num(f_axis{1});

baseline = inputdlg('The baseline limits (in ms)');
baseline = str2num(baseline{1});

test_type = inputdlg('Test Type? Paired t-test = 1 Independent t-test = 0');
test_type = str2num(test_type{1});

name1 = inputdlg('the name of one condition (group)');
name1 = name1{1};

name2 = inputdlg('the name of the other condition (group)');
name2 = name2{1};


%% ��������
[filename1, pathname1, filterindex] = uigetfile('*.mat', 'Pick the mat file of measures of one condition (group)');
[filename2, pathname2, filterindex] = uigetfile('*.mat', 'Pick the mat file of measures of the other condition (group)');
A_fc = importdata(strcat(pathname1,filename1));
B_fc = importdata(strcat(pathname2,filename2));

%% ��ȡ����Ȥ�缫�Ե����ݣ������л��߽���
A_fc_x = squeeze(A_fc(:,:,channel_pair(1,1),channel_pair(1,2),:));
B_fc_x = squeeze(B_fc(:,:,channel_pair(1,1),channel_pair(1,2),:));
clear A_fc B_fc

A_fc_x = A_fc_x - repmat(mean(A_fc_x(:,time_axis >= baseline(1,1) & time_axis <= baseline(1,2),:),2),[1 length(time_axis) 1]);
B_fc_x = B_fc_x - repmat(mean(B_fc_x(:,time_axis >= baseline(1,1) & time_axis <= baseline(1,2),:),2),[1 length(time_axis) 1]);

%% �����ʱƵ�����t���飬��ʹ��FDR����
if test_type == 1 % paired
    for f = 1:size(A_fc_x,1)
        for t = 1:size(A_fc_x,2)
            [~,p(f,t)] = ttest(squeeze(A_fc_x(f,t,:)),squeeze(B_fc_x(f,t,:)));
        end
    end
else  % indpendent
    for f = 1:size(A_fc_x,1)
        for t = 1:size(A_fc_x,2)
            [~,p(f,t)] = ttest2(squeeze(A_fc_x(f,t,:)),squeeze(B_fc_x(f,t,:)));
        end
    end
end    
    
[~, p_masked] = fdr(p,0.05);

%% ��ͼ
A_fc_x_mean = mean(A_fc_x,3);
B_fc_x_mean = mean(B_fc_x,3);

figure;
subplot(2,2,1); imagesc(time_axis,f_axis,A_fc_x_mean);axis xy; colorbar;title(name1)
subplot(2,2,2); imagesc(time_axis,f_axis,B_fc_x_mean);axis xy; colorbar;title(name2)
subplot(2,2,3); imagesc(time_axis,f_axis,p);axis xy; colorbar;title('uncorrected p value')
subplot(2,2,4); imagesc(time_axis,f_axis,p_masked);axis xy; colorbar;title('significant area after FDR')

%pairs(idx_left,:)
% pairs(idx_right,:)



A_fc_x_mean = mean(A_fc_x,3);
B_fc_x_mean = mean(B_fc_x,3);

% --- �ؼ��޸ģ�Ϊp < 0.05����NaN���־��� ---

% 1. ΪδУ����pֵ�������־���ֻ����p < 0.05��ֵ��������ΪNaN��
p_uncorrected_plot = p;
p_uncorrected_plot(p >= 0.05) = NaN; 

% 2. ΪFDRУ�����pֵ�������־���
% �ϸ�������Ҫ��Ҳʹ��p < 0.05�������֣�
%p_FDR_plot = p_masked;
%p_FDR_plot(p >= 0.05) = NaN; 

% **�����顿**
% �����ϣ��p_maskedͼ����ʾ'p_masked'��'����'�����򣬶�ͨ��'p_masked'�Ѿ�ֻ�������������pֵ
% ��ô�������߼������ֿ����ǣ�
% p_FDR_plot = p_masked;
% p_FDR_plot(p_masked >= 0.05) = NaN; % ������p_masked�в������Ĳ��� (���p_masked�а���pֵ)
% 
% ���ߣ����p_masked��һ����Ԫ���� (�� 0/1)������Ҫ�������֣�ֱ�ӻ��Ƽ��ɡ�
% �����ұ������ϸ�����Ҫ������ַ�ʽ��

figure;
subplot(2,2,1); imagesc(time_axis,f_axis,A_fc_x_mean);axis xy; colorbar;title(name1)
subplot(2,2,2); imagesc(time_axis,f_axis,B_fc_x_mean);axis xy; colorbar;title(name2)
% ʹ���µ����־���
subplot(2,2,3); imagesc(time_axis,f_axis,p_uncorrected_plot);axis xy; colorbar;title('uncorrected p value (p < 0.05 only)')
% ʹ���µ����־���
subplot(2,2,4); imagesc(time_axis,f_axis,p_FDR_plot);axis xy; colorbar;title('significant area after FDR (p < 0.05 only)')





%ѡ�缫�Ժ��ҵ�����ʱƵ�������ȡ����

t_roi = inputdlg('The temporal limits of TF-ROI (in ms)');
t_roi = str2num(t_roi{1});

f_roi = inputdlg('The spectral limits of TF-ROI (in Hz)');
f_roi = str2num(f_roi{1});


A_fc_roxi = squeeze(mean(mean(A_fc_x(f_axis >= f_roi(1,1) & f_axis <= f_roi(1,2),time_axis >= t_roi(1,1) & time_axis <= t_roi(1,2),:,:,:),1),2));
B_fc_roxi = squeeze(mean(mean(B_fc_x(f_axis >= f_roi(1,1) & f_axis <= f_roi(1,2),time_axis >= t_roi(1,1) & time_axis <= t_roi(1,2),:,:,:),1),2));


%%
%% ��6���Խ������ͳ�Ʒ���  ֻ�Ƽ������Ƚ�

%��ѡʱƵ����ȷ�������缫

%��ֱ��ʹ��tf_roi_fc_test.m

%Ҳ��ʹ���������



clc;clear all;close all

%% ָ�������Ϣ
time_axis = inputdlg('The time axis of TFR (in ms)');
time_axis = str2num(time_axis{1});

f_axis = inputdlg('The frequency axis of TFR (in Hz)');
f_axis = str2num(f_axis{1});

baseline = inputdlg('The baseline limits (in ms)');
baseline = str2num(baseline{1});

t_roi = inputdlg('The temporal limits of TF-ROI (in ms)');
t_roi = str2num(t_roi{1});

f_roi = inputdlg('The spectral limits of TF-ROI (in Hz)');
f_roi = str2num(f_roi{1});

test_type = inputdlg('Test Type? Paired t-test = 1 Independent t-test = 0');
test_type = str2num(test_type{1});

name1 = inputdlg('the name of one condition (group)');
name1 = name1{1};

name2 = inputdlg('the name of the other condition (group)');
name2 = name2{1};

[filename3, pathname3, filterindex] = uigetfile('*.set', 'Pick an eeglab file used to computed the measures');

%% �������ݣ�����ȡTF-ROI��ָ�꣬����TF-ROI���ݱ�����current folder��
[filename1, pathname1, filterindex] = uigetfile('*.mat', 'Pick the mat file of measures of one condition (group)');
[filename2, pathname2, filterindex] = uigetfile('*.mat', 'Pick the mat file of measures of the other condition (group)');
A_fc = importdata(strcat(pathname1,filename1));
B_fc = importdata(strcat(pathname2,filename2));

A_fc = A_fc - repmat(mean(A_fc(:,time_axis >= baseline(1,1) & time_axis <= baseline(1,2),:,:,:),2),[1 length(time_axis) 1 1 1]);
B_fc = B_fc - repmat(mean(B_fc(:,time_axis >= baseline(1,1) & time_axis <= baseline(1,2),:,:,:),2),[1 length(time_axis) 1 1 1]);

A_fc_roi = squeeze(mean(mean(A_fc(f_axis >= f_roi(1,1) & f_axis <= f_roi(1,2),time_axis >= t_roi(1,1) & time_axis <= t_roi(1,2),:,:,:),1),2));
B_fc_roi = squeeze(mean(mean(B_fc(f_axis >= f_roi(1,1) & f_axis <= f_roi(1,2),time_axis >= t_roi(1,1) & time_axis <= t_roi(1,2),:,:,:),1),2));

fc_roi{1,1} = A_fc_roi;
fc_roi{1,2} = B_fc_roi;
save fc_roi.mat fc_roi  t_roi f_roi

%% ����ĳ��eeglab��ʽ�Ե����ݣ���ȡ��ʵ���Ե����ݵ������Ϣ
EEG = pop_loadset('filename',filename3,'filepath',pathname3);
chanlocs = EEG.chanlocs;
channel_number = size(EEG.data,1);
fc_number = channel_number*(channel_number-1)/2;

%% ȷ��ÿ���缫�Էֱ��Ӧ�������缫
pairs = nchoosek(1:channel_number,2); 
pairs = pairs(:,[2 1]); % ��PS_computation������������ǣ�����ֵ��������ֵ��������pairs�ǵ�һ��С�ڵڶ�����ֵ������ֵС�������֣��������ｫ��һ�к͵ڶ����û� 

%% �ֱ������������������𣩵���ƽ����ͨͼ
A_avg = mean(A_fc_roi,3);
A_avg2 = zeros(fc_number,1);
for i = 1:fc_number
    A_avg2(i,1) = A_avg(pairs(i,1),pairs(i,2));
end  

B_avg = mean(B_fc_roi,3);
B_avg2 = zeros(fc_number,1);
for i = 1:fc_number
    B_avg2(i,1) = B_avg(pairs(i,1),pairs(i,2));
end  

AB_avg2 = [A_avg2; B_avg2];

ds.chanPairs = pairs;  
ds.connectStrength = A_avg2; 
ds.connectStrengthLimits = [min(AB_avg2) max(AB_avg2)];
figure;title(name1);topoplot_connect(ds, chanlocs); 

ds.chanPairs = pairs;  
ds.connectStrength = B_avg2; 
ds.connectStrengthLimits = [min(AB_avg2) max(AB_avg2)];
figure;title(name2);topoplot_connect(ds, chanlocs);  
%% ��������������𣩼�ͳ�Ʒ���
for  i = 1:fc_number   
    A_fc_roi_2(i,:) = A_fc_roi(pairs(i,1),pairs(i,2),:);
    B_fc_roi_2(i,:) = B_fc_roi(pairs(i,1),pairs(i,2),:);
end

if test_type == 1 % paired
    for i = 1:fc_number 
       [~,p_right(i,1)] = ttest(A_fc_roi_2(i,:),B_fc_roi_2(i,:),0.05,'right');
       [~,p_left(i,1)] = ttest(A_fc_roi_2(i,:),B_fc_roi_2(i,:),0.05,'left');
    end
else  % indpendent
    for i = 1:fc_number 
       [~,p_right(i,1)] = ttest2(A_fc_roi_2(i,:),B_fc_roi_2(i,:),0.05,'right');
       [~,p_left(i,1)] = ttest2(A_fc_roi_2(i,:),B_fc_roi_2(i,:),0.05,'left');
    end
end
[~,h_right] = fdr(p_right,0.05);
idx_right = find(h_right == 1);
[~,h_left] = fdr(p_left,0.05);
idx_left = find(h_left == 1);

ds.chanPairs = pairs(idx_right,:); 
ds.connectStrength = h_right(idx_right,1);
ds.connectStrengthLimits = [0 1];
figure;title(strcat('significant pairs: ',name1,'>',name2));topoplot_connect(ds, chanlocs);

ds.chanPairs = pairs(idx_left,:); 
ds.connectStrength = h_left(idx_left,1);
ds.connectStrengthLimits = [0 1];
figure;title(strcat('significant pairs: ',name1,'<',name2));topoplot_connect(ds, chanlocs);



%%
%% 

%ѡ�缫�Ժ��ҵ�����ʱƵ�������ȡ����

t_roi = inputdlg('The temporal limits of TF-ROI (in ms)');
t_roi = str2num(t_roi{1});

f_roi = inputdlg('The spectral limits of TF-ROI (in Hz)');
f_roi = str2num(f_roi{1});


A_fc_roxi = squeeze(mean(mean(A_fc_x(f_axis >= f_roi(1,1) & f_axis <= f_roi(1,2),time_axis >= t_roi(1,1) & time_axis <= t_roi(1,2),:,:,:),1),2));
B_fc_roxi = squeeze(mean(mean(B_fc_x(f_axis >= f_roi(1,1) & f_axis <= f_roi(1,2),time_axis >= t_roi(1,1) & time_axis <= t_roi(1,2),:,:,:),1),2));


%%ѡʱƵ�����ҵ������缫����ȡ����


channel_pair = inputdlg('The channel pair to be test (e.g., [6 5])');
channel_pair = str2num(channel_pair{1});


A_fc_xi = squeeze(A_fc_roi(channel_pair(1,1),channel_pair(1,2),:));
B_fc_xi = squeeze(B_fc_roi(channel_pair(1,1),channel_pair(1,2),:));


%% ============================================================
%% ׷��ģ�飺�����������ӵ� BrainNet Viewer (3D ��ͼ)
%% ============================================================
% �߼�������������� idx_right (A>B) �� idx_left (A<B) ����Ϊ .edge �ļ�
% �������Ϳ����� BrainNet Viewer �ﻭ���Ǹ���SR��񡱵ĸ߼�ͼ�ˡ�

choice_bnv = questdlg('ͳ����ɡ��Ƿ񵼳������������� 3D �������ͼ (BrainNet Viewer)?', ...
    '3D ����', 'Yes', 'No', 'Yes');

if strcmp(choice_bnv, 'Yes')
    % 1. ׼������·��
    save_dir = uigetdir(pwd, 'ѡ�񱣴� .node �� .edge �ļ����ļ���');
    if isequal(save_dir, 0), return; end
    
    % ---------------------------------------------------------
    % A. ���� .node �ļ� (�缫����)
    % ---------------------------------------------------------
    % BrainNet Viewer ��Ҫ MNI ���ꡣ������� .set �Ǳ�׼�缫�����ǽ��н���ת����
    nChans = length(chanlocs);
    node_data = zeros(nChans, 6); % [X, Y, Z, Color, Size, Label]
    scale_factor = 85; % �Ŵ�ϵ����ȷ���缫����Ƥ�����
    
    for i = 1:nChans
        % ����ת�� (������ -> �ѿ��� -> �Ŵ�)
        if isfield(chanlocs, 'X') && ~isempty(chanlocs(i).X)
            x=chanlocs(i).X; y=chanlocs(i).Y; z=chanlocs(i).Z;
        elseif isfield(chanlocs, 'theta')
            theta = chanlocs(i).theta; radius = chanlocs(i).radius;
            % ��ͶӰ����
            phi = radius * pi/2; 
            x = sin(phi) * cos(theta * pi/180);
            y = sin(phi) * sin(theta * pi/180);
            z = cos(phi);
        else
            x=0; y=0; z=0;
        end
        
        vec = [x, y, z];
        if norm(vec)>0, vec = vec/norm(vec)*scale_factor; end
        
        node_data(i, 1:3) = vec;
        node_data(i, 4) = 1; % Ĭ����ɫ
        node_data(i, 5) = 2; % Ĭ�ϴ�С
    end
    
    % ���� .node �ļ�
    node_file = fullfile(save_dir, 'Electrodes.node');
    fid = fopen(node_file, 'w');
    for i = 1:nChans
        fprintf(fid, '%.4f %.4f %.4f %d %d %s\n', ...
            node_data(i,1), node_data(i,2), node_data(i,3), ...
            node_data(i,4), node_data(i,5), chanlocs(i).labels);
    end
    fclose(fid);
    
    % ---------------------------------------------------------
    % B. ���� .edge �ļ� (�������Ӿ���)
    % ---------------------------------------------------------
    % �������������ļ���
    % 1. Pos_Diff.edge (��Ӧ A > B, �� idx_right)
    % 2. Neg_Diff.edge (��Ӧ A < B, �� idx_left)
    
    % --- ���� A > B (Right Tail) ---
    matrix_right = zeros(nChans, nChans);
    if ~isempty(idx_right)
        % idx_right ������pairs���к�
        sig_pairs_right = pairs(idx_right, :); 
        % ��Ӧ�� tֵ �� pֵ (������ 1 ��ʾ����������������� connectStrength)
        % �������ʾ����ǿ�ȣ������� mean difference
        for k = 1:length(idx_right)
            p1 = sig_pairs_right(k, 1);
            p2 = sig_pairs_right(k, 2);
            % BrainNet Viewer �����ǶԳƵ�
            matrix_right(p1, p2) = 1; 
            matrix_right(p2, p1) = 1; 
        end
    end
    dlmwrite(fullfile(save_dir, [name1 '_GT_' name2 '.edge']), matrix_right, 'delimiter', '\t');
    
    % --- ���� A < B (Left Tail) ---
    matrix_left = zeros(nChans, nChans);
    if ~isempty(idx_left)
        sig_pairs_left = pairs(idx_left, :); 
        for k = 1:length(idx_left)
            p1 = sig_pairs_left(k, 1);
            p2 = sig_pairs_left(k, 2);
            matrix_left(p1, p2) = 1;
            matrix_left(p2, p1) = 1;
        end
    end
    dlmwrite(fullfile(save_dir, [name1 '_LT_' name2 '.edge']), matrix_left, 'delimiter', '\t');
    
    msgbox({'BrainNet Viewer �ļ�������!'; ...
            ['1. Node�ļ�: Electrodes.node']; ...
            ['2. Edge�ļ�(A>B): ' name1 '_GT_' name2 '.edge']; ...
            ['3. Edge�ļ�(A<B): ' name1 '_LT_' name2 '.edge']; ...
            ''; ...
            '��� BrainNet Viewer ���� ICBM152 ģ��������ļ����л�ͼ��'}, ...
            '���');
end
































%%  ��Ϣ̬����



%%



%%

%% ��һ��Ԥ����
%% ���ݸ�ʽ��ʽת��
clc;clear;
%% Specify Basic information of different groups     ��cdt�ļ�ת����set�ļ���������Ϊ500��ȥ�����õ缫�����档
group_dir = 'path_to_your_data';     % �˴�·����Ҫ����Ϊ�Լ����ļ�Ŀ¼
group_files = dir([group_dir, filesep, '*.dap']);  %filesep��\����˼
for i=1:length(group_files)
    subj_fn = group_files(i).name;
    EEG = loadcurry(strcat(group_dir, filesep, subj_fn), 'CurryLocations', 'False');    %����ԭʼ����
    EEG = pop_resample( EEG, 250);   %������
    EEG = pop_select( EEG, 'rmchannel',{'HEO','VEO','EKG','EMG','TRIGGER','CB1','CB2'});
    EEG = pop_reref( EEG, []); %ȫ�Բο�
    %EEG = pop_reref( EEG, [33 43] );%˫����ͻ�ο�
     %��Ϣ̬�ֶ� �ֳ�����һ��
    EEG = eeg_regepochs(EEG, 'recurrence', 2, 'limits',[0 2], 'rmbase',NaN);
    %����ֳ�m��һ�Σ�
    %EEG = eeg_regepochs(EEG, 'recurrence', m, 'limits',[0 m], 'rmbase',NaN);
    %���»���eeglab�Ӵ�
    eeglab redraw
    EEG = pop_eegfiltnew(EEG, 'locutoff',0.1,'plotfreqz',1);
    EEG = pop_eegfiltnew(EEG, 'hicutoff',40,'plotfreqz',1);
    EEG = pop_eegfiltnew(EEG, 'locutoff',48,'hicutoff',52,'revfilt',1,'plotfreqz',1); 
    EEG = pop_saveset( EEG, 'filename',strcat(group_files(i).name(1:end-4), '.set'), 'filepath',strcat(group_dir, filesep, '_step1'));   %ע����Ҫ�����д���֮ǰ���ļ�Ŀ¼�½�һ��_resam_remch���ļ��У�������ͬ
end

%%�ļ���_step1���ڴ�����Ǿ�������������ͨ�˲����ݲ��˲���ȥ���޹ص缫֮������ݣ�
%�ļ���_preica���ڴ�����Ǿ���������޻��λ���֮������ݣ�
%�ļ���_ica���ڴ����������ICA֮������ݣ�
%�ļ���_rm_ica���ڴ��ICLabel�Զ�ȥ��α���Լ�ȫ��ƽ���زο�֮������ݣ�
%ĸ�ļ���demo�´��ԭʼ�ɼ������ݡ�
%% �˵���������eeglab����_step1���ļ��������޻������Σ����л������ֵ�����л�����ɾ�������֮�󱣴����ݵ�_preica�ļ��У��ļ����Ʊ��ֲ��䡣

%%�ļ����룺��eeglab�����File->load existing dataset��
%�ļ����棺���File->save current dataset as��

%%   ���� ICA
clc;clear;
group1_dir = 'path_to_your_data';     % �˴�·����Ҫ����Ϊ�Լ����ļ�Ŀ¼
group1_dir1 = 'path_to_your_data';     % �˴�·����Ҫ����Ϊ�Լ����ļ�Ŀ¼
group1_files = dir([group1_dir1, filesep, '*.set']);  %filesep��\����˼
for i=1:length(group1_files)
    subj_fn = group1_files(i).name;
    EEG = pop_loadset('filename',strcat(subj_fn(1:end-4), '.set'), 'filepath', strcat(group1_dir, filesep, '_preica')); %��������
    EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'interrupt','on');   % ��ICA
    EEG = pop_saveset( EEG, 'filename',strcat(group1_files(i).name(1:end-4), '.set'), 'filepath',strcat(group1_dir, filesep, '_ica'));  %��������
    
end
%% ʹ��ICLabel�Զ�ȥ��ICA�ɷ�
clc;clear;
group1_dir ='path_to_your_data';        % �˴�·����Ҫ����Ϊ�Լ����ļ�Ŀ¼
group1_dir2 = 'path_to_your_data'; 
group1_files = dir([group1_dir2, filesep, '*.set']);  %filesep��\����˼
for i=1:length(group1_files)
    subj_fn = group1_files(i).name;
    EEG = pop_loadset('filename',strcat(subj_fn(1:end-4), '.set'), 'filepath', group1_dir2);
    EEG = pop_iclabel(EEG, 'default');
    EEG = pop_icflag(EEG, [NaN NaN;0.5 1;0.5 1;NaN NaN;NaN NaN;NaN NaN;NaN NaN]); % ���α���ɷ֡���������Զ����趨��ֵ������ΪBrain, Muscle, Eye, Heart, Line Noise, Channel Noise, Other.
    EEG = pop_subcomp( EEG, [], 0)   %ȥ������α���ɷ�
    %EEG = pop_reref( EEG, []);    %ȫ��ƽ���زο�
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',strcat(group1_files(i).name(1:end-4), '.set'), 'filepath',strcat(group1_dir, filesep, '_rm_ica')); 
end

%% �ֶ�������� ȷ��Ԥ����֮��������Ǹɾ���


%%
%%
%%
%% ��2���֣����ݷ�Ϊ�����ļ��������ţ����й�������
%%
clear;clc
%�����������ڵ�·��
data_path = 'path_to_your_data';
%�������ݵı���·��
save_path = 'path_to_your_data';
%���������ڵ�·������Ϊ����·��
cd(data_path)
%ɸѡ��ǰ·�������е�vhdr��β���ļ�
files = dir('*.set');
%��ȡ�ļ���
fn = {files.name};
%����ÿ������
for i = 1:length(fn)
    %����set��ʽ������
    EEG = pop_loadset('filename',fn{i},'filepath',data_path);
    EEG = eeg_checkset( EEG ); 
    EEG = pop_select( EEG,'nochannel',{'AF3' 'AF4' 'F7' 'F8' 'FT7' 'FT8' 'T7' 'T8' 'TP7' 'TP8' 'P7' 'P8' 'PO7' 'PO8'});
    EEG = eeg_checkset( EEG );
    EEG = pop_resample( EEG, 250); 
    %EEG = eeg_regepochs(EEG, 'recurrence', 2, 'limits',[0 2], 'rmbase',NaN);
    %EEG = pop_rmbase( EEG, [-1000 0] ,[]);
    %����set��ʽ������
    %EEG.subject = fn{1,i}(1:4);
    %EEG.group = '0';
    %EEG.condition = '22';
    %original_subject_name = fn{i}(1:4); 
    % ����EEG.subjectΪ�µ�����
    %EEG.subject = [original_subject_name '_' EEG.group '_' EEG.condition];
    EEG = pop_saveset( EEG, 'filename',fn{i},'filepath',save_path);
end
%%
%% ������ csdת���������ɲ�����
%compute_CSDdata.m����ÿ�����ݵ�CSD����ѡ����
%���д˽ű�ǰ����Ҫ��CSDtoolbox���ؽ�MATLAB·�������д˽ű��󣬻ᵯ�������Ի�����ֱ��������Ի���ָ������ת����EEG�������ڵ��ļ��С��͡�ת�������ɵ�CSD�������ڵ�·����



%%
%% ���Ĳ������ݽ��м���

%%compute_rest_fc.m
% �˽ű������ڼ��㾲Ϣ̬EEG�Ĺ�����ָͨ�꣬����coherence��PLV��PLI����
%


clear all; close all; clc

%% ָ�������Ϣ
Data_Dir = uigetdir([],'Path of the EEG datasets'); %ָ�����������Ե��������ڵ�·����������Ԥ�����ķֶ��Ե����ݡ�Ҳ������CSDת����ķֶ����ݣ�����һ��������һ�����һ���ļ��У�ÿ���ļ��зֱ����иô���
Output_Dir = uigetdir([],'Path to store the measures');%����õ���ָ���ļ������·��

band = inputdlg('the limits of band');%ָ����Ҫ������Ƶ�ʵķ�Χ����λ��Hz��
band = str2num(band{1}); %��band�������ַ�ת��Ϊ��ֵ��ע�⣺���б��ű�ʱMATLAB��·���в�Ӧ����HERMES������

bandname = inputdlg('the name of the band you computed'); % ָ��ǰ����������Ƶ�ε�����
bandname = bandname{1};

prefix = inputdlg('the prefix of computed measures');%ָ�������ָ���ļ���ǰ׺���������������������
prefix = prefix{1};

%% ��ȡEEG����·���а�����set�ļ�
Dir_Data = dir(fullfile(Data_Dir,'*.set')); 
FileNames = {Dir_Data.name};

%% compute 1st measure (coherence)
for subj = 1:numel(FileNames)  
    EEG = pop_loadset('filename',FileNames{1,subj},'filepath',Data_Dir); % ����ĳ�����Ե�����
    EEG = eeg_checkset( EEG );
    N = EEG.pnts; %ÿ�εĳ��ȣ�������
    SampleRate = EEG.srate;%ȡ����
    NFFT = 2^nextpow2(N);%����ÿ�γ��ȵġ���С��2��N�η�
    Freq = SampleRate/2*linspace(0,1,NFFT/2+1);%Ƶ���� 
    for chan = 1:size(EEG.data,1)
        for epochs = 1:size(EEG.data,3)
            ffts(:,chan,epochs) = fft(hanning(N).*squeeze(EEG.data(chan,:,epochs))',NFFT);% �Ըñ���ÿ���缫��ÿ���ֶν���FFT            
        end
    end
    for x = 1:size(EEG.data,1)
            for y = 1:size(EEG.data,1)
                fx = squeeze(ffts(:,x,:));
                Pxx = fx.*conj(fx)/N;
                MeanPx = mean(Pxx,2); % ����coherenceʱ��x�缫��power
                fy = squeeze(ffts(:,y,:));
                Pyy = fy.*conj(fy)/N; 
                MeanPy = mean(Pyy,2); % ����coherenceʱ��y�缫��power
                Pxy = fx.*conj(fy)/N;
                MeanPxy = mean(Pxy,2); %% Sxy�����������缫�Ľ�����
                C = (abs(MeanPxy).^2)./(MeanPx.*MeanPy); % ���
                coh(:,x,y,subj) = C; % coherence��Ƶ��*�缫*�缫*����
            end
    end
    clear ffts
end
coh = coh(1:NFFT/2 + 1,:,:,:);  

idx = dsearchn(Freq', band'); %ȷ��Ƶ����������Freq��Ƶ�����е�λ��
coh = squeeze(mean(coh(idx(1,1):idx(2,1),:,:,:),1));%����ĳ��Ƶ�ε�ƽ��coh

save(strcat(Output_Dir,'\',prefix,'_',bandname,'_coh.mat'),'coh', '-v7.3');

%% 2nd and 3rd measures (phase-locking value and phase lag index)   
for subj = 1:numel(FileNames)
    EEG = pop_loadset('filename',FileNames{1,subj},'filepath',Data_Dir);
    EEG = eeg_checkset( EEG );
    eeg_filtered = eegfilt(reshape(EEG.data, [size(EEG.data,1) size(EEG.data,2)*size(EEG.data,3)]),...
                   EEG.srate,band(1,1),band(1,2),0,3*fix(EEG.srate/band(1,1)),0,'fir1',0); % ����������ݽ��д�ͨ�˲���ʹ��eegfilt���д�ͨ�˲�ʱ��Ҫ���ֶ��������±�Ϊ�������� 
         
    for channels = 1:size(EEG.data,1)
        band_phase(channels,:) = angle(hilbert(eeg_filtered(channels,:))); %����缫����Hilbert�任������ȡ��λ
    end    
    perc10w =  floor(size(band_phase,2)*0.1);% ȷ�����ݳ���10%�Ƕ��ٸ�������
    band_phase = band_phase(:,perc10w+1:end-perc10w); %��Hilbert�任��������β��λ���㲻׼ȷ����ȥ��ǰ10%�ͺ�10%���������λ
    epoch_num = floor(size(band_phase,2)/size(EEG.data,2)); % ȷ��ʣ������������ת��Ϊ�ֶ����ݣ����Էֳɶ��ٶ�
    band_phase = band_phase(:,1:epoch_num*size(EEG.data,2)); % ���ݿ��ԷֳɵĶ�������ȡ����
    band_phase = reshape(band_phase,[size(EEG.data,1) size(EEG.data,2) epoch_num]);% ����������ת��Ϊ��ά���缫*������*�ֶ�
    
    for x = 1:size(band_phase,1)
         for y = 1:size(band_phase,1)
             for epochs = 1:size(band_phase,3)
                 x_phase = squeeze(band_phase(x,:,epochs)); % ��ȡ�缫���е�һ���缫��ĳ���ֶε���λ
                 y_phase = squeeze(band_phase(y,:,epochs)); % ��ȡ�缫���еڶ����缫��ĳ���ֶε���λ
                 rp = x_phase - y_phase; % ���������缫��ĳ���ֶε���λ��
                 %%% PLV
                 sub_plv(x,y,epochs) = abs(sum(exp(1i*rp))/length(rp)); % ����ĳ������ĳ���缫����ĳ���ֶε�PLV
                 %%% PLI
                 sub_pli(x,y,epochs) = abs(mean(sign((abs(rp)- pi).*rp))); % ����ĳ������ĳ���缫����ĳ���ֶε�PLI
             end
         end
    end
   pli(:,:,subj) = mean(sub_pli,3); % �Ըñ��Ը����ֶε�PLI����ƽ��ֵ��pli����ά���ǵ缫*�缫*����
   plv(:,:,subj) = mean(sub_plv,3); % �Ըñ��Ը����ֶε�PLV����ƽ��ֵ��plv����ά���ǵ缫*�缫*����
   clear band_phase sub_pli sub_plv
end



save(strcat(Output_Dir,'\',prefix,'_',bandname,'_plv.mat'),'plv', '-v7.3');
save(strcat(Output_Dir,'\',prefix,'_',bandname,'_pli.mat'),'pli', '-v7.3');

%% compute 4th measure (weighted phase lag index)
for subj = 1:numel(FileNames)
    EEG = pop_loadset('filename',FileNames{1,subj},'filepath',Data_Dir);
    EEG = eeg_checkset( EEG );
    eeg_filtered = eegfilt(reshape(EEG.data, [size(EEG.data,1) size(EEG.data,2)*size(EEG.data,3)]),...
                   EEG.srate,band(1,1),band(1,2),0,3*fix(EEG.srate/band(1,1)),0,'fir1',0); 
    
    for channels = 1:size(EEG.data,1)
        band_hilbert(channels,:) = hilbert(eeg_filtered(channels,:)); % �Ըñ���ÿ��ͨ�����źŽ���hilbert�任���õ������źţ�a + b*i��������wPLI����Ҫ��ȡ��λ
    end 
    perc10w =  floor(size(band_hilbert,2)*0.1);
    band_hilbert = band_hilbert(:,perc10w+1:end-perc10w);
    epoch_num = floor(size(band_hilbert,2)/size(EEG.data,2));
    band_hilbert = band_hilbert(:,1:epoch_num*size(EEG.data,2));
    band_hilbert = reshape(band_hilbert,[size(EEG.data,1) size(EEG.data,2) epoch_num]); 
    
    for x = 1:size(band_hilbert,1)
         for y = 1:size(band_hilbert,1)
             for epochs = 1:size(band_hilbert,3)
                 x_hilbert = band_hilbert(x,:,epochs);
                 y_hilbert = band_hilbert(y,:,epochs);
                 crossspec = x_hilbert.* conj(y_hilbert); % ������
                 crossspec_imag = imag(crossspec); % �����׵��鲿
                 sub_wpli(x,y,epochs) = abs(mean(crossspec_imag))/mean(abs(crossspec_imag)); % ����ĳ������ĳ���缫����ĳ���ֶε�wPLI
             end
         end
    end
    wpli(:,:,subj) = mean(sub_wpli,3); 
    clear band_hilbert sub_wpli
end
wpli(isnan(wpli)) = 0;% ��Ϊ�Խ����ϵ�wpli��ֵΪnan���ʽ����������
save(strcat(Output_Dir,'\',prefix,'_',bandname,'_wpli.mat'),'wpli', '-v7.3');  












%%
clc;clear all;close all

%% ָ�������Ϣ
[filename1, pathname1, filterindex] = uigetfile('*.mat', 'Pick the mat file of measures of one condition (group)'); % ָ����һ������������ָ���ļ�
[filename2, pathname2, filterindex] = uigetfile('*.mat', 'Pick the mat file of measures of the other condition (group)'); % ָ���ڶ�������������ָ���ļ�

[filename3, pathname3, filterindex] = uigetfile('*.set', 'Pick an eeglab file used to computed the measures'); % ��ΪҪ���Ƶ缫��缫������ͼ�������Ҫ֪���缫������Ϣ�������������ĳ�����Ե��Ե�����

test_type = inputdlg('Test Type? Paired t-test = 1 Independent t-test = 0'); % ָ�����������
test_type = str2num(test_type{1});

name1 = inputdlg('the name of one condition (group)'); % ָ����һ������������������
name1 = name1{1};

name2 = inputdlg('the name of the other condition (group)'); % ָ���ڶ�������������������
name2 = name2{1};

%% ����ĳ��eeglab��ʽ�Ե����ݣ���ȡ��ʵ���Ե����ݵ������Ϣ
EEG = pop_loadset('filename',filename3,'filepath',pathname3);
chanlocs = EEG.chanlocs; % ��ͨ��λ����Ϣ������ chanlocs������
channel_number = size(EEG.data,1); % ��ȡ�缫��Ŀ
fc_number = channel_number*(channel_number-1)/2; % ��ȡ�缫�Ե���Ŀ

%% ȷ��ÿ���缫�Էֱ��Ӧ�������缫
pairs = nchoosek(1:channel_number,2); 

%% ������������������𣩵�ָ��
A_fc  = importdata(strcat(pathname1,filename1)); % �����һ���������ߵ�һ�鱻�ԵĹ������Ӿ���ά���ǵ缫*�缫*����
B_fc  = importdata(strcat(pathname2,filename2)); % ����ڶ����������ߵڶ��鱻�ԵĹ������Ӿ���ά���ǵ缫*�缫*����

%% �ֱ������������������𣩵���ƽ����ͨͼ
A_fc_avg = mean(A_fc,3); % �Ա���ά����ƽ�����õ���һ����������������ƽ��ˮƽ���Ӿ���ά���ǵ缫*�缫
A_fc_avg2 = zeros(fc_number,1); % �Ա���ά����ƽ�����õ��ڶ�����������������ƽ��ˮƽ���Ӿ���ά���ǵ缫*�缫
for i = 1:fc_number
    A_fc_avg2(i,1) = A_fc_avg(pairs(i,1),pairs(i,2));% ��ȡ���Ӿ��������ǵĹ�������ֵ
end  

B_fc_avg = mean(B_fc,3);
B_fc_avg2 = zeros(fc_number,1);
for i = 1:fc_number
    B_fc_avg2(i,1) = B_fc_avg(pairs(i,1),pairs(i,2));
end  

AB_fc_avg2 = [A_fc_avg2; B_fc_avg2];

ds.chanPairs = pairs;  
ds.connectStrength = A_fc_avg2; 
ds.connectStrengthLimits = [min(AB_fc_avg2) max(AB_fc_avg2)];
figure;title(name1);topoplot_connect(ds, chanlocs); 

ds.chanPairs = pairs;  
ds.connectStrength = B_fc_avg2; 
ds.connectStrengthLimits = [min(AB_fc_avg2) max(AB_fc_avg2)];
figure;title(name2);topoplot_connect(ds, chanlocs); 

%% ��������������𣩼�ͳ�Ʒ���
for  i = 1:fc_number   
    A_fc_2(i,:) = A_fc(pairs(i,1),pairs(i,2),:);
    B_fc_2(i,:) = B_fc(pairs(i,1),pairs(i,2),:);
end

if test_type == 1 % paired
    for i = 1:fc_number 
       [~,p_right(i,1)] = ttest(A_fc_2(i,:),B_fc_2(i,:),0.05,'right');
       [~,p_left(i,1)] = ttest(A_fc_2(i,:),B_fc_2(i,:),0.05,'left');
    end
else  % indpendent
    for i = 1:fc_number 
       [~,p_right(i,1)] = ttest2(A_fc_2(i,:),B_fc_2(i,:),0.05,'right');
       [~,p_left(i,1)] = ttest2(A_fc_2(i,:),B_fc_2(i,:),0.05,'left');
    end
end
[~,h_right] = fdr(p_right,0.05);
idx_right = find(h_right == 1);
[~,h_left] = fdr(p_left,0.05);
idx_left = find(h_left == 1);

ds.chanPairs = pairs(idx_right,:); 
ds.connectStrength = h_right(idx_right,1);
ds.connectStrengthLimits = [0 1];
figure;title(strcat('significant pairs: ',name1,'>',name2));
topoplot_connect(ds, chanlocs);

ds.chanPairs = pairs(idx_left,:); 
ds.connectStrength = h_left(idx_left,1);
ds.connectStrengthLimits = [0 1];
figure;title(strcat('significant pairs: ',name1,'<',name2));
topoplot_connect(ds, chanlocs);


%%
%%
%% Ѱ�������缫

pairs(idx_right,:)
pairs(idx_left,:)


% һ��Ҫ������Ǹ���ǰ
channel_pairw = inputdlg('The channel pair to be test (e.g., [6 5])');
channel_pairw = str2num(channel_pairw{1});



% ����˵������ A_fc_roi ��Ϊ A_fc���ҵڶ���Ӧ����ȡ B_fc ������
A_fc_xtj = squeeze(A_fc(channel_pairw(1,1),channel_pairw(1,2),:));
B_fc_xtj = squeeze(B_fc(channel_pairw(1,1),channel_pairw(1,2),:));
