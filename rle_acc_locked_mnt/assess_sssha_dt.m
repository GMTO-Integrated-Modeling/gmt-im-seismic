%
% RLE SSSHA Assessment Script
%
% Sep 2024
%

%% Script flags
rle_id = 2;
reload_model = true; %false; %
plot_sssha_acc = true;
compare_sa = true;

useLogLog = false;
% Create function handle based on useLogLog
if useLogLog, compplotFunc = @loglog;
else, compplotFunc = @semilogx;
end

dt_file = sprintf('./data/model-20250306_2131-RLE0%d_.parquet',rle_id);
% dt_file = sprintf('./data/model-20250306_2132-RLE0%d_.parquet',rle_id);
dt_file_label = replace(extractBetween(dt_file,"data/","-RLE"),'_','\_');

% PSD settings (Welch algorithm)
psd.M = 1024;
psd.nFFT = 2^14;
psd.detrend = 'none';

%% Simulation data
%%
try
    parquetINFO = parquetinfo(dt_file);
    sssha_data = parquetread(dt_file,"SampleRate",1e3,...
        "SelectedVariableNames",parquetINFO.VariableNames);
    % "OSSPayloads6D";"OSS00GroundAcc";"OSSHardpointD";"OSSM1Lcl";"MountEncoders"
catch
    warning('Unable to run parquetread(). Try Matlab 2022b, or later.');
end

t = seconds(sssha_data.Time);
Ts = diff(t(1:2));

%% Acceleration and Mount ENC plots
%%
t_range = [0.005,max(t)];%[];%[4.19,4.208];%[];%
if(isempty(t_range))
    t_idx  = [1, length(t)];
else
    t_idx = [find(t >= t_range(1),1,"first"),find(t < t_range(2),1,"last")];
end

if (any(contains(parquetINFO.VariableNames,"OSS00GroundAcc")) && 1)
    gnd_acc = reshape(cell2mat(sssha_data.OSS00GroundAcc),3,[]);
    
    % GND_ACC plot
    acc_ps = zeros(psd.nFFT/2+1, 3);
    for ik = 1:3
        [acc_ps(:,ik),freqP] = utils.pwelch(gnd_acc(ik,:)'...
            ,2*psd.M,[],psd.nFFT,1/Ts,'onesided',psd.detrend);
    end
    if(plot_sssha_acc)
        figure(rle_id)
        set(gcf,'position',[123   230   400   400])
        subplot(2,1,1)
        plot(t(t_idx(1):t_idx(2)),gnd_acc(:, t_idx(1):t_idx(2))');
        set(gca,'ColorOrderIndex',1); hold on;
        plot(t([t_idx(1),t_idx(2)]), kron([1 1],mean(gnd_acc,2)), '--');
        ylabel('Accelerations (m/s^2)');
        xlabel('Time (s)');
        legend('H1','H2','V');
        title(sprintf('%s - RLE0%d',dt_file_label{1},rle_id))
        grid on; axis tight; hold off

        subplot(2,1,2)
        semilogx(freqP, acc_ps);
        xlabel('Frequency (Hz)');
        ylabel('GND Acc ((m/s^2)^2/Hz)');
        grid on; axis tight;
    end

    gnd_num_x = cumsum(Ts*cumsum(Ts*gnd_acc(:, t_idx(1):t_idx(2))'));
    if (any(contains(parquetINFO.VariableNames,"OSS00Ground6D")) && 1)        
        gnd_D = reshape(cell2mat(sssha_data.OSS00Ground6D),6,[]);
        account4lmk = false;
        if(account4lmk)
            gnd_num_x = gnd_num_x - 8.5294e9/2.1605e12*...
                cumsum(Ts*cumsum(Ts*gnd_D(1:3, t_idx(1):t_idx(2))'));
            warning("Accounting for the large mass spring on the GND position calculation!");
        end
    else
        gnd_D = [];
    end

    % GND motion verification plot
    figure(10+rle_id)
    set(gcf,'position',[423   250   740   400])
    plot(t(t_idx(1):t_idx(2)),gnd_num_x, '-');
    hold on;
    legend_str = {'\int\int H1 acc','\int\int H2 acc','\int\int V acc'};
        
    if(~isempty(gnd_D))
        set(gca,'ColorOrderIndex',1); 
        plot(t(t_idx(1):t_idx(2)), gnd_D(1:3, t_idx(1):t_idx(2))',':','Linewidth',2);
        legend_str = [legend_str(:)',{'GND x motion'},{'GND y motion'},{'GND z motion'}];
%         plot(t(t_idx(1):t_idx(2)), gnd_D(4:6, t_idx(1):t_idx(2))', '-'); % ZERO        
    end
    if(any(contains(parquetINFO.VariableNames,"Pier6D")) && 1)
        pier_D = reshape(cell2mat(sssha_data.Pier6D),12,[]);
        plot(t(t_idx(1):t_idx(2)), pier_D(1:3, t_idx(1):t_idx(2))', '-.');
        legend_str = [legend_str(:)',{'Pier(B) X^{\rightarrow}'},...
            {'Pier(B) Y^{\rightarrow}'},{'Pier(B) Z^{\rightarrow}'}];
        set(gca,'ColorOrderIndex',4);
        plot(t(t_idx(1):t_idx(2)), pier_D((1:3)+6, t_idx(1):t_idx(2))',...
            '--','Linewidth',1.5);
        legend_str = [legend_str(:)',{'Pier(T) X^{\rightarrow}'},...
            {'Pier(T) Y^{\rightarrow}'},{'Pier(T) Z^{\rightarrow}'}];
    end
    
    ylabel('Motion (m)');
    xlabel('Time (s)');
    legend(legend_str,'Location','Southeast');
    title(sprintf('%s - RLE0%d',dt_file_label{1},rle_id))
    grid on; axis tight; hold off;
end


%% Pier ACC or [Motion plots]
%%
if(any(contains(parquetINFO.VariableNames,"Pier6D")) && 1)
    pier_D = reshape(cell2mat(sssha_data.Pier6D),12,[]);
    pier_acc_dt = 1/Ts*diff(1/Ts*diff(pier_D(1:3,t_idx(1)-2:t_idx(2))'));

    if(compare_sa)                
        dt_file_ = sprintf('./data/model-20250306_2132-RLE0%d_.parquet', rle_id); %#ok<*UNRCH> 
        sssha_data = parquetread(dt_file_,"SampleRate",1e3,...
            "SelectedVariableNames","Pier6D");
        pier_D_ = reshape(cell2mat(sssha_data.Pier6D),12,[]);
        pier_acc_dt_ = 1/Ts*diff(1/Ts*diff(pier_D_(1:3,t_idx(1)-2:t_idx(2))'));
        dt_file_label_ = replace(extractBetween(dt_file_,"data/","-RLE"),'_','\_');

        pier_b_Txyz_ps_ = zeros(psd.nFFT/2+1, 3);
        for iu=1:3
            [pier_b_Txyz_ps_(:,iu),freqP] = utils.pwelch(pier_acc_dt_(:,iu),...
                2*psd.M,[],psd.nFFT,1/Ts,'onesided',psd.detrend);
        end
        % ---
        s3_path = "/home/rromano/mnt";        
%         model_label_id = '20230817_1808';
        model_label_id = '20240408_1535';
        model_folder = dir(fullfile(s3_path,sprintf('*%s*largeMass',model_label_id)));
        model_path = model_folder.name;

        dt_fname_prot = "model-%s-RLE0%d_wiM1c_wiSGMC_1.parquet";
        dt_file__ = fullfile(s3_path, model_path,...
            sprintf(dt_fname_prot,model_label_id,rle_id));
        fprintf("Post-processing data from \n%s\n",dt_file__);

        sssha_data = parquetread(dt_file__,"SampleRate",1e3,...
            "SelectedVariableNames","Pier6D");
        pier_D__ = reshape(cell2mat(sssha_data.Pier6D),12,[]);
        pier_acc_dt__ = 1/Ts*diff(1/Ts*diff(pier_D__(1:3,t_idx(1)-2:t_idx(2))'));
        dt_file_label__ = replace(model_label_id,'_','\_');

        pier_b_Txyz_ps__ = zeros(psd.nFFT/2+1, 3);
        for iu=1:3
            [pier_b_Txyz_ps__(:,iu),freqP] = utils.pwelch(pier_acc_dt__(:,iu),...
                2*psd.M,[],psd.nFFT,1/Ts,'onesided',psd.detrend);
        end
        % ---
%         model_label_id = '20240408_1535';
%         model_folder = dir(fullfile(s3_path,sprintf('*%s*largeMass',model_label_id)));
%         model_path = model_folder.name;

%         dt_file___ = fullfile(s3_path, model_path,...
%             sprintf(dt_fname_prot,model_label_id,rle_id));
        dt_file___ = sprintf('./data/model-20250313_0917-RLE0%d_.parquet', rle_id); %#ok<*UNRCH> 
        sssha_data = parquetread(dt_file___,"SampleRate",1e3,...
            "SelectedVariableNames","Pier6D");
        fprintf("Post-processing data from \n%s\n",dt_file___);

%         sssha_data = parquetread(dt_file___,"SampleRate",1e3,...
%             "SelectedVariableNames","Pier6D");
        pier_D___ = reshape(cell2mat(sssha_data.Pier6D),12,[]);
        pier_acc_dt___ = 1/Ts*diff(1/Ts*diff(pier_D___(1:3,t_idx(1)-2:t_idx(2))'));
        dt_file_label___ = replace(extractBetween(dt_file___,"data/","-RLE"),'_','\_');
%         replace(model_label_id,'_','\_');

        pier_b_Txyz_ps___ = zeros(psd.nFFT/2+1, 3);
        for iu=1:3
            [pier_b_Txyz_ps___(:,iu),freqP] = utils.pwelch(pier_acc_dt___(:,iu),...
                2*psd.M,[],psd.nFFT,1/Ts,'onesided',psd.detrend);
        end
    end

    plot_ddot = true;%false;
    %     pier_t = pier_D(7:12,t_idx(1):t_idx(2))';

    if(plot_ddot)
        fid_offset = 100;%200;
        pier_out_dt = pier_acc_dt;
        y_label1 = 'Pier Acc time response (%s/s^2)';
        y_label2 = 'Pier Acc PSD (%s^2/s^4/Hz)';
        plot_title_str = "Bottom Pier Node Acceleration";
%         pier_out_dt = 1/Ts*diff(1/Ts*diff(...
%             pier_D(1:6,t_idx(1)-2:t_idx(2))' - gnd_D(:, t_idx(1)-2:t_idx(2))'));
%         y_label1 = 'Relative Acc time response (%s/s^2)';
%         y_label2 = 'Relative Acc PSD (%s^2/s^4/Hz)';
%         plot_title_str = "Pier(B)-GND Relative Acceleration";

    else
        fid_offset = 400;
%         pier_out_dt = pier_b_x;
%         y_label1 = 'Pier motion time response (%s)';
%         y_label2 = 'Pier motion PSD (%s^2/Hz)';

%         y_label1 = '(T-B) motion time response (%s)';
%         y_label2 = '(T-B) motion PSD (%s^2/Hz)';
%         plot_title_str = "Relative pier motion (Top-bottom)";
%         pier_out_dt = pier_D(7:12,t_idx(1):t_idx(2))' - pier_D(1:6,t_idx(1):t_idx(2))';

        y_label1 = '(B-GND) motion time response (%s)';
        y_label2 = '(B-GND) motion PSD (%s^2/Hz)';
        plot_title_str = "Relative bottom pier motion (Bottom Pier - GND)";
        pier_out_dt = pier_D(1:6,t_idx(1):t_idx(2))';% - gnd_D(:, t_idx(1):t_idx(2))';
    end

    pier_b_Txyz_ps = zeros(psd.nFFT/2+1, 3);
%     pier_b_Rxyz_ps = zeros(psd.nFFT/2+1, 3);
    for iu=1:3
        [pier_b_Txyz_ps(:,iu),freqP] = utils.pwelch(pier_out_dt(:,iu),...
            2*psd.M,[],psd.nFFT,1/Ts,'onesided',psd.detrend);
%         [pier_b_Rxyz_ps(:,iu),~] = utils.pwelch(pier_out_dt(:,iu+3),...
%             2*psd.M,[],psd.nFFT,1/Ts,'onesided',psd.detrend);
    end

    legend_str = {'X^{\rightarrow}','Y^{\rightarrow}','Z^{\rightarrow}'};
    % TRANSLATIONS
    figure(1000-fid_offset+rle_id)
    set(gcf,'position',[123   80   740   400])
    subplot(2,1,1)
    plot(t(t_idx(1):t_idx(2)),pier_out_dt(:,1:3));
    ylabel(sprintf(y_label1,'m'));
    xlabel('Time (s)'); grid on; axis tight;
    legend(legend_str);    
    title(plot_title_str)
    subplot(2,1,2)
    semilogx(freqP, pier_b_Txyz_ps);
    xlabel('Frequency (Hz)'); ylabel(sprintf(y_label2,'m')); grid on; axis tight;

% ROTATIONS    
%     figure(000)
%     set(gcf,'position',[423   150   740   400])
%     subplot(2,1,1)
%     plot(t(t_idx(1):t_idx(2)),pier_out_dt(:,4:6));
%     ylabel(sprintf(y_label1,'rad'));
%     xlabel('Time (s)'); grid on; axis tight;
%     legend(legend_str);
%     title(plot_title_str)
%     subplot(2,1,2)
%     semilogx(freqP, pier_b_Rxyz_ps);
%     xlabel('Frequency (Hz)'); ylabel(sprintf(y_label2,'rad')); grid on; axis tight;
end

%% RLE Spectral Acceleration
%%

% [SRS_STRUCT] = compute_response_spectra (EQ_STRUCT, ZETA)
% B. Smith, 18 Sept 2024
zeta = .02; % payload damping ratio
STEPS = 4; % number of frequency steps per FWHM, 4 yields <3.6% scalloping
FMIN = 1; %  minimum frequency [Hz]

fs = 1/Ts; % sampling frequency
dT = 1/fs;
fmax = fs/2;
q = 1/(2*zeta); % resonance Q factor, Q = f/FWHM

nFreqs = round(log(fmax/FMIN)/log(1+1/(q*STEPS)));
fSRS = logspace(log10(FMIN),log10(fmax),nFreqs); % frequency vector [Hz]

% preallocate result
sa_data = zeros(nFreqs, min(size(pier_acc_dt))+1);
if(compare_sa)
    sa_data_ = zeros(nFreqs, min(size(pier_acc_dt_))+1);
    sa_data__ = zeros(nFreqs, min(size(pier_acc_dt__))+1);
    sa_data___ = zeros(nFreqs, min(size(pier_acc_dt___))+1);
    gnd_sa_data = zeros(nFreqs, min(size(pier_acc_dt))+1);
end
sa_data(:,1) = fSRS;
for j = 2:size(sa_data,2)
    sa_data(:,j) = SpectralA04(pier_acc_dt(:,j-1), fSRS, dT, zeta);
    if(compare_sa)
        sa_data_(:,j) = SpectralA04(pier_acc_dt_(:,j-1), fSRS, dT, zeta);
        sa_data__(:,j) = SpectralA04(pier_acc_dt__(:,j-1), fSRS, dT, zeta);
        sa_data___(:,j) = SpectralA04(pier_acc_dt___(:,j-1), fSRS, dT, zeta);
        gnd_sa_data(:,j) = SpectralA04(gnd_acc(j-1,:)', fSRS, dT, zeta);
    end
end

%%
figure(2000-fid_offset+rle_id)
set(gcf,'position',[423   150   740   400])
ylabel_str = ["H1","H2","V"];
for ik = 1:3    
    subplot(3,1,ik)
    plot(sa_data(:,1),sa_data(:,ik+1),'--','Linewidth',1.5)
    hold on;
    if(compare_sa)
        plot(sa_data(:,1),sa_data_(:,ik+1),'-');
        plot(sa_data(:,1),sa_data__(:,ik+1),'-');
        plot(sa_data(:,1),sa_data___(:,ik+1),'-');
        plot(sa_data(:,1),gnd_sa_data(:,ik+1),'-.',...
            'LineWidth',1.5,'Color',[.3 .3 .3]);
    end
    xlim([0, 100])
    grid on; ylabel(ylabel_str{ik}+" (m/s^2)"); hold off;
end
xlabel("Frequency (Hz)")
subplot(3,1,1);
title(sprintf("RLE%d - SA response (payload damping ratio:%g)",...
    rle_id, zeta));
if(compare_sa)
    legend(dt_file_label{1}, dt_file_label_{1},...
        dt_file_label__, dt_file_label___{1}, 'GND Acc');
end

%%
figure(3000-fid_offset+rle_id)
set(gcf,'position',[423   150   740   400])
ylabel_str = ["H1 PSD","H2 PSD","V PSD"];
for ik = 1:3    
    subplot(3,1,ik)
    compplotFunc(freqP, pier_b_Txyz_ps(:,ik),'--','Linewidth',1.5);
    hold on;
    if(compare_sa)
        compplotFunc(freqP, pier_b_Txyz_ps_(:,ik),'-','Linewidth',1);
        compplotFunc(freqP, pier_b_Txyz_ps__(:,ik),'-','Linewidth',1);
        compplotFunc(freqP, pier_b_Txyz_ps___(:,ik),'-','Linewidth',1);
        compplotFunc(freqP, acc_ps(:,ik),'-.',...
            'LineWidth',1.5,'Color',[.3 .3 .3]);
    end
    
    xlim([0.5, 30])
    grid on; ylabel(ylabel_str{ik}+" (m^2/s^4/Hz)"); hold off;
end
xlabel("Frequency (Hz)")
subplot(3,1,1);
title(sprintf("RLE%d - Bottom pier Acc PSD", rle_id));
if(compare_sa)
    legend([dt_file_label{1},'- K_{lat}=3e9N/m - LR'],...
        [dt_file_label_{1},'- K_{lat}=0.96e9N/m - LR'],...
        [dt_file_label__,'- K_{lat}=0.96e9N/m - CL'],...
        [dt_file_label___{1},'- K_{lat}=0.96e9N/m - CL'], 'GND Acc',...
        'Location','southwest');
end

%% Spectral Acceleration response
%%
function [Sa] = SpectralA04(InputSignal,f,dt,d)
%   SpectralA04(InputSignal,f,dt,d)
%   delivers a vector of the spectral acceleration response of the signals
%   represented in the InputSignal variable.
%   f, corresponds to the frequency vector at which the Spectral [Hz]
%   acceleration is evaluated.
%   dt corresponds to the time step of the input signal [s]
%   d corresponds to the damping ratio at which the spectral response need
%   to be evaluated
%
% AO4 - updated to max abs
    ns=min(size(InputSignal)); % number of signals
    Sa=zeros(length(f),ns); % preallocate result variable
    for i1=1:length(f)

      % continuous transfer function for harmonic oscillator
      filt1=tf([f(i1)*4*pi*d (f(i1)*2*pi)^2],[1 d*4*f(i1)*pi (f(i1)*2*pi)^2]);

      % convert continuous TF to z domain
      % FOH seems to work best per from SpectralA04_check.m
%     filt1d=c2d(filt1,dt,'zoh');
%     filt1d=c2d(filt1,dt,'prewarp',f(i1)*2*pi);
      filt1d=c2d(filt1,dt,'foh');
%     temp=get(filt1d);  %lists structure fieldnames
%      b=filt1d.num{1}; % extract coeffs (Octave compatible)
%      a=filt1d.den{1};

      [b, a] = tfdata(filt1d,'v'); % extract coeffs (Octave/Matlab compatible)

      for i2=1:ns  % filter signals
        Sa(i1,i2)=max(abs(filter(b,a,InputSignal(:,i2))));
      end

    end
end










