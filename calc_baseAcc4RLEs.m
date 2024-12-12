%
% calc_baseAcc4RLE.m
%
% Description: the script calculates the base acceleration time series
% using simulation data generated with multiple structural models a priori.
%
% R. A. Romano - December 2024.
%

% - Mount the s3 volume using:
% mount-s3 --cache ~/.s3-cache/ gmto.im.grim ~/mnt
% - Update s3_path variable below


%% General Settings
%%
s3_path = "/home/rromano/mnt";
dt_fname_prot = "model-%s-RLE0%d_wiM1c_wiSGMC_1.parquet";
rle_cases = (1:7);
struct_models = [...
    "20240408_1535";...
    "20241003_1800";...
    "20230817_1808";...
    "20241021_1535"];

%% Numerical calculation of the base accelerations
%%
% Preallocate result
pier_rle = cell(7,1);

tic
for k_model = 1:numel(struct_models)
    model_label_id = struct_models(k_model);
    model_folder = dir(fullfile(s3_path,...
        sprintf('*%s*largeMass',model_label_id)));
    model_path = model_folder.name;

    %  Loop over the RLE cases
    parfor i_rle = 1:numel(rle_cases)
        sssha_data = [];
        dt_file = fullfile(s3_path, model_path,...
            sprintf(dt_fname_prot,model_label_id,rle_cases(i_rle)));
        fprintf("Post-processing data from %s\n",dt_file);

        try
            parquetINFO = parquetinfo(dt_file);
            sssha_data = parquetread(dt_file,"SampleRate",1e3,...
                "SelectedVariableNames",parquetINFO.VariableNames);
        catch
            warning('Unable to run parquetread(). Try Matlab 2022b, or later.');
        end

        t = seconds(sssha_data.Time);
        dT = diff(t(1:2));
        pier_D = reshape(cell2mat(sssha_data.Pier6D),12,[])';

        % Numerical derivative to calculate the pier node acceleration
        pier_rle{i_rle}.acc = 1/dT*diff(1/dT*diff(pier_D(:,1:3)));
        % Time vector
        pier_rle{i_rle}.time = t;
        % Base displacements
        pier_rle{i_rle}.disp = pier_D(:,1:3);
    end

    pp_dt_fname = sprintf('grsim-%s-RLE_Acc',model_label_id);
    save(pp_dt_fname,'pier_rle');
    fprintf("Base accelerations for model %s saved in %s.\n",...
        model_label_id,pp_dt_fname)
end
toc

%%
