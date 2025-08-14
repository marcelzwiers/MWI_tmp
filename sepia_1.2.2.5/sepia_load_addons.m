%% sepia_load_addons
%
% Description: a script provides add-on ability in SEPIA
% names
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 27 June 2020
% Date modified: 12 June 2021
% Date modified: 3 August 2022 (v1.1)
%
% DO NOT change the variable name
% DO NOT change the order of the entities, add a new one at the end instead
%
%% find addons in these directories
addons_dir              = fullfile(SEPIA_HOME,'addons');
addons_unwrap_dir       = fullfile(addons_dir,'phase_unwrap');
addons_echo_combine_dir	= fullfile(addons_dir,'echo_combine');
addons_bfr_dir          = fullfile(addons_dir,'bfr');
addons_qsm_dir          = fullfile(addons_dir,'qsm');
addons_swismwi_dir   	= fullfile(addons_dir,'swi_smwi');

%% Phase unwrapping addons
listing = dir(addons_unwrap_dir);
addons  = struct();

for klist = 3:length(listing)
    if listing(klist).isdir 
        curr_dir = fullfile(addons_unwrap_dir,listing(klist).name);
        if exist(fullfile(curr_dir,'addon_config.m'),'file')
            addons = read_config(fullfile(curr_dir,'addon_config.m'), addons);
            methodUnwrapName{end+1}        = addons.method;
            wrapper_Unwrap_function{end+1} = addons.wrapper_function;
            gui_unwrap_exclusion{end+1}    = addons.gui_exclude_voxel;
        end
    end

end

%% Echo phase combinationaddons
listing = dir(addons_echo_combine_dir);

for klist = 3:length(listing)
    if listing(klist).isdir 
        curr_dir = fullfile(addons_echo_combine_dir,listing(klist).name);
        if exist(fullfile(curr_dir,'addon_config.m'),'file')
            addons = read_config(fullfile(curr_dir,'addon_config.m'), addons);
            methodEchoCombineName{end+1}        = addons.method;
            wrapper_EchoCombine_function{end+1} = addons.wrapper_function;
            if ~isempty(addons.gui_method_panel)
                function_EchoCombine_method_panel{end+1} = addons.gui_method_panel;
            end
            if ~isempty(addons.config_function)
                config_EchoCombine_function{end+1} = addons.config_function;
            end
        end
    end

end

%% BFR addons
listing = dir(addons_bfr_dir);

for klist = 3:length(listing)
    if listing(klist).isdir 
        curr_dir = fullfile(addons_bfr_dir,listing(klist).name);
        if exist(fullfile(curr_dir,'addon_config.m'),'file')
            addons = read_config(fullfile(curr_dir,'addon_config.m'), addons);
            methodBFRName{end+1}        = addons.method;
            wrapper_BFR_function{end+1} = addons.wrapper_function;
            if ~isempty(addons.gui_method_panel)
                function_BFR_method_panel{end+1} = addons.gui_method_panel;
            end
            if ~isempty(addons.config_function)
                config_BFR_function{end+1} = addons.config_function;
            end
        end
    end

end

%% QSM addons
listing = dir(addons_qsm_dir);

for klist = 3:length(listing)
    if listing(klist).isdir 
        curr_dir = fullfile(addons_qsm_dir,listing(klist).name);
        if exist(fullfile(curr_dir,'addon_config.m'),'file')
            addons = read_config(fullfile(curr_dir,'addon_config.m'), addons);
            methodQSMName{end+1}        = addons.method;
            wrapper_QSM_function{end+1} = addons.wrapper_function;
            if ~isempty(addons.gui_method_panel)
                function_QSM_method_panel{end+1} = addons.gui_method_panel;
            end
            if ~isempty(addons.config_function)
                config_QSM_function{end+1} = addons.config_function;
            end
        end
    end

end

%% SWI/SMWI addons
listing = dir(addons_swismwi_dir);

for klist = 3:length(listing)
    if listing(klist).isdir 
        curr_dir = fullfile(addons_swismwi_dir,listing(klist).name);
        if exist(fullfile(curr_dir,'addon_config.m'),'file')
            addons = read_config(fullfile(curr_dir,'addon_config.m'), addons);
            methodSWISMWIName{end+1}        = addons.method;
            wrapper_SWISMWI_function{end+1} = addons.wrapper_function;
            if ~isempty(addons.gui_method_panel)
                function_SWISMWI_method_panel{end+1} = addons.gui_method_panel;
            end
            if ~isempty(addons.config_function)
                config_SWISMWI_function{end+1} = addons.config_function;
            end
        end
    end

end


function addons = read_config(config_file, addons)
% Read the addon configuration file and return the addons struct

if ~isdeployed
    run(config_file)   % This will populate the 'addons' struct
else
    % For deployed applications, we need to read the file directly
    try
        [p, f]  = fileparts(config_file);
        addons_ = jsondecode(fileread(fullfile(p, f + ".json")));
        for fieldname = fieldnames(addons_)'
            addons.(fieldname{1}) = addons_.(fieldname{1});
        end
    catch exception
        error('Failed to read addon JSON-configuration: %s\nRun addon_config2json before compiling the runtime', exception.message);
    end
end

end