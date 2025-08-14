function addon_config2json()
% Convert SEPIA addon configuration m-files to JSON format
% to enable reading the data during compiled runtime

files = dir(fullfile('**', 'addon_config.m'));
for configfile = fullfile({files.folder}, {files.name})
    % Load config script
    clear addons
    run(configfile{1});    % This populates 'addons' struct

    % Convert struct to JSON and save to a JSON file
    config = jsonencode(addons);
    fid = fopen(fullfile(fileparts(configfile{1}), 'addon_config.json'), 'w');
    fprintf(fid, '%s', config);
    fclose(fid);
end
