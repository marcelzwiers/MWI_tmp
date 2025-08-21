function Hdr = spm_write_vol_gz(Hdr, data, fname)
% FUNCTION Hdr = SPM_WRITE_VOL_GZ(Hdr, data, fname)
%
% A wrapper around SPM_WRITE_VOL that removes the pinfo field, writes a
% .nii or a .nii.gz file. If fname is provided, it will override the file
% name in Hdr.fname.
%__________________________________________________________________________
%   SPM_WRITE_VOL
%
%   Write an image volume to disk, setting scales and offsets as appropriate
%   FORMAT V = spm_write_vol(V,Y)
%   V (input)  - a structure containing image volume information (see spm_vol)
%   Y          - a one, two or three dimensional matrix containing the image voxels
%   V (output) - data structure after modification for writing.
%
%   Note that if there is no 'pinfo' field, then SPM will figure out the
%   max and min values from the data and use these to automatically determine
%   scalefactors.  If 'pinfo' exists, then the scalefactor in this is used.

if nargin > 2 && ~isempty(fname)
    Hdr.fname = fname;
end

[fname, ~, Ext] = myfileparts(Hdr.fname);
Hdr             = rmfield(Hdr, 'pinfo');
[~,~]           = mkdir(fname);
switch Ext
    case '.nii.gz'
        Hdr.fname = spm_file(Hdr.fname, 'ext','');
        Hdr       = spm_write_vol(Hdr, data);
        gzip(Hdr.fname)
        delete(Hdr.fname)
    case '.nii'
        Hdr = spm_write_vol(Hdr, data);
    otherwise
        error('Unknown file extenstion %s in %s', Ext, FileName)
end
