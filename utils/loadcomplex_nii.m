
function Compleximage=loadcomplex_nii(filemagnitude,filephase,varargin);
%Magnitude=loadavw_nii(filemagnitude);
%FieldMap=loadavw_nii(filephase);
if nargin==2
Magnitude=load_untouch_nii(filemagnitude);
FieldMap=load_untouch_nii(filephase);
else
Magnitude=load_untouch_nii(filemagnitude,varargin{1});
FieldMap=load_untouch_nii(filephase,varargin{1});
end;    
Compleximage=FieldMap;
if max(max(max(max((FieldMap.img)))))>4090
    if max(max(max(max(-(FieldMap.img)))))>4090
        Compleximage.img=single(Magnitude.img).*exp(i*single(FieldMap.img)/4096*pi);
    else
        Compleximage.img=single(Magnitude.img).*exp(i*single(FieldMap.img)/4096*2*pi);
    end;
else
    Compleximage.img=single(Magnitude.img).*exp(i*single(FieldMap.img));
end
Compleximage.hdr.dime.bitpix=8;
Compleximage.hdr.dime.datatype=8;