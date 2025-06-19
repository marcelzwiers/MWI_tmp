function imab_overlay(greyimage,redimage,varargin);
% it an image with multiple slices with a backgroung image in gray and the overlay in red
% the fist input is the background 3D matrix, the second input is the 3D image
% to appear on the top
% usage
% imab_overlay(greyimage,redimage); computes colorbar automatically
% imab_overlay(greyimage,redimage,[IntensityRangeGray],[IntensityRangeRed]);
% imab_overlay(greyimage(:,:,slices),redimage(:,:,slices),[IntensityRangeGray],[IntensityRangeRed]);

if nargin==2
    scaleg=[prctile(greyimage(:),2) prctile(greyimage(:),98)];
    
    scaler=[max(redimage(:))/5 max(redimage(:))/3];
else
    if nargin>=3
        if ~isempty(varargin{1})
            scaleg=varargin{1}
        else
            scaleg=[prctile(greyimage(:),2) prctile(greyimage(:),98)];
            
        end
    end
    if nargin>=4
        if ~isempty(varargin{2})
            scaler=varargin{2}
        else
            scaler=[max(redimage(:))/5 max(redimage(:))/3];
            
        end
    end
    
end
greyimage=makemosaic(greyimage);

redimage=makemosaic(redimage);

data=(greyimage-scaleg(1))/(scaleg(2)-scaleg(1));
data(data<0)=0;data(data>1)=1;
%blanking positive values
data(redimage>scaler(1))=0;
%blanking negative values
data(redimage<-scaler(1))=0;
datargb=repmat(data,[1 1 3]);

% changes the red channel
data(redimage>scaler(1))=(redimage(redimage>scaler(1))-scaler(1))/(scaler(2)-scaler(1));
data(data>1)=1;
datargb(:,:,1)=data;


% changes the blue channel
data=(greyimage-scaleg(1))/(scaleg(2)-scaleg(1));
data(data<0)=0;data(data>1)=1;
data(redimage>scaler(1))=0;
data(redimage<-scaler(1))=0;
data(redimage<-scaler(1))=(abs(redimage(redimage<-scaler(1)))-scaler(1))/(scaler(2)-scaler(1));
data(data>1)=1;
datargb(:,:,3)=data;



% datargb = flipdim(permute(datargb,[2 1 3 4]),1);
datargb = (permute(datargb,[2 1 3 4]));
image(datargb);axis equal tight;set(gca,'Xtick',[],'YTick',[])


function imall=makemosaic(im,MaxN);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function imall=makemosaic(im,MaxN);
% Make a mosaic image for display the image with "show.m"
% i.e., the 3D image [im] transforms to a mosaic 2D image [imall]
% If [im] is 4D, [im(:,:,:,1)] will be used
% NOTE : First, singleton dimensions will be removed; 64x64x1x20 -> 3D
% Input :
%   [im] : 3D or 4D image
%   [MaxN](option): The number of colons, Default is 5
% Output :
%   [imall]: mosaic 2D image
% Usages,
% imall=makemosaic(im,MaxN);
% imall=makemosaic(im);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Examples,
% let [im] is a matrix of 64x64x20
% imall=makemosaic(im,10);
% [imall] is a 2x10 image of 64x64, size(imall)= 128 x 640
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) Jaemin Shin
% jaemins@gatech.edu
% 01/30/06
% updated
% 02/28/06 : bug fixed

if exist('MaxN','var') == 0
    MaxN = 5;
end
im = squeeze(im);
dim = size(im);
if length(dim) < 2;
    error('Input is 1D or 2D signal')
elseif length(dim) ==4
    im = squeeze(im(:,:,:,1));
    disp('4D : TimePoint 1 was used')
elseif length(dim) > 4
    error('5D or Higher dimension does not support')
end
Nrow = ceil(dim(3)/MaxN);
Rcol = mod(MaxN - mod(dim(3),MaxN),MaxN);

if dim(3) <= MaxN
    imall = reshape(im,[dim(1) dim(2)*dim(3)]);
    imall = [imall,zeros(dim(1),dim(2)*Rcol)];
else
    imall = reshape(im(:,:,1:MaxN),[dim(1) dim(2)*MaxN]);
    for ii=2:Nrow-1 % bug fixed
        temp = reshape(im(:,:,(ii-1)*MaxN+1:ii*MaxN),[dim(1) dim(2)*MaxN]);
        imall = cat(1,imall,temp);
    end
    temp = reshape(im(:,:,(Nrow-1)*MaxN+1:end),[dim(1) dim(2)*(MaxN-Rcol)]);
    temp = [temp,zeros(dim(1),dim(2)*Rcol)];
    imall = cat(1,imall,temp);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = imcat(data, outdim, indim, updown)
% function out = imcat(data, outdim, indim, updown)
% concatenate data along indim dimension and add it to dimension outdim

if nargin < 4  updown = 1;    end
if nargin < 3  indim = 3;     end

if updown
    switch indim
        case 3
            temp = data(:,:,1,:);
            if size(data,3)> 1
                for i = 2:size(data,3)
                    temp = cat(outdim,temp, data(:,:,i,:));
                end
            end
            
        case 4
            temp = data(:,:,:,1);
            if size(data,4)> 1
                for i = 2:size(data,4)
                    temp = cat(outdim,temp,data(:,:,:,i));
                end
            end
            
        case 5
            temp = data(:,:,:,:,1);
            if size(data,5)> 1
                for i = 2:size(data,5)
                    temp = cat(outdim,temp,data(:,:,:,:,i));
                end
            end
            
    end
    
else
    switch indim
        case 3
            temp = data(:,:,end,:);
            if size(data,3)> 1
                for i = size(data,3)-1:-1:1
                    temp = cat(outdim,temp, data(:,:,i,:));
                end
            end
            
        case 4
            temp = data(:,:,:,end);
            if size(data,4)> 1
                for i = size(data,4)-1:-1:1
                    temp = cat(outdim,temp,data(:,:,:,i));
                end
            end
            
        case 5
            temp = data(:,:,:,:,end);
            if size(data,5)> 1
                for i = size(data,5)-1:-1:1
                    temp = cat(outdim,temp,data(:,:,:,:,i));
                end
            end
            
    end
end

axis off
set(gca,'box','off')
out = temp;