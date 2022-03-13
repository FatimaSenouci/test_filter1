function varargout = cancer(varargin)
% CANCER MATLAB code for cancer.fig
%      CANCER, by itself, creates a new CANCER or raises the existing
%      singleton*.
%
%      H = CANCER returns the handle to a new CANCER or the handle to
%      the existing singleton*.
%
%      CANCER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CANCER.M with the given input arguments.
%
%      CANCER('Property','Value',...) creates a new CANCER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cancer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cancer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cancer

% Last Modified by GUIDE v2.5 19-Sep-2021 17:19:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cancer_OpeningFcn, ...
                   'gui_OutputFcn',  @cancer_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
function J=regiongrowing(I,x,y,reg_maxdist)
% This function performs "region growing" in an image from a specified
% seedpoint (x,y)
%
% J = regiongrowing(I,x,y,t) 
% 
% I : input image 
% J : logical output image of region
% x,y : the position of the seedpoint (if not given uses function getpts)
% t : maximum intensity distance (defaults to 0.2)
%
% The region is iteratively grown by comparing all unallocated neighbouring pixels to the region. 
% The difference between a pixel's intensity value and the region's mean, 
% is used as a measure of similarity. The pixel with the smallest difference 
% measured this way is allocated to the respective region. 
% This process stops when the intensity difference between region mean and
% new pixel become larger than a certain treshold (t)
%
% Example:
%
% I = im2double(imread('medtest.png'));
% x=198; y=359;
% J = regiongrowing(I,x,y,0.2); 
% figure, imshow(I+J);
%
% Author: D. Kroon, University of Twente

if(exist('reg_maxdist','var')==0), reg_maxdist=0.2; end
if(exist('y','var')==0), figure, imshow(I,[]); [y,x]=getpts; y=round(y(1)); x=round(x(1)); end

J = zeros(size(I)); % Output 
Isizes = size(I); % Dimensions of input image

reg_mean = I(x,y); % The mean of the segmented region
reg_size = 1; % Number of pixels in region

% Free memory to store neighbours of the (segmented) region
neg_free = 10000; neg_pos=0;
neg_list = zeros(neg_free,3); 

pixdist=0; % Distance of the region newest pixel to the regio mean

% Neighbor locations (footprint)
neigb=[-1 0; 1 0; 0 -1;0 1];

% Start regiogrowing until distance between regio and posible new pixels become
% higher than a certain treshold
while(pixdist<reg_maxdist&&reg_size<numel(I))

    % Add new neighbors pixels
    for j=1:4,
        % Calculate the neighbour coordinate
        xn = x +neigb(j,1); yn = y +neigb(j,2);
        
        % Check if neighbour is inside or outside the image
        ins=(xn>=1)&&(yn>=1)&&(xn<=Isizes(1))&&(yn<=Isizes(2));
        
        % Add neighbor if inside and not already part of the segmented area
        if(ins&&(J(xn,yn)==0)) 
                neg_pos = neg_pos+1;
                neg_list(neg_pos,:) = [xn yn I(xn,yn)]; J(xn,yn)=1;
        end
    end

    % Add a new block of free memory
    if(neg_pos+10>neg_free), neg_free=neg_free+10000; neg_list((neg_pos+1):neg_free,:)=0; end
    
    % Add pixel with intensity nearest to the mean of the region, to the region
    dist = abs(neg_list(1:neg_pos,3)-reg_mean);
    [pixdist, index] = min(dist);
    J(x,y)=2; reg_size=reg_size+1;
    
    % Calculate the new mean of the region
    reg_mean= (reg_mean*reg_size + neg_list(index,3))/(reg_size+1);
    
    % Save the x and y coordinates of the pixel (for the neighbour add proccess)
    x = neg_list(index,1); y = neg_list(index,2);
    
    % Remove the pixel from the neighbour (check) list
    neg_list(index,:)=neg_list(neg_pos,:); neg_pos=neg_pos-1;
end

% Return the segmented area as logical matrix
J=J>1;

% End initialization code
% --- Executes just before untitled is made visible.

% --- Executes just before cancer is made visible.

function cancer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cancer (see VARARGIN)

% Choose default command line output for cancer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes cancer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = cancer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function [black white]=countBW(i)
  black=length(i(i==0));
  white=length(i(i==1));

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global img1 img2
[path, nofile] = imgetfile();
if nofile()
    msgbox(sprintf('Image not found !!!'),'Error','warning');
    return
end
img1 = imread(path);
img1 = im2double(img1);
img2 = img1;
axes(handles.axes1);
imshow(img1);

title('\fontsize{20}|-\color[rgb]{1,0,1}Tumor')

% --- Executes on button press in pushbutton2.


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global img1  J I

axes(handles.axes2)
axes(handles.axes3)


%clear all ; clc ; close all ; 
I = im2double(img1);
%I = im2double(imread('C:/medtest.png'));
%figure ; imshow(I );
J = regiongrowing(I);
%figure ; imshow(I+J);
axes(handles.axes2);
imshow(I); title('\fontsize{15}|-\color[rgb]{1,1,1} andometrial ');
axes(handles.axes3)
imshow(I+J); title('\fontsize{15}|-\color[rgb]{1,1,1} Tumor Alone ');





%***************************




% --- Executes on button press in pushbutton3.

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global J tumor_label tumor
%axes(handles.axes4)
axes(handles.axes4)
%axes(handles.axes6)



bw=im2bw(J,0.7);
label=bwlabel(J);
stats=regionprops(label,'Solidity','Area'); %the function regionprops to mesur the solidity and the area where the solidity of the tumor is mor then then andometrial  
density=[stats.Solidity];
area=[stats.Area];
high_dense_area=density>0.5; % we can use 0.2 to detect tumor in early stage  
max_area=max(area(high_dense_area));
tumor_label=find(area==max_area);
tumor=ismember(label,tumor_label);
%*************median filter *********************
if size(J,3)==3
   J=rgb2gray(J);
end
tumor = medfilt2(tumor);
%******************structur element*********************** 
se=strel('square',5);
tumor=imdilate(tumor,se); % delatation of image 
% tumor = imerode(tumor,se); % erosion gives less quality then 

%tumor = imgradient(tumor,'prewitt');% use the gradient "perwit" 
[B,L]=bwboundaries(tumor,'noholes');
 
%axes(handles.axes4);
%imshow(J); title('\fontsize{15}|-\color[rgb]{1,1,1} andometrial ');
axes(handles.axes4)
imshow(tumor); title('\fontsize{15}|-\color[rgb]{1,1,1} Tumor Alone ');
%axes(handles.axes6);
%imshow(J); title('\fontsize{15}|-\color[rgb]{1,1,1}Detected Tumor');
props = regionprops('table',tumor,'all');




hold on
for i=1:length(B)
    plot(B{i}(:,2),B{i}(:,1), 'y' ,'linewidth',1.45);
    
    
    
end

%title('Detected Tumor');
hold off;

%tumor = regionprops(e,'Area','BoundingBox'); 
%area_valut= [tumor.Area]
%set(handles.area2 , 'string',areavalut);


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handleguide to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global tumor 
g = regionprops(tumor,'Area'); 
area_valut= [g.Area]
set(handles.area2 , 'string',area_valut);
%area =bwarea(tumor) % this is another method to mesur the area
%------------------------- count white pixel--------------------- 
[b w]=countBW(tumor)

set(handles.whitePixel , 'string',w);
%set(handles.blackPixel , 'string',b);
%----------%compute the entropy---------------------------- 
en = entropy(tumor);
set(handles.entropy2 , 'string',en);
%set(handles.whitePixel , 'string',meanG);
%the mean 
meanval = mean2(tumor)
set(handles.mean , 'string',meanval);
%std 
stdval = std2(tumor)
set(handles.stander , 'string',stdval);
% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)global img1
global img1
axes(handles.axes6)
if size(img1,3)==3
    img1=rgb2gray(img1);
end
K = medfilt2(img1);
axes(handles.axes6);
imshow(K);title('\fontsize{20}|-\color[rgb]{1,0,1}Med Filter');
img1=K;

%.............

 


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global img1
axes(handles.axes7)
if size(img1,3)==3
    img1=rgb2gray(img1);
end

K = wiener2(img1);
axes(handles.axes7);
imshow(K);
title('\fontsize{20}|-\color[rgb]{1,0,1}weiner Filter');


% --- Executes on button press in watershed.
% function watershed_Callback(hObject, eventdata, handles)
% % hObject    handle to watershed (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% global img1 I4
% axes(handles.axes8)
% axes(handles.axes9)
% axes(handles.axes10)
% 
% I = rgb2gray(img1);
% %imshow(I)
% 
% text(732,501,'Image courtesy of Corel(R)',...
%      'FontSize',7,'HorizontalAlignment','right')
%  
%  hy = fspecial('sobel');
% hx = hy';
% Iy = imfilter(double(I), hy, 'replicate');
% Ix = imfilter(double(I), hx, 'replicate');
% gradmag = sqrt(Ix.^2 + Iy.^2);
% %figure
% %imshow(gradmag,[]), title('Gradient magnitude (gradmag)')
% %*********
% L = watershed(gradmag);
% Lrgb = label2rgb(L);
% %figure, imshow(Lrgb), title('Watershed transform of gradient magnitude (Lrgb)')
% %======
% se = strel('disk', 20);
% Io = imopen(I, se);
% %figure
% %imshow(Io), title('Opening (Io)')
% %============
% Ie = imerode(I, se);
% Iobr = imreconstruct(Ie, I);
% %figure
% %imshow(Iobr), title('Opening-by-reconstruction (Iobr)')
% %============
% Ioc = imclose(Io, se);
% %figure
% %imshow(Ioc), title('Opening-closing (Ioc)')
% %============
% Iobrd = imdilate(Iobr, se);
% Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
% Iobrcbr = imcomplement(Iobrcbr);
% %figure
% %imshow(Iobrcbr), title('Opening-closing by reconstruction (Iobrcbr)')
% %==========
% fgm = imregionalmax(Iobrcbr);
% %figure
% %imshow(fgm), title('Regional maxima of opening-closing by reconstruction (fgm)')
% %.............
% I2 = I;
% I2(fgm) = 255;
% %figure
% %imshow(I2), title('Regional maxima superimposed on original image (I2)')
% %==========
% se2 = strel(ones(5,5));
% fgm2 = imclose(fgm, se2);
% fgm3 = imerode(fgm2, se2);
% %=============
% fgm4 = bwareaopen(fgm3, 20);
% I3 = I;
% I3(fgm4) = 255;
% %figure
% %imshow(I3)
% title('Modified regional maxima superimposed on original image (fgm4)')
% %==============
% % bw = imbinarize(Iobrcbr);
% % figure
% % imshow(bw), title('Thresholded opening-closing by reconstruction (bw)')
% %==============
% D = bwdist(I3);
% DL = watershed(D);
% bgm = DL == 0;
% %figure
% %imshow(bgm), title('Watershed ridge lines (bgm)')
% %===============
% gradmag2 = imimposemin(gradmag, bgm | fgm4);
% %==============
% L = watershed(gradmag2);
% %===========
% I4 = I;
% I4(imdilate(L == 0, ones(3, 3)) | bgm | fgm4) = 255;
%  %figure
%  axes(handles.axes8)
% imshow(I4)
% title('Markers and object boundaries superimposed on original image (I4)')
% %=================
% Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
% %figure
% axes(handles.axes9)
% imshow(Lrgb)
% title('Colored watershed label matrix (Lrgb)')
% %=================
% figure
% imshow(I)
% 
% hold on
% %axes(handles.axes10)
% himage = imshow(Lrgb);
% 
%  himage.AlphaData = 0.3;
% title('Lrgb superimposed transparently on original image')
