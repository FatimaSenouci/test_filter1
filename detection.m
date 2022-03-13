function varargout = detection(varargin)
% DETECTION MATLAB code for detection.fig
%      DETECTION, by itself, creates a new DETECTION or raises the existing
%      singleton*.
%
%      H = DETECTION returns the handle to a new DETECTION or the handle to
%      the existing singleton*.
%
%      DETECTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DETECTION.M with the given input arguments.
%
%      DETECTION('Property','Value',...) creates a new DETECTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before detection_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to detection_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help detection

% Last Modified by GUIDE v2.5 09-Aug-2021 00:29:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @detection_OpeningFcn, ...
                   'gui_OutputFcn',  @detection_OutputFcn, ...
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


% --- Executes just before detection is made visible.
function detection_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to detection (see VARARGIN)

% Choose default command line output for detection
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes detection wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = detection_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


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
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global img1
axes(handles.axes2)
if size(img1,3)==3
    img1=rgb2gray(img1);
end
K = medfilt2(img1);
axes(handles.axes2); 
imshow(K);title('\fontsize{20}|-\color[rgb]{1,0,1}Med Filter');


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global img1  
axes(handles.axes4)
%img1 = ones(3,3);

%structuring element 
se1 = strel('disk',11) ; 
ero = imerode(img1,se1); 

%K = ero;
axes(handles.axes4);
imshow(ero);title('\fontsize{20}|-\color[rgb]{1,0,1}erosion ');
img1 = ero  ;


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global img1
axes(handles.axes5)
%structuring element 
se2 = strel('squar',11); 
del = imdilate(img1,se2);



axes(handles.axes5);
imshow(del);title('\fontsize{20}|-\color[rgb]{1,0,1} delatation ');

img1= del ;

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global img1
axes(handles.axes6)
axes(handles.axes7)


[Gmag, Gdir] = imgradient(img1,'prewitt');
axes(handles.axes6);
imshow(Gmag); title('\fontsize{15}|-\color[rgb]{1,0,1}Gradient Magnitude, Gmag ');
axes(handles.axes7)
imshow(Gdir); title('\fontsize{15}|-\color[rgb]{1,0,1}Gradient direction, Gdir ');

%imshowpair(Gmag, Gdir, 'montage');
%title('\fontsize{20}|-\color[rgb]{1,0,1}Gradient Magnitude, Gmag ');
%title('Gradient Magnitude, Gmag (left), and Gradient Direction, Gdir (right), using Prewitt method')
%axes(handles.axes6);
%imshow(Gmag);title('\fontsize{20}|-\color[rgb]{1,0,1}Med Filter');


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global img1
axes(handles.axes8)
axes(handles.axes9)


[Gx,Gy] = imgradientxy(img1);
axes(handles.axes8);
imshow(Gx); title('\fontsize{15}|-\color[rgb]{1,0,1}, Gmag ');
axes(handles.axes9)
imshow(Gy); title('\fontsize{15}|-\color[rgb]{1,0,1} Gy ');


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global img1 tumor_label
axes(handles.axes13)
axes(handles.axes14)
axes(handles.axes15)

%img=imread('C:/cancer.jpg');
bw=im2bw(img1,0.7);
label=bwlabel(bw);
stats=regionprops(label,'Solidity','Area');
density=[stats.Solidity];
area=[stats.Area];
high_dense_area=density>0.5;
max_area=max(area(high_dense_area));
tumor_label=find(area==max_area);
tumor=ismember(label,tumor_label);
%structur element 
se=strel('square',5);
tumor=imdilate(tumor,se); % delatation of image 
%figure(2);
%subplot(1,3,1);
%imshow(img,[]);
%title('Brain');

%subplot(1,3,2);
%imshow(tumor,[]);
%title('Tumor Alone');
[B,L]=bwboundaries(tumor,'noholes');

axes(handles.axes13);
imshow(img1); title('\fontsize{15}|-\color[rgb]{1,0,1}, andometrial ');
axes(handles.axes14)
imshow(tumor); title('\fontsize{15}|-\color[rgb]{1,0,1} Tumor Alone ');
axes(handles.axes15);
imshow(img1); title('\fontsize{15}|-\color[rgb]{1,0,1},Detected Tumor');




%subplot(1,3,3);
%imshow(img,[]);
hold on
for i=1:length(B)
    plot(B{i}(:,2),B{i}(:,1), 'y' ,'linewidth',1.45);
    
    
    
end

%title('Detected Tumor');
hold off;

%axes(handles.axes13);
%imshow(img1); title('\fontsize{10}|-\color[rgb]{1,0,1}, andometrial ');
%axes(handles.axes14)
%imshow(tumor); title('\fontsize{20}|-\color[rgb]{1,0,1} Tumor Alone ');
%axes(handles.axes15);
%imshow(img1); title('\fontsize{10}|-\color[rgb]{1,0,1},Detected Tumor');
