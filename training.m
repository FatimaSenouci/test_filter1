outputFolder =fullfile('Data')
rootFolder = fullfile(outputFolder,'DataEndo'); % load data 
categories = {'endomettreSAIN','endomettreMalade'}; %we are choosing 2 categories 
imds = imageDatastore(fullfile(rootFolder,categories),'LabelSource','foldernames','readfcn',@im_pro); %create image data stor to manage data so the images in categories are now in imds 