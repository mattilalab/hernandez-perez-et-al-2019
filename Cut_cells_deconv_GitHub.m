% Andreas Bruckbauer 11/11/2018
% cluster analysis for Mattila lab
% manually cuts single cell data and saves as ome.tiff
% No background subtraction and Gaussian filtering
% This is for deconvolved files 

clear all;
close all;

% main_dir = 'the path of the folder with your images';
% save_dir = 'the path where you want to save your images';
% experiment = 'name of your exp';
% regex = 'cmle.tif';

main_dir = '';
save_dir = '';
experiment = '';
regex = 'cmle.tif';


display_adjust1 = 0.2;        % for display only
display_adjust2 = 0.1;       % for display only
cell_size = 10;            % approx size of cells to determine square region

% output text file
results_file = [main_dir '\' experiment '_results.txt'];
fid1 = fopen(results_file,'w');

f = get_file_list(main_dir,regex);

file_counter = 0;
for i = 1:size(f,1)
    if (exist([main_dir '\' f(i).name]) == 2)
        fprintf(1,'%s\n',f(i).name);
        fprintf(fid1,'%s\n',f(i).name);            
        file_counter = file_counter+1;
    end
end

% open file for channel 1
index = 1;
filename = strcat(main_dir,'\',f(index).name, f(index).ext);
fprintf(1,'\nfile %i: \t %s\n',index,filename);
fprintf(fid1,'\nfile %i: \t %s\n',index,filename);
name = f(index).name;

data = bfopen(filename);
omeMeta = data{1, 4};
pixel_size = omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER);
pixel_size = double(pixel_size);
pixelSizeZ = omeMeta.getPixelsPhysicalSizeZ(0).value(ome.units.UNITS.MICROMETER);
pixel_sizeZ = double(pixelSizeZ);

series1 = data{1,1};
number_planes = size(series1,1);
for i = 1:round(number_planes)
    im(:,:,i) = double(series1{i, 1});
end
max_proj_image = max(im, [], 3);
max_proj_image_norm = max_proj_image/max(max(max_proj_image));
figure('Position',[50 50 800 800]);
imagesc(imadjust(max_proj_image_norm,[0 display_adjust1], [0 1]));
colormap gray;
axis image;
title(['Ch1 maximum intensity projection']);
hold on;
pause(0.1);
savename = [main_dir '\' name '_ch1_max_proj.png'];
set(gcf,'PaperPositionMode','auto');
print('-dpng',savename);

% open file for channel 2
index = 2;
filename = strcat(main_dir,'\',f(index).name, f(index).ext);
fprintf(1,'\nfile %i: \t %s\n',index,filename);
fprintf(fid1,'\nfile %i: \t %s\n',index,filename);
data2 = bfopen(filename);
series2 = data2{1,1};
number_planes = size(series2,1);
for i = 1:number_planes
    im2(:,:,i) = double(series2{i, 1});
end
avg_proj_image = mean(im2,3);
avg_proj_image_norm = avg_proj_image/max(max(avg_proj_image));
figure('Position',[870 50 800 800]);
imagesc(imadjust(avg_proj_image_norm,[0 display_adjust2], [0 1]));
colormap gray;
axis image;
title(['Ch2 average intensity projection']);
hold on;
real_number_planes = round(number_planes);

% select cells
square_size = round(sqrt(2)*cell_size/pixel_size);
count_param = 1;
number_regions = 0;
M = [];

while (count_param == 1)
   choice = questdlg('Do you want to add a region?','','Yes');
   switch choice
   case 'Yes'
       number_regions = number_regions +1;
       h = imrect(gca, [400 400 square_size square_size]);
       setResizable(h,0)
       position = wait(h);
       M = [M; position];
       rectangle('Position',position,'LineWidth',2,'EdgeColor','w');
   case 'No'
      count_param = 0;
   case 'Cancel'
      count_param = 0;
   end 
end 

% number_regions = size(M,1);
% for k = 1: number_regions
%     position = M(k,:);
%      rectangle('Position',position,'LineWidth',2,'EdgeColor','w');
% end

savename = [save_dir '\' name '_ch2_avg_proj.png'];
set(gcf,'PaperPositionMode','auto');
print('-dpng',savename);


%%
for k = 1:number_regions

    clear RGB;
    position = M(k,:);
    % crop
    for i = 1:real_number_planes
        cropped_region_ch1(:,:,i) = imcrop(im(:,:,i),round(position));
        img_ch1_filter(:,:,i) = cropped_region_ch1(:,:,i); 
    end
    img_ch1_max_proj_image = max(cropped_region_ch1, [], 3);
    img_ch1_max_proj_norm = img_ch1_max_proj_image/max(max(img_ch1_max_proj_image));
        
    for i = 1:real_number_planes
        cropped_region_ch2(:,:,i) = imcrop(im2(:,:,i),round(position));
        img_ch2_filter(:,:,i) = cropped_region_ch2(:,:,i); 
    end
    img_ch2_max_proj_image = max(img_ch2_filter, [], 3);
    img_ch2_max_proj_norm = img_ch2_max_proj_image/max(max(img_ch2_max_proj_image));
    
    img_width = size(img_ch1_max_proj_norm,1);
    RGB = zeros(img_width,img_width,3);
    RGB(:,:,1) = img_ch1_max_proj_norm/0.8;
    RGB(:,:,2) = img_ch2_max_proj_norm/0.8;
            
    figure;
    imagesc(RGB);
    axis image;
    pause(0.1);
    title(['cell ' num2str(k)]);
    savename = [main_dir '\' name '_cell_' num2str(k) '_RGB.png'];
    set(gcf,'PaperPositionMode','auto');
    print('-dpng',savename);
    close;
    
    % save as OME tiff
    save_data(:,:,:,1) = uint16(img_ch1_filter);
    save_data(:,:,:,2) = uint16(img_ch2_filter);
    metadata = createMinimalOMEXMLMetadata(save_data,'XYZCT');
    pixelSize = ome.units.quantity.Length(java.lang.Double(pixel_size), ome.units.UNITS.MICROMETER);
    metadata.setPixelsPhysicalSizeX(pixelSize,0);
    metadata.setPixelsPhysicalSizeY(pixelSize,0);
    pixelSizeZ = ome.units.quantity.Length(java.lang.Double(pixel_sizeZ), ome.units.UNITS.MICROMETER);
    metadata.setPixelsPhysicalSizeZ(pixelSizeZ,0);
    save_name = [save_dir '\' name '_cell_' num2str(k) '.ome.tiff'];
    bfsave(save_data, save_name, 'metadata', metadata);
end

%%

fclose(fid1);

clear im;
clear im2;
clear omeMeta;
clear metadata;
clear pixelSize;
clear pixelSizeZ;
clear data;
clear data2;
clear series1;
clear series2;

%  save matlab file
matlabFile = [main_dir '\' experiment '.mat'];
save(matlabFile);   
   
fprintf('\nfinished!\n');
