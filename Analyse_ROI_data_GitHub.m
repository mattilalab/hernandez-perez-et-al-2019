% Andreas Bruckbauer 05/01/2019
% cluster analysis for Mattila lab
% analyses cut cells directly (no Icy)
% 3D distance between vesicles and mtoc
% uses thresh_tool from Brandon Kuczenski, see copyright notice


clear all;
close all;

% data_dir = 'path where your cropped images are';
% experiment = 'name of your experiment';

data_dir = '';
experiment = '';

display_adjust1 = 0.5;        % for display only Ch1
display_adjust2 = 0.2;        % for display only Ch2
var_threshold_1 = 0.7;        % for cluster detection, 1 = Otsu threshold
var_threshold_2 = 1.0;        % for MTOC detection, 1 = Otsu threshold
distance_threshold = 2;       % cutoff for central cluster counting in micrometer
show2D = 0;                   % show 3D results as 2D stack

% output text file
results_file = [data_dir '\' experiment '_results.txt'];
fid1 = fopen(results_file,'w');

f = get_file_list(data_dir,'ome.tiff');

file_counter = 0;
for i = 1:size(f,2)
    if (exist([data_dir '\' f(i).name f(i).ext]) == 2)
        fprintf(1,'%s\n',f(i).name);
        fprintf(fid1,'%s\n',f(i).name);            
        file_counter = file_counter+1;
    end
end

number_files = file_counter;
fprintf(1,'\n number of files: %i\n',number_files);
fprintf(fid1,'\n number of files: %i\n',number_files);

all_distances = [];

file_counter = 0;
for f_index = 1:number_files
    
    total_intensity1 = [];
    total_intensity2 = [];
    
    filename = [data_dir,'\',f(f_index).name, f(f_index).ext];
        
    data = bfopen(filename);
    omeMeta = data{1, 4};
    pixel_size = omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER);
    pixel_size = double(pixel_size);
    pixelSizeZ = omeMeta.getPixelsPhysicalSizeZ(0).value(ome.units.UNITS.MICROMETER);
    pixel_sizeZ = double(pixelSizeZ);
    
    series1 = data{1,1};
    number_planes = size(series1,1);

    for i = round(number_planes/2)+1:number_planes
        im2(:,:,i-round(number_planes/2)) = double(series1{i, 1});
    end
    avg_proj_image = max(im2,[], 3);
    image_size = size(avg_proj_image,1);
    max_ch2 = max(max(avg_proj_image));
    avg_proj_image_norm = avg_proj_image/max_ch2;
    avg_proj_image_norm_adj = imadjust(avg_proj_image_norm,[0 display_adjust2], [0 1]);
 
    % Ch1
    for i = 1:round(number_planes/2)
        im(:,:,i) = double(series1{i, 1});
    end
    max_proj_image = max(im, [], 3);
    
    max_ch1 = max(max(max_proj_image));
    max_proj_image_norm = max_proj_image/max_ch1;
    max_proj_image_norm_adj = imadjust(max_proj_image_norm,[0 display_adjust1], [0 1]);
    
    % plot RGB image
    img_width = size(max_proj_image,1);
    RGB = zeros(img_width,img_width,3);
    RGB(:,:,1) = max_proj_image_norm_adj;
    RGB(:,:,2) = avg_proj_image_norm_adj;
    figure('Position',[300, 300, 400, 400]);
    imagesc(RGB);
    axis image;
    hold on;
    pause(0.1);
    
    choice = questdlg('Analyse this cell?','','Yes');
    switch choice
        case 'Yes'
            
            file_counter = file_counter +1;
            fprintf(1,'\nfile %i: \t %s\n',f_index,filename);
            fprintf(fid1,'\nfile %i: \t %s\n\n',f_index,filename);
            close(gcf);
  
            % cluster analysis
            level= var_threshold_1 * graythresh(max_proj_image_norm);
            level = thresh_tool_mod(max_proj_image_norm,[],level);
            for i = 1:round(number_planes/2)
                image_norm(:,:,i) = im(:,:,i)/ max_ch1;
                mask(:,:,i) = imbinarize(image_norm(:,:,i),level);
                B{i} = bwboundaries(mask(:,:,i),'noholes');  % for 2D representation
            end

            % 3D stats
            stats = regionprops(mask,im,'Area','Centroid','Area','MeanIntensity');
            number_cluster1 = size(stats,1);
            if number_cluster1 > 0 
                for i = 1:number_cluster1
                    total_intensity1(i) = stats(i).MeanIntensity*stats(i).Area;
                end
            end

            % show results in 2D
            if show2D == 1
                figure;
                hold on;
                colormap gray;
                axis image;
                % plot bounderies
                for j = 1:round(number_planes/2)
                    display_image = imadjust(image_norm(:,:,j),[0 display_adjust1], [0 1]);
                    imagesc(display_image);
                    for i=1:length(B{j})
                        boundary = B{j}{i};
                        plot(boundary(:,2),boundary(:,1),'g','LineWidth',2);
                    end
                    pause(0.1);
                end
                close;
            end

            % MTOC, Ch2
            level2= var_threshold_2 * graythresh(avg_proj_image_norm);
            for i = 1:round(number_planes/2)
                image_norm2(:,:,i) = im2(:,:,i)/ max_ch2;
                mask2(:,:,i) = imbinarize(image_norm2(:,:,i),level2);
                B2{i} = bwboundaries(mask2(:,:,i),'noholes');  % for 2D representation
            end

            % 3D stats MTOC
            stats2 = regionprops(mask2,im2,'Area','Centroid','Area','MeanIntensity'); 
            number_cluster2 = size(stats2,1);
            if number_cluster2 > 0 
                for i = 1:number_cluster2
                    total_intensity2(i) = stats2(i).MeanIntensity*stats2(i).Area;
                end
            end

            I = find(total_intensity2 == max(total_intensity2));
            mtoc_centroid = stats2(I).Centroid;

            x_mtoc = mtoc_centroid(1);
            y_mtoc = mtoc_centroid(2);
            z_mtoc = mtoc_centroid(3);

            % plot RGB image
            img_width = size(max_proj_image,1);
            RGB = zeros(img_width,img_width,3);
            RGB(:,:,1) = max_proj_image_norm_adj;
            RGB(:,:,2) = avg_proj_image_norm_adj;
            figure;
            imagesc(RGB);
            axis image;
            hold on;
            pause(0.1);

            % plot mtoc
            plot(x_mtoc,y_mtoc ,'ocy','MarkerSize',10); % change size no 4  

            % center of cell
            e2 = imellipse(gca,[round(img_width/6) round(img_width/6) 100 100]);
            setFixedAspectRatioMode(e2,1);
            position2 = wait(e2);
            cell_x = mean(position2(:,1));
            cell_y = mean(position2(:,2));
            radius = sqrt((position2(1,1)-cell_x)^2+(position2(1,2)-cell_y)^2);

            % calculate distances and plot cluster
            weighted_distance = [];
            cluster_mtoc_distance = [];
            cluster_centre_distance = [];
            number_clusters = 0;
            for k = 1:size(stats,1)
                x_IgG = stats(k).Centroid(1);
                y_IgG = stats(k).Centroid(2);
                z_IgG = stats(k).Centroid(3);
                Total_Int_IgG = stats(k).MeanIntensity*stats(k).Area;
                Mean_Int_IgG = stats(k).MeanIntensity;
                Area_IgG = stats(k).Area;

                % calculate 2D distance to centre;
                cluster_centre_distance(k) = sqrt((x_IgG-cell_x)^2+(y_IgG-cell_y)^2);
                if cluster_centre_distance(k) <= radius

                    number_clusters = number_clusters +1;

                    % calculate distance
                    cluster_mtoc_distance(number_clusters) = sqrt((pixel_size*(x_IgG-x_mtoc))^2 ...
                    +(pixel_size*(y_IgG-y_mtoc))^2+(pixel_sizeZ*(z_IgG-z_mtoc))^2);

                    weighted_distance(number_clusters) = cluster_mtoc_distance(number_clusters)*Total_Int_IgG;

                    total_intensity(number_clusters) = Total_Int_IgG;

                    all_distances = [all_distances; f_index, number_clusters, cluster_mtoc_distance(number_clusters),...
                         Area_IgG, Mean_Int_IgG];

                    % plot IgG cluster
                    plot(x_IgG,y_IgG ,'or','MarkerSize',10); % change size no 4

                    % plot line to show distance
                    plot([x_mtoc x_IgG],[y_mtoc y_IgG]','cy-');
                end
            end
            weighted_cluster_distance(file_counter) = sum(weighted_distance)/sum(total_intensity);
            avg_distance(file_counter) = mean(cluster_mtoc_distance); 
            I = find(cluster_mtoc_distance <= distance_threshold);
            number_central_clusters(file_counter) = size(I,2);

            title(['cell: ' num2str(file_counter) ' dist: ' num2str(avg_distance(file_counter),'%.3g') ' weight. dist: ' num2str(weighted_cluster_distance(file_counter),'%.3g')]);
            pause(0.2);
            savename = [data_dir '\' f(file_counter).name '_dist_RGB.png'];
            set(gcf,'PaperPositionMode','auto');
            print('-dpng',savename);

            close(gcf);
            
            case 'No'
                 close(gcf);
            case 'cancel'
                close(gcf);
     end 
end
      
avg_distance_overall = mean(avg_distance);
fprintf(1,'\nAverage distance over all clusters and cells (um): %f\n',avg_distance_overall);
fprintf(fid1,'\nAverage distance over all clusters and cells (um): %f\n',avg_distance_overall);

avg_weighted_distance_overall = mean(weighted_cluster_distance);
fprintf(1,'Average intensity weighted distance over all clusters and cells (um): %f\n',avg_weighted_distance_overall);
fprintf(fid1,'Average intensity weighted distance over all clusters and cells (um): %f\n',avg_weighted_distance_overall);

avg_number_central_clusters = mean(number_central_clusters);
fprintf(1,'Average number central clusters over all cells: %f\n',avg_number_central_clusters);
fprintf(fid1,'Average number central clusters over all cells: %f\n',avg_number_central_clusters);

% create file with all results
savename = [data_dir '\' experiment '_all_values.txt'];
fid2 = fopen(savename,'w');
fprintf(fid2,'filename, file_number, cluster_number, cluster_mtoc_distance, volume_cluster, mean_int_cluster, total_int_cluster, total_int_all_clusters_cell\n');
for i =1:size(all_distances,1)
    % find all clusters from same cell
    I = find(all_distances(:,1) ==  all_distances(i,1));
    cell_int_all_cluster = sum(all_distances(I,4).*all_distances(I,5));
    fprintf(fid2,'%s, %i, %i, %f, %i, %f, %f, %f\n', f(all_distances(i,1)).name, all_distances(i,1), ...
        all_distances(i,2), all_distances(i,3), all_distances(i,4), all_distances(i,5),...
        all_distances(i,4).*all_distances(i,5),cell_int_all_cluster);
end
fclose(fid2);


%% generate plots

figure;
h1 = histogram(all_distances(:,3),'Normalization','count');
h1.NumBins = 25;
xlabel('Distance MTOC / um');
ylabel('counts');
title(experiment,'Interpreter', 'none');
pause(0.1);
savename = [data_dir '\'  experiment '_histogram_distance.png'];
set(gcf,'PaperPositionMode','auto');
print('-dpng',savename);

figure;
plot(all_distances(:,3),all_distances(:,4),'bo','MarkerSize',3,'MarkerFaceColor','b');
xlabel('Distance / um');
ylabel('Volume / pixel');
title(experiment,'Interpreter', 'none');
savename = [data_dir '\'  experiment '_Distance_Volume.png'];
set(gcf,'PaperPositionMode','auto');
print('-dpng',savename);

figure;
plot(all_distances(:,3),all_distances(:,5),'bo','MarkerSize',3,'MarkerFaceColor','b');
xlabel('Distance / um');
ylabel('Cluster mean intensity / a.u');
title(experiment,'Interpreter', 'none');
savename = [data_dir '\'  experiment '_Distance_Intensity.png'];
set(gcf,'PaperPositionMode','auto');
print('-dpng',savename);

figure;
plot(all_distances(:,3),all_distances(:,4).*all_distances(:,5),'bo','MarkerSize',3,'MarkerFaceColor','b');
xlabel('Distance / um');
ylabel('Cluster total intensity / a.u');
title(experiment,'Interpreter', 'none');
savename = [data_dir '\'  experiment '_Distance_Total_Intensity.png'];
set(gcf,'PaperPositionMode','auto');
print('-dpng',savename);

fclose(fid1);

clear omeMeta;
clear metadata;
clear pixelSizeZ

%  save matlab file
matlabFile = [data_dir '\' experiment '.mat'];
save(matlabFile);   
   
fprintf('\nfinished!\n');


