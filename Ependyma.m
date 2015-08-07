%% Ependyma is the main function. Calculates angles of basal body orientation vectors (BBOV) (Figure 2)

% 1. Create ROIs of cell outlines for every image
% 2. Save ROIs as zip files (names have to correspond to tiff files)
% 3. Run the script

%% Additional functions used and included:
% combinecsv - combines csv files from a folder into a single cell aray
% cell2double - converts cell array of double into an array of doubles
% conv16_8 - converts 16-bit image into 8-bit and scales based on brightest pixel
% drawfreehandzones - overlays freehand masks generated in ImageJ
% drawzones - overlays line masks generated in ImageJ
% exclude_str_cell - excludes elements of cell array, that match a criterium
% findstr_cell - outputs indices of cell of elements of cell array, that match criterium
% findstr_cell_partial - outputs indices of cell of elements of cell array, that partially match criterium
% histc_full - modified histc to include the elements outside of bin range

%% Functions used and NOT included. Can be downloaded from FileExchange:
% export_fig						http://www.mathworks.com/matlabcentral/fileexchange/23629-export-fig
% Circular Statistics Toolbox       http://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox--directional-statistics-
% arrow                             http://www.mathworks.com/matlabcentral/fileexchange/278-arrow-m
% CSVIMPORT                         http://www.mathworks.com/matlabcentral/fileexchange/23573-csvimport
% linepoints						http://www.mathworks.com/matlabcentral/fileexchange/1853-linepoints
% rdir                              http://www.mathworks.com/matlabcentral/fileexchange/47125-rdir-m
% ReadImageJROI 					http://www.mathworks.com/matlabcentral/fileexchange/32479-load---import-imagej-roi-files


%% Ependyma
function Ependyma(top_folder,varargin)
% Ependyma - calculates angles of basal body orientation vectors (BBOV) 
% top_folder - top folder with subfolder containing annotated tiff files
% bb_size_thresh - maximum allowed ratio of BB to cell area
% all_data_sets - names of the data-sets
% header - excel header

if nargin>1
    bb_size_thresh=varargin{1};
else
    bb_size_thresh=0.2;
end

if nargin>2
    all_data_sets=varargin{2};
else
    all_data_sets={'WT','HET','KD'};
end

if nargin>3
    header=varargin{3};
else
    header={'','Vector length/cell length (absolute distribution)','Vector length/cell length (relative distribution)','Vector projection length/cell length (absolute distribution)','Vector projection length/cell length (relative distribution)'};
end

%% list folders
folders=rdir([top_folder,'\*']);
dircheck=[folders.isdir];
folders={folders.name};
folders=folders(dircheck);
[selectedFolders,ok]=listdlg('ListString',folders,'SelectionMode','multiple'); %selection is an array of double

if ok
   selectedFolders=folders(selectedFolders)';
else
    msgbox('No folders were selected','Exiting...');
    return
end

    for data_dir=1:size(selectedFolders,1)
        input_dir=selectedFolders{data_dir,1};    
        output_dir=[input_dir,'zoned\'];
        mkdir(output_dir);
        files=dir([input_dir,'\*.tiff']);
        files={files.name};
    %% Cropping Regions
        for f=1:size(files,2)
            imfile=[input_dir,files{1,f}];
            zipname=strrep(imfile,'.tiff','.zip');
            if exist(zipname,'file')==0
                continue
            end
            unzip(zipname,input_dir);
            roi_files=dir([input_dir,'*.roi']);
            roi_files={roi_files.name};
            image_info = imfinfo(imfile);
            num_images = numel(image_info);
            input_images=uint16([]);
            for p=1:num_images
                im=imread(imfile,p);
                input_images=cat(3,input_images,im);
            end

            ROI=cell(size(roi_files,2),2);
            for r=1:size(roi_files,2)
                [sROI] = ReadImageJROI([input_dir,roi_files{1,r}]);
                ROI{r,1}=sROI.mnCoordinates; 
                ROI{r,2}=sROI.vnRectBounds;
                delete([input_dir,roi_files{1,r}]);
            end

            if strcmpi(sROI.strType,'Freehand')
                im_array  = drawfreehandzones( ROI, input_images );
            else
                im_array  = drawzones( ROI, input_images);
            end
            imwrite(im_array{1,3},[output_dir,files{1,f}],'tif');
             for z=2:size(im_array,1)
                    imwrite(im_array{z,3},[output_dir,files{1,f}],'tif','WriteMode','append');
             end
        end

    %% Segmentation
        files=dir([input_dir,'zoned\*.tiff']);
        files={files.name};
        files2segment_ind=exclude_str_cell('_segmented.tif',files);
        files=files(1,files2segment_ind);
        fid=fopen([input_dir,'detailed_',num2str(data_dir),'.csv'],'w');
        fprintf(fid,'File,Cell,BB Area,Area Ratio,X(Cell),Y(Cell),X(BB),Y(BB),Angle,Length,Corrected Length(cos),Cells long axis length,Cell Orientation,Basal body position\n');
        fid_sum=fopen([input_dir,'summary_',num2str(data_dir),'.csv'],'w');
        fprintf(fid_sum,'File,Number of cells,Mean Angle,Std Angle,Cell length,Distally-positioned basal body,Proximally-positioned basal body\n');

        for i=1:size(files,2)
            angle_array=[]; cell_size_arr=[]; bb_position_array={};
            curr_im=imfinfo([input_dir,'zoned\',files{1,i}]);
            for fr=1:numel(curr_im)
                im=imread([input_dir,'zoned\',files{1,i}],'Index',fr); 
                im_8=conv16_8(im);
                borders_celll=bwperim(im);
                BWdfill_cell = imfill(borders_celll, 'holes');
                cell_detected=bwconncomp(BWdfill_cell);
                cell_props=regionprops(cell_detected,'Area','Centroid','MajorAxisLength','Orientation');
                %cell_props=regionprops(cell_detected,'all');            
                area=[cell_props.Area];
                [~,ind]=sort(area,'descend');
                cell_centre=cell_props(ind(1)).Centroid; 
                cell_centre=fliplr(cell_centre);
                cell_length=cell_props(ind(1)).MajorAxisLength; 
                cell_size_arr=cat(1,cell_size_arr,cell_length);
                cell_orient=cell_props(ind(1)).Orientation;
                BWdfill_cell = ismember(labelmatrix(cell_detected),ind(1));

    %% Rotating the mask to identify proximal/distal parts in the image            
                mask_rotated=imrotate(BWdfill_cell,-cell_orient,'loose');  
                cell_props_rotated=regionprops(mask_rotated,'Area','BoundingBox');
                area=[cell_props_rotated.Area];
                [~,ind]=sort(area,'descend');
                cell_area=area(ind(1));

    %% Create cell's sub-regions            
                thirds=[ceil(cell_props_rotated(ind(1)).BoundingBox(1)),round(cell_props_rotated(ind(1),1).BoundingBox(1)+cell_props_rotated(ind(1),1).BoundingBox(3)*3/8),round(cell_props_rotated(ind(1),1).BoundingBox(1)+cell_props_rotated(ind(1),1).BoundingBox(3)*5/8),floor(cell_props_rotated(ind(1),1).BoundingBox(1)+cell_props_rotated(ind(1),1).BoundingBox(3))];
                sub_cell_masks=cell(1,3);
                for thr=1:3
                    masked=mask_rotated;
                    if thr==1 && thirds(1)<2                   
                        masked(:,thirds(thr+1):size(mask_rotated,2))=0;
                    else
                        masked(:,thirds(thr+1):size(mask_rotated,2))=0;
                        masked(:,1:thirds(thr))=0;
                    end
                    if thr==3 && thirds(4)==size(mask_rotated,2)
                        masked(:,1:thirds(thr))=0;
                    else
                        masked(:,thirds(thr+1):size(mask_rotated,2))=0;
                        masked(:,1:thirds(thr))=0;
                    end                
                    sub_cell_masks{1,thr}=masked;
                end
                cell_parts={'Distal','Proximal','Distal'};            
                im_median=medfilt2(im);
                lvl=graythresh(im_median);
                BW=im2bw(im_median,lvl);
                BWdfill = imfill(BW, 'holes');
                BWclosed=imclose(BWdfill,strel('disk',6));
                objects=bwconncomp(BWclosed);                   
                props=regionprops(objects,'Area','Centroid');
                area=[props.Area];
                [~,ind]=sort(area,'descend');
                bb_area=area(ind(1));
                area_ratio=bb_area/cell_area;
                if area_ratio>bb_size_thresh
                    if lvl>0
                        lvl=lvl*1.3;
                    else
                       lvl=0.001; 
                    end
                    while area_ratio>bb_size_thresh
                        BW=im2bw(im_median,lvl);
                        BWdfill = imfill(BW, 'holes');
                        BWclosed=imclose(BWdfill,strel('disk',6));
                        objects=bwconncomp(BWclosed);    
                        props=regionprops(objects,'Area','Centroid');
                        area=[props.Area];
                        [~,ind]=sort(area,'descend');
                        bb_area=area(ind(1));
                        area_ratio=bb_area/cell_area;
                        lvl=lvl*1.3;
                    end
                end            
                bb_centre=props(ind(1)).Centroid;
                bb_centre=fliplr(bb_centre);
                bb = ismember(labelmatrix(objects),ind(1));
                bb = imfill(bb, 'holes');
                imshow(bb);figure(gcf);
                bb_mask_rotated=imrotate(bb,-cell_orient,'loose');
                olap_percent=zeros(1,3);
                for thr=1:3
                    olap=and(sub_cell_masks{1,thr},bb_mask_rotated);
                    olap_area=size(find(olap>0),1);
                    olap_percent(1,thr)=olap_area/bb_area;
                end
                [~,max_ind]=max(olap_percent);
                bb_position=cell_parts{1,max_ind};
                bb_position_array=cat(1,bb_position_array,bb_position);            
    %% Color output generation            
                borders=bwperim(bb);
                borders_8=255*uint8(borders);
                green=cat(3,borders_8,im_8);
                green=max(green,[],3);
                red=im_8;
                if max_ind==2
                    segmented=cat(3,red,green,im_8);
                else
                    segmented=cat(3,red,green,im_8);
                end
                imshow(segmented);figure(gcf);
                arrow(fliplr(round(cell_centre)),fliplr(round(bb_centre)),3,'Width',0.5, 'EdgeColor', 'r', 'FaceColor', 'r');
                export_fig('temp.tif','-m4');
                segmented=imread('temp.tif');
                segmented=imresize(segmented,size(red,1)/size(segmented,1));
                delete('temp.tif');

                len=pdist2(cell_centre,bb_centre);
                y=-(bb_centre(1)-cell_centre(1));
                x=(bb_centre(2)-cell_centre(2));            
                angle=atan2d(y,x);
                angle_array=cat(1,angle_array,angle);
                angle_deg=atand(y/x);
                dif_angle=degtorad(cell_orient-angle_deg);
                corrected_length=abs(cos(dif_angle))*len;
                fprintf(fid,'%s,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%s\n',files{1,i},fr,area(ind(1)),area_ratio,cell_centre(2),cell_centre(1),bb_centre(2),bb_centre(1),angle,len,corrected_length,cell_length,cell_orient,bb_position);
                if fr>1
                    imwrite(segmented,[input_dir,'zoned\',strrep(files{1,i},'.tif','_segmented.tif')],'tif','WriteMode','append'); 
                else
                    imwrite(segmented,[input_dir,'zoned\',strrep(files{1,i},'.tif','_segmented.tif')],'tif'); 
                end
            end
            angle_array = circ_ang2rad(angle_array);
            mean_angle=round(wrapTo360(radtodeg(circ_mean(angle_array))));
            std_angle=radtodeg(circ_std(angle_array));
            dist_bb=findstr_cell('Distal',bb_position_array);
            prox_bb=findstr_cell('Proximal',bb_position_array);
            fprintf(fid_sum,'%s,%d,%f,%f,%f,%d,%d\n',files{1,i},size(angle_array,1),mean_angle,std_angle,mean(cell_size_arr),size(dist_bb,1),size(prox_bb,1));
        end

    end

    dta=combinecsv(top_folder,'*detailed_*.csv',[9:12,14]);         %scanning folder for csv files and extracting angles, Length, Corrected Length(cos), Cells long axis length
    all_files=dta{1,1}(:,1);
    dist_prox_all=dta{1,1}(:,2);
    summary={'Type','Cell Number','Circular Standard Deviation (mean)','Circular Standard Deviation (stdev)','Number of proximally-positioned bb','Number of distally-positioned bb'};
    for ds=1:numel(all_data_sets)

        dataset=all_data_sets{1,ds};
        files_set_ind=findstr_cell_partial(dataset,all_files');
        if isempty(files_set_ind)
           continue; 
        end
        files=all_files(files_set_ind);
        dist_prox=dist_prox_all(files_set_ind,:);
        prox_ind=findstr_cell('Proximal',dist_prox);
        %input_folder=[top_folder,dataset,'\'];
        out_folder=top_folder;
        %file2load='i:\Wrana lab\Ira\IF108_1607_EMB1-7\EMB1_KD_detailed.csv';
        un_files=unique(files);
        angles=cell(size(un_files,1),6);
        for i=1:size(un_files,1)
            ind=findstr_cell(un_files{i,1},all_files);
            angles{i,1}=dta{1,2}(ind,1);                        %all angles
            ang_rad=circ_ang2rad(angles{i,1});
            angles{i,2}=circ_rad2ang(circ_mean(ang_rad));       %mean angle for a given field
            angles{i,3}=circ_rad2ang(circ_std(ang_rad));        %std angle for a given field
            ang_diff=90-angles{i,2};
            corrected_angles=wrapTo360(angles{i,1}+ang_diff);
            angles{i,4}=corrected_angles;                       %corrected angles
            angles{1,5}=cat(1,angles{1,5},corrected_angles);    %combined corrected angles
        end
        angles{1,6}=(2*dta{1,2}(:,2))./dta{1,2}(:,4);       %calculating ratio of distance to cilia from the centre of mass to long cell's axis
        angles{1,7}=(2*dta{1,2}(:,3))./dta{1,2}(:,4);       %calculating ratio of angle-corrected distance to cilia from the centre of mass to long cell's axis
        detailed_out=cat(2,un_files,angles(:,3));
        all_std=cell2double(angles(:,3));
        mean_std=mean(all_std);
        std_std=std(all_std);
        disp([num2str(mean_std),' +/- ',num2str(std_std)]);
        disp(size(angles{1,5},1));
        curr_summary=cat(2,all_data_sets(1,ds),num2cell(numel(files_set_ind)),num2cell(mean_std),num2cell(std_std),num2cell(numel(prox_ind)),num2cell(numel(files_set_ind)-numel(prox_ind)));
        summary=cat(1,summary,curr_summary);
        ang_rad=wrapToPi(circ_ang2rad(angles{1,5}));
        edg=-pi:pi/8:pi;
        graph_id=rose(ang_rad,edg);
        x = get(graph_id,'Xdata');
        y = get(graph_id,'Ydata');
        g=patch(x,y,'y','FaceColor','b');
        saveas(graph_id,[out_folder,dataset,'_graph.jpg']);
        hHiddenText = findall(gca,'type','text');
        f_col=repmat({'off'},numel(hHiddenText),1);
        set(hHiddenText,{'Visible'},f_col);
        saveas(graph_id,[out_folder,dataset,'_normalized_angles_notext.jpg']);

        %close(gcf);
        rat_edg=0:0.1:1;
        ratio_binned=histc_full(angles{1,6},rat_edg);
        ratio_binned_rel=ratio_binned./sum(ratio_binned);
        ratio_binned2=histc_full(angles{1,7},rat_edg);
        ratio_binned_rel2=ratio_binned2./sum(ratio_binned2);
        header{1,1}=all_data_sets{1,ds}; 
        xls_out=num2cell(cat(2,rat_edg',ratio_binned,ratio_binned_rel,ratio_binned2,ratio_binned_rel2));
        xls_out=cat(1,header,xls_out);

        xlswrite([out_folder,'results.xlsx'],xls_out,[all_data_sets{1,ds},' Binned']);
        xlswrite([out_folder,'results.xlsx'],detailed_out,[all_data_sets{1,ds},' Detailed']);
    end
    xlswrite([out_folder,'results.xlsx'],summary,'Summary');
    fclose all;
end


%% combinecsv
function all_data = combinecsv( input_folder,fileFilter,columns)
% combinecsv combines all csv-files from a folder into a single cell aray
% input_folder - folder with csv files
% fileFilter - wildcard filter to extract files
% columns - indices of columns to extract
% Assumes that first colum contains text data (file names) and there are
% no more text columns
%{1,1} - file names 
%{1,2} - selected numerical data

    files=dir([input_folder,'\',fileFilter]);
    files={files.name};
    all_data=cell(1,2);
    for i=1:size(files,2)
       cur_data=csvimport([input_folder,files{1,i}]);
       num_cols=[];txt_cols=[];
       for col=1:size(cur_data,2)
            if(ischar(cur_data{2,col}));
                txt_cols=cat(1,txt_cols,col);
            else
                num_cols=cat(1,num_cols,col);
            end       
       end
       data=cell2double(cur_data(2:size(cur_data,1),num_cols));
       %textdata=cur_data(2:size(cur_data,1),txt_cols);

       txt_needed=intersect(txt_cols,columns);   
       if isempty(txt_needed)
           num_cols=columns;
           txt_cols=1;
       else
           txt_cols=cat(2,1,txt_needed);
           num_cols = setdiff(columns,txt_needed)-1;
       end

       all_data{1,1}=cat(1,all_data{1,1},cur_data(2:size(cur_data,1),txt_cols));
       all_data{1,2}=cat(1,all_data{1,2},data(:,num_cols));
    end

end

%% cell2double
function [ Y ] = cell2double( X )
%X - Cell array with numbers
%Y - array of doubles

    Y=zeros(size(X,1),size(X,2));

    for c=1:size(X,2)
       for r=1:size(X,1)
           if ischar(X{r,c})
               Y(r,c)=str2double(X{r,c});
           elseif isempty(X{r,c})
               Y(r,c)=nan;
           else
               Y(r,c)=X{r,c}(1);
           end
       end
    end

end

%% Conversion of 16-bit to 8-bit image
function [ i ] = conv16_8( I )
%conv16_8 converts 16-bit image into 8-bit and scales it
%   I - 16-bit image
%   i - 8-bit image

           doub_image=double(I);
           doub_image=doub_image/max(max(doub_image))*255;
           i=uint8(doub_image);
           
end

%% drawfreehandzones
function [ im_array ] = drawfreehandzones( coord, I)
%drawfreehandzones draws lines on an image, based on coordinates from coordinates cell
%array (z,1). z - number of lines that were drawn;
%output is an array with "zoned" images (column 2), as well a mask (column
%1)
%I - input image. can be a multi-page
%
    %m=size(I,1); n=size(I,2); 
    im_array=cell(size(coord,1),3);
    
    for z=1:size(coord,1)   %create lines separating zones
        curr_coord=coord{z,1};
        im_array{z,1} = roipoly(I, curr_coord(:,1)', curr_coord(:,2)');

        for c=1:size(I,3)
           im_array{z,3} = cat(3, im_array{z,3}, I(:,:,c).*uint16(im_array{z,1}));
        end
        
        for pos=1:2
            if coord{z,2}(1,pos)==0
                coord{z,2}(1,pos)=1;
            end
        end
                
        im_array{z,2} = im_array{z,3}((coord{z,2}(1,1)):(coord{z,2}(1,3)),(coord{z,2}(1,2)):(coord{z,2}(1,4)),:);
    end

end

%% drawzones
function [ im_array ] = drawzones( coord, I)
%drawzones draws lines on an image, based on coordinates from coordinates cell
%array (z,1). z - number of lines that were drawn;
%output is an array 
% Column 1: masks (logical images)
% Column 2: masked images and cropped to bounding box
% Column 3: masked uncropped images

%I - input image. can be a multi-page
%%
m=size(I,1); n=size(I,2); 
im=uint16(zeros(m,n));
lines=[];
    for z=1:size(coord,1)   %create lines separating zones
        curr_coord=coord{z,1};
        [~,order]=sort(curr_coord(:,1),1);
        curr_coord=curr_coord(order,:);
        curr_coord(1,1)=1; curr_coord(size(curr_coord,1),1)=n;  %correct the lines, so it begins from the edge of an image
        curr_line=curr_coord(1,:);
        for segm=2:size(curr_coord,1)
            curr_segm=linepoints(curr_coord(segm-1,:),curr_coord(segm,:));
            curr_line=cat(1,curr_line,curr_segm');
        end
        [~,order]=unique(curr_line(:,1));
        lines=cat(2,lines,curr_line(order,2));
    end
    
    %sort the lines so they are arranged
    positions=min(lines,[],1);
    [~,order]=sort(positions,'ascend');
    lines=lines(:,order);

    for col=1:size(lines,1)
        curr_line=zeros(m,1);
        for z=1:(size(lines,2))
            if z<size(lines,2)
                curr_line(lines(col,z):lines(col,z+1),:)=z;
            else
                curr_line(lines(col,z):m,:)=z;
            end
        end
        im(:,col)=curr_line;
    end
    im_array=cell(size(lines,2)+1,2);

    for z=0:size(lines,2)
        blank_im=uint16(zeros(m,n));
        blank_im(im==z)=1;
        im_array{z+1,1}=blank_im;
        blank_im=repmat(blank_im,[1 1 size(I,3)]);
        im_array{z+1,2}=blank_im.*I;
    end
end

%% exclude_str_cell
function indices = exclude_str_cell( str, cell_array )
%str - string to search for
%cell_array - cell array where the seach is done
    indices=[];
    if size(cell_array,1)==1
        if size(cell_array,2)>1
            cell_array=cell_array';
        end
    end

    for r=1:size(cell_array,1)
        if findstr(str,cell_array{r,1})
        else
            indices=cat(1,indices,r);
        end
    end

end

%% findstr_cell
function indices = findstr_cell( str, cell_array )
%str - string to search for
%cell_array - cell array where the seach is done

    indices=[];
    for r=1:size(cell_array,1)
        if strcmpi(str,cell_array{r,1})
            indices=cat(1,indices,r);
        end
    end
end

%% findstr_cell_partial
function indices = findstr_cell_partial( str, cell_array )
%str - string to search for
%cell_array - cell array where the seach is done
    indices=[];
    if size(cell_array,1)==1
        if size(cell_array,2)>1
            cell_array=cell_array';
        end
    end

    for r=1:size(cell_array,1)
        if strfind(cell_array{r,1},str)
            indices=cat(1,indices,r);
        end
    end

end

%% histc_full
function [ bins,below,above ] = histc_full( x,edges )
%hitsc_full calculates hist_c and adds out of range values as additional
%output

    bins=histc( x,edges );
    below=length(find(x<min(edges)));
    above=length(find(x>max(edges)));

end
