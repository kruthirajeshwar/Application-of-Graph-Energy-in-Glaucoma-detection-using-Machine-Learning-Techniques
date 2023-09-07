clc;
clear all;
close all;

all_files = dir('D:\Kruthi\Research\Retinal Scanning\DRIVE_Augumented');
sorted_files = natsortfiles(all_files);

% Label_file = readtable('corrected_dataset.csv','NumHeaderLines',1); % skips the first three rows of data
%Labels = zeros(height(Label_file));
%Data_array = zeros(height(Label_file),12);
label = 0;
for image = 2:length(sorted_files)
    imagepath = fullfile(sorted_files(image).folder,sorted_files(image).name);

    %folderpath = 'D:\Kruthi\Research\Retinal Scanning\DRIVE';
    %addpath(folderpath)
    %addpath('C:\Users\Kamya\Documents\Research\Retinal scanning\Biometric image database')
    %addpath 'C:\Users\Kamya\Documents\Research\Retinal scanning'
    %I = imread('C:\Users\Kamya\Documents\Research\Retinal scanning\Biometric image database\01.tif');
    %testfile = '\Train\Glaucoma_Negative\002.jpg';
    %datapath = fullfile(folderpath,testfile)
    %I = imread(datapath);
    I = imread(imagepath);
    Img = im2gray(I);
    %figure, imshow(I);
    
    % Detecting the mask
    eyeMask=imbinarize(Img, graythresh(Img));
    eyeMask=imfill(eyeMask, 'holes');
    eyeMask=bwareaopen(eyeMask,100);
    
    %Applying the CLAHE Algorithm
    J = adapthisteq(Img,'numTiles',[64 64],'nBins',256);
    %figure, imshow(J)

    %Removing noise
    J(~eyeMask)=0;
    
    %Smoothening using an Gaussian filter
    JF = imgaussfilt(J,20);
    %h = fspecial('average', [9 9]);
    %JF = imfilter(J, h);
    %figure, imshow(JF)
    
    %Extracting the contrasted segments
    Z = imsubtract(JF, J);
    %figure, imshow(Z)
    % Otsu Thresholding
    level = graythresh(Z);
    seg_value=0.02;
    BW = imbinarize(Z, level);
    %figure, imshow(BW)
    
    %Morphological opening
    BW2 = bwareaopen(BW, 200);
    
    BW2(~imerode(eyeMask, strel('disk', 6)))=0;
    %figure, imshow(BW2);
    
    BW2=imdilate(BW2,strel('disk',30));
    %figure, imshow(BW2)
    
    se = strel('disk',1);
    %BW2=medfilt2(BW2,[2 2]);
    BW2=imopen(BW2,se);
    %figure, imshow(BW2)
    
    % BW2=bwmorph(BW2,'skel',Inf);
    % skelImg=bwmorph(BW2,'spur',1);
    skelImg = bwskel(BW2);
    %figure, imshow(skelImg)
    
    branchImg=detectBranchpoints(skelImg);
    [row column]=find(branchImg);
    branchPts=[row column];
    
    i=1;
    brancharray=zeros(height(skelImg),width(skelImg));
    crossoverarray=zeros(height(skelImg),width(skelImg));
    
    for r=2:height(skelImg)-1 
        for c=2:width(skelImg)-1 
            if(branchImg(r,c)==1)
                     neighbor(1) = skelImg(r-1,c-1); 
                     %Upper_left. r = row, c = column. 
                     neighbor(2) = skelImg(r-1,c); 
                     %Upper_middle. r = row, c = column. 
                     neighbor(3) = skelImg(r-1,c+1); 
                     %Upper_right. r = row, c = column. 
                     neighbor(4) = skelImg(r,c-1); 
                     neighbor(5) = skelImg(r,c);
                     %left. r = row, c = column. 
                     neighbor(6) = skelImg(r,c+1); 
                     %right. r = row, c = column. 
                     neighbor(7) = skelImg(r+1,c+1); 
                     %Lowerleft. r = row, c = column. 
                     neighbor(8) = skelImg(r+1,c); 
                     %lower_middle. r = row, c = column. 
                     neighbor(9) = skelImg(r+1,c-1); 
                     %Lower_left. r = row, c = column. 
                     checkpt(i)=sum(neighbor(:) == 1);
                     if(checkpt(i)==4)
                         brancharray(r,c)=1;
                     else
                         crossoverarray(r,c)=1;
                     end
                     i=i+1;
            end
        end
    end
    
    %% Further segregation of cross and branch points
    number_of_crossandbPts=sum(sum(crossoverarray(:)==1),sum(brancharray(:)==1)); 
    nodeIDbc=cell(1,number_of_crossandbPts);
    node=0;
    
    for r=1:height(skelImg) 
        for c=1:width(skelImg)
            if(brancharray(r,c)==1)
                node=node+1;
                nodeIDbc{node}=[node,r,c,1];
            elseif(crossoverarray(r,c)==1)
                node=node+1;
                nodeIDbc{node}=[node,r,c,2]; 
            end
        end
    end
    
    ncrossover=sum(sum(crossoverarray(:)==1),sum(brancharray(:)==1));
    nodecrossover=cell(1,ncrossover);
    ccounter=1;
    i=1;
    j=0;
    % disp(length(nodeIDbc))
    while(i<=length(nodeIDbc))
        if(nodeIDbc{i}(4)==1)
            x1=nodeIDbc{i}(2);
            y1=nodeIDbc{i}(3);
            j=1;
            while(j<=length(nodeIDbc))
                if(nodeIDbc{j}(4)==1 && i~=j)
                    x2=nodeIDbc{j}(2);
                    y2=nodeIDbc{j}(3);
                    dist=sqrt(((x2-x1)^2)+((y2-y1)^2));
                    if(dist<13)
                        xcross=round((x1+x2)/2);
                        ycross=round((y1+y2)/2);
                        skelImg(xcross,ycross)=1;
                        nodecrossover{ccounter}=[ccounter, xcross, ycross];
                        ccounter=ccounter+1;
                        if (i>j)
                            nodeIDbc(i)=[];
                            nodeIDbc(j)=[];
                        else
                            nodeIDbc(i)=[];
                            nodeIDbc(j-1)=[];
                        end
    
                    end
                end
                j=j+1;
    %             disp(j)
            end
        end
        i=i+1;
    end
    
    
    for r=1:height(skelImg) 
        for c=1:width(skelImg)
            if(crossoverarray(r,c)==1)
                ccounter=ccounter+1;
                nodecrossover{ccounter}=[ccounter,r,c];
            end
        end
    end
    
     nodecrossover=nodecrossover(~cellfun('isempty',nodecrossover));
    
    for i=1:length(nodeIDbc)
        if(nodeIDbc{i}(1)~=i)
            nodeIDbc{i}(1)=i;
            r1(i,1)=nodeIDbc{i}(2);
            c1(i,1)=nodeIDbc{i}(3);
        else
            r1(i,1)=nodeIDbc{i}(2);
            c1(i,1)=nodeIDbc{i}(3);
        end
    end
    
    rc1=[r1,c1];
    
    for i=1:length(nodecrossover)
        if(nodecrossover{i}(1)~=i)
            nodecrossover{i}(1)=i;
            r2(i,1)=nodecrossover{i}(2);
            c2(i,1)=nodecrossover{i}(3);
        else
            r2(i,1)=nodecrossover{i}(2);
            c2(i,1)=nodecrossover{i}(3);
        end
    end
    rc2=[r2,c2];
    %nodeIDbc-all branch points
    %nodecrossover- crossover points
    
    [labeledImage, numberOfRegions] = bwlabel(skelImg);
    G=graph;    %crossover
    H=graph;    %branch points
    
    %branch graph
    for i=1:length(rc1)
        for j=i+1:length(rc1)
            [img_height, img_width] = size(labeledImage);
            if rc1(i,1) > img_height || rc1(i,2) > img_width || rc1(j,1) > img_height || rc1(j,2) > img_width
                break
            else
                if(labeledImage(rc1(i,1), rc1(i,2)) == labeledImage(rc1(j,1), rc1(j,2)) && j~=(length(rc1)-1) && i~=(length(rc1)-1))
                    for m=1:length(nodeIDbc)
                        if(nodeIDbc{m}(2)==rc1(i,1)&&nodeIDbc{m}(3)==rc1(i,2))
                           h1=nodeIDbc{m}(1);
                           x1=nodeIDbc{m}(2);
                           y1=nodeIDbc{m}(3);
                        elseif(nodeIDbc{m}(2)==rc1(j,1)&&nodeIDbc{m}(3)==rc1(j,2))
                            h2=nodeIDbc{m}(1);
                            x2=nodeIDbc{m}(2);
                            y2=nodeIDbc{m}(3);
                        end
                    end
                    Weight=sqrt(((x2-x1)^2)+((y2-y1)^2));
                    H=addedge(H,h1,h2,Weight);
                end
            end
        end

        H=simplify(H);

        BA=adjacency(H,'weighted');

        
        for i=1:length(rc2)
            for j=i+1:length(rc2)
                [img_height, img_width] = size(labeledImage);
                if rc2(i,1) > img_height || rc2(i,2) > img_width || rc2(j,1) > img_height || rc2(j,2) > img_width
                    break
                else
                    if(labeledImage(rc2(i,1), rc2(i,2)) == labeledImage(rc2(j,1), rc2(j,2)))
                        for m=1:length(nodecrossover)
                            if(nodecrossover{m}(2)==rc2(i,1)&&nodecrossover{m}(3)==rc2(i,2))
                                h1=nodecrossover{m}(1);
                                x1=nodecrossover{m}(2);
                                y1=nodecrossover{m}(3);
                            elseif(nodecrossover{m}(2)==rc2(j,1)&&nodecrossover{m}(3)==rc2(j,2))
                                h2=nodecrossover{m}(1);
                                x2=nodecrossover{m}(2);
                                y2=nodecrossover{m}(3);
                            end
                        end
                        Weight=sqrt(((x2-x1)^2)+((y2-y1)^2));
                        G=addedge(G,h1,h2,Weight);
                    end
                end
            end
        end
    end

    G=simplify(G);
    CA=full(adjacency(G,'weighted'));
    
    CA( ~any(CA,2), : ) = [];  %rows
    CA( :, ~any(CA,1) ) = [];  %columns
    
    BA( ~any(BA,2), : ) = [];  %rows
    BA( :, ~any(BA,1) ) = [];  %columns
    
    max_value = max(CA,[],'all');
    CA = CA./max_value;
    CA = 1-CA;
    fuzzy_graph_CA = graph(CA);
    fuzzy_graph_CA = simplify(fuzzy_graph_CA);
    fuzzy_CA = adjacency(fuzzy_graph_CA,'weighted');

    max_value = max(BA,[],'all');
    BA = BA./max_value;
    BA = 1-BA;
    fuzzy_graph_BA = graph(BA);
    fuzzy_graph_BA = simplify(fuzzy_graph_BA);
    fuzzy_BA = adjacency(fuzzy_graph_BA,'weighted');

    vector_fea_Branch = energyfeatures(fuzzy_BA);
    
    vector_fea_Cross = energyfeatures(fuzzy_CA);

%     Labels = Label_file.Var2(image-2);
    if (rem((image-3),9)==0)
        label = label+1;
    end

    feature_vector = [vector_fea_Branch, vector_fea_Cross, label];

    writematrix(feature_vector,'Recognition_dataset.xlsx','WriteMode','append');

    disp(image)

end

