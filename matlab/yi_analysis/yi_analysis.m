function yi_analysis(data_struct)


%loop through mice
for c=1:1:numel(data_struct)
    data=data_struct(c);
    roidir=data.roidir;
    imsize=data.params.imsize;
    %load ROIs, same for all conditions
    fnames=dir(roidir);
    count=1;
    for i=1:1:numel(fnames)
        if ~isempty(strfind(fnames(i).name,'roi'))
            data.data.roi{count}.name=fnames(i).name;
            
            %read roi
            roidata=ReadImageJROI([roidir filesep data.data.roi{count}.name]);
            data.data.roi{count}.corners=roidata.mnCoordinates;
            data.data.roi{count}.mask=poly2mask(data.data.roi{count}.corners(:,1),data.data.roi{count}.corners(:,2),imsize(1),imsize(2));
            data.data.roi{count}.data=[];
            count=count+1;
        end
    end
    
    %loop through trial types / stimulus conditions
    %figure for results
    %nconditions=data.params.trialTypes;
    nconditions=numel(data.datadir);
    for v=1:1:nconditions
        nTrials=data.params.nTrials;
        framesPerTrial=data.params.framesPerTrial;
        datadir=[data.maindatadir filesep data.datadir{v}];
        
        %load on data frame by frame
        fnames=dir([datadir filesep '*.tiff']);
        
        if numel(fnames) < framesPerTrial
            continue
        end

        %loop through tifs on this trial
        %correction for additional frames at beginning to ignore
        startframe=mod(numel(fnames),nTrials*framesPerTrial)+1;
        count=1;
        im=zeros(imsize(1),imsize(2),numel(startframe:numel(fnames)));

        %image loading loop
        for i=startframe:numel(fnames)
            if ~isempty(strfind(fnames(i).name,'tiff'))
                im(:,:,count)=imread([fnames(i).folder filesep fnames(i).name]);
            end
            
            %extract raw mean f for frame
            for k=1:1:numel(data.data.roi)
                data.data.roi{k}.frame(count)=count;
                data.data.roi{k}.meanF(count)=sum(sum(data.data.roi{k}.mask.*double(im(:,:,count))))./(sum(sum(data.data.roi{k}.mask)));
            end
            count=count+1;
        end

        %optional kalman filter on raw stack
        if data.params.kalman==1
            try
                im=Kalman_Stack_Filter(im);
            catch
                %keyboard
            end
        end

        %%optional post filtering stages loop
        for i=1:1:size(im,3)
            if data.params.bksub_each==1
                im(:,:,i)=im(:,:,i)-rollingImBkgrnd(im(:,:,i),50);
            end
            if data.params.norm==1
                tt=im(:,:,i);
                im(:,:,i)=im(:,:,i)./mean(tt(:));
            end
        end
   
        %data extraction loop for ROI on filtered data
        for i=1:1:size(im,3)
            for k=1:1:numel(data.data.roi)
                data.data.roi{k}.meanF_filter(i)=sum(sum(data.data.roi{k}.mask.*double(im(:,:,i))))./(sum(sum(data.data.roi{k}.mask)));
            end
        end
        
        %generate mean and max  image across trials
        im2=reshape(im,[imsize(1) imsize(2) framesPerTrial nTrials]);

        %max for renormaliziopn
        maxx=max(im2(:));
       
        %pre image is average of first 10 frames
        preim=im2(:,:,1:data.params.stimFrame,:);
        preim_mean=mean(preim,4);
        preim_mean16=preim_mean.*((2^16)./maxx);
        preim_mean16=mean(preim_mean,3);
        
        %post images,  each frame after, we take all of them
        for bb=1:1:(data.params.stimFrame+1:size(im2,3))
            postim{bb}=im2(:,:,data.params.stimFrame+1:data.params.stimFrame+bb,:);
            postim_mean{bb}=mean(postim{bb},4);
            
            %compute df images
            df{bb}=mean(postim_mean{bb},3)-mean(preim_mean,3);
            dfup{bb}=df{bb};
            dfdn{bb}=df{bb};
            dfup{bb}(dfup{bb}<0)=0;
            dfdn{bb}(dfup{bb}>0)=0;
            dfdn{bb}=dfdn{bb}*-1;
            
            %rescale to 16-bit the mean
            dfdn16{bb}=dfdn{bb}.*((2^16)./maxx);
            dfup16{bb}=dfup{bb}.*((2^16)./maxx);
            postim_mean16{bb}=mean(postim_mean{bb},3);
            postim_mean16{bb}=postim_mean16{bb}.*((2^16)./maxx);
        end
        
        %create a max dfup, and dfdwn image for mask a la Yi
        dfdn16_max=max(reshape(cell2mat(dfdn16),[imsize numel(dfdn16)]),[],3);      
        dfup16_max=max(reshape(cell2mat(dfup16),[imsize numel(dfup16)]),[],3);   
        
        %create threshold image, and masks
        percentileThresh=0.25;
        dfdn16_max_thrsh=dfdn16_max; dfdn16_max_thrsh(find(dfdn16_max_thrsh<(max(dfdn16_max_thrsh(:)).*percentileThresh)))=0; dfdn16_max_thrsh=medfilt2(dfdn16_max_thrsh,[3 3]);
        dfup16_max_thrsh=dfup16_max; dfup16_max_thrsh(find(dfup16_max_thrsh<(max(dfup16_max_thrsh(:)).*percentileThresh)))=0; dfup16_max_thrsh=medfilt2(dfup16_max_thrsh,[3 3]);
        dfdn16_max_msk=dfdn16_max_thrsh;dfdn16_max_msk(find(dfdn16_max_msk>0))=1;
        dfup16_max_msk=dfup16_max_thrsh;dfup16_max_msk(find(dfup16_max_msk>0))=1;
        
        % df calculation for ROIs - use data from kalman filtered stacks,
        % meanF_filter
        figure(1000);imagesc(medfilt2(dfdn16_max.*-1+dfup16_max,[3 3])),colormap('gray'),colorbar; %for plotting all rois
        winLength=framesPerTrial*2;
        quartileFrac=0.2;
        for k=1:1:numel(data.data.roi)
            fZero=medfilt1(data.data.roi{k}.meanF_filter,winLength);
            data.data.roi{k}.fZero=fZero;
            data.data.roi{k}.df_F=(data.data.roi{k}.meanF_filter-data.data.roi{k}.fZero)./data.data.roi{k}.meanF_filter;
        end
        
        %df and dFF calculation for up and down masks
        tmp=reshape(sum(sum(repmat(dfdn16_max_msk,[1 1 size(im,3)]).*im)),size(im,3),1); 
        dndF=tmp./sum(dfdn16_max_msk(:));
        fZero=medfilt1(dndF,winLength);
        dndF_F=(dndF-fZero)./dndF;
        
        tmp=reshape(sum(sum(repmat(dfup16_max_msk,[1 1 size(im,3)]).*im)),size(im,3),1); 
        updF=tmp./sum(dfup16_max_msk(:));    
        fZero=medfilt1(updF,winLength);
        updF_F=(updF-fZero)./updF;        
        
        %plot for the up and down mask
        %%create figure name
        figure(1001+v)
        %mouse name
        mn=strrep(data.name,'_',' ');
        %condition name
        cn=[num2str(data.MPa{v}) 'MPa, ' num2str(data.DC{v}) 'DC'];
        hand=gcf; hand.Name=[mn ' ' cn];
        
        %plot down responses, and their mean
        subplot(2,4,1),plot(reshape(dndF_F,framesPerTrial,nTrials),'k'), axis([2 60 -0.4 0.6])
        hold on
        subplot(2,4,1),plot(mean(reshape(dndF_F,framesPerTrial,nTrials),2),'r','Linewidth',3)
        hold off
        %down image
        subplot(2,4,2),imagesc(dfdn16_max_thrsh),colormap('gray')
        
        %plot up responses, and their mean
        subplot(2,4,5),plot(reshape(updF_F,framesPerTrial,nTrials),'k'), axis([2 60 -0.4 0.6])
        hold on
        subplot(2,4,5),plot(mean(reshape(updF_F,framesPerTrial,nTrials),2),'g','Linewidth',3)
        hold off
        %up image
        subplot(2,4,6),imagesc(dfup16_max_thrsh),colormap('gray');
        
        %plot combined up and down image
        subplot(2,4,3),imagesc(medfilt2(dfdn16_max,[3 3])),colormap('gray'); %for plotting all rois
        subplot(2,4,7),imagesc(medfilt2(dfup16_max,[3 3])),colormap('gray'); %for plotting all rois
        
        %overlay rois that are overlapping on this
        for k=1:1:numel(data.data.roi)
            %test for rois that are part of the up/down masks
            upids{v}=[];
            dnids{v}=[];
            % perc overlap betewrrn roi and mask
            % if greather frac_thresh, its overlapped
            frac_thresh=0.05;
            frac=sum(sum((dfup16_max_msk+data.data.roi{k}.mask)>=2))./sum(data.data.roi{k}.mask(:));
            if frac>=frac_thresh
                upids{v}=[upids{v} k];
            end
            frac=sum(sum((dfdn16_max_msk+data.data.roi{k}.mask)>=2))./sum(data.data.roi{k}.mask(:));
            if frac>=frac_thresh
                dnids{v}=[dnids{v} k];
            end
            if any(dnids{v}==k)
                subplot(2,4,3), hold on; line(data.data.roi{k}.corners(:,1),data.data.roi{k}.corners(:,2),'color','r','lineWidth',1)
                %and line data from roi
                subplot(2,4,4),hold on; h=plot(mean(reshape(data.data.roi{k}.df_F,framesPerTrial,nTrials),2),'r','lineWidth',1.5); axis([2 60 -0.4 0.6])
                uistack(h,'top')
                hold off
            elseif any(upids{v}==k)
                subplot(2,4,7), hold on; line(data.data.roi{k}.corners(:,1),data.data.roi{k}.corners(:,2),'color','g','lineWidth',1)
                %and line data from roi
                subplot(2,4,8),hold on; h=plot(mean(reshape(data.data.roi{k}.df_F,framesPerTrial,nTrials),2),'g','lineWidth',1.5); axis([2 60 -0.4 0.6])
                uistack(h,'top')
                hold off;
            else
                subplot(2,4,3), hold on; line(data.data.roi{k}.corners(:,1),data.data.roi{k}.corners(:,2),'color','cyan','lineWidth',0.5)
                hold off
                subplot(2,4,7), hold on; line(data.data.roi{k}.corners(:,1),data.data.roi{k}.corners(:,2),'color','cyan','lineWidth',0.5)
                hold off
                subplot(2,4,4),hold on; plot(mean(reshape(data.data.roi{k}.df_F,framesPerTrial,nTrials),2),'color',[0.8 0.8 0.8],'lineWidth',0.5), axis([2 60 -0.4 0.6])
                hold off
                subplot(2,4,8),hold on; plot(mean(reshape(data.data.roi{k}.df_F,framesPerTrial,nTrials),2),'color',[0.8 0.8 0.8],'lineWidth',0.5), axis([2 60 -0.4 0.6])
                hold off
            end
        end
        hold off

        %extract peaks from up and down , and plot
        meantrace=mean(reshape(dndF_F,framesPerTrial,nTrials),2);
        data.mask.downdata_dFF{v}=min(meantrace(10:30));
        
        meantrace=mean(reshape(updF_F,framesPerTrial,nTrials),2);
        data.mask.updata_dFF{v}=max(meantrace(10:30));


        %save figure of data
        savefig(gcf,[data.maindatadir filesep [mn '_' cn] '.fig'])
        %keyboard
    end %end condition/trial type loop
    %save data
    save([data.maindatadir filesep data.name '_data.mat'],'data');
    %keyboard
end %end animal loop