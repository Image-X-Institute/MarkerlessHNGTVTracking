function Projs2Pics(ProjFile,OutputDir)
%Convert a his file projection into a jpg

    close('all')

    %[projHeader,projImg] = MhaRead('TurkeyData\img_1.3.46.423632.135285.1575434474.15\Proj.mha');

    %CTDir = 'Z:\2RESEARCH\2_ProjectData\RemoveTheMask\CTData\HNSCC\HNSCC-01-0001\12-05-1998-RT SIMULATION-43582\Recon\DRR.mha';
    if exist(ProjFile,'file')
        [ctHeader,ctImg] = MhaRead(ProjFile);
    else
        [filedir,name,ext] = fileparts(ProjFile);
        if exist(fullfile(filedir,[name,'1',ext]),'file')
            [ctHeader,ctImg] = MhaRead(fullfile(filedir,[name,'1',ext]));            
        elseif exist(fullfile(filedir,[name,'_1',ext]),'file')
            [ctHeader,ctImg] = MhaRead(fullfile(filedir,[name,'_1',ext])); 
        end
    end
    % figure()
    % ax1 = gca;
    % fig1 = gcf;

    figure()
    ax2 = gca;
    fig2 = gcf;

    %baseFile = 'Z:\2RESEARCH\2_ProjectData\RemoveTheMask\CTData\HNSCC\HNSCC-01-0001\12-05-1998-RT SIMULATION-43582\Recon\ProjPics';

    if ~exist(OutputDir,'dir')
        mkdir(OutputDir)
    end

    %TurkeyProjFile = 'TurkeyData/ProjPic';

    %CTProjFile = 'TurkeyData/ProjPic';

    for i = 1:ctHeader.Dimensions(3)
    %for i = 361:ctHeader.Dimensions(3)
        %if mod(i,3) == 0

            %imshow(projImg(:,:,i),[],'parent',ax1,'border','tight')

            %saveas(fig1,num2str(i,[TurkeyProjFile,'/UnWarpedTurkey%0.3d.jpg']))
            %saveas(fig1,fullfile(baseFile,num2str(i,'UnWarpedTurkey%0.3d.jpg')))

            imshow(ctImg(:,:,i),[],'parent',ax2,'border','tight')

            ctImLoop = ctImg(:,:,i);
            if median(ctImg(:,:,1)) > 0
                ctImLoop(ctImLoop < 0) = 0;
            else
                ctImLoop = ctImLoop - min(ctImLoop(:));
            end
            %saveas(fig2,num2str(i,[CTProjFile,'/CTTurkey%0.3d.jpg']))    
            %saveas(fig2,fullfile(OutputDir,num2str(i,'Proj%0.3d.jpg')))
            %saveas(fig2,fullfile(OutputDir,num2str(i,'DRR_%04d.jpg')))
            imLoop = mat2gray(ctImLoop);
        %end
        %imwrite(projImg(:,:,i),num2str(i,[TurkeyProjFile,'/WarpedTurkey%3d.jpg']))
        imwrite(imLoop,fullfile(OutputDir,num2str(i,'DRR_%04d.jp2')), 'Mode', 'lossless', 'CompressionRatio', 1)
    end
    
    clear ctImg
    
    if ~exist(ProjFile,'file')
        if exist(fullfile(filedir,[name,'2',ext]),'file')
            [ctHeader2,ctImg2] = MhaRead(fullfile(filedir,[name,'2',ext]));
        elseif exist(fullfile(filedir,[name,'_2',ext]),'file')
            [ctHeader2,ctImg2] = MhaRead(fullfile(filedir,[name,'_2',ext]));
        end
        for i = 1:ctHeader2.Dimensions(3)
        %for i = 361:ctHeader.Dimensions(3)
            %if mod(i,3) == 0

                %imshow(projImg(:,:,i),[],'parent',ax1,'border','tight')

                %saveas(fig1,num2str(i,[TurkeyProjFile,'/UnWarpedTurkey%0.3d.jpg']))
                %saveas(fig1,fullfile(baseFile,num2str(i,'UnWarpedTurkey%0.3d.jpg')))

                imshow(ctImg2(:,:,i),[],'parent',ax2,'border','tight')

                ctImLoop = ctImg2(:,:,i);
                if median(ctImg2(:,:,1)) > 0
                    ctImLoop(ctImLoop < 0) = 0;
                else
                    ctImLoop = ctImLoop - min(ctImLoop(:));
                end               
                imLoop = mat2gray(ctImLoop);
                %saveas(fig2,num2str(i,[CTProjFile,'/CTTurkey%0.3d.jpg']))    
                %saveas(fig2,fullfile(OutputDir,num2str(i,'Proj%0.3d.jpg')))
                %saveas(fig2,fullfile(OutputDir,num2str(i,'DRR_%04d.jpg')))

            %end
            %imwrite(projImg(:,:,i),num2str(i,[TurkeyProjFile,'/WarpedTurkey%3d.jpg']))
            imwrite(imLoop,fullfile(OutputDir,num2str(i+ctHeader.Dimensions(3),'DRR_%04d.jp2')), 'Mode', 'lossless', 'CompressionRatio', 1)
        end    
    end
end