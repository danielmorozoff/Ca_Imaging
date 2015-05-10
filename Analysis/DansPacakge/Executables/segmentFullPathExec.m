function segmentFullPathExec(dirName)
animalName= 's10';
date = '2015_03_26';
number_runs = 1; % will be variable by exp

fName = strcat(animalName,'_',date,'_main');
imaging_planes = 4;
number_trials = 16; % will be variable by run

fullExperimentVar = {'runNumber';'planeNumber';'trialNumber';'cellData'};

finalOut{1,1} = 'animalName';
finalOut{1,2} = animalName;
finalOut{2,1} = 'data';

runPlaneStack = {}; %This is a runsxplane stack with i,j entry being a trialStack
for (r =1:number_runs)
    fprintf('New run %d\n',r);
    for(i=1:imaging_planes)
        fprintf('New plane %d\n',i);
        cnt = 0;
        trialStack = {};
        centers = [];
        radii = [];
        sI =[];
        for(j= 1:number_trials)
            
            fprintf('New trial %d\n',j);
            if     j<10
                j_display = strcat('00',int2str(j));
            elseif j<100
                j_display = strcat('0',int2str(j));
            else
                j_display = int2str(j);
            end
            imageName = strcat('Image_Registration_4_',fName,'_',j_display,'.tif');
            path = strcat(dirName,'/run',int2str(r),'/fov_','0100',int2str(i),'/fluo_batch_out/');
            
            segmentPath  =strcat(path,'session_maxproj_chan_01');
            imagePath = strcat(path,imageName);
            
            fprintf('Segment Image:%s \n',segmentPath)
            fprintf('Stack to process:%s \n',imagePath)
            [cell_signals,centers,radii,sI] = segment(segmentPath,imagePath,centers,radii,sI,cnt);
            cnt = cnt+1;
        %     plot_traces();
        
        
        fullExperimentVar{1,2} = r;
        fullExperimentVar{2,2} = i;
        fullExperimentVar{3,2} = j;
        fullExperimentVar{4,2} = cell_signals;
        
        trialStack{j,1} =  fullExperimentVar;
        end
        runPlaneStack{r,i} = trialStack;
    end
end
finalOut{2,2} = runPlaneStack;

outputFileName = strcat('segmentationData_',animalName,number_runs,'_',date);
save(outputFileName,'finalOut');
end
