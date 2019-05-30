clear

% --------------------------
exBeDDE_root = '.\exBeDDE'; % path to exBeDDE

% one of {'VI', 'RI', 'FADE', 'VSI', 'FSIM', 'FSIMc', 'GMSD', 'PSNR', 'SSIM'}
IQA_metric = 'VI';

% true for testing on foggy groups; 
% false for testing on defogged groups; 
% if you want to test our RI, set it to false.
fog_group_only = true;

% --------------------------

cityNames = [{'beijing'},{'changsha'}, {'chengdu'},{'hangzhou'},{'hefei'}, ...
              {'hongkong'},{'lanzhou'},{'nanchang'},{'shanghai'},{'shenyang'}, ...
              {'tianjing'},{'wuhan'}];

obj_root = fullfile(exBeDDE_root, 'statistics');

exclude_subfolders = {'.','..','gt','mask','fog'};


objScoreMat = load(fullfile(obj_root, ['fog_',IQA_metric,'_eval.mat']));
objScoreMat = objScoreMat.eval_result;
objFiles = objScoreMat.files;
objScores = objScoreMat.scores;

seltIQAScores = []; % ordered as subjective scores
seltSubjScores = []; % ordered as subjective scores

cityNum = length(cityNames);

PLCCs = zeros(cityNum,1);SRCCs = zeros(cityNum,1);KRCCs = zeros(cityNum,1);
MAEs = zeros(cityNum,1);RMSs = zeros(cityNum,1);
imgNums = zeros(cityNum,1);

for city_id = 1:cityNum
    curCityName = cityNames{city_id};
    
    if fog_group_only
        subfolders = {'fog'};
    else
        subfoldInfo = dir(fullfile(exBeDDE_root,curCityName));
        subfolders = cell(0);
        for info_id = 1:length(subfoldInfo)
            curSubfoldName = subfoldInfo(info_id).name;
            if subfoldInfo(info_id).isdir && ~any(strcmp(curSubfoldName,exclude_subfolders))
                subfolders = [subfolders; {curSubfoldName}];
            end
        end
    end
    
    subFoldNum = length(subfolders);
    subPLCCs = zeros(subFoldNum,1);subSRCCs = zeros(subFoldNum,1);subKRCCs = zeros(subFoldNum,1);
    subMAEs = zeros(subFoldNum,1);subRMSs = zeros(subFoldNum,1);
    for subf_id = 1:subFoldNum
        curSubFold = subfolders{subf_id};

        curFold = fullfile(exBeDDE_root,curCityName,curSubFold);
        
        if strcmp(curSubFold, 'fog')
            scoreMatPath = fullfile(curFold, [curCityName,'_',curSubFold,'_scores.mat']);
        else
            scoreMatPath = fullfile(curFold, [curSubFold,'_scores.mat']);
        end
        curScoreMat = load(scoreMatPath);
        curScoreMat = curScoreMat.imageScores;

        curImNames = curScoreMat.image_names;
        curImScores = 1.0001 - curScoreMat.scores;

        curImNum = length(curImNames);
        curIQAScores = zeros(curImNum,1);

        for im_id = 1:curImNum
            curImageName = curImNames{im_id};

            if strcmp(curSubFold, 'fog')
                file_id = find(strcmp(curImageName, objFiles) == 1);
                curIQAScores(im_id) = objScores(file_id);
            else
                %for dehazed images
                tempStrs = split(curImageName(1:end-4), '_');
                dehazedMeth = tempStrs{end};
                dezObjScoreMat = load(fullfile(obj_root, [dehazedMeth,'_',IQA_metric,'_eval.mat']));
                dezObjScoreMat = dezObjScoreMat.eval_result;
                dezObjFiles = dezObjScoreMat.files;
                dezObjScores = dezObjScoreMat.scores;

                file_id = find(strcmp(curImageName, dezObjFiles) == 1);
                curIQAScores(im_id) = dezObjScores(file_id);
            end
        end

    %     nlMapFun = 'y~b1*(1/2-1/(1+exp(b2*(x-b3))))+b4*x+b5';
        nlMapFun = @(b,x)b(1).*(1/2-1./(1+exp(b(2).*(x-b(3)))))+b(4).*x+b(5);
        mdl = fitnlm(curIQAScores, curImScores, nlMapFun, [1,1,1,1,1]);
        coefs = mdl.Coefficients{:,1};

        curSRCC = corr(curImScores, curIQAScores, 'type', 'spearman');
        curKRCC = corr(curImScores, curIQAScores, 'type', 'kendall');

        curIQAScores_adj = nlMapFun(coefs, curIQAScores);
        curPLCC = corr(curImScores, curIQAScores_adj, 'type', 'Pearson');

        curMAE = mean(abs(curIQAScores_adj-curImScores));
        curRMS = sqrt( mean((curIQAScores_adj-curImScores).^2) );
    
        
        imgNums(city_id) = imgNums(city_id) + curImNum;
        subPLCCs(subf_id) = curPLCC;
        subSRCCs(subf_id) = curSRCC;
        subKRCCs(subf_id) = curKRCC;
        subMAEs(subf_id) = curMAE;
        subRMSs(subf_id) = curRMS;
    end 
    
    fprintf('- current city: %s, image num: %d) \n',curCityName,imgNums(city_id));
    fprintf('SRCC: ');fprintf('%f ', subSRCCs);fprintf('\n');
    fprintf('KRCC: ');fprintf('%f ', subKRCCs);fprintf('\n');
    fprintf('PLCC: ');fprintf('%f ', subPLCCs);fprintf('\n');
    fprintf('MAE: ');fprintf('%f ', subMAEs);fprintf('\n');
    fprintf('RMS: ');fprintf('%f ', subRMSs);fprintf('\n');
    
    PLCCs(city_id) = mean(subPLCCs);
    SRCCs(city_id) = mean(subSRCCs);
    KRCCs(city_id) = mean(subKRCCs);
    MAEs(city_id) = mean(subMAEs);
    RMSs(city_id) = mean(subRMSs);
end

% avgPLCC = sum(imgNums.*PLCCs/sum(imgNums));
% avgSRCC = sum(imgNums.*SRCCs/sum(imgNums));
% avgKRCC = sum(imgNums.*KRCCs/sum(imgNums));
% avgMAE = sum(imgNums.*MAEs/sum(imgNums));
% avgRMS = sum(imgNums.*RMSs/sum(imgNums));
avgPLCC = mean(PLCCs);
avgSRCC = mean(SRCCs);
avgKRCC = mean(KRCCs);
avgMAE = mean(MAEs);
avgRMS = mean(RMSs);

fprintf('avg SRCC: %f \navg KRCC: %f \navg PLCC: %f \navg MAE: %f \navg RMS: %f\n', ...
                avgSRCC,avgKRCC,avgPLCC, avgMAE,avgRMS)


