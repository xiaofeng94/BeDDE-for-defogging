clear
% --------------------------
BeDDE_root = '.\BeDDE'; % set it to BeDDE root.

method_name = 'fog'; % set it to the folder for our results
eval_method = 'VI'; % 'VI', 'RI'

% --------------------------


% allFolderInfo = dir(image_root);
% cityFolders = cell(0);
% 
% for fold_id = 1:length(allFolderInfo)
%     folder = allFolderInfo(fold_id).name;
%     isDir = allFolderInfo(fold_id).isdir;
%     if isDir && strcmp(folder, '.') == 0 && strcmp(folder, '..') == 0
%         cityFolders = [cityFolders; {folder}];
%     end
% end
cityFolders = [{'beijing'},{'changsha'}, {'chengdu'},{'hangzhou'},{'hefei'}, ...
                  {'hongkong'},{'lanzhou'},{'nanchang'},{'shanghai'},{'shenyang'}, ...
                  {'tianjing'},{'wuhan'},{'chongqing'},{'guangdong'},{'guilin'},...
                  {'hohhot'},{'kunming'},{'nanjing'},{'sian'},{'taiwan'},...
                  {'xining'},{'yinchuan'},{'zhengzhou'}];

cityNum = length(cityFolders);

allScores = [];
allImages = cell(0);
for city_id = 1:cityNum
    cityName = cityFolders{city_id,1};
    
    clearImgFolder = fullfile(BeDDE_root,cityName,'gt');
    maskFolder = fullfile(BeDDE_root,cityName,'mask');
    imageFolder = fullfile(BeDDE_root,cityName,method_name);
    
    clearImgPath = fullfile(clearImgFolder,[cityName,'_clear.png']);
    if exist(clearImgPath, 'file')
        clearImg = imread(clearImgPath);
    else
        disp(['No clear image for [', cityName, ']']);
        continue
    end
    
    imageFileInfo = dir(fullfile(imageFolder, '*.png'));
    disp(['----- ', cityName])
    for img_id = 1:length(imageFileInfo)
        defogImgName = imageFileInfo(img_id).name;
        defogImgPath = fullfile(imageFolder, defogImgName);
        
%         if strcmp(method_name, 'fog')
%             curMaskPath = fullfile(maskFolder, [defogImgName(1:end-4), '_mask.mat']);
%         else
%             curMaskPath = fullfile(maskFolder, [defogImgName(1:end-length(method_name)-4), 'mask.mat']);
%         end 
        imgNameStrs = split(defogImgName(1:end-4), '_');
        curMaskPath = fullfile(maskFolder, sprintf('%s_%s_mask.mat', imgNameStrs{1}, imgNameStrs{2}));

        defogImg = imread(defogImgPath);
        
        if any(size(defogImg) ~= size(clearImg))
            defogImg = imresize(defogImg, [size(clearImg,1), size(clearImg,2)]);      
        end
        
        curMask = load(curMaskPath);
        curMask = curMask.mask;
        
        curScore = getIQAScoreBy(clearImg, defogImg, curMask, eval_method);
 
        allScores = [allScores; curScore];
        allImages = [allImages; {defogImgName}];
        
        disp([defogImgName, ': ', num2str(curScore)]);
    end
end

totalImgNum = length(allScores);
disp(['test image count: ', num2str(totalImgNum)]);
disp(['score std error: ', num2str(sqrt(var(allScores)))]);
disp(['average score: ', num2str(sum(allScores)/totalImgNum)]);

% save evaluation results to <image_root>/statistics/<defogging_methond>_<eval_method>_eval.mat
result_save_root = fullfile(BeDDE_root, './statistics');
if ~exist(result_save_root, 'dir')
    mkdir(result_save_root);
end

eval_result.scores = allScores;
eval_result.files = allImages;
save(fullfile(result_save_root, [method_name,'_',eval_method,'_eval.mat']), 'eval_result');
