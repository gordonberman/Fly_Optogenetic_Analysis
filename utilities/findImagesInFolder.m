function images = findImagesInFolder(folder,fileType,frontConstraint)

    if ~strcmp(folder(end),slashVal())
        folder = [folder slashVal()];
    end

    if nargin > 2 && ~isempty(frontConstraint)
        a = dir([folder frontConstraint '*' fileType ]);
    else
       a = dir([folder '*' fileType]);
    end
    
    images = {a(:).name}';
    images = strcat(folder,images);
    
    if min(size(images)) == 0
        images = {};
    end
    
