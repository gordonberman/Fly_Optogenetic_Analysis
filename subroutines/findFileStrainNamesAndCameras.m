function fileData = findFileStrainNamesAndCameras(files,led_files)

    N = length(files);
    L = length(led_files);
    fileData.N = N;
    fileData.files = files;
    fileData.led_files = led_files;
    
    cameraNumbers = zeros(N,1);
    strains = cell(N,1);
    filming_sessions = cell(N,1);
    isTroy = false(N,1);
    
    for i=1:N
        
        [~,name,~] = fileparts(files{i});
        
        idx = strfind(name,'Cam');
        idx2 = find(name == '_');
        q = idx2(find(idx2>idx(1),1,'first'));
        
        strains{i} = name(1:idx2(1)-1);
        if strcmp(strains{i},'troy')
            strains{i} = name(1:idx2(2)-1);
            isTroy(i) = true;
        end
        if strcmp(strains{i},'trial2')
            strains{i} = name(idx2(1)+1:idx2(2)-1);
        end
        
        cameraNumbers(i) = str2double(name((idx(1)+3):(q-1)));
        filming_sessions{i} = name(1:idx(1)-2);
        
    end
    
    
    fileData.strains = strains;
    fileData.cameraNumbers = cameraNumbers;
    fileData.filming_sessions = filming_sessions;
    
    [fileData.strainNames,~,fileData.strainNumbers] = unique(strains);
    [fileData.filming_session_names,~,fileData.filming_session_numbers] = unique(filming_sessions);
    
    fileData.isControl = cameraNumbers < 7;
    fileData.isControl(isTroy) = false;
        
    led_filming_sessions = cell(L,1);
    led_filming_session_numbers = zeros(L,1);
    for i=1:L
        
        [~,name,~] = fileparts(led_files{i});
        idx = strfind(name,'Frames');
        led_filming_sessions{i} = name(1:idx(1)-2);
                
        for j=1:length(fileData.filming_session_names)
            if strcmp(led_filming_sessions{i},fileData.filming_session_names{j})
                led_filming_session_numbers(i) = j;
                break;
            end
        end
        
    end
    fileData.led_filming_sessions = led_filming_sessions;
    fileData.led_filming_session_numbers = led_filming_session_numbers;
    

    
    
    
    
    
    
    