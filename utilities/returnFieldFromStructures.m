function out = returnFieldFromStructures(structs,field) 

    N = length(structs);
    
    if ischar(structs{1}.(field)) || iscell(structs{1}.(field))
        
        out = cell(N,1);
        for i=1:N
            a = structs{i}.(field);
            if length(a) > 1
                out{i} = a;
            else
                out{i} = a{1};
            end
        end
        
    else
        
        s = size(structs{1}.(field));
        
        if sum(s > 1) < 2
            L = length(structs{1}.(field));
            out = zeros(N,L);
            for i=1:N;
                a = structs{i}.(field);
                if ~isempty(a)
                    out(i,:) = structs{i}.(field);
                end
            end
        else
            out = cell(N,1);
            for i=1:N
                out{i} = structs{i}.(field);
            end
        end
        
    end