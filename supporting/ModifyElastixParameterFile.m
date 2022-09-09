function ModifyElastixParameterFile(FileName,ParamName,ParamVal, varargin)
%Modify the parameters in the elastix parameter file
    
    %% Inputs

    %FileName - Elastix Parameter File Name
    %ParamName - The parameter to be modified

    %% Read Elastix File
    
    fid = fopen(FileName); 
    if fid > 0
        try
            A = fscanf(fid,'%c');
            fclose(fid);
        catch
            fclose(fid);
            error(['Problem reading file',FileName])
        end        
    else
        error(['Elastix File ',FileName,' could not be opened'])
    end
      
    
    %% Write parameter to Elastix File
    
    BracketOpen = strfind(A,'(');
    BracketClose = strfind(A,')');      
    
    for i = 1:numel(BracketOpen)
        strParam = A(BracketOpen(i)+1:BracketClose(i)-1);
        C = strsplit(strParam);
        if strcmp(C{1},ParamName)
            C{2} = ParamVal;
            A2 = [A(1:BracketOpen(i)),C{1},' ',C{2},A(BracketClose(i):end)];
        end
    end
    
    %% Rewrite elastix Parameter File
    fid = fopen(FileName,'w'); 
    try
        fprintf(fid,'%c',A2);
    catch
        disp(['Problem writing to file',FileName])
    end
    fclose(fid);
end