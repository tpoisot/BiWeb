classdef Reading
   
    methods(Static)
       
        function bipartite_web = CREATE_FROM_MATRIX_WITH_LABELS(filename)
           
            fid = fopen(filename); 
            if(fid==-1)
                error('The file could not be open. Check the name of the file.');
            end
                
            row_labels = {}; col_labels = {};
            lread = fgetl(fid);
            idx = find(lread == '"');
            
            for i = 1:2:length(idx)-1
                col_labels{(i+1)/2} = lread(idx(i)+1:idx(i+1)-1);
            end
            
            lread = fgetl(fid);
            i = 1;
            matrix = [];
            while(lread~=-1)
                idx = find(lread == '"');
                row_labels{i} = lread(idx(1)+1:idx(2)-1);
                matrix(i,:) = str2num(lread(idx(2)+1:end));
                i = i + 1;
                lread = fgetl(fid);
            end
            
            bipartite_web = Bipartite(matrix);
            bipartite_web.row_labels = row_labels;
            bipartite_web.col_labels = col_labels;
            [paths namefile ext] = fileparts(filename);
            bipartite_web.name = namefile;
        end
        
    end
    
end