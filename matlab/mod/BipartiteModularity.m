classdef BipartiteModularity < handle
   
    properties
       
        matrix               = [];  %Bipartite adjacency matrix
        rr                   = [];  %Red Nodes (rows) Communities matrix. Size = n_rows*CommunityQuantity
        tt                   = [];  %Blue Nodes (columns) Communities. Size = n_cols*CommunityQuantity
        n_rows               = 0;   %Number of rows
        n_cols               = 0;   %Number of columns
        n_edges              = 0;
        bb                   = [];  %Original - Null.
        index_rows           = [];  %Register of the swaps in Rows.
        index_cols           = [];  %Register of the swaps in Cols.
        trials               = 10;
        Qb                   = 0;
        Qr                   = 0;
        N                    = 0;
        row_modules          = [];
        col_modules          = [];
        done                 = 0;
        
    end
    
    methods(Abstract)
        
        obj = Detect(obj);
        
    end
    
    methods
       
        function obj = SortModules(obj)
           
            if(obj.n_rows > obj.n_cols)
                [~,idx] = sort(sum(obj.rr),'descend');
            else
                [~,idx] = sort(sum(obj.tt),'descend');
            end
                        
            obj.index_rows = 1:obj.n_rows;
            obj.index_cols = 1:obj.n_cols;
            
            obj.rr = obj.rr(:,idx);
            obj.tt = obj.tt(:,idx);
            
            sorted_matrix = obj.matrix;
            
            row_global = []; col_global = [];
            for i = 1:obj.N
               
                row_i = find(obj.rr(:,i));
                col_i = find(obj.tt(:,i));
                
                [~,row_loc] = sort(sum(sorted_matrix(row_i,col_i),2),'descend');
                [~,col_loc] = sort(sum(sorted_matrix(row_i,col_i),1),'ascend');
                
                row_global = [row_global; row_i(row_loc)];
                col_global = [col_global; col_i(col_loc)];
                
            end
           
            col_global = flipud(col_global);
            obj.rr = obj.rr(row_global,:);
            obj.tt = obj.tt(col_global,:);
            
            obj.index_rows = row_global;
            obj.index_cols = col_global;
            
            obj.row_modules = obj.row_modules(row_global);
            obj.col_modules = obj.col_modules(col_global);
            
        end
        
    end
    
end