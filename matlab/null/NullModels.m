classdef NullModels < handle
    
    methods (Access = private)
    %private so that you can't instatiate.
        function out = NullModels

        end
    end 
   
    methods(Static)
        
        function matrix = NULL_1(MatrixOrSizeRows, n_cols, ones)

            if(nargin == 1)

                [n_rows n_cols] = size(MatrixOrSizeRows);
                ones = sum(sum(MatrixOrSizeRows>0));

            else 
                 n_rows = MatrixOrSizeRows;
            end

            hasNull = 1;
            while(hasNull)
                tmpr=rand(n_rows*n_cols,1);
                [tmpy, tmps]=sort(tmpr);
                x=zeros(n_rows,n_cols);
                x(tmps(1:ones))=1;
                hasNull = 0;%sum(sum(x,1)==0)+sum(sum(x,2)==0);
            end

            matrix = x;
        end
        
        function matrix = NULL_2(n_rows, n_cols, probRows,probCols)
            
            if(nargin ==1)
                matrix = n_rows;
                [n_rows n_cols] = size(matrix);
                probCols = sum(matrix)/n_rows;
                probRows = sum(matrix,2)/n_cols;
            end
            
            if(n_rows ~= length(probRows))
                error('The row size is different from the assigned probabilities');
            end
            if(n_cols ~= length(probCols))
                error('The column size is different from the assigned probabilities');
            end
            
            ss = size(probRows);
            if(ss(1) == 1); probRows = probRows'; end;
            
            pr = repmat(probRows,1,n_cols);
            pc = repmat(probCols,n_rows,1);
                      
            mat = (pr+pc)/2;
            
            hasNull = 1;
            while(hasNull)
                matrix = rand(n_rows,n_cols) <= mat;
                hasNull = 0;% = sum(sum(matrix,1)==0) + sum(sum(matrix,2)==0);
            end

        end
        
        function matrix = NULL_3ROW(n_rows, n_cols, probRows)
            %Return a random network in wich colums filling according to
            %probCols
            
            if(nargin ==1)
                matrix = n_rows;
                [n_rows n_cols] = size(matrix);
                probRows = sum(matrix,2)/n_cols;
            end
            
            if(n_rows ~= length(probRows))
                error('The row size is different from the assigned probabilities');
            end
            
            hasNull = 1;
            while(hasNull)
                matrix = rand(n_rows, n_cols) < repmat(probRows,1,n_cols);
                hasNull = sum(sum(matrix,1)==0) + sum(sum(matrix,2)==0);
            end            
        end
        
        function matrix = NULL_3COL(n_rows, n_cols, probCols)

            if(nargin ==1)
                matrix = n_rows;
                [n_rows n_cols] = size(matrix);
                probCols = sum(matrix)/n_rows;
            end
            
            if(n_cols ~= length(probCols))
                error('The column size is different from the assigned probabilities');
            end
            
            hasNull = 1;
            while(hasNull)
                matrix = rand(n_rows, n_cols) < repmat(probCols,n_rows,1);
                hasNull = sum(sum(matrix,1)==0) + sum(sum(matrix,2)==0);
            end   
        end
        
        function rmatrices = NULL_MODEL(adjmatrix,model,replicates)
           
            RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
            
            adjmatrix = adjmatrix ~= 0;
            
            n_edges = sum(adjmatrix(:));
            [n_rows n_cols] = size(adjmatrix);
            p_cols = sum(adjmatrix)/n_rows;
            p_rows = sum(adjmatrix,2)/n_cols;
            p = n_edges / (n_rows * n_cols);
            
            nullmodel = 1;
            if(strcmp(func2str(model),'NullModels.NULL_1')); nullmodel = 1; end;
            if(strcmp(func2str(model),'NullModels.NULL_2')); nullmodel = 2; end;
            if(strcmp(func2str(model),'NullModels.NULL_3ROW')); nullmodel = 3; end;
            if(strcmp(func2str(model),'NullModels.NULL_3COL')); nullmodel = 4; end;
            
            rmatrices = cell(1,replicates);
            for i = 1:replicates
                %fprintf('iteration %i\n',i);
                switch nullmodel
                    case 1
                        rmatrices{i} = NullModels.NULL_1(n_rows,n_cols,n_edges);
                    case 2
                        rmatrices{i} = NullModels.NULL_2(n_rows,n_cols,p_rows,p_cols);
                    case 3
                        rmatrices{i} = NullModels.NULL_3ROW(n_rows,n_cols,p_rows);
                    case 4
                        rmatrices{i} = NullModels.NULL_3COL(n_rows,n_cols,p_cols);
                end
            end
            
        end
        
    end
end