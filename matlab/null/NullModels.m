classdef NullModels < handle
    
    methods (Access = private)
    %private so that you can't instatiate.
        function out = NullModels

        end
    end 
   
    methods(Static)
        
        function matrix = NULL_ER(MatrixOrSizeRows, sizeCols, ones)

            if(nargin == 1)

                [sizeRows sizeCols] = size(MatrixOrSizeRows);
                ones = sum(sum(MatrixOrSizeRows>0));

            else 
                 sizeRows = MatrixOrSizeRows;
            end

            tmpr=rand(sizeRows*sizeCols,1);
            [tmpy, tmps]=sort(tmpr);
            x=zeros(sizeRows,sizeCols);
            x(tmps(1:ones))=1;

            matrix = x;
        end
        
        function matrix = NULL_CE(sizeRows, sizeCols, probRows,probCols)
            
            if(nargin ==1)
                matrix = sizeRows;
                [sizeRows sizeCols] = size(matrix);
                probCols = sum(matrix)/sizeRows;
                probRows = sum(matrix,2)/sizeCols;
            end
            
            if(sizeRows ~= length(probRows))
                error('The row size is different from the assigned probabilities');
            end
            if(sizeCols ~= length(probCols))
                error('The column size is different from the assigned probabilities');
            end
            
            ss = size(probRows);
            if(ss(1) == 1); probRows = probRows'; end;
            
            pr = repmat(probRows,1,sizeCols);
            pc = repmat(probCols,sizeRows,1);
                      
            mat = (pr+pc)/2;
            matrix = rand(sizeRows,sizeCols) <= mat;

        end
        
        function matrix = NULL_ROW(sizeRows, sizeCols, probRows)
            %Return a random network in wich colums filling according to
            %probCols
            
            if(nargin ==1)
                matrix = sizeRows;
                [sizeRows sizeCols] = size(matrix);
                probRows = sum(matrix,2)/sizeCols;
            end
            
            if(sizeRows ~= length(probRows))
                error('The row size is different from the assigned probabilities');
            end
            
            matrix = rand(sizeRows, sizeCols) < repmat(probRows,1,sizeCols);
                        
        end
        
        function matrix = NULL_COL(sizeRows, sizeCols, probCols)

            if(nargin ==1)
                matrix = sizeRows;
                [sizeRows sizeCols] = size(matrix);
                probCols = sum(matrix)/sizeRows;
            end
            
            if(sizeCols ~= length(probCols))
                error('The column size is different from the assigned probabilities');
            end
            
            matrix = rand(sizeRows, sizeCols) < repmat(probCols,sizeRows,1);
            
        end
        
    end
end