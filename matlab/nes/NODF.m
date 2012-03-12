%==========================================================================
% Name        : NODF.m
% Author      : Cesar Flores
% Created     : 23/Jan/2012
% Updated     : 23/Jan/2012
% Description : Represents the NODF algorithm
%==========================================================================

classdef NODF < handle

    properties(GetAccess = 'public', SetAccess = 'private')
        matrix      = {};
        n_rows      = [];
        n_cols      = [];
        index_rows  = [];
        index_cols  = [];
        nodf        = [];
        nodf_rows   = [];
        nodf_cols   = [];
    end
    
    %CONSTRUCTOR AND MAIN PROCEDURE ALGORITHM
    methods

        function obj = NODF(bipmatrix)
            obj.matrix = bipmatrix>0; %Normalize the matrix
            [obj.n_rows obj.n_cols] = size(bipmatrix);  %Number of Rows
            obj.index_rows = 1:obj.n_rows;
            obj.index_cols = 1:obj.n_cols;
            
            obj.CalculateNestedness();
        end
        
        function obj = CalculateNestedness(obj)
           
            m = obj.n_rows;
            n = obj.n_cols;
            dm = m*(m-1)/2;
            dn = n*(n-1)/2;
            denom = dm + dn;
            
            obj.SortMatrix();
            obj.CalculateNpaired();
            obj.nodf = obj.nodf / (denom);
            obj.nodf_rows = obj.nodf_rows / (dm);
            obj.nodf_cols = obj.nodf_cols / (dn);
               
        end
    end

    methods
       
        function obj = CalculateNpaired(obj)

            sumrows = sum(obj.matrix,2);
            sumcols = sum(obj.matrix,1);
            
            obj.nodf = 0;
            obj.nodf_rows = 0;
            obj.nodf_cols = 0;
            %Fill for rows
            for i = 1:obj.n_rows
                for j = i+1:obj.n_rows
                    if( sumrows(j) < sumrows(i) && sumrows(j) > 0 )
                        obj.nodf_rows = obj.nodf_rows + 100*sum(obj.matrix(i, obj.matrix(j,:)==1))/sumrows(j);
                    end
                end
            end
            
            %Fill for columns
            for k = 1:obj.n_cols
                for l = k+1:obj.n_cols
                    if( sumcols(l) < sumcols(k) && sumcols(l) > 0 )
                        obj.nodf_cols = obj.nodf_cols + 100*sum(obj.matrix(obj.matrix(:,l)==1,k))/sumcols(l);
                    end
                end
            end
            
            obj.nodf = obj.nodf_cols + obj.nodf_rows;
        end
        
        function obj = SortMatrix(obj)
            
            [values obj.index_rows] = sort(sum(obj.matrix,2),'descend');
            [values obj.index_cols] = sort(sum(obj.matrix,1),'descend');
        
            obj.matrix = obj.matrix(obj.index_rows, obj.index_cols);
            
        end
        
    end
end
