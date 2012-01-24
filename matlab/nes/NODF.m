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
        nodf_up     = [];
        nodf_low    = [];
    end
    
    %CONSTRUCTOR AND MAIN PROCEDURE ALGORITHM
    methods

        function obj = NODF(bipmatrix)
            obj.matrix = bipmatrix; %Normalize the matrix
            [obj.n_rows obj.n_cols] = size(bipmatrix);  %Number of Rows
            obj.index_rows = 1:obj.n_rows;
            obj.index_cols = 1:obj.n_cols;
            
            obj.CalculateNestedness();
        end
        
        function obj = CalculateNestedness(obj)
           
            m = obj.n_rows;
            n = obj.n_cols;
            denom = n*(n-1)/2 + m*(m-1)/2;
            
            obj.SortMatrix();
            obj.CalculateNpaired();
            obj.nodf = obj.nodf / (100*denom);
            obj.nodf_up = obj.nodf_up / (100*denom);
            obj.nodf_low = obj.nodf_low / (100*denom);
               
        end
    end

    methods
       
        function obj = CalculateNpaired(obj)

            sumrows = sum(obj.matrix,2);
            sumcols = sum(obj.matrix,1);
            
            obj.nodf = 0;
            obj.nodf_up = 0;
            obj.nodf_low = 0;
            %Fill for rows
            for i = 1:obj.n_rows
                for j = i+1:obj.n_cols
                    if( sumrows(j) < sumrows(i) && sumrows(j) > 0 )
                        obj.nodf_up = obj.nodf_up + 100*sum(obj.matrix(i, obj.matrix(j,:)==1))/sumrows(j);
                    end
                end
            end
            
            %Fill for columns
            for k = 1:obj.n_cols
                for l = k+1:obj.n_cols
                    if( sumcols(l) < sumcols(k) && sumcols(l) > 0 )
                        obj.nodf_low = obj.nodf_low + 100*sum(obj.matrix(obj.matrix(:,l)==1,k))/sumcols(l);
                    end
                end
            end
            
            obj.nodf = obj.nodf_low + obj.nodf_up;
        end
        
        function obj = SortMatrix(obj)
            
            [values obj.index_rows] = sort(sum(obj.matrix,2),'descend');
            [values obj.index_cols] = sort(sum(obj.matrix,1),'descend');
        
            obj.matrix = obj.matrix(obj.index_rows, obj.index_cols);
            
        end
        
    end
end
