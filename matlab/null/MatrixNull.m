classdef MatrixNull < handle
    
    methods (Access = private)
    %private so that you can't instatiate.
        function out = MatrixNull

        end
    end 
   
    methods(Static)
        
        function [matrix nested_sites modular_sites] = NESTED_TO_MODULAR(n_size, n_links)
           
            n = n_size;
            matrix = flipud(tril(ones(n)));
            nested_sites = [];
            modular_sites = [];
            for i = 1:n/2
                %matrix( (i-1)*n + n/2 + 1 : n*i-i+1) = 1;
                %matrix( n/2*n + (i-1)*n + 1 : n/2*n + (i-1)*n + 1 + n/2 - i ) = 1;
                %matrix( n/2*n + (i-1)*n + n/2 : n/2*n + (i-1)*n + n) = 1;
                nested_sites = [nested_sites (i-1)*n + n/2 + 1 : n*i-i+1];
                nested_sites = [nested_sites n/2*n + (i-1)*n + 1 : n/2*n + (i-1)*n + 1 + n/2 - i];
                modular_sites = [modular_sites n/2*n + (i-1)*n + n/2 + 1: n/2*n + (i-1)*n + n ];
            end
            
            n_nest = length(nested_sites); n_modul = length(modular_sites);
            
            random_nest_switched = randsample(nested_sites,min(n_nest,n_links),0);
            random_modul_switched = randsample(modular_sites,min(n_modul,n_links),0);
            
            matrix(random_nest_switched) = 0;
            matrix(random_modul_switched) = 1;
            
            imagesc(matrix);
            
        end
        
        function matrix = SORT_MATRIX(matrix)
           
            [a ir] = sort(sum(matrix,2),'descend');
            [a ic] = sort(sum(matrix,1),'descend');
           
            matrix = matrix(ir,ic);
            
        end
        
        function matrix = RANDOM_SORT(matrix)
           
            [n_rows n_cols] = size(matrix);
            
            matrix = matrix(randperm(n_rows),randperm(n_cols));
            
        end
             
    end
end