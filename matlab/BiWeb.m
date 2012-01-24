%==========================================================================
% Name        : BiWeb.m
% Author      : Cesar Flores
% Created     : 23/Jan/2012
% Updated     : 23/Jan/2012
% Description : Represent a complete data structure of a bipartite network.
%==========================================================================

classdef BiWeb

    properties
        matrix                = [];
        modules               = {};
        nestedness            = {};
        num_edges             = 0; %Number of edges
        fill                  = 0; %Fill of the matrix
        n_rows                = 0; %Number of rows
        n_cols                = 0; %Number of columns
    end
    
    methods
        
        function obj = BiWeb(matrix)
            %NetworkBipartite(matrix): Creates bipartite network from an
            %specific bipartite adjacency matrix (teams vs actors)
            %
            % Input arguments -
            %   matrix          = This is the adjacency matrix. The files
            %                     will represent the teams and the columns
            %                     the actors
            % Output          -
            %   obj             = An object of the class NetworkBipartite
            %
            % Example:
            %   A = stats.phage_host;
            %   netbip = NetworkBipartite(A);
            if(nargin == 0)
               error('You need to specify the matrix bipartite network');
            end
            
            [obj.n_rows obj.n_cols] = size(matrix);
            
            obj.matrix = matrix > 0;
            obj.fill = sum(sum(obj.matrix>0))/(obj.n_rows*obj.n_cols);
            obj.num_edges = sum(sum(obj.matrix));
            
            obj.nestedness = NODF(obj.matrix);
            obj.modules = LPBrim(obj.matrix);
            
        end
        
    end
       
end

