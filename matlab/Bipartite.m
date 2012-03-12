%==========================================================================
% Name        : BiWeb.m
% Author      : Cesar Flores
% Created     : 23/Jan/2012
% Updated     : 23/Jan/2012
% Description : Represent a complete data structure of a bipartite network.
%==========================================================================

classdef Bipartite

    properties
        webmatrix             = [];%Interaction matrix (not neccesearly a binary matrix).
        adjacency             = [];%Adjacency matrix (binary matrix)
        modules               = {};
        nestedness            = {};
        num_edges             = 0; %Number of edges
        connectance           = 0; %Fill of the webmatrix
        n_rows                = 0; %Number of rows
        n_cols                = 0; %Number of columns
        size_webmatrix        = 0;
        specificity           = [];
        row_degrees           = [];
        col_degrees           = [];
        rr                    = [];
        ssi                   = [];
        bperf                 = [];
    end
    
    methods
        
        function obj = Bipartite(web)
            %NetworkBipartite(webmatrix): Creates bipartite network from an
            %specific bipartite adjacency webmffffffffffffffsatrix (teams vs actors)
            %
            % Input arguments -
            %   webmatrix       = This is the adjacency webmatrix. The files
            %                     will represent the teams and the columns
            %                     the actors
            % Output          -
            %   obj             = An object of the class NetworkBipartite
            %
            % Example:
            %   A = stats.phage_host;
            %   netbip = NetworkBipartite(A);
            if(nargin == 0)
               error('You need to specify the webmatrix bipartite network');
            end
            
            obj.webmatrix = web;
            
            %General Info
            [obj.n_rows obj.n_cols] = size(web);
            obj.size_webmatrix = obj.n_rows * obj.n_cols;
            
            %Connectance
            obj.adjacency = obj.webmatrix > 0;
            obj.num_edges = sum(obj.adjacency(:));
            obj.connectance = sum(obj.adjacency(:))/numel(obj.adjacency);
            
            %Specifity and all
            obj.row_degrees = sum(obj.adjacency,2);
            obj.col_degrees = sum(obj.adjacency,1)';
            obj.specificity = SpeFunc.SPECIFICITY(obj.webmatrix);
            obj.rr = SpeFunc.RESOURCE_RANGE(obj.webmatrix);
            obj.ssi = SpeFunc.SPECIES_SPECIFICITY_INDEX(obj.webmatrix);
            obj.bperf = max(obj.webmatrix)';
            
            
            obj.nestedness = NODF(obj.webmatrix);
            obj.modules = LPBrim(obj.webmatrix);
            
        end
        
    end
       
end

