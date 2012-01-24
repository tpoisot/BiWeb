%==========================================================================
% Name        : LPBrim.m
% Author      : Cesar Flores
% Created     : 23/Jan/2012
% Updated     : 23/Jan/2012
% Description : Represents the LP&BRIM algorithm 
%==========================================================================
%              

classdef LPBrim < handle

    properties%(GetAccess = 'private', SetAccess = 'private')
        matrix               = [];  %Bipartite adjacency matrix
        rr                   = [];  %Red Nodes (rows) Communities matrix. Size = n_rows*CommunityQuantity
        tt                   = [];  %Blue Nodes (columns) Communities. Size = n_cols*CommunityQuantity
        n_rows               = 0;   %Number of rows
        n_cols               = 0;   %Number of columns
        n_edges              = 0;
        PP                   = [];  %Null model matrix
        bb                   = [];  %Original - Null.
        index_rows           = [];  %Register of the swaps in Rows.
        index_cols           = [];  %Register of the swaps in Cols.
        red_labels           = 0;
        blue_labels          = 0;
        trials               = 10;
        Qb                   = 0;
        N                    = 0;
    end
    
    %CONSTRUCTOR AND MAIN PROCEDURES ALGORITHM
    methods
        
        function obj = LPBrim(bipmatrix)
            
            obj.matrix = bipmatrix > 0;
            
            [obj.n_rows obj.n_cols] = size(obj.matrix);
            
            obj.trials = 10;
            
            if all(all(obj.matrix == 0))
                obj.bb = obj.matrix;
            else
                coldeg = sum(obj.matrix, 1);
                rowdeg = sum(obj.matrix, 2);
                obj.n_edges = sum(rowdeg);
                obj.bb = obj.matrix - (1/obj.n_edges) * rowdeg * coldeg;
            end
            
            obj.CalculateModularity();
            
        end
        
        
        function obj = CalculateModularity(obj)
           
            obj.LP();
            
            obj.BRIM();
            
            obj.CleanCommunities();
            
        end
            
    end
    
    methods
       
        function obj = LP(obj)
            
            obj.red_labels = 1:obj.n_rows;
            obj.blue_labels = 1:obj.n_cols;
            cmax = max(obj.n_rows,obj.n_cols);
            obj.rr = zeros(obj.n_rows,cmax);
            obj.tt = zeros(obj.n_cols,cmax);
            obj.rr(1:obj.n_rows,1:obj.n_rows) = eye(obj.n_rows);
            obj.tt(1:obj.n_cols,1:obj.n_cols) = eye(obj.n_cols);
            
            obj.n_edges = sum(sum(obj.matrix));
            
            qmax = obj.CalculateQValue();
            qmaxglobal = qmax;
            ttMaxGlobal = obj.tt;
            rrMaxGlobal = obj.rr;
            for i = 1:obj.trials
            
                cmax = max(obj.n_rows,obj.n_cols);
                obj.red_labels = 1:obj.n_rows;
                obj.blue_labels = 1:obj.n_cols;
                obj.rr = zeros(obj.n_rows,cmax);
                obj.tt = zeros(obj.n_cols,cmax);
                obj.rr(1:obj.n_rows,1:obj.n_rows) = eye(obj.n_rows);
                obj.tt(1:obj.n_cols,1:obj.n_cols) = eye(obj.n_cols);
                qmax = obj.CalculateQValue();
                rrMax = obj.rr;
                ttMax = obj.tt;

                while(1)

                    %Propagate red labels to blue labels
                    for i = 1:obj.n_cols
                        neighbors = obj.matrix(:,i);
                        labels = obj.red_labels(neighbors);
                        uniq = unique(labels);
                        nuniq = length(uniq);
                        rperm = randperm(nuniq);
                        uniq = uniq(rperm);
                        count = arrayfun(@(x) sum(x==labels), uniq);
                        [maxv index] = max(count);
                        obj.blue_labels(i) = uniq(index);
                    end

                    uniq = unique(obj.blue_labels);
                    nbluec = length(uniq);
                    newlab = arrayfun(@(x) find(uniq==x), obj.blue_labels);
                    obj.blue_labels = newlab;

                    %Propagate blue labels to red labels
                    for i = 1:obj.n_rows
                        neighbors = obj.matrix(i,:);
                        labels = obj.blue_labels(neighbors);
                        uniq = unique(labels);
                        nuniq = length(uniq);
                        rperm = randperm(nuniq);
                        uniq = uniq(rperm);
                        count = arrayfun(@(x) sum(x==labels), uniq);
                        [maxv index] = max(count);
                        obj.red_labels(i) = uniq(index);
                    end

                    uniq = unique(obj.red_labels);
                    nredc = length(uniq);
                    newlab = arrayfun(@(x) find(uniq==x), obj.red_labels);
                    obj.red_labels = newlab;

                    cmax = max(nbluec, nredc);
                    indxblue = sub2ind([obj.n_cols cmax], 1:obj.n_cols, obj.blue_labels);
                    indxred = sub2ind([obj.n_rows cmax], 1:obj.n_rows, obj.red_labels);

                    obj.rr = zeros(obj.n_rows,cmax);
                    obj.tt = zeros(obj.n_cols,cmax);
                    obj.rr(indxred) = 1;
                    obj.tt(indxblue) = 1;

                    q = obj.CalculateQValue();

                    if(q > qmax) 
                        qmax = q;
                        rrMax = obj.rr;
                        ttMax = obj.tt;
                    else
                        break;
                    end

                end
                
                if(qmax > qmaxglobal)
                    qmaxglobal = qmax;
                    rrMaxGlobal = rrMax;
                    ttMaxGlobal = ttMax;
                end
            end
            obj.Qb = qmaxglobal;
            obj.rr = rrMaxGlobal;
            obj.tt = ttMaxGlobal;
        end
        
        function obj = BRIM(obj)
           
            inducingBlueFlag = 1;
            obj.n_edges = sum(obj.matrix(:));
            
            preQ = obj.CalculateQValue();
            while(1)
               
                if(inducingBlueFlag)
                    obj.AssignBlueNodes();
                else
                    obj.AssignRedNodes();
                end
                
                newQ = obj.CalculateQValue();
                
                if(newQ <= preQ)
                    break;
                else
                    preQ = newQ;
                    inducingBlueFlag = ~inducingBlueFlag;
                end
                
            end
        end
        
        function q = CalculateQValue(obj)
            
            q = trace(obj.rr' * obj.bb * obj.tt) / obj.n_edges;
            
        end
        
        function obj = AssignRedNodes(obj)
            rr = obj.bb*obj.tt;
            [maxx,maxi]=max(rr');
            tmps=size(rr);
            obj.rr=zeros(obj.n_rows,tmps(2));
            idx=sub2ind(size(obj.rr),1:obj.n_rows,maxi);
            obj.rr(idx)=1;
        end
        
        function obj = AssignBlueNodes(obj)
            
            rr = obj.bb'*obj.rr;
            [maxx,maxi]=max(rr');
            tmps=size(rr);
            obj.tt=zeros(obj.n_cols,tmps(2));
            idx=sub2ind(size(obj.tt),1:obj.n_cols,maxi);
            obj.tt(idx)=1;
        end
        
        
        function obj = CleanCommunities(obj)
           
            %CLEAN THE COMMUNITIES
            com1 = find(any(obj.tt));
            com2 = find(any(obj.rr));
            
            if(length(com1) >= length(com2)) %Not sure if this validation is necesarry.
                obj.rr = obj.rr(:,com1);
                obj.tt = obj.tt(:,com1);
            else
                obj.rr = obj.rr(:,com2);
                obj.tt = obj.tt(:,com2);
            end
            
            obj.N = size(obj.rr,2);
            
        end
        
    end
    
end