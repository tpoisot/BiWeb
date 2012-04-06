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
        name                  = {};
        tests                 = {};
        nodf_strict           = 1;
        row_labels            = {};
        col_labels            = {};
    end
    
    methods
        
        function obj = Bipartite(web,namebip)
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
               error('You need to specify a double matrix or a txt file in matrix format');
            end
            
            if(nargin == 1 && isa(web,'double')); namebip = 'No name'; end;
            
            if(nargin == 1 && isa(web,'char'))
                [paths namefile ext] = fileparts(web);
                namebip = namefile;
                web = dlmread(web,' ');
            end
            
            emptyrows = sum(web,2)==0;
            emptycols = sum(web,1)==0;
            
            web(emptyrows,:) = [];
            web(:,emptycols) = [];
            
            obj.name = namebip;
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
            obj.tests = Test(obj);
            
            for i = 1:obj.n_rows; obj.row_labels{i} = sprintf('row_%03i',i); end;
            for i = 1:obj.n_cols; obj.col_labels{i} = sprintf('col_%03i',i); end;
            
        end
        
        function obj = SpeciesLevel(obj, toscreen, tofile, filename)
           
            aindex = obj.name==' ';
            namebip = obj.name;
            namebip(aindex) = '_';
            
            if(nargin == 1)
                toscreen = 1;
                tofile = 1;
                filename = [namebip, '_sp.txt'];
            elseif(nargin == 2)
                tofile = 1;
                filename = [namebip, '_sp.txt'];
            elseif (nargin == 3)
                filename = [namebip, '_sp.txt'];
            end
            
            header = '';
            colwitdth = 9;
            formatstring = ['%-',num2str(colwitdth),'s'];
            formatint = ['%-',num2str(colwitdth),'d'];
            formatfloat = ['%-',num2str(colwitdth),'.3f'];
            estring = '----';
            
            header = [header, sprintf(repmat(formatstring,1,6), 'Specie','Level','Degree','SPE','SSI','RR')];
            if(obj.modules.done); header = [header, sprintf(repmat(formatstring,1,1), 'Module')]; end;
            header = [header, '\n'];
            
            info = '';
            for i = 1:obj.n_rows
                aindex = obj.row_labels{i}==' ';
                label_name = obj.row_labels{i};
                label_name(aindex) = '_';
                info = [info, sprintf(formatstring,label_name(1:min(colwitdth,length(label_name))))];
                info = [info, sprintf(formatstring, 'top')];
                info = [info, sprintf(formatint, obj.row_degrees(i))];
                info = [info, sprintf(formatfloat, obj.specificity(i))];
                info = [info, sprintf(formatfloat, obj.ssi(i))];
                info = [info, sprintf(formatfloat, obj.rr(i))];
                if(obj.modules.done); info = [info, sprintf(formatint, obj.modules.row_modules(i))]; end;
                info = [info, '\n'];
            end
            
            for i = 1:obj.n_cols
                aindex = obj.col_labels{i}==' ';
                label_name = obj.col_labels{i};
                label_name(aindex) = '_';
                info = [info, sprintf(formatstring,label_name(1:min(colwitdth,length(label_name))))];
                info = [info, sprintf(formatstring, 'bottom')];
                info = [info, sprintf(formatint, obj.col_degrees(i))];
                info = [info, sprintf(repmat(formatstring,1,3),estring,estring,estring)];
                if(obj.modules.done); info = [info, sprintf(formatint, obj.modules.col_modules(i))]; end;
                info = [info, '\n'];
            end
            
            if(toscreen)
                fprintf(header);
                fprintf(info);
            end
            
            if(tofile)
                fid = fopen([namebip,'-sp.txt'],'w');
                fprintf(fid,header);
                fprintf(fid,info);
                fclose(fid);
            end
            
        end
        
        function obj = NetworkLevelSingle(obj, toscreen, tofile, filename)
            %To complete, maybe have to delete at the end of the day
            if(nargin == 1)
                toscreen = 1;
                tofile = 1;
                filename = [obj.name, '_networklevel.txt'];
            elseif(nargin == 2)
                tofile = 1;
                filename = [obj.name, '_networklevel.txt'];
            elseif (nargin == 3)
                filename = [obj.name, '_networklevel.txt'];
            end
            format = '%5.2f';
            formatint = '%5d';
            fileexist = exist(filename,'file');
            header = '';
            if(~fileexist)
                header = [header, '\nGeneral properties\n'];
                header = [header,   '------------------\n'];
                header = [header, '\tNetwork Name:      ', obj.name,'\n'];
                header = [header, '\tRows:              ', num2str(obj.n_rows,formatint),'\n'];
                header = [header, '\tCols:              ', num2str(obj.n_cols,formatint),'\n'];
                header = [header, '\tSize:              ', num2str(obj.size_webmatrix,formatint),'\n'];
                header = [header, '\tInteractions:      ', num2str(obj.num_edges,formatint),'\n'];
                header = [header, '\tConnectance:       ', num2str(obj.connectance,format),'\n'];
                header = [header, '\t<Specificity>:     ', num2str(mean(obj.specificity),format),'\n'];
                header = [header, '\t<Res. Range.>:     ', num2str(mean(obj.rr),format),'\n'];
                header = [header, '\tNODF:              ', num2str(obj.nestedness.nodf,format),'\n'];
                header = [header, '\tNODF Rows:         ', num2str(obj.nestedness.nodf_rows,format),'\n'];
                header = [header, '\tNODF Cols:         ', num2str(obj.nestedness.nodf_cols,format),'\n'];
                
            end
            
            if(toscreen)
                fprintf(header);
            end
            
        end
        
        function obj = NetworkLevel(obj, toscreen, tofile)
            
            filename = 'networklevel.txt';
            
            if(nargin == 1)
                toscreen = 1;
                tofile = 1;
            else(nargin == 2)
                tofile = 1;
            end
            
            header = '';
            colwitdth = 9;
            formatstring = ['%-',num2str(colwitdth),'s'];
            formatint = ['%-',num2str(colwitdth),'d'];
            formatfloat = ['%-',num2str(colwitdth),'.3f'];
            estring = '----';
            
            header = [header, sprintf(repmat(formatstring,1,5), 'Name','Conn','Size','n_rows','n_cols')];
            header = [header, sprintf(repmat(formatstring,1,5), 'nodf','nodf_r','nodf_c','aspe','arr')];
            header = [header, sprintf(repmat(formatstring,1,5), 'resp','inc','mN','mQb','mQr')];
            header = [header, sprintf(repmat(formatstring,1,5), 'n_sim','n_pval','nicLow','nicUp','nco_sim')];
            header = [header, sprintf(repmat(formatstring,1,4), 'nco_pval','ncoicLow','ncoicUp','nro_sim')];
            header = [header, sprintf(repmat(formatstring,1,3), 'nro_pval','nroicLow','nroicUp')];
            header = [header, sprintf(repmat(formatstring,1,4), 'mQr_sim','mQr_val','mQrIclow','mQrIcup')];
            header = [header, sprintf(repmat(formatstring,1,4), 'mQb_sim','mQb_val','mQbIclow','mQbIcup')];
            header = [header, sprintf(repmat(formatstring,1,3), 'tReps','null_name')];
            header = [header, '\n'];
            %header = [header, '-----------------------------------------------------------------\n'];
            
            [R I] = SpeFunc.IR(obj.webmatrix);
            aindex = obj.name==' ';
            namebip = obj.name;
            namebip(aindex) = '_';
            info = sprintf(formatstring,namebip(1:min(colwitdth-1,length(namebip))));
            info = [info, sprintf(formatfloat, obj.connectance)];
            info = [info, sprintf(formatint, obj.size_webmatrix)];
            info = [info, sprintf(formatint, obj.n_rows)];
            info = [info, sprintf(formatint, obj.n_cols)];
            info = [info, sprintf(formatfloat, obj.nestedness.nodf)];
            info = [info, sprintf(formatfloat, obj.nestedness.nodf_rows)];
            info = [info, sprintf(formatfloat, obj.nestedness.nodf_cols)];
            info = [info, sprintf(formatfloat, mean(obj.specificity))];
            info = [info, sprintf(formatfloat, mean(obj.rr))];
            info = [info, sprintf(formatfloat, R)];
            info = [info, sprintf(formatfloat, I)];
            
            if(obj.modules.done)
                info = [info, sprintf(formatint, obj.modules.N)];
                info = [info, sprintf(formatfloat, obj.modules.Qb)];
                info = [info, sprintf(formatfloat, obj.modules.Qr)];
            else
                info = [info, sprintf(repmat(formatstring,1,3),estring,estring,estring)];
            end
            
            if(obj.tests.nest_done)
                info = [info, sprintf(formatfloat, obj.tests.devnest(3))];
                info = [info, sprintf(formatfloat, obj.tests.devnest(2))];
                info = [info, sprintf(formatfloat, obj.tests.devnest(4))];
                info = [info, sprintf(formatfloat, obj.tests.devnest(5))];
                info = [info, sprintf(formatfloat, obj.tests.devnest_cols(3))];
                info = [info, sprintf(formatfloat, obj.tests.devnest_cols(2))];
                info = [info, sprintf(formatfloat, obj.tests.devnest_cols(4))];
                info = [info, sprintf(formatfloat, obj.tests.devnest_cols(5))];
                info = [info, sprintf(formatfloat, obj.tests.devnest_rows(3))];
                info = [info, sprintf(formatfloat, obj.tests.devnest_rows(2))];
                info = [info, sprintf(formatfloat, obj.tests.devnest_rows(4))];
                info = [info, sprintf(formatfloat, obj.tests.devnest_rows(5))];
            else
                info = [info, sprintf(repmat(formatstring,1,12),estring,estring,estring,estring,estring,estring, ...
                    estring,estring,estring,estring,estring,estring)];
            end
            
            if(obj.tests.modul_done)
                info = [info, sprintf(formatfloat, obj.tests.dev_qr(3))];
                info = [info, sprintf(formatfloat, obj.tests.dev_qr(2))];
                info = [info, sprintf(formatfloat, obj.tests.dev_qr(4))];
                info = [info, sprintf(formatfloat, obj.tests.dev_qr(5))];

                info = [info, sprintf(formatfloat, obj.tests.dev_qb(3))];
                info = [info, sprintf(formatfloat, obj.tests.dev_qb(2))];
                info = [info, sprintf(formatfloat, obj.tests.dev_qb(4))];
                info = [info, sprintf(formatfloat, obj.tests.dev_qb(5))];
            else
                info = [info, sprintf(repmat(formatstring,1,8),estring,estring,estring, ... 
                    estring,estring,estring,estring,estring)];
            end
            
            if(obj.tests.modul_done || obj.tests.nest_done)
                info = [info, sprintf(formatint, obj.tests.replicates)];
                info = [info, sprintf(formatstring, func2str(obj.tests.model))];
            else
                info = [info, sprintf(repmat(formatstring,1,2),estring,estring)];
            end
            
            info = [info, '\n'];
            
            if(toscreen)
                fprintf(header);
                fprintf(info);
            end
            
            if(tofile)
                if(~exist(filename,'file')) 
                    fid = fopen(filename,'w');
                    fprintf(fid, header);
                else
                    fid = fopen(filename,'a');
                end
                fprintf(fid,info);
            end
            
        end
        
    end
       
end

