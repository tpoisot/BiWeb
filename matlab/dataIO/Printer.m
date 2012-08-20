classdef Printer < handle
    
    properties
        
        bipweb = {};
        
    end
    
    methods
        
        function obj = Printer(webbip)
            
            obj.bipweb = webbip;
            
        end
        
        function obj = SpeciesLevel(obj, toscreen, tofile)
           
            webbip = obj.bipweb;
            
            aindex = webbip.name==' ';
            namebip = webbip.name;
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
            if(webbip.modules.done); header = [header, sprintf(repmat(formatstring,1,1), 'Module')]; end;
            header = [header, '\n'];
            
            info = '';
            for i = 1:webbip.n_rows
                aindex = webbip.row_labels{i}==' ';
                label_name = webbip.row_labels{i};
                label_name(aindex) = '_';
                info = [info, sprintf(formatstring,label_name(1:min(colwitdth,length(label_name))))];
                info = [info, sprintf(formatstring, 'top')];
                info = [info, sprintf(formatint, webbip.row_degrees(i))];
                info = [info, sprintf(formatfloat, webbip.specificity(i))];
                info = [info, sprintf(formatfloat, webbip.ssi(i))];
                info = [info, sprintf(formatfloat, webbip.rr(i))];
                if(webbip.modules.done); info = [info, sprintf(formatint, webbip.modules.row_modules(i))]; end;
                info = [info, '\n'];
            end
            
            for i = 1:webbip.n_cols
                aindex = webbip.col_labels{i}==' ';
                label_name = webbip.col_labels{i};
                label_name(aindex) = '_';
                info = [info, sprintf(formatstring,label_name(1:min(colwitdth,length(label_name))))];
                info = [info, sprintf(formatstring, 'bottom')];
                info = [info, sprintf(formatint, webbip.col_degrees(i))];
                info = [info, sprintf(repmat(formatstring,1,3),estring,estring,estring)];
                if(webbip.modules.done); info = [info, sprintf(formatint, webbip.modules.col_modules(i))]; end;
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
            webbip = obj.bipweb;
            
            if(nargin == 1)
                toscreen = 1;
                tofile = 1;
                filename = [webbip.name, '_networklevel.txt'];
            elseif(nargin == 2)
                tofile = 1;
                filename = [webbip.name, '_networklevel.txt'];
            elseif (nargin == 3)
                filename = [webbip.name, '_networklevel.txt'];
            end
            format = '%5.2f';
            formatint = '%5d';
            fileexist = exist(filename,'file');
            header = '';
            if(~fileexist)
                header = [header, '\nGeneral properties\n'];
                header = [header,   '------------------\n'];
                header = [header, '\tNetwork Name:      ', webbip.name,'\n'];
                header = [header, '\tRows:              ', num2str(webbip.n_rows,formatint),'\n'];
                header = [header, '\tCols:              ', num2str(webbip.n_cols,formatint),'\n'];
                header = [header, '\tSize:              ', num2str(webbip.size_webmatrix,formatint),'\n'];
                header = [header, '\tInteractions:      ', num2str(webbip.num_edges,formatint),'\n'];
                header = [header, '\tConnectance:       ', num2str(webbip.connectance,format),'\n'];
                header = [header, '\t<Specificity>:     ', num2str(mean(webbip.specificity),format),'\n'];
                header = [header, '\t<Res. Range.>:     ', num2str(mean(webbip.rr),format),'\n'];
                header = [header, '\tNODF:              ', num2str(webbip.nestedness.nodf,format),'\n'];
                header = [header, '\tNODF Rows:         ', num2str(webbip.nestedness.nodf_rows,format),'\n'];
                header = [header, '\tNODF Cols:         ', num2str(webbip.nestedness.nodf_cols,format),'\n'];                
            end
            
            if(toscreen)
                fprintf(header);
            end
            
        end
        
        function obj = NetworkLevel(obj, toscreen, tofile)
            webbip = obj.bipweb;
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
            header = [header, sprintf(repmat(formatstring,1,6), 'n_sim','n_pval','nicLow','nicUp','n_z','nco_sim')];
            header = [header, sprintf(repmat(formatstring,1,5), 'nco_pval','ncoicLow','ncoicUp','nco_z','nro_sim')];
            header = [header, sprintf(repmat(formatstring,1,4), 'nro_pval','nroicLow','nroicUp','nro_z')];
            header = [header, sprintf(repmat(formatstring,1,4), 'mQr_sim','mQr_val','mQrIclow','mQrIcup','mQr_z')];
            header = [header, sprintf(repmat(formatstring,1,4), 'mQb_sim','mQb_val','mQbIclow','mQbIcup','mQb_z')];
            header = [header, sprintf(repmat(formatstring,1,3), 'tReps','null_name')];
            header = [header, '\n'];
            %header = [header, '-----------------------------------------------------------------\n'];
            
            [R I] = SpeFunc.IR(webbip.webmatrix);
            aindex = webbip.name==' ';
            namebip = webbip.name;
            namebip(aindex) = '_';
            info = sprintf(formatstring,namebip(1:min(colwitdth-1,length(namebip))));
            info = [info, sprintf(formatfloat, webbip.connectance)];
            info = [info, sprintf(formatint, webbip.size_webmatrix)];
            info = [info, sprintf(formatint, webbip.n_rows)];
            info = [info, sprintf(formatint, webbip.n_cols)];
            info = [info, sprintf(formatfloat, webbip.nestedness.nodf)];
            info = [info, sprintf(formatfloat, webbip.nestedness.nodf_rows)];
            info = [info, sprintf(formatfloat, webbip.nestedness.nodf_cols)];
            info = [info, sprintf(formatfloat, mean(webbip.specificity))];
            info = [info, sprintf(formatfloat, mean(webbip.rr))];
            info = [info, sprintf(formatfloat, R)];
            info = [info, sprintf(formatfloat, I)];
            
            if(webbip.modules.done)
                info = [info, sprintf(formatint, webbip.modules.N)];
                info = [info, sprintf(formatfloat, webbip.modules.Qb)];
                info = [info, sprintf(formatfloat, webbip.modules.Qr)];
            else
                info = [info, sprintf(repmat(formatstring,1,3),estring,estring,estring)];
            end
            
            if(webbip.tests.nest_done)
                info = [info, sprintf(formatfloat, webbip.tests.nestvals.ci(1))];
                info = [info, sprintf(formatfloat, webbip.tests.nestvals.p)];
                info = [info, sprintf(formatfloat, webbip.tests.nestvals.ci(2))];
                info = [info, sprintf(formatfloat, webbip.tests.nestvals.ci(3))];
                info = [info, sprintf(formatfloat, webbip.tests.nestvals.zscore)];
                info = [info, sprintf(formatfloat, webbip.tests.nestvals_cols.ci(1))];
                info = [info, sprintf(formatfloat, webbip.tests.nestvals_cols.p)];
                info = [info, sprintf(formatfloat, webbip.tests.nestvals_cols.ci(2))];
                info = [info, sprintf(formatfloat, webbip.tests.nestvals_cols.ci(3))];
                info = [info, sprintf(formatfloat, webbip.tests.nestvals_cols.zscore)];
                info = [info, sprintf(formatfloat, webbip.tests.nestvals_rows.ci(1))];
                info = [info, sprintf(formatfloat, webbip.tests.nestvals_rows.p)];
                info = [info, sprintf(formatfloat, webbip.tests.nestvals_rows.ci(2))];
                info = [info, sprintf(formatfloat, webbip.tests.nestvals_rows.ci(3))];
                info = [info, sprintf(formatfloat, webbip.tests.nestvals_rows.zscore)];
            else
                info = [info, sprintf(repmat(formatstring,1,12),estring,estring,estring,estring,estring,estring, ...
                    estring,estring,estring,estring,estring,estring)];
            end
            
            if(webbip.tests.modul_done)
                info = [info, sprintf(formatfloat, webbip.tests.qr_vals.ci(1))];
                info = [info, sprintf(formatfloat, webbip.tests.qr_vals.p)];
                info = [info, sprintf(formatfloat, webbip.tests.qr_vals.ci(2))];
                info = [info, sprintf(formatfloat, webbip.tests.qr_vals.ci(3))];
                info = [info, sprintf(formatfloat, webbip.tests.qr_vals.zscore)];

                info = [info, sprintf(formatfloat, webbip.tests.qb_vals.ci(1))];
                info = [info, sprintf(formatfloat, webbip.tests.qb_vals.p)];
                info = [info, sprintf(formatfloat, webbip.tests.qb_vals.ci(2))];
                info = [info, sprintf(formatfloat, webbip.tests.qb_vals.ci(3))];
                info = [info, sprintf(formatfloat, webbip.tests.qb_vals.zscore)];
            else
                info = [info, sprintf(repmat(formatstring,1,8),estring,estring,estring, ... 
                    estring,estring,estring,estring,estring)];
            end
            
            if(webbip.tests.modul_done || webbip.tests.nest_done)
                info = [info, sprintf(formatint, webbip.tests.replicates)];
                info = [info, sprintf(formatstring, func2str(webbip.tests.model))];
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