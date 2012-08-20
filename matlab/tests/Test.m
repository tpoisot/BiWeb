classdef Test < handle

    properties(GetAccess = 'public', SetAccess = 'private')
        bipweb       = [];
        nulls        = {};
        nestvals      = [];
        nestvals_rows = [];
        nestvals_cols = [];
        qb_vals       = [];
        qr_vals       = [];
        nest_row_contrib = [];
        nest_col_contrib = []; 
        model        = @NullModels.NULL_CE;
        replicates   = 100;
        modul_done   = 0;
        nest_done    = 0;
        nest_contrib_done = 0;
    end
    
    %CONSTRUCTOR AND MAIN PROCEDURE ALGORITHM
    methods

        function obj = Test(webbip)
            obj.bipweb = webbip;
        end
        
        function obj = DoNulls(obj,nullmodel,replic)
           
            obj.modul_done = 0;
            obj.nest_done = 0;
            
            if(nargin == 1)
                obj.model = @NullModels.NULL_1;
                obj.replicates = 100;
            elseif(nargin == 2)
                obj.model = nullmodel;
                obj.replicates = 100;
            else
                obj.model = nullmodel;
                obj.replicates = replic;
            end
            
            obj.nulls = NullModels.NULL_MODEL(obj.bipweb.adjacency,obj.model,obj.replicates);
            
        end
        
        function obj = Nestedness(obj)
           
            if(isempty(obj.nulls))
                obj.DoNulls();
            end
            [obj.nestvals obj.nestvals_rows obj.nestvals_cols] = Test.GET_DEV_NEST(obj.bipweb,obj.nulls);
            obj.nest_done = 1;
        end
        
        function obj = Modularity(obj)
           
            if(isempty(obj.nulls))
                obj.DoNulls();
            end
            if(obj.bipweb.modules.done == 0)
                obj.bipweb.modules.Detect(100);
            end
            
            [obj.qb_vals obj.qr_vals] = Test.GET_DEV_MODUL(obj.bipweb, obj.nulls);
            obj.modul_done = 1;
        end
        
        function obj = NestednessContributions(obj)
           
            if(isempty(obj.nulls))
                obj.DoNulls();
            end
            [obj.nest_row_contrib obj.nest_col_contrib] = Test.GET_NEST_CONTRIBUTIONS(obj.bipweb.adjacency,obj.nulls);
            obj.nest_contrib_done = 1;
        end
        
    end
    
    methods(Static)
        
        function [out_b out_r] = GET_DEV_MODUL(webbip,rmatrices)
           
            wQr = webbip.modules.Qr;
            wQb = webbip.modules.Qb;
            n = length(rmatrices);
            
            Qb_random = zeros(n,1);
            Qr_random = zeros(n,1);
            
            for i = 1:n
                lp = LPBrim(rmatrices{i});
                lp.Detect(10);
                Qb_random(i) = lp.Qb;
                Qr_random(i) = lp.Qr;    
            end
            
            [hb pb cib] = ttest(Qb_random, wQb);
            [hr pr cir] = ttest(Qr_random, wQr);
            
            cib = [sum(cib)/2; cib];
            cir = [sum(cir)/2; cir];
            
            z_qb = (webbip.modules.Qb - mean(Qb_random))/std(Qb_random);
            z_qr = (webbip.modules.Qr - mean(Qr_random))/std(Qr_random);
            
            out_b.Qb = wQb; out_b.p = pb; out_b.ci = cib; out_b.zscore = z_qb;
            out_r.Qr = wQr; out_r.p = pr; out_r.ci = cir; out_r.zscore = z_qr;
            
        end
        
        function [out out_row out_col] = GET_DEV_NEST(webbip,rmatrices) 
            
            n = length(rmatrices);
            expect = zeros(n,1);
            expect_row = zeros(n,1);
            expect_col = zeros(n,1);
            
            for i = 1:n
                Nodf = NODF(rmatrices{i},webbip.nodf_strict);
                expect(i) = Nodf.nodf;
                expect_row(i) = Nodf.nodf_rows;
                expect_col(i) = Nodf.nodf_cols;
            end
            
            [h p ci] = ttest(expect, webbip.nestedness.nodf);
            [h p_row ci_row] = ttest(expect_row, webbip.nestedness.nodf_rows);
            [h p_col ci_col] = ttest(expect_col, webbip.nestedness.nodf_cols);
            
            
            ci = [sum(ci)/2; ci];
            ci_row = [sum(ci_row)/2; ci_row];
            ci_col = [sum(ci_col)/2; ci_col];
            z_nest = (webbip.nestedness.nodf - mean(expect))/std(expect);
            z_nest_row = (webbip.nestedness.nodf_rows - mean(expect_row))/std(expect_row);
            z_nest_col = (webbip.nestedness.nodf_cols - mean(expect_col))/std(expect_col);
            
            out.nodf = webbip.nestedness.nodf; out.p = p; out.ci = ci; out.zscore = z_nest;
            out_row.nodf = webbip.nestedness.nodf_rows; out_row.p = p_row; out_row.ci = ci_row; out_row.zscore = z_nest_row;
            out_col.nodf = webbip.nestedness.nodf_cols; out_col.p = p_col; out_col.ci = ci_col; out_col.zscore = z_nest_col;
            
        end
        
        function [c_rows c_cols] = GET_NEST_CONTRIBUTIONS(matrix, rmatrices)
            
            nodf = NODF(matrix);
            N = nodf.nodf;
            [n_rows n_cols] = size(matrix);
            nn = length(rmatrices);
            row_contrib = zeros(nn,n_rows);
            col_contrib = zeros(nn,n_cols);
            orig_matrix = matrix;
            for i = 1:n_rows
                temp_row = matrix(i,:);
                for j = 1:nn
                    matrix(i,:) = rmatrices{j}(i,:);
                    nodf = NODF(matrix,0);
                    row_contrib(j,i) = nodf.nodf;
                    if(isnan(nodf.nodf))
                        continue;
                    end
                end
                matrix(i,:) = temp_row;
            end
            
            for i = 1:n_cols
                temp_col = matrix(:,i);
                for j = 1:nn
                    matrix(:,i) = rmatrices{j}(:,i);
                    nodf = NODF(matrix,0);
                    col_contrib(j,i) = nodf.nodf;
                end
                matrix(:,i) = temp_col;
            end
            std_rows = std(row_contrib);
            std_cols = std(col_contrib);
            mean_rows = mean(row_contrib);
            mean_cols = mean(col_contrib);
            
            c_rows = (N - mean_rows)./std_rows;
            c_cols = (N - mean_cols)./std_cols;
            
        end
    end

end
