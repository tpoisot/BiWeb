classdef Test < handle

    properties(GetAccess = 'public', SetAccess = 'private')
        bipweb       = [];
        nulls        = {};
        devnest      = [];
        devnest_rows = [];
        devnest_cols = [];
        dev_qb        = [];
        dev_qr        = [];
        model        = @NullModels.NULL_CE;
        replicates   = 100;
        modul_done   = 0;
        nest_done   = 0;
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
            [obj.devnest obj.devnest_rows obj.devnest_cols] = Test.GET_DEV_NEST(obj.bipweb,obj.nulls);
            obj.nest_done = 1;
        end
        
        function obj = Modularity(obj)
           
            if(isempty(obj.nulls))
                obj.DoNulls();
            end
            if(obj.bipweb.modules.done == 0)
                obj.bipweb.modules.Detect(100);
            end
            
            [obj.dev_qb obj.dev_qr] = Test.GET_DEV_MODUL(obj.bipweb, obj.nulls);
            obj.modul_done = 1;
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
            
            out_b = [wQb; pb; cib];
            out_r = [wQr; pr; cir];
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
            
            out = [webbip.nestedness.nodf; p; ci];
            out_row = [webbip.nestedness.nodf_rows; p_row; ci_row];
            out_col = [webbip.nestedness.nodf_cols; p_col; ci_col];
            
        end
    end

end
