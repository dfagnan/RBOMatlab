clear all;
addpath '..'
test_scripts = {@init_params_fn_NBTx_A @init_params_fn_NBTx_B @init_params_fn_AEA_A1_guarantee @init_params_fn_DDT @init_params_fn_DDT_geometric};
expected = [0.091 0.106 0.089 0.134 0.134];
NSIMUS = 10000;
for i =1:length(expected)
    chance = 2;
    while (chance >0)
        sprintf('Running Test %d\n',i)
        params = test_scripts{i}();
        params.simu.NSIMUS=NSIMUS;
        results=simulate_CFs_fn(1,params);

        HORIZON    = length(results.cash(1,:));
        NGUARPERS  = size(results.guarantee_draw,2);
        NYEARS     = HORIZON/2;

        ROE_raw        = results.ROE;
        ROE_annualized = ((1+ROE_raw).^(1/NYEARS))-1;
        ROE_mean       = mean(ROE_raw);
        ROE_mean_a     = mean(ROE_annualized);
        ROE_std_a      = std(ROE_annualized);
        sprintf('Test gap: %f',abs(ROE_mean_a - expected(i))/(ROE_std_a/sqrt(params.simu.NSIMUS)))
        if abs(ROE_mean_a - expected(i))/(ROE_std_a/sqrt(params.simu.NSIMUS))>1.645
            chance = chance - 1;
            keyboard
            assert(chance > 0, 'Test Failed!');
        else
            chance = 0;
        end
        sprintf('Test %d complete.\n',i)
    end
end

display('All tests complete!')