function [] = test_admm(prob_name)

[At,b,c,K] = fromsdpa(prob_name);
At =At';

for i = 0:2
    params.use_indirect = i;
    for j = 1:2
      parCoLO.domain = j; parCoLO.range = 0; parCoLO.EQorLMI = j;
      admm_r(At,b,c,K,[],parCoLO);
    end
end
