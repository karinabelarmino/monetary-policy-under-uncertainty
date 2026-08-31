Em modelos bayesianos, o resultado começa antes da estimação: começa na prior.

No artigo que publicamos na *Central Bank Review*, o BVAR baseline combinou uma Minnesota hierárquica com priors dummy de soma dos coeficientes e raiz unitária.

No novo release do projeto, mantive exatamente a mesma base, as cinco variáveis, duas defasagens e a identificação por restrições de sinal. Mudei apenas a especificação da prior.

Comparei cinco alternativas em 36 previsões de um passo à frente, de julho de 2019 a junho de 2022.

O resultado não produziu uma “vencedora” em todos os critérios — e essa é justamente a parte interessante:

- A Minnesota com média de passeio aleatório e sem priors dummy teve o melhor desempenho conjunto em MAE e CRPS.
- A Minnesota com média zero apresentou o menor RMSE.
- Em relação ao baseline publicado, a primeira alternativa reduziu o CRPS padronizado em 3,8% nesse período.

Isso não torna a prior original “errada”. Mostra que uma escolha adequada ao objetivo estrutural de um artigo pode não ser a melhor para previsão fora da amostra em um período que inclui a pandemia e o início de um forte ciclo de aperto monetário.

Ao mesmo tempo, as principais respostas ao impulso permaneceram muito próximas entre as duas especificações. A prior alterou mais o desempenho preditivo do que a interpretação econômica central do modelo.

Além dos resultados, deixei no repositório a avaliação fora da amostra, os diagnósticos de convergência e estabilidade, a comparação das respostas ao impulso e um modo rápido para validar a execução.

Projeto completo: https://github.com/karinabelarmino/monetary-policy-under-uncertainty/tree/main/extensions/prior-sensitivity

Em séries macroeconômicas, você costuma avaliar priors por desempenho preditivo, plausibilidade estrutural ou por uma combinação dos dois?

#BayesianEconometrics #Econometrics #TimeSeries #DataScience #RStats
