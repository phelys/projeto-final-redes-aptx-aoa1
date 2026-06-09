# Redes gênicas e a perda de APTX (AOA1): reprogramação do transcriptoma e o módulo das ataxias de reparo de DNA

Projeto final da disciplina **Biologia de Sistemas (PPGAB-UFPR/UTFPR)**.
Estudo de **biologia de sistemas** sobre o gene **APTX (aprataxina)** — causador da **Ataxia com Apraxia
Oculomotora tipo 1 (AOA1)** — usando transcriptoma de células humanas **APTX-knockout vs. selvagem**,
e situando o APTX no **módulo funcional das ataxias de reparo de DNA**.

> **Pergunta:** como a perda de APTX (gene da AOA1) reorganiza a rede de co-expressão gênica, e que
> processos (reparo de DNA, resposta imune inata) são afetados? O APTX e seus parceiros formam um
> módulo coerente na rede de interações?

**Tema de disciplina alinhado à área de pesquisa do autor (ataxias).** A tese de doutorado é sobre
**outra** ataxia (SCA3, análise oculomotora por vídeo); aqui exploramos uma faceta molecular
complementar da neurologia das ataxias — a AOA1 e o reparo de DNA.

---

## Dados (GEO)

| Papel | Accession | Material | Desenho | Plataforma |
|---|---|---|---|---|
| Primário | **GSE245766** | **Microglia humana (HMC3)** | **2×2:** genótipo (APTX-KO / WT) × estímulo imune (IS / NS), 3 réplicas → **n=12** | RNA-seq |

- **Humano** (linhagem celular HMC3) — não é modelo murino.
- Desenho **fatorial 2×2** → permite estudar efeito do **genótipo**, do **estímulo** e da **interação**,
  além de **co-expressão diferencial** entre as 4 condições.
- *Caveat:* é **modelo celular** (microglia), não tecido de paciente; n pequeno (3/condição).
- Dado de iPSC com mutação APTX existe, mas é **"sob solicitação"** (não usado).

---

## Estrutura

```
projeto-final-redes-aptx-aoa1/
├── R/
│   ├── 00_setup.R            # pacotes (CRAN + Bioconductor; org.Hs.eg.db, STRINGdb)
│   ├── config.R              # accession, organismo, desenho 2×2, genes-semente
│   ├── 01_download.R         # baixa GSE245766 (GEOquery)
│   ├── 02_preprocess.R       # contagens → metadados (genótipo×estímulo), filtro, vst
│   ├── 03_diff_expression.R  # DESeq2 fatorial: efeito do KO + interação (baseline)
│   ├── 04_wgcna.R            # rede de co-expressão + módulos
│   ├── 05_diff_coexpression.R# CoDiNA: rede APTX-KO × WT
│   ├── 06_enrichment.R       # GO/KEGG (humano): reparo de DNA, imunidade inata
│   ├── 07_classification.R   # eigengenes → classificador KO×WT (AUC, F1, LOO-CV)
│   ├── 08_network_medicine.R # módulo das ataxias de reparo de DNA (STRING/igraph)
│   └── run_all.R             # pipeline em ordem
├── data/{raw,processed}/     # baixados/processados (NÃO versionar)
├── results/{tables,figures}/ # saídas
└── docs/                     # relatório (artigo) + declaração de uso de IA
```

---

## Como executar

```r
source("R/00_setup.R")   # uma vez (instala dependências)
source("R/run_all.R")    # baixa do GEO e roda tudo
```

Requisitos: **R ≥ 4.2**, internet (download do GEO e do STRING na 1ª execução).

> ⚠️ **Antes da análise:** confira o mapeamento amostra→condição (genótipo/estímulo) nos metadados do
> GEO — `02_preprocess.R` traz um `TODO` no ponto exato.

---

## Métodos (alinhados à disciplina)

- **Co-expressão** — WGCNA (módulos, *eigengenes*, *hubs*).
- **Co-expressão diferencial** — CoDiNA (APTX-KO × WT).
- **Expressão diferencial** — DESeq2 fatorial (baseline / método de referência).
- **Medicina de Redes** — STRINGdb + `igraph`: o **módulo de reparo de DNA** (APTX, SETX, PNKP, ATM,
  TDP1, XRCC1, XRCC4, …); proximidade, centralidade e conectividade.
- **Interpretação** — enriquecimento **GO/KEGG** (clusterProfiler, humano).
- **Validação quantitativa** — classificador com **AUC, precisão, recall, F1, matriz de confusão**
  (LOO-CV), comparado a um **baseline de DEGs**.

### Sobre o n pequeno (M << N)
n = 12 (3/condição): foco em **módulos** (não arestas isoladas), *bootstrap*/permutação e regularização.
O desenho fatorial ajuda — há réplicas e contraste controlado.

---

## Por que APTX / AOA1 (relevância)
- **APTX (aprataxina)** repara **quebras de fita simples de DNA** (resolve intermediários abortivos de
  ligação 5′-adenilados); parceira de **XRCC1/XRCC4**.
- Papel **recém-descrito na imunidade inata** (cGAS-STING, RIG-I/MAVS) — eixo novo e rico para redes.
- Pertence ao grupo das **ataxias recessivas de reparo de DNA** (com SETX/AOA2, PNKP/AOA4, ATM/AT),
  um módulo funcional natural para Medicina de Redes.

## Reprodutibilidade
*Seed* fixa (`config.R`); `sessionInfo` salvo em `results/`. Dados brutos não versionados.

## Licença / uso
Projeto acadêmico. Dados de terceiros sob os termos do GEO/STRING e das publicações originais.
