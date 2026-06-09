# [Título] Reprogramação da rede gênica pela perda de APTX (AOA1): co-expressão, reparo de DNA e imunidade inata

**Autor:** Pedro Henrique Ramos Prado · **Disciplina:** Biologia de Sistemas (PPGAB) · **Data:** 2026

> Esqueleto em formato de artigo (escolher a revista-alvo e ajustar às normas).
> Métricas e figuras são geradas pelo pipeline em `R/` e ficam em `results/`.

---

## Resumo / Abstract
*(Objetivo → dados → métodos → principais resultados/métricas → conclusão. PT + EN.)*
**Palavras-chave:** AOA1; APTX/aprataxina; reparo de DNA; redes de co-expressão; WGCNA; medicina de redes.

## 1. Introdução
- AOA1 e o gene **APTX (aprataxina)**: reparo de **quebras de fita simples de DNA**; parceria com XRCC1/XRCC4.
- Papel emergente do APTX na **imunidade inata** (cGAS-STING, RIG-I/MAVS).
- As **ataxias recessivas de reparo de DNA** como módulo funcional (APTX, SETX, PNKP, ATM, TDP1).
- Biologia de sistemas como lente; objetivo e hipóteses.
- *(Nota de contexto: a tese de doutorado do autor trata de outra ataxia — SCA3, fenótipo oculomotor por vídeo; este trabalho explora uma faceta molecular da neurologia das ataxias.)*

## 2. Materiais e Métodos
- **Dado:** GSE245766 — microglia humana HMC3, **2×2** (APTX-KO/WT × estímulo imune), n=12. Modelo celular (limitação declarada).
- **Pré-processamento:** filtro (`filterByExpr`), normalização (`vst`).
- **Expressão diferencial:** DESeq2 fatorial — efeito do genótipo, do estímulo e **interação** (baseline/referência).
- **Co-expressão:** WGCNA (β por scale-free, módulos, eigengenes; localização do APTX e genes-semente).
- **Co-expressão diferencial:** CoDiNA (APTX-KO × WT).
- **Medicina de Redes:** STRINGdb + igraph — módulo de reparo de DNA, centralidade, integração com DEGs.
- **Interpretação:** enriquecimento GO/KEGG (humano).
- **Validação quantitativa:** classificador (eigengenes vs DEGs), LOO-CV; AUC/precisão/recall/F1.
- **Código:** repositório GitHub [inserir link]; *seed*, `sessionInfo`.

## 3. Resultados e Discussão
- DEGs (efeito do KO; do estímulo; interação) + volcano.
- Módulos WGCNA e associação com genótipo/estímulo; onde cai o APTX.
- Rede diferencial (CoDiNA): genes que mudam de conectividade com a perda de APTX.
- Enriquecimento: reparo de DNA, **resposta imune inata / sensoriamento de DNA citosólico**.
- Medicina de Redes: APTX como hub no módulo de reparo de DNA; quais parceiros são DEGs.
- **Métricas** (tabela): rede × baseline; matriz de confusão; ROC.
- Discussão à luz da biologia da AOA1 (reparo de DNA ↔ neurodegeneração ↔ neuroinflamação).

## 4. Conclusão
- Contribuições; dificuldades (n pequeno, modelo celular); síntese; o que os achados sugerem sobre o papel de rede do APTX.

## 5. Referências
*(APTX/aprataxina e reparo de SSB; APTX e imunidade inata; WGCNA — Langfelder & Horvath 2008; CoDiNA — Gysi et al.; STRINGdb; clusterProfiler — Yu et al.; GSE245766.)*

## 6. Declaração de Uso de Inteligência Artificial
Ver `docs/declaracao-uso-ia.md`.
