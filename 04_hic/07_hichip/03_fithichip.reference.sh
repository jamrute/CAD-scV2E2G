#!/usr/bin/env bash
# =============================================================================
# 04_hic/hichip/03_fithichip.sh — FitHiChIP v11.0, 5 kb, peak-to-all, default params
# Status  : SCAFFOLD — fill fithichip_config_template.txt per cell type
# =============================================================================
set -euo pipefail
FITHICHIP=${FITHICHIP:?path to FitHiChIP repo}
for ct in HCASMC HCAEC; do
  bash "${FITHICHIP}/FitHiChIP_HiCPro.sh" -C "config_${ct}.txt"
done
