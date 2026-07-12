#!/bin/bash
# WISP auto inversion for the 2024-07-19 Mw 7.4 San Pedro de Atacama event (us7000n05d).
# Mirrors the 2026 Calama run: header fix, then ffm auto model on both nodal planes.
source ~/miniforge3/etc/profile.d/conda.sh
conda activate ff-env
cd /home/msseo/works/17.Venezuela_2026/event_2024_spda || exit 1
for f in SAC_PZs_*.sac; do [ -e "$f" ] && mv "$f" "${f%.sac}"; done
python ../scripts/wisp_fix_headers.py . us7000n05d_cmt_CMT > wisp_prep.log 2>&1 || { echo "PREP FAILED" >> wisp_auto.log; exit 1; }
ffm model run . auto_model -g us7000n05d_cmt_CMT -d . -t body -t surf > wisp_auto.log 2>&1
echo WISP2024DONE >> wisp_auto.log
