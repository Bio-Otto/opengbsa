#!/usr/bin/env bash
# Runs both arms (OpenGBSA and Amber MMPBSA.py) over all five VRC01 systems,
# strictly sequentially so the two never contend for CPU.
#
# These are VALIDATION runs, not timing runs, so OpenGBSA is allowed its normal
# multi-worker CPU path -- verified to give bit-identical results to the
# single-worker path on this system (10 frames, both -47.72 +/- 5.19 kcal/mol).
# Pinning to one worker and one OpenMM thread costs ~34 s/frame versus
# ~6.7 s/frame unpinned, i.e. 14 h instead of 3 h for five systems, for no
# change in the answer.
set -uo pipefail
cd "$(dirname "$0")"

SYSTEMS=(wt g54a d99a n100aa w100ba)
# AMBERHOME: an AmberTools installation (e.g. a conda env with ambertools);
# `opengbsa` must be on PATH (e.g. the opengbsa conda env activated).
AMBER="${AMBERHOME:?set AMBERHOME to an AmberTools installation}"
export AMBERHOME="$AMBER"
command -v opengbsa >/dev/null || { echo "opengbsa not on PATH" >&2; exit 1; }

mkdir -p logs

for s in "${SYSTEMS[@]}"; do
    echo "############ $s ############"

    # --- Amber arm -------------------------------------------------------
    out="amber_run/$s/${s}_300_out.dat"
    if grep -q "DELTA TOTAL" "$out" 2>/dev/null; then
        echo "  amber: already done"
    else
        echo "  amber: running..."
        mkdir -p "amber_run/$s"
        cat > "amber_run/$s/mmpbsa_300.in" <<EOF
gp120/VRC01 $s, OBC2 (igb=5), 0.15 M salt, 300 frames
&general
   startframe=1, endframe=300, interval=1,
   verbose=2, keep_files=0,
/
&gb
   igb=5, saltcon=0.15,
/
EOF
        ( cd "amber_run/$s" && rm -f _MMPBSA_* && \
          "$AMBER/bin/MMPBSA.py" -O -i mmpbsa_300.in \
            -cp "../../amber/$s/${s}_complex.prmtop" \
            -rp "../../amber/$s/${s}_receptor.prmtop" \
            -lp "../../amber/$s/${s}_ligand.prmtop" \
            -y "${s}_dry_300.nc" -o "${s}_300_out.dat" ) \
          > "logs/amber_${s}.log" 2>&1
        grep -E "DELTA TOTAL" "$out" 2>/dev/null || echo "  amber FAILED (see logs/amber_${s}.log)"
    fi

    # --- OpenGBSA arm ----------------------------------------------------
    if grep -q "Mean Binding Energy" "logs/opengbsa_${s}.log" 2>/dev/null; then
        echo "  opengbsa: already done"
    else
        echo "  opengbsa: running..."
        opengbsa "vrc01_${s}_config.yaml" > "logs/opengbsa_${s}.log" 2>&1
        grep -iE "Mean binding energy" "logs/opengbsa_${s}.log" | tail -1 \
          || echo "  opengbsa FAILED (see logs/opengbsa_${s}.log)"
    fi
done

echo
echo "############ SUMMARY ############"
printf "%-9s %12s %12s %9s\n" system amber opengbsa diff
for s in "${SYSTEMS[@]}"; do
    a=$(grep -E "^DELTA TOTAL" "amber_run/$s/${s}_300_out.dat" 2>/dev/null | awk '{print $3}')
    o=$(grep -oE "Mean binding energy \(Enthalpy\): -?[0-9.]+" "logs/opengbsa_${s}.log" 2>/dev/null | awk '{print $NF}')
    if [ -n "$a" ] && [ -n "$o" ]; then
        printf "%-9s %12s %12s %9s\n" "$s" "$a" "$o" \
               "$(awk -v x="$a" -v y="$o" 'BEGIN{printf "%.2f", (x-y<0?y-x:x-y)}')"
    else
        printf "%-9s %12s %12s %9s\n" "$s" "${a:-?}" "${o:-?}" "-"
    fi
done
