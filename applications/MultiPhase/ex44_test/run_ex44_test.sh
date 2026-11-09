#!/usr/bin/env bash
set -euo pipefail

# ------------------ EDIT ONLY THIS SECTION ------------------
ELEMENTS=(0 1 3 4)       # 0=HEX, 1=TET, 3=QUAD, 4=TRI
PSI="ball"          # ball or nonconvex
CURVED=1              # 0=straight, 1=analytical curved BALL
LEVELS=2              # 0 through LEVELS
CURVATURE=0.40         # only used for curved BALL
NONCONVEX_TRI_TAU=0.09  # TRI7 shape (0,0.25)
NONCONVEX_TET_TAU=0.09  # TET15 shape (0,0.25)
NONCONVEX_HEX_A=4       # HEX27 coefficients in [-1,1]^3 reference coordinates
NONCONVEX_HEX_B=4
NONCONVEX_HEX_C=160     # nonconvex if c/(a*b)>3
NONCONVEX_HEX_DIGITS=50 # high-precision HEX analytical surface reference
NONCONVEX_QUAD_A=4      # QUAD9 x^2 coefficient (strict interior: >4)
NONCONVEX_QUAD_B=4      # QUAD9 y^2 coefficient (strict interior: >4)
NONCONVEX_QUAD_C=160    # QUAD9 x^2*y^2 coefficient
MODE="moment"          # moment, algoim, compare
DEGREE=10
ALGOIM_ORDER=9
ALGOIM_SUBDIVISIONS=0
VTK=1
EXECUTABLE="./MultiPhase_ex44_test"
# ------------------------------------------------------------

cd "$(dirname "${BASH_SOURCE[0]}")"
element_names=()
for element in "${ELEMENTS[@]}"; do
  case "$element" in
    0) element_names+=(hex) ;;
    1) element_names+=(tet) ;;
    3) element_names+=(quad) ;;
    4) element_names+=(tri) ;;
    *) echo "Unsupported element ID: $element" >&2; exit 1 ;;
  esac
done
if ((${#element_names[@]} == 0)); then
  echo "ELEMENTS must contain at least one element ID" >&2
  exit 1
fi
elements_csv="$(IFS=,; echo "${element_names[*]}")"

cmd=("$EXECUTABLE"
  --elements "$elements_csv"
  --psi "$PSI"
  --curved "$CURVED"
  --levels "$LEVELS"
  --curvature "$CURVATURE"
  --mode "$MODE"
  --degree "$DEGREE"
  --algoim-order "$ALGOIM_ORDER"
  --algoim-subdivisions "$ALGOIM_SUBDIVISIONS"
  --vtk "$VTK")

if [[ ",$PSI," == *",nonconvex,"* ]]; then
  cmd+=(--nonconvex-quad-a "$NONCONVEX_QUAD_A"
        --nonconvex-quad-b "$NONCONVEX_QUAD_B"
        --nonconvex-quad-c "$NONCONVEX_QUAD_C"
        --nonconvex-tri-tau "$NONCONVEX_TRI_TAU"
        --nonconvex-tet-tau "$NONCONVEX_TET_TAU"
        --nonconvex-hex-a "$NONCONVEX_HEX_A"
        --nonconvex-hex-b "$NONCONVEX_HEX_B"
        --nonconvex-hex-c "$NONCONVEX_HEX_C"
        --nonconvex-hex-digits "$NONCONVEX_HEX_DIGITS")
fi

"${cmd[@]}"
