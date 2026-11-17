#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# USER CONFIGURATION
# ============================================================
# ------------------------------------------------------------
# MPI
# ------------------------------------------------------------

NPROC=8
NPRINT=1000
NREINIT=1

# ------------------------------------------------------------
# Problem
# ------------------------------------------------------------

DIM=3
SIMULATION="par"
SYMMETRIC=0
AXISYMMETRIC=0

SURFACE_PRESSURE_DIAGNOSTIC="${SURFACE_PRESSURE_DIAGNOSTIC:-0}"
SURFACE_STANDALONE_MATCH="${SURFACE_STANDALONE_MATCH:-0}"

# ------------------------------------------------------------
# Grid configurations
# ------------------------------------------------------------

UNIFORM_LEVELS=(1)
ADAPTIVE_LEVELS=(0)
LEVEL_OFFSET_PSI=0
LEVEL_OFFSET=0

# ------------------------------------------------------------
# Time discretization
# ------------------------------------------------------------

NSTEPS=(1000)

# ------------------------------------------------------------
# Cutfem integration
# ------------------------------------------------------------

CUTFEM_INTEGRATION="quadratic"     # linear | quadratic
SURFACE_WEIGHTS="${SURFACE_WEIGHTS:-3}" # 1: legacy scalar weight; 3: component RHS weights
# Standalone parity is an opt-in diagnostic control; production defaults are unchanged.
if [[ "${SURFACE_STANDALONE_MATCH}" == "1" ]]; then
  if [[ "${SURFACE_WEIGHTS}" == "1" ]]; then DEFAULT_MOMENT_ORDER=17; else DEFAULT_MOMENT_ORDER=10; fi
  DEFAULT_MOMENT_TOL=1e-12
  DEFAULT_FIT_PRECISION=quad
  DEFAULT_GK_RULE=15/31
  DEFAULT_PROJECT_Q1=1
  if [[ "${SURFACE_WEIGHTS}" == "3" ]]; then DEFAULT_NORMAL_FIRST=1; else DEFAULT_NORMAL_FIRST=0; fi
  DEFAULT_UNCONSTRAINED=1
  DEFAULT_STRICT_WARNINGS=1
  DEFAULT_EXTENDED_BUDGET=1
  DEFAULT_DIAGNOSTIC_ONLY=1
else
  DEFAULT_MOMENT_ORDER=10
  DEFAULT_MOMENT_TOL=1e-8
  DEFAULT_FIT_PRECISION=long-double
  DEFAULT_GK_RULE=7/15
  DEFAULT_PROJECT_Q1=0
  DEFAULT_NORMAL_FIRST=0
  DEFAULT_UNCONSTRAINED=0
  DEFAULT_STRICT_WARNINGS=0
  DEFAULT_EXTENDED_BUDGET=0
  DEFAULT_DIAGNOSTIC_ONLY=0
fi
SURFACE_MOMENT_ORDER="${SURFACE_MOMENT_ORDER:-$DEFAULT_MOMENT_ORDER}"
SURFACE_MOMENT_TOL="${SURFACE_MOMENT_TOL:-$DEFAULT_MOMENT_TOL}"
SURFACE_FIT_PRECISION="${SURFACE_FIT_PRECISION:-$DEFAULT_FIT_PRECISION}"
SURFACE_GK_RULE="${SURFACE_GK_RULE:-$DEFAULT_GK_RULE}"
SURFACE_PROJECT_Q1="${SURFACE_PROJECT_Q1:-$DEFAULT_PROJECT_Q1}"
SURFACE_NORMAL_FIRST="${SURFACE_NORMAL_FIRST:-$DEFAULT_NORMAL_FIRST}"
SURFACE_UNCONSTRAINED_PRESSURE="${SURFACE_UNCONSTRAINED_PRESSURE:-$DEFAULT_UNCONSTRAINED}"
SURFACE_STRICT_WARNINGS="${SURFACE_STRICT_WARNINGS:-$DEFAULT_STRICT_WARNINGS}"
SURFACE_EXTENDED_BUDGET="${SURFACE_EXTENDED_BUDGET:-$DEFAULT_EXTENDED_BUDGET}"
SURFACE_PRESSURE_ONLY="${SURFACE_PRESSURE_ONLY:-$DEFAULT_DIAGNOSTIC_ONLY}"

# ------------------------------------------------------------
# Velocity multigrid
# ------------------------------------------------------------

MULTIGRID_TYPE="mg"

# ------------------------------------------------------------
# Executable / options
# ------------------------------------------------------------

EXECUTABLE="./MultiPhase_ex46"
OPTIONS_FILE="./input/FS_solver.options"

# ------------------------------------------------------------
# PETSc configuration
# ------------------------------------------------------------

KSP_MONITOR=false
KSP_VIEW=false
OPTIONS_LEFT=true


# ============================================================
# Basic checks
# ============================================================

SOURCE_SCRIPT_DIR="$(
  dirname "$(readlink -f "${BASH_SOURCE[0]}")"
)"

for variable in SURFACE_STANDALONE_MATCH SURFACE_PROJECT_Q1 SURFACE_NORMAL_FIRST SURFACE_UNCONSTRAINED_PRESSURE SURFACE_STRICT_WARNINGS SURFACE_EXTENDED_BUDGET SURFACE_PRESSURE_ONLY; do
  value="${!variable}"
  if [[ "$value" != 0 && "$value" != 1 ]]; then
    echo "ERROR: $variable must be 0 or 1."
    exit 1
  fi
done

if [[ "${SURFACE_STANDALONE_MATCH}" == "1" && "${SURFACE_PRESSURE_DIAGNOSTIC}" != "1" ]]; then
  echo "ERROR: SURFACE_STANDALONE_MATCH requires SURFACE_PRESSURE_DIAGNOSTIC=1."
  exit 1
fi

if [[ "${SURFACE_FIT_PRECISION}" != "long-double" && "${SURFACE_FIT_PRECISION}" != "quad" ]]; then
  echo "ERROR: SURFACE_FIT_PRECISION must be long-double or quad."
  exit 1
fi

if [[ "${SURFACE_PRESSURE_DIAGNOSTIC}" != "0" &&
      "${SURFACE_PRESSURE_DIAGNOSTIC}" != "1" ]]; then
  echo "ERROR: SURFACE_PRESSURE_DIAGNOSTIC must be 0 or 1."
  exit 1
fi

if [[ "${SURFACE_PRESSURE_DIAGNOSTIC}" == "1" &&
      "${SIMULATION}" != "par" ]]; then
  echo "ERROR: Surface-pressure diagnostic is only for par."
  exit 1
fi

if [[ "${SURFACE_WEIGHTS}" != "1" && "${SURFACE_WEIGHTS}" != "3" ]]; then
  echo "ERROR: SURFACE_WEIGHTS must be 1 or 3."
  exit 1
fi

if [[ "${MULTIGRID_TYPE}" != "hmg" &&
      "${MULTIGRID_TYPE}" != "mg" ]]; then
  echo "ERROR: MULTIGRID_TYPE must be 'hmg' or 'mg'."
  exit 1
fi

if (( NPROC < 1 )); then
  echo "ERROR: NPROC must be >= 1."
  exit 1
fi

if (( LEVEL_OFFSET < 0 )); then
  echo "ERROR: LEVEL_OFFSET must be >= 0."
  exit 1
fi

if (( SYMMETRIC == 1 && AXISYMMETRIC == 1 )); then
  echo "ERROR: SYMMETRIC and AXISYMMETRIC cannot both be 1."
  exit 1
fi

if [[ "${CUTFEM_INTEGRATION}" != "linear" &&
      "${CUTFEM_INTEGRATION}" != "quadratic" ]]; then
  echo "ERROR: CUTFEM_INTEGRATION must be linear or quadratic."
  exit 1
fi


# ============================================================
# Generate FS_solver.options
# ============================================================

generate_options_file() {

  local uniform_levels="$1"
  local adaptive_levels="$2"
  local level_offset="$3"

  local nlevels

  nlevels=$((uniform_levels + adaptive_levels - level_offset))

  if (( nlevels < 1 )); then

    echo "ERROR: Number of solver levels must be >= 1." >&2
    echo "       N = uniform_levels + adaptive_levels - level_offset" >&2
    echo "       N = ${uniform_levels} + ${adaptive_levels} - ${level_offset}" >&2

    return 1
  fi

  mkdir -p "$(dirname "${OPTIONS_FILE}")"

  : > "${OPTIONS_FILE}"


  # ==========================================================
  # Level 0
  # ==========================================================

  cat >> "${OPTIONS_FILE}" <<EOF_OPTIONS
# ============================================================
# Level 0
# ============================================================

# Velocity multigrid type:
#   hmg -> AMG
#   mg  -> GMG

-ksp_type fgmres
-ksp_rtol 1.e-8
-ksp_atol 1.e-10
-ksp_divtol 1.e+50
-ksp_max_it 200

-level-0ksp_max_it 200

-level-0ksp_norm_type unpreconditioned
# -level-0ksp_monitor

-level-0fieldsplit_0_pc_type ${MULTIGRID_TYPE}
-level-0fieldsplit_0_sub_ksp_type preonly
-level-0fieldsplit_0_sub_pc_type lu
-level-0fieldsplit_0_sub_factor_mat_solver_type mumps


-level-0fieldsplit_1_pc_type ${MULTIGRID_TYPE}
-level-0fieldsplit_1_sub_ksp_type preonly
-level-0fieldsplit_1_sub_pc_type lu
-level-0fieldsplit_1_sub_factor_mat_solver_type mumps

EOF_OPTIONS


  # ==========================================================
  # Outer multigrid levels >= 1
  # ==========================================================

  for ((i=1; i<nlevels; ++i)); do

    cat >> "${OPTIONS_FILE}" <<EOF_OPTIONS
# ============================================================
# Level ${i}
# ============================================================

# ------------------------------------------------------------
# FieldSplit level solver
# ------------------------------------------------------------

-level-${i}ksp_norm_type unpreconditioned
# -level-${i}ksp_view_final_residual
# -level-${i}ksp_monitor

-level-${i}ksp_type gmres
-level-${i}ksp_rtol 1.e-8
-level-${i}ksp_atol 1.e-10
-level-${i}ksp_divtol 1.e+50
-level-${i}ksp_max_it 200

-level-${i}pc_fieldsplit_schur_fact_type upper
-level-${i}pc_fieldsplit_schur_precondition selfp

# ------------------------------------------------------------
# Velocity
# ------------------------------------------------------------

-level-${i}fieldsplit_0_ksp_type preonly
-level-${i}fieldsplit_0_pc_type ${MULTIGRID_TYPE}
-level-${i}fieldsplit_0_mg_levels_ksp_type chebyshev
-level-${i}fieldsplit_0_mg_levels_ksp_max_it 4
-level-${i}fieldsplit_0_mg_levels_pc_type jacobi

# ------------------------------------------------------------
# Pressure
# ------------------------------------------------------------

-level-${i}fieldsplit_1_ksp_type preonly
-level-${i}fieldsplit_1_pc_type ${MULTIGRID_TYPE}
-level-${i}fieldsplit_1_mg_levels_ksp_type chebyshev
-level-${i}fieldsplit_1_mg_levels_ksp_max_it 4
-level-${i}fieldsplit_1_mg_levels_pc_type jacobi

EOF_OPTIONS

  done

  printf '%s\n' "${nlevels}"
}


# ============================================================
# PETSc command-line options
# ============================================================

PETSC_OPTIONS=()

if [[ "${KSP_MONITOR}" == true ]]; then
  PETSC_OPTIONS+=(
    -ksp_monitor_true_residual
  )
fi

if [[ "${KSP_VIEW}" == true ]]; then
  PETSC_OPTIONS+=(
    -ksp_view
  )
fi

if [[ "${OPTIONS_LEFT}" == true ]]; then
  PETSC_OPTIONS+=(
    -options_left
  )
fi

if (( NPROC > 1 )); then
  PETSC_OPTIONS+=(
    -matptap_via hypre
  )
fi


# ============================================================
# Number of runs
# ============================================================

TOTAL_RUNS=$(( \
  ${#UNIFORM_LEVELS[@]} * \
  ${#ADAPTIVE_LEVELS[@]} * \
  ${#NSTEPS[@]} \
))

RUN_INDEX=0


# ============================================================
# Run all requested configurations
# ============================================================

for uniform_levels in "${UNIFORM_LEVELS[@]}"; do

  for adaptive_levels in "${ADAPTIVE_LEVELS[@]}"; do

    # --------------------------------------------------------
    # Generate PETSc options for this grid
    # --------------------------------------------------------

    NLEVELS=$(
      generate_options_file \
        "${uniform_levels}" \
        "${adaptive_levels}" \
        "${LEVEL_OFFSET}"
    )


    # --------------------------------------------------------
    # Run all requested time resolutions on this grid
    # --------------------------------------------------------

    for nsteps in "${NSTEPS[@]}"; do

      RUN_INDEX=$((RUN_INDEX + 1))

      echo
      echo "============================================================"
      echo "ex46 run ${RUN_INDEX}/${TOTAL_RUNS}"
      echo "============================================================"
      echo "MPI ranks          : ${NPROC}"
      echo "Uniform levels     : ${uniform_levels}"
      echo "Adaptive levels    : ${adaptive_levels}"
      echo "Level offset       : ${LEVEL_OFFSET}"
      echo "Level offset psi   : ${LEVEL_OFFSET_PSI}"
      echo "CutFEM integration : ${CUTFEM_INTEGRATION}"
      echo "Surface RHS weights: ${SURFACE_WEIGHTS}"
      echo "Standalone parity  : ${SURFACE_STANDALONE_MATCH}"
      echo "Moment order/tol   : ${SURFACE_MOMENT_ORDER} / ${SURFACE_MOMENT_TOL}"
      echo "Fitting/GK         : ${SURFACE_FIT_PRECISION} / ${SURFACE_GK_RULE}"
      echo "Project Q1         : ${SURFACE_PROJECT_Q1}"
      echo "Normal first       : ${SURFACE_NORMAL_FIRST}"
      echo "Unconstrained P    : ${SURFACE_UNCONSTRAINED_PRESSURE}"
      echo "Solver levels      : ${NLEVELS}"
      echo "Number of steps    : ${nsteps}"
      echo "Dimension          : ${DIM}"
      echo "Simulation         : ${SIMULATION}"
      echo "Velocity MG        : ${MULTIGRID_TYPE}"
      echo "Options file       : ${OPTIONS_FILE}"
      echo "============================================================"
      echo

      DIAG_PREFIX="levelsetdiagnostic/surface_pressure/par_U${uniform_levels}_A${adaptive_levels}_NS${nsteps}_W${SURFACE_WEIGHTS}_RUN${RUN_INDEX}"

      if [[ "${SURFACE_PRESSURE_DIAGNOSTIC}" == "1" ]]; then
        mkdir -p "$(dirname "${DIAG_PREFIX}")"
      fi

      if [[ "${SURFACE_EXTENDED_BUDGET}" == "1" ]]; then
        export FEMUS_CUTFEM_UNLIMITED=1
      else
        unset FEMUS_CUTFEM_UNLIMITED
      fi

      mpirun -n "${NPROC}" \
        "${EXECUTABLE}" \
        "${PETSC_OPTIONS[@]}" \
        -options_file "${OPTIONS_FILE}" \
        --uniform-levels "${uniform_levels}" \
        --adaptive-levels "${adaptive_levels}" \
        --level-offset "${LEVEL_OFFSET}" \
        --level-offset-psi "${LEVEL_OFFSET_PSI}" \
        --cutfem-integration "${CUTFEM_INTEGRATION}" \
        --surface-moment-order "${SURFACE_MOMENT_ORDER}" \
        --surface-weights "${SURFACE_WEIGHTS}" \
        --surface-moment-tolerance "${SURFACE_MOMENT_TOL}" \
        --surface-fit-precision "${SURFACE_FIT_PRECISION}" \
        --gk-regular-2d "${SURFACE_GK_RULE}" \
        --gk-tangent-2d "${SURFACE_GK_RULE}" \
        --gk-regular-3d "${SURFACE_GK_RULE}" \
        --gk-tangent-3d "${SURFACE_GK_RULE}" \
        --surface-pressure-project-q1 "${SURFACE_PROJECT_Q1}" \
        --surface-pressure-normal-first "${SURFACE_NORMAL_FIRST}" \
        --surface-pressure-unconstrained "${SURFACE_UNCONSTRAINED_PRESSURE}" \
        --surface-pressure-strict-warnings "${SURFACE_STRICT_WARNINGS}" \
        --surface-pressure-only "${SURFACE_PRESSURE_ONLY}" \
        --nsteps "${nsteps}" \
        --dim "${DIM}" \
        --simulation "${SIMULATION}" \
        --symmetric "${SYMMETRIC}" \
        --axisymmetric "${AXISYMMETRIC}" \
        --nprint "${NPRINT}" \
        --nreinit "${NREINIT}" \
        --surface-pressure-diagnostic "${SURFACE_PRESSURE_DIAGNOSTIC}" \
        --surface-pressure-prefix "${DIAG_PREFIX}" \

      if [[ "${SURFACE_PRESSURE_DIAGNOSTIC}" == "1" ]]; then
        ANALYSIS_ARGS=()
        if [[ "${SURFACE_UNCONSTRAINED_PRESSURE}" == "1" ]]; then
          ANALYSIS_ARGS+=(--unconstrained-pressure)
        fi
        OPENBLAS_NUM_THREADS=1 python3 \
          "${SOURCE_SCRIPT_DIR}/diagnostic/solve_surface_pressure_range.py" \
          "${DIAG_PREFIX}" "${ANALYSIS_ARGS[@]}"
      fi

      echo

    done
  done
done
