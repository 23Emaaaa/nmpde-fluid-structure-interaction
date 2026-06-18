#!/usr/bin/env bash

set -eu

if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
  echo "Usage: bash scripts/compare_solution_runs.sh RUN_A.solution.csv RUN_B.solution.csv [REPORT.txt]" >&2
  exit 1
fi

run_a="$1"
run_b="$2"
report_path="${3:-}"

awk -F, '
function abs(x) { return x < 0 ? -x : x }
function sqrt_or_zero(x) { return x > 0 ? sqrt(x) : 0 }

BEGIN {
  tol = 1e-12
}

NR == FNR {
  if (FNR == 1) next
  row = FNR
  material[row] = $1 + 0
  x[row] = $2 + 0
  y[row] = $3 + 0
  measure[row] = $4 + 0
  vx_a[row] = $5 + 0
  vy_a[row] = $6 + 0
  p_a[row] = $7 + 0
  dx_a[row] = $8 + 0
  dy_a[row] = $9 + 0
  next
}

FNR == 1 { next }

{
  row = FNR
  if (!(row in material)) {
    print "Mismatch: second file has more rows than the first one." > "/dev/stderr"
    exit 1
  }

  if (($1 + 0) != material[row]) {
    printf("Mismatch at row %d: material %s != %s\n", row, $1, material[row]) > "/dev/stderr"
    exit 1
  }

  if (abs(($2 + 0) - x[row]) > tol || abs(($3 + 0) - y[row]) > tol) {
    printf("Mismatch at row %d: cell center (%s,%s) != (%s,%s)\n",
           row, $2, $3, x[row], y[row]) > "/dev/stderr"
    exit 1
  }

  if (abs(($4 + 0) - measure[row]) > tol) {
    printf("Mismatch at row %d: measure %s != %s\n", row, $4, measure[row]) > "/dev/stderr"
    exit 1
  }

  if (material[row] == 12) {
    dvx = ($5 + 0) - vx_a[row]
    dvy = ($6 + 0) - vy_a[row]
    dp = ($7 + 0) - p_a[row]

    dvel_norm = sqrt(dvx * dvx + dvy * dvy)
    vel_diff_sq += measure[row] * dvel_norm * dvel_norm
    vel_ref_sq += measure[row] * (vx_a[row] * vx_a[row] + vy_a[row] * vy_a[row])
    if (dvel_norm > vel_linf)
      vel_linf = dvel_norm

    pres_diff_sq += measure[row] * dp * dp
    pres_ref_sq += measure[row] * p_a[row] * p_a[row]
    if (abs(dp) > pres_linf)
      pres_linf = abs(dp)

    total_diff_sq += measure[row] * (dvel_norm * dvel_norm + dp * dp)
    total_ref_sq += measure[row] * (vx_a[row] * vx_a[row] + vy_a[row] * vy_a[row] + p_a[row] * p_a[row])
    if (dvel_norm > total_linf)
      total_linf = dvel_norm
    if (abs(dp) > total_linf)
      total_linf = abs(dp)
  }
  else if (material[row] == 13) {
    ddx = ($8 + 0) - dx_a[row]
    ddy = ($9 + 0) - dy_a[row]

    ddisp_norm = sqrt(ddx * ddx + ddy * ddy)
    disp_diff_sq += measure[row] * ddisp_norm * ddisp_norm
    disp_ref_sq += measure[row] * (dx_a[row] * dx_a[row] + dy_a[row] * dy_a[row])
    if (ddisp_norm > disp_linf)
      disp_linf = ddisp_norm

    total_diff_sq += measure[row] * ddisp_norm * ddisp_norm
    total_ref_sq += measure[row] * (dx_a[row] * dx_a[row] + dy_a[row] * dy_a[row])
    if (ddisp_norm > total_linf)
      total_linf = ddisp_norm
  }

  seen[row] = 1
}

END {
  for (row in material)
    if (!(row in seen)) {
      print "Mismatch: second file has fewer rows than the first one." > "/dev/stderr"
      exit 1
    }

  print "Consistency report"
  print "reference_file=" ARGV[1]
  print "comparison_file=" ARGV[2]
  printf("velocity_l2_error=%.17g\n", sqrt_or_zero(vel_diff_sq))
  printf("velocity_linf_error=%.17g\n", vel_linf + 0)
  printf("velocity_relative_l2_error=%.17g\n",
         vel_ref_sq > 0 ? sqrt_or_zero(vel_diff_sq) / sqrt(vel_ref_sq) : 0)
  printf("pressure_l2_error=%.17g\n", sqrt_or_zero(pres_diff_sq))
  printf("pressure_linf_error=%.17g\n", pres_linf + 0)
  printf("pressure_relative_l2_error=%.17g\n",
         pres_ref_sq > 0 ? sqrt_or_zero(pres_diff_sq) / sqrt(pres_ref_sq) : 0)
  printf("displacement_l2_error=%.17g\n", sqrt_or_zero(disp_diff_sq))
  printf("displacement_linf_error=%.17g\n", disp_linf + 0)
  printf("displacement_relative_l2_error=%.17g\n",
         disp_ref_sq > 0 ? sqrt_or_zero(disp_diff_sq) / sqrt(disp_ref_sq) : 0)
  printf("total_l2_error=%.17g\n", sqrt_or_zero(total_diff_sq))
  printf("total_linf_error=%.17g\n", total_linf + 0)
  printf("total_relative_l2_error=%.17g\n",
         total_ref_sq > 0 ? sqrt_or_zero(total_diff_sq) / sqrt(total_ref_sq) : 0)
}
' "$run_a" "$run_b" | {
  if [ -n "$report_path" ]; then
    tee "$report_path"
  else
    cat
  fi
}
