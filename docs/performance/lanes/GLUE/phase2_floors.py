#!/usr/bin/env python3
"""Analytic work/break-even calculations only. Never executes or times DFT.

72^3/128 occupied bands are the documented Si64 fixture model, not extracted
from the relayed e947a35 results (raw results were not supplied).
"""
import json
import math
from pathlib import Path


def calculate():
    n, occupied, maps = 72**3, 128, 15
    pairs = occupied // 2
    dense_flops = 5*n*math.log2(n)
    # Three axis passes, each reads and writes complex values once. Logical
    # model only: caches, vendor algorithms and transpose scratch differ.
    axis_bytes = 3*32*n
    # Serial split contribution: clear (8), read/write += (16), fold read (8)
    # per band/point. Direct accumulation reads/writes rho (16) instead.
    saved_bytes = 16*occupied*n*maps
    return {
        'scope': 'analytic fixture model; no measured hardware floor or speed claim',
        'grid_points': n, 'occupied_bands': occupied, 'paired_inverse_ffts_per_map': pairs,
        'assumed_density_maps': maps, 'paired_inverse_ffts_15_maps': pairs*maps,
        'dense_complex_fft_flops_per_pair_estimate': dense_flops,
        'dense_fft_flops_15_maps_estimate': dense_flops*pairs*maps,
        'logical_axis_pass_bytes_15_maps': axis_bytes*pairs*maps,
        'kept_psi_write_bytes_15_maps': 16*n*pairs*maps,
        'split_contribution_bytes_15_maps': 32*occupied*n*maps,
        'direct_accumulation_saved_logical_bytes_15_maps': saved_bytes,
        'density_seconds_per_map_normalized': 6.3032/maps,
        'density_seconds_per_pair_normalized': 6.3032/(maps*pairs),
        'dense_fft_effective_gflops_s_if_all_density_time_is_fft': dense_flops*pairs*maps/6.3032/1e9,
        'fusion_effective_bandwidth_gb_s_for_0_8s_break_even': saved_bytes/.8/1e9,
        'density_fraction_to_remove_serial': .8/6.3032,
        'density_fraction_to_remove_mpi8': .2/.7331,
        'energy_fraction_to_remove_serial': .8/.9897,
        'energy_fraction_to_remove_mpi8': .2/.2865,
        'xc_asymmetry_proxy_seconds_if_spread_over_max': .1559*.374,
        'verification_potential_reuse_serial_optimistic_seconds': .7749*2/15,
        'verification_potential_reuse_mpi8_optimistic_seconds': .1226*2/15,
        'omitting_two_verification_density_maps_serial_invalid_ceiling_seconds': 6.3032*2/15,
    }


if __name__ == '__main__':
    path = Path(__file__).with_name('phase2-floors.json')
    path.write_text(json.dumps(calculate(), indent=2, allow_nan=False)+'\n')
    print('Wrote analytic counts and break-even values; no benchmark executed.')
