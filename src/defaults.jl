export defaults


"Dictionary of default arguments."
const defaults = Dict{ASCIIString, Any}(
    "out_dir"            => "exprose_out",
    "n_strucs"           => 250,
    "bound_weight"       => 0.3,
    "other_ratio"        => 20.0,
    "out_prefix"         => "out",
    "iters_per_atom"     => 60000,
    "discard_ratio"      => 1.5,
    "box_size"           => 100.0,
    "cyc_step_ratio"     => 50.0,
    "ens_align_cutoff"   => 0.001,
    "align_cycles"       => 5,
    "align_cutoff"       => 0.5,
    "mod_points"         => 120,
    "mod_intra_cutoff"   => 0.0,
    "mod_intra_tolerance"  => 0.1,
    "mod_inter_cutoff"   => 7.0,
    "mod_inter_tolerance"  => 0.1,
    "mod_min_bound_dist" => 0.001,
)
