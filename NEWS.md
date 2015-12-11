# brainGraph 0.49.0

## Minor changes
* Added argument checking for most functions
* `count_interlobar` no longer takes `atlas.dt` as an argument
* `dti_create_mats` now accepts argument `P` for "number of samples"
* `permute.group` no longer takes `atlas.dt` as an argument
* `vertex_attr_dt` is now essentially a wrapper for `igraph`'s function
    `as_data_frame`
