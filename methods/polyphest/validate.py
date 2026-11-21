def extra_checks(dataset_dir, params) -> list[str]:
    if "filter_strategy" not in params:
        return ["polyphest: params.filter_strategy is required (see methods/polyphest/examples/)"]
    return []
