def extra_checks(dataset_dir, params) -> list[str]:
    if params.get("mode") not in ("iter_first", "original"):
        return [f"grampa: params.mode must be iter_first or original, got {params.get('mode')!r}"]
    return []
