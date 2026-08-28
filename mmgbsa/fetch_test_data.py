"""
Downloads and verifies external test/validation datasets from Zenodo.

Mirrors the pattern used by GROMACS's `regressiontests` (a separate,
version-tagged package fetched on demand) and Amber's own test data
distribution, but via Pooch (https://www.fatiando.org/pooch/) -- the
standard Python-ecosystem tool for this (used by SciPy, scikit-image,
napari, MNE) -- rather than a bespoke downloader, since Pooch already
handles Zenodo DOI resolution, SHA256 verification, resumable downloads,
and archive extraction.

This module is intentionally decoupled from any specific dataset: all
dataset/file declarations live in `mmgbsa/test_data_manifest.py`. Adding a
new Zenodo record later means editing only that file.
"""
import logging
from pathlib import Path

log = logging.getLogger(__name__)

# Repository root: mmgbsa/fetch_test_data.py -> mmgbsa/ -> repo root
REPO_ROOT = Path(__file__).resolve().parent.parent


def _get_manifest():
    from .test_data_manifest import TEST_DATA_RECORDS
    return TEST_DATA_RECORDS


def _build_pooch(record):
    """Builds a `pooch.Pooch` instance for one manifest record."""
    import pooch

    registry = {
        fname: info["sha256"] for fname, info in record["files"].items()
    }
    missing_hashes = [fname for fname, sha in registry.items() if not sha]
    if missing_hashes:
        raise ValueError(
            f"Dataset '{record['name']}' has files with no sha256 recorded yet "
            f"in mmgbsa/test_data_manifest.py: {missing_hashes}. "
            f"This dataset cannot be fetched until those are filled in."
        )

    return pooch.create(
        path=pooch.os_cache("opengbsa") / record["name"],
        base_url=f"doi:{record['doi']}/",
        registry=registry,
    )


def list_datasets():
    """
    Returns a list of (name, description, n_files, all_hashes_present) for
    every dataset declared in the manifest -- used by the CLI's
    `--list-test-data` to show what's available without downloading anything.
    """
    results = []
    for record in _get_manifest():
        all_present = all(info["sha256"] for info in record["files"].values())
        results.append((
            record["name"],
            record["description"],
            len(record["files"]),
            all_present,
        ))
    return results


def fetch_all(only=None, force=False):
    """
    Downloads (if missing or checksum-stale) and extracts every dataset
    declared in the manifest, or a subset of them.

    Parameters
    ----------
    only : list of str, optional
        Dataset `name`s (from the manifest) to fetch. If None, fetches
        every declared dataset. Use this to avoid pulling the full ~10 GB
        set when only one system's fixtures are needed.
    force : bool, default False
        If True, re-downloads even if a locally cached copy with a matching
        hash already exists (Pooch's own default behavior already skips
        re-downloading valid cached files, so this is rarely needed --
        mainly useful if a fixture was manually deleted/corrupted).

    Returns
    -------
    dict : {dataset_name: {filename: extracted_path_str}}
    """
    import pooch as pooch_lib

    records = _get_manifest()
    if only:
        wanted = set(only)
        records = [r for r in records if r["name"] in wanted]
        found_names = {r["name"] for r in records}
        missing = wanted - found_names
        if missing:
            raise ValueError(f"Unknown dataset name(s) in `only`: {sorted(missing)}. "
                              f"Available: {[r['name'] for r in _get_manifest()]}")

    results = {}
    for record in records:
        log.info(f"Fetching dataset '{record['name']}' ({record['doi']})...")
        p = _build_pooch(record)
        dataset_results = {}
        for fname, info in record["files"].items():
            is_archive = fname.endswith((".tar.gz", ".tgz", ".zip"))
            processor = pooch_lib.Untar() if fname.endswith((".tar.gz", ".tgz")) \
                else pooch_lib.Unzip() if fname.endswith(".zip") \
                else None
            # Pooch's own `extract_dir` on a processor is relative to ITS
            # cache location, not an arbitrary repo path (confirmed via
            # pooch.processors.ExtractorProcessor's own docstring: "extract_dir
            # ... is interpreted as a relative path (relative to the cache
            # location..."). So extraction happens inside Pooch's managed
            # cache (e.g. ~/.cache/opengbsa/<dataset>/), and we symlink the
            # extracted result into the repo's expected `test/data/<name>`
            # location afterward -- this keeps Pooch in charge of avoiding
            # redundant re-downloads across environments/repo clones while
            # still making fixtures show up where existing tests expect them.
            fetch_kwargs = {"progressbar": True}
            if processor is not None:
                fetch_kwargs["processor"] = processor
            if force:
                # Pooch has no direct "force" kwarg on fetch(); the documented
                # way to force a re-fetch is to remove the cached file first.
                cached_path = Path(p.path) / fname
                if cached_path.exists():
                    cached_path.unlink()
            fetched = p.fetch(fname, **fetch_kwargs)
            dataset_results[fname] = fetched

            if is_archive:
                # `fetched` is a list of extracted file paths, all under
                # "<cache>/<archive_name>.untar/<top-level-dir>/..." --
                # confirmed directly (Pooch's Untar extracts into a
                # "{fname}.untar" folder, one level above the archive's own
                # top-level directory). Every archive built by
                # package_test_data.sh tars exactly one named top-level
                # directory (e.g. "6t1h_1151_comp/"), so that directory --
                # two levels up from any extracted file -- is what we
                # symlink into the repo.
                untar_dir = Path(fetched[0])
                while untar_dir.name and not untar_dir.name.endswith(".untar"):
                    untar_dir = untar_dir.parent
                if not untar_dir.name:
                    raise RuntimeError(
                        f"Could not locate the '.untar' extraction root for {fname} "
                        f"among its extracted files: {fetched[:3]}..."
                    )
                # The single top-level directory inside the archive.
                top_level_dirs = {Path(f).relative_to(untar_dir).parts[0] for f in fetched}
                if len(top_level_dirs) != 1:
                    raise RuntimeError(
                        f"Expected exactly one top-level directory in {fname}, found "
                        f"{sorted(top_level_dirs)} -- package_test_data.sh should tar a "
                        f"single named directory per archive."
                    )
                extracted_root = untar_dir / next(iter(top_level_dirs))
                link_target = REPO_ROOT / info["extract_to"] / extracted_root.name
                link_target.parent.mkdir(parents=True, exist_ok=True)
                if link_target.is_symlink() or link_target.exists():
                    if link_target.is_symlink():
                        link_target.unlink()
                    else:
                        log.warning(f"{link_target} already exists and is not a symlink -- "
                                    f"leaving it as-is rather than overwriting real data.")
                        continue
                link_target.symlink_to(extracted_root, target_is_directory=True)
                log.info(f"  {fname} -> {extracted_root} (linked at {link_target})")
            else:
                log.info(f"  {fname} -> {fetched}")
        results[record["name"]] = dataset_results

    return results
