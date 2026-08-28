"""
Declarative registry of external test/validation datasets hosted on Zenodo.

Each entry in `TEST_DATA_RECORDS` describes one Zenodo record and the files
within it that should be fetched and where they extract to, relative to the
repository root. This is the single place to edit when adding a new dataset
-- no other code needs to change. See `mmgbsa/fetch_test_data.py` for the
Pooch-based download/verification logic that consumes this list.

To add a new dataset: append a new dict to `TEST_DATA_RECORDS` with a unique
`name`, the Zenodo `doi`, and a `files` mapping of {filename: {sha256,
extract_to}}. Get each file's sha256 with `sha256sum <file>` before
uploading, or from Zenodo's own file listing API
(`GET /api/records/<id>/files`, which reports an md5 checksum -- if only
md5 is available from Zenodo, compute sha256 yourself from a local copy and
record that instead, since Pooch verifies against whatever hash type you
give it as long as it's consistently sha256 across this file).

Scope note: deliberately kept to ONE representative fixture per
input-format/binding-mode combination rather than every historical test
variant -- e.g. one Amber-native (.prmtop), one GROMACS-native (.tpr), and
one coordinate-mode (plain PDB) protein-ligand system, plus one
protein-protein/peptide (binding_mode=ppi) system, rather than duplicating
near-identical fixtures. This keeps the download small and each fixture's
purpose unambiguous. New *combinations* (e.g. nucleic acid systems, once
supported) get their own new entries here, not variations on existing ones.
"""

TEST_DATA_RECORDS = [
    {
        "name": "core-fixtures",
        "doi": "10.5281/zenodo.22139199",
        "description": (
            "Core unit-test fixtures used by test/unit/: representative "
            "protein-ligand systems per native input format/loading path "
            "(6T1H, Amber .prmtop-native, including larger multi-ligand "
            "variants; 7KHZ, GROMACS .tpr/.gro-native; 6XJ3, plain-PDB "
            "coordinate-mode; 3IF6, a dimer complex), plus one "
            "protein-protein/peptide (binding_mode=ppi) NAMD/CHARMM .psf "
            "fixture (originally sourced from Zenodo 7186684, see "
            "test/data/charmm_namd_test's own provenance notes)."
        ),
        "files": {
            "6t1h_prmtop_test.tar.gz": {
                "sha256": "afe6337f9327246f8008fcbcc649965483c7981ad45d5d2c897cee1fbac41740",
                "extract_to": "test/data",
            },
            "6t1h_1151_comp.tar.gz": {
                "sha256": "b4f4cc182d4e0fb419ca49eff80043c095da29136c449a1f58aff9c74b6778be",
                "extract_to": "test/data",
            },
            "6t1h_6466_comp.tar.gz": {
                "sha256": "64d311a8ac69577947ca49a5db51a90e5b93ed67492f4be0047212d3f5aff18f",
                "extract_to": "test/data",
            },
            "7khz_tpr_test.tar.gz": {
                "sha256": "9603513345d62aca0e6101190a8ba992a35b0e8f482b931426e5a6799888a3be",
                "extract_to": "test/data",
            },
            "7khz_gro_test.tar.gz": {
                "sha256": "8cfb20fdba10b0b7dad95b60d0a9ac5a325499c0f7cb62f114d25d02f01e52e4",
                "extract_to": "test/data",
            },
            "6xj3_pdb_test.tar.gz": {
                "sha256": "d928e3a7ef3f2e783c475b541a41aba2e5cbe883a94a049fa372e478f0bcc9af",
                "extract_to": "test/data",
            },
            "3if6_dimer_test.tar.gz": {
                "sha256": "19741f00ba8c0a8f84ebe76e09161f69d57844a96c5eddb406a7dc1d29d4eeca",
                "extract_to": "test/data",
            },
            "charmm_namd_test.tar.gz": {
                "sha256": "dbdc5614a240858de24211077b9b8d75ee4a38e068aae42e9f03c2a4dd2debf1",
                "extract_to": "test/data",
            },
        },
    },
    # Future datasets (new validation systems, additional force-field
    # coverage, etc.) get appended here as new dict entries -- each is
    # independent, so partial fetches (`fetch_all(only=[...])`) and partial
    # additions to this list never require touching existing entries.
]
