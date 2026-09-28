# FIXTURE: never tag or merge this branch

`test/untagged-provenance-fixture` holds one leaf commit that no tag contains, by construction.
SWATGenX's `scripts/ship_gate/test_engine_source_provenance.py`
(`test_commit_that_exists_but_is_untagged_fails_naming_zero_tags`) uses it as its example of an
EXISTING-but-UNTAGGED commit (Q615). The test used `main`'s tip until `main` was tagged
(`ship/v1.0-51.1-main`, 2026-09-28), which made its premise false. Tagging or merging this commit breaks that test.
