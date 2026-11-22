# Changelog

## 2.1.0

### Changed

- MetaTrawl 2 writes only to databases it created. Schema-v1 databases
  (MetaTrawl 1.x), unversioned registries (MetaTrawl 0.x), and files that
  2.0.x adopted from 0.x (`legacy-adopted`) are read-only. Write commands stop
  with an error before changing the file. Keep writing to them with the release
  that created them, or start a new database for ZipStrain 2 profiles.
- A file is created as a new v2 database only when it is missing or empty.
  Version metadata is written before the registry tables, so an interrupted
  first open cannot leave a file that looks like a 0.x registry.
- `runs list`, `profiles list`, `profiles remaining`, and `status` open the
  database read-only, so they work on every older database.
- `profiles compact-database` accepts only sources that MetaTrawl 2 can write.

### Removed

- `metatrawl database migrate`. Opening a 0.x registry no longer adopts it as
  schema v2.
- `metatrawl database declare-contract`. Older profiles stay
  `legacy-unknown`; they are not relabelled with settings that cannot be
  verified.

### Fixed

- The pair API reads MetaTrawl 0.x databases, which have no
  `sample_provenance` table, when `allow_unknown_provenance=True` is passed.

## 2.0.1

- Imports ZipStrain 2 gene catalogs whose gene IDs exceed the signed 64-bit
  range, and accepts four-column BED files during matrix builds.
