# Moyo

🚧 This project is under construction 🚧

## Goals
-  Find symmetry operations of a given crystal structure, identify its crystallographic group, and symmetrize the given crystal structure
- Well-defined tolerance for finding symmetry operations
- No dependency on existing symmetry-finder packages
- Simplify crystal symmetry algorithm by extensively exploiting the group structures of crystallographic groups

## Non-goals
- Crystallographic groups in other than three dimensions
- Matching two similar crystal structures

## Details

### Module dependency

```
math <- base <- data <- identify <- standardize <- lib
        ^---- search <--------------|
```

## Acknowledgments

We thank Dr. Yusuke Seto for providing the crystallographic database.
