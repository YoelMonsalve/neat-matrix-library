# CHANGELOG

## Current version (?)

### [2022-04-21]

- [X] In `nml_mat_free()`: adding sentences to set all pointers to `NULL`
- [ ] Adding documentation in a 'JavaDoc' style
- [ ] We don't really need a flag `is_square`, also this could change if we plan to add methods for resizing a matrix. It is maybe preferable to have a method `square()`, returning a boolean value.
- [ ] Do we have a type `boolean` incorporated ?

### [2022-04-25]
- [ ] Adding symbolic constants for `R_SUCCESS` (0), `R_FAILURE` (1). Using these constants as macros will make the code more flexible and modifiable.
- [ ] Adding symbolic constants for `TRUE`, `True`, `true`, `FALSE`, `False`,`false`.