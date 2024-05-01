# CHANGELOG



## v0.2.0 (2024-05-01)

### Chore

* chore: upgrade Codecov action ([`97873f4`](https://github.com/kalekundert/macromol_dataframe/commit/97873f4120823724b681cbe98afa1334e459b111))

### Documentation

* docs: write a README file ([`65606c1`](https://github.com/kalekundert/macromol_dataframe/commit/65606c1dbdf313e3952907874ad1972fa1b7bd3e))

### Feature

* feat: automatically import the pymol module ([`28c4be1`](https://github.com/kalekundert/macromol_dataframe/commit/28c4be1e6ec68693afaf89a0efc01fe2cbafee05))

* feat: parse more details about each biological assembly ([`2a9fb23`](https://github.com/kalekundert/macromol_dataframe/commit/2a9fb231c35e01840a08b9976dee09ca77a5567b))

* feat: parse parenthetical assembly operation expressions ([`7f9a2b2`](https://github.com/kalekundert/macromol_dataframe/commit/7f9a2b28904bd8d9b1c1dab1719f1904c8f35633))

* feat: load atoms from PyMOL selections ([`e676d7b`](https://github.com/kalekundert/macromol_dataframe/commit/e676d7b07b57706d6077f81f3584098ef1fe46c8))

* feat: add a function to generate PDB-style paths ([`2ada4d0`](https://github.com/kalekundert/macromol_dataframe/commit/2ada4d079f0ec35f0e12ebfd2c9e665ee1a6379c))

* feat: automatically clip occupancies to the range [0, 1] ([`5afadda`](https://github.com/kalekundert/macromol_dataframe/commit/5afaddaec20aefa3fc99cd0310d6d383d0bd3a1f))

* feat: expose all dataframes parsed from mmCIF ([`b697c9b`](https://github.com/kalekundert/macromol_dataframe/commit/b697c9b83759b275b66fa697a45bf36e33d8515f))

### Fix

* fix: prune hydrogren and water ([`d5f9e2c`](https://github.com/kalekundert/macromol_dataframe/commit/d5f9e2cfc4b8db6282ec054960942b06a65e675a))

### Test

* test: complain about unsupported operation expressions ([`cb791ff`](https://github.com/kalekundert/macromol_dataframe/commit/cb791ffe370eeaf68f24f142c19618a7eb97dae5))


## v0.1.0 (2024-03-27)

### Chore

* chore: configure automated releases ([`364760a`](https://github.com/kalekundert/macromol_dataframe/commit/364760a292103e1ab036e27eaaaeecd9d34ee424))

* chore: fix lint error ([`367d0bd`](https://github.com/kalekundert/macromol_dataframe/commit/367d0bd5337cf0a8835692c8e76d90b80a37458c))

* chore: require python 3.10 ([`9b16b9f`](https://github.com/kalekundert/macromol_dataframe/commit/9b16b9f6451bcb0588f52ca1ee84121c215c4722))

* chore: apply cookiecutter ([`a23cbc2`](https://github.com/kalekundert/macromol_dataframe/commit/a23cbc294c01f15091dee31a6de5a1c1b6bfa3b7))

### Feature

* feat: add a way to load asymmetric units ([`966fe50`](https://github.com/kalekundert/macromol_dataframe/commit/966fe50af9a864d200e405b1541287910ca4d58d))

* feat: round floats when writing mmCIF files ([`99b0b6f`](https://github.com/kalekundert/macromol_dataframe/commit/99b0b6f5fd88e879a48458533ed07aec3ad58d57))

* feat: initial implementation

Most of the code is copied from atompaint, but the mmCIF parser is new. ([`034227d`](https://github.com/kalekundert/macromol_dataframe/commit/034227dff44a6b477b484dc3f828116a5a70393c))

### Refactor

* refactor: make some test helpers more easily available ([`524bf82`](https://github.com/kalekundert/macromol_dataframe/commit/524bf826d8e8eeedbd8b1fbb7cfb9b2b0e3ab85c))
