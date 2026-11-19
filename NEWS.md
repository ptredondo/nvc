# nvc 0.1.0

* Fix integer overflow in `xi.vec()` for large sample sizes by computing with a double-precision sample size.
* Initial private package draft using the revised manuscript functions.
* Export ten public functions; keep KenCoh eigen and core helpers internal.
* Add help pages and smoke tests without changing the supplied estimators.
