# welford

Online algorithm for mean and variance, with support for uneven weights.

This implements the [Welford's online algorithm][welford-wiki] for
computing mean and variance in a single pass.
The initial implementation was based on [this StackOverflow
answer](https://stats.stackexchange.com/a/235151/146964) by
[Tim](https://stats.stackexchange.com/users/35989/tim).

Weights are treated as suggested by [West][west-wiki]. The implementation uses
weights as [frequencies instead of reliabilities][weighted-variance] for
the calculation of variances.

[welford-wiki]: https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Welford's_online_algorithm
[west-wiki]: https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Weighted_incremental_algorithm
[weighted-variance]: https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Weighted_sample_variance

## References

- [Welford, B. P. 1962. Technometrics 4 (3): 419–20.](https://doi.org/10.1080/00401706.1962.10490022)
- [West, D. H. D. 1979. Communications of the ACM 22 (9): 532–35.](https://doi.org/10.1145/359146.359153)

## Related projects

-   [`rolling-stats`](https://crates.io/crates/rolling-stats)
