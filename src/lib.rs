use num_traits::cast;
use num_traits::Num;
use num_traits::NumAssign;
use num_traits::NumCast;
use num_traits::Zero;

/// Online algorithm for mean and variance, with support for uneven weights.
///
/// This implements the [Welford's online algorithm][welford-wiki] for
/// computing mean and variance in a single pass.
/// The initial implementation was based on [this StackOverflow
/// answer](https://stats.stackexchange.com/a/235151/146964) by
/// [Tim](https://stats.stackexchange.com/users/35989/tim).
///
/// Weights are treated as suggested by [West][west-wiki]. The implementation
/// uses weights as [frequencies instead of reliabilities][weighted-variance]
/// for the calculation of variances.
///
/// [welford-wiki]: https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Welford's_online_algorithm
/// [west-wiki]: https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Weighted_incremental_algorithm
/// [weighted-variance]: https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Weighted_sample_variance
///
/// ## References
///
/// - [Welford, B. P. 1962. Technometrics 4 (3): 419–20.](https://doi.org/10.1080/00401706.1962.10490022)
/// - [West, D. H. D. 1979. Communications of the ACM 22 (9): 532–35.](https://doi.org/10.1145/359146.359153)
#[derive(Clone, Copy, Debug)]
pub struct Welford<T, W = usize> {
    mean: Option<T>,
    total: W,
    msq: T,
}

impl<T> Welford<T>
where
    T: Zero,
{
    /// Create a new *unweighted* Welford calculator.
    ///
    /// # Examples
    /// ```
    /// use welford::Welford;
    /// let mut w = Welford::<f32>::new();
    /// ```
    #[must_use]
    pub fn new() -> Self {
        Self {
            mean: None,
            total: 0,
            msq: T::zero(),
        }
    }
}

impl<T> Default for Welford<T>
where
    T: Zero,
{
    /// Alias for `Welford::new()`.
    ///
    /// # Examples
    /// ```
    /// use welford::Welford;
    /// let w = Welford::<f32>::default();
    /// ```
    fn default() -> Self {
        Self::new()
    }
}

impl<T> Welford<T>
where
    T: Copy + Num + NumAssign + NumCast,
{
    /// Add a new sample to the calculator.
    ///
    /// The weight is assumed to be unity.
    ///
    /// # Examples
    /// ```
    /// # use welford::Welford;
    /// let mut w = Welford::new();
    /// w.push(1.0);
    /// ```
    pub fn push(&mut self, value: T) {
        self.push_weighted(value, 1)
    }
}

impl<T, W> Welford<T, W>
where
    T: Zero,
    W: Zero,
{
    /// Create a new *weighted* Welford calculator.
    ///
    /// # Examples
    /// ```
    /// use welford::Welford;
    /// let mut w = Welford::<f32>::with_weights();
    /// ```
    #[must_use]
    pub fn with_weights() -> Self {
        Self {
            mean: None,
            total: W::zero(),
            msq: T::zero(),
        }
    }
}

impl<T, W> Welford<T, W>
where
    T: Copy + Num + NumAssign + NumCast,
    W: Copy + Num + NumAssign + NumCast + PartialOrd,
{
    /// Add a new sample to the calculator with a given weight.
    ///
    /// # Examples
    /// ```
    /// # use welford::Welford;
    /// let mut w = Welford::with_weights();
    /// w.push_weighted(1.0, 2);
    /// ```
    pub fn push_weighted(&mut self, value: T, weight: W) {
        self.total += weight;

        if self.mean.is_none() {
            self.mean = Some(value);
        }

        // self.mean is Some(T) from here on.
        let delta = value - self.mean.unwrap();

        let total = cast(self.total).expect("failed to cast W to T");
        let weighted_delta = delta * cast(weight).expect("failed to cast W to T");

        *self.mean.as_mut().unwrap() += weighted_delta / total;

        let delta2 = value - self.mean.unwrap();
        self.msq += weighted_delta * delta2;
    }

    /// Get the mean of the samples so far.
    ///
    /// # Examples
    /// ```
    /// # use welford::Welford;
    /// let mut w = Welford::new();
    /// w.push(1.0);
    /// w.push(2.0);
    /// assert_eq!(w.mean(), Some(1.5));
    /// ```
    pub fn mean(&self) -> Option<T> {
        self.mean
    }

    /// Get the variance of the samples so far.
    ///
    /// Weights are treated as
    /// [frequencies instead of reliabilities][weighted-variance].
    ///
    /// [weighted-variance]: https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Weighted_sample_variance
    pub fn var(&self) -> Option<T> {
        if self.total > W::one() {
            let total: T = cast(self.total).expect("failed to cast W to T");
            Some(self.msq / (total - T::one()))
        } else {
            None
        }
    }

    /// Merge the contents of another Welford calculator into this one.
    ///
    /// # Examples
    /// ```
    /// # use welford::Welford;
    /// let mut w1 = Welford::new();
    /// w1.push(1.0f32);
    /// w1.push(2.0);
    /// let mut w2 = Welford::new();
    /// w2.push(3.0f32);
    /// w2.push(4.0);
    /// w1.merge(w2);
    /// assert_eq!(w1.mean(), Some(2.5));
    /// ```
    pub fn merge(&mut self, other: Self) {
        let weight = other.total;

        if weight == W::zero() {
            return;
        } else if self.total == W::zero() {
            *self = other;
            return;
        }

        // self.mean is Some(T) from here on since totals have been updated.
        // WARN: Probably unstable, see <https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Parallel_algorithm>.
        let delta = other.mean.unwrap() - self.mean.unwrap();

        let total = self.total + weight;
        let weighted_delta = delta * cast(weight).expect("failed to cast W to T");

        let mean_corr = weighted_delta / cast(total).expect("failed to cast W to T");
        *self.mean.as_mut().unwrap() += mean_corr;

        self.msq +=
            other.msq + delta * cast(self.total).expect("failed to cast W to T") * mean_corr;

        self.total = total;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_welford() {
        let mut w = Welford::default();
        assert_eq!(w.mean(), None);
        assert_eq!(w.var(), None);

        w.push(1.0);
        assert_eq!(w.mean(), Some(1.0));
        assert_eq!(w.var(), None);

        w.push(3.0);
        assert_eq!(w.mean(), Some(2.0));
        assert_eq!(w.var(), Some(2.0));

        w.push(5.0);
        assert_eq!(w.mean(), Some(3.0));
        assert_eq!(w.var(), Some(4.0));
    }

    #[test]
    fn test_weighted_welford() {
        let mut w = Welford::with_weights();
        assert_eq!(w.mean(), None);
        assert_eq!(w.var(), None);

        w.push_weighted(1.0, 3.0);
        assert_eq!(w.mean(), Some(1.0));
        assert_eq!(w.var(), Some(0.0));

        w.push_weighted(3.0, 2.0);
        assert_eq!(w.mean(), Some(1.8));
        assert_eq!(w.var(), Some(1.2));

        w.push_weighted(5.0, 1.0);
        assert_eq!(w.mean(), Some(2.333_333_333_333_333_5));
        assert_eq!(w.var(), Some(2.666_666_666_666_666_5));
    }

    #[test]
    fn test_merge() {
        let mut w1 = Welford::new();
        let mut w2 = Welford::new();

        w1.push(1.0);
        w1.push(3.0);
        w1.push(5.0);
        w1.push(7.0);

        w2.push(2.0);
        w2.push(4.0);
        w2.push(6.0);
        w2.push(8.0);

        w1.merge(w2);
        assert_eq!(w1.mean(), Some(4.5));
        assert_eq!(w1.var(), Some(6.0));
    }

    #[test]
    fn test_weighted_merge() {
        let mut w1 = Welford::with_weights();
        let mut w2 = Welford::with_weights();

        w1.push_weighted(1.0, 4.0);
        w1.push_weighted(3.0, 3.0);
        w1.push_weighted(5.0, 2.0);
        w1.push_weighted(7.0, 1.0);

        w2.push_weighted(2.0, 4.0);
        w2.push_weighted(4.0, 3.0);
        w2.push_weighted(6.0, 2.0);
        w2.push_weighted(8.0, 1.0);

        w1.merge(w2);
        assert_eq!(w1.mean(), Some(3.5));
        assert_eq!(w1.var(), Some(4.473_684_210_526_316));
    }
}
