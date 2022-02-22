use num_traits::{real::Real, FromPrimitive, NumAssign};

/// An·online·calculator·for·both·mean·and·variance.
///
/// References:
/// - https://doi.org/10.1080/00401706.1962.10490022
/// - https://stats.stackexchange.com/a/235151/146964
pub struct Welford<T> {
    n: usize,
    mean: Option<T>,
    msq: T,
}

impl<T> Welford<T>
where
    T: Real + NumAssign + FromPrimitive,
{
    pub fn new() -> Self {
        Self {
            n: 0,
            mean: None,
            msq: T::zero(),
        }
    }

    pub fn push(&mut self, value: T) {
        self.n += 1;
        let count = self.total();

        if self.mean.is_none() {
            self.mean = Some(value);
        }

        // self.mean is Some(T) from here on.
        let delta = value - self.mean.unwrap();

        *self.mean.as_mut().unwrap() += delta / count;

        let delta2 = value - self.mean.unwrap();
        self.msq += delta * delta2;
    }

    fn total(&self) -> T {
        // The only place we could panic.
        T::from_usize(self.n).expect("failed to convert usize to T")
    }

    /// Mean.
    pub fn mean(&self) -> Option<T> {
        self.mean
    }

    /// Sample variance.
    pub fn var(&self) -> Option<T> {
        if self.n > 1 {
            Some(self.msq / (self.total() - T::one()))
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_welford() {
        let mut w = Welford::new();
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
}
