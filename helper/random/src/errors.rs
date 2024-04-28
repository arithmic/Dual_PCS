use core::fmt;
// ================================================================================================
/// Defines errors which can occur when drawing values from a random coin.
#[derive(Debug, PartialEq, Eq)]
pub enum RandomCoinError {
    /// A valid element could not be drawn from the field after the specified number of tries.
    FailedToDrawFieldElement(usize),
    /// The required number of integer values could not be drawn from the specified domain after
    /// the specified number of tries.
    FailedToDrawIntegers(usize, usize, usize),
}
impl fmt::Display for RandomCoinError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::FailedToDrawFieldElement(num_tries) => {
                write!(
                    f,
                    "failed to generate a valid field element after {} tries",
                    num_tries
                )
            }
            Self::FailedToDrawIntegers(num_expected, num_actual, num_tries) => {
                write!(
                    f,
                    "needed to draw {} integers from a domain, but drew only {} after {} tries",
                    num_expected, num_actual, num_tries
                )
            }
        }
    }
}
