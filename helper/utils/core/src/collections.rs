/******************************************************************************
 * Copyright (c) 2022 FOLIUM LABS PRIVATE LIMITED. and its affiliates.        *
 ******************************************************************************/
//! Feature-based re-export of common collection components.
//!
//! When `std` feature is enabled, this module exports collections from the Rust standard library.
//! When `alloc` feature is enabled, same collected are provided without relying on the Rust
//! standard library.
#[cfg(not(feature = "std"))]
pub use alloc::collections::{BTreeMap, BTreeSet};
#[cfg(not(feature = "std"))]
pub use alloc::vec::{self as vec, Vec};
#[cfg(feature = "std")]
pub use std::collections::{BTreeMap, BTreeSet};
#[cfg(feature = "std")]
pub use std::vec::{self as vec, Vec};
