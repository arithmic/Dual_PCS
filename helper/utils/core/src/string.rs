/******************************************************************************
 * Copyright (c) 2022 FOLIUM LABS PRIVATE LIMITED. and its affiliates.        *
 ******************************************************************************/
//! Feature-based re-export of common string components.
//!
//! When `std` feature is enabled, this module exports string components from the Rust standard
//! library. When `alloc` feature is enabled, same components are provided without relying on the
//! Rust standard library.
#[cfg(not(feature = "std"))]
pub use alloc::string::{String, ToString};
#[cfg(feature = "std")]
pub use std::string::{String, ToString};
