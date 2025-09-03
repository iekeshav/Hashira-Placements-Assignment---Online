use std::fs::File;
use std::io::{self, Read};
use serde_json::Value;
use num_bigint::{BigInt, Sign};
use num_rational::BigRational;
use num_traits::{One, Zero, Signed};

/// Convert a single character to its numeric value in the given base.
/// Supports digits ('0'–'9') and letters ('a'–'z' or 'A'–'Z').
fn char_to_digit(c: char) -> u32 {
    match c {
        '0'..='9' => (c as u32) - ('0' as u32),
        'a'..='z' => (c as u32) - ('a' as u32) + 10,
        'A'..='Z' => (c as u32) - ('A' as u32) + 10,
        _ => 0,
    }
}

/// Convert a string representation of a number in a given base into a BigInt.
fn parse_base_value(s: &str, base: u32) -> BigInt {
    let mut result = BigInt::zero();
    let b = BigInt::from(base);
    for ch in s.chars() {
        let digit = BigInt::from(char_to_digit(ch));
        result = result * &b + digit;
    }
    result
}

/// Compute the greatest common divisor of two BigInts.
fn gcd(a: &BigInt, b: &BigInt) -> BigInt {
    let mut x = a.clone();
    let mut y = b.clone();
    // Convert to positive values for gcd computation
    if x.is_negative() { x = -x; }
    if y.is_negative() { y = -y; }
    while !y.is_zero() {
        let r = &x % &y;
        x = y;
        y = r;
    }
    x
}

/// Compute the polynomial of degree k-1 that passes through the given k points using Lagrange interpolation.
/// Returns a vector of BigRational coefficients where index corresponds to the power of x (i.e. coeffs[0] is constant).
fn compute_polynomial(points: &[(BigInt, BigInt)]) -> Vec<BigRational> {
    let k = points.len();
    // The resulting polynomial has degree at most k-1
    let mut poly: Vec<BigRational> = vec![BigRational::from_integer(BigInt::zero()); k];
    for (i, (xi, yi)) in points.iter().enumerate() {
        // Numerator polynomial for basis i starts as [1]
        let mut numer: Vec<BigRational> = vec![BigRational::from_integer(BigInt::one())];
        // Denominator for basis i
        let mut denom = BigRational::from_integer(BigInt::one());
        for (j, (xj, _)) in points.iter().enumerate() {
            if i == j { continue; }
            // Multiply numerator polynomial by (x - xj)
            let mut new_numer = vec![BigRational::from_integer(BigInt::zero()); numer.len() + 1];
            for (deg, coeff) in numer.iter().enumerate() {
                // term for x^deg * (-xj)
                let neg_xj = BigRational::from_integer(-xj.clone());
                new_numer[deg] = new_numer[deg].clone() + coeff * &neg_xj;
                // term for x^deg * x
                new_numer[deg + 1] = new_numer[deg + 1].clone() + coeff.clone();
            }
            numer = new_numer;
            // Update denominator
            let diff = xi - xj;
            denom = denom * BigRational::from_integer(diff);
        }
        // Multiply numerator polynomial by y_i and divide by denom
        let yi_rat = BigRational::from_integer(yi.clone());
        for d in 0..numer.len() {
            numer[d] = numer[d].clone() * &yi_rat / denom.clone();
        }
        // Add this basis polynomial to the result
        if numer.len() > poly.len() {
            poly.resize(numer.len(), BigRational::from_integer(BigInt::zero()));
        }
        for d in 0..numer.len() {
            poly[d] = poly[d].clone() + numer[d].clone();
        }
    }
    poly
}

/// Determine if all polynomial coefficients are positive integers (greater than zero) and return those integers.
/// Returns `Some(Vec<BigInt>)` if all coefficients are positive integers; otherwise returns `None`.
fn are_coeffs_positive_integer(coeffs: &[BigRational]) -> Option<Vec<BigInt>> {
    let mut result = Vec::with_capacity(coeffs.len());
    for coeff in coeffs {
        // Check if denominator is 1
        if !coeff.is_integer() {
            return None;
        }
        let value = coeff.to_integer();
        if value.is_negative() || value.is_zero() {
            return None;
        }
        result.push(value);
    }
    Some(result)
}

/// Attempt to find a valid polynomial from any combination of k points among the n given points.
/// A valid polynomial is defined as one with all positive integer coefficients.
/// Returns the coefficients of the valid polynomial if found.
fn find_valid_polynomial(
    entries: &[(BigInt, BigInt)],
    k: usize,
) -> Option<Vec<BigInt>> {
    let n = entries.len();
    let mut comb: Vec<usize> = (0..k).collect();
    loop {
        // Build points for this combination
        let mut pts: Vec<(BigInt, BigInt)> = Vec::with_capacity(k);
        for &idx in &comb {
            pts.push((entries[idx].0.clone(), entries[idx].1.clone()));
        }
        let poly = compute_polynomial(&pts);
        if let Some(coeffs) = are_coeffs_positive_integer(&poly) {
            return Some(coeffs);
        }
        // Generate next combination lexicographically
        let mut i = k;
        let mut found = false;
        while i > 0 {
            i -= 1;
            if comb[i] != i + n - k {
                comb[i] += 1;
                for j in i+1..k {
                    comb[j] = comb[j-1] + 1;
                }
                found = true;
                break;
            }
        }
        if !found {
            break;
        }
    }
    None
}

fn main() {
    // Read JSON input from stdin
    let mut input = String::new();
    io::stdin().read_to_string(&mut input).expect("Failed to read input");
    let v: Value = serde_json::from_str(&input).expect("Invalid JSON input");
    let keys_obj = v.get("keys").expect("Missing 'keys' field");
    let n = keys_obj.get("n").expect("Missing 'n'").as_u64().expect("Invalid 'n'") as usize;
    let k = keys_obj.get("k").expect("Missing 'k'").as_u64().expect("Invalid 'k'") as usize;
    // Collect (x, y) pairs
    let mut entries: Vec<(BigInt, BigInt)> = Vec::with_capacity(n);
    let obj = v.as_object().expect("JSON must be an object");
    for (key, val) in obj.iter() {
        if key == "keys" { continue; }
        // parse key as x coordinate
        let x = BigInt::parse_bytes(key.as_bytes(), 10).expect("Invalid key for x");
        let base_str = val.get("base").expect("Missing base").as_str().expect("Base must be string");
        let base: u32 = base_str.parse().expect("Invalid base");
        let val_str = val.get("value").expect("Missing value").as_str().expect("Value must be string");
        let y = parse_base_value(val_str, base);
        entries.push((x, y));
    }
    // Sort entries by x ascending
    entries.sort_by(|a, b| a.0.cmp(&b.0));
    // If the exact number of entries isn't equal to provided n, you may still continue.
    let n_entries = entries.len();
    assert!(n_entries >= k, "Not enough points to solve");
    // Attempt to find a valid polynomial with all positive integer coefficients by exploring combinations
    let coeffs_opt = find_valid_polynomial(&entries, k);
    if let Some(coeffs) = coeffs_opt {
        // Output the constant term (coefficient for x^0)
        println!("{}", coeffs[0]);
    } else {
        // If no valid polynomial found, fall back to computing constant via Lagrange on first k points
        let fallback_points: Vec<(BigInt, BigInt)> = entries.iter().take(k).map(|(x, y)| (x.clone(), y.clone())).collect();
        let poly = compute_polynomial(&fallback_points);
        // constant is poly[0]
        let const_r = &poly[0];
        if const_r.is_integer() {
            println!("{}", const_r.to_integer());
        } else {
            // If constant is not integer, print fraction as numerator/denominator
            println!("{}/{}", const_r.numer(), const_r.denom());
        }
    }
}
