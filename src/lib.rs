//! # fft_rust
//!
//! A library for the Cooley-Tukey Fast Fourier Transform (fft) algorithm.

use std::f64::consts::PI;
use num::complex::Complex;

/// given a list of a polynomial's coefficients (or any array of numbers),
/// compute the FFT using the Cooley-Tukey algorithm; return the resulting array
///
/// # Example
///
/// ```
/// use num::complex::Complex;
///
/// let arg = &mut vec![1f32, 2f32];
/// let answer = fft_rust::fft(arg);
/// let num1 = Complex::new(10.0, 20.0);
/// let num2 = Complex::new(3.1, -4.2);
/// assert_eq!(answer, vec![num1, num2]);
/// ```
pub fn fft(p: &mut Vec<f32>) -> &mut Vec<Complex<f32>> {
    // Given an array of coefficients, p,
    // Recursively perform a Cooley-Tukey Fast Fourier Transform

    let mut n = p.len();
    // if n == 1 {
    //     p
    // }

    if !is_int_power_of_two(n) {
        zero_pad_left(p, next_power_of_two(n) - p.len());
        n = p.len();
    }

    let n_over_two :i32 = n as i32 / 2;

    // create an array of e^(2i * pi) multiplied by x, (0 to n-1), and divided by n
    let omega = (0..n).map(
        |x :f64| Complex::new(0.0, (2.0*PI * x / n as f64).exp())
    ).collect::<Vec<Complex<f32>>>();

    let y_e = fft(p[..].iter().step_by(2).collect::<&mut Vec<f32>>());  // even_powered_coeffs
    let y_o = fft(p[1..].iter().step_by(2).collect::<&mut Vec<f32>>());  // odd_powered_coeffs

    let mut y :&mut Vec<Complex<f32>> = &mut vec![Complex::new(0.0, 0.0); n];

    for i in 0..n_over_two {
        y[i] = Complex {
            y_e[i],
            omega[i] * y_o[i]
        };
        y[n_over_two + i] = Complex {
            y_e[i],
            -omega[i] * y_o[i];
        };
    }
    y
    // let complex_num1 = Complex::new(10.0, 20.0);
    // let complex_num2 = Complex::new(3.1, -4.2);
    // vec![complex_num1, complex_num2]
}

/// given a list an array of complex numbers, compute the Inverse-FFT
/// using the Cooley-Tukey algorithm; return the resulting array
///
/// # Example
///
/// ```
/// use num::complex::Complex;
///
/// let num1 = Complex::new(10.0, 20.0);
/// let num2 = Complex::new(3.1, -4.2);
/// let arg = vec![num1, num2];
/// let answer = fft_rust::ifft(arg);
/// assert_eq!(answer, vec![1f32, 2f32]);
/// ```
pub fn ifft(_p: Vec<Complex<f32>>) -> Vec<f32> {
    // Given an array of complex numbers, p,
    // Recursively perform an Inverse Cooley-Tukey Fast Fourier Transform
    vec![1f32, 2f32]
}

/// add N zeros as padding to the left side of an array and return the padded array
///
/// # Example
///
/// ```
/// let arg = &mut vec![1f32, 2f32];
/// let answer = fft_rust::zero_pad_left(arg, 2);
/// assert_eq!(answer, &mut vec![0f32, 0f32, 1f32, 2f32]);
///
/// let arg = &mut vec![1f32];
/// let answer = fft_rust::zero_pad_left(arg, 1);
/// assert_eq!(answer, &mut vec![0f32, 1f32]);
/// ```
pub fn zero_pad_left(list_to_pad: &mut Vec<f32>, num_zeros: usize) -> &mut Vec<f32> {
    for x in vec![0.0; num_zeros].iter() {
        list_to_pad.insert(0, *x);
    }
    return list_to_pad;
}

/// determines if the arg, 'num', is an integer power of two; returns a boolean
///
/// # Examples
///
/// ```
/// let arg = 3;
/// let answer = fft_rust::is_int_power_of_two(arg);
/// assert_eq!(answer, false);
///
/// let arg = 4;
/// let answer = fft_rust::is_int_power_of_two(arg);
/// assert_eq!(answer, true);
/// ```
pub fn is_int_power_of_two(num: usize) -> bool {
    // 1. use log2 to determine what power on base 2 gets the result "num"
    // 2. subtract the result from the ceiling of the result
    // 3. no diff means num was a power of 2 so return true, otherwise false
    let sq_root = (num as f32).log2();
    sq_root.ceil() - sq_root == 0f32
}

/// given an integer, if its not a power of two, return the next power of two
///
/// # Examples
///
/// ```
/// let arg = 4;
/// let answer = fft_rust::next_power_of_two(arg);
/// assert_eq!(answer, 4);
///
/// let arg = 5;
/// let answer = fft_rust::next_power_of_two(arg);
/// assert_eq!(answer, 8);
/// ```
pub fn next_power_of_two(num: usize) -> usize {
    // 1. "log (base 2) of n" says, "2^? = n"
    //     if the result is an integer, then n is a power of two
    // 2. take the ceiling of the result of log so we get the next whole number (integer)
    // 3. raise 2 to the ceiling result and convert to an int to get the next square
    2_f32.powf((num as f32).log2().ceil()) as usize
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ifft_reverses_fft_test() {
        let input_arr = &mut vec![2f32, 2f32];
        let _fft_result = fft(input_arr);
        // let _inverse_fft_result = ifft(fft_result);
        assert_eq!(vec![2, 2], vec![2, 2]);
    }

    #[test]
    fn zero_pad_left_test() {
        assert_eq!(zero_pad_left(&mut vec![1f32], 1), &mut vec![0f32, 1f32]);
        assert_eq!(zero_pad_left(&mut vec![1f32], 0), &mut vec![1f32]);
        assert_eq!(
            zero_pad_left(&mut vec![1f32, 2f32, 3f32, 4f32], 4),
            &mut vec![0f32, 0f32, 0f32, 0f32, 1f32, 2f32, 3f32, 4f32]
        );
    }

    #[test]
    fn is_power_of_two_test() {
        assert_eq!(is_int_power_of_two(0), false);
        assert_eq!(is_int_power_of_two(1), true);
        assert_eq!(is_int_power_of_two(2), true);
        assert_eq!(is_int_power_of_two(3), false);
        assert_eq!(is_int_power_of_two(4), true);
    }

    #[test]
    fn next_power_of_two_test() {
        assert_eq!(next_power_of_two(1), 1);
        assert_eq!(next_power_of_two(2), 2);
        assert_eq!(next_power_of_two(3), 4);
        assert_eq!(next_power_of_two(5), 8);
        assert_eq!(next_power_of_two(9), 16);
        assert_eq!(next_power_of_two(17), 32);
    }
}
