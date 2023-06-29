//! # fft_rust
//!
//! A library for the Cooley-Tukey Fast Fourier Transform (fft) algorithm.

// pub fn fft<T>(p: Vec<T>) -> Vec<T> {
//     /// Given an array of coefficients, p,
//     /// Recursively perform a Cooley-Tukey Fast Fourier Transform
//
//     // let n = p.len();
//     // if n == 1 {
//     //     p
//     // }
//
//     // if not is_power_of_two(n):
//     // p = zero_pad_left(p, next_square(n) - len(p))
//     // n = len(p)
//
//     // n_over_two = int(n / 2)
//     //
//     // w = np.exp(2j * pi * np.arange(n) / n)
//     //
//     // p_e, p_o = p[::2], p[1::2]  # even_powered_coeffs, odd_powered_coeffs
//     // y_e, y_o = fft(p_e), fft(p_o)
//     //
//     // y = [0] * n
//     //
//     // for i in range(n_over_two):
//     // y[i] = y_e[i] + (w[i] * y_o[i])
//     // y[n_over_two + i] = y_e[i] - (w[i] * y_o[i])
//     vec![2, 2]
// }

// pub fn ifft<T>(p: Vec<T>) -> Vec<T> {
//     /// Given an array of coefficients, p,
//     /// Recursively perform an Inverse Cooley-Tukey Fast Fourier Transform
//     vec![2, 2]
// }

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
pub fn is_int_power_of_two(num: i32) -> bool {
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
pub fn next_power_of_two(num: i32) -> i32 {
    // 1. "log (base 2) of n" says, "2^? = n"
    //     if the result is an integer, then n is a power of two
    // 2. take the ceiling of the result of log so we get the next whole number (integer)
    // 3. raise 2 to the ceiling result and convert to an int to get the next square
    2_f32.powf((num as f32).log2().ceil()) as i32
}

#[cfg(test)]
mod tests {
    use super::*;

    // #[test]
    // fn ifft_reverses_fft_test() {
    //     let input_arr = vec![2, 2];
    //     let fft_result = fft(input_arr);
    //     assert_eq!(ifft(fft_result), p);
    // }

    #[test]
    fn is_power_of_two_test() {
        assert_eq!(is_int_power_of_two(-1), false);
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
