//! # fft_rust
//!
//! A library for the Cooley-Tukey Fast Fourier Transform (fft) algorithm.

use std::f64::consts::PI;
use num::complex::Complex;

/// given a list of a polynomial's coefficients (or any array of numbers),
/// compute the FFT using the Cooley-Tukey algorithm; return the resulting array
///
/// # Examples
///
/// ```
/// use num::complex::Complex;
///
/// let arg = &mut vec![Complex::new(2f64, 0f64), Complex::new(2f64, 0f64)];
/// assert_eq!(fft_rust::fft(arg), vec![Complex::new(4f64, 0f64), Complex::new(0f64, 0f64)]);
///
/// let arg = &mut vec![
///     Complex::new(2f64, 0f64), Complex::new(7f64, 0f64),
///     Complex::new(-1f64, 0f64)
/// ];
/// assert_eq!(fft_rust::fft(arg), vec![
///     Complex::new(8f64, 0f64), Complex::new(-7f64, -3f64),
///     Complex::new(6f64, 0f64), Complex::new(-7f64, 3f64)
/// ]);
///
/// let arg = &mut vec![
///     Complex::new(4f64, 0f64), Complex::new(3f64, 0f64),
///     Complex::new(-5f64, 0f64), Complex::new(1f64, 0f64)
/// ];
/// assert_eq!(fft_rust::fft(arg), vec![
///     Complex::new(3f64, 0f64), Complex::new(9f64, -2f64),
///     Complex::new(-5f64, 0f64), Complex::new(9f64, 2f64)
/// ]);
/// ```
pub fn fft(p: &mut Vec<Complex<f64>>) -> Vec<Complex<f64>> {
    // Given an array of coefficients, p,
    // Recursively perform a Cooley-Tukey Fast Fourier Transform

    // if p is length 1, then just return p
    let mut n = p.len();
    if n == 1 {
        return p.iter().map(|x| *x).collect::<Vec<Complex<f64>>>();
    }

    // ensure length of p is a power of two
    if !is_int_power_of_two(n) {
        zero_pad_left(p, next_power_of_two(n) - p.len());
        n = p.len();
    }

    // even_powered_coeffs
    let y_e = fft(
        &mut p[..].iter().step_by(2).map(|x| *x).collect::<Vec<_>>()
    );
    // odd_powered_coeffs
    let y_o = fft(
        &mut p[1..].iter().step_by(2).map(|x| *x).collect::<Vec<_>>()
    );

    butterfly(n, y_e, y_o, false)
}

/// given a list an array of complex numbers, compute the Inverse-FFT
/// using the Cooley-Tukey algorithm; return the resulting array
///
/// # Example
///
/// ```
/// use num::complex::Complex;
///
/// let arg = &mut vec![Complex::new(4f64, 0f64), Complex::new(0f64, 0f64)];
/// assert_eq!(fft_rust::ifft(arg), vec![Complex::new(2f64, 0f64), Complex::new(2f64, 0f64)]);
///
/// let arg = &mut vec![Complex::new(-2f64, 0f64), Complex::new(4f64, 0f64)];
/// assert_eq!(fft_rust::ifft(arg), vec![Complex::new(1f64, 0f64), Complex::new(-3f64, 0f64)]);
///
/// let arg = &mut vec![
///     Complex::new(8f64, 0f64), Complex::new(-7f64, -3f64),
///     Complex::new(6f64, 0f64), Complex::new(-7f64, 3f64)
/// ];
/// assert_eq!(fft_rust::approximate_complex(&mut fft_rust::ifft(arg)), &mut vec![
///     Complex::new(0f64, 0f64), Complex::new(2f64, 0f64),
///     Complex::new(7f64, 0f64), Complex::new(-1f64, 0f64)
/// ]);
///
/// let arg = &mut vec![
///     Complex::new(3f64, 0f64), Complex::new(9f64, -2f64),
///     Complex::new(-5f64, 0f64), Complex::new(9f64, 2f64)
/// ];
/// assert_eq!(fft_rust::approximate_complex(&mut fft_rust::ifft(arg)), &mut vec![
///     Complex::new(4f64, 0f64), Complex::new(3f64, 0f64),
///     Complex::new(-5f64, 0f64), Complex::new(1f64, 0f64)
/// ]);
/// ```
pub fn ifft(p: &mut Vec<Complex<f64>>) -> Vec<Complex<f64>> {
    // Given an array of complex numbers, p,
    // Recursively perform an Inverse Cooley-Tukey Fast Fourier Transform

    // if p is length 1, then just return p
    let mut n = p.len();
    if n == 1 {
        return p.iter().map(|x| *x).collect::<Vec<Complex<f64>>>();
    }

    // ensure length of p is a power of two
    if !is_int_power_of_two(n) {
        zero_pad_left(p, next_power_of_two(n) - p.len());
        n = p.len();
    }

    // even_powered_coeffs
    let y_e = ifft(
        &mut p[..].iter().step_by(2).map(|x| *x).collect::<Vec<_>>()
    );
    // odd_powered_coeffs
    let y_o = ifft(
        &mut p[1..].iter().step_by(2).map(|x| *x).collect::<Vec<_>>()
    );

    butterfly(n, y_e, y_o, true)
}

/// given two lists (one of even and one of odd coefficient outputs)
/// and boolean, inverse, run the butterfly algorithm; return the resulting array
fn butterfly(
    n: usize, y_e: Vec<Complex<f64>>, y_o: Vec<Complex<f64>>, inverse: bool
) -> Vec<Complex<f64>> {
    let sign = (if inverse {1} else {-1}) as f64;

    // create an array of e^(2i * pi) multiplied by k, (0 to n-1), and divided by N
    let omega = (0..n).map(
        |k| Complex::new(0.0, sign * 2.0 * PI * (k as f64)/(n as f64)).exp()
    ).collect::<Vec<Complex<f64>>>();

    let n_over_two = n as i32 / 2;
    let mut y = vec![Complex::new(0.0, 0.0); n];

    if inverse {
        for i in 0..n_over_two {
            let idx_a = i as usize;
            let idx_b = n_over_two as usize + idx_a;

            let temp = omega[idx_a] * y_o[idx_a];
            y[idx_a] = (y_e[idx_a] + temp) / 2f64;
            y[idx_b] = (y_e[idx_a] - temp) / 2f64;
        }
    } else {
        for i in 0..n_over_two {
            let idx_a = i as usize;
            let idx_b = n_over_two as usize + idx_a;

            let temp = omega[idx_a] * y_o[idx_a];
            y[idx_a] = y_e[idx_a] + temp;
            y[idx_b] = y_e[idx_a] - temp;
        }
    }
    y
}

/// round complex numbers that are really close to an integer
///
/// # Example
///
/// ```
/// use num::complex::Complex;
///
/// let arg = &mut vec![Complex::new(1f64, 0.000001f64), Complex::new(2.00000001f64, 0f64)];
/// assert_eq!(fft_rust::approximate_complex(arg), &mut vec![
///     Complex::new(1f64, 0f64), Complex::new(2f64, 0f64)
/// ]);
/// ```
pub fn approximate_complex(p :&mut Vec<Complex<f64>>) -> &mut Vec<Complex<f64>>{
    for mut complex_num in &mut *p {
        complex_num.re = complex_num.re.round();
        complex_num.im = complex_num.im.round();
    }
    p
}


/// do convolution/multiplication, i.e. fft -> multiply point-wise -> ifft
///
/// # Example
///
/// ```
/// use num::complex::Complex;
///
/// let list_a = &mut vec![  // binary 0110 -> decimal 6
///     Complex::new(0f64, 0f64), Complex::new(1f64, 0f64),
///     Complex::new(1f64, 0f64), Complex::new(0f64, 0f64)
/// ];
/// let list_b = &mut vec![  // binary 0010 -> decimal 2
///     Complex::new(0f64, 0f64), Complex::new(0f64, 0f64),
///     Complex::new(1f64, 0f64), Complex::new(0f64, 0f64)
/// ];
/// let output1 = &mut fft_rust::convolve(list_a, list_b);
/// assert_eq!(fft_rust::approximate_complex(output1), &mut vec![  // binary 1100 -> decimal 12
///     Complex::new(1f64, 0f64), Complex::new(1f64, 0f64),
///     Complex::new(0f64, 0f64), Complex::new(0f64, 0f64)
/// ]);
///
/// let list_c = &mut vec![  // decimal 25
///     Complex::new(0f64, 0f64), Complex::new(0f64, 0f64),
///     Complex::new(2f64, 0f64), Complex::new(5f64, 0f64)
/// ];
/// let list_d = &mut vec![  // decimal 25
///     Complex::new(0f64, 0f64), Complex::new(0f64, 0f64),
///     Complex::new(2f64, 0f64), Complex::new(5f64, 0f64)
/// ];
/// list_c.reverse();  // why do i need to reverse this one to get the answer??
/// let output2 = &mut fft_rust::convolve(list_c, list_d);
/// assert_eq!(fft_rust::approximate_complex(output2), &mut vec![  // decimal 400 + 200 + 25 = 625
///     Complex::new(0f64, 0f64), Complex::new(4f64, 0f64),
///     Complex::new(20f64, 0f64), Complex::new(25f64, 0f64)
/// ]);
/// ```
pub fn convolve(
    list1 :&mut Vec<Complex<f64>>, list2 :&mut Vec<Complex<f64>>
) -> Vec<Complex<f64>> {
    let fft_list1 = &mut fft(list1);
    let fft_list2 = &mut fft(list2);
    let conv_output = &mut multiply_pointwise(fft_list1, fft_list2);
    ifft(conv_output)
}

/// point-wise multiplication
///
/// # Example
///
/// ```
/// use num::complex::Complex;
///
/// let list_a = &mut vec![
///     Complex::new(0f64, 0f64), Complex::new(1f64, 0f64),
///     Complex::new(1f64, 0f64), Complex::new(0f64, 0f64)
/// ];
/// let list_b = &mut vec![
///     Complex::new(0f64, 0f64), Complex::new(0f64, 0f64),
///     Complex::new(1f64, 0f64), Complex::new(0f64, 0f64)
/// ];
/// let output1 = &mut fft_rust::multiply_pointwise(list_a, list_b);
/// assert_eq!(output1, &mut vec![
///     Complex::new(0f64, 0f64), Complex::new(0f64, 0f64),
///     Complex::new(1f64, 0f64), Complex::new(0f64, 0f64)
/// ]);
///
/// let list_c = &mut vec![
///     Complex::new(-2f64, 0f64), Complex::new(3f64, 0f64),
///     Complex::new(2f64, 0f64), Complex::new(5f64, 0f64)
/// ];
/// let list_d = &mut vec![
///     Complex::new(-5f64, 0f64), Complex::new(-3f64, 0f64),
///     Complex::new(2f64, 0f64), Complex::new(5f64, 0f64)
/// ];
/// let output2 = &mut fft_rust::multiply_pointwise(list_c, list_d);
/// assert_eq!(output2, &mut vec![
///     Complex::new(10f64, 0f64), Complex::new(-9f64, 0f64),
///     Complex::new(4f64, 0f64), Complex::new(25f64, 0f64)
/// ]);
/// ```
pub fn multiply_pointwise(
    list1 :&mut Vec<Complex<f64>>, list2 :&mut Vec<Complex<f64>>
) -> Vec<Complex<f64>> {
    let list1_len = list1.len();
    let list2_len = list2.len();
    if list1_len != list2_len {
        if list1_len > list2_len {
            zero_pad_left(list2, list1_len - list2_len);
        }
        else {
            zero_pad_left(list1, list2_len - list1_len);
        }
    }

    list1.iter().enumerate().map(|(i, v)| v.conj() * list2[i]).collect::<Vec<_>>()
}


/// add N zeros as padding to the left side of an array and return the padded array
///
/// # Example
///
/// ```
/// use num::complex::Complex;
///
/// let arg = &mut vec![Complex::new(1f64, 0f64), Complex::new(2f64, 0f64)];
/// let answer = fft_rust::zero_pad_left(arg, 2);
/// assert_eq!(answer, &mut vec![
///     Complex::new(0f64, 0f64), Complex::new(0f64, 0f64),
///     Complex::new(1f64, 0f64), Complex::new(2f64, 0f64)
/// ]);
///
/// let arg = &mut vec![Complex::new(1f64, 0f64)];
/// let answer = fft_rust::zero_pad_left(arg, 1);
/// assert_eq!(answer, &mut vec![Complex::new(0f64, 0f64), Complex::new(1f64, 0f64)]);
/// ```
pub fn zero_pad_left(
    list_to_pad: &mut Vec<Complex<f64>>, num_zeros: usize
) -> &mut Vec<Complex<f64>> {
    for x in vec![Complex::new(0f64, 0f64); num_zeros].iter() {
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
    fn fft_test() {
        assert_eq!(
            fft(&mut vec![
                Complex::new(2f64, 0f64), Complex::new(2f64, 0f64)
            ]),
            vec![Complex::new(4f64, 0f64), Complex::new(0f64, 0f64)]
        );
        assert_eq!(
            fft(&mut vec![
                Complex::new(1f64, 0f64), Complex::new(3f64, 0f64)
            ]),
            vec![Complex::new(4f64, 0f64), Complex::new(-2f64, 0f64)]
        );
        assert_eq!(
            fft(&mut vec![
                Complex::new(1f64, 0f64), Complex::new(1f64, 0f64),
                Complex::new(1f64, 0f64), Complex::new(1f64, 0f64)
            ]),
            vec![Complex::new(4f64, 0f64), Complex::new(0f64, 0f64),
                 Complex::new(0f64, 0f64), Complex::new(0f64, 0f64)
            ]
        );
        assert_eq!(
            fft(&mut vec![
                Complex::new(2f64, 0f64), Complex::new(7f64, 0f64),
                Complex::new(-1f64, 0f64)
            ]),
            vec![Complex::new(8f64, 0f64), Complex::new(-7f64, -3f64),
                 Complex::new(6f64, 0f64), Complex::new(-7f64, 3f64)
            ]
        );
        assert_eq!(
            fft(&mut vec![
                Complex::new(4f64, 0f64), Complex::new(3f64, 0f64),
                Complex::new(-5f64, 0f64), Complex::new(1f64, 0f64)
            ]),
            vec![Complex::new(3f64, 0f64), Complex::new(9f64, -2f64),
                 Complex::new(-5f64, 0f64), Complex::new(9f64, 2f64)
            ]
        );
    }

    #[test]
    fn ifft_test() {
        assert_eq!(
            ifft(&mut vec![
                Complex::new(4f64, 0f64), Complex::new(0f64, 0f64)
            ]),
            vec![Complex::new(2f64, 0f64), Complex::new(2f64, 0f64)]
        );
        assert_eq!(
            ifft(&mut vec![
                Complex::new(4f64, 0f64), Complex::new(-2f64, 0f64)
            ]),
            vec![Complex::new(1f64, 0f64), Complex::new(3f64, 0f64)]
        );
        assert_eq!(
            ifft(&mut vec![
                Complex::new(4f64, 0f64), Complex::new(0f64, 0f64),
                Complex::new(0f64, 0f64), Complex::new(0f64, 0f64)
            ]),
            vec![Complex::new(1f64, 0f64), Complex::new(1f64, 0f64),
                 Complex::new(1f64, 0f64), Complex::new(1f64, 0f64)
            ]
        );
    }

    #[test]
    fn ifft_reverses_fft_test() {
        assert_eq!(
            approximate_complex(&mut ifft( &mut fft(&mut vec![
                Complex::new(2f64, 0f64), Complex::new(2f64, 0f64)
            ]))),
            &mut vec![Complex::new(2f64, 0f64), Complex::new(2f64, 0f64)]
        );
        assert_eq!(
            approximate_complex(&mut ifft( &mut fft(&mut vec![
                Complex::new(2f64, 0f64), Complex::new(2f64, 0f64),
                Complex::new(2f64, 0f64)
            ]))),
            &mut vec![Complex::new(0f64, 0f64), Complex::new(2f64, 0f64),
                 Complex::new(2f64, 0f64), Complex::new(2f64, 0f64)
            ]
        );
        assert_eq!(
            approximate_complex(&mut ifft( &mut fft(&mut vec![
                Complex::new(2f64, 0f64), Complex::new(7f64, 0f64),
                Complex::new(-1f64, 0f64)
            ]))),
            &mut vec![Complex::new(0f64, 0f64), Complex::new(2f64, 0f64),
                 Complex::new(7f64, 0f64), Complex::new(-1f64, 0f64)
            ]
        );
    }

    #[test]
    fn convolve_test() {
        let list_a = &mut vec![  // binary 0110 -> decimal 6
            Complex::new(0f64, 0f64), Complex::new(1f64, 0f64),
            Complex::new(1f64, 0f64), Complex::new(0f64, 0f64)
        ];
        let list_b = &mut vec![  // binary 0010 -> decimal 2
            Complex::new(0f64, 0f64), Complex::new(0f64, 0f64),
            Complex::new(1f64, 0f64), Complex::new(0f64, 0f64)
        ];
        let output = &mut convolve(list_a, list_b);
        assert_eq!(approximate_complex(output), &mut vec![  // binary 1100 -> decimal 12
            Complex::new(1f64, 0f64), Complex::new(1f64, 0f64),
            Complex::new(0f64, 0f64), Complex::new(0f64, 0f64)
        ]);

        let list_c = &mut vec![  // decimal 25
            Complex::new(0f64, 0f64), Complex::new(0f64, 0f64),
            Complex::new(2f64, 0f64), Complex::new(5f64, 0f64)
        ];
        let list_d = &mut vec![  // decimal 25
            Complex::new(0f64, 0f64), Complex::new(0f64, 0f64),
            Complex::new(2f64, 0f64), Complex::new(5f64, 0f64)
        ];
        list_c.reverse();  // why do i need to reverse this one to get the answer??
        let output2 = &mut convolve(list_c, list_d);
        assert_eq!(approximate_complex(output2), &mut vec![  // decimal 400 + 200 + 25 = 625
            Complex::new(0f64, 0f64), Complex::new(4f64, 0f64),
            Complex::new(20f64, 0f64), Complex::new(25f64, 0f64)
        ]);
    }

    #[test]
    fn multiply_pointwise_test() {
        let list_a = &mut vec![  // binary 0110 -> decimal 6
            Complex::new(0f64, 0f64), Complex::new(1f64, 0f64),
            Complex::new(1f64, 0f64), Complex::new(0f64, 0f64)
        ];
        let list_b = &mut vec![  // binary 0010 -> decimal 2
            Complex::new(0f64, 0f64), Complex::new(0f64, 0f64),
            Complex::new(1f64, 0f64), Complex::new(0f64, 0f64)
        ];
        let output = &mut multiply_pointwise(list_a, list_b);
        assert_eq!(output, &mut vec![  // binary 0010 -> decimal 2
            Complex::new(0f64, 0f64), Complex::new(0f64, 0f64),
            Complex::new(1f64, 0f64), Complex::new(0f64, 0f64)
        ]);
    }

    #[test]
    fn approximate_complex_test() {
        let arg = &mut vec![Complex::new(1f64, 0.000001f64), Complex::new(2.00000001f64, 0f64)];
        assert_eq!(approximate_complex(arg), &mut vec![
            Complex::new(1f64, 0f64), Complex::new(2f64, 0f64)
        ]);
    }

    #[test]
    fn zero_pad_left_test() {
        assert_eq!(
            zero_pad_left(&mut vec![Complex::new(1f64, 0f64)], 1),
            &mut vec![Complex::new(0f64, 0f64), Complex::new(1f64, 0f64)]
        );
        assert_eq!(
            zero_pad_left(&mut vec![Complex::new(1f64, 0f64)], 0),
            &mut vec![Complex::new(1f64, 0f64)]
        );
        assert_eq!(
            zero_pad_left(&mut vec![
                Complex::new(1f64, 0f64), Complex::new(2f64, 0f64),
                Complex::new(3f64, 0f64), Complex::new(4f64, 0f64)
            ], 4),
            &mut vec![
                Complex::new(0f64, 0f64), Complex::new(0f64, 0f64),
                Complex::new(0f64, 0f64), Complex::new(0f64, 0f64),
                Complex::new(1f64, 0f64), Complex::new(2f64, 0f64),
                Complex::new(3f64, 0f64), Complex::new(4f64, 0f64)
            ]
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
