use image::{EncodableLayout, GrayImage};
use image::io::Reader as ImageReader;
use num::complex::ComplexFloat;
use std::time::Instant;

fn main() {
    let start = Instant::now();

    let img = ImageReader::open("img/baboon.ppm").unwrap().decode().unwrap();
    let gray_img = img.to_luma8();
    // gray_img.save(format!("img/output/gray_img.jpg")).unwrap();

    let fft_img = &mut fft_2d(gray_img);
    let fft_img_shift = &mut fft_shift(fft_img, "fft_img_shift");
    let fft_reverse_img_shift = fft_shift(fft_img_shift, "fft_reverse_img_shift");
    let _ifft_img = &mut ifft_2d(fft_reverse_img_shift);

    let duration = start.elapsed();
    println!("Time elapsed to do 2-D fft and shift of image is: {:?}", duration);
}

fn fft_shift(img: &mut GrayImage, filename: &str) -> GrayImage {
    let (width, height) = img.dimensions();

    let mut shifted_fft_img_20_ln = GrayImage::new(width, height);
    let mut shifted_fft_img = GrayImage::new(width, height);
    for (y1, row) in img.enumerate_rows() {
        for (x1, _, p) in row {
            let x2 = ((x1 + (width) / 2) % width);
            let y2 = (y1 + height / 2) % height;
            shifted_fft_img.put_pixel(x2, y2, *p);
        }
    }
    for (i, pixel) in shifted_fft_img.as_bytes().iter().enumerate() {
        let x = i % width as usize;
        let y = i / width as usize;
        shifted_fft_img_20_ln.put_pixel(x as u32,y as u32, image::Luma([(20.0 * (*pixel as f32).ln()) as u8]));
    }
    shifted_fft_img_20_ln.save(format!("img/output/{filename}.jpg")).unwrap();
    shifted_fft_img
}

fn fft_2d(img: GrayImage) -> GrayImage {
    let (width, height) = img.dimensions();

    let mut row_fft_mat = vec![0.0f64; (width * height) as usize];
    let mut row_fft_img_20_ln = GrayImage::new(width, height);
    for y in 0..height {
        let row_idx = y * width;

        let float_img: Vec<_> = img.as_bytes().to_vec()[row_idx as usize..(row_idx + width) as usize].iter().map(
            |x| *x as f64
        ).collect();

        let mut row_complex = fft_rust::float_array_to_complex(float_img);

        let mut row_fft = fft_rust::fft(&mut row_complex);

        let row_fft_real = fft_rust::complex_array_to_float(&mut row_fft);

        for (i, pixel) in row_fft_real[row_fft_real.len() - width as usize..].iter().enumerate() {
            row_fft_mat[(row_idx + i as u32) as usize] = *pixel;
            // let arr_idx = row_idx + i as u32;
            // println!("{}", format!("row_fft_mat[{arr_idx}] = {pixel};"))
        }

        for (i, pixel) in row_fft_real[row_fft_real.len() - width as usize..].iter().enumerate() {
            row_fft_img_20_ln.put_pixel(i as u32, y, image::Luma([(20.0 * pixel.ln()) as u8]));
        }
    }
    row_fft_img_20_ln.save(format!("img/output/row_fft_img_20_ln.jpg")).unwrap();

    let mut fully_fft_img = GrayImage::new(width, height);
    let mut fully_fft_img_20_ln = GrayImage::new(width, height);
    for x in 0..width {
        // let col_idx = x * width;

        let float_img: Vec<_> = row_fft_mat[x as usize..].iter().step_by(width as usize).map(
            |x| *x as f64
        ).collect();

        let mut col_complex = fft_rust::float_array_to_complex(float_img);

        let mut col_fft = fft_rust::fft(&mut col_complex);

        let col_fft_real = fft_rust::complex_array_to_float(&mut col_fft);

        for (i, pixel) in col_fft_real[col_fft_real.len() - height as usize..].iter().enumerate() {
            fully_fft_img.put_pixel(x,i as u32, image::Luma([*pixel as u8]));
        }

        for (i, pixel) in col_fft_real[col_fft_real.len() - height as usize..].iter().enumerate() {
            fully_fft_img_20_ln.put_pixel(x,i as u32, image::Luma([(20.0 * pixel.ln()) as u8]));
        }
    }
    fully_fft_img_20_ln.save(format!("img/output/fft_img.jpg")).unwrap();
    fully_fft_img
}

fn ifft_2d(img: GrayImage) -> GrayImage {
    let (width, height) = img.dimensions();

    let mut row_ifft_mat = vec![0.0f64; (width * height) as usize];
    let mut row_ifft_img_20_ln = GrayImage::new(width, height);
    for y in 0..height {
        let row_idx = y * width;

        let float_img: Vec<_> = img.as_bytes().to_vec()[row_idx as usize..(row_idx + width) as usize].iter().map(
            |x| *x as f64
        ).collect();

        let mut row_complex = fft_rust::float_array_to_complex(float_img);

        let mut row_ifft = fft_rust::ifft(&mut row_complex);

        let row_ifft_real = fft_rust::complex_array_to_float(&mut row_ifft);

        for (i, pixel) in row_ifft_real[row_ifft_real.len() - width as usize..].iter().enumerate() {
            row_ifft_mat[(row_idx + i as u32) as usize] = *pixel;
            // let arr_idx = row_idx + i as u32;
            // println!("{}", format!("row_fft_mat[{arr_idx}] = {pixel};"))
        }

        for (i, pixel) in row_ifft_real[row_ifft_real.len() - width as usize..].iter().enumerate() {
            row_ifft_img_20_ln.put_pixel(i as u32, y, image::Luma([(20.0 * pixel.ln()) as u8]));
        }
    }
    row_ifft_img_20_ln.save(format!("img/output/row_ifft_img_20_ln.jpg")).unwrap();

    let mut fully_ifft_img = GrayImage::new(width, height);
    let mut fully_ifft_img_20_ln = GrayImage::new(width, height);
    for x in 0..width {
        // let col_idx = x * width;

        let float_img: Vec<_> = row_ifft_mat[x as usize..].iter().step_by(width as usize).map(
            |x| *x as f64
        ).collect();

        let mut col_complex = fft_rust::float_array_to_complex(float_img);

        let mut col_ifft = fft_rust::ifft(&mut col_complex);

        let col_ifft_real = fft_rust::complex_array_to_float(&mut col_ifft);

        for (i, pixel) in col_ifft_real[col_ifft_real.len() - height as usize..].iter().enumerate() {
            fully_ifft_img.put_pixel(x,i as u32, image::Luma([*pixel as u8]));
        }

        for (i, pixel) in col_ifft_real[col_ifft_real.len() - height as usize..].iter().enumerate() {
            fully_ifft_img_20_ln.put_pixel(x,i as u32, image::Luma([(20.0 * pixel.ln()) as u8]));
        }
    }
    fully_ifft_img_20_ln.save(format!("img/output/ifft_img.jpg")).unwrap();
    fully_ifft_img
}