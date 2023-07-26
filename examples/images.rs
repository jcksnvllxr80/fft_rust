use image::{EncodableLayout, GenericImageView};
use image::io::Reader as ImageReader;

fn main() {
    // let img = ImageReader::open("img/lenna.png").unwrap().decode().unwrap();
    // let img = ImageReader::open("img/airplane.bmp").unwrap().decode().unwrap();
    let img = ImageReader::open("img/baboon.ppm").unwrap().decode().unwrap();

    let img_bytes = img.as_bytes();
    let img_dim = img.dimensions();

    let img_r_bytes: Vec<_> = img_bytes[..].iter().step_by(3).collect();
    let img_g_bytes: Vec<_> = img_bytes[1..].iter().step_by(3).collect();
    let img_b_bytes: Vec<_> = img_bytes[2..].iter().step_by(3).collect();

    component_to_grey(String::from("grey_red_img"), img_r_bytes, img_dim);
    component_to_grey(String::from("grey_green_img"), img_g_bytes, img_dim);
    component_to_grey(String::from("grey_blue_img"), img_b_bytes, img_dim);
}

fn component_to_grey(out_file_name: String, img_bytes: Vec<&u8>, img_dim: (u32, u32)) {
    let mut component_img = image::GrayImage::new(img_dim.0, img_dim.1);
    for x in 0..img_dim.0 {
        for y in 0..img_dim.1 {
            component_img.put_pixel(x, y, image::Luma([*img_bytes[(y*img_dim.0 + x) as usize]]));
        }
    }
    component_img.save(format!("img/output/{out_file_name}.jpg")).unwrap();
}