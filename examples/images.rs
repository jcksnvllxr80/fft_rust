use image::GenericImageView;
use image::io::Reader as ImageReader;

fn main() {
    let img = ImageReader::open("img/lenna.png").unwrap().decode().unwrap();
    let sub_img = img.view(128,128,25,25).to_image();
    println!("dimensions -> {:?}", sub_img.dimensions());
    // println!("bytes -> {:?}", sub_img.as_bytes());
    println!("small img -> {:?}", img.view(256,256,2,2).to_image());
    println!("color -> {:?}", img.color());
    println!("pixel (x0, y0) -> {:?}", img.get_pixel(256,256));
    println!("pixel (x1, y0) -> {:?}", img.get_pixel(257,256));
    println!("pixel (x0, y1) -> {:?}", img.get_pixel(256,257));
    println!("pixel (x1, y1) -> {:?}", img.get_pixel(257,257));

    sub_img.save("img/output.jpg").unwrap();
}