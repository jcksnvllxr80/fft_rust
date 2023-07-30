// use std::fs::File;
// use std::io::BufReader;
// use rodio::{Decoder, OutputStream, source::Source};
use std::fs::File;
use std::path::Path;


fn main() {
    // // Get a output stream handle to the default physical sound device
    // let (_stream, stream_handle) = OutputStream::try_default().unwrap();
    // // Load a sound from a file, using a path relative to Cargo.toml
    // let file = BufReader::new(File::open(
    //     "C:/Users/A-A-Ron/Music/royalty free youtube music/birdhouse baby bird timelapse/Cosmic Drift - DivKid.mp3"
    // ).unwrap());
    // // Decode that sound file into a source
    // let source = Decoder::new(file).unwrap();
    // // Play the sound directly on the device
    // let data = source.convert_samples();
    // // stream_handle.play_raw(source.convert_samples());
    //
    // // The sound plays in a separate audio thread,
    // // so we need to keep the main thread alive while it's playing.
    // std::thread::sleep(std::time::Duration::from_secs(30));

    let mut inp_file = File::open(Path::new("audio/song.wav")).unwrap();
    let (header, data) = wav::read(&mut inp_file).unwrap();
    println!("header -> {:?}", header);
    let song_bytes = data.as_sixteen().unwrap();
    println!("song_bytes -> {:?}", song_bytes[441000..441100].iter().step_by(2).map(
        |x| x
    ).collect::<Vec<_>>());
}