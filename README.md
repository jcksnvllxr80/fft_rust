# Cooley-Tukey FFT in Rust

implemented cooley-tukey fft and ifft in Rust with documentation, tests, and examples galore.

## purpose

The reason for creating this is that I wanted to practice rust and get more familiar with the fft algorithm.

## prerequisites

need to have rust installed

## getting started

clone this repo and read the documentation

## documentation

run the following once you have cloned the repo
```
cargo doc --open
```

## future work / improvements

- add parallelism with threads to speed up the algorithm
- pull out the redundant and costly weight (omega) operations that are based on N outside the recursion  
  since they only need to be calculated once for each N
- use the lib by playing back an .mp3 file while graphing its frequency domain
- use the lib by creating a convolution reverb; try it out on some audio files
  - drum hits
  - voice recordings
