use std::f32::consts::PI;
use plotters::prelude::*;
use plotters::style::full_palette::{GREEN_900, PURPLE};
use fft_rust;


const OUT_FILE_NAME: &'static str = "plotters-doc-data/fft.png";
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let root_area = BitMapBackend::new(
        OUT_FILE_NAME, (1024, 768)
    ).into_drawing_area();
    root_area.fill(&WHITE)?;
    let root_area = root_area.titled(
        "Deconstruct waveforms with FFT", ("sans-serif", 60)
    )?;
    let child_drawing_areas = root_area.split_evenly((3, 1));

    // define parameters
    let fs = 1024.0;
    let period = 1.0/fs;  // Sampling period
    let length = 2048;  // Length of signal - needs to be a power of two
    let t : Vec<f32> = (0..(length-1)).map(|x| (x as f32)*period).collect();  // Time vector
    // amplitude and freq for 2 sine waves
    let fq_1 = 50.0;
    let amp_1 = 0.7;
    let fq_2 = 120.0;
    let amp_2 = 1.0;

    let signal1: Vec<f32>  = t.iter().map(|x| amp_1*(2.0*PI*fq_1*x).sin()).collect();
    let signal2: Vec<f32>  = t.iter().map(|x| amp_2*(2.0*PI*fq_2*x).sin()).collect();

    let mut cc = ChartBuilder::on(
        &child_drawing_areas[0]
    )
        .margin(5)
        .set_all_label_area_size(50)
        .caption(
            "green: {freq: 50; amp: 0.7}; purple: {freq: 120; amp: 1.0}",
            ("sans-serif", 40)
        ).build_cartesian_2d(0.0f32..0.3, -1.1f32..1.1f32)?;

    cc.configure_mesh().x_labels(20).y_labels(10).disable_mesh()
        .x_label_formatter(&|v| format!("{:.1}", v))
        .y_label_formatter(&|v| format!("{:.1}", v))
        .draw()?;

    cc.draw_series(LineSeries::new(
        t.iter().map(|x| *x).zip(signal1.iter().map(|y| *y)),
        &GREEN_900,
    ))?
        .label("sinusoid 1")
        .legend(|(x, y)| PathElement::new(
            vec![(x, y), (x + 20, y)],
            &GREEN_900)
        );

    cc.draw_series(LineSeries::new(
        t.iter().map(|x| *x).zip(signal2.iter().map(|y| *y)),
        &PURPLE,
    ))?
        .label("sinusoid 2")
        .legend(|(x, y)| PathElement::new(
            vec![(x, y), (x + 20, y)],
            &PURPLE)
        );


    let combined_signal: Vec<f32> = t.iter().enumerate().map(
        |(i, _)| signal1[i] + signal2[i]
    ).collect();

    let mut cc = ChartBuilder::on(
        &child_drawing_areas[1]
    )
        .margin(5)
        .set_all_label_area_size(50)
        .caption(
            "Combined purple and green waveforms", ("sans-serif", 40)
        ).build_cartesian_2d(0.0f32..0.3, -2.1f32..2.1f32)?;

    cc.configure_mesh().x_labels(20).y_labels(10).disable_mesh()
        .x_label_formatter(&|v| format!("{:.1}", v))
        .y_label_formatter(&|v| format!("{:.1}", v))
        .draw()?;

    cc.draw_series(LineSeries::new(
        t.iter().map(|x| *x).zip(combined_signal.iter().map(|y| *y)),
        &BLUE,
    ))?
        .label("sinusoid")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));


    let y = combined_signal.iter().map(|x| *x as f64).collect();
    let mut y_complex = fft_rust::float_array_to_complex(y);
    let mut y_fft = fft_rust::fft(&mut y_complex);
    // println!("delta length between original signal and the fft: {:?}", (y_fft.len() - length));

    let y_fft_real = fft_rust::complex_array_to_float(&mut y_fft);
    let y_fft_real_scale: Vec<_> = y_fft_real.iter().map(
        |x| *x/(length as f64)
    ).collect();
    let p1_temp: Vec<f64> = y_fft_real_scale[..=(length/2 + 1)].iter().map(|x| *x).collect();
    let p1: Vec<f64>  = p1_temp[..].iter().enumerate().map(
        |(i, p)| if i == 0 {*p} else if i == p1_temp.len() - 1 {*p} else {2.0*p}
    ).collect();
    let f: Vec<f32> = (0..(length/2)).map(|i| (i as f32)*fs/(length as f32)).collect();
    let fft_graph = f.iter().zip(p1.iter());

    let mut cc = ChartBuilder::on(
        &child_drawing_areas[2]
    )
        .margin(5)
        .set_all_label_area_size(20)
        .caption("FFT of Combined Waveform shows freq and amp of each signal",
                 ("sans-serif", 40)
        )
        .build_cartesian_2d(0.0f32..200f32, 0f32..1.1f32)?;

    cc.configure_mesh().x_labels(20).y_labels(10).disable_mesh()
        .x_label_formatter(&|v| format!("{:.1}", v))
        .y_label_formatter(&|v| format!("{:.1}", v))
        .draw()?;

    cc.draw_series(LineSeries::new(
        fft_graph.map(|(x, y)| ((*x as f32), (*y as f32))),
        &BLACK,
    ))?
        .label("FFT")
        .legend(
            |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &GREEN)
        );

    // To avoid the IO failure being ignored silently, we manually call the present function
    root_area.present().expect(
        "Unable to write to file, does the 'plotters-doc-data' dir exists under current dir"
    );
    println!("Result has been saved to {}", OUT_FILE_NAME);
    Ok(())
}

#[test]
fn entry_point() {
    main().unwrap()
}