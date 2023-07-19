use druid::{
    widget::{Flex, Label, Slider},
    AppLauncher, Data, Lens, Widget, WidgetExt, WindowDesc,
};
use plotters::prelude::*;
use plotters_druid::Plot;
use std::f32::consts::PI;


#[derive(Clone, Data, Lens)]
struct State {
    fq_1: f64,
}

fn build_plot_widget() -> impl Widget<State> {
    Plot::new(|_size, data: &State, root| {
        let fq_1 = data.fq_1 as f32;

        // let res = 400;
        let font1 = FontDesc::new(
            FontFamily::SansSerif,
            16.,
            FontStyle::Normal
        );
        let font2 = FontDesc::new(
            FontFamily::SansSerif,
            16.,
            FontStyle::Normal
        );
        let font3 = FontDesc::new(
            FontFamily::SansSerif,
            16.,
            FontStyle::Normal
        );

        let root = root.titled(
            "Time to Frequency Domain with FFT", ("sans-serif", 60)
        ).unwrap();
        let child_drawing_areas = root.split_evenly((3, 1));

        let mut time_chart = ChartBuilder::on(
            &child_drawing_areas[0]
        )
            .x_label_area_size(30)
            .y_label_area_size(30)
            .margin_right(10)
            .build_cartesian_2d(0.0..1_f32, -1.1..1.1_f32)
            .unwrap();

        time_chart
            .configure_mesh()
            .axis_style(&RGBColor(28, 28, 28))
            .x_label_style(font1.clone().with_color(&WHITE))
            .y_label_style(font1.clone().with_color(&WHITE))
            .draw()
            .unwrap();

        let mut combined_time_chart = ChartBuilder::on(
            &child_drawing_areas[1]
        )
            .x_label_area_size(30)
            .y_label_area_size(30)
            .margin_right(10)
            .build_cartesian_2d(0.0..1_f32, -1.8..1.8_f32)
            .unwrap();

        combined_time_chart
            .configure_mesh()
            .axis_style(&RGBColor(28, 28, 28))
            .x_label_style(font1.clone().with_color(&WHITE))
            .y_label_style(font1.clone().with_color(&WHITE))
            .draw()
            .unwrap();

        let mut freq_chart = ChartBuilder::on(
            &child_drawing_areas[2]
        )
            .x_label_area_size(30)
            .y_label_area_size(30)
            .margin_right(10)
            .build_cartesian_2d(0.0..510_f32, -0.1..1.1_f32)
            .unwrap();

        freq_chart
            .configure_mesh()
            .axis_style(&RGBColor(28, 28, 28))
            .x_label_style(font2.clone().with_color(&WHITE))
            .y_label_style(font2.clone().with_color(&WHITE))
            .draw()
            .unwrap();

        let fs = 1024.0;
        let period = 1.0/fs;  // Sampling period
        let length = 2048;  // Length of signal - needs to be a power of two
        let t : Vec<f32> = (0..(length-1)).map(|x| (x as f32)*period).collect();  // Time vector
        let amp_1 = 1.0;
        let amp_2  = 0.7;
        let fq_2 = 70.0;

        // create sigs 1 and 2
        let y_sig_1 : Vec<_> = t.iter().map(
            |x| amp_1*(2.0*PI*fq_1*x).sin() as f64
        ).collect();
        let y_sig_2 : Vec<_> = t.iter().map(
            |x| amp_2*(2.0*PI*fq_2*x).sin() as f64
        ).collect();

        // combine sigs 1 and 2
        let y_sig : Vec<_> = y_sig_1.iter().zip(&y_sig_2).map(
            |(y1, y2)| y1 + y2
        ).collect();
        // zip in t to the each signal
        let signal = t.iter().map(
            |t| *t
        ).zip(y_sig.iter().map(|y| *y as f32));

        let signal_1 = t.iter().map(
            |t| *t
        ).zip(y_sig_1.iter().map(|y| *y as f32));

        let signal_2 = t.iter().map(
            |t| *t
        ).zip(y_sig_2.iter().map(|y| *y as f32));

        let color1 = Palette99::pick(1);
        let color2 = Palette99::pick(2);
        let color3 = Palette99::pick(3);
        let color4 = Palette99::pick(4);

        time_chart
            .draw_series(LineSeries::new(signal_1, &color3))
            .unwrap()
            .label("signal 1")
            .legend(move |(x, y)| {
                PathElement::new(
                    vec![(x, y), (x + 20, y)],
                    ShapeStyle::from(&color3).stroke_width(2),
                )
            });

        time_chart
            .draw_series(LineSeries::new(signal_2, &color4))
            .unwrap()
            .label("signal 2")
            .legend(move |(x, y)| {
                PathElement::new(
                    vec![(x, y), (x + 20, y)],
                    ShapeStyle::from(&color4).stroke_width(2),
                )
            });

        time_chart
            .configure_series_labels()
            .position(SeriesLabelPosition::UpperRight)
            .background_style(&RGBColor(41, 41, 41))
            .border_style(&RGBColor(28, 28, 28))
            .label_font(font1.with_color(&WHITE))
            .draw()
            .unwrap();

        combined_time_chart
            .draw_series(LineSeries::new(signal, &color1))
            .unwrap()
            .label("signal 1 + signal 2")
            .legend(move |(x, y)| {
                PathElement::new(
                    vec![(x, y), (x + 20, y)],
                    ShapeStyle::from(&color1).stroke_width(2),
                )
            });

        combined_time_chart
            .configure_series_labels()
            .position(SeriesLabelPosition::UpperRight)
            .background_style(&RGBColor(41, 41, 41))
            .border_style(&RGBColor(28, 28, 28))
            .label_font(font3.with_color(&WHITE))
            .draw()
            .unwrap();

        let mut y_complex = fft_rust::float_array_to_complex(y_sig);
        let mut y_fft = fft_rust::fft(&mut y_complex);

        let y_fft_real = fft_rust::complex_array_to_float(&mut y_fft);
        let y_fft_real_scale: Vec<_> = y_fft_real.iter().map(
            |x| *x/(length as f64)
        ).collect();
        let p1_temp: Vec<f64> = y_fft_real_scale[..=(length/2 + 1)].iter().map(
            |x| *x
        ).collect();
        let p1: Vec<f64>  = p1_temp[..].iter().enumerate().map(
            |(i, p)| if i == 0 {*p} else if i == p1_temp.len() - 1 {*p} else {2.0*p}
        ).collect();
        let f: Vec<f32> = (0..(length/2)).map(
            |i| (i as f32)*(fs as f32)/(length as f32)
        ).collect();
        let fft = f.iter().map(
            |f| *f
        ).zip(
            p1.iter().map(
                |p| *p as f32
            )
        );

        freq_chart
            .draw_series(LineSeries::new(fft, &color2))
            .unwrap()
            .label("FFT of combined signal")
            .legend(move |(x, y)| {
                PathElement::new(
                    vec![(x, y), (x + 20, y)],
                    ShapeStyle::from(&color2).stroke_width(2),
                )
            });

        freq_chart
            .configure_series_labels()
            .position(SeriesLabelPosition::UpperRight)
            .background_style(&RGBColor(41, 41, 41))
            .border_style(&RGBColor(28, 28, 28))
            .label_font(font2.with_color(&WHITE))
            .draw()
            .unwrap();
    })
}

fn build_slider_widget(name: String, min: f64, max: f64) -> impl Widget<f64> {
    Flex::row()
        .with_child(Label::new(name))
        .with_flex_child(
            Slider::new().with_range(min, max).env_scope(|env, _| {
                // remove the width limit in [`Slider`]
                env.set(druid::theme::WIDE_WIDGET_WIDTH, f64::INFINITY)
            }),
            1.,
        )
        .must_fill_main_axis(true)
        .fix_height(20.)
}

fn build_root_widget() -> impl Widget<State> {
    Flex::column()
        .with_flex_child(build_plot_widget(), 1.)
        .with_spacer(5.)
        .with_child(build_slider_widget(
            "freq_1".to_string(),
            0., 500.
        ).lens(State::fq_1))
        .padding(10.)
}

fn main() {
    let main_window = WindowDesc::new(build_root_widget())
        .title("FFT: time vs freq")
        .window_size((1500.0, 850.0));

    let freq1 = 50.0;
    AppLauncher::with_window(main_window)
        .launch(State { fq_1: freq1 })
        .expect("Failed to launch application");
}
