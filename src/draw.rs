use plotters::prelude::*;

pub fn plot_salmon_runs(
    measured: Vec<f64>,
    calculated: Vec<f64>,
    start_year: i32,
    output_file: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let measured_years: Vec<i32> = (0..measured.len()).map(|i| start_year + i as i32).collect();

    let calculated_years: Vec<i32> = (0..calculated.len())
        .map(|i| start_year + i as i32)
        .collect();

    let end_year = start_year + measured.len().max(calculated.len()) as i32;

    let max_y = measured
        .iter()
        .chain(calculated.iter())
        .cloned()
        .filter(|v| v.is_finite())
        .fold(0.0_f64, f64::max);

    let root = BitMapBackend::new(output_file, (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("Salmon Runs (Measured vs Calculated)", ("sans-serif", 30))
        .margin(30)
        .x_label_area_size(40)
        .y_label_area_size(50)
        .build_cartesian_2d(start_year..end_year, 0.0..(max_y * 1.1 + 1.0))?;

    chart
        .configure_mesh()
        .x_desc("Year")
        .y_desc("Salmon Count")
        .draw()?;

    // Measured
    chart
        .draw_series(LineSeries::new(
            measured_years
                .iter()
                .zip(measured.iter())
                .map(|(x, y)| (*x, *y)),
            &BLUE,
        ))?
        .label("Measured")
        .legend(|(x, y)| PathElement::new([(x, y), (x + 20, y)], &BLUE));

    // Calculated
    chart
        .draw_series(LineSeries::new(
            calculated_years
                .iter()
                .zip(calculated.iter())
                .map(|(x, y)| (*x, *y)),
            &RED,
        ))?
        .label("Calculated")
        .legend(|(x, y)| PathElement::new([(x, y), (x + 20, y)], &RED));

    chart
        .configure_series_labels()
        .border_style(&BLACK)
        .draw()?;

    Ok(())
}