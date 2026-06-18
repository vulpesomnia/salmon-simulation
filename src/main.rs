use ggca::correlation::pearson_correlation;
use ndarray::Array2;
use ndarray::arr2;

mod draw;

fn main() {
    // AINEISTOT \\
    let measured_salmon_run: Vec<f64> = vec![
        31782.0, 17224.0, 23076.0, 59606.0, 53605.0, 101094.0, 58360.0, 98311.0, 40952.0, 47028.0,
        65520.0, 69149.0, 93121.0, 52030.0, 20245.0, 24627.0,
       45043.0, //0.0, 0.0, 0.0, 0.0, 0.0, // future
    ];
    let fishing_tornionjoki: Vec<f64> = vec![
        13400.0, 8600.0, 10100.0, 20000.0, 15200.0, 24200.0, 23200.0, 29200.0, 18000.0, 21100.0,
        20900.0, 23300.0, 30900.0, 19800.0, 10900.0, 10800.0,
        10800.0, 10800.0, 10800.0, 10800.0, 10800.0, 10800.0,
    ];
   /* let mut fishing_sea: Vec<f64> = vec![
        116333.0, 68618.0, 68867.0, 69057.0, 68335.0, 68589.0, 64771.0, 61873.0, 52386.0, 57208.0,
        57144.0, 52907.0, 52561.0, 40341.0, 28121.0, 15901.0,
    ];*/
    // COASTAL FISHING
    let fishing_coastal: Vec<f64> = vec![
        116333.0, 68618.0, 68867.0, 69057.0, 68335.0, 68589.0, 64771.0, 61873.0, 52386.0, 57208.0,
        57144.0, 52907.0, 52561.0, 38140.0, 31515.0, 28071.0, 
    ];
    // OFFSHORE FISHING
    let fishing_offshore: Vec<f64> = vec![
        60464.0, 79088.0, 74400.0, 57522.0, 33420.0, 30875.0, 27976.0, 26576.0, 28271.0, 36249.0,
        34090.0, 23188.0, 34235.0, 10298.0, 5039.0, 4977.0,
    ];
    let mut fishing_sea: Vec<f64> =
        fishing_coastal.iter()
        .zip(fishing_offshore.iter())
        .map(|(a, b)| a + b)
        .collect();
    fishing_sea.append(&mut vec![33048.0, 33048.0, 33048.0, 33048.0, 33048.0, 33048.0, 33048.0]);
    let seal_amounts: Vec<f64> = vec![
        15313.0, 19627.0, 23941.0, 28255.0, 29496.0, 32244.0, 30285.0, 30070.0, 30348.0, 34235.0,
        38121.0, 40075.0, 41936.0, 36761.0, 45798.0, 48711.0, 48087.0, 49882.0, 51678.0, 53473.0,
        55269.0, 57064.0
    ]; 
    // PARAMETRIT \\
    /*
    let alpha = 0.6;
    let beta = 0.5;
    let p_j: f64 = 0.3; // Merellejääntitodennäköisyys
    let spwn: f64 = 600.0; // Kutumäärä naaraslohille
    let p_l: f64 = 0.1; // Poikasten laskeutumistodennäköisyys
    */

    let mut best_spawn = 0.0;
    let mut best_alpha = 0.0;
    let mut best_beta = 0.0;
    let mut best_p_j = 0.0;
    let mut best_p_l = 0.0;
    let mut best_gamma = 0.0;
    
    let mut best_r2 = f64::MAX;
    let mut best_calculated_salmon_run: Vec<f64> = Vec::new();

    let BY_HAND = false;

    for spawn in 0..20 {
        println!("{}%", spawn as f64 / 20.0 * 100.0);
    for a in 0..30 {
        for b in 0..30 {
            for c in 0..60 {
                
            for pj in 0..10 {
                for pl in 0..10 {

                let spwn;
                let alpha;
                let beta;
                let gamma;
                let p_j;
                let p_l;

                if BY_HAND {
                    spwn = 1200.0;//600.0;
                    alpha = 0.8;//0.52;
                    beta = 1.3;//0.5;
                    gamma = 150000.0;//340000.0;
                    p_j = 0.6;//0.2;
                    p_l = 0.3;//0.1;
                } else {
                    spwn = spawn as f64 * 100.0;
                    alpha = a as f64 / 10.0;
                    beta = b as f64 / 10.0;
                    gamma = c as f64 * 10000.0;
                    p_j = pj as f64 / 10.0;
                    p_l = pl as f64 / 10.0;
                }
                    
    // MUUT ARVOT \\
    // Kutulohen selviytymistodennäköisyys
    let sps: f64 = 0.075;
    // Kudusta -> Smolt selviytymistodennäköisyys
    let fs = 0.01;
    // Hylkeiden tuottaman kuolleisuuden suhde meren ja joen välillä.
    let seal_weight_sea = 0.2;
    let seal_weight_river = 0.8;
    // LESLIE-MALLI \\
    let leslie_matrix = arr2(&[
        //FROM: P1, P2, P3, P4, P5; M1, J1; M2, J2; M3, J3; M4+, J4+;
        [
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, spwn, 0.0, spwn, 0.0, spwn, 0.0, spwn,
        ], // TO: P1
        [
            1.0 - p_l,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ], // P2
        [
            0.0,
            1.0 - p_l,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ], // P3
        [
            0.0,
            0.0,
            1.0 - p_l,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ], // P4
        [
            p_l, p_l, p_l, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ], // P5
        [
            0.0, 0.0, 0.0, 0.0, fs, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ], // M1
        [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0 - p_j,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ], // J1
        [
            0.0, 0.0, 0.0, 0.0, 0.0, p_j, sps, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ], // M2
        [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0 - p_j,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ], // J2
        [
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, p_j, sps, 0.0, 0.0, 0.0, 0.0,
        ], // M3
        [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0 - p_j,
            0.0,
            0.0,
            0.0,
        ], // J3
        [
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, p_j, sps, p_j, sps,
        ], // M4+
        [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0 - p_j,
            0.0,
        ], // J4+
    ]);
    let pop_length: usize = leslie_matrix.nrows();
    let sea_indices = vec![5, 7, 9, 11];
    let river_indices = vec![6, 8, 10, 12];


    let start_population: Array2<f64> = calc_start(1.0, 100, &leslie_matrix);
    let total_river_ratio: f64 = river_indices
        .iter()
        .map(|&i| start_population[[i, 0]])
        .sum::<f64>();
    let scale = measured_salmon_run[0] as f64 / (total_river_ratio);
    let start_population = &start_population * scale * beta;
    let mut salmon_amounts = Vec::new(); // 0 = Naaras, 1 = Uros
    for _ in 0..2 {
        salmon_amounts.push(start_population.clone() / 2.0);
    }
    let simulation_years: usize = measured_salmon_run.len();
    let mut calculated_salmon_run: Vec<f64> = vec![0.0; simulation_years];
    for i in 0..simulation_years {
        // Laske kalastus vuodelle i:
        // Tornionjoen lohen osuus Itämerellä
        let sea_amt: f64 = sea_indices.iter().map(|j| salmon_amounts[0][[*j, 0]] + salmon_amounts[1][[*j, 0]]).sum();
        let fishing_sea_tornionjoki_weight: f64 = sea_amt / (sea_amt + gamma);

        let h_sea = fishing_sea[i] * fishing_sea_tornionjoki_weight;
        let h_river = fishing_tornionjoki[i];
        let total_sea_female = sea_indices
            .iter()
            .map(|&k| salmon_amounts[0][[k, 0]])
            .sum::<f64>();
        let total_sea_male = sea_indices
            .iter()
            .map(|&k| salmon_amounts[1][[k, 0]])
            .sum::<f64>();
        let total_river_female = river_indices
            .iter()
            .map(|&k| salmon_amounts[0][[k, 0]])
            .sum::<f64>();
        let total_river_male = river_indices
            .iter()
            .map(|&k| salmon_amounts[1][[k, 0]])
            .sum::<f64>();
        let h_sea_female = h_sea * total_sea_female / (total_sea_female + total_sea_male);
        let h_sea_male = h_sea * total_sea_male / (total_sea_female + total_sea_male);
        let h_river_female = h_river * total_river_female / (total_river_female + total_river_male);
        let h_river_male = h_river * total_river_male / (total_river_female + total_river_male);
        let mut fishing_female: Array2<f64> = Array2::zeros((pop_length, 1));
        let mut fishing_male: Array2<f64> = Array2::zeros((pop_length, 1));
        // Älä ota ensimmäisen vuoden kalastusta mukaan:
        if i != 0 {
            distribute_value(
                h_sea_female,
                &sea_indices,
                &salmon_amounts[0],
                &mut fishing_female,
            );
            distribute_value(
                h_sea_male,
                &sea_indices,
                &salmon_amounts[1],
                &mut fishing_male,
            );
            distribute_value(
                h_river_female,
                &river_indices,
                &salmon_amounts[0],
                &mut fishing_female,
            );
            distribute_value(
                h_river_male,
                &river_indices,
                &salmon_amounts[1],
                &mut fishing_male,
            );
        }
        // Laske hyljekuolleisuus:
        let s_sea = seal_amounts[i] * seal_weight_sea;
        let s_river = seal_amounts[i] * seal_weight_river;
        let seal_mortality_sea_female =
            s_sea * total_sea_female / (total_sea_female + total_sea_male);
        let seal_mortality_sea_male = s_sea * total_sea_male / (total_sea_female + total_sea_male);
        let seal_mortality_river_female =
            s_river * total_river_female / (total_river_female + total_river_male);
        let seal_mortality_river_male =
            s_river * total_river_male / (total_river_female + total_river_male);
        let mut seal_mortality_female: Array2<f64> = Array2::zeros((pop_length, 1));
        let mut seal_mortality_male: Array2<f64> = Array2::zeros((pop_length, 1));
        distribute_value(
            seal_mortality_sea_female,
            &sea_indices,
            &salmon_amounts[0],
            &mut seal_mortality_female,
        );
        distribute_value(
            seal_mortality_sea_male,
            &sea_indices,
            &salmon_amounts[1],
            &mut seal_mortality_male,
        );
        distribute_value(
            seal_mortality_river_female,
            &river_indices,
            &salmon_amounts[0],
            &mut seal_mortality_female,
        );
        distribute_value(
            seal_mortality_river_male,
            &river_indices,
            &salmon_amounts[1],
            &mut seal_mortality_male,
        );
        // Laske seuraava Leslie-mallin askel:
        // Samalla jaetaan kutu naaraiden ja urosten kesken.
        // NAARAS
        let mut female_matrix = leslie_matrix.clone();
        for j in &river_indices {
            female_matrix[[0, *j]] = spwn * 0.5;
        }
        salmon_amounts[0] = female_matrix
            .dot(&(&salmon_amounts[0] - &fishing_female - alpha * &seal_mortality_female));
        // UROS
        let mut male_matrix = leslie_matrix.clone();
        for j in &river_indices {
            male_matrix[[0, *j]] = 0.0;
        }
        salmon_amounts[1] =
            male_matrix.dot(&(&salmon_amounts[1] - &fishing_male - alpha * &seal_mortality_male));
        salmon_amounts[1][[0, 0]] = salmon_amounts[0][[0, 0]];
        // Laske nousulohet vuodelle i:
        for sex in 0..2 {
            for j in &river_indices {
                calculated_salmon_run[i] += salmon_amounts[sex][[*j, 0]];
            }
        }
    }
        let r_2: f64 =  1.0 - r2(&measured_salmon_run, &calculated_salmon_run);//dtw(&measured_salmon_run, &calculated_salmon_run);//1.0 - r2(&measured_salmon_run, &calculated_salmon_run);
        if r_2 < best_r2 {
            best_r2 = r_2;
            best_spawn = spwn;
            best_alpha = alpha;
            best_beta = beta;
            best_gamma = gamma;
            best_p_j = p_j;
            best_p_l = p_l;
            best_calculated_salmon_run = calculated_salmon_run;
        }
                } 
            }
        }
    }
if BY_HAND == true {
    break;
}
}
if BY_HAND == true {
    break;
}
}
    // Tulosta vastaukset
    draw::plot_salmon_runs(
        measured_salmon_run.clone(),
        best_calculated_salmon_run.clone(),
        2009,
        "salmon_runs.png",
    )
    .unwrap();
    println!(
        "Calculated salmon runs (2009-{}): ",
        2009 + best_calculated_salmon_run.len() - 1
    );
    println!("best spawn: {}", best_spawn);
    println!("best alpha: {}", best_alpha);
    println!("best beta: {}", best_beta);
    println!("best gamma: {}", best_gamma);
    println!("best p_j: {}", best_p_j);
    println!("best p_l: {}", best_p_l);
    /*for salmon_run in &best_calculated_salmon_run {
        print!(
            "{} ",
            if salmon_run.floor() < 0.0 {
                0.0
            } else {
                salmon_run.floor()
            }
        );
    }*/
    print!("\n");
        println!(
            "R2: {:.2}",
            best_r2
        );
    
}
// Neliövirhe
fn mse(measured: &Vec<f64>, calculated: &Vec<f64>) -> f64 {
    let mut sum: f64 = 0.0;
    let n: usize = measured.len();
    for i in 0..n {
        sum += (measured[i] - calculated[i]) * (measured[i] - calculated[i]);
    }
    return sum / (n as f64);
}
// Varianssi
fn var(nums: &Vec<f64>) -> f64 {
    let mut ans: f64 = 0.0;
    let n: usize = nums.len();
    let average: f64 = nums.iter().sum::<f64>() / n as f64;
    for i in 0..n {
        ans += (nums[i] - average) * (nums[i] - average);
    }
    return ans / (n as f64);
}
// Selitysaste
fn r2(measured: &Vec<f64>, calculated: &Vec<f64>) -> f64 {
    return 1.0 - mse(measured, calculated) / var(measured);
}
// Laske ikävaiheluokkien jakauman
fn calc_start(total_population: f64, years_to_run: usize, leslie: &Array2<f64>) -> Array2<f64> {
    let mut population: Array2<f64> = Array2::ones((leslie.nrows(), 1));
    for _ in 0..years_to_run {
        population = leslie.dot(&population);
    }
    let total: f64 = (0..leslie.nrows()).map(|k| population[[k, 0]]).sum::<f64>();
    let normalized: Vec<f64> = population
        .iter()
        .map(|&x| x * total_population as f64 / total)
        .collect();
    return Array2::from_shape_vec((leslie.nrows(), 1), normalized).unwrap();
}
// Jaa skalaari johonkin pystyvektoriin ottaen huomiin niissä indekseissä olevat runsausvektorin lohet.
fn distribute_value(
    input: f64,
    indices: &Vec<usize>,
    salmon_amount: &Array2<f64>,
    output: &mut Array2<f64>,
) {
    let total = indices.iter().map(|k| salmon_amount[[*k, 0]]).sum::<f64>();
    for &k in indices {
        let h = input * salmon_amount[[k, 0]] / total;
        output[[k, 0]] = if h <= 0.0 { 0.0 } else { h };
    }
}

fn dtw(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len();
    let m = y.len();

    let mut cost = vec![vec![f64::INFINITY; m + 1]; n + 1];
    cost[0][0] = 0.0;

    for i in 1..=n {
        for j in 1..=m {
            let dist = (x[i - 1] - y[j - 1]).abs();

            cost[i][j] = dist
                + cost[i - 1][j]
                    .min(cost[i][j - 1])
                    .min(cost[i - 1][j - 1]);
        }
    }

    cost[n][m]
}