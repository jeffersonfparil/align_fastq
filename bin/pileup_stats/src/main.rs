use std::env;
use std::fs::File;
use std::io::{self, prelude::*, BufReader, BufWriter};

trait Stats {
    fn min(&self) -> io::Result<f64>;
    fn max(&self) -> io::Result<f64>;
    fn sum(&self) -> io::Result<f64>;
    fn mean(&self) -> io::Result<f64>;
    fn raise(&self, power:f64) -> io::Result<Vec<f64>>;
    fn var(&self) -> io::Result<f64>;
}

impl Stats for Vec<f64> {
    fn min(&self) -> io::Result<f64> {
        let n = self.len();
        let mut out = self[0];
        if n > 1 {
            for i in 1..n {
                if self[i] < out {
                    out = self[i];
                }
            }
        }
        Ok(out)
    }

    fn max(&self) -> io::Result<f64> {
        let n = self.len();
        let mut out = self[0];
        if n > 1 {
            for i in 1..n {
                if self[i] > out {
                    out = self[i];
                }
            }
        }
        Ok(out)
    }

    fn sum(&self) -> io::Result<f64> {
        let n = self.len();
        let mut out = self[0];
        if n > 1 {
            for i in 1..self.len() {
                out += self[i];
            }
        }
        Ok(out)
    }

    fn mean(&self) -> io::Result<f64> {
        let n = self.len() as f64;
        let sum = self.sum().unwrap();
        Ok(sum / n)
    }

    fn raise(&self, power: f64) -> io::Result<Vec<f64>> {
        let n = self.len();
        let mut out: Vec <f64> = Vec::new();
        for i in 0..n {
            out.push(f64::powf(self[i], power));
        }
        Ok(out)
    }

    fn var(&self) -> io::Result<f64> {
        let n = self.len() as f64;
        // println!("{:?}", n);
        let sum_of_squares = self.raise(2.0).unwrap().sum().unwrap();
        let mean_squared = f64::powf(self.mean().unwrap(), 2.0);
        let out = (sum_of_squares / n) - mean_squared;
        Ok(out)
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let fname = &args[1];
    let fname_out = fname.to_owned() + ".stats";
    // println!("Input: {:?} | Output: {:?}", fname, fname_out);

    let file: File = File::open(fname).expect("Input file not found.");
    let reader: BufReader<File> = BufReader::new(file);

    // let file_out: File = File::create_new(fname_out).expect("Output fle already exists. Please delete it first. Thanks a bunch!");
    let file_out : File = File::create(fname_out).unwrap();
    let mut writer = BufWriter::new(file_out);

    let mut i = 0;
    for line in reader.lines() {
        i += 1;
        let l = &line.unwrap();
        let v: Vec<&str> = l.split("\t").collect();
        // println!("{:?}", v);
        let n = v.len();
        let mut coverage: Vec<f64> = Vec::new();
        for j in (3..n).step_by(3) {
            // println!("{}", j);
            let cov = match v[j].to_string().parse::<f64>() {
                Ok(x) => x,
                Err(_) => return (),
            };
            coverage.push(cov);
            // println!("{}", cov);
        }
        // println!("Minimum: {:?}", coverage.min().unwrap());
        // println!("Maximum: {:?}", coverage.max().unwrap());
        // println!("Sum: {:?}", coverage.sum().unwrap());
        // println!("Mean: {:?}", coverage.mean().unwrap());
        // println!("Variance: {:?}", coverage.var().unwrap());
        let stats = vec![v[0].to_string(),
                                 v[1].to_string(),
                                 coverage.min().unwrap().to_string(),
                                 coverage.max().unwrap().to_string(),
                                 coverage.sum().unwrap().to_string(),
                                 coverage.mean().unwrap().to_string(),
                                 coverage.var().unwrap().to_string()].join("\t") + "\n";
        let out = stats.as_bytes();
        writer.write_all(out).expect("Unable to write for some reason...LOL");

        // if i > 10 {
        //     break
        // }
    }
}
