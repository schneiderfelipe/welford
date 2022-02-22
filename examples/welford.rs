use std::{
    error::Error,
    io::{stdin, stdout, Write},
};
use welford::Welford;

fn main() -> Result<(), Box<dyn Error>> {
    let mut w = Welford::new();

    while let Ok(value) = {
        print!("Please enter a real number (anything else to quit): ");
        stdout().flush()?;

        let mut buffer = String::new();
        stdin().read_line(&mut buffer)?;
        buffer.trim().parse::<f64>()
    } {
        w.push(value);
    }

    println!("mean: {:?}", w.mean());
    println!("variance: {:?}", w.var());

    Ok(())
}
