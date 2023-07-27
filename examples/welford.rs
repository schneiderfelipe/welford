use std::error::Error;
use std::io::stdin;
use std::io::stdout;
use std::io::Write;

use welford::Welford;

fn main() -> Result<(), Box<dyn Error>> {
    let mut w = Welford::<f32>::new();

    while let Ok(value) = {
        print!("Please enter a real number (anything else to quit): ");
        stdout().flush()?;

        let mut buffer = String::new();
        stdin().read_line(&mut buffer)?;
        buffer.trim().parse()
    } {
        w.push(value);
    }

    println!("mean: {:?}", w.mean());
    println!("variance: {:?}", w.var());

    Ok(())
}
