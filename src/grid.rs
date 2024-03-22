use std::io::{self, Write};


#[derive(Debug, Clone, Copy)]
pub enum Direction {
    Forward,
    Backward,
    Jump,
    Stand,
    Quit
}


#[derive(Debug, Clone, Copy)]
pub struct Player {
    pub position: usize,
    pub action: Direction
}


impl Player {
    pub fn new() -> Player {
        Player {
            position: 0,
            action: Direction::Stand
        }
    }

    pub fn parse_input(&mut self, input: &str) {
        let proc_input = input.trim().to_lowercase();

        match proc_input.as_str() {
            "w" => self.action = Direction::Jump,
            "d" => self.action = Direction::Forward,
            "a" => self.action = Direction::Backward,
            "q" => self.action = Direction::Quit,
            &_ => self.action = Direction::Stand
        }
    }


    pub fn execute_action(&mut self, grid: &Grid) {
        match self.action {
            Direction::Forward => {
                self.position += 1;
                println!("{}", grid.sv_state.get(self.position).expect("end of vec"));
            },
            Direction::Backward => {
                self.position -= 1;
                println!("{}", grid.sv_state.get(self.position).expect("end of vec"));
            },
            Direction::Jump => {
                println!("{}", grid.sv_state.get(self.position).expect("end of vec"));
            },
            Direction::Stand => {
                println!("{}", grid.sv_state.get(self.position).expect("end of vec"));
            },
            Direction::Quit => {
                println!("{}","-".repeat(80));
                println!(r"
 _                  _                
| |                | |               
| |__  _   _  ___  | |__  _   _  ___ 
| '_ \| | | |/ _ \ | '_ \| | | |/ _ \
| |_) | |_| |  __/ | |_) | |_| |  __/
|_.__/ \__, |\___| |_.__/ \__, |\___|
        __/ |              __/ |     
       |___/              |___/      
                         ");
                println!("{}","-".repeat(80));
            }
        }
    }
}


#[derive(Debug)]
pub struct Grid {
    pub bin_id: Vec<i64>,
    pub sv_state: Vec<String>,
    pub llr_to_ref: Vec<f64>
}


pub fn get_input(player: &mut Player) {
    // Prompt
    println!("");
    print!("> ");
    io::stdout().flush().unwrap();

    let mut input_str = String::new();

    io::stdin()
        .read_line(&mut input_str)
        .expect("Failed to read move");
    println!("");

    // Parse and set the action
    player.parse_input(input_str.as_str());

}
