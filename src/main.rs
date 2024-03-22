pub mod grid;
use std::vec;

use crate::grid::Direction;
use crate::grid::Grid;
use crate::grid::Player;
use crate::grid::get_input;

fn main() {

    // initial message
    println!("{}","-".repeat(80));
    println!(r"
    __  ___                 _         ____                             
   /  |/  /___  _________ _(_)____   / __ \__  ______  ____  ___  _____
  / /|_/ / __ \/ ___/ __ `/ / ___/  / /_/ / / / / __ \/ __ \/ _ \/ ___/
 / /  / / /_/ (__  ) /_/ / / /__   / _, _/ /_/ / / / / / / /  __/ /    
/_/  /_/\____/____/\__,_/_/\___/  /_/ |_|\__,_/_/ /_/_/ /_/\___/_/     
                                                                       
             ");
    println!("{}","-".repeat(80));


    // initialisations
    let grid = Grid{
        bin_id: vec![1, 2, 3, 4],
        sv_state: vec!["ref".to_string(), "amp".to_string(), "del".to_string(), "inv".to_string()],
        llr_to_ref: vec![3.2, 1.6, 2.9, 4.1]
    };
    let mut player = Player::new();
    println!("{}", grid.sv_state.get(player.position).expect("end of vec"));

    while player.position <= grid.bin_id.len() {
        get_input(&mut player);
        player.execute_action(&grid);
        
        match player.action {
            Direction::Quit => break,
            _ => continue
        }
    }
}
